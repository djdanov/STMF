exp.forecast <- function(fname, fnameLT, fy = 2010, lastyear = 2020, maxAge = 100, frcst = 1) {
# pop exposures  
  #fy = min(fy, 2005)
  #fy = 2005

  mx = mx_LC(fname, maxAge=100, fy = 2005, lastyear = lastyear, frcst)
  a = read.LDB(fname, fy = 2005, maxAge = 100, lastyear = lastyear)
  p = a$pop
  y = a$y
  rl = a$rl
  
  obs = readExp(fnameLT, fy = fy, ly = lastyear, maxAge = maxAge)
  if (!is.null(mx)) { 
    Ex = matrix(NA, maxAge + 1, dim(mx)[2] - 1)
    Dx = Ex
    for (i in 1 : dim(Ex)[2]) {
      d1 = mx[, i + 1] * p
      p1 = p - d1
      p1[maxAge] = p1[maxAge] + p1[maxAge + 1]
      p1 = c(a$b[i] - a$imr[i] * a$b[i], p1[1 : maxAge])
      Ex[, i] = (p + p1) / 2
      nd = length(d1)
      d0u = d1[1]
      d1[1 : (nd - 1)] = d1[1 : (nd - 1)] / 2
      d1[1] = d0u * rl  
      d1[nd] = d1[nd] + d1[nd - 1]
      d1[2 : (nd - 1)] = d1[2 : (nd - 1)] + d1[1 : (nd - 2)] 
      d1[1] = d0u * (1 - rl) + a$imr[i] * a$b[i]
      Dx[, i] = d1
      p = p1
    }
    Ex = cbind(obs$Ex, Ex)
    Dx = cbind(obs$Dx, Dx)
  }
  else {
    Ex = obs$Ex
    Dx = obs$Dx
  }
  #write.csv(Ex,file="fout.csv",row.names=FALSE,quote=FALSE, na=".")
  return(structure(list(Ex = Ex, Dx = Dx, y = y)))
}

#Lee-Carter forecast

mx_LC<-function(fname, maxAge=100, fy = 2000, lastyear = 2020, frcst = 1) {
  # extrapolate mx using LC model  

  dta = read.LDB2LC(fname, fy, maxAge)
  maxinterp=lastyear - max(dta$year)
  if (frcst & maxinterp > 0) {
    dta.lca=lca(dta,years = fy : max(dta$year), adjust="e0", interpolate=TRUE)
    dta.fcast=forecast(dta.lca, maxinterp)
    rates=data.frame(cbind(dta$age, dta.fcast$rate$rates))
    names(rates)<-c('age', dta.fcast$year)
  }
  else {
    if (maxinterp <= 0) {
      rates = NULL
    }
    else {
      r = dta$rate$rates[, ncol(dta$rate$rates)]
      rates = NULL
      for (i in 1 : maxinterp) {
        rates = cbind(rates, r)
      }
      rates = data.frame(cbind(dta$age, rates))
    }
  }
  #  write.csv(rates,file=fout,row.names=FALSE,quote=FALSE, na=".")
  return(rates)
}

#reads Lexis file (HMD)
read.LDB2LC<-function(fname, fy=2000, maxAge=100) {
  #reads from Lexis files, returns LC data
  require(demography)
  maxA=130
  
  dta=read.csv(fname, header=FALSE, stringsAsFactors = F, na.strings = ".")
  maxY=max(dta[,1] - 1)
  dta=dta[dta[,1] <= maxY & dta[, 1] >= fy,]
  
  pop=dta[dta[,3] == 2, ]
  pop=pop[,c(1,2,5)]
  years=fy : maxY
  pop=data.frame(pop)
  names(pop)<-c("Year","Age","Pop")
  pop=reshape(pop, v.names="Pop", idvar="Age",timevar="Year", direction="wide")
  
  pop[maxAge+1,]=colSums(pop[(maxAge+1):131,])
  pop=pop[1 : (maxAge+1), ]
  pop=pop[, 2 : dim(pop)[2]]
  
  d1=dta[dta[,3] == 1 & dta[, 2] > 0,]
  d1=d1[,c(1,2,6)]
  d1=data.frame(d1)
  names(d1)<-c("Year","Age","D")
  #d1=d1[d1$Year>=years[1] & d1$Year<maxY,]
  d1=reshape(d1, v.names="D", idvar="Age",timevar="Year", direction="wide")
  d1[maxAge+1,]=colSums(d1[(maxAge+1):130,])
  d1=d1[1 : (maxAge+1), ]
  d1=d1[, 2 : dim(d1)[2]]
  
  d2=dta[dta[,3] == 2, ]
  d2=d2[,c(1,2,6)]
  d2=data.frame(d2)
  names(d2)<-c("Year","Age","D")
  #d2=d2[d2$Year>years[1],]
  d2=reshape(d2, v.names="D", idvar="Age",timevar="Year", direction="wide")
  d2[maxAge+1,]=colSums(d2[(maxAge+1):131,])
  d2=d2[1:(maxAge+1),]
  d2=d2[,2:dim(d2)[2]]
  
  # p1=pop+(d1-d2)/3
  pop[pop == 0] = 0.01
  # p1[1,]=pop[1,];
  rates=(d1+d2)/pop
  
  ages=0:maxAge
  return(demogdata(rates, pop, ages, years, 'mortality', 'lc', 'rates', 0))
}

#reads pop exposures from the HMD LT (csv)

readExp<-function(fname, fy=2000, ly = ly, maxAge=100) {
  # reads pop exposures from LT (.csv)
  dta=read.csv(fname, header=FALSE, stringsAsFactors = F, na.strings = c("",".","NA"), colClasses = c("numeric"))
  dta = dta[dta[, 1] >= fy & dta[, 1] <= ly, ]
  maxY=max(dta[,1])
  
  pop=data.frame(dta[,c(1,3,4)])
  names(pop)<-c("Year","Age","Pop")
  pop=reshape(pop, v.names="Pop", idvar="Age",timevar="Year", direction="wide")
  
  pop[maxAge + 1, ] = colSums(pop[(maxAge + 1) : 111, ])
  pop[maxAge + 1, 1] = 100
  pop = pop[1 : (maxAge + 1), ]
  
  d=data.frame(dta[,c(1,3,5)])
  names(d)<-c("Year","Age","D")
  d=reshape(d, v.names="D", idvar="Age",timevar="Year", direction="wide")
  
  d[maxAge + 1, ] = colSums(d[(maxAge + 1) : 111, ])
  d[maxAge + 1, 1] = 100
  d = d[1 : (maxAge + 1), ]
  
  return(structure(list(Ex = pop, Dx = d)))
  #return(pop)
}

# reads data from the HMD Lexis file
# returns structure: last year pop, IMR, births 
read.LDB<-function(fname, fy = 2000, maxAge = 100, lastyear = 2020) {
  require(Hmisc)
  maxA = 130
  
  dta = read.csv(fname, header = FALSE, stringsAsFactors = F, na.strings = ".")
  maxY = max(dta[, 1])
  dta = dta[dta[, 1] >= fy,]
  
  d0u = sum(dta[dta[, 2] == 0 & dta[, 3] == 2, 6])
  d1l = sum(dta[dta[, 2] == 1 & dta[, 3] == 1, 6])

  rl = d1l/(d1l + d0u)
  
  pop = dta[dta[, 3] == 2 & dta[, 1] == maxY, 5]
  pop[maxAge + 1] = sum(pop[(maxAge + 1) : 131])
  pop = pop[1 : (maxAge + 1)]
  
  p = dta[dta[, 1] < maxY & dta[, 2] == 0 & dta[, 3] == 1, 5]
  d = dta[dta[, 1] < maxY & dta[, 2] == 0 & dta[, 3] == 1, 6]
  years = fy : (maxY - 1)
  imrs = d / p
  #mod_lm = lm(imrs ~ years)
  imr = predict(lm(imrs ~ years), data.frame(years = maxY : lastyear))
  b = predict(lm(p ~ years), data.frame(years = maxY : lastyear))
  return(structure(list(pop = pop, imr = imr, b = b, y = maxY - 1, rl = rl)))
}
