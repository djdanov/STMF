stmf <- function(cntr, frcst = 1) {
  require(stringr)
  #ISL: frcst = 0
  if (!file.exists(paste(".\\STMF\\", cntr, "stmf.csv", sep = ""))) {
    cat(paste(cntr, ": no data\n"))
    return(0)
  }
  checkHMD(cntr)
  cat(cntr)
  cat("\n")
  ages = c(0, 15, 65, 75, 85, 260)
  ascale = c("0-14", "15-64","65-74","75-84","85+", "Total")
  sex = c("b", "m", "f")
  stmf.in = read.csv(paste(".\\STMF\\", cntr, "stmf.csv", sep = ""), stringsAsFactors = F, na.strings = ".")
  stmf.in$Week =as.numeric(stmf.in$Week)
  stmf.in = stmf.in[!is.na(stmf.in$Deaths) & !is.na(stmf.in$Week) & stmf.in$Week <= 53, ]
  stmf.in$Age[stmf.in$Age == "TOT"] = 300
  stmf.in$Age[stmf.in$Age == "UNK"] = 500
  stmf.in$Age = as.numeric(stmf.in$Age)
  stmf.in$AgeInterval[stmf.in$AgeInterval == "+"] = 130
  stmf.in$AgeInterval[is.na(stmf.in$AgeInterval)] = 999
  stmf.in$AgeInterval = as.numeric(stmf.in$AgeInterval)
  stmf.in$Deaths = as.numeric(stmf.in$Deaths)
  
  cat("Year:", unique(stmf.in$Year), "=>", length(unique(stmf.in$Year)), "\n")
  cat("Week:", unique(stmf.in$Week), "=>", length(unique(stmf.in$Week)), "\n")
  cat("Age:", unique(stmf.in$Age), "\n")
  cat("Age interval:", unique(stmf.in$AgeInterval), "\n")
  #cat("Age interval:", unique(stmf.in$AgeInterval), "\n")
  cat("Deaths (non-numeric values):", length(stmf.in$Deaths[is.na(stmf.in$Deaths)]), "\n")
  years = unique(stmf.in$Year)

  d = NULL
  fy = min(years)
  ly = max(years)
  maxAge = 100
  
  Expm = exp.forecast(paste(".\\HMD\\m", cntr,".txt", sep=""), paste(".\\HMD\\", cntr,"mltper_1x1.csv", sep=""), 
                     fy = fy, lastyear = ly, maxAge = maxAge, frcst = frcst)
  Expf = exp.forecast(paste(".\\HMD\\f", cntr,".txt", sep=""), paste(".\\HMD\\", cntr,"fltper_1x1.csv", sep=""), 
                      fy = fy, lastyear = ly, maxAge = maxAge, frcst = frcst) 

  yHMD = Expm$y
  Dxm = Expm$Dx
  yminHMD = as.numeric(str_sub(names(Expm$Ex)[2], -4, -1))

  Dxmagg = exp.agg(Dxm, ages = head(ages, -1))
  Dxf = Expf$Dx
  Dxfagg = exp.agg(Dxf, ages = head(ages, -1))
  Dxb= Expm$Dx + Expf$Dx
  Dxb[, 1] = Dxm[, 1]

  Expm =  exp.agg(Expm$Ex, ages = head(ages, -1))  
  Expf =  exp.agg(Expf$Ex, ages =  head(ages, -1))  
  
  yi = 2
  stmf.in = stmf.in[stmf.in$Age <= 110 | stmf.in$Age == 500, ] # no total
  yH = yminHMD
  for (y in unique(stmf.in$Year)) { #
    if (y < yminHMD) {
      next 
    } 
    else {
      if (y != yH) {browser() }
      yH = yH + 1
    }
    cat("year", y, ":")
    nm = 0
    nf = 0
    nb = 0
    W = max(stmf.in$Week[stmf.in$Year == y])
    for (w in 1 : W) { #unique(stmf.in$Week[stmf.in$Year == y])
      d1 = stmf.in[stmf.in$Year == y & stmf.in$Sex == "m" & stmf.in$Week == w, ]
      if (nrow(d1)) {
        dm = stmf.out(d1, ages, Expm[, yi]/ 52, Dxm[, c(1, yi)], y, yHMD, w, s = "m")
        d = rbind(d, dm)
        splitsex = 0
        nm = nm + 1
      }
      else {
        splitsex = 1
      }
      d1 = stmf.in[stmf.in$Year == y & stmf.in$Sex == "f" & stmf.in$Week == w, ]
      if (nrow(d1)) {
        df = stmf.out(d1, ages, Expf[, yi]/ 52, Dxf[, c(1, yi)], y, yHMD, w, s = "f")
        d = rbind(d, df)
        nf = nf + 1
      }
      d1 = stmf.in[stmf.in$Year == y & stmf.in$Sex == "b" & stmf.in$Week == w, ]
      if (nrow(d1)) {
        nb = nb + 1
        if (splitsex) {
          d1 = stmf.out(d1, ages, (Expf[, yi] + Expm[, yi])/ 52, Dxb[, c(1, yi)], y, yHMD, w, s = "b")
          ds = split.sex(d1, Dxmagg[, yi], Dxfagg[, yi], ages, Expm[, yi]/ 52, Expf[, yi]/ 52)
          d = rbind(d, ds, d1)
        }
        else {
          t = sum(d1$Death) - as.numeric(df[9])- as.numeric(dm[9])
          if (t){
            #dmf = unk.sex(d1, df, dm, Expm[, yi]/ 52, Expf[, yi]/ 52)
            cat("m+f != tot:", y, w,"diff: ", t, "\n")
          }
          else {
            dmf = rbind(dm, df)
          }
          d1 = combine.sex(df, dm, (Expf[, yi] + Expm[, yi])/ 52)
          d = rbind(d, d1)
        }
      }
      else {
        if (!splitsex) {
          d1 = combine.sex(df, dm, (Expf[, yi] + Expm[, yi])/ 52)
          d = rbind(d, d1)
        }
      }
    }
    yi = yi + 1
    cat(nm, nf, nb, "\n")
  }
  d = data.frame(cbind(cntr, d))
  names(d) <- c("CountryCode", "Year", "Week", "Sex", ascale, ascale, "Split", "SplitSex", "Forecast")
  #stmf.sexcheck(d, correct = 1)
  #stmf.graphs(d, paste(".\\Output\\", cntr, sep = ""))
  write.csv(d, paste(".\\Output\\", cntr, "stmfout.csv", sep = ""), row.names = FALSE, quote = FALSE, na=".")
  return(d)  
}

stmf.sexcheck <- function(d, correct = 1) {
  d[, 5 : 16] = lapply(d[, 5 : 16], function(x) as.numeric(as.character(x)))
  for(y in unique(d$Year)) {
    for(w in unique(d$Week[d$Year == y])) {
      df = d[d$Year == y & d$Week == w & d$Sex == "f", 5 : 10]
      dm = d[d$Year == y & d$Week == w & d$Sex == "m", 5 : 10]
      db = d[d$Year == y & d$Week == w & d$Sex == "b", 5 : 10]
      dif = df + dm - db
      if (length(dif[dif > 0])) {
        cat("Error in sex balance (m + f) > b: ", y, w, "\n")
      }
      if (length(dif[dif < 0])) {#
        if (correct) {
          dm1 = dm
          df1 = df
          dm1[dm == 0] = 0.49
          df1[df == 0] = 0.49
          rm = dm1 / (df1 + dm1)
          rf = 1 - rm
          d[d$Year == y & d$Week == w & d$Sex == "m", 5 : 10] = dm + rm * dif
          d[d$Year == y & d$Week == w & d$Sex == "f", 5 : 10] = df + rf * dif
          d[d$Year == y & d$Week == w & d$Sex == "m", 11 : 16] = d[d$Year == y & d$Week == w & d$Sex == "m", 11 : 16] * (dm + rm * dif) / dm  
          d[d$Year == y & d$Week == w & d$Sex == "f", 11 : 16] = d[d$Year == y & d$Week == w & d$Sex == "f", 11 : 16] * (df + rf * dif) / df  
          cat("Warning: (m + f) < b", y, w, "corrected\n")
        }
        else {
          cat("Error in sex balance: (m + f) < b", y, w, "\n")
        }
      }
    }
  }
}

combine.sex <- function(df, dm, Ex) {
  df[, 4 : 17] = lapply(df[, 4 : 17], function(x) as.numeric(as.character(x)))
  dm[, 4 : 17] = lapply(dm[, 4 : 17], function(x) as.numeric(as.character(x)))
  db = df
  db[1, 3] = "b"
  db[1, 4 : 9] = df[1, 4 : 9] + dm[1, 4 : 9]
  db[1, 10 : 15] = db[1, 4 : 9] / c(Ex, sum(Ex))
  return(db)
}


split.sex <- function(d1, Dxm, Dxf, ages, expsm, expsf) {
  d = data.frame(rbind(d1, d1))
  d[, 4 : 17] = lapply(d[, 4 : 17], function(x) as.numeric(as.character(x)))
  d[, 3] = as.character(d[, 3])
  r = Dxm / (Dxf + Dxm)
  d[1, 4 : 8] = r * d[1, 4 : 8]
  d[1, 9] = sum(d[1, 4 : 8])
  d[1, 10 : 14] = d[1, 4 : 8] / expsm
  d[1, 15] = d[1, 9] / sum(expsm)
  d[1, 3] = "m"
  d[2, 3] = "f"
  #females
  d[2, 4 : 9] = d[2, 4 : 9] - d[1, 4 : 9]
  d[2, 10 : 14] = d[2, 4 : 8] / expsf #rates
  d[2, 15] = d[2, 9] / sum(expsf) #CDR
  d[, 17] = 1 #split sex
  return(d)
}

stmf.split <- function(dta, mx, ages) {
  split = 0
  dta = dta[order(dta$Age), ]
  dta1 = NULL
  ai = 2
  for (i in 1 : dim(dta)[1]) {
    if(dta$Age[i] + dta$AgeInterval[i] == ages[ai]) {
        ai = ai + 1
    }
    while(dta$Age[i] + dta$AgeInterval[i] > ages[ai]) {
        da = split1(dta[i, ], ages[ai], mx)
        dta1 = rbind(dta1, da[1, ])
        dta[i, ] = da[2, ]
        ai = ai + 1
        split = 1
    }
    dta1 = rbind(dta1, dta[i, ])
  }

  st = data.frame(cbind(dta1[c("Age", "AgeInterval", "Deaths")], split))
  names(st) = c("Age", "AgeInterval", "Deaths", "Split")
  st = st[order(st$Age), ]
  return(st)
}


split1 <- function(dta, a, d){
  d = d[d[, 1] >= dta$Age & d[, 1] < (dta$Age + dta$AgeInterval),]
  dta1 = dta
  r = sum(d[d[, 1] < a, 2]) / sum(d[, 2])
  dta$Deaths = r * dta$Deaths
  dta1$Deaths = dta1$Deaths - dta$Deaths
  dta$AgeInterval = a - dta$Age
  dta1$Age = a
  dta1$AgeInterval = dta$Age + dta1$AgeInterval - a
  return(rbind(dta, dta1))
}

stmf.agg <- function(dta, ages) {
  st = vector(mode = "numeric", length = length(ages) - 1)
  for (i in 1 : (length(ages) - 1)) {
    st[i] = sum(dta$Deaths[dta$Age >= ages[i] & (dta$Age + dta$AgeInterval) <= ages[i+1]])
  }
  return(st)
}

stmf.unk <- function(dta) {
  unk = dta$Deaths[dta$Age == 500]
  if (length(unk)) {
    dta = dta[dta$Age < 500, ]
    dta$Deaths = dta$Deaths + dta$Deaths / sum(dta$Deaths) * unk
  }
  return(dta)
}

stmf.out <- function(d1, ages, exps1, mx, y, yHMD, w, s) {
 frcst = if (y > yHMD) 1 else 0
 splitsex = 0
 d1 = stmf.unk(d1)
 d1 = stmf.split(d1, mx, ages)
 split = d1$Split[1]
 d1 = stmf.agg(d1[c("Age", "AgeInterval", "Deaths")], ages)
 r1 = d1 / exps1
 d = c(y, w, s, d1, sum(d1), r1, sum(d1) / sum(exps1), split, splitsex, frcst)
 return(data.frame(t(d), stringsAsFactors = F))
}

stmf.graphs <- function(dta, country) {
  dta = dta[-1]
  require(lattice)
  require(gplots)
  pdf(file = paste(country, "checks.pdf", sep = ""), paper='a4r', width = 11, height = 8, onefile = T)
  Sex = c("both sexes", "males", "females")
  agegr = c("Ages 0-14, ", "Ages 15-64, ", "Ages 65-74, ", "Ages 75-84, ", "Ages 85+, ", "All ages, ")
  sex = c("b", "m", "f")
  for (s in 1 : 3) {
    for (i in 4 : 9) {
      stmf.graphs1(dta, sex = sex[s], title = paste(agegr[i - 3], Sex[s], sep = ""), 
                   c = i, ytitle = "death counts")
      stmf.graphs1(dta, sex = sex[s], title = paste(agegr[i - 3], Sex[s], sep = ""), 
                   c = i + 6, ytitle = "death rates")
    }
  }
  for (s in 1 : 3) {
    for (i in 4 : 9) {
      stmf.graphs2(dta, sex = sex[s], title = paste(agegr[i - 3], Sex[s], sep = ""), 
                   c = i, ytitle = "death counts")
      stmf.graphs2(dta, sex = sex[s], title = paste(agegr[i - 3], Sex[s], sep = ""), 
                   c = i + 6, ytitle = "death rates")
    }
  }
  Year = as.numeric(as.character(dta$Year))
  for (s in 1 : 3) {
    for (yi in min(Year) : max(Year)) {
      stmf.graphs3(dta, sex = sex[s], title = paste(Sex[s], ", year ", yi, sep = ""), 
                   c = 4:9, y = yi, ytitle = "death counts")
      stmf.graphs3(dta, sex = sex[s], title = paste(Sex[s],  ", year ", yi, sep = ""), 
                   c = 10:15, y = yi,  ytitle = "death rates")
    }
  }
  dev.off()
}

stmf.graphs1 <- function(dta, sex = "b", title = "Deaths counts", c = 9, ytitle = "deaths") {
  dta = dta[dta$Sex == sex, ]
  dta$Week = as.numeric(as.character(dta$Week))
  dta = data.frame(dta[, c(1, 2, c)])
  names(dta) = c("Year", "Week", "Deaths")
  dta$Deaths = as.numeric(as.character(dta$Deaths))
  years = as.numeric(as.character(dta$Year))
  ncol = length(unique(dta$Year))
  cls = rainbow(ncol, rev = T)
  LWD = vector(mode = "numeric", ncol) + 1.5
  LWD[ncol] = 2
  print( xyplot(Deaths ~ Week, dta, group = Year, type = "l", lty = 1, lwd = LWD,  col = cls,
         grid = T, xlim = c(1, 53), scales=list(x = list(at = c(1,seq(0, 51, by=5)),
                                                         labels = c(1,seq(0, 51, by=5)))),
         key=list(text=list(unique(dta$Year)), space='right', cex=1.2,
                  points=list(pch = 15, col = cls, cex=1.2), columns=1),
         main = title,
         ylab = ytitle, xlab = "week"
        ))

} 

stmf.graphs2 <- function(dta, sex = "b", title = "Deaths counts", c = 11, ytitle = "deaths") {
  
  dta = dta[dta$Sex == sex, ]
  dta = data.frame(dta[, c(1, 2, c)])
  names(dta) = c("Year", "Week", "Deaths")

  dta$Week = as.numeric(as.character(dta$Week))
  dta$Deaths = as.numeric(as.character(dta$Deaths))
  dta$Year = as.numeric(as.character(dta$Year))
  #dta = dta[dta$Year >= 2017 ,]
  ncol = length(unique(dta$Year))
  cls = rainbow(ncol, rev = T)
  LWD = vector(mode = "numeric", ncol) + 1.5
  LWD[ncol] = 2
  
  boxplot2(dta$Deaths ~ dta$Year, col=cls, ylab = ytitle, xlab = "year", varwidth=T, main = title, top = T)
  grid(nx=NA, ny=NULL) #grid over boxplot
  par(new=TRUE)
  boxplot2(dta$Deaths ~ dta$Year, col=cls, ylab = ytitle, xlab = "year", varwidth=T, main = title, top = T, add = TRUE)
} 



stmf.graphs3 <- function(dta, sex = "b", title = "Deaths counts", c = 10:15, y = 2019, ytitle = "deaths") {
  
  d = dta[dta$Sex == sex & dta$Year == y, ]
  d = data.frame(d[, c(2, c)])
  dta = NULL
  age = c(0, 15, 65, 75, 85, 100)
  d$Week = as.numeric(as.character(d$Week))
  ascale = c("0-14", "15-64","65-74","75-84","85+", "Total")
  for (i in 2 : dim(d)[2]) {
    dta = rbind(dta, cbind(age[i - 1], d[, 1], as.numeric(as.character(d[, i]))))
  }
  dta = data.frame(dta)
  names(dta) = c("Age", "Week", "Deaths")
  
  ncol = length(c)
  cls = rainbow(ncol, rev = T)
  LWD = vector(mode = "numeric", ncol) + 1.5
  
  boxplot2(dta$Deaths ~ dta$Age, col=cls, ylab = ytitle, xlab = "age group", main = title, names = ascale, varwidth=T, top = T)
  grid(nx=NA, ny=NULL) #grid over boxplot
  par(new=TRUE)
  boxplot2(dta$Deaths ~ dta$Age, col=cls, ylab = ytitle, xlab = "age group", main = title, names = ascale, varwidth=T, top = T)
} 

checkHMD <- function(cntr) {
  #require(RCurl)
  hmd = "https://hmd:lexis1x1@www.mortality.org/RunResults/"
#LDB females
  fname = paste("f", cntr, ".txt", sep = "")
  fname_local = paste("./HMD/", fname, sep = "")
  if (!file.exists(fname_local)) {
    cat(fname)
    cat(': downloading ... ')
    try(download.file(paste(hmd, cntr, "/LexisDB/", fname, sep=""), fname_local))
    cat('Done\n\r')
  } 
#LDB males
  fname = paste("m", cntr, ".txt", sep = "")
  fname_local = paste("./HMD/", fname, sep = "")
  if (!file.exists(fname_local)) {
    cat(fname)
    cat(': downloading ... ')
    try(download.file(paste(hmd, cntr, "/LexisDB/", fname, sep=""), fname_local))
    cat('Done\n\r')
  } 
#LT females
  fname = paste(cntr, "fltper_1x1.csv", sep = "")
  fname_local = paste("./HMD/", fname, sep = "")
  if (!file.exists(fname_local)) {
    cat(fname)
    cat(': downloading ... ')
    try(download.file(paste(hmd, cntr, "/STATS/fltper_1x1.csv", sep=""), fname_local))
    cat('Done\n\r')
  } 
#LT males  
  fname = paste(cntr, "mltper_1x1.csv", sep = "")
  fname_local = paste("./HMD/", fname, sep = "")
  if (!file.exists(fname_local)) {
    cat(fname)
    cat(': downloading ... ')
    try(download.file(paste(hmd, cntr, "/STATS/mltper_1x1.csv", sep=""), fname_local))
    cat('Done\n\r')
  } 
  return(1)
}

exp.agg <- function(exps, ages = c(0, 15, 65, 75, 85)) {
  #aggregates pop exps by age 
  ages = c(ages, 130)
  ne = dim(exps)[2]
  Ex = matrix(NA, length(ages) - 1, ne - 1)
  for (i in 1 : (length(ages) - 1)) {
    Ex[i, ] = colSums(exps[exps[, 1] >= ages[i] & exps[, 1] < ages[i+1], 2 : ne])
  }
  Ex = data.frame(cbind(head(ages, -1), Ex))
  return(Ex)
}
