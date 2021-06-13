# 12 m mice

file.list <- list.files('../data', full.names = TRUE)

# 3088 - S1 - 16 Oct	2 copies check mouse numbers
# 3207 - S1 - 09 Mar	no trials
# check for data from 3084 and 3104


dat <- lapply(file.list, function(x) {
	a <- read.csv(x, nrow = 2)
	if(ncol(a) == 25) {
		a <- read.csv(x, comment.char = '#')
	} else {
		a <- read.csv(x, skip = 6, comment.char = '#')
	}
	if(is.numeric(a$Podour)){a$Podour <- as.character(a$Podour)}
	a
})
length(dat) == length(file.list)

dat <- do.call(rbind, dat)
dat <- dat[-which(is.na(dat$animal)),]
summary(dat)
str(dat)
dat$animal <- as.factor(dat$animal)
levels(dat$animal)
# write.csv(data.frame(animal = levels(dat$animal), geno = 'wt', sex = 'f', dob = 'yyyy-mm-dd', strain = '5x'), '12m geno.csv', row.names = FALSE) 


## fix errors on data
dat$StudyName[dat$StudyName == 'REV 1'] <- 'REV1'
dat$StudyName[dat$StudyName == 'OP0'] <- 'INTRO'
dat$Podour[dat$Podour %in% c('1PPM EA ', '1ppm EEA', '1PPM EA', '1 PPM EA')] <- '1ppm EA'
dat$Podour[dat$Podour %in% c('0.1 PPM EA', '0.1PPM EA ', '0.1PPM EA', '0.1PM EA')] <- '0.1ppm EA'
dat$Podour[dat$Podour %in% c('0.01PPM EA', '0.01PPM EA ')] <- '0.01ppm EA'
dat$Podour[dat$Podour %in% c('0.001PPM EA', '0.001', '0.OO1PPM EA')] <- '0.001ppm EA'
dat$Podour[dat$Podour %in% c('0.0001PPM EA')] <- '0.0001ppm EA'
dat$Podour[dat$Podour %in% c('0.00001PPM EA')] <- '0.00001ppm EA'
dat$Nodour[dat$Nodour == 'BANK'] <- 'BLANK'
dat$Nodour[dat$Nodour == 'SEARMINT'] <- 'SPEARMINT'
dat$Podour[dat$Podour == 'SEARMINT'] <- 'SPEARMINT'
dat$Podour[dat$Podour == 'ATCHOULI'] <- 'PATCHOULI'

dat <- dat[-which(dat$animal == '2988' & dat$StudyName == 'SEN2'),]

dat$Podour <- as.factor(as.character(dat$Podour))
dat$Nodour <- as.factor(as.character(dat$Nodour))

geno <- read.csv('12m geno.csv')
geno$dob <- as.Date(geno$dob)
dat <- merge(dat, geno)



concList <- 10^-(c(0:5,5))

datout <- data.frame(mouse=1, stage=1, block=1, hit=1, ptrials=1, cr=1, ntrials=1, geno=1, sex=1, strain=1, incomplete = FALSE)
linenum <- 1
# datout$nblock <- 0
for(mouse in levels(factor(dat$animal))) {
  mmdat <- dat[dat$animal == mouse,]
  
  for(stage in levels(factor(mmdat$StudyName))[which(substr(levels(factor(mmdat$StudyName)), 1, 3 ) == 'SEN')] ) {  ## had used: levels(factor(mmdat$StudyName))[-7:-9]
    ldat <- dat[dat$animal == mouse & dat$StudyName == stage & dat$SS == FALSE,]
ldat <- ldat[order(ldat$date,ldat$counter),]
daycount <- 1
addcount <- 0
ldat$nblock <- 0
incomplete <- FALSE
for(day in levels(factor(ldat$date)) ) {
  if(nrow(ldat[ldat$date == day,]) > 0 ) {
    ldat$ncounter[ldat$date == day] <- 1:nrow(ldat[ldat$date == day,]) + addcount
    addcount <- max(ldat$ncounter[ldat$date == day]) + 1
  }
  
}
ldat <- ldat[ldat$ncounter <= 100,]
for(newblock in 1:5) {
  startblock <- ( 20 * (newblock - 1) ) +1
  endblock <- ( 20 * newblock )
  if(startblock > nrow(ldat)) {
  	incomplete <- TRUE
  } else if( endblock > nrow(ldat) ) {
    ldat$nblock[startblock:nrow(ldat)] <- newblock
  } else {
    ldat$nblock[startblock:endblock] <- newblock
  }
}
for(block in levels(factor(ldat$nblock))) {
  d <- ldat[ldat$nblock == block,]
  hit <- length(d$animal[d$TrialType == 'P' & d$correct == TRUE])
  p.trials <- length(d$animal[d$TrialType == 'P'])
  cr <- length(d$animal[d$TrialType == 'N' & d$correct == TRUE])
  n.trials <- length(d$animal[d$TrialType == 'N'])
  
  datout[linenum,]$mouse <- mouse
  datout[linenum,]$stage <- stage
  datout[linenum,]$block <- as.numeric(block)
  datout[linenum,]$hit <- hit
  datout[linenum,]$ptrials <- p.trials
  datout[linenum,]$cr <- cr
  datout[linenum,]$ntrials <- n.trials
  datout[linenum,]$geno <- levels(factor(d$geno))
  datout[linenum,]$sex <- levels(factor(d$sex))
  datout[linenum,]$strain <- levels(factor(d$strain))
  datout[linenum,]$incomplete <- incomplete
  linenum <- linenum + 1
  
}
  }
}

summary(datout)
str(datout)
ndat <- datout[datout$incomplete == FALSE,]
## Create aggregate to get block hit, FA, miss, CR
ndat$FA <- ndat$ntrials - ndat$cr
ndat$miss <- ndat$ptrials - ndat$hit

ndat$hit <- ndat$hit / ndat$ptrials
ndat$cr <- ndat$cr / ndat$ntrials
ndat$FA <- ndat$FA / ndat$ntrials
ndat$miss <- ndat$miss / ndat$ptrials

# For hits and FA replace 0 with (.5 / n) and 1 with (n -.5) / n ; where n is number of trials


ndat$hit[ndat$hit == 1] <- (ndat$ptrials[ndat$hit == 1] - .5) / ndat$ptrials[ndat$hit == 1]
ndat$hit[ndat$hit == 0] <- (.5 / ndat$ptrials[ndat$hit == 0]) 
ndat$FA[ndat$FA == 1] <- (ndat$ntrials[ndat$FA == 1] - .5) / ndat$ntrials[ndat$FA == 1]
ndat$FA[ndat$FA == 0] <- (.5 / ndat$ntrials[ndat$FA == 0])

# calculate d' (sensitivity index) and threshold
ndat$delta <- qnorm(ndat$hit) - qnorm(ndat$FA)
ndat$crit <- -.5 * (qnorm(ndat$hit) + qnorm(ndat$FA))










########
n <- 100000
# 3x
d3 <- aggregate(delta ~ sex * geno, FUN = mean, data = ndat[ndat$strain == '3x',])
d3$ci <- aggregate(delta ~ sex * geno
	, FUN = function(x) {
		tmp <- replicate(n, mean(sample(x, length(x), replace = T)))
		tmp <- sort(tmp)
		tmp[c(.025*n, .975*n)]
	}
	, data = ndat[ndat$strain == '3x',]
)[,3]

c3 <- aggregate(crit ~ sex * geno, FUN = mean, data = ndat[ndat$strain == '3x',])
c3$ci <- aggregate(crit ~ sex * geno
	, FUN = function(x) {
		tmp <- replicate(n, mean(sample(x, length(x), replace = T)))
		tmp <- sort(tmp)
		tmp[c(.025*n, .975*n)]
	}
	, data = ndat[ndat$strain == '3x',]
)[,3]

# 5x
d5 <- aggregate(delta ~ sex * geno, FUN = mean, data = ndat[ndat$strain == '5x',])
d5$ci <- aggregate(delta ~ sex * geno
	, FUN = function(x) {
		tmp <- replicate(n, mean(sample(x, length(x), replace = T)))
		tmp <- sort(tmp)
		tmp[c(.025*n, .975*n)]
	}
	, data = ndat[ndat$strain == '5x',]
)[,3]


c5 <- aggregate(crit ~ sex * geno, FUN = mean, data = ndat[ndat$strain == '5x',])
c5$ci <- aggregate(crit ~ sex * geno
	, FUN = function(x) {
		tmp <- replicate(n, mean(sample(x, length(x), replace = T)))
		tmp <- sort(tmp)
		tmp[c(.025*n, .975*n)]
	}
	, data = ndat[ndat$strain == '5x',]
)[,3]


m3d <- aov(delta ~ sex * geno, data = ndat[ndat$strain == '3x',])
summary(m3d)

m3c <- aov(crit ~ sex * geno, data = ndat[ndat$strain == '3x',])
summary(m3c)

## only males
# m3dm <- aov(delta ~ geno, data = ndat[ndat$strain == '3x' & ndat$sex == 'm',])
# summary(m3dm)

# m3cm <- aov(crit ~ geno, data = ndat[ndat$strain == '3x' & ndat$sex == 'm',])
# summary(m3cm)


m5d <- aov(delta ~ sex * geno, data = ndat[ndat$strain == '5x',])
summary(m5d)

m5c <- aov(crit ~ sex * geno, data = ndat[ndat$strain == '5x',])
summary(m5c)

#####################
## Plots
#####################

save.plots <- F

lnwd <- 2.5 # line width
cmult <- 2 # cex multiplier
## 3x delta
if(save.plots){
	pdf(file = 'sdt plots/3xDelta.pdf'
		, width = 7
		, height = 6
	)
}

par(mar = c(5,5,4,2) + .01
	# , cex = cmult
)
plot(NULL
	, xlim = c(0.75, 2.25)
	, ylim = c(1, 3)
	, ylab = substitute(paste("Sensitivity index (", italic("d'"), ")"  ))
	, xaxt = 'n'
	, bty = 'l'
	, xlab = ''
	, cex.axis = cmult
	, cex.lab = cmult
)
axis(1
	, at = c(1, 2)
	, labels = c('WT', 'TG')
	, cex.axis = 2
)
d3$x <- c(2,2,1,1)
lines(delta ~ x
	, data = d3[d3$sex == 'f',]
	, col = 'red'
	, lwd = lnwd
)
lines(delta ~ x
	, data = d3[d3$sex == 'm',]
	, col = 'blue'
	, lwd = lnwd
)
arrows(x0 = d3$x
	, y0 = d3$ci[,1]
	, x1 = d3$x
	, y1 = d3$ci[,2]
	, length = 1/8
	, angle = 90
	, code = 3
	, col = ifelse(d3$sex == 'm', 'blue', 'red')
	, lwd = lnwd
)
legend('topright'
	, legend = c('3x Females', '3x Males')
	, col = c('red', 'blue')
	, lty = 1
	, bty = 'n'
	, lwd = lnwd
	, cex = cmult
)
if(save.plots){
	dev.off()
}

## 3x crit

if(save.plots){
	pdf(file = 'sdt plots/3xCrit.pdf'
		, width = 7
		, height = 6
	)
}
par(mar = c(5,5,4,2) + .01
	# , cex = cmult
)

plot(NULL, xlim = c(0.75, 2.25)
	, ylim = c(-.8, 0)
	, ylab = "Criterion"
	, xaxt = 'n'
	, bty = 'l'
	, xlab = ''
	, cex.axis = cmult
	, cex.lab = cmult
)
axis(1
	, at = c(1, 2)
	, labels = c('WT', 'TG')
	, cex.axis = 2
)
c3$x <- c(2,2,1,1)
lines(crit ~ x
	, data = c3[c3$sex == 'f',]
	, col = 'red'
	, lwd = lnwd
)
lines(crit ~ x
	, data = c3[c3$sex == 'm',]
	, col = 'blue'
	, lwd = lnwd
)
arrows(x0 = c3$x
	, y0 = c3$ci[,1]
	, x1 = c3$x
	, y1 = c3$ci[,2]
	, length = 1/8
	, angle = 90
	, code = 3
	, col = ifelse(d3$sex == 'm', 'blue', 'red')
	, lwd = lnwd
)
legend('topright'
	, legend = c('3x Females', '3x Males')
	, col = c('red', 'blue')
	, lty = 1
	, bty = 'n'
	, lwd = lnwd
	, cex = cmult
)
if(save.plots){
	dev.off()
}



###########
# 5x plots



## 5x delta

if(save.plots){
	pdf(file = 'sdt plots/5xDelta.pdf'
		, width = 7
		, height = 6
	)
}
par(mar = c(5,5,4,2) + .01
	# , cex = cmult
)

plot(NULL
	, xlim = c(0.75, 2.25)
	, ylim = c(1, 3)
	, ylab = substitute(paste("Sensitivity index (", italic("d'"), ")"  ))
	, xaxt = 'n'
	, bty = 'l'
	, xlab=''
	, cex.axis = cmult
	, cex.lab = cmult
)
axis(1
	, at = c(1, 2)
	, labels = c('WT', 'TG')
	, cex.axis = 2
)
d5$x <- c(2,2,1,1)
lines(delta ~ x
	, data = d5[d5$sex == 'f',]
	, col = 'red'
	, lwd = lnwd
)
lines(delta ~ x
	, data = d5[d5$sex == 'm',]
	, col = 'blue'
	, lwd = lnwd
)
arrows(x0 = d5$x
	, y0 = d5$ci[,1]
	, x1 = d5$x
	, y1 = d5$ci[,2]
	, angle = 90
	, code = 3
	, col = ifelse(d3$sex == 'm', 'blue', 'red')
	, lwd = lnwd
)
legend('topright'
	, legend = c('5x Females', '5x Males')
	, col = c('red', 'blue')
	, lty = 1
	, bty = 'n'
	, lwd = lnwd
	, cex = cmult
)
if(save.plots){
	dev.off()
}



## 5x crit

if(save.plots){
	pdf(file = 'sdt plots/5xCrit.pdf'
		, width = 7
		, height = 6
	)
}
par(mar = c(5,5,4,2) + .01
	# , cex = cmult
)

plot(NULL, xlim = c(0.75, 2.25)
	, ylim = c(-.8, .3)
	, ylab = "Criterion"
	, xaxt = 'n'
	, bty = 'l'
	, xlab = ''
	, cex.axis = cmult
	, cex.lab = cmult
)
axis(1
	, at = c(1, 2)
	, labels = c('WT', 'TG')
	, cex.axis = 2
)
c5$x <- c(2,2,1,1)
lines(crit ~ x
	, data = c5[c5$sex == 'f',]
	, col = 'red'
	, lwd = lnwd
)
lines(crit ~ x
	, data = c5[c5$sex == 'm',]
	, col = 'blue'
	, lwd = lnwd
)
arrows(x0 = c5$x
	, y0 = c5$ci[,1]
	, x1 = c5$x
	, y1 = c5$ci[,2]
	, length = 1/8
	, angle = 90
	, code = 3
	, col = ifelse(d3$sex == 'm', 'blue', 'red')
	, lwd = lnwd
)
legend('topright'
	, legend = c('5x Females', '5x Males')
	, col = c('red', 'blue')
	, lty = 1
	, bty = 'n'
	, lwd = lnwd
	, cex = cmult
)
if(save.plots){
	dev.off()
}



##### subject values
# write.csv(agdat, 'sensitivity data.csv')

#######################



m3d <- aov(delta ~ sex * geno, data = ndat[ndat$strain == '3x',])
summary(m3d)

m3c <- aov(crit ~ sex * geno, data = ndat[ndat$strain == '3x',])
summary(m3c)

m5d <- aov(delta ~ sex * geno, data = ndat[ndat$strain == '5x',])
summary(m5d)

m5c <- aov(crit ~ sex * geno, data = ndat[ndat$strain == '5x',])
summary(m5c)

#############
# 3x
#############

## dee prime

m3d <- aov(delta ~ sex * geno, data = ndat[ndat$strain == '3x',])
summary(m3d)

d3x1 <- aggregate(delta ~ conc  * geno, data = ndat[ndat$strain == '3x',], mean)
d3x1$conc2 <- rep(1:7, 2)
plot(delta ~ conc2, data = d3x1[d3x1$geno == 'wt',], col = 'blue', ty = 'l', ylim = c(0, 6) )
lines(delta ~ conc2, data = d3x1[d3x1$geno == 'tg',], col = 'red')

d3x2 <- aggregate(delta ~ conc  * sex, data = ndat[ndat$strain == '3x',], mean)
d3x2$conc2 <- rep(1:7, 2)
plot(delta ~ conc2, data = d3x2[d3x2$sex == 'm',], col = 'blue', ty = 'l', ylim = c(0, 6) )
lines(delta ~ conc2, data = d3x2[d3x2$sex == 'f',], col = 'red')


## crit

m3c <- aov(crit ~ sex * geno, data = ndat[ndat$strain == '3x',])
summary(m3c)

c3x1 <- aggregate(crit ~ conc  * geno, data = ndat[ndat$strain == '3x',], mean)
c3x1$conc2 <- rep(1:7, 2)
plot(crit ~ conc2, data = c3x1[c3x1$geno == 'wt',], col = 'blue', ty = 'l', ylim = c(-2.5, 0) )
lines(crit ~ conc2, data = c3x1[c3x1$geno == 'tg',], col = 'red')



####################
## SDT image
####################

## 3x no sex effect
d3agg <- aggregate(delta ~ geno, FUN = mean, data = ndat[ndat$strain == '3x',])
c3agg <- aggregate(crit ~ geno, FUN = mean, data = ndat[ndat$strain == '3x',])
lnwd <- 2
# 3xTG

if(save.plots){
	pdf(file = 'sdt plots/3xTGsdt.pdf'
		, width = 7
		, height = 4.5
	)
}
par(mar = c(4,2,2,2) + .01
	, cex = cmult
)

plot(NULL
	, xlim = c(-3, 4.7)
	, ylim = c(0, .52)
	, ylab = ''
	, xlab = ''
	, main = '3xTG'
	, bty = 'n'
	, yaxt = 'n'
)
abline(h=0)

curve(dnorm(x)	
	, -5
	, 5
	, n = 5000
	, add = T
	, lwd = lnwd
)
curve(dnorm(x
		, mean = d3agg$delta[d3agg$geno == 'tg']	
	)
	, -5
	, 5
	, add = T
	, lty = 'dashed'
	, lwd = lnwd
)
segments(x0 = d3agg$delta[d3agg$geno == 'tg'] / 2 + c3agg$crit[d3agg$geno == 'tg']
	, y0 = 0
	, x1 = d3agg$delta[d3agg$geno == 'tg'] / 2 + c3agg$crit[d3agg$geno == 'tg']
	, y1 = dnorm(d3agg$delta[d3agg$geno == 'tg'] / 2 + c3agg$crit[d3agg$geno == 'tg'])
	, lwd = lnwd
)
arrows(x0 = 0
	, y0 = .45
	, x1 = d3agg$delta[d3agg$geno == 'tg']
	, y1 = .45
	, length = .17
	, angle = 90
	, code = 3
	, lwd = lnwd
)
text(x = d3agg$delta[d3agg$geno == 'tg'] / 2
	, y = .475
	, 'd\''
	, font = 3
	, cex = .5
)
# text(x = c3agg$crit[d3agg$geno == 'tg']
	# , y = .45
	# , 'c'
	# , font = 3
# )
# abline(h=.45)
if(save.plots){
	dev.off()
}



# 3xWT

if(save.plots){
	pdf(file = 'sdt plots/3xWTsdt.pdf'
		, width = 7
		, height = 4.5
	)
}
par(mar = c(4,2,2,2) + .01
	, cex = cmult
)

plot(NULL
	, xlim = c(-3, 4.7)
	, ylim = c(0, .52)
	, ylab = ''
	, xlab = ''
	, main = '3xWT'
	, bty = 'n'
	, yaxt = 'n'
)
abline(h=0)

curve(dnorm(x)	
	, -5
	, 5
	, n = 5000
	, add = T
	, lwd = lnwd
)
curve(dnorm(x
		, mean = d3agg$delta[d3agg$geno == 'wt']	
	)
	, -5
	, 5
	, add = T
	, lty = 'dashed'
	, lwd = lnwd
)
segments(x0 = d3agg$delta[d3agg$geno == 'wt'] / 2 + c3agg$crit[d3agg$geno == 'wt']
	, y0 = 0
	, x1 = d3agg$delta[d3agg$geno == 'wt'] / 2 + c3agg$crit[d3agg$geno == 'wt']
	, y1 = dnorm(d3agg$delta[d3agg$geno == 'wt'] / 2 + c3agg$crit[d3agg$geno == 'wt'])
	, lwd = lnwd

)
arrows(x0 = 0
	, y0 = .45
	, x1 = d3agg$delta[d3agg$geno == 'wt']
	, y1 = .45
	, length = .17
	, angle = 90
	, code = 3
	, lwd = lnwd
)
text(x = d3agg$delta[d3agg$geno == 'wt'] / 2
	, y = .475
	, 'd\''
	, font = 3
	, cex = .5
)
# # text(x = c3agg$crit[d3agg$geno == 'wt']
	# , y = .45
	# , 'c'
	# , font = 3
# )
# abline(h=.45)
if(save.plots){
	dev.off()
}



## 5x 
d5agg <- aggregate(delta ~ geno, FUN = mean, data = ndat[ndat$strain == '5x',])
c5agg <- aggregate(crit ~ geno, FUN = mean, data = ndat[ndat$strain == '5x',])

# 5xTG

if(save.plots){
	pdf(file = 'sdt plots/5xTGsdt.pdf'
		, width = 7
		, height = 4.5
	)
}
par(mar = c(4,2,2,2) + .01
	, cex = cmult
)

plot(NULL
	, xlim = c(-3, 4.7)
	, ylim = c(0, .52)
	, ylab = ''
	, xlab = ''
	, main = '5xTG'
	, bty = 'n'
	, yaxt = 'n'
)
abline(h=0)

curve(dnorm(x)	
	, -5
	, 5
	, n = 5000
	, add = T
	, lwd = lnwd
)
curve(dnorm(x
		, mean = d5agg$delta[d5agg$geno == 'tg']	
	)
	, -5
	, 5
	, add = T
	, lty = 'dashed'
	, lwd = lnwd
)
segments(x0 = d5agg$delta[d5agg$geno == 'tg'] / 2 + c5agg$crit[d5agg$geno == 'tg']
	, y0 = 0
	, x1 = d5agg$delta[d5agg$geno == 'tg'] / 2 + c5agg$crit[d5agg$geno == 'tg']
	, y1 = dnorm(d5agg$delta[d3agg$geno == 'tg'] / 2 + c5agg$crit[d3agg$geno == 'tg'])
	, lwd = lnwd

)
arrows(x0 = 0
	, y0 = .45
	, x1 = d5agg$delta[d5agg$geno == 'tg']
	, y1 = .45
	, length = .17
	, angle = 90
	, code = 3
	, lwd = lnwd
)
text(x = d5agg$delta[d5agg$geno == 'tg'] / 2
	, y = .475
	, 'd\''
	, font = 3
	, cex = .5
)
# text(x = c5agg$crit[d5agg$geno == 'tg']
	# , y = .45
	# , 'c'
	# , font = 3
# )
# abline(h=.45)
if(save.plots){
	dev.off()
}


# 5xWT

if(save.plots){
	pdf(file = 'sdt plots/5xWTsdt.pdf'
		, width = 7
		, height = 4.5
	)
}
par(mar = c(4,2,2,2) + .01
	, cex = cmult
)

plot(NULL
	, xlim = c(-3, 4.7)
	, ylim = c(0, .52)
	, ylab = ''
	, xlab = ''
	, main = '5xWT'
	, bty = 'n'
	, yaxt = 'n'
)
abline(h=0)

curve(dnorm(x)	
	, -5
	, 5
	, n = 5000
	, add = T
	, lwd = lnwd
)
curve(dnorm(x
		, mean = d5agg$delta[d5agg$geno == 'wt']	
	)
	, -5
	, 5
	, add = T
	, lty = 'dashed'
	, lwd = lnwd
)
segments(x0 = d5agg$delta[d5agg$geno == 'wt'] / 2 + c5agg$crit[d5agg$geno == 'wt']
	, y0 = 0
	, x1 = d5agg$delta[d5agg$geno == 'wt'] / 2 + c5agg$crit[d5agg$geno == 'wt']
	, y1 = dnorm(d5agg$delta[d3agg$geno == 'wt'] / 2 + c5agg$crit[d3agg$geno == 'wt'])
	, lwd = lnwd
)
arrows(x0 = 0
	, y0 = .45
	, x1 = d5agg$delta[d5agg$geno == 'wt']
	, y1 = .45
	, length = .17
	, angle = 90
	, code = 3
	, lwd = lnwd
)
text(x = d5agg$delta[d5agg$geno == 'wt'] / 2
	, y = .475
	, 'd\''
	, font = 3
	, cex = .5
)
# text(x = c5agg$crit[d5agg$geno == 'wt']
	# , y = .45
	# , 'c'
	# , font = 3
# )
# abline(h=.45)
if(save.plots){
	dev.off()
}