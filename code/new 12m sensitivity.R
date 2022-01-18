#### load packages ####
library(lme4)
library(nlme)
library(Hmisc)
library(tidyverse)
library(lubridate)
library(MuMIn)
library(PerformanceAnalytics)
library(here)
library(cowplot)

#### read in data ####

# create list of data files
fileList <- list.files(here("data", "12m data"), full.names = TRUE)

# read in data files from list
dat <- lapply(fileList, function(x) {
  a <- read.csv(x, nrow = 2)
  if (ncol(a) == 25) {
    a <- read.csv(x, comment.char = "#")
  } else {
    a <- read.csv(x, skip = 6, comment.char = "#")
  }
  if (is.numeric(a$Podour)) {
    a$Podour <- as.character(a$Podour)
  }
  a
})
length(dat) == length(fileList)

# combine files and convert to tibble
dat <- bind_rows(dat) %>%
  as_tibble(dat) %>%
  filter(!is.na(animal))
dat

# fix typos in data
dat$StudyName[dat$StudyName == "REV 1"] <- "REV1"
dat$StudyName[dat$StudyName == "OP0"] <- "INTRO"
dat$Podour[dat$Podour %in% c("1PPM EA ", "1ppm EEA", "1PPM EA", "1 PPM EA")] <- "1ppm EA"
dat$Podour[dat$Podour %in% c("0.1 PPM EA", "0.1PPM EA ", "0.1PPM EA", "0.1PM EA")] <- "0.1ppm EA"
dat$Podour[dat$Podour %in% c("0.01PPM EA", "0.01PPM EA ")] <- "0.01ppm EA"
dat$Podour[dat$Podour %in% c("0.001PPM EA", "0.001", "0.OO1PPM EA")] <- "0.001ppm EA"
dat$Podour[dat$Podour %in% c("0.0001PPM EA")] <- "0.0001ppm EA"
dat$Podour[dat$Podour %in% c("0.00001PPM EA")] <- "0.00001ppm EA"
dat$Nodour[dat$Nodour == "BANK"] <- "BLANK"
dat$Nodour[dat$Nodour == "SEARMINT"] <- "SPEARMINT"
dat$Podour[dat$Podour == "SEARMINT"] <- "SPEARMINT"
dat$Podour[dat$Podour == "ATCHOULI"] <- "PATCHOULI"
dat$Podour[dat$Podour == "0.00001ppm EA" & dat$StudyName == "SEN5"] <- "0.0001ppm EA"
dat <- dat[-which(dat$animal == "2988" & dat$StudyName == "SEN2"), ]

# make factors
dat$animal <- as.factor(dat$animal)
dat$Podour <- as.factor(as.character(dat$Podour))
dat$Nodour <- as.factor(as.character(dat$Nodour))

# add genotypes and strains
geno <- read_csv(here("data", "12m geno.csv"),
  col_types = "fffDf"
)
dat <- dat %>%
  left_join(geno)

table(dat$experimenter)

#### process data ####

# intro data
dat5xIntro <- dat %>%
  filter(
    strain == "5x",
    StudyName == "INTRO",
    SS == FALSE
  ) %>%
  group_by(animal, Podour, block, geno, sex)

# OP1 data
dat5xOP1 <- dat %>%
  filter(
    strain == "5x",
    StudyName == "OP1",
    SS == FALSE
  ) %>%
  group_by(animal, Podour, block, geno, sex)

# rev data
dat5xRev <- dat %>%
  filter(
    strain == "5x",
    StudyName == "REV1",
    SS == FALSE
  ) %>%
  group_by(animal, Podour, block, geno, sex)

# create list of concentrations
concList <- 10^-(c(0:5, 5))

# get 5xFAD sensitivity data, drop training and short sample trials
dat5x <- dat %>%
  filter(
    strain == "5x",
    StudyName %in% paste0("SEN", 1:7),
    SS == FALSE
  ) %>%
  group_by(animal, StudyName, Podour, block, geno, sex)

# summarize 5x data
dat5xSum <- dat5x %>%
  summarise(
    "acc" = sum(correct) / length(correct),
    "crit" = acc >= .85
  )
#### Intro analysis ####

# number of errors
introErrors <- dat5xIntro %>%
  group_by(animal, geno, sex) %>%
  summarise("errors" = sum(correct == FALSE))

mfIntroErrors <- lme(errors ~ geno * sex,
  random = ~ 1 | animal,
  data = introErrors,
  method = "ML"
)

mNullIntroErrors <- lme(errors ~ 1,
  random = ~ 1 | animal,
  data = introErrors,
  method = "ML"
)

options(na.action = "na.fail")
(T1 <- now())
introErrorsTable <- dredge(mfIntroErrors)
now() - T1
introErrorsTable

#### OP1 analysis ####

# number of errors
op1Errors <- dat5xOP1 %>%
  group_by(animal, geno, sex) %>%
  summarise("errors" = sum(correct == FALSE))

mfOp1Errors <- lme(errors ~ geno * sex,
  random = ~ 1 | animal,
  data = op1Errors,
  method = "ML"
)

mNullOp1Errors <- lme(errors ~ 1,
  random = ~ 1 | animal,
  data = op1Errors,
  method = "ML"
)

options(na.action = "na.fail")
(T1 <- now())
op1ErrorsTable <- dredge(mfOp1Errors)
now() - T1
op1ErrorsTable

#### Rev analysis ####

# number of errors
revErrors <- dat5xRev %>%
  group_by(animal, geno, sex) %>%
  summarise("errors" = sum(correct == FALSE))

mfRevErrors <- lme(errors ~ geno * sex,
  random = ~ 1 | animal,
  data = revErrors,
  method = "ML"
)

mNullRevErrors <- lme(errors ~ 1,
  random = ~ 1 | animal,
  data = revErrors,
  method = "ML"
)

options(na.action = "na.fail")
(T1 <- now())
revErrorsTable <- dredge(mfRevErrors)
now() - T1
revErrorsTable


#### sensitivity ####
#### model mice hiting crit ####

# did mice reach crit at least once per odour
dat5xCrit <- dat5xSum %>%
  group_by(animal, Podour, geno, sex) %>%
  summarise("Crit" = any(crit)) %>%
  mutate(
    "nChar" = nchar(as.character(Podour)),
    "conc" = as.numeric(substr(Podour, 1, nChar - 6))
  )

# model using conc as num
mfCrit <- glmer(Crit ~ conc * geno * sex + (1 | animal),
  data = dat5xCrit,
  family = binomial
)

mNullCrit <- glmer(Crit ~ (1 | animal),
  data = dat5xCrit,
  family = binomial
)

options(na.action = "na.fail")
(T1 <- now())
critTable <- dredge(mfCrit)
now() - T1
critTable

mSelectedCrit <- glmer(Crit ~ conc * sex + (1 | animal),
  data = dat5xCrit,
  family = binomial
)

anova(mSelectedCrit, mNullCrit, test = "Chisq")
summary(mSelectedCrit)
confSelectedCrit <- confint(mSelectedCrit)
confSelectedCrit

# avgCrit <- model.avg(critTable
#                      , subset = delta < 2)
# confint(avgCrit)
#
# summary(avgCrit)


#### plot mice hitting crit ####

plotP85 <- dat5xCrit %>%
  group_by(geno, sex, conc) %>%
  summarise(
    "crit" = mean(Crit),
    "ciUpper" = smean.cl.boot(Crit)[3],
    "ciLower" = smean.cl.boot(Crit)[2]
  ) %>%
  ungroup() %>%
  mutate(
    "xPos" = 1 + log10(1 / conc),
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5XFAD",
      "B6SJL"
    )
  ) %>%
  ggplot(aes(
    x = xPos,
    y = crit,
    shape = paste(Genotype, sex)
  )) +
  # geom_errorbar(aes(ymin = ciLower
  #                   , ymax = ciUpper)
  #               , linetype = 1
  #               , width = .2
  #               , position = position_dodge(width = .5)) +
  geom_point(
    size = 2,
    position = position_dodge(width = .4)
  ) +
  theme_classic() +
  scale_shape_manual(values = c(1, 2, 16, 17)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = concList[1:6]
  ) +
  labs(
    y = "Proportion reaching 85%",
    x = "Odour Concentration (ppm)",
    linetype = "Genotype",
    shape = "Genotype",
    colour = "Sex"
  ) +
  theme(
    legend.position = c(.9, .2),
    legend.title = element_blank()
  )

#### plot acc on best block ####
plotBestBlock <- dat5xSum %>%
  mutate(
      acc = acc * 100,
    "nChar" = nchar(as.character(Podour)),
    "conc" = as.numeric(substr(Podour, 1, nChar - 6)),
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5XFAD",
      "B6SJL"
    )
  ) %>%
  group_by(animal, conc, Genotype, sex) %>%
  summarise("best" = max(acc)) %>%
  mutate("xPos" = 1 + log10(1 / conc)) %>%
  ggplot(aes(
    x = xPos,
    y = best,
    shape = paste(Genotype, sex)
  )) +
  geom_point(
    size = 2,
    position = position_dodge(width = .4)
  ) +
  geom_hline(
    yintercept = 85,
    linetype = "dashed"
  ) +
  theme_classic() +
  scale_shape_manual(values = c(1, 2, 16, 17)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = concList[1:6]
  ) +
  labs(
    y = "Accuracy (%)",
    x = "Odour Concentration (ppm)",
    colour = "Sex"
  ) +
  theme(
    legend.position = c(.9, .2),
    legend.title = element_blank()
  )

#### model correct responses ####
# model correct responses for each trial
dat5xConc <- dat5x %>%
  ungroup() %>%
  mutate(
    "nChar" = nchar(as.character(Podour)),
    "conc" = as.numeric(substr(Podour, 1, nChar - 6)),
    "block" = ifelse(StudyName == "SEN7",
      block + 5,
      block
    )
  )
# start with full model, this part will be slow
# mfCorrect <- glmer(correct ~ conc * block * geno * sex + (1 | animal),
#   data = dat5xConc,
#   family = binomial,
#   control = glmerControl(
#     optimizer = "bobyqa",
#     optCtrl = list(maxfun = 2e5)
#   )
# )
#
# mNullCorrect <- glmer(correct ~ +(1 | animal),
#   data = dat5xConc,
#   family = binomial,
#   control = glmerControl(
#     optimizer = "bobyqa",
#     optCtrl = list(maxfun = 2e5)
#   )
# )
#
# # takes ~ 3.3h to run
# options(na.action = "na.fail")
# (T1 <- now())
# correctTable <- dredge(mfCorrect)
# now() - T1
# correctTable
#
# # save(correctTable, file = 'correctTable.RData')
#
# anova(mfCorrect, mNullCorrect, test = "Chisq")
# summary(mfCorrect)
# confmfCorrect <- confint(mfCorrect)

# save(confmfCorrect, file = 'correctConfInt.RData')


# model.avg(correctTable, subset = delta < 2)

#### model Accuracy ####

acc5xdat <- dat5xConc %>%
  group_by(animal, sex, geno, block, conc) %>%
  # filter(conc == concList[7]
  #        , block == 10) %>%
  summarise("acc" = mean(correct)) %>%
  mutate("concLog" = log10(conc))

mfAcc <- lme(acc ~ sex * geno * block * conc,
  random = ~ 1 | animal,
  data = acc5xdat,
  method = "ML"
)

accTable <- dredge(mfAcc)
accTable

mAccSelected <- lme(acc ~ (sex + geno + block + conc)^2 - conc:geno,
  random = ~ 1 | animal,
  data = acc5xdat,
  method = "ML"
)
summary(mAccSelected)
intervals(mAccSelected)

mNullAcc <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = acc5xdat,
  method = "ML"
)

anova(mAccSelected, mNullAcc)
AICc(mNullAcc)

aggregate(acc ~ geno,
  FUN = mean,
  data = acc5xdat
)
aggregate(acc ~ geno,
  FUN = sd,
  data = acc5xdat
)

aggregate(acc ~ sex,
  FUN = mean,
  data = acc5xdat
)
aggregate(acc ~ sex,
  FUN = sd,
  data = acc5xdat
)

aggregate(acc ~ conc,
  FUN = mean,
  data = acc5xdat
)
aggregate(acc ~ conc,
  FUN = sd,
  data = acc5xdat
)

aggregate(acc ~ block,
  FUN = mean,
  data = acc5xdat
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5xdat
)

# by sex
mAccSelectedMales <- lme(acc ~ (geno + block + conc)^2 - conc:geno,
  random = ~ 1 | animal,
  data = acc5xdat[acc5xdat$sex == "m", ],
  method = "ML"
)
summary(mAccSelectedMales)
aggregate(acc ~ geno,
  FUN = mean,
  data = acc5xdat[acc5xdat$sex == "m", ]
)

mAccSelectedFemales <- lme(acc ~ (geno + block + conc)^2 - conc:geno,
  random = ~ 1 | animal,
  data = acc5xdat[acc5xdat$sex == "f", ],
  method = "ML"
)
summary(mAccSelectedFemales)
aggregate(acc ~ geno,
  FUN = mean,
  data = acc5xdat[acc5xdat$sex == "f", ]
)

# by geno
mAccSelectedWt <- lme(acc ~ (sex + block + conc)^2,
  random = ~ 1 | animal,
  data = acc5xdat[acc5xdat$geno == "wt", ],
  method = "ML"
)
summary(mAccSelectedWt)

mAccSelectedTg <- lme(acc ~ (sex + block + conc)^2,
  random = ~ 1 | animal,
  data = acc5xdat[acc5xdat$geno == "tg", ],
  method = "ML"
)
summary(mAccSelectedTg)

aggregate(acc ~ sex,
  data = acc5xdat[acc5xdat$geno == "tg", ],
  FUN = mean
)
aggregate(acc ~ sex,
  data = acc5xdat[acc5xdat$geno == "wt", ],
  FUN = mean
)

aggregate(acc ~ sex,
  data = acc5xdat[acc5xdat$geno == "tg", ],
  FUN = sd
)
aggregate(acc ~ sex,
  data = acc5xdat[acc5xdat$geno == "wt", ],
  FUN = sd
)

# geno and sex
mAccSelectedWtM <- lme(acc ~ (block + conc)^2,
  random = ~ 1 | animal,
  data = acc5xdat[acc5xdat$geno == "wt" & acc5xdat$sex == "m", ],
  method = "ML"
)
summary(mAccSelectedWtM)

mAccSelectedTgM <- lme(acc ~ (block + conc)^2,
  random = ~ 1 | animal,
  data = acc5xdat[acc5xdat$geno == "tg" & acc5xdat$sex == "m", ],
  method = "ML"
)
summary(mAccSelectedTgM)

mAccSelectedWtF <- lme(acc ~ (block + conc)^2,
  random = ~ 1 | animal,
  data = acc5xdat[acc5xdat$geno == "wt" & acc5xdat$sex == "f", ],
  method = "ML"
)
summary(mAccSelectedWtF)

mAccSelectedTgF <- lme(acc ~ (block + conc)^2,
  random = ~ 1 | animal,
  data = acc5xdat[acc5xdat$geno == "tg" & acc5xdat$sex == "f", ],
  method = "ML"
)
summary(mAccSelectedTgF)

# block effect at first and last conc

aggregate(acc ~ block,
  FUN = mean,
  data = acc5xdat[acc5xdat$conc == 1, ]
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5xdat[acc5xdat$conc == 1, ]
)

aggregate(acc ~ block,
  FUN = mean,
  data = acc5xdat[acc5xdat$conc == 0.00001, ]
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5xdat[acc5xdat$conc == 0.00001, ]
)

# effect of odour conc for each sex

aggregate(acc ~ conc,
  FUN = mean,
  data = acc5xdat[acc5xdat$sex == "f", ]
)
aggregate(acc ~ conc,
  FUN = sd,
  data = acc5xdat[acc5xdat$sex == "f", ]
)

aggregate(acc ~ conc,
  FUN = mean,
  data = acc5xdat[acc5xdat$sex == "m", ]
)
aggregate(acc ~ conc,
  FUN = sd,
  data = acc5xdat[acc5xdat$sex == "m", ]
)

# effect of block for each sex

aggregate(acc ~ block,
  FUN = mean,
  data = acc5xdat[acc5xdat$sex == "f", ]
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5xdat[acc5xdat$sex == "f", ]
)

aggregate(acc ~ block,
  FUN = mean,
  data = acc5xdat[acc5xdat$sex == "m", ]
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5xdat[acc5xdat$sex == "m", ]
)

# effect of block for each geno

aggregate(acc ~ block,
  FUN = mean,
  data = acc5xdat[acc5xdat$geno == "wt", ]
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5xdat[acc5xdat$geno == "wt", ]
)

aggregate(acc ~ block,
  FUN = mean,
  data = acc5xdat[acc5xdat$geno == "tg", ]
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5xdat[acc5xdat$geno == "tg", ]
)


#### Accuracy on first odour conc ####

accHighDat <- acc5xdat %>%
  filter(conc == concList[1])

mfAccHigh <- lme(acc ~ sex * geno * block,
  random = ~ 1 | animal,
  data = accHighDat,
  method = "ML"
)

mNullAccHigh <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = accHighDat,
  method = "ML"
)

accHighTable <- dredge(mfAccHigh)
accHighTable

mAccHighSelected <- lme(acc ~ sex * block,
  random = ~ 1 | animal,
  data = accHighDat,
  method = "ML"
)
summary(mAccHighSelected)
intervals(mAccHighSelected)
anova(mAccHighSelected, mNullAccHigh)

aggregate(acc ~ block,
  FUN = mean,
  data = accHighDat
)
aggregate(acc ~ block,
  FUN = sd,
  data = accHighDat
)

aggregate(acc ~ sex,
  FUN = mean,
  data = accHighDat
)
aggregate(acc ~ sex,
  FUN = sd,
  data = accHighDat
)

aggregate(acc ~ sex * block,
  FUN = mean,
  data = accHighDat
)
aggregate(acc ~ sex * block,
  FUN = sd,
  data = accHighDat
)

#### Acc on 2nd conc ####

acc2Dat <- acc5xdat %>%
  filter(conc == concList[2])

mfAcc2 <- lme(acc ~ sex * geno * block,
  random = ~ 1 | animal,
  data = acc2Dat,
  method = "ML"
)

mNullAcc2 <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = acc2Dat,
  method = "ML"
)

acc2Table <- dredge(mfAcc2)
acc2Table

mAcc2Selected <- lme(acc ~ sex + geno + block + block:geno + geno:sex,
  random = ~ 1 | animal,
  data = acc2Dat,
  method = "ML"
)

anova(mAcc2Selected, mNullAcc2)

aggregate(acc ~ block,
  FUN = mean,
  data = acc2Dat
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc2Dat
)

aggregate(acc ~ geno,
  FUN = mean,
  data = acc2Dat
)
aggregate(acc ~ geno,
  FUN = sd,
  data = acc2Dat
)

aggregate(acc ~ sex,
  FUN = mean,
  data = acc2Dat
)
aggregate(acc ~ sex,
  FUN = sd,
  data = acc2Dat
)

aggregate(acc ~ geno * block,
  FUN = mean,
  data = acc2Dat
)
aggregate(acc ~ geno * block,
  FUN = sd,
  data = acc2Dat
)

aggregate(acc ~ sex * geno,
  FUN = mean,
  data = acc2Dat
)
aggregate(acc ~ sex * geno,
  FUN = sd,
  data = acc2Dat
)

#### Acc on 3rd conc ####

acc3Dat <- acc5xdat %>%
  filter(conc == concList[3])

mfAcc3 <- lme(acc ~ sex * geno * block,
  random = ~ 1 | animal,
  data = acc3Dat,
  method = "ML"
)

mNullAcc3 <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = acc3Dat,
  method = "ML"
)

acc3Table <- dredge(mfAcc3)
acc3Table

mfAcc3Selected <- lme(acc ~ block,
  random = ~ 1 | animal,
  data = acc3Dat,
  method = "ML"
)

anova(mfAcc3Selected, mNullAcc3)

#### Acc on 4th conc ####

acc4Dat <- acc5xdat %>%
  filter(conc == concList[4])

mfAcc4 <- lme(acc ~ sex * geno * block,
  random = ~ 1 | animal,
  data = acc4Dat,
  method = "ML"
)

mNullAcc4 <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = acc4Dat,
  method = "ML"
)

acc4Table <- dredge(mfAcc4)
acc4Table

mfAcc4Selected <- lme(acc ~ sex + block,
  random = ~ 1 | animal,
  data = acc4Dat,
  method = "ML"
)

anova(mfAcc4Selected, mNullAcc4)

aggregate(acc ~ block,
  FUN = mean,
  data = acc4Dat
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc4Dat
)

aggregate(acc ~ sex,
  FUN = mean,
  data = acc4Dat
)
aggregate(acc ~ sex,
  FUN = sd,
  data = acc4Dat
)

aggregate(acc ~ sex * block,
  FUN = mean,
  data = acc4Dat
)
aggregate(acc ~ sex * block,
  FUN = sd,
  data = acc4Dat
)

#### Acc on 5th conc ####

acc5Dat <- acc5xdat %>%
  filter(conc == concList[5])

mfAcc5 <- lme(acc ~ sex * geno * block,
  random = ~ 1 | animal,
  data = acc5Dat,
  method = "ML"
)

mNullAcc5 <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = acc5Dat,
  method = "ML"
)

acc5Table <- dredge(mfAcc5)
acc5Table

mfAcc5Selected <- lme(acc ~ sex + block,
  random = ~ 1 | animal,
  data = acc5Dat,
  method = "ML"
)

anova(mfAcc5Selected, mNullAcc5)

aggregate(acc ~ block,
  FUN = mean,
  data = acc5Dat
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5Dat
)

aggregate(acc ~ sex,
  FUN = mean,
  data = acc5Dat
)
aggregate(acc ~ sex,
  FUN = sd,
  data = acc5Dat
)

aggregate(acc ~ sex * block,
  FUN = mean,
  data = acc5Dat
)
aggregate(acc ~ sex * block,
  FUN = sd,
  data = acc5Dat
)

#### Accuracy on last odour conc, first day ####

acc5.1Dat <- acc5xdat %>%
  filter(
    conc == concList[7],
    block < 6
  )

mfAcc5.1 <- lme(acc ~ sex * geno * block,
  random = ~ 1 | animal,
  data = acc5.1Dat,
  method = "ML"
)

mNullAcc5.1 <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = acc5.1Dat,
  method = "ML"
)

acc5.1Table <- dredge(mfAcc5.1)
acc5.1Table

mAcc5.1Selected <- lme(acc ~ block,
  random = ~ 1 | animal,
  data = acc5.1Dat,
  method = "ML"
)

anova(mAcc5.1Selected, mNullAcc5.1)

aggregate(acc ~ block,
  FUN = mean,
  data = acc5.1Dat
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5.1Dat
)


#### Accuracy on last odour conc, second day ####

acc5.2Dat <- acc5xdat %>%
  filter(
    conc == concList[7],
    block > 5
  )

mfAcc5.2 <- lme(acc ~ sex * geno * block,
  random = ~ 1 | animal,
  data = acc5.2Dat,
  method = "ML"
)

mNullAcc5.2 <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = acc5.2Dat,
  method = "ML"
)

acc5.2Table <- dredge(mfAcc5.2)
acc5.2Table

mAcc5.2Selected <- lme(acc ~ block,
  random = ~ 1 | animal,
  data = acc5.1Dat,
  method = "ML"
)

anova(mAcc5.1Selected, mNullAcc5.1)

aggregate(acc ~ block,
  FUN = mean,
  data = acc5.1Dat
)
aggregate(acc ~ block,
  FUN = sd,
  data = acc5.1Dat
)

#### Accuracy on last odour conc ####

accLowestDat <- acc5xdat %>%
  filter(conc == concList[7])

mfAccLowest <- lme(acc ~ sex * geno * block,
  random = ~ 1 | animal,
  data = accLowestDat,
  method = "ML"
)

mNullAccLowest <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = accLowestDat,
  method = "ML"
)

accLowestTable <- dredge(mfAccLowest)
accLowestTable

mAccLowestSelected1 <- lme(acc ~ block,
  random = ~ 1 | animal,
  data = accLowestDat,
  method = "ML"
)

mAccLowestSelected2 <- lme(acc ~ block * sex,
  random = ~ 1 | animal,
  data = accLowestDat,
  method = "ML"
)

summary(mAccLowestSelected1)
intervals(mAccLowestSelected1)
anova(mAccLowestSelected1, mNullAccLowest)

summary(mAccLowestSelected2)
intervals(mAccLowestSelected2)

# by sex

mAccLowestSelectedMales <- lme(acc ~ block,
  random = ~ 1 | animal,
  data = accLowestDat[accLowestDat$sex == "m", ],
  method = "ML"
)

mAccLowestSelectedFemales <- lme(acc ~ block,
  random = ~ 1 | animal,
  data = accLowestDat[accLowestDat$sex == "f", ],
  method = "ML"
)
summary(mAccLowestSelectedMales)
summary(mAccLowestSelectedFemales)

#### learning to learn ####

accLtoLDat <- acc5xdat %>%
  filter(block %in% c(1, 6)) %>%
  mutate("day" = ifelse(block == 6,
    -concLog + 1,
    -concLog
  ))

mfLtoL <- lme(acc ~ sex * geno * day,
  random = ~ 1 | animal,
  data = accLtoLDat,
  method = "ML"
)

mNullLtoL <- lme(acc ~ 1,
  random = ~ 1 | animal,
  data = accLtoLDat,
  method = "ML"
)

accLtoLTable <- dredge(mfLtoL)
accLtoLTable

mLtoLSelected <- lme(acc ~ sex + day,
  random = ~ 1 | animal,
  data = accLtoLDat,
  method = "ML"
)

anova(mLtoLSelected, mNullLtoL)

aggregate(acc ~ day,
  FUN = mean,
  data = accLtoLDat
)
aggregate(acc ~ day,
  FUN = sd,
  data = accLtoLDat
)

aggregate(acc ~ sex,
  FUN = mean,
  data = accLtoLDat
)
aggregate(acc ~ sex,
  FUN = sd,
  data = accLtoLDat
)

#### plot accuracy ####
datCorrPlot <- dat5xConc %>%
  group_by(block, geno, sex, conc) %>%
  summarise(
    "acc" = mean(correct),
    "ciUpper" = smean.cl.boot(correct)[3],
    "ciLower" = smean.cl.boot(correct)[2]
  ) %>%
  ungroup() %>%
  mutate(
    "xPos" = 5 * log10(1 / conc) + block + log10(1 / conc),
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  )

axisBreaks <- c((5 * log10(1 / concList) + 3 + log10(1 / concList))[1:5], 35.5)

#### plot 12m accuracy per conc ####

pointSize <- 2

test_multi_plot_acc <-
    datCorrPlot %>%
    mutate(acc = 100 * acc,
           ciUpper = 100 * ciUpper,
           ciLower = 100 * ciLower,
           conc = conc %>% as.factor %>% fct_inseq() %>% fct_rev()) %>%
    ggplot(aes(x = xPos,
               y = acc,
               shape = paste(Genotype, sex),
               linetype = paste(Genotype, sex))) +
    geom_point(size = pointSize) +
    geom_line() +
    geom_errorbar(aes(
        ymin = ciLower,
        ymax = ciUpper
    ),
    linetype = 1,
    width = .5
    ) +
    facet_wrap(vars(conc),
               scales = "free_x") +
    theme_classic() +
    scale_shape_manual(values = c(17, 16, 2, 1)) + # c(1, 2, 16, 17)) +
    scale_linetype_manual(values = c("dashed", "dashed", "solid", "solid")) +
    scale_y_continuous(breaks = seq(40, 100, by = 10)) +
    scale_x_continuous(
        breaks = axisBreaks,
        labels = c(
            format(concList[1],
                   nsmall = 0,
                   scientific = FALSE
            ),
            format(concList[2],
                   nsmall = 0,
                   scientific = FALSE
            ),
            format(concList[3],
                   nsmall = 0,
                   scientific = FALSE
            ),
            format(concList[4],
                   nsmall = 0,
                   scientific = FALSE
            ),
            format(concList[5],
                   nsmall = 0,
                   scientific = FALSE
            ),
            format(concList[6],
                   nsmall = 0,
                   scientific = FALSE
            )
        )
    ) +
    labs(
        y = "Accuracy (%)",
        x = "Odour Concentration (ppm)",
        linetype = "Genotype",
        shape = "Genotype",
        colour = "Sex"
    ) +
    theme(
        legend.position = c(.9, .7),
        legend.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )



plotAcc12m <- ggplot(
  data = filter(
    datCorrPlot,
    conc == 1e0
  ) %>%
      mutate(acc = 100 * acc,
             ciUpper = 100 * ciUpper,
             ciLower = 100 * ciLower),
  aes(
    x = xPos,
    y = acc,
    shape = paste(Genotype, sex),
    linetype = paste(Genotype, sex)
  )
) +
  geom_point(size = pointSize) +
  geom_line() +
  geom_errorbar(aes(
    ymin = ciLower,
    ymax = ciUpper
  ),
  linetype = 1,
  width = .5
  ) +
  geom_point(
    data = filter(
      datCorrPlot,
      conc == 1e-1
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    size = pointSize
  ) +
  geom_line(data = filter(
    datCorrPlot,
    conc == 1e-1
  ) %>%
      mutate(acc = 100 * acc,
             ciUpper = 100 * ciUpper,
             ciLower = 100 * ciLower)) +
  geom_errorbar(
    data = filter(
      datCorrPlot,
      conc == 1e-1
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    aes(
      ymin = ciLower,
      ymax = ciUpper
    ),
    linetype = 1,
    width = .5
  ) +
  geom_point(
    data = filter(
      datCorrPlot,
      conc == 1e-2
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    size = pointSize
  ) +
  geom_line(data = filter(
    datCorrPlot,
    conc == 1e-2
  ) %>%
      mutate(acc = 100 * acc,
             ciUpper = 100 * ciUpper,
             ciLower = 100 * ciLower)) +
  geom_errorbar(
    data = filter(
      datCorrPlot,
      conc == 1e-2
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    aes(
      ymin = ciLower,
      ymax = ciUpper
    ),
    linetype = 1,
    width = .5
  ) +
  geom_point(
    data = filter(
      datCorrPlot,
      conc == 1e-3
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    size = pointSize
  ) +
  geom_line(data = filter(
    datCorrPlot,
    conc == 1e-3
  ) %>%
      mutate(acc = 100 * acc,
             ciUpper = 100 * ciUpper,
             ciLower = 100 * ciLower)) +
  geom_errorbar(
    data = filter(
      datCorrPlot,
      conc == 1e-3
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    aes(
      ymin = ciLower,
      ymax = ciUpper
    ),
    linetype = 1,
    width = .5
  ) +
  geom_point(
    data = filter(
      datCorrPlot,
      conc == 1e-4
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    size = pointSize
  ) +
  geom_line(data = filter(
    datCorrPlot,
    conc == 1e-4
  ) %>%
      mutate(acc = 100 * acc,
             ciUpper = 100 * ciUpper,
             ciLower = 100 * ciLower)) +
  geom_errorbar(
    data = filter(
      datCorrPlot,
      conc == 1e-4
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    aes(
      ymin = ciLower,
      ymax = ciUpper
    ),
    linetype = 1,
    width = .5
  ) +
  geom_point(
    data = filter(
      datCorrPlot,
      conc == 1e-5,
      block < 6
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    size = pointSize
  ) +
  geom_line(data = filter(
    datCorrPlot,
    conc == 1e-5,
    block < 6
  ) %>%
      mutate(acc = 100 * acc,
             ciUpper = 100 * ciUpper,
             ciLower = 100 * ciLower)) +
  geom_errorbar(
    data = filter(
      datCorrPlot,
      conc == 1e-5,
      block < 6
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    aes(
      ymin = ciLower,
      ymax = ciUpper
    ),
    linetype = 1,
    width = .5
  ) +
  geom_point(
    data = filter(
      datCorrPlot,
      conc == 1e-5,
      block >= 6
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    size = pointSize
  ) +
  geom_line(data = filter(
    datCorrPlot,
    conc == 1e-5,
    block >= 6
  ) %>%
      mutate(acc = 100 * acc,
             ciUpper = 100 * ciUpper,
             ciLower = 100 * ciLower)) +
  geom_errorbar(
    data = filter(
      datCorrPlot,
      conc == 1e-5,
      block >= 6
    ) %>%
        mutate(acc = 100 * acc,
               ciUpper = 100 * ciUpper,
               ciLower = 100 * ciLower),
    aes(
      ymin = ciLower,
      ymax = ciUpper
    ),
    linetype = 1,
    width = .5
  ) +
  # geom_hline(aes(yintercept = .85)
  #            , linetype = 'dotted') +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) + # c(1, 2, 16, 17)) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid", "solid")) +
  scale_y_continuous(breaks = seq(40, 100, by = 10)) +
  scale_x_continuous(
    breaks = axisBreaks,
    labels = c(
      format(concList[1],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[2],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[3],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[4],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[5],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[6],
        nsmall = 0,
        scientific = FALSE
      )
    )
  ) +
  labs(
    y = "Accuracy (%)",
    x = "Odour Concentration (ppm)",
    linetype = "Genotype",
    shape = "Genotype",
    colour = "Sex"
  ) +
  theme(
    legend.position = c(.9, .2),
    legend.title = element_blank()
  )



#### model correct for lowest conc ####
dat5xLowest <- dat5xConc %>%
  filter(conc == 0.00001)
mfLowest <- glmer(correct ~ block * geno * sex + (1 | animal),
  data = dat5xLowest,
  family = binomial
)

mNullLowest <- glmer(correct ~ (1 | animal),
  data = dat5xLowest,
  family = binomial
)


options(na.action = "na.fail")
(T1 <- now())
lowestTable <- dredge(mfLowest)
now() - T1
lowestTable

mSelectedLowest <- glmer(correct ~ block + geno + sex + block:geno + block:sex + (1 | animal),
  data = dat5xLowest,
  family = binomial
)

anova(mSelectedLowest, mNullLowest, test = "Chisq")
summary(mSelectedLowest)
# confint(mSelectedLowest)
# 2.5 %      97.5 %
#   .sig01        0.94530523  1.74915815
# (Intercept)   1.80788476  3.92085305
# block        -0.02387467  0.10005092
# genotg       -1.63290092  0.59093018
# sexf         -2.06241662  0.20541797
# block:genotg -0.13831398 -0.01071236
# block:sexf    0.02722849  0.15369207

#### model errors per conc ####

dat5xErrors <- dat5x %>%
  group_by(animal, Podour, geno, sex) %>%
  summarise("errors" = sum(correct == FALSE)) %>%
  mutate(
    "nChar" = nchar(as.character(Podour)),
    "conc" = as.numeric(substr(Podour, 1, nChar - 6))
  )

mfErrors <- lme(errors ~ geno * sex * conc,
  random = ~ 1 | animal,
  data = dat5xErrors,
  method = "ML"
)

mNullErrors <- lme(errors ~ 1,
  random = ~ 1 | animal,
  data = dat5xErrors,
  method = "ML"
)

options(na.action = "na.fail")
(T1 <- now())
errorsTable <- dredge(mfErrors)
now() - T1
errorsTable

mSelectedErrors1 <- lme(errors ~ sex * conc,
  random = ~ 1 | animal,
  data = dat5xErrors,
  method = "ML"
)

mSelectedErrors2 <- lme(errors ~ sex + conc + geno + conc:sex + geno:sex,
  random = ~ 1 | animal,
  data = dat5xErrors,
  method = "ML"
)

anova(mSelectedErrors1, mNullErrors)
summary(mSelectedErrors1)

anova(mSelectedErrors2, mNullErrors)
summary(mSelectedErrors2)


#### highest cons errors ####


dat5xHighErrors <- dat5x %>%
  filter(StudyName == "SEN1") %>%
  group_by(animal, geno, sex) %>%
  summarise("errors" = sum(correct == FALSE))

mfHighErrors <- lme(errors ~ geno * sex,
  random = ~ 1 | animal,
  data = dat5xHighErrors,
  method = "ML"
)

mNullHighErrors <- lme(errors ~ 1,
  random = ~ 1 | animal,
  data = dat5xHighErrors,
  method = "ML"
)

options(na.action = "na.fail")
(T1 <- now())
errorsHighTable <- dredge(mfHighErrors)
now() - T1
errorsHighTable

mfHighMalesErrors <- lme(errors ~ geno,
  random = ~ 1 | animal,
  data = dat5xHighErrors[dat5xHighErrors$sex == "m", ],
  method = "ML"
)
summary(mfHighMalesErrors)

mfHighFemalesErrors <- lme(errors ~ geno,
  random = ~ 1 | animal,
  data = dat5xHighErrors[dat5xHighErrors$sex == "f", ],
  method = "ML"
)
summary(mfHighFemalesErrors)

#### corr errors on intro with sen1 ####
cor.test(introErrors$errors, dat5xHighErrors$errors)


##### block 5 of each day ####

block5dat <- dat5x %>%
  filter(block %in% c(5, 10)) %>%
  mutate("day" = as.numeric(substr(StudyName, 4, 4))) %>%
  group_by(animal, day, geno, sex) %>%
  summarise("acc" = mean(correct))

mfBlock5 <- lme(acc ~ sex * geno * day,
  random = ~ 1 | animal,
  data = block5dat,
  method = "ML"
)

block5Table <- dredge(mfBlock5)
block5Table


mfBlock5Males <- lme(acc ~ geno * day,
  random = ~ 1 | animal,
  data = block5dat[block5dat$sex == "m", ],
  method = "ML"
)
summary(mfBlock5Males)

mfBlock5Females <- lme(acc ~ geno * day,
  random = ~ 1 | animal,
  data = block5dat[block5dat$sex == "f", ],
  method = "ML"
)
summary(mfBlock5Females)

##### Best block of each day ####

bestBlockdat <- dat5x %>%
  mutate("day" = as.numeric(substr(StudyName, 4, 4))) %>%
  group_by(animal, block, day, geno, sex) %>%
  summarise("acc" = mean(correct)) %>%
  group_by(animal, day, geno, sex) %>%
  summarise("maxAcc" = max(acc)) # %>%
# mutate('xPos' = 5 * (day - 1) + day)

mfBest <- lme(maxAcc ~ sex * geno * day,
  random = ~ 1 | animal,
  data = bestBlockdat,
  method = "ML"
)

bestTable <- dredge(mfBest)
bestTable
#
# mBestSelect <- lme(maxAcc ~ sex + day
#                    , random = ~1|animal
#                    , data = bestBlockdat
#                    , method = 'ML')
#
# mBestNull <- lme(maxAcc ~ 1
#                    , random = ~1|animal
#                    , data = bestBlockdat
#                    , method = 'ML')
#
# anova(mBestSelect, mBestNull)
# AICc(mBestSelect, mBestNull)
#
# aggregate(maxAcc ~ day,
#           data = bestBlockdat,
#           FUN = mean)
#
# aggregate(maxAcc ~ day,
#           data = bestBlockdat,
#           FUN = sd)


mSexBest <- lme(maxAcc ~ sex,
  random = ~ 1 | animal,
  data = bestBlockdat,
  method = "ML"
)
summary(mSexBest)

plotBestAcc12m <- bestBlockdat %>%
  mutate(
      maxAcc = 100 * maxAcc,
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "geno" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  ggplot(aes(
    x = day,
    y = jitter(maxAcc,
      factor = 4
    ),
    shape = paste(geno, sex)
  )) +
  geom_point(position = position_dodge(width = .5)) +
  geom_hline(
    yintercept = 85,
    linetype = "dashed"
  ) +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) + #  c(1, 2, 16, 17)) +
  scale_x_continuous(
    breaks = c(1:5, 6.5),
    labels = c(
      format(concList[1],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[2],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[3],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[4],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[5],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[6],
        nsmall = 0,
        scientific = FALSE
      )
    )
  ) +
  labs(
    y = "Max accuracy (%)",
    x = "Odour Concentration (ppm)",
    linetype = "Genotype",
    shape = "Genotype",
    colour = "Sex"
  ) +
  theme(
    legend.position = c(.8, .2),
    legend.title = element_blank()
  ) +
  coord_cartesian(ylim = c(20, 100))

plot85Pr12m <- bestBlockdat %>%
  mutate(
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "geno" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  group_by(day, geno, sex) %>%
  summarise(pr85 = mean(maxAcc >= .85)) %>%
  ggplot(aes(
    x = day,
    y = pr85,
    shape = paste(geno, sex)
  )) +
  geom_point(position = position_dodge(width = .5)) +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) +
  scale_x_continuous(
    breaks = c(1:5, 6.5),
    labels = c(
      format(concList[1],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[2],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[3],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[4],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[5],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[6],
        nsmall = 0,
        scientific = FALSE
      )
    )
  ) +
  labs(
    y = "Proportion reaching 85%",
    x = "Odour Concentration (ppm)",
    linetype = "Genotype",
    shape = "Genotype",
    colour = "Sex"
  ) +
  theme(
    legend.position = c(.8, .2),
    legend.title = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1))


thresholdDat <- bestBlockdat %>%
  mutate(
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "geno" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  group_by(day, geno, sex) %>%
  mutate(pr85 = maxAcc >= .85)

mf85 <- glmer(pr85 ~ sex * geno * day + (1 | animal),
  family = binomial,
  data = thresholdDat
)

table85 <- dredge(mf85)
table85

mSelect85 <- glmer(pr85 ~ day + (1 | animal),
  family = binomial,
  data = thresholdDat
)

mNull85 <- glmer(pr85 ~ (1 | animal),
  family = binomial,
  data = thresholdDat
)

anova(mSelect85, mNull85)
AICc(mSelect85, mNull85)

aggregate(pr85 ~ day,
  data = thresholdDat,
  FUN = mean
)

aggregate(pr85 ~ day,
  data = thresholdDat,
  FUN = sd
)

binomSD <- function(x) {
  sqrt(length(x) * mean(x) * (1 - mean(x)))
}
aggregate(pr85 ~ day,
  data = thresholdDat,
  FUN = binomSD
)

#### look at just S- trials ####

sMinusdat <- dat5x %>%
  filter(TrialType == "N") %>%
  mutate("day" = as.numeric(substr(StudyName, 4, 4))) %>%
  group_by(animal, block, day, geno, sex) %>%
  summarise(
    "fa" = mean(1 - correct),
    acc = mean(correct)
  )

mfSMinus <- lme(fa ~ sex * block * day * geno,
  random = ~ 1 | animal,
  data = sMinusdat,
  method = "ML"
)

sMinusTable <- dredge(mfSMinus)
sMinusTable

mSMinusSelect <- lme(fa ~ sex + block + day + geno + block:day + day:geno + day:sex,
  random = ~ 1 | animal,
  data = sMinusdat,
  method = "ML"
)

summary(mSMinusSelect)

mSMinusNull <- lme(fa ~ 1,
  random = ~ 1 | animal,
  data = sMinusdat,
  method = "ML"
)

anova(mSMinusSelect, mSMinusNull)
AICc(mSMinusSelect, mSMinusNull)

# by geno
sMinusdat %>%
  group_by(geno) %>%
  summarise(
    mn = mean(fa),
    stdDv = sd(fa)
  )

# tg mice

sMinusTg <- sMinusdat %>%
  filter(geno == "tg")

mfSMinusTg <- lme(fa ~ sex * block * day,
  random = ~ 1 | animal,
  data = sMinusTg,
  method = "ML"
)

dredge(mfSMinusTg)

# wt mice

sMinusWt <- sMinusdat %>%
  filter(geno == "wt")

mfSMinusWt <- lme(fa ~ sex * block * day,
  random = ~ 1 | animal,
  data = sMinusWt,
  method = "ML"
)

dredge(mfSMinusWt)

# by day
sMinusdat %>%
  group_by(day) %>%
  summarise(
    mn = mean(fa),
    stdDv = sd(fa)
  )

sMinusdat %>%
  group_by(sex) %>%
  summarise(
    mn = mean(fa),
    stdDv = sd(fa)
  )

sMinusdat %>%
  group_by(geno) %>%
  summarise(
    mn = mean(fa),
    stdDv = sd(fa)
  )

# day 1

sMinusD1 <- sMinusdat %>%
  filter(day == 1)

mfSMinusD1 <- lme(fa ~ sex * block * geno,
  random = ~ 1 | animal,
  data = sMinusD1,
  method = "ML"
)

dredge(mfSMinusD1)

mSMinusD1Select <- lme(fa ~ sex * block,
  random = ~ 1 | animal,
  data = sMinusD1,
  method = "ML"
)

mSMinusD1Null <- lme(fa ~ 1,
  random = ~ 1 | animal,
  data = sMinusD1,
  method = "ML"
)

anova(mSMinusD1Select, mSMinusD1Null)
AICc(mSMinusD1Select, mSMinusD1Null)

sMinusD1 %>%
  group_by(block) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

sMinusD1 %>%
  filter(sex == "m") %>%
  group_by(block) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )
sMinusD1 %>%
  filter(sex == "f") %>%
  group_by(block) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

sMinusD1 %>%
  group_by(sex) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

# day 2

sMinusD2 <- sMinusdat %>%
  filter(day == 2)

mfSMinusD2 <- lme(fa ~ sex * block * geno,
  random = ~ 1 | animal,
  data = sMinusD2,
  method = "ML"
)

dredge(mfSMinusD2)

mSMinusD2Null <- lme(fa ~ 1,
  random = ~ 1 | animal,
  data = sMinusD2,
  method = "ML"
)

anova(mfSMinusD2, mSMinusD2Null)
AICc(mfSMinusD2, mSMinusD2Null)

sMinusD2 %>%
  group_by(geno) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

sMinusD2 %>%
  group_by(block) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

sMinusD2 %>%
  filter(sex == "m") %>%
  group_by(block) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )
sMinusD2 %>%
  filter(sex == "f") %>%
  group_by(block) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

sMinusD2 %>%
  group_by(sex) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

# day 3

sMinusD3 <- sMinusdat %>%
  filter(day == 3)

mfSMinusD3 <- lme(fa ~ sex * block * geno,
  random = ~ 1 | animal,
  data = sMinusD3,
  method = "ML"
)

dredge(mfSMinusD3)




# day 4

sMinusD4 <- sMinusdat %>%
  filter(day == 4)

mfSMinusD4 <- lme(fa ~ sex * block * geno,
  random = ~ 1 | animal,
  data = sMinusD4,
  method = "ML"
)

dredge(mfSMinusD4)

mSMinusD4Select <- lme(fa ~ block + geno,
  random = ~ 1 | animal,
  data = sMinusD4,
  method = "ML"
)

mSMinusD4Null <- lme(fa ~ 1,
  random = ~ 1 | animal,
  data = sMinusD4,
  method = "ML"
)

anova(mSMinusD4Select, mSMinusD4Null)
AICc(mSMinusD4Select, mSMinusD4Null)

sMinusD4 %>%
  group_by(block) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

sMinusD4 %>%
  group_by(sex) %>%
  summarize(
    famean = mean(fa),
    fasd = sd(fa)
  )

# day 5

sMinusD5 <- sMinusdat %>%
  filter(day == 5)

mfSMinusD5 <- lme(fa ~ sex * block * geno,
  random = ~ 1 | animal,
  data = sMinusD5,
  method = "ML"
)

dredge(mfSMinusD5)

# day 6

sMinusD6 <- sMinusdat %>%
  filter(day == 6)

mfSMinusD6 <- lme(fa ~ sex * block * geno,
  random = ~ 1 | animal,
  data = sMinusD6,
  method = "ML"
)

dredge(mfSMinusD6)

mSMinusD6Select <- lme(fa ~ block,
  random = ~ 1 | animal,
  data = sMinusD6,
  method = "ML"
)

mSMinusD6Null <- lme(fa ~ 1,
  random = ~ 1 | animal,
  data = sMinusD6,
  method = "ML"
)

anova(mSMinusD6Select, mSMinusD6Null)
AICc(mSMinusD6Select, mSMinusD6Null)

sMinusD6 %>%
  group_by(block) %>%
  summarise(
    famean = mean(fa),
    fasd = sd(fa)
  )

# day 7

sMinusD7 <- sMinusdat %>%
  filter(day == 7)

mfSMinusD7 <- lme(fa ~ sex * block * geno,
  random = ~ 1 | animal,
  data = sMinusD7,
  method = "ML"
)

dredge(mfSMinusD7)

mSMinusD7Select <- lme(fa ~ sex + block,
  random = ~ 1 | animal,
  data = sMinusD7,
  method = "ML"
)

mSMinusD7Null <- lme(fa ~ 1,
  random = ~ 1 | animal,
  data = sMinusD7,
  method = "ML"
)

anova(mSMinusD7Select, mSMinusD7Null)
AICc(mSMinusD7Select, mSMinusD7Null)

sMinusD7 %>%
  group_by(block) %>%
  summarise(
    famean = mean(fa),
    fasd = sd(fa)
  )

sMinusD7 %>%
  group_by(sex) %>%
  summarise(
    famean = mean(fa),
    fasd = sd(fa)
  )


sMinusSdat <- sMinusdat %>%
  group_by(block, day, geno, sex) %>%
  summarise("fa" = mean_cl_boot(1 - acc)) %>%
  mutate(
    "xPos" = 5 * (day - 1) + block,
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "geno" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  )

##### Plot False Alarms ####
plotFa12m <- ggplot(
  data = filter(
    sMinusSdat,
    day == 1
  ) %>%
    mutate(fa = 100 * fa),
  aes(
    x = xPos,
    y = fa$y,
    ymin = fa$ymin,
    ymax = fa$ymax,
    shape = paste(geno, sex),
    linetype = paste(geno, sex)
  )
) +
  geom_point(size = pointSize) +
  geom_line() +
  geom_errorbar(
    data = sMinusSdat %>%
        mutate(fa = 100 * fa),
    linetype = 1,
    width = .5
  ) +
  geom_point(
    data = filter(
      sMinusSdat,
      day == 2
    ) %>%
        mutate(fa = 100 * fa),
    size = pointSize
  ) +
  geom_line(data = filter(
    sMinusSdat,
    day == 2
  ) %>%
      mutate(fa = 100 * fa)) +
  geom_point(
    data = filter(
      sMinusSdat,
      day == 3
    ) %>%
        mutate(fa = 100 * fa),
    size = pointSize
  ) +
  geom_line(data = filter(
    sMinusSdat,
    day == 3
  ) %>%
      mutate(fa = 100 * fa)) +
  geom_point(
    data = filter(
      sMinusSdat,
      day == 4
    ) %>%
        mutate(fa = 100 * fa),
    size = pointSize
  ) +
  geom_line(data = filter(
    sMinusSdat,
    day == 4
  ) %>%
      mutate(fa = 100 * fa)) +
  geom_point(
    data = filter(
      sMinusSdat,
      day == 5
    ) %>%
        mutate(fa = 100 * fa),
    size = pointSize
  ) +
  geom_line(data = filter(
    sMinusSdat,
    day == 5
  ) %>%
      mutate(fa = 100 * fa)) +
  geom_point(
    data = filter(
      sMinusSdat,
      day == 6
    ) %>%
        mutate(fa = 100 * fa),
    size = pointSize
  ) +
  geom_line(data = filter(
    sMinusSdat,
    day == 6
  ) %>%
      mutate(fa = 100 * fa)) +
  geom_point(
    data = filter(
      sMinusSdat,
      day == 7
    ) %>%
        mutate(fa = 100 * fa),
    size = pointSize
  ) +
  geom_line(data = filter(
    sMinusSdat,
    day == 7
  ) %>%
      mutate(fa = 100 * fa)) +
  # geom_hline(aes(yintercept = .85)
  #            , linetype = 'dotted') +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) + # c(1, 2, 16, 17)) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid", "solid")) +
  scale_x_continuous(
    breaks = c(seq(3, 23, by = 5), 30.5) # axisBreaks
    , labels = c(
      format(concList[1],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[2],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[3],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[4],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[5],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[6],
        nsmall = 0,
        scientific = FALSE
      )
    )
  ) +
  labs(
    y = "False alarm rate (%)",
    x = "Odour Concentration (ppm)",
    linetype = "Genotype",
    shape = "Genotype",
    colour = "Sex"
  ) +
  theme(
    legend.position = c(.9, .85),
    legend.title = element_blank(),
    legend.background = element_blank()
  )

##### FAs as % of total errors ####

faDat <- dat5x %>%
  mutate("day" = as.numeric(substr(StudyName, 4, 4))) %>%
  group_by(animal, block, day, geno, sex) %>%
  summarise(
    "fa" = sum(!correct[TrialType == "N"]),
    errors = sum(1 - correct)
  ) %>%
  mutate(percentFA = fa / errors)

#### short samples across days ####

ssDat <- dat %>%
  filter(
    StudyName %in% paste0("SEN", 1:7),
    SS == TRUE
  ) %>%
  mutate("day" = as.numeric(substr(StudyName, 4, 4))) %>%
  group_by(animal, block, day, geno, sex) %>%
  summarise("N" = n()) # %>%
# group_by(block, day, geno, sex) %>%
# summarise('N' = mean(N)) %>%
# mutate('xPos' = 5 * (day - 1) + block)

ssmf <- lme(N ~ sex * block * day * geno,
  random = ~ 1 | animal,
  data = ssDat,
  method = "ML"
)
dredge((ssmf))

ssmSelected <- lme(N ~ sex + block + day + geno + block:day + block:sex + day:geno + block:day:geno,
  random = ~ 1 | animal,
  data = ssDat,
  method = "ML"
)
summary(ssmSelected)

ggplot(
  data = filter(
    ssDat,
    day == 1
  ),
  aes(
    x = xPos,
    y = N,
    shape = paste(geno, sex),
    linetype = paste(geno, sex)
  )
) +
  geom_point(size = pointSize) +
  geom_line() +
  geom_point(
    data = filter(
      ssDat,
      day == 2
    ),
    size = pointSize
  ) +
  geom_line(data = filter(
    ssDat,
    day == 2
  )) +
  geom_point(
    data = filter(
      ssDat,
      day == 3
    ),
    size = pointSize
  ) +
  geom_line(data = filter(
    ssDat,
    day == 3
  )) +
  geom_point(
    data = filter(
      ssDat,
      day == 4
    ),
    size = pointSize
  ) +
  geom_line(data = filter(
    ssDat,
    day == 4
  )) +
  geom_point(
    data = filter(
      ssDat,
      day == 5
    ),
    size = pointSize
  ) +
  geom_line(data = filter(
    ssDat,
    day == 5
  )) +
  geom_point(
    data = filter(
      ssDat,
      day == 6
    ),
    size = pointSize
  ) +
  geom_line(data = filter(
    ssDat,
    day == 6
  )) +
  geom_point(
    data = filter(
      ssDat,
      day == 7
    ),
    size = pointSize
  ) +
  geom_line(data = filter(
    ssDat,
    day == 7
  )) +
  geom_hline(aes(yintercept = .85),
    linetype = "dotted"
  ) +
  theme_classic() +
  scale_shape_manual(values = c(1, 2, 16, 17)) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid", "solid")) +
  scale_x_continuous(
    breaks = axisBreaks,
    labels = concList[1:6]
  ) +
  labs(
    y = "Short Sample (n)",
    x = "Odour Concentration (ppm)",
    linetype = "Genotype",
    shape = "Genotype",
    colour = "Sex"
  ) +
  theme(
    legend.position = "none" # c(.9, .2)
    , legend.title = element_blank()
  )

#### time to complete ####

timeDat <- dat5x %>%
  mutate("day" = as.numeric(substr(StudyName, 4, 4))) %>%
  group_by(animal, day, geno, sex) %>%
  summarise("length" = max(time))

mfTime <- lme(length ~ sex * day * geno,
  random = ~ 1 | animal,
  data = timeDat,
  method = "ML"
)
dredge(mfTime)


timeDatS <- timeDat %>%
  group_by(day, geno, sex) %>%
  summarise("length" = mean(length))

timeDatS %>%
  ggplot(aes(
    x = day,
    y = length,
    shape = paste(geno, sex)
  )) +
  geom_point()

#### sensitivity index for last day ####

senDat <- dat5x %>%
  filter(StudyName %in% "SEN7") %>%
  group_by(animal, geno, sex) %>%
  summarise(
    "hits" = sum(correct[TrialType == "P"]),
    "falseAlarms" = sum(!correct[TrialType == "N"]),
    "pTrials" = length(correct[TrialType == "P"]),
    "nTrials" = length(correct[TrialType == "N"]),
    "hitRate" = ifelse(hits != pTrials,
      hits / pTrials,
      hits / (pTrials + 1)
    ),
    "faRate" = ifelse(falseAlarms != 0,
      falseAlarms / nTrials,
      (falseAlarms + 1) / nTrials
    ),
    "deePrime" = qnorm(hitRate) - qnorm(faRate)
  )

mfDeePrime <- lm(deePrime ~ geno * sex,
  data = senDat
)

options(na.action = "na.fail")
deePrimeTable <- dredge(mfDeePrime)
deePrimeTable

plotDeePrimBox12m <- senDat %>%
  mutate(
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  ggplot(aes(
    x = paste(Genotype, sex),
    y = deePrime
  )) +
  geom_boxplot(fill = "grey") +
  theme_classic() +
  labs(
    x = "",
    y = expression(paste("Sensitivity Index (", italic("d'"), ")"))
  )

#### sensitivity index all data ####

senAllDat <- dat5x %>%
  mutate(
    "nChar" = nchar(as.character(Podour)),
    "conc" = as.numeric(substr(Podour, 1, nChar - 6))
  ) %>%
  # filter(block == 5) %>%
  group_by(animal, geno, sex, conc) %>%
  summarise(
    "hits" = sum(correct[TrialType == "P"]),
    "falseAlarms" = sum(!correct[TrialType == "N"]),
    "pTrials" = length(correct[TrialType == "P"]),
    "nTrials" = length(correct[TrialType == "N"]),
    "hitRate" = ifelse(hits != pTrials,
      hits / pTrials,
      hits / (pTrials + 1)
    ),
    "hitRate" = ifelse(hitRate == 0,
      1 / (pTrials + 1),
      hitRate
    ),
    "faRate" = ifelse(falseAlarms != 0,
      falseAlarms / nTrials,
      (falseAlarms + 1) / nTrials
    ),
    "deePrime" = qnorm(hitRate) - qnorm(faRate)
  )

senAllmf <- lme(deePrime ~ geno * sex * conc,
  data = senAllDat,
  random = ~ 1 | animal,
  method = "ML"
)

options(na.action = "na.fail")
allDeePrimeTable <- dredge(senAllmf)
allDeePrimeTable

senAllSelected <- lme(deePrime ~ sex + geno + conc + sex:conc + geno:sex,
  data = senAllDat,
  random = ~ 1 | animal,
  method = "ML"
)
senAllNull <- lme(deePrime ~ 1,
  data = senAllDat,
  random = ~ 1 | animal,
  method = "ML"
)
anova(senAllSelected, senAllNull)

senAllDat %>%
  group_by(geno) %>%
  summarise(
    mn = mean(deePrime),
    stndev = sd(deePrime)
  )

senAllDat %>%
  group_by(sex) %>%
  summarise(
    mn = mean(deePrime),
    stndev = sd(deePrime)
  )

#### plot 12m d' by conc ####

sen12mDat <- dat5x %>%
  mutate(
    "day" = as.numeric(substr(StudyName, 4, 4)),
    "nChar" = nchar(as.character(Podour)),
    "conc" = as.numeric(substr(Podour, 1, nChar - 6))
  ) %>%
  # filter(block == 5) %>%
  group_by(animal, geno, sex, day) %>%
  summarise(
    "hits" = sum(correct[TrialType == "P"]),
    "falseAlarms" = sum(!correct[TrialType == "N"]),
    "pTrials" = length(correct[TrialType == "P"]),
    "nTrials" = length(correct[TrialType == "N"]),
    "hitRate" = ifelse(hits != pTrials,
      hits / pTrials,
      hits / (pTrials + 1)
    ),
    "hitRate" = ifelse(hitRate == 0,
      1 / (pTrials + 1),
      hitRate
    ),
    "faRate" = ifelse(falseAlarms != 0,
      falseAlarms / nTrials,
      (falseAlarms + 1) / nTrials
    ),
    "deePrime" = qnorm(hitRate) - qnorm(faRate)
  )

plotDeePrime12m <- sen12mDat %>%
  ungroup() %>%
  mutate(
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  group_by(geno, sex, day) %>%
  # summarise('meanDeePrime' = mean(deePrime)) %>%
  ggplot(aes(
    x = day,
    y = deePrime,
    shape = paste(Genotype, sex)
  )) +
  stat_summary(
    fun.data = "mean_cl_boot",
    position = position_dodge(width = .25)
  ) +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) +
  scale_x_continuous(
    breaks = c(1:5, 6.5),
    labels = c(
      format(concList[1],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[2],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[3],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[4],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[5],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[6],
        nsmall = 0,
        scientific = FALSE
      )
    )
  ) +
  labs(
    y = expression(paste("Sensitivity Index (", italic("d'"), ")")),
    x = "Odour Concentration (ppm)",
    shape = ""
  ) +
  theme(
    legend.position = c(.9, .2),
    legend.title = element_blank()
  )

plotDeePrimeBox12m <- sen12mDat %>%
  ungroup() %>%
  mutate(
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  group_by(Genotype, sex) %>%
  ggplot(aes(
    x = paste(Genotype, sex),
    y = deePrime
  )) +
  geom_boxplot(fill = "grey") +
  theme_classic() +
  labs(
    x = "",
    y = expression(paste("Sensitivity Index (", italic("d'"), ")"))
  )

# model last 4 concentrations

sen4Dat <- senAllDat %>%
  filter(conc < .1)

sen4mf <- lme(deePrime ~ geno * sex * conc,
  data = sen4Dat,
  random = ~ 1 | animal,
  method = "ML"
)

options(na.action = "na.fail")
sen4DeePrimeTable <- dredge(sen4mf)
sen4DeePrimeTable

sen4Selected <- lme(deePrime ~ sex,
  data = sen4Dat,
  random = ~ 1 | animal,
  method = "ML"
)

sen4Null <- lme(deePrime ~ 1,
  data = sen4Dat,
  random = ~ 1 | animal,
  method = "ML"
)

anova(sen4Selected, sen4Null)

# model all put 1st conc

sen5Dat <- senAllDat %>%
  filter(conc < 1)

sen5mf <- lme(deePrime ~ geno * sex * conc,
  data = sen5Dat,
  random = ~ 1 | animal,
  method = "ML"
)

dredge(sen5mf)

#### response bias for last day ####

biasDat <- dat5x %>%
  filter(StudyName == "SEN7") %>%
  group_by(animal, geno, sex) %>%
  summarise(
    "hits" = sum(correct[TrialType == "P"]),
    "falseAlarms" = sum(!correct[TrialType == "N"]),
    "pTrials" = length(correct[TrialType == "P"]),
    "nTrials" = length(correct[TrialType == "N"]),
    "hitRate" = ifelse(hits != pTrials,
      hits / pTrials,
      hits / (pTrials + 1)
    ),
    "faRate" = ifelse(falseAlarms != 0,
      falseAlarms / nTrials,
      (falseAlarms + 1) / nTrials
    ),
    "respCrit" = -.5 * (qnorm(hitRate) + qnorm(faRate))
  )

mfRespCrit <- lm(respCrit ~ geno * sex,
  data = biasDat
)
respCritTable <- dredge(mfRespCrit)
respCritTable

respCritSelected <- lm(respCrit ~ sex,
  data = biasDat
)
mNullRespCrit <- lm(respCrit ~ 1,
  data = biasDat
)
summary(respCritSelected)
anova(respCritSelected, mNullRespCrit)
confint(respCritSelected)

aggregate(respCrit ~ sex,
  data = biasDat,
  FUN = mean
)

aggregate(respCrit ~ sex,
  data = biasDat,
  FUN = sd
)

plotResponseBiasBox12m <- biasDat %>%
  mutate(
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  ggplot(aes(
    x = paste(Genotype, sex),
    y = respCrit
  )) +
  geom_boxplot(fill = "grey") +
  theme_classic() +
  labs(
    x = "",
    y = expression(paste("Response Bias (", italic("c"), ")"))
  )

#### response bias for concs ####

biasConcDat <- dat5x %>%
  group_by(animal, geno, sex, StudyName) %>%
  summarise(
    "hits" = sum(correct[TrialType == "P"]),
    "falseAlarms" = sum(!correct[TrialType == "N"]),
    "pTrials" = length(correct[TrialType == "P"]),
    "nTrials" = length(correct[TrialType == "N"]),
    "hitRate" = ifelse(hits != pTrials,
      hits / pTrials,
      hits / (pTrials + 1)
    ),
    "faRate" = ifelse(falseAlarms != 0,
      falseAlarms / nTrials,
      (falseAlarms + 1) / nTrials
    ),
    "respCrit" = -.5 * (qnorm(hitRate) + qnorm(faRate))
  )

mfBiasConc <- lme(respCrit ~ geno * sex * StudyName,
  data = biasConcDat,
  random = ~ 1 | animal,
  method = "ML"
)

respConcTable <- dredge(mfBiasConc)
respConcTable

mBiasConcSelected <- lme(respCrit ~ sex + StudyName,
  data = biasConcDat,
  random = ~ 1 | animal,
  method = "ML"
)

mBiasConcNull <- lme(respCrit ~ 1,
  data = biasConcDat,
  random = ~ 1 | animal,
  method = "ML"
)

anova(mBiasConcSelected, mBiasConcNull)

aggregate(respCrit ~ StudyName,
  FUN = mean,
  data = biasConcDat
)

biasConcDat %>%
  ungroup() %>%
  summarise(mnC = mean(respCrit))

## drop day 1

biasConcM1Dat <- dat5x %>%
  filter(StudyName != "SEN1") %>%
  filter(StudyName != "SEN2") %>%
  group_by(animal, geno, sex, StudyName) %>%
  summarise(
    "hits" = sum(correct[TrialType == "P"]),
    "falseAlarms" = sum(!correct[TrialType == "N"]),
    "pTrials" = length(correct[TrialType == "P"]),
    "nTrials" = length(correct[TrialType == "N"]),
    "hitRate" = ifelse(hits != pTrials,
      hits / pTrials,
      hits / (pTrials + 1)
    ),
    "faRate" = ifelse(falseAlarms != 0,
      falseAlarms / nTrials,
      (falseAlarms + 1) / nTrials
    ),
    "respCrit" = -.5 * (qnorm(hitRate) + qnorm(faRate))
  )

mfBiasM1Conc <- lme(respCrit ~ geno * sex * StudyName,
  data = biasConcM1Dat,
  random = ~ 1 | animal,
  method = "ML"
)

respConcM1Table <- dredge(mfBiasM1Conc)
respConcM1Table

mBiasM1ConcSeleced <- lme(respCrit ~ sex,
  data = biasConcM1Dat,
  random = ~ 1 | animal,
  method = "ML"
)

mBiasM1ConcNull <- lme(respCrit ~ 1,
  data = biasConcM1Dat,
  random = ~ 1 | animal,
  method = "ML"
)

anova(mBiasM1ConcSeleced, mBiasM1ConcNull)

#### plot bias across days ####

plotResponseBiasConcBox12m <- biasConcDat %>%
  ungroup() %>%
  mutate(
    "day" = as.numeric(substr(StudyName, 4, 4)),
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  group_by(geno, sex, day) %>%
  # summarise('meanDeePrime' = mean(deePrime)) %>%
  ggplot(aes(
    x = day - 1,
    y = respCrit,
    shape = paste(Genotype, sex)
  )) +
  stat_summary(
    fun.data = "mean_cl_boot",
    position = position_dodge(width = .25)
  ) +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) +
  scale_x_continuous(
    breaks = c(0:4, 5.5),
    labels = c(
      format(concList[1],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[2],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[3],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[4],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[5],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[6],
        nsmall = 0,
        scientific = FALSE
      )
    )
  ) +
  labs(
    y = expression(paste("Response Bias (", italic("c"), ")")),
    x = "Odour Concentration (ppm)",
    shape = ""
  ) +
  theme(
    legend.position = c(.8, .8),
    legend.title = element_blank()
  )

#### Add frailty data ####

frailtySenDat <- read_csv("analysis/frailty scores.csv",
  col_types = "ffffdd"
) %>%
  left_join(senDat)

cor.test(~ frailtyScore + deePrime,
  data = frailtySenDat
)

frailtyBiasDat <- read_csv("analysis/frailty scores.csv",
  col_types = "ffffdd"
) %>%
  left_join(biasDat)

cor.test(~ frailtyScore + respCrit,
  data = frailtyBiasDat
)

dat5xErrorsTot <- dat5xErrors %>%
  group_by(animal) %>%
  summarise("totalErrors" = sum(errors))

frailtyErrorDat <- read_csv("analysis/frailty scores.csv",
  col_types = "ffffdd"
) %>%
  left_join(dat5xErrorsTot)

cor.test(~ frailtyScore + totalErrors,
  data = frailtyErrorDat
)

t.test(frailtyScore ~ geno, data = frailtyErrorDat)

#### Short Sample correlations ####

ssDat <- dat %>%
  filter(!(StudyName %in% c("INTRO", "OP1", "REV1"))) %>%
  group_by(animal, StudyName, geno, sex) %>%
  summarise(
    ShortSamples = sum(SS),
    HitRate = sum(correct[TrialType == "P" & !SS]) / sum(TrialType == "P" & !SS),
    FARate = sum(!correct[TrialType == "N" & !SS]) / sum(TrialType == "N" & !SS),
    qHitRate = if_else(HitRate == 1,
      .99,
      HitRate
    ),
    qHitRate = if_else(HitRate == 0,
      .01,
      qHitRate
    ),
    qFARate = if_else(FARate == 0,
      .01,
      FARate
    ),
    qFARate = if_else(FARate == 1,
      .99,
      qFARate
    ),
    respBias = -.5 * (qnorm(qHitRate) + qnorm(qFARate)),
    deePrime = qnorm(qHitRate) - qnorm(qFARate)
  ) %>%
  select(-qHitRate, -qFARate)

chart.Correlation.linear <-
  function(R, histogram = TRUE, method = c("pearson", "kendall", "spearman"), ...) { # @author R Development Core Team
    # @author modified by Peter Carl & Marek Lahoda
    # @auther further modified by Kyle Roddick
    # Visualization of a Correlation Matrix. On top the (absolute) value of the correlation plus the result
    # of the cor.test as stars. On botttom, the bivariate scatterplots, with a linear regression fit.
    # On diagonal, the histograms with probability, density and normal density (gaussian) distribution.

    x <- checkData(R, method = "matrix")

    if (missing(method)) method <- method[1] # only use one
    cormeth <- method

    # Published at http://addictedtor.free.fr/graphiques/sources/source_137.R
    panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", method = cormeth, cex.cor, ...) {
      usr <- par("usr")
      on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use = use, method = method) # MG: remove abs here
      txt <- format(c(r, 0.123456789), digits = digits)[1]
      txt <- paste(prefix, txt, sep = "")
      if (missing(cex.cor)) cex <- 0.8 / strwidth(txt)

      test <- cor.test(as.numeric(x), as.numeric(y), method = method)
      # borrowed from printCoefmat
      Signif <- symnum(test$p.value,
        corr = FALSE, na = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols = c("***", "**", "*", ".", " ")
      )
      # MG: add abs here and also include a 30% buffer for small numbers
      text(0.5, 0.5, txt, cex = cex * (abs(r) + .3) / 1.3)
      text(.8, .8, Signif, cex = cex, col = 2)
    }

    # remove method from dotargs
    dotargs <- list(...)
    dotargs$method <- NULL
    rm(method)

    hist.panel <- function(x, ... = NULL) {
      par(new = TRUE)
      hist(x,
        col = "light gray",
        probability = TRUE,
        axes = FALSE,
        main = "",
        breaks = "FD"
      )
      lines(density(x, na.rm = TRUE),
        col = "red",
        lwd = 1
      )
      # adding line representing density of normal distribution with parameters correponding to estimates of mean and standard deviation from the data
      ax.x <- seq(min(x), max(x), 0.1) # ax.x containts points corresponding to data range on x axis
      density.est <- dnorm(ax.x, mean = mean(x), sd = sd(x)) # density corresponding to points stored in vector ax.x
      # lines(ax.x, density.est, col = "blue", lwd = 1, lty = 1)                                # adding line representing density into histogram
      rug(x)
    }

    # Linear regression line fit over points
    reg <- function(x, y, ...) {
      points(x, y, ...)
      abline(lm(y ~ x), col = "red")
    }

    # Draw the chart
    if (histogram) {
      pairs(x, gap = 0, lower.panel = reg, upper.panel = panel.cor, diag.panel = hist.panel)
    } else {
      pairs(x, gap = 0, lower.panel = reg, upper.panel = panel.cor)
    }
  }

png(here("figures", "CorrelationMatrix.png"),
  width = 8.5,
  height = 8.5,
  units = "in",
  res = 300
)

ssDat %>%
  # filter(StudyName != 'SEN1') %>%
  ungroup() %>%
  select(
    `Short Samples` = ShortSamples,
    `Hit Rate` = HitRate,
    `False Alarm Rate` = FARate,
    # `Response Bias` = respBias,
    # `Sensitivity Index` = deePrime
  ) %>%
  chart.Correlation.linear(method = "spearman")

dev.off()

cor.test(ssDat$ShortSamples, ssDat$HitRate,
  method = "spearman"
)

cor.test(ssDat$ShortSamples, ssDat$FARate,
  method = "spearman"
)

cor.test(ssDat$ShortSamples, ssDat$respBias,
  method = "spearman"
)

cor.test(ssDat$ShortSamples, ssDat$deePrime,
  method = "spearman"
)

cor.test(ssDat$HitRate, ssDat$FARate,
  method = "spearman"
)

cor.test(ssDat$HitRate, ssDat$respBias,
  method = "spearman"
)

cor.test(ssDat$HitRate, ssDat$deePrime,
  method = "spearman"
)

cor.test(ssDat$FARate, ssDat$respBias,
  method = "spearman"
)

cor.test(ssDat$FARate, ssDat$deePrime,
  method = "spearman"
)

cor.test(ssDat$respBias, ssDat$deePrime,
  method = "spearman"
)

ssDat %>%
  # filter(StudyName != 'SEN1') %>%
  ungroup() %>%
  select(ShortSamples:deePrime) %>%
  chart.Correlation(method = "pearson")

ssDat %>%
  # filter(StudyName != 'SEN1') %>%
  ungroup() %>%
  select(ShortSamples:deePrime) %>%
  chart.Correlation(method = "kendal")



cor.test(ssDat$ShortSamples, ssDat$HitRate)

ssDat %>%
  ungroup() %>%
  filter(geno == "wt") %>%
  select(ShortSamples:deePrime) %>%
  chart.Correlation()

ssDat %>%
  ungroup() %>%
  filter(geno == "tg") %>%
  select(ShortSamples:deePrime) %>%
  chart.Correlation()

ssDat %>%
  ggplot(aes(
    x = ShortSamples,
    y = HitRate,
    colour = geno
  )) +
  geom_point() +
  geom_smooth(method = "lm")

ssDat %>%
  ggplot(aes(
    x = ShortSamples,
    y = respBias,
    colour = geno
  )) +
  geom_point() +
  geom_smooth(method = "lm")

ssDat %>%
  ggplot(aes(
    x = ShortSamples,
    y = deePrime,
    colour = geno
  )) +
  geom_point() +
  geom_smooth(method = "lm")


##### Short Sample analysis ####

mfSS <- lme(ShortSamples ~ StudyName * geno * sex,
  data = ssDat,
  random = ~ 1 | animal,
  method = "ML"
)

dredge(mfSS)

mSSSelect <- lme(ShortSamples ~ sex,
  data = ssDat,
  random = ~ 1 | animal,
  method = "ML"
)
mSSNull <- lme(ShortSamples ~ 1,
  data = ssDat,
  random = ~ 1 | animal,
  method = "ML"
)
anova(mSSSelect, mSSNull)

ssDat %>%
  group_by(sex) %>%
  summarise(
    ssMean = mean(ShortSamples),
    ssSD = sd(ShortSamples)
  )

ssDat %>%
  ungroup() %>%
  mutate(
    "day" = as.numeric(substr(StudyName, 4, 4)),
    "sex" = ifelse(sex == "m",
      "Male",
      "Female"
    ),
    "Genotype" = ifelse(geno == "tg",
      "5xFAD",
      "B6SJL"
    )
  ) %>%
  group_by(geno, sex, day) %>%
  # summarise('meanDeePrime' = mean(deePrime)) %>%
  ggplot(aes(
    x = day - 1,
    y = ShortSamples,
    shape = paste(Genotype, sex)
  )) +
  stat_summary(
    fun.data = "mean_cl_boot",
    position = position_dodge(width = .25)
  ) +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) +
  scale_x_continuous(
    breaks = c(0:4, 5.5),
    labels = c(
      format(concList[1],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[2],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[3],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[4],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[5],
        nsmall = 0,
        scientific = FALSE
      ),
      format(concList[6],
        nsmall = 0,
        scientific = FALSE
      )
    )
  ) +
  labs(
    y = "Short Samples (n)",
    x = "Odour Concentration (ppm)",
    shape = ""
  ) +
  theme(
    legend.position = c(.8, .8),
    legend.title = element_blank()
  )

ggsave(here("figures", "newNewPlots", "shortSamples.png"),
  width = 8.5,
  height = 5.5
)

#### Combined plots ####

# figure1 <- plot_grid(plotAcc12m, plotBestAcc12m + theme(legend.background = element_blank()), plot85Pr12m,
#   ncol = 1,
#   labels = "AUTO"
# )



ggsave(here("figures", "just12m", "figure1.eps"),
  plotDeePrimeBox12m,
  height = 5.5,
  width = 5.5,
  dpi = 1200
)


ggsave(here("figures", "just12m", "figure2.eps"),
  plotResponseBiasBox12m,
  height = 5.5,
  width = 5.5,
  dpi = 1200
)


# figure3 <- plot_grid(plotDeePrimeBox12m, plotDeePrime12m,
#   ncol = 1,
#   labels = "AUTO"
# )

ggsave(here("figures", "just12m", "figure3.eps"),
  plotAcc12m,
  height = 5.5,
  width = 5.5,
  dpi = 1200
)


# figure4 <- plot_grid(plotResponseBiasBox12m,
#   plotResponseBiasConcBox12m + theme(legend.background = element_blank(), legend.position = c(.9, .14)),
#   ncol = 1,
#   labels = "AUTO"
# )
# ggsave(here("figures", "just12m", "figure4.png"),
#   figure4,
#   height = 8.5,
#   width = 5.5
# )
