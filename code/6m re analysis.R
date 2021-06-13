library(tidyverse)
library(nlme)

#### read in data ####
MSdat <- read_csv( '../data/MSc sensitivity.csv'
                   , col_types = 'fdddd')
MSgeno <- read_csv( '../data/MSc genotypes.csv'
                    , col_types = 'fffff')
dat6m <- MSdat %>% 
  left_join(MSgeno) %>% 
  mutate('geno' = Alz) %>% 
  group_by(Mouse, Consentration, Block, geno, sex) %>% 
  select(-retinal_degeneration, -dysferlin, -Alz)


#### plot acc ####
datCorrPlot6m <- dat6m %>% 
  mutate(block = Block,
         conc = Consentration) %>% 
  group_by(block, geno, sex, conc) %>% 
  summarise('acc' = mean((Hits + Correct_rejections)/20)
            , 'ciUpper' = smean.cl.boot((Hits + Correct_rejections)/20)[3]
            , 'ciLower' = smean.cl.boot((Hits + Correct_rejections)/20)[2]) %>%
  ungroup() %>% 
  mutate('xPos' = 5 * log10(1/conc) + block + log10(1/conc)
         , 'sex' = ifelse(sex == 'm'
                          , 'Male'
                          , 'Female')
         , 'Genotype' = ifelse(geno == 'tg'
                               , '5xFAD'
                               , 'B6SJL'))

axisBreaks <- c((5 * log10(1/concList) + 3 + log10(1/concList))[1:5], 35.5)


pointSize <-  2
plotAcc6m <- ggplot(data = filter(datCorrPlot6m
                                   , conc == 1e0)
                     , aes(x = xPos
                           , y = acc
                           , shape = paste(Genotype, sex)
                           , linetype = paste(Genotype, sex))) +
  geom_point(size = pointSize) +
  geom_line() +
  geom_errorbar(aes(ymin = ciLower
                    , ymax = ciUpper)
                , linetype = 1
                , width = .5) +
  geom_point(data = filter(datCorrPlot6m
                           , conc == 1e-1)
             , size = pointSize) +
  geom_line(data = filter(datCorrPlot6m
                          , conc == 1e-1)) +
  geom_errorbar(data = filter(datCorrPlot6m
                              , conc == 1e-1)
                , aes(ymin = ciLower
                      , ymax = ciUpper)
                , linetype = 1
                , width = .5) +
  geom_point(data = filter(datCorrPlot6m
                           , conc == 1e-2)
             , size = pointSize) +
  geom_line(data = filter(datCorrPlot6m
                          , conc == 1e-2)) +
  geom_errorbar(data = filter(datCorrPlot6m
                              , conc == 1e-2)
                , aes(ymin = ciLower
                      , ymax = ciUpper)
                , linetype = 1
                , width = .5) +
  geom_point(data = filter(datCorrPlot6m
                           , conc == 1e-3)
             , size = pointSize) +
  geom_line(data = filter(datCorrPlot6m
                          , conc == 1e-3)) +
  geom_errorbar(data = filter(datCorrPlot6m
                              , conc == 1e-3)
                , aes(ymin = ciLower
                      , ymax = ciUpper)
                , linetype = 1
                , width = .5) +
  geom_point(data = filter(datCorrPlot6m
                           , conc == 1e-4)
             , size = pointSize) +
  geom_line(data = filter(datCorrPlot6m
                          , conc == 1e-4)) +
  geom_errorbar(data = filter(datCorrPlot6m
                              , conc == 1e-4)
                , aes(ymin = ciLower
                      , ymax = ciUpper)
                , linetype = 1
                , width = .5) +
  geom_point(data = filter(datCorrPlot6m
                           , conc == 1e-5
                           , block < 6)
             , size = pointSize) +
  geom_line(data = filter(datCorrPlot6m
                          , conc == 1e-5
                          , block < 6)) +
  geom_errorbar(data = filter(datCorrPlot6m
                              , conc == 1e-5
                              , block < 6)
                , aes(ymin = ciLower
                      , ymax = ciUpper)
                , linetype = 1
                , width = .5) +
  geom_point(data = filter(datCorrPlot6m
                           , conc == 1e-5
                           , block >= 6)
             , size = pointSize) +
  geom_line(data = filter(datCorrPlot6m
                          , conc == 1e-5
                          , block >= 6)) +
  geom_errorbar(data = filter(datCorrPlot6m
                              , conc == 1e-5
                              , block >= 6)
                , aes(ymin = ciLower
                      , ymax = ciUpper)
                , linetype = 1
                , width = .5) +
  # geom_hline(aes(yintercept = .85)
  #            , linetype = 'dotted') +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) + #c(1, 2, 16, 17)) +
  scale_linetype_manual(values = c('dashed', 'dashed', 'solid', 'solid')) +
  scale_y_continuous(breaks = seq(.4, 1, by = .1)) +
  scale_x_continuous(breaks = axisBreaks
                     , labels = c(format(concList[1]
                                         , nsmall = 0
                                         , scientific = FALSE)
                                  , format(concList[2]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[3]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[4]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[5]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[6]
                                           , nsmall = 0
                                           , scientific = FALSE))) +
  labs(y = 'Accuracy (%)'
       , x = 'Odour Concentration (ppm)'
       , linetype = 'Genotype'
       , shape = 'Genotype'
       , colour = 'Sex') +
  theme(legend.position = c(.9, .2)
        , legend.title = element_blank())


#### Best block ####

bestBlockdat6m <- dat6m %>% 
  mutate(day = -log10(Consentration) + 1,
         day = if_else(Block > 5,
                       7,
                       day)) %>%
  group_by(Mouse, Block, day, geno, sex) %>% 
  summarise('acc' = (Hits + Correct_rejections)/20) %>% 
  group_by(Mouse, day, geno, sex) %>% 
  summarise('maxAcc' = max(acc))

plotBestAcc6m <- bestBlockdat6m %>% 
  mutate('sex' = ifelse(sex == 'm'
                        , 'Male'
                        , 'Female')
         , 'geno' = ifelse(geno == 'tg'
                           , '5xFAD'
                           , 'B6SJL')) %>% 
  ggplot(aes(x = day,
             y = jitter(maxAcc,
                        factor = 2),
             shape = paste(geno, sex))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_hline(yintercept = .85,
             linetype = 'dashed') +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) + # c(1, 2, 16, 17)) +
  scale_x_continuous(breaks = c(1:5, 6.5)
                     , labels = c(format(concList[1]
                                         , nsmall = 0
                                         , scientific = FALSE)
                                  , format(concList[2]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[3]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[4]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[5]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[6]
                                           , nsmall = 0
                                           , scientific = FALSE))) +
  labs(y = 'Max accuracy (%)'
       , x = 'Odour Concentration (ppm)'
       , linetype = 'Genotype'
       , shape = 'Genotype'
       , colour = 'Sex') +
  theme(legend.position = c(.1, .2)
        , legend.title = element_blank()) +
  coord_cartesian(ylim = c(.2, 1))

plot85Pr6m <- bestBlockdat6m %>% 
  mutate('sex' = ifelse(sex == 'm'
                        , 'Male'
                        , 'Female')
         , 'geno' = ifelse(geno == 'tg'
                           , '5xFAD'
                           , 'B6SJL')) %>% 
  group_by(day, geno, sex) %>% 
  summarise(pr85 = mean(maxAcc >= .85)) %>% 
  ggplot(aes(x = day,
             y = pr85,
             shape = paste(geno, sex))) +
  geom_point(position = position_dodge(width = .5)) +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) +
  scale_x_continuous(breaks = c(1:5, 6.5)
                     , labels = c(format(concList[1]
                                         , nsmall = 0
                                         , scientific = FALSE)
                                  , format(concList[2]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[3]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[4]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[5]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[6]
                                           , nsmall = 0
                                           , scientific = FALSE))) +
  labs(y = 'Proportion reaching 85%'
       , x = 'Odour Concentration (ppm)'
       , linetype = 'Genotype'
       , shape = 'Genotype'
       , colour = 'Sex') +
  theme(legend.position = c(.2, .2)
        , legend.title = element_blank()) +
  coord_cartesian(ylim = c(0, 1))



#### sensitivity index last concentration ####
senDat6m <- dat6m %>% 
  filter(Consentration == 1e-5) %>% 
  group_by(Mouse, geno, sex) %>% 
  summarise('hits' = sum(Hits)
            , 'falseAlarms' = sum(10 - Correct_rejections)
            , 'hitRate' = ifelse(hits != 100
                                 , hits / 100
                                 , hits / (100 + 1))
            , 'faRate' = ifelse(falseAlarms != 0
                                , falseAlarms / 100
                                , (falseAlarms + 1) / 100)
            , 'faRate' = ifelse(faRate != 1
                                , faRate
                                , falseAlarms / 101)
            , 'deePrime' = qnorm(hitRate) - qnorm(faRate)
            , 'age' = '6m')

mfDeePrime6m <- lm(deePrime ~ geno * sex
                 , data = senDat6m)

options(na.action = "na.fail")
deePrimeTable6m <- dredge(mfDeePrime6m)
deePrimeTable6m

plotDeePrimBox6m <- senDat6m %>% 
  mutate('sex' = ifelse(sex == 'm'
                        , 'Male'
                        , 'Female')
         , 'Genotype' = ifelse(geno == 'tg'
                               , '5xFAD'
                               , 'B6SJL')) %>% 
  ggplot(aes(x = paste(Genotype, sex)
             , y = deePrime)) +
  geom_boxplot(fill = 'grey') +
  theme_classic() +
  labs(x = ''
       , y = expression(paste('Sensitivity Index (', italic('d\''), ')')))

#### compare d' of 12m and 6m on final conc ####

senDat12m <- senDat %>% 
  ungroup() %>% 
  mutate('Mouse' = animal
         , 'age' = '12m') %>% 
  select(-animal, -pTrials, -nTrials)

senDatAge <- bind_rows(senDat6m, senDat12m)

mfDeePrimeAge <- lm(deePrime ~ geno * sex * age
                   , data = senDatAge)

options(na.action = "na.fail")
deePrimeTableAge <- dredge(mfDeePrimeAge)
deePrimeTableAge

deePrimeAgeSelected <- lm(deePrime ~ age
                          , data = senDatAge)
deePrimeAgeNull <- lm(deePrime ~ 1
                          , data = senDatAge)
anova(deePrimeAgeSelected, deePrimeAgeNull)

aggregate(deePrime ~ age
          , data = senDatAge
          , FUN = mean)

aggregate(deePrime ~ age
          , data = senDatAge
          , FUN = sd)

senDatAge$age <- factor(senDatAge$age, c('6m', '12m'))

senDatAge %>% 
  mutate('sex' = ifelse(sex == 'm'
                        , 'Male'
                        , 'Female')
         , 'Genotype' = ifelse(geno == 'tg'
                               , '5xFAD'
                               , 'B6SJL')) %>% 
  ggplot(aes(x = paste(Genotype, sex)
             , y = deePrime)) +
  geom_boxplot(fill = 'grey') +
  theme_classic() +
  labs(x = ''
       , y = expression(paste('Sensitivity Index (', italic('d\''), ')'))) +
  facet_grid(cols = vars(age))

#### sensitivity all conc ####

senAllDat6m <- dat6m %>% 
  group_by(Mouse, Consentration, geno, sex) %>% 
  # filter(Block == 5) %>% 
  summarise('hits' = sum(Hits)
            , 'falseAlarms' = sum(10 - Correct_rejections)
            , 'hitRate' = ifelse(hits != 100
                                 , hits / 100
                                 , hits / (100 + 1))
            , 'faRate' = ifelse(falseAlarms != 0
                                , falseAlarms / 100
                                , (falseAlarms + 1) / 100)
            , 'faRate' = ifelse(faRate != 1
                                , faRate
                                , falseAlarms / 101)
            , 'deePrime' = qnorm(hitRate) - qnorm(faRate))

senAll6mmf <- lme(deePrime ~ geno * sex * Consentration
                , data = senAllDat6m
                , random = ~1|Mouse
                , method = 'ML')

options(na.action = "na.fail")
allDeePrimeTable6m <- dredge(senAll6mmf)
allDeePrimeTable6m  

allDeePrimeSelected6m <- lme(deePrime ~ sex
                             , data = senAllDat6m
                             , random = ~1|Mouse
                             , method = 'ML')

allDeePrimeNull6m <- lme(deePrime ~ 1
                             , data = senAllDat6m
                             , random = ~1|Mouse
                             , method = 'ML')

anova(allDeePrimeSelected6m, allDeePrimeNull6m)

aggregate(deePrime ~ sex
          , FUN = mean
          , data = senAllDat6m)
aggregate(deePrime ~ sex
          , FUN = sd
          , data = senAllDat6m)


sen6mDat <- dat6m %>% 
  mutate(day = -log10(Consentration) + 1,
         day = if_else(Block > 5,
                       7,
                       day)) %>% 
  group_by(Mouse, day, geno, sex) %>% 
  # filter(Block == 5) %>% 
  summarise('hits' = sum(Hits)
            , 'falseAlarms' = sum(10 - Correct_rejections)
            , 'hitRate' = ifelse(hits != 100
                                 , hits / 100
                                 , hits / (100 + 1))
            , 'faRate' = ifelse(falseAlarms != 0
                                , falseAlarms / 100
                                , (falseAlarms + 1) / 100)
            , 'faRate' = ifelse(faRate != 1
                                , faRate
                                , falseAlarms / 101)
            , 'deePrime' = qnorm(hitRate) - qnorm(faRate))

plotDeePrime6m <- sen6mDat %>% 
  ungroup() %>% 
  mutate('sex' = ifelse(sex == 'm'
                          , 'Male'
                          , 'Female')
         , 'Genotype' = ifelse(geno == 'tg'
                               , '5xFAD'
                               , 'B6SJL')) %>% 
  group_by(geno, sex, day) %>% 
  # summarise('meanDeePrime' = mean(deePrime)) %>% 
  ggplot(aes(x = day
             , y = deePrime
             , shape = paste(Genotype, sex))) +
  stat_summary(fun.data = 'mean_cl_boot'
               , position = position_dodge(width = .25)) +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) +
  scale_x_continuous(breaks = c(1:5, 6.5)
                     , labels = c(format(concList[1]
                                                  , nsmall = 0
                                                  , scientific = FALSE)
                                           , format(concList[2]
                                                    , nsmall = 0
                                                    , scientific = FALSE)
                                           , format(concList[3]
                                                    , nsmall = 0
                                                    , scientific = FALSE)
                                           , format(concList[4]
                                                    , nsmall = 0
                                                    , scientific = FALSE)
                                           , format(concList[5]
                                                    , nsmall = 0
                                                    , scientific = FALSE)
                                           , format(concList[6]
                                                    , nsmall = 0
                                                    , scientific = FALSE))) +
  labs(y = expression(paste('Sensitivity Index (', italic('d\''), ')'))
       , x = 'Odour Concentration (ppm)'
       , shape = '') +
  theme(legend.position = c(.8, .8)
        , legend.title = element_blank())

#### response bias ####

biasDat6m <- dat6m %>% 
  group_by(Mouse, geno, sex) %>% 
  summarise('hits' = sum(Hits)
            , 'falseAlarms' = sum(10 - Correct_rejections)
            , 'hitRate' = ifelse(hits != 350
                                 , hits / 350
                                 , hits / (350 + 1))
            , 'faRate' = ifelse(falseAlarms != 0
                                , falseAlarms / 350
                                , (falseAlarms + 1) / 350)
            , 'respCrit' = -.5 * (qnorm(hitRate) + qnorm(faRate))
            , 'age' = '6m')

mfRespCrit6m <- lm(respCrit ~ geno * sex
                 , data = biasDat6m)
respCritTable6m <- dredge(mfRespCrit6m)  
respCritTable6m

plotResponseBiasBox6m <- biasDat6m %>% 
  mutate('sex' = ifelse(sex == 'm'
                        , 'Male'
                        , 'Female')
         , 'Genotype' = ifelse(geno == 'tg'
                               , '5xFAD'
                               , 'B6SJL')) %>% 
  ggplot(aes(x = paste(Genotype, sex)
             , y = respCrit)) +
  geom_boxplot(fill = 'grey') +
  theme_classic() +
  labs(x = ''
       , y = expression(paste('Response Bias (', italic('c'), ')')))

dat6m %>% 
  group_by(Mouse, Consentration) %>% 
  summarise('hits' = sum(Hits)
            , 'n' = max(Block) * 20
            , 'falseAlarms' = (n/2) - sum(Correct_rejections)
            , 'hitRate' = ifelse(hits != (n/2)
                                 , hits / (n/2)
                                 , hits / ((n/2) + 1))
            , 'faRate' = ifelse(falseAlarms != 0
                                , falseAlarms / (n/2)
                                , (falseAlarms + 1) / (n/2))
            , faRate = if_else(falseAlarms == (n/2),
                               falseAlarms / (n/2 + 1),
                               faRate)
            , 'respCrit' = -.5 * (qnorm(hitRate) + qnorm(faRate))
           ) %>% 
  group_by(Consentration) %>% 
  summarise(mnBias = mean(respCrit))

#### compare bias 6m vs 12m ####

biasDat12m <- biasDat %>% 
  ungroup() %>% 
  mutate('Mouse' = animal
         , 'age' = '12m') %>% 
  select(-animal, -pTrials, -nTrials)

biasDatAll <- bind_rows(biasDat6m, biasDat12m)

mfRespCritAll <- lm(respCrit ~ geno * sex * age
                   , data = biasDatAll)
respCritTableAll <- dredge(mfRespCritAll)  
respCritTableAll

respCritAllSelected <- lm(respCrit ~ sex * age
                          , data = biasDatAll)
respCritAllNull <- lm(respCrit ~ 1
                          , data = biasDatAll)
anova(respCritAllSelected, respCritAllNull)

aggregate(respCrit ~ age
          , data = biasDatAll
          , FUN = mean)
aggregate(respCrit ~ age
          , data = biasDatAll
          , FUN = sd)

biasDatAll$age <- factor(biasDatAll$age,  c('6m', '12m'))

biasDatAll %>%
  mutate('sex' = ifelse(sex == 'm'
                        , 'Male'
                        , 'Female')
         , 'Genotype' = ifelse(geno == 'tg'
                               , '5xFAD'
                               , 'B6SJL')) %>% 
  ggplot(aes(x = paste(Genotype, sex)
             , y = respCrit)) +
  geom_boxplot(fill = 'grey') +
  theme_classic() +
  labs(x = ''
       , y = expression(paste('Response Bias (', italic('c'), ')'))) +
  facet_grid(cols = vars(age))

##### false alarm rate comparison ####

biasDatAll %>% 
  group_by(age) %>% 
  summarise(faMean = mean(faRate),
            faSD = sd(faRate),
            hitMean = mean(hitRate),
            hitSD = sd(hitRate),
            critMean = mean(respCrit),
            critSD = sd(respCrit))

datFa6m <- dat6m %>% 
  group_by(geno, sex, Consentration, Block) %>% 
  summarise(fa = mean_cl_boot((10 - Correct_rejections)/10)) %>% 
  mutate('xPos' = 5 * log10(1/Consentration) + Block + log10(1/Consentration)
                , 'sex' = ifelse(sex == 'm'
                                 , 'Male'
                                 , 'Female')
                , 'Genotype' = ifelse(geno == 'tg'
                                      , '5xFAD'
                                      , 'B6SJL'))

dat6mFa <- dat6m %>% 
  mutate(day = log10(1/Consentration) + 1,
         day = ifelse(Block > 5,
                      day + 1,
                      day),
         Block = if_else(day == 7,
                         Block - 5,
                         Block),
         fa = (10 - Correct_rejections) / 10 )

mfSMinus6m <- lme(fa ~ sex * Block * day * geno
                , random = ~1|Mouse
                , data = dat6mFa
                , method = 'ML')
dredge(mfSMinus6m)

mSelect <- lme(fa ~ Block + day + sex + Block:day 
                  , random = ~1|Mouse
                  , data = dat6mFa
                  , method = 'ML')
mNull <- lme(fa ~ 1
               , random = ~1|Mouse
               , data = dat6mFa
               , method = 'ML')
anova(mSelect, mNull)
AIC(mSelect, mNull)

dat6mFa %>% 
  group_by(sex) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))
dat6mFa %>% 
  group_by(day) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))
dat6mFa %>% 
  group_by(Block) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))

# day 1
dat6mFaD1 <- dat6mFa %>% 
  filter(day == 1)

mfSMinus6mD1 <- lme(fa ~ sex * Block * geno
                  , random = ~1|Mouse
                  , data = dat6mFaD1
                  , method = 'ML')
dredge(mfSMinus6mD1)

mSMinun6mD1Select <- lme(fa ~ Block
                         , random = ~1|Mouse
                         , data = dat6mFaD1
                         , method = 'ML')

mSMinun6mD1Null<- lme(fa ~ 1
                         , random = ~1|Mouse
                         , data = dat6mFaD1
                         , method = 'ML')

anova(mSMinun6mD1Select, mSMinun6mD1Null)
AICc(mSMinun6mD1Select, mSMinun6mD1Null)


dat6mFaD1 %>% 
  group_by(Block) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))

# day 2
dat6mFaD2 <- dat6mFa %>% 
  filter(day == 2)

mfSMinus6mD2 <- lme(fa ~ sex * Block * geno
                    , random = ~1|Mouse
                    , data = dat6mFaD2
                    , method = 'ML')
dredge(mfSMinus6mD2)

mSMinun6mD2Select <- lme(fa ~ Block * sex
                         , random = ~1|Mouse
                         , data = dat6mFaD2
                         , method = 'ML')

mSMinun6mD2Null<- lme(fa ~ 1
                      , random = ~1|Mouse
                      , data = dat6mFaD2
                      , method = 'ML')

anova(mSMinun6mD2Select, mSMinun6mD2Null)
AICc(mSMinun6mD2Select, mSMinun6mD2Null)


dat6mFaD2 %>% 
  group_by(Block) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))
dat6mFaD2 %>% 
  group_by(sex) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))

# day 3
dat6mFaD3 <- dat6mFa %>% 
  filter(day == 3)

mfSMinus6mD3 <- lme(fa ~ sex * Block * geno
                    , random = ~1|Mouse
                    , data = dat6mFaD3
                    , method = 'ML')
dredge(mfSMinus6mD3)

mSMinun6mD3Select <- lme(fa ~ (Block + geno + sex)^2 - geno:sex
                         , random = ~1|Mouse
                         , data = dat6mFaD3
                         , method = 'ML')

mSMinun6mD3Null<- lme(fa ~ 1
                      , random = ~1|Mouse
                      , data = dat6mFaD3
                      , method = 'ML')

anova(mSMinun6mD3Select, mSMinun6mD3Null)
AICc(mSMinun6mD3Select, mSMinun6mD3Null)

dat6mFaD3 %>% 
  group_by(Block) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))
dat6mFaD3 %>% 
  group_by(geno) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))
dat6mFaD3 %>% 
  group_by(sex) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))

# day 4
dat6mFaD4 <- dat6mFa %>% 
  filter(day == 4)

mfSMinus6mD4 <- lme(fa ~ sex * Block * geno
                    , random = ~1|Mouse
                    , data = dat6mFaD4
                    , method = 'ML')
dredge(mfSMinus6mD4)

mSMinun6mD4Select <- lme(fa ~ Block + geno
                         , random = ~1|Mouse
                         , data = dat6mFaD4
                         , method = 'ML')

mSMinun6mD4Null<- lme(fa ~ 1
                      , random = ~1|Mouse
                      , data = dat6mFaD4
                      , method = 'ML')

anova(mSMinun6mD4Select, mSMinun6mD4Null)
AICc(mSMinun6mD4Select, mSMinun6mD4Null)


dat6mFaD4 %>% 
  group_by(Block) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))

dat6mFaD4 %>% 
  group_by(geno) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))


# day 5
dat6mFaD5 <- dat6mFa %>% 
  filter(day == 5)

mfSMinus6mD5 <- lme(fa ~ sex * Block * geno
                    , random = ~1|Mouse
                    , data = dat6mFaD5
                    , method = 'ML')
dredge(mfSMinus6mD5)

mSMinun6mD5Select <- lme(fa ~ Block
                         , random = ~1|Mouse
                         , data = dat6mFaD5
                         , method = 'ML')

mSMinun6mD5Null<- lme(fa ~ 1
                      , random = ~1|Mouse
                      , data = dat6mFaD5
                      , method = 'ML')

anova(mSMinun6mD5Select, mSMinun6mD5Null)
AICc(mSMinun6mD5Select, mSMinun6mD5Null)

dat6mFaD5 %>% 
  group_by(Block) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))

# day 6
dat6mFaD6 <- dat6mFa %>% 
  filter(day == 6)

mfSMinus6mD6 <- lme(fa ~ sex * Block * geno
                    , random = ~1|Mouse
                    , data = dat6mFaD6
                    , method = 'ML')
dredge(mfSMinus6mD6)

mSMinun6mD6Select <- lme(fa ~ Block
                         , random = ~1|Mouse
                         , data = dat6mFaD6
                         , method = 'ML')

mSMinun6mD6Null<- lme(fa ~ 1
                      , random = ~1|Mouse
                      , data = dat6mFaD6
                      , method = 'ML')

anova(mSMinun6mD6Select, mSMinun6mD6Null)
AICc(mSMinun6mD6Select, mSMinun6mD6Null)


# day 7
dat6mFaD7 <- dat6mFa %>% 
  filter(day == 7)

mfSMinus6mD7 <- lme(fa ~ sex * Block * geno
                    , random = ~1|Mouse
                    , data = dat6mFaD7
                    , method = 'ML')
dredge(mfSMinus6mD7)

mSMinun6mD7Select <- lme(fa ~ Block + geno
                         , random = ~1|Mouse
                         , data = dat6mFaD7
                         , method = 'ML')

mSMinun6mD7Null<- lme(fa ~ 1
                      , random = ~1|Mouse
                      , data = dat6mFaD7
                      , method = 'ML')

anova(mSMinun6mD7Select, mSMinun6mD7Null)
AICc(mSMinun6mD7Select, mSMinun6mD7Null)

dat6mFaD7 %>% 
  group_by(Block) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))

dat6mFaD7 %>% 
  group_by(geno) %>% 
  summarise(mn = mean(fa),
            stdDv = sd(fa))

##### False alarm plot ####

plotFa6m <- ggplot(data = filter(datFa6m
                                  , Consentration == 1e0)
                    , aes(x = xPos
                          , y = fa$y
                          , ymin = fa$ymin
                          , ymax - fa$ymax
                          , shape = paste(Genotype, sex)
                          , linetype = paste(Genotype, sex))) +
  geom_point(size = pointSize) +
  geom_line() +
  geom_errorbar(data = datFa6m
                , aes(x = xPos
                      , y = fa$y
                      , ymin = fa$ymin
                      , ymax = fa$ymax)
                , linetype = 1
                , width = .5) +
  geom_point(data = filter(datFa6m
                           , Consentration == 1e-1)
             , size = pointSize) +
  geom_line(data = filter(datFa6m
                          , Consentration == 1e-1)) +
  geom_point(data = filter(datFa6m
                           , Consentration == 1e-2)
             , size = pointSize) +
  geom_line(data = filter(datFa6m
                          , Consentration == 1e-2)) +
  geom_point(data = filter(datFa6m
                           , Consentration == 1e-3)
             , size = pointSize) +
  geom_line(data = filter(datFa6m
                          , Consentration == 1e-3)) +
  geom_point(data = filter(datFa6m
                           , Consentration == 1e-4)
             , size = pointSize) +
  geom_line(data = filter(datFa6m
                          , Consentration == 1e-4)) +
  geom_point(data = filter(datFa6m
                           , Consentration == 1e-5
                           , Block < 6)
             , size = pointSize) +
  geom_line(data = filter(datFa6m
                          , Consentration == 1e-5
                          , Block < 6)) +
  geom_point(data = filter(datFa6m
                           , Consentration == 1e-5
                           , Block >= 6)
             , size = pointSize) +
  geom_line(data = filter(datFa6m
                          , Consentration == 1e-5
                          , Block >= 6)) +
  # geom_hline(aes(yintercept = .85)
  #            , linetype = 'dotted') +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) + #c(1, 2, 16, 17)) +
  scale_linetype_manual(values = c('dashed', 'dashed', 'solid', 'solid')) +
  # scale_y_continuous(breaks = seq(0, 1, by = .2)) +
  scale_x_continuous(breaks = axisBreaks
                     , labels = c(format(concList[1]
                                         , nsmall = 0
                                         , scientific = FALSE)
                                  , format(concList[2]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[3]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[4]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[5]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[6]
                                           , nsmall = 0
                                           , scientific = FALSE))) +
  labs(y = 'False alarm rate (%)'
       , x = 'Odour Concentration (ppm)'
       , linetype = 'Genotype'
       , shape = 'Genotype'
       , colour = 'Sex') +
  theme(legend.position = c(.2, .85)
        , legend.title = element_blank()
        , legend.background = element_blank())




#### bias all concentrations ####

biasConcDat6m <- dat6m %>% 
  group_by(Mouse, Consentration, geno, sex) %>% 
  summarise('hits' = sum(Hits)
            , 'falseAlarms' = sum(10 - Correct_rejections)
            , 'hitRate' = ifelse(hits != 350
                                 , hits / 350
                                 , hits / (350 + 1))
            , 'faRate' = ifelse(falseAlarms != 0
                                , falseAlarms / 350
                                , (falseAlarms + 1) / 350)
            , 'respCrit' = -.5 * (qnorm(hitRate) + qnorm(faRate))
            , 'age' = '6m')

mfBiasConc6m <- lme(respCrit ~ sex * geno * Consentration
    , random = ~1|Mouse
    , data = biasConcDat6m
    , method = 'ML')


dat6m %>% 
  mutate(day = -log10(Consentration) + 1,
         day = if_else(Block > 5,
                       7,
                       day)) %>% 
  group_by(Mouse, Day, geno, sex) %>% 
  summarise('hits' = sum(Hits)
            , 'falseAlarms' = sum(10 - Correct_rejections)
            , 'hitRate' = ifelse(hits != 350
                                 , hits / 350
                                 , hits / (350 + 1))
            , 'faRate' = ifelse(falseAlarms != 0
                                , falseAlarms / 350
                                , (falseAlarms + 1) / 350)
            , 'respCrit' = -.5 * (qnorm(hitRate) + qnorm(faRate))
            , 'age' = '6m')

plotResponseBiasConcBox6m <- dat6m %>% 
  mutate(day = -log10(Consentration) + 1,
         day = if_else(Block > 5,
                       7,
                       day)) %>% 
  group_by(Mouse, day, geno, sex) %>% 
  summarise('hits' = sum(Hits)
            , 'falseAlarms' = sum(10 - Correct_rejections)
            , 'hitRate' = ifelse(hits != 50
                                 , hits / 50
                                 , hits / (50 + 1))
            , 'faRate' = ifelse(falseAlarms != 0
                                , falseAlarms / 50
                                , (falseAlarms + 1) / 50)
            , 'respCrit' = -.5 * (qnorm(hitRate) + qnorm(faRate))
            , 'age' = '6m') %>%  
  ungroup() %>% 
  mutate('sex' = ifelse(sex == 'm'
                          , 'Male'
                          , 'Female')
         , 'Genotype' = ifelse(geno == 'tg'
                               , '5xFAD'
                               , 'B6SJL')) %>% 
  group_by(geno, sex, day) %>% 
  # summarise('meanDeePrime' = mean(deePrime)) %>% 
  ggplot(aes(x = day - 1
             , y = respCrit
             , shape = paste(Genotype, sex))) +
  stat_summary(fun.data = 'mean_cl_boot'
               , position = position_dodge(width = .25)) +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 2, 1)) +
  scale_x_continuous(breaks = c(0:4, 5.5)
                     , labels = c(format(concList[1]
                                         , nsmall = 0
                                         , scientific = FALSE)
                                  , format(concList[2]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[3]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[4]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[5]
                                           , nsmall = 0
                                           , scientific = FALSE)
                                  , format(concList[6]
                                           , nsmall = 0
                                           , scientific = FALSE))) +
  labs(y = expression(paste('Response Bias (', italic('c\''), ')'))
       , x = 'Odour Concentration (ppm)'
       , shape = '') +
  theme(legend.position = c(.8, .2)
        , legend.title = element_blank())


biasConc6mTable <- dredge(mfBiasConc6m)
biasConc6mTable

mBiasConc6mSelected <- lme(respCrit ~ Consentration
                    , random = ~1|Mouse
                    , data = biasConcDat6m
                    , method = 'ML')

mBiasConc6mNull <- lme(respCrit ~ 1
                           , random = ~1|Mouse
                           , data = biasConcDat6m
                           , method = 'ML')

anova(mBiasConc6mSelected, mBiasConc6mNull)

#### learning ####

dat6l <- dat6m %>% 
  mutate('acc' = (Hits + Correct_rejections)/20)

mf6mAcc <- lme(acc ~ sex * geno * Block * Consentration
             , random = ~1|Mouse
             , data = dat6l
             , method = 'ML')

acc6mTable <- dredge(mf6mAcc)
acc6mTable
# acc6mTable[acc6mTable$AICc == AICc(mNull6mAcc),]

mSelected6mAcc <- lme(acc ~ sex + Block * Consentration
                     , random = ~1|Mouse
                     , data = dat6l
                     , method = 'ML')

mNull6mAcc <- lme(acc ~ 1
                      , random = ~1|Mouse
                      , data = dat6l
                      , method = 'ML')

AICc(mSelected6mAcc, mNull6mAcc)
anova(mSelected6mAcc, mNull6mAcc)

aggregate(acc ~ sex
          , data = dat6l
          , FUN = mean)
aggregate(acc ~ sex
          , data = dat6l
          , FUN = sd)

aggregate(acc ~ Block
          , data = dat6l
          , FUN = mean)
aggregate(acc ~ Block
          , data = dat6l
          , FUN = sd)

aggregate(acc ~ Consentration
          , data = dat6l
          , FUN = mean)
aggregate(acc ~ Consentration
          , data = dat6l
          , FUN = sd)

aggregate(acc ~ Block
          , data = dat6l[dat6l$Consentration == 1e-05,]
          , FUN = mean)
aggregate(acc ~ Block
          , data = dat6l[dat6l$Consentration == 1e-05,]
          , FUN = sd)

#### acc 1st conc ####

acc1Dat6m <- dat6l %>% 
  filter(Consentration == concList[1])

mfAcc1.6m <- lme(acc ~ sex * geno * Block
                 , random = ~1|Mouse
                 , data = acc1Dat6m
                 , method = 'ML')

acc1.6mTable <- dredge(mfAcc1.6m)
acc1.6mTable

mAcc1.6mSelected <- lme(acc ~ sex + Block
                 , random = ~1|Mouse
                 , data = acc1Dat6m
                 , method = 'ML')

mAcc1.6mNull <- lme(acc ~ 1
                        , random = ~1|Mouse
                        , data = acc1Dat6m
                        , method = 'ML')

anova(mAcc1.6mSelected, mAcc1.6mNull)

aggregate(acc ~ sex
          , FUN = mean
          , data = acc1Dat6m)
aggregate(acc ~ sex
          , FUN = sd
          , data = acc1Dat6m)

#### acc 2nd conc ####

acc2Dat6m <- dat6l %>% 
  filter(Consentration == concList[2])

mfAcc2.6m <- lme(acc ~ sex * geno * Block
                 , random = ~1|Mouse
                 , data = acc2Dat6m
                 , method = 'ML')

acc2.6mTable <- dredge(mfAcc2.6m)
acc2.6mTable

mAcc2.6mSelected <- lme(acc ~ sex * Block
                 , random = ~1|Mouse
                 , data = acc2Dat6m
                 , method = 'ML')

mAcc2.6mNull <- lme(acc ~ 1
                        , random = ~1|Mouse
                        , data = acc2Dat6m
                        , method = 'ML')

anova(mAcc2.6mSelected, mAcc2.6mNull)

aggregate(acc ~ sex
          , FUN = mean
          , data = acc2Dat6m)
aggregate(acc ~ sex
          , FUN = sd
          , data = acc2Dat6m)

aggregate(acc ~ Block
          , FUN = mean
          , data = acc2Dat6m)
aggregate(acc ~ Block
          , FUN = sd
          , data = acc2Dat6m)

aggregate(acc ~ Block * sex
          , FUN = mean
          , data = acc2Dat6m)
aggregate(acc ~ Block * sex
          , FUN = sd
          , data = acc2Dat6m)


#### acc 3rd conc ####

acc3Dat6m <- dat6l %>% 
  filter(Consentration == concList[3])

mfAcc3.6m <- lme(acc ~ sex * geno * Block
                 , random = ~1|Mouse
                 , data = acc3Dat6m
                 , method = 'ML')

acc3.6mTable <- dredge(mfAcc3.6m)
acc3.6mTable

mAcc3.6mSelected <- lme(acc ~ sex + geno * Block
                 , random = ~1|Mouse
                 , data = acc3Dat6m
                 , method = 'ML')

mAcc3.6mNull <- lme(acc ~ 1
                        , random = ~1|Mouse
                        , data = acc3Dat6m
                        , method = 'ML')

anova(mAcc3.6mSelected, mAcc3.6mNull)

aggregate(acc ~ geno
          , FUN = mean
          , data = acc3Dat6m)
aggregate(acc ~ geno
          , FUN = sd
          , data = acc3Dat6m)

aggregate(acc ~ Block * geno
          , FUN = mean
          , data = acc3Dat6m)
aggregate(acc ~ Block * geno
          , FUN = sd
          , data = acc3Dat6m)

aggregate(acc ~ sex
          , FUN = mean
          , data = acc3Dat6m)
aggregate(acc ~ sex
          , FUN = sd
          , data = acc3Dat6m)

#### acc 4th conc ####

acc4Dat6m <- dat6l %>% 
  filter(Consentration == concList[4])

mfAcc4.6m <- lme(acc ~ sex * geno * Block
                 , random = ~1|Mouse
                 , data = acc4Dat6m
                 , method = 'ML')

acc4.6mTable <- dredge(mfAcc4.6m)
acc4.6mTable

mAcc4.6mSelected <- lme(acc ~ geno + Block
                 , random = ~1|Mouse
                 , data = acc4Dat6m
                 , method = 'ML')

mAcc4.6mNull <- lme(acc ~ 1
                        , random = ~1|Mouse
                        , data = acc4Dat6m
                        , method = 'ML')

anova(mAcc4.6mSelected, mAcc4.6mNull)

aggregate(acc ~ geno
          , FUN = mean
          , data = acc4Dat6m)
aggregate(acc ~ geno
          , FUN = sd
          , data = acc4Dat6m)

aggregate(acc ~ Block
          , FUN = mean
          , data = acc4Dat6m)
aggregate(acc ~ Block
          , FUN = sd
          , data = acc4Dat6m)

#### acc 5th conc ####

acc5Dat6m <- dat6l %>% 
  filter(Consentration == concList[5])

mfAcc5.6m <- lme(acc ~ sex * geno * Block
                 , random = ~1|Mouse
                 , data = acc5Dat6m
                 , method = 'ML')

acc5.6mTable <- dredge(mfAcc5.6m)
acc5.6mTable

mAcc5.6mSelected <- lme(acc ~ Block
                 , random = ~1|Mouse
                 , data = acc5Dat6m
                 , method = 'ML')

mAcc5.6mNull <- lme(acc ~ 1
                        , random = ~1|Mouse
                        , data = acc5Dat6m
                        , method = 'ML')

anova(mAcc5.6mSelected, mAcc5.6mNull)

aggregate(acc ~ Block
          , FUN = mean
          , data = acc5Dat6m)
aggregate(acc ~ Block
          , FUN = sd
          , data = acc5Dat6m)

#### acc lowest conc day 1 ####

acc6Dat6m <- dat6l %>% 
  filter(Consentration == concList[6]
         , Block < 6)

mfAcc6.6m <- lme(acc ~ sex * geno * Block
                 , random = ~1|Mouse
                 , data = acc6Dat6m
                 , method = 'ML')

acc6.6mTable <- dredge(mfAcc6.6m)
acc6.6mTable

mAcc6.6mSelected <- lme(acc ~ Block
                 , random = ~1|Mouse
                 , data = acc6Dat6m
                 , method = 'ML')

mAcc6.6mNull <- lme(acc ~ 1
                        , random = ~1|Mouse
                        , data = acc6Dat6m
                        , method = 'ML')

anova(mAcc6.6mSelected, mAcc6.6mNull)

aggregate(acc ~ Block
          , FUN = mean
          , data = acc6Dat6m)
aggregate(acc ~ Block
          , FUN = sd
          , data = acc6Dat6m)

#### acc lowest conc day 2 ####

acc7Dat6m <- dat6l %>% 
  filter(Consentration == concList[6]
         , Block  >5)

mfAcc7.6m <- lme(acc ~ sex * geno * Block
                 , random = ~1|Mouse
                 , data = acc7Dat6m
                 , method = 'ML')

acc7.6mTable <- dredge(mfAcc7.6m)
acc7.6mTable

mAcc7.6mSelected <- lme(acc ~ geno + Block
                        , random = ~1|Mouse
                        , data = acc7Dat6m
                        , method = 'ML')

mAcc7.6mNull <- lme(acc ~ 1
                    , random = ~1|Mouse
                    , data = acc7Dat6m
                    , method = 'ML')

anova(mAcc7.6mSelected, mAcc7.6mNull)

aggregate(acc ~ geno
          , FUN = mean
          , data = acc7Dat6m)
aggregate(acc ~ geno
          , FUN = sd
          , data = acc7Dat6m)

aggregate(acc ~ Block
          , FUN = mean
          , data = acc7Dat6m)
aggregate(acc ~ Block
          , FUN = sd
          , data = acc7Dat6m)
