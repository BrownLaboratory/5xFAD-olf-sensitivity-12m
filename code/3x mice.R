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

options(na.action = "na.fail")

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
dat3xIntro <- dat %>%
    filter(
        strain == "3x",
        StudyName == "INTRO",
        SS == FALSE
    ) %>%
    group_by(animal, Podour, block, geno, sex)

# OP1 data
dat3xOP1 <- dat %>%
    filter(
        strain == "3x",
        StudyName == "OP1",
        SS == FALSE
    ) %>%
    group_by(animal, Podour, block, geno, sex)

# rev data
dat3xRev <- dat %>%
    filter(
        strain == "3x",
        StudyName == "REV1",
        SS == FALSE
    ) %>%
    group_by(animal, Podour, block, geno, sex)

# create list of concentrations
concList <- 10^-(c(0:5, 5))

# get 3x sensitivity data, drop training and short sample trials
dat3x <- dat %>%
    filter(
        strain == "3x",
        StudyName %in% paste0("SEN", 1:7),
        SS == FALSE,
        !(StudyName == "SEN5" & block > 5)
    ) %>%
    group_by(animal, StudyName, Podour, block, geno, sex)

# summarize 3x data
dat3xSum <- dat3x %>%
    summarise(
        "acc" = sum(correct) / length(correct),
        "crit" = acc >= .85
    )
#### Intro analysis ####

# number of errors
introErrors <- dat3xIntro %>%
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

introErrorsTable <- dredge(mfIntroErrors)

introErrorsTable

#### OP1 analysis ####

# number of errors
op1Errors <- dat3xOP1 %>%
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

op1ErrorsTable <- dredge(mfOp1Errors)
op1ErrorsTable

#### Rev analysis ####

# number of errors
revErrors <- dat3xRev %>%
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

revErrorsTable <- dredge(mfRevErrors)
revErrorsTable

mSelectRevErrors <- get.models(revErrorsTable,
                               delta == 0)[[1]]

anova(mSelectRevErrors, mNullRevErrors)

revErrors %>%
    group_by(sex) %>%
    summarise(mean_errors = mean(errors),
              sd_errors = sd(errors))

#### sensitivity ####
#### model mice hiting crit ####

# did mice reach crit at least once per odour
dat3xCrit <- dat3xSum %>%
    group_by(animal, Podour, geno, sex) %>%
    summarise("Crit" = any(crit)) %>%
    mutate(
        "nChar" = nchar(as.character(Podour)),
        "conc" = as.numeric(substr(Podour, 1, nChar - 6))
    )

# model using conc as num
mfCrit <- glmer(Crit ~ conc * geno * sex + (1 | animal),
                data = dat3xCrit,
                family = binomial
)

mNullCrit <- glmer(Crit ~ (1 | animal),
                   data = dat3xCrit,
                   family = binomial
)

critTable <- dredge(mfCrit)
critTable

mSelectedCrit <- glmer(Crit ~ conc * sex + (1 | animal),
                       data = dat3xCrit,
                       family = binomial
)

anova(mSelectedCrit, mNullCrit, test = "Chisq")
summary(mSelectedCrit)
# confSelectedCrit <- confint(mSelectedCrit)
# confSelectedCrit




#### model Accuracy ####
dat3xConc <- dat3x %>%
    ungroup() %>%
    mutate(
        "nChar" = nchar(as.character(Podour)),
        "conc" = as.numeric(substr(Podour, 1, nChar - 6)),
        "block" = ifelse(StudyName == "SEN7",
                         block + 5,
                         block
        )
    )


acc3xdat <- dat3xConc %>%
    group_by(animal, sex, geno, block, conc) %>%
    # filter(conc == concList[7]
    #        , block == 10) %>%
    summarise("acc" = mean(correct)) %>%
    mutate("concLog" = log10(conc))

mfAcc <- lme(acc ~ sex * geno * block * conc,
             random = ~ 1 | animal,
             data = acc3xdat,
             method = "ML"
)

accTable <- dredge(mfAcc)
accTable

mAccSelected <- get.models(accTable, delta == 0)[[1]]
summary(mAccSelected)
intervals(mAccSelected)

mNullAcc <- lme(acc ~ 1,
                random = ~ 1 | animal,
                data = acc3xdat,
                method = "ML"
)

anova(mAccSelected, mNullAcc)
AICc(mNullAcc)



#### plot accuracy ####
datCorrPlot <- dat3xConc %>%
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
                            "3xTg-AD",
                            "B6129"
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
