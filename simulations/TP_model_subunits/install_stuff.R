pcklist <- c(
    "knitr",
    "latex2exp",
    "pryr",
    "gdata", 
    "reshape2",
    "ggplot2",
    "tidyr", 
    "plyr",  
    "Hmisc", 
    "tidyverse",
    "latex2exp", 
    "FSA",
    "lme4",
    "nlme",
    "lmerTest",
    "ez",
    "multcomp",
    "languageR")

lapply (pcklist,
        install.packages,
        dependencies = TRUE)
