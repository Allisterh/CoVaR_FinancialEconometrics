rm(list=ls())
setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 7")

library(rugarch)
library(Rsolnp)
library(copula)


################################################################################
### Problem 1                                                                ###
################################################################################

source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/PattonCopula.R")

################################################################################
### Problem 2                                                                ###
################################################################################

### Point 1)
# Consider the same couple of assets you choose in the previous exercise set

# load dataset from rugarch package
data(dji30ret)

# save returns in vector (ticker AA and AXP first 2000 observations)
dj_returns <- dji30ret[seq(1, 2000), c("AA", "AXP")]

# Estimate a the Patton’s model with Gaussian copula.
copNorm <- Estimate_Patton(CopType = "norm", mY = dj_returns, M = 10)

# Estimate a the Patton’s model with Student’s t copula.
copT <- Estimate_Patton(CopType = "t", mY = dj_returns, M = 10)

# Compare the filtered copula correlation parameter of the Gaussian and Student’s t copula models.
if(!require(ggplot2)){install.packages('ggplot2')}
if(!require(ggthemes)){install.packages('ggthemes')}
df <- data.frame(
    time = as.Date(row.names(dj_returns)),
    corr = c(copNorm$vCor,
            copT$vCor),
    lab = rep(c("Gaussian","Student's t"),
            each = nrow(dj_returns))
    )
# Plot
ggplot(data = df, aes(x=time, y=corr, col = lab)) + 
                                    geom_line() +
                                    labs(title = paste0("Correlation of Bivariate Elliptical Random Variable")) +
                                    scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month",
                                    date_labels = "%m-%Y") + 
                                    labs(y = "Correlation") +
                                    theme_economist() +
                                    theme(legend.title=element_blank())
setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 7")
ggsave("./img/Correlation of Gaussian Copula and t Copula.pdf")

# The correlation of the Gaussian copula is more stable than the correlation of the Student’s t copula.

# Compare the two copula models using BIC. Which model is selected?
sort(c(BIC_Gaussian = copNorm$BIC, BIC_t = copT$BIC)) # lower BIC values are preferred
# The Student’s t copula model is selected.
