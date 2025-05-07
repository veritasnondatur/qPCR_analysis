#### Scatterplot and correlation analysis of Nalm6 with Pbx1a OE (pCLIG) ###


# BLK #

# Import dataset and construct matrix
library(readxl)
data <- read_excel("qPCR_Nalm6_Pbx1aOE_allruns_analysis.xlsx", 
                   sheet = "data_for_R", col_types = c("text", 
                                                       "numeric", "numeric"))
View(data)


# Assign variables for data to plot
scatterplot <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
PBX1 <- data$hPBX1
BLK <- data$hBLK


# Create plot for TALE-HD
library(ggplot2)
plot(scatterplot, PBX1, type = "b", col = "purple", pch = 8,                      # Add PBX1 values
     main = "Expression of PBX1 and BLK in Nalm6 with(out) Pbx1a OE",    # Change main title & axis labels
     ylab = "Expression rel. to RPLP0 [ddCt]",
     ylim = c(0,0.5),
     xlim = c(0.9, 16))
lines(scatterplot, BLK, type = "b",  col = "darkgreen", pch = 15)            
             

legend("topright", #inset=c(-0.2,0),                                       # Add legend to plot
       legend = c("PBX1", 
                  "BLK" 
                  ),
       col = c("purple", 
               "darkgreen", 
               "black", 
               "grey", 
               "brown"),
       pch = c(8, 
               15
               ))


# Compute correlation coefficient
library("ggpubr")
cor(PBX1, BLK, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, BLK, method=c("pearson"))

ggscatter(data, x = "hPBX1", y = "hBLK", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to RPLP0 [ddCt]", ylab = "BLK Expression rel. to RPLP0 [ddCt]")

############################################################################

# More plots/analysis #

# Import dataset and construct matrix
library(readxl)
data2 <- read_excel("qPCR_Nalm6_Pbx1aOE_allruns_analysis.xlsx", 
                    sheet = "data_for_R_moreTargets", col_types = c("text", 
                                                                    "numeric", "numeric", "numeric", 
                                                                    "numeric", "numeric", "numeric", 
                                                                    "numeric"))
na.omit(data2)
View(data2)


# Assign variables for data to plot
PBX1 <- data2$hPBX1
NOTCH3 <- data2$hNOTCH3
TENM4<- data2$hTENM4
E2F5 <- data2$hE2F5
HPRT1 <- data2$hHPRT1
ETV5 <- data2$hETV5

# Compute correlation coefficient for NOTCH3
library("ggpubr")
cor(PBX1, NOTCH3, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, NOTCH3, method=c("pearson"))

ggscatter(data2, x = "hPBX1", y = "hNOTCH3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to RPLP0 [ddCt]", ylab = "NOTCH3 Expression rel. to RPLP0 [ddCt]")


# Compute correlation coefficient for TENM4
library("ggpubr")
cor(PBX1, TENM4, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, TENM4, method=c("pearson"))

ggscatter(data2, x = "hPBX1", y = "hTENM4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to RPLP0 [ddCt]", ylab = "TENM4 Expression rel. to RPLP0 [ddCt]")


# Compute correlation coefficient for E2F5
library("ggpubr")
cor(PBX1, E2F5, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, E2F5, method=c("pearson"))

ggscatter(data2, x = "hPBX1", y = "hE2F5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to RPLP0 [ddCt]", ylab = "E2F5 Expression rel. to RPLP0 [ddCt]")


# Compute correlation coefficient for HPRT1
library("ggpubr")
cor(PBX1, HPRT1, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, HPRT1, method=c("pearson"))

ggscatter(data2, x = "hPBX1", y = "hHPRT1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to RPLP0 [ddCt]", ylab = "HPRT1 Expression rel. to RPLP0 [ddCt]")


# Compute correlation coefficient for ETV5
library("ggpubr")
cor(PBX1, ETV5, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, ETV5, method=c("pearson"))

ggscatter(data2, x = "hPBX1", y = "hETV5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to RPLP0 [ddCt]", ylab = "ETV5 Expression rel. to RPLP0 [ddCt]")


############################################################################

# More plots/analysis, relative to HPRT1 instead of RPLP0 #

# Import dataset and construct matrix
library(readxl)
data3 <- read_excel("qPCR_Nalm6_Pbx1aOE_allruns_analysis.xlsx", 
                    sheet = "combined_relHPRT1", col_types = c("text", 
                                                               "numeric", "numeric", "numeric", 
                                                               "numeric", "numeric", "numeric"))
na.omit(data3)
View(data3)


# Assign variables for data to plot
PBX1 <- data3$hPBX1
NOTCH3 <- data3$hNOTCH3
TENM4<- data3$hTENM4
E2F5 <- data3$hE2F5
RPLP0 <- data3$hRPLP0
ETV5 <- data3$hETV5

# Compute correlation coefficient for NOTCH3
library("ggpubr")
cor(PBX1, NOTCH3, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, NOTCH3, method=c("pearson"))

ggscatter(data3, x = "hPBX1", y = "hNOTCH3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to HPRT1 [ddCt]", ylab = "NOTCH3 Expression rel. to HPRT1 [ddCt]")


# Compute correlation coefficient for TENM4
library("ggpubr")
cor(PBX1, TENM4, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, TENM4, method=c("pearson"))

ggscatter(data3, x = "hPBX1", y = "hTENM4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to HPRT1 [ddCt]", ylab = "TENM4 Expression rel. to HPRT1 [ddCt]")


# Compute correlation coefficient for E2F5
library("ggpubr")
cor(PBX1, E2F5, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, E2F5, method=c("pearson"))

ggscatter(data3, x = "hPBX1", y = "hE2F5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to HPRT1 [ddCt]", ylab = "E2F5 Expression rel. to HPRT1 [ddCt]")


# Compute correlation coefficient for RPLP0
library("ggpubr")
cor(PBX1, RPLP0, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, RPLP0, method=c("pearson"))

ggscatter(data3, x = "hPBX1", y = "hRPLP0", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to HPRT1 [ddCt]", ylab = "RPLP0 Expression rel. to HPRT1 [ddCt]")


# Compute correlation coefficient for ETV5
library("ggpubr")
cor(PBX1, ETV5, method = c("pearson"))    # use pearson correlation, as this does not assume any kind of distribution (like normal dist) of data
cor.test(PBX1, ETV5, method=c("pearson"))

ggscatter(data3, x = "hPBX1", y = "hETV5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          main = "Pearson Correlation",
          xlab = "PBX1 Expression rel. to HPRT1 [ddCt]", ylab = "ETV5 Expression rel. to HPRT1 [ddCt]")
