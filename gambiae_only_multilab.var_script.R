#How quantified variability was determined from WHO cone bioassay data  

#Imput your directory and dataset
setwd("")
Vesterg.multilab.dat <- read.csv("Anonymised_Multilab_An.gambiae_only_dataset.csv")


#load packages 
library(ggplot2)
library(lme4)
library(lmtest)
library(lsmeans)
library(beeswarm)
library(ggbeeswarm)

#first line is code to check for rows with % mortality values >100 as quality check and second is code to remove any rows with moralities >100%
Vesterg.multilab.dat  <- subset(Vesterg.multilab.dat , perc_Mortality <= 100) 

#Trying to work out the different products (LLIN codes)
table(Vesterg.multilab.dat $LLIN_code)

#different sites
table(Vesterg.multilab.dat $Testing.lab_coded)


df <-data.frame(
  Product=Vesterg.multilab.dat$LLIN_code,
  Test_lab = Vesterg.multilab.dat$Testing.lab_coded
)

table(df)


###############################################################################
#observe spread of data - such as violin plots. geom_points adds raw data points to the plot

lab_colors <- c("Lab A" = "red", "Lab B" = "blue", "Lab C" = "green", "Lab D"= "magenta", 
                "Lab E"= "grey", "Lab F"="cyan", "Lab I"="yellow") 

ggplot(Vesterg.multilab.dat, aes(x = Testing.lab_coded, y = perc_Mortality, color = Testing.lab_coded)) + 
  geom_violin() +
  geom_beeswarm(color = "black", cex = 1.1) +
  labs(title = "Mortality variation between LLINs and Labs",
       x = "Test laboratory",
       y = "Percentage mortality") +
  ylim(0, 100) +
  facet_wrap(~ LLIN_code, nrow = 2) +  # Wrap into two rows
  scale_color_manual(name = "Test Laboratory", values = lab_colors) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic() +  # Removes grid lines and gives a white background
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




##################################################################################

#calculating raw dead numbers from proportion for binomal glmms

#divide by 100 to get proportion not a percentage
Vesterg.multilab.dat  $proportionalmortality <-Vesterg.multilab.dat $perc_Mortality/100
#then multiple totmosquitoes by the proportion
Vesterg.multilab.dat  $Number.of.dead <-Vesterg.multilab.dat $Totmosquitoes * Vesterg.multilab.dat $proportionalmortality

#making sure number of deads are a integer
class(Vesterg.multilab.dat$Number.of.dead)
Vesterg.multilab.dat$Number.of.dead <- as.integer(round(Vesterg.multilab.dat$Number.of.dead))

any(!is.integer(Vesterg.multilab.dat$Number.of.dead))

#bin glmms


#with just test lab as a fixed effect
#Vesterg.multilab.glmm.TL <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~  Testing.lab_coded + (1|LLIN_code) + (1|Netpiece) , data=Vesterg.multilab.dat, family="binomial")
#summary(Vesterg.multilab.glmm.TL)


######################################################################################

#subsetting by different labs for separate glmm analyses


#Lab A
#subsetting 
Vesterg.labA.dat<- subset(Vesterg.multilab.dat, Testing.lab_coded == "Lab A")


# Calculating the overall mean percent mortality
LabA_mean <- mean(Vesterg.labA.dat$perc_Mortality)

# Calculating the overall standard deviation of percent mortality
LabA_sd <- sd(Vesterg.labA.dat$perc_Mortality)

LabA_cv <- (LabA_sd / LabA_mean)

LabA_cv


#glmm
Vesterg.labA.glmm1 <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ 1 +  (Netpiece|LLIN_code) , data=Vesterg.labA.dat, family="binomial")
summary(Vesterg.labA.glmm1 )


#Lab B

Vesterg.labB.dat<- subset(Vesterg.multilab.dat, Testing.lab_coded == "Lab B")

# Calculating the overall mean percent mortality
LabB_mean <- mean(Vesterg.labB.dat$perc_Mortality)

# Calculating the overall standard deviation of percent mortality
LabB_sd <- sd(Vesterg.labB.dat$perc_Mortality)

LabB_cv <- (LabB_sd / LabB_mean)

LabB_cv


Vesterg.labB.glmm <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ 1 +  (1|LLIN_code) + (1|Netpiece) , data=Vesterg.labB.dat, family="binomial")
summary(Vesterg.labB.glmm )


#Lab C

Vesterg.labC.dat<- subset(Vesterg.multilab.dat, Testing.lab_coded == "Lab C")

# Calculating the overall mean percent mortality
LabC_mean <- mean(Vesterg.labC.dat$perc_Mortality)

# Calculating the overall standard deviation of percent mortality
LabC_sd <- sd(Vesterg.labC.dat$perc_Mortality)

LabC_cv <- (LabC_sd / LabC_mean)

LabC_cv


Vesterg.labC.glmm <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ 1 +  (1|LLIN_code) + (1|Netpiece) , data=Vesterg.labC.dat, family="binomial")
summary(Vesterg.labC.glmm )

#Lab D

Vesterg.labD.dat<- subset(Vesterg.multilab.dat, Testing.lab_coded == "Lab D")

# Calculating the overall mean percent mortality
LabD_mean <- mean(Vesterg.labD.dat$perc_Mortality)

# Calculating the overall standard deviation of percent mortality
LabD_sd <- sd(Vesterg.labD.dat$perc_Mortality)

LabD_cv <- (LabD_sd / LabD_mean)

LabD_cv


Vesterg.labD.glmm <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ 1 +  (1|LLIN_code) + (1|Netpiece) , data=Vesterg.labD.dat, family="binomial")
summary(Vesterg.labD.glmm )

#Lab E

Vesterg.labE.dat<- subset(Vesterg.multilab.dat, Testing.lab_coded == "Lab E")


# Calculating the overall mean percent mortality
LabE_mean <- mean(Vesterg.labE.dat$perc_Mortality)

# Calculating the overall standard deviation of percent mortality
LabE_sd <- sd(Vesterg.labE.dat$perc_Mortality)

LabE_cv <- (LabE_sd / LabE_mean)

LabE_cv


Vesterg.labE.glmm <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ 1 +  (1|LLIN_code) + (1|Netpiece) , data=Vesterg.labE.dat, family="binomial")
summary(Vesterg.labE.glmm )


#Lab F

Vesterg.labF.dat<- subset(Vesterg.multilab.dat, Testing.lab_coded == "Lab F")


# Calculating the overall mean percent mortality
LabF_mean <- mean(Vesterg.labF.dat$perc_Mortality)

# Calculating the overall standard deviation of percent mortality
LabF_sd <- sd(Vesterg.labF.dat$perc_Mortality)

LabF_cv <- (LabF_sd / LabF_mean)

LabF_cv

Vesterg.labF.glmm <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ 1 +  (1|LLIN_code) + (1|Netpiece) , data=Vesterg.labF.dat, family="binomial")
summary(Vesterg.labF.glmm )


#Lab I

Vesterg.labI.dat<- subset(Vesterg.multilab.dat, Testing.lab_coded == "Lab I")


# Calculating the overall mean percent mortality
LabI_mean <- mean(Vesterg.labI.dat$perc_Mortality)

# Calculating the overall standard deviation of percent mortality
LabI_sd <- sd(Vesterg.labI.dat$perc_Mortality)

LabI_cv <- (LabI_sd / LabI_mean)

LabI_cv

Vesterg.labI.glmm <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ 1 +  (Netpiece|LLIN_code) , data=Vesterg.labI.dat, family="binomial")
summary(Vesterg.labI.glmm )

hist(Vesterg.labI.dat$perc_Mortality)
###############################################################################

#removing out Lab C from the main dataset 

Vesterg.multilab.dat.new <- subset(Vesterg.multilab.dat, Testing.lab_coded != "Lab C") 

#######################################################################################

# get 0.6 and 0.3 SD Var from summary (rounded)
Vesterg.multilab.glmm.gambiae <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ Testing.lab_coded + (1|LLIN_code) + (1|Netpiece) , data=Vesterg.multilab.dat.new, family="binomial")
summary(Vesterg.multilab.glmm.gambiae)


Vesterg.multilab.glmm.gambiae1 <-glmer(cbind(Number.of.dead,Totmosquitoes-Number.of.dead)~ Testing.lab_coded +  (Netpiece|LLIN_code) , data=Vesterg.multilab.dat.new, family="binomial")
summary(Vesterg.multilab.glmm.gambiae1)

#using lrtest to determine explanatory power differences between the models
lrtest(Vesterg.multilab.glmm.gambiae,Vesterg.multilab.glmm.gambiae1)

