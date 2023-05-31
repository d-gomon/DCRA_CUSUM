rm(list = ls())
library(dplyr)
library(survival)
library(parallel)
load("df.Rdata")
#devtools::load_all("S:/Scripting&Onderwijs/3.staging/packages/dicapackage_2022")

#devtools::install_github("d-gomon/success")
library(success)
#We make 2 dfs, one for overal survival: datsurv
#               one for duration of stay: datstay


#----------------------
#Data processing here
#----------------------

#Use hospital id as unit identifier
#df$unit <- df$id

#Calculate entrytime
#Calculate difference between minimal entry time and patient entry time
#df <- df %>% mutate(entrytime = CalculateInterval(min(datok), datok, units = "days",
#                                                  lowerinclude = 0, upperinclude = 10000))



#Determine arrival rate
hosps <- unique(df$unit)
psivals <- numeric(length(hosps))
for(i in seq_along(hosps)){
  tdat <- subset(df, unit == hosps[i])
  psivals[i] <- nrow(tdat)/(max(tdat$entrytime) - min(tdat$entrytime))
}

#psivals <- arrival_rate(df)

hist(psivals)

casemixvars = c("geslacht", "bmi_cat", "leeftijd_cat", "charlson", "asa_cat", "dubbeltumor",
                "preoperatieve_tumorcomplicaties", "pt_stadium", "m_stadium", "typok_cat",
                "aanvullende_resectie_ivm_metastasen", "aanvullende_resectie_ivm_doorgroei_tumor")


othervars <- c("entrytime", "unit", "survtime", "censorid")

#Retain only the useful columns
#dat <- select(df, union(casemixvars, othervars))

#-------------------
##Overal survival
#-------------------

datsurv <- df

#Make survtime

datsurv$survtime <- df$tijdovl
datsurv$censorid <- ifelse(datsurv$survtime > 90 | is.na(datsurv$survtime), 0, 1)
datsurv$survtime[which(datsurv$survtime > 90)] <- 90
datsurv$survtime[which(is.na(datsurv$survtime))] <- 90
#Make censorid

datsurv <- datsurv[, colnames(datsurv) %in% c(casemixvars, othervars)]
datsurv <- datsurv[!(rowSums(is.na(datsurv)) > 0),]




#--------------------
##Duration of stay
#--------------------

datstay <- df

#Make survtime

datstay$survtime <- df$opnduur
datstay <- datstay[-which(is.na(datstay$survtime)),]

#Make censorid

datstay$censorid <- rep(1, nrow(datstay))

datstay <- datstay[, colnames(datstay) %in% c(casemixvars, othervars)]
datstay <- datstay[!(rowSums(is.na(datstay)) > 0),]


#Need to create: psivals, phmodels, glmmodels for both survtime and duration of stay.
#All simulation results can be done afterwards.
#Construct charts (Bernoulli, BK and CGR) for both survtime and duration of stay.


#---------------------
#Choose time period over which to construct the charts
#---------------------



#Voor nu: beperk data tot patienten tussen 1 jan 2019 en 31-12-2021
datsurv <- subset(datsurv, entrytime >= 2922)
datstay <- subset(datstay, entrytime >= 2922)


#Hoe veel patienten worden in elk ziekenhuis ongeveer elke dag behandeld?
arr_rate_surv <- arrival_rate(datsurv)
arr_rate_stay <- arrival_rate(datstay)



#----------------------
##Construct the charts
#----------------------

#----------------------
#Risk-adjustment models
#----------------------




#Survival
followupsurv <- 90
#Probably 90
exprfitglmsurv <- as.formula(paste0("(survtime <= followupsurv) & (censorid == 1)~",
                                    paste(casemixvars, collapse=" + ")))
glmmodsurv <- glm(exprfitglmsurv, data = datsurv, family = binomial(link = "logit"))

exprfitsurv <- as.formula(paste0("Surv(survtime, censorid) ~ ", paste(casemixvars, collapse=" + ")))
coxmodsurv <- coxph(exprfitsurv, data = datsurv)

#Duration of stay = 2 weken
followupstay <- 21

exprfitglmstay <- as.formula(paste0("(survtime <= followupstay) & (censorid == 1)~",
                                    paste(casemixvars, collapse=" + ")))
glmmodstay <- glm(exprfitglmstay, data = datstay, family = binomial(link = "logit"))

exprfitstay <- as.formula(paste0("Surv(survtime, censorid) ~ ", paste(casemixvars, collapse=" + ")))
coxmodstay <- coxph(exprfitstay, data = datstay)

#Haal originele data weg om geheugen vrij te maken
rm(df)
#Maak geheugen daadwerkelijk vrij
gc()

psi_surv <- arrival_rate(datsurv)
psi_stay <- arrival_rate(datstay)


#----------------------
#Survival times
#----------------------
hosps <- as.numeric(names(psi_surv))
charts_surv <- vector(mode = "list", length = length(arr_rate_surv))

for(i in 1:length(hosps)){
  print(i)
  charts_surv[[i]]$ber <- bernoulli_cusum(data = subset(datsurv, unit == hosps[i]),
                                          glmmod = glmmodsurv, followup = followupsurv,
                                          theta = log(2))
  charts_surv[[i]]$bk <- bk_cusum(data = subset(datsurv, unit == hosps[i]),
                                  coxphmod = coxmodsurv, theta = log(2))
  charts_surv[[i]]$cgr <- cgr_cusum(data = subset(datsurv, unit == hosps[i]),
                                    coxphmod = coxmodsurv)
}

save(charts_surv, file = "charts_surv.Rdata")



#----------------------
#Duration of stay
#----------------------
hosps <- as.numeric(names(psi_stay))

charts_stay <- vector(mode = "list", length = length(arr_rate_stay))

for(i in 1:length(hosps)){
  print(i)
  charts_stay[[i]]$ber <- bernoulli_cusum(data = subset(datstay, unit == hosps[i]),
                                          glmmod = glmmodstay, followup = followupstay,
                                          theta = -log(2))
  charts_stay[[i]]$bk <- bk_cusum(data = subset(datstay, unit == hosps[i]),
                                  coxphmod = coxmodstay, theta = -log(2))
  charts_stay[[i]]$cgr <- cgr_cusum(data = subset(datstay, unit == hosps[i]),
                                    coxphmod = coxmodstay, detection = "lower")
}

save(charts_stay, file = "charts_stay.Rdata")





save(psi_surv, psi_stay, psivals, glmmodsurv, glmmodstay, coxmodsurv, coxmodstay, charts_surv, charts_stay, file = "DICA_charts.Rdata")
