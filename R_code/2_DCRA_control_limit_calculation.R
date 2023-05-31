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

#------------
#Determine control limits (both surv and duration of stay)
#-----------

#--------------
#Survival
#--------------

psi_surv <- arrival_rate(datsurv)
psi_stay <- arrival_rate(datstay)

#Create Matrix to store control limits for each individual hospital.
h_surv <- matrix(NA, nrow = length(psi_surv), ncol = 3)
colnames(h_surv) <- c("Ber", "BK", "CGR")
rownames(h_surv) <- names(psi_surv)


#Stel eerste 4 ziekenhuizen zijn gedaan: dit kun je checken door print(h_surv)
#Verander for loop dan in for(i in 5:71) (71 altijd als laatste, 5 is dan 4 + 1)
for(i in 1:20){
  print(i)
  h_surv[i, 1] <- bernoulli_control_limit(time = 1095, followup = followupsurv, psi = psi_surv[i], n_sim = 200,
                                          glmmod = glmmodsurv, baseline_data = datsurv, theta = log(2), pb = TRUE)$h
  h_surv[i,2] <- bk_control_limit(time = 1095, psi = psi_surv[i], theta = log(2), coxphmod = coxmodsurv, n_sim = 200,
                                  baseline_data = datsurv, pb = TRUE)$h
  h_surv[i,3] <- cgr_control_limit(time = 1095, psi = psi_surv[i], coxphmod = coxmodsurv, n_sim = 200,
                                   baseline_data = datsurv, pb = TRUE, ncores = 5)$h
}
save(h_surv, file = "h_surv1_20.Rdata")
for(i in 21:40){
  print(i)
  h_surv[i, 1] <- bernoulli_control_limit(time = 1095, followup = followupsurv, psi = psi_surv[i], n_sim = 200,
                                          glmmod = glmmodsurv, baseline_data = datsurv, theta = log(2), pb = TRUE)$h
  h_surv[i,2] <- bk_control_limit(time = 1095, psi = psi_surv[i], theta = log(2), coxphmod = coxmodsurv, n_sim = 200,
                                  baseline_data = datsurv, pb = TRUE)$h
  h_surv[i,3] <- cgr_control_limit(time = 1095, psi = psi_surv[i], coxphmod = coxmodsurv, n_sim = 200,
                                   baseline_data = datsurv, pb = TRUE, ncores = 5)$h
}
save(h_surv, file = "h_surv21_40.Rdata")
for(i in 41:60){
  print(i)
  h_surv[i, 1] <- bernoulli_control_limit(time = 1095, followup = followupsurv, psi = psi_surv[i], n_sim = 200,
                                          glmmod = glmmodsurv, baseline_data = datsurv, theta = log(2), pb = TRUE)$h
  h_surv[i,2] <- bk_control_limit(time = 1095, psi = psi_surv[i], theta = log(2), coxphmod = coxmodsurv, n_sim = 200,
                                  baseline_data = datsurv, pb = TRUE)$h
  h_surv[i,3] <- cgr_control_limit(time = 1095, psi = psi_surv[i], coxphmod = coxmodsurv, n_sim = 200,
                                   baseline_data = datsurv, pb = TRUE, ncores = 5)$h
}
save(h_surv, file = "h_surv41_60.Rdata")
for(i in 61:71){
  print(i)
  h_surv[i, 1] <- bernoulli_control_limit(time = 1095, followup = followupsurv, psi = psi_surv[i], n_sim = 200,
                                          glmmod = glmmodsurv, baseline_data = datsurv, theta = log(2), pb = TRUE)$h
  h_surv[i,2] <- bk_control_limit(time = 1095, psi = psi_surv[i], theta = log(2), coxphmod = coxmodsurv, n_sim = 200,
                                  baseline_data = datsurv, pb = TRUE)$h
  h_surv[i,3] <- cgr_control_limit(time = 1095, psi = psi_surv[i], coxphmod = coxmodsurv, n_sim = 200,
                                   baseline_data = datsurv, pb = TRUE, ncores = 5)$h
}
save(h_surv, file = "h_surv61_71.Rdata")
#ALTIJD ALS JE OPSLAAT, CHECK DAT JE NIET DE OUDE FILE OVERSCHRIJFT!!! Maak er 
#bijvoorbeeld file = "control_surv_DICA2.Rdata" van
save(h_surv, psi_surv, followupsurv, file = "control_surv_DICA.Rdata")

#--------------
#Duration of stay
#--------------

psi_stay <- arrival_rate(datstay)

#Create Matrix to store control limits for each individual hospital.
h_stay <- matrix(NA, nrow = length(psi_stay), ncol = 3)
colnames(h_surv) <- c("Ber", "BK", "CGR")
rownames(h_surv) <- names(psi_stay)


#Stel eerste 4 ziekenhuizen zijn gedaan: dit kun je checken door print(h_stay)
#Verander for loop dan in for(i in 5:71) (71 altijd als laatste, 5 is dan 4 + 1)
for(i in 1:20){
  print(i)
  h_stay[i,1] <- bernoulli_control_limit(time = 1095, followup = followupstay, psi = psi_stay[i], n_sim = 200,
                                         glmmod = glmmodstay, baseline_data = datstay, theta = -log(2), pb = TRUE)$h
  h_stay[i,2] <- bk_control_limit(time = 1095, psi = psi_stay[i], theta = -log(2), coxphmod = coxmodstay, n_sim = 200,
                                  baseline_data = datstay, pb = TRUE)$h
  h_stay[i,3] <- cgr_control_limit(time = 1095, psi = psi_stay[i], coxphmod = coxmodstay, n_sim = 200,
                                   baseline_data = datstay, pb = TRUE, ncores = 5, detection = "lower")$h
}
save(h_stay, file = "h_stay1_20.Rdata")

for(i in 21:40){
  print(i)
  h_stay[i,1] <- bernoulli_control_limit(time = 1095, followup = followupstay, psi = psi_stay[i], n_sim = 200,
                                         glmmod = glmmodstay, baseline_data = datstay, theta = -log(2), pb = TRUE)$h
  h_stay[i,2] <- bk_control_limit(time = 1095, psi = psi_stay[i], theta = -log(2), coxphmod = coxmodstay, n_sim = 200,
                                  baseline_data = datstay, pb = TRUE)$h
  h_stay[i,3] <- cgr_control_limit(time = 1095, psi = psi_stay[i], coxphmod = coxmodstay, n_sim = 200,
                                   baseline_data = datstay, pb = TRUE, ncores = 5, detection = "lower")$h
}
save(h_stay, file = "h_stay21_40.Rdata")

for(i in 41:60){
  print(i)
  h_stay[i,1] <- bernoulli_control_limit(time = 1095, followup = followupstay, psi = psi_stay[i], n_sim = 200,
                                         glmmod = glmmodstay, baseline_data = datstay, theta = -log(2), pb = TRUE)$h
  h_stay[i,2] <- bk_control_limit(time = 1095, psi = psi_stay[i], theta = -log(2), coxphmod = coxmodstay, n_sim = 200,
                                  baseline_data = datstay, pb = TRUE)$h
  h_stay[i,3] <- cgr_control_limit(time = 1095, psi = psi_stay[i], coxphmod = coxmodstay, n_sim = 200,
                                   baseline_data = datstay, pb = TRUE, ncores = 5, detection = "lower")$h
}
save(h_stay, file = "h_stay41_60.Rdata")
for(i in 61:71){
  print(i)
  h_stay[i,1] <- bernoulli_control_limit(time = 1095, followup = followupstay, psi = psi_stay[i], n_sim = 200,
                                         glmmod = glmmodstay, baseline_data = datstay, theta = -log(2), pb = TRUE)$h
  h_stay[i,2] <- bk_control_limit(time = 1095, psi = psi_stay[i], theta = -log(2), coxphmod = coxmodstay, n_sim = 200,
                                  baseline_data = datstay, pb = TRUE)$h
  h_stay[i,3] <- cgr_control_limit(time = 1095, psi = psi_stay[i], coxphmod = coxmodstay, n_sim = 200,
                                   baseline_data = datstay, pb = TRUE, ncores = 5, detection = "lower")$h
}
save(h_stay, file = "h_stay61_71.Rdata")

save(h_stay, psi_stay, followupstay, file = "control_stay_DICA.Rdata")

