# script to compare BRM with and without proposal information against
# ARIMA model for accuracy on Grainger expenditure data

library(data.table)
library(lubridate)
library(forecast)
library(janitor)
library(ggplot2)
library(rstanarm)
library(broom)
library(tidyverse)
library(performance)
library(see)
library(bayesplot)
library(bayestestR)
library(ggeffects)
library(BayesPostEst)
library(bayesforecast)

rm(list=ls())

## read in data
setwd('~/Projects/arima_stuff')
exp <- fread('MonthlyCGExp_7252023.csv')
proposals <- fread('Proposal_Details.csv')

## expenditure data
exp <- clean_names(exp)
exp[,cal_year:=year(cg_snapshot_dt)]
exp[,cal_monthn:=month(cg_snapshot_dt)]
exp[,cal_monthd:=month(cg_snapshot_dt,label=TRUE)]

exp <- exp[,sum(exp_ytd),by=c('cg_year','cg_quarter','cg_snapshot_dt','cal_year','cal_monthn','cal_monthd')]
setkeyv(exp,cols=c('cg_year','cg_quarter','cg_snapshot_dt'))
sfy_monthn <- c(1,2,3,4,5,6,7,8,9,10,11,12)
exp<- cbind(exp,rep(sfy_monthn,5))
colnames(exp)[7] <- 'month_exp'
colnames(exp)[8] <- 'sfy_month'


exp[,lag:=shift(month_exp,fill=first(month_exp)),by=cg_year]
exp[,lag:=shift(month_exp,fill=first(month_exp)),by=cg_year]
exp[,diff:=month_exp-lag]

# Impute values for July of each SFY
impute <- exp[diff!=0,.(min(diff),median(diff),max(diff)),by=cg_year]

exp[cg_year == 2019 & diff == 0, `:=`(imputed_diff1, impute[cg_year == 2019, V1])]
exp[cg_year == 2020 & diff == 0, `:=`(imputed_diff1, impute[cg_year == 2020, V1])]
exp[cg_year == 2021 & diff == 0, `:=`(imputed_diff1, impute[cg_year == 2021, V1])]
exp[cg_year == 2022 & diff == 0, `:=`(imputed_diff1, impute[cg_year == 2022, V1])]
exp[cg_year == 2023 & diff == 0, `:=`(imputed_diff1, impute[cg_year == 2023, V1])]

exp[cg_year == 2019 & diff == 0, `:=`(imputed_diff2, impute[cg_year == 2019, V2])]
exp[cg_year == 2020 & diff == 0, `:=`(imputed_diff2, impute[cg_year == 2020, V2])]
exp[cg_year == 2021 & diff == 0, `:=`(imputed_diff2, impute[cg_year == 2021, V2])]
exp[cg_year == 2022 & diff == 0, `:=`(imputed_diff2, impute[cg_year == 2022, V2])]
exp[cg_year == 2023 & diff == 0, `:=`(imputed_diff2, impute[cg_year == 2023, V2])]

exp[cg_year == 2019 & diff == 0, `:=`(imputed_diff3, impute[cg_year == 2019, V3])]
exp[cg_year == 2020 & diff == 0, `:=`(imputed_diff3, impute[cg_year == 2020, V3])]
exp[cg_year == 2021 & diff == 0, `:=`(imputed_diff3, impute[cg_year == 2021, V3])]
exp[cg_year == 2022 & diff == 0, `:=`(imputed_diff3, impute[cg_year == 2022, V3])]
exp[cg_year == 2023 & diff == 0, `:=`(imputed_diff3, impute[cg_year == 2023, V3])]

# fill in missing values with diffs 1, 2, or 3
exp[is.na(imputed_diff1), `:=`(imputed_diff1, diff)]
exp[is.na(imputed_diff2), `:=`(imputed_diff2, diff)]
exp[is.na(imputed_diff3), `:=`(imputed_diff3, diff)]
keycols=c('cg_year','cg_quarter','cg_snapshot_dt')
setkeyv(exp,cols=keycols)

pre_merge_exp <- exp[,.(sfy_year=cg_year,cal_year,cal_monthd,imputed_diff2)]

## proposal data
proposals <- clean_names(proposals)
proposals[,sbmt_year:=year(proposal_submission_date)]
proposals[,sbmt_month:=month(proposal_submission_date,label=TRUE)]
pre_merge_proposals <- proposals[,length(proposal_number),by=c('sbmt_year','sbmt_month')]


# merge expenditures and proposals
merged <- merge(pre_merge_exp,pre_merge_proposals,by.x=c('cal_year','cal_monthd'),by.y=c('sbmt_year','sbmt_month'))
merged$cal_monthd <- factor(merged$cal_monthd, levels=c('Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'))
colnames(merged)[5] <- 'prop_count'

# model building BRM
bm1 <- stan_glm(imputed_diff2 ~ sfy_year + cal_monthd, data=merged)   
bm2 <- stan_glm(imputed_diff2 ~ sfy_year + cal_monthd + prop_count, data=merged) 

pp_check(bm1,"dens_overlay")
pp_check(bm2,"dens_overlay")

p_direction(bm1)
p_direction(bm2)

#### 
# build test set sfy 2023

test <- merged[sfy_year==2023]
df <- cbind(test[,.(cal_monthd,imputed_diff2)],
            exp_mdl=apply(posterior_predict(bm1,merged[sfy_year==2023]),2,median),
            exp_prop_mdl=apply(posterior_predict(bm2,merged[sfy_year==2023]),2,median))

m_df <- reshape2::melt(df,id='cal_monthd')
ggplot(m_df,aes(x=cal_monthd,y=value,group=variable)) + geom_point(aes(color=variable)) + 
  geom_line(aes(color=variable)) + theme_bw()



df[,exp_dev:=100*((imputed_diff2 - exp_mdl)/imputed_diff2)]
df[,exp_prop_dev:=100*((imputed_diff2 - exp_prop_mdl)/imputed_diff2)]

#### model building ARIMA
# training set sfy 2019-2022
arima_train <- exp[cg_year <=2022]

# convert to time series
train_ts <- ts(arima_train$imputed_diff1,start=c(2018,7),end=c(2022,6),frequency=12)
if(sum(is.na(train_ts)) !=0) {
  print("Missing Values")
} else {
  print("No missing values")
}

est_lambda<- BoxCox.lambda(train_ts)
arima_model <- auto.arima(train_ts,lambda=est_lambda,approximation=FALSE)
arima_forecast <- forecast(arima_model, level=c(95),h =1*12)
