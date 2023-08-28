# Predict Grainger sfy 2023 expenditures based on sfy 2019-2022
# Models: (i)  BRM exp 
# (ii) BRM exp + proposal counts 
# (iii) ARIMA model 
# (iv) BSTS variants

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
library(bsts)

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
exp_model <- stan_glm(imputed_diff2 ~ sfy_year + cal_monthd, data=merged)   
exp_prop_model <- stan_glm(imputed_diff2 ~ sfy_year + cal_monthd + prop_count, data=merged) 

pp_check(exp_model,"dens_overlay")
pp_check(exp_prop_model,"dens_overlay")

p_direction(exp_model)
p_direction(exp_prop_model)

#### 

# predict and compare to observed values of sfy 2024
df <- cbind(merged[sfy_year==2023][,.(sfy_year,cal_monthd,imputed_diff2)],
            exp_mdl=apply(posterior_predict(exp_model,merged[sfy_year==2023][,.(sfy_year,cal_monthd)]),2,median),
            exp_prop_mdl=apply(posterior_predict(exp_prop_model,merged[sfy_year==2023][,.(sfy_year,cal_monthd,prop_count)]),2,median))

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

df2 <- cbind(df,arima=as.vector(arima_forecast$mean))

## BSTS approach...local + seasonal
ss <- AddLocalLinearTrend(list(), train_ts)
ss <- AddSeasonal(ss, train_ts, nseasons = 48)
bsts1 <- bsts(train_ts,
               state.specification = ss,
               niter = 1000)
pred1 <- predict(bsts1, horizon = 12)
plot(pred1,plot.original=48)
df3 <- cbind(df2,bsts1=pred1$mean)

# semi-local linear controls uncertainty through stationarity
ss2 <- AddSemilocalLinearTrend(list(), train_ts)
ss2 <- AddSeasonal(ss2, train_ts, nseasons = 48)
bsts2 <- bsts(train_ts, state.specification = ss2, niter = 1000)
pred2 <- predict(bsts2, horizon = 12)

df4 <- cbind(df3,bsts2=pred2$mean)
df4[,ex_pd:=100*(imputed_diff2-exp_mdl)/imputed_diff2]
df4[,ex_pr_pd:=100*(imputed_diff2-exp_prop_mdl)/imputed_diff2]
df4[,arima_pd:=100*(imputed_diff2-arima)/imputed_diff2]
df4[,bsts1_pd:=100*(imputed_diff2-bsts1)/imputed_diff2]
df4[,bsts2_pd:=100*(imputed_diff2-bsts2)/imputed_diff2]

m_df4 <- reshape2::melt(df4,id='cal_monthd')
ggplot(m_df4[variable %in% c('cal_monthd','imputed_diff2','bsts1','bsts2')],aes(x=cal_monthd,y=value,group=variable)) + 
  geom_point(aes(color=variable)) + 
  geom_line(aes(color=variable)) + geom_hline(yintercept =0) + theme_bw()


df5 <- df4[,.(cal_monthd,ex_pd,ex_pr_pd,arima_pd,bsts1_pd,bsts2_pd )]
m_df5 <- reshape2::melt(df5,id='cal_monthd')

ggplot(m_df5,aes(x=cal_monthd,y=value,group=variable)) + geom_point(aes(color=variable)) + 
  geom_line(aes(color=variable)) + geom_hline(yintercept =0) + theme_bw() + 
  ylab("Percent Deviation From Observed")
