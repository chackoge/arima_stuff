library(data.table)
library(lubridate)
library(forecast)
library(janitor)
library(ggplot2)
rm(list=)
x <- fread('MonthlyCGExp_7252023.csv')
x <- clean_names(x)
x[,cal_monthn:=month(cg_snapshot_dt)]
x[,cal_monthd:=month(cg_snapshot_dt,label=TRUE)]

y <- x[,sum(exp_ytd),by=c('cg_year','cg_quarter','cg_snapshot_dt','cal_monthn','cal_monthd')]
setkeyv(y,cols=c('cg_year','cg_quarter','cg_snapshot_dt'))
sfy_monthn <- c(1,2,3,4,5,6,7,8,9,10,11,12)
y <- cbind(y,rep(sfy_monthn,5))

y[,lag:=shift(V1,fill=first(V1)),by=cg_year]
y[,diff:=V1-lag]

# Impute values for July of each SFY
y1 <- y[diff!=0,.(min(diff),median(diff),max(diff)),by=cg_year]

y[cg_year == 2019 & diff == 0, `:=`(imputed_diff1, y1[cg_year == 2019, V1])]
y[cg_year == 2020 & diff == 0, `:=`(imputed_diff1, y1[cg_year == 2020, V1])]
y[cg_year == 2021 & diff == 0, `:=`(imputed_diff1, y1[cg_year == 2021, V1])]
y[cg_year == 2022 & diff == 0, `:=`(imputed_diff1, y1[cg_year == 2022, V1])]
y[cg_year == 2023 & diff == 0, `:=`(imputed_diff1, y1[cg_year == 2023, V1])]

y[cg_year == 2019 & diff == 0, `:=`(imputed_diff2, y1[cg_year == 2019, V2])]
y[cg_year == 2020 & diff == 0, `:=`(imputed_diff2, y1[cg_year == 2020, V2])]
y[cg_year == 2021 & diff == 0, `:=`(imputed_diff2, y1[cg_year == 2021, V2])]
y[cg_year == 2022 & diff == 0, `:=`(imputed_diff2, y1[cg_year == 2022, V2])]
y[cg_year == 2023 & diff == 0, `:=`(imputed_diff2, y1[cg_year == 2023, V2])]

y[cg_year == 2019 & diff == 0, `:=`(imputed_diff3, y1[cg_year == 2019, V3])]
y[cg_year == 2020 & diff == 0, `:=`(imputed_diff3, y1[cg_year == 2020, V3])]
y[cg_year == 2021 & diff == 0, `:=`(imputed_diff3, y1[cg_year == 2021, V3])]
y[cg_year == 2022 & diff == 0, `:=`(imputed_diff3, y1[cg_year == 2022, V3])]
y[cg_year == 2023 & diff == 0, `:=`(imputed_diff3, y1[cg_year == 2023, V3])]

# fill in missing values with diffs 1, 2, or 3
y[is.na(imputed_diff1), `:=`(imputed_diff1, diff)]
y[is.na(imputed_diff2), `:=`(imputed_diff2, diff)]
y[is.na(imputed_diff3), `:=`(imputed_diff3, diff)]
keycols=c('cg_year','cg_quarter','cg_snapshot_dt')
setkeyv(y,cols=keycols)

# training set sfy 2019-2022
train <- y[cg_year <=2022]

# convert to time series
train_ts <- ts(train$imputed_diff1,start=c(2018,7),end=c(2022,6),frequency=12)
if(sum(is.na(train_ts)) !=0) {
  print("Missing Values")
} else {
  print("No missing values")
}

# model 
est_lambda<- BoxCox.lambda(train_ts)
my_model <- auto.arima(train_ts,lambda=est_lambda,approximation=FALSE)

# forecast
my_forecast <- forecast(my_model, level=c(95),h =1*12)
predicted <- as.vector(my_forecast$mean)
lower <- as.vector(my_forecast$lower)
upper <- as.vector(my_forecast$upper)
# fc <- data.frame(forecast=as.matrix(my_forecast$mean))
ob <- y[cg_year > 2022,.(cg_year,year(cg_snapshot_dt),cal_monthd,observed=imputed_diff1)]
df <- cbind(ob,predicted,upper,lower)
melted_df <- melt(df[,.(cal_monthd,observed,predicted,upper,lower)],id.vars="cal_monthd")
# set order of months
melted_df$cal_monthd <- factor(melted_df$cal_monthd,levels=c('Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May',"Jun"))
# plot
p <- ggplot(melted_df,aes(x=cal_monthd,y=log10(value),group=variable)) + geom_point() + geom_line(aes(color=variable))
p1 <- p + theme_bw() + xlab("SFY Month") + ylab("log10 Monthly Exp")
pdf(paste('arima_',today(),'.pdf'))
print(p1)
dev.off()

pdf('arima_forecast.pdf')
plot(my_forecast,main="ARIMA Monthly")
dev.off()

df[,diff_perc:=round(100*(predicted-observed)/observed,1)]
fwrite(df,file=paste0('df_',today(),'.csv'))

## Proposal Data
proposals <- fread('Proposal_Details.csv')
proposals <- clean_names(proposals)
proposals[,sbmt_year:=year(proposal_submission_date)]
proposals[,sbmt_month:=month(proposal_submission_date,label=TRUE)]
proposals[total_requested_dollar <= 300000,budget_type:= 't1']
proposals[total_requested_dollar > 300000,budget_type:= 't2']

summ_proposals <- proposals[,.N,by=c('sfy','sbmt_year','sbmt_month','budget_type')][order(sbmt_year,sbmt_month)]
wide_proposal_counts <- dcast(summ_proposals,sbmt_year + sbmt_month ~ budget_type, value.var="N")

## Expenditure data
x <- fread("MonthlyCGExp_7252023.csv")
x <- clean_names(x)
z <- x[, sum(exp_ytd), c("cg_year", "cg_quarter", "cg_snapshot_dt")]
z <- z[order(cg_year, cg_quarter, cg_snapshot_dt)]
z[, `:=`(lag, shift(V1, fill = first(V1))), by = "cg_year"]
z[, `:=`(diff, V1 - lag)]
z <- cbind(z[order(cg_year, cg_quarter, cg_snapshot_dt)], rep(month, 5))
# remove extra row added by recycling due to extra month
# z <- z[-nrow(z)] # not needed in July2023. Need better strategy though

# set index
keycols <- c("cg_year", "cg_quarter", "cg_snapshot_dt")
setkeyv(z, cols = keycols)
