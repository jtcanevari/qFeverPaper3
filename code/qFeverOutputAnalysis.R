#------------------------------------------------------------#
# Q fever output analysis
#------------------------------------------------------------#
library(lubridate); library(dplyr); library(ggplot2); library(scales)


#unlist simulations
for(i in 1:nsims){
  modelOut[[i]] <- do.call(data.frame, modelOut[[i]])
  modelOut[[i]]$iter <- i}
dat <- do.call(rbind, modelOut)

#flag kidding event with a 1
dat <- dat %>% group_by(iter) %>% mutate(kevent = ave(K, FUN=function(x) c(0, diff(x))))
#flag kidding of an infected doe with a 1
dat <- dat %>% group_by(iter) %>% mutate(kinf.event=ave(KI, FUN=function(x) c(0, diff(x))))
#flag abortions
dat <- dat %>% group_by(iter) %>% mutate(abort = ave(A, FUN=function(x) c(0, diff(x))))

#obtain year
dat$date <- as.Date.numeric(dat$time, origin = dmy("13-04-2014"))
dat$year <- year(dat$date)

#count events
tdat <- plyr::count(dat[dat$kevent == 1,], c("iter","year"))
inf <- plyr::count(dat[dat$kinf.event == 1,], c("iter","year"))
abort <- plyr::count(dat[dat$abort == 1,], c("iter","year"))

#merge kiddings with count of infected does kiddings
tdat <- merge(tdat, inf, by=c("iter","year"), all.x = TRUE )
colnames(tdat)[c(3,4)] <- c("total.kid", "inf.kid")

#merge kiddings with count of abortions
tdat <- merge(tdat, abort, by=c("iter","year"), all.x = TRUE )
colnames(tdat)[5] <- c("abortions")

#calculate shedding incidence risk and abortion risk
tdat$incidence <- round(tdat$inf.kid/tdat$total.kid,3)
tdat$abortion.risk<- round(tdat$abortions / (tdat$abortions + tdat$total.kid),3)
tdat <- subset(tdat, !year==max(year))

#plot shedding by farm and year as dot plot.
# windows();plot(tdat$incidence ~ tdat$year,ylim = c(0,0.2),col='red', ylab='Prevalence', xlab='Time')

# setwd('C:\\Users\\jcanevari\\Documents\\Projects\\PhD\\qFeverPaper3\\docs')

# pdf('outputDutchV03.pdf')
ggplot()+
  geom_point(aes(x=year,y=incidence, colour=factor(iter)),data=tdat)+
  geom_line(aes(x=year,y=incidence, colour = factor(iter)),data=tdat)+
  scale_x_continuous(breaks=2014:2028, labels = 1:15)+
  theme(legend.position = "none")+
  labs(title = "Inicidence of shedding does at kidding",
       subtitle = "Dutch version 50 iters 15 years",
       caption = paste("alhpa = ", params[6], "beta = ", params[7]),
       x = "year", y = "incidence")


ggplot()+
  geom_point(aes(x=year,y=abortion.risk, colour = factor(iter)),data=tdat)+
  geom_line(aes(x=year,y=abortion.risk, colour = factor(iter)),data=tdat)+
  scale_x_continuous(breaks=2014:2028, labels = 1:15)+
  theme(legend.position = "none")+
  labs(title = "Abortion risk over time",
       subtitle = "Dutch version 100 iters 15 years",
       caption = paste("alhpa = ", params[6], "beta = ", params[7]),
       x = "year", y = "risk")
# dev.off()

#donde termina una epidemia?
tdat$active.kinf <- ifelse(tdat$inf.kid>0,1,0)
tdat$active.abo <- ifelse(tdat$abortions>0,1,0)
tdat$active <- ifelse(tdat$active.kinf == 1 | tdat$active.abo == 1,1,0)
surv.dat <- data.frame(tapply(tdat$year[tdat$active == 1], tdat$iter[tdat$active == 1], max))
colnames(surv.dat)<-'year'

surv.dat$start <- 0
surv.dat$stop <- surv.dat$year-2013 #this is your stop for survival
surv.dat$status <- 1
surv.dat$status[surv.dat$year == max(tdat$year)] <-0

library(survival); library(survminer)
km_fit <- survfit(Surv(stop, status) ~ 1, data=surv.dat)

ggsurvplot(km_fit, data = surv.dat, ggtheme = theme_bw(), risk.table = TRUE,xlab='Time (years)',conf.int = FALSE)

#---------------------------------------------------------
# need to know proportion of infected pregnant that are IP, IP2 and JP respectively

#flag kidding of an infected doe with a 1
dat <- dat %>% group_by(iter) %>% mutate(ip.track=ave(IP, FUN=function(x) c(0, diff(x))))
dat$ipKids<-NA
dat$ipKids[dat$ip.track == -1 & dat$kevent == 1] <- 1

#flag kidding of an infected class2 doe with a 1
dat <- dat %>% group_by(iter) %>% mutate(ip2.track=ave(IP2, FUN=function(x) c(0, diff(x))))
dat$ip2Kids<-NA
dat$ip2Kids[dat$ip2.track == -1 & dat$kevent == 1] <- 1

#flag kidding of an persistently infected doe with a 1
dat <- dat %>% group_by(iter) %>% mutate(jp.track=ave(JP, FUN=function(x) c(0, diff(x))))
dat$jpKids<-NA
dat$jpKids[dat$jp.track == -1 & dat$kevent == 1] <- 1

tot.ipkid <- sum(dat$ipKids, na.rm = T)
tot.ip2kid <- sum(dat$ip2Kids, na.rm = T)
tot.jpkid <- sum(dat$jpKids, na.rm = T)
tot.infkid <- sum(tot.ipkid, tot.ip2kid, tot.jpkid)

ip.prop <- tot.ipkid/tot.infkid
ip2.prop <- tot.ip2kid/tot.infkid
jp.prop <- tot.jpkid/tot.infkid


df <- data.frame(
  group = c("IP", "IP2", "JP"),
  value = c(ip.prop, ip2.prop, jp.prop)
)

loc.lab <- df$value/2 + cumsum(df$value)[1:nrow(df)]
lab <- round(df$value,2)

# png('distShedComp.png')
ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  xlab('')+
  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value)), size=5)
  
#