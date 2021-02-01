library(dplyr)
library(readr)
library(ggplot2)
library(scales)
path <- "/Users/Pinar/git/2018-paper-graphene/sigcomm-paper/figures/"

dat <- read_delim("/Users/Pinar/git/2018-paper-graphene/basic_sim/data/ping-pong.csv",
                  delim = '\t',
                  col_names = c('TFPs', 'blkSize', 'bound', 'fraction', 'mempoolSize', 'fprS', 'realFprS', 'fprR', 'a', 'b', 'x_star', 'y_star', 'z', 'count', 'observed_FPs', 'decoded', 'pingpong_decode', 'graphene', 'first_IBLT', 'first_BF', 'second_IBLT', 'second_BF', 'extra', 'IBLT_rows_first', 'IBLT_rows_second', 'compact'),
                  col_types = cols(.default = col_double(), decoded = col_character(), pingpong_decode = col_character()))
dat$de <- dat$decoded =='True'
dat$dePing <- dat$pingpong_decode =='True'
dat$combine <- dat$de | dat$dePing

##################################################################################################
dat2 <- read_delim("/Users/Pinar/Desktop/empirical-dist/short-half.csv",
                  delim = '\t',
                  col_names = c('TFPs', 'blkSize', 'bound', 'fraction', 'mempoolSize', 'fprS', 'fprR', 'a', 'b', 'x_star', 'y_star', 'z', 'count', 'observed_FPs', 'decoded', 'graphene', 'first_IBLT', 'first_BF', 'second_IBLT', 'second_BF', 'extra', 'IBLT_rows_first', 'IBLT_rows_second', 'compact'),
                  col_types = cols(.default = col_double(), decoded = col_character()))
dat2$de <- dat2$decoded =='True'
dat2$dePing <- FALSE
tru <- cbind(dat2[,1:15], pingpong_decode='False', dat2[,16:26])

dat <- rbind(dat, tru)

##################################################################################################
# Analyze how many graphene blocks decoded
r1 <- dat %>% 
  group_by(blkSize, bound, fraction) %>% 
  summarize(avg=mean(de), 
            var=var(de), 
            cnt=n(), 
            ci_max=ifelse(avg == 1, -log(.05)/n(),1-avg+ 1.96*sqrt(avg*(1-avg)/n())), 
            ci_min=ifelse(avg == 1, -log(.05)/n(),1-avg- 1.96*sqrt(avg*(1-avg)/n())), 
            type='without'
            )
dat$combine <- dat$de | dat$dePing
r2 <- dat %>% 
  group_by(blkSize, bound, fraction) %>% 
  summarize(avg=mean(combine), 
            var=var(combine), 
            cnt=n(), 
            ci_max=ifelse(avg == 1,max(0,1-avg+ (-log(.05)/n())),1-avg+ 1.96*sqrt(avg*(1-avg)/n())), 
            ci_min=ifelse(avg == 1,0, 1-avg-1.96*sqrt(avg*(1-avg)/n())),
            ci_min=max(0,ci_min),
            type='with ping-pong'
            )
results <- rbind(r1, r2)

pdf(paste0(path,"decode-rate-graphene-extended.pdf"),w=5,h=4.5)
ggplot(data=results)+# %>% filter(bound == 0.1))+
  theme_gray(16)+
  aes(x=fraction, y=1-avg, color=factor(type), group=factor(type))+
  geom_point(size=1)+
  geom_line()+
  xlab("Fraction of blk receiver has in mempool")+
  ylab("Decode failure probability\n (logscale)")+
  geom_errorbar(aes(ymin=ci_min,ymax=ci_max), size=.4, alpha=.5)+
  scale_x_continuous(breaks=seq(0, 1, by=0.1))+
  scale_y_log10()+
  geom_hline(aes(yintercept=bound), linetype=2, color='black')+
  #facet_wrap(blkSize ~ bound, ncol=3,scales='free')
  facet_grid(blkSize ~ .)+
  theme(legend.position="top")+
  scale_color_discrete(name="")
dev.off()

##################################################################################################
# Analyze each piece of message
results.1 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(first_IBLT), type='first_IBLT')
results.2 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(first_BF), type='first_BF')
results.3 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(second_IBLT), type='second_IBLT')
results.4 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(second_BF), type='second_BF')
results.5 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(extra), type='extra')
results <- bind_rows(results.1, results.2, results.3, results.4,results.5)
#cb <- dat %>% group_by(blkSize, bound, fraction)  %>% summarize(size=mean(compact), type='Compact')#, ci=1.96*sqrt(avg*(1-avg)/n()))

pdf(paste0(path,"size-by-parts.pdf"),w=5,h=4)
ggplot(data=results)+# %>% filter (bound == 0.005))+
  aes(x=fraction, y=(size/1024))+
  theme_gray(18)+
  geom_bar(stat = 'identity', aes(fill = type))+#, position = 'dodge')+
  xlab("Fraction of blk in mempool")+
  ylab("Avg encoding size by parts (KB)")+
  geom_line(aes(y=blkSize/1024*(1-fraction)*ifelse(blkSize==200,1,3)+8*blkSize/1024), linetype=2, color='black')+
  #geom_line(data=cb, aes(color=type), linetype=2, color='black')+
  #facet_grid(bound ~ ., scales='free')+
  scale_x_continuous(breaks=seq(0, 1, by=0.2))+
  scale_y_continuous(breaks=pretty_breaks(3))+
  #labs(fill = "Msg")+
  scale_fill_discrete(name="",
                      breaks=c("extra", "first_BF", "first_IBLT", "second_BF", "second_IBLT"),
                      labels=c("getdata","BF S","IBLT I", "BF R", "IBLT J"))+
  theme(legend.position="top")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  facet_grid( blkSize ~ ., scales='free')
dev.off()

##################################################################################################
#dat[dat$x_star == -1,]$x_star = 0
# results.1 <- dat %>% group_by(fraction, bound) %>% summarize(avg=mean(x_star), var=var(x_star), cnt=n(), x='x*', ci=1.96*sqrt(avg*(1-avg)/n()))
# results.2 <-data.frame(fraction=results.1$fraction, bound=results.1$bound, avg = 1500*results.1$fraction, var=0, cnt=0, x='x', ci=0)
# results.1$val <- results.2$avg - results.1$avg
# results <- bind_rows(results.1, results.2)

r <- dat %>% filter(IBLT_rows_second != 0)
result <- data.frame(fraction=r$fraction, blkSize=r$blkSize, bound=r$bound, result=(r$blkSize*r$fraction) - r$x_star, x_star = r$x_star, actual_x =(r$blkSize*r$fraction), z=r$z)
result$val <- result$result >= 0
final <- result %>% group_by(blkSize, fraction, bound) %>% summarize(avg=mean(val), var=var(val), cnt=n(), ci=1.96*sqrt(avg*(1-avg)/n()))

pdf(paste0(path,"x_star.pdf"),w=5,h=4.5)
ggplot(data=final %>% filter (fraction <= 0.9))+
  aes(x=fraction, y=avg)+
  theme_gray(16)+
  geom_point(size=0.5)+
  geom_line()+
  xlab("Fraction of block in mempool")+
  ylab("Fraction of time \nx* is a lower bound")+
  geom_errorbar(aes(ymin=avg+ci,ymax=avg-ci), width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=seq(0, 0.9, by=0.1))+
  scale_y_continuous(breaks=pretty_breaks(4))+
  geom_hline(aes(yintercept=1-bound), linetype=2, color='red')+
  #facet_wrap(blkSize ~ bound, ncol=3,scales='free')
  facet_grid(blkSize ~ .)
  #geom_hline(yintercept=.999, linetype=2)
dev.off()

##################################################################################################
r <- dat %>% filter(IBLT_rows_second != 0)
result <- data.frame(fraction=r$fraction, blkSize=r$blkSize, bound=r$bound, result=r$y_star - r$observed_FPs, y_star = r$y_star, FPs = r$observed_FPs)
result$val <- result$result >= 0
final <- result %>% group_by(blkSize, fraction, bound) %>% summarize(avg=mean(val), var=var(val), cnt=n(), ci=1.96*sqrt(avg*(1-avg)/n()))

pdf(paste0(path,"y_star.pdf"),w=5,h=4.5)
ggplot(data=final %>% filter (fraction <= 0.90))+#%>% filter (bound == 0.005)  +
  aes(x=fraction, y=avg)+
  theme_gray(16)+
  geom_point(size=0.5)+
  geom_line()+
  xlab("Fraction of block in mempool")+
  ylab("Fraction of time \ny* is an upper bound")+
  geom_errorbar(aes(ymin=avg+ci,ymax=avg-ci), width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=seq(0, 0.9, by=0.1))+
  #scale_y_continuous(breaks=seq(0, 1, by=0.1))+
  geom_hline(aes(yintercept=1-bound), linetype=2, color='red')+
  #facet_wrap(blkSize ~ bound, ncol=3,scales='free')
  facet_grid(blkSize ~ .)
  #geom_hline(yintercept=.999, linetype=2)
dev.off()

##################################################################################################
results.1 <- dat %>% group_by(fraction, bound) %>% summarize(avg=mean(y_star), var=var(y_star), cnt=n(), x='estimated', ci=1.96*sqrt(avg*(1-avg)/n()))
results.2 <- dat %>% group_by(fraction, bound) %>% summarize(avg=mean(observed_FPs), var=var(observed_FPs), cnt=n(), x='true', ci=1.96*sqrt(avg*(1-avg)/n()))
results <- bind_rows(results.1, results.2)

ggplot(data=results %>% filter(fraction <= 0.95))+# %>% filter(bound == 0.1))+
  aes(x=fraction, y=avg, color=factor(x), group=factor(x))+
  geom_point(size=0.25)+
  geom_line()+
  xlab("Fraction of blk receiver has in mempool")+
  ylab("Fraction of Graphene blks decoded")+
  geom_errorbar(aes(ymin=avg+var,ymax=avg-var), width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=seq(0, 1, by=0.1))+
  #geom_hline(aes(yintercept=1-bound), linetype=2, color='red')+
  #facet_wrap(blkSize ~ bound, ncol=3,scales='free')
  facet_grid(bound ~ .)+
  #geom_hline(yintercept=.999, linetype=2)
  theme_gray(14)

##################################################################################################
results.1 <- dat %>% group_by(fraction, bound) %>% summarize(avg=mean(a), var=var(a), cnt=n(), x='estimated', ci=1.96*sqrt(avg*(1-avg)/n()))
results.2 <- dat %>% group_by(fraction, bound) %>% summarize(avg=mean(observed_FPs), var=var(observed_FPs), cnt=n(), x='true', ci=1.96*sqrt(avg*(1-avg)/n()))
results <- bind_rows(results.1, results.2)

ggplot(data=results %>% filter(fraction <= 0.95))+# %>% filter(bound == 0.1))+
  aes(x=fraction, y=avg, color=factor(x), group=factor(x))+
  geom_point(size=0.25)+
  geom_line()+
  xlab("Fraction of blk receiver has in mempool")+
  ylab("Fraction of Graphene blks decoded")+
  geom_errorbar(aes(ymin=avg+var,ymax=avg-var), width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=seq(0, 1, by=0.1))+
  #geom_hline(aes(yintercept=1-bound), linetype=2, color='red')+
  #facet_wrap(blkSize ~ bound, ncol=3,scales='free')
  facet_grid(bound ~ .)+
  #geom_hline(yintercept=.999, linetype=2)
  theme_gray(14)

##################################################################################################
results <- dat %>% group_by(TFPs, blkSize) %>% summarize(avg=mean(observed_FPs), var=var(observed_FPs), cnt=n()) #ci=1.96*sqrt(avg*(1-avg)/n()), )

results <- dat %>% group_by(TFPs, blkSize) %>% summarize(avg=mean(FPR), var=var(FPR), cnt=n()) #ci=1.96*sqrt(avg*(1-avg)/n()), )

ggplot(data=results)+# %>% filter(TFPs == 2000))+
  aes(x=TFPs, y=avg)+
  geom_point()+
  geom_line()+
  xlab("Txns in mempool not in blk")+
  #ylab("Number of txns that falsely pass \n thru sender's BF")+
  ylab("FPR of sender's BF")+
  geom_errorbar(aes(ymin=avg+var,ymax=avg-var), width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=seq(0, 4000, by=500))+
  facet_grid(blkSize ~ .)
##################################################################################################
