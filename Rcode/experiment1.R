library(dplyr)
library(readr)
library(scales)
library(ggplot2)

path <- "/Users/Pinar/git/2018-paper-graphene/sigcomm-paper/figures/"

dat <- read_delim("/Users/Pinar/git/2018-paper-graphene/basic_sim/data/experiment1.csv",
                  delim = '\t',
                  col_names = c('mempoolSize', 'blkSize', 'bound', 'FPR', 'a', 'IBLT_rows', 'TFPs', 'observed_FPs', 'compactBlk', 'decoded', 'graphene', 'first_IBLT', 'first_BF', 'multiple'))
dat$bool <- dat$decoded =='True'

# Analyze how many graphene blocks decoded
results <- dat %>% group_by(multiple, blkSize, bound) %>% summarize(avg=mean(bool), var=var(bool), cnt=n(), ci=1.96*sqrt(avg*(1-avg)/n()))

# Analyze how many graphene blocks DID NOT decode
#results <- dat %>% group_by(mempoolSize, blkSize) %>% summarize(avg=mean(bool == FALSE), var=var(bool), cnt=n()) #ci=1.96*sqrt(avg*(1-avg)/n()), )

pdf(paste0(path,"decode-rate-graphene-original.pdf"),w=5,h=4)
ggplot(data=results)+theme_gray(18)+
  aes(x=multiple, y=1-avg)+
  geom_point(size=0.25)+
  geom_line()+
  xlab("Txns in mempool not in blk \n as a multiple of blk size")+
  ylab("Decode failure probability")+
  geom_errorbar(aes(ymin=1-avg+ci,ymax=1-avg-ci), width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=pretty_breaks(5))+
  scale_y_continuous(breaks=pretty_breaks(3))+
  geom_hline(aes(yintercept=bound), linetype=2, color='red')+
  #facet_wrap(blkSize ~ bound, ncol=3,scales='free')
  facet_grid(blkSize ~ ., scales='free')
  #geom_hline(yintercept=.999, linetype=2)
  #scale_y_continuous(lim=c(.99,NA))
dev.off()

# Analyze avg size of graphene blks
#results1 <- dat %>% group_by(multiple, blkSize, bound) %>% summarize(avg=mean(compactBlk), var=var(compactBlk), cnt=n(), type='Compact')#, ci=1.96*sqrt(avg*(1-avg)/n()))
results1 <- dat %>% group_by(multiple, blkSize, bound) %>% summarize(avg=mean(graphene), var=var(graphene), cnt=n(), type='Graphene')#, ci=1.96*sqrt(avg*(1-avg)/n()))
results2 <- data.frame(multiple=results1$multiple, blkSize=results1$blkSize, bound=results1$bound, avg=8*results1$blkSize, var=0, cnt=0, type='Compact')
results_size <- rbind(data.frame(results1), results2)

pdf(paste0(path,"size-graphene-original.pdf"),w=5,h=4)
ggplot(data=results_size %>% mutate(type=ifelse(type=="Compact","Compact Blocks",type)))+# %>% filter(bound == 0.01))+# %>% filter(blkSize == 1000))+
  theme_grey(18)+
  aes(x=multiple, y=avg/1024, color=factor(type), group=factor(type))+
  geom_point()+
  geom_line()+
  xlab("Txns in mempool not in blk \n as a multiple of blk size")+
  ylab("Avg encoding size (KB)")+
  #geom_errorbar(aes(ymin=avg+var,ymax=avg-var))+#, width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=seq(0, 5, by=0.5))+
  scale_y_continuous(breaks=pretty_breaks(4))+
  facet_grid(blkSize ~ ., scales='free')+
  labs(color = "")+
  theme(legend.position="top",
        legend.box.margin = margin(-10,-10,-10,-10))
dev.off()

##################################################################################
# Analyze each piece of message
results.1 <- dat %>% group_by(TFPs, blkSize, bound) %>% summarize(size=mean(first_IBLT), type='first_IBLT')
results.2 <- dat %>% group_by(TFPs, blkSize, bound) %>% summarize(size=mean(first_BF), type='first_BF')
#results.3 <- dat %>% group_by(TFPs, fractionHas) %>% summarize(size=mean(second_IBLT), type='second_IBLT')
#results.4 <- dat %>% group_by(TFPs, fractionHas) %>% summarize(size=mean(second_BF), type='second_BF')
results.5 <- dat %>% group_by(TFPs, blkSize, bound) %>% summarize(size=mean(extra), type='extra')
results <- bind_rows(results.1, results.2, results.5)#results.3, results.4, 

ggplot(data=results)+
  aes(x=TFPs, y=(size))+
  geom_bar(stat = 'identity', aes(fill = type))+#, position = 'dodge')+
  xlab("Txns in mempool not in blk")+
  ylab("Avg encoding size by parts")+
  facet_grid(blkSize ~ bound, scales='free')+
  scale_x_continuous(breaks=seq(0, 4000, by=2000))

##################################################################################
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

