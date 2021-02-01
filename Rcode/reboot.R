library(dplyr)
library(readr)
library(ggplot2)
path <- "/Users/Pinar/git/2018-paper-graphene/sigcomm-paper/figures/"

dat <- read_delim("/Users/Pinar/git/2018-paper-graphene/basic_sim/data/reboot.csv",
                  delim = '\t',
                  col_names = c('TFPs', 'blkSize', 'bound', 'fraction', 'mempoolSize', 'fprS', 'realFprS', 'fprR', 'a', 'b', 'x_star', 'z', 'count', 'observed_FPs', 'decoded', 'pingpong_decode', 'graphene', 'first_IBLT', 'first_BF', 'second_IBLT', 'second_BF', 'extra', 'IBLT_rows_first', 'IBLT_rows_second', 'compact', 'third_BF', 'fprRSecond'),
                  col_types = cols(.default = col_double(), decoded = col_character(), pingpong_decode = col_character()))
dat$de <- dat$decoded =='True'
dat$dePing <- dat$pingpong_decode =='True'
dat$combine <- dat$de | dat$dePing

##################################################################################################
# Analyze how many graphene blocks decoded
r1 <- dat %>% group_by(blkSize, bound, fraction) %>% summarize(avg=mean(de), var=var(de), cnt=n(), ci=1.96*sqrt(avg*(1-avg)/n()), type='without')
dat$combine <- dat$de | dat$dePing
r2 <- dat %>% group_by(blkSize, bound, fraction) %>% summarize(avg=mean(combine), var=var(combine), cnt=n(), ci=1.96*sqrt(avg*(1-avg)/n()), type='with ping-pong')
results <- rbind(r1, r2)

pdf(paste0(path,"decode-rate-reboot.pdf"),w=5,h=4)
ggplot(data=results)+# %>% filter(bound == 0.1))+
  aes(x=fraction, y=avg, color=factor(type), group=factor(type))+
  geom_point(size=0.25)+
  geom_line()+
  xlab("Fraction of blk receiver has in mempool")+
  ylab("Fraction of Graphene blks decoded")+
  geom_errorbar(aes(ymin=avg+ci,ymax=avg-ci), width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=seq(0, 1, by=0.1))+
  #scale_y_continuous(breaks=pretty_breaks(20))+
  geom_hline(aes(yintercept=1-bound), linetype=2, color='black')+
  #facet_wrap(blkSize ~ bound, ncol=3,scales='free')
  facet_grid(blkSize ~ .)+
  #geom_hline(yintercept=.999, linetype=2)
  theme_gray(14)+
  theme(legend.position="top")+
  scale_color_discrete(name="Algorithm")
dev.off()

##################################################################################################
# Analyze avg size of graphene blks
#results1 <- dat %>% group_by(multiple, blkSize, bound) %>% summarize(avg=mean(compactBlk), var=var(compactBlk), cnt=n(), type='Compact')#, ci=1.96*sqrt(avg*(1-avg)/n()))
results1 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(avg=mean(graphene), var=var(graphene), cnt=n(), type='Graphene')#, ci=1.96*sqrt(avg*(1-avg)/n()))
results2 <- data.frame(blkSize=results1$blkSize, fraction=results1$fraction, bound=results1$bound, avg=results1$blkSize*(1-results1$fraction)*ifelse(results1$blkSize==200,1,3)+8*results1$blkSize, var=0, cnt=0, type='Compact Blocks')
results_size <- rbind(data.frame(results1), results2)

pdf(paste0(path,"size-reboot.pdf"),w=5,h=4)
ggplot(data=results_size)+# %>% filter(bound == 0.01))+# %>% filter(blkSize == 1000))+
  aes(x=fraction, y=avg/1024, color=factor(type), group=factor(type))+
  theme_gray(18)+
  geom_point()+
  geom_line()+
  xlab("Fraction of txns common to mempools")+
  ylab("Avg encoding size (KB)")+
  #geom_errorbar(aes(ymin=avg+var,ymax=avg-var))+#, width=0.025, alpha=0.5)+
  scale_x_continuous(breaks=seq(0, 1, by=0.2))+
  facet_grid(blkSize ~ ., scales='free')+
  labs(color = "")+
  theme(legend.position="top", 
        legend.box.margin = margin(-10,-10,-10,-10))
dev.off()

##################################################################################################
# Analyze each piece of message
results.1 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(first_IBLT), type='first_IBLT')
results.2 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(first_BF), type='first_BF')
results.3 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(second_IBLT), type='second_IBLT')
results.4 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(second_BF), type='second_BF')
results.5 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(third_BF), type='third_BF')
results.6 <- dat %>% group_by(blkSize, fraction, bound) %>% summarize(size=mean(extra), type='extra')
results <- bind_rows(results.1, results.2, results.3, results.4,results.5, results.6)
cb <- dat %>% group_by(blkSize, bound, fraction)  %>% summarize(size=mean(compact), type='Compact')#, ci=1.96*sqrt(avg*(1-avg)/n()))

pdf(paste0(path,"size-by-parts-reboot.pdf"),w=5,h=4)
ggplot(data=results)+# %>% filter (bound == 0.005))+
  aes(x=fraction, y=(size))+
  geom_bar(stat = 'identity', aes(fill = type))+#, position = 'dodge')+
  xlab("Fraction of blk in mempool")+
  ylab("Avg blk size by parts (bytes)")+
  #geom_line(data=cb, aes(color=type), linetype=2, color='black')+
  geom_line(aes(y=blkSize*(1-fraction)*ifelse(blkSize==200,1,3)+8*blkSize), linetype=2, color='black')+
  #facet_grid(bound ~ ., scales='free')+
  scale_x_continuous(breaks=seq(0, 1, by=0.1))+
  #scale_y_continuous(breaks=pretty_breaks(10))+
  theme_gray(14)+
  #labs(fill = "Msg")+
  scale_fill_discrete(name="Msg type",
                      breaks=c("extra", "first_BF", "first_IBLT", "second_BF", "second_IBLT"),
                      labels=c("getdata","BF 1","IBLT 1", "BF 2", "IBLT 2"))+
  theme(legend.position="top")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  facet_grid(blkSize ~ ., scales='free')
dev.off()
