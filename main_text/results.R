library(magrittr)
library(tidyverse)
library(ggpubr)
library(DescTools)


############
#main text

est<-readRDS("est.rds")%>%filter(!grepl("trans",method))
epsilon<-readRDS("epsilon.rds")%>%filter(!grepl("trans",method))
Qhat<-readRDS("Qhat.rds")%>%select(!contains("trans"))
est$method%<>%fct_relabel(~gsub("ee","dml",.))%>%fct_relabel(~gsub("dml_tr","dml_cl",.))

load("setup.RData")

#estimator
ggplot(est,aes(x=method,y=ATT.est))+geom_violin()+
    stat_summary(fun="mean",geom="point",alpha=.3)+stat_summary(fun="median",geom="point",shape=2)+
    geom_hline(yintercept=true.ATT,alpha=.2,linetype=2)+
    facet_wrap(~n,labeller="label_both")+ylab(expression(paste("Estimate of ATT ",theta)))+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))
ggsave("est.pdf",height=2,width=6)
(tab<-est%>%group_by(method,n)%>%summarize(bias=(mean(ATT.est)-true.ATT)*1e4,MSE=mean((ATT.est-true.ATT)^2)*1e4,coverage=mean(CIcover)))
est%>%filter(n==ns[1],method=="dml")%>%mutate(is.neg=EY0.est<0)%>%{BinomCI(x=sum(.$is.neg),n=nrow(.))}


#CI coverage
est%>%group_by(n,method)%>%summarize(coverage=mean(CIcover))


epsilon%>%filter(!grepl("p",method))%>%{table(abs(.$epsilon)>10,xor(.$anycase_train,.$anycase_test),deparse.level=2)}
epsilon%>%filter(grepl("p",method))%>%{table(abs(.$epsilon)>10,!(.$anycase_train | .$anycase_test),deparse.level=2)}
table(abs(epsilon$epsilon)>2,!(epsilon$anycase_train | epsilon$anycase_test),deparse.level=2)

View(est%>%filter(n==ns[1]))
View(epsilon%>%filter(n==ns[1]))
small.n.sim.id<-4
large.n.sim.id<-3

#epsilon
epsilon%>%filter(n==ns[2])%>%{summary(.$epsilon)}
epsilon%>%filter(n==ns[1],sim.id %in% small.n.sim.id)%>%View
p1<-ggplot(epsilon,aes(x=epsilon))+geom_histogram()+
    geom_vline(xintercept=0,alpha=.2,linetype=2)+
    facet_grid(method~n,scales="free")+
    xlab(expression(paste(epsilon,"*")))+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))
p2<-ggplot(epsilon,aes(x=tr(epsilon,-10,10)))+geom_histogram()+
    geom_vline(xintercept=0,alpha=.2,linetype=2)+
    facet_grid(method~n)+
    xlab(expression(paste("Clipped ",epsilon,"*")))+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))
ggarrange(p1,p2,ncol=2,labels="auto")
ggsave("epsilon.pdf",height=3.5,width=8)

#qq-plot of initial Q and Qhat
MRAD<-Qhat%>%pivot_longer(starts_with("tmle"),names_to="method",values_to="Qhat")%>%group_by(n,method,sim.id)%>%summarize(MRAD=mean(abs(Qhat-Q)/Qhat))%>%mutate(method=sub("_Qhat","",method))
MRAD<-left_join(MRAD,epsilon%>%select(n,method,sim.id,anycase_test,anycase_train),by=c("n","method","sim.id"))
p1<-ggplot(MRAD,aes(x=MRAD))+geom_histogram()+
    geom_vline(xintercept=0,alpha=.2,linetype=2)+
    facet_grid(method~n,scales="free")+
    xlab("MRAD")+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))
p2<-ggplot(MRAD,aes(x=tr(MRAD,0,10)))+geom_histogram()+
    geom_vline(xintercept=0,alpha=.2,linetype=2)+
    facet_grid(method~n)+
    xlab("Clipped MRAD")+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))
ggarrange(p1,p2,ncol=2,labels="auto")
ggsave("MRAD.pdf",height=3.5,width=8)

d<-Qhat%>%filter((sim.id %in% small.n.sim.id & n==ns[1]) | (sim.id %in% large.n.sim.id & n==ns[2]))%>%pivot_longer(starts_with("tmle"),names_to="method",values_to="Qhat")%>%rename(`Simulation run #`=sim.id)
d$method%<>%fct_relabel(~gsub("_Qhat","",.))
ggplot(d,aes(x=Q,y=Qhat))+geom_jitter(alpha=.01,,width=.002,height=.0005)+
    geom_abline(slope=1,intercept=0,alpha=.2,linetype=2,linewidth=1.5)+
    facet_wrap(n+`Simulation run #`+method~.,labeller="label_both",scales="free",ncol=4)+
    facet_grid(n~method,labeller="label_both",scales="free")+
    xlab(expression(paste("Initial ",hat(Q)[v])))+
    ylab(expression(paste("Targeted Fluctuation ",hat(Q)[v],"*")))+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))
ggsave("Q.pdf",height=3.5,width=8)




#cherry-picked TMLE performance
epsilon.include<-epsilon%>%filter(n==ns[1])%>%group_by(method,sim.id)%>%summarize(include=all(abs(epsilon)<=10))%>%filter(include)%>%select(!include)
to.include<-unique(epsilon.include$sim.id)
p1<-ggplot(est%>%filter(n==ns[1],(method %in% c("dml","dml_cl")) | (paste(method,sim.id) %in% (paste(epsilon.include$method,epsilon.include$sim.id)))),aes(x=method,y=ATT.est))+geom_violin()+
    stat_summary(fun="mean",geom="point",alpha=.3)+stat_summary(fun="median",geom="point",shape=2)+
    geom_hline(yintercept=true.ATT,alpha=.2,linetype=2)+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))+rremove("ylab")+rremove("xlab")
(tab<-est%>%filter(n==ns[1],sim.id %in% to.include)%>%group_by(method,n)%>%summarize(bias=(mean(ATT.est)-true.ATT)*1e4,MSE=mean((ATT.est-true.ATT)^2)*1e4,coverage=mean(CIcover)))

Q.include<-MRAD%>%filter(n==ns[1])%>%group_by(method,sim.id)%>%summarize(include=all(MRAD<=10))%>%filter(include)%>%mutate(method=sub("_Qhat","",method))
to.include<-unique(Q.include$sim.id)
p2<-ggplot(est%>%filter(n==ns[1],(method %in% c("dml","dml_cl")) | (paste(method,sim.id) %in% (paste(Q.include$method,Q.include$sim.id)))),aes(x=method,y=ATT.est))+geom_violin()+
    stat_summary(fun="mean",geom="point",alpha=.3)+stat_summary(fun="median",geom="point",shape=2)+
    geom_hline(yintercept=true.ATT,alpha=.2,linetype=2)+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))+rremove("ylab")+rremove("xlab")
(tab<-est%>%filter(n==ns[1],sim.id %in% to.include)%>%group_by(method,n)%>%summarize(bias=(mean(ATT.est)-true.ATT)*1e4,MSE=mean((ATT.est-true.ATT)^2)*1e4,coverage=mean(CIcover)))

p<-ggarrange(p1,p2,ncol=2,labels="auto")
annotate_figure(p,left=text_grob(expression(paste("Estimate of ATT ",theta)),rot=90,size=10),bottom=text_grob("method",size=10))
ggsave("check_diagnosis.pdf",height=2,width=6)


##############
#small bounded TMLE
est<-readRDS("est.rds")%>%filter(grepl("trans|ee",method))
epsilon<-readRDS("epsilon.rds")%>%filter(grepl("trans|ee",method))
Qhat<-readRDS("Qhat.rds")%>%select(contains("trans") | contains("ee"))
est$method%<>%fct_relabel(~gsub("ee","dml",.))%>%fct_relabel(~gsub("dml_tr","dml_cl",.))


load("setup.RData")

#estimator
ggplot(est,aes(x=method,y=ATT.est))+geom_violin()+
    stat_summary(fun="mean",geom="point",alpha=.3)+stat_summary(fun="median",geom="point",shape=2)+
    geom_hline(yintercept=true.ATT,alpha=.2,linetype=2)+
    facet_wrap(~n,labeller="label_both")+ylab(expression(paste("Estimate of ATT ",theta)))+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))
ggsave("est_small_bound.pdf",height=3,width=6)
(tab<-est%>%group_by(method,n)%>%summarize(bias=(mean(ATT.est)-true.ATT)*1e4,MSE=mean((ATT.est-true.ATT)^2)*1e4))



#############
#large bounded TMLE
est<-readRDS("large_Q_bound/est.rds")%>%filter(grepl("trans|ee",method))
epsilon<-readRDS("large_Q_bound/epsilon.rds")%>%filter(grepl("trans|ee",method))
Qhat<-readRDS("large_Q_bound/Qhat.rds")%>%select(contains("trans") | contains("ee"))
est$method%<>%fct_relabel(~gsub("ee","dml",.))%>%fct_relabel(~gsub("dml_tr","dml_cl",.))

load("large_Q_bound/setup.RData")

#estimator
ggplot(est,aes(x=method,y=ATT.est))+geom_violin()+
    stat_summary(fun="mean",geom="point",alpha=.3)+stat_summary(fun="median",geom="point",shape=2)+
    geom_hline(yintercept=true.ATT,alpha=.2,linetype=2)+
    facet_wrap(~n,labeller="label_both")+ylab(expression(paste("Estimate of ATT ",theta)))+
    theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5))
ggsave("est_large_bound.pdf",height=3,width=6)
(tab<-est%>%group_by(method,n)%>%summarize(bias=(mean(ATT.est)-true.ATT)*1e4,MSE=mean((ATT.est-true.ATT)^2)*1e4))
