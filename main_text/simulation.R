library(magrittr)
library(dplyr)
library(SuperLearner)

load("setup.RData")

library(batchtools)
# options(batchtools.progress=FALSE)

# reg<-loadRegistry("simulation",writeable=TRUE)
reg<-makeRegistry(file.dir="simulation",seed=942,packages=c("magrittr","dplyr","SuperLearner"))
batchExport(mget(ls()))

batchMap(run.once,args=expand.grid(n=ns,sim.id=seq_len(N)))
submitJobs(getJobTable()[,chunk:=1],resources=list(chunks.as.arrayjobs=TRUE))
waitForJobs()
reduceResults(function(aggr,res) c(aggr,list(res$est)),init=list(),findDone())%>%bind_rows%>%saveRDS("est.rds")
reduceResults(function(aggr,res) c(aggr,list(res$epsilon)),init=list(),findDone())%>%bind_rows%>%saveRDS("epsilon.rds")
reduceResults(function(aggr,res) c(aggr,list(res$Qhat)),init=list(),findDone())%>%bind_rows%>%saveRDS("Qhat.rds")
