library(magrittr)
library(dplyr)
library(SuperLearner)

N<-2e2
ns<-c(3e2,2e3)
logit<-binomial()$linkfun
expit<-binomial()$linkinv

xgboost.learners<-create.Learner("SL.xgboost",tune=list(ntree=c(500,1e3),minobspernode=c(1,5)))
ranger.learners<-create.Learner("SL.ranger",tune=list(num.trees=c(500,1e3)))
nnet.learners<-create.Learner("SL.nnet",tune=list(size=c(25,50,100,200)))
SL.library<-c("SL.glm",xgboost.learners$names,ranger.learners$names,nnet.learners$names)

tr<-function(x,l=-1e4,u=1e4){
    pmax(pmin(x,u),l)
}

Q.true<-function(X1,X2,X3){
    expit(-4.64+(X1+X2+X3)/3)
}
g.true<-function(X1,X2,X3){
    expit(-1.4+.1*X1+.1*X2-.1*X3)
}
set.seed(9543)
d<-tibble(X1=runif(1e5,-1,1),X2=runif(1e5,-1,1),X3=runif(1e5,-1,1),g=g.true(X1,X2,X3),A=rbinom(1e5,1,g),Q=Q.true(X1,X2,X3))
truth<-d%>%filter(A==1)%>%summarize(x=mean(Q))%>%pull(x)
true.ATT<--truth
rm(d)

Q.bound<-.05

run.once<-function(n,sim.id){
    folds<-split(1:n,rep(1:2,each=n/2))
    data<-tibble(X1=runif(n,-1,1),X2=runif(n,-1,1),X3=runif(n,-1,1),A=rbinom(n,1,g.true(X1,X2,X3)),Y=ifelse(A==1,0,rbinom(n,1,Q.true(X1,X2,X3))),Y.trans=Y/Q.bound)
    
    Q<-lapply(folds,function(fold){
        model<-SuperLearner(data[-fold,]%>%filter(A==0)%>%pull(Y),data[-fold,]%>%filter(A==0)%>%select(X1:X3),family=binomial(),newX=data[fold,]%>%select(X1:X3),SL.library=SL.library)
        as.numeric(model$SL.predict)
    })%>%do.call(what=c)
    Q.trans<-lapply(folds,function(fold){
        d<-data[-fold,]
        X<-cbind(1,as.matrix(d%>%filter(A==0)%>%select(starts_with("X"))))
        optim.out<-optim(par=rep(0,4),fn=function(beta){
            pi<- expit(X%*%beta)
            pi[pi==0] <- .Machine$double.neg.eps
            pi[pi==1]<- 1-.Machine$double.neg.eps
            -sum(d$Y.trans[d$A==0]*log(pi)+(1-d$Y.trans[d$A==0])*log(1-pi))
        },gr=function(beta){
            pi<- expit(X%*%beta)
            pi[pi==0] <- .Machine$double.neg.eps
            pi[pi==1]<- 1-.Machine$double.neg.eps
            resid<-d$Y.trans[d$A==0]-pi
            -crossprod(X,resid)
        },method="BFGS")
        as.numeric(expit(cbind(1,as.matrix(data[fold,]%>%select(starts_with("X"))))%*%optim.out$par))
    })%>%do.call(what=c)
    g<-lapply(folds,function(fold){
        model<-SuperLearner(data[-fold,]%>%pull(A),data[-fold,]%>%select(X1:X3),family=binomial(),newX=data[fold,]%>%select(X1:X3),SL.library=SL.library)
        as.numeric(model$SL.predict)
    })%>%do.call(what=c)
    g.truncated<-any(g<.05 | g>.5)
    g<-tr(g,.05,.5)
    logitQ<-tr(logit(Q))
    logitQ.trans<-tr(logit(Q.trans))
    gamma<-lapply(folds,function(fold){
        rep(mean(data$A[fold]),length(fold))
    })%>%do.call(what=c)
    H<-g/(1-g)/gamma
    
    tmle_c.fold<-lapply(folds,function(fold){
        model<-glm(data$Y[fold]~H[fold]+0,family=binomial(),offset=logitQ[fold],subset=data$A[fold]==0)
        epsilon<-coef(model)
        Qhat<-expit(logitQ[fold]+epsilon*H[fold])
        
        list(Qhat=Qhat,epsilon=epsilon,anycase_train=any(data$Y[-fold]==1),anycase_test=any(data$Y[fold]==1))
    })
    tmle_c<-sapply(1:2,function(v){
        sum(data$A[folds[[v]]]*tmle_c.fold[[v]]$Qhat)/sum(data$A[folds[[v]]])
    })%>%mean
    tmle_c_epsilon<-sapply(tmle_c.fold,function(x){
        x$epsilon
    })
    tmle_c_anycase_train<-sapply(tmle_c.fold,function(x){
        x$anycase_train
    })
    tmle_c_anycase_test<-sapply(tmle_c.fold,function(x){
        x$anycase_test
    })
    tmle_c_Qhat<-lapply(tmle_c.fold,function(x){
        x$Qhat
    })%>%do.call(what=c)
    tmle_c_SE<-lapply(1:2,function(v){
        ifelse(data$A[folds[[v]]]==0,H[folds[[v]]]*(data$Y[folds[[v]]]-tmle_c.fold[[v]]$Qhat),(tmle_c.fold[[v]]$Qhat-tmle_c)/gamma[folds[[v]]])
    })%>%do.call(what=c)%>%{sqrt(mean(.^2)/n)} #because Y=0 when A=1, this is also the SE for ATT; CI for psi is equivalent to CI for ATT
    tmle_c_CI<-tmle_c+qnorm(c(.025,.975))*tmle_c_SE
    tmle_c_CIcover<-truth>=tmle_c_CI[1] && truth<=tmle_c_CI[2]
    
    tmle_c_trans.fold<-lapply(folds,function(fold){
        d<-data[fold,]
        optim.out<-optim(par=0,fn=function(ep){
            pi<-expit(logitQ.trans[fold][d$A==0]+ep*H[fold][d$A==0])
            pi[pi==0] <- .Machine$double.neg.eps
            pi[pi==1]<- 1-.Machine$double.neg.eps
            logLike<-sum(d$Y.trans[d$A==0]*log(pi)+(1-d$Y.trans[d$A==0])*log(1-pi))
            if(is.na(logLike) | is.nan(logLike) | is.infinite(logLike)){
                eps<-0
                logLike<-99
            }
            return(-logLike)
        },gr=function(ep){
            pi<-expit(logitQ.trans[fold][d$A==0]+ep*H[fold][d$A==0])
            pi[pi==0] <- .Machine$double.neg.eps
            pi[pi==1]<- 1-.Machine$double.neg.eps
            resid<-d$Y.trans[d$A==0]-pi
            gr<-crossprod(H[fold][d$A==0],resid)
            if( sum(is.na(gr))>0 | sum(is.nan(gr))>0 | sum(is.infinite(gr))>0 ){
                gr<-99
            }
            return(-gr)
        },method="BFGS")
        epsilon<-optim.out$par
        Qhat<-expit(logitQ.trans[fold]+epsilon*H[fold])*Q.bound
        list(Qhat=Qhat,epsilon=epsilon,anycase_train=any(data$Y[-fold]==1),anycase_test=any(data$Y[fold]==1))
    })
    tmle_c_trans<-sapply(1:2,function(v){
        sum(data$A[folds[[v]]]*tmle_c_trans.fold[[v]]$Qhat)/sum(data$A[folds[[v]]])
    })%>%mean
    tmle_c_trans_epsilon<-sapply(tmle_c_trans.fold,function(x){
        x$epsilon
    })
    tmle_c_trans_anycase_train<-sapply(tmle_c_trans.fold,function(x){
        x$anycase_train
    })
    tmle_c_trans_anycase_test<-sapply(tmle_c_trans.fold,function(x){
        x$anycase_test
    })
    tmle_c_trans_Qhat<-lapply(tmle_c_trans.fold,function(x){
        x$Qhat
    })%>%do.call(what=c)
    tmle_c_trans_SE<-lapply(1:2,function(v){
        ifelse(data$A[folds[[v]]]==0,H[folds[[v]]]*(data$Y[folds[[v]]]-tmle_c_trans.fold[[v]]$Qhat),(tmle_c_trans.fold[[v]]$Qhat-tmle_c_trans)/gamma[folds[[v]]])
    })%>%do.call(what=c)%>%{sqrt(mean(.^2)/n)}
    tmle_c_trans_CI<-tmle_c_trans+qnorm(c(.025,.975))*tmle_c_trans_SE
    tmle_c_trans_CIcover<-truth>=tmle_c_trans_CI[1] && truth<=tmle_c_trans_CI[2]
    
    tmle_w.fold<-lapply(folds,function(fold){
        model<-glm(data$Y[fold]~1,family=binomial(),offset=logitQ[fold],weights=H[fold]*(data$A[fold]==0))
        epsilon<-coef(model)
        Qhat<-expit(logitQ[fold]+epsilon)

        list(Qhat=Qhat,epsilon=epsilon,anycase_train=any(data$Y[-fold]==1),anycase_test=any(data$Y[fold]==1))
    })
    tmle_w<-sapply(1:2,function(v){
        sum(data$A[folds[[v]]]*tmle_w.fold[[v]]$Qhat)/sum(data$A[folds[[v]]])
    })%>%mean
    tmle_w_epsilon<-sapply(tmle_w.fold,function(x){
        x$epsilon
    })
    tmle_w_anycase_train<-sapply(tmle_w.fold,function(x){
        x$anycase_train
    })
    tmle_w_anycase_test<-sapply(tmle_w.fold,function(x){
        x$anycase_test
    })
    tmle_w_Qhat<-lapply(tmle_w.fold,function(x){
        x$Qhat
    })%>%do.call(what=c)
    tmle_w_SE<-lapply(1:2,function(v){
        ifelse(data$A[folds[[v]]]==0,H[folds[[v]]]*(data$Y[folds[[v]]]-tmle_w.fold[[v]]$Qhat),(tmle_w.fold[[v]]$Qhat-tmle_w)/gamma[folds[[v]]])
    })%>%do.call(what=c)%>%{sqrt(mean(.^2)/n)}
    tmle_w_CI<-tmle_w+qnorm(c(.025,.975))*tmle_w_SE
    tmle_w_CIcover<-truth>=tmle_w_CI[1] && truth<=tmle_w_CI[2]
    
    
    tmle_w_trans.fold<-lapply(folds,function(fold){
        d<-data[fold,]
        weights<-H[fold][d$A==0]
        optim.out<-optim(par=0,fn=function(ep){
            pi<-expit(logitQ.trans[fold][d$A==0]+ep)
            pi[pi==0] <- .Machine$double.neg.eps
            pi[pi==1]<- 1-.Machine$double.neg.eps
            logLike<-sum(weights*(d$Y.trans[d$A==0]*log(pi)+(1-d$Y.trans[d$A==0])*log(1-pi)))
            if(is.na(logLike) | is.nan(logLike) | is.infinite(logLike)){
                eps<-0
                logLike<-99
            }
            return(-logLike)
        },gr=function(ep){
            pi<-expit(logitQ.trans[fold][d$A==0]+ep)
            pi[pi==0] <- .Machine$double.neg.eps
            pi[pi==1]<- 1-.Machine$double.neg.eps
            resid<-weights*(d$Y.trans[d$A==0]-pi)
            gr<-sum(resid)
            if( sum(is.na(gr))>0 | sum(is.nan(gr))>0 | sum(is.infinite(gr))>0 ){
                gr<-99
            }
            return(-gr)
        },method="BFGS")
        epsilon<-optim.out$par
        Qhat<-expit(logitQ.trans[fold]+epsilon*H[fold])*Q.bound
        list(Qhat=Qhat,epsilon=epsilon,anycase_train=any(data$Y[-fold]==1),anycase_test=any(data$Y[fold]==1))
    })
    tmle_w_trans<-sapply(1:2,function(v){
        sum(data$A[folds[[v]]]*tmle_w_trans.fold[[v]]$Qhat)/sum(data$A[folds[[v]]])
    })%>%mean
    tmle_w_trans_epsilon<-sapply(tmle_w_trans.fold,function(x){
        x$epsilon
    })
    tmle_w_trans_anycase_train<-sapply(tmle_w_trans.fold,function(x){
        x$anycase_train
    })
    tmle_w_trans_anycase_test<-sapply(tmle_w_trans.fold,function(x){
        x$anycase_test
    })
    tmle_w_trans_Qhat<-lapply(tmle_w_trans.fold,function(x){
        x$Qhat
    })%>%do.call(what=c)
    tmle_w_trans_SE<-lapply(1:2,function(v){
        ifelse(data$A[folds[[v]]]==0,H[folds[[v]]]*(data$Y[folds[[v]]]-tmle_w_trans.fold[[v]]$Qhat),(tmle_w_trans.fold[[v]]$Qhat-tmle_w_trans)/gamma[folds[[v]]])
    })%>%do.call(what=c)%>%{sqrt(mean(.^2)/n)}
    tmle_w_trans_CI<-tmle_w_trans+qnorm(c(.025,.975))*tmle_w_trans_SE
    tmle_w_trans_CIcover<-truth>=tmle_w_trans_CI[1] && truth<=tmle_w_trans_CI[2]
    
    
    model<-glm(data$Y~H+0,family=binomial(),offset=logitQ,subset=data$A==0)
    epsilon<-coef(model)
    Qhat<-expit(logitQ+epsilon*H)
    tmle_cp<-sum(data$A*Qhat)/sum(data$A)
    tmle_cp_epsilon<-epsilon
    tmle_cp_anycase_train<-tmle_cp_anycase_test<-any(data$Y==1)
    tmle_cp_Qhat<-Qhat
    tmle_cp_SE<-ifelse(data$A==0,H*(data$Y-tmle_cp_Qhat),(tmle_cp_Qhat-tmle_cp)/gamma)%>%{sqrt(mean(.^2)/n)}
    tmle_cp_CI<-tmle_cp+qnorm(c(.025,.975))*tmle_cp_SE
    tmle_cp_CIcover<-truth>=tmle_cp_CI[1] && truth<=tmle_cp_CI[2]
    
    
    optim.out<-optim(par=0,fn=function(ep){
        pi<-expit(logitQ.trans[data$A==0]+ep*H[data$A==0])
        pi[pi==0] <- .Machine$double.neg.eps
        pi[pi==1]<- 1-.Machine$double.neg.eps
        logLike<-sum(data$Y.trans[data$A==0]*log(pi)+(1-data$Y.trans[data$A==0])*log(1-pi))
        if(is.na(logLike) | is.nan(logLike) | is.infinite(logLike)){
            eps<-0
            logLike<-99
        }
        return(-logLike)
    },gr=function(ep){
        pi<-expit(logitQ.trans[data$A==0]+ep*H[data$A==0])
        pi[pi==0] <- .Machine$double.neg.eps
        pi[pi==1]<- 1-.Machine$double.neg.eps
        resid<-data$Y.trans[data$A==0]-pi
        gr<-crossprod(H[data$A==0],resid)
        if( sum(is.na(gr))>0 | sum(is.nan(gr))>0 | sum(is.infinite(gr))>0 ){
            gr<-99
        }
        return(-gr)
    },method="BFGS")
    epsilon<-optim.out$par
    Qhat<-expit(logitQ.trans+epsilon*H)*Q.bound
    tmle_cp_trans<-sum(data$A*Qhat)/sum(data$A)
    tmle_cp_trans_epsilon<-epsilon
    tmle_cp_trans_anycase_train<-tmle_cp_trans_anycase_test<-any(data$Y==1)
    tmle_cp_trans_Qhat<-Qhat
    tmle_cp_trans_SE<-ifelse(data$A==0,H*(data$Y-tmle_cp_trans_Qhat),(tmle_cp_trans_Qhat-tmle_cp_trans)/gamma)%>%{sqrt(mean(.^2)/n)}
    tmle_cp_trans_CI<-tmle_cp_trans+qnorm(c(.025,.975))*tmle_cp_trans_SE
    tmle_cp_trans_CIcover<-truth>=tmle_cp_trans_CI[1] && truth<=tmle_cp_trans_CI[2]
    
    model<-glm(data$Y~1,family=binomial(),offset=logitQ,weights=H*(data$A==0))
    epsilon<-coef(model)
    Qhat<-expit(logitQ+epsilon)
    tmle_wp<-sum(data$A*Qhat)/sum(data$A)
    tmle_wp_epsilon<-epsilon
    tmle_wp_anycase_train<-tmle_wp_anycase_test<-any(data$Y==1)
    tmle_wp_Qhat<-Qhat
    tmle_wp_SE<-ifelse(data$A==0,H*(data$Y-tmle_wp_Qhat),(tmle_wp_Qhat-tmle_wp)/gamma)%>%{sqrt(mean(.^2)/n)}
    tmle_wp_CI<-tmle_wp+qnorm(c(.025,.975))*tmle_wp_SE
    tmle_wp_CIcover<-truth>=tmle_cp_CI[1] && truth<=tmle_cp_CI[2]
    
    
    weights<-H[data$A==0]
    optim.out<-optim(par=0,fn=function(ep){
        pi<-expit(logitQ.trans[data$A==0]+ep)
        pi[pi==0] <- .Machine$double.neg.eps
        pi[pi==1]<- 1-.Machine$double.neg.eps
        logLike<-sum(weights*(data$Y.trans[data$A==0]*log(pi)+(1-data$Y.trans[data$A==0])*log(1-pi)))
        if(is.na(logLike) | is.nan(logLike) | is.infinite(logLike)){
            eps<-0
            logLike<-99
        }
        return(-logLike)
    },gr=function(ep){
        pi<-expit(logitQ.trans[data$A==0]+ep)
        pi[pi==0] <- .Machine$double.neg.eps
        pi[pi==1]<- 1-.Machine$double.neg.eps
        resid<-weights*(data$Y.trans[data$A==0]-pi)
        gr<-sum(resid)
        if( sum(is.na(gr))>0 | sum(is.nan(gr))>0 | sum(is.infinite(gr))>0 ){
            gr<-99
        }
        return(-gr)
    },method="BFGS")
    epsilon<-optim.out$par
    Qhat<-expit(logitQ.trans+epsilon)*Q.bound
    tmle_wp_trans<-sum(data$A*Qhat)/sum(data$A)
    tmle_wp_trans_epsilon<-epsilon
    tmle_wp_trans_anycase_train<-tmle_wp_trans_anycase_test<-any(data$Y==1)
    tmle_wp_trans_Qhat<-Qhat
    tmle_wp_trans_SE<-ifelse(data$A==0,H*(data$Y-tmle_wp_trans_Qhat),(tmle_wp_trans_Qhat-tmle_wp_trans)/gamma)%>%{sqrt(mean(.^2)/n)}
    tmle_wp_trans_CI<-tmle_wp_trans+qnorm(c(.025,.975))*tmle_wp_trans_SE
    tmle_wp_trans_CIcover<-truth>=tmle_cp_CI[1] && truth<=tmle_cp_CI[2]
    
    
    plug.in.terms<-lapply(folds,function(fold){
        ifelse(data$A[fold]==1,Q[fold]/gamma[fold],0)
    })%>%do.call(what=c)
    ee.correction<-lapply(folds,function(fold){
        ifelse(data$A[fold]==0,(data$Y[fold]-Q[fold])*H[fold],0)
    })%>%do.call(what=c)
    ee<-mean(plug.in.terms)+mean(ee.correction)
    ee.IF<-lapply(folds,function(fold){
        ifelse(data$A[fold]==0,(data$Y[fold]-Q[fold])*H[fold],(Q[fold]-ee)/gamma[fold])
    })%>%do.call(what=c)
    ee_SE<-sqrt(mean(ee.IF^2)/n)
    ee_CI<-ee+qnorm(c(.025,.975))*ee_SE
    ee_CIcover<-truth>=ee_CI[1] && truth<=ee_CI[2]
    
    ee_tr<-tr(ee,0,1)
    ee_tr_CI<-tr(ee_CI,0,1)
    ee_tr_CIcover<-truth>=ee_tr_CI[1] && truth<=ee_tr_CI[2]
    
    anycase<-any(data$Y==1)
    
    EY1<-sum(data$A*data$Y)/sum(data$A)
    
    list(est=tibble(n=n,sim.id=sim.id,
                    method=c("tmle_c","tmle_c_trans","tmle_cp","tmle_cp_trans","tmle_w","tmle_w_trans","tmle_wp","tmle_wp_trans","ee","ee_tr"),
                    EY0.est=c(tmle_c,tmle_c_trans,tmle_cp,tmle_cp_trans,tmle_w,tmle_w_trans,tmle_wp,tmle_wp_trans,ee,ee_tr),
                    ATT.est=EY1-EY0.est,
                    SE=c(tmle_c_SE,tmle_c_trans_SE,tmle_cp_SE,tmle_cp_trans_SE,tmle_w_SE,tmle_w_trans_SE,tmle_wp_SE,tmle_wp_trans_SE,ee_SE,ee_SE),
                    CI.l=c(tmle_c_CI[1],tmle_c_trans_CI[1],tmle_cp_CI[1],tmle_cp_trans_CI[1],tmle_w_CI[1],tmle_w_trans_CI[1],tmle_wp_CI[1],tmle_wp_trans_CI[1],ee_CI[1],ee_tr_CI[1]),
                    CI.u=c(tmle_c_CI[2],tmle_c_trans_CI[2],tmle_cp_CI[2],tmle_cp_trans_CI[2],tmle_w_CI[2],tmle_w_trans_CI[2],tmle_wp_CI[2],tmle_wp_trans_CI[2],ee_CI[2],ee_tr_CI[2]),
                    CIcover=c(tmle_c_CIcover,tmle_c_trans_CIcover,tmle_cp_CIcover,tmle_cp_trans_CIcover,tmle_w_CIcover,tmle_w_trans_CIcover,tmle_wp_CIcover,tmle_wp_trans_CIcover,ee_CIcover,ee_tr_CIcover),
                    g.truncated=g.truncated,anycase=anycase),
         epsilon=tibble(n=n,sim.id=sim.id,
                        method=c(rep(c("tmle_c","tmle_c_trans","tmle_w","tmle_w_trans"),each=2),"tmle_cp","tmle_cp_trans","tmle_wp","tmle_wp_trans"),
                        epsilon=c(tmle_c_epsilon,tmle_c_trans_epsilon,tmle_w_epsilon,tmle_w_trans_epsilon,tmle_cp_epsilon,tmle_cp_trans_epsilon,tmle_wp_epsilon,tmle_wp_trans_epsilon),
                        anycase_train=c(tmle_c_anycase_train,tmle_c_trans_anycase_train,tmle_w_anycase_train,tmle_w_trans_anycase_train,tmle_cp_anycase_train,tmle_cp_trans_anycase_train,tmle_wp_anycase_train,tmle_wp_trans_anycase_train),
                        anycase_test=c(tmle_c_anycase_test,tmle_c_trans_anycase_test,tmle_w_anycase_test,tmle_w_trans_anycase_test,tmle_cp_anycase_test,tmle_cp_trans_anycase_test,tmle_wp_anycase_test,tmle_wp_trans_anycase_test)),
         Qhat=tibble(n=n,sim.id=sim.id,Q=Q,tmle_c_Qhat,tmle_c_trans_Qhat,tmle_cp_Qhat,tmle_cp_trans_Qhat,tmle_w_Qhat,tmle_w_trans_Qhat,tmle_wp_Qhat,tmle_wp_trans_Qhat)
    )
}

save.image("setup.RData")
