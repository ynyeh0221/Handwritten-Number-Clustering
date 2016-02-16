require(MASS)
require(mvtnorm)

# Read handwritten digits data
myData=read.csv("C:\\Users\\Juliana Yeh\\Documents\\semeion.csv",header=FALSE)
# Build data matrix with (thresholded) pixel and label data
myX=data.matrix(myData[,1:256])
myLabel=apply(myData[,257:266],1,function(xx){
  return(which(xx=="1")-1)
})
#dev.new(width=20,height=20)#open a graph window
#par(mai=c(0.05,0.05,0.05,0.05),mfrow=c(2,4))#mai: margin
#mfrow=c(nrows, ncols)# to create a matrix of nrows x ncols plots that are filled in by row
# Examine a (random) few
#indices=sample(dim(myX)[1])
#continueLoop=TRUE
#counter=1
#while(continueLoop){
#  image(t(matrix(myX[indices[counter],],byrow=TRUE,16,16)[16:1,]),col=gray(0:1),axes=FALSE)
#  box()
#  # Should loop continue?
#  cat("\n","Quit? (enter any character)","\n") # prompt
#  myQ=scan(what=character(),n=1)
#  if(length(myQ)>0){
#    continueLoop=FALSE
#  }else{
#    counter=counter+1
#  }
#}
#Do PCA
# Build a low dimensional (principle components) representation
#xBar=colMeans(myX)
#xTilde=t(apply(myX,1,function(xx){
#  return(xx-xBar)
#}))
#apply: Returns a vector or array or list of values obtained by applying a function to margins of an array or matrix.
#apply(X, MARGIN, FUN, ...)
#MARGIN: a vector giving the subscripts which the function will be applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns. Where X has named dimnames, it can be a character vector selecting dimension names.
#MARGIN==1 (BY ROW), MARGIN==2 (BY COLUMN)
#mySVD=svd(xTilde)
#myPCs=xTilde%*%mySVD$v

#function of EM algorithm
EM_mainfunction=function(X, N, iterations)
{
  #Initialization
  prob=matrix(0,n,N)
  temp=cbind(1:n,class)
  prob[temp]=1;
  Nmatrix=apply(prob,2,function(xx){return(sum(xx))})
  pi=Nmatrix/n
  mu=km$centers
  Sigma=rep(list(diag(0,d,d)),N)
  Var=array(0,dim=c(d,d,N))
  for(i in 1:N)
  {
    temp = matrix(0,d,d)
    for(j in 1:n)
    {
      temp=temp+(myX[j,]-mu[i,])%*%t(myX[j,]-mu[i,])*prob[j,i]
    }
    Var[,,i]=temp/Nmatrix[i]
    eig=svd(Var[,,i])
    sig=1/(d-q)*sum(eig$d[(q+1):d])
    if(q!=0)
    {
      Wq=eig$v[,1:q]%*%diag(apply(as.matrix(eig$d[1:q]),1,function(xx){return(sqrt(xx-sig))}),q,q)
      Sigma[[i]]=Wq%*%t(Wq)+sig*diag(1,d,d)
    }
    else
    {
      Sigma[[i]]=sig*diag(1,d,d)
    }
  }
  old_log_like=sum(log(rowSums(prob)))
  parameters<-list(mu=mu,Sigma=Sigma,pi=pi)
  
  log_likelihood=list()
  
  #Repeat the E step and M step
  for(k in 1:iterations)
  {
    print(k)
    
    #E step
    prob=matrix(0,n,N)
    for(j in 1:N)
    {
      #It takes really a long time.
      #for(i in 1:n)
      #{
      #  prob[i,j]=parameters$pi[j]*1/sqrt(det(parameters$Sigma[[j]]))*exp(-1/2*t(X[i,]-parameters$mu[j,])%*%solve(parameters$Sigma[[j]])%*%(X[i,]-parameters$mu[j,]))
      #}
      prob[,j]=pi[j]*dmvnorm(myX, matrix(parameters$mu[j,]), matrix(unlist(parameters$Sigma[[j]]), ncol=d, byrow = TRUE))
    }
    
    log_like=sum(log(rowSums(prob)))
    log_likelihood<-c(log_likelihood,log_like)
    if(abs((log_like-old_log_like)/log_like)<log_threshold)
    {
      break;
    }
    old_log_like=log_like
    #prob=prob/rowSums(prob)
    prob=t(apply(prob,1,function(xx){return(xx/sum(xx))}))
    
    #M step
    Nmatrix=colSums(prob)
    
    #Update mu
    mu=t(prob)%*%myX/Nmatrix
    
    #Update pi
    pi=Nmatrix/n
    
    #Update Sigma
    Matrix=list()
    for(j in 1:N)
    {
      temp=matrix(0,d,d)
      for(i in 1:n)
      {
        temp=temp+1/Nmatrix[j]*prob[i,j]*(X[i,]-mu[j,])%*%t(X[i,]-mu[j,])
      }
      Matrix[[j]]<-temp
    }
    Sigma=rep(list(diag(d)),N)
    for(i in 1:N)
    {
      eig=svd(Matrix[[i]])
      sig=1/(d-q)*sum(eig$d[(q+1):d])
      if(q!=0)
      {
        Wq=eig$v[,1:q]%*%diag(apply(as.matrix(eig$d[1:q]),1,function(xx){return(sqrt(xx-sig))}),q,q)
        Sigma[[i]]=Wq%*%t(Wq)+sig*diag(1,d,d)
      }
      else
      {
        Sigma[[i]]=sig*diag(1,d,d)
      }
    }
    
    parameters<-list(mu=mu,Sigma=Sigma,pi=pi)
  }
  for(i in 1:n)
  {
    class[i]=which.max(prob[i,])
  }
  #Return results
  return(list(log_likelihood=log_likelihood, parameters=parameters))
}

#Update class
#class=apply(myX,1,function(xx){return(which.min(apply(mu,1,function(yy){return(sum((xx-yy)^2))})))})


#Parameters
K=10
q=6
d=dim(myX)[2]
n=dim(myX)[1]
log_threshold=0.00001
km=kmeans(myX, K, iter.max=100, nstart=30)
class=km$cluster

#Do EM algorithm
result<-EM_mainfunction(myX, K, 100)#N=10, iteration=100

#Draw graphes of log_likehood and AIC with q
log_like_list=matrix(unlist(result$log_likelihood),ncol=1,byrow=TRUE)
#final_mu=matrix(unlist(result$parameters$mu), ncol=d, byrow = TRUE)-xBar
par(mai=c(1,1,1,1),cex=0.9)
plot(log_like_list,xlab="Iterations",ylab="log_likelihood", main = paste("Log-Likelihood (q=", q, ")"))
plot(-2*log_like_list+d*q+1-q*(q-1)/2,xlab="Iterations",ylab="AIC", main = paste("AIC (q=", q, ")"))

#Visualization
dev.new(width=6, height=10)
par(mai=c(0,0,0,0),cex=0.8,mfrow=c(K,6))
x_k_new = array(0,dim=c(6,d,K))
for(i in 1:K)
{
  x_k_new[1,1:d,i]=result$parameters$mu[i,]
  x_k_new[2:6,1:d,i]=rmvnorm(n=5,result$parameters$mu[i,],matrix(unlist(result$parameters$Sigma[[i]]), ncol=d, byrow = TRUE))
}
for(i in 1:K)
{
  for(j in 1:6)
  {
    image(t(matrix(x_k_new[j,,i],byrow=TRUE,16,16)[16:1,]),col=gray(0:128/128),axes=FALSE)
    box()
  }
}

#Accuracy
majorityclass=matrix(0,K,1)
correctcount=matrix(0,K,1)
for(numberplus1 in 1:K)
{
  counts=matrix(0,K)
  for(j in 1:K)
  {
    for(i in 1:n)
    {
      if(class[i]==j && myLabel[i]==numberplus1-1)
      {
        counts[j]=counts[j]+1
      }
    }
  }
  correctcount[numberplus1]=counts[which.max(counts)]
  majorityclass[numberplus1]=which.max(counts)#The miscluster rate of each number 0,...,9
}
miscounts=matrix(0,K)
for(i in 1:K)
{
  mcount=0;
  for(k in 1:n)
  {
    if(myLabel[k]==i-1 && class[k]!=majorityclass[i])
    {
      mcount=mcount+1
    }
  }
  miscounts[i]=mcount
}
for(i in 1:K)
{
  #miscounts[i]=100*miscounts[i]/(length(myX[myLabel==i-1])/d)
  miscounts[i]=(1-correctcount[i]/(length(myX[myLabel==i-1])/d))*100
}
totalmis=sum(miscounts)/K
dev.new(width=6, height=10)
par(mai=c(0.3,0.3,0.3,0.3),mgp=c(2.4, 0.8, 0), las=1,cex=0.7)
barplot(t(miscounts),names.arg=c("0","1","2","3","4","5","6","7","8","9"),xlab="hand-written digits", ylab="miscategorization
        rate %", main=paste("Overall mis-categorization rate=", round(totalmis, digits=3), "% (q=",q,")"))
box()


