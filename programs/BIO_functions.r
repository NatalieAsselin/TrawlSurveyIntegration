##Functions adapted from https://github.com/LobsterScience


sdmTMBcv_tntpreds <- function(object){
	#training and test predictions
  
	x  = object$data
	k = length(object$models)
	ou = list()
	for(i in 1:k){
		g = object$models[[i]]
		d = subset(x,!(cv_fold %in% i))
		e = subset(x,cv_fold %in% i)
		w <- predict(g, newdata = d,  
		             offset = if (!is.null(offset)) 
		               d$log.effort
		             else rep(0, nrow(d))) 
		d$pred = g$family$linkinv(w$est) ##same as exp but specific to the link function used
		d$tt = 'train'
		e$pred = e$cv_predicted
		e$tt = 'test'
		d = rbind(d,e)
		ou[[i]] = d
		}
	return(ou)
}

pbias<- function(x,y){
  sum(y-x)/sum(x) * 100  
}

mae<- function(x,y){
  sum(abs(x-y))/length(x)
}

rmse = function(x,y){
  sqrt((sum((y-x)^2))/length(x))
  
}

cv_SpaceTimeFolds = function(dat,idCol='IDS',k_folds=5){
  dat$fold_id = NA
  x = split(dat,f=dat[,idCol])
  o = list()
  s = 1:k_folds
  for(i in 1:length(x)){
    m = x[[i]]
    n = nrow(m)
    f=rep(s,each=floor(n/k_folds))
    if(length(f)!=n) f=c(f,s[1:(n-length(f))])
    m$fold_id = sample(f,size=n)
    o[[i]] = m
  }
  return(do.call(rbind,o))
}
