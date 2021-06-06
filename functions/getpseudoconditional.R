#-----------------------------------------------------------------------------------------------------------------------------------
# get pseudo conditional survival probabilities 
# t is the survival time
# d is the censoring indicator
# qt is a vector of time points that are used to divide the time interval
# output has subject id (id) and time points (s) and pseudo conditional survival probabilities (pseudost) for subject=id and at time s
#------------------------------------------------------------------------------------------------------------------------------------
getPseudoConditional <- function(t, d, qt){
  #browser()
  s <- c(0, qt)  
  n=length(t)
  ns=length(s)-1  # the number of intervals
  D <- do.call(cbind, lapply(1:ns, function(j)  (s[j] < t) * (t <= s[j+1]) * (d == 1)))
  R <- do.call(cbind, lapply(1:ns, function(j) ifelse(s[j] < t, 1, 0)))
  Delta<-do.call(cbind, lapply(1:ns, function(j) pmin(t,s[j+1])-s[j]))
  
  # format into long formate
  dd.tmp=cbind.data.frame(id=rep(1:n,ns),s=rep(c(0,qt[-length(qt)]), each=n), y=c(R*Delta),d=c(D))
  
  dd=dd.tmp[dd.tmp$y>0,]
  pseudost=rep(NA, nrow(dd))
  for (j in 1:ns){
    index= (dd$s==s[j])
    dds=dd[index,]
    pseudost[index]=pseudosurv(time=dds$y, event=dds$d, tmax=s[j+1]-s[j])$pseudo
    print(j)
  }
  dd$pseudost=pseudost  
  
  return(dd[,c(1,2,5)])
}


#--------------------------------------------------------------------------------------
# There are two example neural network models used in the paper.
# You should tunning the hyperparamters to find the best neural network for your own study.
#-------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# Example of a neural network with two hidden layers implemented in R keras in the paper
#-----------------------------------------------------------------------------------------
pseudoDNN.train <- function(x_train, y_train){
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 8, kernel_regularizer = regularizer_l2(0.0001), activation = "tanh",
                input_shape = dim(x_train)[[2]]) %>%
    layer_dense(units = 4, kernel_regularizer = regularizer_l2(0.01),
                activation = "tanh") %>%
    layer_dense(units = 1, activation='sigmoid')
  
  model %>% compile(
    optimizer = optimizer_adam(lr = 0.0025),
    loss = "mse",
    metrics = c("mae")
  )
  model %>% fit(x_train, y_train,
                epochs = 30, batch_size = 64,
                verbose = 0)
  
  model
}

#----------------------------------------------------------------------------------------
# Another example of a neural network with one hidden layer implemented in R keras in the paper
#-----------------------------------------------------------------------------------------
pseudoDNN.train <- function(x_train, y_train){
  # use selu instead of relu for some studies
  model <- keras_model_sequential() %>%
    layer_dense(units=16,  activation = "selu",bias_initializer = initializer_constant(0.0),
                input_shape = dim(x_train)[[2]]) %>%
    layer_dropout(rate = 0.2) %>%
    layer_dense(units = 1, activation='sigmoid')
  
  model %>% compile(
    optimizer = optimizer_rmsprop(lr = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  model %>% fit(x_train, y_train,
                epochs = 1000, batch_size =256,
                verbose = 0)
  model
}

#----------------------------------------------------------------------------
#prediction based on a keras model
#----------------------------------------------------------------------------
pseudoDNN.predict <- function(model, x_test){
  ypred <- model %>% predict(x_test)
  ypred
}

