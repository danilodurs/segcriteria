#' log-transformation function for segmentation
#'
#' @param rna.data : a vector of counts from RNA-sequencing
#'
#' @return a vector of log-transformed data
#'
#'
#' @examples
#' log.data <- log.transform(dataset1)
#' plot(dataset1, type="l")
#' plot(log.data, type="l")
log.transform <- function(rna.data)
{
  log.data <- log(rna.data+1)
  n <- length(rna.data)
  ##transformation "HALL"
  x   <- log.data
  wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
  mat <- wei %*% t(x)
  mat[2, -n] = mat[2, -1]
  mat[3, -c(n-1, n)] = mat[3, -c(1, 2)]
  mat[4, -c(n-2, n-1, n)] = mat[4, -c(1, 2, 3)]
  est.sd <- sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3))
  return(log.data/est.sd)
}



#' Segment mean per data
#'
#' @param tau a vector of changepoints
#' @param rna.data a vector of counts
#'
#' @return a vector of means per data
#'
#' @examples
#' library(robseg)
#' log.trans <- log.transform(dataset1)
#' seg_rob  <- Rob_seg.std(x = log.trans,  loss = "Outlier", lambda = beta, lthreshold=3)
#' tau       <- seg_rob$t.est
#' av.line   <- mean.seg(rna.data = data, tau = tau)
#' plot(av.line)
#'
mean.seg <- function(tau,rna.data)
{

  tau      <- c(0,tau)
  all.mean <- NULL
  for ( i in 2:length(tau))#### pour chaque segment
  {
    seg      <- rna.data[(tau[i-1]+1):tau[i]]####on recupere les donnees du segments
    mymean   <- mean(seg) #####on fait la moyenne du segment
    meanseg  <- rep(mymean,length(seg)) ##average for each data of the segment
    all.mean <- c(all.mean,meanseg)
  }
  cat(".")
  return(all.mean)
}

#' Computing the MSE from a pair of datasets, depending on a given penalty
#'
#' @param rna.data1 a vector of data
#' @param rna.data2 a vector of data, same length as rna.data1
#' @param beta a value of penalty
#'
#' @return the value of MSE from two datasets by a given penalty
#'
#' @examples
#' log.dataset1 <- log.transform(dataset1)
#' log.dataset2 <- log.transform(dataset2)
#' mse.val <- mse.func(rna.data1=log.dataset1,rna.data2=log.dataset2,beta=25*log(length(dataset1)))
#'
#'
mse.func <- function(rna.data1,rna.data2,beta)
{
  n          <- length(rna.data1)
  seg_rob1   <- Rob_seg.std(x = rna.data1,  loss = "Outlier", lambda = beta, lthreshold=3)
  listoftau1 <- seg_rob1$t.est
  seg_rob2   <- Rob_seg.std(x = rna.data2,  loss = "Outlier", lambda = beta, lthreshold=3)
  listoftau2 <- seg_rob2$t.est
  meanset1   <- mean.seg( tau = seg_rob1$t.est , rna.data = rna.data1 )
  meanset2   <- mean.seg( tau = seg_rob2$t.est , rna.data = rna.data2 )
  mse        <- (sum((meanset1-meanset2)^2)/n)
  cat('*')
  return(mse)
}


#' Segmentation selection on MSE criterion
#'
#' @param list.rna.data a list of vector of dataset
#' @param penalty_range a vector of length 2, respectively defining the minimum penalty and the maximum penalty
#' @return A list with a vector of penalties and a matrix of mse
#' @examples
#' l.d1    <- log.transform(dataset1)
#' l.d2    <- log.transform(dataset2)
#' l.d3    <- log.transform(dataset3)
#' l.data  <- list(l.d1,l.d2,l.d3)
#' mse.res <- mse.penalties(list.rna.data=l.data, penalty_range=c(10,100))
#'
#'
mse.penalties <- function(list.rna.data , penalty_range)
{
  crops          <- NULL
  n              <- length(list.rna.data[[1]])
  for(rna.data in list.rna.data)
  {
    crop  <- CROPS.RFPOP(data = rna.data , min_pen = penalty_range[1]*log(n) , max_pen = penalty_range[2]*log(n))
    crops <- c(crops,crop[[1]][2,])
    #save all intermediate penalties from all datasets
  }
  intermediate.penalties  <- sort(crops)
  doublons       <- which(duplicated(intermediate.penalties))
  intermediate.penalties  <- intermediate.penalties[-doublons]
  ##calcul des taus estimés et des mse
  dataset.pairs <- combn(list.rna.data,2,simplify = F)
  mat.results <- NULL
  for(dsp in dataset.pairs)
  {
    row.mse <- NULL
    for(beta in intermediate.penalties)
    {
      row.mse <- c(row.mse,mse.func(rna.data1 = dsp[[1]],rna.data2 = dsp[[2]],beta = beta))
      cat('+')
    }
    mat.results <- rbind(mat.results,row.mse)
    cat('|')
  }
  i              <- 1
  row.mse <- NULL

  row.mse <- colMeans(mat.results)
  #while(i <= length(mat.results[2,]))
  #{
  #  row.mse <- c(row.mse,mean(mat.results[,i]))
  #  i <- i+1
  #}
  mat.results <- rbind(mat.results,row.mse)
  cat('|')
  list.results <- list(intermediate.penalties,mat.results)

  class(list.results) <- "MSE"

  return(list.results)

}

#' NID computed from two changepoint datasets
#'
#' @param rna.data1 a vector of data
#' @param rna.data2 a vector of data
#' @param beta penalty value
#'
#' @return the NID from the 2 segmentations, a decimal value between 0 and 1.
#'
#' @examples
#' log.dataset1 <- log.transform(dataset1)
#' log.dataset2 <- log.transform(dataset2)
#' nid.val <- nid.func(rna.data1=log.dataset1,rna.data2=log.dataset2,beta=25*log(length(dataset1)))
#'
nid.func <- function(rna.data1,rna.data2,beta)
{
  seg_rob1      <- Rob_seg.std(x = rna.data1,  loss = "Outlier", lambda = beta, lthreshold=3)
  seg_rob2      <- Rob_seg.std(x = rna.data2,  loss = "Outlier", lambda = beta, lthreshold=3)
  seg1          <- seg_rob1$t.est
  seg2          <- seg_rob2$t.est
  seg1          <- c(0,seg1)
  seg2          <- c(0,seg2)
  class.data.1 <- NULL
  class.data.2 <- NULL
  for(i in 2:length(seg1))##attribut l'indice du segment à chaque position
  {
    class.data.1 <- c(class.data.1,rep(i-1,seg1[i]-seg1[i-1]))
  }
  for(j in 2:length(seg2))##attribut l'indice du segment à chaque position
  {
    class.data.2 <- c(class.data.2,rep(j-1,seg2[j]-seg2[j-1]))
  }
  nid.out12 <- NID(class.data.1,class.data.2)
  cat("|")
  return(nid.out12)
}

#' Segmentation penalty criterion by NID
#'
#' @param penalty_range a vector of length 2 respectively defining the minimum penalty and the maximum penalty
#' @param list.rna.data a list of datasets with the same length
#'
#' @return A list with a vector of penalties and a matrix of nid
#' @examples
#' l.d1    <- log.transform(dataset1)
#' l.d2    <- log.transform(dataset2)
#' l.d3    <- log.transform(dataset3)
#' l.data  <- list(l.d1,l.d2,l.d3)
#' nid.res <- nid.penalties(list.rna.data=l.data, penalty_range=c(10,100))
#'
nid.penalties <- function(list.rna.data, penalty_range)
{
  n              <- length(list.rna.data[[1]])
  crops          <- NULL
  for(rna.data in list.rna.data)
  {
    crop  <- CROPS.RFPOP(data = rna.data , min_pen = penalty_range[1]*log(n) , max_pen = penalty_range[2]*log(n))
    crops <- c(crops,crop[[1]][2,])
    #save all intermediate penalties from all datasets
  }
  intermediate.penalties  <- sort(crops)
  doublons       <- which(duplicated(intermediate.penalties))
  intermediate.penalties  <- intermediate.penalties[-doublons]
  ##calcul des taus estimés
  dataset.pairs <- combn(list.rna.data,2,simplify = F)
  mat.results <- NULL
  for(dsp in dataset.pairs)
  {
    row.nid <- NULL
    for(beta in intermediate.penalties)
    {
      row.nid <- c(row.nid,nid.func(rna.data1 = dsp[[1]],rna.data2 = dsp[[2]],beta = beta))
      cat('+')
    }
    mat.results <- rbind(mat.results,row.nid)
    cat('|')
  }
  i              <- 1
  row.nid <- NULL
  while(i <= length(mat.results[2,]))
  {
    row.nid <- c(row.nid,mean(mat.results[,i]))
    i <- i+1
  }
  mat.results <- rbind(mat.results,row.nid)
  cat('|')
  list.results <- list(intermediate.penalties,mat.results)
  class(list.results) <- "NID"
  return(list.results)
}


#' Average count per segment for each replicat
#'
#' @param list.tau A list of changepoints, for each dataset, for a given penalty.
#' @param list.rna.data A list of dataset (RNA countings)
#'
#' @return a matrix of average counts per segmentation per replica. The last line is the actual average counts per segment.
#'
#' @examples
#' library(robseg)
#' log.d1        <- log.transform(dataset1)
#' log.d2        <- log.transform(dataset2)
#' my.list.data  <- list(log.d1,log.d2)
#' rob_seg1      <- Rob_seg.std( x = log.d1 , loss = "Outlier", lambda = 25*log(length(dataset1)), lthreshold=3)
#' rob_seg2      <- Rob_seg.std( x = log.d2 , loss = "Outlier", lambda = 25*log(length(dataset1)), lthreshold=3)
#' my.list.tau   <- list(rob_seg1$t.est,rob_seg2$t.est)
#' cps      <- counts.per.seg(list.tau = my.list.tau,list.rna.data = my.list.data)
counts.per.seg <- function(list.tau,list.rna.data)
{
  ##"Union" of all changepoins
  vect.tau <- unlist(list.tau)
  vect.tau <- sort(vect.tau[!duplicated(vect.tau)])
  average.counts    <- NULL
  matrix.of.results <- NULL
  #Adding the minimum changepoint
  vect.tau <- c(0,vect.tau)
  #For each segment
  i<- 1
  #For each dataset
  for(dataset in list.rna.data)
  {
    row.segment <- NULL
    while(i < length(vect.tau)){
      #Defining a segment( a vector of data delimited by taus)
      curr.segment  <- dataset[(vect.tau[i]+1):vect.tau[i+1]]
      row.segment <- c(row.segment, mean(curr.segment))
      i <- i+1
    }
    #for each dataset, adding the set of segments
    matrix.of.results <- rbind(matrix.of.results, row.segment)
  }

  #Adding the row of results
  matrix.of.results <- rbind(cmatrix.of.results, colMeans(matrix.of.results))
  return(matrix.of.results)
}

#' Plotting MSE object
#'
#' @param mse : the MSE objet returned by the function mse.penalties
#'
#' @return the graph according to the parameters entered in the function
#'
#' @examples
#' mse.beta <- mse.penalties(list(log.transform(dataset1),log.transform(dataset2),log.transform(dataset3)), penalty_range = c(15,75) )
#' plot(mse.beta)
#'
plot.MSE <- function(mses)
{
  plot(y=mses[[2]][4,], x=mses[[1]], type = "s",col="orange",xlab = "penalty value",ylab = "MSE")
}

#' Plotting NID object
#'
#' @param nid : the NID objet returned by the function nid.penalties
#'
#' @return the graph according to the parameters entered in the function
#'
#' @examples
#' nid.beta <- nid.penalties(rna.data1 = dataset1 , rna.data2 = dataset2, penalty_range = c(15,75) )
#' plot(nid.beta)
#'
plot.NID <- function(nids)
{
  plot(y=nids[[2]][4,], x=nids[[1]], type = "s",col="orange",xlab = "penalty value",ylab = "NID")
}
