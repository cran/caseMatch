## This is the function from case.match_sent to Braun_5apr2013.R
## copied directly.  Then I added weights to do weighted mahalanobis
## matching.  Weighting formulat from: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3760213/

###########################################
## Wrap the case-matcher as a function
###########################################

case.match <- function(data, id.var, case.N=2, 
                        distance="mahalanobis",
                        design.type="most similar",
                        match.case=NULL,
                        number.of.matches.to.return=1,
                        treatment.var=NULL,
                        leaveout.vars=NULL,
                        max.variance=FALSE,
                        variance.tolerance=.1,
                        max.spread=FALSE,
                        varweights=NULL){
       ## "design.type" can be "most similar" or "most different".
       ## "id.var" must be the name of the id variable in quotes
       ## "distance" can be "euclidean", "standardized", or "mahalanobis"

    ## Error messages
  if (case.N<2) {
      stop("case.N must be an integer 2 or larger", call. = FALSE) }

    ## Warns you if you have too many treatment variables specified
  if(length(treatment.var)>1){
    stop("Only one treatment variable can be specified.")}

    ## warns you if the treatment variable is missing
  if(max(names(data)==treatment.var, na.omit=T)==F) {
      stop("Treatment variable is missing from the dataset.")}

    ## warns you if one or more of the leaveout vars isn't in the dataset
  if(class(leaveout.vars)=="character"){
    for(i in 1:length(leaveout.vars)){
      if(max(names(data)==leaveout.vars[i], na.omit=T)==F) {
        stop("At least 1 leaveout.vars variable is missing from the dataset.") 
      }
    }
  }  
  
  if(max.variance==T & max.spread==T){
    stop("Cannot specify max.variance==TRUE & max.spread==TRUE.  Pick one or the other.")}





    ## setup stuff

    ## pull out the unit names
  unit.names <- as.vector(data[,which(names(data)==id.var)])

    ## this first part deals with it if it is numeric
  if(class(leaveout.vars)=="numeric") { 
    leaveout <- leaveout.vars
  }

    ## if it is character
      ## this uses charmatch() to find the column numbers of the 
      ## user specified variable names
  if(class(leaveout.vars)=="character"){
    leaveout <- charmatch(leaveout.vars,names(data))
  }

    ## pull out a matrix with just the variables to use for matching
    ## if there are no variables to leave out and no treatment var
  if(is.null(leaveout.vars)==T & is.null(treatment.var)==T){
    X <- data[,-which(names(data)==id.var)]}
    ## if there are vars to leave out but no treatment var
  if(is.null(leaveout.vars)==F & is.null(treatment.var)==T){
    X <- data[,-c(which(names(data)==id.var), leaveout)]}
    ## if there are no vars to leave out but one treatment var
  if(is.null(leaveout.vars)==T & is.null(treatment.var)==F){
    X <- data[,-c(which(names(data)==id.var), 
                  which(names(data)==treatment.var))]
    treat <- data[,which(names(data)==treatment.var)]
  }
    ## if there are vars to leave out and one treatment var
  if(is.null(leaveout.vars)==F & is.null(treatment.var)==F){
    X <- data[,-c(which(names(data)==id.var),
                  leaveout,
                  which(names(data)==treatment.var))]
    treat <- data[,which(names(data)==treatment.var)]
  }

  ## make sure the variable weights match the dim of X
  if(!is.null(varweights)){
    if(length(varweights) != ncol(X)){stop("Number of variable weights does not match number of variables.")}
  }
  ## make sure that variable weight names match the names of X
  if(!is.null(varweights)){
    if(!is.null(names(varweights))){
      if(sum(names(varweights) != colnames(X))>0){stop("Names of variable weights do not match names of variables.")}
    }
  }
  

    ## Sets N as the overall number of observations
  N <- nrow(data)
    ## creates a new id variable to use throughout the algorith
  id <- seq(1,N,1)
    ## create the covariance matrix of the whole data for later
  if(distance=="mahalanobis") {covX <- cov(X)}
  if(distance=="euclidean"){covX <- diag(ncol(X))}
    ## standardize the data by its standard dev. if specified.
  if(distance=="standardized") {
    stddev <- apply(X,MARGIN=2,FUN=sd)
    X <- t(apply(X,MARGIN=1,FUN=function(x){x/stddev}))
  }

    ## other messages
  if(is.null(match.case)==T) {
    print(paste("There are", choose(N, case.N), "possible case combinations"))
  }

  if(is.null(match.case)==F) {
    print(paste("There are", choose((N-1), (case.N-1)), "possible case combinations"))
  }

    ## set the memory to the maximum allowed to allow space for the permutations
  #max.mem <- round(memory.limit(), 2)
  #memory.limit(max.mem)
  #round(memory.limit(), 2)

    ## if there is a specific case to match:
  if(is.null(match.case)==F) { 
      ## Issue a warning if the name of the unit to match doesn't exist
    if(length(grep(match.case,unit.names, fixed=T))==0) {
      stop("match.case does not match any of the unit IDs", call. = FALSE) }
      ## This pulls out the id of the unit to match
    match.case.id <- id[grep(match.case,unit.names, fixed=T)]
      ## make a new id var that has the unit to match having id=1
    id2 <- rep(NA,length(id))
    id2[id==match.case.id] <- 1
    id2[id!=match.case.id] <- seq(2,N,1)
    id <- id2
    match.case.id <- 1
      ## I reorder the dataset so that the match.case is first
    X <- X[order(id),]
      ## reorder the unit names and the id and treat variables as well
    unit.names <- unit.names[order(id)]
    if(is.null(treatment.var)==F){treat <- treat[order(id)]}
    id <- sort(id)
  
      ## If there is no treat var, I do all the possible combinations:
    if(is.null(treatment.var)==T){
        ## combn() lists all possible combinations
        ## each column is a combination
        ## a different set of combinations -- only those that include the match.case
      combin <- t(combn(id[-match.case.id],(case.N-1)))
      combinations <- cbind(rep(match.case.id,nrow(combin)), combin)
    }
    ## If there is a treat var, I do only the combinations of the other treatments:
    if(is.null(treatment.var)==F){
        ## a different set of combinations -- only those that include the match.case
        ## using only the units with different treatments
      treat.value <- treat[id==match.case.id]
      combin <- t(combn(id[-unique(c(match.case.id,which(treat==treat.value)))],(case.N-1)))  ## note, I already took out 
      combinations <- cbind(rep(match.case.id,nrow(combin)), combin)
    }
  }

     ## If there isn't a case to match, then do this:
  if(is.null(match.case)==T) { 
    combinations <- t(combn(id,case.N))
  }

    
    ## make a holder matrix to hold the distances between each unit
  holder <- matrix(NA,nrow(X),nrow(X))

    ## weights for variables
  if(is.null(varweights)){varweights <- rep(1,ncol(X))}
  W <- diag(varweights)
  invCovW <- W %*% solve(covX) %*% W
    ## in the genmatch paper in Review of Econ. and Stat., Jas has a different version of weights that gives a different answer!
  #jasW <- t(chol(covX)) %*% W %*% (chol(covX))

    ## fill the holder with the distances between each unit
  for(i in 1:(nrow(X)-1)){
    vec1 <- X[i,]  ## The first vector
    for(j in (i+1):nrow(X)){  ## compare to evey row below the current row
      vec2 <- X[j,]
         ## Euclidean distance (using an identity matrix)
      if(distance=="euclidean" | distance=="standardized"){
        #holder[j,i] <- mahalanobis(as.numeric(vec1),as.numeric(vec2),diag(length(vec1)))
        ## as of 26jan2014 I allow weights.  This means that the code for euclidean 
        ## is now the same as for mahal (with different "invCovW".  I could combine.
        holder[j,i] <- mahalanobis(as.numeric(vec1),as.numeric(vec2),invCovW,inverted = T)
      }
         ## Mahalanobis distance
      if(distance=="mahalanobis"){
        ## As of 26jan2014 I allow weights for variables
        #holder[j,i] <- mahalanobis(as.numeric(vec1),as.numeric(vec2),covX)
        holder[j,i] <- mahalanobis(as.numeric(vec1),as.numeric(vec2),invCovW,inverted = T)
        ## Jas' weights would be 
        #holder[j,i] <- mahalanobis(as.numeric(vec1),as.numeric(vec2),invCovW,inverted = T)
      }
    }
  }

    ## I calculate the distances differently depending on whether
    ##  there is a particular case to match.  If there is no case
    ##  to match, then I calculate the sum of distances between
    ##  all possible combinations of the units.  If there IS a case
    ##  to match, then I calculate just the difference between the
    ##  the match.case and every other case.
  if(is.null(match.case)==T) {
      ## Make holder2 as a matrix with the distances of 
      ## every possible combination of the units.
      ## "holder2" has case.N choose 2 columns because I want
      ##   to sum the distance of every possible pairwise combination
      ##   of the units within each ordering
    holder2 <- matrix(NA,nrow(combinations),(choose(case.N,2)))

      ## make a matrix of all possible pairings of the case.N units 
    pairings <- t(combn(case.N,2))

      ## name the columns of the holder
    colnames(holder2)<- paste(pairings[,1],"-",pairings[,2], sep="")

      ## fill holder2
    for(k in 1:(choose(case.N,2))){
      for(i in 1:nrow(combinations)){
        holder2[i,k] <- holder[combinations[i,pairings[k,2]],combinations[i,pairings[k,1]]]
      }
    }
  }  ## ends "if" statement

  if(is.null(match.case)==F) { 
    holder2 <- matrix(NA,nrow(combinations),(case.N-1))
      ## make a matrix of all possible pairings of the case.N units 
    pairings <- cbind(rep(1,length(seq(2,case.N,1))), seq(2,case.N,1))

      ## name the columns of the holder
    colnames(holder2)<- paste(pairings[,1],"-",pairings[,2], sep="")

      ## fill holder2
    for(k in 1:ncol(holder2)){
      for(i in 1:nrow(combinations)){
        holder2[i,k] <- holder[combinations[i,pairings[k,2]],1]  ## column is 1 because I made the match.case.id=1
      }
    }
  }

    ## This sums the distances to find the overall mahalanobis distance
  total.distance <- apply(holder2, MARGIN=1, FUN=sum)

  
    ## I need to attach the unit names to it, but I'm not sure how to 
    ##  do it without a forloop.
    ## Without a loop, it looks like :  distances <- cbind(total.distance, combinations)
    ##  but this only gives the internally coded id numbers 1-N.
  distances <- data.frame(total.distance)
  for(i in 1:case.N){
    distances <- cbind(distances, unit.names[combinations[,i]])
  }

    ## Sort by the overall distance
  distances.sort <- distances[order(distances[,1]),]
  colnames(distances.sort) <- c("distances", rep("unit id",case.N))
  
    ## sort holder2 by the overall distance
  holder2.sort <- data.frame(holder2[order(distances[,1]),])
  colnames(holder2.sort)<- paste(pairings[,1],"-",pairings[,2], sep="")



    ## THIS SECTION MAXIMIZES THE VARIANCE ON THE TREAT VAR
    ## Note that maximizing the variance IS NOT always what we want
    ## var(c(1,1,4,4)) > var(c(1,2,3,4))

    ## If max.variance=T, pull out the treatment variable 
    ##   and create a vector of the variances of the variable
  if(max.variance==TRUE){
    variances <- rep(NA,nrow(distances.sort))
    for(i in 1:nrow(distances.sort)){
      current.row <- distances.sort[i,]
        ## I think I can just pull out these rows because I havent' reshuffled the data
      units.matched <- unlist(apply(current.row, MARGIN=2, FUN=function(x){which(unit.names==x)}))
      values.to.check.variance <- data[units.matched, colnames(data)==treatment.var]
      variances[i] <- var(values.to.check.variance)
    }
  }

    ## trims the holders by the variances
  if(max.variance==TRUE){
      ## If there are enough matches that have the highest possible variance
      ##  I just drop all the matches that don't have the max variance
    if(length(variances[variances==max(variances,na.rm=T)]) > number.of.matches.to.return){
      distances.sort.trim <- distances.sort[variances==max(variances,na.rm=T),]
      treat.variance <- variances[variances==max(variances,na.rm=T)]
    }
      ## If there are too few matches that have the maximum variance,
      ##  I take the top "variance.tolerance" proportion of the observations
    if(length(variances[variances==max(variances,na.rm=T)]) < number.of.matches.to.return){
      variances >= quantile(variances,probs=variance.tolerance)
      distances.sort.trim <- distances.sort[variances >= quantile(variances,probs=(1-variance.tolerance)),]
      treat.variance <- variances[variances >= quantile(variances,probs=(1-variance.tolerance))]
    }
    distances.sort <- cbind(distances.sort.trim,treat.variance)
  }



    ## THIS SECTION MAXIMIZES THE SPREAD ON THE TREAT VAR
    ## Note that maximizing the variance IS NOT always what we want
    ## we want c(1,2,3,4)) instead of c(1,1,4,4)

    ## If max.spread=T, pull out the treatment variable 
    ##   and create a vector of the variances of the variable
  if(max.spread==TRUE){

      ## This is a function that calculates the overall spread
      ##   of a vector -- the sum of the distances between the ordered components
    spreadem <- function(vec, overall.min=min(vec), overall.max=min(vec)){
      vec <- sort(vec)
      spread.dist <- c()
      for(k in 1:(length(vec)-1)){
        spread.dist <- c(spread.dist, vec[k+1]-vec[k])
      }
      return(sum(spread.dist))
    }

      ## Then create a holder for the spreads
    spread <- rep(NA,nrow(distances.sort))
    for(i in 1:nrow(distances.sort)){
      current.row <- distances.sort[i,]
        ## I think I can just pull out these rows because I haven't reshuffled the data
      units.matched <- unlist(apply(current.row, MARGIN=2, FUN=function(x){which(unit.names==x)}))
      values.to.check.spread <- sort(data[units.matched, colnames(data)==treatment.var])
        ## calculates the difference between order elements of the vector using spreadem()
      spread[i] <- spreadem(values.to.check.spread)
    }
  }

    ## trims the holders by the spreads
  if(max.spread==TRUE){
      ## If there are enough matches that have the highest possible spread
      ##  I just drop all the matches that don't have the minimum difference
      ##  from the ideal spread
    if(length(spread[spread==max(spread,na.rm=T)]) > number.of.matches.to.return){
      distances.sort.trim <- distances.sort[spread==max(spread,na.rm=T),]
      treat.spread <- spread[spread==max(spread,na.rm=T)]
    }
      ## If there are too few matches that have the maximum variance,
      ##  I take the top "variance.tolerance" proportion of the observations
    if(length(spread[spread==max(spread,na.rm=T)]) < number.of.matches.to.return){
      distances.sort.trim <- distances.sort[spread >= quantile(spread,probs=(1-variance.tolerance)),]
      treat.spread <- spread[spread >= quantile(spread,probs=(1-variance.tolerance))]
    }
    distances.sort <- cbind(distances.sort.trim,treat.spread)
  }






    ## PUT THE FINAL TABLE TOGETHER
    ## Do different things if it is "most similar" or "most different"
    ## Most similar
  if(design.type=="most similar"){
      ##  the object that lists the closest matches
    matches <- distances.sort[1:number.of.matches.to.return,]
    rownames(matches) <- seq(1,number.of.matches.to.return,1)
      ## put together a list with the distances between each pair in the matched set
      ## Have to do something different with case.N=2 
    if(case.N>=3) {
      match.dist <- holder2.sort[1:number.of.matches.to.return,]
      } else {
      match.dist <- data.frame(as.matrix(holder2.sort[1:number.of.matches.to.return,]))
      colnames(match.dist) <- colnames(holder2.sort)
    }
    match.units <- matches[,-1]
    case.distances <- c()
    for(i in 1:nrow(match.dist)){
      case.distances[[i]] <- data.frame(match.dist[i,])
      colnames(case.distances[[i]])<- paste(as.matrix(match.units[i,pairings[,1]]),"-",as.matrix(match.units[i,pairings[,2]]), sep="")
    }
    return(list("cases"=matches, "case.distances"=case.distances))
  }

    ## most different
  if(design.type=="most different"){
      ## the object with the most different matches
    matches <- distances.sort[nrow(distances.sort):(nrow(distances.sort)-number.of.matches.to.return + 1),]
    rownames(matches) <- rev(seq((nrow(distances.sort)-number.of.matches.to.return + 1),nrow(distances.sort),1))
      ## put together a list with the distances between each pair in the matched set 
      ## Need to do something special for case.N=2
    if(case.N>=3) {
      match.dist <- holder2.sort[nrow(distances.sort):(nrow(distances.sort)-number.of.matches.to.return + 1),]
      } else {
      match.dist <- data.frame(as.matrix(holder2.sort[nrow(distances.sort):(nrow(distances.sort)-number.of.matches.to.return + 1),]))
      colnames(match.dist) <- colnames(holder2.sort)
    }
    match.units <- matches[,-1]
    case.distances <- c()
    for(i in 1:nrow(match.dist)){
      case.distances[[i]] <- data.frame(match.dist[i,])
      colnames(case.distances[[i]])<- paste(as.matrix(match.units[i,pairings[,1]]),"-",as.matrix(match.units[i,pairings[,2]]), sep="")
    }
    return(list("cases"=matches, "case.distances"=case.distances))
  }




}   ## END OF FUNCTION

