
case.match <- function(data, ## the data set
                        id.var, ## a variable in the data set, specifying unit ids. Must be the name of the id variable in quotes
                        case.N=2, ## the number of cases the function should return
                        distance="mahalanobis", ## the distance metric to use: "mahalanobis", "euclidean", or "standardized"
                        design.type="most similar", ## specify "most similar" or "most different"
                        match.case=NULL, ## the id of a specific case to be matched
                        greedy.match="pareto", ## what matches to report: "pareto", "all", "greedy"
                        number.of.matches.to.return=1, ## how many possible pairings should the function return?
                        treatment.var=NULL, ## the name of the treatment variable, if any
                        outcome.var=NULL, ## the name of the outcome variable, if any
                        leaveout.vars=NULL, ## specify variables by name to leave out of the matching
                        max.variance=FALSE, ## maximize variance of the treatment variable?
                        max.variance.outcome=FALSE, ## maximize variance of the outcome variable?
                        variance.tolerance=.1,
                        max.spread=FALSE, ## maximize "spread" on the treatment, which tries to evenly space cases
                        max.spread.outcome=FALSE, ## maximize "spread" on the outcome, which tries to evenly space cases
                        varweights=NULL ## Optional vector of variable weights
                        ){
  
  ## Make sure data is coming in as a data frame, not a tibble (in response to a bug-fix request)
  data <- as.data.frame(data)
  ## Error messages
  if (case.N<2) {
      stop("case.N must be an integer 2 or larger", call. = FALSE) }
  ## stops if length(match.case)>1 (the user is trying to find matches for two or more specific cases)
  if(length(match.case)>1){
    stop("match.case should be a string vector length = 1", call. = FALSE) }
  ## Warns you if you have too many treatment variables specified
  if(length(treatment.var)>1){
    stop("Only one treatment variable can be specified.")}
  ## warns you if the treatment variable is missing
  if(max(names(data)==treatment.var, na.omit=T)==F) {
      stop("Treatment variable is missing from the dataset.")}
  ## warns you if one or more of the leaveout vars isn't in the dataset
  if(inherits(leaveout.vars, "character")){
    for(i in 1:length(leaveout.vars)){
      if(max(names(data)==leaveout.vars[i], na.omit=T)==F) {
        stop("At least 1 leaveout.vars variable is missing from the dataset.") 
      }
    }
  }  
  ## warning about specifying both max.variance and max.spread
  if(max.variance==T & max.spread==T){
    stop("Cannot specify max.variance==TRUE & max.spread==TRUE.  Pick one or the other.")}
  ## Warns you if you have too many outcome variables specified
  if(length(outcome.var)>1){
    stop("Only one outcome variable can be specified.")}
  ## warns you if the outcome variable is missing
  if(max(names(data)==outcome.var, na.omit=T)==F) {
    stop("Outcome variable is missing from the dataset.")}
  ## warns about max.variance.outcome==T & max.spread.outcome==T
  if(max.variance.outcome==T & max.spread.outcome==T){
    stop("Cannot specify max.variance.outcome==TRUE & max.spread.outcome==TRUE.  Pick one or the other.")}
  ## warns about is.null(treatment.var)==T & max.variance==T
  if(is.null(treatment.var)==T & max.variance==T){
    stop("Cannot specify max.variance==TRUE without specifying treatment variable")}
  ## warns about is.null(treatment.var)==T & max.spread==T
  if(is.null(treatment.var)==T & max.spread==T){
    stop("Cannot specify max.spread==TRUE without specifying treatment variable")}
  ## warns about is.null(outcome.var)==T & max.variance.outcome==T
  if(is.null(outcome.var)==T & max.variance.outcome==T){
    stop("Cannot specify max.variance.outcome==TRUE without specifying outcome variable")}
  ## warns about is.null(outcome.var)==T & max.spread.outcome==T
  if(is.null(outcome.var)==T & max.spread.outcome==T){
    stop("Cannot specify max.spread.outcome==TRUE without specifying outcome variable")}
  if(design.type !="most similar" & design.type !="most different"){
    stop("Design type is not correctly specified")}
  if(greedy.match=="greedy" & design.type=="most different" & is.null(match.case)==T & case.N>2){
    stop("Most different design with greedy matching only allows case.N = 2. Try greedy.match = 'pareto'")}
  
  
  ## Preliminary data manipulation

  ## If leaveout variables are numeric...
  if(inherits(leaveout.vars,"numeric")) { 
    leaveout <- leaveout.vars
  }
  ## If leaveout vector is character, use charmatch() to find the column numbers of the user specified variable names
  if(inherits(leaveout.vars,"character")){
    leaveout <- charmatch(leaveout.vars,names(data))
  }
  ## Print which variables are being matched
  if(is.null(leaveout.vars)==F ){
    print(paste0("Matching variables: ", paste(names(data[,-c(leaveout),drop=F]),collapse=", ")))
  } else { 
    print(paste0("Matching variables: ", paste(names(data),collapse=", ")))
  }

  ## Take out NAs in the data.  Only keep observations that are complete (information for matching, treatment, etc.)
  ## Code differs slightly, depending on if leaveout.vars is specified
  if(is.null(leaveout.vars)==F){
    ## Print message if there are missing variables
    if(ifelse(sum(as.numeric(is.na(data[,-c(leaveout),drop=F])))>0, TRUE, FALSE)==TRUE){
      print ("Some observations are missing data on the specified variables.  These observations have been removed.")
    }
    missing<-as.vector(row.names(na.omit(data[,-c(leaveout),drop=F])))
    data<-data[which(rownames(data) %in% missing),,drop=F]
  } else {
    # Print message if there are missing variables
    if(ifelse(sum(as.numeric(is.na(data)))>0, TRUE, FALSE)==TRUE){
      print ("Some observations are missing data on the specified variables.  These observations have been removed.")
    }
    missing<-as.vector(row.names(na.omit(data)))
    data<-data[which(rownames(data) %in% missing),,drop=F]
  }
  
  ## pull out the unit names
  unit.names <- as.vector(data[,which(names(data)==id.var)])
  
  ## pull out a matrix with just the variables to use for matching
  ## if there are no variables to leave out and no treatment or outcome var
  if(is.null(leaveout.vars)==T & is.null(treatment.var)==T & is.null(outcome.var)==T){
    X <- data[,-which(names(data)==id.var),drop=F]
  }
  ## if there are vars to leave out but no treatment or outcome var
  if(is.null(leaveout.vars)==F & is.null(treatment.var)==T & is.null(outcome.var)==T){
    X <- data[,-c(which(names(data)==id.var), leaveout),drop=F]
  }
  ## if there are vars to leave out and one treatment var 
  if(is.null(leaveout.vars)==F & is.null(treatment.var)==F & is.null(outcome.var)==T){
    X <- data[,-c(which(names(data)==id.var),
                  leaveout,
                  which(names(data)==treatment.var)),drop=F]
    treat <- data[,which(names(data)==treatment.var)]
  } 
  ## if there are no vars to leave out but one treatment var
  if(is.null(leaveout.vars)==T & is.null(treatment.var)==F & is.null(outcome.var)==T){
    X <- data[,-c(which(names(data)==id.var), 
                  which(names(data)==treatment.var)),drop=F]
    treat <- data[,which(names(data)==treatment.var)]
  }
  ## if there are no vars to leave out but one outcome var
  if(is.null(leaveout.vars)==T & is.null(treatment.var)==T & is.null(outcome.var)==F){
    X <- data[,-c(which(names(data)==id.var), 
                  which(names(data)==outcome.var)),drop=F]
    outcome <- data[,which(names(data)==outcome.var)]
  }
  ## if there are no vars to leave out but one treatment var and outcome var
  if(is.null(leaveout.vars)==T & is.null(treatment.var)==F & is.null(outcome.var)==F){
    X <- data[,-c(which(names(data)==id.var), 
                  which(names(data)==treatment.var),
                  which(names(data)==outcome.var)),drop=F]
    treat <- data[,which(names(data)==treatment.var)]
    outcome <- data[,which(names(data)==outcome.var)]
  }
  ## if there are vars to leave out and one outcome var
  if(is.null(leaveout.vars)==F & is.null(treatment.var)==T & is.null(outcome.var)==F){
    X <- data[,-c(which(names(data)==id.var), 
                  leaveout,
                  which(names(data)==outcome.var)),drop=F]
    outcome <- data[,which(names(data)==outcome.var)]
  }
  ## if there are vars to leave out and one treatment var and one outcome var
  if(is.null(leaveout.vars)==F & is.null(treatment.var)==F & is.null(outcome.var)==F){
    X <- data[,-c(which(names(data)==id.var),
                  leaveout,
                  which(names(data)==treatment.var),
                  which(names(data)==outcome.var)),drop=F]
    treat <- data[,which(names(data)==treatment.var)]
    outcome <- data[,which(names(data)==outcome.var)]
  }
  
  ## Warnings about mismatch in variable weights, now that X is tidied up
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
  ## creates a new id variable to use throughout the algorithm
  id <- seq(1,N,1)
  ## create the covariance matrix of the whole data
  if(distance=="mahalanobis") {covX <- cov(X)}
  if(distance=="euclidean"){covX <- diag(ncol(X))}
  ## standardize the data by its standard dev. if specified.
  if(distance=="standardized") {
    stddev <- apply(X,MARGIN=2,FUN=sd)
    X <- as.data.frame(t(apply(X,MARGIN=1,FUN=function(x){x/stddev})))
    covX <- diag(ncol(X))
  }
  ## Warning messages about missing covX
  if(exists("covX")==F) {
    stop("Matching distance not specified properly. Check spelling")}
  ## Check if matrix can be inverted
  f <- function(m){ "matrix" %in% class(try(solve(m),silent=T)) }
  if(f(covX)==F) {
    stop("Matching variables are too collinear to calculate the mahalanobis distance.  
         Either retry with euclidean distance or a different combination of matching variables.")}
  
  ## Printe messages about the number of matches
  if(is.null(match.case)==T){
    print(paste("There are",N,"choose",case.N,"=", choose(N, case.N), "possible case combinations"))
  }
  if(is.null(match.case)==F){
    print(paste("There are",(N-1),"choose",(case.N-1),"=", choose((N-1), (case.N-1)), "possible case combinations"))
  }

  ## Prep if there is a specific case to match:
  if(is.null(match.case)==F){ 
    ## Issue a warning if the name of the unit to match doesn't exist
    if(length(grep(match.case,unit.names, fixed=T))==0) {
      stop("match.case does not match any of the unit IDs", call. = FALSE) }
    ## This pulls out the id of the unit to match
    match.case.id <- id[unit.names==match.case]
    ## Issue a warning if the there are too many units to match 
    if(length(match.case.id)>1) {
      stop("match.case matches more than one unit ID", call. = FALSE) }
    ## make a new id var that has the unit to match having id=1
    id2 <- rep(NA,length(id))
    id2[id==match.case.id] <- 1
    id2[id!=match.case.id] <- seq(2,N,1)
    id <- id2
    match.case.id <- 1
    ## I reorder the dataset so that the match.case is first
    X <- X[order(id),,drop=F]
    ## reorder the unit names and the id and treat variables as well
    unit.names <- unit.names[order(id)]
    if(is.null(treatment.var)==F){treat <- treat[order(id)]}
    if(is.null(outcome.var)==F){outcome <- outcome[order(id)]}
    id <- sort(id)
    
    ## If there is no treat or outcome var, I do all the possible combinations:
    if(is.null(treatment.var)==T & is.null(outcome.var)==T){
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
      print ("Matching only on units that have a different treatment value")
      treat.value <- treat[id==match.case.id]
      combin <- t(combn(id[-unique(c(match.case.id,which(treat==treat.value)))],(case.N-1)))  ## note, I already took out 
      combinations <- cbind(rep(match.case.id,nrow(combin)), combin)
    }
    
    ## If there is a outcome variance specified, I do only the combinations of the other outcomes:
    if(is.null(outcome.var)==F){
       ## a different set of combinations -- only those that include the match.case
       ## using only the units with different outcomes
       print ("Matching only on units that have a different outcome value")
       outcome.value <- outcome[id==match.case.id]
       combin <- t(combn(id[-unique(c(match.case.id,which(outcome==outcome.value)))],(case.N-1)))  ## note, I already took out 
       combinations <- cbind(rep(match.case.id,nrow(combin)), combin)
     }
  } ## end of if(is.null(match.case)==F)

  ## If there isn't a case to match, then do this:
  if(is.null(match.case)==T) { 
    combinations <- t(combn(id,case.N))
  }

  ## Calculate distances
  ## make a holder matrix to hold the distances between each unit
  holder <- matrix(NA,nrow(X),nrow(X))
  ## weights for variables
  if(is.null(varweights)){varweights <- rep(1,ncol(X))}
  W <- diag(varweights)
  invCovW <- W %*% solve(covX) %*% W
  rownames(invCovW) <- colnames(invCovW) <- names(X)
  ## in the genmatch paper in Review of Econ. and Stat., Jas has a different version of weights that gives a different answer!
  #jasW <- t(chol(covX)) %*% W %*% (chol(covX))

  ## Fast distance calculation code
  ## code from C:\Users\Richard Nielsen\Desktop\Papers\pscore paradox\paradoxMdisc\paradoxFuns_10aug2010.R
  ## myMH function from Ben Hansen at http://www.stat.lsa.umich.edu/~bbh/optmatch/doc/mahalanobisMatching.pdf (broken link)
  ## Now at http://www.mit.edu/~rnielsen/mahalDistanceCode.txt
  myMH <- function(Tnms, Cnms, inv.cov, data) {
    stopifnot(!is.null(dimnames(inv.cov)[[1]]),# dim(inv.cov)[1] > 1,
    all.equal(dimnames(inv.cov)[[1]], dimnames(inv.cov)[[2]]),
    all(dimnames(inv.cov)[[1]] %in% names(data)))
    covars <- dimnames(inv.cov)[[1]]
    xdiffs <- as.matrix(data[Tnms, covars])
    xdiffs <- xdiffs - as.matrix(data[Cnms, covars])
    rowSums((xdiffs %*% inv.cov) * xdiffs)
  }
  ## only calculate distances for the case to be matched, if one is specified
  if(is.null(match.case)){
    holder <- outer(rownames(X), rownames(X), FUN = myMH, inv.cov = invCovW, data = X)
  } else {
    holder <- outer(match.case.id, rownames(X), FUN = myMH, inv.cov = invCovW, data = X)
    holder <- t(holder)
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

    ## fill holder2 (faster thanks to solution from jblood94 here: https://stackoverflow.com/questions/76374499/can-a-for-loop-be-avoided-when-comparing-matrices)
    holder2 <- matrix(
      holder[matrix(combinations[,pairings[,2:1]], ncol = 2)],
      nrow(combinations), nrow(pairings),
      dimnames = list(NULL, paste(pairings[,1], pairings[,2], sep = "-"))
    )
  }  ## ends if(is.null(match.case)==T)

  if(is.null(match.case)==F) { 
    holder2 <- matrix(NA,nrow(combinations),(case.N-1))
    ## make a matrix of all possible pairings of the case.N units 
    pairings <- cbind(rep(1,length(seq(2,case.N,1))), seq(2,case.N,1))

    ## name the columns of the holder
    colnames(holder2)<- paste(pairings[,1],"-",pairings[,2], sep="")

    ## fill holder2
    holder2 <- matrix(
      holder[matrix(combinations[,pairings[,2:1]], ncol = 2)],
      nrow(combinations), nrow(pairings),
      dimnames = list(NULL, paste(pairings[,1], pairings[,2], sep = "-"))
    )
  }  ## ends is.null(match.case)==F

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
  
  ## THIS SECTION MAXIMIZES THE VARIANCE ON THE TREAT VAR
  ## Note that maximizing the variance IS NOT always what we want
  ## var(c(1,1,4,4)) > var(c(1,2,3,4))

  ## If max.variance=T, pull out the treatment variable 
  ##   and create a vector of the variances of the variable
  if(max.variance==TRUE){
    tvals <- apply(combinations,MARGIN=2,function(x){data[[treatment.var]][x]})
    variances <- apply(tvals,MARGIN=1,var)
    
    ## Sort by the overall distance
    distances.sort <- distances[order(distances[,1]),]
    variances <- variances[order(distances[,1])]
    colnames(distances.sort) <- c("distances", paste(rep("unit id",case.N),1:case.N))
    ## sort holder2 by the overall distance
    holder2.sort <- data.frame(holder2[order(distances[,1]),])
    colnames(holder2.sort)<- paste(pairings[,1],"-",pairings[,2], sep="")

    ## trims the holders by the variances
    
    ## If there are enough matches that have the highest possible variance
    ##  I just drop all the matches that don't have the max variance
    if(length(variances[which(variances==max(variances,na.rm=T))]) > number.of.matches.to.return){
      distances.sort.trim <- distances.sort[which(variances==max(variances,na.rm=T)),]
      treat.variance <- variances[which(variances==max(variances,na.rm=T))]
    }
    ## If there are too few matches that have the maximum variance,
    ##  I take the top "variance.tolerance" proportion of the observations
    if(length(variances[which(variances==max(variances,na.rm=T))]) < number.of.matches.to.return){
      distances.sort.trim <- distances.sort[which(variances >= quantile(variances,probs=(1-variance.tolerance), na.rm=TRUE)),]
      treat.variance <- variances[which(variances >= quantile(variances,probs=1-variance.tolerance, na.rm=TRUE))]
    }
    distances.sort <- cbind(distances.sort.trim,treat.variance)
  } ## end if(max.variance==TRUE)
  

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
    } ## end spreadem()
    
    ## Then create a holder for the spreads
    tvals <- apply(combinations,MARGIN=2,function(x){data[[treatment.var]][x]})
    spread <- apply(tvals,MARGIN=1,spreadem)
  
    ## Sort by the overall distance
    distances.sort <- distances[order(distances[,1]),]
    spread <- spread[order(distances[,1])]
    colnames(distances.sort) <- c("distances", paste(rep("unit id",case.N),1:case.N))
    ## sort holder2 by the overall distance
    holder2.sort <- data.frame(holder2[order(distances[,1]),])
    colnames(holder2.sort)<- paste(pairings[,1],"-",pairings[,2], sep="")

    ## trims the holders by the spreads
    
    ## If there are enough matches that have the highest possible spread
    ##  I just drop all the matches that don't have the minimum difference
    ##  from the ideal spread
    if(length(spread[which(spread==max(spread,na.rm=T))]) > number.of.matches.to.return){
      distances.sort.trim <- distances.sort[which(spread==max(spread,na.rm=T)),]
      treat.spread <- spread[which(spread==max(spread,na.rm=T))]
    }
    ## If there are too few matches that have the maximum variance,
    ##  I take the top "variance.tolerance" proportion of the observations
    if(length(spread[which(spread==max(spread,na.rm=T))]) < number.of.matches.to.return){
      distances.sort.trim <- distances.sort[which(spread >= quantile(spread,probs=(1-variance.tolerance), na.rm=TRUE)),]
      treat.spread <- spread[which(spread >= quantile(spread,probs=(1-variance.tolerance), na.rm=TRUE))]
    }
    distances.sort <- cbind(distances.sort.trim,treat.spread)
  }

  ## Maximize variance for outcome variable
  if(max.variance.outcome==TRUE){
    ovals <- apply(combinations,MARGIN=2,function(x){data[[outcome.var]][x]})
    variances <- apply(ovals,MARGIN=1,var)
    
    ## Sort by the overall distance
    distances.sort <- distances[order(distances[,1]),]
    variances <- variances[order(distances[,1])]
    colnames(distances.sort) <- c("distances", paste(rep("unit id",case.N),1:case.N))
    ## sort holder2 by the overall distance
    holder2.sort <- data.frame(holder2[order(distances[,1]),])
    colnames(holder2.sort)<- paste(pairings[,1],"-",pairings[,2], sep="")
  
    ## trims the holders by the variances
    ## If there are enough matches that have the highest possible variance
    ##  I just drop all the matches that don't have the max variance
    if(length(variances[which(variances==max(variances,na.rm=T))]) > number.of.matches.to.return){
      distances.sort.trim <- distances.sort[which(variances==max(variances,na.rm=T)),]
      outcome.variance <- variances[which(variances==max(variances,na.rm=T))]
    }
    
    ## If there are too few matches that have the maximum variance,
    ##  I take the top "variance.tolerance" proportion of the observations
    if(length(variances[which(variances==max(variances,na.rm=T))]) < number.of.matches.to.return){
      distances.sort.trim <- distances.sort[which(variances >= quantile(variances,probs=(1-variance.tolerance), na.rm=TRUE)),]
      outcome.variance <- variances[which(variances >= quantile(variances,probs=1-variance.tolerance, na.rm=TRUE))]
    }
    distances.sort <- cbind(distances.sort.trim,outcome.variance)
  }
  
  ## THIS SECTION MAXIMIZES THE SPREAD ON THE OUTCOME VAR
  ## Note that maximizing the variance IS NOT always what we want
  ## we want c(1,2,3,4)) instead of c(1,1,4,4)
  
  ## If max.spread.outcome=T, pull out the outcome variable 
  ##   and create a vector of the variances of the variable
  if(max.spread.outcome==TRUE){
    
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
    ovals <- apply(combinations,MARGIN=2,function(x){data[[outcome.var]][x]})
    spread <- apply(ovals,MARGIN=1,spreadem)
    
    ## Sort by the overall distance
    distances.sort <- distances[order(distances[,1]),]
    spread <- spread[order(distances[,1])]
    colnames(distances.sort) <- c("distances", paste(rep("unit id",case.N),1:case.N))
    ## sort holder2 by the overall distance
    holder2.sort <- data.frame(holder2[order(distances[,1]),])
    colnames(holder2.sort)<- paste(pairings[,1],"-",pairings[,2], sep="")
    
    ## trims the holders by the spreads
    ## If there are enough matches that have the highest possible spread
    ##  I just drop all the matches that don't have the minimum difference
    ##  from the ideal spread
    if(length(spread[which(spread==max(spread,na.rm=T))]) > number.of.matches.to.return){
      distances.sort.trim <- distances.sort[which(spread==max(spread,na.rm=T)),]
      outcome.spread <- spread[which(spread==max(spread,na.rm=T))]
    }
    ## If there are too few matches that have the maximum variance,
    ##  I take the top "variance.tolerance" proportion of the observations
    if(length(spread[which(spread==max(spread,na.rm=T))]) < number.of.matches.to.return){
      distances.sort.trim <- distances.sort[which(spread >= quantile(spread,probs=(1-variance.tolerance), na.rm=TRUE)),]
      outcome.spread <- spread[which(spread >= quantile(spread,probs=(1-variance.tolerance), na.rm=TRUE))]
    }
    distances.sort <- cbind(distances.sort.trim,outcome.spread)
  }

  # Added 09/2016
  # Options to get unique matches only  
  # Start with pareto matches
  # Most similar
  
  if(!exists("distances.sort")){
    ## Sort by the overall distance
    distances.sort <- distances[order(distances[,1]),]
    colnames(distances.sort) <- c("distances", paste(rep("unit id",case.N),1:case.N))
    ## sort holder2 by the overall distance
    holder2.sort <- data.frame(holder2[order(distances[,1]),])
    colnames(holder2.sort)<- paste(pairings[,1],"-",pairings[,2], sep="")
  }
  
  if(greedy.match=="pareto" & design.type=="most similar" & is.null(match.case)==T & case.N<3){
    
    # take out matches where both units have a better match
    sub_matches<-c()
    for (i in 1:length(data[,which(names(data)==id.var)])){
      sub_matches[i]=min(distances.sort[distances.sort[,c(2)]==data[,which(names(data)==id.var)][i] | distances.sort[,c(3)]==data[,which(names(data)==id.var)][i],1]) 
    }
    
    is.na(sub_matches) <- do.call(cbind,lapply(sub_matches, is.infinite))
    sub_matches<-na.omit(sub_matches)
    sub_matches<-sort(unique(sub_matches))
    test_final=distances.sort[distances.sort[,1] %in% sub_matches,]
    
    if(nrow(test_final) > number.of.matches.to.return){
      distances.sort=test_final  
    }
    ## If there are too few matches 
    if(nrow(test_final) < number.of.matches.to.return){
      print("For more matches, change greedy.match to 'all'")
      distances.sort <- na.omit(test_final)
    }
  }

  # Most different
  
  if(greedy.match=="pareto" & design.type=="most different" & is.null(match.case)==T){
    distances.sort2<-distances.sort[order(-distances.sort[,1]),]
    sub_matches<-c()
    for (i in 1:length(data[,which(names(data)==id.var)])){
      # Original line that I think has a mistake because it assumes only case.N=2
      #sub_matches[i]=max(distances.sort2[distances.sort2[,c(2)]==data[,which(names(data)==id.var)][i] | distances.sort2[,c(3)]==data[,which(names(data)==id.var)][i],1]) 
      # A fixed version that mirrors the original code for case.N=3
      #sub_matches[i]=max(distances.sort2[distances.sort2[,c(2)]==data[,which(names(data)==id.var)][i] | distances.sort2[,c(3)]==data[,which(names(data)==id.var)][i] | distances.sort2[,c(4)]==data[,which(names(data)==id.var)][i],1]) 
      ## A fixed version that's more flexible
      idx <- data[,which(names(data)==id.var)][i]
      sub_matches[i]=max(distances.sort2[rowSums(apply(distances.sort2[,2:(2+case.N-1)],MARGIN=2,FUN=function(x){x==idx}))>0, 1])
    }
    tmp <- sub_matches
    sub_matches==tmp
    is.na(sub_matches) <- do.call(cbind,lapply(sub_matches, is.infinite))
    sub_matches<-na.omit(sub_matches)
    sub_matches<-sort(unique(sub_matches))
    test_final=distances.sort[distances.sort[,1] %in% sub_matches,]
    
    if(nrow(test_final) > number.of.matches.to.return){
      distances.sort=test_final  
    }
    ## If there are too few matches 
    if(nrow(test_final) < number.of.matches.to.return){
      print("For more matches, change greedy.match to 'all'")
      distances.sort <- test_final
    }
  }
  
  # Greedy matches
  
  if(greedy.match=="greedy" & design.type=="most similar" &is.null(match.case)==T){
    # Create Matrix of distances
    distance.temp<-distances.sort[,c(2,3)]
    m <- matrix(unlist(distance.temp), ncol=ncol(distance.temp))
    
    for (i in 1:length(data[,which(names(data)==id.var)])){
    
      caserows=which(m==data[,which(names(data)==id.var)][i], arr.ind=TRUE)
      caserows2<-caserows[order(caserows[,1])]
      row.keep=caserows2[-1]
      row.keep2<-na.omit(row.keep)
      if (length(row.keep2)>0){
        m=m[-row.keep2,]
      }
    }
    test_final2=distances.sort[distances.sort[,2] %in% m[,1],]
    
    test_final3<-test_final2[match(unique(test_final2[,2]), test_final2[,2]),]
    
    if(nrow(test_final3) > number.of.matches.to.return){
      distances.sort=test_final3  
    }
    ## If there are too few matches 
    if(nrow(test_final3) < number.of.matches.to.return){
      print("For more matches, change greedy.match to 'pareto' or 'all'")
      distances.sort <- na.omit(test_final3)
    }
  }
  
  # Most different
  
  if(greedy.match=="greedy" & design.type=="most different" & is.null(match.case)==T){
    distance.temp<-distances.sort[,c(2,3)][order(-distances.sort[,1]),]
    distance.matrix <- matrix(unlist(distance.temp), ncol=ncol(distance.temp))
    
    for (i in 1:length(data[,which(names(data)==id.var)])){
      print(i)
      caserows=which(distance.matrix==data[,which(names(data)==id.var)][i], arr.ind=TRUE)
      caserows2<-caserows[order(caserows[,1])]
      row.keep=caserows2[-1]
      row.keep2<-na.omit(row.keep)
      if (length(row.keep2)>0){
        distance.matrix=distance.matrix[-row.keep2,]
      }
    }
    test_final2=distances.sort[distances.sort[,2] %in% m[,1],]    
    test_final2<-test_final2[order(-test_final2[,1]),]
    test_final3<-test_final2[match(unique(test_final2[,2]), test_final2[,2]),]
    test_final3<-test_final3[order(test_final3[,1]),]
  
      if(nrow(test_final3) > number.of.matches.to.return){
      distances.sort=test_final3  
    }
    ## If there are too few matches 
    if(nrow(test_final3) < number.of.matches.to.return){
      print("For more matches, change greedy.match to 'pareto' or 'all'")
      distances.sort <- test_final3
    }
  }
  
  if(greedy.match=="all"){
    distances.sort=distances.sort
  }
  
  
  
  ## Gather the results
  ## Do different things if it is "most similar" or "most different"
  ## Most similar
  if(design.type=="most similar"){
    ##  the object that lists the closest matches
    matches <- distances.sort[1:number.of.matches.to.return,]
    rownames(matches) <- seq(1,number.of.matches.to.return,1)
    ## put together a list with the distances between each pair in the matched set
    ## Have to do something different with case.N=2 
    if(case.N>=3) {
      holder2.sort<-holder2.sort[which(rownames(holder2.sort) %in% rownames(distances.sort)),]
      match.dist <- holder2.sort[1:number.of.matches.to.return,]
      } else {
      match.dist <- data.frame(as.matrix(distances.sort[1:number.of.matches.to.return, 1]))
      colnames(match.dist) <- colnames(holder2.sort)
    }
    match.units <- matches[,-1]
    case.distances <- c()
    for(i in 1:nrow(match.dist)){
      case.distances[[i]] <- data.frame(match.dist[i,])
      colnames(case.distances[[i]])<- paste(as.matrix(match.units[i,pairings[,1]]),"-",as.matrix(match.units[i,pairings[,2]]), sep="")
    }
    print(head(matches[1:number.of.matches.to.return, 1: as.numeric(case.N+1)]))
    return(list("cases"=matches, "case.distances"=case.distances))  
  }

  ## most different
  if(design.type=="most different"){
    ## the object with the most different matches
    if(number.of.matches.to.return>nrow(distances.sort)){
      print("To produce more matches, change the matching specification")
      number.of.matches.to.return=nrow(distances.sort)
      }
    
    matches <- distances.sort[nrow(distances.sort):(nrow(distances.sort)-number.of.matches.to.return + 1),]
    rownames(matches) <- rev(seq((nrow(distances.sort)-number.of.matches.to.return + 1),nrow(distances.sort),1))
    ## put together a list with the distances between each pair in the matched set 
    ## Need to do something special for case.N=2
    if(case.N>=3) {
      holder2.sort<-holder2.sort[which(rownames(holder2.sort) %in% rownames(distances.sort)),]
      match.dist <- holder2.sort[nrow(distances.sort):(nrow(distances.sort)-number.of.matches.to.return + 1),]
      } else {
      match.dist <- data.frame(as.matrix(distances.sort[nrow(distances.sort):(nrow(distances.sort)-number.of.matches.to.return + 1), 1]))
      colnames(match.dist) <- colnames(holder2.sort)
    }
    match.units <- matches[,-1]
    case.distances <- c()
    for(i in 1:nrow(match.dist)){
      case.distances[[i]] <- data.frame(match.dist[i,])
      colnames(case.distances[[i]])<- paste(as.matrix(match.units[i,pairings[,1]]),"-",as.matrix(match.units[i,pairings[,2]]), sep="")
    }
    print(head(matches[1:number.of.matches.to.return, 1: as.numeric(case.N+1)]))
    return(list("cases"=matches, "case.distances"=case.distances))
  }
} ## end case.match()

