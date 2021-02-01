# SLIP MODEL 
# CURRENT VERSION : 
# Addition of fixation parameter; parameter to filter frameshift indel sequences 
# * Capitalized variables are global, initialized by setup() function *
source("~/vindels/2_within-host/slip-model-utils.r")


prior <- function(param){
  # some of these are weakly informative and can be narrowed
  prior.pe <- dunif(param[1], min=1e-7, max=1e-3, log=T)  
  prior.ps <- dunif(param[2], min=0.55, max=0.99, log=T) 
  prior.rate <- dlnorm(param[3],meanlog=log(0.00001),sdlog=1.1,log=T) #dunif(param[3], min=1e-7, max=1e-3, log=T)
  prior.fix <- dbeta(param[4], shape1=3, shape2=27, log=T) 
  
  return(prior.pe + prior.ps + prior.rate + prior.fix)
}

posterior <- function(param, slip, llh){
  if (any(param > 1)){
    return(log(0))
  }else{
    return(c(prior(param) , likelihood(param, slip, llh)))
  }
}

# Function responsible for proposing changes to the slip.list
changeSlip <- function(slip.vector, sd=0.5){
  index <- getSlipLocations(slip.vector)    # convert the slip list into indices
  
  tab <- table(index[[1]])
  
  # choose a slip event to change
  toEdit <- sample(length(index[[1]]),1)
  
  stack <- tab[as.character(index[[1]][toEdit])]
  if (stack > 1){   # have a set probability of moving the entire stack instead of individual events 
    
    all.idx <- which(index[[1]] == as.numeric(names(stack[1])))
    rnum <- runif(1)
    if (rnum < 0.9){
      # rewrite the toEdit variable with all indices in the stack
      toEdit <- all.idx 
    }else if(rnum > 0.9 && rnum < 0.95){
      if (runif(1) < 0.5){
        toEdit <- sample(all.idx, floor(length(all.idx) * 0.5))
      }else{
        toEdit <- sample(all.idx, ceiling(length(all.idx) * 0.5))
      }
    } 
  }
  
  # this is to ensure that the proposed change is never outside the slip region
  proposal <- index[[1]][toEdit] + delta(sd)
  while(any(proposal <= 0) || any(proposal > length(slip.vector))){
    proposal <- index[[1]][toEdit] + delta(sd)
  }
  # save the change to the slip list
  index[[1]][toEdit] <- proposal
  
  # save the whole list
  getSlipVector(index[[1]],index[[2]])
}

proposalFunction <- function(param, slip_current, llh_current){
  num <- runif(1)
  s2p <- 0.96     # "slip-to-parameter": dictates how often slips are changed compared to parameters
  
  if (num > s2p){  # CHANGE PARAMETERS 
    num <- runif(1)
    
    # weighting variables of each parameter (weight 4 is the remaining complement; 1-w3-w2-w1)
    w1 <- 0.2
    w2 <- 0.3
    w3 <- 0.3
    
    if (num < w1){         # 0.2 density
      param[1] <- rlnorm(1,meanlog=log(param[1]),sdlog=0.1)
      llh_proposed <- llh_current  # stays the same
      
    }else if(num > w1 && num < w2+w1){    # 0.3 density
      param[2] <- rlnorm(1,meanlog=log(param[2]),sdlog=0.02)
      llh_proposed <- llh_current   # stays the same
      
    }else if(num > w1+w2 && num < w1+w2+w3){     # 0.3 density
      param[3] <- rlnorm(1,meanlog=log(param[3]),sdlog=0.08)
      llh_proposed <- seqllh(param[3], slip_current)  # recalcuate using the new rate
      
    }else{        # 0.2 density
      param[4] <- rlnorm(1,meanlog=log(param[4]),sdlog=0.008)
      llh_proposed <- llh_current
    }
    slip_proposed <- slip_current    # always stays the same
  
  }else{   # CHANGE SLIP LIST
    # randomly choose a sequence to edit
    rand <- sample(length(SLIP.IDX),1)
    changed <- SLIP.IDX[rand]          # this is the location at which the slip.list was changed
    # convert it to indices
    slip_proposed <- slip_current
    slip_proposed[[changed]] <- changeSlip(slip_current[[changed]]) # this is the new slip.list with a single change
    
    # recalculate SINGLE pairwise llh at position = changed
    llh_proposed <- llh_current
    new.tip <- getTip(INDELS$tip[changed], slip_proposed[[changed]])
    llh_proposed[changed] <- pairllh(anc.seqs[changed], new.tip, param[3], BRANCHES[changed])
  }
  return(list(param=param, slip=slip_proposed, llh=llh_proposed))
}

runMCMC <- function(startvalue, iterations, runno, notes){
  # timing
  start.time <- proc.time()
  
  # initialize the chain
  chain <- array(dim = c(iterations+1,4))
  chain[1,] <- startvalue
  
  # start the slip list and the llh list
  slip_current <- SLIP.LIST
  llh_current <- seqllh(startvalue[3], SLIP.LIST)
  
  # keep a logfile up to date
  logfile <- file(paste0("~/PycharmProjects/hiv-withinhost/15_modeling/slip-", 
                         as.character(runno), "-",substr(gsub("[\\ :-]","",Sys.time()), 7, 12),
                         ".csv"), "w")
  notes <- gsub("^", "#",notes)
  notes <- gsub("\n", "\n#", notes)
  write(notes, file=logfile)
  write("p(Enter), p(Stay), Rate, Fix, Likelihood, Posterior, Slip-changed, Accept, Time", file=logfile, append=T)
  
  for (i in 1:iterations){
    # calculate posterior of current position
    p.current <- posterior(chain[i,], slip_current, llh_current)
    
    # generate new proposal 
    proposal <- proposalFunction(chain[i,], slip_current, llh_current)
    
    # calculate posterior of the new proposal (parameters, slip_proposed, llh_proposed)
    p.next <- posterior(proposal[[1]], proposal[[2]], proposal[[3]])
    #print(paste("Current:", p.current, "Next:", p.next, sep=" "))
    
    # --- For checking if the slip list has changed ---
    # s.change <- any(unname(unlist(slip_current))!=unname(unlist(proposal[[2]])))
    # print(paste0("Sliplist change proposed: ", s.change))
    
    # Metropolis-Hastings: Calculate the proportion of the proposal to the current position
    prop <- exp(sum(p.next) - sum(p.current))
    
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
      chain[i+1,] <- proposal[[1]]
      slip_current <- proposal[[2]]
      llh_current <- proposal[[3]]
      llh <- p.next[2]
      #print("Accept")
      accept <- T
    
    # if the proportion is less than the random uniform sample, REJECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
      #print("Reject")
      accept <- F
      llh <- p.current[2]
    }
    
    # Save the slip list every  iterations
    if (i %% 100 == 0){
      print(paste(c("STATE",i,":", chain[i,], sum(p.current)), collapse=" "))
      write(paste(c(chain[i,], llh, sum(p.current), as.numeric(accept), (proc.time() - start.time)[[3]]), collapse=",") , file=logfile, append=T)
    }
    
    # Save the slip list every 50000 iterations
    if (i %% 50000 == 0){
      slip_current <<- slip_current
      llh_current <<- llh_current
      slip <- unlist(lapply(slip_current, function(x){
        paste(x, collapse="")
      }))
      writeLines(slip, con=paste0("~/PycharmProjects/hiv-withinhost/15_modeling/list-",as.character(runno),"-",i,'.csv'))
      writeLines(as.character(llh_current), con=paste0("~/PycharmProjects/hiv-withinhost/15_modeling/llh-",as.character(runno),"-",i,'.csv'))
    }
  }
  close(logfile)
  return(list(chain=chain, slip=slip_current))
  
}
