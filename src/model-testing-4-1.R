setwd("~/indel-modeling/src")
source('slip-model-utils.r')
source("slippage-model-4-1.r")

nt <- c("A", "C", "G", "T")
lens <- c(75, 126, 105,  90 , 33)   # median nucleotide lengths of the five variable loops
f <- c(0.4158261, 0.1649874, 0.1761463, 0.2423612)
names(f) <- nt

# Generates a nucleotide sequence based on defined input frequencies 
genSeq <- function(len, freq=f){

  seq <- sapply(1:len, function(x){
    num <- runif(1)
    if (num < freq[1]){
      "A"
    }else if (num > freq[1] && num < (freq[2]+freq[1])){
      "C"
    }else if (num > (freq[2]+freq[1]) && num < (freq[3]+freq[2]+freq[1])){
      "G"
    }else{
      "T"
    }
  })
  paste(seq, collapse="")
}

# Simulates a pair of sequences 
  # Draws nucleotides based on constant underlying values
  # Creates nucleotides based on the rate parameter
  # Creates indel events based on the p.enter and p.stay parameters
  # Filters out frameshift-length indels based on the "fix" parameter
simPair <- function(param){
  p.enter <- param[1]
  p.stay <- param[2]
  rate <- param[3]
  fix <- param[4]
  
  vlen <- lens[sample(1:5, 1)]
  anc <- genSeq(vlen)
  
  # pick a branch rate value from a distribution 
  branch <- rlnorm(1, meanlog=4, sdlog=1.05)
  
  # ----SUBSTITUTIONS--------
  # add in substitutions based on the substitution rate value (param 3) 
  
  anc.chars <- strsplit(anc, "")[[1]]
  tmat <- getMat(rate, branch)
  
  # generate the tip sequence based on the transition rate matrix 
  tip <- paste(sapply(1:nchar(anc), function(n){
    rand <- runif(1)
    probs <- tmat[anc.chars[n],]
    
    if (rand < probs[1]){
      "A"
    }else if (rand > probs[1] && rand < (probs[1]+probs[2])){
      "C"
    }else if (rand > (probs[1]+probs[2]) && rand < (probs[1]+probs[2]+probs[3])){
      "G"
    }else{
      "T"
    }
  }), collapse="")
  
  # ------ INDELS -------
  # determine the number of insertions that occur 
  count <- sum(rbinom(1,vlen,p.enter))
  noFilter <- 0
  if (count > 0){
    
    # generate a vector of insertion lengths 
    lens <- c()
    for (n in 1:count){
      len <- 0
      p.exit <- 0
      while(p.exit < p.stay){
        len <- len + 1
        p.exit <- runif(1)
      }
      lens[n] <- len
    }
    noFilter <- length(lens)
    # apply a boolean filter to remove 91% of non-3 insertions 
    toRemove <- sapply(1:length(lens), function(x){
      if (lens[x] %% 3 != 0){
        runif(1) < fix
      }else{
        T
      }
    })
    lens <- lens[toRemove]
    
    # add any indels that make it through the filtering 
    if (length(lens) > 0){
      for (n in lens){
        # add the length and choose a random location for the slip event 
        idx <- sample(1:vlen,1)
        
        # generate the sequence to insert / delete 
        indel <- genSeq(n)
        
        # generate the tip and ancestor sequence by adding / removing sequence 
        tip <- insert(tip, indel, idx)
        anc <- insert(anc, rep("-",n), idx)
      }
    }
  }
  #print("finished indels")
  return(list(tip=tip,anc=anc,branch=branch, count=noFilter))
  # to estimate the ACTUAL rate of indels, you need to make failures occur after p.enter has been selected
  # algorithm:
  # probability of enter is chosen
  # draw a number from a poisson process
  # if the number is %% 3 == 0: 
  # keep it 100 percent
  # else if the number if %%3 != 0:
  # there's a low probability that it will be kept (penalty)
  
}

# Wrapper function for "simPair" above
# Runs simPair for set number of iterations and returns a dataframe with 
# all relevant information
simSeqs <- function(iter, param){
  all.seqs <- sapply(1:iter, function(n){
    if (n %% 1000 == 0 ){
      print(n)
    }
    pair <- simPair(param)
    # VALUE 1 = Tip, VALUE 2 = Ancestor, VALUE 3 = Branch length
    return(c(pair[[1]], pair[[2]], pair[[3]], pair[[4]]))
  })
  
  
  insertions <- as.data.frame(t(all.seqs), stringsAsFactors = F)
  colnames(insertions) <- c("tip", "anc", "branch", 'no.filter')
  
  # Generate the length and position columns
  data <- t(unname(sapply(insertions$anc, function(x){
    # returns c(length, position) of insertion events
    gaps <- gregexpr("-",x)[[1]]
    if (length(gaps) == 1 && gaps == -1){
      return (c(NA, NA))
    }else{
      return (c(length(gaps), max(gaps)))
    }
  })))
  
  insertions$len <- data[,1]
  insertions$pos <- data[,2]
  return(insertions)
}


#  ---- Start MCMC ----
param <- c(0.00052,0.88, 0.00001, 0.09)      # Load "true" parameter values into a vector
insertions <- simSeqs(25000, param)          # Simulate based on these values
setup(insertions$tip, insertions$anc, insertions$len, insertions$pos, insertions$branch, T)
startvalue <- c(0.0008, 0.65, 0.000001, 0.25)
notes <- "Test #22
Same as 21, just with modified priors and proposal
truevalues:(0.00052, 0.88, 0.00001, 0.09)
startvalues:(0.0008, 0.65, 0.000001, 0.25)
priors: all uninformative, uniform, broad; except beta on fixation
shuffle: on
"        # Notes will get appended to the top of output files
chain <- runMCMC(startvalue, 200000, '22-test', notes)

