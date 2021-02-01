# Replication Slippage Model - Utils 
# File containing supporting functions for the slippage model 
# 
require(expm)
require(stringr)
require(parallel)

# Function to convert ancestral sequences into vectors of slip events
createSlips <- function(anc, len, pos){
  # start out with a base vector containing nchar number of zeros 
  base <- rep(0,nchar(anc)) 
  # remove the gap characters from the ancestral sequence 
  anc <- gsub("-","",anc)

  # if there is no insertion, return the vector of zeros 
  if (len == 0 || is.na(len)){
    return (base)
    
  # if there is an insertion, include the slip count at the appropriate position
  }else{
    pos <- as.numeric(pos) - (len-1)   # pos marks the end of the insertion; adjust by subtracting length-1 
    base[pos] <- len
    return (base)
  }
}

# Receives a vector of sequences and returns a vector of c(A,C,G,T) nucleotide 
# proportions found in the data
estimateFreq <- function(seqs){
  nt <- c("A", "C", "G", "T")
  output <- c()
  for (n in 1:4){
    counts <- sum(unname(sapply(seqs,function(x) str_count(x, nt[n]))))
    output[n] <- counts / sum(unname(sapply(seqs, nchar)))
  }
  output
}

# AIM TO REMOVE -- DOES NOT PROVIDE USE MAKE DEPRECATED 
# Function to convert vector of slip locations c(4,4,4,4,4,4)
# to slip vectors c(0,0,0,6,0,0)
getSlipVector <- function(locs, length){
  
  vect <- rep(0,length)
  if (length(locs) == 0){
    return (vect)
  }else{
    tab <- table(locs)
    vect <- replace(vect, as.numeric(names(tab)), as.numeric(tab))
    return(vect)
  }
}

# Function to convert vector of slip vectors c(0,0,0,6,0,0) 
# to slip locations c(4,4,4,4,4,4)
getSlipLocations <- function(slip){
  # used to go from c(0,0,0,3,0,0,0) to c(4,4,4)
  nonzeros <- which(slip!=0)
  if (length(nonzeros) == 0){
    return (list(loc=c(),len=length(slip)))
  }else{
    return (list(loc=rep(nonzeros, slip[nonzeros]),len=length(slip)))
  }
}

# Function to determine the new tip sequence from the old one using the slip vector
getTip <- function(oldtip, slip){
  
  nonzeros <- which(slip != 0)
  
  # Return immediately if no slip events
  if (length(nonzeros) == 0){
    return(oldtip)
  
  # Modify the tip sequence if one or more slip events
  }else{
    toCopy <- rep(T, nchar(oldtip))     # default vector of T values 
    tip.chars <- strsplit(oldtip, "")[[1]]   # nucleotide sequence (split up)
    s.loc <- getSlipLocations(slip)[[1]]     # returns the locations of all slip events 
    tab <- table(s.loc)
    
    # for every distinct slip event specified in the slip vector, 
    for (slip in 1:length(tab)){
      start <- as.numeric(names(tab)[slip])
      stop <- start + (tab[[slip]] - 1)
      
      # If either stop or start exceeds the boundary
      if (start > nchar(oldtip) || stop > nchar(oldtip)){
        start <- nchar(oldtip) - (stop - start)
        stop <- nchar(oldtip)
      }
      # re-adjustment of positions 
      # perform this action on all indices after the first
      if (slip > 1 && !toCopy[start]){
        adjust <- min(which(toCopy[start:length(toCopy)])) - 1
        start <- start + adjust
        stop <- stop + adjust
      }
      #print(start)
      #print(stop)
      toCopy[start:stop] <- F
    }
    #tip.chars[toCopy] <- "-"
    return(paste0(tip.chars[toCopy],collapse=""))
  }
}

# Function for drawing from discrete normal values
# Used for moving rearranging slip positions
delta <- function(sd=2){
  x <- rnorm(1,mean=0,sd=sd)
  if (x < 0){
    x <- abs(x)
    x <- ceiling(x)
    x <- -x
  }else{
    x <- ceiling(x)
  }
  x
}

# Returns the F81 transition probability matrix
# Depends on global variable FREQ c(A,C,G,T)
getMat <- function(rate, branch){
  mat <- matrix(rep(FREQ, each=4), nrow=4, ncol=4,dimnames=list(nt,nt))  # set up the matrix
  
  inv.freq <- sapply(1:4, function(x) 1-FREQ[x])    # save the inverse of each frequency
  
  mat <- mat * rate         # multiply by rate 
  
  diag(mat) <- sapply(1:4, function(x) -(rate*inv.freq[x]))  # adjust the diagonals so that every row sums to 0
  
  mat <- branch * mat    # multiply by branch length
  
  # exponentiate and return
  expm(mat)
}

pairllh <- function(anc, newtip, rate, branch){
  tmat <- getMat(rate,branch)
  achars <- strsplit(anc, "")[[1]]
  tchars <- strsplit(newtip, "")[[1]]

  # Determines the F81 log-likelihood of observing this pair of sequences
  result <- mapply(function(tchar,achar){
    # initializes simple row vector with 1 for the given nucleotide, 0 for the others
    tip.llh <- matrix(as.numeric(tchar == nt), nrow=4,ncol=1,dimnames=list(nt))
    
    # finalize the calculation for tip likelihood
    # dot product
    llh <- tmat %*% tip.llh    # tip.llh = c(1,0,0,0) for A
    llh <- llh * f
    if (!achar %in% nt){
      print(achar)
    }
    # likelihood given the exact nucleotide (state) that we see in the ancestor
    final.llh <- llh[achar,]
    return(log(final.llh))
  }, achars, tchars)
  
  sum(result)
}

# likelihood of the entire slip.list
# HIGHEST COMPLEXITY (very slow) -- only calculated when: 
  # a) rate parameter is changed
  # b) before the MCMC starts for the first iteration
seqllh <- function(rate, slip.list, cores=8){
  tip.seqs <- unname(mcmapply(getTip, indels$tip, slip.list, mc.cores=cores))
  rate <- rep(rate, length(anc.seqs))
  unname(mcmapply(pairllh, anc.seqs, tip.seqs, rate, branches, mc.cores=cores))
}

# Custom function to process insertions
insert <- function(str, indel,pos ){
  vect <- strsplit(str,"")[[1]]
  if (pos < 1){
    return (NA)
  }
  
  if (pos == 1){    # Add before the first nucleotide
    return (paste(c(indel, vect[1:length(vect)]),collapse=""))
  }else if (pos == length(vect)+1){   # Add after the last nucleotide
    return (paste(c(vect[1:length(vect)], indel),collapse=""))
  }else{
    return (paste(c(vect[1:pos-1], indel, vect[pos:length(vect)]),collapse=""))
  }
}

# Custom function to process deletions 
delete <- function(str, indel, pos){
  vect <- rep(T, nchar(str))
  
  if (pos > (nchar(str))){
    return (NA)
  }
  
  start <- pos - nchar(indel) + 1
  end <- pos
  
  vect[start:end] <- F
  
  chars <- strsplit(str, "")[[1]]
  
  return(paste(chars[vect],collapse=""))
}
