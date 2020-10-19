#!/usr/bin/env R

# Utilities used by scripts.

#' Get index blocks, with remainder handling
#'
#' Get index blocks, with remainder handling. This is used by scripts
#' that process large data objects processively, in chunks/blocks of 
#' sample indices. Note the default block size of 50 is recommended for
#' processes using entire array's worth of data (e.g. from an RGChannel set,
#' getRed, getGreen, getMeth, getUnmeth, etc.), as test showed diminished
#' performance with greater or lesser samples on typical compute systems.
#'
#' @param slength Length of the index vector (integer)
#' @param bsize The maximum block size (integer, default = 50).
#' @returns List of index blocks.
#' @export
getblocks <- function(slength, bsize = 50){
  iv <- list()
  if(slength < bsize){
    iv[[1]] <- seq(1, slength, 1)
  } else{
    sc <- 1; ec <- sc + bsize - 1
    nblocks <- slength %/% bsize
    for(b in 1:nblocks){
      iv[[b]] <- seq(sc, ec, 1)
      sc <- ec + 1; ec <- ec + bsize
    }
    # add final indices
    if(nblocks < (slength/bsize)){
      iv[[length(iv) + 1]] <- seq(sc, slength, 1)
    }
  }
  return(iv)
}
