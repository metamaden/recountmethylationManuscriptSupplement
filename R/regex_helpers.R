#!/usr/bin/env R

# Helper functions for regular expression term matching and mapping

#' Get regex patterns for a string, or the affirmative match for the 
#' given string
#'
#' Gets the regex pattern matching the affirmative of a given string `v`.
#' Does progressive capitalization on values separated by spaces
#' for each value in v, appends flanking '.*' (matches any char)
#' for each value in v, appends "|" OR conditional separator
#'
#' @param v Character string or vector of such character strings. 
#' @returns Regex patterns for detecting matching strings from metadata.
#' @seealso get_pstr_neg, get_filt, appendvar
#' @export
get_pstr <- function(v){
  rs <- ""
  for(ci in seq(length(v))){
    c <- v[ci]
    if(ci == 1){rs = paste0(".*", c, ".*")} else{
      rs = paste(c(rs, paste0(".*", c, ".*")), collapse = "|")
    }
    uv <- unlist(strsplit(c, " ")) # num space-sep units
    # for each unit use lower- and uppercase
    uvstr <- c(paste(uv, collapse = " ")); uvl <- list(uv)
    for(i in seq(length(uv))){
      uvi <- c()
      for(ui in seq(length(uv))){
        chari <- uv[ui]
        if(ui <= i){
          if(nchar(chari)>1){
            ssi <- paste0(toupper(substr(chari, 1, 1)),
                         substr(chari, 2, nchar(chari)))
          } else{
            ssi <- paste0(toupper(substr(chari, 1, 1)))
          }
        }
        else{ssi <- chari}
        uvi <- c(uvi, ssi)
      }; uvl[[length(uvl)+1]] <- uvi
    }
    # append to new str
    for(si in 1:length(uvl)){
      s <- uvl[[si]]
      if(length(uv) > 1){
        if(!si==1){
          # space sep
          rs <- paste(c(rs, paste0(".*", paste(s, collapse = " "), ".*")), collapse = "|")
        }
        # underline sep
        rs <- paste(c(rs, paste0(".*", paste(s, collapse = "_"), ".*")), collapse = "|")
        # dash sep
        rs <- paste(c(rs, paste0(".*", paste(s, collapse = "-"), ".*")), collapse = "|")
      } else{
        if(!si==1){
          rs <- paste(c(rs, paste0(".*", s, ".*")), collapse = "|")
        }
      }
    }
  }
  return(rs)
}

#' Get the negation regex patterns for a given regex pattern
#'
#' Does progressive capitalization on values separated by spaces
#' for each value in v, appends flanking '.*' (matches any char)
#' for each value in v, appends "|" OR conditional separator
#'
#' @param pstr Regex patterns for matching, or output from `get_pstr()`.
#' @returns Regex for negation of indicated patterns.
#' @seealso get_pstr, get_filt, appendvar
#' @export
get_pstr_neg <- function(pstr){
  pstrg = gsub("\\.\\*", "_", pstr); pstrg = gsub("\\|", "", pstrg)
  uv = unlist(strsplit(pstrg, "_")); uv = uv[!uv==""]
  for(ui in 1:length(uv)){
    s = uv[ui]
    if(ui == 1){
      ns = paste(paste0(".*non-", s, ".*"), paste0(".*Non-", s, ".*"), paste0(".*non ", s, ".*"),
                 paste0(".*Non ", s, ".*"), paste0(".*not ", s, ".*"), paste0(".*Not ", s, ".*"),
                 paste0(".*", s, "-free.*"), paste0(".*", s, "-Free.*"), paste0(".*", s, " free.*"),
                 paste0(".*", s, " Free.*"), sep = "|")
    } else{
      ns = paste(ns, paste0(".*non-", s, ".*"), paste0(".*Non-", s, ".*"), paste0(".*non ", s, ".*"),
                 paste0(".*Non ", s, ".*"), paste0(".*not ", s, ".*"), paste0(".*Not ", s, ".*"),
                 paste0(".*", s, "-free.*"), paste0(".*", s, "-Free.*"), paste0(".*", s, " free.*"),
                 paste0(".*", s, " Free.*"), sep = "|")
    }
  }
  return(ns)
}

#' Get the outcome of pattern matching a string with regex patterns
#'
#' Does pattern matching for regex patterns constructed from a character
#' string or vector of such character strings. 
#'
#' @param v Character string or vector of such strings.
#' @param m Preprocessed metadata used for pattern matching/lookups (data.frame, mdpre).
#' @param filtrel Logical symbol joining each regex pattern (default "|").
#' @param ntfilt Regex pattern corresponding to negative lookup filter (default NULL).
#' @param ptfilt Regex pattern corresponding to positive lookup filter (default NULL).
#' @returns The result of assessing a regex match on a metadata variable.
#' @seealso appendvar
#' @export
get_filt <- function(v, m = mdpre, filtrel = "|", ntfilt = NULL, ptfilt = NULL,
                     varl = c("gsm_title", "sample_type", "disease_state", 
                              "anatomic_location", "misc")){
  # positive lookup filter
  filtl <- grepl(v, m[,varl[1]])
  if(!is.null(ptfilt)){
    message("Using positive lookup filter with ptfilt: ", ptfilt)
    filtl <- filtl & grepl(get_pstr(ptfilt), m[,varl[1]])
  }
  # negative lookup filter
  if(!is.null(ntfilt)){
    message("Using negative lookup filter with ntfilt: ", ntfilt)
    nfiltv <- get_pstr_neg(ntfilt)
    filtl <- filtl & !grepl(nfiltv, m[,varl[1]])
  }
  # proceed if additional vars specified
  if(length(varl)>1){
    if(!filtrel %in% c("|", "&")){
      message("Please provide a valid filter relation symbol.")
      return(NULL)
    } else{
      for(vi in varl[2:length(varl)]){
        if(filtrel == "|"){
          filtl <- filtl | grepl(v, m[,vi])
          if(!is.null(ntfilt)){filtl <- filtl & !grepl(nfiltv, m[,vi])}
        } else if(filtrel == "&"){
          filtl <- filtl & grepl(v, m[,vi])
          if(!is.null(ntfilt)){filtl <- filtl & !grepl(nfiltv, m[,vi])}
        }
      }
    }
  }
  return(filtl)
}

#' Append mapped terms to a metadata variable
#'
#' Appends new variable data to a metadata variable, preserving current
#' variable terms.
#'
#' @param var Variable in m for which to append new terms.
#' @param val Character string to append to matched samples.
#' @param filtv Vector of boolean values identifying samples for which to
#'  append the new term.
#' @param m The postprocessed metadata for which to append the new terms
#'  (data.frame, mdpost).
#' @returns The result of appending new terms to specified var.
#' @seealso get_filt, get_pstr_neg, get_pstr
#' @export
appendvar <- function(var, val, filtv, m = mdpost){
  varr <- m[, var]
  # get composite filter
  filti <- !grepl(paste0("(^|;)", val, "(;|$)"), varr)
  compfilt <- filti & filtv
  # assess filter results
  if(length(compfilt[compfilt]) == 0){
    message("No unique values to append. Returning var unchanged.")
    return(varr)
  } else{
    varr[compfilt] <- ifelse(varr[compfilt] == "NA", val,
                            paste(varr[compfilt], val, sep = ";")
    )
    message("Appended n = ", length(varr[compfilt]), " values")
    return(varr)
  }
  return(NULL)
}

#' Match and combine 2 datasets
#'
#' Orders 2 datasets on corresponding variables, then appends an NA
#' matrix as needed before binding on columns. This is mainly used to 
#' combine different types of sample metadata on common/shared GSM IDs.
#' Note columns should be matchable after calling `as.character()` for
#' each.
#'
#' @param d1 First set to combine (matrix).
#' @param d2 Second set to combine (matrix).
#' @param ci1 Column index in first dataset to match on
#' @param ci2 Column index in second dataset to match on
#' @returns Union of the 2 datasets, including NA values where appropriate.
#' @export
match1to2 <- function(d1, d2, ci1 = 2, ci2 = 1){
  id.all <- unique(c(d1[, ci1], d2[, ci2]))
  id1 <- id.all[!id.all %in% d1[, ci1]]
  id2 <- id.all[!id.all %in% d2[, ci2]]
  # append na slices as necessary
  if(length(id1) > 0){
    nav <- rep(rep("NA", length(id1)), ncol(d1) - 1)
    mna <- matrix(c(id1, nav), nrow = length(id1), ncol = ncol(d1))
    d1 <- rbind(d1, mna)
  }
  if(length(id2) > 0){
    nav <- rep(rep("NA", length(id2)), ncol(d2) - 1)
    mna <- matrix(c(id2, nav), nrow = length(id2), ncol = ncol(d2))
    d2 <- rbind(d2, mna)
  }
  # reorder and assign title var
  match.id1 <- match(as.character(d1[, ci1]), as.character(d2[, ci2]))
  order.id1 <- order(match.id1); d1 <- d1[order.id1,]
  match.id2 <- match(as.character(d2[, ci2]), as.character(d1[, ci1]))
  order.id2 <- order(match.id2); d2 <- d2[order.id2,]
  cond <- identical(as.character(d2[, ci2]), as.character(d1[, ci1]))
  if(cond){return(cbind(d1, d2))} else{
    stop("there was an issue matching the final datasets...")
  }
  return(NULL)
}