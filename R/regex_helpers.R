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
  rs = ""
  for(ci in 1:length(v)){
    c = v[ci]
    if(ci == 1){
      rs = paste0(".*", c, ".*")
    } else{
      rs = paste(c(rs, paste0(".*", c, ".*")), collapse = "|")
    }
    uv = unlist(strsplit(c, " ")) # num space-sep units
    # for each unit use lower- and uppercase
    uvstr = c(paste(uv, collapse = " "))
    uvl = list(uv)
    for(i in 1:length(uv)){
      uvi = c()
      for(ui in 1:length(uv)){
        chari = uv[ui]
        if(ui <= i){
          if(nchar(chari)>1){
            ssi = paste0(toupper(substr(chari, 1, 1)),
                         substr(chari, 2, nchar(chari)))
          } else{
            ssi = paste0(toupper(substr(chari, 1, 1)))
          }
        }
        else{
          ssi = chari
        }
        uvi = c(uvi, ssi)
      }
      uvl[[length(uvl)+1]] = uvi
    }
    # append to new str
    for(si in 1:length(uvl)){
      s = uvl[[si]]
      if(length(uv) > 1){
        if(!si==1){
          # space sep
          rs = paste(c(rs, paste0(".*", paste(s, collapse = " "), ".*")), collapse = "|")
        }
        # underline sep
        rs = paste(c(rs, paste0(".*", paste(s, collapse = "_"), ".*")), collapse = "|")
        # dash sep
        rs = paste(c(rs, paste0(".*", paste(s, collapse = "-"), ".*")), collapse = "|")
      } else{
        if(!si==1){
          rs = paste(c(rs, paste0(".*", s, ".*")), collapse = "|")
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
#' @param m Metadata used for pattern matching (data.frame).
#' @param filtrel Logical symbol joining each regex pattern (default "|").
#' @param nfilt Whether to also use negative lookup
#' @param ntfilt Specific character strings of terms to filter on, or for which
#' to use the negative lookup/regex negation (default "").
#' @param ptfilt Additional positive/affirmative match terms to filter (default "").
#' @returns The result of assessing a regex match on a metadata variable.
#' @seealso appendvar
#' @export
get_filt <- function(v, m = md, filtrel = "|", nfilt = FALSE, ntfilt = "", ptfilt = "",
                     varl = c("gsm_title", "sample_type", "disease_state", 
                              "anatomic_location", "misc")){
  if(!filtrel %in% c("|", "&")){
    message("Please provide a valid filter relation symbol.")
    return(NULL)
  }
  # positive match filter
  if(ptfilt == ""){
    filtl = grepl(v, m[,colnames(m)==varl[1]])
  } else{
    filtl = grepl(get_pst(ptfilt), m[,colnames(m)==varl[1]])
    filtl = grepl(v, m[,colnames(m)==varl[1]])
  }
  # negative match filter
  if(nfilt){
    message("Using negative lookup/exclusion filter...")
    nfiltv = get_pstr_neg(v)
    filtl = filtl & !grepl(nfiltv, m[,colnames(m)==varl[1]])
  }
  # term filter
  if(!ntfilt == ""){
    message("Using term lookup filter...")
    filtl = filtl & !grepl(get_pstr(ntfilt), m[,colnames(m)==varl[1]])
  }
  # proceed if additional vars specified
  if(length(varl)>1){
    for(vi in varl[2:length(varl)]){
      if(filtrel == "|"){
        filtl = filtl | grepl(v, m[,colnames(m)==vi])
        if(nfilt){
          filtl = filtl & !grepl(nfiltv, m[,colnames(m)==vi])
        }
      }
      if(filtrel == "&"){
        filtl = filtl | grepl(v, m[,colnames(m)==vi])
        if(nfilt){
          filtl = filtl & !grepl(nfiltv, m[,colnames(m)==vi])
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
#' @param var The variable in the metadata to assess.
#' @param val Character string for regex matching.
#' @param filtv Vector of character strings to filter on.
#' @param m Metadata for term lookup (data.frame).
#' @returns The result of assessing a regex match on a metadata variable.
#' @seealso get_filt, get_pstr_neg, get_pstr
#' @export
appendvar <- function(var, val, filtv, m){
  varr = m[, colnames(m) == var]
  # get composite filter
  filti = !grepl(paste0("(^|;)", val, "(;|$)"), varr); compfilt = filti & filtv
  # assess filter results
  if(length(compfilt[compfilt]) == 0){
    message("No unique values to append. Returning var unchanged.")
    return(varr)
  } else{
    varr[compfilt] = ifelse(varr[compfilt] == "NA", val,
                            paste(varr[compfilt], val, sep = ";")
    )
    message("Appended n = ", length(varr[compfilt]), " values")
    return(varr)
  }
  return(NULL)
}
