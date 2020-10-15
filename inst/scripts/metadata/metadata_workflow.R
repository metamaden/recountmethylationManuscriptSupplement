#!/usr/bin/env R

# Author: Sean Maden
# 
# This is the main sample metadata workflow script. It describes the
# steps to mine, map, and learn sample metadata. This script assumes 
# no SOFT files have been downloaded, and shows the process from data
# download through to metadata mapping and DNAm model predictions.
#
# Note this script manages a full workflow to mine, map, and learn
# sample metadata from GSE SOFT files and DNAm.
# 
# The first step, obtaining the GSE SOFT files, is performed using 
# the `recount-methylation-server`, which has a number of dependencies,
# including the Entrez Utilities tools, Python3, the celery Python 
# library, MongoDB, and SQLite.

#------------------------------
# assign entity names and paths
#------------------------------

# server script paths
server.dir <- "recount-methylation-server"
serverpy.path <- file.path(server.dir, "src", "server.py")
processsoft.path <- file.path(server.dir, "src", "process_soft")

# supplement package path
pkg.name <- "recountmethylationManuscriptSupplement"
md.scripts.dir <- system.file("scripts", "metadata", package = pkg.name)
tables.scripts.dir <- system.file("scripts", "tables", package = pkg.name)

# new directories
soft.path <- "gse_soft"
gsmjson.path <- "gsm_json"
gsmjsonfilt.path <- paste0(gsmjson.path, "_ filt")
gsescripts.path <- "gse_tables"

#------------------------
# download gse soft files
#------------------------
# run the server/download manager
runserver <- system(paste0("python3 ", serverpy.path))

#-----------------
# extract gsm json
#-----------------
# get gsm metadata as json files

# from a python3 session:

# import os
# import sys
# sys.path.insert(0, os.path.join("recount-methylation-server","src"))
# import process_soft
# 
# process_soft.expand_soft()
# process_soft.extract_gsm_soft()
# process_soft.gsm_soft2json()

# get filtered json files
source(file.path(md.scripts.dir, "jsonfilt.R"))

# get sample titles from final datasets
source(file.path(md.scripts.dir, "get_gsm_titles.R"))

#-------------------------------------
# run metasra-pipeline and map results
#-------------------------------------
# from a python3 session:
# process_soft.msrap_screens(nscreensi=200, nmaxscreens=35, qcprint=True)

# 

#-------------------------
# get metadata tables list
#-------------------------
source(file.path(md.scripts.dir, "make_gse_annolist.R"))

#------------------------------
# process soft-derived metadata
#------------------------------
# preprocessing
source(file.path(md.scripts.dir, "metadata_preprocessing.R"))

# do postprocessing
source(file.path(md.scripts.dir, "metadata_postprocessing.R"))

#----------------
# append new data 
#----------------
# append data 
source(file.path(tables.scripts.dir, "tableS1.R"))

#----------------
# save and return
#----------------
