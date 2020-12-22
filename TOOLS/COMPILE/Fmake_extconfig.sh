#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
# ===============
# Fmake_extconfig.sh
# ===============
# ---------------
# Make the directories for an externally supported configuration
# and base the initial versions on the nearest equivalent from the
# reference set (as named in uspcfg.txt)
# ---------------
# SYNOPSIS
# ========
# ::
#  $ Fmake_extconfig.sh
# DESCRIPTION
# ===========
# - Make the config directory 
# - Create repositories needed :
#  - EXP00 for namelist
#  - MY_SRC for user sources
#  - BLD for compilation 
# EXAMPLES
# ========
# ::
#  $ ./Fmake_extconfig.sh CONFIG_NAME REF_CONFIG_NAME 
# TODO
# ====
# option debug
# EVOLUTIONS
# ==========
# $Id: Fmake_extconfig.sh 3715 2012-11-28 16:06:02Z acc $
#   * creation
#-
\mkdir  ${1}
\mkdir  ${1}/EXP00
\mkdir  ${1}/MY_SRC
\cp -R  ${2}/cpp_${2}.fcm ${1}/cpp_${1}.fcm
\cp -R  ${2}/EXP00/*namelist* ${1}/EXP00/.
\cp -R  ${2}/EXP00/*.xml ${1}/EXP00/.
[ -f ${2}/EXP00/AGRIF_FixedGrids.in ] &&  \cp -R  ${2}/EXP00/AGRIF_FixedGrids.in ${1}/EXP00/.
[ -d    ${2}/MY_SRC ]                 &&  \cp  ${2}/MY_SRC/* ${1}/MY_SRC/. 2> /dev/null
