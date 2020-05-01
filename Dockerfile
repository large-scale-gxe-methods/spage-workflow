FROM rocker/r-ubuntu:20.04

MAINTAINER Kenny Westerman <kewesterman@mgh.harvard.edu>

RUN apt-get update && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
RUN R \
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_version('SPAtest', version='3.0.0')"
RUN R -e "devtools::install_github('WenjianBi/SPAGE', dep=F)"

RUN apt-get update && apt-get install -y dstat atop
