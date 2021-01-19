FROM wzhou88/saige:0.38 

MAINTAINER Kenny Westerman <kewesterman@mgh.harvard.edu>

#RUN apt-get update && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
#RUN R -e "install.packages('devtools')"
RUN R -e "options(repos=structure(c(CRAN='https://cloud.r-project.org/'))); devtools::install_version('SPAtest', version='3.0.0')"
RUN R -e "devtools::install_github('WenjianBi/SPAGE', dep=F)"
RUN R -e "devtools::install_github('https://github.com/duytpm16/bgenR', dep=F)"

COPY SPAGE.R /SPAGE.R

#RUN wget http://bitbucket.org/gavinband/bgen/get/master.tar.gz \
#	&& tar xvzf master.tar.gz \
#	&& cd gavinband-bgen-44fcabbc5c38 \
#	&& ./waf configure \
#	&& ./waf 
#ENV BGENIX /gavinband-bgen-44fcabbc5c38/build/apps/bgenix

RUN R -e "options(repos=structure(c(CRAN='https://cloud.r-project.org/'))); install.packages(c('dplyr', 'readr'))"
COPY format_spage_phenos.R format_spage_output.R /
RUN apt-get update && apt-get install tzdata

RUN apt-get update && apt-get install -y dstat atop
