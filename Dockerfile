# get shiny serves plus tidyverse packages image
FROM rocker/rstudio

#system libraries of general use
RUN apt-get update && apt-get install -y \
    software-properties-common \
    libssl-dev \
    build-essential 

#RUN apt-add-repository ppa:marutter/rrutter3.5
#RUN apt-add-repository ppa:marutter/c2d4u3.5
#RUN apt-add-repository ppa:zarquon42/statismo-develop
#RUN sudo apt update --allow-insecure-repositories
#RUN sudo apt install statismo-dev 

# install R packages required
#RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
#RUN R -e "devtools::install_github("zarquon42b/RvtkStatismo",ref="develop"))"

#build will fail with with small RAM allotments -- set aside at least 4GB 


