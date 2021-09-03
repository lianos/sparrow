FROM bioconductor/bioconductor_docker:devel
LABEL maintainer="Steve Lianoglou (slianoglou@gmail.com)"

# Newer versions of bioconductor_docker image sets options::repos to a
# frozen RSPM date. This is good for reproducibility at large, but I like
# to live on the edge/latest.
RUN sudo echo 'options(repos = c(CRAN = "https://packagemanager.rstudio.com/all/__linux__/bionic/latest"))' >> /usr/local/lib/R/etc/Rprofile.site

# Required for textshaping library
RUN sudo apt-get update -y
RUN sudo apt-get install -y \
  libfribidi-dev \
  libmagick++-dev \
  libv8-dev \
  && sudo rm -rf \
    /var/lib/apt/lists/ \
    /tmp/downloaded_packages/ \
    /tmp/*.rds

ADD packages.txt /tmp/
ADD install.R /tmp/

# Installation of 'foundational' compiled packages / language bridged packages
# can be wonky if we install the prepackaged stuff, so we (re)compile this
# by source as many of the downstream packages we end up relying on uses
# header libraries from here. So also this related note about pre-built packages
# and RSPM.
#   https://community.rstudio.com/t/package-manager-and-travis/73816
RUN Rscript -e "\
  install.packages(c(\
      'Rcpp', \
      'cpp11', \
      'igraph', \
      'magick', \
      'mclust', \
      'V8'), \
    repos = 'https://cran.rstudio.com', \
    type = 'source')"

# Use the default RSPM mojo to install the rest of the packages
RUN R -f /tmp/install.R

# Now that we've got the base packages here, we will want to update them
# every now and again without reinstalling all of them again, so we can
# do that here. Run a no-op command to flip this on/of/rerun
# RUN pwd
RUN Rscript -e "\
  options(repos = 'https://cran.rstudio.com/'); \
  BiocManager::install(update = TRUE, ask = FALSE)" \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e "\
  BiocManager::install(c(\
    'lianos/sparrow', \
    'lianos/sparrow.shiny'))" \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Init command for s6-overlay
CMD ["/init"]
