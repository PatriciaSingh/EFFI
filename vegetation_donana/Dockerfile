FROM rocker/r-ver:4.1.2

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    gdal-bin \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install additional R packages
RUN R -e "install.packages(c('ows4R', 'sf', 'terra', 'dplyr', 'ggplot2', 'readr', 'jagsUI', 'coda'), repos='http://cran.rstudio.com/')"

# Copy the R scripts and data files into the container
COPY . /home/rstudio/vegetation_donana

# Set the working directory
WORKDIR /home/rstudio/vegetation_donana

# Set the default command to run an R script
CMD ["Rscript", "main_analyses_donana.r"]
