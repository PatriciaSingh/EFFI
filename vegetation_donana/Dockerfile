FROM rocker/r-ver:4.1.2

# Copy the R scripts and data files into the container
COPY . /home/rstudio/vegetation_donana

# Set the working directory
WORKDIR /home/rstudio/vegetation_donana

RUN apt-get update && apt-get install -y \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libcurl4-openssl-dev \
    libsqlite3-dev \
    libudunits2-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*



# Set the default command to run an R script
RUN Rscript install_libraries.R && echo "script1.R has finished running"
RUN Rscript main_analyses_donana.R
# CMD ["Rscript", "main_analyses_donana.r"]
