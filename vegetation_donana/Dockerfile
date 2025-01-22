FROM eco4cast/rocker-neon4cast

# Copy the R scripts and data files into the container
COPY . /home/rstudio/vegetation_donana

# Set the working directory
WORKDIR /home/rstudio/vegetation_donana


# Set the default command to run an R script
RUN Rscript install_libraries.R && echo "script1.R has finished running"
RUN Rscript main_analyses_donana.R
# CMD ["Rscript", "main_analyses_donana.r"]
