# Base R Shiny image
FROM rocker/shiny

USER root

# Make a directory in the container
RUN mkdir /home/rnaCARD_shiny

# Install R dependencies
RUN Rscript -e "install.packages(c('ggplot2', 'readr', 'shinyWidgets', 'shinyjs', 'future', 'multiprocess', 'promises', 'shinycssloaders', 'DT'), repos='http://cran.us.r-project.org/')"

RUN apt-get update && apt-get install -y libcurl4-openssl-dev libssl-dev

RUN Rscript -e "install.packages(c('ggthemes', 'plotly'), repos='http://cran.us.r-project.org/')"

RUN Rscript -e "install.packages(c('future'), repos='http://cran.us.r-project.org/')"


# Copy the Shiny app code
COPY rnaCARD_shiny  /home/rnaCARD_shiny

# Expose the application port
EXPOSE 8180

# Run the R Shiny app
CMD ["R", "-e", "shiny::runApp('/home/rnaCARD_shiny', host = getOption('shiny.host', '0.0.0.0'), port=8180)"]
