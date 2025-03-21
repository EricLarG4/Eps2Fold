# Use the latest Rocker Shiny image
FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages required by the app
RUN R -e "install.packages(c('tidyverse', 'readxl', 'data.table', 'magrittr', 'shiny', 'bslib', 'bsicons', 'thematic', 'ggrepel', 'ggthemes', 'ggtext', 'patchwork', 'plotly', 'DT', 'FactoMineR', 'cluster', 'FactoInvestigate'), dependencies = TRUE)"

# Copy the app files to the container
WORKDIR /srv/shiny-server
COPY . .

# Ensure correct permissions
RUN chmod -R 755 /srv/shiny-server

# Expose the Shiny port
EXPOSE 3838

# Start the Shiny app
CMD ["/usr/bin/shiny-server"]

#end
