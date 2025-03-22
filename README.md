![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)


![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/ericlarg4/eps2fold/docker-publish.yml?style=for-the-badge)
![GitHub last commit](https://img.shields.io/github/last-commit/ericlarg4/eps2fold?style=for-the-badge)


![GitHub License](https://img.shields.io/github/license/ericlarg4/eps2fold?style=for-the-badge)

# &epsilon;2Fold

## Launching the app

The app can be used on a web server (go to section 1), or run locally through:

* your own web browser (not install required; section 2)
* an R interpreter (R install necessary, GUI optional; section 3)
* in a Docker container. (docker install necessary, GUI optional; section 4).

All dependencies will be installed automatically where necessary.

### 1. Run the app online

Acess the web app on the [Connect Posit Cloud server](https://ericlarg4-eps2fold.share.connect.posit.cloud/).

### 2. Run the app locally in your web browser

Launch the app in your browser by opening the following web page: [https://ericlarg4.github.io/Eps2Fold/](https://ericlarg4.github.io/Eps2Fold/).

This method is slower, particularly in Firefox - [for some reason](https://github.com/posit-dev/shinylive/issues/191).

### 3. Run the app directly in R

*Requires to have R installed: [https://www.r-project.org/](https://www.r-project.org/)*

*The use of a GUI like [Rstudio](https://posit.co/download/rstudio-desktop/) or [Positron](https://positron.posit.co/) is advised.*

Clone/download the repository then launch `app.R` through the R console or a GUI (see below)

On first use, dependencies will automatically be installed, which may take some time.

#### 3.1. GUI

Open `app.R` then click on "Run App" or "Run Shiny App".

#### 3.2. Terminal

The following instructions work for Windows users. Other OS may need adjustments.

1. If R was not added in PATH, locate the installation path, typically `& "C:\Program Files\R\R-X.Y.Z\"`, where X.Y.Z is the R version number.
2. Navigate to the app directory, i.e. Eps2Fold if it was cloned from GitHub
3. Launch the app in your default web browser with: `& "R_PATH\bin\x64\Rscript.exe" -e "shiny::runApp(launch.browser = TRUE)"`, where `R_PATH` is the path determined at step 1.  For instance: `& "C:\Program Files\R\R-4.4.3\bin\x64\Rscript.exe" -e "shiny::runApp(launch.browser = TRUE)"`.

### 4. Run the app in a Docker container

*Requires to have installed Docker: [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/).*

*Docker desktop may be installed, but this is not required.*

*Does not require to have R installed.*

For the first use or to update the app, pull the image, otherwise proceed directly to 4.2.

#### 4.1 Pulling the image

In a terminal, pull the Docker image: `docker pull ghcr.io/ericlarg4/eps2fold:main`

#### 4.2 Running the app

Proceed to either 2.2.1 or 2.2.2 if you prefer to use docker desktop.

##### 4.2.1 Terminal

1. Launch the app: `docker run --rm -p 3838:3838 ghcr.io/ericlarg4/eps2fold:main`
2. Open your browser and go to [http://localhost:3838/](http://localhost:3838/)
3. To stop and remove the process, simply type ctrl+C from the terminal.

##### 4.2.2 Docker desktop program

1. In the "Images" tab, run the image by clicking on the "play" button, then indicate *3838* as *Host port* in the optional settings.
2. In the Containers/Apps tab, click on the `Open in browser` icon.
3. Open your browser and go to [http://localhost:3838/](http://localhost:3838/)
4. To stop the process, click on the stop button
5. Optionally, to remove the container (the image will remained available), click on the trash button.

If you want to stop and remove the container once you're done with the app, click on the remove icon.

## Use of the app

See the wiki.
