![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)

![CodeFactor Grade](https://img.shields.io/codefactor/grade/github/ericlarg4/eps2fold?style=for-the-badge)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/ericlarg4/eps2fold/docker-publish.yml?style=for-the-badge)
![GitHub last commit](https://img.shields.io/github/last-commit/ericlarg4/eps2fold?style=for-the-badge)

![GitHub License](https://img.shields.io/github/license/ericlarg4/eps2fold?style=for-the-badge)

# &epsilon;2Fold

## Use of the app
Please refer to the [Wiki](https://github.com/EricLarG4/Eps2Fold/wiki).

## Launching the app

The app can be used on a [Web server](#1-run-the-app-online), or run locally through:

* [Your own web browser](#2-run-the-app-locally-in-your-web-browser) (not install required)
* [An R interpreter](#3-run-the-app-directly-in-r) (R install necessary, GUI optional)
* [In a Docker container](#4-run-the-app-in-a-docker-container) (docker install necessary, GUI optional.

All dependencies will be installed automatically where necessary.

### 1. Run the app online

Access the web app on the [Connect Posit Cloud server](https://ericlarg4-eps2fold.share.connect.posit.cloud/).

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
3. Launch the app in your default web browser with:

```PowerShell
& "R_PATH\bin\x64\Rscript.exe" -e "shiny::runApp(launch.browser = TRUE)"
```

where `R_PATH` is the path determined at step 1. For example:

```PowerShell
& "C:\Program Files\R\R-4.4.3\bin\x64\Rscript.exe" -e "shiny::runApp(launch.browser = TRUE)"
```

### 4. Run the app in a Docker container

*Requires to have installed Docker: [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/).*

*Docker desktop may be installed, but this is not required.*

*Does not require to have R installed.*

For the first use or to update the app, pull the image, otherwise proceed directly to [4.2](https://github.com/EricLarG4/Eps2Fold/edit/main/README.md#42-running-the-app).

#### 4.1 Pulling the image

In a terminal, pull the Docker image:

```PowerShell
docker pull ghcr.io/ericlarg4/eps2fold:main
```

#### 4.2 Running the app

Proceed to either [4.2.1](https://github.com/EricLarG4/Eps2Fold/edit/main/README.md#421-terminal) for terminal-based launch, or [4.2.2](https://github.com/EricLarG4/Eps2Fold/edit/main/README.md#422-docker-desktop-program) if you prefer to use docker desktop.

##### 4.2.1 Terminal

1. Launch the app:

```PowerShell
docker run --rm -p 3838:3838 ghcr.io/ericlarg4/eps2fold:main
```
  
2. Open your browser and go to [http://localhost:3838/](http://localhost:3838/)
3. To stop and remove the process, simply type ctrl+C from the terminal.

##### 4.2.2 Docker desktop program

1. In the "Images" tab, run the image by clicking on the "play" button

<img width="948" height="91" alt="image" src="https://github.com/user-attachments/assets/5b4df8f0-caed-4832-9bcc-e7eef788a677" />

2. indicate *3838* as *Host port* in the optional settings and click on *Run*

<img width="596" height="573" alt="image" src="https://github.com/user-attachments/assets/f1ba0d12-fce8-4142-8ebf-fb047f428d48" />

3. Click on 3838:3838 on the container page that opened (you can also find it on the *Containers* page

<img width="1062" height="347" alt="image" src="https://github.com/user-attachments/assets/948ea6ed-1082-4962-b171-561bffa86945" />

4. To stop the process, click on the stop button
5. Optionally, to remove the container (the image will remained available), click on the trash button.

## Use of the app

See the wiki.

## Changelog

### Version 1.2.0

Eps2Fold can now work in CD or Eps2Fold-only modes for user predictions.
Eps2Fold deals with empty or missing user data sheet by displaying the training set output only.
