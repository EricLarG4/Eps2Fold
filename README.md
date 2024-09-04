# G4_IDS

## Lauching the app

The app can be used natively in R, or in a Docker container. Either way, there is no need to install extra software/packages/dependencies appart from R or Docker themselves.

It is also possible to test the app in the cloud wihtout installing anything. It only requires a Docker account.

### 1. Run the app directly in R

*Requires to have R installed: https://www.r-project.org/* 

*A GUI like Rstudio is not necessary but can be used.*

Clone/download the repository (the Docker folder is not necessary) then proceeds through the R console or a GUI.

#### 1.1. GUI

In RStudio open ui.R and/or server.R then click on Run App.

#### 1.2. R Console

- Launch the R console.

*Windows users: R.exe is located somewhere like: `C:\Program Files\R\R-4.0.5\bin`*

- Set the working directory to where the repository has been saved: `setwd("Path_to_IDS_folder")`, where `Path_to_IDS_folder` is the path to the repository.

*Windows users: make sure to use / and not \ in the folder path*

- Load and run the app with:

```
source("ui.R")
source("server.R")
shiny::shinyApp(ui, server)
```
or open app.R in RStudio and click on the Run button.

### 2. Run the app in a Docker container

*Does not require to have R installed*

*Requires to have installed Docker: https://docs.docker.com/get-docker/*

For the first use or to update the app, pull the image, otherwise proceed directly to 2.2.

#### 2.1 Pulling the image

In a terminal, pull the Docker image: `docker pull ericlarg4/shiny-ids-app:release` 

#### 2.2 Running the app

Proceed to either 2.2.1 or 2.2.2 if you prefer to use a GUI.

##### 2.2.1 Terminal

- Launch the app: `docker run --rm -p 3838:3838 ericlarg4/shiny-ids-app:release`

- Open your browser and go to `http://localhost:3838/`

If you want to stop and remove the container once you are done with the app, use:

`docker ps -q --filter ancestor=ericlarg4/shiny-ids-app:release | xargs docker stop`

##### 2.2.2 Docker desktop program

- In the image tab, run the image using 3838 as Local Host in the optional settings.

- In the Containers/Apps tab, click on the `Open in browser` icon.

If you want to stop and remove the container once you're done with the app, click on the remove icon.

### 3. Test the app in a Docker playground

*Requires a free Docker account: https://hub.docker.com/signup*

- Go to https://labs.play-with-docker.com/.

- Click on start then login.

- Add a new instance.

- Pull the image with: `docker pull ericlarg4/shiny-ids-app:release`.

- Run the image with: `docker run -p 3838:3838 ericlarg4/shiny-ids-app:release`

- At the top of the screen, click on the 3838 button (next to open port).

## Use of the app

See the wiki.




