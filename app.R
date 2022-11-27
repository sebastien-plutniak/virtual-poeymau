library(shiny)

# Run app:
source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)
