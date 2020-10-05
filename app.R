library(shiny)
library(shinythemes)
library(reshape2)
library(ggplot2)
library(spatstat)
library(dplyr)
library(car)
library(mgcv)
library(rgl)
library(plotly)

get.surface.model <- function(df, layer.name){
  layer.df <- df[df$couche == layer.name, c("x_rand", "y_rand", "z_rand") ]
  # calcul du GAM:
  fit <- gam(z_rand ~ s(x_rand, y_rand), data = layer.df)
  # ajout des valeurs d'altitude prédites:
  layer.df$pred <- predict(fit)
  # création d'un tableau pour le samplage:
  x <- seq(min(layer.df$x_rand), max(layer.df$x_rand), len = 100)
  y <- seq(min(layer.df$y_rand), max(layer.df$y_rand), len = 100)
  plot.df <- expand.grid(x_rand=x, y_rand=y)
  
  plot.df$predict <- predict(fit, newdata = plot.df) 
  z <- dcast(plot.df, x_rand~y_rand, value.var="predict")
  # préparation des valeurs d'output:
  x.values <- z[,1]
  y.values <- as.numeric(colnames(z[,-1]))
  z.matrix <- as.matrix(z[-1]) * -1
  col <- as.character(unique(df[df$couche == layer.name,]$couche.col))
  col <- gsub("[0-9]*", "", col)
  
  list(z.matrix = z.matrix, x = x.values, y = y.values, color = col)
}


ui <- shinyUI(fluidPage(
  theme = shinytheme("flatly"),  # slate  flatly
  
  fluidRow(
    column(2,
           h3("Virtual Poeymaü"),
           p(
             a("S. Plutniak", href="https://sebastien-plutniak.github.io",  target="_blank"),
             "/",
             a("PAVO", href="https://pavo.hypotheses.org/",  target="_blank"),
             "project"
           ),
           checkboxGroupInput("localisation", h4("Localisation method"),
                              choices = list("Point" = "point",
                                             "Volume" = "volume"),
                              selected = c("point", "volume")),
           
           checkboxGroupInput("objects", h4("Remains class"),
                              choices = list("Stone tools" = "lithique taillé",
                                             "Fauna" = "faune",
                                             "Lithic" = "lithique non silex",
                                             "Core" = "nucléus",
                                             "Bone tools" = "industrie osseuse",
                                             "Pottery" = "céramique"
                              ),
                              selected = NULL),
           h4("Plot options"),
           # checkboxInput("surface", "Surfaces", value = F)
           checkboxInput("burial", "Show burial", value = F),
           checkboxInput("surface", "Show surfaces", value = F),
           numericInput("point.size", "Point size", 3, min = 1, max=5, width="50%")
    ),
    # column(7,
    column(10,
           
           plotlyOutput("plot3d",  width = 900, height = 650),
           # rglwidgetOutput("plot3d",  width = 800, height = 700),
    )
    # column(3,
    #        br(), br(), br(), br(), br(), br(), br(),
    #        # imageOutput("plotLegend")
    # )
  ),
  fluidRow(
    column(5),
    column(4,
           h4("Remains class and localisation method"),
           tableOutput("classLocalStats")
    ),
    column(3,
           h4("Remains by layer"),
           tableOutput("layersStats")
    )
  )
))





# DEFINE SERVER  ----    
server <- function(input, output) { 
  
  df <- read.csv("data/20209999_PCR_PAVO_Poeymau_coordonnees_gen.csv",
                 header = T, stringsAsFactors = F)
  df$id <- 1:nrow(df)
  
  # suppression des lignes aux données manquantes
  df <- df[ ! is.na(df$couche),]
  df <- df[ ! is.na(df$carre_x),]
  
  
  
  # préparation des localisation horizontales ####
  
  # création d'un identifiant de carré
  df$square <- paste(df$carre_y, df$carre_x, sep = "")
  
  # attribution des min / max de la couche dans le carré ou l'objet se trouve ####
  obs.minmax <- group_by(df, couche, square) %>% summarize(
    x.min.obs = min(x_min, na.rm=T),
    x.max.obs = max(x_max, na.rm=T),
    y.min.obs = min(y_min, na.rm=T),
    y.max.obs = max(y_max, na.rm=T))
  # remplacement des valeurs manquantes par les min et max (0 et 100)
  obs.minmax[ abs(obs.minmax$x.min.obs) == Inf, ]$x.min.obs <- 0
  obs.minmax[ abs(obs.minmax$y.min.obs) == Inf, ]$y.min.obs <- 0
  obs.minmax[ abs(obs.minmax$x.max.obs) == Inf, ]$x.max.obs <- 100
  obs.minmax[ abs(obs.minmax$y.max.obs) == Inf, ]$y.max.obs <- 100
  # ajout des valeurs au tableau principal:
  df <- merge(df, obs.minmax, c("couche", "square"), all.x=T)
  # identifiant des objets sans x ou y min / max
  no.xmin <- is.na(df$x_min)
  no.ymin <- is.na(df$y_min)
  no.xmax <- is.na(df$x_min) & is.na(df$x_max)
  no.ymax <- is.na(df$y_min) & is.na(df$y_max)
  # affectation des valeurs min max observées:
  df[no.xmin,]$x_min <- df[no.xmin,]$x.min.obs
  df[no.xmax,]$x_max <- df[no.xmax,]$x.max.obs
  df[no.ymin,]$y_min <- df[no.ymin,]$y.min.obs
  df[no.ymax,]$y_max <- df[no.ymax,]$y.max.obs
  
  # conversion des coordonnées xy relative à chaque carré en des coordonnées 
  # générales au site. On ajoute une pondération à chaque carré en x et en y:
  df$x_correction <- factor(df$carre_x,
                            levels = 11:0,
                            labels = seq(0, 1100, 100) )
  df$x_correction <- as.numeric(as.character(df$x_correction))
  
  df$y_correction <- factor(df$carre_y,
                            levels = c("Y", "Z", LETTERS[1:9]),
                            labels = seq(0, 1000, 100) )
  df$y_correction <- as.numeric(as.character(df$y_correction))
  
  df$carre_y <- factor(df$carre_y, levels = c("Y", "Z", "A", "B", "C", "D", "E"))
  df$carre_x <- factor(df$carre_x, levels = c(4:0))
  
  # application des correctifs:
  df$x_min <- df$x_min + df$x_correction
  df$x_max <- df$x_max + df$x_correction
  df$y_min <- df$y_min + df$y_correction
  df$y_max <- df$y_max + df$y_correction
  
  
  # controle:
  summary(df[df$square=="B2", ]$x_min)
  summary(df[df$carre_y == "D",]$y_min)
  summary(df[df$carre_x == "2",]$x_min)
  
  # préparation des localisations verticales ####
  couches <- c("CS", "CT", "FSH", "CI", "CI ou FIH", "FIH", "foyer", "BS", "BS jaune", "CPE", "CPE ou BI", "CN", "BI-CN", "BI")
  df <- df[ df$couche %in% couches,]
  
  
  # réordonnancement des couches
  df$couche <- factor(df$couche, levels = couches, labels = couches)
  df$couche.col <- factor(df$couche,
                          levels = couches,
                          labels = c(CS="gold", CT="deepskyblue4", FSH="orangered2", 
                                     CI="darkgreen", 'CI ou FSH'="grey30", FIH="tan4", 
                                     foyer="grey31", BS="gold1", 'BS jaune'="gold3", 
                                     CPE="deepskyblue3", 'CPE ou BI' = "grey32", 
                                     CN = "grey10",  'BI-CN'="grey33", BI="orangered3")  )
  
  
  
  # attribution des min / max des altitudes dans le carré ou l'objet se trouve ####
  df[ ( ! is.na(df$z_min)) & is.na(df$z_max), ]$z_max <- df[ ( ! is.na(df$z_min)) & is.na(df$z_max), ]$z_min
  
  # identification des zmin/max par carré et par couche:
  obs.z.minmax <- group_by(df, couche, square) %>% summarize(
    z.min.obs = min(z_min, na.rm=T),
    z.max.obs = max(z_max, na.rm=T))
  
  # identification (complémentaire) des zmin/max par couche:
  obs.z.minmax.couche <- group_by(df, couche) %>% summarize(
    z.min.obs.couche = min(z_min, na.rm=T),
    z.max.obs.couche = max(z_max, na.rm=T))
  # fusion:
  obs.z.minmax <- merge(obs.z.minmax, obs.z.minmax.couche, by = "couche", all.x=T)
  # correction des valeurs infinies:
  sel <- is.infinite(obs.z.minmax$z.min.obs)
  obs.z.minmax[sel, ]$z.min.obs <- obs.z.minmax[sel, ]$z.min.obs.couche
  obs.z.minmax[sel, ]$z.max.obs <- obs.z.minmax[sel, ]$z.max.obs.couche
  
  # ajout des valeurs au tableau principal
  df <- merge(df, obs.z.minmax, c("couche", "square"), all.x=T)
  # identifiant des objets sans x ou y min / max
  no.zmin <- is.na(df$z_min)
  no.zmax <- is.na(df$z_min) & is.na(df$z_max)
  # affectation des valeurs min max observées:
  df[no.zmin,]$z_min <- df[no.zmin,]$z.min.obs
  df[no.zmax,]$z_max <- df[no.zmax,]$z.max.obs
  
  
  # Gestion des incertitudes ####
  # — attribution de coordonnées moyennes pour les objets localisés par volume
  # df$z_moy <- apply(df[, 10:11], 1, mean)
  # df$x_moy <- apply(df[, 12:13], 1, mean)
  # df$y_moy <- apply(df[, 14:15], 1, mean)
  
  # — attribution de coordonnées aléatoires pour les objets localisés par volume
  df$z_rand <- df$z_min
  df$x_rand <- df$x_min
  df$y_rand <- df$y_min
  
  df[which(df[, "z_min"] != df[, "z_max"]), ]$z_rand <- 
    apply(df[which(df[, "z_min"] != df[, "z_max"]), 11:12], 1,
          function(x) sample(x[1]:x[2], 1) )
  df[which(df[, "x_min"] != df[, "x_max"]), ]$x_rand <- 
    apply(df[which(df[, "x_min"] != df[, "x_max"]), 13:14], 1,
          function(x) sample(x[1]:x[2], 1) )
  df[which(df[, "y_min"] != df[, "y_max"]), ]$y_rand <- 
    apply(df[which(df[, "y_min"] != df[, "y_max"]), 15:16], 1,
          function(x) sample(x[1]:x[2], 1) )
  
  # — marquage des modes de localisation ####
  df$localisation_mode <- "point"
  df[which(df[, "z_min"] != df[, "z_max"]),]$localisation_mode <- "volume"
  df[which(df[, "x_min"] != df[, "x_max"]),]$localisation_mode <- "volume"
  df[which(df[, "y_min"] != df[, "y_max"]),]$localisation_mode <- "volume"
  
  # types d'objets ####
  df$objet_type <- tolower(df$objet_type)
  
  # subset
  df.sub <- df[, c("id", "square", "carre_x", "carre_y", "square", "z_rand", "couche", "couche.col", "x_rand", "y_rand",
                   "localisation_mode", "objet_type")]
  
  df.sub <- df.sub[df.sub$couche %in% c("BS", "CS", "CN", "CI ou FIH", "CPE", "FSH",  "CT", "BI", "CI", "FIH"),]
  df.sub <- droplevels(df.sub)
  
  output$plotOuput <- renderPlot({  # plot des diagrammes ----
    
    # type de localisation
    if(input$point) {selection <- "point"}
    if(input$volume) {selection <- "volume"}
    if(input$point & input$volume) {selection <- c("point", "volume")}
    
    df.sub <- df.sub[df.sub$localisation_mode %in% selection, ]
    
    # type d'objets
    if( ! is.null(input$objects)){
      df.sub <- df.sub[df.sub$objet_type %in% input$objects, ]
    }
    
    # ggplot(df.sub,
    #        aes(x = x_rand, y = - z_rand, group = couche, color = couche)) +
    #   theme_light(base_size = 8) +
    #   geom_vline(xintercept = seq(700, 1200, 100), color= "grey50") +
    #   geom_point(aes(shape = localisation_mode),
    #              size = .6,
    #              show.legend = c(color=T, fill=F, shape=T),
    #              fill="grey80" ) +
    #   geom_smooth(size=.5, fill = "grey20") +
    #   scale_x_continuous("Coupe frontale", breaks = seq(750, 1150, 100),
    #                      labels = c(4:0))  +
    #   scale_y_continuous("profondeur (cm)", breaks = - seq(0, 700, 100)) +
    #   scale_color_manual(values = levels(df.sub$couche.col) ) +
    #   scale_shape_manual("Localisation par :", values = c(21, 23)) +
    #   coord_fixed()
  })
  
  
  
  # output$plot3d <- renderRglwidget({ # plot du nuage 3D scatter ----
  output$plot3d <- renderPlotly({ # plot  3D plotly ----
    
    # type de localisation
    df.sub <- df.sub[df.sub$localisation_mode %in% input$localisation, ]
    
    # type d'objets
    if( ! is.null(input$objects)){
      df.sub <- df.sub[df.sub$objet_type %in% input$objects, ]
    }
    
    
    output$classLocalStats <- renderTable({
      stats.df <- table(df.sub$objet_type, df.sub$localisation_mode)
      
      if(nrow(stats.df) > 1 & ncol(stats.df) > 1){
        stats.df <- as.matrix(stats.df)
        stats.df <- stats.df[order(stats.df[,1], decreasing = T), ]
        stats.df <- rbind(stats.df, sum = apply(stats.df, 2, sum))
      }
      
      if(ncol(stats.df) == 1 & nrow(stats.df) > 1){
        stats.df <- stats.df[order(stats.df[,1], decreasing = T), ]
        stats.df <- c(stats.df, sum = sum(stats.df))
        stats.df <- as.data.frame(stats.df)
        colnames(stats.df) <- input$localisation
      } else {
        stats.df <- as.data.frame.matrix(stats.df)
      }
      stats.df
    }, rownames = T)
    
    output$layersStats <- renderTable({
      stats.df <- table(df.sub$couche)
      stats.df <- data.frame(stats.df)
      colnames(stats.df) <- c("Layer", "n")
      stats.df
    }, rownames = F)
    
    
    # output$plotLegend <- renderPlot({  # plot de la légende ----
    #   colors.df <- unique(df.sub[, c("couche", "couche.col")])
    #   colors.df <- colors.df[order(match(colors.df$couche, levels(colors.df$couche))),]
    #   colors.df$y <- nrow(colors.df):1
    # 
    #   ggplot(colors.df) +
    #     theme_void() +
    #     geom_point(aes(y = y, color = couche), x = 0, size= 5, show.legend = F)  +
    #     geom_text(aes(y = y, label = couche) , x = 0.05, size=5, hjust = "inward", show.legend = F) +
    #     scale_color_manual(values =  as.character(factor(colors.df$couche.col)))
    # })
    
    fig <- plot_ly(df.sub, x = ~x_rand * -1, y = ~y_rand * -1, z = ~z_rand * -1,
                   color = ~couche,
                   sizes = input$point.size,
                   colors = as.character(levels(df.sub$couche.col)),
                   text = ~paste('id:', id,
                                 '<br>Square:', square,
                                 '<br>Localisation:', localisation_mode,
                                 '<br>Class:', objet_type)
    )
    # ajout des points
    fig <- fig %>% add_markers(size = input$point.size)
    
    # localisation sépulture
    # summary(df[df$square=="B1" & df$couche =="CI",]$y_min)
    # xmin = 1000, 1100
    # ymin = 300, 400
    # zmin =190, 230
    if(input$burial){
      fig <- fig %>% add_mesh(x =  - c(1000, 1000, 1100, 1100, 1000, 1000, 1100, 1100),
                              y =   - c(300, 400, 400, 300, 300, 400, 400, 300), 
                              z = - c(190,   190,  190,  190,  230,  230,  230,  230), 
                              i = c(7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2), 
                              j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3), 
                              k = c(0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6), 
                              intensity = 1, color = 1, colors = rgb(.8, .8, .8), showscale=F,
                              text = "Square: B1<br>Approximate location of the burial")
    }
    
    # ajout des surfaces ####
    if(input$surface){
      # sélection des surfaces à calculer :
      couches <- table(df.sub$couche) 
      couches <- names(couches[couches > 100])
      
      # calcul des surfaces:
      surf.list <- lapply(couches, get.surface.model, df=df.sub)
      
      i <- 1
      while(i <= length(surf.list)){
        fig <- fig %>% add_surface(z=surf.list[[i]]$z.matrix,
                                   x = surf.list[[i]]$x * -1,
                                   y = surf.list[[i]]$y * -1, 
                                   inherit = F,
                                   colorscale = list(c(0, 1), c("black", surf.list[[i]]$col )),
                                   opacity = .7, showscale=F)
        i <- i + 1
      }
    }
    
    # paramètrage de la visualisation ####
    fig %>% layout(scene = list(
      # camera = list(eye = list(x=-1.25, y=-2, z=1.25)),
      # aspectmode = 'manual',
      xaxis = list(title = 'X',
                   range =  -c(1200, 700),
                   tickmode = "array",
                   tickvals = -seq(750, 1200, 50),
                   ticktext = c(rbind(levels(df.sub$carre_x), ""))
      ),
      yaxis = list(title = 'Y',
                   range = - c(700, 0),
                   tickmode = "array",
                   tickvals = -seq(50, 700, 50),
                   ticktext = c(rbind(levels(df.sub$carre_y), ""))
      ),
      zaxis = list(title = 'Depth (m)',
                   range = c(-800, 0),
                   tickmode = "array",
                   tickvals = - seq(100, 700, 100),
                   ticktext =  seq(1, 7, 1)
      )
    ))
    
    
    # rgl.open(useNULL=T)
    # scatter3d(x = df.sub$x_rand, y = - df.sub$z_rand, z =  df.sub$y_rand,
    #           # point.col = df.sub$couche.col,
    #           sphere.size=.1,
    #           groups = df.sub$couche,
    #           surface = input$surface,  fit = "smooth", parallel = F,
    #           surface.col= as.character(levels(df.sub$couche.col)),
    #           grid = F, residuals = F,
    #           # revolutions = 2, speed = 2,
    #           # ellipsoid = TRUE,
    #           axis.scales = T,
    #           axis.ticks = F,
    #           xlab = "X", ylab = "Depth (cm)",
    #           zlab = "Y")
    # 
    # rglwidget()
    
    
  })
  
  
  # fin du serveur
}

# Run app:
shinyApp(ui = ui, server = server)


