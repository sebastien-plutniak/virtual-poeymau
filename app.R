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
           checkboxInput("c14", "Show C14 samples", value = F),
           checkboxInput("hearth", "Highlight hearths", value = F),
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
    column(2),
    column(10, 
      h4("Object details (type an id)")
           )
  ),
  fluidRow(
    column(2),
    column(2,
      numericInput("id", "id", value = 1, min = 1, max=16000, width="50%")
      ),
    column(8,
      uiOutput("id.table")
    )
  ),
  br(),
  
  
  fluidRow(
    column(2),
    column(4,
           h4("Remains class and localisation method"),
           tableOutput("classLocalStats")
    ),
    column(4,
           h4("Remains by layer"),
           tableOutput("layersStats")
    ),
    column(2)
  )
))





# DEFINE SERVER  ----    
server <- function(input, output) { 
  
  df <- read.csv("data/20209999_PCR_PAVO_Poeymau_coordonnees_gen.csv",
                 header = T, stringsAsFactors = F)
  
  # suppression des lignes aux données manquantes
  df <- df[ ! is.na(df$couche),]
  df <- df[ ! is.na(df$carre_x),]
  
  # suppression des localisations douteuses
  df <- df[ ! df$localisation_douteuse == "oui",]
  
  # préparation des localisation horizontales ####
  
  # création d'un identifiant de carré
  df$square <- paste(df$carre_y, df$carre_x, sep = "")
  
  # attribution des x et y min / max de la couche dans le carré ou l'objet se trouve ####
  obs.minmax <- group_by(df, couche, square) %>% summarize(
    x.min.obs = min(x_min, na.rm=T),
    x.max.obs = max(x_max, na.rm=T),
    y.min.obs = min(y_min, na.rm=T),
    y.max.obs = max(y_max, na.rm=T))
  # remplacement des valeurs infinies par les min et max (0 et 100)
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
  
  # attribution coordonnées xy pour les objets localisés par sous-carrés ----
  # création tableau de référence des coordonnées des sous-carrés
  subsquare.ref <- data.frame(sous.carre = 1:9, 
                              x_min = c(0, 33, 67),
                              x_max = c(33, 66, 100),
                              y_min = c(0,0,0, 33,33,33, 67,67,67),
                              y_max = c(33,33,33, 66,66,66, 100,100,100))
  # fusion tableau et tableau de références des sous-carrés:
  df <- merge(df, subsquare.ref, by="sous.carre", all.x=T, suffixes = c("", "B"))
  # remplacement des coordonnées 
  df[df$sous.carre %in% 1:9, c("x_min", "x_max", "y_min", "y_max") ] <- 
    df[df$sous.carre %in% 1:9, c("x_minB", "x_maxB", "y_minB", "y_maxB") ]
  
  
  # conversion des coordonnées xy relative à chaque carré en des coordonnées 
  # générales au site. On ajoute une pondération à chaque carré en x et en y:
  df$x_correction <- factor(df$carre_x,
                            levels = 11:-2,
                            labels = seq(0, 1300, 100) )
  df$x_correction <- as.numeric(as.character(df$x_correction))
  
  df$y_correction <- factor(df$carre_y,
                            levels = c("Y", "Z", LETTERS[1:9]),
                            labels = seq(0, 1000, 100) )
  df$y_correction <- as.numeric(as.character(df$y_correction))
  
  df$carre_y <- factor(df$carre_y, levels = c("Y", "Z", "A", "B", "C", "D", "E"))
  df$carre_x <- factor(df$carre_x, levels = c(4:-2))
  
  # application des correctifs:
  df$x_min <- df$x_min + df$x_correction
  df$x_max <- df$x_max + df$x_correction
  df$y_min <- df$y_min + df$y_correction
  df$y_max <- df$y_max + df$y_correction
  
  
  # controle:
  summary(df[df$square == "B2", ]$x_min)
  summary(df[df$carre_y == "D",]$y_min)
  summary(df[df$carre_x == "2",]$x_min)
  
  # préparation des localisations verticales ####
  couches <- c("CS", "CT", "FSH", "CI", "CI ou FIH", "FIH",  "BS",
               "BS jaune", "CPE", "CPE ou BI", "BI")
  df <- df[ df$couche %in% couches,]
  
  # réordonnancement des couches
  df$couche <- factor(df$couche, levels = couches, labels = couches)
  df$couche.col <- factor(df$couche,
                levels = couches,
                labels = c(CS="gold", CT="deepskyblue4", FSH="orangered2", 
                           CI="darkgreen", 'CI ou FSH'="grey30", FIH="tan4", 
                           BS="gold1", 'BS jaune'="gold3", 
                           CPE="deepskyblue3",  
                           'CPE ou BI' = "grey32", BI="orangered3")  )
  
  
  
  # attribution des min / max des altitudes dans le carré ou l'objet se trouve ####
  selection <- ( ! is.na(df$z_min)) & is.na(df$z_max)
  df[selection, ]$z_max <- df[selection, ]$z_min
  
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
  
  
  # — attribution de coordonnées aléatoires pour les objets localisés par volume----
  df$z_rand <- df$z_min
  df$x_rand <- df$x_min
  df$y_rand <- df$y_min
  
  df[which(df[, "z_min"] != df[, "z_max"]), ]$z_rand <- 
    apply(df[which(df[, "z_min"] != df[, "z_max"]), c("z_min", "z_max") ], 1,
          function(x) sample(x[1]:x[2], 1) )
  df[which(df[, "x_min"] != df[, "x_max"]), ]$x_rand <- 
    apply(df[which(df[, "x_min"] != df[, "x_max"]), c("x_min", "x_max")], 1,
          function(x) sample(x[1]:x[2], 1) )
  df[which(df[, "y_min"] != df[, "y_max"]), ]$y_rand <- 
    apply(df[which(df[, "y_min"] != df[, "y_max"]), c("y_min", "y_max")], 1,
          function(x) sample(x[1]:x[2], 1) )
  
  # — marquage des modes de localisation ####
  df$localisation_mode <- "point"
  df[which(df[, "z_min"] != df[, "z_max"]),]$localisation_mode <- "volume"
  df[which(df[, "x_min"] != df[, "x_max"]),]$localisation_mode <- "volume"
  df[which(df[, "y_min"] != df[, "y_max"]),]$localisation_mode <- "volume"
  
  # types d'objets ####
  df$objet_type <- tolower(df$objet_type)
  
  # subset
  df.sub <- df[, c("id", "square", "carre_x", "carre_y", "square", "z_rand",
                   "couche", "couche.col",  "sous.couche", "x_rand", "y_rand",
                   "localisation_mode", "objet_type")]
  
  # actuellement toutes les couches:
  df.sub <- df.sub[df.sub$couche %in% c("CS", "CT", "FSH", "CI", "CI ou FIH", "FIH", "BS", "CPE", "CPE ou BI", "BI"),]
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

    # tableau des classes d'objets  #### 
      output$classLocalStats <- renderTable({
        
      stats.df <- table(df.sub$objet_type, df.sub$localisation_mode)

      if(nrow(stats.df) > 1 & ncol(stats.df) > 1){
        stats.df <- as.matrix(stats.df)
        stats.df <- stats.df[order(stats.df[,1], decreasing = T), ]
        stats.df <- rbind(stats.df, Total = apply(stats.df, 2, sum))
      }

      if(ncol(stats.df) == 1 & nrow(stats.df) > 1){
        stats.df <- stats.df[order(stats.df[,1], decreasing = T), ]
        stats.df <- c(stats.df, Total = sum(stats.df))
        stats.df <- as.data.frame(stats.df)
        colnames(stats.df) <- input$localisation
      } else {
        stats.df <- as.data.frame.matrix(stats.df)
      }
      stats.df
    }, rownames = T, digits=0)
    
    # tableau des couches  #### 
    output$layersStats <- renderTable({
      
        stats.df <- group_by(df.sub, couche, localisation_mode) %>% summarise(n = n())
        stats.df <- dcast(stats.df, couche~localisation_mode, value.var="n")
        stats.df[is.na(stats.df)] <- 0
        rownames(stats.df) <- stats.df[, 1]
        stats.df <- stats.df[, -1]
        stats.df$Total <- apply(stats.df, 1, sum)
        colnames(stats.df) <- c("Point", "Volume", "Total")
        stats.df <- rbind(stats.df,
                          "Total" = apply(stats.df, 2, sum))
        stats.df
      }, rownames = T, digits=0)

    
    # tableau de l'id sélectionné  #### 
    output$id.tab <- renderTable({
      df.tab <- df[df$id == input$id, c("id", "square", "couche", "z_min", "z_max", "objet_texte", "objet_alteration", "objet_type", "objet_matiere")]
      colnames(df.tab) <- c("id", "square", "layer", "z min", "z max", "description", "alteration", "class", "material")
      df.tab
    }, digits=0)
    
    output$id.table <- renderUI({div(style = 'overflow-x: scroll; overflow: auto', 
                  tableOutput('id.tab'))})
    
    
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
   
    # display C14 (temporary approach 20201023)  ----
    df.sub$point.size <- input$point.size
    size.scale <- input$point.size
    if(input$c14){
      tmp.id <- c(14579, 14599, 10741, 2032, 7967)
      df.sub[df.sub$id %in% tmp.id, ]$point.size <- 5
      size.scale <- c(input$point.size, input$point.size * 10)
    }
    # end temporary approach
    
    # highlight hearths ----
    if(input$hearth){
      levels(df.sub$couche) <- c(levels(df.sub$couche), "hearth")
      selection <- df.sub$sous.couche == "foyer" & df.sub$objet_type == ""
      df.sub[selection ,]$couche <- "hearth"
      levels(df.sub$couche.col) <- c(levels(df.sub$couche.col), "black")
    }
    
    fig <- plot_ly(df.sub, x = ~x_rand * -1, y = ~y_rand * -1, z = ~z_rand * -1,
                   color = ~couche,
                   colors = as.character(levels(df.sub$couche.col)),
                   size  = ~point.size,
                   sizes = size.scale,
                   marker = list(symbol = 'square', sizemode = 'diameter'),
                   text = ~paste('id:', id,
                                 '<br>Square:', square,
                                 '<br>Localisation:', localisation_mode,
                                 '<br>Class:', objet_type)
    )
    # ajout des points
    # fig <- fig %>% add_markers(size = input$point.size)
    fig <- fig %>% add_markers()

    # localisation sépulture
      # summary(df[df$square=="B1" & df$couche =="CI",]$y_min)
      # xmin = 1000, 1100
      # ymin = 300, 400
      # zmin =190, 230
    if(input$burial){
      fig <- fig %>% add_mesh(x =  - c(1000, 1000, 1300, 1300, 1000, 1000, 1300, 1300),
                       y =   - c(300, 400, 400, 300, 300, 400, 400, 300), 
                       z = - c(190,   190,  190,  190,  230,  230,  230,  230), 
                       i = c(7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2), 
                       j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3), 
                       k = c(0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6), 
                       intensity = 1, color = 1, colors = rgb(.8, .8, .8), showscale=F,
                       text = paste("Square: B1<br>",
                               "Localisation: volume<br>",
                               "Class: burial", sep=""))
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
      xaxis = list(title = 'X',
                   range =  -c(1400, 700),
                   tickmode = "array",
                   tickvals = -seq(750, 1400, 50),
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
                   tickvals = - seq(0, 750, 50),
                   ticktext =  c("", "", c(rbind(1:7, "")))
      ),
      aspectmode = "manual", 
      aspectratio = list(x = 1, y = 1, z = 800/700)
    ))

    
  })
    
  # fin du serveur
}

# Run app:
shinyApp(ui = ui, server = server)


