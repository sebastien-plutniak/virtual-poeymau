library(ggplot2)
library(car)
library(mgcv) # for gam()
library(rgl)
library(cxhull)

get.cxhull.model <- function(df, layer.name){
  layer.df <- df[df$layer == layer.name, c("xrand", "yrand", "zrand") ]
  layer.df <- unique(layer.df)
  hull.df <- cxhull::cxhull(as.matrix(layer.df), triangulate = T) # compute convex hull
  hull.df <- cxhull::hullMesh(hull.df)                            # extract mesh
  nfaces <- nrow(hull.df$faces)
  hull.df <- - hull.df$vertices              # convert to negative coordinates
  color <- as.character(unique(df[df$layer == layer.name,]$layer.col))
  color <- gsub("[0-9]*", "", color)
  list("hull"=hull.df, "color"=color, "nfaces"=nfaces)
}    

get.surface.model <- function(df, layer.name){
  layer.df <- df[df$layer == layer.name, c("xrand", "yrand", "zrand") ]
  # calcul du GAM:
  fit <- gam(zrand ~ s(xrand, yrand), data = layer.df)
  # ajout des valeurs d'altitude prédites :
  layer.df$pred <- predict(fit)
  # création d'un tableau pour le samplage:
  x <- seq(min(layer.df$xrand), max(layer.df$xrand), len = 100)
  y <- seq(min(layer.df$yrand), max(layer.df$yrand), len = 100)
  plot.df <- expand.grid(xrand=x, yrand=y)
  
  plot.df$predict <- predict(fit, newdata = plot.df) 
  z <- dcast(plot.df, xrand ~ yrand, value.var="predict")
  # préparation des valeurs d'output:
  x.values <- z[,1]
  y.values <- as.numeric(colnames(z[,-1]))
  z.matrix <- as.matrix(z[-1]) * -1
  col <- as.character(unique(df[df$layer == layer.name,]$layer.col))
  col <- gsub("[0-9]*", "", col)
  
  list(z.matrix = z.matrix, x = x.values, y = y.values, color = col)
}


# DEFINE SERVER  ----    
server <- function(input, output) {
  
  # CAVE MAP ####
  cave.grid <- expand.grid(square_x = 11:-2, 
                           square_y = c("Y", "Z", LETTERS[1:9]))
  
  cave.grid$square_x <- factor(cave.grid$square_x, levels = 11:-2)
  
  cave.contour <- data.frame(matrix(c(1.2,2, 2,2,  2.4,2, 1.5,4, 1.3,5, 1.1,6, 1,6.5,
                                      .8,7.5, 1.2,8, 2,8.2, 3,8.4, 4,8.4, 4.4,9, 5,9.5,
                                      6,10.8, 7,8.2,  8,7.2, 8.2,7, 9,5.9, 10,4, 11,2.6,
                                      11.5,2, 12,2),
                                    ncol=2, byrow = T))
  colnames(cave.contour) <- c("x", "y")
  cave.contour  <- cave.contour + .5  
  
  cave.map <- ggplot() +
    theme_minimal(base_size = 11) +
    geom_tile(data = cave.grid, aes(square_x, y=square_y), alpha=0) +
    geom_vline(xintercept =  seq(0.5, 13.5, 1), colour = "grey70" ) +
    geom_hline(yintercept =  seq(0.5, 12.5, 1), colour = "grey70" ) +
    geom_path(data=cave.contour, aes(x = x, y = y) ) +
    coord_fixed() + xlab("") + ylab("") 
  
  cave.sep.map <- cave.map +
    geom_rect(aes(xmin = 10.5, xmax = 13.5 , ymin= 3.5, ymax= 4.5),
              fill="grey", alpha=.8) 
  # end cave map
  
  # preprocessing Poeymau data: ####
  source("data-preprocessing.R")
  
  # selection specific to the Poeymau cave:
  dataset <- reactive({
    df.list <- list()
    if("1950s" %in% input$periods){
      period1950s.sub.df <- period1950s.df[, c("id", "square", "square_x",
                                               "square_y", "square", "zrand",
                                               "layer", "sublayer", "xrand", "yrand",
                                               "localisation_mode", "object_type")]
      period1950s.sub.df <- period1950s.sub.df[period1950s.sub.df$layer %in% c("CS", "CT", "FSH", "CI", "CI ou FIH", "FIH", "BS", "CPE", "CPE ou BI", "BI"),]
      df.list$"1950s" <- period1950s.sub.df
    }
    
    if("1970s" %in% input$periods){
      period1970s.sub.df <- period1970s.df[, c("id", "square", "square_x",
                                               "square_y", "square", "zrand",
                                               "layer", "sublayer", "xrand",
                                               "yrand", "localisation_mode", "object_type")]
      df.list$"1970s" <- period1970s.sub.df
    }
    
    df <- Reduce(rbind, df.list)
    
    # type de localisation:
    df <- df[df$localisation_mode %in% input$localisation, ]
    # type d'objets:
    if( ! is.null(input$objects)){
      df <- df[df$object_type %in% input$objects, ]
    }
    
    df$square_y <- factor(df$square_y,
                              levels = c("Y", "Z", "A", "B", "C", "D", "E", "F"))
    df$square_x <- factor(df$square_x, levels = c(7:-2))
    
    
    # réordonnancement des couches:
    layers <- c("alsh", "CS", "CT", "FSH", "CI", "CI ou FIH", "FIH",  "BS",
                "BS jaune", "CPE", "CPE ou BI", "BI")
    
    df$layer <- factor(df$layer, levels = layers, labels = layers)
    df$layer.col <- factor(df$layer,
                               levels = layers,
                               labels = c(alsh ="blue", CS="gold", CT="deepskyblue4", FSH="orangered2", 
                                          CI="darkgreen", 'CI ou FSH'="grey30", FIH="tan4", 
                                          BS="gold1", 'BS jaune'="gold3", 
                                          CPE="deepskyblue3",  
                                          'CPE ou BI' = "grey32", BI="orangered3")  )
    df
  }) # end of the poeymau-specific definition of the dataset
  
  
  # tableau des classes d'objets ---
  output$classLocalStats <- renderTable({
    dataset <- dataset()
    stats.df <- table(dataset$object_type, dataset$localisation_mode)
    
    if(nrow(stats.df) > 1 & ncol(stats.df) > 1){
      stats.df <- as.matrix(stats.df)
      stats.df <- stats.df[order(stats.df[,1], decreasing = T), ]
      stats.df <- rbind(stats.df, total = apply(stats.df, 2, sum))
      stats.df <- cbind(stats.df, total = apply(stats.df, 1, sum))
    }
    
    if(ncol(stats.df) == 1 & nrow(stats.df) > 1){
      stats.df <- stats.df[order(stats.df[,1], decreasing = T), ]
      stats.df <- c(stats.df, total = sum(stats.df))
      stats.df <- as.data.frame(stats.df)
      colnames(stats.df) <- input$localisation
    } else {
      stats.df <- as.data.frame.matrix(stats.df)
    }
    stats.df
  }, rownames = T, digits=0)
  
  # tableau des couches ----
  output$layersStats <- renderTable({
    dataset <- dataset()
    stats.df <- group_by(dataset, layer, localisation_mode) %>%
      summarise(n = n())
    stats.df <- dcast(stats.df, layer~localisation_mode, value.var="n")
    stats.df[is.na(stats.df)] <- 0
    rownames(stats.df) <- stats.df[, 1]
    
    if(ncol(stats.df) == 3){
      stats.df <- stats.df[, -1]
      stats.df$total <- apply(stats.df, 1, sum)
      colnames(stats.df) <- c("point", "volume", "total")
    }
    if(ncol(stats.df) == 2){
      labels <- stats.df$layer
      loc.type <- colnames(stats.df)[2]
      stats.df <- stats.df[,-1]  
      names(stats.df) <- labels
      stats.df <- as.data.frame(stats.df)
      colnames(stats.df) <- loc.type
      
    }
    stats.df <- rbind(stats.df, "total" = apply(stats.df, 2, sum))
    stats.df
  }, rownames = T, digits=0)
  
  
  # tableau de l'id sélectionné ---
  output$id.tab <- renderTable({
    df.tab <- period1950s.df[period1950s.df$id == input$id, c("id", "square", "layer", "zmin", "zmax", "object_text", "object_alteration", "object_type", "object_material")]
    colnames(df.tab) <- c("id", "square", "layer", "zmin", "zmax", "description", "alteration", "class", "material")
    df.tab
  }, digits=0)
  
  output$id.table <- renderUI({
    div(style = 'overflow-x: scroll; overflow: auto', 
        tableOutput('id.tab'))
  })
  
  
  # plot  3D plotly ----
  output$plot3d <- renderPlotly({ 
    dataset <- dataset()
    dataset$point.size <- input$point.size
    size.scale <- input$point.size
    
    # highlight hearths ----
    if(input$hearth){
      levels(dataset$layer) <- c(levels(dataset$layer), "hearth")
      selection <- dataset$sublayer == "foyer" & dataset$object_type == ""
      dataset[selection ,]$layer <- "hearth"
      levels(dataset$layer.col) <- c(levels(dataset$layer.col), "black")
    }
    
    # Plot initial ----
    fig <- plot_ly(dataset, x = ~xrand * -1, y = ~yrand * -1, z = ~zrand * -1,
                   color = ~layer,
                   colors = as.character(levels(dataset$layer.col)),
                   size  = ~point.size,
                   sizes = size.scale,
                   axis="x2",
                   marker = list(symbol = 'square', sizemode = 'diameter'),
                   text = ~paste('id:', id,
                                 '<br>Square:', square,
                                 '<br>Localisation:', localisation_mode,
                                 '<br>Class:', object_type)
    ) %>%   config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "poeymau3D",
        width = 600, height = 600
      ))
    
    #  — add points ####
    fig <- fig %>% add_markers()
    
    # — add burial ####
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

    
    # — add C14 ----
    if(input$c14){
      c14.id <- c(14579, 14599, 10741, 2032, 7967)
      c14.df <- dataset[dataset$id %in% c14.id, ]
      c14.df[, c("xrand", "yrand", "zrand")] <- - c14.df[, c("xrand", "yrand", "zrand")]
      fig <- fig %>% add_trace(x=~xrand, y=~yrand, z=~zrand, data=c14.df,
                                 name = 'C14', mode = 'markers',
                                 sizes=5,
                                 color = I("black"), inherit = F)
    }
    
    # — add surfaces ####
    if(input$surface){
      # sélection des surfaces à calculer :
      layers <- table(dataset$layer) 
      layers <- names(layers[layers > 100])
      
      # compute regression surfaces:
      surf.list <- lapply(layers, get.surface.model, df=dataset)
      # add traces:
      for(i in  1:length(surf.list)){
        fig <- fig %>% add_surface(z = surf.list[[i]]$z.matrix,
                                   x = surf.list[[i]]$x * -1,
                                   y = surf.list[[i]]$y * -1, 
                                   inherit = F,
                                   colorscale = list(c(0, 1), c("black", surf.list[[i]]$col )),
                                   hoverinfo="skip",
                                   opacity = .7, showscale=F)
      }
    }
    
    # — add convex hull ####
    if(input$cxhull){
      # sélection des enveloppes à calculer :
      layers <- table(dataset$layer) 
      layers <- names(layers[layers > 100])
      
      # calcul des enveloppes:
      mesh.list <- lapply(layers, get.cxhull.model, df=dataset)
      
      for(i in 1:length(mesh.list)){
        fig <- fig %>% add_mesh(x = mesh.list[[i]][[1]][,1],
                                y = mesh.list[[i]][[1]][,2],
                                z = mesh.list[[i]][[1]][,3],
                                facecolor = rep(mesh.list[[i]]$color, mesh.list[[i]]$nfaces),
                                hoverinfo="skip",  showscale=F,
                                opacity = 0.4, alphahull=0 , inherit = F)
      }
    }    
    
    # — add grid and bottom map ----
    coordx <- data.frame(
      id = c(0,0, 1,1, 2,2, 3,3, 4,4, 5,5, 6,6, 7,7, 8,8, 9,9, 10,10,
             11,11, 12,12, 13,13, 14,14),
      x = - c(0,0, 100,100, 200,200, 300,300, 400,400, 500,500, 600,600,
              700,700, 800,800, 900,900, 1000,1000, 1100,1100, 1200,1200,
              1300,1300, 1400,1400),
      y = - rep(c(0, 1200), 15),
      z =  -800)
    coordy <- data.frame(
      id = c(0,0, 1,1, 2,2, 3,3, 4,4, 5,5, 6,6, 7,7, 8,8, 9,9, 10,10,
             11,11, 12,12),
      y = - c(0,0, 100,100, 200,200, 300,300, 400,400, 500,500, 600,600,
              700,700, 800,800, 900,900, 1000,1000, 1100,1100, 1200,1200),
      x = - rep(c(0, 1400), 13),
      z =  -800)
    
    fig <- fig  %>% add_paths(x = ~x,   y = ~y, z = ~z, data=coordx,
                       split = ~id,
                       color = I("grey50"), showlegend=F,
                       hoverinfo="skip",
                       inherit = F)   %>% 
                    add_paths(x = ~x,   y = ~y, z = ~z, data=coordy,
                      split = ~id,
                      color = I("grey50"), showlegend=F,
                      hoverinfo="skip",
                      inherit = F)  %>% 
                    # ajout plan
                    add_paths(x = -cave.contour$x * 100,
                              y = -cave.contour$y * 100,
                              z = rep(-795, nrow(cave.contour)),
                              color = I("black"), showlegend=F,
                              hoverinfo="skip",
                              inherit = F)
    
    # — paramètrage de la visualisation et affichage ####
    fig %>% layout(scene = list(
      # camera = list(eye = list(x=-1.25, y=-2, z=1.25)),
      xaxis = list(title = 'X',
                   tickmode = "array",
                   range =  -c(1400, 100),
                   tickvals = -seq(150, 1400, 50),
                   ticktext = c(rbind(10:-2, " ")),
                   zeroline=F, showline=F
      ),
      yaxis = list(title = 'Y',
                   tickmode = "array",
                   range = - c(1200, 0),
                   tickvals = -seq(50, 1100, 50),
                   ticktext = c(rbind(c("Y", "Z", LETTERS[1:9]), " "))
      ),
      zaxis = list(title = 'Depth (m)',
                   range =  c(-800, 0),
                   tickmode = "array",
                   tickvals = - seq(0, 700, 100),
                   ticktext =   0:7
      ),
      aspectmode = "manual", 
      aspectratio = list(x = 1, y = 1, z = input$ratio *.33)
    ))  #end layout
  }) # end plot3d
  
  
  # Section X ####
  min.max.X <- eventReactive(input$goButtonX, {
    list(coordx = seq(input$sectionXx[1], input$sectionXx[2]),
         coordy = seq(input$sectionXy[1], input$sectionXy[2]))
  })
  
  
  output$sectionXplot <- renderPlotly({
    dataset <- dataset()
    min.max.X <- min.max.X()
    
    sel <- dataset$xrand %in% min.max.X$coordx & dataset$yrand %in% min.max.X$coordy
    sectionX.df <- dataset[sel,  ]
    
    sectionX.df$coordx <- sectionX.df$yrand
    
    sectionX <- plot_ly(sectionX.df, x = ~coordx * -1, y = ~zrand * -1,
                        color = ~layer,
                        colors = as.character(levels(sectionX.df$layer.col)),
                        size  = 1,
                        sizes = c(1,5),
                        marker = list(symbol = 'square', sizemode = 'diameter'),
                        text = ~paste('id:', id,
                                      '<br>Square:', square,
                                      '<br>Localisation:', localisation_mode,
                                      '<br>Class:', object_type)
    )   %>%  
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "PoeymauSectionX",
          width = 600, height = 600
        ))  %>%
      add_markers() %>% 
      layout(xaxis = list(title="Y", 
                          zeroline = FALSE, 
                          range=c(0, -750),
                          tickvals = -seq(50, 750, 50),
                          ticktext = c(rbind(levels(dataset$square_y), ""), "")
      ),
      yaxis = list(title="Depth (m)",
                   zeroline = FALSE,
                   tickvals = - seq(0, 800, 100),
                   scaleanchor="x",
                   ticktext =   0:8,
                   range=c(-800,0)
      )
      )
    
    if(sum(sectionX.df$xrand %in% 1000:1300) > 0){ # ajout sépulture
      sectionX <-
        sectionX %>% 
        add_polygons(x = c(-400, -300, -300, -400, -400),
                     y = c(-190, -190, -230, -230, -190),
                     mode="line",
                     line = list(width=.1, color="grey"),
                     fillcolor = 'burlywood4', opacity = 0.9,
                     marker = list(opacity=0, size=.01),
                     text = paste("Square: B1<br>",
                                  "Localisation: volume<br>",
                                  "Class: burial", sep=""),
                     name = 'Burial', inherit=F)
    }
    
    sectionX
  })# end sectionX
  
  
  
  # Section Y ####
  
  min.max.Y <- eventReactive(input$goButtonY, {
    list(coordx = seq(input$sectionYx[1], input$sectionYx[2]),
         coordy = seq(input$sectionYy[1], input$sectionYy[2]))
  })
  
  output$sectionYplot <- renderPlotly({
    dataset <- dataset()
    min.max.Y <- min.max.Y()
    
    sel <- dataset$yrand %in% min.max.Y$coordy & dataset$xrand %in% min.max.Y$coordx
    
    sectionY.df <- dataset[sel,  ]
    sectionY.df$coordx <- sectionY.df$xrand
    
    sectionY <- plot_ly(sectionY.df, x = ~coordx, y = ~zrand * -1,
                        color = ~layer,
                        colors = as.character(levels(sectionY.df$layer.col)),
                        size  = 1,
                        sizes = c(1,5),
                        marker = list(symbol = 'square', sizemode = 'diameter'),
                        text = ~paste('id:', id,
                                      '<br>Square:', square,
                                      '<br>Localisation:', localisation_mode,
                                      '<br>Class:', object_type)
    )   %>%  
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "PoeymauSectionY",
          width = 600, height = 600
        )) %>%
      add_markers() %>% 
      layout(
        xaxis = list(title="X", 
                     zeroline = FALSE,
                     tickvals = seq(450, 1400, 50),
                     ticktext = c(rbind(levels(dataset$square_x), " ")),
                     range =  c(400, 1400)
        ),
        yaxis = list(title="Depth (m)",
                     zeroline = FALSE,
                     tickvals = - seq(0, 800, 100),
                     scaleanchor="x",
                     ticktext =   0:8,
                     range=c(-800,0)
        )
      )
    
    if(sum(sectionY.df$yrand %in% 300:400) > 0){ # ajout sépulture
      sectionY <-
        sectionY %>%
        add_polygons(x = c(1300, 1000, 1000, 1300, 1300),
                     y = c(-190, -190, -230, -230, -190),
                     mode="line",
                     line = list(width=.1, color="grey"),
                     fillcolor = 'burlywood4', opacity = 0.9,
                     marker = list(opacity=0, size=.01),
                     text = paste("Square: B1<br>",
                                  "Localisation: volume<br>",
                                  "Class: burial", sep=""),
                     name = 'Burial', inherit=F)
    }
    
    sectionY
  }) #end section Y
  
  
  
  # Plan Z ####
  
  planZ <- eventReactive(input$goButtonZ, {
    dataset <- dataset()
    min.max.Z <- seq(input$planZ[1], input$planZ[2])
    sel <- dataset$zrand %in% min.max.Z
    planZ.df <- dataset[sel,  ]
  })
  
  output$planZplot <- renderPlotly({
    dataset <- dataset()
    planZ.df <- planZ()
    
    planZ <- plot_ly(planZ.df, x = ~xrand, y = ~yrand * -1 ,
                     color = ~layer,
                     colors = as.character(levels(planZ.df$layer.col)),
                     size  = 1,
                     sizes = c(1,5),
                     marker = list(symbol = 'square', sizemode = 'diameter'),
                     text = ~paste('id:', id,
                                   '<br>Square:', square,
                                   '<br>Localisation:', localisation_mode,
                                   '<br>Class:', object_type)
    ) %>%
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "poeymauPlan",
          width = 600, height = 600
        )) %>%
      add_markers() %>%
      layout(
        xaxis = list(title="X",
                     zeroline = FALSE,
                     tickvals = seq(450, 1400, 50),
                     ticktext = c(rbind(levels(dataset$square_x), " ")),
                     range =  c(400, 1400)
        ),
        yaxis = list(title="Y",
                     zeroline = FALSE,
                     range = c(0, -750),
                     tickvals = -seq(50, 750, 50),
                     ticktext = c(rbind(levels(dataset$square_y), " "), " "),
                     scaleanchor="x"
        )
      )
    
    if(sum(planZ.df$zrand %in% 300:400) > 0){ # ajout sépulture
      planZ <-
        planZ %>%
        add_polygons(x = c(1300, 1000, 1000, 1300, 1300),
                     y = - c(300, 300, 400, 400, 300),
                     mode="line",
                     line = list(width=.1, color="grey"),
                     fillcolor = 'burlywood4', opacity = 0.9,
                     marker = list(opacity=0, size=.01),
                     text = paste("Square: B1<br>",
                                  "Localisation: volume<br>",
                                  "Class: burial", sep=""),
                     name = 'Burial', inherit=F)
    }
    planZ
    
  })
  
  
  # Density Plan
  density.plan <- eventReactive(input$goButtonZ, ({
    dataset <- dataset()
    planZ.df <- planZ()
    
    density.map <- cave.map 
    
    if(sum(planZ.df$zrand %in% 300:400) > 0){
      density.map <- cave.sep.map  
    }
    
    dens.map <- density.map +
      geom_point(data=planZ.df, aes(x = xrand /  100 + .5,
                                    y = yrand /  100 + .5,
                                    color = layer),
                 size = .1) +
      geom_density2d(data=planZ.df, aes(x = xrand  /  100 + .5,
                                        y = yrand  /  100 + .5),
                     size = .4, color = "red" ) +
      scale_color_viridis_d(option ="inferno", begin=.1, end = .9, alpha=.2)+
      guides(color = guide_legend(override.aes = list(size = 2, alpha=1) ) )  
  }))
  
  output$density.plan <- renderPlot({ density.plan() })
  
  output$download.density.plan <- downloadHandler(
    filename = "density-plan.svg",
    content = function(file) {
      ggsave(file, plot = density.plan())
    }
  )
  
  # cave map ####
  
  output$cave.mapX <- renderPlot({
    cave.sep.map +
      geom_rect(aes(
        ymin= input$sectionXy[1] / 100 + .5,
        ymax= input$sectionXy[2] / 100 + .5,
        xmin= input$sectionXx[1] / 100 + .5,
        xmax= input$sectionXx[2] / 100 + .5),
        fill="red", alpha=.7
      )
  })
  
  output$cave.mapY <- renderPlot({
    cave.sep.map +
      geom_rect(aes( 
        ymin= input$sectionYy[1] / 100 + .5,
        ymax= input$sectionYy[2] / 100 + .5,
        xmin= input$sectionYx[1] / 100 + .5,
        xmax= input$sectionYx[2] / 100 + .5),
        fill="red", alpha=.7)
  })
  
  
 #  excavation history ####
  
  output$cave.map.history <- renderPlot({
    hist.sub <- hist.df[hist.df$year == input$history.date, ]
    
    cave.map +
      geom_tile(data = hist.sub, aes(square_x, y=square_y, fill = state),
                show.legend = F) +
      scale_fill_manual("State:", values = c(rgb(.4, .7, 0, .5), rgb(1,1,1,0)) ) 
  })
  
  
  output$cave.map.history.grid <- renderPlot(
    cave.map +
      geom_tile(data = hist.df, aes(square_x, y=square_y, fill = state),
                show.legend = F)  +
      scale_fill_manual("State:", values = c(rgb(.4, .7, 0, .5), rgb(1,1,1,0)) ) +
      facet_wrap(~year) +
      theme(axis.text.x = element_text(size=.1),
            axis.text.y = element_text(size=.1))
  )
  
  
} # end of server.R
