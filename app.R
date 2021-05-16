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
  layer.df <- df[df$layer == layer.name, c("xrand", "yrand", "zrand") ]
  # calcul du GAM:
  fit <- gam(zrand ~ s(xrand, yrand), data = layer.df)
  # ajout des valeurs d'altitude prédites:
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


ui <- shinyUI(
  fluidPage(
  theme = shinytheme("cosmo"),  # lumen  flatly simplex
  
  sidebarLayout(
    sidebarPanel(
           h4("Virtual Poeymaü"),
           p(
             "Dev. by: ",
             a("S. Plutniak", href="https://sebastien-plutniak.github.io",  target="_blank"),
             br(),
             a("PAVO", href="https://pavo.hypotheses.org/",  target="_blank"),
             "project",
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
           checkboxGroupInput("periods", h4("Periods"),
                              choices = list("1950s" = "1950s",
                                             "1970-80s" = "1970s"),
                              selected = c("1950s")),
    width=2), # end sidebarpanel
    
    mainPanel(
      tabsetPanel(
        tabPanel("Introduction",
          column(12, align="center",
           tags$div(
        HTML("<div style=width:450px;, align=left>
              <br><br>
   <p><b>Welcome to the <i>Virtual Poeymaü</i> application</b></p>
   <p>
     The <a target=blank href=https://fr.wikipedia.org/wiki/Grotte_du_Poeyma%C3%BC>Poeymaü</a>,
     a cave located in the western Pyrenees
     (in <a target=blank, href=https://www.openstreetmap.org/relation/2120173>Arudy</a>, France),
     was excavated by Georges Laplace and his collaborators from the early 1950s to the late 1980s.
     </p>
   <p>
     The site is currently restudied in the framework of the
     <a target=_blank, href=https://pavo.hypotheses.org><i>PAVO</i></a> research project
     (<i>Préhistoire ancienne de la vallée d'Ossau</i>). 
     Archaeological information was extracted from Laplace's excavation field notes related to
     the 1950s campaigns (archives of the <i><a href=https://musee-prehistoire-eyzies.fr>Musée national de Préhistoire</a></i>). 
     About 15,000 objects are documented and located in the space of the cave.
   </p>
   <p>
      The <i>Virtual Poeymaü</i> is a web application to explore this data set. 
      It allows to 
      filter the data, 
      plot them with interactive 3D visualisation, 
      generate and download longitudinal and transverse sections, and 
      generate summary tables.
   </p>
    <hr>
    For more information:
    <ul>
      <li> <b>R code:</b> <a target=blank, href=https://github.com/sebastien-plutniak/virtual-poeymau>github repository</a>.</li>
      <li> <b>reference:</b> Plutniak, S. 2021, “Virtual Poeymaü: a web application to explore the archaeological data from the excavation archives of the Poeymaü cave (France)”,  <i>Zenodo</i>. DOI: <a target=blank, href=https://doi.org/10.5281/zenodo.4443132>10.5281/zenodo.4443132</a>.</li>
      <li> <b>documentation:</b> 
        <ul>
          <li>Plutniak, S. 2021, “<a target=blank, href=https://pavo.hypotheses.org/98>Virtual Poeymaü : visualisation 3D interactive des données archéologiques extraites des archives de fouilles de la grotte du Poeymaü</a> », <i>Préhistoire ancienne de la vallée d'Ossau</i>.</li>
          <li>Plutniak, S. 2020, “Archives des fouilles de la grotte du Poeymaü (1951-1956) : informatisation et spatialisation des restes archéologiques”. Pétillon, J.-M. et B. Marquebielle (dir.), <i>Préhistoire ancienne de la vallée d'Ossau. Paléoenvironnement et sociétés de chasseurs-collecteurs dans le piémont pyrénéen. Projet collectif de recherche. Bilan 2020</i>, DRAC-SRA Nouvelle-Aquitaine, <a href=https://hal.archives-ouvertes.fr/hal-03092989>hal-03092989</a>.</li>
        </ul>

      </li>
    </ul>
    <hr>  
  
              </div>" )
            ) # end div()
                 ) # end column
        ),      #end tabPanel
        tabPanel("3D plot",    
                 fluidRow(
                   column(10,
                          plotlyOutput("plot3d",  width = 900, height = 650),
                          numericInput("id", "id", 1, min=1, max=16000, width="8%"),
                          uiOutput("id.table")
                   ),
                   column(2,
                          h4("3D plot options"),
                          checkboxInput("burial", "Show burial", value = F),
                          checkboxInput("surface", "Show surfaces", value = F),
                          checkboxInput("c14", "Show C14 samples", value = F),
                          checkboxInput("hearth", "Highlight hearths", value = F),
                          numericInput("point.size", "Point size", 3, min = 1, max=5, width="50%"),
                   )
                 )     # end fluid row
                 
        ),      #end tabPanel
      tabPanel("Section X", 
               fluidRow(
                column(10,
                sliderInput("sectionX", "Min/max X", width="100%",
                                    min=700, max = 1400, value = c(900, 950) )),
                column(1, br(), actionButton("goButtonX", "Draw"))
               ),
               fluidRow(
                 column(7,
                    plotlyOutput("sectionXplot", width = 650, height = 500)
                 ),
                 column(3,
                    imageOutput("cave.mapX", width = "250px", height = "250px")
                 )
               ) #end fluidrow
      ), # end tabPanel
      tabPanel("Section Y", 
               fluidRow(
                 column(10,
                        sliderInput("sectionY", "Min/max Y",  width="95%",
                    min=0, max = 700, value = c(300, 350))),
                  column(1, br(), actionButton("goButtonY", "Draw"))
                 ),
               fluidRow(
                 column(7,
                        plotlyOutput("sectionYplot", width = 650, height = 500)
                 ),
                 column(3,
                     imageOutput("cave.mapY", width = "250px", height = "250px")
                 )
               )#end fluidrow
      ), # end tabPanel
      tabPanel("Summary tables", 
         fluidRow(
           column(2),
           column(5,
                  h4("Remains class and localisation method"),
                  tableOutput("classLocalStats")
           ),
           column(5,
                  h4("Remains by layer"),
                  tableOutput("layersStats")
           ),
         ) #end fluidrow
      ) #end tabPanel
      ), # end  tabsetPanel
    width=10) # end mainPanel
  ) #sidebarLayout  
  ) #endfluidPage
  
) #end  shinyUI
          



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
  geom_path(data=cave.contour, aes(x = x, y = y) )+
  coord_fixed() + xlab("") + ylab("") 
# end cave map
  
  
period1950s.df <- read.csv("data/20299999_PCR_PAVO_Poeymau_coordonnees_gen.csv",
               header = T, stringsAsFactors = F)
# suppression des lignes aux données manquantes
period1950s.df <- period1950s.df[ ! is.na(period1950s.df$layer),]
period1950s.df <- period1950s.df[ ! is.na(period1950s.df$square_x),]

# suppression des localisations douteuses
period1950s.df <- period1950s.df[ ! period1950s.df$uncertain_location == "oui",]

# préparation des localisations horizontales ####

# création d'un identifiant de carré
period1950s.df$square <- paste(period1950s.df$square_y, period1950s.df$square_x, sep = "")

# attribution des x et y min / max de la couche dans le carré ou l'objet se trouve ####
obs.minmax <- group_by(period1950s.df, layer, square) %>% summarize(
x.min.obs = min(xmin, na.rm=T),
x.max.obs = max(xmax, na.rm=T),
y.min.obs = min(ymin, na.rm=T),
y.max.obs = max(ymax, na.rm=T))
# remplacement des valeurs infinies par les min et max (0 et 100)
obs.minmax[ abs(obs.minmax$x.min.obs) == Inf, ]$x.min.obs <- 0
obs.minmax[ abs(obs.minmax$y.min.obs) == Inf, ]$y.min.obs <- 0
obs.minmax[ abs(obs.minmax$x.max.obs) == Inf, ]$x.max.obs <- 100
obs.minmax[ abs(obs.minmax$y.max.obs) == Inf, ]$y.max.obs <- 100

# ajout des valeurs au tableau principal:
period1950s.df <- merge(period1950s.df, obs.minmax, c("layer", "square"), all.x=T)
# identifiant des objets sans x ou y min / max
no.xmin <- is.na(period1950s.df$xmin)
no.ymin <- is.na(period1950s.df$ymin)
no.xmax <- is.na(period1950s.df$xmin) & is.na(period1950s.df$xmax)
no.ymax <- is.na(period1950s.df$ymin) & is.na(period1950s.df$ymax)
# affectation des valeurs min max observées:
period1950s.df[no.xmin,]$xmin <- period1950s.df[no.xmin,]$x.min.obs
period1950s.df[no.xmax,]$xmax <- period1950s.df[no.xmax,]$x.max.obs
period1950s.df[no.ymin,]$ymin <- period1950s.df[no.ymin,]$y.min.obs
period1950s.df[no.ymax,]$ymax <- period1950s.df[no.ymax,]$y.max.obs

# attribution coordonnées xy pour les objets localisés par sous-carrés ----
# création tableau de référence des coordonnées des sous-carrés
subsquare.ref <- data.frame(subsquare = 1:9, 
                      xmin = c(0, 33, 67),
                      xmax = c(33, 66, 100),
                      ymin = c(0,0,0, 33,33,33, 67,67,67),
                      ymax = c(33,33,33, 66,66,66, 100,100,100))
# fusion tableau et tableau de références des sous-carrés:
period1950s.df <- merge(period1950s.df, subsquare.ref, by="subsquare", all.x=T, suffixes = c("", "B"))
# remplacement des coordonnées 
period1950s.df[period1950s.df$subsquare %in% 1:9, c("xmin", "xmax", "ymin", "ymax") ] <- 
period1950s.df[period1950s.df$subsquare %in% 1:9, c("xminB", "xmaxB", "yminB", "ymaxB") ]


# conversion des coordonnées xy relative à chaque carré en des coordonnées 
# générales au site. On ajoute une pondération à chaque carré en x et en y:
period1950s.df$xcorrection <- factor(period1950s.df$square_x,
                    levels = 11:-2,
                    labels = seq(0, 1300, 100) )
period1950s.df$xcorrection <- as.numeric(as.character(period1950s.df$xcorrection))

period1950s.df$ycorrection <- factor(period1950s.df$square_y,
                    levels = c("Y", "Z", LETTERS[1:9]),
                    labels = seq(0, 1000, 100) )
period1950s.df$ycorrection <- as.numeric(as.character(period1950s.df$ycorrection))

# application des correctifs:
period1950s.df$xmin <- period1950s.df$xmin + period1950s.df$xcorrection
period1950s.df$xmax <- period1950s.df$xmax + period1950s.df$xcorrection
period1950s.df$ymin <- period1950s.df$ymin + period1950s.df$ycorrection
period1950s.df$ymax <- period1950s.df$ymax + period1950s.df$ycorrection

# controle:
# summary(period1950s.df[period1950s.df$square == "B2", ]$x_min)
# summary(period1950s.df[period1950s.df$carre_y == "D",]$y_min)
# summary(period1950s.df[period1950s.df$carre_x == "2",]$x_min)

# préparation des localisations verticales ####
layers <- c("alsh", "CS", "CT", "FSH", "CI", "CI ou FIH", "FIH",  "BS",
       "BS jaune", "CPE", "CPE ou BI", "BI")
period1950s.df <- period1950s.df[ period1950s.df$layer %in% layers,]

# attribution des min / max des altitudes dans le carré ou l'objet se trouve ####
selection <- ( ! is.na(period1950s.df$zmin)) & is.na(period1950s.df$zmax)
period1950s.df[selection, ]$zmax <- period1950s.df[selection, ]$zmin

# identification des zmin/max par carré et par couche:
obs.z.minmax <- group_by(period1950s.df, layer, square) %>% summarize(
z.min.obs = min(zmin, na.rm=T),
z.max.obs = max(zmax, na.rm=T))

# identification (complémentaire) des zmin/max par couche:
obs.z.minmax.layer <- group_by(period1950s.df, layer) %>% summarize(
z.min.obs.layer = min(zmin, na.rm=T),
z.max.obs.layer = max(zmax, na.rm=T))
# fusion:
obs.z.minmax <- merge(obs.z.minmax, obs.z.minmax.layer,
                by = "layer", all.x=T)
# correction des valeurs infinies:
sel <- is.infinite(obs.z.minmax$z.min.obs)
obs.z.minmax[sel, ]$z.min.obs <- obs.z.minmax[sel, ]$z.min.obs.layer
obs.z.minmax[sel, ]$z.max.obs <- obs.z.minmax[sel, ]$z.max.obs.layer

# ajout des valeurs au tableau principal
period1950s.df <- merge(period1950s.df,
                        obs.z.minmax,
                        c("layer", "square"), all.x=T)
# identifiant des objets sans x ou y min / max
no.zmin <- is.na(period1950s.df$zmin)
no.zmax <- is.na(period1950s.df$zmin) & is.na(period1950s.df$zmax)
# affectation des valeurs min max observées:
period1950s.df[no.zmin,]$zmin <- period1950s.df[no.zmin,]$z.min.obs
period1950s.df[no.zmax,]$zmax <- period1950s.df[no.zmax,]$z.max.obs

# Gestion des incertitudes ####
# — attribution de coordonnées moyennes pour les objets localisés par volume
# period1950s.df$z_moy <- apply(period1950s.df[, 10:11], 1, mean)
# period1950s.df$x_moy <- apply(period1950s.df[, 12:13], 1, mean)
# period1950s.df$y_moy <- apply(period1950s.df[, 14:15], 1, mean)

# — attribution de coordonnées aléatoires pour les objets
# localisés par volume----
period1950s.df$zrand <- period1950s.df$zmin
period1950s.df$xrand <- period1950s.df$xmin
period1950s.df$yrand <- period1950s.df$ymin

period1950s.df[which(period1950s.df[, "zmin"] != period1950s.df[, "zmax"]), ]$zrand <- 
apply(period1950s.df[which(period1950s.df[, "zmin"] != period1950s.df[, "zmax"]), c("zmin", "zmax") ], 1,
  function(x) sample(x[1]:x[2], 1) )
period1950s.df[which(period1950s.df[, "xmin"] != period1950s.df[, "xmax"]), ]$xrand <- 
apply(period1950s.df[which(period1950s.df[, "xmin"] != period1950s.df[, "xmax"]), c("xmin", "xmax")], 1,
  function(x) sample(x[1]:x[2], 1) )
period1950s.df[which(period1950s.df[, "ymin"] != period1950s.df[, "ymax"]), ]$yrand <- 
apply(period1950s.df[which(period1950s.df[, "ymin"] != period1950s.df[, "ymax"]), c("ymin", "ymax")], 1,
  function(x) sample(x[1]:x[2], 1) )

# — marquage des modes de localisation ####
period1950s.df$localisation_mode <- "point"
period1950s.df[which(period1950s.df[, "zmin"] != period1950s.df[, "zmax"]),]$localisation_mode <- "volume"
period1950s.df[which(period1950s.df[, "xmin"] != period1950s.df[, "xmax"]),]$localisation_mode <- "volume"
period1950s.df[which(period1950s.df[, "ymin"] != period1950s.df[, "ymax"]),]$localisation_mode <- "volume"
# types d'objets ####
period1950s.df$object_type <- tolower(period1950s.df$object_type)


# DATA LIVACHE 1970s####
period1970s.df <- read.csv2("data/20210513_Livache1997poeymau.csv")
# sélection des quatre sous-couches:
period1970s.df <- period1970s.df[ period1970s.df$layer %in% 8:11, ] 
period1970s.df <- period1970s.df[ period1970s.df$square_y %in% c("B", "C", "D", "E"), ] # sélection des carrés
    
# ALTITUDES : sampling avec altitude min et max de la couche alsh
# NB : peut être amélioré en définissant des min/max par carrés ou sous-carrés
# NB : CBH = partie supérieure de FIH
# min/max global pour les z: -108 et -60 
thick.layers <- (108 - 60) / 4 # on divise en 4 "tailles" (=décapages)
period1970s.df$zmin <- NA
period1970s.df$zmax <- NA
period1970s.df[which(period1970s.df$layer == 8), ]$zmin <- 60
period1970s.df[which(period1970s.df$layer == 8), ]$zmax <- 60 + thick.layers
period1970s.df[which(period1970s.df$layer == 9), ]$zmin <- 60 + thick.layers
period1970s.df[which(period1970s.df$layer == 9), ]$zmax <- 60 + 2 * thick.layers
period1970s.df[which(period1970s.df$layer == 10), ]$zmin <- 60 + 2 * thick.layers
period1970s.df[which(period1970s.df$layer == 10), ]$zmax <- 60 + 3 * thick.layers
period1970s.df[which(period1970s.df$layer == 11), ]$zmin <- 60 + 3 * thick.layers
period1970s.df[which(period1970s.df$layer == 11), ]$zmax <- 60 + 4 * thick.layers
# sampling des z aléatoires:
period1970s.df$zrand <- apply(period1970s.df[, c("zmin", "zmax")], 1, function(x) sample(seq(x[1], x[2]), 1) )  

# LOCALISATION HORIZONTALE
# 1) Par carrés :
period1970s.df$x <- factor(period1970s.df$square_x,
                           levels = 11:-2,
                           labels = seq(0, 1300, 100) )
period1970s.df$x <- as.numeric(as.character(period1970s.df$x))
period1970s.df$y <- factor(period1970s.df$square_y,
                           levels = c("Y", "Z", LETTERS[1:9]),
                           labels = seq(0, 1000, 100) )
period1970s.df$y <- as.numeric(as.character(period1970s.df$y))

# 2) par sous-carrés (cf. rapport 1977: adoption d'une sous-division des carrés en 4 sous-carrés)
subsquare.ref <- data.frame(subsquare = c(1:4, NA), 
                      sub_x_min = c(0, 50, 0, 50, 0),
                      sub_x_max = c(50, 100, 50, 100, 100),
                      sub_y_min = c(0,0, 50,50, 0),
                      sub_y_max = c(50,50, 100,100, 100))
period1970s.df <- merge(period1970s.df, subsquare.ref,
                        by.x="subsquare", all.x=T)

# sampling des x et y par tirage aléatoire compris dans les min/max des subsquares:
period1970s.df$xrand <- period1970s.df$x + apply(period1970s.df, 1,
              function(i) sample(seq(i["sub_x_min"], i["sub_x_max"]), 1) )
period1970s.df$yrand <- period1970s.df$y + apply(period1970s.df, 1,
              function(i) sample(seq(i["sub_y_min"], i["sub_y_max"]), 1)  )

period1970s.df$layer <- "alsh"
period1970s.df$localisation_mode <- "volume"
period1970s.df$object_type <- "lithique taillé"
period1970s.df$square <- paste(period1970s.df$square_x,
                               period1970s.df$square_y, sep="")
period1970s.df <- period1970s.df
  
  
df.sub <- reactive({
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
  
df.sub <- Reduce(rbind, df.list)
    
    # type de localisation:
    df.sub <- df.sub[df.sub$localisation_mode %in% input$localisation, ]
    # type d'objets:
    if( ! is.null(input$objects)){
      df.sub <- df.sub[df.sub$object_type %in% input$objects, ]
    }
    
    df.sub$square_y <- factor(df.sub$square_y,
                              levels = c("Y", "Z", "A", "B", "C", "D", "E"))
    df.sub$square_x <- factor(df.sub$square_x, levels = c(7:-2))
    
    
    # réordonnancement des couches
    layers <- c("alsh", "CS", "CT", "FSH", "CI", "CI ou FIH", "FIH",  "BS",
                "BS jaune", "CPE", "CPE ou BI", "BI")
    
    df.sub$layer <- factor(df.sub$layer, levels = layers, labels = layers)
    df.sub$layer.col <- factor(df.sub$layer,
                           levels = layers,
                           labels = c(alsh ="blue", CS="gold", CT="deepskyblue4", FSH="orangered2", 
                                      CI="darkgreen", 'CI ou FSH'="grey30", FIH="tan4", 
                                      BS="gold1", 'BS jaune'="gold3", 
                                      CPE="deepskyblue3",  
                                      'CPE ou BI' = "grey32", BI="orangered3")  )
    df.sub
  })
  
  
output$plotOuput <- renderPlot({  # plot des diagrammes ----
    df.sub <- df.sub()
    # type de localisation
    if(input$point) {selection <- "point"}
    if(input$volume) {selection <- "volume"}
    if(input$point & input$volume) {selection <- c("point", "volume")}

    df.sub <- df.sub[df.sub$localisation_mode %in% selection, ]

    # type d'objets
    if( ! is.null(input$objects)){
      df.sub <- df.sub[df.sub$object_type %in% input$objects, ]
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

# tableau des classes d'objets  #### 
output$classLocalStats <- renderTable({
  df.sub <- df.sub()
  stats.df <- table(df.sub$object_type, df.sub$localisation_mode)

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
    df.sub <- df.sub()
    stats.df <- group_by(df.sub, layer, localisation_mode) %>%
      summarise(n = n())
    stats.df <- dcast(stats.df, layer~localisation_mode, value.var="n")
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
    df.sub <- df.sub()
    
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
      levels(df.sub$layer) <- c(levels(df.sub$layer), "hearth")
      selection <- df.sub$sublayer == "foyer" & df.sub$object_type == ""
      df.sub[selection ,]$layer <- "hearth"
      levels(df.sub$layer.col) <- c(levels(df.sub$layer.col), "black")
    }
    
    fig <- plot_ly(df.sub, x = ~xrand * -1, y = ~yrand * -1, z = ~zrand * -1,
                   color = ~layer,
                   colors = as.character(levels(df.sub$layer.col)),
                   size  = ~point.size,
                   sizes = size.scale,
                   marker = list(symbol = 'square', sizemode = 'diameter'),
                   text = ~paste('id:', id,
                                 '<br>Square:', square,
                                 '<br>Localisation:', localisation_mode,
                                 '<br>Class:', object_type)
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
    } %>%   config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "poeymau3D",
        width = 600, height = 600
      ))
    
    # ajout des surfaces ####
    if(input$surface){
      # sélection des surfaces à calculer :
      layers <- table(df.sub$layer) 
      layers <- names(layers[layers > 100])
        
      # calcul des surfaces:
      surf.list <- lapply(layers, get.surface.model, df=df.sub)
      
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
                   # range =  -c(1400, 700),
                   range =  -c(1400, 400),
                   tickmode = "array",
                   tickvals = -seq(450, 1400, 50),
                   ticktext = c(rbind(levels(df.sub$square_x), ""))
      ),
      yaxis = list(title = 'Y',
                   range = - c(700, 0),
                   tickmode = "array",
                   tickvals = -seq(50, 700, 50),
                   ticktext = c(rbind(levels(df.sub$square_y), ""))
      ),
      zaxis = list(title = 'Depth (m)',
                   range =  c(-800, 0),
                   tickmode = "array",
                   tickvals = - seq(0, 750, 50),
                   ticktext =  c("", "", c(rbind(1:7, "")))
      ),
      aspectmode = "manual", 
      aspectratio = list(x = 1, y = 1, z = 800/700)
    ))
    
  })
  

# Section X ####
min.max.X <- eventReactive(input$goButtonX, {
  seq(input$sectionX[1], input$sectionX[2]) 
})
  
  
output$sectionXplot <- renderPlotly({
  df.sub <- df.sub()
  min.max.X <- min.max.X()
  sectionX.df <- df.sub[df.sub$xrand %in% min.max.X,  ]
  sectionX.df$coordx <- sectionX.df$yrand
  
  plot_ly(sectionX.df, x = ~coordx * -1, y = ~zrand * -1,
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
        filename = "sectionX",
        width = 600, height = 600
    ))  %>%
    add_markers() %>% 
    layout(xaxis = list(title="Y", 
                              zeroline = FALSE, 
                              range=c(0, -750),
                              tickvals = -seq(50, 750, 50),
                              ticktext = c(rbind(levels(df.sub$square_y), ""), "")
                              ),
                 yaxis = list(title="Depth (m)",
                              zeroline = FALSE,
                              tickvals = - seq(0, 800, 100),
                              scaleanchor="x",
                              ticktext =   0:8,
                              range=c(-800,0)
                              )
                 )
})


# Section Y ####

min.max.Y <- eventReactive(input$goButtonY, {
  seq(input$sectionY[1], input$sectionY[2])
})

output$sectionYplot <- renderPlotly({
  df.sub <- df.sub()
  min.max.Y <- min.max.Y()
  
  sectionY.df <- df.sub[df.sub$yrand %in% min.max.Y,  ]
  sectionY.df$coordx <- sectionY.df$xrand
  
  plot_ly(sectionY.df, x = ~coordx, y = ~zrand * -1,
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
        filename = "sectionX",
        width = 600, height = 600
      )) %>%
    add_markers() %>% 
    layout(
    xaxis = list(title="X", 
              zeroline = FALSE,
              tickvals = seq(450, 1400, 50),
              ticktext = c(rbind(levels(df.sub$square_x), "")),
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
})


# cave map ####

output$cave.mapX <- renderPlot({
    cave.map +
      geom_rect(aes(ymin=.5, ymax=7.5, 
                    xmin= input$sectionX[1] / 100 + .5,
                    xmax= input$sectionX[2] / 100 + .5),
                fill="red", alpha=.7)
})

output$cave.mapY <- renderPlot({
  df.sub <- df.sub()
  
  cave.map +
    geom_rect(aes(xmin=4.5, xmax=14.5,
                  ymin= input$sectionY[1] /100 + .5,
                  ymax= input$sectionY[2] /100 + .5 ),
              fill="red", alpha=.7)
})

  # fin du serveur
}

# Run app:
shinyApp(ui = ui, server = server)
