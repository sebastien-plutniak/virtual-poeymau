library(shinythemes)
library(shinyBS)
library(plotly) # for plotlyOutput() in ui.R

ui <- shinyUI(
  fluidPage(
    theme = shinytheme("cosmo"),  # lumen  flatly simplex
    
    sidebarLayout(
      sidebarPanel(
        h4("Virtual Poeymaü v0.3"),
        p(
          "Dev. by: ",
          a("S. Plutniak", href="https://sebastien-plutniak.github.io", target="_blank"),
          br(),
          a("PAVO", href="https://pavo.hypotheses.org/", target="_blank"),
          "project",
        ),
        checkboxGroupInput("localisation", h4("Localisation method"),
                           choices = list("Point" = "point",
                                          "Volume" = "volume"),
                           selected = c("point", "volume")),
        
        checkboxGroupInput("objects", h4("Remains class"),
                           choices = list("Undetermined" = "",
                                          "Stone tools" = "lithique taillé",
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
                            HTML("<div style=width:52%;, align=left>
              <br><br>
   <p><b>Welcome to the <i>Virtual Poeymaü</i> application</b></p>
   <p>
     The <a target=blank href=https://fr.wikipedia.org/wiki/Grotte_du_Poeyma%C3%BC>Poeymaü</a>,
     a cave located in the western Pyrenees
     (in <a target=blank, href=https://www.openstreetmap.org/relation/2120173>Arudy</a>, France),
     was excavated by Georges Laplace and his collaborators from the early 1950s to the late 1980s.
     </p>
   <p>
     The site is currently restudied by the participants in the
     <a target=_blank, href=https://pavo.hypotheses.org><i>PAVO</i></a> research project
     (<i>Préhistoire ancienne de la vallée d'Ossau</i>). 
     Archaeological information was extracted from Laplace's excavation field notes related to
     the 1950s campaigns (archives of the <i><a href=https://musee-prehistoire-eyzies.fr>Musée national de Préhistoire</a></i>). 
     About 15,000 objects are documented and located in the cave's space. 
   </p>
   <p>
      The <i>Virtual Poeymaü</i> is a web application to explore this data set. 
      It can  
      filter the data, 
      plot them with interactive 3D visualisation, 
      generate and download 2D plans, longitudinal and transverse sections,  
      generate summary tables, and
      explore the timeline of the excavation from 1951 to 1985.
   </p>
   <p> The remains are located by two methods, called:
      <ul>
          <li>“point”: using the xyz coordinates when they are given in the fieldwork documentation and,</li>
          <li>“volume”: if only ranges of xyz coordinates are given, by sampling random values within these ranges.</li>
      </ul>
   </p>
    <hr>
    For more information:
    <ul>
      <li> <b>Software:</b> 
        <ul>
           <li> <b>R code:</b> <a target=blank, href=https://github.com/sebastien-plutniak/virtual-poeymau>github repository</a>.</li>
           <li> <b>reference:</b> Plutniak, S. 2021, “Virtual Poeymaü: a web application to explore the archaeological data from the excavation archives of the Poeymaü cave (France)”,  <i>Zenodo</i>. DOI: <a target=blank, href=https://doi.org/10.5281/zenodo.4765693>10.5281/zenodo.4765693</a>.</li>
        </ul>
        <li> <b>Documentation:</b> 
        <ul>
          <li>Plutniak, S. 2021, “Virtual Poeymaü : visualisation 3D interactive des données archéologiques extraites des archives de fouilles de la grotte du Poeymaü », <i>Préhistoire ancienne de la vallée d'Ossau</i>, HDL: <a target=blank, href=https://hdl.handle.net/10670/1.juoxn4>10670/1.juoxn4</a>.</li>
          <li>Plutniak, S. 2020, “Archives des fouilles de la grotte du Poeymaü (1951–1956) : informatisation et spatialisation des restes archéologiques”. Pétillon, J.-M. et B. Marquebielle (dir.), <i>Préhistoire ancienne de la vallée d'Ossau. Paléoenvironnement et sociétés de chasseurs-collecteurs dans le piémont pyrénéen. Projet collectif de recherche. Bilan 2020</i>, DRAC-SRA Nouvelle-Aquitaine, <a href=https://hal.archives-ouvertes.fr/hal-03092989>hal-03092989</a>.</li>
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
                            plotlyOutput("plot3d",  width = 800, height = 650),
                            numericInput("id", "id", 1, min=1, max=16000, width="12%"),
                            uiOutput("id.table")
                     ),
                     column(2,
                            br(),
                            h4("3D plot options"),
                            checkboxInput("burial", "Show burial", value = F),
                            checkboxInput("c14", "Show C14 samples", value = F),
                            checkboxInput("hearth", "Highlight hearths", value = F),
                            checkboxInput("surface", "Compute surfaces", value = F),
                            checkboxInput("cxhull", "Compute hull", value = F),
                            numericInput("point.size", "Point size", 2, min=1, max=5,
                                         width="50%"),
                            sliderInput("ratio", "Vertical ratio", width="100%", sep = "",
                                        min=1, max=3, value=2, step=.2)
                     )  # end column
                   )  # end fluid row
          ),      #end tabPanel
          tabPanel("Plan", 
                   fluidRow(
                     column(11,
                            sliderInput("planZ", "Z: min/max",  width="100%",  sep = "",
                                        min=0, max = 800, value = c(300, 350))
                     ),
                     column(1, br(),
                            actionButton("goButtonZ", "Draw"),)
                   ),
                   fluidRow(column(12, align="center",
                                   plotlyOutput("planZplot",   width = "95%"),
                                   h3("Density contour"),
                                   imageOutput("density.plan", width = "100%"),
                                   downloadButton("download.density.plan", "Download contour plot (svg)")
                   ) #end column
                   ) #end fluidrow
          ), # end tabPanel             
          tabPanel("Section X", 
                   fluidRow(
                     column(11,
                            sliderInput("sectionYx", "X: min/max", width="100%",  sep = "",
                                        min=400, max = 1400, value = c(400, 1400)),
                            sliderInput("sectionYy", "Y: min/max",  width="100%",  sep = "",
                                        min=0, max = 800, value = c(550, 600))
                     ),
                     column(1,
                            br(), actionButton("goButtonY", "Draw"))
                   ),
                   fluidRow(
                     column(9,
                            plotlyOutput("sectionYplot", width = 650, height = 500)
                     ),
                     column(3,
                            imageOutput("cave.mapY", width = "250px", height = "250px")
                     )
                   )#end fluidrow
          ), # end tabPanel
          tabPanel("Section Y", 
                   fluidRow(
                     column(11,
                            sliderInput("sectionXx", "X: min/max", width="100%",  sep = "",
                                        min=400, max = 1400, value =  c(1050, 1100)),
                            sliderInput("sectionXy", "Y: min/max",  width="100%",  sep = "",
                                        min=0, max = 800, value = c(0, 800))
                     ),
                     column(1, br(),
                            actionButton("goButtonX", "Draw"),)
                   ),
                   fluidRow(
                     column(9,
                            plotlyOutput("sectionXplot", width = 650, height = 500)
                     ),
                     column(3,
                            imageOutput("cave.mapX", width = "250px", height = "250px")
                     )
                   ) #end fluidrow
          ), # end tabPanel
          tabPanel("Summary tables", 
                   fluidRow(
                     column(2),
                     column(5,
                            h4("Remains class and localisation method"),
                            tableOutput("classLocalStats")
                     ),
                     column(5,
                            h4("Remains by layer and localisation method"),
                            tableOutput("layersStats")
                     ),
                   ) #end fluidrow
          ), #end tabPanel
          tabPanel("Excavation timeline", 
                   sliderInput("history.date", "Year", width="100%",  sep = "",
                               min=1951, max = 1985, value = 1951),
                   fluidRow(
                     column(7, 
                            imageOutput("cave.map.history", width = "100%", height = "500px")),
                     column(5,
                            imageOutput("cave.map.history.grid", width = "100%", height = "400px")),
                   ) #end fluidrow
          ) #end tabPanel
        ), # end  tabsetPanel
        width=10) # end mainPanel
    ) #sidebarLayout  
  ) #endfluidPage
) #end  shinyUI
