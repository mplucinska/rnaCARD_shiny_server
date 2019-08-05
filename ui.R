# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinythemes)
library(shinyBS)

shinyUI(
  fluidPage(
    useShinyjs(),
    navbarPage(theme = shinytheme("yeti"), "rna", position = "fixed-top", collapsible = TRUE,
               fluid = TRUE,
               header=singleton(tags$head(tags$head(tags$script(src="/scripts/fornac.js")),
                                          tags$head(tags$script(src="/js/vis.js")),
                                          tags$head(tags$script(src="/js/sequence-viewer.js")),
                                          tags$head(
                                            tags$link(rel = "stylesheet", type = "text/css", href = "my_theme.css")
                                          ))),
               #tab rnaCARD
               tabPanel("rnaCARD",
                        div(class = "container-fluid main-container",
                            h3(strong("rnaCARD"), class = "text-muted"),
                            hr(),
                            # rnaCARD input panel
                            conditionalPanel(condition = "input.submit == '0' & input.load_example == '0'",
                                             splitLayout(
                                               #paste input
                                               div(style = "margin: 30px; background-color: #ffffff; border-width: 0px; border-color: #ffffff;",
                                                   h5(strong("Paste sequence & structures")),
                                                   verticalLayout(
                                                     textAreaInput("sequence", h6(strong("NUCLEOTIDE SEQUENCE"), class = 'text-primary'), resize = "vertical"),
                                                     textAreaInput("str1", h6(strong("STRUCTURE 1"), class = 'text-primary'), resize = "vertical"),
                                                     textAreaInput("str2", h6(strong("STRUCTURE 2"), class = 'text-primary'), resize = "vertical")
                                                   )
                                               ),
                                               # file input
                                               div(style = "margin: 30px; background-color: #ffffff; border-width: 0px; border-color: #ffffff;",
                                                   h5(strong("Or upload input file")),
                                                   br(),
                                                   fileInput("file1", "", multiple = FALSE)
                                               )
                                             ),
                                             fluidRow(
                                               column(5, ""),
                                               column(2,actionButton("submit", label = "Submit",  class = "btn-primary",  style = "width: 200px;")),
                                               column(1,actionButton("load_example", label = "Example",  class = "btn-link"))
                                             )
                            ), #rnaCARD input panel end
                            
                            #rnaCARD results panel
                            conditionalPanel(condition = "input.submit == '1' || input.load_example == '1'",
                                             h5(strong('What are you looking for?'), class = "text-muted"),
                                             radioGroupButtons(inputId = "mode", 
                                                               label = " ",
                                                               choices = c("similar motifs", "differential motifs"), 
                                                               status = "primary", 
                                                               checkIcon = list(yes = icon("ok", lib = "glyphicon"),
                                                                                no = icon("remove",lib = "glyphicon"))),
                                             br(),
                                             uiOutput("selectID_card"),
                                             br(),
                                             bsCollapse(id = "collapseExample", multiple = TRUE, open = "sequence & structures",
                                                        bsCollapsePanel(h5("sequence & structures"),
                                                                        #fluidRow(column(width = 1, p("SEQUENCE")), column(width = 10, p(textOutput("mot_seq")))),
                                                                        h6("SEQUENCE"),
                                                                        div(textOutput("mot_seq"), style = "font-family: Courier,courier;"),
                                                                        br(),
                                                                        #fluidRow(column(width = 1, p("STRUCTURE 1")), column(width = 10, p(textOutput("mot_str1")))),
                                                                        h6("STRUCTURE 1"),
                                                                        div(textOutput("mot_str1"), style = "font-family: Courier,courier;"),
                                                                        br(),
                                                                        #fluidRow(column(width = 1, p("STRUCTURE 2")), column(width = 10, p(textOutput("mot_str2"))))
                                                                        h6("STRUCTURE 2"),
                                                                        div(textOutput("mot_str2"), style = "font-family: Courier,courier;")
                                                        )
                                             ),
                                             h5(strong('LIST of MOTIFS'), class = "text-primary"),
                                             br(),
                                             DTOutput('table_motifs'),
                                             br(),
                                             #Forna overview
                                             h5(strong(textOutput('mot_num')), class = "text-primary", id = "motif_count"),
                                             br(),
                                             #div(class = 'alert', style = "display: block; max-height: 600px; border-width: 6px; border-color: #823329; border-radius: 15px;",
                                             #    p("SEQUENCE"),
                                             #    p(textOutput("mot_seq")),
                                             #    br(),
                                             #    p("STRUCTURE 1"),
                                             #    p(textOutput("mot_str1")),
                                             #    br(),
                                             #    p("STRUCTURE 2"),
                                             #    textOutput("mot_str2")
                                             #),
                                             splitLayout(h5("#structure1"), h5("#structure2"), id = 'motif_header'),
                                             splitLayout(
                                               div(class = "breadcrumb",id = "rna_m", style = "max-height: 600px; border-radius: 5px; background-color: white; border-width: 4px;"),
                                               div(class = "breadcrumb",id = "rna_m2", style = "max-height: 600px;border-radius: 5px; background-color: white; border-width: 4px;")
                                             ),          
                                             br(),
                                             h5(strong("STRUCTURES OVERVIEW", class = "text-primary")),
                                             br(),
                                             splitLayout(
                                                div(class = "breadcrumb",id = "rna1", style = "max-height: 600px; border-radius: 5px; background-color: white; border-width: 4px;"),
                                                div(class = "breadcrumb",id = "rna2", style = "max-height: 600px; border-radius: 5px; background-color: white; border-width: 4px;")
                                             )
                            )
                        ) #results panel view end
               ), #tab rnaCARD end
               tabPanel("Help",
                h3()
               )#end tab Help
    )# navbarPage end
  )# fluidPage end
)# UI end