# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinythemes)
library(shinyBS)

shinyUI(fluidPage(
  useShinyjs(),
  navbarPage(
    theme = shinytheme("yeti"),
    "rna",
    position = "fixed-top",
    collapsible = TRUE,
    fluid = TRUE,
    header = singleton(tags$head(
      tags$head(tags$script(src = "/scripts/fornac.js")),
      tags$head(tags$script(src = "/js/vis.js")),
      tags$head(tags$script(src = "/js/sequence-viewer.js")),
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "my_theme.css")
      )
    )),
    #tab rnaCARD
    tabPanel(
      "rnaCARD",
      div(
        class = "container-fluid main-container",
        
        # rnaCARD input panel
        conditionalPanel(
          condition = "input.submit == '0' && input.start == '0'",
          style = "background-image: url(g38710.png); background-size: cover;",
          
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),

          
          fluidRow(
            column(6, align="center", offset =  3,
                   h1(strong("rnaCARD"), class = "text-muted")
            )
          ),
          br(),
          br(),
          fluidRow(
            column(6, align="center", offset =  3,
                   h3("Welcome in rnaCARD webserver!")
            )
          ),
          
          fluidRow(
            column(6, align="center", offset =  3,
                   h4("Click below to start comparative analysis of RNA secondary structures.")
            )
          ),

          br(),
          br(),
          fluidRow(
            column(6, align="center", offset =  3,
                   actionButton(
                     "start",
                     label = "Start analysis",
                     class = "btn-success  btn-lg",
                     style = "width: 200px;"
                   )
            )
          ),
        div(style="height: 400px;")
        ),
        
        conditionalPanel(
          condition = "input.submit == '0' && input.start == '1'",
          h3(strong("rnaCARD"), class = "text-muted"),
          hr(),
          splitLayout(
            #paste input
            div(
              style = "margin: 30px; background-color: #ffffff; border-width: 0px; border-color: #ffffff;",
              h5(strong("Paste sequence & structures")),
              verticalLayout(
                textAreaInput(
                  "sequence",
                  h6(strong("NUCLEOTIDE SEQUENCE"), class = 'text-primary'),
                  "TAATGCCTTTGTTTGGCCAAGCTATGTGCAAATATCACAAATTAAAAATTTGGTTAAGCAGTTAGGCTGGACCTAATATTTTAGAAAAACCTAATTTTTTTTGTGGACCCATTTTCGATATTTACTCACAAATGGAATTCAAGGGGAACAACTTCGGTCTCAGCACTTTAATTATTCTTCTCGTTCCCACCTAATTTCGCAATTTATTGTCCTTGACTTCTACCACGAGAAAAAAATTAAGAAAATGCAACGCTGCCCGTGCAGGGTTTTCTGAGCGGGATGAAAAAATCAGACAAATATCCAAGTTATGAGTAATTACTTTGTTGGAAGGAGGGAGCAGAGGATAAGGAAATTCTTAAAACTGTTATGTATATAAAGGAAGAACCATTTCTAGTTATTTCACTTTTTGATACTTGTCAACTATCTTAGTAAAAATACAGAACTCTATAAAGAACCACAGAAAAATCGACAGCAATGACAAGCATTGACATTAACAACTTACAAAATACCTTTCAACAAGCTATGAATATGAGCGGCTCCCCAGGCGCTGT",
                  resize = "vertical"
                ),
                textAreaInput(
                  "str1",
                  h6(strong("STRUCTURE 1"), class = 'text-primary'),
                  "(((((.....(((((((((((...((((.......)))).........))))))))))).((((((.....))))))...((((......)))).(((((((((((..((((((..((.......))..)))))).......((((...........))))......((..((((.(((((.((((((.(((......((((....)))).(((((.((((......)))).....))))).....(((....)))((((.((((....)))).))))(((......)))..........((((.....((((....))))..)))).))).)))))))))))))))..)).(((((......(((.(((((.....((((......)))).(((......)))....(((....))).((((...))))....))))).)))......)))))))))))))))).......(((((.....))))).)))))....((((....(((.(((.....))).))).....))))((((.((...)).)))).",
                  resize = "vertical"
                ),
                textAreaInput(
                  "str2",
                  h6(strong("STRUCTURE 2"), class = 'text-primary'),
                  "(((((.....(((((((((((...((((.......)))).........))))))))))).((((((.....))))))..................(((((((((((..((((((..((.......))..)))))).......((((...........))))......((..((((.(((((.((((((.(((......((((....)))).(((((.((((......)))).....))))).....(((....)))((((.((((....)))).))))(((......)))..........((((.....((((....))))..)))).))).)))))))))))))))..)).(((((......(((.(((((..(((((((((((........))).)))).)))).((((....))))...............))))).)))......)))))))))))))))).......(((((.....))))).)))))....((((....(((.(((.....))).))).....))))((((.((...)).)))).",
                  resize = "vertical"
                )
              ),
              #actionButton("submit", label = "Submit",  class = "btn-primary")
            ),
            # file input
            div(
              style = "margin: 30px; background-color: #ffffff; border-width: 0px; border-color: #ffffff;",
              h5(strong("Or upload input file")),
              br(),
              fileInput("file2", "", multiple = FALSE)
            )
          ),
          fluidRow(column(5, ""), column(
            2,
            actionButton(
              "submit",
              label = "Submit",
              class = "btn-primary",
              style = "width: 200px;"
            )
          ), #column(1,actionButton("load_example", label = "Example",  class = "btn-link")),
          column(3, ""))
        ),
        #rnaCARD input panel end
        
        #rnaCARD results panel
        conditionalPanel(
          condition = "input.submit == '1' && input.start == '1'",
        fluidRow(
          column(
            7,
            h3(strong("rnaCARD"), class = "text-muted"),
            span(textOutput("done"), style = "color:white")
          ),
          column(
            2,
            actionButton("new_analysis", "New analysis", class = "btn-success")
          ),
          
          column(
            3,
            downloadButton("downloadData", "Download results", style =
                             'font-size:90%')
          ),
          column(1, "")
        ),
        hr(),
        br(),
        uiOutput("selectID_card"),
        br(),
        h5(strong('LIST of MOTIFS'), class = "text-primary"),
        br(),
          DTOutput('table_motifs'),
          br(),
          #Forna overview
          h5(strong(textOutput('mot_num')), class = "text-primary"),
          br(),
          
          splitLayout(h5("#structure1"), h5("#structure2")),
          splitLayout(
            div(class = "breadcrumb", id = "rna_m", style = "max-height: 600px; border-radius: 5px; background-color: white; border-width: 4px;"),
            div(class = "breadcrumb", id = "rna_m2", style = "max-height: 600px;border-radius: 5px; background-color: white; border-width: 4px;")
          ),
          br(),
          h5(strong("STRUCTURES OVERVIEW", class = "text-primary")),
          br(),
          splitLayout(
            div(class = "breadcrumb", id = "rna1", style = "max-height: 600px; border-radius: 5px; background-color: white; border-width: 4px;"),
            div(class = "breadcrumb", id = "rna2", style = "max-height: 600px; border-radius: 5px; background-color: white; border-width: 4px;")
          )
          
        )
      ) #results panel view end
    ),
    #tab rnaCARD end
    tabPanel(
      "Help",
      div(
        class = "container-fluid main-container",
        h5("rnaCARD is availble in command-line version on github."),
        br(),
        actionButton(
          inputId = 'ab1',
          label = "Go to github",
          icon = icon("send"),
          onclick = "window.open('https://github.com/mplucinska/rnaCARD')"
        )
      )
    )#end tab Help
  )# navbarPage end
)# fluidPage end
)# UI end