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
                            conditionalPanel(condition = "input.submit == '0'",
                                             splitLayout(
                                               #paste input
                                               div(style = "margin: 30px; background-color: #ffffff; border-width: 0px; border-color: #ffffff;",
                                                   h5(strong("Paste sequence & structures")),
                                                   verticalLayout(
                                                     textAreaInput("sequence", h6(strong("NUCLEOTIDE SEQUENCE"), class = 'text-primary'),  "TAATGCCTTTGTTTGGCCAAGCTATGTGCAAATATCACAAATTAAAAATTTGGTTAAGCAGTTAGGCTGGACCTAATATTTTAGAAAAACCTAATTTTTTTTGTGGACCCATTTTCGATATTTACTCACAAATGGAATTCAAGGGGAACAACTTCGGTCTCAGCACTTTAATTATTCTTCTCGTTCCCACCTAATTTCGCAATTTATTGTCCTTGACTTCTACCACGAGAAAAAAATTAAGAAAATGCAACGCTGCCCGTGCAGGGTTTTCTGAGCGGGATGAAAAAATCAGACAAATATCCAAGTTATGAGTAATTACTTTGTTGGAAGGAGGGAGCAGAGGATAAGGAAATTCTTAAAACTGTTATGTATATAAAGGAAGAACCATTTCTAGTTATTTCACTTTTTGATACTTGTCAACTATCTTAGTAAAAATACAGAACTCTATAAAGAACCACAGAAAAATCGACAGCAATGACAAGCATTGACATTAACAACTTACAAAATACCTTTCAACAAGCTATGAATATGAGCGGCTCCCCAGGCGCTGT", resize = "vertical"),
                                                     textAreaInput("str1", h6(strong("STRUCTURE 1"), class = 'text-primary'), "(((((.....(((((((((((...((((.......)))).........))))))))))).((((((.....))))))...((((......)))).(((((((((((..((((((..((.......))..)))))).......((((...........))))......((..((((.(((((.((((((.(((......((((....)))).(((((.((((......)))).....))))).....(((....)))((((.((((....)))).))))(((......)))..........((((.....((((....))))..)))).))).)))))))))))))))..)).(((((......(((.(((((.....((((......)))).(((......)))....(((....))).((((...))))....))))).)))......)))))))))))))))).......(((((.....))))).)))))....((((....(((.(((.....))).))).....))))((((.((...)).)))).", resize = "vertical"),
                                                     textAreaInput("str2", h6(strong("STRUCTURE 2"), class = 'text-primary'), "(((((.....(((((((((((...((((.......)))).........))))))))))).((((((.....))))))..................(((((((((((..((((((..((.......))..)))))).......((((...........))))......((..((((.(((((.((((((.(((......((((....)))).(((((.((((......)))).....))))).....(((....)))((((.((((....)))).))))(((......)))..........((((.....((((....))))..)))).))).)))))))))))))))..)).(((((......(((.(((((..(((((((((((........))).)))).)))).((((....))))...............))))).)))......)))))))))))))))).......(((((.....))))).)))))....((((....(((.(((.....))).))).....))))((((.((...)).)))).", resize = "vertical")
                                                   ),
                                                   actionButton("submit", label = "Submit",  class = "btn-primary")
                                               ),
                                               # file input
                                               div(style = "margin: 30px; background-color: #ffffff; border-width: 0px; border-color: #ffffff;",
                                                   h5(strong("Or upload input file")),
                                                   br(),
                                                   fileInput("file1", "", multiple = FALSE),
                                                   actionButton("submit_file", label = "Submit",  class = "btn-primary")
                                               )
                                             )
                            ), #rnaCARD input panel end
                            
                            #rnaCARD results panel
                            conditionalPanel(condition = "input.submit == '1'",
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
                                             h5(strong('LIST of MOTIFS'), class = "text-primary"),
                                             br(),
                                             DTOutput('table_motifs'),
                                             br(),
                                             #Forna overview
                                             h5(strong(textOutput('mot_num')), class = "text-primary"),
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
                                             bsCollapse(id = "collapseExample", multiple = TRUE, open = "sequence & structures",
                                                        bsCollapsePanel(h5("sequence & structures"),
                                                                        #fluidRow(column(width = 1, p("SEQUENCE")), column(width = 10, p(textOutput("mot_seq")))),
                                                                        h6("SEQUENCE"),
                                                                        h6(textOutput("mot_seq")),
                                                                        br(),
                                                                        #fluidRow(column(width = 1, p("STRUCTURE 1")), column(width = 10, p(textOutput("mot_str1")))),
                                                                        h6("STRUCTURE 1"),
                                                                        h6(textOutput("mot_str1")),
                                                                        br(),
                                                                        #fluidRow(column(width = 1, p("STRUCTURE 2")), column(width = 10, p(textOutput("mot_str2"))))
                                                                        h6("STRUCTURE 2"),
                                                                        h6(textOutput("mot_str2"))
                                                        )
                                             ),
                                             splitLayout(h5("#structure1"), h5("#structure2")),
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
                        div(class = "container-fluid main-container",
                            h4(strong("About rnaNORM"), class = "text-muted"),
                            br(),
                            br(),
                            p("rnaNORM is a method of signal calculation that eliminate read distribution bias and prevent underestimation of reactivity."),
                            p("Check 'How it works?' section for more info. There is also command-line version availible."),
                            br(),
                            h5(strong("How rnaNORM works?"), class = "text-primary"),
                            hr(),
                            img(src='grafika_rnaNORM_small.png', style="display: block; margin-right: auto; width: 70%; align:left"),
                            br(),
                            tags$ol(
                              tags$li(h5("Read distrinbution is identified by analysis of regression shift towards control sample. First stage of normalization process is calculation of log2 fold change of modified counts with respect to control counts.")), 
                              tags$li(h5("Distribution of log2 fold changes is used for density function estimation. It enables to calculate value of fold change characteristic for background, which next is used as normalization factor. In background signal after correction, normalized control counts are approximate value of counts in modified sample. Only positions above standard deviation are used for reactivity calculation.")), 
                              tags$li(h5("Normalized counts are used for reactivity calculation for each position by substraction of normalized control signal from modified. Final profile is scaled with 2/8 normalization (Deigan et al., 2009)."))
                            ),
                            #hr(),
                            #h5("rnaNORM is method reactivity calculation that eliminate read distribution bias and prevent underestimation of reactivity."),
                            br(),
                            h5(strong("Input format"), class = "text-primary"),
                            hr(),
                            h5("Input consist of 4 tab delimited columns (default order):"),
                            fluidRow(
                              column(1,"id"),
                              column(2, "position"),
                              column(3, "number_of_stops_in_control"),
                              column(2, "number_of_stops_in_modified")
                            ),
                            br(),
                            h5("User can also specify column order in file."),
                            h5("It is possible to upload multiple transcripts."),
                            h5("Script counting stops from chemical probing experiments is availible on github. Input file is indexed BAM file with reads aligned to transcriptome."),
                            br(),
                            h5(strong("Output"), class = "text-primary"),
                            hr(),
                            h5("Transcript of interest can be chosen from drop-down list."),
                            br(),
                            tags$ol(
                              tags$li(h5("Normalized reactivity for each position in transcript visualized on bar plot.")),
                              tags$li(h5("Table with input counts, normalized counts, reactivity and information whether reactivity passed filter or not (P - passed, F - FAILED)")),
                              tags$li(h5("Dot plot representing counts correlation before and after normalization. Each dot on the plot referes to number of counts in control and modified sample in particular position in transcript."))
                            ),
                            br(),
                            h5(strong("Downloading results"), class = "text-primary"),
                            hr(),
                            h5("Button 'Download results for transcript' enables downloading text file with results for selected transcript from the list."),
                            h5("In order to download results for all uploaded transcripts click on 'Calculate all'. It starts calculation for all transcripts. It may take a while. After finishing calculation, button changes to green 'Download all' button."),
                            br(),
                            h5(strong("Download command-line version of rnaNORM from github"), class = "text-primary"),
                            hr(),
                            h5("rnaNORM is availble in command-line version on github."),
                            br(),
                            actionButton(inputId='ab1', label="Go to github", 
                                         icon = icon("send"), 
                                         onclick ="window.open('https://github.com/mplucinska/rnaNORM')")
                        )
               )#end tab Help
    )# navbarPage end
  )# fluidPage end
)# UI end