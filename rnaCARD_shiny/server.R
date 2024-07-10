library(ggplot2)
library(ggthemes)
library(shiny)
library(plotly)
library(DT)
library(zip)
library(shinyalert)

shinyServer(function(input, output, session) {
  observeEvent(input$submit, {
    output$sequence <- renderText({
      input$sequence
    })
    sessions_id <- stringi::stri_rand_strings(1, 5)

    if (!is.null(input$file2)) {
        
        input_file_card <- input$file2$datapath
      
    } else {
      
      input_file_card <- paste0("www/working_dir/", sessions_id, "_card_input")
      
      fileConn<-file(input_file_card)
      writeLines(c(
        ">sequence1",
        paste0("seq\t", input$sequence),
        paste0("str1\t", input$str1),
        paste0("seq\t", input$str2)
      ), fileConn)
      close(fileConn)
      
    }
    a <- system(
        paste0('python3 rnaCARD.py -i ', input_file_card, ' --match --os --prefix www/working_dir/', sessions_id)
    )
          
    if (a == "1"){
      shinyalert("Error", "Incorrect input format. Correct format described in Help section", type = "warning", callbackR = function(x) { session$reload() })
    } else {
      output_card_out_match <- paste0("www/working_dir/", sessions_id, "_matched_whole_transcripts.txt")
      output_card_motifs_match <- paste0("www/working_dir/", sessions_id, "_matched_motifs_out.txt")

      if(file.info(output_card_out_match)$size == 0){
        shinyalert("Error", "Incorrect input format. Correct format described in Help section", type = "warning", callbackR = function(x) { session$reload() })
      } else {
        card_out <-
          read_delim(
            output_card_out_match,
            "\t",
            escape_double = FALSE,
            col_names = FALSE,
            trim_ws = TRUE
          )
        
        card_motifs <-
          read_delim(
            output_card_motifs_match,
            "\t",
            escape_double = FALSE,
            col_names = FALSE,
            trim_ws = TRUE
          )
        
        colnames(card_motifs) <-
          c(
            "ID",
            "motif_number",
            "start position(s)",
            "end position(s)" ,
            "shape",
            "sequence",
            "structure 1",
            "structure 2",
            "s1",
            "e1",
            "s2",
            "e2"
          )
        
        state = reactiveValues(choice = unique(card_out$X1))
        output$selectID_card <- renderUI({
          selectInput("selected_tr", strong(h5("Select transcript ID:")), state$choice, selected = 1)
        })
        
        
        list_files <- list.files("./www/working_dir/", full.names = T)
        selected <- list_files[grep(sessions_id, list_files)]
        
        dir.create(tmp <- tempfile())
        
        file.copy(from=selected, to=tmp, 
                  overwrite = TRUE, recursive = FALSE, 
                  copy.mode = TRUE)
        
        zip_name <- paste0(tmp, "/", sessions_id, "_", "rnaCARD_outputs", ".zip")
        zip::zip(zipfile = zip_name, files = tmp, include_directories = FALSE)
        
        out_name <-  paste0(sessions_id, "_", "rnaCARD_outputs", ".zip")
        output$downloadData <- downloadHandler(
          filename <- function() {
            paste(out_name)
          },
          content <- function(file) {
            output_file_name <- paste(zip_name)
            file.copy(output_file_name, file)
          }
        )
        
        
        observeEvent(input$selected_tr, {
          state$val <- input$selected_tr
          selected_transcript <-
            subset(card_out, card_out$X1 == state$val)
          

          
          draw_overview_structure(
            selected_transcript$X3,
            selected_transcript$X4,
            selected_transcript$X5,
            selected_transcript$X2
          )
          
          selected_motifs <-
            subset(card_motifs, card_motifs$ID == state$val)
          print(selected_motifs)
          create_table(selected_motifs)

          observeEvent(input$table_motifs_rows_selected, {
            if(input$table_motifs_rows_selected > nrow(selected_transcript)){
              sel <- 1
            } else {
              sel <- input$table_motifs_rows_selected 
            }
            if (!is.null(input$table_motifs_rows_selected)) {
              draw_motif(selected_motifs[sel, ], selected_transcript)
            }
          })
        })
        
        output$mot_seq = renderPrint({
          s = subset(card_out, card_out$X1 == state$val)$X2
          if (length(s)) {
            cat(s, sep = ', ')
          }
        })
        
        output$mot_str1 = renderPrint({
          s = subset(card_out, card_out$X1 == state$val)$X3
          if (length(s)) {
            cat(s, sep = ', ')
          }
        })
        
        output$mot_str2 = renderPrint({
          s = subset(card_out, card_out$X1 == state$val)$X3
          if (length(s)) {
            cat(s, sep = ', ')
          }
        })
        
        output$mot_s = renderPrint({
          s = subset(card_out, card_out$X1 == state$val)$X4
          if (length(s)) {
            cat(s, sep = ', ')
          }
        })
        
        output$id_tr = renderPrint({
          s = input$selected_tr
          if (length(s)) {
            cat(s, sep = ', ')
          }
        })
        
        output$mot_num = renderPrint({
          s = input$table_motifs_rows_selected
          if (length(s)) {
            cat('MOTIF #')
            cat(s, sep = ', ')
          }
        })
      }
      }

    

  }, ignoreInit = TRUE)
  
  observeEvent(input$new_analysis, {
    session$reload()
  })
  
  colors <-
    c(
      "#ef9a9a",
      "#F48FB1",
      "#CE93D8",
      "#B39DDB",
      "#9FA8DA",
      "#90CAF9",
      "#81D4FA",
      "#80DEEA",
      "#80CBC4",
      "#A5D6A7",
      "#C5E1A5",
      "#E6EE9C",
      "#FFF59D",
      "#FFE082",
      "#FFCC80",
      "#FFAB91",
      "#BCAAA4",
      "#B0BEC5"
    )
  
  
  draw_motif <- function(motif, transcript) {

    if (grepl(" ", motif$s1, fixed = TRUE)) {
      seq_m1 <-
        paste(
          substr(
            transcript$X2,
            strsplit(motif$s1, " ")[[1]][1],
            strsplit(motif$e1, " ")[[1]][1]
          ),
          "&",
          substr(
            transcript$X2,
            strsplit(motif$s1, " ")[[1]][2],
            strsplit(motif$e1, " ")[[1]][2]
          ),
          sep = ""
        )
      seq_m2 <-
        paste(
          substr(
            transcript$X2,
            strsplit(motif$s2, " ")[[1]][1],
            strsplit(motif$e2, " ")[[1]][1]
          ),
          "&",
          substr(
            transcript$X2,
            strsplit(motif$s2, " ")[[1]][2],
            strsplit(motif$e2, " ")[[1]][2]
          ),
          sep = ""
        )
      col_m1 <-
        changeColors(paste(strsplit(transcript$X5, " ")[[1]][strsplit(motif$s1, " ")[[1]][1]:strsplit(motif$e1, " ")[[1]][1]], collapse = " "))
      col_m2 <-
        changeColors(paste(strsplit(transcript$X5, " ")[[1]][strsplit(motif$s2, " ")[[1]][1]:strsplit(motif$e2, " ")[[1]][1]], collapse = " "))
      str_m1 <-
        paste(
          substr(
            transcript$X3,
            strsplit(motif$s1, " ")[[1]][1],
            strsplit(motif$e1, " ")[[1]][1]
          ),
          "&",
          substr(
            transcript$X3,
            strsplit(motif$s1, " ")[[1]][2],
            strsplit(motif$e1, " ")[[1]][2]
          ),
          sep = ""
        )
      str_m2 <-
        paste(
          substr(
            transcript$X4,
            strsplit(motif$s2, " ")[[1]][1],
            strsplit(motif$e2, " ")[[1]][1]
          ),
          "&",
          substr(
            transcript$X4,
            strsplit(motif$s2, " ")[[1]][2],
            strsplit(motif$e2, " ")[[1]][2]
          ),
          sep = ""
        )
    } else{
      seq_m1 <- substr(transcript$X2, motif$s1, motif$e1)
      seq_m2 <- substr(transcript$X2, motif$s2, motif$e2)
      col_m1 <-
        changeColors(paste(strsplit(transcript$X5, " ")[[1]][motif$s1:motif$e1], collapse = " "))
      col_m2 <-
        changeColors(paste(strsplit(transcript$X5, " ")[[1]][motif$s2:motif$e2], collapse = " "))
      str_m1 <- substr(transcript$X3, motif$s1, motif$e1)
      str_m2 <- substr(transcript$X4, motif$s2, motif$e2)
    }
    
    session$sendCustomMessage("message_motif",
                              list(
                                seq_m = seq_m1,
                                str_m = str_m1,
                                col = col_m1
                              ))
    session$sendCustomMessage("message_motif2",
                              list(
                                seq_m = seq_m2,
                                str_m = str_m2,
                                col = col_m2
                              ))
  }
  
  create_table <- function(selected_motifs) {
    output$table_motifs <-
      renderDT(
        datatable(
          as.data.frame(selected_motifs)[c(2,3,4,5)],
          class = 'table-light',
          options = list(pageLength = 15) ,
          selection = list(mode = 'single', selected = c(1)),
          filter = 'none',
          rownames = FALSE
        ) %>% formatStyle('motif_number', backgroundColor = styleEqual(
          unique(selected_motifs$motif_number), colors[1:nrow(selected_motifs)]
        ))
      )
  }
  
  changeColors <- function(string_col) {
    list_col <- as.list(strsplit(string_col, " ")[[1]])
    string_new_colors <- ""
    j <- 1
    for (i in list_col) {
      if (i != 0) {
        new  <- paste(j, ":", colors[as.numeric(i)], sep = "")
        string_new_colors <- paste(string_new_colors, new)
      }
      j = j + 1
    }
    return(string_new_colors)
  }
  
  draw_overview_structure <- function(str1_, str2_, col, seq)
  {
    #FORNA view
    new_colors <- changeColors(col)
    session$sendCustomMessage(
      "mymessage",
      list(
        sequence = seq,
        str1 = str1_,
        str2 = str2_,
        size = input$char_num,
        color = new_colors
      )
    )
    session$sendCustomMessage(
      "mymessage2",
      list(
        sequence = seq,
        str1 = str1_,
        str2 = str2_,
        size = input$char_num,
        color = new_colors
      )
    )
  }
  
  traceback()
})
