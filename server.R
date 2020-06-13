library(ggplot2)
library(ggthemes)
library(shiny)
library(plotly)
library(DT)

shinyServer(function(input, output, session) {
  
  observeEvent(c(input$submit, input$load_example), {
    results <- parseResults(example = TRUE)
    if (input$submit){
      if (!is.null(input$file1)) {
        df <- read.csv(input$file1$datapath, sep = "\t", header = F)
        write.table(df, "www/working_dir/input_file_rnaCARD.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
        system(
          'python rnaCARD.py -i www/working_dir/input_file_rnaCARD.txt --mismatch --os'
        )
        system(
          'python rnaCARD.py -i www/working_dir/input_file_rnaCARD.txt --match --mismatch --os'
        )
      } else if (input$sequence != "") {
        write(">transcript1", file = "card_input_file.txt")
        write(input$sequence, file = "card_input_file.txt", append = TRUE)
        
        write("seq", file = "card_input_file.txt")
        write(paste("seq", input$sequence, sep = "\t") ,
              file = "card_input_file.txt",
              append = TRUE)
        
        write("str1", file = "card_input_file.txt")
        write(paste("str1", input$str1, sep = "\t") ,
              file = "card_input_file.txt",
              append = TRUE)
        
        write("str2", file = "card_input_file.txt")
        write(paste("str2", input$str2, sep = "\t") ,
              file = "card_input_file.txt",
              append = TRUE)
        
        system(
          'python rnaCARD.py -i www/working_dir/card_input_file.txt --mismatch --os'
        )
        system(
          'python rnaCARD.py -i www/working_dir/card_input_file.txt --match --os'
        )
      }
      results <- parseResults()
      
    }
    
    if (input$load_example) {
      results <- parseResults(example = TRUE)
    }
    
    card_out <- results$card_out
    card_motifs <- results$card_motifs
    card_mismatch_out <- results$card_mismatch_out
    card_mismatch_motifs <- results$card_mismatch_motifs
    
    state = reactiveValues(choice = unique(card_out$X1))
    
    output$selectID_card <- renderUI({
      selectInput("selected_tr", strong(h5("Select transcript ID:")), state$choice, selected = 1)
    })
    
    observeEvent(input$mode, {
      # assign dataframes
      if (input$mode == 'similar motifs') {
        card_out <- card_out
        card_motifs <- card_motifs
      }
      else {
        card_out <- card_mismatch_out
        card_motifs <- card_mismatch_motifs
      }
      
      observeEvent(input$selected_tr, {
        print(input$selected_tr)
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
        
        create_table(selected_motifs)
        
        if (input$mode == 'similar motifs') {
          observeEvent(input$table_motifs_rows_selected, {
            if (!is.null(input$table_motifs_rows_selected)) {
              draw_motif(selected_motifs[input$table_motifs_rows_selected,], selected_transcript)
            }
          })
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
    print(motif)
    print(transcript)
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
      print(col_m1)
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
    print(str_m1)
    print(seq_m1)
    print(col_m1)
    
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
  
  #icons to remove from plotly modebar
  list_buttons <-
    list(
      "autoScale2d" ,
      'Collaborate',
      'toggleSpikelines',
      'lasso2d',
      'pan2d',
      'sendDataToCloud',
      'select2d',
      'zoomIn2d',
      'zoomOut2d',
      'hoverClosestCartesian',
      'hoverCompareCartesian'
    )

  parseResults <- function(example = FALSE){
    # parse and display CARD results (column names etc.)
    if (example == FALSE){
      card_out <-
        read_delim(
          "www/working_dir/matched_whole_transcript.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
      
      card_motifs <-
        read_delim(
          "www/working_dir/matched_motifs_out.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
      
      card_mismatch_out <-
        read_delim(
          "www/working_dir/mismatched_whole_transcripts.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
      
      card_mismatch_motifs <-
        read_delim(
          "www/working_dir/mismatched_motifs_out.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
    } else {
      card_out <-
        read_delim(
          "www/working_dir/example_matched_whole_transcripts.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
      
      card_motifs <-
        read_delim(
          "www/working_dir/example_matched_motifs_out.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
      
      card_mismatch_out <-
        read_delim(
          "www/working_dir/example_mismatched_whole_transcripts.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
      
      card_mismatch_motifs <-
        read_delim(
          "www/working_dir/example_mismatched_motifs_out.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
    }
    
    
    colnames(card_mismatch_motifs) <-
      c(
        "ID",
        "motif_number",
        "start position(s)",
        "end position(s)" ,
        "structure 1",
        "structure 2"
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
    
    results <-
      list(
        card_out = card_out,
        card_motifs = card_motifs,
        card_mismatch_out = card_mismatch_out,
        card_mismatch_motifs = card_mismatch_motifs
      )
    return(results)
  }

  })
