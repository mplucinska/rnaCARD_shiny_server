library(ggplot2)
library(ggthemes)
library(shiny)
library(plotly)
#library(R4RNA)
library(DT)

shinyServer(function(input, output, session) {
  observeEvent(input$submit, {
    output$sequence <- renderText({
      input$sequence
    })
    
    system(
      'python rnaCARD.py -i www/working_dir/structures_2mfe_Fold_homosap_100_500.txt --match --os'
    )
    
    
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
    View(card_motifs)
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
      create_table(selected_motifs)
      observeEvent(input$table_motifs_rows_selected, {
        print(input$table_motifs_rows_selected)
        if (!is.null(input$table_motifs_rows_selected)) {
          draw_motif(selected_motifs[input$table_motifs_rows_selected, ], selected_transcript)
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
  
  observeEvent(input$submit_norm, {
    if (!is.null(input$file2)) {
      df <- read.csv(input$file2$datapath,
                     sep = "\t",
                     header = F)
      download_all(df)
      state = reactiveValues(choice = unique(df$V1))
      output$selectID <- renderUI({
        selectInput("selected", strong(h5("Select transcript ID:")), state$choice  , selected = 1)
      })
    } else{
      dfp <-
        as.data.frame(matrix(unlist(strsplit(
          input$counts, "\n"
        )[[1]])))
      df <-
        data.frame(do.call('rbind', strsplit(as.character(dfp$V1), '\t', fixed =
                                               TRUE)))
      state = reactiveValues(choice = unique(df$X1))
      download_all(df)
      output$selectID <- renderUI({
        selectInput("selected", strong(h5("Select transcript ID:")), state$choice  , selected = 1)
      })
    }
    
    observeEvent(input$selected, {
      state$val <- input$selected
      data_selected_1 <- subset(df, df$V1 == state$val)
      #print(as.numeric(input$ID_col))
      data_selected <-
        cbind(data_selected_1[input$ID_col],
              data_selected_1[input$position_col],
              data_selected_1[input$control_col],
              data_selected_1[input$treated_col])
      #if (sum(data_selected$V3) > 0){
      res <- rnaPRE_results(data_selected)
      normalized_all_data <- res[1]
      table_res <- res[2]
      draw_plots(normalized_all_data, table_res)
      #}
    })
  })
  
  download_all_button <- function() {
    output$download <- downloadHandler(
      filename <- function() {
        paste("aaaa", "csv", sep = ".")
      },
      content = function(file) {
        file.copy("www/working_dir/output_all_rnaPRE.txt", file)
      }
    )
  }
  
  download_all <- function(df) {
    #promise for download all data
    d_all <- future({
      write.table(
        df,
        "www/working_dir/input_all_rnaPRE",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
      )
      system(
        'python new_normalization.py -i www/working_dir/input_all_rnaPRE -o www/working_dir/output_all_rnaPRE.txt '
      )
      normalized_all_data <-
        read_delim(
          "www/working_dir/output_all_rnaPRE.txt",
          "\t",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
    }) %plan% multiprocess
    
    observeEvent(input$calculate, {
      disable("calculate")
      while (!resolved(d_all)) {
        Sys.sleep(1)
      }
      download_all_button()
      output$done <- renderText(1)
    })
  }
  
  observeEvent(input$load_example, {
    output$done <- renderText(0)
    df <- read.csv("test_multiple" , sep = "\t", header = F)
    df <- df[which(df$V1 %in% c('RDN18-1', 'RDN25-1')), ]
    
    
    download_all(df)
    
    state = reactiveValues(choice = unique(df$V1))
    output$selectID <- renderUI({
      selectInput("selected", h6(strong("Select transcript ID:")), state$choice, selected = 1)
    })
    
    observeEvent(input$selected, {
      state$val <- input$selected
      data_selected_1 <- subset(df, df$V1 == state$val)
      data_selected <-
        cbind(data_selected_1[input$ID_col],
              data_selected_1[input$position_col],
              data_selected_1[input$control_col],
              data_selected_1[input$treated_col])
      res <- rnaPRE_results(data_selected)
      normalized_all_data <- res[1]
      table_res <- res[2]
      draw_plots(normalized_all_data, table_res)
      #}
      
    })
  })
  
  observeEvent(input$new_analysis, {
    session$reload()
  })
  
  
  rnaPRE_results <- function(df) {
    write.table(
      df,
      "www/working_dir/input_file_rnaPRE",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
    system(
      'python new_normalization.py -i www/working_dir/input_file_rnaPRE -o www/working_dir/output_file_rnaPRE.txt '
    )
    normalized_all_data <-
      read_delim(
        "www/working_dir/output_file_rnaPRE.txt",
        "\t",
        escape_double = FALSE,
        col_names = FALSE,
        trim_ws = TRUE
      )
    system('rm www/working_dir/output_file_rnaPRE.txt')
    maxy <-
      read_delim(
        "www/working_dir/maxy_output_file_rnaPRE.txt",
        "\t",
        escape_double = FALSE,
        col_names = FALSE,
        trim_ws = TRUE
      )
    output$maxy2 <- renderText(as.numeric(maxy[1, 1]))
    system('rm www/working_dir/maxy_output_file_rnaPRE.txt')
    output$maxy <- renderText(as.numeric(maxy[1, 1]))
    if (maxy == '0') {
      createAlert(
        session,
        "alert",
        "exampleAlert",
        title = "Oops",
        content = "Choose another transcript. Not enough positions with stops to calculate reactivity.",
        append = FALSE,
        style = 'danger'
      )
    }
    print(as.numeric(maxy[1, 1]))
    normalized_all_data$X7[normalized_all_data$X8 == "F"] <- 0
    normalized_all_data$colour[normalized_all_data$X7 < 1.2] <-
      "0.9 - 1.2"
    normalized_all_data$colour[normalized_all_data$X7 < 0.9] <-
      "0.6 - 0.9"
    normalized_all_data$colour[normalized_all_data$X7 < 0.6] <-
      "< 0.6"
    normalized_all_data$colour[normalized_all_data$X7 > 1.2] <-
      "> 1.2"
    write.csv(normalized_all_data[, c(1, 2, 5, 7)],
              "www/working_dir/output_rnaPRE.csv",
              row.names = FALSE)
    table_res <-
      as.data.frame(
        cbind(
          df,
          normalized_all_data$X5 ,
          normalized_all_data$X6 ,
          normalized_all_data$X8
        )
      )
    colnames(table_res) <-
      c(
        "ID",
        "position",
        "counts in control",
        "counts in modified",
        "normalized count in control" ,
        "reactivity",
        "passed filter"
      )
    return(list(normalized_all_data, table_res))
  }
  
  draw_plots <- function(ndata, table_res) {
    ndata <- as.data.frame(ndata)
    #download selscted transcript
    print(unique(ndata$X1))
    output$downloadData <- downloadHandler(filename <- function() {
      paste(unique(ndata$X1), "csv", sep = ".")
    },
    content <- function(file) {
      file.copy("www/working_dir/output_rnaPRE.csv", file)
    })
    
    output$plot_scatter1 <- renderPlotly({
      p <- plot_ly(ndata, x = ~ X3, y = ~ X4)
      ggplotly(p) %>% config(displayModeBar = T,
                             modeBarButtonsToRemove = list_buttons) %>% layout(
                               margin = m,
                               yaxis = list(title = "stops control"),
                               xaxis = list(title = "stops treated", range = c(0, max(
                                 ndata$X3, ndata$X4
                               ))),
                               xaxis = list(title = "")
                             )
    })
    
    output$plot_scatter2 <- renderPlotly({
      p <- plot_ly(ndata,
                   x = ~ X3,
                   y = ~ X5,
                   color = ~ X8)
      ggplotly(p) %>% config(displayModeBar = T,
                             modeBarButtonsToRemove = list_buttons) %>% layout(
                               showlegend = FALSE,
                               margin = m,
                               yaxis = list(title = "normalized stops control"),
                               xaxis = list(title = "stops treated", range = c(0, max(
                                 ndata$X3, ndata$X5
                               )))
                             )
    })
    
    output$plot_histogram <- renderPlotly({
      p <- plot_ly(ndata, x = ~ log)
      ggplotly(p) %>% config(displayModeBar = T,
                             modeBarButtonsToRemove = list_buttons) %>% layout(
                               showlegend = FALSE,
                               margin = m,
                               xaxis = list(title = "log2(FC)")
                             )
    })
    
    output$plot_histogram2 <- renderPlotly({
      p <- plot_ly(ndata, x = ~ X7)
      ggplotly(p) %>% config(displayModeBar = T,
                             modeBarButtonsToRemove = list_buttons) %>% layout(
                               showlegend = FALSE,
                               margin = m,
                               xaxis = list(title = "log2(FC)")
                             )
    })
    
    output$table_results <-
      renderDT(
        as.data.frame(table_res)[-1],
        filter = list(position = 'bottom', clear = FALSE),
        options = list(pageLength = 15),
        rownames = FALSE
      ) #%>% formatStyle(1:5, 'text-align' = 'left')
    
    output$plot <- renderPlotly({
      colour_map = c(
        `0.9 - 1.2` = "#AD4C41",
        `0.6 - 0.9` = "#F2C05C",
        `< 0.6 ` = "#C1BAB6",
        `> 1.2` = "#D0864E"
      )
      #p <- ggplot(normalized_all_data, aes(X2, as.character(X7), fill = colour)) + geom_bar()
      p <-
        plot_ly(
          ndata,
          x = ~ X2,
          y = ~ X7,
          color = ~ colour,
          colors = c("#C1BAB6", "#AD4C41", "#D0864E", "#F2C05C"),
          type = 'bar'
        ) #+ scale_fill_manual(breaks = c("0.9 - 1.2", "0.6 - 0.9", "< 0.6", "> 1.2"), values=c("#AD4C41", "#F2C05C", "#C1BAB6", "#D0864E"))
      ggplotly(p) %>% config(displayModeBar = T,
                             modeBarButtonsToRemove = list_buttons) %>% layout(
                               yaxis = list(fixedrange = T, title = "reactivity"),
                               legend = list(orientation = 'h', y = -0.2),
                               xaxis = list(title = "")
                             )
    })
    
    m <- list(l = 50,
              r = 50,
              b = 50,
              t = 50)
  }
})
