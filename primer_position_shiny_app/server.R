### Load in required packages
library(shiny)
library(rtracklayer)
library(GenomicFeatures)
library(AnnotationDbi)
library(biovizBase)
library(ggbio)
options(shiny.maxRequestSize=500*1024^2)
speciesMap <- available.species()
### Empty variables for later
variables = reactiveValues(q = "", p = "no", txdb = "", c = "" , grl.primer.temp = NULL, chromes = "")
names <- c("Chr","strand","start","start2","end","end2","gene_id","primername") # For csv checking
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
strand.color <- c("magenta", "black", "green") # For Plotting
shinyServer(function(input, output) {
  observeEvent(input$txdb2, {
    if(input$existing == "A.tha"){
      variables$txdb <- loadDb("data/A.tha.TxDb")
      variables$c <- "Chr"
    }
    if(input$existing == "B.nap"){
      variables$txdb <- loadDb("data/B.nap.TxDb")
      variables$c <- "chr"
    }
    if(input$existing == "H.sap"){
      variables$txdb <- loadDb("data/H.sap.TxDb")
      variables$c <- ""
    }
    if(input$existing == "S.lyc"){
      variables$txdb <- loadDb("data/S.lyc.TxDb")
      variables$c <- "SL2.50ch"
    }
    variables$p <- "yes"
    variables$q <- "Database loaded"
  })
  observeEvent(input$txdb1, {
    variables$grl.primer.temp <- NULL
    variables$txdb <- ""
    variables$p <- "no"
    ### Attempt to make TxDb
    if(!is.null(input$file1)){
      variables$txdb <- try(makeTxDbFromGFF(
        file = input$file1$datapath,
        organism = input$organism,
        format = "gff3"))
    } else {
      variables$q <- "Please provide valid GFF file and organism name."
    }
    # Notify user TxDb has been created and sets "c", used later for compatibility between csv and gff
    if(class(variables$txdb) == "TxDb"){
      variables$q <- "Database has been created. Please select primer csv if you have not done so yet."
      variables$p <- "yes"
      variables$c <- input$prefix
    } else {
      variables$q <- "Please provide valid GFF file and organism name."
      variables$txdb <- ""
    }
  })
  observeEvent(input$image, {
    variables$grl.primer.temp <- NULL
    if(variables$p == "yes"){
      # Resets notification message and sets up primer pair variable
      variables$q <- ""
      # Checks for the presence of csv and the neccessary columns in csv
      if(!is.null(input$file2)){
        primer.info <- try(read.csv(input$file2$datapath, colClasses = c("Chr"="character")))
        if(sum(names %in% colnames(primer.info)) != 8){
          variables$q <- "Please provide valid primer csv file. Do you have the columns (primername,Chr,start,end,start2,end2,gene_id,strand)?"
          rm(primer.info)
          return()
        }
        # Check for presence of correct strand info
        if((sum((as.character(primer.info$strand) == "+") | (as.character(primer.info$strand) == "*") | (as.character(primer.info$strand) == "-")) != nrow(primer.info)) 
           | is.na(sum((as.character(primer.info$strand) == "+") | (as.character(primer.info$strand) == "*") | (as.character(primer.info$strand) == "-")))){
          variables$q <- "Strand values need to be one of the following (+ , - , *)"
          return()
        }
        i <- input$pair1
        if(i!=1){
          i <- i*2-1
        }
      } else {
        variables$q <- "Please provide a primer csv file."
        return()
      }
      # Check that primers selected exist
      if(i <= nrow(primer.info) & i>=1){
        primers <- primer.info
        gene_id <-unique(primer.info[i:i+1,"gene_id"])
        # Check that there was 2 primers provide
        if(is.na(gene_id)){
          variables$q <- "Only one primer was provided in the pair."
          return()
        }
        # Check that primers are for the same gene
        if(primers[i,"gene_id"] != primers[i+1,"gene_id"]) {
          variables$q <- "Target gene of primers do no match."
          return()
        }
        # Check primers are on different strands
        if(primers[i,"strand"] == primers[i+1,"strand"]){
          variables$q <- "Both primers are on the same strand."
          return()
        }
        if(((as.character(primers[i,"primername"]) == "") == TRUE) | ((as.character(primers[i+1,"primername"]) == "") == TRUE)){
          variables$q <- "Name of one or both primers is missing."
          return()
        }
      } else {
        variables$q <- "Primer set does not exist."
        return()
      }
      # making primer GRange object 
      primers$primername<-as.character(primers$primername)
      primers$start2 <- as.numeric(as.character(primers$start2))
      primers$end2 <- as.numeric(as.character(primers$end2))
      # add primer info for AT1G01060 (fake primers)
      primer1 <- try(GRanges(paste(variables$c,primers[i,"Chr"],sep=""), IRanges(primers[i,"start"], primers[i,"end"]),strand=as.character(primers[i,"strand"]),tx_id="",tx_name=primers[i,"primername"],gene_id=primers[i,"gene_id"],model="utr"))
      primer2 <- try(GRanges(paste(variables$c,primers[i+1,"Chr"],sep=""), IRanges(primers[i+1,"start"], primers[i+1,"end"]),strand=as.character(primers[i+1,"strand"]),tx_id="",tx_name=primers[i+1,"primername"],gene_id=primers[i+1,"gene_id"],model="utr"))
      # if primer is splitted
      if((class(primer1) == "try-error") | (class(primer2) == "try-error")){
        variables$q <- "start locations need to be less than end locations."
        return()
      }
      if(!is.na(primers[i,"start2"]))  {
        print(primers[i,])
        primer1b <- try(GRanges(paste(variables$c,primers[i,"Chr"],sep=""), IRanges(as.numeric(as.character(primers[i,"start2"])), as.numeric(as.character(primers[i,"end2"]))),strand=as.character(primers[i,"strand"]),tx_id="",tx_name=primers[i,"primername"],gene_id=primers[i,"gene_id"],model="utr"))
        primer1 <- try(c(primer1,primer1b))
        if(class(primer1b) == "try-error"){
          variables$q <- "start locations need to be less than end locations."
          return()
        }
      }
      if(!is.na(primers[i+1,"start2"]))  {
        print(primers[i+1,])
        primer2b <- try(GRanges(paste(variables$c,primers[i+1,"Chr"],sep=""), IRanges(as.numeric(as.character(primers[i+1,"start2"])), as.numeric(as.character(primers[i+1,"end2"]))),strand=as.character(primers[i+1,"strand"]),tx_id="",tx_name=primers[i+1,"primername"],gene_id=primers[i+1,"gene_id"],model="utr"))
        primer2<- try(c(primer2,primer2b))
        if(class(primer2b) == "try-error"){
          variables$q <- "start locations need to be less than end locations."
          return()
        }
      }
      primers.temp<-c(primer1,primer2) #GRange object
      #print(strand(primers.temp))
      print(primers.temp)
      if("ChrNA" %in% seqnames(primers.temp)){
        variables$q <- "Primers need to have valid Chr value"
        return()
      }
      transcripts.selected.temp <- transcripts(variables$txdb,filter=list(gene_id=gene_id))
      if(length(transcripts.selected.temp) == 0){
        transcripts.selected.temp <- transcripts(variables$txdb,filter=list(tx_name=gene_id))
        genes <- as.character(gene_id)
        genes <- rep(genes, length(seqnames(transcripts.selected.temp)))
        holder <- unlist(transcripts.selected.temp)
        holder$gene_id <- genes
        transcripts.selected.temp <- relist(holder,transcripts.selected.temp)
        if(length(transcripts.selected.temp) == 0){
          variables$q <- "gene_id's in primer csv do not match gene_id's in database."
          return()
        }
      }
      try(if(as.character(unique(seqnames(primers.temp))) != as.character(unique(seqnames(transcripts.selected.temp)))){
        variables$q <- "Chromosome of gene and primer are not the same"
        return()
      })
      gr.temp <- try(crunch(variables$txdb,which=transcripts.selected.temp))
        if(class(gr.temp) == "try-error"){
          variables$q <- "Please provide valid primer csv file. Do you have the columns (primername,Chr,start,end,start2,end2,gene_id,strand)?"
          return()
        }
      ## change column to 'model'
      colnames(values(gr.temp))[4] <- "model"
      # into GRL
      gr.primer.temp<-c(gr.temp,primers.temp)
      # into GRL
      variables$grl.primer.temp <- split(gr.primer.temp, gr.primer.temp$tx_name)
      # print 
      print(variables$grl.primer.temp)
      variables$q <- paste("Plot for", gene_id)
    } else {
      variables$q <- "Please make a TxDb first"
    }
  })
  observeEvent(input$image2, {
    variables$grl.primer.temp <- NULL
    if(variables$p == "yes"){
      # Resets notification message and sets up primer pair variable
      variables$q <- ""
      # Checks for the presence of csv and the neccessary columns in csv
      if(!is.null(input$file3)){
        primer.info <- try(read.csv(input$file3$datapath, colClasses = c("Chr"="character")))
        if(sum(names %in% colnames(primer.info)) != 8){
          variables$q <- "Please provide valid primer csv file. Do you have the columns (primername,Chr,start,end,start2,end2,gene_id,strand)?"
          rm(primer.info)
          return()
        }
        # Check for presence of correct strand info
        if((sum((as.character(primer.info$strand) == "+") | (as.character(primer.info$strand) == "*") | (as.character(primer.info$strand) == "-")) != nrow(primer.info)) 
           | is.na(sum((as.character(primer.info$strand) == "+") | (as.character(primer.info$strand) == "*") | (as.character(primer.info$strand) == "-")))){
          variables$q <- "Strand values need to be one of the following (+ , - , *)"
          return()
        }
        i <- input$pair2
        if(i!=1){
          i <- i*2-1
        }
      } else {
        variables$q <- "Please provide a primer csv file."
        return()
      }
      # Check that primers selected exist
      if(i <= nrow(primer.info) & i>=1){
        primers <- primer.info
        gene_id <-unique(primer.info[i:i+1,"gene_id"])
        # Check that there was 2 primers provide
        if(is.na(gene_id)){
          variables$q <- "Only one primer was provided in the pair."
          return()
        }
        # Check that primers are for the same gene
        if(primers[i,"gene_id"] != primers[i+1,"gene_id"]) {
          variables$q <- "Target gene of primers do no match."
          return()
        }
        # Check primers are on different strands
        if(primers[i,"strand"] == primers[i+1,"strand"]){
          variables$q <- "Both primers are on the same strand."
          return()
        }
        if(((as.character(primers[i,"primername"]) == "") == TRUE) | ((as.character(primers[i+1,"primername"]) == "") == TRUE)){
          variables$q <- "Name of one or both primers is missing."
          return()
        }
      } else {
        variables$q <- "Primer set does not exist."
        return()
      }
      # making primer GRange object 
      primers$primername<-as.character(primers$primername)
      primers$start2 <- as.numeric(as.character(primers$start2))
      primers$end2 <- as.numeric(as.character(primers$end2))
      # add primer info for AT1G01060 (fake primers)
      primer1 <- try(GRanges(paste(variables$c,primers[i,"Chr"],sep=""), IRanges(primers[i,"start"], primers[i,"end"]),strand=as.character(primers[i,"strand"]),tx_id="",tx_name=primers[i,"primername"],gene_id=primers[i,"gene_id"],model="utr"))
      primer2 <- try(GRanges(paste(variables$c,primers[i+1,"Chr"],sep=""), IRanges(primers[i+1,"start"], primers[i+1,"end"]),strand=as.character(primers[i+1,"strand"]),tx_id="",tx_name=primers[i+1,"primername"],gene_id=primers[i+1,"gene_id"],model="utr"))
      # if primer is splitted
      if((class(primer1) == "try-error") | (class(primer2) == "try-error")){
        variables$q <- "start locations need to be less than end locations."
        return()
      }
      if(!is.na(primers[i,"start2"]))  {
        print(primers[i,])
        primer1b <- try(GRanges(paste(variables$c,primers[i,"Chr"],sep=""), IRanges(as.numeric(as.character(primers[i,"start2"])), as.numeric(as.character(primers[i,"end2"]))),strand=as.character(primers[i,"strand"]),tx_id="",tx_name=primers[i,"primername"],gene_id=primers[i,"gene_id"],model="utr"))
        primer1 <- try(c(primer1,primer1b))
        if(class(primer1b) == "try-error"){
          variables$q <- "start locations need to be less than end locations."
          return()
        }
      }
      if(!is.na(primers[i+1,"start2"]))  {
        print(primers[i+1,])
        primer2b <- try(GRanges(paste(variables$c,primers[i+1,"Chr"],sep=""), IRanges(as.numeric(as.character(primers[i+1,"start2"])), as.numeric(as.character(primers[i+1,"end2"]))),strand=as.character(primers[i+1,"strand"]),tx_id="",tx_name=primers[i+1,"primername"],gene_id=primers[i+1,"gene_id"],model="utr"))
        primer2<- try(c(primer2,primer2b))
        if(class(primer2b) == "try-error"){
          variables$q <- "start locations need to be less than end locations."
          return()
        }
      }
      primers.temp<-c(primer1,primer2) #GRange object
      #print(strand(primers.temp))
      print(primers.temp)
      if("ChrNA" %in% seqnames(primers.temp)){
        variables$q <- "Primers need to have valid Chr value"
        return()
      }
      transcripts.selected.temp <- transcripts(variables$txdb,filter=list(gene_id=gene_id))
      if(length(transcripts.selected.temp) == 0){
        transcripts.selected.temp <- transcripts(variables$txdb,filter=list(tx_name=gene_id))
        genes <- as.character(gene_id)
        genes <- rep(genes, length(seqnames(transcripts.selected.temp)))
        holder <- unlist(transcripts.selected.temp)
        holder$gene_id <- genes
        transcripts.selected.temp <- relist(holder,transcripts.selected.temp)
        if(length(transcripts.selected.temp) == 0){
          variables$q <- "gene_id's in primer csv do not match gene_id's in database."
          return()
        }
      }
      try(if(as.character(unique(seqnames(primers.temp))) != as.character(unique(seqnames(transcripts.selected.temp)))){
        variables$q <- "Chromosome of gene and primer are not the same"
        return()
      })
      gr.temp <- try(crunch(variables$txdb,which=transcripts.selected.temp))
      if(class(gr.temp) == "try-error"){
        variables$q <- "Please provide valid primer csv file. Do you have the columns (primername,Chr,start,end,start2,end2,gene_id,strand)?"
        return()
      }
      ## change column to 'model'
      colnames(values(gr.temp))[4] <- "model"
      # into GRL
      gr.primer.temp<-c(gr.temp,primers.temp)
      # into GRL
      variables$grl.primer.temp <- split(gr.primer.temp, gr.primer.temp$tx_name)
      # print 
      print(variables$grl.primer.temp)
      variables$q <- paste("Plot for", gene_id)
    } else {
      variables$q <- "Please load a TxDb first"
    }
  })
  output$plot <- renderPlot({
    try(autoplot(variables$grl.primer.temp,aes(type=model,fill=strand)) + 
          scale_fill_manual(values= strand.color) + 
          theme(axis.text.y=element_blank(),panel.grid.major = element_blank(), 
                panel.grid.minor =  element_blank(), 
                panel.background = element_blank(), 
                axis.line.y = element_line(colour = "white"),
                axis.ticks.y=element_line(color="white")))
  })
  output$text <- renderText({
    print(variables$q)
  })

})
