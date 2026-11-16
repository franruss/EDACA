
source("definitions.R")


library(BiocManager)
options(repos = BiocManager::repositories())

shinyServer(function(input, output) {
  
  #Training and Test
    results2 <-  eventReactive(input$run, {
    comparison=NULL
    in_method <- input$method
    in_normtrans <- "tmm-logcpm"
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    in_cross_val <- input$cross_val
    in_cross_val = as.numeric(in_cross_val)
    in_repeats_num <- input$repeats_num
    in_repeats_num = as.numeric(in_repeats_num)
    in_tuneLength <- input$tuneLength
    in_tuneLength = as.numeric(in_tuneLength)
    # in_method <- "svmRadial"
      # by clicking the 'Run' button the code below is executed!
      # setwd("C:/Users/f.russo.ENGIBBC/Desktop/EDACA_paper/")
      cervical  = NULL
      templates = NULL
      if (isTRUE(input$use_example) && is.null(input$file3) && is.null(input$file4)) {
        cervical  = read.delim("data_mrna_seq_read_counts.txt", header = TRUE, sep = "\t", row.names = 1)
        templates = read.delim("design_data_mrna_seq_read_counts.txt", header = FALSE, sep = "\t", row.names = 1)
      }else{ #use an external file
        inFile3 <- input$file3
        inFile4 <- input$file4
        
        path <- inFile3$datapath
        is_zip <- tryCatch(!is.null(unzip(path, list = TRUE)),error = function(e) FALSE)
        if (is_zip) {
          exdir <- file.path(tempdir(), paste0("unzipped_", Sys.getpid()))
          dir.create(exdir, showWarnings = FALSE)
          unzip(path, exdir = exdir)
          txt_file <- list.files(exdir,pattern = "\\.txt$",full.names = TRUE)
          validate(need(length(txt_file) > 0, "ZIP contains no .txt file"))
          print("txt_file[1]")
          print(txt_file[1])
          cervical <- read.table(file = txt_file[1], header = TRUE,sep = "\t",row.names=1)
        } else {cervical = read.table(file = inFile3$datapath, header = TRUE,sep = "\t",row.names=1)}
        
        #cervical = read.table(file = inFile3$datapath, header = TRUE,sep = "\t",row.names=1)
        templates = read.table(file= inFile4$datapath, header = FALSE,sep = "\t",row.names=1)
      }
      ### Se ci sono dei pazienti con valori mancanti o con dei NA in (templates)[2,] cancellarli  da templates e poi da cervical
      valid_ids = !is.na( as.character(templates[2,]))
      templates = templates[,valid_ids]
      cervical = cervical[,gsub("-", ".",((templates)[1,]))]
      # cervical = read.table(file ="data", header = TRUE,sep = "\t",row.names=1)
      # templates = read.table(file ="design_data_mrna_seq_read_counts.txt", header = FALSE,sep = "\t",row.names=1)
      # templates = read.table(file ="TCGA-CHOL_design_gender.txt", header = FALSE,sep = "\t",row.names=1)
      # in_vars_num = 4
      # QUI DEVE ESSE L'UTENTE A SCEGLIERE QUALE SAMPLE e' A E QUALE e' B 
      # Anche qui deve poter decidere l'utente quali sono i pazienti del training e quali del test
      class <- data.frame(condition = factor( as.character(templates[2,]) ))
      set.seed(2128)
      # We do not perform a differential expression analysis to select differentially
      # expressed genes. However, in practice, DE analysis might be performed before
      # fitting classifiers. Here, we selected top  "in_vars_num"  features having the highest
      # gene-wise variances in order to decrease computational cost.
      # Invece di far cosi bisogna prendere le features che fornisce l'utente.
      # Si potrebbe dare questa qui sotto come opzione
      if ( input$selector == 'variance' ){
         vars <- sort(apply(cervical, 1, var, na.rm = TRUE), decreasing = TRUE) #highest gene-wise variances
         data <- cervical[names(vars)[1: in_vars_num], ]
      }else{
        gene_list = as.character(input$predictor_list)
        # gene_list = as.character("MGAT3,MGAT4C,A4GALT,A3GALT2,ABO")
        genes = strsplit(gene_list,",")[[1]]
        #entrez_ids <- mapIds(org.Hs.eg.db, keys = as.character(genes), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
        #entrez_ids <- na.omit(entrez_ids)
        #genes_list <- as.character(entrez_ids)
        ensembl_ids <- mapIds(org.Hs.eg.db,keys = as.character(genes), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
        genes_list <- as.character(ensembl_ids)
        canonical_rows <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", rownames(cervical))
        cervical_clean <- cervical[canonical_rows, ] # remove _PAR_|_HAP|_PATCH
        rownames(cervical_clean) <- sub("\\..*$", "", rownames(cervical_clean)) #remove .35
        data <- cervical[genes_list, ]
      }
      # data <- cervical[as.character(templates[,1]), ]
      # nTest <- ceiling(ncol(data) * 0.3) #30% for test
      # ind <- sample(ncol(data), nTest, FALSE) #pick random numbers to select the train samples
      ind <- which( templates[3,]=="F")  #pick those samples from templates to select the test and train samples
      data.train <- as.matrix(data[ ,-ind] + 1)
      data.test <- as.matrix(data[ ,ind] + 1)
      classtr <- data.frame(condition = class[-ind, ])
      classts <- data.frame(condition = class[ind, ])
      data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,design = formula(~condition))
      data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,design = formula(~condition))
      withProgress(message = 'Running, Please Wait!', detail = 'This may take a while...',  value = 0, {
        # Number of times we'll go through the loop
        n <- 100
        # Define control lists.
        ctrl.continuous <- trainControl   (method = "repeatedcv", number = in_cross_val, repeats = in_repeats_num)
        ctrl.discrete   <- discreteControl(method = "repeatedcv", number = in_cross_val, repeats = in_repeats_num,tuneLength = in_tuneLength)
        ctrl.voom       <- voomControl    (method = "repeatedcv", number = in_cross_val, repeats = in_repeats_num,tuneLength = in_tuneLength)
        discrete_classifiers = c("PLDA", "PLDA2", "NBLDA")
        voom_classifiers = c("voomDLDA", "voomDQDA", "voomNSC")
        fit=NULL
        if (is.element(in_method, discrete_classifiers)){
            print("discrete classifiers")
            fit <- classify(data = data.trainS4, method = in_method, normalize = "cpm", ref = templates[2,3], control = ctrl.discrete)
        }else if(is.element(in_method, voom_classifiers)){
            fit <- classify(data = data.trainS4, method = in_method, preProcessing = in_normtrans, ref = templates[2,3],tuneLength = in_tuneLength, control = ctrl.voom)
        }else{
            fit <- classify(data = data.trainS4, method = in_method, preProcessing = in_normtrans, ref = templates[2,3],tuneLength = in_tuneLength, control = ctrl.continuous)
        }
      #Here we use predict fun for predictions
      pred.svm = predict(fit, data.testS4)
      predicted_classes <- relevel(pred.svm, ref = templates[2,3])
      actual_classes    <- relevel(classts$condition, ref = templates[2,3])
      comparison = rbind(as.character(predicted_classes),as.character(actual_classes))
      comparison = cbind(c("predicted classes","true classes"),comparison)
      colnames(comparison) = c("classes" , colnames(data.test))
      for (i in 1:n) {
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("\n  ",i, "% completed!"))
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.0005)
      }
      })
    comparison_no_first_col = comparison[,-1] #togli la prima colonna perchè sono solo i nomi delle variabili
    trained_fit=NULL
    if(is.element(in_method, voom_classifiers) || is.element(in_method, discrete_classifiers)){
       trained_fit= MLSeq::trained(fit)
    }else{
       trained_fit= MLSeq::trained(fit)[[4]]
    }
    #selected_genes = selectedGenes(fit)
    #print(selected_genes)
    #show(fit)
    #print(trained(fit))
    return(list(comparison = comparison, fit = fit, trained_fit = trained_fit, pred.svm = pred.svm, classts = classts, ref = templates[2,3])) 
  })
    
    output$summaryTrain <- renderPrint({
      req(results2())
      fit = results2()$fit
      trained(fit)   # prints the summary
    })
    
    output$summary <- renderText({
      req(results2())
      comparison = results2()$comparison
      comparison_no_first_col = comparison[,-1] #togli la prima colonna perchè sono solo i nomi delle variabili
      in_method <- input$method
      in_vars_num <- input$vars_num
      in_vars_num = as.numeric(in_vars_num)
      in_cross_val <- input$cross_val
      in_cross_val = as.numeric(in_cross_val)
      in_repeats_num <- input$repeats_num
      in_repeats_num = as.numeric(in_repeats_num)
      in_tuneLength <- input$tuneLength
      in_tuneLength = as.numeric(in_tuneLength)
      predictor=NULL
      if ( input$selector == 'variance' ){
          predictor_info = paste("Number of predictors having the highest variance used: ",in_vars_num, sep="")  
      }else{
          predictor_info = paste("Predictors used: ", as.character(input$predictor_list), sep="")  
      }
      summary_res = paste(paste("You picked the algorithm: ",as.character(in_method), sep=""),
                          paste("Reference class used: ", results2()$ref, sep=""),
                          predictor_info,
                          paste("Cross Validation Fold: ",in_cross_val, sep=""),
                          paste("Number of Repeats: ", in_repeats_num, sep=""),
                          paste("Training Circles: ",in_tuneLength, sep=""),
                          paste("Total Comparisons: ",dim(comparison_no_first_col)[2], sep=""),
                          paste("Right predictions: ", sum(comparison_no_first_col[1,]==comparison_no_first_col[2,]), sep=""),
                          paste("Wrong predictions: ", sum(comparison_no_first_col[1,]!=comparison_no_first_col[2,]), sep=""),
                          paste("Testing Accuracy: ",  round(sum(comparison_no_first_col[1,]==comparison_no_first_col[2,])/dim(comparison_no_first_col)[2],4),sep=""),
                    sep="\n")
      return(summary_res)
    })
  
  result_classification <- reactive({  
    req(results2())
    comparison <- results2()$comparison
    comparison_no_first_col = comparison[,-1] #togli la prima colonna perchè sono solo i nomi delle variabili
    comparison
  })
  
  output$unsupervised_result_classification <- renderDataTable({  
    # C'è un errore! non bisogna farlo per tutte le classi in template. Ma solo per quelle che vengono identificate nel test set del metodo scelto
    inFile3 <- input$file3
    inClusters <- input$clusters
    inInitialConf <- input$initialConf
    mydata = read.table(file = inFile3$datapath, header = TRUE,sep = "\t",row.names=1)
    #mydata = read.table(file = "c", header = TRUE,sep = "\t",row.names=1)
    dim(mydata)
    keep <- rowSums(mydata) >= 10
    mydata <- mydata[keep,]
    dim(mydata)
    #z <- t(mydata)
    #means <- apply(z,2,mean)
    #sds <- apply(z,2,sd)
    #nor <- scale(z,center=means,scale=sds)
    #k2 <- kmeans(nor, centers = inClusters, nstart = inInitialConf)
    #fviz_cluster(k2, data = nor)
    x = mydata
    x = log10(x+1)
    # metodi per clusterizzare: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
    clustering_used_method = "ward.D"
    distance_used_method = "euclidean"
    # metodi per la distanza: "euclidean"
    # ottieni i cluster
    # fit <- hclust(dist(t(x),method="euclidean"),method=used_method)
    # invece di scriverla come il rigo sopra: con due funzioni innestate, e' meglio definire una funzione per volta in modo separato 
    # e poi richiamarle insieme innestandole con la funzione fit
    hclustfunc <- function(x)  hclust(x, method = clustering_used_method) # definisco la funzione hclustfunc da hclust
    distfunc <-   function(x)    dist(x, method = distance_used_method)      # definisco la funzione distfunc da dist
    fit <- hclustfunc(distfunc((t(x)))) #le innesto una dentro l'atra, chiamando prima distfunc
    clusters <- cutree(fit, inClusters) # inClusters e' il numero di gruppi che vuoi ottenere 
    t(clusters)
    #plot(fit)
  })
  
  
  output$download <- downloadHandler(
    filename = function(){"results.txt"}, 
   content = function(fname){
      write.table((rnaseq_analysis(input$file5,input$text1)), fname, col.names = TRUE,row.names = FALSE,quote=FALSE, sep = "\t")
    }
  )
  
  
  output$results5 <- renderDataTable({  
    inFile5 <- input$file5
    question <- input$text1
    if (is.null(inFile5) || is.null(question) ){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
       counts = read.table(inFile5$datapath,row.names=1)
  #    rownames(counts) <- sapply(rownames(counts) , function(x) 
  #      toString(tags$a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",x,"&keywords=",x), x))
  #    )
  #    resulttable = datatable(as.data.frame(counts), escape=FALSE)
      id=rownames(counts)
      (cbind(id,counts))[1:4,]
    }
  })
  
  output$results8 <- renderDataTable({  
    counts = matrix(data=" ",nrow=2,ncol=2)
    colnames(counts) = c(" "," ")
    if (isTRUE(input$use_example) && is.null(input$file3)) {
      counts  = read.delim("data_mrna_seq_read_counts.txt", header = TRUE, sep = "\t", row.names = 1)
    }
    if(isFALSE(input$use_example) && !is.null(input$file3)){
      inFile3 <- input$file3
      
      path <- inFile3$datapath
      is_zip <- tryCatch(!is.null(unzip(path, list = TRUE)),error = function(e) FALSE)
      if (is_zip) {
        exdir <- file.path(tempdir(), paste0("unzipped_", Sys.getpid()))
        dir.create(exdir, showWarnings = FALSE)
        unzip(path, exdir = exdir)
        txt_file <- list.files(exdir,pattern = "\\.txt$",full.names = TRUE)
        validate(need(length(txt_file) > 0, "ZIP contains no .txt file"))
        print("txt_file[1]")
        print(txt_file[1])
        counts <- read.table(file = txt_file[1], header = TRUE,sep = "\t",row.names=1)
      } else {counts = read.table(file = inFile3$datapath, header = TRUE,sep = "\t",row.names=1)}
      
    }
    datatable(counts) 
  })
  
  
  output$results9 <- renderDataTable({  
      counts = matrix(data=" ",nrow=2,ncol=2)
      colnames(counts) = c(" "," ")
      if (isTRUE(input$use_example) && is.null(input$file4)) {
        counts  = read.delim("design_data_mrna_seq_read_counts.txt", header = TRUE, sep = "\t", row.names = 1)
      }
      if(isFALSE(input$use_example) && !is.null(input$file4)){
      inFile4 = input$file4
      counts = read.table(inFile4$datapath,header=TRUE,sep = "\t",row.names=1)
      }
      t(counts)
  })
  
  build_data <- reactive( {
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- as.numeric(input$vars_num2)
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{ 
      cervical = read.table(file = inFile11$datapath, header = TRUE,sep = "\t",row.names=1)
      templates = read.table(file= inFile12$datapath, header = TRUE,sep = "\t",row.names=1)
      # in_vars_num=4
      # layers= c(1,2)
      # templates = read.table(file ="training_test_design_matrix_ermes.txt", header = TRUE,sep = "\t",row.names=1)
      # cervical = read.table(file ="Ion_AmpliSeq_Transcriptome_Expression_ERMES_results.txt", header = TRUE,sep = "\t",row.names=1)
      templates = templates[,colnames(cervical)] #importante per ordinare correttamente 
      class <- as.factor(as.character(templates[1,]))
      # Si potrebbe dare questa qui sotto come opzione
      vars <- sort(apply(cervical, 1, var, na.rm = TRUE), decreasing = TRUE) #highest gene-wise variances
      data <- cervical[names(vars)[1: in_vars_num], ]
      data = t(data)
      class = as.list(class)
      data = data.frame(data,class=unlist(class)) #importante unlist altrimenti neural net non funziona
      data
    }
  })
  
  build_data.train <- reactive( {
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- as.numeric(input$vars_num2)
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    req(build_data())
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{ 
      templates = read.table(file= inFile12$datapath, header = TRUE,sep = "\t",row.names=1)
      cervical = read.table(file = inFile11$datapath, header = TRUE,sep = "\t",row.names=1)
      # in_vars_num=4
      # layers= c(1,2)
      # templates = read.table(file ="training_test_design_matrix_ermes.txt", header = TRUE,sep = "\t",row.names=1)
      # cervical = read.table(file ="Ion_AmpliSeq_Transcriptome_Expression_ERMES_results.txt", header = TRUE,sep = "\t",row.names=1)
      templates = templates[,colnames(cervical)] #importante per ordinare correttamente 
      data <- build_data()
      ind <- which( templates[2,]=="F")  
      data.train <- (data[ -ind,]) #pick those samples from templates to select train samples with "-ind"
      data.train
    }
  })
  
  build_data.test <- reactive( {
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- as.numeric(input$vars_num2)
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    req(build_data())
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{ 
      templates = read.table(file= inFile12$datapath, header = TRUE,sep = "\t",row.names=1)
      cervical = read.table(file = inFile11$datapath, header = TRUE,sep = "\t",row.names=1)
      # in_vars_num=4
      # layers= c(1,2)
      # templates = read.table(file ="training_test_design_matrix_ermes.txt", header = TRUE,sep = "\t",row.names=1)
      # cervical = read.table(file ="Ion_AmpliSeq_Transcriptome_Expression_ERMES_results.txt", header = TRUE,sep = "\t",row.names=1)
      data <- build_data()
      templates = templates[,colnames(cervical)] #importante per ordinare correttamente 
      ind <- which( templates[2,]=="F")  
      data.test <- (data[ ind,]) #pick those samples from templates to select test samples with "ind"
      data.test
    }
  })
  
  build_formula <- reactive( {
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- input$vars_num2
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    req(build_data.train())
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{ 
       data.train <- build_data.train()
       gene_names <- colnames(data.train)[1:(length(colnames(data.train))-1)] #nomi dei geni
       f <- as.formula(paste("class ~", paste(gene_names, collapse = " + ")))
       f
    }
  })
  
  
  build_network <- reactive( {
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- input$vars_num2
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    req(build_formula())
    req(build_data.train())
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{
      f <- build_formula()
      data.train <- build_data.train()
      nn <- neuralnet(f,data = as.data.frame(data.train),hidden =  layers,linear.output = F,lifesign = 'full')
      nn
    }
  })
  
  output$show_network <- renderPlot({  
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- input$vars_num2
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    req(build_network())
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{
      nn <-build_network()
      plotnet(nn,col.hidden = 'darkgreen',     
           col.hidden.synapse = 'darkgreen',
           show.weights = T,
           information = T,
           fill = 'lightblue')
    }
  })
  
  
  output$show_testing_results <- renderDataTable({  
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- input$vars_num2
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    req(build_network())
    req(build_formula())
    req(build_data.test())
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{
      #  pred <- predict(nn, data.test)
      data.test <- build_data.test()
      nn <-build_network()
      pred <- predict(nn, data.test)
      #labels = unique(as.character(data.test$class))
      labels = levels(data.test$class)
      true_values <- (data.test[,dim(data.test)[2]])
      colnames(pred) = labels
      predictions <- (data.frame(max.col(pred))[,1])
      res = rbind(rbind(rownames(data.test),as.character(true_values)),labels[predictions])
      res = t(res)
      colnames(res) = c("ID", "true_values", "predictions")
      table_res = table(res[,2], res[,3])
      table_res = as.data.frame(table_res)
      colnames(table_res) = c("true values", "predictions","frequency")
      #table_res
      t(res)
      
    }
  })
  
  output$show_input <- renderDataTable({  
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- input$vars_num2
    if (is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
      counts = read.table(inFile11$datapath,header=TRUE, row.names=1)
      id=rownames(counts)
      (cbind(id,counts))[1:4,]
    }
  })
  
  
  output$show_design <- renderDataTable({  
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- input$vars_num2
    if (is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
      counts = read.table(inFile12$datapath,header=TRUE,row.names=1)
      id=rownames(counts)
      cbind(id,counts)
    }
  })
  
  output$show_accuracy <- renderDataTable({  
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- input$vars_num2
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    req(build_network())
    req(build_formula())
    req(build_data.test())
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{
      #  pred <- predict(nn, data.test)
      data.test <- build_data.test()
      nn <-build_network()
      pred <- predict(nn, data.test)
      check = as.numeric(data.test$class) == max.col(pred)
      right_predictions = sum(check)
      accuracy = (sum(check)/nrow(data.test))
      
      res = cbind(right_predictions,accuracy)
      res2 =  cbind(nrow(data.test),res)
      colnames(res2) = c("Total Comparisons","Correct Predictions","Accuracy")
      as.data.frame(res2)
    }
  })
  
  output$show_result_summary <- renderDataTable({  
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    in_normtrans = "tmm-logcpm"
    in_method = input$method
    in_cross_val <- input$cross_val
    in_cross_val = as.numeric(in_cross_val)
    in_repeats_num <- input$repeats_num
    in_repeats_num = as.numeric(in_repeats_num)
    in_tuneLength <- input$tuneLength
    in_tuneLength = as.numeric(in_tuneLength)
    if (in_vars_num<=1000){
      comparison = req(result_classification())
      comparison = comparison[,-1]
      right_predictions = sum(comparison[1,]==comparison[2,])
      accuracy = round((right_predictions/dim(comparison)[2]),4)
      res = cbind(right_predictions,accuracy)
      res2 =  cbind(dim(comparison)[2],res)
      res3 =  cbind(res2,in_method)
      #res4 =  cbind(res3,in_normtrans)
      res5 =  cbind(res3,in_cross_val)
      res6 =  cbind(res5,in_repeats_num)
      res7 =  cbind(res6,in_tuneLength)
      res2 =  cbind(res7,in_vars_num)
      colnames(res2) = c("Total Comparisons","Correct Predictions","Testing Accuracy","Method","Number of Folds", "Number of Repeats","Tune Length", "Highest Var Features")
      res2 = as.data.frame(res2)
      res2
    }else{
      res2= "Number of features greater then 500, Please, choose a number less then 500."
      res2
    }
  })
  
  output$show_confu_matrix <- renderDataTable({ 
    res = confu_matrix()
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    if (in_vars_num<=1000){
       result = as.data.frame(res$conf_mat)
       result
    }else{
       result = "Number of features greater then 500, Please, choose a number less then 500."
    }

  })
  
  output$show_complete_results <- renderDataTable({ 
    res = confu_matrix()
    result = res$complete_results[,-1]
    round( result, 4)
  })
  
  confu_matrix <- reactive({
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    res=NULL
    conf_matrix=NULL
    if (in_vars_num<=1000){
    comparison = req(result_classification())
    
    comparison = comparison[,-1]
    comparison_df <- data.frame(True = as.character(comparison[2,]), Predicted = as.character(comparison[1,]))
    lvls <- sort(union(comparison_df$Predicted, comparison_df$True))
    conf_mat <- table(Predicted = factor(comparison_df$Predicted, levels = lvls), True = factor(comparison_df$True, levels = lvls) )
 
    conf_matrix <- conf_mat 
    
    # Create confusion matrix: True (columns), Predicted (rows)
    #print("Confusion Matrix:")
    #print(conf_matrix)
    
    # Identify matching class labels (those that exist in both predicted and true)
    matching_labels <- intersect(rownames(conf_matrix), colnames(conf_matrix))
    
    # Sum of correct predictions = sum of diagonal elements for matching labels
    correct_predictions <- sum(sapply(matching_labels, function(label) conf_matrix[label, label]))
    
    # Total number of predictions = sum of all elements
    total_predictions <- sum(conf_matrix)
    
    # Calculate accuracy
    accuracy <- correct_predictions / total_predictions

    all_classes <- union(rownames(conf_matrix), colnames(conf_matrix))
    
    # Accuracy
    #print('accuracy')
    #print( accuracy)
   
    # Per-class precision and recall
    precision <- diag(conf_mat) / rowSums(conf_mat) # Precisione = TP / (TP + FP) 
    recall <- diag(conf_mat) / colSums(conf_mat) #  Recall = TP / (TP + FN) 
    #Recall is vital when the cost of missing a positive case (a false negative) is high, such as in medical screening for diseases or fraud detection
    
    # Build final results table
    results_table <- data.frame(
      Class = all_classes,
      Precision = round(precision[all_classes], 7),
      Recall = round(recall[all_classes], 7)
    )
    res = results_table
    
    }else{
      results_table= "Number of features greater then 500, Please, choose a number less then 500."
      res = results_table
    }
    list(complete_results = res, conf_mat = conf_matrix)
  })
  
  confu_matrix2 <- reactive({

      comparison = req(result_classification())
      comparison = comparison[,-1]
      comparison_df <- data.frame(True = as.character(comparison[2,]), Predicted = as.character(comparison[1,]))
      lvls <- sort(union(comparison_df$Predicted, comparison_df$True))
      conf_mat <- table(Predicted = factor(comparison_df$Predicted, levels = lvls), True = factor(comparison_df$True, levels = lvls) )
      df = as.data.frame(conf_mat)
      
      expanded <- df[rep(seq_len(nrow(df)), df$Freq), 1:2]
      expanded$True <- factor(expanded$True, levels = lvls)
      expanded$Predicted <- factor(expanded$Predicted, levels = lvls)
      
      # Confusion matrix + metrics
      cm <- caret::confusionMatrix(expanded$Predicted, expanded$True)
      list(by_class_res = cm$byClass, overall_res = cm$overall)

  })
  

  output$show_by_class_results <- renderDataTable({ 
    res = confu_matrix2()
    round(res$by_class_res,4)
  })
  by_class_results <- reactive({ 
    res = req(confu_matrix2())
    round(res$by_class_res,4)
  })
  # Downloadable csv of selected dataset ----
  output$download_by_class_results <- downloadHandler(
    filename = function() {
      paste("class_results_",input$file3,"_",input$file4,"_method_",input$method,"_predictors_",input$vars_num,"_fold_",
            input$cross_val,"_repeats_",input$repeats_num,"_tuneLength_",input$tuneLength,".txt", sep = "")
    },
    content = function(file) {
      write.table(by_class_results(), file, row.names = TRUE, quote= FALSE, sep = "\t")
    }
  )
  
  output$show_overall_results <- renderDataTable({ 
    res = confu_matrix2()
    values = res$overall
    round(as.data.frame(values),4)
  })
  
  output$show_overall_results_plot <- renderPlotly({
    res = req(confu_matrix2())
    p <- plot_ly(y = round(res$overall,4), x= c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue"),  type= 'bar',
                 marker = list(size = 10, color = 'blue',
                               line = list(color = 'blue',width = 2))
    )
    p <- layout(p,title = "Bar Plot", 
                xaxis = list(title = ''), yaxis = list(title = ''))
    p
  })
  
  overall_results <- reactive({ 
    res = req(confu_matrix2())
    round(as.data.frame(res$overall),4)
  })
  
  # Downloadable csv of selected dataset ----
  output$download_overall_results <- downloadHandler(
    filename = function() {
      paste("global_results_",input$file3,"_",input$file4,"_method_",input$method,"_predictors_",input$vars_num,"_fold_",
            input$cross_val,"_repeats_",input$repeats_num,"_tuneLength_",input$tuneLength,".txt", sep = "")
    },
    content = function(file) {
      write.table(overall_results(), file, row.names = TRUE, quote= FALSE, sep = "\t")
    }
  )
  
  output$show_result_classification <- renderDataTable({ 
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    if (in_vars_num<=1000){
    comparison = req(result_classification())
    res=t(comparison)
    colnames(res) = res[1,]
    res = res[-1,]
    }else{
      res= "Number of features greater then 500, Please, choose a number less then 500."
      res
    }
    res
  })
  

  output$show_training_results <- renderDataTable({ 
    req(results2())
    in_method = input$method
    discrete_classifiers = c("PLDA", "PLDA2", "NBLDA")
    voom_classifiers = c("voomDLDA", "voomDQDA", "voomNSC")
    res=NULL
    if(in_method == "voomNSC") {res = round(results2()$trained_fit@finalModel$model$modelRes,4)}
    if((in_method == "voomDLDA")||(in_method == "voomDQDA") || (in_method == "PLDA")){
      res = as.data.frame(round(results2()$trained_fit@finalModel$accuracy,7))
      colnames(res)="Training Accuracy"
    }
    if( !is.element(in_method, discrete_classifiers) && !is.element(in_method, voom_classifiers)){res = round(results2()$trained_fit,7)}
    res
  })
  
  output$cmPlot <- renderPlotly({
    res = confu_matrix()
    df = as.data.frame(res$conf_mat)
    cm <- xtabs(Freq ~ True + Predicted, data = df)
    total <- sum(cm)
    
    # Base heatmap
    p <- plot_ly(
      x = colnames(cm),
      y = rownames(cm),
      z = cm,
      type = "heatmap",
      colors = colorRamp(c("#D6EAF8", "#2E86C1")),
      hoverinfo = "x+y+z"
    )
    
    # Add annotations (count + % of total)
    annotations <- list()
    for (i in 1:nrow(cm)) {
      for (j in 1:ncol(cm)) {
        count <- cm[i, j]
        perc <- round(100 * count / total, 1)  # % of total
        annotations <- c(annotations, list(
          list(
            x = colnames(cm)[j],
            y = rownames(cm)[i],
            text = paste0(count, " (", perc, "%)"),
            showarrow = FALSE,
            font = list(color = "white", size = 12)
          )
        ))
      }
    }
    p <- layout(
      p,
      title = "Confusion Matrix Heatmap",
      xaxis = list(title = "Predicted"),
      yaxis = list(title = "True"),
      annotations = annotations
    )
    p
  })
  
  output$costPlot <- renderPlotly({
    req(results2())
    in_method = input$method
    discrete_classifiers = c("PLDA", "PLDA2", "NBLDA")
    voom_classifiers = c("voomDLDA", "voomDQDA", "voomNSC")
    if(is.element(in_method, voom_classifiers) || is.element(in_method, discrete_classifiers)){
      p <- plot(results2()$fit)
    }else{
    res = trained(results2()$fit)[[4]]
    Training_cicles = (1:dim(res)[1])
    res = cbind(Training_cicles,res)
    p <- plot_ly(data = res, y = ~Accuracy, x = ~Training_cicles,  type= 'scatter', mode= 'lines+markers',
                   marker = list(size = 10, color = 'rgba(255, 182, 193, .9)',
                   line = list(color = 'rgba(152, 0, 0, .8)',width = 2))
                )
    p <- layout(p,title = "Training Accuracy values across the training cicles", 
                xaxis = list(title = 'Training Cicles'), yaxis = list(title = 'Training Accuracy'))
    }
    p
  })
  
  
  output$explore_methods <- renderDataTable({ 
    inFile1 <- input$file1
    inFile2 <- input$file2
    if (is.null(inFile1) || is.null(inFile2) ){return(NULL)}else{
    cervical = read.table(file = inFile1$datapath, header = TRUE,sep = "\t",row.names=1)
    templates = read.table(file= inFile2$datapath, header = FALSE,sep = "\t",row.names=1)
    comparison=NULL
    in_method_vector <- c('cforest','RRF','svmRadial','bagFDA','bayesglm','amdai','sparseLDA','lda','spls','hda','nnet',
                          'pcaNNet','kknn','avNNet','glmnet','ctree','treebag','BstLm','bstSm','rpart',
                          "voomDLDA", "voomDQDA", "voomNSC","PLDA", "PLDA2", "NBLDA","mlp","pam")
    in_normtrans <- "tmm-logcpm"
    in_repeats_num <- 1
    in_vars_num  <- 30
    in_tuneLength <- 1
    in_cross_val <- 4
    temp_results = list(in_method=NULL, in_vars_num=NULL, in_tuneLength=NULL, in_cross_val=NULL, accuracy=NULL)
    results=NULL
    class <- data.frame(condition = factor( as.character(templates[2,]) ))
    set.seed(2128)
    vars <- sort(apply(cervical, 1, var, na.rm = TRUE), decreasing = TRUE) #highest gene-wise variances
    data <- cervical[names(vars)[1: in_vars_num], ]
    ind <- which( templates[3,]=="F")  #pick those samples from templates to select the test and train samples
    data.train <- as.matrix(data[ ,-ind] + 1)
    data.test <- as.matrix(data[ ,ind] + 1)
    classtr <- data.frame(condition = class[-ind, ])
    classts <- data.frame(condition = class[ind, ])
    data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,design = formula(~condition))
    data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,design = formula(~condition))
    #ctrl.continuous <- trainControl   (method = "repeatedcv", number = in_cross_val, repeats = in_repeats_num)
    # Define control lists.
    ctrl.continuous <- trainControl   (method = "repeatedcv", number = in_cross_val, repeats = in_repeats_num)
    ctrl.discrete   <- discreteControl(method = "repeatedcv", number = in_cross_val, repeats = in_repeats_num,tuneLength = in_tuneLength)
    ctrl.voom       <- voomControl    (method = "repeatedcv", number = in_cross_val, repeats = in_repeats_num,tuneLength = in_tuneLength)
    discrete_classifiers = c("PLDA", "PLDA2", "NBLDA")
    voom_classifiers = c("voomDLDA", "voomDQDA", "voomNSC")
    fit=NULL
    for(in_method  in in_method_vector){
      temp_results$in_method = in_method
      temp_results$in_vars_num = in_vars_num
      temp_results$in_tuneLength = in_tuneLength
      temp_results$in_cross_val = in_cross_val
      print(paste("in_method: ", in_method,sep=""))
      if (is.element(in_method, discrete_classifiers)){
        print("discrete classifiers")
        fit <- classify(data = data.trainS4, method = in_method, normalize = "cpm", ref = templates[2,3], control = ctrl.discrete)
      }else if(is.element(in_method, voom_classifiers)){
        fit <- classify(data = data.trainS4, method = in_method, preProcessing = in_normtrans, ref = templates[2,3],tuneLength = in_tuneLength, control = ctrl.voom)
      }else{
        fit <- classify(data = data.trainS4, method = in_method, preProcessing = in_normtrans, ref = templates[2,3],tuneLength = in_tuneLength, control = ctrl.continuous)
      }
      #fit <- classify(data = data.trainS4, method = in_method, preProcessing = in_normtrans, ref = templates[2,1],tuneLength = in_tuneLength, control = ctrl.continuous)
      pred.svm = predict(fit, data.testS4)
      predicted_classes <- relevel(pred.svm, ref = templates[2,3])
      actual_classes    <- relevel(classts$condition, ref = templates[2,3])
      comparison = rbind(as.character(predicted_classes),as.character(actual_classes))
      comparison = cbind(c("predicted classes","true classes"),comparison)
      colnames(comparison) = c("classes" , colnames(data.test))
      comparison_no_first_col = comparison[,-1] #togli la prima colonna perchè sono solo i nomi delle variabili
      #trained_fit= MLSeq::trained(fit)[[4]]
      trained_fit=NULL
      if(is.element(in_method, voom_classifiers) || is.element(in_method, discrete_classifiers)){
        trained_fit= MLSeq::trained(fit)
      }else{
        trained_fit= MLSeq::trained(fit)[[4]]
      }
      comparison = comparison_no_first_col
      comparison_df <- data.frame(True = as.character(comparison[2,]), Predicted = as.character(comparison[1,]))
      lvls <- sort(union(comparison_df$Predicted, comparison_df$True))
      conf_mat <- table(Predicted = factor(comparison_df$Predicted, levels = lvls), True = factor(comparison_df$True, levels = lvls) )
      df = as.data.frame(conf_mat)
      expanded <- df[rep(seq_len(nrow(df)), df$Freq), 1:2]
      expanded$True <- factor(expanded$True, levels = lvls)
      expanded$Predicted <- factor(expanded$Predicted, levels = lvls)
      cm <- caret::confusionMatrix(expanded$Predicted, expanded$True)
      accuracy = round(sum(comparison_no_first_col[1,]==comparison_no_first_col[2,])/dim(comparison_no_first_col)[2],4)
      print(paste("Testing Accuracy: ", accuracy, sep=""))
      temp_results$accuracy = accuracy
      results = rbind(results, temp_results)
      by_class_res = cm$byClass
      overall_res = cm$overall
    }
    res = (as.data.frame(results))
    res=data.frame(unlist(res$in_method), unlist(res$in_vars_num), unlist(res$in_tuneLength), unlist(res$in_cross_val), unlist(res$accuracy))
    res[which( max(res$unlist.res.accuracy.)==res$unlist.res.accuracy. ),]
    sort(res$unlist.res.accuracy.)
    res = res %>% arrange(desc(unlist.res.accuracy.))
    res
    }
  })

})


  
  




