
source("definitions.R")


library(BiocManager)
options(repos = BiocManager::repositories())

shinyServer(function(input, output) {
  
  #Classification
  output$results <- renderTable({  
    inFile <- input$file1
    inFile2 <- input$file2
    if (is.null(inFile) || is.null(inFile2) ){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
      emat = read.csv(inFile$datapath, header=TRUE, sep="\t",row.names=1)
      #emat = read.table(file ="all_entrez.txt", header = TRUE,sep = "\t",row.names=1)
      #templates = read.table(file ="prediction_table.txt", header = TRUE,sep = "\t")
      templates = read.table(inFile2$datapath, header=TRUE, sep="\t")
      counts <- emat
      templates$class = factor(templates$class, levels = unique(templates$class) )  #devo dare i livelli a templates$class altrimenti la funzione sotto non gira
      counts2 <- limma::normalizeQuantiles(log2(counts + 0.25))
      res = ntp_mine(counts2, templates, verbose =TRUE)
      sample = rownames(res)
      res = cbind(sample,res)
      res
    }
  })
  
  #Training and Test
    results2 <- reactive({
    comparison=NULL
    inFile3 <- input$file3
    inFile4 <- input$file4
    in_method <- input$method
    in_normtrans <- input$normtrans
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    # in_method <- "svmRadial"
    if (is.null(inFile3) || is.null(inFile4) ){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
      # setwd("C:/Users/f.russo.ENGIBBC/Desktop/Reti Transomiche/")
      cervical = read.table(file = inFile3$datapath, header = TRUE,sep = "\t",row.names=1)
      templates = read.table(file= inFile4$datapath, header = FALSE,sep = "\t",row.names=1)
      # cervical = read.table(file ="all.txt", header = TRUE,sep = "\t",row.names=1)
      # templates = read.table(file ="training_test.txt", header = FALSE,sep = "\t",row.names=1)
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
      vars <- sort(apply(cervical, 1, var, na.rm = TRUE), decreasing = TRUE) #highest gene-wise variances
      data <- cervical[names(vars)[1: in_vars_num], ]
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
      # Support Vector Machines with Radial Kernel
      # QUI L'UTENTE DEVE SCERGLIERE IL METHOD :
      withProgress(message = 'Running, Please Wait!', detail = 'This may take a while...',  value = 0, {
        # Number of times we'll go through the loop
        n <- 100
        print("You picked the algorithm:")
        print(as.character(in_method))
        print("You selected the top ")
        print(in_vars_num)  
        print("You selected the Norm. plus Trans.")
        print(in_normtrans)  
        print(" features having the highest variances!")
        #using 5-fold cross validation without repeats. (number = 5,repeats = 1)
      fit <- classify(data = data.trainS4, method = in_method, preProcessing = in_normtrans, ref = templates[2,1],
                     control = trainControl(method = "repeatedcv", number = 4,repeats = 1, classProbs = TRUE))
      pred.svm = predictClassify(fit, data.testS4)
      predicted_classes <- relevel(pred.svm, ref = templates[2,1])
      actual_classes    <- relevel(classts$condition, ref = templates[2,1])
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
    }
    comparison_no_first_col = comparison[,-1] #togli la prima colonna perchè sono solo i nomi delle variabili
    print("Right predictions:")
    print(sum(comparison_no_first_col[1,]==comparison_no_first_col[2,]))
    print("Accuracy:")
    print(sum(comparison_no_first_col[1,]==comparison_no_first_col[2,])/dim(comparison_no_first_col)[2])
    comparison 
  })
  
  result_classification <- reactive({  
    req(results2())
    comparison <- results2()
    comparison_no_first_col = comparison[,-1] #togli la prima colonna perchè sono solo i nomi delle variabili
    print("Right predictions:")
    print(sum(comparison_no_first_col[1,]==comparison_no_first_col[2,]))
    print("Accuracy:")
    print(sum(comparison_no_first_col[1,]==comparison_no_first_col[2,])/dim(comparison_no_first_col)[2])
    comparison
  })
  
  output$unsupervised_result_classification <- renderTable({  
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
  
  
  output$results5 <- renderTable({  
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

  
  output$results8 <- renderTable({  
    inFile3 <- input$file3
    if (is.null(inFile3) ){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
      counts = read.table(inFile3$datapath,header=TRUE, row.names=1)
      id=rownames(counts)
      (cbind(id,counts))[1:8,]
    }
  })

  
  output$results9 <- renderTable({  
    inFile3 <- input$file3
    inFile4 <- input$file4
    if (is.null(inFile3) || is.null(inFile4) ){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
      counts = read.table(inFile4$datapath,header=TRUE,row.names=1)
      id=rownames(counts)
      cbind(id,counts)
    }
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
      print(head(data))
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
  
  output$show_training_results <- renderTable({  
    inFile11 <- input$file11
    inFile12 <- input$file12
    layers_words <- input$text10
    in_vars_num <- input$vars_num2
    layers <- as.numeric(strsplit(layers_words," ")[[1]])
    req(build_network())
    req(build_formula())
    req(build_data.train())
    if (is.null(layers) || is.null(inFile11) || is.null(inFile12) || is.null(in_vars_num) || is.null(layers_words)){return(NULL)}else{
      #  pred <- predict(nn, data.train)
      data.train <- build_data.train()
      nn <-build_network()
      
      pred <- predict(nn, data.train)
      #labels = unique(as.character(data.train$class))
      labels = levels(data.train$class)
      true_values <- (data.train[,dim(data.train)[2]])
      colnames(pred) = labels
      predictions <- (data.frame(max.col(pred))[,1])
      res = rbind(rbind(rownames(data.train),as.character(true_values)),labels[predictions])
      res = t(res)
      colnames(res) = c("ID", "true_values", "predictions")
      table_res = table(res[,2], res[,3])
      table_res = as.data.frame(table_res)
      colnames(table_res) = c("true values", "predictions","frequency")
      table_res
    }
  })
  
  output$show_testing_results <- renderTable({  
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
  
  output$show_input <- renderTable({  
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
  
  
  output$show_design <- renderTable({  
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
  
  output$show_accuracy <- renderTable({  
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
  
  output$show_result_summary <- renderTable({  
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    in_normtrans = input$normtrans
    in_method = input$method
    if (in_vars_num<=500){
      comparison = req(result_classification())
      comparison = comparison[,-1]
      right_predictions = sum(comparison[1,]==comparison[2,])
      print('right_predictions')
      print(right_predictions)
      accuracy = round((right_predictions/dim(comparison)[2]),3)
      print('accuracy')
      print(accuracy)
      res = cbind(right_predictions,accuracy)
      res2 =  cbind(dim(comparison)[2],res)
      res3 =  cbind(res2,in_method)
      res4 =  cbind(res3,in_normtrans)
      res5 =  cbind(res4, "4-fold validation without repeats")
      res2 =  cbind(res5,in_vars_num)
      colnames(res2) = c("Total comparisons","Correct predictions","Accuracy","Method chosen", "Normalization and Transformation chosen", "Type of training", "Number of features with highest variances chosen")
      res2 = as.data.frame(res2)
      res2
    }else{
      res2= "Number of features greater then 500, Please, choose a number less then 500."
      res2
    }
  })
  
  output$show_confu_matrix <- renderTable({ 
    res = confu_matrix()
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    result = NULL
    if (in_vars_num<=500){
       result = t(res$conf_mat)
    }else{
       result = "Number of features greater then 500, Please, choose a number less then 500."
    }
    result
  })
  
  output$show_complete_results <- renderTable({ 
    res = confu_matrix()
    res$complete_results
  })
  
  confu_matrix <- reactive({
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    res=NULL
    conf_matrix=NULL
    if (in_vars_num<=500){
    comparison = req(result_classification())
    
    comparison = comparison[,-1]
    comparison_df <- data.frame(True = as.character(comparison[2,]), Predicted = as.character(comparison[1,]))
    conf_mat <- table(Predicted = comparison_df$Predicted, True = comparison_df$True)
    conf_matrix <- conf_mat 
    
    # Create confusion matrix: True (columns), Predicted (rows)
    print("Confusion Matrix:")
    print(conf_matrix)
    
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
    print('accuracy')
    print( accuracy)
   
    # Per-class precision and recall
    precision <- diag(conf_mat) / rowSums(conf_mat)
    recall <- diag(conf_mat) / colSums(conf_mat)
    
    # Build final results table
    results_table <- data.frame(
      Class = all_classes,
      Precision = round(precision[all_classes], 6),
      Recall = round(recall[all_classes], 6)
    )
    
    res = results_table
    }else{
      results_table= "Number of features greater then 500, Please, choose a number less then 500."
      res = results_table
    }
    list(complete_results = res, conf_mat = conf_matrix)
  })
  
  
  output$show_result_classification <- renderTable({ 
    in_vars_num <- input$vars_num
    in_vars_num = as.numeric(in_vars_num)
    if (in_vars_num<=500){
    comparison = req(result_classification())
    comparison
    }else{
      comparison= "Number of features greater then 500, Please, choose a number less then 500."
      comparison
    }
  })
  

})


  
  




