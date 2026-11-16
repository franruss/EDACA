library(DESeq2)
library(earth)
library(markdown)
library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
library(shinythemes)
library(plotly)
library(diagram)
library(gplots)
library(png)
library(DT)
library(knitr)
library(MLSeq)
library(pamr)
library(caret)
library(kernlab)
library(CMScaller)
library(glmnet) 
library(RSNNS)
library(kknn)
library(klaR)
library(LiblineaR)
library(sparseLDA)
library(kerndwd)
library(party)
library(randomForest) 
library(RRF)
library(spls)
library(shinyalert)
library(preprocessCore)
library(stats)
library(keras)
library(mlbench)
library(dplyr)
library(magrittr)
library(tidyverse)
library(neuralnet)
library(BiocManager)
library(NeuralNetTools)
library(dslabs)
library(dplyr)
library(igraph)
library(ggraph)
library(ggrepel)
library(ggplot2)
library(tidyr)
library(arm)
library(bst)
library(hda)
library(adaptDA) #https://cran.r-project.org/src/contrib/Archive/adaptDA/
library(AnnotationDbi)
library(org.Hs.eg.db)
source("definitions.R")

options(shiny.maxRequestSize=70*1024^2)
options(repos = BiocManager::repositories())
# Define UI 
shinyUI(fluidPage(theme = shinytheme("yeti"),
            #shinythemes::themeSelector(),  
            ## Custom style
            tagList(
              # CSS for navbar elements
              tags$style(
                HTML('.navbar {font-size: 17px;}',
                     '.navbar-default .navbar-brand {color: #f2a2ec; font-weight: bold}',
                     '.navbar-default .navbar-brand {color: #f2a2ec;}',
                     '.navbar-default .navbar-nav > li > a {color: #b31414;}')
              ), 
              # CSS for fonts
              tags$style('h3 {font-weight: normal;}',
                         'h4 {font-weight: normal;}',
                         '* {font-family: Ubuntu;}'),
              # CSS for all buttons (class .btn)
              tags$style('.btn {color: #f1f518; 
                      background-color: #0d822c; 
                      border-color: #36447a;}
                .btn:hover {background-color: #734357;}'), 
              # CSS for errors and validation messages
              tags$style('.shiny-output-error-validation {
                 color: #e35300;
                 font-weight: bold;}'),
              # CSS for individual elements: note the '#id' syntax
              tags$style('#map-readme {color: #63071d;
                             background-color: #dce8e8;
                             border-color: #c84407;}
                #map-readme:hover {background-color: #e36a0075;}
                #trend-readme {color: #d1b6db;
                               background-color: #ff770050;
                               border-color: #c84407;}
                #trend-readme:hover {background-color: #e36a0075;}')
            ), # end of style block
            # Note how both styling *and* navbarPage are within fluidPage.
            
  # Application title
  titlePanel("EDACA: Expression DAta Classifier App"),
  navbarPage("A Web App for Classifying RNASeq Expression Data", 
             
 tabPanel("Supervised Machine Learning Algorithms",
          sidebarLayout(
            sidebarPanel(
              helpText("You can use this interface by uploading an expression file in a tab separated format, as the one shown below (we just display the first four lines):"),
              helpText(" "),
              helpText(HTML("Gene&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID27&nbsp;&nbsp;&nbsp;&nbsp;ID33&nbsp;&nbsp;&nbsp;ID68&nbsp;&nbsp;&nbsp;ID70")) ,
              helpText(HTML("A1BG0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;410&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;14&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4")),
              helpText(HTML("A1CF&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;192&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;156&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;63")),
              helpText(HTML("GGACT&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;28&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;23&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;15")),
              helpText(HTML("...  ...  ...")),
              helpText(" "),
              helpText("You also need to upload a design file in a tab separated format with 3 lines only. The first line describes the samples, the second line describes the classification known for each sample and the third line specifies which sample the 
                       algorithm you choose must use in the training process (with letter 'T') and in the testing process (with letter 'F'), as the one shown below:"),
              helpText(" "),
              helpText(HTML("Samples&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbspID27&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID33&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID68&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID70")) ,
              helpText(HTML("type&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cms1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cms2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cms3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;NotClassified")),
              helpText(HTML("train&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;T&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;T&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;T&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F")),
              helpText(" "),
              helpText("Finally, choose you classificator with the radio button from the list below and click on 'Train and Test!' button."),
              helpText("---------------------------------------------------------------------------------- "),
              helpText(" "),
              helpText("YOU CAN TRY THIS SAMPLE FILE !!!"),
              helpText(" "),
              checkboxInput("use_example", "Check this box to try EDACA on an expression and a design sample datasets!", value = FALSE),
              helpText(" "),
              helpText("These sample files are named '' and are available at: "),
              helpText(" "),
              helpText("---------------------------------------------------------------------------------- "),
              fileInput('file3', 'Upload an expression file',accept=c('text/csv','text/comma-separated-values,text/plain','.csv',".txt.gz", ".gz",".zip")),
              fileInput('file4', 'Upload a design file',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
              sliderInput(inputId = "cross_val", label = "Cross validation fold:", min = 2, max = 10, value = 4,width = 200),
              radioButtons('selector','Choose the ay to select predictors: ', c('use predictors with highest variance'='variance', 'use a list of predictors'='pred_list')),
              textInput(inputId = "predictor_list", label = "Write a list of predictors:","ACTR3B,ANLN,BAG1,BCL2,BIRC5,BLVRA,CCNB1,CCNE1,CDC20,CDC6,NUF2,CDH3,CENPF,CEP55,CXXC5,EGFR,ERBB2,ESR1,EXO1,FGFR4,FOXA1,FOXC1,GPR160,GRB7,KIF2C,NDC80,KRT14,KRT17,KRT5,MAPT,MDM2,MELK,MIA,MKI67,MLPH,MMP11,MYBL2,MYC,NAT1,ORC6,PGR,PHGDH,PTTG1,RRM2,SFRP1,SLC39A6,TMEM45B,TYMS,UBE2C,UBE2T"),
              helpText("Here, you can select a number of features (e.g.: genes) having the highest variances accross samples to be considered by the method."),
              sliderInput(inputId = "vars_num", label = "Number of predictors with highest variance:", min = 2, max = 1000, value = 4,width = 600),
              sliderInput(inputId = "repeats_num", label = "Number of repeats:", min = 1, max = 4, value = 1,width = 100),
              sliderInput(inputId = "tuneLength", label = "Training cicles:", min = 2, max = 10, value = 4,width = 200),
              selectInput("method", "Select a classification algorithm: ", c(
                             ' 1) Conditional Inference Random Forest    '='cforest',
                             ' 2) Regularized Random Forest              '='RRF',
                             ' 3) Support Vector Machine                 '='svmRadial',
                             ' 4) Bagged Flexible Discriminant Analysis  '='bagFDA',
                             ' 5) Bayesian generalized linear model      '='bayesglm',
                             ' 6) Adaptive Mixture Discriminant Analysis '='amdai',
                             ' 7) Sparse Linear Discriminant Analysis    '='sparseLDA',
                             ' 8) Linear Discriminant Analysis           '='lda',
                             ' 9) Sparse Partial Least Squares           '='spls',
                             '10) Heteroscedastic Discriminant Analysis  '='hda',
                             '11) Neural Network                         '='nnet',
                             '12) Neural Networks with PCA step          '='pcaNNet',
                             '13) k-Nearest Neighbors Network            '='kknn',
                             '14) Averaged Neural Network                '='avNNet',
                             '15) Generalized Linear Model with Elastic Net Regression' = 'glmnet',
                             '16) Conditional Inference Tree             '='ctree',
                             '17) Boosted Linear Model                   '='BstLm',
                             '18) Boosted Smoothing Spline               '='bstSm',
                             '19) CART (rpart)                           '='rpart',
                             '20) Bagged CART (treebag)                  '='treebag',
                             '21) Multi-Layer Perceptron                 '='mlp',
                             '22) Poisson Linear Discriminant Analysis   '='PLDA',
                             '23) Poisson Linear Discriminant Analysis with Power Transformation'='PLDA2',
                             '24) Negative Binomial Linear Discriminant Analysis '='NBLDA',
                             '25) Voom-based Diagonal Linear Discriminant Analysis' = 'voomDLDA',
                             '26) Voom-based Diagonal Quadratic Discriminant Analysis' = 'voomDQDA',
                             '27) Voom-based Nearest Shrunken Centroids' = 'voomNSC',
                             '28) Nearest Shrunken Centroids' = 'pam'
                             #'AdaBoost.M1, SAMME and Bagging' = 'AdaBag'
                             #'Multi-Layer Perceptron, with multiple layers ' = 'mlpML',
                             #'Naive Bayes' = 'nb',
                             #'blackboost' = 'blackboost',
                            ),
              ),
              actionButton("run", "Train and Test!"),
              tags$head(tags$style(
                HTML(".shiny-notification {background-color:yellow;position:fixed;top: 50%;left: 5%;right: 5%;}"))),
              tags$hr(),
              width = 10,
            ),
            # Show a tabset
            mainPanel(
              tabsetPanel(type = "tabs", tabPanel("Your expression file", dataTableOutput('results8'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Your design file", dataTableOutput('results9'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Summary of the Training process", verbatimTextOutput("summaryTrain"))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Training results (with Training Accuracy)", dataTableOutput('show_training_results'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Training Accuracy plot", plotlyOutput('costPlot'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Summary of the Testing process", verbatimTextOutput("summary"))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Test results", dataTableOutput('show_result_classification'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Test result table", dataTableOutput('show_result_summary'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Confusion matrix heatmap", plotlyOutput('cmPlot'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Global results (with Testing Accuracy) table", dataTableOutput('show_overall_results'))),
              downloadButton("download_overall_results", "Download"),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Global results (with Testing Accuracy) bar plot", plotlyOutput('show_overall_results_plot'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Class Results", dataTableOutput('show_by_class_results'))),
              downloadButton("download_by_class_results", "Download"),
              tags$hr(),
              tags$hr(),
            )
          )
 ),
 
 tabPanel("Look for the Best Machine Learning Tools for your Data",
          sidebarLayout(
            sidebarPanel(
              helpText("You can use this interface by uploading an expression file in a tab separated format, as the one shown below (we just display the first four lines):"),
              helpText(" "),
              helpText(HTML("Gene&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID27&nbsp;&nbsp;&nbsp;&nbsp;ID33&nbsp;&nbsp;&nbsp;ID68&nbsp;&nbsp;&nbsp;ID70")) ,
              helpText(HTML("A1BG0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;410&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;14&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4")),
              helpText(HTML("A1CF&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;192&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;156&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;63")),
              helpText(HTML("GGACT&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;28&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;23&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;15")),
              helpText(HTML("...  ...  ...")),
              helpText(" "),
              helpText("You also need to upload a design file in a tab separated format with 3 lines only. The first line describes the samples, the second line describes the classification known for each sample and the third line specifies which sample the 
                       algorithm you choose must use in the training process (with letter 'T') and in the testing process (with letter 'F'), as the one shown below:"),
              helpText(" "),
              helpText(HTML("Samples&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbspID27&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID33&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID68&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID70")) ,
              helpText(HTML("type&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cms1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cms2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cms3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;NotClassified")),
              helpText(HTML("train&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;T&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;T&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;T&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F")),
              helpText(" "),
              helpText("Finally, choose you classificator with the radio button from the list below and click on 'Explore All Tools!' button."),
              fileInput('file1', 'Upload an expression file',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
              fileInput('file2', 'Upload a design file',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
              actionButton("run2", "Explore All Tools!"),
              tags$head(tags$style(
                HTML(".shiny-notification {background-color:yellow;position:fixed;top: 50%;left: 5%;right: 5%;}"))),
              tags$hr(),
              width = 10,
            ),
            # Show a tabset
            mainPanel(
              tabsetPanel(type = "tabs", tabPanel("Results", dataTableOutput('explore_methods'))),
              tags$hr(),
            )
          )
 ),
)))