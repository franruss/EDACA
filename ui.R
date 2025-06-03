library(markdown)
library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
library(plotly)
library(diagram)
library(gplots)
library(shinythemes)
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
library(DESeq2)
source("definitions.R")

options(shiny.maxRequestSize=70*1024^2)
options(repos = BiocManager::repositories())
# Define UI 
shinyUI(fluidPage(theme = shinytheme("flatly"),
            #shinythemes::themeSelector(),  
  # Application title
  titlePanel("EDACA: Expression DAta Classifier App"),
  navbarPage("A Web App for Classifying RNASeq Expression Data", 
             
 tabPanel("Supervised Machine Learning Algorithms",
          sidebarLayout(
            sidebarPanel(
              helpText("You can use this interface by uploading an expression file in a tab separetad format, as the one shown below (we just display the first four lines):"),
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
              fileInput('file3', 'Choose an expression file',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
              fileInput('file4', 'Choose a training and test design matrix',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
              helpText("Here, you can select a number of features (e.g.: genes) having the highest variances accross samples to be considered by the method in order to decrease computational cost."),
              numericInput('vars_num', 'Number of features with highest variances to be considered (MAX allowed is 500)',"20"),
              radioButtons('method', 'Select a method: ',
                           c(
                             'Conditional Random Forest    '='cforest',
                             'Regularized Random Forest    '='RRF',
                             'Support Vector Machine       '='svmRadial',
                             'Adaptive Mixture Discriminant Analysis' = 'amdai',
                             'Sparse Linear Discriminant Analysis' = 'sparseLDA',
                             'Linear Discriminant Analysis '='lda',
                             'Voom-based Diagonal Quadratic Discriminant Analysis' = 'voomDQDA',
                             'Sparse Partial Least Squares '='spls',
                             'Neural Network               '='nnet',
                             'Neural Networks with Feature Extraction' = 'pcaNNet',
                             'k-Nearest Neighbors Network' = 'kknn',
                             'Averaged Neural Network' = 'avNNet',
                             'Multi-Layer Perceptron' = 'mlp',
                             'Multi-Layer Perceptron, with multiple layers ' = 'mlpML',
                             'Generalized Linear Model with Elastic Net Regression' = 'glmnet',
                             'Conditional Inference Tree' = 'ctree',
                             'Bagged CART' = 'treebag',
                             'Naive Bayes' = 'nb'
                            ),
              ),
              radioButtons('normtrans', 'Select a Normalization-Transformation combination: ',
                           c('Deseq Median Ratio Normalization + Variance Stabilizing Transformation'      = 'deseq-vst',
                             #'Deseq Median Ratio Normalization + Regularized Logarithmic Transformation'   = 'deseq-rlog', #computational expensive
                             #'Deseq Median Ratio Normalization + Log of Counts-Per-Million Transformation'='deseq-logcpm', #errors
                             'Trimmed Mean of M values (TMM) Norm. + Log of Counts-Per-Million Transformation '='tmm-logcpm'
                           ),
              ),
              submitButton("Train and Test!"),
              tags$head(tags$style(
                HTML(".shiny-notification {background-color:yellow;position:fixed;top: 50%;left: 5%;right: 5%;}"))),
              tags$hr(),
              width = 10,
            ),
            # Show a tabset
            mainPanel(
              tabsetPanel(type = "tabs", tabPanel("Your expression file", tableOutput('results8'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Your design matrix", tableOutput('results9'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Test results", tableOutput('show_result_classification'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Result summary", tableOutput('show_result_summary'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Confusion matrix", tableOutput('show_confu_matrix'))),
              tags$hr(),
              tabsetPanel(type = "tabs", tabPanel("Complete results", tableOutput('show_complete_results'))),
            )
          )
 ),
)))