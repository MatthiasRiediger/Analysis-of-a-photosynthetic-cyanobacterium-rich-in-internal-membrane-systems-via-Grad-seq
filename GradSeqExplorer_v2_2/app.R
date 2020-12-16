#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(shinyWidgets)
library(ggplot2)
library(gplots)
library(cowplot)
library(gridExtra)
library(grid)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)
library(scales)
library(DT)


#setwd("C:/Users/matth/OneDrive/Dokumente/HessLab_Data_PhD_Backup_130320/projects/3_GradSeq_MS/ShinyApp/Shiny/App/GradSeqExplorer_v2_karsten")
#setwd("U:/documents/Data/projects/3_GradSeq_MS/ShinyApp/Shiny/App/GradSeqExplorer_v2_karsten")

constProt <- read.csv("const_prot_norm.csv", sep=";")

mbgd_genomes <- read.csv("mbgd_genomes.csv", sep=";")

Shiny_dt <- read.csv("Shiny_datatable_v2_1.csv", sep=";")
add_info <- read.csv("add_info.csv", sep=",")

KEGG <- read.csv("KEGG.csv", sep=";")
KEGG <- merge(KEGG, Shiny_dt, by=c("locustag"))

mbgd <- read.csv("cons_mbgd.csv", sep=";")
mbgd <- mbgd  %>%  separate_rows(u_syn) %>% arrange(-as.numeric(Cyanobacteria))
mbgd <- cbind(mbgd[4], mbgd[7:50])
mbgd <- mbgd[!duplicated(mbgd$u_syn),]
mbgd_merge <- droplevels(subset(mbgd, Conservation != ""))
mbgd_merge <- cbind(mbgd_merge[1], mbgd_merge[43:44])
names(mbgd_merge)[1] <- "locustag"
mbgd <-  reshape2::melt(mbgd, id.vars=c("u_syn", "Proteobacteria", "Arabidopsis", "Cyanobacteria", "Conservation", "Percentage"))
names(mbgd)[1] <- c("syn_locustag")
KEGG <- merge(KEGG, mbgd_merge, by=c("locustag"),all.x = TRUE)
levels(KEGG$Conservation) <- c(levels(KEGG$Conservation), "unknown") 
KEGG <- KEGG %>% replace_na(list(Cyanobacteria = 1, Conservation = "unknown")) 
KEGG <- KEGG[!duplicated(KEGG),] 


SVM <- read.csv("SVM_mbgd.csv", sep=";")

syn_svm <- SVM %>%  filter(org == "u_syn") %>% arrange(-Ortholog_clustersize) %>% 
  distinct(locus_tags, .keep_all = TRUE)
syn_svm <- cbind(syn_svm[,1:2], syn_svm[6])
  names(syn_svm)[1:3] <- c("locustag", "Ortholog_ClusterID", "syn_svm")
unique_syn_svm <- syn_svm$Ortholog_ClusterID[!duplicated(syn_svm$Ortholog_ClusterID)]
svm_input <- SVM %>% filter(Ortholog_ClusterID %in% unique_syn_svm)
svm_input <- merge(svm_input, syn_svm, by=c("Ortholog_ClusterID"), all = TRUE) 
svm_median <- aggregate(svm_input[6], list(svm_input$locustag), median, na.rm=TRUE)
  names(svm_median) <- c("locustag", "SVM.score.Median")
svm_input <- merge(svm_input, svm_median, by=c("locustag"), all = TRUE) 
svm <- droplevels(subset(svm_input, org == "u_syn"))
svm <- cbind(svm[1], svm[,8:9])
  names(svm) <- c("locustag", "SVM.score", "SVM.score.Median")
svm <- svm[!duplicated(svm),]

KEGG <- merge(KEGG, svm, by=c("locustag"),all.x = TRUE)
KEGG <- KEGG[!duplicated(KEGG),] 



#syn_svm <- SVM %>% filter(locus_tags %in% gene_selected) %>% 
 # distinct(locus_tags, .keep_all = TRUE)





#Categories <- c("Category/Summary", "Category/Detail", "Pathway/Complex", "Pathway", "Pathway/Overview")

locustag <- as.character(Shiny_dt$Gene_ID_IndMode)

AnnotationListBar <- as.character(KEGG$MainCategory)
AnnotationListBar <- sort(AnnotationListBar[!duplicated(AnnotationListBar)])

AnnotationListBarComplex <- as.character(as.character(Shiny_dt$Pathway.Complex))
AnnotationListBarComplex <- sort(AnnotationListBarComplex[!duplicated(AnnotationListBarComplex)])



ui <- navbarPage("Grad-Seq Analysis of Synechocystis sp. PCC 6803", id="nav",
                 tabPanel("Introduction",
                          sidebarLayout( position = "left",
                                         sidebarPanel(
                                            p("If you are using this Shiny App for your own research please cite:", style = "font-family: 'times'; font-size: 16pt"),
                                           br(),  
                                           p("Matthias Riediger, Philipp Spaet, Raphael Bilger, Karsten Voigt, Boris Macek, Wolfgang R Hess,
                                              Analysis of a Photosynthetic Cyanobacterium Rich in Internal Membrane Systems via Gradient Profiling by Sequencing (Grad-seq), The Plant Cell, koaa017, ", 
                                             span(a("https://doi.org/10.1093/plcell/koaa017", href="https://doi.org/10.1093/plcell/koaa017")), style = "font-family: 'times'; font-size: 16pt"),
                                           br(),
                                           p("and visit our", span(a("Cyanolab Homepage", href="http://www.cyanolab.de/")), "for the latest updates of our research.", style = "font-family: 'times'; font-size: 16pt"),
                                           br(),
                                           br(),
                                           strong("Abstract", style = "font-family: 'times';font-size: 32pt"),
                                            p("Although regulatory small RNAs have been reported in photosynthetic cyanobacteria, the lack of clear RNA chaperones involved in their regulation 
                                              poses a conundrum. Here, we analyzed the full complement of cellular RNAs and proteins using gradient profiling by sequencing (Grad-seq) in",  
                                              em("Synechocystis"),"6803. Complexes with overlapping subunits such as the CpcG1-type versus the CpcL-type phycobilisomes or the PsaK1 versus PsaK2 
                                              photosystem I pre(complexes) could be distinguished, supporting the high quality of this approach. Clustering of the in-gradient distribution 
                                              profiles followed by several additional criteria yielded a short list of potential RNA chaperones that include a YlxR homolog and a cyanobacterial 
                                              homolog of the KhpA/B complex. The data suggest previously undetected complexes between accessory proteins and CRISPR-Cas systems, such as a Csx1-Csm6 
                                              ribonucleolytic defense complex. Moreover, the exclusive association of either RpoZ or 6S RNA with the core RNA polymerase complex and the existence of 
                                              a reservoir of inactive sigma-antisigma complexes is suggested. The", em("Synechocystis"),"Grad-seq resource provides a comprehensive resource for the 
                                              unctional assignment of RNA-protein complexes and multisubunit protein complexes in a photosynthetic organism.", style = "font-family: 'times'; font-size: 16pt")
                                          ),                             
                                          mainPanel(
                                            strong("Features", style = "font-family: 'times';font-size: 32pt"),
                                            br(),
                                            strong("General Description", style = "font-family: 'times'; font-size: 24pt"),
                                            p("Customizable visualizations for the", span(em("Synechocystis")), " sp. PCC6803 proteome and transcriptome, 
                                            according to the separation of RNA-protein and multi-protein complexes in a sucrose density gradient. The dataset has been subjected to a hierarchical clustering approach and all detected 
                                            proteins and RNAs were assigned to one of 17 distinct clusters, based on their individual sedimentation characteristics. These
                                            visualizations are also available as download in '.png' format (300 dpi).", style = "font-family: 'times'; font-size: 16pt"),
                                            br(),
                                            strong("Grad-Seq Explorer", style = "font-family: 'times'; font-size: 24pt"),
                                            p("Heatmap visualizations of sedimentation profiles ", span(strong("('Grad-Seq Heatmap')."))," The options can be set to search for 'locus tags/gene names', 'functional categories'
                                            (in hierarchical search order of ", span(a("KEGG BRITE", href="https://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+gn:T00004")),
                                            " -- database - category - pathway - complex/type) or comparison between selected locus tags 
                                            and selected functional categories. Multiple filter options are available for search within 'functional categories'.
                                            Further visualization tools are available for proteins only. The ", span(strong("'Phylogeny Heatmap'")), "shows the phylogenetic occurrence of selected proteins within selected cyanobacteria,",  
                                            span(em("A.thaliana, E.coli")), " or ", span(em("S.enterica")), "based on domclust of the", span(a("Microbial Genome Database (MBGD)", href="http://mbgd.genome.ad.jp/")),
                                            " and all selected genomes, given in the supplementary information tab", span(strong("(Table S1)."))," The", span(strong("'SVM score Boxplot'")), "represents the SVM score for the prediction of RNA-binding proteins by",
                                            span(a("RNApred", href="http://crdd.osdd.net/raghava/rnapred/")),"and all selected", span(em(" Synechocystis ")), "proteins (red dot) and their respective orthologues (boxplot), determined by MBGD. For details of the color classification refer to ",
                                            span(strong(" Figure S1 ")), "in the supplementary information tab.
                                            A tabular output of the selected proteins and RNAs with all information is available as well.",style = "font-family: 'times'; font-size: 16pt"),
                                            br(),
                                            strong("Gradient Composition", style = "font-family: 'times'; font-size: 24pt"),
                                            p("Barplots of proteins and RNA types by cluster assignment and peak fractions.",style = "font-family: 'times'; font-size: 16pt"),
                                            br(),
                                            strong("Protein Categories", style = "font-family: 'times'; font-size: 24pt"),
                                            p("Barplots of proteins of selected 'functional categories' by cluster assignment and peak fractions.", 
                                            "(in hierarchical search order of ", span(a("(KEGG BRITE", href="https://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+gn:T00004")),
                                            " -- database - category - pathway - complex/type)", style = "font-family: 'times'; font-size: 16pt")
                                            
                                            
                                          )
                                        
                        )
                 ),
                 
                 tabPanel("Tutorial",
                          
                          mainPanel(
                           
                            img(src='Tutorial1.jpg', height="100%", width="100%", align = "center"),
                            img(src='Tutorial2.jpg', height="100%", width="100%", align = "center"),
                            img(src='Tutorial3.jpg', height="100%", width="100%", align = "center"),
                            img(src='Tutorial4.jpg', height="100%", width="100%", align = "center"),
                            img(src='Tutorial5.jpg', height="100%", width="100%", align = "center"),
                            img(src='Tutorial6.jpg', height="100%", width="100%", align = "center")
                            
                            
                          )
                 ), #tabpanel
                 
   
                 tabPanel("Grad-Seq Explorer",
                          sidebarLayout( position = "left",
                                  sidebarPanel(
                                    
                                           selectInput(inputId = "height", label = "Please adjust Plot Height: ",
                                                       choices = c(1:10*10), selected=20),
                                           
                                           radioButtons(inputId = "NormChoice",  choices = c("Constant Volume", "Constant Protein"), label="Normalization (Measurements & Clustering performed with Constant Volume)", selected = c("Constant Volume"), inline=TRUE),
                                           
                                           
                                           radioButtons(inputId = "HeatmapChoice",  choices = c("Locustags", "Pathways", "Both"), label="Heatmap Selection", selected = c("Locustags"), inline=TRUE),
                                           
            
                                           conditionalPanel(condition = "input.HeatmapChoice == 'Both' ",
                                                    sliderInput(inputId = "heightRatio", label = "Please adjust Plot Height ratio: ",
                                                                        min = 0.1, max = 10, value = 1.5,  step = 0.1 )
                                                    
                                           ),                                    
                                           conditionalPanel(condition = "input.HeatmapChoice == 'Locustags' || input.HeatmapChoice == 'Both'",
                                                            
                                                            pickerInput(inputId = "Gene_ID", label = "Please choose Gene Identifier: ",
                                                                        choices = unique(as.character(Shiny_dt$Gene_ID_IndMode)), 
                                                                        selected=c("psaA/slr1834 (Protein) - TU958", "psaA/slr1834 (RNA) - TU958", "slr1834-3UTR (RNA) - TU958",
                                                                                   "slr1834-5UTR (RNA) - TU958", "TU958 (slr1834)"), 
                                                                        options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple = T)

                                           ),
                                           conditionalPanel(condition = "input.HeatmapChoice == 'Pathways' || input.HeatmapChoice == 'Both' ",
                                                            
                                                            pickerInput(inputId = "database", label = "Please select Functional Categories: ",
                                                                        choices = unique(as.character(KEGG$database)), selected=c("syn00194 Photosynthesis Proteins"), 
                                                                        options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                                            pickerInput(inputId = "category", choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T), 
                                                            pickerInput(inputId = "pathway",  choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                                            pickerInput(inputId = "complex",  choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                                            
                                                            pickerInput(inputId = "Feature", label = "Please select Feature Type(s): ",
                                                                        choices =  c("Protein", "3UTR", "5UTR", "mRNA", "sRNA", "tRNA", "rRNA", "aRNA", "crRNA", "Transposase", "transcript"),
                                                                        selected = "Protein", options = list(`actions-box` = TRUE),multiple = T,
                                                                        choicesOpt = list(style = c("background: lightgrey;",
                                                                                                    "background: rgba(150, 0, 100, 0.3);",
                                                                                                    "background: rgba(100, 0, 50, 0.3);",
                                                                                                    "background: rgba(200, 0, 150, 0.1);",
                                                                                                    "background: rgba(255, 0, 0, 0.5);",
                                                                                                    "background: rgba(0, 255, 255, 0.3);",
                                                                                                    "background: rgba(200, 0, 225, 0.5);",
                                                                                                    "background: rgba(0, 0, 255, 0.5);", 
                                                                                                    "background: rgba(255, 205, 0, 0.5);", 
                                                                                                    "background: rgba(0, 255, 0, 0.5);", 
                                                                                                    "background: rgba(255, 150, 0, 0.5);")) ),
                                                            
                                                            selectInput(inputId = "HeatmapOrder", label = "Please choose Heatmap Order: ", 
                                                                        choices= c("Cluster", "Feature Type", "Category", "Locustag", "Abundance"), selected="Cluster"),
                                                            
                                                            pickerInput(inputId = "Cluster", label = "Please select Cluster: ", 
                                                                        choices =  c(0:17), selected= c(1:17), options = list(`actions-box` = TRUE),multiple = T),
                                                            pickerInput(inputId = "Maximum", label = "Please select maximum Peak Fraction(s): ",
                                                                        choices = c(1:13,14.5,16.5,18), selected = c(1:13,14.5,16.5,18), options = list(`actions-box` = TRUE),multiple = T),
                                                            pickerInput(inputId = "Peak_no", label = "Please select number of Peaks: ",
                                                                        choices = c(1:6), selected = c(1:6), options = list(`actions-box` = TRUE),multiple = T),
                                                            
                                                            pickerInput(inputId = "Peaks", label = "Please select Peak Fraction Combinations(s): ",
                                                                        choices = c(1:13,14.5,16.5,18), selected = c(1:13,14.5,16.5,18), options = list(`actions-box` = TRUE),multiple = T),
                                                            
                                                            pickerInput(inputId = "Location", label = "Please select Location(s): ",
                                                                        choices = c("Chr", "pSYSA", "pSYSG", "pSYSM", "pSYSX", "pCA2.4"), 
                                                                        selected = c("Chr", "pSYSA", "pSYSG", "pSYSM", "pSYSX", "pCA2.4"), options = list(`actions-box` = TRUE),multiple = T),
                                                            
                                                            sliderInput(inputId = "AA_LENGTH", label = "Protein length [AA] (Protein only)",
                                                                        min = 38, max = 4200, value = c(38,4200),  step = 1 ),
                                                            
                                                            pickerInput(inputId = "Phylogeny", label = "Please select Phylogenetic groups (Protein only)", 
                                                                        choices = c("other Cyanobacteria", "other Cyanobacteria and A.thaliana", "other Cyanobacteria and S.enterica/E.coli", "All", "unknown"), 
                                                                        selected = c("other Cyanobacteria", "other Cyanobacteria and A.thaliana", "other Cyanobacteria and S.enterica/E.coli", "All", "unknown"),
                                                                        options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                                            
                                                            #pickerInput(inputId = "Phylogeny_Detail", label = "Please select Organisms (Protein only)", 
                                                            #            choices = levels(as.factor(mbgd$variable)), selected = levels(as.factor(mbgd$variable)),
                                                            #            options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                                            
                                                            
                                                            sliderInput(inputId = "Cyanobacteria", label = "Conservation within selected Cyanobacteria (Protein only)",
                                                                        min = 0, max = 1, value = c(0,1),  step = 0.05 ),
                                                            
                                                            radioButtons(inputId = "Median", label = "Please select SVM score classification: ",
                                                                        choices = c("Synechocystis 6803", "Median"), selected = c("Synechocystis 6803"), inline=TRUE),
                                                            
                                                            
                                                            sliderInput(inputId = "SVM_score", label = "Synechocystis SVM score (Protein only)",
                                                                        min = -16, max = 6, value = c(-16,6),  step = 0.1 ),
                                                            
                                                            sliderInput(inputId = "Median_SVM_score", label = "Median SVM score (Protein only)",
                                                                        min = -16, max = 6, value = c(-16,6),  step = 0.1 )
                                                            
                                                            
                                           ),
                                           sliderInput(inputId = "Plot_width", label = "Please select Plot width for download:",
                                                       min = 1, max = 100, value = 20,  step = 5 ),
                                           
                                          downloadButton("downloadPlot", "Download Grad-Seq Heatmap"),
                                          downloadButton("downloadPlot2", "Download Phylogeny Heatmap"),
                                          downloadButton("downloadPlot3", "Download SVM score Boxplot")

                                        ), #sidebarPanel
                                        mainPanel( 
                                          tabsetPanel(type = "tabs",  
                                            tabPanel("Grad-Seq Heatmap", 
                                            uiOutput("GradSeqHeatmap_summaryUI") 
                                          ),
                                          tabPanel("Phylogeny Heatmap (Protein only)", 
                                             uiOutput("GradSeqPhylogeny_summaryUI") 
                                          ),
                                          tabPanel("SVM score Boxplot (Protein only)", 
                                            uiOutput("GradSeqSVM_summaryUI") 
                                          ),
                                          tabPanel("Tabular Output", 
                                                   uiOutput(outputId = "GradSeqTable_summaryUI")
                                          ) #tabPanel
                                        ) #tabsetPanel
                                        
                                      ) #mainPanel
                          ) #sidebarLayout

                          
                 ), # tabPanel
                 tabPanel("Gradient Composition",
                          sidebarLayout( position = "left",      
                                         sidebarPanel( 
                                           selectInput(inputId = "heightBar", label = "Please adjust Plot Height: ",
                                                       choices = c(1:10*10), selected=40),
                                           selectInput(inputId = "Bar", label = "Please set Bar Chart Type: ",
                                                       choices = c("Absolute Counts", "Relative Counts"), selected="Absolute Counts"),
                                           selectInput(inputId = "logtrans", label = "Scale: Normal / Log10", 
                                                       choices = c("Normal", "Log10"), selected = "Log10"),
                                           
                                          pickerInput(inputId = "BarFeature", label = "Please select Feature Type(s): ",
                                                      choices =  c("Protein", "3UTR", "5UTR", "mRNA", "sRNA", "tRNA", "rRNA", "aRNA", "crRNA", "Transposase"),
                                                      selected = c("Protein", "3UTR", "5UTR", "mRNA", "sRNA", "tRNA", "rRNA", "aRNA", "crRNA", "Transposase"),
                                                      options = list(`actions-box` = TRUE),multiple = T,
                                                      choicesOpt = list(style = c("background: lightgrey;",
                                                                                  "background: rgba(150, 0, 100, 0.3);",
                                                                                  "background: rgba(100, 0, 50, 0.3);",
                                                                                  "background: rgba(200, 0, 150, 0.1);",
                                                                                  "background: rgba(255, 0, 0, 0.5);",
                                                                                  "background: rgba(0, 255, 255, 0.3);",
                                                                                  "background: rgba(200, 0, 225, 0.5);",
                                                                                  "background: rgba(0, 0, 255, 0.5);", 
                                                                                  "background: rgba(255, 205, 0, 0.5);", 
                                                                                  "background: rgba(0, 255, 0, 0.5);")) ),
                                          
                                           pickerInput(inputId = "BarCluster", label = "Please select Cluster: ", 
                                                       choices =  c(0:17), selected= c(1:17), options = list(`actions-box` = TRUE),multiple = T),
                                           pickerInput(inputId = "BarMaximum", label = "Please select maximum Peak Fraction(s): ",
                                                       choices = c(1:13,14.5,16.5,18), selected = c(1:13,14.5,16.5,18), options = list(`actions-box` = TRUE),multiple = T),
                                           pickerInput(inputId = "BarPeak_no", label = "Please select number of Peaks: ",
                                                       choices = c(1:6), selected = c(1:6), options = list(`actions-box` = TRUE),multiple = T),
                                           
                                           pickerInput(inputId = "BarPeaks", label = "Please select Peak Fraction Combinations(s): ",
                                                       choices = c(1:13,14.5,16.5,18), selected = c(1:13,14.5,16.5,18), options = list(`actions-box` = TRUE),multiple = T),
                                           
                                           pickerInput(inputId = "BarLocation", label = "Please select Location(s): ",
                                                       choices = c("Chr", "pSYSA", "pSYSG", "pSYSM", "pSYSX", "pCA2.4"), 
                                                       selected = c("Chr", "pSYSA", "pSYSG", "pSYSM", "pSYSX", "pCA2.4"), options = list(`actions-box` = TRUE),multiple = T),

                                          sliderInput(inputId = "BarPlot_width", label = "Please select Plot width for download:",
                                                                     min = 1, max = 100, value = 20,  step = 5 ),
                                          downloadButton("downloadBarPlot", "Download Barplot of Gradient Composition")
                                                                            
                                         ), #sidebarPanel
                                         mainPanel( 
                                           uiOutput("GradSeqBarplot_summaryUI") 
                                         ) #mainPanel
                          ) #sidebarLayout
                 ), # tabPanel
                 tabPanel("Protein Categories",
                          sidebarLayout( position = "left",      
                                         sidebarPanel( 
                                           selectInput(inputId = "heightBar2", label = "Please adjust Plot Height: ",
                                                       choices = c(1:10*10), selected=40),
                                           selectInput(inputId = "Bar2", label = "Please set Bar Chart Type: ",
                                                       choices = c("Absolute Counts", "Relative Counts"), selected="Absolute Counts"),
                                           selectInput(inputId = "logtrans2", label = "Scale: Normal / Log10", 
                                                       choices = c("Normal", "Log10"), selected = "Log10"),
                              
                                           pickerInput(inputId = "database2", label = "Please select Functional Categories: ",
                                                       choices = unique(as.character(KEGG$database)), selected=c("syn00194 Photosynthesis Proteins"), 
                                                       options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                           pickerInput(inputId = "category2", choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T), 
                                           pickerInput(inputId = "pathway2",  choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                           pickerInput(inputId = "complex2",  choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                           
                                           pickerInput(inputId = "Bar2Cluster", label = "Please select Cluster: ", 
                                                       choices =  c(0:17), selected= c(1:17), options = list(`actions-box` = TRUE),multiple = T),
                                           pickerInput(inputId = "Bar2Maximum", label = "Please select maximum Peak Fraction(s): ",
                                                       choices = c(1:13,14.5,16.5,18), selected = c(1:13,14.5,16.5,18), options = list(`actions-box` = TRUE),multiple = T),
                                           pickerInput(inputId = "Bar2Peak_no", label = "Please select number of Peaks: ",
                                                       choices = c(1:6), selected = c(1:6), options = list(`actions-box` = TRUE),multiple = T),
                                           
                                           pickerInput(inputId = "Bar2Peaks", label = "Please select Peak Fraction Combinations(s): ",
                                                       choices = c(1:13,14.5,16.5,18), selected = c(1:13,14.5,16.5,18), options = list(`actions-box` = TRUE),multiple = T),
                                           
                                           pickerInput(inputId = "Bar2Location", label = "Please select Location(s): ",
                                                       choices = c("Chr", "pSYSA", "pSYSG", "pSYSM", "pSYSX", "pCA2.4"), 
                                                       selected = c("Chr", "pSYSA", "pSYSG", "pSYSM", "pSYSX", "pCA2.4"), options = list(`actions-box` = TRUE),multiple = T),
                                           
                                           
                                           pickerInput(inputId = "Bar2Phylogeny", label = "Please select Phylogenetic groups", 
                                                       choices = c("other Cyanobacteria", "other Cyanobacteria and A.thaliana", "other Cyanobacteria and S.enterica/E.coli", "All", "unknown"), 
                                                       selected = c("other Cyanobacteria", "other Cyanobacteria and A.thaliana", "other Cyanobacteria and S.enterica/E.coli", "All", "unknown"),
                                                       options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple = T),
                                           
                                           sliderInput(inputId = "Bar2Cyanobacteria", label = "Conservation within selected Cyanobacteria",
                                                       min = 0, max = 1, value = c(0,1),  step = 0.05 ),
                                           
                                           sliderInput(inputId = "Bar2SVM_score", label = "SVM score",
                                                       min = -16, max = 6, value = c(-16,6),  step = 0.1 ),
                                           
                                           sliderInput(inputId = "Bar2Plot_width", label = "Please select Plot width for download:",
                                                       min = 1, max = 100, value = 20,  step = 5 ),
                                           
                                           downloadButton("downloadBar2Plot", "Download Databse Barplot"),
                                           downloadButton("downloadBar2_cat_Plot", "Download Category Barplot"),
                                           downloadButton("downloadBar2_path_Plot", "Download Pathway Barplot"),
                                           downloadButton("downloadBar2_compl_Plot", "Download Complex / Type Barplot"),
                                           downloadButton("downloadBar2_cons_Plot", "Download Phylogenetic group Barplot"),
                                           downloadButton("downloadBar2_loc_Plot", "Download Location Barplot")
                                           
                                           
                                         ), #sidebarPanel
                                         mainPanel(
                                           tabsetPanel(type = "tabs",  
                                              tabPanel("KEGG Database",             
                                                uiOutput("GradSeqBar2plot_summaryUI") 
                                              ),
                                              tabPanel("Category",             
                                                       uiOutput("GradSeqBar2plot_cat_summaryUI") 
                                              ),
                                              tabPanel("Pathway",             
                                                       uiOutput("GradSeqBar2plot_path_summaryUI") 
                                              ),
                                              tabPanel("Complex / Type",             
                                                       uiOutput("GradSeqBar2plot_compl_summaryUI") 
                                              ),
                                              tabPanel("Phylogenetic group",
                                                       uiOutput("GradSeqBar2plot_cons_summaryUI")
                                              ),
                                              tabPanel("Location",
                                                       uiOutput("GradSeqBar2plot_loc_summaryUI")
                                              )
                                           ) # tabsetpanel
                                         ) #mainPanel
                          ) #sidebarLayout
                 ), 
                 tabPanel("Supplementary Information",
                          
                          mainPanel(
                               p(span(strong("Table S1: ")), "All genomes which were used in this study are shown in this table, based on the genbank/refeq files listed in the ", 
                              span(a("MBGD genomelist.", href="http://mbgd.genome.ad.jp/htbin/genomelist")), style = "font-family: 'times'; font-size: 16pt"),
                              br(),
                               uiOutput("GradSeqTable_mbgd_genomesUI"),
                              br(),
                              br(),
                              p(span(strong("Figure S1: ")), "SVM scores for proteins of ", span(em("Synechocystis ")), "sp. PCC6803, grouped by 
                              Gene Ontology Terms (Nucleotide binding, DNA binding, RNA binding) or 'No binding' if no such term existed. 
                              A SVM score < 0.49 is reached by the majority of proteins with no reported binding to nucleic acids, while almost all
                              proteins with reported RNA binding properties exceeded that threshold, showing the validity of the algorithm.
                              Therefore, all proteins of the dataset have been categorized into three groups in the SVM score Boxplots. Those proteins which fail
                              to pass a SVM score > 0.49 (grey group, probably no RNA binding), those proteins which are above that threshold but are still in the range of the upper quartile 
                              of the DNA binding group (yellow group, probably RNA binding) and those proteins which surpass the upper quartile of the DNA binding group with an SVM score >1.07 
                              (green group, RNA binding very likely).", style = "font-family: 'times'; font-size: 16pt"),
                              br(),
                              img(src='Figure_S2.jpg', height="20%", width="50%", align = "center")
                              
                          )
                )#tabpanel
                                           
)#navbarPage



server <- function(input, output, session){
  # Description:
  

 
  # Grad-Seq Explorer Tab outputs:
  HeatmapSize <- reactive({ req(input$height)
    as.numeric(input$height)
  })
  HeatmapHeight <- reactive( (20 * HeatmapSize()) )
  
  HeatmapRatio <- reactive({ req(input$heightRatio)
    as.numeric(input$heightRatio)
  })
  HeatmapHeightRatio <- reactive( (HeatmapRatio()) )
  
  Heatmap_Width <- reactive({ req(input$Plot_width)
    as.numeric(input$Plot_width)
  })
  HeatmapWidth <- reactive( (Heatmap_Width()) )
  
  Barplot_Width <- reactive({ req(input$BarPlot_width)
    as.numeric(input$BarPlot_width)
  })
  BarplotWidth <- reactive( (Barplot_Width()) )
  
  Bar2plot_Width <- reactive({ req(input$Bar2Plot_width)
    as.numeric(input$Bar2Plot_width)
  })
  Bar2plotWidth <- reactive( (Bar2plot_Width()) )
  
  

  
  database <- reactive({
    req(input$database)
    filter(KEGG, database %in% input$database)  
  })
  observeEvent(database(), { 
    updatePickerInput(session, "category", choices = unique(as.character(database()$Category))) 
  })
  category <- reactive({    req(input$category)
    filter(database(), Category %in% input$category)
  })
  observeEvent(category(), { 
    updatePickerInput(session, "pathway", choices = unique(as.character(category()$Pathway)))
  })
  pathway <- reactive({    req(input$pathway)
    filter(category(), Pathway %in% input$pathway)
  })
  observeEvent(pathway(), {
    updatePickerInput(session, "complex", choices = unique(as.character(pathway()$Complex.Type)))
  })
  complex <- reactive({    req(input$complex)
    filter(pathway(), Complex.Type %in% input$complex)
  })
  
  
  GradSeqHeatmap_Ind <- reactive({
    
    Gene_ID_in <- input$Gene_ID
    
    
    # Grad-Seq Heatmap Data.frame:
    data <- Shiny_dt
    
    if(input$NormChoice == "Constant Volume") {
    data[,6:21] <- as.data.frame( t( scale(t(data[, 6:21]), center = TRUE, scale = TRUE) ) )
    }
    
    data <- data %>% filter( Gene_ID_IndMode %in% Gene_ID_in ) 
    
    
    data <- cbind(data[,1:28], data[36])
    data <-  reshape2::melt(data, id.vars=c("Cluster", "Reproducibility",  "Feature.Type", "Feature.Subtype", "Amino.Acid.Sequence.Length",
                                  "Log2.Abundance.", "Maximum.Fraction", "Average.Fraction", "No.of.peaks", "Peak.fractions", "Average.Spearmann.Correlation.Coefficient.between.Replicates", 
                                  "PearsonCorrCoeff_to_ME", "Gene_ID_IndMode"))
    
    if(input$NormChoice == "Constant Protein") {
      data <- merge(data, constProt, by=c("variable", "Feature.Type") )
      data$value <- data$value * data$Factor
      data <- dcast(data, Cluster+Reproducibility+Feature.Type+Feature.Subtype+Amino.Acid.Sequence.Length+
                        Log2.Abundance.+Maximum.Fraction+Average.Fraction+No.of.peaks+Peak.fractions+Average.Spearmann.Correlation.Coefficient.between.Replicates+ 
                        PearsonCorrCoeff_to_ME+Gene_ID_IndMode ~ variable, FUN=sum)
      data[,14:29] <- as.data.frame( t( scale(t(data[,14:29]), center = TRUE, scale = TRUE) ) )
    data <-  reshape2::melt(data, id.vars=c("Cluster", "Reproducibility",  "Feature.Type", "Feature.Subtype", "Amino.Acid.Sequence.Length",
                                  "Log2.Abundance.", "Maximum.Fraction", "Average.Fraction", "No.of.peaks", "Peak.fractions", "Average.Spearmann.Correlation.Coefficient.between.Replicates", 
                                  "PearsonCorrCoeff_to_ME", "Gene_ID_IndMode"))
    }
    
    data_dedup <- data[!duplicated(data$Gene_ID_IndMode),]
    
     
    
    heat_ind <- ggplot(data, aes(x = variable, y = Gene_ID_IndMode, fill = as.factor(Feature.Subtype) , alpha = as.numeric(value)) ) + 
      geom_tile() +
      #scale_fill_gradient2("z-score",low = "white", mid="white", high="black", midpoint = 0) +
      scale_alpha_continuous(range = c(0, 1), limits=c(-0.25,3), oob = squish) +
      scale_fill_manual("Feature Type",values = c("3UTR" = "pink3", "5UTR" = "pink4", "aRNA" ="blue", "crRNA" = "gold2", 
                                                  "mRNA"="pink", "sRNA" = "red", "Protein" = "black", "rRNA" ="purple", 
                                                  "Transposase" = "green2", "tRNA" = "turquoise1", "transcript" ="orange"), guide = guide_colorbar(barwidth = 0.8, barheight = 18)) +
      scale_y_discrete(name = "Locustags",
                       limits=data_dedup$Gene_ID_IndMode[ order( -as.numeric(data_dedup$Cluster), -as.numeric(data_dedup$Maximum.Fraction), -as.numeric(data_dedup$Average.Fraction) ) ] ) +
      scale_x_discrete(name = "Fractions") +
      #geom_segment( aes(x=17, xend=17, y = 0, yend = 200  ),color="black", size=8) +
      theme_bw() 
    if(input$HeatmapChoice == "Locustags") {
      heat_ind + labs(alpha = "z-score (<0 = white)   ", fill="Feature Type") + guides(alpha = guide_legend(label.position = "bottom"), fill = guide_legend(label.position = "right")) + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing.y = unit(-0.1, "lines"), 
                       legend.background = element_rect(linetype = 2, size = 0.5, colour = 1), legend.position = "bottom") 
    } else if(input$HeatmapChoice == "Both") {
    heat_ind + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing.y = unit(-0.1, "lines"), 
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1), legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank() ) 
    }  
        
  })
  GradSeqHeatmap <- reactive({
    
    req(input$complex)
    
    Cluster_in <- input$Cluster
    Feature_in <- input$Feature
    Max_Peak <- input$Maximum
    No_Peak <- input$Peak_no
    Peak_in <- input$Peaks
    
    Location_in <- input$Location
    Phylogeny_in <- input$Phylogeny
    Cyanos_in_low <- input$Cyanobacteria[1]
    Cyanos_in_high <- input$Cyanobacteria[2]
    
    SVM_low <- input$SVM_score[1]
    SVM_high <- input$SVM_score[2]
    
    SVM_median_low <- input$Median_SVM_score[1]
    SVM_median_high <- input$Median_SVM_score[2]
    
    AA_length_low <- input$AA_LENGTH[1]
    AA_length_high <- input$AA_LENGTH[2]
    
    #phylogeny_detail <- input$Phylogeny_Detail
    
     
    
    # Grad-Seq Heatmap Data.frame:
    data <- complex()
    
    
    #data <- droplevels(subset(KEGG, Complex.Type == "Photosystem II (P680 chlorophyll a)" & Feature.Type == "Protein"))
    #data<- KEGG
    
    if(input$NormChoice == "Constant Volume") {
    data[,19:34] <- as.data.frame( t( scale(t(data[, 19:34]), center = TRUE, scale = TRUE) ) )
    }
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    #phylo_org <- mbgd %>% filter(value > 0, variable %in% phylogeny_detail)
    #phylo_org <- mbgd$syn_locustag
    
    protein_data <- data %>% filter(Feature.Type == "Protein",
                                    Cluster %in% Cluster_in,
                                    Feature.Subtype %in% Feature_in,
                                    Maximum.Fraction %in% Max_Peak,
                                    No.of.peaks %in% No_Peak,
                                    location %in% Location_in,
                                    Conservation %in% Phylogeny_in,
                                    Cyanobacteria >= Cyanos_in_low,
                                    Cyanobacteria <= Cyanos_in_high,
                                    SVM.score >= SVM_low  | SVM.score == "NA",
                                    SVM.score <= SVM_high  | SVM.score == "NA",
                                    SVM.score.Median >= SVM_median_low,
                                    SVM.score.Median <= SVM_median_high,
                                    Amino.Acid.Sequence.Length >= AA_length_low,
                                    Amino.Acid.Sequence.Length <= AA_length_high) #, locustag %in% phylo_org)
    
    RNA_data <- data %>% filter(Feature.Type != "Protein",
                                Cluster %in% Cluster_in,
                                Feature.Subtype %in% Feature_in,
                                Maximum.Fraction %in% Max_Peak,
                                No.of.peaks %in% No_Peak,
                                location %in% Location_in) 
    data <- rbind(protein_data, RNA_data)
    
    
  
    data <- data[,1:48]
    data <-  reshape2::melt(data, id.vars=c("locustag","genename","annotation","KEGG.Orthology","database","MainCategory","MainPathway",
                                  "Category","undefined","Pathway","Complex.Type","Details","Summary","Cluster","Reproducibility",
                                  "Log2.Abundance.","Average.Spearmann.Correlation.Coefficient.between.Replicates",
                                  "PearsonCorrCoeff_to_ME","Maximum.Fraction","Average.Fraction","No.of.peaks","Peak.fractions",
                                  "Feature.Type","Feature.Subtype","Amino.Acid.Sequence.Length","locustag_ID","genename_1","TU",
                                  "annotation..cyanobase.","GO.NAME","GOT_DNA_RNA_Nucleotide","Gene_ID_IndMode"))     
    
    if(input$NormChoice == "Constant Protein") {
      data <- merge(data, constProt, by=c("variable", "Feature.Type") )
      data$value <- data$value * data$Factor
      data <- dcast(data, locustag+genename+annotation+KEGG.Orthology+database+MainCategory+MainPathway+
                        Category+undefined+Pathway+Complex.Type+Details+Summary+Cluster+Reproducibility+
                        Log2.Abundance.+Average.Spearmann.Correlation.Coefficient.between.Replicates+
                        PearsonCorrCoeff_to_ME+Maximum.Fraction+Average.Fraction+No.of.peaks+Peak.fractions+
                        Feature.Type+Feature.Subtype+Amino.Acid.Sequence.Length+locustag_ID+genename_1+TU+
                        annotation..cyanobase.+GO.NAME+GOT_DNA_RNA_Nucleotide+Gene_ID_IndMode ~ variable, FUN=sum)
      data[,33:48] <- as.data.frame( t( scale(t(data[,33:48]), center = TRUE, scale = TRUE) ) )
      data <-   reshape2::melt(data, id.vars=c("locustag","genename","annotation","KEGG.Orthology","database","MainCategory","MainPathway",
                                                   "Category","undefined","Pathway","Complex.Type","Details","Summary","Cluster","Reproducibility",
                                                   "Log2.Abundance.","Average.Spearmann.Correlation.Coefficient.between.Replicates",
                                                   "PearsonCorrCoeff_to_ME","Maximum.Fraction","Average.Fraction","No.of.peaks","Peak.fractions",
                                                   "Feature.Type","Feature.Subtype","Amino.Acid.Sequence.Length","locustag_ID","genename_1","TU",
                                                   "annotation..cyanobase.","GO.NAME","GOT_DNA_RNA_Nucleotide","Gene_ID_IndMode"))  
    }
    
    
    data_dedup <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    HeatOrder <- switch(input$HeatmapOrder,
                        "Cluster" = -as.numeric(data_dedup$Cluster), 
                        "Feature Type"  = data_dedup$Feature.Subtype, 
                        "Category" = data_dedup$Summary,
                        "Locustag" = data_dedup$locustag_ID,
                        "Abundance" = data_dedup$Log2.Abundance.)
    
    
    ggplot(data, aes(x = variable, y = Gene_ID_IndMode, fill = as.factor(Feature.Subtype) , alpha = as.numeric(value)) ) + 
      geom_tile() +
      #scale_fill_gradient2("z-score",low = "white", mid="white", high="black", midpoint = 0) +
      scale_alpha_continuous(range = c(0, 1), limits=c(-0.25,3), oob = squish ) +
      scale_fill_manual("Feature Type",values = c("3UTR" = "pink3", "5UTR" = "pink4", "aRNA" ="blue", "crRNA" = "gold2", 
                                             "mRNA"="pink", "sRNA" = "red", "Protein" = "black", "rRNA" ="purple", 
                                             "Transposase" = "green2", "tRNA" = "turquoise1", "transcript" ="orange")) +
      scale_y_discrete(name = "Pathways",
                       limits=data_dedup$Gene_ID_IndMode[ order(HeatOrder, -as.numeric(data_dedup$Cluster), -as.numeric(data_dedup$Maximum.Fraction), -as.numeric(data_dedup$Average.Fraction) ) ] ) +
      scale_x_discrete(name = "Fractions") +
      #geom_segment( aes(x=17, xend=17, y = 0, yend = 200  ),color="black", size=8) +
      theme_bw() + labs(alpha = "z-score (<0 = white)   ", fill="Feature Type") + guides(alpha = guide_legend(label.position = "bottom"), fill = guide_legend(label.position = "right")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing.y = unit(-0.1, "lines"), #legend.position = c(0.8,0.8),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1), legend.position = "bottom")  #+ eval(call)
    
    
    
    
    
    
    
    
  })
  output$GradSeqHeatmap_summary <- renderPlot({ 
    if(input$HeatmapChoice == "Locustags"){
      GradSeqHeatmap_Ind()
    }  
    else if(input$HeatmapChoice == "Pathways"){
      GradSeqHeatmap()
    }  
    else if(input$HeatmapChoice == "Both"){
      plot_grid(GradSeqHeatmap_Ind(), GradSeqHeatmap(), ncol=1, align="v",  rel_heights = c(1,HeatmapHeightRatio()) )
    }  
  })
  output$GradSeqHeatmap_summaryUI <- renderUI({
    plotOutput("GradSeqHeatmap_summary", height = HeatmapHeight())
  })
  
  
  GradSeqPhylogeny_Ind <- reactive({
    
    Gene_ID_in <- input$Gene_ID
    data <- Shiny_dt
    data <- data %>% filter( Gene_ID_IndMode %in% Gene_ID_in ) 
    
    data <- cbind(data[,1:5], data[,22:36])
    data <- droplevels(subset(data, Feature.Type == "Protein"))
    
    cons_order <- droplevels(subset(mbgd, syn_locustag == "Order"))
    
    names(mbgd)[1] <- "locustag"
    mbgd_in <- merge(data, mbgd, by=c("locustag"), all.x = TRUE)
    mbgd_in <- mbgd_in[!duplicated(mbgd_in),]
    
    data_dedup <- mbgd_in[!duplicated(mbgd_in$Gene_ID_IndMode),]
    
    
    
    HeatOrder_x <- -as.numeric(cons_order$value)
    
    
    heat_ind <- ggplot(mbgd_in, aes(x = variable, y = Gene_ID_IndMode, fill = as.factor(value))) + 
      geom_tile( alpha=0.7 ) +
      #scale_fill_gradient2("Conservation",low = "grey50", mid="yellow4", high="green", midpoint = 0.3) +
      #scale_fill_manual(values =c("0" = "red", "0.33333" = "yellow1", "0.5" ="yellow3", "0.6" = "greenyellow", "0.75" = "greenyellow", "0.8" = "greenyellow", "1" = "green"))+
      scale_fill_manual(values = c("0" = "red","0.2" = "yellow1",  "0.25" = "yellow1", "0.33" = "yellow1", 
                                   "0.4" = "yellow2", "0.5" = "yellow2", 
                                   "0.6" = "yellow3", "0.67" = "yellow3", 
                                   "0.75" = "greenyellow", "0.8" = "greenyellow", 
                                   "1" = "green", "-1" = "grey50" ) ) +
      scale_y_discrete(name = "Locustags",
                       limits=data_dedup$Gene_ID_IndMode[ order( -as.numeric(data_dedup$Cluster), -as.numeric(data_dedup$Maximum.Fraction), -as.numeric(data_dedup$Average.Fraction) ) ] ) +
      scale_x_discrete(name = "Organisms", 
                       limits=cons_order$variable[ order(-HeatOrder_x, as.numeric(cons_order$value))  ] ) +
      #geom_segment( aes(x=17, xend=17, y = 0, yend = 200  ),color="black", size=8) +
      theme_bw() 
    
    if(input$HeatmapChoice == "Locustags") {
    heat_ind  + labs(fill="Relative protein conservation (-1 = NA)    " ) + guides(fill = guide_legend(label.position = "bottom", nrow = 1) ) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing.y = unit(-0.1, "lines"), #legend.position = c(0.8,0.8),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1), legend.position = "bottom", axis.text.x = element_text(angle = 45, size=8, hjust = 1) )
    } else if(input$HeatmapChoice == "Both") {      
    heat_ind + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing.y = unit(-0.1, "lines"), #legend.position = c(0.8,0.8),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1), legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank() )  
      
    }
    
  })
  GradSeqPhylogeny <- reactive({
    
    req(input$complex)
    
    Cluster_in <- input$Cluster
    Feature_in <- input$Feature
    Max_Peak <- input$Maximum
    No_Peak <- input$Peak_no
    Peak_in <- input$Peaks
    
    Location_in <- input$Location
    Phylogeny_in <- input$Phylogeny
    Cyanos_in_low <- input$Cyanobacteria[1]
    Cyanos_in_high <- input$Cyanobacteria[2]
    
    SVM_low <- input$SVM_score[1]
    SVM_high <- input$SVM_score[2]
    
    SVM_median_low <- input$Median_SVM_score[1]
    SVM_median_high <- input$Median_SVM_score[2]
    
    AA_length_low <- input$AA_LENGTH[1]
    AA_length_high <- input$AA_LENGTH[2]
    #phylogeny_detail <- input$Phylogeny_Detail
    
    
    # Grad-Seq Heatmap Data.frame:
    data <- complex()
    #data <- KEGG
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    #phylo_org <- mbgd %>% filter(value > 0, variable %in% phylogeny_detail)
    #phylo_org <- mbgd$syn_locustag
    
    protein_data <- data %>% filter(Feature.Type == "Protein",
                                    Cluster %in% Cluster_in,
                                    Feature.Subtype %in% Feature_in,
                                    Maximum.Fraction %in% Max_Peak,
                                    No.of.peaks %in% No_Peak,
                                    location %in% Location_in,
                                    Conservation %in% Phylogeny_in,
                                    Cyanobacteria >= Cyanos_in_low,
                                    Cyanobacteria <= Cyanos_in_high,
                                    SVM.score >= SVM_low  | SVM.score == "NA",
                                    SVM.score <= SVM_high  | SVM.score == "NA",
                                    SVM.score.Median >= SVM_median_low,
                                    SVM.score.Median <= SVM_median_high,
                                    Amino.Acid.Sequence.Length >= AA_length_low,
                                    Amino.Acid.Sequence.Length <= AA_length_high) #, locustag %in% phylo_org) 
    
    data <- protein_data
    
 
    data <- data[,1:48]
    data <- cbind(data[1:18], data[35:48])
    data <- droplevels(subset(data, Feature.Type == "Protein"))
    
    cons_order <- droplevels(subset(mbgd, syn_locustag == "Order"))
    
    names(mbgd)[1] <- "locustag"
    mbgd_in <- merge(data, mbgd, by=c("locustag"), all.x = TRUE)
    mbgd_in <- mbgd_in[!duplicated(mbgd_in),]

    data_dedup <- mbgd_in[!duplicated(mbgd_in$Gene_ID_IndMode),]
    
    
    
    HeatOrder_x <- -as.numeric(cons_order$value)
    
    HeatOrder <- switch(input$HeatmapOrder,
                        "Cluster" = -as.numeric(data_dedup$Cluster), 
                        "Feature Type"  = data_dedup$Feature.Subtype, 
                        "Category" = data_dedup$Summary,
                        "Locustag" = data_dedup$locustag_ID,
                        "Abundance" = data_dedup$Log2.Abundance.)
   
 
    ggplot(mbgd_in, aes(x = variable, y = Gene_ID_IndMode, fill = as.factor(value))) + 
      geom_tile( alpha=0.7 ) +
      #scale_fill_gradient2("Conservation",low = "grey50", mid="yellow4", high="green", midpoint = 0.3) +
      #scale_fill_manual(values =c("0" = "red", "0.33333" = "yellow1", "0.5" ="yellow3", "0.6" = "greenyellow", "0.75" = "greenyellow", "0.8" = "greenyellow", "1" = "green"))+
      scale_fill_manual(values = c("0" = "red","0.2" = "yellow1",  "0.25" = "yellow1", "0.33" = "yellow1", 
                                   "0.4" = "yellow2", "0.5" = "yellow2", 
                                   "0.6" = "yellow3", "0.67" = "yellow3", 
                                   "0.75" = "greenyellow", "0.8" = "greenyellow", 
                                   "1" = "green", "-1" = "grey50" ) ) +
      scale_y_discrete(name = "Pathways",
                       limits=data_dedup$Gene_ID_IndMode[ order(HeatOrder, -as.numeric(data_dedup$Cluster), -as.numeric(data_dedup$Maximum.Fraction), -as.numeric(data_dedup$Average.Fraction) ) ] ) +
      scale_x_discrete(name = "Organisms", 
                       limits=cons_order$variable[ order(-HeatOrder_x, as.numeric(cons_order$value))  ] ) +
      #geom_segment( aes(x=17, xend=17, y = 0, yend = 200  ),color="black", size=8) +
      theme_bw() + labs(fill="Relative protein conservation (-1 = NA)    " ) + guides(fill = guide_legend(label.position = "bottom", nrow = 1) ) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing.y = unit(-0.1, "lines"), #legend.position = c(0.8,0.8),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1), legend.position = "bottom", axis.text.x = element_text(angle = 45, size=8, hjust = 1))  

    
  })
  output$GradSeqPhylogeny_summary <- renderPlot({ 
    if(input$HeatmapChoice == "Locustags"){
      GradSeqPhylogeny_Ind()
    }  
    else if(input$HeatmapChoice == "Pathways"){
      GradSeqPhylogeny()
    }  
    else if(input$HeatmapChoice == "Both"){
      plot_grid(GradSeqPhylogeny_Ind(), GradSeqPhylogeny(), ncol=1, align="v",  rel_heights = c(1,HeatmapHeightRatio()) )
    }  
  })
  output$GradSeqPhylogeny_summaryUI <- renderUI({
    plotOutput("GradSeqPhylogeny_summary", height = HeatmapHeight())
  })
  
  
  grob <- grobTree(textGrob("Synechocystis 6803", x=0.1,  y=-0.1, hjust=0,
                            gp=gpar(col="red", fontsize=12, fontface="italic")))
  
  GradSeqSVM_Ind <- reactive({
    
   
    Gene_ID_in <- input$Gene_ID
    data <- Shiny_dt
    data <- data %>% filter( Gene_ID_IndMode %in% Gene_ID_in ) 
    
    data <- cbind(data[,1:5], data[,22:36])
    data <- droplevels(subset(data, Feature.Type == "Protein"))
    
    
    syn_svm <- SVM %>% filter(locus_tags %in% data$locustag) %>% 
      arrange(-Ortholog_clustersize) %>% 
      distinct(locus_tags, .keep_all = TRUE)
    syn_svm <- cbind(syn_svm[,1:2], syn_svm[6])
    names(syn_svm)[1:3] <- c("locustag", "Ortholog_ClusterID", "syn_svm")
    unique_syn_svm <- syn_svm$Ortholog_ClusterID[!duplicated(syn_svm$Ortholog_ClusterID)]
    svm_input <- SVM %>% filter(Ortholog_ClusterID %in% unique_syn_svm)
    svm_input <- merge(svm_input, syn_svm, by=c("Ortholog_ClusterID"), all = TRUE) 
    
    svm_median <- aggregate(svm_input[6], list(svm_input$locustag), median, na.rm=TRUE)
    names(svm_median) <- c("locustag", "SVM.score.Median")
    svm_input <- merge(svm_input, svm_median, by=c("locustag"), all = TRUE) 
    
    svm_input <- merge(svm_input, data, by=c("locustag"))
    svm_input <- droplevels(subset(svm_input, Feature.Type == "Protein"))
    
    if(input$Median == "Synechocystis 6803") {
    svm_input = within(svm_input, {
      svm_group = ifelse(syn_svm <= 0.49, "no_RNA_binding", ifelse(syn_svm > 1.07, "likely_RNA_binding", "potential_RNA_binding" ) )
    })
    }
    if(input$Median == "Median") {
      svm_input = within(svm_input, {
        svm_group = ifelse(SVM.score.Median <= 0.49, "no_RNA_binding", ifelse(SVM.score.Median > 1.07, "likely_RNA_binding", "potential_RNA_binding" ) )
      })
    }
    
    data_dedup <- svm_input[!duplicated(svm_input$Gene_ID_IndMode),]
     
    svm_ind <- ggplot(data = svm_input, aes(x = Gene_ID_IndMode, y = SVM.score, fill = svm_group))+
      geom_boxplot()+
      scale_fill_manual(values = c("no_RNA_binding" = "grey70", "likely_RNA_binding" = "green2", "potential_RNA_binding" = "yellow") )+
      geom_point( aes(x=Gene_ID_IndMode, y=syn_svm), color="red", size=2, shape=19) + 
      #annotation_custom(grob)+
      scale_x_discrete(name = "Locustags",
                       limits=data_dedup$Gene_ID_IndMode[ order(-as.numeric(data_dedup$Cluster), -as.numeric(data_dedup$Maximum.Fraction), -as.numeric(data_dedup$Average.Fraction) ) ] ) +
      geom_hline(yintercept = 0.49, linetype = "dashed", color = "yellow4", size= 1.5)+
      geom_hline(yintercept = 1.07, linetype = "dashed", color = "green2", size= 1.5)+
      scale_y_continuous(name = "SVM score")+
      coord_flip()+
      theme_bw()
    
      if(input$HeatmapChoice == "Locustags") {
      svm_ind + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position = "bottom")
      } else if(input$HeatmapChoice == "Both") {
      svm_ind + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position = "none", axis.title.x=element_blank() )
      }
    
  })
  GradSeqSVM <- reactive({
    
    req(input$complex)
    
    Cluster_in <- input$Cluster
    Feature_in <- input$Feature
    Max_Peak <- input$Maximum
    No_Peak <- input$Peak_no
    Peak_in <- input$Peaks
    
    Location_in <- input$Location
    Phylogeny_in <- input$Phylogeny
    Cyanos_in_low <- input$Cyanobacteria[1]
    Cyanos_in_high <- input$Cyanobacteria[2]
    
    SVM_low <- input$SVM_score[1]
    SVM_high <- input$SVM_score[2]
    
    SVM_median_low <- input$Median_SVM_score[1]
    SVM_median_high <- input$Median_SVM_score[2]
    
    AA_length_low <- input$AA_LENGTH[1]
    AA_length_high <- input$AA_LENGTH[2]
    #phylogeny_detail <- input$Phylogeny_Detail
    
    
    # Grad-Seq Heatmap Data.frame:
    data <- complex()
    #data <- KEGG
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    #phylo_org <- mbgd %>% filter(value > 0, variable %in% phylogeny_detail)
    #phylo_org <- mbgd$syn_locustag
    
    protein_data <- data %>% filter(Feature.Type == "Protein",
                                    Cluster %in% Cluster_in,
                                    Feature.Subtype %in% Feature_in,
                                    Maximum.Fraction %in% Max_Peak,
                                    No.of.peaks %in% No_Peak,
                                    location %in% Location_in,
                                    Conservation %in% Phylogeny_in,
                                    Cyanobacteria >= Cyanos_in_low,
                                    Cyanobacteria <= Cyanos_in_high,
                                    SVM.score >= SVM_low  | SVM.score == "NA",
                                    SVM.score <= SVM_high  | SVM.score == "NA",
                                    SVM.score.Median >= SVM_median_low,
                                    SVM.score.Median <= SVM_median_high,
                                    Amino.Acid.Sequence.Length >= AA_length_low,
                                    Amino.Acid.Sequence.Length <= AA_length_high )#, locustag %in% phylo_org) 
    
    data <- protein_data
    
    
       
    
    
    data <- data[,1:48]
    data <- cbind(data[1:18], data[35:48])
    data <- droplevels(subset(data, Feature.Type == "Protein"))
    
    syn_svm <- SVM %>% filter(locus_tags %in% data$locustag) %>% 
      arrange(-Ortholog_clustersize) %>% 
      distinct(locus_tags, .keep_all = TRUE)
    syn_svm <- cbind(syn_svm[,1:2], syn_svm[6])
    names(syn_svm)[1:3] <- c("locustag", "Ortholog_ClusterID", "syn_svm")
    unique_syn_svm <- syn_svm$Ortholog_ClusterID[!duplicated(syn_svm$Ortholog_ClusterID)]
    svm_input <- SVM %>% filter(Ortholog_ClusterID %in% unique_syn_svm)
    svm_input <- merge(svm_input, syn_svm, by=c("Ortholog_ClusterID"), all = TRUE) 
    
    svm_median <- aggregate(svm_input[6], list(svm_input$locustag), median, na.rm=TRUE)
    names(svm_median) <- c("locustag", "SVM.score.Median")
    svm_input <- merge(svm_input, svm_median, by=c("locustag"), all = TRUE) 
    
    svm_input <- merge(svm_input, data, by=c("locustag"))
    svm_input <- droplevels(subset(svm_input, Feature.Type == "Protein"))
    
    if(input$Median == "Synechocystis 6803") {
      svm_input = within(svm_input, {
        svm_group = ifelse(syn_svm <= 0.49, "no_RNA_binding", ifelse(syn_svm > 1.07, "likely_RNA_binding", "potential_RNA_binding" ) )
      })
    }
    if(input$Median == "Median") {
      svm_input = within(svm_input, {
        svm_group = ifelse(SVM.score.Median <= 0.49, "no_RNA_binding", ifelse(SVM.score.Median > 1.07, "likely_RNA_binding", "potential_RNA_binding" ) )
      })
    }
    
    
    
    data_dedup <- svm_input[!duplicated(svm_input$Gene_ID_IndMode),]
    
    HeatOrder <- switch(input$HeatmapOrder,
                        "Cluster" = -as.numeric(data_dedup$Cluster), 
                        "Feature Type"  = data_dedup$Feature.Subtype, 
                        "Category" = data_dedup$Summary,
                        "Locustag" = data_dedup$locustag_ID,
                        "Abundance" = data_dedup$Log2.Abundance.)
    
    ggplot(data = svm_input, aes(x=Gene_ID_IndMode, y=SVM.score, fill = svm_group))+
      geom_boxplot()+
      scale_fill_manual(values = c("no_RNA_binding" = "grey70", "likely_RNA_binding" = "green2", "potential_RNA_binding" = "yellow") )+
      geom_point( aes(x=Gene_ID_IndMode, y=syn_svm), color="red", size=2, shape=19) + 
      #annotation_custom(grob)+
      scale_x_discrete(name = "Pathways",
                       limits=data_dedup$Gene_ID_IndMode[ order(HeatOrder, -as.numeric(data_dedup$Cluster), -as.numeric(data_dedup$Maximum.Fraction), -as.numeric(data_dedup$Average.Fraction) ) ] ) +
      geom_hline(yintercept = 0.49, linetype = "dashed", color = "yellow4", size= 1.5)+
      geom_hline(yintercept = 1.07, linetype = "dashed", color = "green2", size= 1.5)+
      scale_y_continuous(name = "SVM score")+
      coord_flip()+
      theme_bw()+
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position = "bottom")
    
  })
  output$GradSeqSVM_summary <- renderPlot({ 
    if(input$HeatmapChoice == "Locustags"){
      GradSeqSVM_Ind()
    }  
    else if(input$HeatmapChoice == "Pathways"){
      GradSeqSVM()
    }  
    else if(input$HeatmapChoice == "Both"){
      plot_grid(GradSeqSVM_Ind(), GradSeqSVM(), ncol=1, align="v",  rel_heights = c(1,HeatmapHeightRatio()) )
    }  
  })
  output$GradSeqSVM_summaryUI <- renderUI({
    plotOutput("GradSeqSVM_summary", height = HeatmapHeight())
  })
  
  
  output$GradSeqTable_summary <- DT::renderDataTable({ 
    if(input$HeatmapChoice == "Locustags"){
      Gene_ID_in <- input$Gene_ID
      data_Ind <- KEGG
      data_Ind <- data_Ind %>% filter( Gene_ID_IndMode %in% Gene_ID_in ) 
      data_Ind <- data_Ind[order(data_Ind$Gene_ID_IndMode),]
      
      # merge the categories into one
      database <- cbind(data_Ind[48], data_Ind[5])
      database <- database[!duplicated(database),]
      database <- ddply(database, .(Gene_ID_IndMode), summarise, database = list(as.character(database)))
      
      category <- cbind(data_Ind[48], data_Ind[8])
      category <- category[!duplicated(category),]
      category <- ddply(category, .(Gene_ID_IndMode), summarise, category = list(as.character(Category)))
      
      pathway <- cbind(data_Ind[48], data_Ind[10])
      pathway <- pathway[!duplicated(pathway),]
      pathway <- ddply(pathway, .(Gene_ID_IndMode), summarise, pathway = list(as.character(Pathway)))
      
      complex <- cbind(data_Ind[48], data_Ind[11])
      complex <- complex[!duplicated(complex),]
      complex <- ddply(complex, .(Gene_ID_IndMode), summarise, complex = list(as.character(Complex.Type)))
      
      details <- cbind(data_Ind[48], data_Ind[12])
      details <- details[!duplicated(details),]
      details <- ddply(details, .(Gene_ID_IndMode), summarise, details = list(as.character(Details)))
      
      data_Ind <- data_Ind[!duplicated(data_Ind$Gene_ID_IndMode),]
     
      data_Ind <- dplyr::left_join(data_Ind, database, by="Gene_ID_IndMode")
      data_Ind <- dplyr::left_join(data_Ind, category, by="Gene_ID_IndMode")
      data_Ind <- dplyr::left_join(data_Ind, pathway, by="Gene_ID_IndMode")
      data_Ind <- dplyr::left_join(data_Ind, complex, by="Gene_ID_IndMode")
      data_Ind <- dplyr::left_join(data_Ind, details, by="Gene_ID_IndMode")
      data_Ind$database.x <- data_Ind$database.y;  names(data_Ind)[5] <- "database" 
      data_Ind$Category <- data_Ind$category
      data_Ind$Pathway <- data_Ind$pathway
      data_Ind$Complex.Type <- data_Ind$complex
      data_Ind$Details <- data_Ind$details
      
      
      
      
      data_Ind <- cbind(data_Ind[,1:18], data_Ind[,35:52])
      df_Ind <- as.data.frame(cbind(data_Ind[,26:27],data_Ind[24],data_Ind[,28:29],data_Ind[33],data_Ind[,14:16],data_Ind[,19:22],data_Ind[,4:5],data_Ind[8],data_Ind[,10:12],
                                    data_Ind[31],data_Ind[25],data_Ind[,34:36]))

            
      if(nrow(df_Ind[which(df_Ind$Feature.Subtype != "Protein"),]) > 0) {
        df_Ind_RNA <- df_Ind %>% filter(Feature.Subtype != "Protein") %>% mutate(SVM.score = "NA")
        df_Ind_Protein <- df_Ind %>% filter(Feature.Subtype == "Protein") 
        df_Ind <- rbind(df_Ind_Protein, df_Ind_RNA)
      }
         
      if(nrow(df_Ind[which(df_Ind$Conservation == "unknown"),]) > 0) {
        df_Ind_unknown <- df_Ind %>% filter(Conservation == "unknown") %>% mutate(Cyanobacteria = "NA")
        df_Ind_known <- df_Ind %>% filter(Conservation != "unknown") 
        df_Ind <- rbind(df_Ind_known, df_Ind_unknown) 
      }
      
      names(df_Ind) <- c("Locustag", "Genename", "Feature Type", "TU", "Annotation", "Location", "Cluster", "Reproducibility", "Log2(Abundance)", "Maximum Fraction", 
                     "Average Fraction", "No of Peaks", "Peak Fractions", "KEGG Orthology", "Database", "Category", "Pathway", "Complex/Type", "Details",
                     "Gene Ontology Nucleic acid binding", "AA Length", "Relative Conservation in Cyanobacteria", "Phylogenetic Group", "SVM Score")
      # Data.frame output:  
      df_Ind <- as_tibble(df_Ind) 
    }  
    else if(input$HeatmapChoice == "Pathways"){
      req(input$complex)
      
      Cluster_in <- input$Cluster
      Feature_in <- input$Feature
      Max_Peak <- input$Maximum
      No_Peak <- input$Peak_no
      Peak_in <- input$Peaks
      
      Location_in <- input$Location
      Phylogeny_in <- input$Phylogeny
      Cyanos_in_low <- input$Cyanobacteria[1]
      Cyanos_in_high <- input$Cyanobacteria[2]
      
      SVM_low <- input$SVM_score[1]
      SVM_high <- input$SVM_score[2]
      
      SVM_median_low <- input$Median_SVM_score[1]
      SVM_median_high <- input$Median_SVM_score[2]
      
      
      data <- complex()
      #data<- KEGG
      data[,19:34] <- as.data.frame( t( scale(t(data[, 19:34]), center = TRUE, scale = TRUE) ) )
      
      data <- transform(data, All_peaks = Peak.fractions)
      data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
      data <- data[!duplicated(data$Gene_ID_IndMode),]
      
      protein_data <- data %>% filter(Feature.Type == "Protein",
                                      Cluster %in% Cluster_in,
                                      Feature.Subtype %in% Feature_in,
                                      Maximum.Fraction %in% Max_Peak,
                                      No.of.peaks %in% No_Peak,
                                      location %in% Location_in,
                                      Conservation %in% Phylogeny_in,
                                      Cyanobacteria >= Cyanos_in_low,
                                      Cyanobacteria <= Cyanos_in_high,
                                      SVM.score >= SVM_low  | SVM.score == "NA",
                                      SVM.score <= SVM_high  | SVM.score == "NA",
                                      SVM.score.Median >= SVM_median_low,
                                      SVM.score.Median <= SVM_median_high) 
      
      RNA_data <- data %>% filter(Feature.Type != "Protein",
                                  Cluster %in% Cluster_in,
                                  Feature.Subtype %in% Feature_in,
                                  Maximum.Fraction %in% Max_Peak,
                                  No.of.peaks %in% No_Peak) 
      data <- rbind(protein_data, RNA_data)
      
      
      # merge the categories into one
      database <- cbind(data[48], data[5])
      database <- database[!duplicated(database),]
      database <- ddply(database, .(Gene_ID_IndMode), summarise, database = list(as.character(database)))
      
      category <- cbind(data[48], data[8])
      category <- category[!duplicated(category),]
      category <- ddply(category, .(Gene_ID_IndMode), summarise, category = list(as.character(Category)))
      
      pathway <- cbind(data[48], data[10])
      pathway <- pathway[!duplicated(pathway),]
      pathway <- ddply(pathway, .(Gene_ID_IndMode), summarise, pathway = list(as.character(Pathway)))
      
      complex <- cbind(data[48], data[11])
      complex <- complex[!duplicated(complex),]
      complex <- ddply(complex, .(Gene_ID_IndMode), summarise, complex = list(as.character(Complex.Type)))
      
      details <- cbind(data[48], data[12])
      details <- details[!duplicated(details),]
      details <- ddply(details, .(Gene_ID_IndMode), summarise, details = list(as.character(Details)))
      
      data <- data[!duplicated(data$Gene_ID_IndMode),]
      
      data <- dplyr::left_join(data, database, by="Gene_ID_IndMode")
      data <- dplyr::left_join(data, category, by="Gene_ID_IndMode")
      data <- dplyr::left_join(data, pathway, by="Gene_ID_IndMode")
      data <- dplyr::left_join(data, complex, by="Gene_ID_IndMode")
      data <- dplyr::left_join(data, details, by="Gene_ID_IndMode")
      data$database.x <- data$database.y;  names(data)[5] <- "database" 
      data$Category <- data$category
      data$Pathway <- data$pathway
      data$Complex.Type <- data$complex
      data$Details <- data$details
      
      
      data <- cbind(data[,1:18], data[,35:52])
      df <- as.data.frame(cbind(data[,26:27],data[24],data[,28:29],data[33],data[,14:16],data[,19:22],data[,4:5],data[8],data[,10:12],data[31],data[25],data[,34:36]))
      
      if(nrow(df[which(df$Feature.Subtype != "Protein"),]) > 0) {
        df_RNA <- df %>% filter(Feature.Subtype != "Protein") %>% mutate(SVM.score = "NA")
        df_Protein <- df %>% filter(Feature.Subtype == "Protein") 
        df <- rbind(df_Protein, df_RNA)
      }
      if(nrow(df[which(df$Conservation == "unknown"),]) > 0) {
        df_unknown <- df %>% filter(Conservation == "unknown") %>% mutate(Cyanobacteria = "NA")
        df_known <- df %>% filter(Conservation != "unknown") 
        df <- rbind(df_known, df_unknown) 
      }
      
      
      names(df) <- c("Locustag", "Genename", "Feature Type", "TU", "Annotation", "Location", "Cluster", "Reproducibility", "Log2(Abundance)", "Maximum Fraction", 
                     "Average Fraction", "No of Peaks", "Peak Fractions", "KEGG Orthology", "Database", "Category", "Pathway", "Complex/Type", "Details",
                     "Gene Ontology Nucleic acid binding", "AA Length", "Relative Conservation in Cyanobacteria", "Phylogenetic Group", "SVM Score")
      # Data.frame output:  
      df <- as_tibble(df)  
      
    }  
    else if(input$HeatmapChoice == "Both"){
      Gene_ID_in <- input$Gene_ID
      data_Ind <- KEGG
      data_Ind <- data_Ind %>% filter( Gene_ID_IndMode %in% Gene_ID_in ) 
      
      # merge the categories into one
      database <- cbind(data_Ind[48], data_Ind[5])
      database <- database[!duplicated(database),]
      database <- ddply(database, .(Gene_ID_IndMode), summarise, database = list(as.character(database)))
      
      category <- cbind(data_Ind[48], data_Ind[8])
      category <- category[!duplicated(category),]
      category <- ddply(category, .(Gene_ID_IndMode), summarise, category = list(as.character(Category)))
      
      pathway <- cbind(data_Ind[48], data_Ind[10])
      pathway <- pathway[!duplicated(pathway),]
      pathway <- ddply(pathway, .(Gene_ID_IndMode), summarise, pathway = list(as.character(Pathway)))
      
      complex <- cbind(data_Ind[48], data_Ind[11])
      complex <- complex[!duplicated(complex),]
      complex <- ddply(complex, .(Gene_ID_IndMode), summarise, complex = list(as.character(Complex.Type)))
      
      details <- cbind(data_Ind[48], data_Ind[12])
      details <- details[!duplicated(details),]
      details <- ddply(details, .(Gene_ID_IndMode), summarise, details = list(as.character(Details)))
      
      data_Ind <- data_Ind[!duplicated(data_Ind$Gene_ID_IndMode),]
      
      data_Ind <- dplyr::left_join(data_Ind, database, by="Gene_ID_IndMode")
      data_Ind <- dplyr::left_join(data_Ind, category, by="Gene_ID_IndMode")
      data_Ind <- dplyr::left_join(data_Ind, pathway, by="Gene_ID_IndMode")
      data_Ind <- dplyr::left_join(data_Ind, complex, by="Gene_ID_IndMode")
      data_Ind <- dplyr::left_join(data_Ind, details, by="Gene_ID_IndMode")
      data_Ind$database.x <- data_Ind$database.y;  names(data_Ind)[5] <- "database" 
      data_Ind$Category <- data_Ind$category
      data_Ind$Pathway <- data_Ind$pathway
      data_Ind$Complex.Type <- data_Ind$complex
      data_Ind$Details <- data_Ind$details
      
      
      data_Ind <- cbind(data_Ind[,1:18], data_Ind[,35:52])
      df_Ind <- as.data.frame(cbind(data_Ind[,26:27],data_Ind[24],data_Ind[,28:29],data_Ind[33],data_Ind[,14:16],data_Ind[,19:22],data_Ind[,4:5],data_Ind[8],data_Ind[,10:12],
                                    data_Ind[31],data_Ind[25],data_Ind[,34:36]))
      
      if(nrow(df_Ind[which(df_Ind$Feature.Subtype != "Protein"),]) > 0) {
        df_Ind_RNA <- df_Ind %>% filter(Feature.Subtype != "Protein") %>% mutate(SVM.score = "NA")
        df_Ind_Protein <- df_Ind %>% filter(Feature.Subtype == "Protein") 
        df_Ind <- rbind(df_Ind_Protein, df_Ind_RNA)
      }
      
      if(nrow(df_Ind[which(df_Ind$Conservation == "unknown"),]) > 0) {
        df_Ind_unknown <- df_Ind %>% filter(Conservation == "unknown") %>% mutate(Cyanobacteria = "NA")
        df_Ind_known <- df_Ind %>% filter(Conservation != "unknown") 
        df_Ind <- rbind(df_Ind_known, df_Ind_unknown) 
      }
      
      
      names(df_Ind) <- c("Locustag", "Genename", "Feature Type", "TU", "Annotation", "Location", "Cluster", "Reproducibility", "Log2(Abundance)", "Maximum Fraction", 
                     "Average Fraction", "No of Peaks", "Peak Fractions", "KEGG Orthology", "Database", "Category", "Pathway", "Complex/Type", "Details", 
                     "Gene Ontology Nucleic acid binding", "AA Length", "Relative Conservation in Cyanobacteria", "Phylogenetic Group", "SVM Score")
      # Data.frame output:  
      df_Ind <- as_tibble(df_Ind)  
      
      
      
      req(input$complex)
      
      Cluster_in <- input$Cluster
      Feature_in <- input$Feature
      Max_Peak <- input$Maximum
      No_Peak <- input$Peak_no
      Peak_in <- input$Peaks
      
      Location_in <- input$Location
      Phylogeny_in <- input$Phylogeny
      Cyanos_in_low <- input$Cyanobacteria[1]
      Cyanos_in_high <- input$Cyanobacteria[2]
      
      SVM_low <- input$SVM_score[1]
      SVM_high <- input$SVM_score[2]
      
      SVM_median_low <- input$Median_SVM_score[1]
      SVM_median_high <- input$Median_SVM_score[2]
      
      
      data <- complex()
      #data<- KEGG
      data[,19:34] <- as.data.frame( t( scale(t(data[, 19:34]), center = TRUE, scale = TRUE) ) )
      
      data <- transform(data, All_peaks = Peak.fractions)
      data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
      data <- data[!duplicated(data$Gene_ID_IndMode),]
      
      protein_data <- data %>% filter(Feature.Type == "Protein",
                                      Cluster %in% Cluster_in,
                                      Feature.Subtype %in% Feature_in,
                                      Maximum.Fraction %in% Max_Peak,
                                      No.of.peaks %in% No_Peak,
                                      location %in% Location_in,
                                      Conservation %in% Phylogeny_in,
                                      Cyanobacteria >= Cyanos_in_low,
                                      Cyanobacteria <= Cyanos_in_high,
                                      SVM.score >= SVM_low  | SVM.score == "NA",
                                      SVM.score <= SVM_high  | SVM.score == "NA",
                                      SVM.score.Median >= SVM_median_low,
                                      SVM.score.Median <= SVM_median_high) 
      
      RNA_data <- data %>% filter(Feature.Type != "Protein",
                                  Cluster %in% Cluster_in,
                                  Feature.Subtype %in% Feature_in,
                                  Maximum.Fraction %in% Max_Peak,
                                  No.of.peaks %in% No_Peak) 
      data <- rbind(protein_data, RNA_data)
      
      
      # merge the categories into one
      database <- cbind(data[48], data[5])
      database <- database[!duplicated(database),]
      database <- ddply(database, .(Gene_ID_IndMode), summarise, database = list(as.character(database)))
      
      category <- cbind(data[48], data[8])
      category <- category[!duplicated(category),]
      category <- ddply(category, .(Gene_ID_IndMode), summarise, category = list(as.character(Category)))
      
      pathway <- cbind(data[48], data[10])
      pathway <- pathway[!duplicated(pathway),]
      pathway <- ddply(pathway, .(Gene_ID_IndMode), summarise, pathway = list(as.character(Pathway)))
      
      complex <- cbind(data[48], data[11])
      complex <- complex[!duplicated(complex),]
      complex <- ddply(complex, .(Gene_ID_IndMode), summarise, complex = list(as.character(Complex.Type)))
      
      details <- cbind(data[48], data[12])
      details <- details[!duplicated(details),]
      details <- ddply(details, .(Gene_ID_IndMode), summarise, details = list(as.character(Details)))
      
      data <- data[!duplicated(data$Gene_ID_IndMode),]
      
      data <- dplyr::left_join(data, database, by="Gene_ID_IndMode")
      data <- dplyr::left_join(data, category, by="Gene_ID_IndMode")
      data <- dplyr::left_join(data, pathway, by="Gene_ID_IndMode")
      data <- dplyr::left_join(data, complex, by="Gene_ID_IndMode")
      data <- dplyr::left_join(data, details, by="Gene_ID_IndMode")
      data$database.x <- data$database.y;  names(data)[5] <- "database" 
      data$Category <- data$category
      data$Pathway <- data$pathway
      data$Complex.Type <- data$complex
      data$Details <- data$details
      
      
      data <- cbind(data[,1:18], data[,35:52])
      df <- as.data.frame(cbind(data[,26:27],data[24],data[,28:29],data[33],data[,14:16],data[,19:22],data[,4:5],data[8],data[,10:12],data[31],data[25],data[,34:36]))
      
       
      if(nrow(df[which(df$Feature.Subtype != "Protein"),]) > 0) {
        df_RNA <- df %>% filter(Feature.Subtype != "Protein") %>% mutate(SVM.score = "NA")
        df_Protein <- df %>% filter(Feature.Subtype == "Protein") 
        df <- rbind(df_Protein, df_RNA)
      }
      if(nrow(df[which(df$Conservation == "unknown"),]) > 0) {
        df_unknown <- df %>% filter(Conservation == "unknown") %>% mutate(Cyanobacteria = "NA")
        df_known <- df %>% filter(Conservation != "unknown") 
        df <- rbind(df_known, df_unknown) 
      }
      
      names(df) <- c("Locustag", "Genename", "Feature Type", "TU", "Annotation", "Location", "Cluster", "Reproducibility", "Log2(Abundance)", "Maximum Fraction", 
                     "Average Fraction", "No of Peaks", "Peak Fractions", "KEGG Orthology", "Database", "Category", "Pathway", "Complex/Type", "Details",
                     "Gene Ontology Nucleic acid binding", "AA Length", "Relative Conservation in Cyanobacteria", "Phylogenetic Group", "SVM Score")
      # Data.frame output:  
      df <- as_tibble(df)  
      df_out <- rbind(df_Ind, df)
     }  
  })
  output$GradSeqTable_summaryUI <- renderUI({
    DT::dataTableOutput("GradSeqTable_summary")
  })
  
 
  output$downloadPlot <- downloadHandler(  filename = function() { paste("GradSeq_Heatmap_",input$HeatmapChoice, '.pdf', sep='') },   content = function(file) {
      if(input$HeatmapChoice == "Locustags"){
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = GradSeqHeatmap_Ind(), device = device)
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(GradSeqHeatmap_Ind())
        dev.off()
       }  
      else if(input$HeatmapChoice == "Pathways"){
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(GradSeqHeatmap())
        dev.off()
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = GradSeqHeatmap(), device = device)
      }  
      else if(input$HeatmapChoice == "Both"){
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(plot_grid(GradSeqHeatmap_Ind(), GradSeqHeatmap(), ncol=1, align="v",  rel_heights = c(1,4) ))
        dev.off()
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = plot_grid(GradSeqHeatmap_Ind(), GradSeqHeatmap(), ncol=1, align="v",  rel_heights = c(1,4) ), device = device)
      }  
    } )
  output$downloadPlot2 <- downloadHandler( filename = function() { paste("Phylogeny_Heatmap_",input$HeatmapChoice, '.pdf', sep='') }, content = function(file) {
      if(input$HeatmapChoice == "Locustags"){
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(GradSeqPhylogeny_Ind())
        dev.off()
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = GradSeqPhylogeny_Ind(), device = device)
      }  
      else if(input$HeatmapChoice == "Pathways"){
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(GradSeqPhylogeny())
        dev.off()
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = GradSeqPhylogeny(), device = device)
      }  
      else if(input$HeatmapChoice == "Both"){
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(plot_grid(GradSeqPhylogeny_Ind(), GradSeqPhylogeny(), ncol=1, align="v",  rel_heights = c(1,4) ))
        dev.off()
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = plot_grid(GradSeqPhylogeny_Ind(), GradSeqPhylogeny(), ncol=1, align="v",  rel_heights = c(1,4) ), device = device)
      }  
    } )
  output$downloadPlot3 <- downloadHandler( filename = function() { paste("SVMscore_Boxplot_",input$HeatmapChoice, '.pdf', sep='') },  content = function(file) {
      if(input$HeatmapChoice == "Locustags"){
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(GradSeqSVM_Ind())
        dev.off()
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = GradSeqSVM_Ind(), device = device)
      }  
      else if(input$HeatmapChoice == "Pathways"){
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(GradSeqSVM())
        dev.off()
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = GradSeqSVM(), device = device)
      }  
      else if(input$HeatmapChoice == "Both"){
        pdf(file, width = HeatmapWidth(), height = HeatmapHeight()/50)
        print(plot_grid(GradSeqSVM_Ind(), GradSeqSVM(), ncol=1, align="v",  rel_heights = c(1,4) ) )
        dev.off()
        #device <- function(..., width, height) grDevices::png(..., width = HeatmapWidth(), height = HeatmapHeight()/50, res = 300, units = "cm")
        #ggsave(file, plot = plot_grid(GradSeqSVM_Ind(), GradSeqSVM(), ncol=1, align="v",  rel_heights = c(1,4) ), device = device)
      }  
    } )
  

  
   # Grad-Seq Statistics Tab outputs
  BarplotSize <- reactive({ req(input$heightBar)
    as.numeric(input$heightBar)
  })
  BarplotHeight <- reactive( (20 * BarplotSize()) )
  
  
  GradSeqBarplot <- reactive({
    
    Cluster_in <- input$BarCluster
    Feature_in <- input$BarFeature
    Max_Peak <- input$BarMaximum
    No_Peak <- input$BarPeak_no
    Peak_in <- input$BarPeaks
    
    Location_in <- input$BarLocation
    
    #Phylogeny_in <- input$BarPhylogeny
    #Cyanos_in_low <- input$BarCyanobacteria[1]
    #Cyanos_in_high <- input$BarCyanobacteria[2]
    
    #SVM_low <- input$BarSVM_score[1]
    #SVM_high <- input$BarSVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG
    data_dedup <- data[!duplicated(data$Cluster),]  %>% filter(Cluster >=1)
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                                    Feature.Subtype %in% Feature_in,
                                    Maximum.Fraction %in% Max_Peak,
                                    No.of.peaks %in% No_Peak,
                                    location %in% Location_in) #,
                                    #Conservation %in% Phylogeny_in,
                                    #Cyanobacteria >= Cyanos_in_low,
                                    #Cyanobacteria <= Cyanos_in_high,
                                    #SVM.score >= SVM_low  | SVM.score == "NA",
                                    #SVM.score <= SVM_high  | SVM.score == "NA") 
    

    data <- cbind(data[,1:18], data[,35:54])
    
    
    BarType <- switch(input$Bar,
                      "Absolute Counts" = "stack",
                      "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$BarOrder,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Cluster), fill = as.factor(Feature.Subtype) ) ) + 
      #geom_bar(color = "black", position = "fill", alpha=0.7 ) +
      geom_bar(color= "black", position = BarType, alpha=0.7) +
      scale_fill_manual("Feature Type",values = c("3UTR" = "pink3", "5UTR" = "pink4", "aRNA" ="blue", "crRNA" = "gold2", 
                                                  "mRNA"="pink", "sRNA" = "red", "Protein" = "black", "rRNA" ="purple", 
                                                  "Transposase" = "green2", "tRNA" = "turquoise1", "transcript" ="orange")) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Cluster", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"),
                       limits=data_dedup$Cluster[ order(as.numeric(data_dedup$Cluster))]) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
        legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
        legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans == "Log10" & input$Bar == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans == "Normal" & input$Bar == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    
})
  GradSeqBarplot_Max <- reactive({
    
    Cluster_in <- input$BarCluster
    Feature_in <- input$BarFeature
    Max_Peak <- input$BarMaximum
    No_Peak <- input$BarPeak_no
    Peak_in <- input$BarPeaks
    
    Location_in <- input$BarLocation
    
    #Phylogeny_in <- input$BarPhylogeny
    #Cyanos_in_low <- input$BarCyanobacteria[1]
    #Cyanos_in_high <- input$BarCyanobacteria[2]
    
    #SVM_low <- input$BarSVM_score[1]
    #SVM_high <- input$BarSVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG  %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    data_dedup <- data[!duplicated(data$Maximum.Fraction),] 
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype %in% Feature_in,
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in) #,
    #Conservation %in% Phylogeny_in,
    #Cyanobacteria >= Cyanos_in_low,
    #Cyanobacteria <= Cyanos_in_high,
    #SVM.score >= SVM_low  | SVM.score == "NA",
    #SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    BarType <- switch(input$Bar,
                      "Absolute Counts" = "stack",
                      "Relative Counts" = "fill"  )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Maximum.Fraction), fill = as.factor(Feature.Subtype) ) ) + 
      #geom_bar(color = "black", position = "fill", alpha=0.7 ) +
      geom_bar(color= "black", position = BarType, alpha=0.7) +
      scale_fill_manual("Feature Type",values = c("3UTR" = "pink3", "5UTR" = "pink4", "aRNA" ="blue", "crRNA" = "gold2", 
                                                  "mRNA"="pink", "sRNA" = "red", "Protein" = "black", "rRNA" ="purple", 
                                                  "Transposase" = "green2", "tRNA" = "turquoise1", "transcript" ="orange")) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Maximum Fraction", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$Maximum.Fraction[ order(as.numeric(data_dedup$Maximum.Fraction))]) ) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans == "Log10" & input$Bar == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans == "Normal" & input$Bar == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    
  })
  GradSeqBarplot_Peak <- reactive({
    
    Cluster_in <- input$BarCluster
    Feature_in <- input$BarFeature
    Max_Peak <- input$BarMaximum
    No_Peak <- input$BarPeak_no
    Peak_in <- input$BarPeaks
    
    Location_in <- input$BarLocation
    
    #Phylogeny_in <- input$BarPhylogeny
    #Cyanos_in_low <- input$BarCyanobacteria[1]
    #Cyanos_in_high <- input$BarCyanobacteria[2]
    
    #SVM_low <- input$BarSVM_score[1]
    #SVM_high <- input$BarSVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>% filter(Cluster >=1)
    data_dedup <- data[!duplicated(data$All_peaks),]  
    data <- data  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data[c(48,54)]),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype %in% Feature_in,
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in) #,
    #Conservation %in% Phylogeny_in,
    #Cyanobacteria >= Cyanos_in_low,
    #Cyanobacteria <= Cyanos_in_high,
    #SVM.score >= SVM_low  | SVM.score == "NA",
    #SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    BarType <- switch(input$Bar,
                      "Absolute Counts" = "stack",
                      "Relative Counts" = "fill"  )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(All_peaks), fill = as.factor(Feature.Subtype) ) ) + 
      #geom_bar(color = "black", position = "fill", alpha=0.7 ) +
      geom_bar(color= "black", position = BarType, alpha=0.7) +
      scale_fill_manual("Feature Type",values = c("3UTR" = "pink3", "5UTR" = "pink4", "aRNA" ="blue", "crRNA" = "gold2", 
                                                  "mRNA"="pink", "sRNA" = "red", "Protein" = "black", "rRNA" ="purple", 
                                                  "Transposase" = "green2", "tRNA" = "turquoise1", "transcript" ="orange")) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "All Peak Fractions", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$All_peaks[ order(as.numeric(data_dedup$All_peaks))]) ) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="bottom", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans == "Log10" & input$Bar == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans == "Normal" & input$Bar == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    
  })
  output$GradSeqBarplot_summary <- renderPlot({
    plot_grid( GradSeqBarplot(), GradSeqBarplot_Max(), GradSeqBarplot_Peak(), ncol=1, align="v",  rel_heights = c(1,1,1.2) )
   
  })
  output$GradSeqBarplot_summaryUI <- renderUI({
    plotOutput("GradSeqBarplot_summary", height = BarplotHeight())
  })
 
 
  output$downloadBarPlot <- downloadHandler( filename = function() { "Gradient_Composition_Barplot.png" }, content = function(file) {  device <- function(..., width, height) grDevices::png(..., width = BarplotWidth(), height = BarplotHeight()/50, res = 300, units = "cm")
        ggsave(file, plot = plot_grid( GradSeqBarplot(), GradSeqBarplot_Max(), GradSeqBarplot_Peak(), ncol=1, align="v",  rel_heights = c(1,1,1.5) ), device = device)
     } )

  
  
    database2 <- reactive({
    req(input$database2)
    filter(KEGG, database %in% input$database2)  
  })
    observeEvent(database2(), { 
      updatePickerInput(session, "category2", choices = unique(as.character(database2()$Category))) 
    })
    category2 <- reactive({    req(input$category2)
      filter(database2(), Category %in% input$category2)
    })
    observeEvent(category2(), { 
      updatePickerInput(session, "pathway2", choices = unique(as.character(category2()$Pathway)))
    })
    pathway2 <- reactive({    req(input$pathway2)
      filter(category2(), Pathway %in% input$pathway2)
    })
    observeEvent(pathway2(), {
      updatePickerInput(session, "complex2", choices = unique(as.character(pathway2()$Complex.Type)))
    })
    complex2 <- reactive({    req(input$complex2)
      filter(pathway2(), Complex.Type %in% input$complex2)
    })
  
  
  Bar2plotSize <- reactive({ req(input$heightBar2)
    as.numeric(input$heightBar2)
  })
  Bar2plotHeight <- reactive( (20 * Bar2plotSize()) )
  
  GradSeqBar2plot <- reactive({
    
    req(input$database2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- database2()
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    data_dedup <- KEGG[!duplicated(as.numeric(KEGG$Cluster)),]  %>% filter(Cluster >=1)
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                      "Absolute Counts" = "stack",
                      "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Cluster), fill = as.factor(database)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Cluster", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"),
                       limits=as.character(data_dedup$Cluster[ order(as.numeric(data_dedup$Cluster))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
  
    
  })
  GradSeqBar2plot_Max <- reactive({
    
    req(input$database2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- database2() %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    KEGG <- KEGG %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    data_dedup <- KEGG[!duplicated(KEGG$Maximum.Fraction),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Maximum.Fraction), fill = as.factor(database)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Maximum Fraction", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$Maximum.Fraction[ order(as.numeric(data_dedup$Maximum.Fraction))]) ) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Peak <- reactive({
    
    req(input$database2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- database2()
    
    data_dedup <- transform(KEGG, All_peaks = Peak.fractions)
    data_dedup <- data_dedup  %>%  separate_rows(All_peaks) %>% filter(Cluster >=1, All_peaks != "NA") 
    data_dedup <- data_dedup[!duplicated(data_dedup$All_peaks),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data[c(48,54)]),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(All_peaks), fill = as.factor(database)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "All Peak Fractions", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$All_peaks[ order(as.numeric(data_dedup$All_peaks))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="bottom", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  
  
  GradSeqBar2plot_cat <- reactive({
    
    req(input$category2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- category2()
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    data_dedup <- KEGG[!duplicated(as.numeric(KEGG$Cluster)),]  %>% filter(Cluster >=1)
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Cluster), fill = as.factor(Category)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Cluster", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"),
                       limits=as.character(data_dedup$Cluster[ order(as.numeric(data_dedup$Cluster))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Max_cat <- reactive({
    
    req(input$category2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- category2() %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    KEGG <- KEGG %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    data_dedup <- KEGG[!duplicated(KEGG$Maximum.Fraction),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Maximum.Fraction), fill = as.factor(Category)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Maximum Fraction", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$Maximum.Fraction[ order(as.numeric(data_dedup$Maximum.Fraction))]) ) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Peak_cat <- reactive({
    
    req(input$category2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- category2()
    
    data_dedup <- transform(KEGG, All_peaks = Peak.fractions)
    data_dedup <- data_dedup  %>%  separate_rows(All_peaks) %>% filter(Cluster >=1, All_peaks != "NA") 
    data_dedup <- data_dedup[!duplicated(data_dedup$All_peaks),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data[c(48,54)]),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(All_peaks), fill = as.factor(Category)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "All Peak Fractions", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$All_peaks[ order(as.numeric(data_dedup$All_peaks))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="bottom", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
 
  GradSeqBar2plot_path <- reactive({
    
    req(input$pathway2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- pathway2()
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    data_dedup <- KEGG[!duplicated(as.numeric(KEGG$Cluster)),]  %>% filter(Cluster >=1)
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Cluster), fill = as.factor(Pathway)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Cluster", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"),
                       limits=as.character(data_dedup$Cluster[ order(as.numeric(data_dedup$Cluster))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Max_path <- reactive({
    
    req(input$pathway2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- pathway2() %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    KEGG <- KEGG %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    data_dedup <- KEGG[!duplicated(KEGG$Maximum.Fraction),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Maximum.Fraction), fill = as.factor(Pathway)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Maximum Fraction", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$Maximum.Fraction[ order(as.numeric(data_dedup$Maximum.Fraction))]) ) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Peak_path <- reactive({
    
    req(input$pathway2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- pathway2()
    
    data_dedup <- transform(KEGG, All_peaks = Peak.fractions)
    data_dedup <- data_dedup  %>%  separate_rows(All_peaks) %>% filter(Cluster >=1, All_peaks != "NA") 
    data_dedup <- data_dedup[!duplicated(data_dedup$All_peaks),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data[c(48,54)]),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(All_peaks), fill = as.factor(Pathway)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "All Peak Fractions", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$All_peaks[ order(as.numeric(data_dedup$All_peaks))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="bottom", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
 
  GradSeqBar2plot_compl <- reactive({
    
    req(input$complex2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- complex2()
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    data_dedup <- KEGG[!duplicated(as.numeric(KEGG$Cluster)),]  %>% filter(Cluster >=1)
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Cluster), fill = as.factor(Complex.Type)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Cluster", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"),
                       limits=as.character(data_dedup$Cluster[ order(as.numeric(data_dedup$Cluster))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Max_compl <- reactive({
    
    req(input$complex2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- complex2() %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    KEGG <- KEGG %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    data_dedup <- KEGG[!duplicated(KEGG$Maximum.Fraction),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Maximum.Fraction), fill = as.factor(Complex.Type)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Maximum Fraction", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$Maximum.Fraction[ order(as.numeric(data_dedup$Maximum.Fraction))]) ) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Peak_compl <- reactive({
    
    req(input$complex2)
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- complex2()
    
    data_dedup <- transform(KEGG, All_peaks = Peak.fractions)
    data_dedup <- data_dedup  %>%  separate_rows(All_peaks) %>% filter(Cluster >=1, All_peaks != "NA") 
    data_dedup <- data_dedup[!duplicated(data_dedup$All_peaks),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data[c(48,54)]),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(All_peaks), fill = as.factor(Complex.Type)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "All Peak Fractions", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$All_peaks[ order(as.numeric(data_dedup$All_peaks))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="bottom", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  
  GradSeqBar2plot_cons <- reactive({
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    data_dedup <- KEGG[!duplicated(as.numeric(KEGG$Cluster)),]  %>% filter(Cluster >=1)
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Cluster), fill = as.factor(Conservation)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Cluster", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"),
                       limits=as.character(data_dedup$Cluster[ order(as.numeric(data_dedup$Cluster))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Max_cons <- reactive({
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    KEGG <- KEGG %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    data_dedup <- KEGG[!duplicated(KEGG$Maximum.Fraction),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Maximum.Fraction), fill = as.factor(Conservation)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Maximum Fraction", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$Maximum.Fraction[ order(as.numeric(data_dedup$Maximum.Fraction))]) ) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Peak_cons <- reactive({
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG
    
    data_dedup <- transform(KEGG, All_peaks = Peak.fractions)
    data_dedup <- data_dedup  %>%  separate_rows(All_peaks) %>% filter(Cluster >=1, All_peaks != "NA") 
    data_dedup <- data_dedup[!duplicated(data_dedup$All_peaks),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data[c(48,54)]),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(All_peaks), fill = as.factor(Conservation)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "All Peak Fractions", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$All_peaks[ order(as.numeric(data_dedup$All_peaks))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="bottom", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  
  GradSeqBar2plot_loc <- reactive({
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    data_dedup <- KEGG[!duplicated(as.numeric(KEGG$Cluster)),]  %>% filter(Cluster >=1)
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Cluster), fill = as.factor(location)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Cluster", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"),
                       limits=as.character(data_dedup$Cluster[ order(as.numeric(data_dedup$Cluster))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Max_loc <- reactive({
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    #data <- droplevels(subset(KEGG, database == "syn00001 KEGG Orthology"))
    KEGG <- KEGG %>% filter(Cluster >=1, Maximum.Fraction != "NA")
    data_dedup <- KEGG[!duplicated(KEGG$Maximum.Fraction),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data$Gene_ID_IndMode),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(Maximum.Fraction), fill = as.factor(location)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "Maximum Fraction", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$Maximum.Fraction[ order(as.numeric(data_dedup$Maximum.Fraction))]) ) +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="none", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  GradSeqBar2plot_Peak_loc <- reactive({
    
    Cluster_in <- input$Bar2Cluster
    Max_Peak <- input$Bar2Maximum
    No_Peak <- input$Bar2Peak_no
    Peak_in <- input$Bar2Peaks
    
    Location_in <- input$Bar2Location
    
    Phylogeny_in <- input$Bar2Phylogeny
    Cyanos_in_low <- input$Bar2Cyanobacteria[1]
    Cyanos_in_high <- input$Bar2Cyanobacteria[2]
    
    SVM_low <- input$Bar2SVM_score[1]
    SVM_high <- input$Bar2SVM_score[2]
    
    
    
    
    # Grad-Seq Heatmap Data.frame:
    data<- KEGG
    
    data_dedup <- transform(KEGG, All_peaks = Peak.fractions)
    data_dedup <- data_dedup  %>%  separate_rows(All_peaks) %>% filter(Cluster >=1, All_peaks != "NA") 
    data_dedup <- data_dedup[!duplicated(data_dedup$All_peaks),]  
    
    
    data <- transform(data, All_peaks = Peak.fractions)
    data <- data  %>%  separate_rows(All_peaks)  %>%  filter(All_peaks %in% Peak_in)
    data <- data[!duplicated(data[c(48,54)]),]
    
    
    data <- data %>% filter(Cluster %in% Cluster_in,
                            Feature.Subtype == "Protein",
                            Maximum.Fraction %in% Max_Peak,
                            No.of.peaks %in% No_Peak,
                            location %in% Location_in,
                            Conservation %in% Phylogeny_in,
                            Cyanobacteria >= Cyanos_in_low,
                            Cyanobacteria <= Cyanos_in_high,
                            SVM.score >= SVM_low  | SVM.score == "NA",
                            SVM.score <= SVM_high  | SVM.score == "NA") 
    
    
    data <- cbind(data[,1:18], data[,35:54])
    
    
    Bar2Type <- switch(input$Bar2,
                       "Absolute Counts" = "stack",
                       "Relative Counts" = "fill"  )
    #Cluster_MaxPeak <- switch(input$Bar2Order,
    #                          "Cluster" = as.factor(data$Cluster), 
    #                          "Maximum Fraction"  = as.factor(data$Maximum.Fraction) )
    
    
    
    barplot <- ggplot(data, aes(x = as.character(All_peaks), fill = as.factor(location)) ) + 
      #geom_bar(color = "black", position = "stack") +
      #geom_bar(color = "black", position = "dodge", stat = "identity") +
      geom_bar(color= "black", position = Bar2Type, alpha=0.7) +
      scale_y_continuous(name = "Counts (absolute / relative)") +
      scale_x_discrete(name = "All Peak Fractions", breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14.5", "16.5", "18"),
                       limits=as.character( data_dedup$All_peaks[ order(as.numeric(data_dedup$All_peaks))]) )  +
      theme_bw() +
      theme(panel.grid=element_blank(), plot.title=element_text(size=14, face="bold"), legend.title = element_blank(),
            legend.position="bottom", panel.spacing.y = unit(-0.5, "lines"),
            legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
    
    if(input$logtrans2 == "Log10" & input$Bar2 == "Absolute Counts" )  {
      barplot + scale_y_log10() + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else if (input$logtrans2 == "Normal" & input$Bar2 == "Absolute Counts" ) {
      barplot + geom_hline(yintercept = 10, linetype = "dashed", color="grey50")+
        geom_hline(yintercept = 0, linetype = "solid")
    } else {
      barplot
    }
    
    #ggplotly(plotly_build(barplot))
    
    
  })
  
  
  output$GradSeqBar2plot_summary <- renderPlot({
    plot_grid( GradSeqBar2plot(), GradSeqBar2plot_Max(), GradSeqBar2plot_Peak(), ncol=1, align="v",  rel_heights = c(1,1,1.2) )
  })
  output$GradSeqBar2plot_summaryUI <- renderUI({
    plotOutput("GradSeqBar2plot_summary", height = Bar2plotHeight())
  })
  
  output$GradSeqBar2plot_cat_summary <- renderPlot({
    plot_grid( GradSeqBar2plot_cat(), GradSeqBar2plot_Max_cat(), GradSeqBar2plot_Peak_cat(), ncol=1, align="v",  rel_heights = c(1,1,1.2) )
  })
  output$GradSeqBar2plot_cat_summaryUI <- renderUI({
    plotOutput("GradSeqBar2plot_cat_summary", height = Bar2plotHeight())
  })
 
  output$GradSeqBar2plot_path_summary <- renderPlot({
    plot_grid( GradSeqBar2plot_path(), GradSeqBar2plot_Max_path(), GradSeqBar2plot_Peak_path(), ncol=1, align="v",  rel_heights = c(1,1,1.2) )
  })
  output$GradSeqBar2plot_path_summaryUI <- renderUI({
    plotOutput("GradSeqBar2plot_path_summary", height = Bar2plotHeight())
  })
  
  output$GradSeqBar2plot_compl_summary <- renderPlot({
    plot_grid( GradSeqBar2plot_compl(), GradSeqBar2plot_Max_compl(), GradSeqBar2plot_Peak_compl(), ncol=1, align="v",  rel_heights = c(1,1,1.2) )
  })
  output$GradSeqBar2plot_compl_summaryUI <- renderUI({
    plotOutput("GradSeqBar2plot_compl_summary", height = Bar2plotHeight())
  })
  
  output$GradSeqBar2plot_cons_summary <- renderPlot({
    plot_grid( GradSeqBar2plot_cons(), GradSeqBar2plot_Max_cons(), GradSeqBar2plot_Peak_cons(), ncol=1, align="v",  rel_heights = c(1,1,1.2) )
  })
  output$GradSeqBar2plot_cons_summaryUI <- renderUI({
    plotOutput("GradSeqBar2plot_cons_summary", height = Bar2plotHeight())
  })
  
  output$GradSeqBar2plot_loc_summary <- renderPlot({
    plot_grid( GradSeqBar2plot_loc(), GradSeqBar2plot_Max_loc(), GradSeqBar2plot_Peak_loc(), ncol=1, align="v",  rel_heights = c(1,1,1.2) )
  })
  output$GradSeqBar2plot_loc_summaryUI <- renderUI({
    plotOutput("GradSeqBar2plot_loc_summary", height = Bar2plotHeight())
  })
  
  
  output$downloadBar2Plot <- downloadHandler( filename = function() { "Protein_Composition_database_Barplot.png" }, content = function(file) {  device <- function(..., width, height) grDevices::png(..., width = Bar2plotWidth(), height = Bar2plotHeight()/50, res = 300, units = "cm")
  ggsave(file, plot =  plot_grid( GradSeqBar2plot(), GradSeqBar2plot_Max(), GradSeqBar2plot_Peak(), ncol=1, align="v",  rel_heights = c(1,1,1.5) ), device = device)
  } )
  output$downloadBar2_cat_Plot <- downloadHandler( filename = function() { "Protein_Composition_category_Barplot.png" }, content = function(file) {  device <- function(..., width, height) grDevices::png(..., width = Bar2plotWidth(), height = Bar2plotHeight()/50, res = 300, units = "cm")
  ggsave(file, plot =  plot_grid( GradSeqBar2plot_cat(), GradSeqBar2plot_Max_cat(), GradSeqBar2plot_Peak_cat(), ncol=1, align="v",  rel_heights = c(1,1,1.5) ), device = device)
  } )
  output$downloadBar2_path_Plot <- downloadHandler( filename = function() { "Protein_Composition_pathway_Barplot.png" }, content = function(file) {  device <- function(..., width, height) grDevices::png(..., width = Bar2plotWidth(), height = Bar2plotHeight()/50, res = 300, units = "cm")
  ggsave(file, plot =  plot_grid( GradSeqBar2plot_path(), GradSeqBar2plot_Max_path(), GradSeqBar2plot_Peak_path(), ncol=1, align="v",  rel_heights = c(1,1,1.5) ), device = device)
  } )
  output$downloadBar2_compl_Plot <- downloadHandler( filename = function() { "Protein_Composition_complex_Barplot.png" }, content = function(file) {  device <- function(..., width, height) grDevices::png(..., width = Bar2plotWidth(), height = Bar2plotHeight()/50, res = 300, units = "cm")
  ggsave(file, plot =  plot_grid( GradSeqBar2plot_compl(), GradSeqBar2plot_Max_compl(), GradSeqBar2plot_Peak_compl(), ncol=1, align="v",  rel_heights = c(1,1,1.5) ), device = device)
  } )
  output$downloadBar2_cons_Plot <- downloadHandler( filename = function() { "Protein_Composition_phylogeny_Barplot.png" }, content = function(file) {  device <- function(..., width, height) grDevices::png(..., width = Bar2plotWidth(), height = Bar2plotHeight()/50, res = 300, units = "cm")
  ggsave(file, plot =  plot_grid( GradSeqBar2plot_cons(), GradSeqBar2plot_Max_cons(), GradSeqBar2plot_Peak_cons(), ncol=1, align="v",  rel_heights = c(1,1,1.5) ), device = device)
  } )
  output$downloadBar2_loc_Plot <- downloadHandler( filename = function() { "Protein_Composition_location_Barplot.png" }, content = function(file) {  device <- function(..., width, height) grDevices::png(..., width = Bar2plotWidth(), height = Bar2plotHeight()/50, res = 300, units = "cm")
  ggsave(file, plot =  plot_grid( GradSeqBar2plot_loc(), GradSeqBar2plot_Max_loc(), GradSeqBar2plot_Peak_loc(), ncol=1, align="v",  rel_heights = c(1,1,1.5) ), device = device)
  } )
  
  output$GradSeqTable_mbgd_genomes <- DT::renderDataTable({ 
    df <- mbgd_genomes
    names(df) <- c("order in shinyApp", "abbrev. in shinyApp", "organism", "abbrev. in MBGD")
    df <- as_tibble(df)
 
  })
  output$GradSeqTable_mbgd_genomesUI <- renderUI({
    DT::dataTableOutput("GradSeqTable_mbgd_genomes")
  })
  

}


shinyApp(ui = ui, server=server)  

