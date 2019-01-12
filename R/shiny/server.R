#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(RISmed)
library(plotly)
library(qdap)
library(wordcloud2)
library(shinycssloaders)

function(input, output) {
  glob_values <- reactiveValues(
    search_query=NULL,
    records=NULL,
    pubmed_data=NULL,
    count=NULL,
    ord=NULL,
    auths_table=NULL,
    title=NULL,
    journal_table=NULL)
  reactiveValues.reset <-function(){
    search_query=NULL
    records=NULL
    pubmed_data=NULL
    count=NULL
    ord=NULL
    auths_table=NULL
    title=NULL
    journal_table=NULL
  }
  observeEvent(input$do,{
    reactiveValues.reset()
    glob_values$search_query <- EUtilsSummary(c(input$keyword),
                                              db="pubmed", 
                                              retmax=500,
                                              datetype='pdat',
                                              mindate=input$from, 
                                              maxdate=input$to)
    
    glob_values$records<- EUtilsGet(glob_values$search_query)
    glob_values$pubmed_data <- data.frame('Title'=ArticleTitle(glob_values$records),
                 'Year'=YearAccepted(glob_values$records),
                 'Journal'=ISOAbbreviation(glob_values$records),
                 'Abstract'=AbstractText(glob_values$records))
    glob_values$count <- as.data.frame(table(YearPubmed(EUtilsGet(glob_values$search_query))))
    names(glob_values$count) <-c ("Year", "Counts")
    
    abstractsOnly<-as.character(glob_values$pubmed_data$Abstract)
    abstractsOnly<-paste(abstractsOnly, sep="", collapse="")
    abstractsOnly<-as.vector(abstractsOnly)
    abstractsOnly<-strip(abstractsOnly)
    stsp<-rm_stopwords(abstractsOnly, stopwords = qdapDictionaries::Top100Words)
    glob_values$ord<-as.data.frame(table(stsp))
    glob_values$ord<-glob_values$ord[order(glob_values$ord$Freq, decreasing=TRUE),]
    
    auths<-Author(glob_values$records)
    Last<-sapply(auths, function(x)paste(x$LastName))
    Fore<-sapply(auths, function(x)paste(x$ForeName))
    name <- paste(unlist(Fore), unlist(Last), sep = " ")
    auths2<-as.data.frame(sort(table(unlist(name)), dec=TRUE))
    auths2 <- as.data.frame(auths2)
    colnames(auths2) <- c("Author","Freq")
    glob_values$auths_table <- auths2
    
    TitleOnly<-as.character(glob_values$pubmed_data$Title)
    TitleOnly<-paste(TitleOnly, sep="", collapse="")
    TitleOnly<-as.vector(TitleOnly)
    TitleOnly<-strip(TitleOnly)
    sts<-rm_stopwords(TitleOnly, stopwords = qdapDictionaries::Top100Words)
    glob_values$title<-as.data.frame(table(sts))
    glob_values$title<-glob_values$title[order(glob_values$title$Freq, decreasing=TRUE),]
    
    journalsOnly <- glob_values$pubmed_data$Journal
    glob_values$journal_table <- as.data.frame(sort(table(unlist(journalsOnly)), dec=TRUE))
  })

  output$summary <- renderText({
    recount <- attr(glob_values$search_query,"count")
    paste("We found", recount, "records.")
  })

  output$barPlot <- renderPlotly(
    if (! is.null(glob_values$count)){
      plot_ly(glob_values$count, x =~Year, y =~Counts, fillcolor=~Year) 
    }
  )

  output$my_wc<-renderWordcloud2(
    if (! is.null(glob_values$ord)){
      wordcloud2(glob_values$ord, size = 2, minRotation = -pi/6, maxRotation = -pi/6,
                 rotateRatio = 1)})
  
  output$author_table <- DT::renderDataTable({if (! is.null(glob_values$auths_table))glob_values$auths_table},
                                             rownames= FALSE,options = list(
    pageLength = 10,
    lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
    scrollX = TRUE,
    fixedHeader = TRUE,
    fixedColumns = TRUE ,
    deferRender = TRUE
  ),
  escape = FALSE
  )
  
  output$my_wc2<-renderWordcloud2(if (! is.null(glob_values$title)){
    wordcloud2(glob_values$title, size = 2, minRotation = -pi/6, maxRotation = -pi/6,
                                            rotateRatio = 1)})
  
  output$journal_table <- DT::renderDataTable({if (!is.null(glob_values$journal_table))glob_values$journal_table},
                                             rownames= FALSE,options = list(
                                               pageLength = 10,
                                               lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                                               scrollX = TRUE,
                                               fixedHeader = TRUE,
                                               fixedColumns = TRUE ,
                                               deferRender = TRUE
                                             ),
                                             escape = FALSE
  )
  
  output$downloadData <- downloadHandler(
    filename = function(){paste(input$keyword, '.csv', sep='')},
    content = function(file) {
      write.csv(glob_values$pubmed_data, file)
    }
  )
  
  output$result_table <- DT::renderDataTable({if (!is.null(glob_values$pubmed_data))glob_values$pubmed_data},
                                              rownames= FALSE,options = list(
                                                pageLength = 10,
                                                lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                                                scrollX = TRUE,
                                                fixedHeader = TRUE,
                                                fixedColumns = TRUE ,
                                                deferRender = TRUE
                                              ),
                                              escape = FALSE
  )
}