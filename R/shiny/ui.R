#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
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
library(dashboardthemes)

logo <- shinyDashboardLogoDIY(
  
  boldText = "PubMed"
  ,mainText = "Search"
  ,textSize = 16
  ,badgeText = "BETA"
  ,badgeTextColor = "white"
  ,badgeTextSize = 2
  ,badgeBackColor = "#40E0D0"
  ,badgeBorderRadius = 3
  
)

dashboardPage(
  dashboardHeader(title = logo),
  dashboardSidebar(width=250,
                   sidebarMenu(
                     menuItem("Home Page",icon=icon("home"),tabName="home_page"),
                     menuItem("Keyword Search",icon=icon("search"),tabName="search")
                     )
  ),
  dashboardBody(
    shinyDashboardThemes(
      theme = "poor_mans_flatly"
    ),
    tabItems(
      # First tab content
      tabItem(tabName = "home_page",
              h4("Search the PubMed database at NCBI using the RISmed package in R."),
              h4("Choose from the menu on the left."),
              h4("** We have two brand new upgrades to our Keyword Search! Now you can see the top ten journal titles and ten most recent article titles for your keyword. See the 'JOURNALS' and 'ARTICLES' tabs."),
              h4("You can still compare the most frequent words for two different timespans to see how the content of the abstracts have changed. Use this feature to track the changing interests of an author, or investigate what led to a spike in interest in a topic. See the 'COMPARISON' tab."),
              h4("Remember to be patient if you're performing a big search."),
              h4("'Scholar Indices' rates the impact of an author based on different metrics."),
              h4("More updates coming soon...")
      ),
      
      # Second tab content
      tabItem(tabName = "search",
              tabsetPanel(
                tabPanel("HOME",
                         h4("Below, type a keyword to search PubMed and find documents that contain that word in the text."),
                         h4("You can even type multiple words. You can search authors, topics, any acronym, etc."),
                         h4("Specify the start and end dates you'd like to search, using the format YYYY."),
                         h4("Then click 'RUN' and scroll through the tabs to see the results."),
                         textInput('keyword','Keyword(s)', value='lncRNA'),
                         textInput('from','From', value='2013'),
                         textInput('to','To', value='2018'),
                         actionButton("do", "RUN", icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                         helpText("If your keyword is something like 'cancer' or 'research', you might bump up against the search limit. Try being more specific or searching fewer years."),
                         textOutput("summary")  %>% withSpinner(),
                         DT::dataTableOutput('result_table') %>% withSpinner()
                         ),
                tabPanel("BARPLOT",
                         h4("A barplot of publications containing your keyword per year."),
                         plotlyOutput("barPlot")  %>% withSpinner()
                         ),
                tabPanel("WORDCLOUD",
                         h4("A wordcloud of the abstracts for your keyword."),
                         wordcloud2Output("my_wc", width = "100%", height = "400px") %>% withSpinner()
                         ),
                tabPanel("AUTHORS",
                         h4("A table of authors with publicatons containing your keyword."),
                         DT::dataTableOutput('author_table') %>% withSpinner()),
                tabPanel("JOURNAL TITLES",
                         h4("A table of journal titles of publicatons containing your keyword."),
                         DT::dataTableOutput('journal_table') %>% withSpinner()),
                tabPanel("ARTICLE TITLES",
                         h4("A table of recent article titles of publicatons containing your keyword."),
                         wordcloud2Output("my_wc2", width = "100%", height = "400px")%>% withSpinner()),
                tabPanel("DOWNLOAD",
                         h4("Download the csv file."),
                         downloadButton('downloadData', 'Download'))
              )
      )
    )
  )
)
