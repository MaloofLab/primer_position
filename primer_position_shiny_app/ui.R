#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Primer Position Viewer"),
  
  helpText(
    p("This app produces an image of your primers' locations relative to the genes they are designed for.
      You can use one of the installed organism databases or create your own. When creating your own database, any standard gff and nomenclature will work. 
      The Primer Information CSV requires 8 filled in columns, (primername, gene_id, Chr, start, end, start2, end2, strand), do not include prefixes in Chr value.
      In order to produce an image you must first load/make a database. After the database has been loaded/made, you will get a notfication message. 
      From there you can proceed to display primers. Primer set is the pair of forward and reverse primers in your csv file from top to bottom.  
      Below is an example of a valid CSV when read into a data frame for primer set 1."),
    img(src = "example.png")
  ),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      radioButtons("choice1",
                   "Use existing organism data or input own data?",
                   c("Use existing","Input own data")),
      
      conditionalPanel(
        condition = 'input.choice1 == "Input own data"',
        
        fileInput('file1', 'GFF File',
                accept=c('text/gff',
                         '.gff',
                         '.gff2',
                         '.gff3',
                         '.gff.gz',
                         '.gff2.gz',
                         '.gff3.gz')),
        
        textInput('organism', label = "Name of Organism", placeholder = "ie. Arabidopsis thaliana"),
        
        textInput('prefix', label = "Chromosome prefix in GFF file", placeholder = "Chr"),
        
        fileInput('file2', 'Primer Information CSV',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      
      
        numericInput('pair1', label = "Primer Set", value = 1, min = 1),
      
        actionButton("txdb1", "Make Database", style='padding:4px'),
      
        actionButton("image", "Display Primers", style='padding:4px')),
      
      conditionalPanel(
        condition = 'input.choice1 == "Use existing"',
        
        selectInput("existing",
                    label = "Choose a species ",
                    choices = list("Arabidopsis thaliana" = "A.tha",
                                   "Brassica napus" = "B.nap",
                                   "Homo sapiens" = "H.sap",
                                   "Solanum lycopersicum" = "S.lyc")),
        
        fileInput('file3', 'Primer Information CSV',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        
        numericInput('pair2', label = "Primer Set", value = 1, min = 1),
        
        actionButton("txdb2", "Load Database", style='padding:4px'),
        
        actionButton("image2", "Display Primers", style='padding:4px')
      )
    ),
    # Show a plot of the generated distribution
    mainPanel(
      textOutput('text'),
      plotOutput('plot')
    )
  )
))
