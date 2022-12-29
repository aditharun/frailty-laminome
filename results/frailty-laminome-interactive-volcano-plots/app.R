library(shiny)
library(tidyverse)
library(readxl)
library(plotly)

heart <- read_excel("shiny-data/heart-analysis.xlsx")
brain <- read_excel("shiny-data/brain-analysis.xlsx")
muscle <- read_excel("shiny-data/muscle-analysis.xlsx")

tissuelist <- c("Skeletal Muscle", "Heart", "Brain")

# Define UI for application that draws a histogram
ui <- fluidPage(
   sidebarLayout(
      sidebarPanel(
        selectInput("tissue", "Tissue", choices = tissuelist),width="auto"),
      mainPanel(
         plotlyOutput('distplot')
      )
  )
)


server <- function(input, output,session) {

  


   output$distplot <- renderPlotly({

      if (input$tissue=="Brain"){
        analysis <- brain
      }

      if (input$tissue=="Skeletal Muscle"){
        analysis <- muscle
      }

      if (input$tissue == "Heart"){
        analysis <- heart
      }


      color <- function(x){
        x <- x %>% mutate(expression=ifelse(`p-value` < 0.05 ,'p-value < 0.05','p-value \U2265 0.05'))
        colnames(x)[1] <- "fold.change"
        x
      }

    analysis <- color(analysis)

    plot_volcano <- function(x, sample, show_legend){

      textsize <- 15
      titlesize <- 18


      cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

      if (show_legend){
        leg <- "bottom"
      } else{
        leg <- "none"
      }

      xlims <- x$log2.foldchange %>% quantile(c(0.05, 0.95)) %>% unname() %>% round()

      ggplot(data = x,
           aes(x = log2.foldchange,
               y = -log10(`p-value`),
               colour=expression)) +
        geom_point(alpha=0.4, size=3.5) +
        scale_color_manual(values=c("orange", "blue"))+
        xlim(xlims) +
        geom_vline(xintercept=c(-.05,.05),lty=2,col="grey70",lwd=0.3) +
        geom_hline(yintercept = 1.301,lty=4,col="grey70",lwd=0.3) +
        xlab("log2 fold change (KO/WT)")+
             ylab("-log10 p-value")+
             ggtitle(sample) +
                 theme_bw()+
                 theme(plot.title = element_text(hjust = 0.5, size=titlesize),
                       legend.position=leg,
                       legend.title = element_blank()) +   theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank() # get rid of minor grid
      ) + theme(
      axis.title.x = element_text(size = textsize),
      axis.text.x = element_text(size = textsize-1),
      axis.title.y = element_text(size = textsize),
      axis.text.y= element_text(size=textsize-1),
      legend.text = element_text(size=textsize-1)
      )

  }

  xlims <- analysis$log2.foldchange %>% quantile(c(0.05, 0.95)) %>% unname()


  p.interactive <- plot_ly(
  analysis, x = ~log2.foldchange, y = ~-log10( `p-value`),
  # Hover text:
  text = ~paste0(toupper(gene_names), "\n", Protein),
  color = ~expression, colors=c("orange", "blue")) %>% layout(title=input$tissue, yaxis=list(title="-log10(p-value)"), xaxis=list(title="Log2 Fold Change (KO/WT)", range=xlims) )
      
   })

}

# Run the application 
shinyApp(ui = ui, server = server)