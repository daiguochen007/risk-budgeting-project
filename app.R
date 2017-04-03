library(shiny)
library(pheatmap)

ret_all_data<- ret_xts(all_data)

ui<- fluidPage(
  tags$h1("Dynamic covariance matrix plot"),
  tags$h6("By Guochen Dai"),
  wellPanel(
    sliderInput(inputId = "sampledays",label = "choose sample days",min = 30, max= 500,step = 10,value = 90),
    sliderInput(inputId = "startdate",label = "choose sample start date",min = index(ret_all_data)[1], max = index(ret_all_data)[nrow(ret_all_data)-500],value=index(ret_all_data)[1])
  ),
  plotOutput(outputId = "htmap")
  
)


server<- function(input,output){
  output$htmap<- renderPlot({
    startnum<- min(which(index(ret_all_data)>=input$startdate))
    title<- paste("Cov matrix with data from",index(ret_all_data)[startnum],"to",index(ret_all_data)[startnum+input$sampledays])
    pheatmap(cov(ret_all_data[startnum:(startnum+input$sampledays),])*252,cluster_rows = F,cluster_cols = F, main=title)
    })
}

shinyApp(ui=ui,server = server)
