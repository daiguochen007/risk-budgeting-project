library(shiny)
library(pheatmap)
library(plot3D)

ret_all_data<- ret_xts(all_data)

ui<- fluidPage(
  tags$h1("Dynamic covariance matrix plot"),
  tags$h6("By Guochen Dai"),
  
  wellPanel(
    sliderInput(inputId = "daterange",label = "choose sample start date",
                 min = index(ret_all_data)[1], max = index(ret_all_data)[nrow(ret_all_data)],
                 value=c(index(ret_all_data)[1],index(ret_all_data)[252]),animate = T),
    textOutput("textdays")
  ),

  plotOutput("plt3d"),
  plotOutput(outputId = "htmap")
  
)

bks=seq(-0.05,0.4,0.025)
server<- function(input,output){
  startnum<- reactive({min(which(index(ret_all_data)>=input$daterange[1]))})
  endnum<- reactive({max(which(index(ret_all_data)<=input$daterange[2]))})
  
  output$textdays<- renderText({paste("count:",endnum()-startnum()+1,"days in sample")})
  
  output$htmap<- renderPlot({
    title<- paste("cov matrix (data:",index(ret_all_data)[startnum()],"to",index(ret_all_data)[endnum()],")")
    cov<- cov(ret_all_data[startnum():endnum(),])*252
    pheatmap(cov,cluster_rows = F,cluster_cols = F,main=title)
    })
  
  output$plt3d<- renderPlot({   
    title<- paste("3D cov matrix (data:",index(ret_all_data)[startnum()],"to",index(ret_all_data)[endnum()],")")
    cov<- cov(ret_all_data[startnum():endnum(),])*252
    persp3D(1:nrow(cov),1:ncol(cov),cov,zlim=c(-0.05,0.7),breaks=bks,main=title)
  })
}

shinyApp(ui=ui,server = server)
