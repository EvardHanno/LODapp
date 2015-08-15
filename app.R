# LoD teaching tool

## Lisa veel: CCa joon + koht kus Ã¼tleb, kui palju piigi stdevist on alla CCa
## Hpeak-i muutes ei arvutata uusi andmeid - kõrgus muutub, aga uusi andmeid ei anta.
## Defineeri krom tekitamise funktsioon eraldi ning siis kasuta seda reactive() kohtades


library(shiny)
library(ggplot2)

ui <- fluidPage(

  
  titlePanel("Title"),
  
    sidebarLayout(
      sidebarPanel(
        # Data point frequency
        sliderInput("points", label = "Number of datapoints on chromatogram per unit of time:", min = 1, max = 5, value = 2.5, step = 0.5),
      
        # mean peak height
        sliderInput("Hpeak", label = "Mean height of peak:", min = 0, max = 15, value = 2, step = 0.2),
        
        # noise stdev
        sliderInput("noise", label = "Standard devation of noise:", min = 0.1, max = 10, value = 1),
        
        # Homo- or heteroscedastic
        selectInput("selectSced", label = "Choose if system is homo- or heteroscedastic:", choices = c(Homoscedastic = "homosc", Heteroscedastic = "heterosc")),
      
        # If homoscedastic then show text explaning stdev of peak = stdev of noise
        conditionalPanel(
          condition = "input.selectSced == 'homosc'",
          textOutput("homoscText")
        ),
        
        # If heteroscedastic then show slider for stdev slope
        conditionalPanel(
          condition = "input.selectSced == 'heterosc'",
          sliderInput("peakHvar", label = "Slope/strength of heteroscedasticity:", min = 0, max = 5, value = 1, step = 0.25)
        ),
      width = 3),

      ## Main panel outputs
  
      mainPanel(
        fluidRow(
          column( 6, 
            plotOutput("krom")
          ),
          
          column( 6,
            plotOutput("graph")
          )
      )  
    )
  )
)
  
server <- function(input, output) {
  
  # Create reactive function of baseline
  baseline <- reactive({
    
    x = seq(-25, 25, 1/input$points)
    
  })
  
  # create reactive function of noise
  noise <- reactive({
  
    noise = rnorm(length(baseline()), sd = input$noise)
    
  })
  
  # Create reactive function of chromatogam
  krom <- reactive({
   
    stdev = 1
    peak = input$Hpeak/0.3183099 * (1/(stdev*(2*pi)^1/2))*exp(-((baseline()-0)^2/(2*stdev^2)))
    
    krom = noise() + peak
    
  })
  
  sced <- reactive({
    if(input$selectSced == "homosc"){
      sced <- input$noise
    }
    else{
      sced <- input$noise + input$peakHvar * input$Hpeak/10
    }
  })
  
  
  # Render of chromatogram
  output$krom <- renderPlot({
    
    qplot(baseline(), krom(), geom = "line", ylim = c(min(krom())-2*sd(noise()), max(krom())+4*sced()+5))
    
  })
  

  
  # Render of graph
  output$graph <- renderPlot({
   
    y <- seq(min(noise())-2*sd(noise()), max(krom())+4*sced()+5, 0.01)
    
    Dist.noise <- (1/(sd(noise())*(2*pi)^1/2))*exp(-((y-mean(noise()))^2/(2*sd(noise())^2)))
    
    Dist.peak <- (1/(sced()*(2*pi)^1/2))*exp(-((y-input$Hpeak)^2/(2*sced()^2)))
    
    
    df <- data.frame(y = y, Dist.noise = Dist.noise, Dist.peak = Dist.peak)
    
    ggplot(data = df) + 
      geom_line(aes(y = Dist.noise, x = y), colour = "red") +
      geom_line(aes(y = Dist.peak, x = y), colour = "blue") + 
      geom_line(aes(x = y, y = rep(0, length(y))), colour = "black") +
      xlim(min(krom())-2*sd(noise()), max(krom())+4*sced()+5) + 
      coord_flip()
    
    
  })
  
  # Create the text that is necessary if Homoscedastic data is chosen
  output$homoscText <- renderText({
    print("In case of homoscedastic data the standard deviation of the peak is set to be equal to the standard deviation of the noise.")
  })
  
  
}

shinyApp(ui = ui, server = server)

