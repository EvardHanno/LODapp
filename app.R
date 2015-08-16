# LoD teaching tool

## S/N values can be added to the chromatorgam

## Warnings about removed rows containing missing values come from plot being smaller then lines
## at the edges. This can be ignored. This error however seems to cause the app to recalculate everything
## several times when the inputs are changed.

library(shiny)
library(ggplot2)

ui <- fluidPage(

  
  titlePanel("Simulation of LC-MS/MS chromatogram and estimation of the resulting Limit of Detection"),
  tags$p("You can find more information about the app and help develop it further",
         tags$a(href = "https://github.com/EvardHanno/LODapp", "in GitHub"),"." ),
  
    sidebarLayout(
      sidebarPanel(
        # Data point frequency
        sliderInput("points", label = "Number of datapoints on chromatogram per unit of time:", min = 1, max = 5, value = 2.5, step = 0.5),
      
        # mean peak height
        sliderInput("Hpeak", label = "Mean height of peak:", min = 0, max = 15, value = 8, step = 0.2),
        
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
         column( 12,
             "Here there will be some text that is explaining the app - what it containes and how to use it."
         ), 
         
         column( 6, 
            plotOutput("krom")
          ),
          
          column( 6,
            plotOutput("graph"),
            tableOutput("pvalue")
          )
      )  
    )
  )
)
  
server <- function(input, output) {
  
    
  # Create reactive function of baseline
  baseline <- reactive({
    
    x <- seq(-25, 25, 1/input$points)
    
  })
  
  
  # create reactive function of noise
  noise <- reactive({
    
    for.reactivity <- c(input$Hpeak, input$peakHvar) # Added here so that noise would be recalculated if they change
    noise <- rnorm(length(baseline()), sd = input$noise)
    
  })
  
  
  # Create reactive function of chromatogam
  krom <- reactive({
   
    stdev <- 1 # This is the st deviation of the signal peak -> higher values mean wider peak
               # DO NOT CHANGE! Figure 2 calculations will become incorrect.
    
    peak <- input$Hpeak/0.3183099 * (1/(stdev*(2*pi)^1/2))*exp(-((baseline()-0)^2/(2*stdev^2)))
    
    krom <- noise() + peak
    
  })
 
 
  # Create reactiv function to create homo- or heteroscedasticity for peak height graph
  # NB! Figure 2 only depends on input values and does not depend on randomness of measurement results
  # Good because: p value is always the same with same input values
  # Bad because: not a simulated picture + in reality we never know the exact population sd and mean +
  # + because the number of measurements is not an input value
  sced <- reactive({
    if(input$selectSced == "homosc"){
      sced <- input$noise 
                          
    }
    else{
      # ERROR: peakHvar is concrete value and therefore stdev of Hpeak does not depend on
      # random measurements -> peakHvar should be multiplied with rnorm() or rt(),
      # divided by sqrt(n) [n = nr of samples; new input], multiplied by qt(0.95, df).
      # 
      # Division by sqrt(n) should be made only if LoD of method with n measurements is estimated (?).
      # Because we are looking at only one chromatorgam then this is not the case here (?).
      # 
      # Should Dist.noise and Dist.peak be calculated from t-dist not norm dist funktion?
      # Should stdev of Dist.nose be divided by sqrt(n) also? Wont peak dissapear into noise?
      # ???
      sced <- input$noise + input$peakHvar * input$Hpeak/10
    }
  })
  
  
  
  # Render of chromatogram
  output$krom <- renderPlot({
    
    df <- data.frame(baseline = baseline(), krom = krom())  
    
    ggplot(data = df) +
        geom_line(aes(x = baseline, y = krom)) +
        ylim(min(krom())-2*sd(noise()), max(krom())+4*sced()+5) +
        ggtitle("Figure 1")
    
    # qplot(baseline(), krom(), geom = "line", title = "Figure 1", ylim = c(min(krom())-2*sd(noise()), max(krom())+4*sced()+5))
    
  })
  
  # Render of Figure 2
  output$graph <- renderPlot({
   
    y <- seq(min(noise())-2*sd(noise()), max(krom())+4*sced()+5, 0.01)
    
    Dist.noise <- ( 1/(input$noise * (2 * pi)^1/2)) * exp(-( (y - 0)^2 / (2 * input$noise^2) ))
    Dist.peak <- ( 1/(sced() * (2 * pi)^1/2)) * exp(-( (y - input$Hpeak)^2 / (2 * sced()^2) ))
    

    # For creating the Critical Limit (CCa) line
    height <- seq(0, max(Dist.noise), max(Dist.noise)/10)
    Critical.limit <- qnorm(0.95) * input$noise # calculation of critical limit value
    pos <- rep(Critical.limit, length(height))
    
    
    df <- data.frame(y = y, Dist.noise = Dist.noise, Dist.peak = Dist.peak) # ggplot wants data in dataframes
    df.crit <- data.frame(pos = pos, height = height)
    
    
    # ggplot to create the Graph
    ggplot(data = df) + 
      geom_line(aes(y = Dist.noise, x = y), colour = "red") +
      geom_line(aes(y = Dist.peak, x = y), colour = "blue") + 
      geom_line(aes(x = y, y = rep(0, length(y))), colour = "black") +
      xlim(min(krom())-2*sd(noise()), max(krom())+4*sced()+5) + 
      geom_line(data = df.crit, aes(x = pos, y = height), colour = "black", linetype = "dashed", size = 1) +
      coord_flip() + 
      theme(legend.position = "none") +
      ggtitle("Figure 2")
    
    
  })
  
  # Create the text that is necessary if Homoscedastic data is chosen
  output$homoscText <- renderText({
    print("In case of homoscedastic data the standard deviation of the peak is set to be equal to the standard deviation of the noise.")
  })
  
  # Create the p value output table
  output$pvalue <- renderTable({
      
      Critical.limit <- qnorm(0.95) * input$noise # calculation of critical limit value
      
      y <- seq(min(noise())-2*sd(noise()), max(krom())+4*sced()+5, 0.01)
      Dist.peak <- (1/(sced()*(2*pi)^1/2))*exp(-((y-input$Hpeak)^2/(2*sced()^2)))
      
      df <- data.frame(y = y, Dist.peak = Dist.peak)
      
      vecBelow <- c(sum(df$Dist.peak[df$y < Critical.limit]))
      vecAll <- c(sum(df$Dist.peak))
      
      
      data.frame(Below_LC = vecBelow, All = vecAll, LC = Critical.limit, percent = vecBelow/vecAll*100)
      
  })
  
}



shinyApp(ui = ui, server = server)

