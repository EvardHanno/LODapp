# LoD teaching tool

## Warnings about removed rows containing missing values come from plot being smaller then lines
## at the edges. This can be ignored. This error however seems to cause the app to recalculate everything
## several times when the inputs are changed.

library(shiny)
library(ggplot2)

ui <- fluidPage(

  
  titlePanel(
      fluidRow( 
          column(9,
            h3("Simulation of LC-MS/MS chromatogram and estimation of the resulting Limit of Detection")
          ),
          column(3,
            img(src = "TY_logo.png", height = 30, width = 257.73)
          )
      ),
      
             windowTitle = "LODapp"),
         
         tags$p("More information about the app (including the script) can be found",
         tags$a(href = "https://github.com/EvardHanno/LODapp", "in GitHub"),
         ". Any suggestions to improve the app are always welcome." ),
         
    sidebarLayout(
      sidebarPanel(
        # Data point frequency
        sliderInput("points", label = "Number of datapoints on chromatogram per unit of time:", min = 2, max = 5, value = 2.5, step = 0.25),
      
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
          sliderInput("peakHvar", label = "Slope/strength of heteroscedasticity:", min = 0, max = 4, value = 1, step = 0.25)
        ),
      width = 3),

      ## Main panel outputs
  
      mainPanel(
        tabsetPanel(
            tabPanel("Figures",
            fluidRow(
             column( 12,
               p(""),
                     
               p("In this app an LC-MS/MS chromatogram is simulated (Figure 1) and the variance of noise
                and peaks are shown (Figure 2). The user can change different parameters of the chromatogram.
                New Figures are created interacively."),
               
               p("It must be noted that any comparsions to real chromatograms should not be made for estimating LoD-s.
                The data here is simulated and makes assumptions (e.g. normality) and therefore will not correspond accurately to
                real life measurements."),
               
               p("A more thorough explanation can be found in the Description tab.")
               
            ), 
             
             column( 6, 
                plotOutput("krom"),
                tableOutput("SN"),
                
                p("In the chromatorgam the height of the peak is set to have the 
                mean value of the distribution of peak heights (blue line in Figure 2)."),
                
                p("Signal-to-noise ratios (S/N) estimated from peak-to-peak noise (SN.ptp) and by 
                  standard deviation of the noise (SN.rms) can be found in the table above.")
                
                
              ),
              
              column( 6,
                plotOutput("graph"),
                tableOutput("pvalue"),
                
                p("Figure 2 shows the distribution of height values of the whole population. The red line represents the
                  distribution of noise and the blue line represents the distribution of peak heights."),
                  
                p("The dashed line represents the 95% quantile of the noise distribution. This means that 95% of the 
                  data points (excluding the peak) are below this line."),
                
                p("The table shows the whole area value of the peak distribution (Area.all),
                  the decision threshold (LC), 
                  the peak distribution area below the LC (Area.BelowLC),
                  and from these values the p value is calculated which shows the percent of area below")
                
              )
            )
            ),
            
            tabPanel("Description",
                
                h4("Figure 1."),
                
                p("Noise in the chromatogram is created by the random number generator in R."),
                
                p("Signal-to-noise ratios (S/N) calculated for the simulated peak by different approaches can be found 
                  in the table below Figure 1. Signal strength is taken as the height of the peak for both cases.
                  Noise is estimated with the following procedures:"),
                
                p("In the S/N peak-to-peak approach (SN.ptp) N is taken as the height difference between the 
                  lowest and highest points of noise."),
                p("In the S/N root mean square approach (SN.rms) N is calculated as the standard deviation
                  of the noise."),
                  
                p("The comparison between the S/N values and the p-values (to estimate LOD) in the table below Figure 2
                  can not be compared to real experimental situations due to assuptions made in the simulations."),
                     
                
                h4("Figure 2."),
                  
                p("Note that usually in LC analysis areas not heights of peaks are used for data interpretation. This is important because
                  taking an area of a peak uses information from more than one data point. However, here when calculating the p values
                  (and the S/N values) only the height of one data point is used.   
                  The result is that for a specific set of parameters visually the peak may not seem separated from the noise but the S/N value and
                  the p value suggest the peak is at LOD. This is especally so when the number of datapoints is set to a low value.
                  However, when a high number of datapoints is set the peak can seem visually clearly present although the 
                  S/N and p values are still poor."),

                p("Moreover, generally in analytical chemistry 
                  the accurate parameters of populations are not known and can only be estimated from measured data.")
                     
                     
            )
        )
    )
  )
)
  
server <- function(input, output) {
  
    
  # Create reactive function of baseline
  baseline <- reactive({
    
    x <- seq(0, 50, 1/input$points)
    
  })
  
  
  # create reactive function of noise
  noise <- reactive({
    
    for.reactivity <- c(input$Hpeak, input$peakHvar) # Added here so that noise would be recalculated if they change
    noise <- rnorm(length(baseline()), sd = input$noise, mean = 50)
    
  })
  
  
  # Create reactive function of chromatogam
  krom <- reactive({
   
    stdev <- 1 # This is the st deviation of the signal peak -> higher values mean wider peak
               # DO NOT CHANGE! Figure 2 calculations will become incorrect.
    
    peak <- input$Hpeak/0.3183099 * (1/(stdev*(2*pi)^1/2))*exp(-((baseline()-25)^2/(2*stdev^2)))
    
    krom <- noise() + peak
    
  })
 
 
  # Create reactiv function to create homo- or heteroscedasticity for peak height graph
  # NB! Figure 2 only depends on input values and does not depend on randomness of measurement results
  # Good because: p value is always the same with same input values
  # Bad because: not a simulated picture, in reality we never know the exact population sd and mean +
  # + because the number of measurements is not an input value
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
    
    df <- data.frame(baseline = baseline(), krom = krom())  
    
    ggplot(data = df) +
        geom_line(aes(x = baseline, y = krom)) +
        ylim(min(krom())-2*sd(noise()), max(krom())+4*sced()+5) +
        ggtitle("Figure 1") + 
        xlab("Time") + 
        ylab("Intensity")
    
    # qplot(baseline(), krom(), geom = "line", title = "Figure 1", ylim = c(min(krom())-2*sd(noise()), max(krom())+4*sced()+5))
    
  })
  
  # Render of Figure 2
  output$graph <- renderPlot({
   
    y <- seq(min(noise())-2*sd(noise()), max(krom())+4*sced()+5, 0.01)
    
    Dist.noise <- ( 1/(input$noise * (2 * pi)^1/2)) * exp(-( (y - 50)^2 / (2 * input$noise^2) ))
    Dist.peak <- ( 1/(sced() * (2 * pi)^1/2)) * exp(-( (y - (input$Hpeak + 50))^2 / (2 * sced()^2) ))
    

    # For creating the Critical Limit (CCa) line
    height <- seq(0, max(Dist.noise), max(Dist.noise)/10)
    Critical.limit <- (qnorm(0.95) * input$noise) + 50 # calculation of critical limit value
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
      ggtitle("Figure 2") + 
      xlab("Height values") + 
      ylab("Probability")
    
  })
  
  # Create the text that is necessary if Homoscedastic data is chosen
  output$homoscText <- renderText({
    print("In case of homoscedastic data the standard deviation of the peak is set to be equal to the standard deviation of the noise.")
  })
  
  # Create the p value output table
  output$pvalue <- renderTable({
      
      # test the quantile() function for this!
      Critical.limit <- (qnorm(0.95) * input$noise) + 50 # calculation of critical limit value
      
      y <- seq(min(noise())-2*sd(noise()), max(krom())+4*sced()+5, 0.01)
      Dist.peak <- ( 1/(sced() * (2 * pi)^1/2)) * exp(-( (y - (input$Hpeak + 50))^2 / (2 * sced()^2) ))
      
      df <- data.frame(y = y, Dist.peak = Dist.peak)
      
      vecBelow <- c(sum(df$Dist.peak[df$y < Critical.limit]))
      vecAll <- c(sum(df$Dist.peak))
      
      
      data.frame(Area.all = vecAll, LC = Critical.limit, Area.belowLC = vecBelow, p.value = vecBelow/vecAll*100)
      
  })
  
  
  # Create the S/N value
  output$SN <- renderTable({
      
      krom <- krom()
      x <- baseline()
      
      # arvutame punktide tiheduse (vajalik backgroundi andmepunktide määramisel)
      points <- 1/ input$points
      
      # eraldame kromatogrammist backgroundi osa
      kat1 <- c(  (abs((min(x) - 0)/points) + 1) : (abs((min(x) - 15)/points) + 1) )
      kat2 <- c( (abs((min(x) - 35)/points) + 1) : (abs((min(x) - 50)/points) + 1) )
      
      background = c( krom[kat1], krom[kat2] )

      n1 = sd(background)
      n2 = max(background) - min(background)
      
      s = max(krom) - abs(mean(background))

      data.frame(SN.ptp = s/n2, SN.rms = s/n1)
      
      
  })
  
  
}



shinyApp(ui = ui, server = server)

