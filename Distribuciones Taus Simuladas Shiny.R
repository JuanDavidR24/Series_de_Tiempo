#install.packages("shiny")
library(shiny)

ui <- fluidPage(
  
  titlePanel(strong("Simulación Pruebas Tau's")),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4(strong("Definición de Parámetros")),
      
      sliderInput(inputId = "nObs.input",
                  label = "Tamaño de la muestra:",
                  value = 250,
                  min = 0,
                  max = 1000),
      
      sliderInput(inputId = "itr.input",
                   label = "Numero de Iteraciones",
                   value = 200,
                   min = 100,
                   max = 1000),
      br(),
      h5("By Juan David Rincón, 2022")
    ),
    
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  tabPanel(title = strong("Tau Tau"), 
                           plotOutput(outputId = "tautau.output"),
                           textOutput(outputId = "TTprct10.output"),
                           textOutput(outputId = "TTprct5.output"),
                           textOutput(outputId = "TTprct1.output")),
                  tabPanel(title = strong("Tau Mu"),
                           plotOutput(outputId = "taumu.output"),
                           textOutput(outputId = "TMprct10.output"),
                           textOutput(outputId = "TMprct5.output"),
                           textOutput(outputId = "TMprct1.output")),
                  tabPanel(title = strong("Tau"),
                           plotOutput(outputId = "tau.output"),
                           textOutput(outputId = "Tprct10.output"),
                           textOutput(outputId = "Tprct5.output"),
                           textOutput(outputId = "Tprct1.output")))
      
    )
  ),
)

server <- function(input, output, session){
  
  # Funciones.
  tau.tau.f <- function(N, itr, y0, alpha){

    vectorTaus <- c()

    for(p in 1:itr){

      at_sim <- rnorm(N, 0, 1)

      vectorYs <- c(y0)
      for(i in 2:N){
        vectorYs[i] <- vectorYs[i-1] + at_sim[i]
      }

      Nquit <- floor(N*0.2)

      vectorYs <- vectorYs[-(1:Nquit)]

      N.new <- N - Nquit

      Y <- as.matrix(vectorYs[-1]-vectorYs[-N.new])
      X <- as.matrix(cbind(1, c(1:(N.new-1)), vectorYs[-N.new]))

      B_gorro <- solve(t(X)%*%X)%*%t(X)%*%Y

      e_gorro <- Y - (X%*%B_gorro)

      sigma2_gorro <- as.numeric((t(e_gorro)%*%e_gorro)/(N.new-3))

      varcovB <- sigma2_gorro*solve(t(X)%*%X)

      sdBetas <- sqrt(diag(varcovB))

      tau <- B_gorro[3]/(sdBetas[3])

      vectorTaus[p] <- tau
    }

    VC <- quantile(vectorTaus, alpha)

    return(list(VC=VC, taus=vectorTaus))
  }
  tau.mu.f <- function(N, itr, y0, alpha){

    vectorTaus <- c()

    for(p in 1:itr){

      at_sim <- rnorm(N, 0, 1)

      vectorYs <- c(y0)
      for(i in 2:N){
        vectorYs[i] <- vectorYs[i-1] + at_sim[i]
      }

      Nquit <- floor(N*0.2)

      vectorYs <- vectorYs[-(1:Nquit)]

      N.new <- N - Nquit

      Y <- as.matrix(vectorYs[-1]-vectorYs[-N.new])
      X <- as.matrix(cbind(1, vectorYs[-N.new]))

      B_gorro <- solve(t(X)%*%X)%*%t(X)%*%Y

      e_gorro <- Y - (X%*%B_gorro)

      sigma2_gorro <- as.numeric((t(e_gorro)%*%e_gorro)/(N.new-2))

      varcovB <- sigma2_gorro*solve(t(X)%*%X)

      sdBetas <- sqrt(diag(varcovB))

      tau <- B_gorro[2]/(sdBetas[2])

      vectorTaus[p] <- tau
    }

    VC <- quantile(vectorTaus, alpha)

    return(list(VC=VC, taus=vectorTaus))
  }
  tau.f <- function(N, itr, y0, alpha){

    vectorTaus <- c()

    for(p in 1:itr){

      at_sim <- rnorm(N, 0, 1)

      vectorYs <- c(y0)
      for(i in 2:N){
        vectorYs[i] <- vectorYs[i-1] + at_sim[i]
      }

      Nquit <- floor(N*0.2)

      vectorYs <- vectorYs[-(1:Nquit)]

      N.new <- N - Nquit

      Y <- as.matrix(vectorYs[-1]-vectorYs[-N.new])
      X <- as.matrix(vectorYs[-N.new])

      B_gorro <- solve(t(X)%*%X)%*%t(X)%*%Y

      e_gorro <- Y - (X%*%B_gorro)

      sigma2_gorro <- as.numeric((t(e_gorro)%*%e_gorro)/(N.new-1))

      varcovB <- sigma2_gorro*solve(t(X)%*%X)

      sdBetas <- sqrt(diag(varcovB))

      tau <- B_gorro[1]/(sdBetas[1])

      vectorTaus[p] <- tau
    }

    VC <- quantile(vectorTaus, alpha)

    return(list(VC=VC, taus=vectorTaus))
  }
  
  # Inputs.
  Nobs <- reactive(input$nObs.input)
  Itr <- reactive(input$itr.input)
  
  # Plots.
  plotTauTau <- reactive(tau.tau.f(Nobs(), Itr(), 0, 0.05)$taus)
  
  # Values.
  
  # Tau Tau.
  value5TauTau <- reactive(round(as.numeric(tau.tau.f(Nobs(), Itr(), 0, 0.05)$VC), 4))
  value10TauTau <- reactive(round(as.numeric(tau.tau.f(Nobs(), Itr(), 0, 0.1)$VC), 4))
  value1TauTau <- reactive(round(as.numeric(tau.tau.f(Nobs(), Itr(), 0, 0.01)$VC), 4))
  
  # Tau Mu.
  value5TauMu <- reactive(round(as.numeric(tau.mu.f(Nobs(), Itr(), 0, 0.05)$VC), 4))
  value10TauMu <- reactive(round(as.numeric(tau.mu.f(Nobs(), Itr(), 0, 0.1)$VC), 4))
  value1TauMu <- reactive(round(as.numeric(tau.mu.f(Nobs(), Itr(), 0, 0.01)$VC), 4))
  
  # Tau.
  value5Tau <- reactive(round(as.numeric(tau.f(Nobs(), Itr(), 0, 0.05)$VC), 4))
  value10Tau <- reactive(round(as.numeric(tau.f(Nobs(), Itr(), 0, 0.1)$VC), 4))
  value1Tau <- reactive(round(as.numeric(tau.f(Nobs(), Itr(), 0, 0.01)$VC), 4))
  
  # Outputs Value.
  
  # Tau Tau.
  output$TTprct5.output <- renderText(paste0("Valor Critico 5%: ", value5TauTau()))
  output$TTprct10.output <- renderText(paste0("Valor Critico 10%: ", value10TauTau()))
  output$TTprct1.output <- renderText(paste0("Valor Critico 1%: ", value1TauTau()))

  # Tau Mu.
  output$TMprct5.output <- renderText(paste0("Valor Critico 5%: ", value5TauMu()))
  output$TMprct10.output <- renderText(paste0("Valor Critico 10%: ", value10TauMu()))
  output$TMprct1.output <- renderText(paste0("Valor Critico 1%: ", value1TauMu()))
  
  # Tau.
  output$Tprct5.output <- renderText(paste0("Valor Critico 5%: ", value5Tau()))
  output$Tprct10.output <- renderText(paste0("Valor Critico 10%: ", value10Tau()))
  output$Tprct1.output <- renderText(paste0("Valor Critico 1%: ", value1Tau()))
  
  # Outpluts Plots.
  
  # Tau Tau.
  output$tautau.output <- renderPlot({
    plot(density(tau.tau.f(Nobs(), Itr(), 0, 0.05)$tau),
         main = "Distribución Tau Tau Simulada")
  })
  
  # Tau Mu.
  output$taumu.output <- renderPlot({
    plot(density(tau.mu.f(Nobs(), Itr(), 0, 0.05)$tau),
         main = "Distribución Tau Mu Simulada")
  })
  
  # Tau.
  output$tau.output <- renderPlot({
    plot(density(tau.f(Nobs(), Itr(), 0, 0.05)$tau),
         main = "Distribución Tau Simulada")
  })
}

shinyApp(ui, server)
