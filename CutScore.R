
library(plyr)
library(ggplot2)
library(tidyr)
library(shiny)
library(shinyjs)
library(pracma)
library(stringr)

ui = fluidPage(
  
 titlePanel("The Cut Score Operating Function (aka Grabovsky Curves)"),
  
 shinyjs::useShinyjs(),

 div(id = "myapp",
      
  h4("Please choose a method of entering examinee ability levels."),

  fluidRow(
  column(4, checkboxInput("toggleAdvanced1", "Simulate Examinees", FALSE)),
  column(4, checkboxInput("toggleAdvanced2", "Upload Examinee Theta .csv File", FALSE))
  ),
      
  shinyjs::hidden(
        div(id = "drawExaminees",
  fluidRow(
  column(4, numericInput("num_obs", "Number of examinees:", 1000)),
  column(4, sliderInput("mability", "Mean Examinee Ability:",
                                 min = -4, max = 4, value = 0.0, step= 0.01)),
  column(4, sliderInput("sdability", "Standard Deviation of Examinee Ability:",
                                 min = 0, max = 2, value = 1.0, step= 0.01))))
  ),
   
  shinyjs::hidden(
        div(id = "csvExaminees",
        fileInput('file1', 'Choose CSV File', accept='.csv'))
  ),

  h4("Please choose a method of entering examination item difficulties."),

  fluidRow(
  column(4, checkboxInput("toggleAdvanced3", "Simulate Examination Items", FALSE)),
  column(4, checkboxInput("toggleAdvanced4", "Upload Exam Difficulty .csv File", FALSE))
  ),
        
  shinyjs::hidden(
        div(id = "drawExamItems",
  fluidRow(
  column(4, numericInput("num_items", "Number of items:", 280)),
  column(4, sliderInput("itemdiff_mean", "Mean Item Difficulty:",
                                 min = -4, max = 4, value = 0.0, step= 0.01)),
  column(4, sliderInput("itemdiff_sd", "Standard Deviation of Item Difficulties:",
                                 min = 0, max = 2, value = 0.2, step= 0.01))))
  ),
      
  shinyjs::hidden(
        div(id = "csvDifficulties", 
        fileInput('file2', 'Choose CSV File', accept= '.csv'))
  ),

  fluidRow(
  column(4, sliderInput("reliability", "Exam Reliability:",
                           min=0, max=1, value=0.7, step = 0.01)),
  column(4, sliderInput("tstar", "Proposed Cut Score:",
                           min = -4, max = 4, value = 1.0, step= 0.01)),
  column(2, radioButtons("cs_meth", "Cut Score Method:",
                            c("Weighted" = "wcsf",
                              "Conditional" = "ccsf",
                              "Both" = "bothcsf"))),
  column(2, numericInput("wcec", "Weight", 1, min = 0.1, max = 10, step = 0.1))
  ),
  
  fluidRow(
  textOutput("textMin")
  
  ),
  
  
  fluidRow(
  column(4, checkboxInput("viewPlots", "View Plots"))
  ),
  
  shinyjs::hidden(
       div(id = "plotAgree",
  fluidRow(
  column(6,
  p("Please note that it takes a bit longer for the output to refresh
     while viewing plots.  Thank you for your patience. Please press the 
    \'Generate Cut Score Plots\' to proceed.")
  ),
  column(4,
       actionButton("genPlots", "Generate Cut Score Plots"))))
  ),

  fluidRow(
  column(4, actionButton("calcMin", "Calculate Minima")),
  column(4, actionButton("reset", "Reset form"))
  ),
  
  fluidRow(
    plotOutput("distPlot")
  ),

tags$head(tags$style("#textMin{font-size: 24px;}"))
))

server = function(input, output) {

  shinyjs::onclick("toggleAdvanced1",
                   c(shinyjs::toggle(id = "drawExaminees", anim = TRUE),
                     shinyjs::toggle("toggleAdvanced2")))

  shinyjs::onclick("toggleAdvanced2",
                   c(shinyjs::toggle(id = "csvExaminees", anim = TRUE),
                     shinyjs::toggle("toggleAdvanced1")))
  
  shinyjs::onclick("toggleAdvanced3",
                   c(shinyjs::toggle(id = "drawExamItems", anim = TRUE),
                     shinyjs::toggle("toggleAdvanced4")))
  
  shinyjs::onclick("toggleAdvanced4",
                   c(shinyjs::toggle(id = "csvDifficulties", anim = TRUE),
                     shinyjs::toggle("toggleAdvanced3")))  
  
  shinyjs::onclick("viewPlots",
                   shinyjs::toggle(id = "plotAgree", anim = TRUE))
  
  values <- reactiveValues()
  
  observe(values$num_itemsR <- as.numeric(input$num_items))
  observe(values$itemdiff_meanR <- as.numeric(input$itemdiff_mean))
  observe(values$itemdiff_sd <- as.numeric(input$itemdiff_sd))
  
  observe(values$num_obsR <- as.numeric(input$num_obs))
  observe(values$mabilityR <- as.numeric(input$mability))
  observe(values$sdabilityR <- as.numeric(input$sdability))

  observe(values$reliabilityR <- as.numeric(input$reliability))
  observe(values$tstarR <- as.numeric(input$tstar))
  
  observe(values$csmethod <- input$cs_meth)
  observe(values$ceweight <- input$wcec)
  
  observeEvent(input$calcMin, {
    output$textMin <- renderText({
      
      item_difficulties <- rnorm(values$num_itemsR,
                                 values$itemdiff_meanR,
                                 values$itemdiff_sd)
      
      examinee_pool <- rnorm(values$num_obsR,
                             values$mabilityR,
                             values$sdabilityR)
      
      if(values$csmethod == "wcsf"){
        
        f <- function(c, 
                      diffs,
                      reliab, 
                      examinees,
                      taustar){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- FN_result$value
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- FP_result$value
          
          CEC <- FP + FN
        }
        
        xmin <- optimize(f,
                         lower = 0,
                         upper = 1,
                         tol = 0.0001,
                         diffs = item_difficulties,
                         reliab = values$reliabilityR,
                         examinees = examinee_pool,
                         taustar = values$tstarR)
        
        best_cut <- round(xmin$minimum,4)
        
        minx <- best_cut - values$crangeR
        maxx <- best_cut + values$crangeR
        
        ggplot_result <- paste0("The weighted cut score with weight ", values$ceweight, " is ",
                                best_cut, ".")
      }
      
      if(values$csmethod == "ccsf"){
        
        ptfail <- pnorm(values$tstarR, values$mabilityR, values$sdabilityR, lower.tail = TRUE)
        ptpass  <- pnorm(values$tstarR, values$mabilityR, values$sdabilityR, lower.tail = FALSE)
        
        con_f <- function(c, 
                          diffs, 
                          reliab,
                          examinees,
                          taustar,
                          ptf,
                          ptp){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- (FN_result$value / ptp)
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- (FP_result$value / ptf)
          
          CEC <- FP + FN
        }
        
        
        conf_xmin1 <- optimize(con_f, 
                               lower = 0, 
                               upper = 1, 
                               tol = 0.0001,
                               diffs = item_difficulties,
                               reliab = values$reliabilityR,
                               examinees = examinee_pool,
                               taustar = values$tstarR,
                               ptf = ptfail,
                               ptp = ptpass)
        
        cond_min1 <- round(conf_xmin1$minimum, 4)
        
        minx <- cond_min1 - values$crangeR
        maxx <- cond_min1 + values$crangeR
        
        
        ggplot_result <- paste0("The conditional cut score is ", cond_min1, ".")
      }
      
      if(values$csmethod == "bothcsf"){
        
        ptfail <- pnorm(values$tstarR, values$mabilityR, values$sdabilityR, lower.tail = TRUE)
        ptpass  <- pnorm(values$tstarR, values$mabilityR, values$sdabilityR, lower.tail = FALSE)
        
        un_f <- function(c, 
                         diffs, 
                         reliab,
                         examinees,
                         taustar){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- FN_result$value
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- FP_result$value
          
          CEC <- FP + FN
        }
        
        con_f <- function(c, 
                          diffs, 
                          reliab,
                          examinees,
                          taustar,
                          ptf,
                          ptp){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- (FN_result$value / ptp)
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- (FP_result$value / ptf)
          
          CEC <- FP + FN
        }
        
        
        unf_xmin1 <- optimize(un_f, 
                              lower = 0, 
                              upper = 1, 
                              tol = 0.0001,
                              diffs = item_difficulties,
                              reliab = values$reliabilityR,
                              examinees = examinee_pool,
                              taustar = values$tstarR)
        
        uncond_min1 <- round(unf_xmin1$minimum, 4)
        
        un_minx <- uncond_min1 - values$crangeR
        un_maxx <- uncond_min1 + values$crangeR
        
        conf_xmin1 <- optimize(con_f, 
                               lower = 0, 
                               upper = 1, 
                               tol = 0.0001,
                               diffs = item_difficulties,
                               reliab = values$reliabilityR,
                               examinees = examinee_pool,
                               taustar = values$tstarR,
                               ptf = ptfail,
                               ptp = ptpass)
        
        cond_min1 <- round(conf_xmin1$minimum, 4)
        
        con_minx <- uncond_min1 - values$crangeR
        con_maxx <- uncond_min1 + values$crangeR
        
        ggplot_result <- paste0("The conditional cut score is ", cond_min1,
                                ", and the weighted cut score with weight ", values$ceweight, " is ",
                                uncond_min1, ".")
      }
      
      
      
      ggplot_result
      
    })
      })
  
  observeEvent(input$genPlots, {
    output$distPlot <- renderPlot({
    
      item_difficulties <- rnorm(values$num_itemsR,
                                 values$itemdiff_meanR,
                                 values$itemdiff_sd)
      
      examinee_pool <- rnorm(values$num_obsR,
                             values$mabilityR,
                             values$sdabilityR)
      
      if(values$csmethod == "wcsf"){
        
        f <- function(c, 
                      diffs,
                      reliab, 
                      examinees,
                      taustar){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- FN_result$value
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- FP_result$value
          
          CEC <- FP + FN
        }
        
        xmin <- optimize(f,
                         lower = 0,
                         upper = 1,
                         tol = 0.0001,
                         diffs = item_difficulties,
                         reliab = values$reliabilityR,
                         examinees = examinee_pool,
                         taustar = values$tstarR)
        
        best_cut <- xmin$minimum
        
        minx <- 0.5
        maxx <- 1
        
        c_vec <- seq(minx, maxx, 0.01)
        
        calculate_FPFN <- function(c, 
                                   diffs, 
                                   reliab,
                                   examinees,
                                   taustar){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- FN_result$value
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- FP_result$value
          
          CEC <- FP + FN
          
          rdf <- data.frame(c, FN, FP, CEC)
          
          return(rdf)
        }
        
        result1 <- ldply(c_vec, 
                         calculate_FPFN,
                         diffs = item_difficulties,
                         reliab = values$reliabilityR,
                         examinees = examinee_pool,
                         taustar = values$tstarR)
        
        result2 <- melt(result1, id.vars = "c")
        names(result2) <- c("c", "Statistic", "value")
        
        
        ggplot_result <- ggplot(result2, aes(x = c, y = value, color = Statistic)) + 
          geom_line(size = 2) +
          facet_grid(Statistic ~ .) + 
          theme(legend.position="none") + 
          geom_vline(xintercept=best_cut) + 
          ggtitle("The Cut Score Operating Function",
                  subtitle = paste0("Suggested Weighted Cut Score is ", round(best_cut, 4), ".")) + 
          theme(plot.title = element_text(hjust = 0.5)) + 
          xlim(minx, maxx) + 
          ylab("") +
          xlab("Percent of Items Correct on Exam")
      }
      
      if(values$csmethod == "ccsf"){
        
        ptfail <- pnorm(values$tstarR, values$mabilityR, values$sdabilityR, lower.tail = TRUE)
        ptpass  <- pnorm(values$tstarR, values$mabilityR, values$sdabilityR, lower.tail = FALSE)
        
        con_f <- function(c, 
                          diffs, 
                          reliab,
                          examinees,
                          taustar,
                          ptf,
                          ptp){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- (FN_result$value / ptp)
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- (FP_result$value / ptf)
          
          CEC <- FP + FN
        }
        
        
        conf_xmin1 <- optimize(con_f, 
                               lower = 0, 
                               upper = 1, 
                               tol = 0.0001,
                               diffs = item_difficulties,
                               reliab = values$reliabilityR,
                               examinees = examinee_pool,
                               taustar = values$tstarR,
                               ptf = ptfail,
                               ptp = ptpass)
        
        cond_min1 <- round(conf_xmin1$minimum, 4)
        
        minx <- 0.5
        maxx <- 1
        
        c_vec <- seq(minx, maxx, 0.01)
        
        calculate_FPFN <-function(c, 
                                  diffs, 
                                  reliab,
                                  examinees,
                                  taustar, 
                                  ptf, 
                                  ptp){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FNc <- FN_result$value / ptp
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FPc <- FP_result$value / ptf
          
          CECc <- FPc + FNc
          
          
          rdf <- data.frame(c, FNc, FPc, CECc)
          
          return(rdf)
        }
        
        
        result1 <- ldply(c_vec, 
                         calculate_FPFN,
                         diffs = item_difficulties,
                         reliab = values$reliabilityR,
                         examinees = examinee_pool,
                         taustar = values$tstarR,
                         ptf = ptfail,
                         ptp = ptpass)
        
        result2 <- melt(result1, id.vars = "c")
        names(result2) <- c("c", "Statistic", "value")
        
        
        
        ggplot_result <- ggplot(result2, aes(x = c, y = value, color = Statistic)) + 
          geom_line(size = 2) +
          facet_grid(Statistic ~ .) + 
          theme(legend.position="none") + 
          geom_vline(xintercept=cond_min1) + 
          ggtitle("The Cut Score Operating Function",
                  subtitle = paste0("Suggested Conditional Cut Score is ", round(cond_min1, 4), ".")) + 
          theme(plot.title = element_text(hjust = 0.5)) + 
          xlim(minx, maxx) + 
          ylab("") +
          xlab("Percent of Items Correct on Exam")
      }
      
      if(values$csmethod == "bothcsf"){
        
        ptfail <- pnorm(values$tstarR, values$mabilityR, values$sdabilityR, lower.tail = TRUE)
        ptpass  <- pnorm(values$tstarR, values$mabilityR, values$sdabilityR, lower.tail = FALSE)
        
        un_f <- function(c, 
                         diffs, 
                         reliab,
                         examinees,
                         taustar){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- FN_result$value
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- FP_result$value
          
          CEC <- FP + FN
        }
        
        con_f <- function(c, 
                          diffs, 
                          reliab,
                          examinees,
                          taustar,
                          ptf,
                          ptp){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- (FN_result$value / ptp)
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- (FP_result$value / ptf)
          
          CEC <- FP + FN
        }
        
        
        unf_xmin1 <- optimize(un_f, 
                              lower = 0, 
                              upper = 1, 
                              tol = 0.0001,
                              diffs = item_difficulties,
                              reliab = values$reliabilityR,
                              examinees = examinee_pool,
                              taustar = values$tstarR)
        
        uncond_min1 <- round(unf_xmin1$minimum, 4)
        
        un_minx <- 0.5
        un_maxx <- 1
        
        un_c_vec <- seq(un_minx, un_maxx, 0.01)
        
        
        calculate_FPFN <- function(c, 
                                   diffs, 
                                   reliab,
                                   examinees,
                                   taustar){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- FN_result$value
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- FP_result$value
          
          CEC <- FP + FN
          
          rdf <- data.frame(c, FN, FP, CEC)
          
          return(rdf)
        }
        
        unresult <- ldply(un_c_vec, 
                          calculate_FPFN,
                          diffs = item_difficulties,
                          reliab = values$reliabilityR,
                          examinees = examinee_pool,
                          taustar = values$tstarR)
        
        us1 <- melt(unresult, id.vars = "c")
        names(us1) <- c("c", "Statistic", "value")
        us1$Mode <- "Weighted"
        us1$Minima <- uncond_min1
        
        conf_xmin1 <- optimize(con_f, 
                               lower = 0, 
                               upper = 1, 
                               tol = 0.0001,
                               diffs = item_difficulties,
                               reliab = values$reliabilityR,
                               examinees = examinee_pool,
                               taustar = values$tstarR,
                               ptf = ptfail,
                               ptp = ptpass)
        
        cond_min1 <- round(conf_xmin1$minimum, 4)
        
        con_minx <- 0.5
        con_maxx <- 1
        
        con_c_vec <- seq(con_minx, con_maxx, 0.01)
        
        calculate_FPFN_con <-function(c, 
                                      diffs, 
                                      reliab,
                                      examinees,
                                      taustar, 
                                      ptf, 
                                      ptp){
          
          reliability <- reliab
          
          m_ability <- mean(examinees)
          theta_var <- var(examinees)
          v_ability <- theta_var - ((1 - reliability)/reliability)
          sigma_ability <- sqrt(v_ability)
          
          itemdiffs <- diffs
          
          # Using the above information we can now calculate the probability
          # of a false negative with the iteration-specific "c" value.
          prob_FN <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * PHI * phi
          }
          
          # Vectorizing the above function in order to pass it to "integrate"
          vFN <- Vectorize(prob_FN, "tau")
          FN_result <- integrate(vFN, lower = taustar, upper = 4)
          FN <- FN_result$value / ptp
          
          
          prob_FP <- function(tau){
            
            workdf <- data.frame(itemdiffs)
            workdf$ptb <- 1 / (1 + exp(itemdiffs - tau))
            workdf$mult <- workdf$ptb*(1-workdf$ptb)
            p_bar <- mean(workdf$ptb)
            SNtau <- sqrt(sum(workdf$mult))
            
            PHI_x <- nrow(workdf) * (c - p_bar) / SNtau
            
            phi_x <- (tau - m_ability) / sigma_ability
            
            PHI <- (1/2) * (1 + erf(PHI_x / sqrt(2)))
            
            phi <- (1 / (sqrt(2*pi))) * exp(-phi_x^2 / 2)
            
            f <- (1/sigma_ability) * (1-PHI) * phi
          }
          
          vFP <- Vectorize(prob_FP, "tau")
          FP_result <- integrate(vFP, lower = -4, upper = taustar)
          FP <- FP_result$value / ptf
          
          CEC <- FP + FN
          
          
          
          rdf <- data.frame(c, FN, FP, CEC)
          
          return(rdf)
        }
        
        conresult <- ldply(con_c_vec, 
                           calculate_FPFN_con,
                           diffs = item_difficulties,
                           reliab = values$reliabilityR,
                           examinees = examinee_pool,
                           taustar = values$tstarR,
                           ptf = ptfail,
                           ptp = ptpass)
        
        cs1 <- melt(conresult, id.vars = "c")
        names(cs1) <- c("c", "Statistic", "value")
        cs1$Mode <- "Conditional"
        cs1$Minima <- cond_min1
        
        ts1 <- rbind(us1, cs1)
        
        ggplot_result <- ggplot(ts1, aes(x = c, y = value, color = Statistic)) + geom_line(size = 2) + 
          facet_grid(Statistic ~ Mode, scales = "free") + theme(legend.position="none") + 
          ggtitle("The Cut Score Operating Function",
                  subtitle = paste0("Suggested Cut Score are ", round(cond_min1, 4), " (conditional) and ", round(uncond_min1, 4), " (weighted).")) + 
          theme(plot.title = element_text(hjust = 0.5)) + 
          ylab("Percent of Examinees Misclassified") +
          xlab("Percent of Items Correct on Exam") + 
          geom_vline(data = ts1, aes(xintercept = as.numeric(Minima)))
      }
      
      ggplot_result
      
    })})
  
  observeEvent(input$reset, {
    shinyjs::reset("myapp")
    "drawExaminees" = FALSE
    shinyjs::show("toggleAdvanced1")
    shinyjs::hide("drawExaminees")
    
    "csvExaminees" = FALSE
    shinyjs::show("toggleAdvanced2")
    shinyjs::hide("csvExaminees")
    
    "drawExamItems" = FALSE
    shinyjs::show("toggleAdvanced3")
    shinyjs::hide("drawExamItems")
    
    "csvDifficulties" = FALSE
    shinyjs::show("toggleAdvanced4")
    shinyjs::hide("csvDifficulties")
    
    "viewPlots" = FALSE
    shinyjs::hide("plotAgree")

  })

}

shinyApp(ui, server)