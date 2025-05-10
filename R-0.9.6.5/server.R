mlemurServer <- function(input, output, session) {
  rv <- reactiveValues()
  rv$BatchFirstRun <- FALSE
  rv$BatchDatasets <- NA
  
  #### Rate Validate Data ####
  rv$SettingsRateUserInput <- callModule(settingsPlating, id = "SettingsRate")
  rv$CountsRateUserInput <- callModule(countsPlating, id = "CountsRate", userSettings=reactive(rv$SettingsRateUserInput()), stack_cols=FALSE)
  loadDatasetWrapper("CountsRate", rv = rv, SettingsPlatingID = "SettingsRate")
  
  #### Rate Calculate ####
  observeEvent(input$calculate, {
    userData <- rv$CountsRateUserInput()
    if (userData$errors == TRUE) {
      output$errorBarRate <- renderText("There are problems with your data that need to be resolved. If you have trouble, check Help.")
      shinyjs::hide("advanced")
      shinyjs::enable("calculate")
      shinyjs::runjs("$('#calculate').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    } else {
      if (userData$setPerPlate == 1) outputData <- calc.rate.int(userData)[1:14] else outputData <- calc.rate.per.plate.int(userData)[1:14]
      results <- data.frame("Calculations" = outputData,
                            row.names = c("&epsilon;", "N<sub>t</sub>", "CV", "&rho;", "&lambda;", "d", ifelse(userData$setPerPlate == 1, "m<sub>p</sub>", "&mu;<sub>p</sub>"), "&phi;",
                                          "m", "m<sup>95&percnt;&ndash;</sup>", "m<sup>95&percnt;&plus;</sup>", "&mu;", "&mu;<sup>95&percnt;&ndash;</sup>", "&mu;<sup>95&percnt;&plus;</sup>"),
                            check.names = FALSE)
      output$tableRate <- reactable::renderReactable({
        reactable::reactable(results,
                             rownames = TRUE,
                             sortable = FALSE,
                             pagination = FALSE,
                             defaultColDef = reactable::colDef(html = TRUE))
      })
      rv$RatetoClipboard <- outputData
      output$errorBarRate <- renderText("")
      shinyjs::show("advanced")
      shinyjs::enable("calculate")
      shinyjs::runjs("$('#calculate').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    }
    gc()
  })
  
  #### Rate Copy to Clipboard ####
  observeEvent(input$clip, {
    clipr::write_clip(paste(as.character(rv$RatetoClipboard), collapse = "\n"), allow_non_interactive = TRUE)
  })
  
  #### Rate Add to Report ####
  rv$Ratereport <- data.frame(
    "Parameters" = c("eff", "Nt", "CV", "Fitness", "Lag", "Death rate", "Residual m", "Rel size of inoc", "m", "m_low", "m_up", "mu", "mu_low", "mu_up")
  )
  
  observeEvent(input$appendToReport, {
    addedSet <- data.frame(c(rv$RatetoClipboard), check.names = FALSE)
    names(addedSet) <- c(input$datasetName)
    rv$Ratereport <- cbind(rv$Ratereport, addedSet)
    shinyFeedback::showToast(type = "success", message = "Appended succesfully", .options = list(positionClass = "toast-bottom-right", progressBar = FALSE, timeOut = 1500))
  })
  
  #### Rate Download Report ####
  output$dl <- downloadHandler(
    filename = function() {
      paste("Report-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(rv$Ratereport, path = file)
    }
  )
  
  #### Rate Erase ####
  observeEvent(input$erase, {
    callModule(countsPlatingUpdate, "CountsRate", usePreset = 0)
  })
  
  #### Rate Load Sample ####
  observeEvent(input$sample, {
    callModule(countsPlatingUpdate, "CountsRate", usePreset = 1)
  })
  
  #### P value Validate Data ####
  rv$SettingsPvalUserInput <- callModule(settingsPlating, id = "SettingsPval")
  rv$CountsStrain1UserInput <- callModule(countsPlating, id = "CountsStrain1", userSettings=reactive(rv$SettingsPvalUserInput()), stack_cols=TRUE)
  rv$CountsStrain2UserInput <- callModule(countsPlating, id = "CountsStrain2", userSettings=reactive(rv$SettingsPvalUserInput()), stack_cols=TRUE)
  loadDatasetWrapper("CountsStrain1", rv = rv, SettingsPlatingID = "SettingsPval")
  loadDatasetWrapper("CountsStrain2", rv = rv, SettingsPlatingID = "SettingsPval")
  
  #### P value Calculate ####
  observeEvent(input$calculate2, {
    userDataStrain1 <- rv$CountsStrain1UserInput()
    userDataStrain2 <- rv$CountsStrain2UserInput()
    if (userDataStrain1$errors == TRUE || userDataStrain2$errors == TRUE) {
      output$errorBarPvalue <- renderText("There are problems with your data that need to be resolved. If you have trouble, check Help.")
      shinyjs::hide("advanced2")
      shinyjs::enable("calculate2")
      shinyjs::runjs("$('#calculate2').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    } else {
      if (userDataStrain1$setPerPlate == 1) Strain1 <- calc.rate.int(userDataStrain1)[1:14] else Strain1 <- calc.rate.per.plate.int(userDataStrain1)[1:14]
      if (userDataStrain2$setPerPlate == 1) Strain2 <- calc.rate.int(userDataStrain2)[1:14] else Strain2 <- calc.rate.per.plate.int(userDataStrain2)[1:14]
      results2 <- data.frame("Strain 1" = Strain1[1:14],
                             "Strain 2" = Strain2[1:14],
                             row.names = c("&epsilon;", "N<sub>t</sub>", "CV", "&rho;", "&lambda;", "d", ifelse(userDataStrain1$setPerPlate == 1, "m<sub>p</sub>", "&mu;<sub>p</sub>"), "&phi;",
                                           "m", "m<sup>95&percnt;&ndash;</sup>", "m<sup>95&percnt;&plus;</sup>", "&mu;", "&mu;<sup>95&percnt;&ndash;</sup>", "&mu;<sup>95&percnt;&plus;</sup>"),
                             check.names = FALSE)
      output$tablePvalue <- reactable::renderReactable({
        reactable::reactable(
          results2,
          rownames = TRUE,
          sortable = FALSE,
          wrap = FALSE,
          pagination = FALSE,
          defaultColDef = reactable::colDef(html = TRUE))
      })
      if (userDataStrain1$setPerPlate == 1) statistics <- calc.pval.int(userDataStrain1, userDataStrain2, as.numeric(Strain1[15]), as.numeric(Strain2[15])) else statistics <- calc.pval.per.plate.int(userDataStrain1, userDataStrain2)
      pvalue <- statistics[1]
      output$pvalueinfo <- renderText({
        if (is.na(pvalue) == TRUE | is.nan(pvalue) == TRUE) {
          Outp <- paste("impossible to calculate using provided datasets.")
        }
        else if (pvalue >= 0.05) {
          Outp <- paste( signif(pvalue, 4), "&mdash; the difference is not significant at &alpha;&equals;0.05.")
        }
        else {
          Outp <- paste( signif(pvalue, 4), "&mdash; the difference is significant at &alpha;&equals;0.05.")
        }
        paste("<em>", "P", "</em> value is ", Outp)
      })
      rv$Strain1ForPower <- as.character(c(length(userDataStrain1$CountsSelective), Strain1[12], Strain1[4], Strain1[1], Strain1[2], ifelse(is.na(suppressWarnings(as.numeric(Strain1[3]))), 0, Strain1[3]), Strain1[5], Strain1[7], Strain1[6], ifelse(is.numeric(userDataStrain1$Inoculum), userDataStrain1$Inoculum, 0), userDataStrain1$useLagFitness))
      rv$Strain2ForPower <- as.character(c(length(userDataStrain2$CountsSelective), Strain2[12], Strain2[4], Strain2[1], Strain2[2], ifelse(is.na(suppressWarnings(as.numeric(Strain2[3]))), 0, Strain2[3]), Strain2[5], Strain2[7], Strain2[6], ifelse(is.numeric(userDataStrain2$Inoculum), userDataStrain2$Inoculum, 0), userDataStrain2$useLagFitness))
      output$errorBarPvalue <- renderText("")
      shinyjs::show("advanced2")
      shinyjs::enable("calculate2")
      shinyjs::runjs("$('#calculate2').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    }
    gc()
  })
  
  #### P value Erase ####
  observeEvent(input$erase2, {
    callModule(countsPlatingUpdate, "CountsStrain1", usePreset = 0)
    callModule(countsPlatingUpdate, "CountsStrain2", usePreset = 0)
  })
  
  #### P value Load Sample ####
  observeEvent(input$sample2, {
    callModule(countsPlatingUpdate, "CountsStrain1", usePreset = 1)
    callModule(countsPlatingUpdate, "CountsStrain2", usePreset = 2)
  })
  
  #### P value Send for power analysis ####
  observeEvent(input$sendToPower, {
    shinyWidgets::updateAwesomeRadio(session = session, inputId = "sampleSizeOrPower", selected = 0)
    callModule(PowerModuleUpdate, "PowerStrain1", usePreset = rv$Strain1ForPower)
    callModule(PowerModuleUpdate, "PowerStrain2", usePreset = rv$Strain2ForPower)
    shinyjs::runjs("SwitchTabToPower();")
  })
  
  #### Corrector Data Validation ####
  observe({
    rv$enteredPvalues <- numerise(input$enteredPvalues)
    rv$selectedMethod <- input$correctionMethod
    
    rv$errorCorrector <- 0
    
    if (PvalueValidator(rv$enteredPvalues, rv$selectedMethod) != "") {
      rv$errorCorrector <- 1
      textInputError(inputId = "enteredPvalues", text = paste(PvalueValidator(rv$enteredPvalues, rv$selectedMethod)))
    } else {
      shinyFeedback::hideFeedback(inputId = "enteredPvalues")
    }
    
  })
  
  #### Corrector Error Bar ####
  observeEvent(input$calculate3, {
    if (rv$errorCorrector == 1) {
      output$errorBarCorrector <- renderText("There are problems with your data that need to be resolved. If you have trouble, check Help.")
      shinyjs::hide("advanced3")
      #### Corrector Calculator ####
    } else {
      output$errorBarCorrector <- renderText("")
      correctedPvalues <- p.adjust(p=rv$enteredPvalues, method=rv$selectedMethod)
      results3 <-
        data.frame(
          "Entered <em>P</em> values" = as.character(rv$enteredPvalues),
          "Corrected <em>P</em> values" = as.character(correctedPvalues),
          check.names = FALSE
        )
      rv$oldPvalues <- c("Entered P values", as.character(rv$enteredPvalues))
      rv$newPvalues <- c("Corrected P values", as.character(correctedPvalues))
      output$tableCorrector <- reactable::renderReactable({
        reactable::reactable(
          results3,
          sortable = FALSE,
          wrap = TRUE,
          defaultColDef = reactable::colDef(html = TRUE)
        )
      })
      shinyjs::show("advanced3")
      rv$PvaltoClipboard <- paste(rv$oldPvalues, rv$newPvalues, sep = "\t", collapse = "\n")
    }
  })
  
  #### Corrector Clipboard ####
  observeEvent(input$clip2, {
    clipr::write_clip(paste(as.character(rv$PvaltoClipboard), collapse = "\n"), allow_non_interactive = TRUE)
  })
  
  #### Corrector Eraser ####
  observeEvent(input$erase3, {
    updateTextAreaInput(session, inputId = "enteredPvalues", value = "")
  })
  
  #### Corrector Sample ####
  observeEvent(input$sample3, {
    updateTextAreaInput(session, inputId = "enteredPvalues", value = "0.064\n0.00007\n0.002")
  })
  
  #### Fold Disable per-plate ####
  shinyjs::disable(id = "SettingsFold-setPerPlate")
  
  #### Fold Validate Data ####
  rv$SettingsFoldUserInput <- callModule(settingsPlating, id = "SettingsFold")
  rv$CountsStrain1FoldUserInput <- callModule(countsPlating, id = "CountsStrain1Fold", userSettings=reactive(rv$SettingsFoldUserInput()), stack_cols=TRUE)
  rv$CountsStrain2FoldUserInput <- callModule(countsPlating, id = "CountsStrain2Fold", userSettings=reactive(rv$SettingsFoldUserInput()), stack_cols=TRUE)
  rv$CountsStrain3FoldUserInput <- callModule(countsPlating, id = "CountsStrain3Fold", userSettings=reactive(rv$SettingsFoldUserInput()), stack_cols=TRUE)
  rv$CountsStrain4FoldUserInput <- callModule(countsPlating, id = "CountsStrain4Fold", userSettings=reactive(rv$SettingsFoldUserInput()), stack_cols=TRUE)
  rv$CountsStrain5FoldUserInput <- callModule(countsPlating, id = "CountsStrain5Fold", userSettings=reactive(rv$SettingsFoldUserInput()), stack_cols=TRUE)
  rv$CountsStrain6FoldUserInput <- callModule(countsPlating, id = "CountsStrain6Fold", userSettings=reactive(rv$SettingsFoldUserInput()), stack_cols=TRUE)
  loadDatasetWrapper("CountsStrain1Fold", rv = rv, SettingsPlatingID = "SettingsFold")
  loadDatasetWrapper("CountsStrain2Fold", rv = rv, SettingsPlatingID = "SettingsFold")
  loadDatasetWrapper("CountsStrain3Fold", rv = rv, SettingsPlatingID = "SettingsFold")
  loadDatasetWrapper("CountsStrain4Fold", rv = rv, SettingsPlatingID = "SettingsFold")
  loadDatasetWrapper("CountsStrain5Fold", rv = rv, SettingsPlatingID = "SettingsFold")
  loadDatasetWrapper("CountsStrain6Fold", rv = rv, SettingsPlatingID = "SettingsFold")
  
  observe({
    if (!paste(input$FoldEquation, collapse = "") %in% c("Strain 1/Strain 2", "Strain 1-Strain 2", "(Strain 1/Strain 2)/(Strain 3/Strain 4)", "(Strain 1-Strain 3)/(Strain 2-Strain 3)", "(Strain 1+Strain 2)/Strain 3")) updateSelectInput(inputId = "FoldUsePreset", selected = 0)
    rv$FoldEquationVector <- gsub(pattern = "Strain ", replacement = "X", x = input$FoldEquation)
    rv$FoldEquationText <- paste(rv$FoldEquationVector, collapse = "")
    rv$FoldEquationState <- tryCatch((function(X1, X2, X3, X4, X5, X6){eval(parse(text=rv$FoldEquationText))})(1,2,3,4,5,6), error = function(e) NA)
    if (is.null(need(rv$FoldEquationState, message = F))) {
      shinyjs::runjs("FoldHideFeedback()")
    } else {
      shinyjs::runjs(paste("FoldShowFeedback('", "Incorrect equation.", "')", sep = ""))
    }
    
    rv$FoldUsedDatasets <- sapply(c("X1", "X2", "X3", "X4", "X5", "X6"), function(x) {length(grep(x, rv$FoldEquationVector))})
    updateSliderInput(inputId = "FoldNumberOfStrains", value = max(2, max(which(rv$FoldUsedDatasets > 0))))
  })
  
  #### Fold Equation Settings ####
  observeEvent(input$FoldUsePreset, {
    if (input$FoldUsePreset == 1) {
      shinyjqui::updateOrderInput(session = session, inputId = "FoldEquation", item_class = "danger", items = c("Strain 1", "/", "Strain 2"))
    } else if (input$FoldUsePreset == 2) {
      shinyjqui::updateOrderInput(session = session, inputId = "FoldEquation", item_class = "danger", items = c("Strain 1", "-", "Strain 2"))
    } else if (input$FoldUsePreset == 3) {
      shinyjqui::updateOrderInput(session = session, inputId = "FoldEquation", item_class = "danger", items = c("(", "Strain 1", "/", "Strain 2", ")", "/", "(", "Strain 3", "/", "Strain 4", ")"))
    } else if (input$FoldUsePreset == 4) {
      shinyjqui::updateOrderInput(session = session, inputId = "FoldEquation", item_class = "danger", items = c("(", "Strain 1", "-", "Strain 3", ")", "/", "(", "Strain 2", "-", "Strain 3", ")"))
    } else if (input$FoldUsePreset == 5) {
      shinyjqui::updateOrderInput(session = session, inputId = "FoldEquation", item_class = "danger", items = c("(", "Strain 1", "+", "Strain 2", ")", "/", "Strain 3"))
    }
  })
  
  observeEvent(input$FoldClear, {
    shinyjqui::updateOrderInput(session = session, inputId = "FoldEquation", item_class = "danger", items = NULL)
    shinyjs::runjs("FoldRemoveElements()")
    updateSelectInput(inputId = "FoldUsePreset", selected = 0)
  })
  
  observe({
    rv$FoldNumberOfStrains <- input$FoldNumberOfStrains
    if (rv$FoldNumberOfStrains == 6){
      for (i in 3:6){
        eval(parse(text = paste("shinyjs::show('FoldColumn", i, "')", sep = "")))
      }
    } else if (rv$FoldNumberOfStrains == 2){
      for (i in 3:6){
        eval(parse(text = paste("shinyjs::hide('FoldColumn", i, "')", sep = "")))
      }
    } else {
      for (i in 3:rv$FoldNumberOfStrains){
        eval(parse(text = paste("shinyjs::show('FoldColumn", i, "')", sep = "")))
      }
      for (i in (rv$FoldNumberOfStrains+1):6){
        eval(parse(text = paste("shinyjs::hide('FoldColumn", i, "')", sep = "")))
      }
    }
  })
  
  #### Fold Calculate ####
  observeEvent(input$calculate5, {
    userDataStrain1 <- rv$CountsStrain1FoldUserInput()
    userDataStrain2 <- rv$CountsStrain2FoldUserInput()
    userDataStrain3 <- rv$CountsStrain3FoldUserInput()
    userDataStrain4 <- rv$CountsStrain4FoldUserInput()
    userDataStrain5 <- rv$CountsStrain5FoldUserInput()
    userDataStrain6 <- rv$CountsStrain6FoldUserInput()
    if (any(c(c(userDataStrain1$errors, userDataStrain2$errors, userDataStrain3$errors, userDataStrain4$errors, userDataStrain5$errors, userDataStrain6$errors)[1:rv$FoldNumberOfStrains], !is.null(need(rv$FoldEquationState, message = F))))) {
      output$errorBarFold <- renderText("There are problems with your data that need to be resolved. If you have trouble, check Help.")
      shinyjs::hide("advanced5")
      shinyjs::enable("calculate5")
      shinyjs::runjs("$('#calculate5').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    } else {
      foldExpr <- gsub(pattern = "X", replacement = "mu", x = rv$FoldEquationText)
      foldExpr <- gsub(pattern = "([1-6])", replacement = "\\[\\1\\]", x = foldExpr)
      foldExpr <- Ryacas::tex(Ryacas::yac_symbol(paste("f==", foldExpr, sep = "")))
      output$funFold <- renderUI({
        withMathJax(
          paste("$$", foldExpr, "$$")
        )
      })
      outputData <- unlist(calc.fold.int(userDataStrain1, userDataStrain2, userDataStrain3,
                                         userDataStrain4, userDataStrain5, userDataStrain6, rv$FoldEquationText))
      
      results5 <- data.frame(
        "Calculations" = outputData,
        row.names = c("f", "f<sup>95&percnt;&ndash;</sup>", "f<sup>95&percnt;&plus;</sup>",
                      "&mu;<sub>1</sub>", "&mu;<sub>2</sub>", "&mu;<sub>3</sub>",
                      "&mu;<sub>4</sub>", "&mu;<sub>5</sub>", "&mu;<sub>6</sub>")[1:length(outputData)],
        check.names = FALSE
      )
      output$tableFold <- reactable::renderReactable({
        reactable::reactable(results5,
                             rownames = TRUE,
                             sortable = FALSE,
                             pagination = FALSE,
                             defaultColDef = reactable::colDef(html = TRUE))
      })
      rv$FoldtoClipboard <- outputData
      output$errorBarFold <- renderText("")
      shinyjs::show("advanced5")
      shinyjs::enable("calculate5")
      shinyjs::runjs("$('#calculate5').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    }
    gc()
  })
  
  #### Fold Copy to Clipboard ####
  observeEvent(input$clipFold, {
    clipr::write_clip(paste(as.character(rv$FoldtoClipboard), collapse = "\n"), allow_non_interactive = TRUE)
  })
  
  #### Fold Erase ####
  observeEvent(input$erase5, {
    for (item in c("CountsStrain1Fold", "CountsStrain2Fold", "CountsStrain3Fold",
                   "CountsStrain4Fold", "CountsStrain5Fold", "CountsStrain6Fold")[1:input$FoldNumberOfStrains]) {
      callModule(countsPlatingUpdate, item, usePreset = 0)
    }
  })
  
  #### Fold Load Sample ####
  observeEvent(input$sample5, {
    for (item in c("CountsStrain1Fold", "CountsStrain2Fold", "CountsStrain3Fold",
                   "CountsStrain4Fold", "CountsStrain5Fold", "CountsStrain6Fold")[1:input$FoldNumberOfStrains]) {
      callModule(countsPlatingUpdate, item, usePreset = (1+((which(c("CountsStrain1Fold", "CountsStrain2Fold", "CountsStrain3Fold",
                                                                     "CountsStrain4Fold", "CountsStrain5Fold", "CountsStrain6Fold") == item)+1)%%2)))
    }
  })
  
  #### Power Validate Data ####
  rv$PowerStrain1UserInput <- callModule(PowerModule, id = "PowerStrain1", sampleSizeOrPower=reactive(input$sampleSizeOrPower))
  rv$PowerStrain2UserInput <- callModule(PowerModule, id = "PowerStrain2", sampleSizeOrPower=reactive(input$sampleSizeOrPower))
  observe({
    if (input$sampleSizeOrPower == 1) {
      rv$PowerValue <- numerise(input$PowerValue)
    }
    if (PowerValidator(rv$PowerValue) != "") {
      rv$PowerValueErrorsDetected <- TRUE
      textInputError(inputId = "PowerValue", text = paste(PowerValidator(rv$PowerValue)))
    } else {
      rv$PowerValueErrorsDetected <- FALSE
      shinyFeedback::hideFeedback(inputId = "PowerValue")
    }
  })
  
  #### Power Calculate ####
  observeEvent(input$calculate6, {
    PowerStrain1UserInput <- rv$PowerStrain1UserInput()
    PowerStrain2UserInput <- rv$PowerStrain2UserInput()
    if (PowerStrain1UserInput$errors || PowerStrain2UserInput$errors || rv$PowerValueErrorsDetected) {
      output$errorBarPower <- renderText("There are problems with your data that need to be resolved. If you have trouble, check Help.")
      shinyjs::hide("advanced6")
      shinyjs::enable("calculate6")
      shinyjs::runjs("$('#calculate6').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    } else {
      if (input$sampleSizeOrPower == 1) {
        outputData <- tryCatch(sample.size(power = rv$PowerValue, rate1 = PowerStrain1UserInput$MutationRate, rate2 = PowerStrain2UserInput$MutationRate,
                                           Nt1 = PowerStrain1UserInput$MeanCells, Nt2 = PowerStrain2UserInput$MeanCells,
                                           e1 = PowerStrain1UserInput$PlatingEfficiency, e2 = PowerStrain2UserInput$PlatingEfficiency,
                                           w1 = PowerStrain1UserInput$Fitness, w2 = PowerStrain2UserInput$Fitness,
                                           lag1 = PowerStrain1UserInput$Lag, lag2 = PowerStrain2UserInput$Lag,
                                           death1 = PowerStrain1UserInput$Death, death2 = PowerStrain2UserInput$Death,
                                           phi1 = PowerStrain1UserInput$Inoculum/PowerStrain1UserInput$MeanCells, phi2 = PowerStrain2UserInput$Inoculum/PowerStrain2UserInput$MeanCells,
                                           cv1 = PowerStrain1UserInput$CV, cv2 = PowerStrain2UserInput$CV,
                                           poisson1 = PowerStrain1UserInput$Residual, poisson2 = PowerStrain2UserInput$Residual),
                               error = function(err) NA)
        output$powerinfo <- renderText({
          if (is.na(outputData) == TRUE | is.nan(outputData) == TRUE) {
            Outp <- paste("impossible to calculate due to error")
          } else {
            Outp <- paste(signif(as.vector(outputData), 4))
          }
          paste("The minimum required sample size to achieve power ", rv$PowerValue, " is ", Outp, ".", sep = "")
        })
        results6 <- data.frame("Strain 1" = c(signif(PowerStrain1UserInput$MutationRate, 4), format(PowerStrain1UserInput$MeanCells, format = "e", digits = 3), PowerStrain1UserInput$PlatingEfficiency,
                                              signif(PowerStrain1UserInput$CV, 4), PowerStrain1UserInput$Fitness, PowerStrain1UserInput$Lag,
                                              PowerStrain1UserInput$Death, PowerStrain1UserInput$Residual, format(PowerStrain1UserInput$Inoculum/PowerStrain1UserInput$MeanCells, format = "e", digits = 3)),
                               "Strain 2" = c(signif(PowerStrain2UserInput$MutationRate, 4), format(PowerStrain2UserInput$MeanCells, format = "e", digits = 3), PowerStrain2UserInput$PlatingEfficiency,
                                              signif(PowerStrain2UserInput$CV, 4), PowerStrain2UserInput$Fitness, PowerStrain2UserInput$Lag,
                                              PowerStrain2UserInput$Death, PowerStrain2UserInput$Residual, format(PowerStrain2UserInput$Inoculum/PowerStrain2UserInput$MeanCells, format = "e", digits = 3)),
                               row.names = c("&mu;", "N<sub>t</sub>", "&epsilon;", "CV", "&rho;", "&lambda;", "d", "m<sub>p</sub>", "&phi;"),
                               check.names = FALSE)
      } else {
        outputData <- tryCatch(power.est(n1 = PowerStrain1UserInput$SampleSize, n2 = PowerStrain2UserInput$SampleSize,
                                         rate1 = PowerStrain1UserInput$MutationRate, rate2 = PowerStrain2UserInput$MutationRate,
                                         Nt1 = PowerStrain1UserInput$MeanCells, Nt2 = PowerStrain2UserInput$MeanCells,
                                         e1 = PowerStrain1UserInput$PlatingEfficiency, e2 = PowerStrain2UserInput$PlatingEfficiency,
                                         w1 = PowerStrain1UserInput$Fitness, w2 = PowerStrain2UserInput$Fitness,
                                         lag1 = PowerStrain1UserInput$Lag, lag2 = PowerStrain2UserInput$Lag,
                                         death1 = PowerStrain1UserInput$Death, death2 = PowerStrain2UserInput$Death,
                                         phi1 = PowerStrain1UserInput$Inoculum/PowerStrain1UserInput$MeanCells, phi2 = PowerStrain2UserInput$Inoculum/PowerStrain2UserInput$MeanCells,
                                         cv1 = PowerStrain1UserInput$CV, cv2 = PowerStrain2UserInput$CV,
                                         poisson1 = PowerStrain1UserInput$Residual, poisson2 = PowerStrain2UserInput$Residual),
                               error = function(err) NA)
        output$powerinfo <- renderText({
          if (is.na(outputData) == TRUE | is.nan(outputData) == TRUE) {
            Outp <- paste("impossible to calculate due to error")
          } else {
            Outp <- paste(signif(as.vector(outputData), 4))
          }
          paste("The calculated power of likelihood ratio test given specified sample sizes is ", Outp, ".", sep = "")
        })
        results6 <- data.frame("Strain 1" = c(PowerStrain1UserInput$SampleSize, signif(PowerStrain1UserInput$MutationRate, 4), format(PowerStrain1UserInput$MeanCells, format = "e", digits = 3), PowerStrain1UserInput$PlatingEfficiency,
                                              signif(PowerStrain1UserInput$CV, 4), PowerStrain1UserInput$Fitness, PowerStrain1UserInput$Lag,
                                              PowerStrain1UserInput$Death, PowerStrain1UserInput$Residual, format(PowerStrain1UserInput$Inoculum/PowerStrain1UserInput$MeanCells, format = "e", digits = 3)),
                               "Strain 2" = c(PowerStrain2UserInput$SampleSize, signif(PowerStrain2UserInput$MutationRate, 4), format(PowerStrain2UserInput$MeanCells, format = "e", digits = 3), PowerStrain2UserInput$PlatingEfficiency,
                                              signif(PowerStrain2UserInput$CV, 4), PowerStrain2UserInput$Fitness, PowerStrain2UserInput$Lag,
                                              PowerStrain2UserInput$Death, PowerStrain2UserInput$Residual, format(PowerStrain2UserInput$Inoculum/PowerStrain2UserInput$MeanCells, format = "e", digits = 3)),
                               row.names = c("n", "&mu;", "N<sub>t</sub>", "&epsilon;", "CV", "&rho;", "&lambda;", "d", "m<sub>p</sub>", "&phi;"),
                               check.names = FALSE)
      }
      output$tablePower <- reactable::renderReactable({
        reactable::reactable(
          results6,
          rownames = TRUE,
          sortable = FALSE,
          wrap = FALSE,
          pagination = FALSE,
          defaultColDef = reactable::colDef(html = TRUE))
      })
      output$errorBarPower <- renderText("")
      shinyjs::show("advanced6")
      shinyjs::enable("calculate6")
      shinyjs::runjs("$('#calculate6').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    }
    gc()
  })
  
  #### Power Erase ####
  observeEvent(input$erase6, {
    callModule(PowerModuleUpdate, "PowerStrain1", usePreset = 0)
    callModule(PowerModuleUpdate, "PowerStrain2", usePreset = 0)
  })
  
  ####  Power Load Sample ####
  observeEvent(input$sample6, {
    callModule(PowerModuleUpdate, "PowerStrain1", usePreset = 1)
    callModule(PowerModuleUpdate, "PowerStrain2", usePreset = 2)
  })
  
  #### BatchCalc file upload ####
  observeEvent(input$userData, {
    output$batchInfo <- renderText("")
    shinyjs::hide("onDataChecked")
    shinyjs::hide("ifDataCorrect")
    shinyjs::hide("onOutputGenerated")
    shinyjs::runjs("$('a[data-value=\"Experiment parameters\"]').click()")
    
    userInput <- tryCatch(readxl::read_excel(input$userData$datapath, sheet = 1, col_types = "text"), error = function(err) {NULL})
    if (is.null(userInput)) {
      output$errorBarFile <- renderText("<p class=\"text-danger\">Invalid file. Please make sure you uploaded the correct file.</p>")
      shinyjs::hide("onDataUploaded")
      return(NULL)
    }
    sheet_number <- length(readxl::excel_sheets(input$userData$datapath))
    if (sheet_number != 1) {
      text1 <- "This file has more than 1 sheets. Only first sheet was loaded.\n"
    } else {text1 <- ""}
    userInput <- userInput[!apply(is.na(userInput), 1, all),!apply(is.na(userInput), 2, all)]
    text2 <- ""
    if (is.na(any(as.character(userInput[[1]])=="CountsSelective", na.rm = TRUE))) {
      output$errorBarFile <- renderText("<p class=\"text-danger\">Invalid file. Please make sure you uploaded the correct file.</p>")
      shinyjs::hide("onDataUploaded")
      return(NULL)
    } else if (!any(as.character(userInput[[1]])=="CountsSelective", na.rm = TRUE)) {
      text2 <- "CountsSelective was not found in the first column, therefore mlemur assumed that you provided only colony counts on selective medium."
      userInput <- as.data.frame(cbind(rep("CountsSelective", nrow(userInput)), as.matrix(userInput)))
    }
    for (i in 2:nrow(userInput)) {
      if (is.na(userInput[[1]][i-1])) {userInput[[1]][i-1] <- "NoName"}
      if (is.na(userInput[[1]][i])) {userInput[[1]][i] <- userInput[[1]][i-1]}
    }
    rv$DataSelective <- subset(userInput, subset = userInput[[1]]=='CountsSelective')
    rv$DataNonselective <- subset(userInput, subset = userInput[[1]]=='CountsNonselective')
    rv$DataPlating <- subset(userInput, subset = userInput[[1]]!='CountsSelective' & userInput[[1]]!='CountsNonselective')
    output$errorBarFile <- renderText(paste(text1, text2, sep = ""))
    shinyjs::show("onDataUploaded")
  })
  
  #### BatchCalc data checker ####
  observeEvent(input$dataCheck, {
    
    rv$BatchFirstRun <- TRUE
    
    validation <- comboValidator(rv$DataPlating, rv$DataSelective, rv$DataNonselective)
    
    for (item in c("Fitness", "Lag", "Residual", "Death", "Inoculum", "CV")) {
      shinyWidgets::updateAwesomeCheckbox(session = session,
                                          inputId = paste("choices", "AllStrains", item, sep = "_"),
                                          value = FALSE)
    }
    
    if (is.null(validation)==FALSE) {
      
      DataPlating <- validation[[1]]
      DataSelective <- validation[[2]]
      DataNonselective <- validation[[3]]
      
      rv$DataPlatingForCalc <- validation[[4]]
      rv$DataSelectiveForCalc <- validation[[5]]
      rv$DataNonselectiveForCalc <- validation[[6]]
      rv$InitialSettings <- t(apply(as.matrix(rv$DataPlatingForCalc)[c(2,3,4,5,6,13),,drop = FALSE], MARGIN=c(1,2), FUN = function(x) {if (isTruthy(x)) {x>0} else FALSE}))
      warningList <- validation[[7]]
      errorList <- validation[[8]]
      
      adjustTabNames(validation[[9]], validation[[10]], validation[[11]])
      
      rv$CVFromCounts <- validation[[12]]
      rv$PvaluePossible <- validation[[13]]
      
      if (warningList != "" & errorList != "") {
        output$batchInfo <- renderText({paste("<p class=\"text-danger\">There are problems with your data that need to be resolved. Look for <i class='fa fa-exclamation-triangle' style = 'color:#b94a48'></i> red triangles for more details:</p>", errorList, "\n",
                                              "<p>Additionally, mlemur has generated the following non-error informations and warnings:</p>", warningList, 
                                              "<p class=\"text-success\">Valid datasets, if any, were added to the Load dataset menu in tabs: Rate, <i>P</i> value, Fold.</p>", sep = "")})
        shinyjs::hide("ifDataCorrect")
      } else if (errorList != "") {
        output$batchInfo <- renderText({paste("<p class=\"text-danger\">There are problems with your data that need to be resolved. Look for red triangles for more details:</p>", errorList,  
                                              "<p class=\"text-success\">Valid datasets, if any, were added to the Load dataset menu in tabs: Rate, <i>P</i> value, Fold.</p>", sep = "")})
        shinyjs::hide("ifDataCorrect")
      } else {
        if (warningList != "") {
          output$batchInfo <- renderText({paste("<p class=\"text-success\">Your data appears to be OK.</p>", "\n",
                                                "<p>Additionally, mlemur has generated the following non-error informations and warnings:</p>", warningList,  
                                                "<p class=\"text-success\">Valid datasets, if any, were added to the Load dataset menu in tabs: Rate, <i>P</i> value, Fold.</p>", sep = "")})
        } else {
          output$batchInfo <- renderText({paste("<p class=\"text-success\">Your data appears to be OK.</p>",
                                                "<p class=\"text-success\">Valid datasets, if any, were added to the Load dataset menu in tabs: Rate, <i>P</i> value, Fold.</p>", sep = "")})
        }
        shinyjs::show("ifDataCorrect")
        
        choiceslist <- list()
        
        for (i in 1:nrow(rv$InitialSettings)) {
          rowchoices <- list()
          for (item in colnames(rv$InitialSettings)) {
            rowchoices <- c(rowchoices, tagList(div(style = "display:inline-block;vertical-align:top;width:60px;padding-left:20px;",
                                                    shinyWidgets::awesomeCheckbox(inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), item, sep = "_"),
                                                                                  label = "",
                                                                                  status = "primary",
                                                                                  value = rv$InitialSettings[i,item]))))
          }
          choiceslist <- c(choiceslist, tagList(
            div(style = "display:block;vertical-align:top;width:100%;",
                div(style = "display:inline-block;vertical-align:top;width:calc(100% - 390px);padding-bottom:15px;", HTML(paste(rownames(rv$InitialSettings)[i]))),
                tagList(rowchoices))))
        }
        
        output$BatchOptionsPval <- renderUI({
          tagList(
            shinyjs::disabled(div(class = "checkboxdiv",
                                  shinyWidgets::awesomeCheckbox("BatchCalculateRates",
                                                                HTML("Mutation rates are always calculated."),
                                                                status = "primary",
                                                                value = TRUE))),
            if (rv$PvaluePossible)
            {tagList(div(class = "checkboxdiv",
                         shinyWidgets::awesomeCheckbox("BatchCalculatePvalues",
                                                       HTML("Calculate <em>P</em> values"),
                                                       status = "primary")),
                     conditionalPanel(condition = "input.BatchCalculatePvalues",
                                      div(class = "checkboxdiv",
                                          shinyWidgets::awesomeCheckbox("batchApplyCorrection",
                                                                        HTML("Apply a <em>P</em> value correction"),
                                                                        status = "primary")),
                                      conditionalPanel(condition = "input.batchApplyCorrection",
                                                       div(class = "checkboxdiv",
                                                           style = "padding-left:27px; padding-bottom:5px;",
                                                           shinyWidgets::awesomeRadio("batchCorrectionMethod",
                                                                                      label = HTML(paste("Select method of correction:", infoTooltip("A more detailed description of available correction methods can be found in Help."))),
                                                                                      choices = c("Bonferroni" = "bonferroni", "Bonferroni-Holm" = "holm", "Benjamini-Hochberg" = "BH"),
                                                                                      selected = "bonferroni",
                                                                                      width = "100%") )
                                      )
                     )
            )}
          )
        })
        output$BatchOptionsChoices <- renderUI({
          tagList(choiceslist)
        })
      }
      
      output$platingTable <- reactable::renderReactable({defaultTable(DataPlating, "Input", TRUE)})
      output$selectiveTable <- reactable::renderReactable({defaultTable(DataSelective, "No")})
      output$nonselectiveTable <- reactable::renderReactable({defaultTable(DataNonselective, "No")})
      
      for (item in colnames(rv$InitialSettings)) {
        if (all(rv$InitialSettings[,item])) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", "AllStrains", item, sep = "_"),
                                              value = TRUE)
        }
      }
      
      rv$PlatingList <- rv$DataPlatingForCalc
      rv$PlatingList <- lapply(rv$PlatingList, as.list)
      for (item in colnames(rv$DataPlatingForCalc)) {
        eval(parse(text = paste("names(rv$PlatingList$`", item, "`) <- rownames(rv$DataPlatingForCalc)", sep = "")))
        eval(parse(text = paste("rv$PlatingList$`", item, "`$CountsSelective <- rv$DataSelectiveForCalc$`", item, "`", sep = "")))
        eval(parse(text = paste("rv$PlatingList$`", item, "`$CountsNonselective <- rv$DataNonselectiveForCalc$`", item, "`", sep = "")))
        eval(parse(text = paste("rv$PlatingList$`", item, "`$model <- TRUE", sep = "")))
        eval(parse(text = paste("rv$PlatingList$`", item, "`$setSel <- ifelse(is.na(rv$PlatingList$`", item, "`$PlatingEfficiency), 1, 2)", sep = "")))
        eval(parse(text = paste("rv$PlatingList$`", item, "`$setNsel <- ifelse(is.na(rv$PlatingList$`", item, "`$MeanCells), 1, 2)", sep = "")))
        eval(parse(text = paste("if (rv$InitialSettings['", item, "','Lag']) {rv$PlatingList$`", item, "`$useLagFitness <- 1} else if (rv$InitialSettings['", item, "','Fitness']) {rv$PlatingList$`", item, "`$useLagFitness <- 2} else {rv$PlatingList$`", item, "`$useLagFitness <- 0}", sep = "")))
        eval(parse(text = paste("rv$PlatingList$`", item, "`$setCV <- 1", sep = "")))
      }
      rv$BatchDatasets <- rv$PlatingList
      names(rv$BatchDatasets) <- as.vector(sapply(names(rv$BatchDatasets), function(x) paste("Batch-", x, sep="")))
      
      for (item in colnames(rv$DataPlatingForCalc)) {
        eval(parse(text = paste("rv$PlatingList$`", item, "`$setCV <- !rv$CVFromCounts[which(colnames(rv$DataPlatingForCalc)=='", item, "')]", sep = "")))
      }
      
      shinyjs::show("onDataChecked")
      
    } else {
      shinyjs::hide("ifDataCorrect")
      output$batchInfo <- renderText("<p class=\"text-danger\">There are no data in this file.</p>")
    }
    
    shinyFeedback::resetLoadingButton("dataCheck")
    
  })
  
  #### BatchCalc toggle observers ####
  observe({
    lapply(X = rownames(rv$InitialSettings), FUN = function(x) {
      observeEvent(input[[paste("choices", gsub(pattern = " ", replacement = "_", x), "Fitness", sep = "_")]], {
        if (input[[paste("choices", gsub(pattern = " ", replacement = "_", x), "Fitness", sep = "_")]]) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", x), "Lag", sep = "_"),
                                              value = FALSE)
        }
      })
    })
    
    lapply(X = rownames(rv$InitialSettings), FUN = function(x) {
      observeEvent(input[[paste("choices", gsub(pattern = " ", replacement = "_", x), "Lag", sep = "_")]], {
        if (input[[paste("choices", gsub(pattern = " ", replacement = "_", x), "Lag", sep = "_")]]) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", x), "Fitness", sep = "_"),
                                              value = FALSE)
        }
      })
    })
  })
  
  observeEvent(input$choices_AllStrains_Fitness, {
    if (rv$BatchFirstRun) {
      if (input$choices_AllStrains_Fitness) {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Fitness", sep = "_"),
                                              value = TRUE)
        }
        shinyWidgets::updateAwesomeCheckbox(session = session,
                                            inputId = "choices_AllStrains_Lag",
                                            value = FALSE)
      } else {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Fitness", sep = "_"),
                                              value = FALSE)
        }
      }
    }
  })
  
  observeEvent(input$choices_AllStrains_Lag, {
    if (rv$BatchFirstRun) {
      if (input$choices_AllStrains_Lag) {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Lag", sep = "_"),
                                              value = TRUE)
        }
        shinyWidgets::updateAwesomeCheckbox(session = session,
                                            inputId = "choices_AllStrains_Fitness",
                                            value = FALSE)
      } else {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Lag", sep = "_"),
                                              value = FALSE)
        }
      }
    }
  })
  
  observeEvent(input$choices_AllStrains_Residual, {
    if (rv$BatchFirstRun) {
      if (input$choices_AllStrains_Residual) {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Residual", sep = "_"),
                                              value = TRUE)
        }
      } else {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Residual", sep = "_"),
                                              value = FALSE)
        }
      }
    }
  })
  
  observeEvent(input$choices_AllStrains_Death, {
    if (rv$BatchFirstRun) {
      if (input$choices_AllStrains_Death) {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Death", sep = "_"),
                                              value = TRUE)
        }
      } else {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Death", sep = "_"),
                                              value = FALSE)
        }
      }
    }
  })
  
  observeEvent(input$choices_AllStrains_Inoculum, {
    if (rv$BatchFirstRun) {
      if (input$choices_AllStrains_Inoculum) {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Inoculum", sep = "_"),
                                              value = TRUE)
        }
      } else {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "Inoculum", sep = "_"),
                                              value = FALSE)
        }
      }
    }
  })
  
  observeEvent(input$choices_AllStrains_CV, {
    if (rv$BatchFirstRun) {
      if (input$choices_AllStrains_CV) {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "CV", sep = "_"),
                                              value = TRUE)
        }
      } else {
        for (i in 1:nrow(rv$InitialSettings)) {
          shinyWidgets::updateAwesomeCheckbox(session = session,
                                              inputId = paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), "CV", sep = "_"),
                                              value = FALSE)
        }
      }
    }
  })
  
  #### BatchCalc calculations ####
  observeEvent(input$calculate4, {

    rv$FinalSettings <- rv$InitialSettings
    for (i in 1:nrow(rv$InitialSettings)) {
      for (item in colnames(rv$InitialSettings)) {
        rv$FinalSettings[i,item] <- input[[paste("choices", gsub(pattern = " ", replacement = "_", rownames(rv$InitialSettings)[i]), item, sep = "_")]]
      }
    }
    rv$PlatingListCopy <- rv$PlatingList
    rv$PlatingListCopy <- removeParameters(rv$PlatingListCopy, rv$FinalSettings)
    
    DataRate <- lapply(rv$PlatingListCopy, calc.rate.int)
    filetosave <- rv$PlatingListCopy
    save(filetosave, file = "latest.RData")
    NoOfMutations <- as.numeric(lapply(DataRate, function(x) {x[15]}))
    DataRate <- as.data.frame(lapply(DataRate, function(x) {x[-15]}), check.names = FALSE)
    
    rownames(DataRate) <- c("&epsilon;", "N<sub>t</sub>", "CV", "&rho;", "&lambda;", "d", "m<sub>p</sub>", "&phi;",
                            "m", "m<sup>95&percnt;&ndash;</sup>", "m<sup>95&percnt;&plus;</sup>", "&mu;", "&mu;<sup>95&percnt;&ndash;</sup>", "&mu;<sup>95&percnt;&plus;</sup>")
    output$batchRateResults <- reactable::renderReactable({defaultTable(DataRate, "Parameters")})
    rv$RateOutput <- data.frame("Parameters" = c("eff", "Nt", "CV", "Fitness", "Lag", "Death", "Residual", "Inoculum", "m", "m_low", "m_up", "mu", "mu_low", "mu_up"), DataRate, check.names = FALSE)

    if (input$BatchCalculatePvalues) {
      
      table.length <- length(DataRate)
      xarg <- rep(1:(table.length-1), (table.length-1):1)
      yarg <- unlist(sapply(2:table.length, function(x) x:table.length))
      
      proxy <- mapply(function(x,y) {calc.pval.int(rv$PlatingListCopy[[x]], rv$PlatingListCopy[[y]], NoOfMutations[[x]], NoOfMutations[[y]])}, x = xarg, y = yarg)

      DataPvalue <- matrix(NA,nrow=table.length,ncol=table.length,dimnames = list(names(DataRate),names(DataRate)))
      DataPvalue[lower.tri(DataPvalue)] <- proxy
      
      if (input$batchApplyCorrection) {
        selectedMethod <- input$batchCorrectionMethod
        DataPvalue[] <- p.adjust(DataPvalue, method=selectedMethod)
      }
      
      DataPvalue[upper.tri(DataPvalue)] <- t(DataPvalue)[upper.tri(DataPvalue)] # do it after correction to avoid doubling; must use transpose to correctly fill second triangle
      
      output$BatchPvalueResults <- reactable::renderReactable({
        DataPvalueR <- matrix(NA,nrow=table.length,ncol=table.length,dimnames = list(names(DataRate),names(DataRate)))
        m <- mapply(function(x,y) {colourPvalue(DataPvalue[x,y])}, x = xarg, y = yarg)
        DataPvalueR[lower.tri(DataPvalueR)] <- m
        DataPvalueR[upper.tri(DataPvalueR)] <- t(DataPvalueR)[upper.tri(DataPvalueR)]
        rownames(DataPvalueR) <- sapply(rownames(DataPvalue), function(x) paste("<b>", x, "</b>", sep = ""))
        defaultTable(DataPvalueR, "Strain")
      })
      
      rv$PvalueOutput <- data.frame("Strain" = colnames(DataPvalue), as.data.frame(DataPvalue, check.names = FALSE))
      shinyjs::show("BatchPvalue")
    } else {
      shinyjs::hide("BatchPvalue")
      shinyjs::hide("BatchEffSize")
      rv$PvalueOutput <- NULL
    }
    
    shinyjs::show("onOutputGenerated")
    
    shinyjs::enable("calculate4")
    shinyjs::runjs("$('#calculate4').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    gc()
  })
  
  #### BatchCalc file export ####
  output$dlBatch <- downloadHandler(
    filename = function() {
      paste("Report-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      FileList <- list("Rates" = rv$RateOutput)
      if (is.null(rv$PvalueOutput)==FALSE) {
        FileList$`P values` <- rv$PvalueOutput
      }
      writexl::write_xlsx(FileList, path = file)
    }
  )
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
}