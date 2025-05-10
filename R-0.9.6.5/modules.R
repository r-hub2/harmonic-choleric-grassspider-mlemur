settingsPlatingUI <- function(id, defaultSettings) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 12,
        shinyWidgets::awesomeRadio(
          inputId = ns("setPerPlate"),
          label = HTML(paste("Method of estimating mutation rate:", infoTooltip("If the first option is selected, first, the average number of mutations will be calculated using selective medium data, and then this value will be divided by an estimate of culture size to obtain an estimate of mutation rate. If the second option is selected, data for selective and non-selective medium will be treated as paired data for each culture and mutation rate will be estimated directly. In this case, the same number of counts on selective and non-selective medium must be provided."))),
          choices = c("Use an estimate of average culture size from a limited number of non-selective plates" = 1, "Use pairs of counts on both media for each culture" = 2),
          selected = 1,
          inline = TRUE
        )
      )
    )
  )
}

settingsPlating <- function(input, output, session) {
  rv <- reactiveValues()
  
  observe({
    rv$value <- list("setPerPlate" = input$setPerPlate)
  })
  return(reactive(rv$value))
}

countsPlatingUI <- function(id, currentTab, stack_cols=FALSE, usePreset, defaultSettings) {
  ns <- NS(id)
  number <- ifelse(stack_cols, 12, 6)
  if (usePreset==1) {useValues <- c("500", "100", "1", "1\n2\n3", "100", "1000000", "100\n200\n300", "0.9", "0.2", "1e+9", "0", "0", "0", "0", "0", "0")}
  else if (usePreset) {useValues <- c("500", "150", "1", "7\n8\n9", "100", "1000000", "90\n100\n110", "1", "0.3", "5e+8", "0", "0", "0", "0", "0", "0")}
  tagList(
    shinyFeedback::useShinyFeedback(),
    fluidRow(
      column(
        width = 12,
        shinyWidgets::awesomeRadio(
          inputId = ns("setSelective"),
          label = HTML(paste("Plating efficiency:", infoTooltip("Choose whether you want to calculate the fraction of culture plated on selective medium (plating efficiency) using appropriate volumes, or already have a pre-calculated value."))),
          choices = c("Calculate using provided volumes and dilutions" = 1, "Use a pre-calculated value" = 2),
          selected = defaultSettings[2],
          inline = TRUE
        ),
        shinyWidgets::awesomeRadio(
          inputId = ns("setNonselective"),
          label = HTML(paste("Average number of cells in the culture:", infoTooltip("Choose whether you want to calculate average number of cells in culture (culture size) using colony counts and appropriate volumes, or already have a pre-calculated value."))),
          choices = c("Calculate using the colony counts and volumes" = 1, "Use a pre-calculated value" = 2),
          selected = defaultSettings[3],
          inline = TRUE
        ),
        shinyWidgets::awesomeCheckbox(
          inputId = ns("setModel"),
          label = HTML(paste("Use advanced experiment parameters", infoTooltip("A more detailed description of available mutant models can be found in Help."))),
          value = FALSE
        )
      )
    ),
    conditionalPanel(
      condition="input.setModel",
      ns = NS(id),
      fluidRow(
        column(
          width = 12,
          shinyWidgets::awesomeRadio(
            inputId = ns("useLagFitness"),
            label = HTML(paste("Specify phenotypic lag or fitness:", infoTooltip("Choose whether you want to input phenotypic lag or mutant fitness."))),
            choices = c("Don't use" = 0, "Specify phenotypic lag" = 1, "Specify mutant fitness" = 2),
            selected = 0,
            inline = TRUE
          )
        ),
        column(
          width = number,
          conditionalPanel(
            condition="input.useLagFitness==1",
            ns = NS(id),
            textInput(
              inputId = ns("Lag"),
              label = HTML(paste("Phenotypic lag (generations):", infoTooltip("Please provide a non-negative number (&ge;&nbsp;0)."))),
              value = useValues[12],
              width = "100%"
            )
          ),
          conditionalPanel(
            condition="input.useLagFitness==2",
            ns = NS(id),
            textInput(
              inputId = ns("Fitness"),
              label = HTML(paste("Mutant relative fitness:", infoTooltip("Please provide a positive number. Fitness &lt;&nbsp;1 means that mutants grow <u>slower</u> than non-mutants, while &gt;&nbsp;1 means that mutants grow faster. Fitness equal to 1 is equivalent to L-C model.<br><b>Note: underflow may be encountered when when both plating efficiency and fitness are much smaller than 0.1.</b>"))),
              value = useValues[8],
              width = "100%"
            )
          )
        ),
        column(
          width = 12,
          textInput(
            inputId = ns("Death"),
            label = HTML(paste("Death rate:", infoTooltip("Please provide a non-negative number that is a fraction of the birth rate (&ge;&nbsp;0 but &lt;&nbsp;1)."))),
            value = useValues[14],
            width = "100%"
          )
        ),
        column(
          width = number,
          textInput(
            inputId = ns("Residual"),
            label = div(id = paste(id, "Residual", "labeltext", sep = "-"), HTML(paste("Residual mutations:", infoTooltip("Please provide a non-negative number. Floating-point numbers will be rounded down to nearest integer. Residual mutations apply only to the portion of culture that was plated.")))),
            value = useValues[13],
            width = "100%"
          )
        ),
        column(
          width = number,
          textInput(
            inputId = ns("Inoculum"),
            label = HTML(paste("Size of the inoculum (number of cells):", infoTooltip("Please provide a non-negative number. Floating-point numbers will be rounded down to nearest integer. Inoculum size cannot exceed the total number of cells in grown culture."))),
            value = useValues[16],
            width = "100%"
          )
        ),
        column(
          width = 12,
          conditionalPanel(
            condition = "input.setPerPlate==1",
            ns = NS(currentTab),
            div(
              style="display:inline-block; vertical-align:top;",
              shinyWidgets::awesomeRadio(
                inputId = ns("setCV"),
                label = HTML(paste("Coefficient of variation:", infoTooltip("Use if you want to account for the variation of the final number of cells in each culture. If using custom value, please provide a positive number. Scientific notation can be used, e.g. 2e-2. Values smaller than 1e-6 will be set to 0, and values bigger than 10 will be set to 10. Leave 0 if you do not want to use the coefficient of variation. When only one colony count or only the final number of cells is provided, it will default to 0 unless supplied with another value."))),
                choices = c("Calculate using the colony counts" = 0, "Use custom value:" = 1),
                selected = 1,
                inline = TRUE
              )
            ),
            div(
              style = "display:inline-block; vertical-align:bottom;",
              textInput(
                inputId = ns("CV"),
                label = NULL,
                value = useValues[11],
                width = "100%",
              )
            )
          )
        )
      ),
      hr()
    ),
    conditionalPanel(
      condition="input.setSelective==1 || input.setNonselective==1",
      ns = NS(id),
      textInput(
        inputId = ns("VolumeTotal"),
        label = HTML(paste("Total culture volume (&mu;l):", infoTooltip("Please provide a positive number."))),
        value = useValues[1],
        width = "100%"
      )
    ),
    fluidRow(
      column(
        width = number,
        h4("Selective medium:"),
        conditionalPanel(
          condition="input.setSelective==1",
          ns = NS(id),
          textInput(
            inputId = ns("VolumeSelective"),
            label = HTML(paste("Volume plated on selective medium (&mu;l):", infoTooltip("Please provide a positive number. Please note that this value, divided by the dilution factor, cannot exceed total culture volume."))),
            value = useValues[2],
            width = "100%"
          ),
          textInput(
            inputId = ns("DilutionSelective"),
            label = HTML(paste("Dilution factor:", infoTooltip("Please provide a positive number &gt;&nbsp;1 if culture was diluted prior plating. If culture was concentrated, use a number &lt;&nbsp;1. If culture was neither diluted nor concentrated, type 1."))),"",
            value = useValues[3],
            width = "100%"
          )
        ),
        conditionalPanel(
          condition="input.setSelective==2",
          ns = NS(id),
          textInput(
            inputId = ns("PlatingEfficiency"),
            label = HTML(paste("Plating efficiency:", infoTooltip("Please provide a positive number so that 0&nbsp;&lt;&nbsp;&epsilon;&nbsp;&le;&nbsp;1."))),
            value = useValues[9],
            width = "100%"
          ),
        ),
        textAreaInput(
          inputId = ns("CountsSelective"),
          label = HTML(paste("Colony counts:", infoTooltip("Zeros and positive numbers are accepted. At least two numbers must be provided, and at least one must be bigger than 0. Paste one number under the other. Do not separate the numbers with a comma, semicolon, or any other character."))),
          value = useValues[4],
          rows = 10,
          width = "100%",
          resize = "vertical"
        ),
        br()
      ),
      column(
        width = number,
        h4("Non-selective medium:"),
        conditionalPanel(
          condition = "input.setNonselective==1",
          ns = NS(id),
          textInput(
            inputId = ns("VolumeNonselective"),
            label = HTML(paste("Volume plated on non-selective medium (&mu;l):", infoTooltip("Please provide a positive number. Please note that this value, divided by the dilution factor, cannot exceed total culture volume."))),
            value = useValues[5],
            width = "100%"
          ),
          textInput(
            inputId = ns("DilutionNonselective"),
            label = HTML(paste("Dilution factor:", infoTooltip("Please provide a positive number &gt;&nbsp;1 if culture was diluted prior plating. If culture was concentrated, use a number &lt;&nbsp;1. If culture was neither diluted nor concentrated, type 1."))),"",
            value = useValues[6],
            width = "100%"
          ),
          textAreaInput(
            inputId = ns("CountsNonselective"),
            label = HTML(paste("Colony counts:", infoTooltip("Zeros and positive numbers are accepted. At least one number must be provided if you have chosen not to account for the variation of the culture size (L-C or M-K models), and at least two if chosen otherwise (B<sup>0</sup> model). At least one colony count must be bigger than 0. Paste one number under the other. Do not separate the numbers with a comma, semicolon, or any other character."))),
            value = useValues[7],
            rows = 10,
            width = "100%",
            resize = "vertical"
          )
        ),
        conditionalPanel(
          condition = "input.setNonselective==2",
          ns = NS(id),
          textInput(
            inputId = ns("MeanCells"),
            label = HTML(paste("Average number of cells in culture:", infoTooltip("Please provide a positive number. Scientific notation can be used, e.g. 4e9."))),
            value = useValues[10],
            width = "100%"
          )
        )
      ),
      hr(),
      column(
        width = 12,
        selectInput(inputId = ns("LoadDataset"),
                    label = HTML(paste("Load dataset:", infoTooltip("Here you can copy data from other form in the app, or use one of datasets described previously in the literature."))),
                    choices = c("-", "Rate-Strain", "Pvalue-Strain1", "Pvalue-Strain2",
                                "Fold-Strain1", "Fold-Strain2", "Fold-Strain3",
                                "Fold-Strain4", "Fold-Strain5", "Fold-Strain6",
                                "Barna-glucose", "Barna-maltose", "Crane-1", "Crane-2", "Demerec", "Foster-1994",
                                "Vaisman-2021-1", "Vaisman-2021-2"),
                    selected = "-")
      )
    )
  )
}

countsPlating <- function(input, output, session, userSettings, stack_cols) {
  ReactValue <- reactiveValues()
  
  observeEvent(input$useLagFitness,{
    if (!stack_cols) {
      if (input$useLagFitness!=0) {
        shinyjs::runjs("document.getElementById('CountsRate-Death-label').parentElement.parentElement.setAttribute('class','col-sm-6')")
      } else {
        shinyjs::runjs("document.getElementById('CountsRate-Death-label').parentElement.parentElement.setAttribute('class','col-sm-12')")
      }
    }
  })
  
  observe({
    
    if (input$setCV==0) {
      shinyjs::disable(id = "CV")
    } else {
      shinyjs::enable(id = "CV")
    }
    
    ReactValue$setNonselective <- input$setNonselective
    ReactValue$setSelective <- input$setSelective
    ReactValue$setModel <- input$setModel
    ReactValue$setPerPlate <- userSettings()$setPerPlate
    ReactValue$setCV <- input$setCV
    ReactValue$useLagFitness <- input$useLagFitness
    
    if (ReactValue$setPerPlate==2) {
      if (input$setNonselective==2) {
        shinyFeedback::showToast(type = "warning", message = "'Average number of cells in the culture' must be 'Calculate using the colony counts and volumes' when 'Use pairs of counts on both media for each culture' is selected.", .options = list(positionClass = "toast-bottom-right", progressBar = FALSE, timeOut = 3000))
      }
      shinyWidgets::updateAwesomeRadio(session = session, inputId = "setNonselective", selected = 1)
      shinyjs::disable(id = "CV")
    }
    
    moduleID <- session$ns("proxy")
    moduleID <- substr(moduleID, 1, nchar(moduleID)-6)
    if (ReactValue$setPerPlate == 2) {
      shinyjs::runjs(paste('document.getElementById("', paste(moduleID, 'Residual', 'labeltext', sep = '-'), '").innerHTML = "', paste('Rate of residual mutations:', infoTooltip('Please provide a non-negative decimal number smaller than 1.')), '";', sep = ''))
    } else {
      shinyjs::runjs(paste('document.getElementById("', paste(moduleID, 'Residual', 'labeltext', sep = '-'), '").innerHTML = "',  paste('Residual mutations:', infoTooltip('Please provide a non-negative number. Floating-point numbers will be rounded down to nearest integer. Residual mutations apply only to the portion of culture that was plated.')), '";', sep = ''))
    }
    
    if (ReactValue$setNonselective==1 | ReactValue$setSelective==1) {
      ReactValue$VolumeTotal <- numerise(input$VolumeTotal)
    } else {
      ReactValue$VolumeTotal <- NA
    }
    if (ReactValue$setSelective==1) {
      ReactValue$VolumeSelective <- numerise(input$VolumeSelective)
      ReactValue$DilutionSelective <- numerise(input$DilutionSelective)
      ReactValue$PlatingEfficiency <- NA
    } else if (ReactValue$setSelective==2) {
      ReactValue$VolumeSelective <- NA
      ReactValue$DilutionSelective <- NA
      ReactValue$PlatingEfficiency <- numerise(input$PlatingEfficiency)
    }
    ReactValue$CountsSelective <- numerise(input$CountsSelective)
    if (ReactValue$setModel) {
      if (ReactValue$useLagFitness == 1) {
        ReactValue$Lag <- numerise(input$Lag)
      } else {
        ReactValue$Lag <- NA
      }
      if (ReactValue$useLagFitness == 2) {
        ReactValue$Fitness <- numerise(input$Fitness)
      } else {
        ReactValue$Fitness <- NA
      }
      ReactValue$Death <- numerise(input$Death)
      ReactValue$Residual <- numerise(input$Residual)
      ReactValue$Inoculum <- numerise(input$Inoculum)
      if (ReactValue$setCV==1) {
        ReactValue$CV <- numerise(input$CV)
      }
    } else {
      ReactValue$Fitness <- NA
      ReactValue$Lag <- NA
      ReactValue$Residual <- NA
      ReactValue$Death <- NA
      ReactValue$Inoculum <- NA
      ReactValue$CV <- NA
      shinyFeedback::hideFeedback(inputId = "CV")
    }
    if (ReactValue$setNonselective==1) {
      ReactValue$VolumeNonselective <- numerise(input$VolumeNonselective)
      ReactValue$DilutionNonselective <- numerise(input$DilutionNonselective)
      ReactValue$CountsNonselective <- numerise(input$CountsNonselective)
      ReactValue$MeanCells <- NA
    } else {
      ReactValue$VolumeNonselective <- NA
      ReactValue$DilutionNonselective <- NA
      ReactValue$CountsNonselective <- NA
      ReactValue$MeanCells <- numerise(input$MeanCells)
    }
    
    ReactValue$errorsDetected <- FALSE
    
    if (ReactValue$setNonselective == 1 || ReactValue$setSelective == 1) {
      if (ValueValidator(ReactValue$VolumeTotal) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "VolumeTotal", text = paste(ValueValidator(ReactValue$VolumeTotal)))
      } else {
        shinyFeedback::hideFeedback(inputId = "VolumeTotal")
      }
    }
    
    if (ReactValue$setModel) {
      
      if (ReactValue$useLagFitness == 1) {
        if (NonNegValueValidator(ReactValue$Lag) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "Lag", text = paste(NonNegValueValidator(ReactValue$Lag)))
        } else {
          shinyFeedback::hideFeedback(inputId = "Lag")
        }
      }
      
      if (ReactValue$useLagFitness == 2) {
        if (ValueValidator(ReactValue$Fitness) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "Fitness", text = paste(ValueValidator(ReactValue$Fitness)))
        } else {
          shinyFeedback::hideFeedback(inputId = "Fitness")
        }
      }
      
      if (ReactValue$setPerPlate == 1) {
        if (NonNegValueValidator(ReactValue$Residual) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "Residual", text = paste(NonNegValueValidator(ReactValue$Residual)))
        } else {
          shinyFeedback::hideFeedback(inputId = "Residual")
        }
      } else {
        if (DeathValidator(ReactValue$Residual) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "Residual", text = paste(DeathValidator(ReactValue$Residual)))
        } else {
          shinyFeedback::hideFeedback(inputId = "Residual")
        }
      }
      
      if (NonNegValueValidator(ReactValue$Inoculum) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "Inoculum", text = paste(NonNegValueValidator(ReactValue$Inoculum)))
      } else {
        shinyFeedback::hideFeedback(inputId = "Inoculum")
      }
      
      if (DeathValidator(ReactValue$Death) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "Death", text = paste(DeathValidator(ReactValue$Death)))
      } else {
        shinyFeedback::hideFeedback(inputId = "Death")
      }
      
      if ((ReactValue$setCV == 1) && (ReactValue$setPerPlate == 1)) {
        if (NonNegValueValidator(ReactValue$CV) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "CV", text = paste(NonNegValueValidator(ReactValue$CV)))
        } else {
          shinyFeedback::hideFeedback(inputId = "CV")
        }
      }
    }
    
    if (ReactValue$setSelective == 1) {
      if (ValueValidator(ReactValue$VolumeSelective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "VolumeSelective", text = paste(ValueValidator(ReactValue$VolumeSelective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "VolumeSelective")
      }
      
      if (ValueValidator(ReactValue$DilutionSelective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "DilutionSelective", text = paste(ValueValidator(ReactValue$DilutionSelective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "DilutionSelective")
      }
      
    } else if (ReactValue$setSelective == 2) {
      if (PlatEffValidator(ReactValue$PlatingEfficiency) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "PlatingEfficiency", text = paste(PlatEffValidator(ReactValue$PlatingEfficiency)))
      } else {
        shinyFeedback::hideFeedback(inputId = "PlatingEfficiency")
      }
    }
    
    if (SelectiveValidator(ReactValue$CountsSelective) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "CountsSelective", text = paste(SelectiveValidator(ReactValue$CountsSelective)))
    } else {
      shinyFeedback::hideFeedback(inputId = "CountsSelective")
    }
    
    if (ReactValue$setNonselective==1) {
      if (ValueValidator(ReactValue$VolumeNonselective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "VolumeNonselective", text = paste(ValueValidator(ReactValue$VolumeNonselective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "VolumeNonselective")
      }
      
      if (ValueValidator(ReactValue$DilutionNonselective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "DilutionNonselective", text = paste(ValueValidator(ReactValue$DilutionNonselective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "DilutionNonselective")
      }
      
      if (ReactValue$setPerPlate == 1) {
        if (NonselectiveValidator(ReactValue$CountsNonselective) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "CountsNonselective", text = paste(NonselectiveValidator(ReactValue$CountsNonselective)))
        } else {
          shinyFeedback::hideFeedback(inputId = "CountsNonselective")
        }
      } else {
        if (NonselectivePerPlateValidator(ReactValue$CountsNonselective, ReactValue$CountsSelective) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "CountsNonselective", text = paste(NonselectivePerPlateValidator(ReactValue$CountsNonselective, ReactValue$CountsSelective)))
        } else {
          shinyFeedback::hideFeedback(inputId = "CountsNonselective")
        }
      }
      
    } else if (ReactValue$setNonselective==2) {
      if (ValueValidator(ReactValue$MeanCells) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "MeanCells", text = paste(ValueValidator(ReactValue$MeanCells)))
      } else {
        shinyFeedback::hideFeedback(inputId = "MeanCells")
      }
    }
    
    if (ReactValue$setSelective==1) {
      if (paste(ValueValidator(ReactValue$VolumeTotal), ValueValidator(ReactValue$VolumeSelective), ValueValidator(ReactValue$DilutionSelective), sep = "") == "") {
        if (PlatingValidator(ReactValue$VolumeSelective, ReactValue$DilutionSelective, ReactValue$VolumeTotal) != "") {
          ReactValue$errorsDetected <- TRUE
          ReactValue$FeedbackAmountSelective <- TRUE
          textInputError(inputId = "VolumeSelective", text = paste(PlatingValidator(ReactValue$VolumeSelective, ReactValue$DilutionSelective, ReactValue$VolumeTotal)))
          textInputError(inputId = "DilutionSelective", text = "")
          textInputError(inputId = "VolumeTotal", text = "")
        } else {
          ReactValue$FeedbackAmountSelective <- FALSE
          shinyFeedback::hideFeedback(inputId = "VolumeSelective")
          shinyFeedback::hideFeedback(inputId = "DilutionSelective")
          shinyFeedback::hideFeedback(inputId = "VolumeTotal")
        }
      }
    }
    
    if (ReactValue$setNonselective==1) {
      if (paste(ValueValidator(ReactValue$VolumeTotal), ValueValidator(ReactValue$VolumeNonselective), ValueValidator(ReactValue$DilutionNonselective), sep = "") == "") {
        if (PlatingValidator(ReactValue$VolumeNonselective, ReactValue$DilutionNonselective, ReactValue$VolumeTotal) != "") {
          ReactValue$errorsDetected <- TRUE
          ReactValue$FeedbackAmountNonselective <- TRUE
          textInputError(inputId = "VolumeNonselective", text = paste(PlatingValidator(ReactValue$VolumeNonselective, ReactValue$DilutionNonselective, ReactValue$VolumeTotal)))
          textInputError(inputId = "DilutionNonselective", text = "")
          textInputError(inputId = "VolumeTotal", text = "")
        } else {
          ReactValue$FeedbackAmountNonselective <- FALSE
          shinyFeedback::hideFeedback(inputId = "VolumeNonselective")
          shinyFeedback::hideFeedback(inputId = "DilutionNonselective")
          if (!ReactValue$FeedbackAmountSelective) shinyFeedback::hideFeedback(inputId = "VolumeTotal")
        }
      }
    }
    
    if (ReactValue$setModel) {
      if (ReactValue$setNonselective == 2) {
        if (paste(NonNegValueValidator(ReactValue$Inoculum), ValueValidator(ReactValue$MeanCells), sep = "") == "") {
          if (ReactValue$Inoculum >= ReactValue$MeanCells) {
            ReactValue$errorsDetected <- TRUE
            textInputError(inputId = "Inoculum", text = paste("Size of the inoculum must be smaller than mean culture size."))
            textInputError(inputId = "MeanCells", text = "")
          } else {
            shinyFeedback::hideFeedback(inputId = "Inoculum")
            shinyFeedback::hideFeedback(inputId = "MeanCells")
          }
        }
      } else if (ReactValue$setNonselective == 1) {
        if (ReactValue$setPerPlate == 1) {
          if (paste(NonNegValueValidator(ReactValue$Inoculum), NonselectiveValidator(ReactValue$CountsNonselective), ValueValidator(ReactValue$VolumeTotal), ValueValidator(ReactValue$VolumeNonselective), ValueValidator(ReactValue$DilutionNonselective), sep = "") == "") {
            if (ReactValue$Inoculum >= (mean(ReactValue$CountsNonselective) * ReactValue$DilutionNonselective / ReactValue$VolumeNonselective * ReactValue$VolumeTotal)) {
              ReactValue$errorsDetected <- TRUE
              textInputError(inputId = "Inoculum", text = paste("Size of the inoculum must be smaller than mean culture size."))
              textInputError(inputId = "CountsNonselective", text = "")
              textInputError(inputId = "DilutionNonselective", text = "")
              textInputError(inputId = "VolumeNonselective", text = "")
              textInputError(inputId = "VolumeTotal", text = "")
            } else {
              shinyFeedback::hideFeedback(inputId = "Inoculum")
              shinyFeedback::hideFeedback(inputId = "CountsNonselective")
              if (!ReactValue$FeedbackAmountNonselective) shinyFeedback::hideFeedback(inputId = "DilutionNonselective")
              if (!ReactValue$FeedbackAmountNonselective) shinyFeedback::hideFeedback(inputId = "VolumeNonselective")
              if (!ReactValue$FeedbackAmountNonselective && !ReactValue$FeedbackAmountSelective) shinyFeedback::hideFeedback(inputId = "VolumeTotal")
            }
          }
        } else {
          if (paste(NonNegValueValidator(ReactValue$Inoculum), NonselectivePerPlateValidator(ReactValue$CountsNonselective, ReactValue$CountsSelective), ValueValidator(ReactValue$VolumeTotal), ValueValidator(ReactValue$VolumeNonselective), ValueValidator(ReactValue$DilutionNonselective), sep = "") == "") {
            if (any(ReactValue$Inoculum >= (ReactValue$CountsNonselective * ReactValue$DilutionNonselective / ReactValue$VolumeNonselective * ReactValue$VolumeTotal))) {
              ReactValue$errorsDetected <- TRUE
              textInputError(inputId = "Inoculum", text = paste("Size of the inoculum must be smaller than each culture size."))
              textInputError(inputId = "CountsNonselective", text = "")
              textInputError(inputId = "DilutionNonselective", text = "")
              textInputError(inputId = "VolumeNonselective", text = "")
              textInputError(inputId = "VolumeTotal", text = "")
            } else {
              shinyFeedback::hideFeedback(inputId = "Inoculum")
              shinyFeedback::hideFeedback(inputId = "CountsNonselective")
              if (!ReactValue$FeedbackAmountNonselective) shinyFeedback::hideFeedback(inputId = "DilutionNonselective")
              if (!ReactValue$FeedbackAmountNonselective) shinyFeedback::hideFeedback(inputId = "VolumeNonselective")
              if (!ReactValue$FeedbackAmountNonselective && !ReactValue$FeedbackAmountSelective) shinyFeedback::hideFeedback(inputId = "VolumeTotal")
            }
          }
        }
      }
    }
  })
  return(reactive(list("VolumeTotal" = ReactValue$VolumeTotal, "Fitness" = ReactValue$Fitness, "VolumeSelective" = ReactValue$VolumeSelective, "DilutionSelective" = ReactValue$DilutionSelective,
                       "PlatingEfficiency" = ReactValue$PlatingEfficiency, "CountsSelective" = ReactValue$CountsSelective, "VolumeNonselective" = ReactValue$VolumeNonselective,
                       "DilutionNonselective" = ReactValue$DilutionNonselective, "CountsNonselective" = ReactValue$CountsNonselective, "MeanCells" = ReactValue$MeanCells, "CV" = ReactValue$CV,
                       "Lag" = ReactValue$Lag, "Residual" = ReactValue$Residual, "Death" = ReactValue$Death, "Inoculum" = ReactValue$Inoculum,
                       "model" = ReactValue$setModel, "errors" = ReactValue$errorsDetected, "setSel" = ReactValue$setSelective, "setNsel" = ReactValue$setNonselective, "setCV" = ReactValue$setCV, "setPerPlate" = ReactValue$setPerPlate,
                       "useLagFitness" = ReactValue$useLagFitness
  )))
}

countsPlatingUpdate <- function(input, output, session, usePreset) {
  inputsVec <- c("VolumeTotal", "VolumeSelective", "DilutionSelective", "CountsSelective", "VolumeNonselective", "DilutionNonselective",
                 "CountsNonselective", "Fitness", "PlatingEfficiency", "MeanCells", "CV", "Lag", "Residual", "Death", "Inoculum")
  if (usePreset==0) {useValues <- rep("", 16)}
  else if (usePreset==1) {useValues <- c("500", "100", "1", "1\n2\n3", "100", "1000000", "100\n200\n300", "0.9", "0.2", "1e+9", "0", "0", "0", "0", "0")}
  else if (usePreset==2) {useValues <- c("500", "150", "1", "7\n8\n9", "100", "1000000", "90\n100\n110", "1", "0.3", "5e+8", "0", "0", "0", "0", "0")}
  else {return(NULL)}
  for (i in 1:16) {
    updateTextInput(session, inputId = inputsVec[i], value = useValues[i])
  }
}

countsPlatingLoadDataset <- function(input, output, session, SettingsPlatingID = NULL, CountsRateUserInput, CountsStrain1UserInput, CountsStrain2UserInput,
                                     CountsStrain1FoldUserInput, CountsStrain2FoldUserInput, CountsStrain3FoldUserInput,
                                     CountsStrain4FoldUserInput, CountsStrain5FoldUserInput, CountsStrain6FoldUserInput,
                                     BatchCalcUserInput) {
  ReactValue <- reactiveValues()
  
  observe({
    updateSelectInput(session = session, inputId = "LoadDataset",
                      choices = c("-", "Rate-Strain", "Pvalue-Strain1", "Pvalue-Strain2",
                                  "Fold-Strain1", "Fold-Strain2", "Fold-Strain3",
                                  "Fold-Strain4", "Fold-Strain5", "Fold-Strain6",
                                  "Barna-glucose", "Barna-maltose", "Crane-1", "Crane-2", "Demerec", "Foster-1994",
                                  "Vaisman-2021-1", "Vaisman-2021-2", names(BatchCalcUserInput())))
  })
  
  observeEvent(input$LoadDataset, {
    ReactValue$CountsRateUserInput <- CountsRateUserInput()
    ReactValue$CountsStrain1UserInput <- CountsStrain1UserInput()
    ReactValue$CountsStrain2UserInput <- CountsStrain2UserInput()
    ReactValue$CountsStrain1FoldUserInput <- CountsStrain1FoldUserInput()
    ReactValue$CountsStrain2FoldUserInput <- CountsStrain2FoldUserInput()
    ReactValue$CountsStrain3FoldUserInput <- CountsStrain3FoldUserInput()
    ReactValue$CountsStrain4FoldUserInput <- CountsStrain4FoldUserInput()
    ReactValue$CountsStrain5FoldUserInput <- CountsStrain5FoldUserInput()
    ReactValue$CountsStrain6FoldUserInput <- CountsStrain6FoldUserInput()
    ReactValue$BatchCalcUserInput <- BatchCalcUserInput()
    dataset <- input$LoadDataset
    
    if (!is.null(dataset)) {
      if (dataset == "-") {return()}
      if (dataset == "Crane-1") {choices <- list("setModel" = FALSE,
                                                 "setSelective" = 1,
                                                 "setNonselective" = 2,
                                                 "CountsSelective" = "121\n129\n146\n173\n181\n185\n193\n207\n222\n241\n287",
                                                 "VolumeTotal" = "2000",
                                                 "VolumeSelective" = "200",
                                                 "DilutionSelective" = "1",
                                                 "MeanCells" = "3.6e9")}
      else if (dataset == "Crane-2") {choices <- list("setModel" = FALSE,
                                                      "setSelective" = 1,
                                                      "setNonselective" = 2,
                                                      "CountsSelective" = "82\n107\n133\n144\n154\n166\n224\n224\n234\n165\n308",
                                                      "VolumeTotal" = "2000",
                                                      "VolumeSelective" = "200",
                                                      "DilutionSelective" = "1",
                                                      "MeanCells" = "3.9e9")}
      else if (dataset == "Demerec") {choices <- list("setModel" = FALSE,
                                                      "setSelective" = 2,
                                                      "setNonselective" = 2,
                                                      "CountsSelective" = "33\n18\n839\n47\n13\n126\n48\n80\n9\n71\n196\n66\n28\n17\n27\n37\n126\n33\n12\n44\n28\n67\n730\n168\n44\n50\n583\n23\n17\n24",
                                                      "PlatingEfficiency" = "1",
                                                      "MeanCells" = "1.9e8")}
      else if (dataset == "Barna-glucose") {choices <- list("setModel" = TRUE,
                                                            "setSelective" = 2,
                                                            "setNonselective" = 2,
                                                            "CountsSelective" = "0\n0\n0\n0\n0\n0\n0\n0\n0\n1\n0\n0\n25\n0\n28\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n54\n0\n0\n21\n1\n2\n0\n0\n0\n1\n14\n1\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n12\n0\n0\n0\n0\n0\n1\n0\n1\n1\n0\n0\n3\n1\n0\n0\n0\n0\n0\n1\n0\n0\n0\n0\n0\n0\n0\n0\n0\n9\n0\n8\n0\n0\n10\n0\n0\n0\n0\n1\n0\n0\n0\n0\n0\n0\n177\n0\n2\n1\n9\n0",
                                                            "PlatingEfficiency" = "1",
                                                            "MeanCells" = "4.2e5",
                                                            "setCV" = "0",
                                                            "useLagFitness" = "1",
                                                            "CV" = "0.119",
                                                            "Inoculum" = "925",
                                                            "Residual" = "0",
                                                            "Lag"= "1.7",
                                                            "Poisson" = "0",
                                                            "Death" = "0")}
      else if (dataset == "Barna-maltose") {choices <- list("setModel" = TRUE,
                                                            "setSelective" = 2,
                                                            "setNonselective" = 2,
                                                            "CountsSelective" = "2\n0\n0\n2\n4\n6\n0\n0\n87\n3\n2\n0\n2\n1\n0\n4\n0\n1\n5\n0\n0\n0\n0\n2\n292\n1\n4\n0\n1\n3\n1\n0\n3\n4\n0\n11\n0\n1\n4\n2\n1\n2\n0\n1\n2\n0\n0\n1\n1\n8\n139\n1\n0\n0\n0\n2\n1\n0\n2\n1\n0\n2\n1\n2\n2\n1\n1\n1\n0\n1\n1\n2\n103\n1\n1\n19\n55\n2\n0\n0\n0\n3\n0\n1\n11\n0\n0\n2\n4\n4\n1\n1\n4\n0\n0\n0\n0\n51\n0\n0",
                                                            "PlatingEfficiency" = "1",
                                                            "setCV" = "0",
                                                            "useLagFitness" = "0",
                                                            "MeanCells" = "4.5e5",
                                                            "CV" = "0.267",
                                                            "Inoculum" = "1230",
                                                            "Residual" = "0",
                                                            "Lag"= "0",
                                                            "Poisson" = "0",
                                                            "Death" = "0")}
      else if (dataset == "Foster-1994") {choices <- list("setModel" = TRUE,
                                                          "setSelective" = 2,
                                                          "setNonselective" = 2,
                                                          "CountsSelective" = "20\n25\n16\n11\n22\n24\n8\n9\n58\n25\n23\n10\n11\n8\n17\n9\n29\n22\n41\n10\n9\n12\n41\n16\n54\n15\n18\n9\n23\n18\n26\n25\n8\n22\n53\n6\n24\n11\n15\n30\n11\n11\n23\n9\n13\n17\n34\n21\n9\n9\n6\n22\n11\n220\n8\n56\n24\n32\n395\n18",
                                                          "PlatingEfficiency" = "1",
                                                          "setCV" = "0",
                                                          "useLagFitness" = "0",
                                                          "MeanCells" = "6.16e8",
                                                          "CV" = "0",
                                                          "Inoculum" = "0",
                                                          "Residual" = "7",
                                                          "Lag"= "0",
                                                          "Poisson" = "0",
                                                          "Death" = "0")}
      else if (dataset == "Vaisman-2021-1") {choices <- list("setModel" = FALSE,
                                                             "setSelective" = 1,
                                                             "setNonselective" = 1,
                                                             "CountsSelective" = "42\n38\n42\n28\n48\n34\n44\n79\n9\n36\n21\n32\n29\n20\n26\n33\n13\n30\n30\n24\n33\n35\n25\n14",
                                                             "CountsNonselective" = "214\n203\n210\n151\n208\n166\n187\n209\n146\n194\n160\n215\n185\n142\n168\n174\n178\n149\n177\n189\n200\n196\n183\n152",
                                                             "VolumeTotal" = "2000",
                                                             "VolumeSelective" = "100",
                                                             "DilutionSelective" = "1",
                                                             "VolumeNonselective" = "100",
                                                             "DilutionNonselective" = "1e6")}
      else if (dataset == "Vaisman-2021-2") {choices <- list("setModel" = FALSE,
                                                             "setSelective" = 1,
                                                             "setNonselective" = 1,
                                                             "CountsSelective" = "478\n92\n97\n254\n390\n141\n639\n439\n233\n409\n639\n313\n405\n310\n317\n350\n473\n343\n291\n82\n361\n457\n350\n317\n350\n361\n353\n285\n350\n332\n428\n387\n169\n237\n328\n478\n365\n296\n364\n409\n849\n1337\n306",
                                                             "CountsNonselective" = "120\n95\n101\n135\n169\n126\n102\n113\n91\n183\n128\n119\n184\n97\n169\n109\n88\n122\n77\n64\n120\n119\n98\n95\n130\n120\n77\n93\n83\n123\n183\n93\n203\n158\n120\n81\n108\n108\n144\n145\n183\n96\n158",
                                                             "VolumeTotal" = "2000",
                                                             "VolumeSelective" = "100",
                                                             "DilutionSelective" = "10",
                                                             "VolumeNonselective" = "100",
                                                             "DilutionNonselective" = "1e6")}
      else if (dataset == "Rate-Strain") {choices <- ReactValue$CountsRateUserInput}
      else if (dataset == "Pvalue-Strain1") {choices <- ReactValue$CountsStrain1UserInput}
      else if (dataset == "Pvalue-Strain2") {choices <- ReactValue$CountsStrain2UserInput}
      else if (dataset == "Fold-Strain1") {choices <- ReactValue$CountsStrain1FoldUserInput}
      else if (dataset == "Fold-Strain2") {choices <- ReactValue$CountsStrain2FoldUserInput}
      else if (dataset == "Fold-Strain3") {choices <- ReactValue$CountsStrain3FoldUserInput}
      else if (dataset == "Fold-Strain4") {choices <- ReactValue$CountsStrain4FoldUserInput}
      else if (dataset == "Fold-Strain5") {choices <- ReactValue$CountsStrain5FoldUserInput}
      else if (dataset == "Fold-Strain6") {choices <- ReactValue$CountsStrain6FoldUserInput}
      else if (!any(is.na(ReactValue$BatchCalcUserInput))) {if (dataset %in% names(ReactValue$BatchCalcUserInput)) {choices <- ReactValue$BatchCalcUserInput[[dataset]]}}
      else {choices <- NULL}
      
      suppressWarnings({
        if (!is.null(choices)) {
          if (!is.null(SettingsPlatingID)) {
            command <- paste("document.getElementById('", SettingsPlatingID, "-setPerPlate1').click();", sep = "")
            shinyjs::runjs(command)
          }
          for (name in names(choices)) {
            if (name %in% c("VolumeTotal", "VolumeSelective", "DilutionSelective", "VolumeNonselective", "DilutionNonselective", "CountsSelective",
                            "Fitness", "PlatingEfficiency", "MeanCells", "CV", "Lag", "Residual", "Death", "Inoculum", "CountsNonselective") && all(!is.na(choices[[name]]))) {
              if (length(choices[[name]]) > 1) {choice <- paste(as.character(choices[[name]]), collapse = "\n")} else {choice <- as.character(choices[[name]])}
              updateTextInput(session = session, inputId = name, value = choice)
            }
            else if (name %in% c("useLagFitness", "setCV")) {shinyWidgets::updateAwesomeRadio(session = session, inputId = name, selected = as.numeric(choices[[name]]))}
            else if (name %in% c("setSelective", "setSel")) {shinyWidgets::updateAwesomeRadio(session = session, inputId = "setSelective", selected = as.numeric(choices[[name]]))}
            else if (name %in% c("setNonselective", "setNsel")) {shinyWidgets::updateAwesomeRadio(session = session, inputId = "setNonselective", selected = as.numeric(choices[[name]]))}
            else if (name %in% c("setModel", "model")) {shinyWidgets::updateAwesomeCheckbox(session = session, inputId = "setModel", value = as.logical(choices[[name]]))}
          }
        }
      })
    }
  })
}

PowerModuleUI <- function(id, usePreset) {
  ns <- NS(id)
  if (usePreset==1) {useValues <- c("20", "1e-9", "1", "0.2", "1e9", "0", "0", "0", "0", "0")}
  else if (usePreset==2) {useValues <- c("20", "2e-9", "1", "0.5", "1.5e9", "0", "0", "0", "0", "0")}
  tagList(
    shinyFeedback::useShinyFeedback(),
    fluidRow(
      column(
        width = 12,
        conditionalPanel(
          condition = "input.sampleSizeOrPower == 0",
          textInput(
            inputId = ns("SampleSize"),
            label = HTML(paste("Sample size:", infoTooltip("Please provide a positive integer."))),
            value = useValues[1],
            width = "100%"
          )
        ),
        textInput(
          inputId = ns("MutationRate"),
          label = HTML(paste("Mutation rate:", infoTooltip("Please provide a positive number so that 0&nbsp;&lt;&nbsp;&mu;&nbsp;&le;&nbsp;1. Scientific notation can be used, e.g. 2.5e-9."))),
          value = useValues[2],
          width = "100%"
        ),
        textInput(
          inputId = ns("MeanCells"),
          label = HTML(paste("Average number of cells in culture:", infoTooltip("Please provide a positive number. Scientific notation can be used, e.g. 4e9."))),
          value = useValues[5],
          width = "100%"
        ),
        hr(),
        textInput(
          inputId = ns("PlatingEfficiency"),
          label = HTML(paste("Plating efficiency:", infoTooltip("Please provide a positive number so that 0&nbsp;&lt;&nbsp;&epsilon;&nbsp;&le;&nbsp;1."))),
          value = useValues[4],
          width = "100%"
        ),
        shinyWidgets::awesomeRadio(
          inputId = ns("useLagFitness"),
          label = HTML(paste("Specify phenotypic lag or fitness:", infoTooltip("Choose whether you want to input phenotypic lag or mutant fitness."))),
          choices = c("Don't use" = 0, "Specify phenotypic lag" = 1, "Specify mutant fitness" = 2),
          selected = 0,
          inline = TRUE
        ),
        conditionalPanel(
          condition="input.useLagFitness==1",
          ns = NS(id),
          textInput(
            inputId = ns("Lag"),
            label = HTML(paste("Phenotypic lag (generations):", infoTooltip("Please provide a positive number."))),
            value = useValues[7],
            width = "100%"
          )
        ),
        conditionalPanel(
          condition="input.useLagFitness==2",
          ns = NS(id),
          textInput(
            inputId = ns("Fitness"),
            label = HTML(paste("Mutant relative fitness:", infoTooltip("Please provide a positive number. Fitness &lt;&nbsp;1 means that mutants grow <u>slower</u> than non-mutants, while &gt;&nbsp;1 means that mutants grow faster. Fitness equal to 1 is equivalent to L-C model.<br><b>Note: underflow may be encountered when when both plating efficiency and fitness are much smaller than 0.1.</b>"))),
            value = useValues[3],
            width = "100%"
          )
        ),
        textInput(
          inputId = ns("Death"),
          label = HTML(paste("Death rate:", infoTooltip("Please provide a non-negative number that is a fraction of the birth rate (&ge;&nbsp;0 but &lt;&nbsp;1)."))),
          value = useValues[9],
          width = "100%"
        ),
        textInput(
          inputId = ns("Residual"),
          label = div(id = paste(id, "Residual", "labeltext", sep = "-"), HTML(paste("Residual mutations:", infoTooltip("Please provide a non-negative number. Floating-point numbers will be rounded down to nearest integer. Residual mutations apply only to the portion of culture that was plated.")))),
          value = useValues[8],
          width = "100%"
        ),
        textInput(
          inputId = ns("Inoculum"),
          label = HTML(paste("Size of the inoculum (number of cells):", infoTooltip("Please provide a non-negative number. Floating-point numbers will be rounded down to nearest integer. Inoculum size cannot exceed the total number of cells in grown culture."))),
          value = useValues[10],
          width = "100%"
        ),
        textInput(
          inputId = ns("CV"),
          label = HTML(paste("Coefficient of variation:", infoTooltip("Please provide a positive number. Scientific notation can be used, e.g. 2e-2. Values smaller than 1e-6 will be set to 0, and values bigger than 10 will be set to 10. Leave 0 if you do not want to use the coefficient of variation."))),
          value = useValues[6],
          width = "100%",
        )
      )
    )
  )
}

PowerModule <- function(input, output, session, sampleSizeOrPower) {
  ReactValue <- reactiveValues()
  
  observe({
    if (sampleSizeOrPower() == 0) {
      ReactValue$SampleSize <- numerise(input$SampleSize)
    } else {
      ReactValue$SampleSize <- NA
    }
    
    ReactValue$useLagFitness <- input$useLagFitness
    
    if (ReactValue$useLagFitness == 1) {
      ReactValue$Lag <- numerise(input$Lag)
    } else {
      ReactValue$Lag <- 0
    }
    if (ReactValue$useLagFitness == 2) {
      ReactValue$Fitness <- numerise(input$Fitness)
    } else {
      ReactValue$Fitness <- 1
    }
    ReactValue$Death <- numerise(input$Death)
    ReactValue$Residual <- numerise(input$Residual)
    ReactValue$Inoculum <- numerise(input$Inoculum)
    ReactValue$CV <- numerise(input$CV)
    ReactValue$MeanCells <- numerise(input$MeanCells)
    ReactValue$MutationRate <- numerise(input$MutationRate)
    ReactValue$PlatingEfficiency <- numerise(input$PlatingEfficiency)
    
    ReactValue$errorsDetected <- FALSE
    
    if (sampleSizeOrPower() == 0) {
      if (PosValueValidator(ReactValue$SampleSize) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "SampleSize", text = paste(PosValueValidator(ReactValue$SampleSize)))
      } else {
        shinyFeedback::hideFeedback(inputId = "SampleSize")
      }
    }
    
    if (ReactValue$useLagFitness == 1) {
      if (NonNegValueValidator(ReactValue$Lag) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "Lag", text = paste(NonNegValueValidator(ReactValue$Lag)))
      } else {
        shinyFeedback::hideFeedback(inputId = "Lag")
      }
    }
    
    if (ReactValue$useLagFitness == 2) {
      if (ValueValidator(ReactValue$Fitness) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "Fitness", text = paste(ValueValidator(ReactValue$Fitness)))
      } else {
        shinyFeedback::hideFeedback(inputId = "Fitness")
      }
    }
    
    if (PowerValidator(ReactValue$MutationRate) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "MutationRate", text = paste(PowerValidator(ReactValue$MutationRate)))
    } else {
      shinyFeedback::hideFeedback(inputId = "MutationRate")
    }
    
    if (NonNegValueValidator(ReactValue$Residual) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "Residual", text = paste(NonNegValueValidator(ReactValue$Residual)))
    } else {
      shinyFeedback::hideFeedback(inputId = "Residual")
    }
    
    if (NonNegValueValidator(ReactValue$Inoculum) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "Inoculum", text = paste(NonNegValueValidator(ReactValue$Inoculum)))
    } else {
      shinyFeedback::hideFeedback(inputId = "Inoculum")
    }
    
    if (DeathValidator(ReactValue$Death) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "Death", text = paste(DeathValidator(ReactValue$Death)))
    } else {
      shinyFeedback::hideFeedback(inputId = "Death")
    }
    
    if (NonNegValueValidator(ReactValue$CV) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "CV", text = paste(NonNegValueValidator(ReactValue$CV)))
    } else {
      shinyFeedback::hideFeedback(inputId = "CV")
    }
    
    if (PlatEffValidator(ReactValue$PlatingEfficiency) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "PlatingEfficiency", text = paste(PlatEffValidator(ReactValue$PlatingEfficiency)))
    } else {
      shinyFeedback::hideFeedback(inputId = "PlatingEfficiency")
    }
    
    if (ValueValidator(ReactValue$MeanCells) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "MeanCells", text = paste(ValueValidator(ReactValue$MeanCells)))
    } else {
      shinyFeedback::hideFeedback(inputId = "MeanCells")
    }
    
    if (paste(NonNegValueValidator(ReactValue$Inoculum), ValueValidator(ReactValue$MeanCells), sep = "") == "") {
      if (ReactValue$Inoculum >= ReactValue$MeanCells) {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "Inoculum", text = paste("Size of the inoculum must be smaller than mean culture size."))
        textInputError(inputId = "MeanCells", text = "")
      } else {
        shinyFeedback::hideFeedback(inputId = "Inoculum")
        shinyFeedback::hideFeedback(inputId = "MeanCells")
      }
    }
  })
  return(reactive(list("MutationRate" = ReactValue$MutationRate, "Fitness" = ReactValue$Fitness, "PlatingEfficiency" = ReactValue$PlatingEfficiency,
                       "MeanCells" = ReactValue$MeanCells, "CV" = ReactValue$CV, "Lag" = ReactValue$Lag, "Residual" = ReactValue$Residual,
                       "Death" = ReactValue$Death, "Inoculum" = ReactValue$Inoculum, "SampleSize" = ReactValue$SampleSize, "errors" = ReactValue$errorsDetected)))
}

PowerModuleUpdate <- function(input, output, session, usePreset) {
  inputsVec <- c("SampleSize", "MutationRate", "Fitness", "PlatingEfficiency", "MeanCells", "CV", "Lag", "Residual", "Death", "Inoculum")
  if (length(usePreset)==11) {useValues <- usePreset; shinyWidgets::updateAwesomeRadio(inputId = "useLagFitness", selected = as.numeric(usePreset[11]))}
  else if (usePreset==0) {useValues <- rep("", 10)}
  else if (usePreset==1) {useValues <- c("20", "1e-9", "1", "0.2", "1e9", "0", "0", "0", "0", "0")}
  else if (usePreset==2) {useValues <- c("20", "2e-9", "1", "0.5", "1.5e9", "0", "0", "0", "0", "0")}
  else {return(NULL)}
  for (i in 1:10) {
    updateTextInput(session, inputId = inputsVec[i], value = useValues[i])
  }
}

abbrevsUI <- function(id, alpha=TRUE) {
  ns <- NS(id)
  tagList(
    h5("Abbreviations:"),
    h6(HTML("&epsilon; &mdash; plating efficiency")),
    h6(HTML("N<sub>t</sub> &mdash; total number of cells in the culture")),
    h6(HTML("CV &mdash; coefficient of variation of N<sub>t</sub>")),
    h6(HTML("&rho; &mdash; relative fitness of the mutant cells")),
    h6(HTML("&lambda; &mdash; mean phenotypic lag (generations)")),
    h6(HTML("d &mdash; death rate (as a fraction of growth rate)")),
    h6(HTML("m<sub>p</sub> &mdash; number of residual mutations")),
    h6(HTML("&phi; &mdash; size of inoculum relative to final culture size")),
    h6(HTML("m &mdash; number of mutations per culture")),
    h6(HTML("m<sup>95&percnt;&ndash;</sup> &mdash; lower limit of 95&percnt; CI for m")),
    h6(HTML("m<sup>95&percnt;&plus;</sup> &mdash; upper limit of 95&percnt; CI for m")),
    h6(HTML("&mu; &mdash; mutation rate per cell per generation")),
    h6(HTML("&mu;<sup>95&percnt;&ndash;</sup> &mdash; lower limit of 95&percnt; CI for &mu;")),
    h6(HTML("&mu;<sup>95&percnt;&plus;</sup> &mdash; upper limit of 95&percnt; CI for &mu;")),
    if (alpha==TRUE) {h6(HTML("&alpha; &mdash; significance level"))}
  )
}