### helper functions for mlemur GUI ###

# Numeriser
numerise <- function(x) {
  y <- strsplit(
    gsub("\\s+\n|\\s+$|\n\\s+", "\n", x),
    "\n"
  )
  y <- lapply(y, function(x) x[nchar(x) >= 1])
  y <- suppressWarnings(as.numeric(unlist(y)))
  return(y)
}

# Input updater
loadDatasetWrapper <- function(id, rv, SettingsPlatingID) {
  callModule(countsPlatingLoadDataset,
             id = id,
             SettingsPlatingID = SettingsPlatingID,
             CountsRateUserInput = reactive(rv$CountsRateUserInput()),
             CountsStrain1UserInput = reactive(rv$CountsStrain1UserInput()),
             CountsStrain2UserInput = reactive(rv$CountsStrain2UserInput()),
             CountsStrain1FoldUserInput = reactive(rv$CountsStrain1FoldUserInput()),
             CountsStrain2FoldUserInput = reactive(rv$CountsStrain2FoldUserInput()),
             CountsStrain3FoldUserInput = reactive(rv$CountsStrain3FoldUserInput()),
             CountsStrain4FoldUserInput = reactive(rv$CountsStrain4FoldUserInput()),
             CountsStrain5FoldUserInput = reactive(rv$CountsStrain5FoldUserInput()),
             CountsStrain6FoldUserInput = reactive(rv$CountsStrain6FoldUserInput()),
             BatchCalcUserInput = reactive(rv$BatchDatasets)
  )
}

# Warning Tooltip
warningTooltip <- function(text) {
  text <- paste(text)
  output <- paste("<i class='fa fa-exclamation-triangle warning-tooltip' data-toggle='tooltip' style = 'color:#b94a48' data-placement='top' data-container='body' title='", text, "'></i>", sep = "")
  return(output)
}

# Info Toolip
infoTooltip <- function(text) {
  text <- paste(text)
  output <- paste("<i class='fa fa-question-circle info-tooltip' data-toggle='tooltip' style = 'color:#89a9cb' data-placement='top' data-container='body' title='", text, "'></i>", sep = "")
  return(output)
}

# Soft Warning Toolip
softWarningTooltip <- function(text) {
  text <- paste(text)
  output <- paste("<i class='fa fa-exclamation-circle info-tooltip' data-toggle='tooltip' style = 'color:#cccccc' data-placement='top' data-container='body' title='", text, "'></i>", sep = "")
  return(output)
}

# Shiny Feedback customisation
textInputError <- function(inputId, text) {
  shinyFeedback::hideFeedback(inputId = inputId)
  shinyFeedback::showFeedback(inputId = inputId,
                              text = text,
                              color = "#b94a48",
                              icon = shiny::icon("exclamation-triangle", lib="font-awesome"),
                              session = shiny::getDefaultReactiveDomain())
}

# Validator of Single Value inputs (Volume, Dilution, Mean Number)
ValueValidator <- function(x) {
  info <- paste("")
  if (length(x) == 0) {
    info <- paste("Please provide a positive value.\n")
  }
  else if (is.na(x) == TRUE) {
    info <- paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else if (x <= 0) {
    info <- paste("A positive number is required.\n")
  }
  return(info)
}

# Validator of Single Value non-negative inputs (Volume, Dilution, Mean Number)
NonNegValueValidator <- function(x) {
  info <- paste("")
  if (length(x) == 0) {
    info <- paste("Please provide a non-negative value.\n")
  }
  else if (is.na(x) == TRUE) {
    info <- paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else if (x < 0) {
    info <- paste("A non-negative number is required.\n")
  }
  return(info)
}

# Validator of Single Value positive inputs (Sample size)
PosValueValidator <- function(x) {
  info <- paste("")
  if (length(x) == 0) {
    info <- paste("Please provide a positive value.\n")
  }
  else if (is.na(x) == TRUE) {
    info <- paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else if (x <= 0) {
    info <- paste("A positive number is required.\n")
  }
  return(info)
}

# Validator of Death parameter
DeathValidator <- function(x) {
  info <- paste("")
  if (length(x) == 0) {
    info <- paste("Please provide a non-negative value.\n")
  }
  else if (is.na(x) == TRUE) {
    info <- paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else if (x < 0) {
    info <- paste("A non-negative number is required.\n")
  }
  else if (x >= 1) {
    info <- paste("The value must be smaller than unity.\n")
  }
  return(info)
}

# Validator of Power parameter
PowerValidator <- function(x) {
  info <- paste("")
  if (length(x) == 0) {
    info <- paste("Please provide a positive value.\n")
  }
  else if (is.na(x) == TRUE) {
    info <- paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else if (x <= 0) {
    info <- paste("A positive number is required.\n")
  }
  else if (x >= 1) {
    info <- paste("The value must be smaller than unity.\n")
  }
  return(info)
}

# Validator of Plating efficiency
PlatEffValidator <- function(x) {
  info <- paste("")
  if (length(x) == 0) {
    info <- paste("Please provide a positive value.\n")
  }
  else if (is.na(x) == TRUE) {
    info <- paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else if (x <= 0) {
    info <- paste("A positive number is required.\n")
  }
  else if (x > 1) {
    info <- paste("The value must be no bigger than unity.\n")
  }
  return(info)
}

# Validator of Colony Counts on Selective Medium
SelectiveValidator <- function(x) {
  info <- paste("")
  if (length(x) == 0) {
    info <- paste("This cell cannot be empty.\n")
  }
  else if (any(is.na(x)) == TRUE) {
    info <- paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else {
    if (any(x<0)) {
      info1 <- paste("A positive number is required.\n")
    } else {
      info1 <- paste("")
    }
    if (length(x) <= 1) {
      info2 <- paste("Please provide at least two positive values.\n")
    } else {
      info2 <- paste("")
    }
    if (sum(x) <= 0) {
      info3 <- paste("At least one colony count must be bigger than 0.\n")
    } else {
      info3 <- paste("")
    }
    info <- paste(info1, info2, info3, sep = "")
  }
  return(info)
}

# Validator of Colony Counts on Nonselective Medium
NonselectiveValidator <- function(x) {
  info <- paste("")
  if (length(x) == 0) {
    info <-
      paste("This cell cannot be empty.\n")
  }
  else if (any(is.na(x)) == TRUE) {
    info <-
      paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else {
    if (any(x< 0)) {
      info1 <-
        paste("A positive number is required.\n")
    } else {
      info1 <- paste("")
    }
    if (sum(x) <= 0) {
      info3 <-
        paste("At least one colony count must be bigger than 0.\n")
    } else {
      info3 <- paste("")
    }
    info <- paste(info1, info3, sep = "")
  }
  return(info)
}

# Validator of Colony Counts on Nonselective Medium - per plate version
NonselectivePerPlateValidator <- function(x, y) {
  info <- paste("")
  if (length(x) == 0) {
    info <- paste("This cell cannot be empty.\n")
  } else if (length(x) != length(y)) {
    info <- paste("Number of counts on non-selective medium must be the same as on selective medium if 'Use pairs of counts on both media for each culture' is selected.\n")
  }
  else if (any(is.na(x)) == TRUE) {
    info <- paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else {
    if (any(x< 0)) {
      info1 <- paste("A positive number is required.\n")
    } else {
      info1 <- paste("")
    }
    if (sum(x) <= 0) {
      info3 <- paste("At least one colony count must be bigger than 0.\n")
    } else {
      info3 <- paste("")
    }
    info <- paste(info1, info3, sep = "")
  }
  return(info)
}

# Validator of Plating Efficiency
PlatingValidator <- function(x, y, z) {
  info <- paste("")
  if (grepl("[a-zA-Z]",
            paste(ValueValidator(x), ValueValidator(y), ValueValidator(z))) == TRUE) {
    info <- paste("")
  } else if (x / y > z) {
    info <-
      paste("The amount plated (",
            as.character(x),
            "&divide;",
            as.character(y),
            ") cannot exceed ",
            as.character(z),
            ".\n",
            sep = "")
  }
  return(info)
}

# Validator of p-values
PvalueValidator <- function(x, y) {
  info <- paste("")
  if (is.null(y) == TRUE) {
    info <-
      paste("Correction method not selected.")
  }
  else if (length(x) == 0) {
    info <-
      paste("This cell cannot be empty.\n")
  }
  else if (any(is.na(x)) == TRUE) {
    info <-
      paste("Non-numeric characters (letters, commas, colons, etc.) are not allowed.\n")
  }
  else {
    if (length(x) <= 1) {
      info1 <-
        paste("At least two <em>P</em> values must be provided for correction.")
    } else {
      info1 <- paste("")
    }
    if (any(x<0)) {
      info2 <-
        paste("<em>P</em> value must be equal to or bigger than 0.")
    } else {
      info2 <- paste("")
    }
    if (any(x>1)) {
      info3 <-
        paste("<em>P</em> value cannot exceed 1.")
    } else {
      info3 <- paste("")
    }
    info <- paste(info1, info2, info3, sep = "")
  }
  return(info)
}

# Tab Names in batch calculator
adjustTabNames <- function(checkPlating, checkSelective, checkNonselective) {
  if (checkPlating != 0) {
    shinyjs::runjs("
                     $('#onDataChecked ul:first li:eq(0) a').html('Experiment parameters <i class=\"fa fa-circle\" style=\"color:#b94a48\"></i>')
                     ")
  } else {
    shinyjs::runjs("
                     $('#onDataChecked ul:first li:eq(0) a').html('Experiment parameters')
                     ")
  }
  if (checkSelective != 0) {
    shinyjs::runjs("
                     $('#onDataChecked ul:first li:eq(1) a').html('Counts on selective medium <i class=\"fa fa-circle\" style=\"color:#b94a48\"></i>')
                     ")
  } else {
    shinyjs::runjs("
                     $('#onDataChecked ul:first li:eq(1) a').html('Counts on selective medium')
                     ")
  }
  if (checkNonselective != 0) {
    shinyjs::runjs("
                     $('#onDataChecked ul:first li:eq(2) a').html('Counts on non-selective medium <i class=\"fa fa-circle\" style=\"color:#b94a48\"></i>')
                     ")
  } else {
    shinyjs::runjs("
                     $('#onDataChecked ul:first li:eq(2) a').html('Counts on non-selective medium')
                     ")
  }
}

# Default settings for Reactable
defaultTable <- function(data, name, first.col.as.rownames=FALSE) {
  if (first.col.as.rownames==FALSE) {data <- cbind(rownames(data), data)}
  colnames(data)[1] <- name
  text <- paste('reactable::reactable(data,
            sortable = FALSE, rownames = FALSE, pagination = FALSE, outlined = TRUE, wrap = TRUE, bordered = TRUE, compact = TRUE, defaultColDef = reactable::colDef(align = "center", html = TRUE),
            columns = list(', name, ' = reactable::colDef(
              style = list(position = "sticky", left = 0, background = "#fff", zIndex = 1),
              headerStyle = list(position = "sticky", left = 0, background = "#fff", zIndex = 1),
              class = "sticky left-col-1",
              headerClass = "sticky left-col-1"
            )
            )
  )', sep = "")
  eval(parse(text = text))
}

colourPvalue <- function(data) {
  if (is.finite(data)) {
    if (data == 0) {
      data2 <- signif(data, digits = 1)
    } else if (data < 0.01) {
      data2 <- formatC(data, format = "e", digits = 2)
    } else {
      data2 <- formatC(data, format = "f", digits = 3)
    }
    if (as.numeric(data) > 0.05) {
      data2 <- paste("<font color = '#b1b1b1'>", data2, "</font>", sep = "")
    }
  }
  return(data2)
}

colourEffsize <- function(data) {
  if (is.finite(data)) {
    if (data == 0) {
      data2 <- signif(data, digits = 1)
    } else if (data < 0.01) {
      data2 <- formatC(data, format = "e", digits = 2)
    } else {
      data2 <- formatC(data, format = "f", digits = 3)
    }
    if (as.numeric(data) < 0.1) {
      data2 <- paste("<font color = '#b1b1b1'>", data2, "</font>", sep = "")
    }
  }
  return(data2)
}

# Combo validator for Data in Batch calculator
comboValidator <- function(DataPlating, DataSelective, DataNonselective) {
  items <- c("VolumeTotal", "Fitness", "Lag", "Residual", "Death", "Inoculum",
             "VolumeSelective", "DilutionSelective", "PlatingEfficiency",
             "VolumeNonselective", "DilutionNonselective", "MeanCells", "CV")
  
  # Removing empty columns
  colsToRemove <- NULL
  for (i in 2:length(DataPlating)) {
    if (length(DataPlating[[i]])==0) {colPlatingEmpty <- TRUE}
    else if (all(is.na(DataPlating[i]))==TRUE) {colPlatingEmpty <- TRUE}
    else {colPlatingEmpty <- FALSE}
    if (length(DataSelective[[i]])==0) {colSelectiveEmpty <- TRUE}
    else if (all(is.na(DataSelective[i]))==TRUE) {colSelectiveEmpty <- TRUE}
    else {colSelectiveEmpty <- FALSE}
    if (length(DataNonselective[[i]])==0) {colNonselectiveEmpty <- TRUE}
    else if (all(is.na(DataNonselective[i]))==TRUE) {colNonselectiveEmpty <- TRUE}
    else {colNonselectiveEmpty <- FALSE}
    if (sum(c(colPlatingEmpty, colSelectiveEmpty, colNonselectiveEmpty))==3) {
      colsToRemove <- c(colsToRemove, i)
    }
  }
  
  if(!is.null(colsToRemove)) {
    DataPlating <- DataPlating[-colsToRemove]
    DataSelective <- DataSelective[-colsToRemove]
    DataNonselective <- DataNonselective[-colsToRemove]
  }
  
  colnames(DataPlating)[1] <- "Parameters"
  colnames(DataSelective[1]) <- "No"
  colnames(DataNonselective[1]) <- "No"
  DataPlating2 <- data.frame(matrix(data = numeric(), nrow = 13, ncol = length(DataPlating[-1]), dimnames = list(items, colnames(DataPlating[-1]))))
  DataSelective2 <- as.list(matrix(data = numeric(), nrow = 1, ncol = length(DataSelective[-1]), dimnames = list(NULL, colnames(DataSelective[-1]))))
  DataNonselective2 <- as.list(matrix(data = numeric(), nrow = 1, ncol = length(DataNonselective[-1]), dimnames = list(NULL, colnames(DataNonselective[-1]))))
  
  names(DataPlating2) <- colnames(DataPlating[-1])
  names(DataSelective2) <- colnames(DataSelective[-1])
  names(DataNonselective2) <- colnames(DataNonselective[-1])
  
  allnames <- names(DataPlating2)
  
  UnrecognizedRows <- 0
  SameRowAnotherTime <- 0
  RowsEmptyOmmited <- 0
  CellNonnumeric <- 0
  CellNonpositive <- 0
  ColumnPlatingTooBigSel <- 0
  ColumnPlatingTooBigNonsel <- 0
  ColumnTooMuchData <- 0
  ColumnInoculumTooBig <- 0
  DataMissingSel <- 0
  DataMissingNonsel <- 0
  checkNonselNumeric <- 0
  checkNonselNonNeg <- 0
  checkNonselZeros <- 0
  checkNonselColEmpty <- 0
  checkSelNumeric <- 0
  checkSelNonNeg <- 0
  checkSel2Counts <- 0
  checkSelZeros <- 0
  checkSelColEmpty <- 0
  colsToSkip <- vector(mode="numeric")
  colsWithErrors <- vector(mode="numeric")
  MeanCellsToUnity <- 0
  PlatEffToUnity <- 0
  PvaluePossible <- TRUE
  CVFromCounts <- rep(FALSE, length(DataPlating)-1)
  NoGoodColumns <- FALSE
  
  parameterNumbers <- list()
  if (length(DataPlating[[1]]) > 0) {
    rowsList <- 1:length(DataPlating[[1]])
  } else {
    rowsList <- NULL
  }
  
  # Mine parameters from user input
  for (item in items) {
    parameterNumbers[[item]] <- which(sapply(DataPlating[1], tolower)==tolower(item))
  }
  
  # Checking for unknown parameters
  if (!is.null(rowsList)) {
    for (i in rowsList) {
      if (!(i %in% as.vector(unlist(parameterNumbers)))) {
        DataPlating[i,1] <- paste(DataPlating[i,1], softWarningTooltip("mlemur does not recognize this parameter name. Check Example to see a list of permitted parameter names."), sep=" ")
        UnrecognizedRows <- UnrecognizedRows + 1
      }
    }
  }
  
  # Checking for redundant rows
  for (i in 1:13) {
    if (length(parameterNumbers[[i]])>1) {
      vector <- parameterNumbers[[i]]
      for (x in vector[-1]) {
        DataPlating[[x,1]] <- paste(DataPlating[x,1], softWarningTooltip("This row was skipped because there is another with the same name."), sep=" ")
        SameRowAnotherTime <- SameRowAnotherTime + 1
      }
      parameterNumbers[[i]] <- vector[1]
    }
  }
  
  ParameterExists <- rbind(DataPlating2, rep(NA, length(DataPlating2)), rep(NA, length(DataPlating2)))
  rownames(ParameterExists)[14] <- "CountsSelective"
  rownames(ParameterExists)[15] <- "CountsNonselective"
  
  if (length(DataPlating)>=2) {
    
    if(length(DataPlating)==2) {PvaluePossible <- FALSE}
    
    # Checking data for each column
    for (i in 2:length(DataPlating)) {
      for (item in items) {
        if (length(parameterNumbers[[item]]) == 0) {ParameterExists[item,i-1] <- FALSE}
        else if (is.na(DataPlating[[parameterNumbers[[item]],i]])) {ParameterExists[item,i-1] <- FALSE}
        else {ParameterExists[item,i-1] <- TRUE}
      }
      if (length(DataNonselective[[i]])==0) {ParameterExists["CountsNonselective",i-1] <- FALSE}
      else if (all(is.na(DataNonselective[[i]])==TRUE)) {ParameterExists["CountsNonselective",i-1] <- FALSE}
      else {ParameterExists["CountsNonselective",i-1] <- TRUE}
      if (length(DataSelective[[i]])==0) {ParameterExists["CountsSelective",i-1] <- FALSE}
      else if (all(is.na(DataSelective[[i]])==TRUE)) {ParameterExists["CountsSelective",i-1] <- FALSE}
      else {ParameterExists["CountsSelective",i-1] <- TRUE}
      
      TooMuchData <- 0
      
      # Checking if Mean Cells should be set to 1
      if (((all(c(ParameterExists["MeanCells",i-1], ParameterExists["VolumeNonselective",i-1], ParameterExists["DilutionNonselective",i-1], ParameterExists["CountsNonselective",i-1])==FALSE) & ParameterExists["PlatingEfficiency",i-1]==FALSE) | (all(c(ParameterExists["MeanCells",i-1], ParameterExists["VolumeTotal",i-1], ParameterExists["VolumeNonselective",i-1], ParameterExists["DilutionNonselective",i-1], ParameterExists["CountsNonselective",i-1])==FALSE) & ParameterExists["PlatingEfficiency",i-1]==TRUE)) & ((all(c(ParameterExists["MeanCells",i-1], ParameterExists["VolumeNonselective",i-1], ParameterExists["DilutionNonselective",i-1], ParameterExists["CountsNonselective",i-1], ParameterExists["PlatingEfficiency",i-1], ParameterExists["VolumeSelective",i-1], ParameterExists["DilutionSelective",i-1])==FALSE) & ParameterExists["VolumeTotal",i-1]==TRUE)==FALSE)) {
        DataPlating2["MeanCells",i-1] <- 1
        if (length(which(DataPlating[1]=="MeanCells"))==0) {
          DataPlating <- as.data.frame(rbind(as.matrix(DataPlating), c("MeanCells", rep(NA, length(DataPlating[-1])))))
        }
        DataPlating[which(DataPlating[1]=="MeanCells")[1],i] <- paste("1", softWarningTooltip("Mean culture size was set to unity because other value could not be found or calculated."), sep=" ")
        MeanCellsToUnity <- MeanCellsToUnity + 1
      } else if (ParameterExists["MeanCells",i-1]==FALSE & any(c(ParameterExists["VolumeTotal",i-1], ParameterExists["VolumeNonselective",i-1], ParameterExists["DilutionNonselective",i-1], ParameterExists["CountsNonselective",i-1])==FALSE)) {
        if (ParameterExists["VolumeTotal",i-1]==FALSE) {
          if (length(which(DataPlating[1]=="VolumeTotal"))==0) {
            DataPlating <- as.data.frame(rbind(as.matrix(DataPlating), c("VolumeTotal", rep(NA, length(DataPlating[-1])))))
          }
          DataPlating[which(DataPlating[1]=="VolumeTotal")[1],i] <- paste(softWarningTooltip("Cannot estimate mean culture size because this value is missing. The column will be skipped."), sep=" ")
          DataMissingNonsel <- DataMissingNonsel + 1
          colsToSkip <- c(colsToSkip, i-1)
        }
        if (ParameterExists["VolumeNonselective",i-1]==FALSE) {
          if (length(which(DataPlating[1]=="VolumeNonselective"))==0) {
            DataPlating <- as.data.frame(rbind(as.matrix(DataPlating), c("VolumeNonselective", rep(NA, length(DataPlating[-1])))))
          }
          DataPlating[which(DataPlating[1]=="VolumeNonselective")[1],i] <- paste(softWarningTooltip("Cannot estimate mean culture size because this value is missing. The column will be skipped."), sep=" ")
          DataMissingNonsel <- DataMissingNonsel + 1
          colsToSkip <- c(colsToSkip, i-1)
        }
        if (ParameterExists["DilutionNonselective",i-1]==FALSE) {
          if (length(which(DataPlating[1]=="DilutionNonselective"))==0) {
            DataPlating <- as.data.frame(rbind(as.matrix(DataPlating), c("DilutionNonselective", rep(NA, length(DataPlating[-1])))))
          }
          DataPlating[which(DataPlating[1]=="DilutionNonselective")[1],i] <- paste(softWarningTooltip("Cannot estimate mean culture size because this value is missing. The column will be skipped."), sep=" ")
          DataMissingNonsel <- DataMissingNonsel + 1
          colsToSkip <- c(colsToSkip, i-1)
        }
        if (ParameterExists["CountsNonselective",i-1]==FALSE) {
          checkNonselColEmpty <- checkNonselColEmpty + 1
          colsToSkip <- c(colsToSkip, i-1)
        }
      }
      
      # Checking if Plating efficiency should be set to 1
      if (((all(c(ParameterExists["PlatingEfficiency",i-1], ParameterExists["VolumeSelective",i-1], ParameterExists["DilutionSelective",i-1])==FALSE) & ParameterExists["MeanCells",i-1]==FALSE) | (all(c(ParameterExists["PlatingEfficiency",i-1], ParameterExists["VolumeTotal",i-1], ParameterExists["VolumeSelective",i-1], ParameterExists["DilutionSelective",i-1])==FALSE) & ParameterExists["MeanCells",i-1]==TRUE)) & ((all(c(ParameterExists["MeanCells",i-1], ParameterExists["VolumeNonselective",i-1], ParameterExists["DilutionNonselective",i-1], ParameterExists["CountsNonselective",i-1], ParameterExists["PlatingEfficiency",i-1], ParameterExists["VolumeSelective",i-1], ParameterExists["DilutionSelective",i-1])==FALSE) & ParameterExists["VolumeTotal",i-1]==TRUE)==FALSE)) {
        DataPlating2["PlatingEfficiency",i-1] <- 1
        if (length(which(DataPlating[1]=="PlatingEfficiency"))==0) {
          DataPlating <- as.data.frame(rbind(as.matrix(DataPlating), c("PlatingEfficiency", rep(NA, length(DataPlating[-1])))))
        }
        DataPlating[which(DataPlating[1]=="PlatingEfficiency")[1],i] <- paste("1", softWarningTooltip("Plating efficiency was set to unity because other value could not be found or calculated."), sep=" ")
        PlatEffToUnity <- PlatEffToUnity + 1
      } else if (ParameterExists["PlatingEfficiency",i-1]==FALSE & any(c(ParameterExists["VolumeTotal",i-1], ParameterExists["VolumeSelective",i-1], ParameterExists["DilutionSelective",i-1])==FALSE)) {
        if (ParameterExists["VolumeTotal",i-1]==FALSE) {
          if (length(which(DataPlating[1]=="VolumeTotal"))==0) {
            DataPlating <- as.data.frame(rbind(as.matrix(DataPlating), c("VolumeTotal", rep(NA, length(DataPlating[-1])))))
          }
          DataPlating[which(DataPlating[1]=="VolumeTotal")[1],i] <- paste(softWarningTooltip("Cannot estimate plating efficiency because this value is missing. The column will be skipped."), sep=" ")
          DataMissingSel <- DataMissingSel + 1
          colsToSkip <- c(colsToSkip, i-1)
        }
        if (ParameterExists["VolumeSelective",i-1]==FALSE) {
          if (length(which(DataPlating[1]=="VolumeSelective"))==0) {
            DataPlating <- as.data.frame(rbind(as.matrix(DataPlating), c("VolumeSelective", rep(NA, length(DataPlating[-1])))))
          }
          DataPlating[which(DataPlating[1]=="VolumeSelective")[1],i] <- paste(softWarningTooltip("Cannot estimate plating efficiency because this value is missing. The column will be skipped."), sep=" ")
          DataMissingSel <- DataMissingSel + 1
          colsToSkip <- c(colsToSkip, i-1)
        }
        if (ParameterExists["DilutionSelective",i-1]==FALSE) {
          if (length(which(DataPlating[1]=="DilutionSelective"))==0) {
            DataPlating <- as.data.frame(rbind(as.matrix(DataPlating), c("DilutionSelective", rep(NA, length(DataPlating[-1])))))
          }
          DataPlating[which(DataPlating[1]=="DilutionSelective")[1],i] <- paste(softWarningTooltip("Cannot estimate plating efficiency because this value is missing. The column will be skipped."), sep=" ")
          DataMissingSel <- DataMissingSel + 1
          colsToSkip <- c(colsToSkip, i-1)
        }
      }
      
      # Checking colony counts on selective medium
      if (ParameterExists["CountsSelective",i-1]==FALSE) {
        colnames(DataSelective)[i] <- paste(colnames(DataSelective)[i], softWarningTooltip("This column will be skipped because no counts were provided."), sep = " ")
        checkSelColEmpty <- checkSelColEmpty + 1
        colsToSkip <- c(colsToSkip, i-1)
      } else if (ParameterExists["CountsSelective",i-1] == TRUE) {
        for (x in 1:nrow(DataSelective)) {
          if (is.na(DataSelective[x,i]) == FALSE ) {
            if (is.na(as.numeric(DataSelective[x,i])) == TRUE) {
              DataSelective[x,i] <- paste(DataSelective[x,i], warningTooltip("Non-numeric characters detected."), sep = " ")
              checkSelNumeric <- checkSelNumeric + 1
              colsWithErrors <- c(colsWithErrors, i-1)
            } else if (as.numeric(DataSelective[x,i]) < 0) {
              DataSelective[x,i] <- paste(DataSelective[x,i], warningTooltip("A negative (&ge;0) value is required."), sep = " ")
              checkSelNonNeg <- checkSelNonNeg + 1
              colsWithErrors <- c(colsWithErrors, i-1)
            } else {
              DataSelective2[[i-1]] <- append(DataSelective2[[i-1]], as.numeric(DataSelective[x,i]))
            }
          }
        }
        DataSelective2[[i-1]] <- DataSelective2[[i-1]][-1]
        if (length(DataSelective2[[i-1]]) == 0) {
          colnames(DataSelective)[i] <- paste(colnames(DataSelective)[i], softWarningTooltip("No valid colony counts. Note: only correct inputs (numeric and &ge;0) were considered. The column will be skipped."), sep = " ")
          checkSelColEmpty <- checkSelColEmpty + 1
          colsToSkip <- c(colsToSkip, i-1)
        } else if (sum(DataSelective2[[i-1]]) <= 0) {
          colnames(DataSelective)[i] <- paste(colnames(DataSelective)[i], warningTooltip("At least 1 colony count must be bigger than 0. Note: only correct inputs (numeric and &ge;0) were considered."), sep = " ")
          checkSelZeros <- checkSelZeros + 1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (length(DataSelective2[[i-1]]) <= 1) {
          colnames(DataSelective)[i] <- paste(colnames(DataSelective)[i], warningTooltip("At least two colony counts must be provided. Note: only correct inputs (numeric and &ge;0) were considered."), sep = " ")
          checkSel2Counts <- checkSel2Counts + 1
          colsWithErrors <- c(colsWithErrors, i-1)
        }
      }
      
      # Checking for coefficient of variation
      if (ParameterExists["CV",i-1]==TRUE) {
        if (is.na(as.numeric(DataPlating[parameterNumbers$CV,i])) == TRUE) {
          DataPlating[parameterNumbers$CV,i] <- paste(DataPlating[parameterNumbers$CV,i], warningTooltip("Non-numeric characters detected."), sep = " ")
          CellNonnumeric <- CellNonnumeric + 1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (as.numeric(DataPlating[parameterNumbers$CV,i]) < 0) {
          DataPlating[parameterNumbers$CV,i] <- paste(DataPlating[parameterNumbers$CV,i], warningTooltip("Coefficient of variation must be a non-negative value (&ge;0)."), sep = " ")
          CellNonpositive <- CellNonpositive + 1
          colsWithErrors <- c(colsWithErrors, i-1)
        }  else {
          if (as.numeric(DataPlating[parameterNumbers$CV,i]) < 1e-5 && as.numeric(DataPlating[parameterNumbers$CV,i]) > 0) {
            DataPlating[parameterNumbers$CV,i] <- paste(DataPlating[parameterNumbers$CV,i], softWarningTooltip("Coefficient of variation is small and will be set to 1e-5 for calculations."), sep = " ")
          } else if (as.numeric(DataPlating[parameterNumbers$CV,i]) > 10) {
            DataPlating[parameterNumbers$CV,i] <- paste(DataPlating[parameterNumbers$CV,i], softWarningTooltip("Coefficient of variation is big and will be set to 10 for calculations."), sep = " ")
          }
          DataPlating2["CV",i-1] <- as.numeric(DataPlating[parameterNumbers$CV,i])
        }
      } else {
        DataPlating2["CV",i-1] <- 0
        CVFromCounts[i-1] <- TRUE
      }
      
      # Checking if Fitness should be evaluated
      if (ParameterExists["Fitness",i-1]==TRUE) {
        if (is.na(as.numeric(DataPlating[parameterNumbers$Fitness,i])) == TRUE) {
          DataPlating[parameterNumbers$Fitness,i] <- paste(DataPlating[parameterNumbers$Fitness,i], warningTooltip("Non-numeric characters detected."), sep = " ")
          CellNonnumeric <- CellNonnumeric+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (as.numeric(DataPlating[parameterNumbers$Fitness,i]) <= 0) {
          DataPlating[parameterNumbers$Fitness,i] <- paste(DataPlating[parameterNumbers$Fitness,i], warningTooltip("Fitness must be bigger than 0."), sep = " ")
          CellNonpositive <- CellNonpositive+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else {
          DataPlating2["Fitness",i-1] <- as.numeric(DataPlating[parameterNumbers$Fitness,i])
        }
        if (ParameterExists["Lag",i-1]==TRUE) {
          DataPlating[parameterNumbers$Fitness,i] <- paste(DataPlating[parameterNumbers$Fitness,i], warningTooltip("Fitness cannot be used together with Phenotypic lag."), sep = " ")
          TooMuchData <- 1
          colsWithErrors <- c(colsWithErrors, i-1)
        }
      }
      
      # Checking if Death should be evaluated
      if (ParameterExists["Death",i-1]==TRUE) {
        if (is.na(as.numeric(DataPlating[parameterNumbers$Death,i])) == TRUE) {
          DataPlating[parameterNumbers$Death,i] <- paste(DataPlating[parameterNumbers$Death,i], warningTooltip("Non-numeric characters detected."), sep = " ")
          CellNonnumeric <- CellNonnumeric+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (as.numeric(DataPlating[parameterNumbers$Death,i]) < 0 | as.numeric(DataPlating[parameterNumbers$Death,i]) >= 1) {
          DataPlating[parameterNumbers$Death,i] <- paste(DataPlating[parameterNumbers$Death,i], warningTooltip("Death rate must not smaller than 0 (&ge;0) and smaller than 1 (&lt;1)."), sep = " ")
          CellNonpositive <- CellNonpositive+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else {
          DataPlating2["Death",i-1] <- as.numeric(DataPlating[parameterNumbers$Death,i])
        }
      } else {
        DataPlating2["Death",i-1] <- 0
      }
      
      # Checking if Inoculum should be evaluated
      if (ParameterExists["Inoculum",i-1]==TRUE) {
        if (is.na(as.numeric(DataPlating[parameterNumbers$Inoculum,i])) == TRUE) {
          DataPlating[parameterNumbers$Inoculum,i] <- paste(DataPlating[parameterNumbers$Inoculum,i], warningTooltip("Non-numeric characters detected."), sep = " ")
          CellNonnumeric <- CellNonnumeric+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (as.numeric(DataPlating[parameterNumbers$Inoculum,i]) < 0) {
          DataPlating[parameterNumbers$Inoculum,i] <- paste(DataPlating[parameterNumbers$Inoculum,i], warningTooltip("Inoculum must be not smaller than 0."), sep = " ")
          CellNonpositive <- CellNonpositive+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else {
          DataPlating2["Inoculum",i-1] <- as.numeric(DataPlating[parameterNumbers$Inoculum,i])
        }
      } else {
        DataPlating2["Inoculum",i-1] <- 0
      }
      
      # Checking if Lag should be evaluated
      if (ParameterExists["Lag",i-1]==TRUE) {
        if (is.na(as.numeric(DataPlating[parameterNumbers$Lag,i])) == TRUE) {
          DataPlating[parameterNumbers$Lag,i] <- paste(DataPlating[parameterNumbers$Lag,i], warningTooltip("Non-numeric characters detected."), sep = " ")
          CellNonnumeric <- CellNonnumeric+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (as.numeric(DataPlating[parameterNumbers$Lag,i]) < 0) {
          DataPlating[parameterNumbers$Lag,i] <- paste(DataPlating[parameterNumbers$Lag,i], warningTooltip("Phenotypic lag must be not smaller than 0."), sep = " ")
          CellNonpositive <- CellNonpositive+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else {
          DataPlating2["Lag",i-1] <- as.numeric(DataPlating[parameterNumbers$Lag,i])
        }
        if (ParameterExists["Fitness",i-1]==TRUE) {
          DataPlating[parameterNumbers$Lag,i] <- paste(DataPlating[parameterNumbers$Lag,i], warningTooltip("Phenotypic lag cannot be used together with Fitness."), sep = " ")
          TooMuchData <- 1
          colsWithErrors <- c(colsWithErrors, i-1)
        }
      } else {
        DataPlating2["Lag",i-1] <- 0
      }
      
      # Checking if Residual should be evaluated
      if (ParameterExists["Residual",i-1]==TRUE) {
        if (is.na(as.numeric(DataPlating[parameterNumbers$Residual,i])) == TRUE) {
          DataPlating[parameterNumbers$Residual,i] <- paste(DataPlating[parameterNumbers$Residual,i], warningTooltip("Non-numeric characters detected."), sep = " ")
          CellNonnumeric <- CellNonnumeric+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (as.numeric(DataPlating[parameterNumbers$Residual,i]) < 0) {
          DataPlating[parameterNumbers$Residual,i] <- paste(DataPlating[parameterNumbers$Residual,i], warningTooltip("Residual must be not smaller than 0."), sep = " ")
          CellNonpositive <- CellNonpositive+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else {
          DataPlating2["Residual",i-1] <- as.numeric(DataPlating[parameterNumbers$Residual,i])
        }
      } else {
        DataPlating2["Residual",i-1] <- 0
      }
      
      # Checking if total volume should be evaluated
      if (ParameterExists["MeanCells",i-1]==FALSE | ParameterExists["PlatingEfficiency",i-1]==FALSE) {
        if (ParameterExists["VolumeTotal",i-1]==TRUE) {
          if (is.na(as.numeric(DataPlating[parameterNumbers$VolumeTotal,i])) == TRUE) {
            DataPlating[parameterNumbers$VolumeTotal,i] <- paste(DataPlating[parameterNumbers$VolumeTotal,i], warningTooltip("Non-numeric characters detected."), sep = " ")
            CellNonnumeric <- CellNonnumeric+1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else if (as.numeric(DataPlating[parameterNumbers$VolumeTotal,i]) <= 0) {
            DataPlating[parameterNumbers$VolumeTotal,i] <- paste(DataPlating[parameterNumbers$VolumeTotal,i], warningTooltip("Total culture volume must be a positive value (&gt;0)."), sep = " ")
            CellNonpositive <- CellNonpositive+1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else {
            DataPlating2["VolumeTotal",i-1] <- as.numeric(DataPlating[parameterNumbers$VolumeTotal,i])
          }
        }
      } else if (ParameterExists["MeanCells",i-1]==TRUE & ParameterExists["PlatingEfficiency",i-1]==TRUE) {
        if (length(as.vector(unlist(parameterNumbers[1]))) > 0) {
          for (x in as.vector(unlist(parameterNumbers[1]))) {
            DataPlating[x,i] <- paste(DataPlating[x,i], softWarningTooltip("This value was skipped because Plating efficiency and Mean Cells values are defined for this strain."), sep=" ")
          }
        }
      }
      # Checking if mean cells should be evaluated
      if (ParameterExists["MeanCells",i-1]==TRUE) {
        if (is.na(as.numeric(DataPlating[parameterNumbers$MeanCells,i])) == TRUE) {
          DataPlating[parameterNumbers$MeanCells,i] <- paste(DataPlating[parameterNumbers$MeanCells,i], warningTooltip("Non-numeric characters detected."), sep = " ")
          CellNonnumeric <- CellNonnumeric+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (as.numeric(DataPlating[parameterNumbers$MeanCells,i]) <= 0) {
          DataPlating[parameterNumbers$MeanCells,i] <- paste(DataPlating[parameterNumbers$MeanCells,i], warningTooltip("Mean number of cells must be a positive value (&gt;0)."), sep = " ")
          CellNonpositive <- CellNonpositive+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else {
          DataPlating2["MeanCells",i-1] <- as.numeric(DataPlating[parameterNumbers$MeanCells,i])
        }
        if (length(as.vector(na.exclude(DataNonselective[[i]]))) > 0) {
          colnames(DataNonselective)[i] <- paste(colnames(DataNonselective)[i], softWarningTooltip("Counts for this strain were skipped because Mean Cells value is defined."), sep=" ")
        }
        if (length(as.vector(unlist(parameterNumbers[c(10,11)]))) > 0) {
          for (x in as.vector(unlist(parameterNumbers[c(10,11)]))) {
            DataPlating[x,i] <- paste(DataPlating[x,i], softWarningTooltip("This value was skipped because Mean Cells value is defined for this strain."), sep=" ")
          }
        }
      } else {
        # Checking volume non-selective
        if (ParameterExists["VolumeNonselective",i-1]==TRUE) {
          if (is.na(as.numeric(DataPlating[parameterNumbers$VolumeNonselective,i])) == TRUE) {
            DataPlating[parameterNumbers$VolumeNonselective,i] <- paste(DataPlating[parameterNumbers$VolumeNonselective,i], warningTooltip("Non-numeric characters detected."), sep = " ")
            CellNonnumeric <- CellNonnumeric + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else if (as.numeric(DataPlating[parameterNumbers$VolumeNonselective,i]) <= 0) {
            DataPlating[parameterNumbers$VolumeNonselective,i] <- paste(DataPlating[parameterNumbers$VolumeNonselective,i], warningTooltip("Volume plated on non-selective medium must be a positive value (&gt;0)."), sep = " ")
            CellNonpositive <- CellNonpositive + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else {
            DataPlating2["VolumeNonselective",i-1] <- as.numeric(DataPlating[parameterNumbers$VolumeNonselective,i])
          }
        }
        # Checking dilution non-selective
        if (ParameterExists["DilutionNonselective",i-1]==TRUE) {
          if (is.na(as.numeric(DataPlating[parameterNumbers$DilutionNonselective,i])) == TRUE) {
            DataPlating[parameterNumbers$DilutionNonselective,i] <- paste(DataPlating[parameterNumbers$DilutionNonselective,i], warningTooltip("Non-numeric characters detected."), sep = " ")
            CellNonnumeric <- CellNonnumeric + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else if (as.numeric(DataPlating[parameterNumbers$DilutionNonselective,i]) <= 0) {
            DataPlating[parameterNumbers$DilutionNonselective,i] <- paste(DataPlating[parameterNumbers$DilutionNonselective,i], warningTooltip("Dilution factor on non-selective medium must be a positive value (&gt;0)."), sep = " ")
            CellNonpositive <- CellNonpositive + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else {
            DataPlating2["DilutionNonselective",i-1] <- as.numeric(DataPlating[parameterNumbers$DilutionNonselective,i])
          }
        }
        # Checking colony counts on non-selective medium
        if (ParameterExists["CountsNonselective",i-1] == TRUE) {
          for (x in 1:nrow(DataNonselective)) {
            if (is.na(DataNonselective[x,i]) == FALSE ) {
              if (is.na(as.numeric(DataNonselective[x,i])) == TRUE) {
                DataNonselective[x,i] <- paste(DataNonselective[x,i], warningTooltip("Non-numeric characters detected."), sep = " ")
                checkNonselNumeric <- checkNonselNumeric + 1
                colsWithErrors <- c(colsWithErrors, i-1)
              } else if (as.numeric(DataNonselective[x,i]) < 0) {
                DataNonselective[x,i] <- paste(DataNonselective[x,i], warningTooltip("A non-negative (&ge;0) value is required."), sep = " ")
                checkNonselNonNeg <- checkNonselNonNeg + 1
                colsWithErrors <- c(colsWithErrors, i-1)
              } else {
                DataNonselective2[[i-1]] <- append(DataNonselective2[[i-1]], as.numeric(DataNonselective[x,i]))
              }
            }
          }
          DataNonselective2[[i-1]] <- DataNonselective2[[i-1]][-1]
          if (length(DataNonselective2[[i-1]]) == 0) {
            checkNonselColEmpty <- checkNonselColEmpty + 1
            colsToSkip <- c(colsToSkip, i-1)
          } else if (sum(DataNonselective2[[i-1]]) <= 0) {
            colnames(DataNonselective)[i] <- paste(colnames(DataNonselective)[i], warningTooltip("At least 1 colony count must be bigger than 0. Note: only correct inputs (numeric and &ge;0) were considered."), sep = " ")
            checkNonselZeros <- checkNonselZeros + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else if (length(DataNonselective2[[i-1]]) <= 1) {
            if (CVFromCounts[i-1]) {
              colnames(DataNonselective)[i] <- paste(colnames(DataNonselective)[i], softWarningTooltip("Not enough colony counts to estimate CV, therefore CV will not be calculated. Note: only correct inputs (numeric and &ge;0) were considered."), sep = " ")
            }
            CVFromCounts[i-1] <- FALSE
          }
        } else {
          CVFromCounts[i-1] <- FALSE
        }
        # Checking plating fraction
        if (!any(is.na(c(DataPlating2["DilutionNonselective",i-1], DataPlating2["VolumeNonselective",i-1], DataPlating2["VolumeTotal",i-1])))) {
          if ((DataPlating2["VolumeNonselective",i-1] / DataPlating2["DilutionNonselective",i-1] / DataPlating2["VolumeTotal",i-1]) > 1) {
            DataPlating[parameterNumbers$DilutionNonselective,i] <- paste(DataPlating[parameterNumbers$DilutionNonselective,i], warningTooltip("Amount plated on non-selective medium (Volume/Dilution) exceeds Total culture volume."), sep = " ")
            DataPlating[parameterNumbers$VolumeNonselective,i] <- paste(DataPlating[parameterNumbers$VolumeNonselective,i], warningTooltip("Amount plated on non-selective medium (Volume/Dilution) exceeds Total culture volume."), sep = " ")
            ColumnPlatingTooBigNonsel <- ColumnPlatingTooBigNonsel + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          }
        }
      }
      # Checking if plating efficiency should be evaluated
      if (ParameterExists["PlatingEfficiency",i-1]==TRUE) {
        if (is.na(as.numeric(DataPlating[parameterNumbers$PlatingEfficiency,i])) == TRUE) {
          DataPlating[parameterNumbers$PlatingEfficiency,i] <- paste(DataPlating[parameterNumbers$PlatingEfficiency,i], warningTooltip("Non-numeric characters detected."), sep = " ")
          CellNonnumeric <- CellNonnumeric+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else if (as.numeric(DataPlating[parameterNumbers$PlatingEfficiency,i]) <= 0 | as.numeric(DataPlating[parameterNumbers$PlatingEfficiency,i]) > 1) {
          DataPlating[parameterNumbers$PlatingEfficiency,i] <- paste(DataPlating[parameterNumbers$PlatingEfficiency,i], warningTooltip("Plating efficiency must be between 0 and 1."), sep = " ")
          CellNonpositive <- CellNonpositive+1
          colsWithErrors <- c(colsWithErrors, i-1)
        } else {
          DataPlating2["PlatingEfficiency",i-1] <- as.numeric(DataPlating[parameterNumbers$PlatingEfficiency,i])
        }
        if (length(as.vector(unlist(parameterNumbers[c(7,8)]))) > 0) {
          for (x in as.vector(unlist(parameterNumbers[c(7,8)]))) {
            DataPlating[x,i] <- paste(DataPlating[x,i], softWarningTooltip("This value was skipped because Plating efficiency value is defined for this strain."), sep=" ")
          }
        }
      } else {
        # Checking volume selective
        if (ParameterExists["VolumeSelective",i-1]==TRUE) {
          if (is.na(as.numeric(DataPlating[parameterNumbers$VolumeSelective,i])) == TRUE) {
            DataPlating[parameterNumbers$VolumeSelective,i] <- paste(DataPlating[parameterNumbers$VolumeSelective,i], warningTooltip("Non-numeric characters detected."), sep = " ")
            CellNonnumeric <- CellNonnumeric + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else if (as.numeric(DataPlating[parameterNumbers$VolumeSelective,i]) <= 0) {
            DataPlating[parameterNumbers$VolumeSelective,i] <- paste(DataPlating[parameterNumbers$VolumeSelective,i], warningTooltip("Volume plated on selective medium must be a positive value (&gt;0)."), sep = " ")
            CellNonpositive <- CellNonpositive + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else {
            DataPlating2["VolumeSelective",i-1] <- as.numeric(DataPlating[parameterNumbers$VolumeSelective,i])
          }
        }
        # Checking dilution selective
        if (ParameterExists["DilutionSelective",i-1]==TRUE) {
          if (is.na(as.numeric(DataPlating[parameterNumbers$DilutionSelective,i])) == TRUE) {
            DataPlating[parameterNumbers$DilutionSelective,i] <- paste(DataPlating[parameterNumbers$DilutionSelective,i], warningTooltip("Non-numeric characters detected."), sep = " ")
            CellNonnumeric <- CellNonnumeric + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else if (as.numeric(DataPlating[parameterNumbers$DilutionSelective,i]) <= 0) {
            DataPlating[parameterNumbers$DilutionSelective,i] <- paste(DataPlating[parameterNumbers$DilutionSelective,i], warningTooltip("Dilution factor on selective medium must be a positive value (&gt;0)."), sep = " ")
            CellNonpositive <- CellNonpositive + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          } else {
            DataPlating2["DilutionSelective",i-1] <- as.numeric(DataPlating[parameterNumbers$DilutionSelective,i])
          }
        }
        # Checking plating fraction
        if (!any(is.na(c(DataPlating2["DilutionSelective",i-1], DataPlating2["VolumeSelective",i-1], DataPlating2["VolumeTotal",i-1])))) {
          if ((DataPlating2["VolumeSelective",i-1] / DataPlating2["DilutionSelective",i-1] / DataPlating2["VolumeTotal",i-1]) > 1) {
            DataPlating[parameterNumbers$DilutionSelective,i] <- paste(DataPlating[parameterNumbers$DilutionSelective,i], warningTooltip("Amount plated on selective medium (Volume/Dilution) exceeds Total culture volume."), sep = " ")
            DataPlating[parameterNumbers$VolumeSelective,i] <- paste(DataPlating[parameterNumbers$VolumeSelective,i], warningTooltip("Amount plated on selective medium (Volume/Dilution) exceeds Total culture volume."), sep = " ")
            ColumnPlatingTooBigSel <- ColumnPlatingTooBigSel + 1
            colsWithErrors <- c(colsWithErrors, i-1)
          }
        }
        # Checking if phenotypic lag and fitness have been specified together
        if (TooMuchData == 1) {
          ColumnTooMuchData <- ColumnTooMuchData + 1
        }
      }
      
      # Checking if inoculum size does not exceed culture size
      if (!any(is.na(c(DataPlating2["DilutionNonselective",i-1], DataPlating2["VolumeNonselective",i-1], DataPlating2["VolumeTotal",i-1], DataPlating2["Inoculum",i-1], DataNonselective2[[i-1]])))) {
        if (DataPlating2["Inoculum",i-1] >= (mean(DataNonselective2[[i-1]]) * DataPlating2["DilutionNonselective",i-1] / DataPlating2["VolumeNonselective",i-1] * DataPlating2["VolumeTotal",i-1])) {
          DataPlating[parameterNumbers$Inoculum,i] <- paste(DataPlating[parameterNumbers$Inoculum,i], warningTooltip("Inoculum is bigger than mean culture size."), sep = " ")
          ColumnInoculumTooBig <- ColumnInoculumTooBig + 1
          colsWithErrors <- c(colsWithErrors, i-1)
        }
      } else if (!any(is.na(c(DataPlating2["Inoculum",i-1], DataPlating2["MeanCells",i-1])))) {
        if (DataPlating2["Inoculum",i-1] >= DataPlating2["MeanCells",i-1]) {
          DataPlating[parameterNumbers$Inoculum,i] <- paste(DataPlating[parameterNumbers$Inoculum,i], warningTooltip("Inoculum is bigger than mean culture size."), sep = " ")
          ColumnInoculumTooBig <- ColumnInoculumTooBig + 1
          colsWithErrors <- c(colsWithErrors, i-1)
        }
      }
    }
    
    if (length(colsToSkip) > 0 || length(colsWithErrors) > 0) {
      badCols <- sort(unique(c(colsToSkip, colsWithErrors)))
      
      DataPlating2 <- DataPlating2[,-badCols,drop = FALSE]
      DataSelective2 <- DataSelective2[-badCols]
      DataNonselective2 <- DataNonselective2[-badCols]
      
      if (length(DataPlating2)==0) NoGoodColumns <- TRUE
    }
    
    colnamesToSkip <- vector(mode = "character")
    if (length(colsToSkip) > 0) {
      colsToSkip <- sort(unique(colsToSkip))
      
      for (i in colsToSkip) {
        colnames(DataPlating)[i+1] <- paste(colnames(DataPlating)[i+1], softWarningTooltip("This column will be skipped due to missing data."), sep = " ")
        colnames(DataSelective)[i+1] <- paste(colnames(DataSelective)[i+1], softWarningTooltip("This column will be skipped due to missing data."), sep = " ")
        colnames(DataNonselective)[i+1] <- paste(colnames(DataNonselective)[i+1], softWarningTooltip("This column will be skipped due to missing data."), sep = " ")
      }
      
      colnamesToSkip <- allnames[colsToSkip]
    }
    
    colnamesWithErrors <- vector(mode = "character")
    colnamesWithErrorsButSkipped <- vector(mode = "character")
    colsWithErrorsButSkipped <- vector(mode = "numeric")
    if (length(colsWithErrors) > 0) {
      colsWithErrors <- sort(unique(colsWithErrors))
      
      if (length(colsToSkip) > 0) {
        colsWithErrorsButSkipped <- colsWithErrors[which(colsWithErrors %in% colsToSkip)]
        if (length(which(colsWithErrors %in% colsToSkip)) > 0) colsWithErrors <- colsWithErrors[-which(colsWithErrors %in% colsToSkip)]
      }
      
      colnamesWithErrors <- allnames[colsWithErrors]
      
      if (length(colsWithErrorsButSkipped) > 0) {
        colnamesWithErrorsButSkipped <- allnames[colsWithErrorsButSkipped]
      }
    }
    
    errorList <- ""
    warningList <- ""
    
    warningList <- paste(ifelse(!is.null(colsToRemove), paste("Columns", paste(colsToRemove, collapse = ", "), "were removed because they were empty.\n", sep = " "), ""),
                         ifelse(length(colsToSkip) > 0, paste("Columns", paste(colnamesToSkip, collapse = ", "), "will be skipped due to missing data.\n", sep = " "), ""),
                         ifelse(PvaluePossible==FALSE, "<em>P</em> values cannot be calculated because there is only one valid dataset.\n", ""),
                         ifelse(MeanCellsToUnity != 0, paste("In", MeanCellsToUnity, "dataset(s), Total culture size was set to 1 because other data were not provided.\n", sep = " "), ""),
                         ifelse(PlatEffToUnity != 0, paste("In", PlatEffToUnity, "dataset(s), plating efficiency was set to 1 because other data were not provided.\n", sep = " "), ""),
                         ifelse(UnrecognizedRows != 0, paste("The are", UnrecognizedRows, "unrecognized input row(s) in Experiment parameters data.\n", sep = " "), ""),
                         ifelse(SameRowAnotherTime != 0, paste(SameRowAnotherTime, "row(s) in Experiment parameters are duplicated and therefore were skipped.\n", sep = " "), ""),
                         ifelse(RowsEmptyOmmited != 0, paste(RowsEmptyOmmited, "row(s) in Experiment parameters are empty and therefore were skipped.\n", sep = " "), ""),
                         sep = "")
    
    errorList <- paste(ifelse(length(colsWithErrors) > 0, paste("<i class='fa fa-exclamation-triangle' style = 'color:#b94a48'></i> <font color = '#b94a48'><b>Columns", paste(colnamesWithErrors, collapse = ", "), "contain fatal errors that need to be resolved before continuing.</b></font>\n", sep = " "), ""),
                       ifelse(length(colsWithErrorsButSkipped) > 0, paste("<i class='fa fa-exclamation-triangle' style = 'color:#b94a48'></i> Columns", paste(colnamesWithErrorsButSkipped, collapse = ", "), "contain errors but the column was skipped due to missing data.\n", sep = " "), ""),
                       ifelse(NoGoodColumns, paste("<i class='fa fa-exclamation-triangle' style = 'color:#b94a48'></i> <font color = '#b94a48'><b>It appears there are no usable datasets in your file.</b></font>\n", sep = " "), ""),
                       ifelse(CellNonnumeric != 0, paste(CellNonnumeric, "input(s) in Experiment parameters contain(s) non-numerical characters.\n", sep = " "), ""),
                       ifelse(CellNonpositive != 0, paste(CellNonpositive, "input(s) in Experiment parameters contain(s) values out of bounds.\n", sep = " "), ""),
                       ifelse(ColumnPlatingTooBigNonsel != 0, paste("In", ColumnPlatingTooBigNonsel, "dataset(s) in Counts on non-selective medium, plating efficiency is bigger than 1.\n", sep = " "), ""),
                       ifelse(checkNonselNumeric != 0, paste(checkNonselNumeric, "count(s) on non-selective medium contain(s) non-numerical characters.\n", sep = " "), ""),
                       ifelse(checkNonselNonNeg != 0, paste(checkNonselNonNeg, "count(s) on non-selective medium are smaller than 0.\n", sep = " "), ""),
                       ifelse(checkNonselZeros != 0, paste(checkNonselZeros, "dataset(s) in Counts on non-selective medium do not contain any positive value.\n", sep = " "), ""),
                       ifelse(ColumnPlatingTooBigSel != 0, paste("In", ColumnPlatingTooBigSel, "dataset(s) in Counts on selective medium, plating efficiency is bigger than 1.\n", sep = " "), ""),
                       ifelse(ColumnTooMuchData != 0, paste("In", ColumnTooMuchData, "dataset(s), both Phenotypic lag and Fitness have been specified. Phenotypic lag cannot be used together with Fitness.\n", sep = " "), ""),
                       ifelse(ColumnInoculumTooBig != 0, paste("In", ColumnInoculumTooBig, "dataset(s), Inoculum is bigger than average culture size.\n", sep = " "), ""),
                       ifelse(checkSelNumeric != 0, paste(checkSelNumeric, "count(s) on selective medium contain(s) numerical characters.\n", sep = " "), ""),
                       ifelse(checkSelNonNeg != 0, paste(checkSelNonNeg, "count(s) on selective medium are smaller than 0.\n", sep = " "), ""),
                       ifelse(checkSelZeros != 0, paste(checkSelZeros, "dataset(s) in Counts on selective medium do not contain any positive value.\n", sep = " "), ""),
                       ifelse(checkSel2Counts != 0, paste(checkSel2Counts, "dataset(s) in Counts on selective medium do not contain only one count.\n", sep = " "), ""),
                       sep = "")
    
    checkMsgPlating <- sum(MeanCellsToUnity, PlatEffToUnity, UnrecognizedRows, SameRowAnotherTime, RowsEmptyOmmited, CellNonnumeric, CellNonpositive)
    checkMsgSelective <- sum(DataMissingSel, ColumnPlatingTooBigSel, checkSelNumeric, checkSelNonNeg, checkSelZeros, checkSel2Counts, checkSelColEmpty)
    checkMsgNonselective <- sum(DataMissingNonsel, ColumnPlatingTooBigNonsel, checkNonselNumeric, checkNonselNonNeg, checkNonselZeros, checkNonselColEmpty)
    
    DataSelective <- DataSelective[-1]
    DataNonselective <- DataNonselective[-1]
    
    output <- list(DataPlating, DataSelective, DataNonselective, DataPlating2, DataSelective2, DataNonselective2, warningList, errorList, checkMsgPlating, checkMsgSelective, checkMsgNonselective, CVFromCounts, PvaluePossible)
    
    return(output)
    
  } else {return(NULL)}
  
}

removeParameters <- function(PlatingList, FinalSettings) {
  for (i in 1:nrow(FinalSettings)) {
    for (item in colnames(FinalSettings)[-6]) {
      number <- which(names(PlatingList[[i]])==item)
      if (!FinalSettings[i,item]) {PlatingList[[i]][[number]] <- NA}
    }
    if (!FinalSettings[i,6]) {PlatingList[[i]][["model"]] <- FALSE}
    if (PlatingList[[i]][["CV"]]!=0) {PlatingList[[i]][["setCV"]] <- TRUE}
  }
  return(PlatingList)
}