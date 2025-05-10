mlemurUI <- function(request) {
  fluidPage(
    shinyjs::useShinyjs(),
    shinyFeedback::useShinyFeedback(),
    withMathJax(),
    tags$head(
      tags$link(rel= "stylesheet", type = "text/css", href = "www/bootstrap.css"),
      tags$script(src = "www/javascript.js"),
      tags$link(rel="shortcut icon", href="www/mlemur.ico"),
      tags$div(HTML("<script type='text/x-mathjax-config'>
                                      MathJax.Hub.Config({
                                      tex2jax: {inlineMath: [['$','$']]}
                                      });
                                     </script>
                                    "))
    ),
    a(name="top"),
    navbarPage(fluid = FALSE,
               id = "navigation",
               title = div(img(src = ("www/mlemur.svg"), style="margin-top:-11px; padding-right:10px; padding-bottom:0px;", class="unselectable", draggable="false", dragstart="false;", height = 45)),
               windowTitle="mlemur: MLE Mutation Rate Calculator",
               footer = column(12,
                               hr(),
                               HTML("<center>mlemur: MLE Mutation Rate Calculator v0.9.6.5 | 2023 | GPL-2</center>"),
                               HTML("<center>This program is distributed in the hope that it will be useful,
                                    but WITHOUT ANY WARRANTY; without even the implied warranty of
                                    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
                                    GNU General Public License for more details.</center>"),
                               HTML("<center>Some code used in this program was inspired by rSalvador by Qi Zheng, distributed under GPL-2 licence. See LICENSE.MD for more details.</center>"),
                               br(),
                               br()
               ),
               collapsible = TRUE,
               #### Mutation Rate Tab ####
               tabPanel(
                 "Rate",
                 titlePanel("Mutation rate"),
                 hr(),
                 fluidRow(
                   column(
                     7,
                     settingsPlatingUI("SettingsRate", c(0, 1, 1)),
                     hr(),
                     countsPlatingUI("CountsRate", "SettingsRate", FALSE, usePreset = 1, defaultSettings = c(0, 1, 1)),
                     fluidRow(
                       column(
                         4,
                         shinyFeedback::loadingButton(
                           "calculate",
                           label = HTML("<i class = \"fas fa-calculator\"></i> Calculate"),
                           style = "width:100%;",
                           loadingLabel = "Calculating\U2026"
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("erase",
                                      label = "Clear all cells",
                                      width = "100%",
                                      icon = icon("trash", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("sample",
                                      label = "Load sample",
                                      width = "100%",
                                      icon = icon("sync", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       )
                     )
                   ),
                   column(
                     5,
                     h3("Results"),
                     htmlOutput(outputId = "errorBarRate"),
                     tags$style(type = "text/css", "#errorBarRate {white-space: pre-wrap; line-height:18px;}"),
                     shinyjs::hidden(
                       div(
                         id = "advanced",
                         reactable::reactableOutput("tableRate"),
                         br(),
                         fluidRow(
                           column(8,
                                  offset = 2,
                                  div(actionButton(inputId = "clip",
                                                   label = "Copy to Clipboard",
                                                   style = "text-align:center",
                                                   icon = icon("clipboard"),
                                                   width = "100%"))
                           )
                         ),
                         br(),
                         hr(),
                         fluidRow(
                           column(8,
                                  offset = 2,
                                  textInput(
                                    inputId = "datasetName",
                                    label = "Dataset name:",
                                    value = "Strain 1",
                                    width = "100%"
                                  ),
                                  actionButton("appendToReport",
                                               label = "Add to report",
                                               width = "100%",
                                               icon = icon("plus-square", lib = "font-awesome")
                                  ),
                                  br(),
                                  br(),
                                  div(downloadButton("dl", "Export report to XLSX"),
                                      style = "text-align:center")
                           )
                         ),
                         br(),
                         hr(),
                         abbrevsUI("abbrevsrate", FALSE)
                       )
                     )
                   )
                 )
               ),
               #### p-value Tab ####
               tabPanel(
                 title = HTML("<em>P</em> value"),
                 titlePanel(div(HTML(
                   "<em>P</em> value"
                 ))),
                 hr(),
                 fluidRow(
                   column(
                     7,
                     settingsPlatingUI("SettingsPval", c(0, 1, 2)),
                     hr(),
                     fluidRow(
                       column(
                         6,
                         h4(HTML("<b>Strain 1</b>")),
                         countsPlatingUI("CountsStrain1", "SettingsPval", TRUE, usePreset = 1, defaultSettings = c(0, 1, 2)),
                       ),
                       column(
                         6,
                         h4(HTML("<b>Strain 2</b>")),
                         countsPlatingUI("CountsStrain2", "SettingsPval", TRUE, usePreset = 2, defaultSettings = c(0, 1, 2)),
                       )
                     ),
                     br(),
                     fluidRow(
                       column(
                         4,
                         shinyFeedback::loadingButton(
                           "calculate2",
                           label = HTML("<i class = \"fas fa-calculator\"></i> Calculate"),
                           style = "width:100%;",
                           loadingLabel = "Calculating\U2026"
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("erase2",
                                      label = "Clear all cells",
                                      width = "100%",
                                      icon = icon("trash", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("sample2",
                                      label = "Load sample",
                                      width = "100%",
                                      icon = icon("sync", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       )
                     )
                   ),
                   column(5,
                          h3("Results"),
                          htmlOutput(outputId = "errorBarPvalue"),
                          tags$style(type = "text/css", "#errorBarPvalue {white-space: pre-wrap; line-height:18px;}"),
                          shinyjs::hidden(
                            div(
                              id = "advanced2",
                              reactable::reactableOutput("tablePvalue"),
                              hr(),
                              htmlOutput("pvalueinfo"),
                              br(),
                              fluidRow(
                                column(8,
                                       offset = 2,
                                       div(actionButton(inputId = "sendToPower",
                                                        label = "Send for power analysis",
                                                        style = "text-align:center",
                                                        icon = icon("paper-plane"),
                                                        width = "100%"))
                                )
                              ),
                              br(),
                              hr(),
                              abbrevsUI("abbrevspval", TRUE)
                            )
                          ))
                 )
               ),
               #### Pvalue Correction Tab ####
               tabPanel(
                 title = HTML("<em>P</em> correction"),
                 titlePanel(div(HTML(
                   "<em>P</em> value correction"
                 ))),
                 hr(),
                 fluidRow(
                   column(
                     7,
                     fluidRow(
                       column(
                         6,
                         shinyWidgets::awesomeRadio(
                           inputId = "correctionMethod",
                           label = "Select method of correction:",
                           choices = c("Bonferroni" = "bonferroni", "Bonferroni-Holm" = "holm", "Benjamini-Hochberg" = "BH"),
                           selected = "bonferroni"
                         ),
                         h6()
                       ),
                       column(
                         6,
                         textAreaInput(
                           inputId = "enteredPvalues",
                           label = HTML(paste("Paste <em>P</em> values here:", infoTooltip("Inputs must be within 0 and 1. Scientific notation can be used as needed. Paste one value under the other. Do not separate the values with a comma, semicolon, or any other character."))),
                           value = "0.064\n0.00007\n0.002",
                           rows = 10,
                           width = "100%",
                           resize = "vertical"
                         ),
                         br()
                       )
                     ),
                     fluidRow(
                       column(
                         4,
                         actionButton(
                           "calculate3",
                           label = "Calculate",
                           width = "100%",
                           class = "btn btn-primary",
                           icon = icon("calculator", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("erase3",
                                      label = "Clear all cells",
                                      width = "100%",
                                      icon = icon("trash", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("sample3",
                                      label = "Load sample",
                                      width = "100%",
                                      icon = icon("sync", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       )
                     )
                   ),
                   column(5,
                          h3("Results"),
                          htmlOutput(outputId = "errorBarCorrector"),
                          tags$style(type = "text/css", "#errorBarCorrector {white-space:pre-wrap; line-height:18px;}"),
                          shinyjs::hidden(
                            div(
                              id = "advanced3",
                              reactable::reactableOutput("tableCorrector"),
                              br(),
                              fluidRow(
                                column(8,
                                       offset = 2,
                                       div(actionButton(inputId = "clip2",
                                                        label = "Copy to Clipboard",
                                                        style = "text-align:center",
                                                        icon = icon("clipboard"),
                                                        width = "100%"))
                                )
                              )
                            )
                          ))
                 )
               ),
               #### Batch Tab ####
               tabPanel("Batch",
                        titlePanel(div(HTML(
                          "Batch mutation rate and <em>P</em> value"
                        ))),
                        hr(),
                        "XLS(X) files only! For instructions, check Help, or see example.",
                        a(href="www/example.xlsx", HTML("<i class=\"fa fa-download\" aria-hidden=\"true\"></i> Download example with comments"), download=NA, target="_blank"),
                        a(href="www/template.xlsx", HTML("<i class=\"fa fa-download\" aria-hidden=\"true\"></i> Download empty template to fill"), download=NA, target="_blank"),
                        fluidRow(
                          column(4,
                                 offset = 4,
                                 h3("Step 1: Upload file"),
                                 fileInput(inputId = "userData", width = "100%", label = "Select XLS/XLSX file", accept = c("application/vnd.ms-excel", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", ".xls", ".xlsx"), multiple = FALSE),
                                 htmlOutput("errorBarFile"),
                                 shinyjs::hidden(
                                   div(
                                     id = "onDataUploaded",
                                     h3("Step 2: Load data"),
                                     shinyFeedback::loadingButton(
                                       "dataCheck",
                                       class = "btn btn-default",
                                       label = "Load & Check data",
                                       style = "width:100%",
                                       loadingLabel = "Checking\U2026"
                                     )
                                   )
                                 )
                          )
                        ),
                        br(),
                        htmlOutput(outputId = "batchInfo"),
                        tags$style(type = "text/css", "#batchInfo {white-space: pre-wrap; line-height:18px;}"),
                        br(),
                        shinyjs::hidden(
                          div(
                            id = "onDataChecked",
                            tabsetPanel(
                              tabPanel(title = "Experiment parameters",
                                       id = "selectorPlating",
                                       reactable::reactableOutput("platingTable")
                              ),
                              tabPanel(title = "Counts on selective medium",
                                       id = "selectorSel",
                                       reactable::reactableOutput("selectiveTable")
                              ),
                              tabPanel(title = "Counts on non-selective medium",
                                       id = "selectorNsel",
                                       reactable::reactableOutput("nonselectiveTable")
                              )
                            ),
                            shinyjs::hidden(
                              div(
                                id = "ifDataCorrect",
                                fluidRow(
                                  column(4,
                                         offset = 4,
                                         h3("Step 3: Select options"),
                                         uiOutput("BatchOptionsPval", inline = FALSE),
                                         br()
                                  ),
                                  column(8,
                                         offset = 2,
                                         tagList(
                                           div(class = "panel panel-info",
                                               div(class="panel-heading", h5(HTML("<font color = '#3a87ad'>mlemur tried to select the additional parameters based on your data. You can change your choice for each strain using checkboxes below:</font>"))),
                                               div(class="panel-body",
                                                   style="padding-bottom: 0; overflow-y: scroll;",
                                                   div(style = "display:block;vertical-align:top;width:100%;",
                                                       div(style = "display:inline-block;vertical-align:top;width:calc(100% - 390px);", HTML("")),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;text-align:center;height:40px;", HTML("Fitness")),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;text-align:center;height:40px;", HTML("Phen. lag")),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;text-align:center;height:40px;", HTML("Resid. mut.")),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;text-align:center;height:40px;", HTML("Death rate")),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;text-align:center;height:40px;", HTML("Inoculum size")),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;text-align:center;height:40px;", HTML("Coef. of var."))),
                                                   div(style = "display:block;vertical-align:top;width:100%;background:#f5f5f5;padding-top:15px;margin-bottom:15px;",
                                                       div(style = "display:inline-block;vertical-align:top;width:calc(100% - 390px);padding-bottom:15px;", HTML("<b>Toggle all</b>")),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;padding-left:20px;",
                                                           shinyWidgets::awesomeCheckbox(inputId = "choices_AllStrains_Fitness",
                                                                                         label = "",
                                                                                         status = "primary",
                                                                                         value = FALSE)),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;padding-left:20px;",
                                                           shinyWidgets::awesomeCheckbox(inputId = "choices_AllStrains_Lag",
                                                                                         label = "",
                                                                                         status = "primary",
                                                                                         value = FALSE)),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;padding-left:20px;",
                                                           shinyWidgets::awesomeCheckbox(inputId = "choices_AllStrains_Residual",
                                                                                         label = "",
                                                                                         status = "primary",
                                                                                         value = FALSE)),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;padding-left:20px;",
                                                           shinyWidgets::awesomeCheckbox(inputId = "choices_AllStrains_Death",
                                                                                         label = "",
                                                                                         status = "primary",
                                                                                         value = FALSE)),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;padding-left:20px;",
                                                           shinyWidgets::awesomeCheckbox(inputId = "choices_AllStrains_Inoculum",
                                                                                         label = "",
                                                                                         status = "primary",
                                                                                         value = FALSE)),
                                                       div(style = "display:inline-block;vertical-align:top;width:60px;padding-left:20px;",
                                                           shinyWidgets::awesomeCheckbox(inputId = "choices_AllStrains_CV",
                                                                                         label = "",
                                                                                         status = "primary",
                                                                                         value = FALSE)),
                                                   )),
                                               div(class="panel-body",
                                                   style="max-height: 250px; overflow-y: scroll; padding-top: 0; padding-bottom: 0;",
                                                   uiOutput("BatchOptionsChoices", inline = FALSE))))
                                  ),
                                  column(4,
                                         offset = 4,
                                         h3("Step 4: Calculate"),
                                         shinyFeedback::loadingButton(
                                           "calculate4",
                                           label = HTML("<i class = \"fas fa-calculator\"></i> Calculate"),
                                           style = "width:100%;",
                                           loadingLabel = "Calculating\U2026"
                                         ),
                                         br()
                                  )
                                )
                              )
                            ),
                            shinyjs::hidden(
                              div(
                                id = "onOutputGenerated",
                                br(),
                                textOutput("infoBatch"),
                                HTML("<h4>Mutation Rates</h4>"),
                                reactable::reactableOutput("batchRateResults"),
                                br(),
                                shinyjs::hidden(
                                  div(
                                    id = "BatchPvalue",
                                    HTML("<h4><i>P</i> values</h4>"),
                                    reactable::reactableOutput("BatchPvalueResults"),
                                    br()
                                  )
                                ),
                                shinyjs::hidden(
                                  div(
                                    id = "BatchEffSize",
                                    HTML("<h4>Effect sizes</h4>"),
                                    reactable::reactableOutput("BatchEffSizeResults"),
                                    br()
                                  )
                                ),
                                fluidRow(
                                  column(4,
                                         offset = 4,
                                         h3("Step 5: Download results"),
                                         div(downloadButton("dlBatch", "Export the results to XLSX"), style =
                                               "text-align:center; width:100%")
                                  )
                                )
                              )
                            )
                          )
                        ),
                        br()
               ),
               #### Fold Tab ####
               tabPanel(
                 title = HTML("Fold"),
                 titlePanel(div(HTML(
                   "Confidence intervals for fold or other function of data"
                 ))),
                 hr(),
                 fluidRow(
                   column(
                     7,
                     fluidRow(
                       column(
                         9,
                         selectInput(inputId = "FoldUsePreset", label = "Use one of the pre-defined functions",
                                     choices = c(list("-" = 0,
                                                      "Strain 1 / Strain 2" = 1,
                                                      "Strain 1 - Strain 2" = 2,
                                                      "(Strain 1 / Strain 2) / (Strain 3 / Strain 4)" = 3,
                                                      "(Strain 1 - Strain 3) / (Strain 2 - Strain 3)" = 4,
                                                      "(Strain 1 + Strain 2) / Strain 3" = 5)),
                                     width = "100%", selected = 1)
                       ),
                       column(
                         3,
                         HTML('<label class = "control-label"></label>'),
                         actionButton("FoldClear", label = "Clear equation", icon = icon("trash", lib = "font-awesome"), style = "height: 100%;")
                       )
                     ),
                     br(),
                     div(
                       class = "panel-body bg-alt",
                       fluidRow(
                         column(
                           12,
                           shinyjqui::orderInput(inputId = 'FoldEquation', label = "Or drag and drop items to the equation field", items = c("Strain 1", "/", "Strain 2"), item_class = "danger", style = "min-height: 60px; width: 100%; padding: 8px 12px; font-size: 14px; line-height: 1.42857143; color: #3D441E; background-color: #f9f9f9; background-image: none; border: 1px solid #cccccc; border-radius: 10px; box-shadow: inset 0 3px 10px rgba(0, 0, 0, 0.075); -webkit-box-shadow: inset 0 3px 10px rgba(0, 0, 0, 0.075); overflow-x: auto; overflow-y: hidden; white-space: nowrap;", placeholder = 'Drag items here...'),
                         )
                       ),
                       shinyjqui::orderInput(inputId = 'FoldStrains', label = 'Strains', items = sapply(1:6, function(x) paste("Strain", x)),
                                             as_source = TRUE, connect = 'FoldEquation', item_class = "danger", width = "100%"),
                       shinyjqui::orderInput(inputId = 'FoldOperators', label = 'Operators', items = c("+", "-", "*", "/", "(", ")"),
                                             as_source = TRUE, connect = 'FoldEquation',  item_class = "danger", width = "100%")
                     ),
                     br(),
                     shinyjs::hidden(
                       sliderInput(inputId = "FoldNumberOfStrains", label = "Number of strains", min = 2, max = 6, value = 2, step = 1)
                     ),
                     hr(),
                     settingsPlatingUI("SettingsFold", c(0, 1, 2)),
                     hr(),
                     fluidRow(
                       style = "display: flex!important; overflow-x: auto;",
                       column(
                         6,
                         style = "min-width: 300px;",
                         h4(HTML("<b>Strain 1</b>")),
                         countsPlatingUI("CountsStrain1Fold", "SettingsFold", TRUE, usePreset = 1, defaultSettings = c(0, 1, 2)),
                       ),
                       column(
                         6,
                         style = "min-width: 300px;",
                         h4(HTML("<b>Strain 2</b>")),
                         countsPlatingUI("CountsStrain2Fold", "SettingsFold", TRUE, usePreset = 2, defaultSettings = c(0, 1, 2)),
                       ),
                       shinyjs::hidden(
                         column(
                           6,
                           id = "FoldColumn3",
                           style = "min-width: 300px;",
                           h4(HTML("<b>Strain 3</b>")),
                           countsPlatingUI("CountsStrain3Fold", "SettingsFold", TRUE, usePreset = 1, defaultSettings = c(0, 1, 2)),
                         )
                       ),
                       shinyjs::hidden(
                         column(
                           6,
                           id = "FoldColumn4",
                           style = "min-width: 300px;",
                           h4(HTML("<b>Strain 4</b>")),
                           countsPlatingUI("CountsStrain4Fold", "SettingsFold", TRUE, usePreset = 2, defaultSettings = c(0, 1, 2)),
                         )
                       ),
                       shinyjs::hidden(
                         column(
                           6,
                           id = "FoldColumn5",
                           style = "min-width: 300px;",
                           h4(HTML("<b>Strain 5</b>")),
                           countsPlatingUI("CountsStrain5Fold", "SettingsFold", TRUE, usePreset = 2, defaultSettings = c(0, 1, 2)),
                         )
                       ),
                       shinyjs::hidden(
                         column(
                           6,
                           id = "FoldColumn6",
                           style = "min-width: 300px;",
                           h4(HTML("<b>Strain 6</b>")),
                           countsPlatingUI("CountsStrain6Fold", "SettingsFold", TRUE, usePreset = 2, defaultSettings = c(0, 1, 2)),
                         )
                       )
                     ),
                     br(),
                     fluidRow(
                       column(
                         4,
                         shinyFeedback::loadingButton(
                           "calculate5",
                           label = HTML("<i class = \"fas fa-calculator\"></i> Calculate"),
                           style = "width:100%;",
                           loadingLabel = "Calculating\U2026"
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("erase5",
                                      label = "Clear all cells",
                                      width = "100%",
                                      icon = icon("trash", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("sample5",
                                      label = "Load sample",
                                      width = "100%",
                                      icon = icon("sync", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       )
                     )),
                   column(5,
                          h3("Results"),
                          htmlOutput(outputId = "errorBarFold"),
                          tags$style(type = "text/css", "#errorBarFold {white-space: pre-wrap; line-height:18px;}"),
                          shinyjs::hidden(
                            div(
                              id = "advanced5",
                              uiOutput("funFold"),
                              reactable::reactableOutput("tableFold"),
                              br(),
                              fluidRow(
                                column(8,
                                       offset = 2,
                                       div(actionButton(inputId = "clipFold",
                                                        label = "Copy to Clipboard",
                                                        style = "text-align:center",
                                                        icon = icon("clipboard"),
                                                        width = "100%"))
                                )
                              )
                            )
                          )
                   )
                 )
               ),
               #### Power Tab ####
               tabPanel(
                 title = HTML("Power"),
                 titlePanel(div(HTML(
                   "Power analysis and sample size determination"
                 ))),
                 hr(),
                 fluidRow(
                   column(
                     7,
                     shinyWidgets::awesomeRadio(inputId = "sampleSizeOrPower", label = HTML(paste("Calculate sample size or power:", infoTooltip("Choose whether you want to calculate the required sample size to achieve specified power, or statistical power given specified sample sizes."))), choices = c("Sample size" = 1, "Power" = 0), selected = 1, inline = T),
                     conditionalPanel(
                       condition = "input.sampleSizeOrPower == 1",
                       textInput(inputId = "PowerValue", label = HTML(paste("Power:", infoTooltip("Please provide a decimal number bigger than 0 but smaller than 1."))), value = "0.8")
                     ),
                     fluidRow(
                       column(
                         6,
                         h4(HTML("<b>Strain 1</b>")),
                         PowerModuleUI("PowerStrain1", usePreset = 1),
                       ),
                       column(
                         6,
                         h4(HTML("<b>Strain 2</b>")),
                         PowerModuleUI("PowerStrain2", usePreset = 2),
                       )
                     ),
                     br(),
                     fluidRow(
                       column(
                         4,
                         shinyFeedback::loadingButton(
                           "calculate6",
                           label = HTML("<i class = \"fas fa-calculator\"></i> Calculate"),
                           style = "width:100%;",
                           loadingLabel = "Calculating\U2026"
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("erase6",
                                      label = "Clear all cells",
                                      width = "100%",
                                      icon = icon("trash", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("sample6",
                                      label = "Load sample",
                                      width = "100%",
                                      icon = icon("sync", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       )
                     )
                   ),
                   column(5,
                          h3("Results"),
                          htmlOutput(outputId = "errorBarPower"),
                          tags$style(type = "text/css", "#errorBarPower {white-space: pre-wrap; line-height:18px;}"),
                          shinyjs::hidden(
                            div(
                              id = "advanced6",
                              reactable::reactableOutput("tablePower"),
                              hr(),
                              htmlOutput("powerinfo"),
                              br(),
                              hr(),
                              abbrevsUI("abbrevspval", TRUE)
                            )
                          ))
                 )
               )
    )
  )
}