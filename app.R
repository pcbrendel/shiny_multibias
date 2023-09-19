# Simultaneous Multi-bias Analysis Shiny App
# Paul Brendel

library(shiny)
library(shinythemes)
library(multibias)
library(tidyverse)
library(doParallel)

# constants ----
param_style <- "display: inline-block;vertical-align:top; width: 120px;"

# UI ----
ui <- navbarPage(
  "Simultaneous Multi-bias Analysis",
  theme = shinytheme("flatly"),
  tabPanel(
    "Instructions",
    tags$style(
      HTML(".section-header {
              font-size: 24px; /* Adjust the font size as needed */
              font-weight: bold;
              margin-bottom: 10px; /* Add spacing below the header */
            }")
    ),
    HTML("<p class='section-header'>Instructions</p>"),
    tags$b("Simultaneously adjust an effect estimate for
           any combination of uncontrolled confounding, exposure 
           misclassification, and selection bias."),
    br(),
    br(),
    "The application is organized into 3 steps:",
    tags$ol(
      tags$li("Upload your data as a text or csv file."),
      tags$li("Identify the variables in your data corresponding to the 
              exposure, outcome, and confounder(s)."),
      tags$li("Identify the biases that are present and input the quantities of 
              the appropriate bias parameters. These parameters will correspond 
              to the bias models that appear at the top of the screen.")
    ),
    "After providing these inputs, press the action button to display the
     bias-adjusted exposure-outcome odds ratio, confidence interval, and
      histogram of results.",
    br(),
    br(),
    "Please note the following:",
    tags$ul(
      tags$li("Exposure and outcome variables must both be binary."),
      tags$li("The selected number of observed confounders cannot exceed three. 
              If this restriction is problematic, consider creating a 
              propensity score.")
    ),
    "Notation:",
    tags$ul(
      tags$li("X = exposure"),
      tags$li("Y = outcome"),
      tags$li("C = known confounder(s)"),
      tags$li("U = unknown confounder"),
      tags$li("X* = misclassified exposure"),
      tags$li("S = selection")
    ),
    br(),
    HTML("<p class='section-header'>Citations</p>"),
    "Please use the following citation when using this website:",
    br(),
    tags$ul(
      tags$li(
        HTML(
          "Brendel PB, Torres AZ, Arah OA, Simultaneous adjustment of
           uncontrolled confounding, selection bias and misclassification <br> 
           in multiple-bias modelling, International Journal of Epidemiology,
           Volume 52, Issue 4, Pages 1220–1230 
           <a href='https://doi.org/10.1093/ije/dyad001' target='_blank'>Link</a>"
        )
      )
    ),
    br(),
    tags$b("Created by Paul Brendel with guidance from Onyebuchi Arah."),
    br(),
    HTML("Please submit any bug reports by opening an issue on <a href='https://github.com/pcbrendel/shiny_multibias/issues' target='_blank'>Github</a>."),
    br(),
    br(),
    h5(
      "Built with",
      img(src = "https://www.rstudio.com/wp-content/uploads/2014/04/shiny.png",
          height = "30px"),
      "by",
      img(src = paste0("https://www.rstudio.com/wp-content/uploads/2018/10/",
                       "RStudio-Logo.png"),
          height = "30px"),
      "."
    )
  ),
  tabPanel(
    "Upload Data",
    sidebarLayout(
      sidebarPanel(
        fileInput(
          "file",
          "Choose File to Upload",
          accept = c("text/csv",
                     "text/comma-separated-values,text/plain",
                     ".csv")
        ),
        checkboxInput("header", "Header", TRUE),
        radioButtons(
          "sep",
          "Separator",
          choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
          selected = ","
        ),
        radioButtons("quote", "Quote",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = '"'),
        radioButtons("disp", "Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head")
      ),
      mainPanel(tableOutput("contents"))
    )
  ),
  tabPanel(
    "Identify Variables",
    sidebarLayout(
      sidebarPanel(
        uiOutput("outcome_select"),
        uiOutput("exposure_select"),
        uiOutput("confounder_select")
      ),
      mainPanel(
        h4("Observed Outcome Regression Model:"),
        verbatimTextOutput("model"),
        htmlOutput("observed")
      )
    )
  ),
  tabPanel(
    "Bias Parameters & Results",
    sidebarLayout(
      sidebarPanel(
        checkboxInput("prob_bias_param", "Probabilistic Bias Parameters?"),
        radioButtons(
          "method",
          "Which Biases?",
          choices = list(
            "Uncontrolled Confounding and Selection Bias" = 1,
            "Uncontrolled Confounding and Exposure Misclassification" = 2,
            "Exposure Misclassification and Selection Bias" = 3,
            "All Three Biases" = 4
          ),
          selected = 1
        ),
        tags$head(
          tags$style(
            type = "text/css",
            "#inline label{ display: table-cell; text-align: center; vertical-align: middle; }
             #inline .form-group { display: table-row;}
             #inline {width: 125px;}"
          )
        ),
        uiOutput("parameters_select")
      ),
      mainPanel(
        h4("Bias Models:"),
        htmlOutput("bias_models"),
        hr(),
        htmlOutput("biased"),
        hr(),
        wellPanel(
          conditionalPanel(
            condition = "input.method == '1'",
            htmlOutput("final1a"),
            htmlOutput("final1b"),
            br(),
            plotOutput("final1c")
          ),
          conditionalPanel(
            condition = "input.method == '2'",
            htmlOutput("final2a"),
            htmlOutput("final2b"),
            br(),
            plotOutput("final2c")
          ),
          conditionalPanel(
            condition = "input.method == '3'",
            htmlOutput("final3a"),
            htmlOutput("final3b"),
            br(),
            plotOutput("final3c")
          ),
          conditionalPanel(
            condition = "input.method == '4'",
            htmlOutput("final4a"),
            htmlOutput("final4b"),
            br(),
            plotOutput("final4c")
          )
        )
      )
    )
  )
)

# SERVER ----

server <- function(input, output) {

  infile <- reactive({
    if (is.null(input$file)) {
      return(NULL)
    } else {
      input$file
    }
  })

  data <- reactive({
    if (is.null(infile())) {
      return(NULL)
    } else {
      read.csv(infile()$datapath)
    }
  })

  output$contents <- renderTable({
    req(input$file)
    df <- read.csv(input$file$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    if (input$disp == "head") {
      return(head(df))
    } else {
      return(df)
    }
  })

  output$exposure_select <- renderUI({
    selectInput("x_select", "Select Binary Exposure Variable (X)",
                choices = as.list(names(data())),
                multiple = FALSE, selected = names(data())[1])
  })

  output$outcome_select <- renderUI({
    selectInput("y_select", "Select Binary Outcome Variable (Y)",
                choices = as.list(names(data())),
                multiple = FALSE, selected = names(data())[2])
  })

  c_names <- reactive({
    drop_vars <- c(input$x_select, input$y_select)
    as.list(names(select(data(), -one_of(drop_vars))))
  })

  output$confounder_select <- renderUI({
    req(input$x_select)
    checkboxGroupInput("c_select", "Select Confounder(s) (C)",
                       choices = c_names())
  })

  form <- reactive({
    ifelse(is.null(input$c_select),
      sprintf("%s~%s", input$y_select, input$x_select),
      sprintf(
        "%s~%s",
        input$y_select,
        paste(
          input$x_select,
          paste0(input$c_select, collapse = "+"),
          sep = "+"
        )
      )
    )
  })

  logreg <- reactive({
    glm(
      as.formula(form()),
      family = binomial(link = "logit"),
      data = data()
    )
  })

  biased_or <- reactive({round(exp(coef(logreg())[2]), 2)})
  biased_or_ci <- reactive(
    c(round(exp(summary(logreg())$coef[2, 1] +
                  summary(logreg())$coef[2, 2] * qnorm(.025)), 2),
      round(exp(summary(logreg())$coef[2, 1] +
                  summary(logreg())$coef[2, 2] * qnorm(.975)), 2))
  )

  output$model <- renderPrint({ # if statement with blank
    req(input$file)
    if (is.null(input$x_select))
      return(NULL)
    print(form())
    print(summary(logreg()))
  })

  output$observed <- renderUI({
    req(input$file)
    if (is.null(input$x_select))
      return(NULL)
    HTML(paste(
      "<b>Observed OR<sub>YX</sub>:</b> ",
      biased_or(),
      ", <b>95% CI:</b> (",
      biased_or_ci()[1],
      ",",
      biased_or_ci()[2],
      ")",
      sep = ""
    ))
  })

  output$biased <- renderUI({
    req(input$file)
    HTML(paste(
      "<b>Observed OR<sub>YX</sub>:</b> ",
      biased_or(),
      ", <b>95% CI:</b> (",
      biased_or_ci()[1],
      ",",
      biased_or_ci()[2],
      ")",
      sep = ""
    ))
  })

  output$parameters_select <- renderUI({

    c_length <- length(input$c_select)

    # UC, SEL ----

    if (input$prob_bias_param == FALSE & input$method == 1) {
      return(
        list(
          sliderInput("level", "Confidence level",
                      min = 0, max = 1, value = .95),
          tags$div(
            id = "inline",
            numericInput(
              "u1_0",
              HTML("&alpha;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "u1_xstar",
              HTML("&alpha;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "u1_y",
              HTML("&alpha;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          if (c_length == 1) {
            tags$div(
              id = "inline",
              numericInput(
                "u1_c",
                HTML("&alpha;<sub>3</sub>:"),
                value = 1,
                step = .1
              )
            )
          },
          if (c_length == 2) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "u1_c1",
                  HTML("&alpha;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "u1_c2",
                  HTML("&alpha;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          if (c_length == 3) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "u1_c1",
                  HTML("&alpha;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "u1_c2",
                  HTML("&alpha;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "u1_c3",
                  HTML("&alpha;<sub>5</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          tags$div(
            id = "inline",
            numericInput(
              "s1_0",
              HTML("&beta;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "s1_xstar",
              HTML("&beta;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "s1_y",
              HTML("&beta;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          br(),
          numericInput("bs_samples", "Number of bootstrap samples:",
                       value = 1000, min = 1, max = 1000),
          br(),
          actionButton("go", "Adjust for biases")
        )
      )
    }

    if (input$prob_bias_param == TRUE & input$method == 1) {
      return(
        list(
          sliderInput("level", "Confidence level",
                      min = 0, max = 1, value = .95),
          selectInput(
            "distribution",
            HTML("Distribution:"),
            choices = list("Normal", "Uniform"),
            selected = input$distribution
          ),

          if (input$distribution == "Normal" && !is.null(input$distribution)
              && !is.null(input$prob_bias_param) && !is.null(input$method)) {
            list(
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_mean",
                  HTML("&alpha;<sub>0</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_sd",
                  HTML("&alpha;<sub>0</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_xstar_mean",
                  HTML("&alpha;<sub>1</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_xstar_sd",
                  HTML("&alpha;<sub>1</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_mean",
                  HTML("&alpha;<sub>2</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_sd",
                  HTML("&alpha;<sub>2</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c_mean",
                      HTML("&alpha;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c_sd",
                      HTML("&alpha;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c1_mean",
                      HTML("&alpha;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c1_sd",
                      HTML("&alpha;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c2_mean",
                      HTML("&alpha;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c2_sd",
                      HTML("&alpha;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c1_mean",
                      HTML("&alpha;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c1_sd",
                      HTML("&alpha;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c2_mean",
                      HTML("&alpha;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c2_sd",
                      HTML("&alpha;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c3_mean",
                      HTML("&alpha;<sub>5</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c3_sd",
                      HTML("&alpha;<sub>5</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_mean",
                  HTML("&beta;<sub>0</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_sd",
                  HTML("&beta;<sub>0</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_mean",
                  HTML("&beta;<sub>1</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_sd",
                  HTML("&beta;<sub>1</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_mean",
                  HTML("&beta;<sub>2</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_sd",
                  HTML("&beta;<sub>2</sub> sd:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },

          if (input$distribution == "Uniform" && !is.null(input$distribution) &&
                !is.null(input$prob_bias_param) && !is.null(input$method)) {
            list(
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_min",
                  HTML("&alpha;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_max",
                  HTML("&alpha;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_xstar_min",
                  HTML("&alpha;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_xstar_max",
                  HTML("&alpha;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_min",
                  HTML("&alpha;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_max",
                  HTML("&alpha;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c_min",
                      HTML("&alpha;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c_max",
                      HTML("&alpha;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c1_min",
                      HTML("&alpha;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c1_max",
                      HTML("&alpha;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c2_min",
                      HTML("&alpha;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c2_max",
                      HTML("&alpha;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c1_min",
                      HTML("&alpha;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c1_max",
                      HTML("&alpha;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c2_min",
                      HTML("&alpha;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c2_max",
                      HTML("&alpha;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c3_min",
                      HTML("&alpha;<sub>5</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "u1_c3_max",
                      HTML("&alpha;<sub>5</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_min",
                  HTML("&beta;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_max",
                  HTML("&beta;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_min",
                  HTML("&beta;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_max",
                  HTML("&beta;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_min",
                  HTML("&beta;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_max",
                  HTML("&beta;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },

          br(),
          numericInput("bs_samples", "Number of bootstrap samples:",
                       value = 1000, min = 1, max = 1000),
          br(),
          actionButton("go", "Adjust for biases")
        )
      )
    }

    # UC, EMC ----

    if (input$prob_bias_param == FALSE & input$method == 2) {
      return(
        list(
          sliderInput("level", "Confidence level",
                      min = 0, max = 1, value = .95),
          tags$div(
            id = "inline",
            numericInput(
              "u1_0",
              HTML("&alpha;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "u1_x",
              HTML("&alpha;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "u1_y",
              HTML("&alpha;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "x1_0",
              HTML("&delta;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "x1_xstar",
              HTML("&delta;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "x1_y",
              HTML("&delta;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          if (c_length == 1) {
            tags$div(
              id = "inline",
              numericInput(
                "x1_c",
                HTML("&delta;<sub>3</sub>:"),
                value = 1,
                step = .1
              )
            )
          },
          if (c_length == 2) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c1",
                  HTML("&delta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c2",
                  HTML("&delta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          if (c_length == 3) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c1",
                  HTML("&delta;<sub>3</sub>:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c2",
                  HTML("&delta;<sub>4</sub>:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c3",
                  HTML("&delta;<sub>5</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          br(),
          numericInput("bs_samples", "Number of bootstrap samples:",
                       value = 1000, min = 1, max = 1000),
          br(),
          actionButton("go", "Adjust for biases")
        )
      )
    }
    if (input$prob_bias_param == TRUE & input$method == 2) {
      return(
        list(
          sliderInput("level", "Confidence level",
                      min = 0, max = 1, value = .95),
          selectInput(
            "distribution",
            HTML("Distribution:"),
            choices = list("Normal", "Uniform"),
            selected = input$distribution
          ),
          if (input$distribution == "Normal" && !is.null(input$distribution)
              && !is.null(input$prob_bias_param) && !is.null(input$method)) {
            list(
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_mean",
                  HTML("&alpha;<sub>0</sub> mean:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_sd",
                  HTML("&alpha;<sub>0</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_x_mean",
                  HTML("&alpha;<sub>1</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_x_sd",
                  HTML("&alpha;<sub>1</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_mean",
                  HTML("&alpha;<sub>2</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_sd",
                  HTML("&alpha;<sub>2</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_mean",
                  HTML("&delta;<sub>0</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_sd",
                  HTML("&delta;<sub>0</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_mean",
                  HTML("&delta;<sub>1</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_sd",
                  HTML("&delta;<sub>1</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_mean",
                  HTML("&delta;<sub>2</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_sd",
                  HTML("&delta;<sub>2</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_mean",
                      HTML("&delta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_sd",
                      HTML("&delta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_mean",
                      HTML("&delta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_sd",
                      HTML("&delta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_mean",
                      HTML("&delta;<sub>5</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_sd",
                      HTML("&delta;<sub>5</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              }
            )
          },

          if (input$distribution == "Uniform" && !is.null(input$distribution)
              && !is.null(input$prob_bias_param) && !is.null(input$method)) {
            list(
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_min",
                  HTML("&alpha;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_max",
                  HTML("&alpha;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_x_min",
                  HTML("&alpha;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_x_max",
                  HTML("&alpha;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_min",
                  HTML("&alpha;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_max",
                  HTML("&alpha;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_min",
                  HTML("&delta;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_max",
                  HTML("&delta;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_min",
                  HTML("&delta;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_max",
                  HTML("&delta;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_min",
                  HTML("&delta;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_max",
                  HTML("&delta;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_min",
                      HTML("&delta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_max",
                      HTML("&delta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_min",
                      HTML("&delta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_max",
                      HTML("&delta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_min",
                      HTML("&delta;<sub>5</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_max",
                      HTML("&delta;<sub>5</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              }
            )
          },
          br(),
          numericInput("bs_samples", "Number of bootstrap samples:",
                       value = 1000, min = 1, max = 1000),
          br(),
          actionButton("go", "Adjust for biases")
        )
      )
    }

    # EMC, Sel ----

    if (input$prob_bias_param == FALSE & input$method == 3) {
      return(
        list(
          sliderInput("level", "Confidence level",
                      min = 0, max = 1, value = .95),
          tags$div(
            id = "inline",
            numericInput(
              "x1_0",
              HTML("&delta;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "x1_xstar",
              HTML("&delta;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "x1_y",
              HTML("&delta;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          if (c_length == 1) {
            tags$div(
              id = "inline",
              numericInput(
                "x1_c",
                HTML("&delta;<sub>3</sub>:"),
                value = 1,
                step = .1
              )
            )
          },
          if (c_length == 2) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c1",
                  HTML("&delta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c2",
                  HTML("&delta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          if (c_length == 3) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c1",
                  HTML("&delta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c2",
                  HTML("&delta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c3",
                  HTML("&delta;<sub>5</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          tags$div(
            id = "inline",
            numericInput(
              "s1_0",
              HTML("&beta;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "s1_xstar",
              HTML("&beta;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "s1_y",
              HTML("&beta;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          if (c_length == 1) {
            tags$div(
              id = "inline",
              numericInput(
                "s1_c",
                HTML("&beta;<sub>3</sub>:"),
                value = 1,
                step = .1
              )
            )
          },
          if (c_length == 2) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c1",
                  HTML("&beta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c2",
                  HTML("&beta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          if (c_length == 3) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c1",
                  HTML("&beta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c2",
                  HTML("&beta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c3",
                  HTML("&beta;<sub>5</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          br(),
          numericInput("bs_samples", "Number of bootstrap samples:",
                       value = 1000, min = 1, max = 1000),
          br(),
          actionButton("go", "Adjust for biases")
        )
      )
    }

    if (input$prob_bias_param == TRUE & input$method == 3) {
      return(
        list(
          sliderInput("level", "Confidence level",
                      min = 0, max = 1, value = .95),
          selectInput("distribution", HTML("Distribution:"),
                      choices = list("Normal", "Uniform"),
                      selected = input$distribution),
          if (input$distribution == "Normal" && !is.null(input$distribution) &&
                !is.null(input$prob_bias_param) && !is.null(input$method)) {
            list(
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_mean",
                  HTML("&delta;<sub>0</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_sd",
                  HTML("&delta;<sub>0</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_mean",
                  HTML("&delta;<sub>1</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_sd",
                  HTML("&delta;<sub>1</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_mean",
                  HTML("&delta;<sub>2</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_sd",
                  HTML("&delta;<sub>2</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_mean",
                      HTML("&delta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_sd",
                      HTML("&delta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_mean",
                      HTML("&delta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_sd",
                      HTML("&delta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_mean",
                      HTML("&delta;<sub>5</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_sd",
                      HTML("&delta;<sub>5</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_mean",
                  HTML("&beta;<sub>0</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_sd",
                  HTML("&beta;<sub>0</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_mean",
                  HTML("&beta;<sub>1</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_sd",
                  HTML("&beta;<sub>1</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_mean",
                  HTML("&beta;<sub>2</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_sd",
                  HTML("&beta;<sub>2</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c_mean",
                      HTML("&beta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c_sd",
                      HTML("&beta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_mean",
                      HTML("&beta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_sd",
                      HTML("&beta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_mean",
                      HTML("&beta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_sd",
                      HTML("&beta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_mean",
                      HTML("&beta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_sd",
                      HTML("&beta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_mean",
                      HTML("&beta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_sd", HTML("&beta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c3_mean",
                      HTML("&beta;<sub>5</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c3_sd",
                      HTML("&beta;<sub>5</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              }
            )
          },
          if (input$distribution == "Uniform" && !is.null(input$distribution) &&
                !is.null(input$prob_bias_param) && !is.null(input$method)) {
            list(
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_min",
                  HTML("&delta;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_max",
                  HTML("&delta;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_min",
                  HTML("&delta;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_max",
                  HTML("&delta;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_min",
                  HTML("&delta;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_max",
                  HTML("&delta;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_min",
                      HTML("&delta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_max",
                      HTML("&delta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_min",
                      HTML("&delta;<sub>4</sub> min:"),
                      value = 1, step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_max",
                      HTML("&delta;<sub>4</sub> max:"),
                      value = 1, step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_min",
                      HTML("&delta;<sub>5</sub> min:"),
                      value = 1, step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_max",
                      HTML("&delta;<sub>5</sub> max:"),
                      value = 1, step = .1
                    )
                  )
                )
              },
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_min",
                  HTML("&beta;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_max",
                  HTML("&beta;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_min",
                  HTML("&beta;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_max",
                  HTML("&beta;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_min",
                  HTML("&beta;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_max",
                  HTML("&beta;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c_min",
                      HTML("&beta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c_max",
                      HTML("&beta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_min",
                      HTML("&beta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_max",
                      HTML("&beta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_min",
                      HTML("&beta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_max",
                      HTML("&beta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_min",
                      HTML("&beta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_max",
                      HTML("&beta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_min",
                      HTML("&beta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_max",
                      HTML("&beta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c3_min",
                      HTML("&beta;<sub>5</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c3_max",
                      HTML("&beta;<sub>5</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              }
            )
          },
          br(),
          numericInput("bs_samples", "Number of bootstrap samples:",
                       value = 1000, min = 1, max = 1000),
          br(),
          actionButton("go", "Adjust for biases")
        )
      )
    }

    # UC, EMC, Sel ----

    if (input$prob_bias_param == FALSE & input$method == 4) {
      return(
        list(
          sliderInput("level", "Confidence level",
                      min = 0, max = 1, value = .95),
          tags$div(
            id = "inline",
            numericInput(
              "u1_0",
              HTML("&alpha;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "u1_x",
              HTML("&alpha;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "u1_y",
              HTML("&alpha;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "x1_0",
              HTML("&delta;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "x1_xstar",
              HTML("&delta;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "x1_y",
              HTML("&delta;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          if (c_length == 1) {
            tags$div(
              id = "inline",
              numericInput(
                "x1_c",
                HTML("&delta;<sub>3</sub>:"),
                value = 1,
                step = .1
              )
            )
          },
          if (c_length == 2) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c1",
                  HTML("&delta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c2",
                  HTML("&delta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          if (c_length == 3) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c1",
                  HTML("&delta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c2",
                  HTML("&delta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "x1_c3",
                  HTML("&delta;<sub>5</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          tags$div(
            id = "inline",
            numericInput(
              "s1_0",
              HTML("&beta;<sub>0</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "s1_xstar",
              HTML("&beta;<sub>1</sub>:"),
              value = 1,
              step = .1
            )
          ),
          tags$div(
            id = "inline",
            numericInput(
              "s1_y",
              HTML("&beta;<sub>2</sub>:"),
              value = 1,
              step = .1
            )
          ),
          if (c_length == 1) {
            tags$div(
              id = "inline",
              numericInput(
                "s1_c",
                HTML("&beta;<sub>3</sub>:"),
                value = 1,
                step = .1
              )
            )
          },
          if (c_length == 2) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c1",
                  HTML("&beta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c2",
                  HTML("&beta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          if (c_length == 3) {
            list(
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c1",
                  HTML("&beta;<sub>3</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c2",
                  HTML("&beta;<sub>4</sub>:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                id = "inline",
                numericInput(
                  "s1_c3",
                  HTML("&beta;<sub>5</sub>:"),
                  value = 1,
                  step = .1
                )
              )
            )
          },
          br(),
          numericInput("bs_samples", "Number of bootstrap samples:",
                       value = 1000, min = 1, max = 1000),
          br(),
          actionButton("go", "Adjust for biases")
        )
      )
    }

    if (input$prob_bias_param == TRUE & input$method == 4) {
      return(
        list(
          sliderInput("level", "Confidence level",
                      min = 0, max = 1, value = .95),
          selectInput("distribution", HTML("Distribution:"),
                      choices = list("Normal", "Uniform"),
                      selected = input$distribution),
          if (input$distribution == "Normal" && !is.null(input$distribution)
              && !is.null(input$prob_bias_param) && !is.null(input$method)) {
            list(
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_mean",
                  HTML("&alpha;<sub>0</sub> mean:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_sd",
                  HTML("&alpha;<sub>0</sub> sd:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_x_mean",
                  HTML("&alpha;<sub>1</sub> mean:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_x_sd",
                  HTML("&alpha;<sub>1</sub> sd:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_mean",
                  HTML("&alpha;<sub>2</sub> mean:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_sd",
                  HTML("&alpha;<sub>2</sub> sd:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_mean",
                  HTML("&delta;<sub>0</sub> mean:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_sd",
                  HTML("&delta;<sub>0</sub> sd:"),
                  value = 1, step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_mean",
                  HTML("&delta;<sub>1</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_sd",
                  HTML("&delta;<sub>1</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_mean",
                  HTML("&delta;<sub>2</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_sd",
                  HTML("&delta;<sub>2</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_mean",
                      HTML("&delta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_sd",
                      HTML("&delta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_mean",
                      HTML("&delta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_sd",
                      HTML("&delta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_mean",
                      HTML("&delta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_sd",
                      HTML("&delta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_mean",
                      HTML("&delta;<sub>5</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_sd",
                      HTML("&delta;<sub>5</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_mean",
                  HTML("&beta;<sub>0</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_sd",
                  HTML("&beta;<sub>0</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_mean",
                  HTML("&beta;<sub>1</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_sd",
                  HTML("&beta;<sub>1</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_mean",
                  HTML("&beta;<sub>2</sub> mean:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_sd",
                  HTML("&beta;<sub>2</sub> sd:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c_mean",
                      HTML("&beta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c_sd",
                      HTML("&beta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_mean",
                      HTML("&beta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_sd",
                      HTML("&beta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_mean",
                      HTML("&beta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_sd",
                      HTML("&beta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_mean",
                      HTML("&beta;<sub>3</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_sd",
                      HTML("&beta;<sub>3</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_mean",
                      HTML("&beta;<sub>4</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_sd",
                      HTML("&beta;<sub>4</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c3_mean",
                      HTML("&beta;<sub>5</sub> mean:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c3_sd",
                      HTML("&beta;<sub>5</sub> sd:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              }
            )
          },

          if (input$distribution == "Uniform" && !is.null(input$distribution)
              && !is.null(input$prob_bias_param) && !is.null(input$method)) {
            list(
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_min",
                  HTML("&alpha;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_0_max",
                  HTML("&alpha;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_x_min",
                  HTML("&alpha;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_x_max",
                  HTML("&alpha;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_min",
                  HTML("&alpha;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "u1_y_max",
                  HTML("&alpha;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_min",
                  HTML("&delta;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_0_max",
                  HTML("&delta;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_min",
                  HTML("&delta;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_xstar_max",
                  HTML("&delta;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_min",
                  HTML("&delta;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "x1_y_max",
                  HTML("&delta;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1, step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_min",
                      HTML("&delta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_max",
                      HTML("&delta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_min",
                      HTML("&delta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c1_max",
                      HTML("&delta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_min",
                      HTML("&delta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c2_max",
                      HTML("&delta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_min",
                      HTML("&delta;<sub>5</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "x1_c3_max",
                      HTML("&delta;<sub>5</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_min",
                  HTML("&beta;<sub>0</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_0_max",
                  HTML("&beta;<sub>0</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_min",
                  HTML("&beta;<sub>1</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_xstar_max",
                  HTML("&beta;<sub>1</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_min",
                  HTML("&beta;<sub>2</sub> min:"),
                  value = 1,
                  step = .1
                )
              ),
              tags$div(
                style = param_style,
                numericInput(
                  "s1_y_max",
                  HTML("&beta;<sub>2</sub> max:"),
                  value = 1,
                  step = .1
                )
              ),
              if (c_length == 1) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c_min",
                      HTML("&beta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c_max",
                      HTML("&beta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 2) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_min",
                      HTML("&beta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_max",
                      HTML("&beta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_min",
                      HTML("&beta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_max",
                      HTML("&beta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              },
              if (c_length == 3) {
                list(
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_min",
                      HTML("&beta;<sub>3</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c1_max",
                      HTML("&beta;<sub>3</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_min",
                      HTML("&beta;<sub>4</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c2_max",
                      HTML("&beta;<sub>4</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c3_min",
                      HTML("&beta;<sub>5</sub> min:"),
                      value = 1,
                      step = .1
                    )
                  ),
                  tags$div(
                    style = param_style,
                    numericInput(
                      "s1_c3_max",
                      HTML("&beta;<sub>5</sub> max:"),
                      value = 1,
                      step = .1
                    )
                  )
                )
              }
            )
          },
          br(),
          numericInput("bs_samples", "Number of bootstrap samples:",
                       value = 1000, min = 1, max = 1000),
          br(),
          actionButton("go", "Adjust for biases")
        )
      )
    }
  })

  output$bias_models <- renderUI({
    j <- "where j = 1:number of confounders"
    if (input$method == 1) {
      x <- HTML("logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X +
                 &alpha;<sub>2</sub>Y + &alpha;<sub>2+j</sub>C<sub>j</sub>")
      y <- HTML("logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X +
                 &beta;<sub>2</sub>Y")
      return(HTML(paste(x, y, j, sep = "<br/>")))
    }
    if (input$method == 2) {
      x <- HTML("logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X +
                 &alpha;<sub>2</sub>Y")
      y <- HTML("logit(P(X=1) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* +
                 &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>")
      return(HTML(paste(x, y, j, sep = "<br/>")))
    }
    if (input$method == 3) {
      x <- HTML("logit(P(X=1) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* +
                 &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>")
      y <- HTML("logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* +
                 &beta;<sub>2</sub>Y + &beta;<sub>2+j</sub>C<sub>j</sub>")
      return(HTML(paste(x, y, j, sep = "<br/>")))
    }
    if (input$method == 4) {
      x <- HTML("logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X +
                 &alpha;<sub>2</sub>Y")
      y <- HTML("logit(P(X=1) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* +
                 &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>")
      z <- HTML("logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* +
                 &beta;<sub>2</sub>Y + &beta;<sub>2+j</sub>C<sub>j</sub>")
      return(HTML(paste(x, y, z, j, sep = "<br/>")))
    }
  })

  # UC, Sel ----

  # check if this is unnecessary
  pu1_vals <- reactive({
    c_length <- length(input$c_select)
    if (input$prob_bias_param == FALSE) {
      if (c_length == 0) {
        return(c(input$u1_0, input$u1_xstar, input$u1_y))
      }
      if (c_length == 1) {
        return(c(input$u1_0, input$u1_xstar, input$u1_y, input$u1_c))
      }
      if (c_length == 2) {
        return(c(input$u1_0, input$u1_xstar, input$u1_y,
                 input$u1_c1, input$u1_c2))
      }
      if (c_length == 3) {
        return(c(input$u1_0, input$u1_xstar, input$u1_y,
                 input$u1_c1, input$u1_c2, input$u1_c3))
      }
    }
  })

  est_uc_sel <- eventReactive(input$go, {

    showNotification("Results are loading. The loading time will vary
                      depending on the size of the data and the number of
                      bootstrap samples.")
    no_cores <- 3
    registerDoParallel(cores = no_cores)
    cl <- makeCluster(no_cores)
    nreps <- input$bs_samples

    if (input$prob_bias_param == FALSE) {

      est <- vector(length = nreps)
      df <- data()
      pu1 <- pu1_vals()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]
        isolate({
          est[i] <- adjust_uc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = pu1,
            s_model_coefs = c(
              input$s1_0,
              input$s1_xstar,
              input$s1_y
            ),
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }

    if (input$prob_bias_param == TRUE & input$distribution == "Normal") {

      c_length <- length(input$c_select)
      est <- vector(length = nreps)
      df <- data()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_uc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = if (c_length == 0) {
              c(rnorm(1, input$u1_0_mean, input$u1_0_sd),
                rnorm(1, input$u1_xstar_mean, input$u1_xstar_sd),
                rnorm(1, input$u1_y_mean, input$u1_y_sd)
              )
            } else if (c_length == 1) {
              c(rnorm(1, input$u1_0_mean, input$u1_0_sd),
                rnorm(1, input$u1_xstar_mean, input$u1_xstar_sd),
                rnorm(1, input$u1_y_mean, input$u1_y_sd),
                rnorm(1, input$u1_c_mean, input$u1_c_sd)
              )
            } else if (c_length == 2) {
              c(rnorm(1, input$u1_0_mean, input$u1_0_sd),
                rnorm(1, input$u1_xstar_mean, input$u1_xstar_sd),
                rnorm(1, input$u1_y_mean, input$u1_y_sd),
                rnorm(1, input$u1_c1_mean, input$u1_c1_sd),
                rnorm(1, input$u1_c2_mean, input$u1_c2_sd)
              )
            } else if (c_length == 3) {
              c(rnorm(1, input$u1_0_mean, input$u1_0_sd),
                rnorm(1, input$u1_xstar_mean, input$u1_xstar_sd),
                rnorm(1, input$u1_y_mean, input$u1_y_sd),
                rnorm(1, input$u1_c1_mean, input$u1_c1_sd),
                rnorm(1, input$u1_c2_mean, input$u1_c2_sd),
                rnorm(1, input$u1_c3_mean, input$u1_c3_sd)
              )
            },
            s_model_coefs = c(
              rnorm(1, input$s1_0_mean, input$s1_0_sd),
              rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
              rnorm(1, input$s1_y_mean, input$s1_y_sd)
            ),
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }

    if (input$prob_bias_param == TRUE & input$distribution == "Uniform") {

      c_length <- length(input$c_select)
      est <- vector(length = nreps)
      df <- data()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_uc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = if (c_length == 0) {
              c(runif(1, input$u1_0_min, input$u1_0_max),
                runif(1, input$u1_xstar_min, input$u1_xstar_max),
                runif(1, input$u1_y_min, input$u1_y_max))
            } else if (c_length == 1) {
              c(runif(1, input$u1_0_min, input$u1_0_max),
                runif(1, input$u1_xstar_min, input$u1_xstar_max),
                runif(1, input$u1_y_min, input$u1_y_max),
                runif(1, input$u1_c_min, input$u1_c_max)
              )
            } else if (c_length == 2) {
              c(runif(1, input$u1_0_min, input$u1_0_max),
                runif(1, input$u1_xstar_min, input$u1_xstar_max),
                runif(1, input$u1_y_min, input$u1_y_max),
                runif(1, input$u1_c1_min, input$u1_c1_max),
                runif(1, input$u1_c2_min, input$u1_c2_max)
              )
            } else if (c_length == 3) {
              c(runif(1, input$u1_0_min, input$u1_0_max),
                runif(1, input$u1_xstar_min, input$u1_xstar_max),
                runif(1, input$u1_y_min, input$u1_y_max),
                runif(1, input$u1_c1_min, input$u1_c1_max),
                runif(1, input$u1_c2_min, input$u1_c2_max),
                runif(1, input$u1_c3_min, input$u1_c3_max)
              )
            },
            s_model_coefs = c(
              runif(1, input$s1_0_min, input$s1_0_max),
              runif(1, input$s1_xstar_min, input$s1_xstar_max),
              runif(1, input$s1_y_min, input$s1_y_max)
            ),
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }
  }, ignoreInit = TRUE)

  output$final1a <- renderUI({
    req(input$file)
    HTML(
      paste(
        "<b>Bias-adjusted OR<sub>YX</sub>:</b>",
        round(median(est_uc_sel()), 2)
      )
    )
  })

  output$final1b <- renderUI({
    req(input$file)
    HTML(
      paste(
        "<b>Bias-adjusted OR<sub>YX</sub> confidence interval:</b> (",
        round(
          quantile(
            est_uc_sel(),
            (1 - input$level) - ((1 - input$level) / 2)
          ),
          2
        ),
        ", ",
        round(
          quantile(
            est_uc_sel(),
            input$level + ((1 - input$level) / 2)
          ),
          2
        ),
        ")",
        sep = ""
      )
    )
  })

  output$final1c <- renderPlot({
    req(input$file)
    hist(
      est_uc_sel(),
      main = "Distribution of bias-adjusted estimates from bootstrap samples",
      xlab = expression('OR'["YX"]),
      col = "red"
    )
  })

  # UC, EMC ----

  px1_vals <- reactive({
    c_length <- length(input$c_select)
    if (c_length == 0) {
      return(c(input$x1_0, input$x1_xstar, input$x1_y))
    }
    if (c_length == 1) {
      return(c(input$x1_0, input$x1_xstar, input$x1_y, input$x1_c))
    }
    if (c_length == 2) {
      return(c(input$x1_0, input$x1_xstar, input$x1_y,
               input$x1_c1, input$x1_c2))
    }
    if (c_length == 3) {
      return(c(input$x1_0, input$x1_xstar, input$x1_y,
               input$x1_c1, input$x1_c2, input$x1_c3))
    }
  })

  est_uc_emc <- eventReactive(input$go, {

    showNotification("Results are loading. The loading time will vary
                      depending on the size of the data and the number of
                      bootstrap samples.")
    no_cores <- 3
    registerDoParallel(cores = no_cores)
    cl <- makeCluster(no_cores)
    nreps <- input$bs_samples

    if (input$prob_bias_param == FALSE) {

      est <- vector(length = nreps)
      df <- data()
      px1 <- px1_vals()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]
        isolate({
          est[i] <- adjust_uc_emc(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = c(input$u1_0, input$u1_x, input$u1_y),
            x_model_coefs = px1,
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }

    if (input$prob_bias_param == TRUE & input$distribution == "Normal") {

      c_length <- length(input$c_select)
      est <- vector(length = nreps)
      df <- data()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_uc_emc(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = c(
              rnorm(1, input$u1_0_mean, input$u1_0_sd),
              rnorm(1, input$u1_x_mean, input$u1_x_sd),
              rnorm(1, input$u1_y_mean, input$u1_y_sd)
            ),
            x_model_coefs = if (c_length == 0) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd))
            } else if (c_length == 1) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c_mean, input$x1_c_sd))
            } else if (c_length == 2) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c1_mean, input$x1_c1_sd),
                rnorm(1, input$x1_c2_mean, input$x1_c2_sd))
            } else if (c_length == 3) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c1_mean, input$x1_c1_sd),
                rnorm(1, input$x1_c2_mean, input$x1_c2_sd),
                rnorm(1, input$x1_c3_mean, input$x1_c3_sd))
            },
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }

    if (input$prob_bias_param == TRUE & input$distribution == "Uniform") {

      c_length <- length(input$c_select)
      est <- vector(length = nreps)
      df <- data()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_uc_emc(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = c(
              rnorm(1, input$u1_0_min, input$u1_0_max),
              rnorm(1, input$u1_x_min, input$u1_x_max),
              rnorm(1, input$u1_y_min, input$u1_y_max)
            ),
            x_model_coefs = if (c_length == 0) {
              c(rnorm(1, input$x1_0_min, input$x1_0_max),
                rnorm(1, input$x1_xstar_min, input$x1_xstar_max),
                rnorm(1, input$x1_y_min, input$x1_y_max))
            } else if (c_length == 1) {
              c(rnorm(1, input$x1_0_min, input$x1_0_max),
                rnorm(1, input$x1_xstar_min, input$x1_xstar_max),
                rnorm(1, input$x1_y_min, input$x1_y_max),
                rnorm(1, input$x1_c_min, input$x1_c_max))
            } else if (c_length == 2) {
              c(rnorm(1, input$x1_0_min, input$x1_0_max),
                rnorm(1, input$x1_xstar_min, input$x1_xstar_max),
                rnorm(1, input$x1_y_min, input$x1_y_max),
                rnorm(1, input$x1_c1_min, input$x1_c1_max),
                rnorm(1, input$x1_c2_min, input$x1_c2_max))
            } else if (c_length == 3) {
              c(rnorm(1, input$x1_0_min, input$x1_0_max),
                rnorm(1, input$x1_xstar_min, input$x1_xstar_max),
                rnorm(1, input$x1_y_min, input$x1_y_max),
                rnorm(1, input$x1_c1_min, input$x1_c1_max),
                rnorm(1, input$x1_c2_min, input$x1_c2_max),
                rnorm(1, input$x1_c3_min, input$x1_c3_max))
            },
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }
  }, ignoreInit = TRUE)

  output$final2a <- renderUI({
    req(input$file)
    HTML(
      paste(
        "<b>Bias-adjusted OR<sub>YX</sub>:</b>",
        round(median(est_uc_emc()), 2)
      )
    )
  })

  output$final2b <- renderUI({
    req(input$file)
    HTML(
      paste(
        "<b>Bias-adjusted OR<sub>YX</sub> confidence interval:</b> (",
        round(
          quantile(
            est_uc_emc(),
            (1 - input$level) - ((1 - input$level) / 2)
          ),
          2
        ),
        ", ",
        round(
          quantile(
            est_uc_emc(),
            input$level + ((1 - input$level) / 2)
          ),
          2
        ),
        ")",
        sep = ""
      )
    )
  })

  output$final2c <- renderPlot({
    req(input$file)
    hist(
      est_uc_emc(),
      main = "Distribution of bias-adjusted estimates from bootstrap samples",
      xlab = expression('OR'["YX"]),
      col = "red"
    )
  })

  # EMC, Sel ----

  ps1_vals <- reactive({
    c_length <- length(input$c_select)
    if (c_length == 0) {
      return(c(input$s1_0, input$s1_xstar, input$s1_y))
    }
    if (c_length == 1) {
      return(c(input$s1_0, input$s1_xstar, input$s1_y, input$s1_c))
    }
    if (c_length == 2) {
      return(c(input$s1_0, input$s1_xstar, input$s1_y,
               input$s1_c1, input$s1_c2))
    }
    if (c_length == 3) {
      return(c(input$s1_0, input$s1_xstar, input$s1_y,
               input$s1_c1, input$s1_c2, input$s1_c3))
    }
  })

  est_emc_sel <- eventReactive(input$go, {

    showNotification("Results are loading. The loading time will vary
                      depending on the size of the data and the number of
                      bootstrap samples.")
    no_cores <- 3
    registerDoParallel(cores = no_cores)
    cl <- makeCluster(no_cores)
    nreps <- input$bs_samples

    if (input$prob_bias_param == FALSE) {

      est <- vector(length = nreps)
      df <- data()
      px1 <- px1_vals()
      ps1 <- ps1_vals()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]
        isolate({
          est[i] <- adjust_emc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            x_model_coefs = px1,
            s_model_coefs = ps1,
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }

    if (input$prob_bias_param == TRUE & input$distribution == "Normal") {

      c_length <- length(input$c_select)
      est <- vector(length = nreps)
      df <- data()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_emc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            x_model_coefs = if (c_length == 0) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd))
            } else if (c_length == 1) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c_mean, input$x1_c_sd))
            } else if (c_length == 2) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c1_mean, input$x1_c1_sd),
                rnorm(1, input$x1_c2_mean, input$x1_c2_sd))
            } else if (c_length == 3) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c1_mean, input$x1_c1_sd),
                rnorm(1, input$x1_c2_mean, input$x1_c2_sd),
                rnorm(1, input$x1_c3_mean, input$x1_c3_sd))
            },
            s_model_coefs = if (c_length == 0) {
              c(rnorm(1, input$s1_0_mean, input$s1_0_sd),
                rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
                rnorm(1, input$s1_y_mean, input$s1_y_sd))
            } else if (c_length == 1) {
              c(rnorm(1, input$s1_0_mean, input$s1_0_sd),
                rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
                rnorm(1, input$s1_y_mean, input$s1_y_sd),
                rnorm(1, input$s1_c_mean, input$s1_c_sd))
            } else if (c_length == 2) {
              c(rnorm(1, input$s1_0_mean, input$s1_0_sd),
                rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
                rnorm(1, input$s1_y_mean, input$s1_y_sd),
                rnorm(1, input$s1_c1_mean, input$s1_c1_sd),
                rnorm(1, input$s1_c2_mean, input$s1_c2_sd))
            } else if (c_length == 3) {
              c(rnorm(1, input$s1_0_mean, input$s1_0_sd),
                rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
                rnorm(1, input$s1_y_mean, input$s1_y_sd),
                rnorm(1, input$s1_c1_mean, input$s1_c1_sd),
                rnorm(1, input$s1_c2_mean, input$s1_c2_sd),
                rnorm(1, input$s1_c3_mean, input$s1_c3_sd))
            },
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }

    if (input$prob_bias_param == TRUE & input$distribution == "Uniform") {

      c_length <- length(input$c_select)
      est <- vector(length = nreps)
      df <- data()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_emc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            x_model_coefs = if (c_length == 0) {
              c(runif(1, input$x1_0_min, input$x1_0_max),
                runif(1, input$x1_xstar_min, input$x1_xstar_max),
                runif(1, input$x1_y_min, input$x1_y_max))
            } else if (c_length == 1) {
              c(runif(1, input$x1_0_min, input$x1_0_max),
                runif(1, input$x1_xstar_min, input$x1_xstar_max),
                runif(1, input$x1_y_min, input$x1_y_max),
                runif(1, input$x1_c_min, input$x1_c_max))
            } else if (c_length == 2) {
              c(runif(1, input$x1_0_min, input$x1_0_max),
                runif(1, input$x1_xstar_min, input$x1_xstar_max),
                runif(1, input$x1_y_min, input$x1_y_max),
                runif(1, input$x1_c1_min, input$x1_c1_max),
                runif(1, input$x1_c2_min, input$x1_c2_max))
            } else if (c_length == 3) {
              c(runif(1, input$x1_0_min, input$x1_0_max),
                runif(1, input$x1_xstar_min, input$x1_xstar_max),
                runif(1, input$x1_y_min, input$x1_y_max),
                runif(1, input$x1_c1_min, input$x1_c1_max),
                runif(1, input$x1_c2_min, input$x1_c2_max),
                runif(1, input$x1_c3_min, input$x1_c3_max))
            },
            s_model_coefs = if (c_length == 0) {
              c(runif(1, input$s1_0_min, input$s1_0_max),
                runif(1, input$s1_xstar_min, input$s1_xstar_max),
                runif(1, input$s1_y_min, input$s1_y_max))
            } else if (c_length == 1) {
              c(runif(1, input$s1_0_min, input$s1_0_max),
                runif(1, input$s1_xstar_min, input$s1_xstar_max),
                runif(1, input$s1_y_min, input$s1_y_max),
                runif(1, input$s1_c_min, input$s1_c_max))
            } else if (c_length == 2) {
              c(runif(1, input$s1_0_min, input$s1_0_max),
                runif(1, input$s1_xstar_min, input$s1_xstar_max),
                runif(1, input$s1_y_min, input$s1_y_max),
                runif(1, input$s1_c1_min, input$s1_c1_max),
                runif(1, input$s1_c2_min, input$s1_c2_max))
            } else if (c_length == 3) {
              c(runif(1, input$s1_0_min, input$s1_0_max),
                runif(1, input$s1_xstar_min, input$s1_xstar_max),
                runif(1, input$s1_y_min, input$s1_y_max),
                runif(1, input$s1_c1_min, input$s1_c1_max),
                runif(1, input$s1_c2_min, input$s1_c2_max),
                runif(1, input$s1_c3_min, input$s1_c3_max))
            },
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }
  }, ignoreInit = TRUE)

  output$final3a <- renderUI({
    req(input$file)
    HTML(
      paste(
        "<b>Bias-adjusted OR<sub>YX</sub>:</b>",
        round(median(est_emc_sel()), 2)
      )
    )
  })

  output$final3b <- renderUI({
    req(input$file)
    HTML(
      paste(
        "<b>Bias-adjusted OR<sub>YX</sub> confidence interval:</b> (",
        round(
          quantile(
            est_emc_sel(),
            (1 - input$level) - ((1 - input$level) / 2)
          ),
          2
        ),
        ", ",
        round(
          quantile(
            est_emc_sel(), input$level + ((1 - input$level) / 2)
          ),
          2
        ),
        ")",
        sep = ""
      )
    )
  })

  output$final3c <- renderPlot({
    req(input$file)
    hist(
      est_emc_sel(),
      main = "Distribution of bias-adjusted estimates from bootstrap samples",
      xlab = expression('OR'["YX"]),
      col = "red"
    )
  })

  # UC, EMC, Sel ----

  est_uc_emc_sel <- eventReactive(input$go, {

    showNotification("Results are loading. The loading time will vary
                      depending on the size of the data and the number of
                      bootstrap samples.")
    no_cores <- 3
    registerDoParallel(cores = no_cores)
    cl <- makeCluster(no_cores)
    nreps <- input$bs_samples

    if (input$prob_bias_param == FALSE) {

      est <- vector(length = nreps)
      df <- data()
      px1 <- px1_vals()
      ps1 <- ps1_vals()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_uc_emc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = c(
              input$u1_0,
              input$u1_x,
              input$u1_y
            ),
            x_model_coefs = px1,
            s_model_coefs = ps1,
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }

    if (input$prob_bias_param == TRUE & input$distribution == "Normal") {

      c_length <- length(input$c_select)
      est <- vector(length = nreps)
      df <- data()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_uc_emc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = c(
              rnorm(1, input$u1_0_mean, input$u1_0_sd),
              rnorm(1, input$u1_x_mean, input$u1_x_sd),
              rnorm(1, input$u1_y_mean, input$u1_y_sd)
            ),
            x_model_coefs = if (c_length == 0) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd))
            } else if (c_length == 1) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c_mean, input$x1_c_sd))
            } else if (c_length == 2) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c1_mean, input$x1_c1_sd),
                rnorm(1, input$x1_c2_mean, input$x1_c2_sd))
            } else if (c_length == 3) {
              c(rnorm(1, input$x1_0_mean, input$x1_0_sd),
                rnorm(1, input$x1_xstar_mean, input$x1_xstar_sd),
                rnorm(1, input$x1_y_mean, input$x1_y_sd),
                rnorm(1, input$x1_c1_mean, input$x1_c1_sd),
                rnorm(1, input$x1_c2_mean, input$x1_c2_sd),
                rnorm(1, input$x1_c3_mean, input$x1_c3_sd))
            },
            s_model_coefs = if (c_length == 0) {
              c(rnorm(1, input$s1_0_mean, input$s1_0_sd),
                rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
                rnorm(1, input$s1_y_mean, input$s1_y_sd))
            } else if (c_length == 1) {
              c(rnorm(1, input$s1_0_mean, input$s1_0_sd),
                rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
                rnorm(1, input$s1_y_mean, input$s1_y_sd),
                rnorm(1, input$s1_c_mean, input$s1_c_sd))
            } else if (c_length == 2) {
              c(rnorm(1, input$s1_0_mean, input$s1_0_sd),
                rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
                rnorm(1, input$s1_y_mean, input$s1_y_sd),
                rnorm(1, input$s1_c1_mean, input$s1_c1_sd),
                rnorm(1, input$s1_c2_mean, input$s1_c2_sd))
            } else if (c_length == 3) { 
              c(rnorm(1, input$s1_0_mean, input$s1_0_sd),
                rnorm(1, input$s1_xstar_mean, input$s1_xstar_sd),
                rnorm(1, input$s1_y_mean, input$s1_y_sd),
                rnorm(1, input$s1_c1_mean, input$s1_c1_sd),
                rnorm(1, input$s1_c2_mean, input$s1_c2_sd),
                rnorm(1, input$s1_c3_mean, input$s1_c3_sd))
            },
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }

    if (input$prob_bias_param == TRUE & input$distribution == "Uniform") {

      c_length <- length(input$c_select)
      est <- vector(length = nreps)
      df <- data()

      or <- foreach(
        i = 1:nreps, .combine = c, .packages = c("dplyr", "multibias"),
        .export = c("input", "isolate")
      ) %dopar% {

        bdf <- df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), ]

        isolate({
          est[i] <- adjust_uc_emc_sel(
            data = bdf,
            exposure = input$x_select,
            outcome = input$y_select,
            confounders = input$c_select,
            u_model_coefs = c(
              rnorm(1, input$u1_0_min, input$u1_0_max),
              rnorm(1, input$u1_x_min, input$u1_x_max),
              rnorm(1, input$u1_y_min, input$u1_y_max)
            ),
            x_model_coefs = if (c_length == 0) {
              c(runif(1, input$x1_0_min, input$x1_0_max),
                runif(1, input$x1_xstar_min, input$x1_xstar_max),
                runif(1, input$x1_y_min, input$x1_y_max))
            } else if (c_length == 1) {
              c(runif(1, input$x1_0_min, input$x1_0_max),
                runif(1, input$x1_xstar_min, input$x1_xstar_max),
                runif(1, input$x1_y_min, input$x1_y_max),
                runif(1, input$x1_c_min, input$x1_c_max))
            } else if (c_length == 2) {
              c(runif(1, input$x1_0_min, input$x1_0_max),
                runif(1, input$x1_xstar_min, input$x1_xstar_max),
                runif(1, input$x1_y_min, input$x1_y_max),
                runif(1, input$x1_c1_min, input$x1_c1_max),
                runif(1, input$x1_c2_min, input$x1_c2_max))
            } else if (c_length == 3) {
              c(runif(1, input$x1_0_min, input$x1_0_max),
                runif(1, input$x1_xstar_min, input$x1_xstar_max),
                runif(1, input$x1_y_min, input$x1_y_max),
                runif(1, input$x1_c1_min, input$x1_c1_max),
                runif(1, input$x1_c2_min, input$x1_c2_max),
                runif(1, input$x1_c3_min, input$x1_c3_max))
            },
            s_model_coefs = if (c_length == 0) {
              c(runif(1, input$s1_0_min, input$s1_0_max),
                runif(1, input$s1_xstar_min, input$s1_xstar_max),
                runif(1, input$s1_y_min, input$s1_y_max))
            } else if (c_length == 1) {
              c(runif(1, input$s1_0_min, input$s1_0_max),
                runif(1, input$s1_xstar_min, input$s1_xstar_max),
                runif(1, input$s1_y_min, input$s1_y_max),
                runif(1, input$s1_c_min, input$s1_c_max))
            } else if (c_length == 2) {
              c(runif(1, input$s1_0_min, input$s1_0_max),
                runif(1, input$s1_xstar_min, input$s1_xstar_max),
                runif(1, input$s1_y_min, input$s1_y_max),
                runif(1, input$s1_c1_min, input$s1_c1_max),
                runif(1, input$s1_c2_min, input$s1_c2_max))
            } else if (c_length == 3) {
              c(runif(1, input$s1_0_min, input$s1_0_max),
                runif(1, input$s1_xstar_min, input$s1_xstar_max),
                runif(1, input$s1_y_min, input$s1_y_max),
                runif(1, input$s1_c1_min, input$s1_c1_max),
                runif(1, input$s1_c2_min, input$s1_c2_max),
                runif(1, input$s1_c3_min, input$s1_c3_max))
            },
            level = input$level
          )$estimate
        })
      }
      stopCluster(cl)
      return(or)
    }
  }, ignoreInit = TRUE)

  output$final4a <- renderUI({
    req(input$file)
    HTML(
      paste(
        "<b>Bias-adjusted OR<sub>YX</sub>:</b>",
        round(median(est_uc_emc_sel()), 2)
      )
    )
  })

  output$final4b <- renderUI({
    req(input$file)
    HTML(
      paste(
        "<b>Bias-adjusted OR<sub>YX</sub> confidence interval:</b> (",
        round(
          quantile(
            est_uc_emc_sel(),
            (1 - input$level) - ((1 - input$level) / 2)
          ),
          2
        ),
        ", ",
        round(
          quantile(
            est_uc_emc_sel(),
            input$level + ((1 - input$level) / 2)
          ),
          2
        ),
        ")",
        sep = ""
      )
    )
  })

  output$final4c <- renderPlot({
    req(input$file)
    hist(
      est_uc_emc_sel(),
      main = "Distribution of bias-adjusted estimates from bootstrap samples",
      xlab = expression('OR'["YX"]),
      col = "red"
    )
  })

}

# Run the application
shinyApp(ui = ui, server = server)