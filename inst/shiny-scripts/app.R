#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)

# Define UI for application
ui <- fluidPage(
    useShinyjs(),

    # Application title
    titlePanel("Identifying Consistent Cancer Cell Lines"),

    # Sidebar
    sidebarLayout(
        # sidebar panel for inputs
        sidebarPanel(
        # STEP 1: upload PharmacoSets
        tags$p("Provided two to three PharmacoSets (either downloaded from Orcestra or
               from the PharmacoGx R package onto your local computer as .rds files),
               intersect them."),

        fileInput(inputId = "pset1",
                  NULL,
                  buttonLabel = "Upload PSet 1",
                  multiple = FALSE,
                  accept = c(".rds")),

        fileInput(inputId = "pset2",
                  NULL,
                  buttonLabel = "Upload PSet 2",
                  multiple = FALSE,
                  accept = c(".rds")),

        fileInput(inputId = "pset3",
                  NULL,
                  buttonLabel = "Upload PSet 3",
                  multiple = FALSE,
                  accept = c(".rds")),

        actionButton(inputId = "upload",
                     label = "Upload PSets!"),

        # STEP 2: intersect psets
        shinyjs::hidden(wellPanel(
            id = "panel_int",
            tags$p("Select categories on which the PharmacoSets should be intersected."),

            checkboxGroupInput(inputId = "intersectOn",
                               label = "Intersect on:",
                               choices = c("drugs", "cell.lines", "concentrations"),
                               selected = c("drugs")),

            actionButton(inputId = "intersect",
                         label = "Intersect PSets!")
        )),

        # STEP 3: select coefficients of interest
        shinyjs::hidden(wellPanel(
            id = "panel_sel",
            checkboxGroupInput(inputId = "coefs",
                                label = "Correlation Coefficients of Interest:",
                                choices = c("pearson", "spearman", "kendall")),

            # TO-DO: add text box explaining what pearson, spearman, kendall are

            # STEP 4: select drugs of interest
            checkboxGroupInput(inputId = "drugs",
                            label = "Drugs of Interest:",
                            choices = NULL),

            # STEP 5: select sensitivity measures of interest
            checkboxGroupInput(inputId = "sens",
                               label = "Sensitivity Measures of Interest:",
                               choices = NULL),

            checkboxInput(inputId = "pval",
                          label = "Include P-Values?"),

            # STEP 6: compute cell line correlations
            actionButton(inputId = "computeCors",
                         label = "Compute Correlations!")

        )),
        # display cell line correlations

        # STEP 7: get consistent cell lines

        # STEP 8: get drug correlations

        # STEP 9: compute concordance


        # STEP 10: plot(s)


    ), # end of sidebar panel

    # main panel for displaying outputs
    mainPanel(
        shinyjs::hidden(wellPanel(
            id = "plotChoices",
            # choose which sens measure to plot
            checkboxGroupInput(inputId = "plotSens",
                               label = "Sensitivity Measure to plot:",
                               choices = NULL),
            # choose which coef to plot
            checkboxGroupInput(inputId = "plotCoefs",
                               label = "Correlation Coefficient to plot:",
                               choices = NULL),
            textInput(inputId = "plotTitle",
                      label = "Title for Plot:"),
            #actionButton(inputId = "showPlot",
            #             label = "Update Plot!")
            plotOutput("cellPlot")
        )),
    ) # end of main panel
)
)

# Define server logic
options(shiny.maxRequestSize=200*1024^2) # increase max file input size
# pset_list <- list()
server <- function(input, output) {
    # STEP 1: save uploaded psets as reactives

    # TO-DO: check if any files are uploaded before allowing to click upload button
    pset_list <- eventReactive(input$upload, { # wait until file is uploaded
        showModal(modalDialog("Uploading PSets", footer = NULL))
        if (!is.null(input$pset1) && !is.null(input$pset2)) {
            ext <- tools::file_ext(input$pset1$name)
            validate(need(ext == "rds", "Invalid pset1 file. Please upload a .rds file"))

            ext <- tools::file_ext(input$pset2$name)
            validate(need(ext == "rds", "Invalid pset2 file. Please upload a .rds file"))

            pset_list <- list()
            pset1 <- readRDS(input$pset1$datapath)
            pset2 <- readRDS(input$pset2$datapath)
            pset_list <- c(pset1, pset2)
        }

        if (!is.null(input$pset3)) {
            ext <- tools::file_ext(input$pset1$name)
            validate(need(ext == "rds", "Invalid file. Please upload a .rds file"))
            pset3 <- readRDS(input$pset3$datapath)
            pset_list <- c(pset_list, pset3)
        }
        removeModal()
        return(pset_list)
    })

    observeEvent(pset_list(),
        shinyjs::showElement(id = "panel_int")
    )
    # TO-DO: check if uploaded files are PharmacoSets before intersecting
    # TO-DO: handle intersectPSet errors (check what they can intersect on?)
    intersected <- eventReactive(input$intersect, { #intersect once button is clicked
        showModal(modalDialog("Intersecting PSets", footer = NULL))
        intersected <- PharmacoGx::intersectPSet(pSets = pset_list(),
                                  intersectOn = input$intersectOn,
                                  verbose = TRUE
                                )
        removeModal()
        shinyjs::showElement(id = "panel_sel")
        return(intersected)
    })

    # use observeEvent to perform an action in response to an event
    # use eventReactive to create a calculated value that only updates in response to an event

    # update the list of choices in the input$cells checkboxes
    #observeEvent(intersected(), {
    #    choices <- intersected()[[1]]@cell$cellid
    #    updateCheckboxGroupInput(inputId = "cells",
    #                      choices = choices,
    #                      selected = choices)
    #})

    # update the list of choices in the input$sens checkboxes
    observeEvent(intersected(), {
        choices <- intersectSensMeasures(intersected())
        updateCheckboxGroupInput(inputId = "sens",
                                 choices = choices)
    })

    # update the list of choices in the input$drugs checkboxes
    observeEvent(intersected(), {
        choices <- intersected()[[1]]@drug$drugid
        updateCheckboxGroupInput(inputId = "drugs",
                          choices = choices,
                          selected = choices)
    })

    #cells <- reactive({ # get cell lines chosen
    #    req(input$cells)
    #    return(input$cells)
    #})

    coefs <- reactive({ # get coefs chosen
        input$coefs
    })

    sens <- reactive({ # get sens chosen
        input$sens
    })

    drugs <- reactive({ # get drugs chosen
        input$drugs
    })

    pval <- reactive({
        input$pval #returns TRUE if checked, FALSE otherwise
    })

    # STEP 6: compute cell line correlations
    cors <- eventReactive(input$computeCors, { # compute once button is clicked
        showModal("Computing Correlations")
        cors <- computeCellLineCorrelation(pSet = intersected(),
                                           coefs = coefs(),
                                           sensMeasures = sens(),
                                           pval = pval())
        removeModal()
        shinyjs::showElement(id = "plotChoices")
        return(cors)
    })
    # TO-DO: handle computeCellLineCorrelation message returns (not enough obs to compute corrs)

    # update the list of choices in the input$drugs checkboxes
    observeEvent(cors(), {
        updateCheckboxGroupInput(inputId = "plotSens",
                                 choices = sens())
        updateCheckboxGroupInput(inputId = "plotCoefs",
                                 choices = coefs())
    })

    #observeEvent(input$showPlot, {
    #    req(input$plotSens, input$plotCoefs)
    #    shinyjs::showElement(id = "plot1")
    #})

    plotSens <- reactive({
        input$plotSens
    })

    plotCoefs <- reactive({
        input$plotCoefs
    })

    plotTitle <- reactive({
        input$plotTitle
    })

    output$cellPlot <- renderPlot({
       plotCorrelations(correlations = cors(),
                        sensMeasure = plotSens(),
                        coefficient = plotCoefs(),
                        title = plotTitle())
    })

}





# Run the application
shinyApp(ui = ui, server = server)

#[END]
