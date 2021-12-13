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

        # STEP 3: select common cell lines of interest
        shinyjs::hidden(wellPanel(
            id = "panel_sel",
            checkboxGroupInput(inputId = "cells",
                                label = "Cell Lines of Interest:",
                                choices = NULL),

        # STEP 4: select drugs of interest
            checkboxGroupInput(inputId = "drugs",
                            label = "Drugs of Interest:",
                            choices = NULL)
        )),

        # STEP 5: compute correlations
        #actionButton(inputId = "computeCors",
        #             label = "Compute Correlations!"),

        # STEP 6: display plot(s)
    ), # end of sidebar panel

    # main panel for displaying outputs
    mainPanel(

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
    observeEvent(intersected(), {
        choices <- intersected()[[1]]@cell$cellid
        updateCheckboxGroupInput(inputId = "cells",
                          choices = choices,
                          selected = choices)
    })

    # update the list of choices in the input$drugs checkboxes
    observeEvent(intersected(), {
        choices <- intersected()[[1]]@drug$drugid
        updateCheckboxGroupInput(inputId = "drugs",
                          choices = choices,
                          selected = choices)
    })

    cells <- reactive({ # get cell lines chosen
        req(input$cells)
        return(input$cells)
    })

    drugs <- reactive({ # get drugs chosen
        req(input$drugs)
        return(input$drugs)
    })


}





# Run the application
shinyApp(ui = ui, server = server)

#[END]
