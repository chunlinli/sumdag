
library(shiny)

# a button: add edges
# two lists: from, to
# this updates d_mat
# a button: testing

# use a specific data
# p-value adjustment slidebar (scale??) ??
# a button: testing each single edges



# Define UI for application
ui <- fluidPage(
    # App title ----
    titlePanel("Inferring cardiovascular-related
                protein-protein interaction network"),
    fluidRow(
        column(
            3,
            h3("Buttons"),
            actionButton("action", "Action"),
            br(),
            br(),
            submitButton("Submit")
        ),
        column(
            3,
            h3("Single checkbox"),
            checkboxInput("checkbox", "Choice A", value = TRUE)
        ),
        column(
            3,
            checkboxGroupInput("checkGroup",
                h3("Checkbox group"),
                choices = list(
                    "Choice 1" = 1,
                    "Choice 2" = 2,
                    "Choice 3" = 3
                ),
                selected = 1
            )
        ),
        column(
            3,
            dateInput("date",
                h3("Date input"),
                value = "2014-01-01"
            )
        )
    ),
    fluidRow(
        column(
            3,
            dateRangeInput("dates", h3("Date range"))
        ),
        column(
            3,
            fileInput("file", h3("File input"))
        ),
        column(
            3,
            h3("Help text"),
            helpText(
                "Note: help text isn't a true widget,",
                "but it provides an easy way to add text to",
                "accompany other widgets."
            )
        ),
        column(
            3,
            numericInput("num",
                h3("Numeric input"),
                value = 1
            )
        )
    ),
    fluidRow(
        column(
            3,
            radioButtons("radio", h3("Radio buttons"),
                choices = list(
                    "Choice 1" = 1, "Choice 2" = 2,
                    "Choice 3" = 3
                ), selected = 1
            )
        ),
        column(
            3,
            selectInput("select", h3("Select box"),
                choices = list(
                    "Choice 1" = 1, "Choice 2" = 2,
                    "Choice 3" = 3
                ), selected = 1
            )
        ),
        column(
            3,
            sliderInput("slider1", h3("Sliders"),
                min = 0, max = 100, value = 50
            ),
            sliderInput("slider2", "",
                min = 0, max = 100, value = c(25, 75)
            )
        ),
        column(
            3,
            textInput("text", h3("Text input"),
                value = "Enter text..."
            )
        )
    ),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(

            # Input: Slider for the number of bins ----
            sliderInput(
                inputId = "bins",
                label = "Number of bins:",
                min = 1,
                max = 50,
                value = 30
            ),
            h2("Installation"),
            p("Shiny is available on CRAN, so you can install it in the usual way from your R console:"),
            code('install.packages("shiny")'),
            br(),
            br(),
            br(),
            br(),
            br(),
            "Shiny is a product of "
        ),

        # Main panel for displaying outputs ----
        mainPanel(

            # Output: Histogram ----
            plotOutput(outputId = "distPlot")
        )
        # mainPanel(
        #     p("p creates a paragraph of text."),
        #     p("A new p() command starts a new paragraph. Supply a style attribute to change the format of the entire paragraph.", style = "font-family: 'times'; font-si16pt"),
        #     strong("strong() makes bold text."),
        #     em("em() creates italicized (i.e, emphasized) text."),
        #     br(),
        #     code("code displays your text similar to computer code"),
        #     div("div creates segments of text with a similar style. This division of text is all blue because I passed the argument 'style = color:blue' to div", style = "color:blue"),
        #     br(),
        #     p(
        #         "span does the same thing as div, but it works with",
        #         span("groups of words", style = "color:blue"),
        #         "that appear inside a paragraph."
        #     )
        # )
    )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {

    # Histogram of the Old Faithful Geyser Data ----
    # with requested number of bins
    # This expression that generates a histogram is wrapped in a call
    # to renderPlot to indicate that:
    #
    # 1. It is "reactive" and therefore should be automatically
    #    re-executed when inputs (input$bins) change
    # 2. Its output type is a plot
    output$distPlot <- renderPlot({
        x <- faithful$waiting
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        hist(x,
            breaks = bins, col = "#007bc2", border = "white",
            xlab = "Waiting time to next eruption (in mins)",
            main = "Histogram of waiting times"
        )
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)



ui <- fluidPage(
  selectInput("dataset", "Dataset", c("diamonds", "rock", "pressure", "cars")),
  conditionalPanel( condition = "output.nrows",
                    checkboxInput("headonly", "Only use first 1000 rows"))
)
server <- function(input, output, session) {
  datasetInput <- reactive({
    switch(input$dataset,
           "rock" = rock,
           "pressure" = pressure,
           "cars" = cars)
  })

  output$nrows <- reactive({
    nrow(datasetInput())
  })

  outputOptions(output, "nrows", suspendWhenHidden = FALSE)
}

shinyApp(ui, server)
