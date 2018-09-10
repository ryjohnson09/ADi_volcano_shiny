###############################
### Volcano Plots Shiny App ###
###############################

library(shiny)
library(tidyverse)
library(plotly)
source("IVTT_stats.R", local = TRUE)
source("PurProt_stats.R", local = TRUE)


###############
### Choices ###
###############

ADi_datasets <- c("ETEC_IgG_IVTT", "ETEC_IgA_IVTT", 
                  "PanEC_IgG_IVTT", "PanEC_IgA_IVTT", 
                  "Purified_Protein_ETEC_IgG", "Purified_Protein_ETEC_IgA")

matched_choices <- list("All Samples" = "all_samples",
                        "Matched Samples" = "matched_samples")

visit_choices <- c("Visit 1 vs 4",
                   "Visit 1 vs 5",
                   "Acute vs Conv")

detection_choices <- c("Culture", "Taq", "Both", "Either")


# Pathogen choices

treat <- suppressWarnings(suppressMessages(read_csv("TrEAT_Clinical_Metadata_tidy.csv")))

culture_choices <- c("All", colnames(treat)[grep(pattern = "_culture$", colnames(treat))])

taq_choices <- c("All", colnames(treat)[grep(pattern = "_taq$", colnames(treat))])

both_choices <- c("All", colnames(treat)[grep(pattern = "_both$", colnames(treat))])

either_choices <- c("All", colnames(treat)[grep(pattern = "_either$", colnames(treat))])

## UI ------------------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Exploring the ADi Datasets"),
  
  # Dataset choice
  selectInput('dataset', 'Choose a dataset:', choices = ADi_datasets),
  
  
  # Visit
  radioButtons('visit', 'Visit:', choices = visit_choices, selected = "Visit 1 vs 5", inline = TRUE),
  helpText("Select patient samples from specified visits"),
  
  
  # Matched samples?
  radioButtons("matched", label = "Matched or All samples",
               choices = matched_choices, inline = TRUE, selected = "matched_samples"),
  helpText("Matched = Only patients that provided samples for all selected visits"),
  
  
  # Pathogen selection
  selectInput("path_detection", "Pathogen Detection Method", choices = detection_choices, selected = "Both"),
  helpText("How are samples determined to be positive for pathogen"),
  # Pathogen Choices
  uiOutput("secondSelection"),
  helpText("If", code("All"), ", then all samples included."),
  
 
  mainPanel(
    plotlyOutput("plot", width = "150%", height = "100%"),
    htmlOutput("legend"),
    htmlOutput("text"))
   
)


## Server --------------------------------------------------------------------------------

server <- function(input, output){
  
  ####################
  ### Read in Data ###
  ####################
  
  new_dataset <- reactive({
    adi_data <- read_csv(ifelse(input$dataset == "ETEC_IgG_IVTT", "ETEC_IgG_IVTT_RawData_tidy.csv",
             ifelse(input$dataset == "ETEC_IgA_IVTT", "ETEC_IgA_IVTT_RawData_tidy.csv",
             ifelse(input$dataset == "PanEC_IgG_IVTT", "PanEC_IgG_IVTT_RawData_tidy.csv",
             ifelse(input$dataset == "PanEC_IgA_IVTT", "PanEC_IgA_IVTT_RawData_tidy.csv",
             ifelse(input$dataset == "Purified_Protein_ETEC_IgG", "IgG_PurifiedProtein_RawData_tidy.csv",
             ifelse(input$dataset == "Purified_Protein_ETEC_IgA", "IgA_PurifiedProtein_RawData_tidy.csv",
             stopApp("Invalid Data Set"))))))))

    # Filter out LOP and PLA samples
    filter(adi_data, !Treatment %in% c("LOP", "PLA"))
  })
  

  
  ############################
  ### Filter by Visit Type ###
  ############################
  visit_filter <- reactive({
    if(input$visit == "Visit 1 vs 4"){
      filter(new_dataset(), visit %in% c(1,4))

    } else if(input$visit == "Visit 1 vs 5"){
      filter(new_dataset(), visit %in% c(1,5))

    } else if(input$visit == "Acute vs Conv"){
      filter(new_dataset(), visit_type %in% c("acute", "convalescent"))

    } else {
      stopApp("Invalid Visit Number")
    }
  })


  ##################################
  ### Filter for matched samples ###
  ##################################

  matched_dataset <- reactive({
    if(input$matched == "matched_samples" & input$visit != "Acute vs Conv"){

      # Get length of visit number
      visit_numbers <- length(unique(visit_filter()$visit))

      # Get counts of patient IDs with visit number
      patient_visit_counts <- visit_filter() %>%
        select(Patients, visit) %>%
        unique() %>%
        group_by(Patients) %>%
        summarise(n =n()) %>%

        # Then filter patients that have that number
        filter(n == visit_numbers) %>%
        pull(Patients)

      # Filter the visit_filter dataset to only
      # inlucde those isolates
      visit_filter() %>%
        filter(Patients %in% patient_visit_counts)

    } else {
        visit_filter()
      }
    })



  #######################################
  ### Render pathogen detection list ###
  #######################################
  output$secondSelection <- renderUI({
    if(input$path_detection == "Culture"){
      selectInput("pathogens", "Pathogens", choices = culture_choices)
    } else if (input$path_detection == "Taq"){
      selectInput("pathogens", "Pathogens", choices = taq_choices)
    } else if (input$path_detection == "Either"){
      selectInput("pathogens", "Pathogens", choices = either_choices)
    } else {
      selectInput("pathogens", "Pathogens", choices = both_choices)
    }
  })


  ##########################
  ### Filter by Pathogen ###
  ##########################
  treat_pathogen_filtered <- reactive({
    
    # Include all patients
    if (input$pathogens == "All"){
      matched_dataset()
    } else {
      # Select for patients with specified pathogen
      pathogen <- sym(input$pathogens)
      matched_dataset() %>%
        filter(!!pathogen == "yes")
    }
  })


  
  stats_dataset <- reactive({

    # Set Groupings
    if (input$visit == "Visit 1 vs 4"){
      group_column <- "visit"
      group_variables <- c(1,4)
    } else if (input$visit == "Visit 1 vs 5"){
      group_column <- "visit"
      group_variables <- c(1,5)
    } else if (input$visit == "Acute vs Conv"){
      group_column <- "visit_type"
      group_variables <- c("acute", "convalescent")
    } else {
      stopApp("Invalid groupings")
    }

    if (startsWith(input$dataset, "Pur")) {
      PurProt_stats(treat_pathogen_filtered(),
                    group_column = group_column,
                    group_variables = group_variables)
    } else {
      IVTT_stats(treat_pathogen_filtered(),
                 group_column = group_column,
                 group_variables = group_variables)
    }
  })



  output$plot <- renderPlotly({
    ggplot(data = stats_dataset(), aes(x = mean_diff, y = -log(p_value))) +
      geom_point(aes(text = ID_Description),size = 3, alpha = 0.7) +
      geom_hline(yintercept = -log(0.05)) +
      geom_hline(yintercept = -log(0.1), linetype = "dashed") +
      ylab("-log(p_value)") +
      xlab("Mean Difference") +
      ggtitle(input$dataset) +
      theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0)
      )
  })
  
  # Text under graph
  output$legend <- renderUI({
    HTML(paste("Proteins above black line have a p-value < 0.05", 
               "Proteins above dashed line have a p-value < 0.1" , sep="<br/>"))
  })
  
  
  output$text <- renderUI({ 
    if(input$visit == "Visit 1 vs 4"){
      visit_1_num <- unique(stats_dataset()$count.1)
      visit_2_num <- unique(stats_dataset()$count.2)
      group_1 <- "Visit 1"
      group_2 <- "Visit 4"
    } else if(input$visit == "Visit 1 vs 5") {
      visit_1_num <- unique(stats_dataset()$count.1)
      visit_2_num <- unique(stats_dataset()$count.2)
      group_1 <- "Visit 1"
      group_2 <- "Visit 5"
    } else if(input$visit == "Acute vs Conv") {
      visit_1_num <- unique(stats_dataset()$count.1)
      visit_2_num <- unique(stats_dataset()$count.2)
      group_1 <- "Acute"
      group_2 <- "Convalescent"
    } else {
      stopApp("Invalid text below graph")}
    
    str1 <- paste0("# of Patients in ", group_1, ": ", visit_1_num)
    str2 <- paste0("# of Patients in ", group_2, ": ", visit_2_num)
    
    HTML(paste(str1, str2, sep = "<br/>"))
  })
  
  

 }

shinyApp(server = server, ui = ui)