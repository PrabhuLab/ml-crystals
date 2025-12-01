# --- Interactive Analysis App for crystract ---
# Filename: app.R

# 1. SETUP: Load all necessary libraries
library(shiny)
library(shinythemes)
library(crystract)
library(data.table)
library(ggplot2)
library(plotly)
library(DT)

# 2. ONE-TIME ANALYSIS: This code runs only ONCE when the app starts.
# It performs the heavy lifting and stores the results in memory.
# ==============================================================================
cat("--- Performing one-time analysis of all CIF files. Please wait... ---\n")

# --- Timing Start ---
app_startup_time <- system.time({
  # --- File Loading ---
  cif_dir <- "/Users/donngo/repos/ml-crystals/packages/crystract/workflow/ClathrateI_all"
  all_cif_paths <- list.files(path = cif_dir, pattern = "\\.cif$", full.names = TRUE)

  # --- Full Analysis ---
  analysis_results <- analyze_cif_files(all_cif_paths)

  # --- Data Processing ---
  target_wyckoff_symbols <- c("6c", "16i", "24k")

  process_file_results <- function(i) {
    if (is.null(analysis_results$unit_cell_metrics[[i]])) return(NULL)
    bonds <- analysis_results$bonded_pairs[[i]]; coords <- analysis_results$atomic_coordinates[[i]]; metrics <- analysis_results$unit_cell_metrics[[i]]
    filtered_bonds <- filter_by_wyckoff_symbol(data_table = bonds, atomic_coordinates = coords, atom_col = "Atom1", wyckoff_symbols = target_wyckoff_symbols)
    avg_distance <- if (nrow(filtered_bonds) > 0) mean(filtered_bonds$Distance) else NA_real_
    lattice_a <- metrics$`_cell_length_a`
    return(data.table(file_name = basename(analysis_results$file_name[i]), lattice_parameter_a = lattice_a, average_distance = avg_distance))
  }

  plot_data <- rbindlist(lapply(1:nrow(analysis_results), process_file_results), use.names = TRUE, fill = TRUE)
  plot_data <- na.omit(plot_data)
}) # --- Timing End ---

cat("--- Analysis complete. App is ready. ---\n")
print(app_startup_time)
# ==============================================================================


# 3. USER INTERFACE (UI): Defines the layout and appearance of the app.
# ==============================================================================
ui <- fluidPage(
  theme = shinytheme("journal"),
  titlePanel("Interactive Clathrate Structure Analysis"),

  sidebarLayout(
    sidebarPanel(
      h4("Exploration Controls"),
      p("The full analysis is complete. Use the controls below to explore the results."),

      # Dropdown menu to select a file
      selectInput(
        "selected_file",
        "Select a File to View Detailed Data:",
        choices = sort(basename(analysis_results$file_name)),
        selected = sort(basename(analysis_results$file_name))[1]
      ),

      hr(),
      h5("About"),
      p(paste("This app analyzed", length(all_cif_paths), "CIF files in", round(app_startup_time['elapsed'], 1), "seconds."))
    ),

    mainPanel(
      # The main output is organized into tabs
      tabsetPanel(
        type = "tabs",
        tabPanel("Main Plot", plotlyOutput("main_plot", height = "600px")),
        tabPanel("Processed Data", DTOutput("processed_data_table")),
        tabPanel("Detailed File Analysis",
                 h4("Showing detailed tables for:"),
                 strong(textOutput("detail_title")),
                 hr(),
                 # UI elements to hold the nested tables for the selected file
                 h5("Unit Cell Metrics"),
                 DTOutput("metrics_table"),
                 h5("Atomic Coordinates"),
                 DTOutput("coords_table"),
                 h5("Bonded Pairs"),
                 DTOutput("bonds_table")
        )
      )
    )
  )
)
# ==============================================================================


# 4. SERVER: Defines the logic of the app.
# It tells R how to build the outputs based on user inputs.
# ==============================================================================
server <- function(input, output) {

  # --- Main Plot ---
  output$main_plot <- renderPlotly({
    p <- ggplot(plot_data, aes(x = lattice_parameter_a, y = average_distance, text = file_name)) +
      geom_point(alpha = 0.7, size = 2.5, color = "#0072B2") +
      geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linetype = "dashed") +
      labs(title = "Average Interatomic Distance vs. Lattice Parameter 'a'", x = "Lattice Parameter a (Å)", y = "Average Interatomic Distance (Å)") +
      theme_bw(base_size = 14)
    ggplotly(p, tooltip = c("x", "y", "text"))
  })

  # --- Processed Data Table ---
  output$processed_data_table <- renderDT({
    datatable(plot_data, extensions = 'Buttons', options = list(dom = 'Bfrtip', buttons = c('copy', 'csv')), rownames = FALSE)
  })

  # --- Detailed Analysis Tables (Reactive) ---
  # This section reacts to the user's selection in the dropdown menu.

  # Find the row index corresponding to the selected file
  selected_file_index <- eventReactive(input$selected_file, {
    which(basename(analysis_results$file_name) == input$selected_file)
  })

  # Title for the detail section
  output$detail_title <- renderText({
    input$selected_file
  })

  # Render the metrics table for the selected file
  output$metrics_table <- renderDT({
    req(selected_file_index()) # Require an index to be found
    df <- analysis_results$unit_cell_metrics[[selected_file_index()]]
    if (!is.null(df)) datatable(df, options = list(dom='t'), rownames=FALSE)
  })

  # Render the coordinates table for the selected file
  output$coords_table <- renderDT({
    req(selected_file_index())
    df <- analysis_results$atomic_coordinates[[selected_file_index()]]
    if (!is.null(df)) datatable(df, extensions = 'Buttons', options = list(dom='Bfrtip', buttons='csv'), rownames=FALSE)
  })

  # Render the bonds table for the selected file
  output$bonds_table <- renderDT({
    req(selected_file_index())
    df <- analysis_results$bonded_pairs[[selected_file_index()]]
    if (!is.null(df)) datatable(df, extensions = 'Buttons', options = list(dom='Bfrtip', buttons='csv'), rownames=FALSE)
  })

}
# ==============================================================================


# 5. RUN APP: Combines the UI and Server into a running application.
# ==============================================================================
shinyApp(ui = ui, server = server)
