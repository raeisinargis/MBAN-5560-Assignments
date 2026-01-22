library(shiny)
library(ggplot2)
library(MASS)
library(shinythemes)

# UI definition
ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  tags$head(
    tags$style(HTML("
      .well { background-color: #ffffff; }
      .nav-tabs { margin-bottom: 20px; }
      .control-label { font-weight: 500; }
      .shiny-output-error { color: #ff0000; }
      .btn { margin: 5px; }
    "))
  ),
  titlePanel("Dynamic Nonparametric Smoothing Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("data_type", "Select Data Type",
                  choices = c("Sine Wave" = "sine",
                              "Linear Trend" = "linear",
                              "Quadratic" = "quadratic",
                              "Volatile Series" = "volatile"),
                  selected = "sine"),
      
      numericInput("n_points", "Number of Points",
                   value = 300, min = 50, max = 1000),
      
      numericInput("noise_level", "Noise Level",
                   value = 0.25, min = 0, max = 2, step = 0.1),
      
      sliderInput("window_size", "Window Size",
                  min = 3, max = 51, value = 3, step = 2),
      
      selectInput("smooth_method", "Smoothing Method",
                  choices = c("Mean" = "mean",
                              "Median" = "median",
                              "Trimmed Mean" = "trimmed",
                              "Local Regression" = "loess",
                              "Kernel" = "kernel"),
                  selected = "mean"),
      
      conditionalPanel(
        condition = "input.smooth_method == 'trimmed'",
        sliderInput("trim_percent", "Trim Percentage",
                    min = 0, max = 0.4, value = 0.1, step = 0.05)
      ),
      
      conditionalPanel(
        condition = "input.smooth_method == 'loess'",
        sliderInput("degree", "Polynomial Degree",
                    min = 0, max = 2, value = 1, step = 1)
      ),
      
      conditionalPanel(
        condition = "input.smooth_method == 'kernel'",
        selectInput("kernel_type", "Kernel Type",
                    choices = c("Gaussian" = "gaussian",
                                "Epanechnikov" = "epanechnikov",
                                "Triangular" = "triangular"),
                    selected = "gaussian"),
        sliderInput("bandwidth", "Bandwidth Multiplier",
                    min = 0.2, max = 2, value = 1, step = 0.1)
      ),
      
      sliderInput("animation_speed", "Animation Speed (ms)",
                  min = 50, max = 1000, value = 100, step = 50),
      
      actionButton("start", "Start/Stop Animation"),
      actionButton("reset", "Reset"),
      
      verbatimTextOutput("stats"),
      
      # Add attribution and copyright
      tags$div(
        style = "margin-top: 30px; padding-top: 20px; border-top: 1px solid #e5e5e5;",
        tags$p(
          style = "color: #666; font-size: 0.9em; margin-bottom: 5px;",
          "Dr. Aydede - SMU",
          tags$br(),
          "Prepared for MBAN Students"
        ),
        tags$p(
          style = "color: #666; font-size: 0.8em; font-style: italic;",
          paste("© Copyright", format(Sys.Date(), "%Y"), "- All Rights Reserved")
        )
      )
    ),
    
    mainPanel(
      plotOutput("smoothing_plot", height = "600px"),
      wellPanel(
        uiOutput("method_explanation")
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Define kernel functions
  kernel_functions <- list(
    gaussian = function(x) dnorm(x, sd = 1),
    epanechnikov = function(x) ifelse(abs(x) <= 1, 3/4 * (1 - x^2), 0),
    triangular = function(x) ifelse(abs(x) <= 1, 1 - abs(x), 0)
  )
  
  # Initialize reactive values
  rv <- reactiveValues(
    current_pos = 1,
    is_running = FALSE
  )
  
  # Generate data based on inputs
  data <- reactive({
    n <- input$n_points
    set.seed(1)
    x <- sort(runif(n) * 2 * pi)
    
    y <- switch(input$data_type,
                "sine" = sin(x),
                "linear" = x/pi - 1,
                "quadratic" = (x/pi - 1)^2 - 1,
                "volatile" = {
                  base <- cumsum(rnorm(n, 0, 0.3))
                  spikes <- numeric(n)
                  spike_positions <- seq(20, n, by = round(n/10))
                  spikes[spike_positions] <- rnorm(length(spike_positions), 0, 2)
                  shifts <- cumsum(sample(c(0, 1, -1), n, prob = c(0.98, 0.01, 0.01), replace = TRUE))
                  base + spikes + shifts
                })
    
    y <- y + rnorm(n) * input$noise_level
    data.frame(time = 1:n, y = y)
  })
  
  # Store smoothed values
  smoothed_data <- reactiveVal(data.frame(time = numeric(0), y = numeric(0)))
  
  # Calculate smoothed value based on method
  calculate_smooth <- function(window_data) {
    values <- window_data$y
    times <- window_data$time
    
    switch(input$smooth_method,
           "mean" = mean(values),
           "median" = median(values),
           "trimmed" = mean(values, trim = input$trim_percent),
           "kernel" = {
             kernel_fn <- kernel_functions[[input$kernel_type]]
             center_time <- mean(times)
             h <- input$bandwidth * IQR(times)/1.34
             x_scaled <- (times - center_time) / h
             weights <- kernel_fn(x_scaled)
             weights <- weights / sum(weights)
             sum(values * weights)
           },
           "loess" = {
             if (length(values) >= 4) {
               if (input$degree == 0) {
                 # Use kernel smoothing for degree 0
                 kernel_fn <- kernel_functions[["gaussian"]]
                 center_time <- mean(times)
                 h <- input$window_size/4
                 x_scaled <- (times - center_time) / h
                 weights <- kernel_fn(x_scaled)
                 weights <- weights / sum(weights)
                 sum(values * weights)
               } else {
                 # Use robust regression for degree 1 and 2
                 if (input$degree == 1) {
                   fit <- try(MASS::rlm(values ~ times), silent = TRUE)
                   if (inherits(fit, "try-error")) {
                     fit <- lm(values ~ times)
                   }
                 } else {
                   fit <- lm(values ~ poly(times, 2))
                 }
                 predict(fit, newdata = data.frame(times = mean(times)))
               }
             } else {
               mean(values)
             }
           })
  }
  
  # Handle animation timing
  observe({
    if (rv$is_running) {
      invalidateLater(input$animation_speed)
      
      isolate({
        current_data <- data()
        window <- input$window_size
        
        if (rv$current_pos <= (nrow(current_data) - window)) {
          window_data <- current_data[rv$current_pos:(rv$current_pos + window), ]
          smooth_val <- calculate_smooth(window_data)
          
          new_point <- data.frame(
            time = rv$current_pos + (window/2 - 1),
            y = smooth_val
          )
          
          current_smoothed <- smoothed_data()
          smoothed_data(rbind(current_smoothed, new_point))
          
          rv$current_pos <- rv$current_pos + 1
        } else {
          rv$is_running <- FALSE
        }
      })
    }
  })
  
  # Handle start/stop button
  observeEvent(input$start, {
    rv$is_running <- !rv$is_running
  })
  
  # Handle reset button
  observeEvent(input$reset, {
    rv$current_pos <- 1
    rv$is_running <- FALSE
    smoothed_data(data.frame(time = numeric(0), y = numeric(0)))
  })
  
  # Render the plot
  output$smoothing_plot <- renderPlot({
    current_data <- data()
    current_smoothed <- smoothed_data()
    
    p <- ggplot(current_data, aes(x = time, y = y)) +
      geom_point(alpha = 0.5) +
      theme_minimal() +
      labs(title = "Dynamic Nonparametric Smoothing",
           x = "Time", y = "Value")
    
    if (rv$current_pos <= (nrow(current_data) - input$window_size)) {
      window_start <- rv$current_pos
      window_end <- rv$current_pos + input$window_size
      window_data <- current_data[window_start:window_end, ]
      
      # Add vertical window lines
      p <- p + 
        geom_vline(xintercept = window_start, color = "green", alpha = 0.5) +
        geom_vline(xintercept = window_end, color = "green", alpha = 0.5)
      
      if (input$smooth_method == "kernel") {
        # Calculate kernel weights
        kernel_fn <- kernel_functions[[input$kernel_type]]
        mid_point <- (window_start + window_end) / 2
        h <- input$bandwidth * IQR(window_data$time)/1.34
        
        # Create smooth kernel curve
        curve_points <- seq(window_start, window_end, length.out = 100)
        x_scaled <- (curve_points - mid_point) / h
        kernel_curve <- data.frame(
          time = curve_points,
          height = kernel_fn(x_scaled)
        )
        
        # Calculate weights and weighted mean
        x_scaled_data <- (window_data$time - mid_point) / h
        weights <- kernel_fn(x_scaled_data)
        weights <- weights / sum(weights)
        weighted_mean <- sum(window_data$y * weights)
        
        # Scale kernel curve for visualization
        y_range <- diff(range(window_data$y))
        scale_factor <- y_range / max(kernel_curve$height) * 0.3
        kernel_curve$height <- weighted_mean + kernel_curve$height * scale_factor
        
        # Show kernel shape and weighted mean
        p <- p + 
          geom_line(data = kernel_curve,
                    aes(x = time, y = height, group = 1),
                    color = "blue", linewidth = 1) +
          geom_hline(yintercept = weighted_mean,
                     color = "blue", linetype = "dashed") +
          geom_point(data = data.frame(
            time = window_data$time,
            y = window_data$y,
            weight = weights
          ),
          aes(x = time, y = y, size = weight),
          color = "blue", alpha = 0.6) +
          scale_size_continuous(range = c(1, 4), guide = "none")
      } else if (input$smooth_method == "loess") {
        if (input$degree == 0) {
          # For degree 0, show horizontal line
          kernel_fn <- kernel_functions[["gaussian"]]
          center_time <- mean(window_data$time)
          h <- input$window_size/4
          x_scaled <- (window_data$time - center_time) / h
          weights <- kernel_fn(x_scaled)
          weights <- weights / sum(weights)
          weighted_mean <- sum(window_data$y * weights)
          
          p <- p + geom_segment(data = data.frame(x = window_start, xend = window_end,
                                                  y = weighted_mean, yend = weighted_mean),
                                aes(x = x, y = y, xend = xend, yend = yend),
                                color = "blue", linewidth = 1)
        } else if (input$degree == 1) {
          # For degree 1, use robust linear regression
          tryCatch({
            fit <- MASS::rlm(y ~ time, data = window_data)
            y_start <- predict(fit, newdata = data.frame(time = window_start))
            y_end <- predict(fit, newdata = data.frame(time = window_end))
            p <- p + geom_segment(data = data.frame(x = window_start, xend = window_end,
                                                    y = y_start, yend = y_end),
                                  aes(x = x, y = y, xend = xend, yend = yend),
                                  color = "blue", linewidth = 1)
          }, error = function(e) {
            fit <- lm(y ~ time, data = window_data)
            y_start <- predict(fit, newdata = data.frame(time = window_start))
            y_end <- predict(fit, newdata = data.frame(time = window_end))
            p <- p + geom_segment(data = data.frame(x = window_start, xend = window_end,
                                                    y = y_start, yend = y_end),
                                  aes(x = x, y = y, xend = xend, yend = yend),
                                  color = "blue", linewidth = 1)
          })
        } else {
          # For degree 2, use polynomial regression
          fit <- lm(y ~ poly(time, 2), data = window_data)
          pred_times <- seq(window_start, window_end, length.out = 50)
          fitted_df <- data.frame(
            time = pred_times,
            y = predict(fit, newdata = data.frame(time = pred_times))
          )
          p <- p + geom_line(data = fitted_df, aes(x = time, y = y, group = 1),
                             color = "blue", linewidth = 1)
        }
      } else {
        # For other methods, show horizontal line at smoothed value
        smooth_val <- calculate_smooth(window_data)
        p <- p + geom_segment(data = data.frame(x = window_start, xend = window_end,
                                                y = smooth_val, yend = smooth_val),
                              aes(x = x, y = y, xend = xend, yend = yend),
                              color = "blue", linewidth = 1)
      }
    }
    
    # Add smoothed line
    if (nrow(current_smoothed) > 0) {
      p <- p + geom_line(data = current_smoothed, 
                         aes(x = time, y = y, group = 1),
                         color = "red", linewidth = 1)
    }
    
    p
  })
  
  # Add explanations for each method
  output$method_explanation <- renderUI({
    method <- input$smooth_method
    
    explanation <- switch(method,
                          "mean" = HTML("
        <h4>Moving Average</h4>
        <p>The moving average method calculates the simple arithmetic mean of all points within the window. 
        Each point in the window receives equal weight in the calculation. This method is straightforward but 
        can be sensitive to outliers and may not handle sharp transitions well.</p>
        <p>The blue horizontal line shows the current window's average, which becomes part of the red smoothed line.</p>
      "),
                          "median" = HTML("
        <h4>Moving Median</h4>
        <p>The moving median takes the middle value of all points within the window. This method is more 
        robust to outliers than the moving average, making it useful for data with occasional extreme values 
        or noise. However, it may produce a less smooth result than other methods.</p>
        <p>The blue horizontal line shows the current window's median value.</p>
      "),
                          "trimmed" = HTML("
        <h4>Trimmed Mean</h4>
        <p>The trimmed mean removes a specified percentage of extreme values from both ends of the data 
        within the window before calculating the average. This combines the smoothness of the mean with 
        some of the robustness of the median. The trim percentage controls how many extreme values are 
        excluded from the calculation.</p>
        <p>The blue horizontal line shows the trimmed mean for the current window.</p>
      "),
                          "loess" = HTML(sprintf("
        <h4>Local Regression (Degree %d)</h4>
        <p>Local regression fits a polynomial of specified degree to the points within each window. 
        Degree 0 fits a constant (similar to weighted mean), degree 1 fits a straight line, and 
        degree 2 fits a quadratic curve. This method can capture more complex patterns in the data.</p>
        <p>The blue line shows the fitted polynomial for the current window. Notice how it adapts to 
        local patterns in the data.</p>
      ", input$degree)),
                          "kernel" = HTML(sprintf("
        <h4>Kernel Smoothing (%s)</h4>
        <p>Kernel smoothing calculates a weighted average where points closer to the center of the window 
        receive higher weights than points near the edges. The blue curve shows the kernel function that 
        determines these weights - points under the peak receive the most weight. The size of the blue 
        points indicates how much each value contributes to the final average.</p>
        <p>The bandwidth parameter controls how quickly the weights decay with distance: larger values 
        create smoother results but may miss local features, while smaller values track the data more 
        closely but may be more sensitive to noise.</p>
        <p>In this visualization, the dashed blue line shows the weighted average that becomes part of 
        the final smoothed result (red line). The blue curve above shows the shape of the kernel 
        function, indicating how much weight each point receives based on its distance from the center.</p>
      ", input$kernel_type))
    )
    
    tagList(
      hr(),
      h3("Method Explanation"),
      explanation
    )
  })
  
  # Update statistics display
  output$stats <- renderText({
    method_info <- if(input$smooth_method == "kernel") {
      sprintf("%s kernel (bandwidth %.1f)", input$kernel_type, input$bandwidth)
    } else if(input$smooth_method == "loess") {
      sprintf("%s (degree %d)", input$smooth_method, input$degree)
    } else {
      input$smooth_method
    }
    
    sprintf("Position: %d/%d\nWindow size: %d\nMethod: %s",
            rv$current_pos,
            nrow(data()),
            input$window_size,
            method_info)
  })
}

# Run the application
shinyApp(ui = ui, server = server)