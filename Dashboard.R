#
library(here)
here::i_am("Dashboard.R")
message("Project root: ", here::here())
library(shiny)
library(leaflet)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(lme4)
library(merTools)
library(here)
library(patchwork)
library(htmltools)
library(scales)
library(stringr)
library(conflicted)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("recode", "dplyr")
conflicted::conflict_prefer("lmer", "lme4")



load(here::here("Data", "Models.RData"))
combined_data <- read.csv(here::here("Data", "combined_labour.csv"), stringsAsFactors = FALSE)
rownames(combined_data) <- combined_data$X; combined_data$X <- NULL

combined_data <- combined_data %>%
  mutate(
    host_country   = factor(sub("_.*", "", rownames(.))),
    home_country   = sub(".*_", "", rownames(.)),
    regional_ident = factor(regional_ident)
  )

# Predictions
combined_data$pred_full   <- predict(random_slope_model,            newdata = combined_data, allow.new.levels = TRUE)
combined_data$pred_male   <- predict(random_slope_model_men_only,   newdata = combined_data, allow.new.levels = TRUE)
combined_data$pred_female <- predict(random_slope_model_women_only, newdata = combined_data, allow.new.levels = TRUE)
combined_data$obs_full    <- combined_data$ptunemp
combined_data$obs_male    <- combined_data$pmunemp
combined_data$obs_female  <- combined_data$pfunemp

# Long format
pred_long <- combined_data %>%
  select(host_country, home_country,
         pred_full, pred_male, pred_female,
         obs_full,  obs_male,  obs_female) %>%
  pivot_longer(cols = starts_with(c("pred","obs")),
               names_to = c("type","sample"),
               names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(discrepancy = pred - obs,
         sample = recode(sample,
                         full="Full Sample", male="Male Only", female="Female Only"))

# ---- Country name recode ----
ne_recode <- c(
  "United States"                       = "United States of America",
  "USA"                                 = "United States of America",
  "UK"                                  = "United Kingdom",
  "Czech Republic"                      = "Czechia",
  "Cape Verde"                          = "Cabo Verde",
  "Swaziland"                           = "Eswatini",
  "Ivory Coast"                         = "Côte d'Ivoire",
  "Cote d'Ivoire"                       = "Côte d'Ivoire",
  "Côte d’Ivoire"                       = "Côte d'Ivoire",
  "Congo (Kinshasa)"                    = "Democratic Republic of the Congo",
  "Congo, Dem. Rep."                    = "Democratic Republic of the Congo",
  "Congo (Brazzaville)"                 = "Republic of the Congo",
  "Congo, Rep."                         = "Republic of the Congo",
  "Burma"                               = "Myanmar",
  "Macedonia"                           = "North Macedonia",
  "Macedonia, FYR"                      = "North Macedonia",
  "Russian Federation"                  = "Russia",
  "Kyrgyz Republic"                     = "Kyrgyzstan",
  "Bahamas"                             = "The Bahamas",
  "Gambia"                              = "The Gambia",
  "Micronesia"                          = "Federated States of Micronesia",
  "Bolivia (Plurinational State of)"    = "Bolivia",
  "Venezuela (Bolivarian Republic of)"  = "Venezuela"
)

pred_long <- pred_long %>%
  mutate(
    host_country    = as.character(host_country),
    home_country_ne = dplyr::recode(home_country, !!!ne_recode, .default = home_country)
  )

# ---- World map ----
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(4326) %>%
  select(admin, geometry) %>%
  rename(ne_admin = admin) %>%
  filter(ne_admin != "Antarctica")

# =========================
# CDLP helper
# =========================
compute_slope_df <- function(fitted_model) {
  fx  <- fixef(fitted_model)["cult_fst_z"]
  re  <- ranef(fitted_model, condVar = TRUE)$host_country
  seA <- attr(ranef(fitted_model, condVar = TRUE)[["host_country"]], "postVar")
  se  <- sqrt(vapply(seq_len(dim(seA)[3]), function(i) seA[2,2,i], numeric(1)))
  tibble(host_country = rownames(re),
         slope = fx + re[,"cult_fst_z"],
         se = se) |>
    mutate(lower = slope - 1.96*se, upper = slope + 1.96*se,
           significant = !(lower <= 0 & upper >= 0),
           sig_color = ifelse(significant, "significant","not significant"))
}

# ---- Helper: Add "the" before United Kingdom ----
display_host_name <- function(country_name) {
  if (country_name == "United Kingdom") {
    return("the United Kingdom")
  } else {
    return(country_name)
  }
}


create_country_plot_sample <- function(country_name, fitted_model, outcome_var, title_suffix = "") {
  pts <- combined_data |> filter(host_country == country_name)
  ann <- compute_slope_df(fitted_model) |> filter(host_country == country_name)
  newdat <- expand.grid(
    cult_fst_z = seq(min(combined_data$cult_fst_z, na.rm=TRUE),
                     max(combined_data$cult_fst_z, na.rm=TRUE), length.out=100),
    host_country=country_name, pop_size_z=0, educ_level_z=0, GDPpc_z=0,
    tt_z=0, gen_dist_z=0, regional_ident="Africa", stringsAsFactors=FALSE
  )
  pi <- predictInterval(merMod=fitted_model, newdata=newdat, level=0.95, n.sims=500,
                        stat="median", type="linear.prediction", include.resid.var=FALSE)
  newdat <- bind_cols(newdat, pi); newdat$sig_color <- ann$sig_color
  
  ggplot() +
    geom_ribbon(data=newdat, aes(x=cult_fst_z, ymin=lwr, ymax=upr),
                fill="grey80", alpha=0.5) +
    geom_line(data=newdat, aes(x=cult_fst_z, y=fit, colour=sig_color), linewidth=1.2) +
    geom_point(data=pts, aes(x=cult_fst_z, y=.data[[outcome_var]]),
               colour="grey60", size=2, alpha=0.6) +
    ggrepel::geom_text_repel(data=pts, aes(x=cult_fst_z, y=.data[[outcome_var]], label=home_country),
                             size=2.5, max.overlaps=Inf, show.legend=FALSE) +
    geom_hline(yintercept=0, linetype="dashed", colour="black", linewidth=0.8, alpha=0.7) +
    geom_text(data=ann, aes(x=-Inf, y=Inf, label=paste0("β = ", round(slope, 2))),
              hjust=-0.1, vjust=1.1, size=5, fontface="bold",
              colour = ifelse(ann$significant, "#d7301f", "#2b8cbe")) +
    scale_colour_manual(values = c(significant = "#d7301f", `not significant` = "#2b8cbe")) +
    labs(
      title = paste0("Cultural Distance Labor Penalty for Migrants in ", display_host_name(country_name), title_suffix),
      x = "Cultural Distance (standardized)",
      y = "Migrant Employment Gap (p.p)",
      caption = paste(
        "Red = significant at p < .05; Blue = not significant.",
        "Shaded band shows 95% CI of model predictions.",
        "Controls: population size, education level, GDP per capita, migrant population size, regional identity, genetic distance.",
        "Dashed line = native-born employment rate (gap = 0).",
        sep = "\n"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title    = element_text(size = 12, face = "bold"),
      axis.text     = element_text(size = 11),
      legend.position = "none",
      plot.caption  = element_text(size = 9, hjust = 0.5)
    )
}

# =========================
# UI
# =========================
ui <- fluidPage(
  div(style="background:linear-gradient(to right,#2b8cbe,#1b3b73); color:white; text-align:center;
             padding:15px 10px; border-radius:6px; margin-bottom:25px;",
      h1("Migrant Employment Gaps Across European Destination Countries by Home Country",
         style="margin:0; font-weight:700; font-size:24px;")),
  # Global polish
  tags$style(HTML("
  body {
    background-color: #f5f5f5;   /* light grey background */
    font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
    color: #222;
  }

  .well, .panel, .main-panel {
    background-color: #ffffff;    /* keep panels white for contrast */
    border-radius: 8px;
    box-shadow: 0 1px 6px rgba(0,0,0,0.08);
    padding: 10px;
  }

  .leaflet {
    border-radius: 8px;
    overflow: hidden;
    box-shadow: 0 2px 10px rgba(0,0,0,0.06);
  }

  .btn, .btn-default, .btn-primary {
    background-color: #1b3b73;
    color: #fff;
    border: none;
    border-radius: 6px;
  }

  .btn:hover { background-color: #204390; }
")),
  sidebarLayout(
    sidebarPanel(width=3,
                 selectInput(
                   "host",
                   "Host country:",
                   choices = sort(unique(pred_long$host_country)),
                   selected = "United Kingdom"
                 ),
                 selectInput("sample","Sample:", choices=c("Full Sample","Male Only","Female Only"), selected="Full Sample"),
                 selectInput("view","Show:",
                             choices=c("Predicted"="pred","Observed"="obs",
                                       "Difference between Predicted and Observed"="discrepancy"),
                             selected="pred"),
                 checkboxInput("show_labels","Overlay numbers on map",FALSE),
                 checkboxInput("show_caption","Show caption note",TRUE),
                 checkboxInput("show_cdlp","Show Cultural Distance Labor Penalty plot",TRUE),
                 downloadButton("download_map","Download PNG")
    ),
    mainPanel(width=9,
              h3(textOutput("mapTitleText"), style="text-align:center; margin-bottom:4px;"),
              div(textOutput("mapSubtitleText"), style="text-align:center; color:#555;"),
              leafletOutput("mapLeaf", height="600px"),
              # Map caption
              div(textOutput("mapCaptionText"),
                  style="text-align:center; color:#666; font-size:12px; margin-top:8px;"),
              # RMSE only for discrepancy
              conditionalPanel(
                condition = "input.view == 'discrepancy'",
                div(style = "text-align:center; font-weight:600; margin-top:6px;",
                    textOutput("rmseText"))
              ),
              # Optional CDLP plot + caption
              conditionalPanel(
                condition = "input.show_cdlp == true",
                plotOutput("cdlpPlot", height="600px"),
                div(textOutput("cdlpCaptionText"),
                    style="text-align:center; color:#666; font-size:12px; margin-top:8px;")
              ),
              
              # ---- Dashboard footer ----
              div(
                style = "text-align:center; color:#888; font-size:11px; margin-top:25px; margin-bottom:4px;",
                HTML("Dashboard assembled by <b>Harry Quantrill</b> for <i>The Cultural Distance Labor Penalty: Cultural Distance is Associated With Migrant Employment Gaps Across Europe</i>.")
              ),
              div(
                style = "text-align:center; color:#888; font-size:11px; margin-bottom:10px;",
                HTML("Employment data were taken from Eurostat's 2021 EU database, available for download here: 
        <a href='https://ec.europa.eu/CensusHub2/query.do?step=selectHyperCube&qhc=false' target='_blank'>
        https://ec.europa.eu/CensusHub2/query.do?step=selectHyperCube&qhc=false</a>.")
              )
              
    )
  )
)

# =========================
# SERVER
# =========================
server <- function(input, output, session) {
  pal_blue <- "#2b8cbe"
  pal_red  <- "#d7301f"
  na_grey  <- "#d9d9d9"  # lighter grey for NA
  
  view_label <- reactive({
    switch(input$view,
           pred="Predicted Migrant Employment Gap in ",
           obs="Observed Migrant Employment Gap in ",
           discrepancy="Difference (Predicted − Observed) in Migrant Employment Gap in ")
  })
  output$mapTitleText <- renderText({
    paste0(view_label(), display_host_name(input$host))
  })
  output$mapSubtitleText <- renderText({ paste0(input$sample, " by Home Country") })
  
  data_sel <- reactive({
    pred_long %>% filter(host_country == input$host, sample == input$sample)
  })
  
  # Build map data & palette
  map_data <- reactive({
    df <- data_sel()
    m <- world %>%
      left_join(df, by=c("ne_admin"="home_country_ne")) %>%
      mutate(fillvar = case_when(
        input$view=="pred"~pred,
        input$view=="obs"~obs,
        input$view=="discrepancy"~discrepancy,
        TRUE~NA_real_))
    
    # symmetric limits so 0 is white across all views
    rng <- suppressWarnings(max(abs(m$fillvar), na.rm=TRUE)); if(!is.finite(rng)||rng==0) rng<-1
    pal <- colorNumeric(c(pal_blue, "white", pal_red), domain = c(-rng, rng), na.color = na_grey)
    
    m$tooltip <- sprintf("<strong>Home:</strong> %s<br><strong>Gap (p.p.):</strong> %s",
                         m$ne_admin,
                         ifelse(is.na(m$fillvar),"NA",sprintf("%.1f",m$fillvar)))
    list(sf = m, pal = pal)
  })
  
  output$mapLeaf <- renderLeaflet({
    leaflet(options = leafletOptions(zoomControl = TRUE, minZoom = 1.2)) %>%
      addProviderTiles(providers$CartoDB.VoyagerNoLabels) %>%
      setView(lng=10, lat=35, zoom=2)
  })
  
  observe({
    md <- map_data()
    m  <- md$sf
    pal <- md$pal
    
    # centroids for number overlay (exclude NA)
    m_valid <- m %>% filter(!is.na(fillvar))
    ct <- NULL
    if (nrow(m_valid) > 0) {
      cent_sf <- suppressWarnings(st_point_on_surface(m_valid))
      coords  <- as.data.frame(st_coordinates(cent_sf))
      names(coords) <- c("lon","lat")
      ct <- cbind(st_drop_geometry(m_valid), coords)
    }
    
    vals_non_na <- m$fillvar[!is.na(m$fillvar)]
    
    leafletProxy("mapLeaf", data = m) %>%
      clearShapes() %>% clearControls() %>% clearMarkers() %>%
      addPolygons(
        fillColor   = ~pal(fillvar),
        weight      = 0.5, color = "#444",
        fillOpacity = 0.9,
        label       = lapply(m$tooltip, HTML),
        labelOptions = labelOptions(style = list("font-size"="11px"))
      ) %>%
      # main gradient legend WITHOUT NA
      addLegend(
        position = "bottomright",
        pal = pal, values = vals_non_na,
        title = HTML("Gap (p.p.)"),
        labFormat = labelFormat(digits = 1),
        opacity = 0.9, na.label = NULL
      ) %>%
      # separate NA swatch so it never overlaps the gradient
      addLegend(
        position = "bottomright",
        colors = na_grey, labels = "No data", opacity = 0.9
      )
    
    # number overlay (non-NA only)
    if (isTRUE(input$show_labels) && !is.null(ct) && nrow(ct) > 0) {
      leafletProxy("mapLeaf") %>%
        addLabelOnlyMarkers(
          data = ct, lng = ~lon, lat = ~lat,
          label = ~sprintf("%.1f", fillvar),
          labelOptions = labelOptions(
            noHide = TRUE, direction = "center", textOnly = TRUE,
            style = list(
              "color" = "#111111",
              "font-weight" = "500",  # lighter than before
              "font-size" = "10px",   # slightly larger, still tidy
              "text-shadow" = "0 0 2px white, 0 0 4px white"
            )
          )
        )
    }
  })
  
  # Map caption
  output$mapCaptionText <- renderText({
    if(!isTRUE(input$show_caption)) return(NULL)
    if(input$view=="discrepancy")
      "Note: Positive values mean the model overestimates the gap."
    else if(input$view=="obs")
      "Migrant employment gaps refer to the difference in the employment rates of migrants from a given home country and the employment rate of the native-born population. A gap of zero therefore indicates parity with the native-born population."
    else
      "Estimates come from models using cultural and genetic distance between home and host country, home country population size, education level, GDP per capita, continent, and migrant population size."
  })
  
  # RMSE (discrepancy only)
  output$rmseText <- renderText({
    req(input$view == "discrepancy")
    df <- data_sel()
    if (nrow(df) == 0) return(NULL)
    rmse <- sqrt(mean((df$pred - df$obs)^2, na.rm = TRUE))
    paste0("RMSE (Predicted − Observed) — ", input$host, " / ", input$sample, ": ",
           sprintf('%.2f', rmse), " p.p.")
  })
  
  # CDLP plot
  output$cdlpPlot <- renderPlot({
    req(input$show_cdlp)
    if(input$sample=="Full Sample")
      create_country_plot_sample(input$host, random_slope_model, "ptunemp")
    else if(input$sample=="Male Only")
      create_country_plot_sample(input$host, random_slope_model_men_only, "pmunemp", " (Male Only)")
    else
      create_country_plot_sample(input$host, random_slope_model_women_only, "pfunemp", " (Female Only)")
  })
  
  # =========================
  # DOWNLOAD PNG
  # =========================
  output$download_map <- downloadHandler(
    filename = function() {
      paste0(
        "map_", gsub(" ", "_", tolower(input$host)), "_",
        gsub(" ", "_", tolower(input$sample)),
        if (isTRUE(input$show_cdlp)) "_with_cdlp", ".png"
      )
    },
    content = function(file) {
      pal_blue <- "#2b8cbe"; pal_red <- "#d7301f"; na_grey <- "#d9d9d9"
      
      # rebuild the mapped sf exactly like the live map
      df <- pred_long %>% filter(host_country == input$host, sample == input$sample)
      m <- world %>%
        left_join(df, by = c("ne_admin" = "home_country_ne")) %>%
        mutate(
          fillvar = case_when(
            input$view == "pred"        ~ pred,
            input$view == "obs"         ~ obs,
            input$view == "discrepancy" ~ discrepancy,
            TRUE                        ~ NA_real_
          )
        )
      
      rng <- suppressWarnings(max(abs(m$fillvar), na.rm = TRUE)); if (!is.finite(rng) || rng == 0) rng <- 1
      
      # caption text (add RMSE for discrepancy)
      cap_text <- NULL
      if (isTRUE(input$show_caption)) {
        if (input$view == "discrepancy") {
          rmse_val <- sqrt(mean((df$pred - df$obs)^2, na.rm = TRUE))
          cap_text <- sprintf("Note: Positive values mean the model overestimates the gap. — RMSE: %.2f p.p.", rmse_val)
        } else if (input$view == "obs") {
          cap_text <- paste(
            "Migrant employment gaps refer to the difference in the employment rates of migrants",
            "from a given home country and the employment rate of the native-born population.",
            "A gap of zero therefore indicates parity with the native-born population.",
            sep = "\n"
          )
        } else {
          cap_text <- paste(
            "Estimates come from models using cultural and genetic distance between home and host country,",
            "home country population size, average education level, GDP per capita, continent,",
            "and host country migrant population size.",
            sep = "\n"
          )
        }
      }
      
      p_map <- ggplot(m) +
        geom_sf(aes(fill = fillvar), color = "#666666", linewidth = 0.3) +  # cleaner borders
        scale_fill_gradient2(
          low = pal_blue, mid = "white", high = pal_red,
          midpoint = 0, limits = c(-rng, rng),
          na.value = na_grey, name = "Gap (p.p.)"
        ) +
        labs(
          title    = paste0(view_label(), input$host),
          subtitle = paste0(input$sample, " by Home Country"),
          caption  = cap_text
        ) +
        theme_minimal(base_size = 13) +  # match CDLP base size
        theme(
          axis.text     = element_blank(),
          axis.ticks    = element_blank(),
          axis.title    = element_blank(),
          panel.grid    = element_blank(),
          plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          legend.title  = element_text(size = 10),
          legend.text   = element_text(size = 9),
          plot.caption  = element_text(size = 9, hjust = 0.5)
        )
      
      # optional numbers on PNG (non-NA only)
      if (isTRUE(input$show_labels)) {
        m_valid <- m %>% filter(!is.na(fillvar))
        if (nrow(m_valid) > 0) {
          cent_sf <- suppressWarnings(st_point_on_surface(m_valid))
          coords  <- as.data.frame(st_coordinates(cent_sf))
          names(coords) <- c("lon", "lat")
          ct <- cbind(st_drop_geometry(m_valid), coords)
          label_col <- "black"
          p_map <- p_map +
            geom_text(
              data = ct,
              aes(lon, lat, label = sprintf('%.1f', fillvar)),
              size = 1.8, color = label_col, inherit.aes = FALSE
            )
        }
      }
      
      if (isTRUE(input$show_cdlp)) {
        p_cdlp <- if (input$sample == "Full Sample") {
          create_country_plot_sample(input$host, random_slope_model, "ptunemp")
        } else if (input$sample == "Male Only") {
          create_country_plot_sample(input$host, random_slope_model_men_only, "pmunemp", " (Male Only)")
        } else {
          create_country_plot_sample(input$host, random_slope_model_women_only, "pfunemp", " (Female Only)")
        }
        combined <- (p_map / p_cdlp) + plot_layout(heights = c(0.55, 0.45))
        ggsave(file, combined, width = 12, height = 14, dpi = 300)
      } else {
        ggsave(file, p_map, width = 11, height = 7.5, dpi = 300)
      }
    }
  )
}

# =========================
# RUN
# =========================
shinyApp(ui, server)
