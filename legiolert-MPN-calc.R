# app.R

library(shiny)
library(bslib)
library(DT)


# -------------------------------------------------------------------------
# Legiolert MPN calculator
#
# Statistical approach:
# Jarvis B, Wilrich C, Wilrich P-T. (2010).
# Reconsideration of the derivation of Most Probable Numbers, their
# standard deviations, confidence bounds and rarity values.
# Journal of Applied Microbiology 109(5), 1660-1667.
# doi:10.1111/j.1365-2672.2010.04792.x
#
# Quanti-Tray/Legiolert:
# 90 small wells and 6 large wells.
# Well volumes were back-calculated by modifying the small well volume
# until calculated MPNs matched the official Legiolert MPN table.
# -------------------------------------------------------------------------


mpn_calc <- function(d, x) {
  
  # Well volumes (mL)
  w <- c(0.198, 13.69666667)
  n <- c(90, 6)
  
  if (length(x) != 2 || any(!is.finite(x))) {
    stop("Enter positive counts for the 90 small and 6 large wells.")
  }
  
  if (any(x < 0) ||
      x[1] > 90 ||
      x[2] > 6 ||
      any(x != floor(x))) {
    stop("Positive-well counts must be whole numbers within the tray limits.")
  }
  
  if (!is.finite(d) || d <= 0) {
    stop("Dilution factor must be greater than 0.")
  }
  
  
  # -----------------------------------------------------------------------
  # ALL WELLS NEGATIVE
  #
  # Jarvis:
  # MPN estimate = 0
  #
  # Upper 95% limit:
  # ln(40) / sum(d_i * w_i * n_i)
  #
  # This is then converted to MPN / 100 mL.
  # -----------------------------------------------------------------------
  
  if (sum(x) == 0) {
    
    upper95 <- log(40) / sum(d * w * n)
    upper95_100 <- upper95 * 100
    
    return(data.frame(
      MPN_100mL = "<1",
      log10_MPN = NA_real_,
      Lower95_100mL = 0,
      Upper95_100mL = upper95_100,
      SD_log10_MPN = NA_real_,
      Rarity_Index = 1,
      Status = paste0(
        "No positive wells detected; MPN estimate = <1 without accounting for dilution. ",
        "Upper 95% limit shown."
      )
    ))
  }
  
  
  # -----------------------------------------------------------------------
  # ALL WELLS POSITIVE
  #
  # No finite MLE exists when every well is positive.
  # -----------------------------------------------------------------------
  
  if (all(x == n)) {
    
    lower95 <- NA_real_
    
    # Lower one-sided 95% limit:
    # solve likelihood probability = 0.025
    lower_eq <- function(mu) {
      sum(
        n * log1p(-exp(-d * w * mu))
      ) - log(0.025)
    }
    
    lower_mu <- tryCatch(
      uniroot(
        lower_eq,
        interval = c(1e-12, 1e10)
      )$root,
      error = function(e) NA_real_
    )
    
    lower95_100 <- ifelse(
      is.finite(lower_mu),
      lower_mu * 100,
      NA_real_
    )
    
    return(data.frame(
      MPN_100mL = Inf,
      log10_MPN = Inf,
      Lower95_100mL = lower95_100,
      Upper95_100mL = Inf,
      SD_log10_MPN = NA_real_,
      Rarity_Index = 1,
      Status = "All wells positive; no finite MPN estimate."
    ))
  }
  
  
  # -----------------------------------------------------------------------
  # MPN ESTIMATE
  # -----------------------------------------------------------------------
  
  mpn_eq <- function(mu) {
    
    sum(
      (x * d * w) /
        (1 - exp(-d * w * mu)) -
        n * d * w
    )
  }
  
  
  lower <- 1e-10
  upper <- 1e10
  
  root <- tryCatch(
    uniroot(
      mpn_eq,
      interval = c(lower, upper)
    )$root,
    error = function(e) NA_real_
  )
  
  if (!is.finite(root)) {
    stop("The MPN could not be estimated for this result.")
  }
  
  mu_hat <- root
  
  
  # -----------------------------------------------------------------------
  # VARIANCE / SD OF LOG10 MPN
  # -----------------------------------------------------------------------
  
  var_ln_mu <- 1 / (
    mu_hat^2 *
      sum(
        x * (d * w)^2 *
          exp(-d * w * mu_hat) /
          (1 - exp(-d * w * mu_hat))^2
      )
  )
  
  sd_ln_mu <- sqrt(var_ln_mu)
  
  sd_log10_mpn <- log10(exp(1)) * sd_ln_mu
  
  
  # -----------------------------------------------------------------------
  # CONVERT TO MPN / 100 mL
  # -----------------------------------------------------------------------
  
  mpn_100 <- mu_hat * 100
  
  log10_mpn <- log10(mpn_100)
  
  
  # -----------------------------------------------------------------------
  # APPROXIMATE 95% CONFIDENCE LIMITS
  #
  # Jarvis:
  # Lower = MPN * exp(-2 * SD_ln)
  # Upper = MPN * exp(+2 * SD_ln)
  # -----------------------------------------------------------------------
  
  lower95 <- mpn_100 * exp(-2 * sd_ln_mu)
  upper95 <- mpn_100 * exp(2 * sd_ln_mu)
  
  
  # -----------------------------------------------------------------------
  # RARITY INDEX
  #
  # r = likelihood of observed result /
  #     maximum likelihood over all possible result combinations
  #
  # Calculated on the log scale for numerical stability.
  # -----------------------------------------------------------------------
  
  p <- 1 - exp(-d * w * mu_hat)
  q <- exp(-d * w * mu_hat)
  
  # Log likelihood of observed result
  logL_obs <- sum(
    lchoose(n, x) +
      x * log(p) +
      (n - x) * log(q)
  )
  
  # Find maximum likelihood among all possible positive-well combinations.
  logL_max <- -Inf
  
  for (x1 in 0:n[1]) {
    for (x2 in 0:n[2]) {
      
      xx <- c(x1, x2)
      
      logL <- sum(
        lchoose(n, xx) +
          xx * log(p) +
          (n - xx) * log(q)
      )
      
      if (logL > logL_max) {
        logL_max <- logL
      }
    }
  }
  
  rarity_index <- exp(logL_obs - logL_max)
  
  
  # -----------------------------------------------------------------------
  # STATUS
  # -----------------------------------------------------------------------
  
  status <- if (mpn_100 > 2272.6) {
    ">2,272.6 MPN/100 mL (above standard undiluted tray range)"
  } else {
    ""
  }
  
  
  data.frame(
    MPN_100mL = mpn_100,
    log10_MPN = log10_mpn,
    Lower95_100mL = lower95,
    Upper95_100mL = upper95,
    SD_log10_MPN = sd_log10_mpn,
    Rarity_Index = rarity_index,
    Status = status
  )
}


# -------------------------------------------------------------------------
# Formatting
# -------------------------------------------------------------------------

fmt_mpn <- function(x) {
  
  if (is.na(x)) return("—")
  
  if (is.infinite(x)) return("∞")
  
  if (x == 0) return("0")
  
  if (x < 10) {
    return(formatC(x, format = "f", digits = 2))
  }
  
  if (x < 100) {
    return(formatC(x, format = "f", digits = 1))
  }
  
  formatC(
    x,
    format = "f",
    digits = 0,
    big.mark = ","
  )
}


fmt_decimal <- function(x, digits = 3) {
  
  if (is.na(x)) return("—")
  
  if (is.infinite(x)) return("∞")
  
  formatC(
    x,
    format = "f",
    digits = digits
  )
}


fmt_rarity <- function(x) {
  
  if (is.na(x)) return("—")
  
  if (x == 0) return("<0.0001")
  
  formatC(
    x,
    format = "f",
    digits = 4
  )
}



# -------------------------------------------------------------------------
# UI
# -------------------------------------------------------------------------

ui <- page_fillable(
  
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#111111",
    font_scale = 0.95
  ),
  
  tags$head(
    
    tags$style(HTML("

       body { background: #f7f7f5; color: #171717; }
      .app-shell { max-width: 900px; margin: 0 auto; padding: 42px 24px 36px; }
      .eyebrow { font-size: 0.75rem; letter-spacing: .14em; text-transform: uppercase;
                 color: #777; font-weight: 600; margin-bottom: 8px; }
      h1 { font-weight: 600; letter-spacing: -.035em; margin-bottom: 8px; }
      .subtitle { color: #666; margin-bottom: 32px; }
      .card-clean { background: white; border: 1px solid #e6e6e3; border-radius: 14px;
                    padding: 24px; margin-bottom: 16px; box-shadow: 0 1px 2px rgba(0,0,0,.02); }
      .section-label { font-size: .76rem; font-weight: 700; letter-spacing: .08em;
                       text-transform: uppercase; color: #777; margin-bottom: 14px; }
      .result-number { font-size: 2.35rem; font-weight: 650; letter-spacing: -.04em; }
      .result-unit { color: #777; font-size: .9rem; margin-top: -4px; }
      .ci-text { color: #444; font-size: .95rem; }
      .status { border-top: 1px solid #eee; padding-top: 14px; margin-top: 18px;
                color: #666; font-size: .86rem; }
      .btn-dark { background:#171717; border-color:#171717; border-radius:9px; }
      .btn-dark:hover { background:#333; border-color:#333; }
      .form-control, .form-select { border-radius: 9px; border-color:#d9d9d5; }
      .form-control:focus, .form-select:focus { border-color:#999; box-shadow: 0 0 0 .15rem rgba(0,0,0,.07); }
      .help-note { color:#777; font-size:.82rem; line-height:1.5; }
      .footer { color:#888; font-size:.78rem; line-height:1.55; padding-top: 8px; }
      .shiny-output-error-validation { color: #a33; }

      /* Results */

      .results-card {
        background: #171717;
        color: white;
        border-radius: 11px;
        padding: 17px 20px;
        margin-bottom: 10px;
      }

      .results-card .section-label {
        color: #aaa;
      }

      .main-result {
        font-size: 2rem;
        font-weight: 650;
        letter-spacing: -.045em;
        line-height: 1;
      }

      .main-unit {
        color: #aaa;
        font-size: .72rem;
        margin-top: 4px;
      }

      .result-stat {
        border-left: 1px solid #3a3a3a;
        padding-left: 18px;
      }

      .result-stat-label {
        color: #aaa;
        font-size: .64rem;
        text-transform: uppercase;
        letter-spacing: .08em;
        margin-bottom: 3px;
      }

      .result-stat-value {
        font-size: 1rem;
        font-weight: 600;
      }

      .ci-label {
        color: #aaa;
        font-size: .64rem;
        text-transform: uppercase;
        letter-spacing: .08em;
        margin-bottom: 2px;
      }

      .ci-value {
        font-size: .85rem;
      }

      .status {
        color: #aaa;
        font-size: .68rem;
        margin-top: 9px;
        padding-top: 8px;
        border-top: 1px solid #333;
      }

     

    "))
  ),
  
  
  div(
    class = "app-shell",
    
    # ---------------------------------------------------------------------
    # HEADER
    # ---------------------------------------------------------------------
    
      h1("Legiolert MPN Calculator"),
    
    
    # ---------------------------------------------------------------------
    # INPUTS
    # ---------------------------------------------------------------------
    
    div(
      class = "card-clean",
      
      div(
        class = "section-label",
        "Tray result"
      ),
      
      fluidRow(
        
        column(
          width = 3,
          
          numericInput(
            "small_positive",
            "Positive small wells",
            value = 0,
            min = 0,
            max = 90,
            step = 1
          ),
          
          div(
            class = "help-note",
            "0–90 positive"
          )
        ),
        
        column(
          width = 3,
          
          numericInput(
            "large_positive",
            "Positive large wells",
            value = 0,
            min = 0,
            max = 6,
            step = 1
          ),
          
          div(
            class = "help-note",
            "0–6 positive"
          )
        ),
        
        column(
          width = 3,
          
          numericInput(
            "dilution",
            "Dilution factor",
            value = 1,
            min = 0.000001,
            step = 0.1
          ),
          
          div(
            class = "help-note",
            "1 = undiluted; 0.1 = 1:10, ect."
          )
        ),
        
        column(
          width = 3,
          
          div(
            style = "padding-top: 29px;",
            
            actionButton(
              "calculate",
              "Calculate MPN",
              class = "btn-dark",
              width = "100%"
            )
          )
        )
      )
    ),
    
    
    # ---------------------------------------------------------------------
    # RESULTS
    # ---------------------------------------------------------------------
    
    uiOutput("result_card"),
    
    
    # ---------------------------------------------------------------------
    # BOTTOM ROW
    # ---------------------------------------------------------------------
    
    fluidRow(
      class = "bottom-row",
      
      # Method / equation
      column(
        width = 12,
        
        div(
          class = "card-clean",
          
          div(
            class = "section-label",
            "Method & equations"
          ),
          # div(
          #   class = "equation-placeholder"
          #   ,
          #   tags$img(
          #     src = "Jarvis equations.png",
          #     class= "equation-image",
          #     alt ="Jarvis equations"
          #   )),
          div(
            class = "method-text",
            
            p(
              style = "margin-bottom: 7px;",
              
              "Calculation uses the MPN formulas ",
              "described by ",
              
              tags$a(
                "Jarvis, Wilrich & Wilrich (2010)",
                href = "https://doi.org/10.1111/j.1365-2672.2010.04792.x",
                target = "_blank"
              ),
              
              ". The confidence interval shown here use the approximate +/- SD interval on the natural log-scale,",
              "which differs slighitly from the IDEXX Water MPN Generator Software results. Well volumes were",
              "back-calculated by modifying the small well volumes until calculated MPNs matched the official MPN table:",
              tags$a(
                "Quanti-Tray/Legiolert MPN Table",
                href = "https://www.idexx.com/files/quanti-tray-legiolert-mpn-table.pdf",
                target = "_blank"
              ),
            ),
              
           
            
            
            div(
              style = "margin-top: 7px;",
              
              "Reference: Jarvis B, Wilrich C, Wilrich P-T. (2010). ",
              tags$em(
                "Journal of Applied Microbiology"
              ),
              " 109(5):1660–1667. ",
              "doi:10.1111/j.1365-2672.2010.04792.x."
            )
          )
        )
      ),
      
      
    #   # Download / reporting note
    #   column(
    #     width = 5,
    #     
    #     div(
    #       class = "card-clean bottom-card",
    #       
    #       div(
    #         class = "section-label",
    #    
    #       ),
    #       
    #       div(
    #         class = "method-text",
    # 
    #         downloadButton(
    #           "download",
    #           "Download result",
    #           class = "btn-dark",
    #           width = "100%"
    #         )
    #         
    #       )
    #     )
    #   )
    # ),
    ),
    
    div(
      class = "footer",
      
      HTML(
        "For research and laboratory data analysis. Verify the applicable ",
        "Legiolert method, dilution, reporting limits, and QC requirements ",
        "for your use case. ChatGPT assisted in the coding of this app."
      )
    )
  )
)


# -------------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------------

server <- function(input, output, session) {
  
  
  result <- eventReactive(
    input$calculate,
    
    {
      
      validate(
        
        need(
          input$small_positive >= 0 &&
            input$small_positive <= 90,
          "Small-well positives must be between 0 and 90."
        ),
        
        need(
          input$large_positive >= 0 &&
            input$large_positive <= 6,
          "Large-well positives must be between 0 and 6."
        ),
        
        need(
          input$dilution > 0,
          "Dilution factor must be greater than 0."
        )
      )
      
      
      tryCatch(
        
        mpn_calc(
          d = input$dilution,
          x = c(
            input$small_positive,
            input$large_positive
          )
        ),
        
        error = function(e) {
          
          validate(
            need(
              FALSE,
              conditionMessage(e)
            )
          )
        }
      )
      
    },
    
    ignoreInit = FALSE
  )
  
  
  # -----------------------------------------------------------------------
  # RESULT CARD
  # -----------------------------------------------------------------------
  
  output$result_card <- renderUI({
    
    r <- result()
    
    
    if (is.infinite(r$MPN_100mL)) {
      
      div(
        class = "results-card",
        
        div(
          class = "section-label",
          "Most Probable Number"
        ),
        
        fluidRow(
          
          column(
            width = 3,
            
            div(
              class = "main-result",
              "∞"
            ),
            
            div(
              class = "main-unit",
              "MPN / 100 mL"
            )
          ),
          
          column(
            width = 2,
            
            div(
              class = "result-stat",
              
              div(
                class = "result-stat-label",
                "log10 MPN"
              ),
              
              div(
                class = "result-stat-value",
                "∞"
              )
            )
          ),
          
          column(
            width = 3,
            
            div(
              class = "result-stat",
              
              div(
                class = "result-stat-label",
                "95% confidence limit"
              ),
              
              div(
                class = "ci-value",
                paste0(
                  "Lower: ",
                  fmt_mpn(r$Lower95_100mL)
                )
              ),
              
              div(
                class = "ci-value",
                paste0(
                  "Upper: ",
                  fmt_mpn(r$Upper95_100mL)
                )
              )
            )
          ),
          
          column(
            width = 2,
            
            div(
              class = "result-stat",
              
              div(
                class = "result-stat-label",
                "Rarity Index"
              ),
              
              div(
                class = "result-stat-value",
                fmt_rarity(r$Rarity_Index)
              )
            )
          ),
            # Download / reporting note
            column(
              width = 2,


                # div(
                #   class = "section-label",
                # 
                # ),

                div(
                  class = "method-text",

                  downloadButton(
                    "download",
                    "Download",
                    class = "btn-light",
                    width = "100%"
                  ))),
        ),
        
        div(
          class = "status",
          r$Status
        )
      )
      
      
    } else {
      
      div(
        class = "results-card",
        
        div(
          class = "section-label",
          "Estimated concentration"
        ),
        
        fluidRow(
          
          column(
            width = 4,
            
            div(
              class = "main-result",
              fmt_mpn(r$MPN_100mL)
            ),
            
            div(
              class = "main-unit",
              "MPN / 100 mL"
            )
          ),
          
          column(
            width = 2,
            
            div(
              class = "result-stat",
              
              div(
                class = "result-stat-label",
                "log10 MPN"
              ),
              
              div(
                class = "result-stat-value",
                fmt_decimal(
                  r$log10_MPN,
                  3
                )
              )
            )
          ),
          
          column(
            width = 2,
            
            div(
              class = "result-stat",
              
              div(
                class = "result-stat-label",
                "95% confidence limits"
              ),
              
              div(
                class = "ci-value",
                paste0(
                  "Lower: ",
                  fmt_mpn(r$Lower95_100mL)
                )
              ),
              
              div(
                class = "ci-value",
                paste0(
                  "Upper: ",
                  fmt_mpn(r$Upper95_100mL)
                )
              )
            )
          ),
          
          column(
            width = 2,
            
            div(
              class = "result-stat",
              
              div(
                class = "result-stat-label",
                "SD log10 MPN"
              ),
              
              div(
                class = "result-stat-value",
                fmt_decimal(
                  r$SD_log10_MPN,
                  3
                )
              ),
              
              div(
                class = "result-stat-label",
                style = "margin-top: 8px;",
                "Rarity Index"
              ),
              
              div(
                class = "result-stat-value",
                fmt_rarity(r$Rarity_Index)
              )
            )
          ),
          column(
            width = 2,
            
            
            # div(
            #   class = "section-label",
            # 
            # ),
            
            div(
              class = "method-text",
              
              downloadButton(
                "download",
                "Download",
                class = "btn-light",
                width = "100%"
              ))),
        ),
        
        div(
          class = "status",
          r$Status
        )
      )
    }
  })
  
  
  # -----------------------------------------------------------------------
  # DOWNLOAD
  # -----------------------------------------------------------------------
  
  output$download <- downloadHandler(
    
    filename = function() {
      
      paste0(
        "legiolert_mpn_",
        Sys.Date(),
        ".csv"
      )
    },
    
    content = function(file) {
      
      r <- result()
      
      out <- data.frame(
        
        small_wells_positive =
          input$small_positive,
        
        large_wells_positive =
          input$large_positive,
        
        dilution_factor =
          input$dilution,
        
        MPN_100mL =
          r$MPN_100mL,
        
        log10_MPN =
          r$log10_MPN,
        
        Lower95_100mL =
          r$Lower95_100mL,
        
        Upper95_100mL =
          r$Upper95_100mL,
        
        SD_log10_MPN =
          r$SD_log10_MPN,
        
        Rarity_Index =
          r$Rarity_Index,
        
        status =
          r$Status
      )
      
      write.csv(
        out,
        file,
        row.names = FALSE
      )
    }
  )
}


shinyApp(
  ui,
  server
)


