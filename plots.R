## Tooltips ====
font_large <- list(family='Spectral', size = 16, color = '#000000')
font_medium <- list(family='Spectral', size = 14, color = '#000000')
font_small <- list(family='Spectral', size = 10, color = '#000000')

tooltip <- list(
  bordercolor = "transparent",
  font = list(color = "#ffffff", family = font_large),
  align = 'center'
)

#scales::show_col(sypalette$hex)
male.color <- sypalette$hex[9]
female.color <- sypalette$hex[1]
other.color <- '#d9d9d9'
#' Diagnosis Date Histogram
#'
#' Function for Producing Histogram
#' @param df a dataframe object
#' @param xcol column to group and count by; used for xaxis
#' @param roll logical. flag to include rolling mean line. Defaults to `TRUE`
#' @param roll_period integer. period to use for rolling mean. Defaults to `14`
plot.histogram <-
  function(df,
           xcol = diagnosisdate,
           tooltip_desc = 'are aged',
           roll = TRUE,
           roll_period = 14) {
    xcol = enquo(xcol)
    
    plot.data <- df %>%
      count_patients(!!xcol, sourcename, row_total = FALSE) %>%
      arrange(!!xcol, sourcename) %>%
      ungroup %>%
      mutate(roll_mean = zoo::rollmean(patients, roll_period, na.pad = TRUE))
    
    plot <- plot.data %>%
      plot_ly(
        x = xcol,
        y = ~ patients,
        color = ~ sourcename,
        colors = rev(sort(sypalette$hex)),
        type = 'bar',
        marker = list(#color = sypalette$hex[[10]],
                      opacity = 0.90,
                      font = font_medium),
        text = ~ abs(patients),
        textposition = 'inside',
        hoverlabel = tooltip,
        hovertemplate = paste0(
          " <span style='font-size:1.2em'><span style='font-size:1.6em'><b>%{y:.0f}</b></span> patients", 
          glue::glue(tooltip_desc), 
          "<b>%{x}</b></span><extra></extra> "
        )
      ) %>%
      layout(
        barmode = 'stack',
        showlegend = TRUE,
        uniformtext = list(minsize = 18, mode = 'hide'),
        legend = list(
          font = syhelpr::font_medium,
          orientation = "h",
          xanchor = "center",
          x = 0.55,
          y = -0.045
        ),
        yaxis = list(title = 'patients',
                     yax.hist),
        xaxis = list(title = '',
                     xax.hist),
        margin = list(
          l = 5,
          r = 10,
          b = 20,
          t = 35,
          pad = 4
        )
      ) %>% config(displaylogo = FALSE)
    
    
    if (!is.factor(plot.data %>% pull(!!xcol))) {
      plot <- plot %>%
        layout(xaxis = list(
          range = c(
            plot.data %>% summarise(min = min(!!xcol) - 1) %>% pull,
            plot.data %>% summarise(max = max(!!xcol) + 1) %>% pull
          )
        ))
    }
    
    if (roll) {
      plot <- plot %>% add_trace(
        x =  xcol,
        y =  ~ roll_mean,
        text = xcol,
        name = paste(roll_period, "Rolling Average"),
        text = roll_period,
        type = 'scatter',
        mode = 'lines+markers',
        marker = list(opacity = 0, color = sypalette$hex[[1]]),
        line = list(
          shape = "spline",
          color = sypalette$hex[[1]],
          width = 5,
          dash = "dot"
        ),
        hoverinfo = text,
        hovertemplate = paste0(
          " <span style='font-size:1.2em'><span style='font-size:1.6em'><b>%{y:.0f}</b></span> ",
          roll_period,
          " rolling average on <b>%{x}</b></span><extra></extra> "
        )
      )
    }
    return(plot)
  }




#' Sex Date Histogram
#'
#' Function for Producing Histogram
#' @param df a dataframe object with
#' @param xcol column to group and count by; used for xaxis
#' @param roll logical. flag to include rolling mean line. Defaults to `TRUE`
#' @param roll_period integer. period to use for rolling mean. Defaults to `14`
plot.age_pyramid <-
  function(df,
           age_col = age_group,
           sex_col = sex) {
    sex_col = enquo(sex_col)
    age_col = enquo(age_col)
    
    female <-
      df %>% filter(!!sex_col == 'Female') %>% count_patients(!!age_col, row_total = FALSE)
    male <-
      df %>% filter(!!sex_col == 'Male') %>% count_patients(!!age_col, row_total = FALSE)
    male.color <- sypalette$hex[1]
    female.color <- sypalette$hex[9]
    other.color <- '#d9d9d9'
    
    plot_ly(
      data = female,
      y = age_col,
      x = ~ patients * -1,
      type = 'bar',
      marker = list(color = female.color),
      name = 'Female',
      text = ~ abs(patients),
      textposition = 'auto',
      hoverlabel = tooltip,
      hovertemplate = ~ paste(
        ' <span style="font-size: 1.5em; color: white"><span style="font-size: 1.25em"><b>%{text}</b></span>',
        'patients aged %{y}</b></span> '
      )
    ) %>%
      add_trace(
        data = male,
        y = age_col,
        x = ~ patients,
        type = 'bar',
        marker = list(color = male.color),
        name = 'Male'
      ) %>%
      layout(
        barmode = "overlay",
        showlegend = TRUE,
        uniformtext = list(minsize = 18, mode = 'hide'),
        legend = list(font =syhelpr::font_large,
                      x = 0.85,
                      y = 0.95),
        yaxis = yax.pyr,
        xaxis = xax.pyr,
        margin = list(
          l = 5,
          r = 10,
          b = 10,
          t = 10,
          pad = 4
        )
      )  %>%
      config(displaylogo = FALSE)
  }




#' Sex Date Histogram
#'
#' Function for Producing Histogramplot.age_pyramid <-
plot.sex_pie <- function(df, sex_col = sex) {
  sex_col = enquo(sex_col)
  
  plot.data <- df %>%
    mutate(!!sex_col := ifelse(sex %in% c('Male', 'Female'), sex, 'Unknown')) %>%
    count_patients(!!sex_col, row_total = F) %>%
    mutate(pct = patients / sum(patients))
  male.color <- sypalette$hex[1]
  female.color <- sypalette$hex[9]
  other.color <- '#d9d9d9'
  
  plot.data %>%
    plot_ly(
      labels = sex_col,
      values = ~ pct,
      text = ~ abs(patients),
      type = 'pie',
      textposition = 'inside',
      textinfo = 'label+percent',
      marker = list(
        colors = c(female.color, male.color, other.color),
        line = list(color = '#FFFFFF', width = 4)
      ),
      insidetextfont = list(color = '#FFFFFF'),
      hovertemplate = ~ glue::glue(
        ' <span style="font-size: 1.5em; color: white"><span style="font-size: 1.25em"><b>{abs(patients)}</b></span>',
        ' {ifelse(abs(patients)<2, \'patient\', \'patients\')}<extra>{sex}</extra> '
      ),
      hoverlabel = tooltip,
      hole = 0.5
    ) %>%
    layout(
      uniformtext = list(minsize = 18, mode = 'hide'),
      barmode = "overlay",
      showlegend = FALSE,
      legend = list(x = .75,
                    y = 0.95),
      margin = list(
        l = 5,
        r = 10,
        b = 10,
        t = 60,
        pad = 4
      )
    )  %>% config(displaylogo = FALSE)
}

#' Pie Chart
#'
#' Function for Producing Histogramplot.age_pyramid <-
# plot.pie <- function(df, col = sex, pal = syhelpr::sypalette$hex) {
#   col_ = enquo(col)
#   
#   plot.data <- df %>%
#     mutate(pie_col := !!col_) %>% 
#     count_patients(pie_col, row_total = F) %>%
#     mutate(pct = patients / sum(patients))
#  
#   plot.data %>%
#     plot_ly(
#       labels = ~ pie_col,
#       values = ~ pct,
#       text = ~ abs(patients),
#       type = 'pie',
#       textposition = 'inside',
#       textinfo = 'label+percent',
#       marker = list(
#         colors = pal,
#         line = list(color = '#FFFFFF', width = 4)
#       ),
#       insidetextfont = list(color = '#FFFFFF', size = 10),
#       hovertemplate = ~ glue::glue(
#         ' <span style="font-size: 1.5em; color: white"><span style="font-size: 1.25em"><b>{abs(patients)}</b></span>',
#         ' {ifelse(abs(patients)<2, \'patient\', \'patients\')}<extra>{pie_col}</extra> '
#       ),
#       hoverlabel = tooltip,
#       hole = 0.5
#     ) %>%
#     layout(
#       uniformtext = list(minsize = 12, mode = 'hide'),
#       barmode = "overlay",
#       showlegend = FALSE,
#       legend = list(x = .75,
#                     y = 0.95),
#       margin = list(
#         l = 5,
#         r = 10,
#         b = 10,
#         t = 60,
#         pad = 4
#       )
#     )  %>% config(displaylogo = FALSE)
# }

plot.pie <- function(df, col = sex, hole = 0.65, outline = 8, pal = sypalette$hex, scale = 1) {
  col_ = enquo(col)
  base_font <- 1.5 * scale
  
  plot.data <- df %>%
    mutate(group = replace_na(!!col_, FALSE)) %>% 
    count_patients(group, row_total = F) %>%
    mutate(pct = patients / sum(patients),
           pct_lab = paste0(round(100 * pct, 1), '%')) %>% 
    mutate(group = factor(group, levels = c(TRUE, FALSE))) %>% 
    arrange(group)
  
  plot.data %>%
    plot_ly(
      labels = ~ group,
      values = ~ pct,
      type = 'pie',
      rotation = ifelse(max(plot.data$pct, na.rm = T) >= 0.8, 90, 0),
      text = ~paste0('<span style="font-weight: 200; font-size:', base_font*.65, 'em;">', 
                     '<span style="font-weight: 200; font-size:', base_font, 'em;">', 
                     pct_lab, 
                     '<br></span>',
                     group, 
                     '</span>'),
      textposition = 'outside',
      textinfo = 'text',
      marker = list(
        colors = pal,
        line = list(color = '#FFFFFF', width = outline)
      ),
      hovertemplate = ~ glue::glue(
        ' <span style="font-size: 1.5em; color: white"><span style="font-size: 1.25em"><b>{format(patients, big.mark=",")}</b></span>',
        ' {ifelse(abs(patients)<2, \'patient\', \'patients\')}<extra>{group}</extra> '
      ),
      hoverlabel = tooltip,
      hole = hole
    ) %>%
    layout(
      uniformtext = list(minsize = 14, mode = 'hide'),
      font = font_medium,
      showlegend = FALSE,
      legend = list(x = .75,
                    y = 0.95),
      margin = list(l=60,r=60,t=60,b=60,pad=20)
    )  %>% config(displaylogo = FALSE)
}

#' Heatmap
#'
#' Function for Producing a Heatmap
plot.heatmap <- function(df) {
  heat.palette <-
    colorRampPalette(c("#e9f9f9", sypalette$hex[1], theme_blue))
  heatmap.tooltip <- list(
    bgcolor = heat.palette(20)[15],
    color = '#FFF',
    font = list(color = "#ffffff", family = "Roboto, sans-serif"),
    bordercolor = 'transparent'
  )
  
  heatmap.data <- cohort %>%
    count_patients(cancer_type, age_group, row_total = FALSE) %>%
    group_by(cancer_type) %>%
    mutate(pct = patients / sum(patients))
  
  plot_ly(
    data = heatmap.data,
    y = ~ cancer_type,
    x = ~ age_group,
    z = ~ patients,
    colors =  heat.palette(20),
    opacity = 0.90,
    text = ~ glue::glue(
      ' <span style="font-size: 1.5em"><span style="font-size: 1.25em"><b>{patients}</span> {cancer_type}</b> patients aged <b>{age_group}</b></span> '
    ),
    hoverinfo = 'text',
    hoverlabel = heatmap.tooltip,
    type = 'heatmap'
  ) %>%
    add_trace(
      data = heatmap.data,
      y = ~ cancer_type,
      x = ~ age_group,
      z = ~ pct,
      colors =  heat.palette(20),
      opacity = 0.90,
      text = ~ glue::glue(
        ' <span style="font-size: 1.5em"><span style="font-size: 1.25em"><b>{patients}</span> {cancer_type}</b> patients aged <b>{age_group}</b></span> '
      ),
      hoverinfo = 'text',
      hoverlabel = heatmap.tooltip,
      type = 'heatmap',
      visible = FALSE
    ) %>%
    layout(
      uniformtext = list(minsize = 18, mode = 'hide'),
      updatemenus = list(
        list(
          active = -1,
          type = 'buttons',
          direction = "right",
          xanchor = 'center',
          yanchor = "top",
          pad = list('r' = 0, 't' = 2, 'b' = 4),
          x = 0.5,
          y = 1.09,
          buttons = list(
            list(
              label = "Overall",
              method = "update",
              args = list(list(visible = c(TRUE, FALSE)))
            ),
            
            list(
              label = "By Cancer Type",
              method = "update",
              args = list(list(visible = c(FALSE, TRUE)))
            )
          )
        )
      ),
      yaxis = list(
        showticklabels = TRUE,
        title = '',
        showgrid = FALSE,
        fixedrange = TRUE,
        categoryarray = rev(unique(heatmap.data$cancer_type)),
        categoryorder = "array"
      ),
      xaxis = xax.hist,
      hovermode = 'closest',
      showlegend = FALSE,
      margin = list(
        l = 5,
        r = 10,
        b = 10,
        t = 40,
        pad = 4
      )
    )  %>% config(displaylogo = FALSE)
}
