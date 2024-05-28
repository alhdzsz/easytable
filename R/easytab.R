
easy_table <- function(model_list,
                       csv = NULL,
                       robust.se = F,
                       control.var = NULL,
                       margins = F,
                       highlight = F) {

  # Dependencies
  if (!requireNamespace("marginaleffects", quietly = TRUE)) {
    install.packages("marginaleffects")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (!requireNamespace("flextable", quietly = TRUE)) {
    install.packages("flextable")
  }
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    install.packages("lmtest")
  }
  if (!requireNamespace("broom", quietly = TRUE)) {
    install.packages("broom")
  }
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    install.packages("sandwich")
  }

  require(broom)
  require(dplyr)
  require(lmtest)
  require(flextable)
  require(sandwich)
  require(marginaleffects)

  # Error messages
  if (!is.list(model_list) || is.null(names(model_list))) {
    stop("Input must be a named list of models. Each element of the list must be a statistical model object.")
  }

  # Function
  mnames <- names(model_list)

  parse_model <- function(model) {
    model_name <- deparse(substitute(model))

    if(robust.se == T & margins == F){
      m <- lmtest::coeftest(model, vcov = sandwich::vcovHC(model, type = "HC")) %>% tidy()
    }
    if(robust.se == F & margins == T){
      m <- marginaleffects::marginaleffects(model) %>% tidy()
    }
    if(robust.se == T & margins == T){
      m1 <- lmtest::coeftest(model, vcov = sandwich::vcovHC(model, type = "HC")) %>%
        tidy() %>%
        filter(term != "(Intercept)") %>%
        select(-estimate)
      m2 <- marginaleffects::marginaleffects(model) %>%
        tidy() %>%
        select(term, estimate)
      m <- left_join(m1,m2)
    }
    if(robust.se == F & margins == F){
      m <- model %>% tidy()
    }

    mod.df <- m %>%
      select(term, estimate, std.error, p.value) %>%
      mutate(
        significance = case_when(
          p.value < 0.01 ~ "***",
          p.value >= 0.01 & p.value <= 0.05 ~ "** ",
          p.value > 0.05 & p.value <= 0.1  ~ "*  ",
          TRUE ~ "   "
        )) %>%
      mutate(across(where(is.numeric), round, digits = 2)) %>%
      mutate(estimate = paste0(estimate, " ",
                               significance, "\n", "(",
                               std.error, ")")) %>%
      select(term, estimate) %>%
      rename(model = estimate)

    measures <-
      data.frame(
        term = c("N", "R sq.", "Adj. R sq.", "AIC"),
        model = c(
          stats::nobs(model),
          ifelse(!is.null(summary(model)$r.squared),summary(model)$r.squared,NA),
          ifelse(!is.null(summary(model)$adj.r.squared),summary(model)$adj.r.squared,NA),
          ifelse(!is.null(summary(model)$aic),summary(model)$aic,NA)
        )
      )

    measures$model <- round(measures$model, 2)

    measures$model <- as.character(measures$model)

    mod.df <- bind_rows(mod.df, measures)

    return(mod.df)

  }

  mtable <- parse_model(model_list[[1]])
  names(mtable)[2] <- mnames[1]

  if(length(model_list) > 1){
    for(i in 2:length(model_list)){
      mod <- model_list[[i]]
      mtableadd <- parse_model(mod)
      names(mtableadd)[2] <- mnames[i]
      mtable <- full_join(mtable,mtableadd)
    }
  }

  if(!is.null(control.var)){
    for(var in control.var){
      mtable$term <- gsub(sprintf("^factor\\(%s\\).*$", var), var, mtable$term)
      mtable$term <- gsub(sprintf("^log\\(%s\\).*$", var), var, mtable$term)
      mtable$term <- gsub(sprintf("^%s.*$", var), var, mtable$term)
    }
  }

  for (col in names(mtable)[-1]) {
    replace_indices <- mtable$term %in% control.var & !is.na(mtable[[col]])
    mtable[replace_indices, col] <- "Y"
  }

  # This is a slight problem!
  for (col in names(mtable)[-1]) {
    mtable <- mtable[!(duplicated(mtable$term) & mtable[[col]] == "Y"), ]
  }

  mtable[is.na(mtable)] <- ""
  mtable <- mtable[order(apply(mtable[, -1], 1, function(row) sum(row == "Y")), decreasing = F),]

  mmes <- mtable %>% filter(term %in% c("N", "R sq.", "Adj. R sq.", "AIC"))
  mmes <- mmes[rowSums(mmes[-1] != "") > 0, ]
  mtable <- mtable %>% filter(!term %in% c("N", "R sq.", "Adj. R sq.", "AIC"))
  mtable <- bind_rows(mtable, mmes)

  ft <- flextable(mtable) %>%
    add_footer_lines("Significance: ***p < .01; **p < .05; *p < .1 ") %>%
    hline(j = 1:ncol(mtable), i = nrow(mtable)-nrow(mmes))

  if(robust.se == T & margins == F){
    ft <- ft %>%
      add_footer_lines("Note: Robust Standard Errors")
  }
  if(robust.se == F & margins == T){
    ft <- ft %>%
      add_footer_lines("Note: Average Marginal Effects (AME)")
  }
  if(robust.se == T & margins == T){
    ft <- ft %>%
      add_footer_lines("Note: Marginal Effects and Robust Standard Errors")
  }

  if(highlight){
    for(i in 2:ncol(mtable)) {
      p_values <- grepl("\\*", mtable[[i]])
      ft <- ft %>%
        bg(j = i,
           i = p_values,
           part = "body", bg = "#e6ffe6")
    }
    for(i in 2:ncol(mtable)) {
      p_values <- grepl("-\\d+(\\.\\d+)? \\*", mtable[[i]])
      ft <- ft %>%
        bg(j = i,
           i = p_values,
           part = "body", bg = "#ffcccc")
    }
  }

  if(!is.null(csv)){
    write.csv(
      mtable,
      file = paste0(csv, ".csv"),
      row.names = F,
    )
  }
  return(ft)
}

