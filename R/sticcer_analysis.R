source(file = "R/functions.R")
library(knitr)
library(kableExtra)

# Add * for sig and ~ small, ~~ meduim, and ~~~ large effect sizes
sigBetweenTwoValues_improved_timing <- function(side, normal, reduced, pvalue, effectsize,
                                         flipEval = FALSE) {
  value <- NULL
  
  ef = ifelse(effectsize == "large", "$^{\\ast}$", "")
  
  
  if(normal >= reduced & pvalue <= 0.05) {
    value <- paste(ef,"$\\APLup$\\textbf{", format(round(reduced, 2), nsmall = 2), "}", sep="")
  } else if(normal <= reduced & pvalue <= 0.05) {
    value <- paste(ef,"$\\APLdown$\\textbf{", format(round(reduced, 2), nsmall = 2), "}", sep="")
  } else {
    value <- as.character(format(round(reduced, 2), nsmall = 2))
  }
  
  return(value)
}

# Add * for sig and ~ small, ~~ meduim, and ~~~ large effect sizes
sigBetweenTwoValues_improved <- function(side, normal, reduced, pvalue, effectsize,
                                         flipEval = FALSE) {
  value <- NULL
  
  ef = ifelse(effectsize == "large", "$^{\\ast}$", "")
  
  if(normal >= reduced & pvalue <= 0.05) {
    value <- paste(ef,"$\\APLdown$\\textbf{", format(round(reduced, 1), nsmall = 1), "}", sep="")
  } else if(normal <= reduced & pvalue <= 0.05) {
    value <- paste(ef,"$\\APLup$\\textbf{", format(round(reduced, 1), nsmall = 1), "}", sep="")
  } else {
    value <- as.character(format(round(reduced, 1), nsmall = 1))
  }
  
  return(value)
}

sticcer_side_analysis_improved <- function(df, func = "mean") {
  
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduedwith = character(),
    metric = character(),
    reduce = character(),
    pvalue = character(), 
    effectsize = character(),
    sig = character()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducedusing <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        original <- df %>% filter(casestudy == case, 
                                  dbms == db, 
                                  datagenerator == gen, 
                                  reducedwith == "sticcer",
                                  fullreduce == TRUE)
        mutationanalysistime.median.normal <- eval(call(func, original$mutationanalysistime/1000))
        new_case <- case
        x <- c(db, new_case, gen, "sticcer", "mutationanalysistime",
               format(round(mutationanalysistime.median.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        originalresultstime.median.normal <- eval(call(func, original$originalresultstime/1000))
        x <- c(db, new_case, gen, "sticcer", "originalresultstime",
               format(round(originalresultstime.median.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        reductionTime.median.normal <- eval(call(func, original$reductionTime/1000))
        x <- c(db, new_case, gen, "sticcer", "reductionTime",
               format(round(reductionTime.median.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        # Total Time
        total_time.normal <- reductionTime.median.normal + mutationanalysistime.median.normal
        x <- c(db, new_case, gen, "sticcer", "totalTime",
               format(round(total_time.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        mutationScore.median.normal <- eval(call(func, original$mutationScore))
        x <- c(db, new_case, gen, "sticcer", "mutationScore",
               format(round(mutationScore.median.normal, 1), nsmall = 1),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        for (r in reducedusing) {
          if (r != "sticcer") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r)
            # mutationanalysistime
            mutationanalysistime.p <- wilcox.test(original$mutationanalysistime,
                                                  reduced$mutationanalysistime)$p.value
            mutationanalysistime.p <- ifelse(is.na(mutationanalysistime.p), 1, mutationanalysistime.p)
            
            mutationanalysistime.e <- effectsize_accurate(original$mutationanalysistime,
                                                          reduced$mutationanalysistime)$size
            
            mutationanalysistime.median.normal <- eval(call(func, original$mutationanalysistime/1000))
            mutationanalysistime.median.reduce <- eval(call(func, reduced$mutationanalysistime/1000))
            
            # originalresultstime
            originalresultstime.p <- wilcox.test(original$originalresultstime,
                                                  reduced$originalresultstime)$p.value
            originalresultstime.p <- ifelse(is.na(originalresultstime.p), 1, originalresultstime.p)
            
            originalresultstime.e <- effectsize_accurate(original$originalresultstime,
                                                          reduced$originalresultstime)$size
            
            originalresultstime.median.normal <- eval(call(func, original$originalresultstime/1000))
            originalresultstime.median.reduce <- eval(call(func, reduced$originalresultstime/1000))
            
            # reductionTime
            reductionTime.p <- wilcox.test(original$reductionTime,
                                                 reduced$reductionTime)$p.value
            reductionTime.p <- ifelse(is.na(reductionTime.p), 1, reductionTime.p)
            
            reductionTime.e <- effectsize_accurate(original$reductionTime,
                                                         reduced$reductionTime)$size
            
            #reductionTime.median.normal <- eval(call(func, original$reductionTime/1000))
            reductionTime.median.reduce <- eval(call(func, reduced$reductionTime/1000))
            #browser()
            # Total Time
            #*****
            total_time.p <- NULL
            total_time.e <- NULL
            if (r == "original") {
              sticcer_total_time <- original$reductionTime + original$mutationanalysistime
              total_time.p <- wilcox.test(reduced$mutationanalysistime,
                                          sticcer_total_time)$p.value
              total_time.p <- ifelse(is.na(reductionTime.p), 1, reductionTime.p)
              
              total_time.e <- effectsize_accurate(reduced$mutationanalysistime,
                                                  sticcer_total_time)$size
            }
            
            total_time.normal <- reductionTime.median.normal + mutationanalysistime.median.normal
            total_time.reduced <- reductionTime.median.reduce + mutationanalysistime.median.reduce
            
            # Mutation
            mutationScore.p <- wilcox.test(original$mutationScore,
                                           reduced$mutationScore)$p.value
            mutationScore.p <- ifelse(is.na(mutationScore.p), 1, mutationScore.p)
            
            mutationScore.e <- effectsize_accurate(original$mutationScore,
                                                   reduced$mutationScore)$size
            
            mutationScore.median.normal <- eval(call(func, original$mutationScore))
            mutationScore.median.reduce <- eval(call(func, reduced$mutationScore))
            
            # Reduction time
            #reductionTime <- (original %>% arrange(randomseed))$testgenerationtime - (reduced %>% arrange(randomseed))$testgenerationtime
            
            new_case <- case

            x <- c(db, new_case, gen, r, "mutationanalysistime",
                   sigBetweenTwoValues_improved_timing("R", mutationanalysistime.median.normal, 
                                               mutationanalysistime.median.reduce, mutationanalysistime.p,
                                               mutationanalysistime.e,
                                               flipEval = TRUE),
                   mutationanalysistime.p, mutationanalysistime.e,
                   ifelse(mutationanalysistime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))

            
            x <- c(db, new_case, gen, r, "originalresultstime",
                   sigBetweenTwoValues_improved_timing("R", originalresultstime.median.normal, 
                                                       originalresultstime.median.reduce, originalresultstime.p,
                                                       originalresultstime.e,
                                                       flipEval = TRUE),
                   originalresultstime.p, originalresultstime.e,
                   ifelse(originalresultstime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "totalTime",
                   format(round(total_time.reduced, 2), nsmall = 2),
                   "","","")
            results <- rbind(results, as.data.frame(t(x)))

            x <- c(db, new_case, gen, r, "reductionTime",
                   sigBetweenTwoValues_improved_timing("R", reductionTime.median.normal,
                                                       reductionTime.median.reduce, reductionTime.p,
                                                       reductionTime.e,
                                                       flipEval = TRUE),
                   reductionTime.p, reductionTime.e,
                   ifelse(reductionTime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "mutationScore",
                   sigBetweenTwoValues_improved("R", mutationScore.median.normal, 
                                       mutationScore.median.reduce, mutationScore.p,
                                       mutationScore.e),
                   mutationScore.p, mutationScore.e,
                   ifelse(mutationScore.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
 
          }
        }
      }
    }
  }
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducedusing", "metric", 
                         "reduce", "pvalue", "effectsize", "sig")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}


ots_side_analysis_improved <- function(df, func = "mean") {
  
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduedwith = character(),
    metric = character(),
    reduce = character(),
    pvalue = character(), 
    effectsize = character(),
    sig = character()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducedusing <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        original <- df %>% filter(casestudy == case, 
                                  dbms == db, 
                                  datagenerator == gen, 
                                  reducedwith == "original",
                                  fullreduce == FALSE)
        mutationanalysistime.median.normal <- eval(call(func, original$mutationanalysistime/1000))
        new_case <- case
        x <- c(db, new_case, gen, "original", "mutationanalysistime",
               format(round(mutationanalysistime.median.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        originalresultstime.median.normal <- eval(call(func, original$originalresultstime/1000))
        x <- c(db, new_case, gen, "original", "originalresultstime",
               format(round(originalresultstime.median.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        reductionTime.median.normal <- eval(call(func, original$reductionTime/1000))
        x <- c(db, new_case, gen, "original", "reductionTime",
               format(round(reductionTime.median.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        # Total Time
        total_time.normal <- reductionTime.median.normal + mutationanalysistime.median.normal
        x <- c(db, new_case, gen, "original", "totalTime",
               format(round(total_time.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        mutationScore.median.normal <- eval(call(func, original$mutationScore))
        x <- c(db, new_case, gen, "original", "mutationScore",
               format(round(mutationScore.median.normal, 1), nsmall = 1),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        for (r in reducedusing) {
          if (r != "original") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r)

            # mutationanalysistime
            mutationanalysistime.p <- wilcox.test(original$mutationanalysistime,
                                                  reduced$mutationanalysistime)$p.value
            mutationanalysistime.p <- ifelse(is.na(mutationanalysistime.p), 1, mutationanalysistime.p)
            
            mutationanalysistime.e <- effectsize_accurate(original$mutationanalysistime,
                                                          reduced$mutationanalysistime)$size
            
            mutationanalysistime.median.normal <- eval(call(func, original$mutationanalysistime/1000))
            mutationanalysistime.median.reduce <- eval(call(func, reduced$mutationanalysistime/1000))
            
            # originalresultstime
            originalresultstime.p <- wilcox.test(original$originalresultstime,
                                                 reduced$originalresultstime)$p.value
            originalresultstime.p <- ifelse(is.na(originalresultstime.p), 1, originalresultstime.p)
            
            originalresultstime.e <- effectsize_accurate(original$originalresultstime,
                                                         reduced$originalresultstime)$size
            
            originalresultstime.median.normal <- eval(call(func, original$originalresultstime/1000))
            originalresultstime.median.reduce <- eval(call(func, reduced$originalresultstime/1000))
            
            # reductionTime
            reductionTime.p <- wilcox.test(original$reductionTime,
                                           reduced$reductionTime)$p.value
            reductionTime.p <- ifelse(is.na(reductionTime.p), 1, reductionTime.p)
            
            reductionTime.e <- effectsize_accurate(original$reductionTime,
                                                   reduced$reductionTime)$size
            
            #reductionTime.median.normal <- eval(call(func, original$reductionTime/1000))
            reductionTime.median.reduce <- eval(call(func, reduced$reductionTime/1000))
            #browser()
            # Total Time
            #*****
            total_time.p <- NULL
            total_time.e <- NULL
            if (r == "original") {
              sticcer_total_time <- original$reductionTime + original$mutationanalysistime
              total_time.p <- wilcox.test(reduced$mutationanalysistime,
                                          sticcer_total_time)$p.value
              total_time.p <- ifelse(is.na(reductionTime.p), 1, reductionTime.p)
              
              total_time.e <- effectsize_accurate(reduced$mutationanalysistime,
                                                  sticcer_total_time)$size
            }
            
            total_time.normal <- reductionTime.median.normal + mutationanalysistime.median.normal
            total_time.reduced <- reductionTime.median.reduce + mutationanalysistime.median.reduce
            
            # Mutation
            mutationScore.p <- wilcox.test(original$mutationScore,
                                           reduced$mutationScore)$p.value
            mutationScore.p <- ifelse(is.na(mutationScore.p), 1, mutationScore.p)
            
            mutationScore.e <- effectsize_accurate(original$mutationScore,
                                                   reduced$mutationScore)$size
            
            mutationScore.median.normal <- eval(call(func, original$mutationScore))
            mutationScore.median.reduce <- eval(call(func, reduced$mutationScore))
            
            # Reduction time
            #reductionTime <- (original %>% arrange(randomseed))$testgenerationtime - (reduced %>% arrange(randomseed))$testgenerationtime
            
            new_case <- case
            
            x <- c(db, new_case, gen, r, "mutationanalysistime",
                   sigBetweenTwoValues_improved_timing("R", mutationanalysistime.median.normal, 
                                                       mutationanalysistime.median.reduce, mutationanalysistime.p,
                                                       mutationanalysistime.e,
                                                       flipEval = TRUE),
                   mutationanalysistime.p, mutationanalysistime.e,
                   ifelse(mutationanalysistime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "originalresultstime",
                   sigBetweenTwoValues_improved_timing("R", originalresultstime.median.normal, 
                                                       originalresultstime.median.reduce, originalresultstime.p,
                                                       originalresultstime.e,
                                                       flipEval = TRUE),
                   originalresultstime.p, originalresultstime.e,
                   ifelse(originalresultstime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "totalTime",
                   format(round(total_time.reduced, 2), nsmall = 2),
                   "","","")
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "reductionTime",
                   sigBetweenTwoValues_improved_timing("R", reductionTime.median.normal,
                                                       reductionTime.median.reduce, reductionTime.p,
                                                       reductionTime.e,
                                                       flipEval = TRUE),
                   reductionTime.p, reductionTime.e,
                   ifelse(reductionTime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))

            x <- c(db, new_case, gen, r, "mutationScore",
                   sigBetweenTwoValues_improved("R", mutationScore.median.normal, 
                                                mutationScore.median.reduce, mutationScore.p,
                                                mutationScore.e),
                   mutationScore.p, mutationScore.e,
                   ifelse(mutationScore.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
          }
        }
      }
    }
  }
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducedusing", "metric", 
                         "reduce", "pvalue", "effectsize", "sig")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}



greedy_vs_hgs_reduction_time <- function(df, func = "mean") {
  
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduedwith = character(),
    metric = character(),
    reduce = character(),
    pvalue = character(), 
    effectsize = character(),
    sig = character()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducedusing <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        original <- df %>% filter(casestudy == case, 
                                  dbms == db, 
                                  datagenerator == gen, 
                                  reducedwith == "greedy",
                                  fullreduce == TRUE)
        new_case <- case
        reductionTime.median.normal <- eval(call(func, original$reductionTime/1000))
        x <- c(db, new_case, gen, "sticcer", "reductionTime",
               format(round(reductionTime.median.normal, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        for (r in reducedusing) {
          if (r == "HGS") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r)
            # reductionTime
            reductionTime.p <- wilcox.test(original$reductionTime,
                                           reduced$reductionTime)$p.value
            reductionTime.p <- ifelse(is.na(reductionTime.p), 1, reductionTime.p)
            
            reductionTime.e <- effectsize_accurate(original$reductionTime,
                                                   reduced$reductionTime)$size
            
            #reductionTime.median.normal <- eval(call(func, original$reductionTime/1000))
            reductionTime.median.reduce <- eval(call(func, reduced$reductionTime/1000))
            x <- c(db, new_case, gen, r, "reductionTime",
                   sigBetweenTwoValues_improved_timing("R", reductionTime.median.normal,
                                                       reductionTime.median.reduce, reductionTime.p,
                                                       reductionTime.e,
                                                       flipEval = TRUE),
                   reductionTime.p, reductionTime.e,
                   ifelse(reductionTime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
          }
        }
      }
    }
  }
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducedusing", "metric", 
                         "reduce", "pvalue", "effectsize", "sig")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}


# Add * for sig and ~ small, ~~ meduim, and ~~~ large effect sizes
sigBetweenTwoValues_total_times <- function(side, normal, reduced, pvalue, effectsize,
                                                flipEval = FALSE) {
  value <- NULL
  
  ef = ifelse(effectsize == "large", "$^{\\ast}$", "")
  
  
  if(normal >= reduced & pvalue <= 0.05) {
    value <- paste(ef,"$\\APLdown$\\textbf{", format(round(reduced, 2), nsmall = 2), "}", sep="")
  } else if(normal <= reduced & pvalue <= 0.05) {
    value <- paste(ef,"$\\APLup$\\textbf{", format(round(reduced, 2), nsmall = 2), "}", sep="")
  } else {
    value <- as.character(format(round(reduced, 2), nsmall = 2))
  }
  
  return(value)
}

ots_vs_sticcer_total_time <- function(df, func = "mean") {
  
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduedwith = character(),
    metric = character(),
    reduce = character(),
    pvalue = character(), 
    effectsize = character(),
    sig = character()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducedusing <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        sticcer <- df %>% filter(casestudy == case, 
                                  dbms == db, 
                                  datagenerator == gen, 
                                  reducedwith == "sticcer",
                                  fullreduce == TRUE)
        mutationanalysistime.median.sticcer <- eval(call(func, sticcer$mutationanalysistime/1000))
        new_case <- case
        x <- c(db, new_case, gen, "sticcer", "mutationanalysistime",
               format(round(mutationanalysistime.median.sticcer, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        reductionTime.median.sticcer <- eval(call(func, sticcer$reductionTime/1000))
        x <- c(db, new_case, gen, "sticcer", "reductionTime",
               format(round(reductionTime.median.sticcer, 2), nsmall = 2),
               "","","")
        results <- rbind(results, as.data.frame(t(x)))
        
        # Total Time
        total_time_sticcer <- reductionTime.median.sticcer + mutationanalysistime.median.sticcer
        # x <- c(db, new_case, gen, "sticcer", "totalTime",
        #        format(round(total_time_sticcer, 2), nsmall = 2),
        #        "","","")
        # results <- rbind(results, as.data.frame(t(x)))

        for (r in reducedusing) {
          if (r == "original") {
            original <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r)

            mutationanalysistime.median.original <- eval(call(func, original$mutationanalysistime/1000))
            #browser()
            # Total Time
            #*****

            sticcer_total_time <- sticcer$reductionTime + sticcer$mutationanalysistime
            
            total_time.p <- wilcox.test(original$mutationanalysistime,
                                        sticcer_total_time)$p.value
            total_time.p <- ifelse(is.na(total_time.p), 1, total_time.p)
            
            total_time.e <- effectsize_accurate(original$mutationanalysistime,
                                                sticcer_total_time)$size
            
            #total_time_sticcer <- reductionTime.median.normal + mutationanalysistime.median.normal
            #total_time.original <- reductionTime.median.reduce + mutationanalysistime.median.reduce
            
          
            new_case <- case
            
            # Total Time
            total_time_sticcer <- reductionTime.median.sticcer + mutationanalysistime.median.sticcer
            x <- c(db, new_case, gen, "sticcer", "totalTime",
                   sigBetweenTwoValues_total_times("R", mutationanalysistime.median.original,
                                                       total_time_sticcer,
                                                       total_time.p,
                                                       total_time.e,
                                                       flipEval = TRUE),
                   total_time.p, total_time.e,
                   ifelse(total_time.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "totalTime",
                   format(round(mutationanalysistime.median.original, 2), nsmall = 2),
                   "","","")
            
            
            results <- rbind(results, as.data.frame(t(x)))

          }
        }
      }
    }
  }
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducedusing", "metric", 
                         "reduce", "pvalue", "effectsize", "sig")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}

