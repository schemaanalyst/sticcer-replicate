library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
require(tidyr)
require(purrr)
source(file = "R/effectsize.R")

read_mutationanalysis_files <- function(data_path) {
  files_one <- dir(path = data_path, pattern = "*mutationanalysis.dat")
  
  data_one <- files_one %>%
    map(function(x) read_csv(file.path(data_path, x))) %>%
    reduce(rbind)
  
  fulldata <- data_one %>% dplyr::mutate(casestudy = as.character(gsub("parsedcasestudy.","",casestudy)))
  
  return(fulldata)
}

read_mutanttiming_files <- function(data_path) {
  reductions <- c('HGS', 'random', 'sticcer', 'original', 'additionalGreedy')
  fulldata <- NULL
  for (r in reductions) {
    if (r != "original") {
      files_one <- NULL
      pattern_one <- paste(r,  ".*-mutanttiming.dat", sep="")
      files_one <- dir(path = data_path, pattern = pattern_one)
      data_one <- files_one %>%
        map(function(x) read_csv(file.path(data_path, x))) %>%
        reduce(rbind)
      data_one <- data_one %>% mutate(reduction = r)
      fulldata <- rbind(fulldata, data_one)
    } else {
      pattern_one <- paste("*-testingSchemas-mutanttiming.dat", sep="")
      files_one <- dir(path = data_path, pattern = pattern_one)
      if (is_empty(files_one) == FALSE) {
        data_one <- files_one %>%
          map(function(x) read_csv(file.path(data_path, x))) %>%
          reduce(rbind)
        
        data_one <- data_one %>% mutate(reduction = r)
        fulldata <- rbind(fulldata, data_one)
      }
      
      pattern_two <- paste("*-allminusitrust-mutanttiming.dat", sep="")
      files_two <- dir(path = data_path, pattern = pattern_two)
      if (is_empty(files_two) == FALSE) {
        data_two <- files_two %>%
          map(function(x) read_csv(file.path(data_path, x))) %>%
          reduce(rbind)
        
        data_two <- data_two %>% mutate(reduction = r)
        fulldata <- rbind(fulldata, data_two)
      }
      
      pattern_three <- paste("*-iTrust-mutanttiming.dat", sep="")
      files_three <- dir(path = data_path, pattern = pattern_three)
      if (is_empty(files_two) == FALSE) {
        data_three <- files_three %>%
          map(function(x) read_csv(file.path(data_path, x))) %>%
          reduce(rbind)
        
        data_three <- data_three %>% mutate(reduction = r)
        fulldata <- rbind(fulldata, data_three)
      }
    }
  }
  fulldata <- fulldata %>% dplyr::mutate(schema = as.character(gsub("parsedcasestudy.","",schema)))
  
  # Changes the column name from schema to casestudy
  colnames(fulldata)[which(names(fulldata) == "schema")] <- "casestudy"
  colnames(fulldata)[which(names(fulldata) == "generator")] <- "datagenerator"
  colnames(fulldata)[which(names(fulldata) == "reduction")] <- "reducedwith"
  
  return(fulldata)
}

sortSchema <- function(df) {
  #browser()
  df <- df %>% mutate_if(is.factor, as.character)
  df <- df %>% mutate(schema = if_else(schema == "Iso3166", "Isoiii", schema))
  df <- df %>% mutate(schema = if_else(schema == "IsoFlav_R2", "IsoFlav", schema))
  df <- df %>% mutate(schema = if_else(schema == "NistDML181", "NistDMLi", schema))
  df <- df %>% mutate(schema = if_else(schema == "NistDML182", "NistDMLii", schema))
  df <- df %>% mutate(schema = if_else(schema == "NistDML183", "NistDMLiii", schema))
  df <- df %>% mutate(schema = if_else(schema == "NistXTS748", "NistXTSEight", schema))
  df <- df %>% mutate(schema = if_else(schema == "NistXTS749", "NistXTSNine", schema))
  df <- df %>% arrange(schema)
  df <- df %>% mutate(schema = paste("\\", schema, "ForTable", sep = ""))
  
  return(df)
}

### Reduction Effectivness 
reductionEffectiveness <- function(df, func = "mean") {
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduction = character(),
    tests = character(),
    inserts = character()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducewith <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        original <- df %>% filter(casestudy == case, 
                                  dbms == db, 
                                  datagenerator == gen, 
                                  fullreduce == FALSE)
        for (r in reducewith) {
          if (r != "original") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r, 
                                     fullreduce == TRUE)
            # Test Suite Reduction Effectiveness 
            tests.eff <- ((1 - eval(call(func, reduced$tests))/eval(call(func, original$tests))) * 100)
            
            # INSERTs Reduction Effectiveness 
            inserts.eff <- ((1 - eval(call(func, reduced$finalinsertscounter))/eval(call(func, original$finalinsertscounter))) * 100)
            
            new_case <- case
            if (r != "sticcer") {
              random <- df %>% filter(casestudy == case, 
                                       dbms == db, 
                                       datagenerator == gen, 
                                       reducedwith == "sticcer", 
                                       fullreduce == TRUE)
              
              tests.p <- wilcox.test(random$tests,
                                     reduced$tests)$p.value
              
              tests.p <- ifelse(is.na(tests.p), 1, tests.p)
              
              tests.e <- effectsize_accurate(random$tests,
                                             reduced$tests)$size
              #browser()
              # ef.tests = ifelse(tests.e == "large",
              #                   "$^{\\ast\\ast\\ast}$",
              #                   ifelse(tests.e == "medium",
              #                          "$^{\\ast\\ast}$",
              #                          ifelse(tests.e == "small","$^{\\ast}$", "")
              #                   )
              # )
              
              ef.tests = ifelse(tests.e == "large", "$^{\\ast}$","")
              
              if(eval(call(func, reduced$tests)) >= eval(call(func, random$tests)) & tests.p <= 0.05) {
                tests.eff <- paste(ef.tests,"$\\APLdown$\\textbf{", round(tests.eff, 0), "}", sep="")
              } else if (eval(call(func, reduced$tests)) <= eval(call(func, random$tests)) & tests.p <= 0.05)  {
                tests.eff <- paste(ef.tests,"$\\APLup$\\textbf{", round(tests.eff, 0), "}", sep="")
              } else {
                tests.eff <- round(tests.eff, 0)
              }
              
              inserts.p <- wilcox.test(random$finalinsertscounter,
                                     reduced$finalinsertscounter)$p.value
              inserts.p <- ifelse(is.na(inserts.p), 1, inserts.p)
              inserts.e <- effectsize_accurate(random$finalinsertscounter,
                                             reduced$finalinsertscounter)$size
              # ef.inserts = ifelse(inserts.e == "large",
              #                     "$^{\\ast\\ast\\ast}$",
              #                     ifelse(inserts.e == "medium",
              #                            "$^{\\ast\\ast}$",
              #                            ifelse(inserts.e == "small","$^{\\ast}$", "")
              #                     )
              # )
              ef.inserts = ifelse(inserts.e == "large", "$^{\\ast}$","")
              
              if(eval(call(func, reduced$finalinsertscounter)) >= eval(call(func, random$finalinsertscounter)) & inserts.p <= 0.05) {
                inserts.eff <- paste(ef.tests,"$\\APLdown$\\textbf{", round(inserts.eff, 0), "}", sep="")
              } else if (eval(call(func, reduced$finalinsertscounter)) <= eval(call(func, random$finalinsertscounter)) & inserts.p <= 0.05)  {
                inserts.eff <- paste(ef.tests,"$\\APLup$\\textbf{", round(inserts.eff, 0), "}", sep="")
              } else {
                inserts.eff <- round(inserts.eff, 0)
              }
              
              # tests.eff <- paste(tests.eff, 
              #                    "\\% (",round(eval(call(func, reduced$tests), 0)), "/", 
              #                    round(eval(call(func, original$tests)), 0), ")" ,sep="")
              # inserts.eff <- paste(inserts.eff, 
              #                      "\\% (",
              #                      round(eval(call(func, reduced$finalinsertscounter), 0)), "/", 
              #                      round(eval(call(func, original$finalinsertscounter)), 0), ")", sep="")
              tests.eff <- paste(tests.eff, "\\%",sep="")
              inserts.eff <- paste(inserts.eff, "\\%", sep="")
              
            } else {
              tests.eff <- paste(round(tests.eff, 0), "\\%",sep="")
              inserts.eff <- paste(round(inserts.eff, 0), "\\%", sep="")
            }

            
            #tests.eff <- round(tests.eff, 0)
            #inserts.eff <- round(inserts.eff, 0)
            
            # tests.eff <- paste(round(tests.eff, 0), "\\% (", round(mean(reduced$tests),0), 
            #                    "/", round(mean(original$tests),0), ")", sep="")
            # inserts.eff <- paste(round(inserts.eff, 0), "\\% (", round(mean(reduced$finalinsertscounter),0), 
            #                      "/", round(mean(original$finalinsertscounter),0), ")", sep="")
            
            x <- c(db, new_case, gen, r, tests.eff, inserts.eff)
            results <- rbind(results, as.data.frame(t(x)))
          }
        }
      }
    }
  }
  
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducewith", "tests.eff", "inserts.eff")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}

### Reduction Effectivness 
reductionEffectivenessBrackets <- function(df, func = "mean") {
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduction = character(),
    tests = character(),
    inserts = character()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducewith <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        original <- df %>% filter(casestudy == case, 
                                  dbms == db, 
                                  datagenerator == gen, 
                                  fullreduce == FALSE)
        for (r in reducewith) {
          if (r != "original") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r, 
                                     fullreduce == TRUE)
            # Test Suite Reduction Effectiveness 
            tests.eff <- ((1 - eval(call(func, reduced$tests))/eval(call(func, original$tests))) * 100)
            
            # INSERTs Reduction Effectiveness 
            inserts.eff <- ((1 - eval(call(func, reduced$finalinsertscounter))/eval(call(func, original$finalinsertscounter))) * 100)
            
            new_case <- case
            if (r != "sticcer") {
              random <- df %>% filter(casestudy == case, 
                                      dbms == db, 
                                      datagenerator == gen, 
                                      reducedwith == "sticcer", 
                                      fullreduce == TRUE)
              tests.eff <- paste("(",round(eval(call(func, reduced$tests), 0)), "/", 
                                 round(eval(call(func, original$tests)), 0), ")" ,sep="")
              inserts.eff <- paste("(",
                                   round(eval(call(func, reduced$finalinsertscounter), 0)), "/", 
                                   round(eval(call(func, original$finalinsertscounter)), 0), ")", sep="")
              
            } else {
              tests.eff <- paste("(",round(eval(call(func, reduced$tests), 0)), "/", 
                                 round(eval(call(func, original$tests)), 0), ")" ,sep="")
              inserts.eff <- paste("(",round(eval(call(func, reduced$finalinsertscounter), 0)), "/", 
                                   round(eval(call(func, original$finalinsertscounter)), 0), ")", sep="")
            }

            x <- c(db, new_case, gen, r, tests.eff, inserts.eff)
            results <- rbind(results, as.data.frame(t(x)))
          }
        }
      }
    }
  }
  
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducewith", "tests.eff", "inserts.eff")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}


### Reduction Effectivness 
reductionEffectivenessRaw <- function(df, func = "mean") {
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduction = character(),
    tests = integer(),
    inserts = integer()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducewith <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        original <- df %>% filter(casestudy == case, 
                                  dbms == db, 
                                  datagenerator == gen, 
                                  fullreduce == FALSE)
        for (r in reducewith) {
          if (r != "original") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r, 
                                     fullreduce == TRUE)
            # Test Suite Reduction Effectiveness 
            tests.eff <- ((1 - eval(call(func, reduced$tests))/eval(call(func, original$tests))) * 100)
            
            # INSERTs Reduction Effectiveness 
            inserts.eff <- ((1 - eval(call(func, reduced$finalinsertscounter))/eval(call(func, original$finalinsertscounter))) * 100)

            # for latex purposes
            new_case <- case
            #new_case <- case
            tests.eff <- round(tests.eff, 0)
            inserts.eff <- round(inserts.eff, 0)
            
            #tests.eff <- round(tests.eff, 0)
            #inserts.eff <- round(inserts.eff, 0)
            
            # tests.eff <- paste(round(tests.eff, 0), "\\% (", round(mean(reduced$tests),0), 
            #                    "/", round(mean(original$tests),0), ")", sep="")
            # inserts.eff <- paste(round(inserts.eff, 0), "\\% (", round(mean(reduced$finalinsertscounter),0), 
            #                      "/", round(mean(original$finalinsertscounter),0), ")", sep="")
            
            x <- c(db, new_case, gen, r, tests.eff, inserts.eff)
            results <- rbind(results, as.data.frame(t(x)))
          }
        }
      }
    }
  }
  
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducewith", "tests.eff", "inserts.eff")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}

doAcount <- function(df) {
  a <- df %>% group_by(dbms, datagenerator, fullreduce, casestudy, reducedwith) %>% tally() %>% filter(n < 30)
  return(a)
}

# Add * for sig and ~ small, ~~ meduim, and ~~~ large effect sizes
sigBetweenTwoValues <- function(side, normal, reduced, pvalue, effectsize,
                                flipEval = FALSE) {
  value <- NULL
  
  ef = ifelse(effectsize == "large",
              "$^{\\ast\\ast\\ast}$",
              ifelse(effectsize == "medium",
                     "$^{\\ast\\ast}$",
                     ifelse(effectsize == "small","$^{\\ast}$", "")
              )
  )
  
  if(side == "L" & flipEval == FALSE) {
    if(normal >= reduced & pvalue <= 0.05) {
      value <- paste(ef,"\\textbf{", round(normal, 0), "}", sep="")
    } else {
      value <- as.character(round(normal, 0))
    }
  }
  
  if(side == "L" & flipEval == TRUE) {
    if(normal <= reduced & pvalue <= 0.05) {
      value <- paste(ef,"\\textbf{", round(normal, 0), "}", sep="")
    } else {
      value <- as.character(round(normal, 0))
    }
  }
  
  if(side == "R" & flipEval == FALSE) {
    if(normal <= reduced & pvalue <= 0.05) {
      value <- paste(ef,"\\textbf{", round(reduced, 0), "}", sep="")
    } else {
      value <- as.character(round(reduced, 0))
    }
  }
  
  if(side == "R" & flipEval == TRUE) {
    if(normal >= reduced & pvalue <= 0.05) {
      value <- paste(ef,"\\textbf{", round(reduced, 0), "}", sep="")
    } else {
      value <- as.character(round(reduced, 0))
    }
  }
  
  return(value)
}

# Add * for sig and ~ small, ~~ meduim, and ~~~ large effect sizes
sigBetweenTwoValuesOnePoint <- function(side, normal, reduced, pvalue, effectsize,
                                        flipEval = FALSE) {
  value <- NULL
  
  ef = ifelse(effectsize == "large",
              "$^{\\ast\\ast\\ast}$",
              ifelse(effectsize == "medium",
                     "$^{\\ast\\ast}$",
                     ifelse(effectsize == "small","$^{\\ast}$", "")
              )
  )
  
  if(side == "L" & flipEval == FALSE) {
    if(normal >= reduced & pvalue <= 0.05) {
      value <- paste(ef,"\\textbf{", round(normal, 1), "}", sep="")
    } else {
      value <- as.character(round(normal, 1))
    }
  }
  
  if(side == "L" & flipEval == TRUE) {
    if(normal <= reduced & pvalue <= 0.05) {
      value <- paste(ef,"\\textbf{", round(normal, 1), "}", sep="")
    } else {
      value <- as.character(round(normal, 1))
    }
  }
  
  if(side == "R" & flipEval == FALSE) {
    if(normal <= reduced & pvalue <= 0.05) {
      value <- paste(ef,"\\textbf{", round(reduced, 1), "}", sep="")
    } else {
      value <- as.character(round(reduced, 1))
    }
  }
  
  if(side == "R" & flipEval == TRUE) {
    if(normal >= reduced & pvalue <= 0.05) {
      value <- paste(ef,"\\textbf{", round(reduced, 1), "}", sep="")
    } else {
      value <- as.character(round(reduced, 1))
    }
  }
  
  return(value)
}


full_analysis <- function(df, func = "mean") {
  
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduedwith = character(),
    metric = character(),
    normal = character(),
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
                                  fullreduce == FALSE)
        for (r in reducedusing) {
          if (r != "original") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r,
                                     fullreduce == TRUE)
            #browser()
            if (length(original$coverage) < 30 || length(reduced$numCoveredReqReduced) < 30) {
              browser()
            }
            # Coverage
            reduced_coverage <- (reduced %>% 
                                   select(numCoveredReqReduced, numCoveredRequirments, coverage) %>% 
                                   mutate(finalcoverage = numCoveredReqReduced/numCoveredRequirments * 100) %>%
                                   mutate(finalcoverage = ifelse(finalcoverage == 100 & finalcoverage != coverage, coverage, finalcoverage)))$finalcoverage
            coverage.p <- wilcox.test(original$coverage,
                                      reduced_coverage)$p.value
            coverage.p <- ifelse(is.na(coverage.p), 1, coverage.p)
            
            coverage.e <- effectsize_accurate(original$coverage,
                                              reduced_coverage)$size
            
            coverage.median.normal <- eval(call(func, original$coverage))
            coverage.median.reduce <- eval(call(func, reduced_coverage))
            
            # Timetaking
            timetaken.p <- wilcox.test(original$timetaken,
                                       reduced$timetaken)$p.value
            timetaken.p <- ifelse(is.na(timetaken.p), 1, timetaken.p)
            
            timetaken.e <- effectsize_accurate(original$timetaken,
                                               reduced$timetaken)$size
            
            timetaken.median.normal <- eval(call(func, original$timetaken/1000))
            timetaken.median.reduce <- eval(call(func, reduced$timetaken/1000))
            
            # testgenerationtime
            testgenerationtime.p <- wilcox.test(original$testgenerationtime,
                                                reduced$testgenerationtime)$p.value
            testgenerationtime.p <- ifelse(is.na(testgenerationtime.p), 1, testgenerationtime.p)
            
            testgenerationtime.e <- effectsize_accurate(original$testgenerationtime,
                                                        reduced$testgenerationtime)$size
            
            testgenerationtime.median.normal <- eval(call(func, original$testgenerationtime/1000))
            testgenerationtime.median.reduce <- eval(call(func, reduced$testgenerationtime/1000))
            
            # mutationanalysistime
            mutationanalysistime.p <- wilcox.test(original$mutationanalysistime,
                                                  reduced$mutationanalysistime)$p.value
            mutationanalysistime.p <- ifelse(is.na(mutationanalysistime.p), 1, mutationanalysistime.p)
            
            mutationanalysistime.e <- effectsize_accurate(original$mutationanalysistime,
                                                          reduced$mutationanalysistime)$size
            
            mutationanalysistime.median.normal <- eval(call(func, original$mutationanalysistime/1000))
            mutationanalysistime.median.reduce <- eval(call(func, reduced$mutationanalysistime/1000))
            
            # Mutation
            mutationScore.p <- wilcox.test(original$mutationScore,
                                           reduced$mutationScore)$p.value
            mutationScore.p <- ifelse(is.na(mutationScore.p), 1, mutationScore.p)
            
            mutationScore.e <- effectsize_accurate(original$mutationScore,
                                                   reduced$mutationScore)$size
            
            mutationScore.median.normal <- eval(call(func, original$mutationScore))
            mutationScore.median.reduce <- eval(call(func, reduced$mutationScore))
            
            # Test Suite Size
            tests.p <- wilcox.test(original$tests,
                                   reduced$tests)$p.value
            tests.p <- ifelse(is.na(tests.p), 1, tests.p)
            
            tests.e <- effectsize_accurate(original$tests,
                                           reduced$tests)$size
            
            tests.median.normal <- eval(call(func, original$tests))
            tests.median.reduce <- eval(call(func, reduced$tests))
            
            # Number of INSERTs
            finalinsertscounter.p <- wilcox.test(original$finalinsertscounter,
                                                 reduced$finalinsertscounter)$p.value
            finalinsertscounter.p <- ifelse(is.na(finalinsertscounter.p), 1, finalinsertscounter.p)
            
            finalinsertscounter.e <- effectsize_accurate(original$finalinsertscounter,
                                                         reduced$finalinsertscounter)$size
            
            finalinsertscounter.median.normal <- eval(call(func, original$finalinsertscounter))
            finalinsertscounter.median.reduce <- eval(call(func, reduced$finalinsertscounter))
            
            new_case <- case
            #new_case <- case
            x <- c(db, new_case, gen, r, "coverage", 
                   sigBetweenTwoValues("L", coverage.median.normal, 
                                       coverage.median.reduce, coverage.p,
                                       coverage.e,
                                       flipEval = FALSE),
                   sigBetweenTwoValues("R", coverage.median.normal, 
                                       coverage.median.reduce, coverage.p,
                                       coverage.e,
                                       flipEval = FALSE), 
                   coverage.p, coverage.e,
                   ifelse(coverage.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "timetaken", 
                   sigBetweenTwoValues("L", timetaken.median.normal, 
                                       timetaken.median.reduce, timetaken.p,
                                       timetaken.e,
                                       flipEval = TRUE),
                   sigBetweenTwoValues("R", timetaken.median.normal, 
                                       timetaken.median.reduce, timetaken.p,
                                       timetaken.e,
                                       flipEval = TRUE),
                   timetaken.p, timetaken.e,
                   ifelse(timetaken.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "testgenerationtime", 
                   sigBetweenTwoValuesOnePoint("L", testgenerationtime.median.normal, 
                                               testgenerationtime.median.reduce, testgenerationtime.p,
                                               testgenerationtime.e,
                                               flipEval = TRUE),
                   sigBetweenTwoValuesOnePoint("R", testgenerationtime.median.normal, 
                                               testgenerationtime.median.reduce, testgenerationtime.p,
                                               testgenerationtime.e,
                                               flipEval = TRUE),
                   testgenerationtime.p, testgenerationtime.e,
                   ifelse(testgenerationtime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "mutationanalysistime", 
                   sigBetweenTwoValuesOnePoint("L", mutationanalysistime.median.normal, 
                                               mutationanalysistime.median.reduce, mutationanalysistime.p,
                                               mutationanalysistime.e,
                                               flipEval = TRUE),
                   sigBetweenTwoValuesOnePoint("R", mutationanalysistime.median.normal, 
                                               mutationanalysistime.median.reduce, mutationanalysistime.p,
                                               mutationanalysistime.e,
                                               flipEval = TRUE),
                   mutationanalysistime.p, mutationanalysistime.e,
                   ifelse(mutationanalysistime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "mutationScore", 
                   sigBetweenTwoValues("L", mutationScore.median.normal, 
                                       mutationScore.median.reduce, mutationScore.p,
                                       mutationScore.e),
                   sigBetweenTwoValues("R", mutationScore.median.normal, 
                                       mutationScore.median.reduce, mutationScore.p,
                                       mutationScore.e),
                   mutationScore.p, mutationScore.e,
                   ifelse(mutationScore.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "tests", 
                   sigBetweenTwoValues("L", tests.median.normal, 
                                       tests.median.reduce, tests.p,
                                       tests.e,
                                       flipEval = TRUE),
                   sigBetweenTwoValues("R", tests.median.normal, 
                                       tests.median.reduce, tests.p,
                                       tests.e,
                                       flipEval = TRUE),
                   tests.p, tests.e,
                   ifelse(tests.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "finalinsertscounter", 
                   sigBetweenTwoValues("L", finalinsertscounter.median.normal, 
                                       finalinsertscounter.median.reduce, 
                                       finalinsertscounter.p,
                                       finalinsertscounter.e,
                                       flipEval = TRUE),
                   sigBetweenTwoValues("R", finalinsertscounter.median.normal, 
                                       finalinsertscounter.median.reduce, 
                                       finalinsertscounter.p,
                                       finalinsertscounter.e,
                                       flipEval = TRUE),
                   finalinsertscounter.p, finalinsertscounter.e,
                   ifelse(finalinsertscounter.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
          }
        }
      }
    }
  }
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducedusing", "metric", 
                         "normal", "reduce", "pvalue", "effectsize",
                         "sig")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}


sticcer_side_analysis <- function(df, func = "mean") {
  
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduedwith = character(),
    metric = character(),
    normal = character(),
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
        for (r in reducedusing) {
          if (r != "sticcer") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r)
            if (length(original$coverage) < 30 || length(reduced$numCoveredReqReduced) < 30) {
              browser()
            }
            # Coverage
            reduced_coverage <- (reduced %>% 
                                   select(numCoveredReqReduced, numCoveredRequirments, coverage) %>% 
                                   mutate(finalcoverage = numCoveredReqReduced/numCoveredRequirments * 100) %>%
                                   mutate(finalcoverage = ifelse(finalcoverage == 100 & finalcoverage != coverage, coverage, finalcoverage)))$finalcoverage
            coverage.p <- wilcox.test(original$coverage,
                                      reduced_coverage)$p.value
            coverage.p <- ifelse(is.na(coverage.p), 1, coverage.p)
            
            coverage.e <- effectsize_accurate(original$coverage,
                                              reduced_coverage)$size
            
            coverage.median.normal <- eval(call(func, original$coverage))
            coverage.median.reduce <- eval(call(func, reduced_coverage))
            
            # Timetaking
            timetaken.p <- wilcox.test(original$timetaken,
                                       reduced$timetaken)$p.value
            timetaken.p <- ifelse(is.na(timetaken.p), 1, timetaken.p)
            
            timetaken.e <- effectsize_accurate(original$timetaken,
                                               reduced$timetaken)$size
            
            timetaken.median.normal <- eval(call(func, original$timetaken/1000))
            timetaken.median.reduce <- eval(call(func, reduced$timetaken/1000))
            
            # testgenerationtime
            testgenerationtime.p <- wilcox.test(original$testgenerationtime,
                                                reduced$testgenerationtime)$p.value
            testgenerationtime.p <- ifelse(is.na(testgenerationtime.p), 1, testgenerationtime.p)
            
            testgenerationtime.e <- effectsize_accurate(original$testgenerationtime,
                                                        reduced$testgenerationtime)$size
            
            testgenerationtime.median.normal <- eval(call(func, original$testgenerationtime/1000))
            testgenerationtime.median.reduce <- eval(call(func, reduced$testgenerationtime/1000))
            
            # mutationanalysistime
            mutationanalysistime.p <- wilcox.test(original$mutationanalysistime,
                                                  reduced$mutationanalysistime)$p.value
            mutationanalysistime.p <- ifelse(is.na(mutationanalysistime.p), 1, mutationanalysistime.p)
            
            mutationanalysistime.e <- effectsize_accurate(original$mutationanalysistime,
                                                          reduced$mutationanalysistime)$size
            
            mutationanalysistime.median.normal <- eval(call(func, original$mutationanalysistime/1000))
            mutationanalysistime.median.reduce <- eval(call(func, reduced$mutationanalysistime/1000))
            
            # Mutation
            mutationScore.p <- wilcox.test(original$mutationScore,
                                           reduced$mutationScore)$p.value
            mutationScore.p <- ifelse(is.na(mutationScore.p), 1, mutationScore.p)
            
            mutationScore.e <- effectsize_accurate(original$mutationScore,
                                                   reduced$mutationScore)$size
            
            mutationScore.median.normal <- eval(call(func, original$mutationScore))
            mutationScore.median.reduce <- eval(call(func, reduced$mutationScore))
            
            # Test Suite Size
            tests.p <- wilcox.test(original$tests,
                                   reduced$tests)$p.value
            tests.p <- ifelse(is.na(tests.p), 1, tests.p)
            
            tests.e <- effectsize_accurate(original$tests,
                                           reduced$tests)$size
            
            tests.median.normal <- eval(call(func, original$tests))
            tests.median.reduce <- eval(call(func, reduced$tests))
            
            # Number of INSERTs
            finalinsertscounter.p <- wilcox.test(original$finalinsertscounter,
                                                 reduced$finalinsertscounter)$p.value
            finalinsertscounter.p <- ifelse(is.na(finalinsertscounter.p), 1, finalinsertscounter.p)
            
            finalinsertscounter.e <- effectsize_accurate(original$finalinsertscounter,
                                                         reduced$finalinsertscounter)$size
            
            finalinsertscounter.median.normal <- eval(call(func, original$finalinsertscounter))
            finalinsertscounter.median.reduce <- eval(call(func, reduced$finalinsertscounter))
            
            new_case <- case
            #new_case <- case
            x <- c(db, new_case, gen, r, "coverage", 
                   sigBetweenTwoValues("L", coverage.median.normal, 
                                       coverage.median.reduce, coverage.p,
                                       coverage.e,
                                       flipEval = FALSE),
                   sigBetweenTwoValues("R", coverage.median.normal, 
                                       coverage.median.reduce, coverage.p,
                                       coverage.e,
                                       flipEval = FALSE), 
                   coverage.p, coverage.e,
                   ifelse(coverage.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "timetaken", 
                   sigBetweenTwoValues("L", timetaken.median.normal, 
                                       timetaken.median.reduce, timetaken.p,
                                       timetaken.e,
                                       flipEval = TRUE),
                   sigBetweenTwoValues("R", timetaken.median.normal, 
                                       timetaken.median.reduce, timetaken.p,
                                       timetaken.e,
                                       flipEval = TRUE),
                   timetaken.p, timetaken.e,
                   ifelse(timetaken.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "testgenerationtime", 
                   sigBetweenTwoValuesOnePoint("L", testgenerationtime.median.normal, 
                                               testgenerationtime.median.reduce, testgenerationtime.p,
                                               testgenerationtime.e,
                                               flipEval = TRUE),
                   sigBetweenTwoValuesOnePoint("R", testgenerationtime.median.normal, 
                                               testgenerationtime.median.reduce, testgenerationtime.p,
                                               testgenerationtime.e,
                                               flipEval = TRUE),
                   testgenerationtime.p, testgenerationtime.e,
                   ifelse(testgenerationtime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "mutationanalysistime", 
                   sigBetweenTwoValuesOnePoint("L", mutationanalysistime.median.normal, 
                                               mutationanalysistime.median.reduce, mutationanalysistime.p,
                                               mutationanalysistime.e,
                                               flipEval = TRUE),
                   sigBetweenTwoValuesOnePoint("R", mutationanalysistime.median.normal, 
                                               mutationanalysistime.median.reduce, mutationanalysistime.p,
                                               mutationanalysistime.e,
                                               flipEval = TRUE),
                   mutationanalysistime.p, mutationanalysistime.e,
                   ifelse(mutationanalysistime.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "mutationScore",
                   sigBetweenTwoValues("L", mutationScore.median.normal, 
                                       mutationScore.median.reduce, mutationScore.p,
                                       mutationScore.e),
                   sigBetweenTwoValues("R", mutationScore.median.normal, 
                                       mutationScore.median.reduce, mutationScore.p,
                                       mutationScore.e),
                   mutationScore.p, mutationScore.e,
                   ifelse(mutationScore.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "tests", 
                   sigBetweenTwoValues("L", tests.median.normal, 
                                       tests.median.reduce, tests.p,
                                       tests.e,
                                       flipEval = TRUE),
                   sigBetweenTwoValues("R", tests.median.normal, 
                                       tests.median.reduce, tests.p,
                                       tests.e,
                                       flipEval = TRUE),
                   tests.p, tests.e,
                   ifelse(tests.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
            
            x <- c(db, new_case, gen, r, "finalinsertscounter", 
                   sigBetweenTwoValues("L", finalinsertscounter.median.normal, 
                                       finalinsertscounter.median.reduce, 
                                       finalinsertscounter.p,
                                       finalinsertscounter.e,
                                       flipEval = TRUE),
                   sigBetweenTwoValues("R", finalinsertscounter.median.normal, 
                                       finalinsertscounter.median.reduce, 
                                       finalinsertscounter.p,
                                       finalinsertscounter.e,
                                       flipEval = TRUE),
                   finalinsertscounter.p, finalinsertscounter.e,
                   ifelse(finalinsertscounter.p <= 0.05, "*", ""))
            results <- rbind(results, as.data.frame(t(x)))
          }
        }
      }
    }
  }
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducedusing", "metric", 
                         "mediansticcer", "reduce", "pvalue", "effectsize",
                         "sig")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}


reductionEffectivenessOTS <- function(df, func = "mean") {
  results <- data.frame(
    dbms = character(),
    schema = character(), 
    datagenerator = character(),
    reduction = character(),
    tests = character(),
    inserts = character()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducewith <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        original <- df %>% filter(casestudy == case, 
                                  dbms == db, 
                                  datagenerator == gen, 
                                  fullreduce == FALSE)
        for (r in reducewith) {
          if (r != "original") {
            reduced <- df %>% filter(casestudy == case, 
                                     dbms == db, 
                                     datagenerator == gen, 
                                     reducedwith == r, 
                                     fullreduce == TRUE)
            # Test Suite Reduction Effectiveness 
            tests.eff <- ((1 - eval(call(func, reduced$tests))/eval(call(func, original$tests))) * 100)
            
            # INSERTs Reduction Effectiveness 
            inserts.eff <- ((1 - eval(call(func, reduced$finalinsertscounter))/eval(call(func, original$finalinsertscounter))) * 100)
            
            new_case <- case
            if (r != "original") {
              random <- original
              
              tests.p <- wilcox.test(random$tests,
                                     reduced$tests)$p.value
              tests.e <- effectsize_accurate(random$tests,
                                             reduced$tests)$size
              
              ef.tests = ifelse(tests.e == "large",
                                "$^{\\ast\\ast\\ast}$",
                                ifelse(tests.e == "medium",
                                       "$^{\\ast\\ast}$",
                                       ifelse(tests.e == "small","$^{\\ast}$", "")
                                )
              )
              
              if(eval(call(func, original$tests)) >= eval(call(func, random$tests)) & tests.p <= 0.05) {
                tests.eff <- paste(ef.tests,"\\textbf{", round(tests.eff, 0), "}", sep="")
              } else if (eval(call(func, original$tests)) <= eval(call(func, random$tests)) & tests.p <= 0.05)  {
                tests.eff <- paste(ef.tests,"\\textit{", round(tests.eff, 0), "}", sep="")
              } else {
                tests.eff <- round(tests.eff, 0)
              }
              
              inserts.p <- wilcox.test(random$finalinsertscounter,
                                       reduced$finalinsertscounter)$p.value
              inserts.e <- effectsize_accurate(random$finalinsertscounter,
                                               reduced$finalinsertscounter)$size
              ef.inserts = ifelse(inserts.e == "large",
                                  "$^{\\ast\\ast\\ast}$",
                                  ifelse(inserts.e == "medium",
                                         "$^{\\ast\\ast}$",
                                         ifelse(inserts.e == "small","$^{\\ast}$", "")
                                  )
              )
              
              if(eval(call(func, original$finalinsertscounter)) >= eval(call(func, random$finalinsertscounter)) & inserts.p <= 0.05) {
                inserts.eff <- paste(ef.tests,"\\textbf{", round(inserts.eff, 0), "}", sep="")
              } else if (eval(call(func, original$finalinsertscounter)) <= eval(call(func, random$finalinsertscounter)) & inserts.p <= 0.05)  {
                inserts.eff <- paste(ef.tests,"\\textit{", round(inserts.eff, 0), "}", sep="")
              } else {
                inserts.eff <- round(inserts.eff, 0)
              }
              
              tests.eff <- paste(tests.eff, 
                                 "\\% (",round(eval(call(func, reduced$tests), 0)), "/", 
                                 round(eval(call(func, original$tests)), 0), ")" ,sep="")
              inserts.eff <- paste(inserts.eff, 
                                   "\\% (",
                                   round(eval(call(func, reduced$finalinsertscounter), 0)), "/", 
                                   round(eval(call(func, original$finalinsertscounter)), 0), ")", sep="")
              
            } else {
              tests.eff <- paste(round(tests.eff, 0), 
                                 "\\% (",round(eval(call(func, reduced$tests), 0)), "/", 
                                 round(eval(call(func, original$tests)), 0), ")" ,sep="")
              inserts.eff <- paste(round(inserts.eff, 0), "\\% (",
                                   round(eval(call(func, reduced$finalinsertscounter), 0)), "/", 
                                   round(eval(call(func, original$finalinsertscounter)), 0), ")", sep="")
            }
            
            
            #tests.eff <- round(tests.eff, 0)
            #inserts.eff <- round(inserts.eff, 0)
            
            # tests.eff <- paste(round(tests.eff, 0), "\\% (", round(mean(reduced$tests),0), 
            #                    "/", round(mean(original$tests),0), ")", sep="")
            # inserts.eff <- paste(round(inserts.eff, 0), "\\% (", round(mean(reduced$finalinsertscounter),0), 
            #                      "/", round(mean(original$finalinsertscounter),0), ")", sep="")
            
            x <- c(db, new_case, gen, r, tests.eff, inserts.eff)
            results <- rbind(results, as.data.frame(t(x)))
          }
        }
      }
    }
  }
  
  colnames(results) <- c("dbms", "schema", "datagenerator", "reducewith", "tests.eff", "inserts.eff")
  
  return(results)
  #return(print(xtable::xtable(results), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
}


# Reduction Time
reduction_time_analysis <- function(df, func = "mean") {
  
  results <- data.frame(
    dbms = character(),
    casestudy = character(), 
    datagenerator = character(),
    reduedwith = character(),
    randomseed = character(),
    reductionTime = character()
  )
  
  dbs <- as.vector(dplyr::distinct(df, dbms))[[1]]
  generators <- as.vector(dplyr::distinct(df, datagenerator))[[1]]
  cases <- as.vector(dplyr::distinct(df, casestudy))[[1]]
  reducedusing <- as.vector(dplyr::distinct(df, reducedwith))[[1]]
  randomseeds <- as.vector(dplyr::distinct(df, randomseed))[[1]]
  #browser()
  for(db in dbs) {
    for(gen in generators) {
      for(case in cases) {
        for(seed in randomseeds) {
          original <- df %>% filter(casestudy == case, 
                                    dbms == db, 
                                    datagenerator == gen, 
                                    randomseed == seed,
                                    #reducedwith == "original",
                                    fullreduce == FALSE)
          for (r in reducedusing) {
            if (r == "original") {
              reductionTimeCalc <- 0
              new_case <- case
              x <- c(db, new_case, gen, r, seed, reductionTimeCalc)
              results <- rbind(results, as.data.frame(t(x)))
            } else {
              reduced <- df %>% filter(casestudy == case, 
                                       dbms == db, 
                                       datagenerator == gen,
                                       randomseed == seed,
                                       reducedwith == r)
              
              # Reduction time
              reductionTimeCalc <- (reduced %>% arrange(randomseed))$testgenerationtime - (original %>% arrange(randomseed))$testgenerationtime
              new_case <- case
              x <- c(db, new_case, gen, r, seed, reductionTimeCalc)
              results <- rbind(results, as.data.frame(t(x)))
            }
          }
        }
      }
    }
  }
  colnames(results) <- c("dbms", "casestudy", "datagenerator", "reducedwith", "randomseed", "reductionTime")
  
  results$reducedwith <- as.character(results$reducedwith)
  results$dbms <- as.character(results$dbms)
  results$casestudy <- as.character(results$casestudy)
  results$datagenerator <- as.character(results$datagenerator)
  results$randomseed <- as.numeric(results$randomseed)
  #results$reductionTime <- as.numeric(results$reductionTime)
  results$reductionTime <- as.numeric(as.character(results$reductionTime))
  
  
  return(results)
}