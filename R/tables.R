source(file = "R/main.R")

library(knitr)
library(kableExtra)

################## HEADERS ###########################################
reductionHeaders <- c("Schema",
                  "\\randomR", "\\Greedy", "\\HGS", "\\sticcer",
                  "\\randomR", "\\Greedy", "\\HGS", "\\sticcer")

reductionHeadersLevelTwo <- c("Schema",
                      "\\\\sticcer" = 2, "\\\\randomR" = 2, "\\\\Greedy" = 2, "\\\\HGS" = 2,
                      "\\\\sticcer" = 2, "\\\\randomR" = 2, "\\\\Greedy" = 2, "\\\\HGS" = 2)

reductionTabs <- c(" ",
                  "\\\\sticcer" = 2, "\\\\Greedy" = 2,
                  "\\\\HGS" = 2, "\\\\randomR" = 2,
                  "\\\\sticcer" = 2, "\\\\Greedy" = 2,
                  "\\\\HGS" = 2, "\\\\randomR" = 2)

generatorsTab <- c(" ", "\\\\AVMD" = 8, "\\\\domino" = 8)

originalReducedTab <- c("Schema",
                        "OTS", "RTS",
                        "OTS", "RTS",
                        "OTS", "RTS",
                        "OTS", "RTS",
                        "OTS", "RTS",
                        "OTS", "RTS",
                        "OTS", "RTS",
                        "OTS", "RTS")

################ Means ###############################
###  reduction effectiveness with min, max, avg
reductioneff <- reductionEffectivenessRaw(mutationanalysis, func = "median")
reductioneff$reducewith <- as.character(reductioneff$reducewith)
reductioneff <- reductioneff %>% mutate(reducewith = ifelse(reducewith == "simpleGreedy" | reducewith == "additionalGreedy", "greedy", reducewith))

reductionTestsTable <- reductioneff %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducewith")) %>%
  filter(variable == "tests.eff") %>%
  dcast(schema ~ dbms + datagenerator + reducewith + variable, value.var = "value") %>%
  sortSchema()

reductionInsertsTable <- reductioneff %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducewith")) %>%
  filter(variable == "inserts.eff") %>%
  dcast(schema ~ dbms + datagenerator + reducewith + variable, value.var = "value") %>%
  sortSchema()

################# DO MIN, MEAN, MAX #########################################
library(psych)
reductionTableDescribe <- describe(reductionTestsTable) %>% select(min, mean, max)
reductionTableDescribe <- as.data.frame(reductionTableDescribe)
reductionTableDescribe <- cbind(rownames(reductionTableDescribe), reductionTableDescribe)
rownames(reductionTableDescribe) <- NULL
colnames(reductionTableDescribe) <- c("var", "min", "mean", "max")
reductionTableDescribe <- reductionTableDescribe %>% dplyr::mutate(var = as.character(gsub("\\*","",var))) %>% filter(var != "schema")
reductionTableDescribe <- reductionTableDescribe %>% melt(id.vars = c("var")) %>% dcast(variable ~ var)
reductionTableDescribe[,-1] <- round(reductionTableDescribe[,-1],0)

reductionTableDescribeINS <- describe(reductionInsertsTable) %>% select(min, mean, max)
reductionTableDescribeINS <- as.data.frame(reductionTableDescribeINS)
reductionTableDescribeINS <- cbind(rownames(reductionTableDescribeINS), reductionTableDescribeINS)
rownames(reductionTableDescribeINS) <- NULL
colnames(reductionTableDescribeINS) <- c("var", "min", "mean", "max")
reductionTableDescribeINS <- reductionTableDescribeINS %>% dplyr::mutate(var = as.character(gsub("\\*","",var))) %>% filter(var != "schema")
reductionTableDescribeINS <- reductionTableDescribeINS %>% melt(id.vars = c("var")) %>% dcast(variable ~ var)
reductionTableDescribeINS[,-1] <- round(reductionTableDescribeINS[,-1],0)

############# GET ORIGINAL Values #################################################
reductioneff <- reductionEffectiveness(mutationanalysis, func = "median")
#reductioneff <- reductionEffectivenessOTS(mutationanalysis, func = "median")
reductioneff$reducewith <- as.character(reductioneff$reducewith)
reductioneff <- reductioneff %>% mutate(reducewith = ifelse(reducewith == "simpleGreedy" | reducewith == "additionalGreedy", "greedy", reducewith))

reductionTestsTable <- reductioneff %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducewith")) %>%
  filter(variable == "tests.eff") %>%
  dcast(schema ~ dbms + datagenerator + reducewith + variable, value.var = "value") %>%
  sortSchema()

colnames(reductionTableDescribe)[1] <- "schema"
reductionTableDescribe[,-1] <- lapply(reductionTableDescribe[,-1], function(x) paste(x, "\\%", sep = ""))
reductionTestsTable <- rbind(reductionTestsTable, reductionTableDescribe)
#browser()
reductionTestsTable <- reductionTestsTable %>% mutate(schema = replace(schema, schema == "min", "Minimum"))
reductionTestsTable <- reductionTestsTable %>% mutate(schema = replace(schema, schema == "max", "Maximum"))
reductionTestsTable <- reductionTestsTable %>% mutate(schema = replace(schema, schema == "mean", "Average"))

reductionEffBrackets <- reductionEffectivenessBrackets(mutationanalysis, func = "median")
redTestTable <- reductionEffBrackets %>% rename(tests.effect = tests.eff, inserts.effect = inserts.eff) %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducewith")) %>%
  filter(variable == "tests.effect") %>%
  dcast(schema ~ dbms + datagenerator + reducewith + variable, value.var = "value") %>%
  sortSchema()

reductionTestsTable <- left_join(reductionTestsTable, redTestTable) %>%
  select("schema",
         "SQLite_avsDefaults_sticcer_tests.eff", "SQLite_avsDefaults_sticcer_tests.effect",
         "SQLite_avsDefaults_random_tests.eff", "SQLite_avsDefaults_random_tests.effect",
         "SQLite_avsDefaults_greedy_tests.eff", "SQLite_avsDefaults_greedy_tests.effect",
         "SQLite_avsDefaults_HGS_tests.eff", "SQLite_avsDefaults_HGS_tests.effect",
         "SQLite_dominoRandom_sticcer_tests.eff", "SQLite_dominoRandom_sticcer_tests.effect",
         "SQLite_dominoRandom_random_tests.eff","SQLite_dominoRandom_random_tests.effect",
         "SQLite_dominoRandom_greedy_tests.eff", "SQLite_dominoRandom_greedy_tests.effect",
         "SQLite_dominoRandom_HGS_tests.eff", "SQLite_dominoRandom_HGS_tests.effect")

#reductionTestsTable <- reductionTestsTable[c(1,5,3,4,2,9,7,8,6)]

kable(reductionTestsTable, format = "latex", booktabs = T, label = "testReduction",
      align=rep('r', length(reductionTestsTable[,1])),
      caption = "Median Reduction Effectiveness of Tests. Higher is better.",
      col.names = NA,
      escape = FALSE, linesep = "") %>%
  row_spec(0, align = "c") %>%
  kable_styling(latex_options = c("striped", "scale_down"), full_width = F) %>%
  add_header_above(reductionHeadersLevelTwo, escape = FALSE)  %>%
  add_header_above(c(" ", "\\\\AVMD" = 8, "\\\\domino" = 8), escape = FALSE)  %>%
  cat(., file = "texTables/testReductionTableWithSummary-median.tex")

rm(reductionTestsTable)
rm(reductionTableDescribe)

##### INSERTS #####################################
reductionInsertsTable <- reductioneff %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducewith")) %>%
  filter(variable == "inserts.eff") %>%
  dcast(schema ~ dbms + datagenerator + reducewith + variable, value.var = "value") %>%
  sortSchema()

colnames(reductionTableDescribeINS)[1] <- "schema"
reductionTableDescribeINS[,-1] <- lapply(reductionTableDescribeINS[,-1], function(x) paste(x, "\\%", sep = ""))
reductionInsertsTable <- rbind(reductionInsertsTable, reductionTableDescribeINS)

reductionInsertsTable <- reductionInsertsTable %>% mutate(schema = replace(schema, schema == "min", "Minimum"))
reductionInsertsTable <- reductionInsertsTable %>% mutate(schema = replace(schema, schema == "max", "Maximum"))
reductionInsertsTable <- reductionInsertsTable %>% mutate(schema = replace(schema, schema == "mean", "Average"))

#reductionInsertsTable <- reductionInsertsTable[c(1,5,3,4,2,9,7,8,6)]

reductionEffBrackets <- reductionEffectivenessBrackets(mutationanalysis, func = "median")
redTestTable <- reductionEffBrackets %>% rename(tests.effect = tests.eff, inserts.effect = inserts.eff) %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducewith")) %>%
  filter(variable == "inserts.effect") %>%
  dcast(schema ~ dbms + datagenerator + reducewith + variable, value.var = "value") %>%
  sortSchema()

reductionInsertsTable <- left_join(reductionInsertsTable, redTestTable) %>%
  select("schema",
         "SQLite_avsDefaults_sticcer_inserts.eff", "SQLite_avsDefaults_sticcer_inserts.effect",
         "SQLite_avsDefaults_random_inserts.eff", "SQLite_avsDefaults_random_inserts.effect",
         "SQLite_avsDefaults_greedy_inserts.eff", "SQLite_avsDefaults_greedy_inserts.effect",
         "SQLite_avsDefaults_HGS_inserts.eff", "SQLite_avsDefaults_HGS_inserts.effect",
         "SQLite_dominoRandom_sticcer_inserts.eff", "SQLite_dominoRandom_sticcer_inserts.effect",
         "SQLite_dominoRandom_random_inserts.eff","SQLite_dominoRandom_random_inserts.effect",
         "SQLite_dominoRandom_greedy_inserts.eff", "SQLite_dominoRandom_greedy_inserts.effect",
         "SQLite_dominoRandom_HGS_inserts.eff", "SQLite_dominoRandom_HGS_inserts.effect")

kable(reductionInsertsTable, format = "latex", booktabs = T, label = "insertReduction",
      align=rep('r', length(reductionInsertsTable[,1])),
      caption = "Median Reduction Effectiveness of INSERTs. Higher is better.",
      col.names = NA,
      escape = FALSE, linesep = "") %>%
  row_spec(0, align = "c") %>%
  kable_styling(latex_options = c("striped", "scale_down"), full_width = F) %>%
  add_header_above(reductionHeadersLevelTwo, escape = FALSE)  %>%
  add_header_above(c(" ", "\\\\AVMD" = 8, "\\\\domino" = 8), escape = FALSE)  %>%
  cat(., file = "texTables/insertsReductionTableWithSummary-median.tex")

rm(reductionInsertsTable)
rm(reductionTableDescribeINS)
rm(reductioneff)

# ########### Do other analysis comparing with STICCER
source(file = "R/sticcer_analysis.R")
comb_resultsAnalysis <- ots_side_analysis_improved(mutationanalysis, func = "median")

mutationsticcerTable <- comb_resultsAnalysis %>%
  select(dbms, schema, reducedusing, metric, datagenerator, reduce) %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducedusing", "metric")) %>%
  filter(metric == "mutationScore") %>%
  dcast(schema ~ dbms + datagenerator + reducedusing + variable, value.var = "value") %>%
  sortSchema()
#browser()
mutationsticcerTable <- mutationsticcerTable[c(1,2,5,6,3,4,7,10,11,8,9)]
#browser()
kable(mutationsticcerTable %>% filter_all(any_vars(str_detect(., '\\$'))), format = "latex", booktabs = T, label = "mutationScore",
      align=rep('r', length(mutationsticcerTable[,1])),
      caption = "Medain Mutation Analysis. A bold number denote significance, and $\\ast$, $\\ast\\ast$, $\\ast\\ast\\ast$ denotes small, medium, and large effect sizes respectively. The $+$ denote that \\sticcer is worse and $-$ denote \\sticcer is better.",
      col.names = c("Schemas",
                    "OTS", "\\sticcer", "\\randomR", "\\Greedy", "\\HGS",
                    "OTS", "\\sticcer", "\\randomR", "\\Greedy", "\\HGS"),
      escape = FALSE, linesep = "") %>%
  row_spec(0, align = "c") %>%
  kable_styling(latex_options = c("striped", "scale_down"), full_width = F) %>%
  add_header_above(c(" ", "\\\\AVMD" = 5, "\\\\domino" = 5), escape = FALSE) %>%
  cat(., file = "texTables/sticcerMutationTable-median.tex")

rm(mutationsticcerTable)

comb_resultsAnalysis <- sticcer_side_analysis_improved(mutationanalysis, func = "median")

mutationanalysistimesticcerTable <- comb_resultsAnalysis %>%
  select(dbms, schema, reducedusing, metric, datagenerator, reduce) %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducedusing", "metric")) %>%
  filter(metric == "mutationanalysistime") %>%
  dcast(schema ~ dbms + datagenerator + reducedusing + variable, value.var = "value") %>%
  sortSchema()

# mutationanalysistimesticcerTable <- mutationanalysistimesticcerTable[c(1,2,5,6,3,4,7,10,11,8,9)]

mutationanalysistimesticcerTable <- mutationanalysistimesticcerTable[c(1,2,6,3,4,7,11,8,9)]

kable(mutationanalysistimesticcerTable, format = "latex", booktabs = T, label = "sticcerMutationTiming",
      align=rep('r', length(mutationanalysistimesticcerTable[,1])),
      caption = "Median Mutation Analysis Timing In Seconds. A bold number denote significance, and $\\ast$, $\\ast\\ast$, $\\ast\\ast\\ast$ denotes small, medium, and large effect sizes respectively. The $+$ denote that \\sticcer is worse and $-$ denote \\sticcer is better.",
      col.names = c("Schemas",
                    "\\sticcer", "\\randomR", "\\Greedy", "\\HGS",
                    "\\sticcer", "\\randomR", "\\Greedy", "\\HGS"),
      escape = FALSE, linesep = "") %>%
  row_spec(0, align = "c") %>%
  kable_styling(latex_options = c("striped", "scale_down"), full_width = F) %>%
  add_header_above(c(" ", "\\\\AVMD" = 5, "\\\\domino" = 5), escape = FALSE) %>%
  cat(., file = "texTables/sticcerMutationTimingTable-medain.tex")

rm(mutationanalysistimesticcerTable)

originalresultstimesticcerTable <- comb_resultsAnalysis %>%
  select(dbms, schema, reducedusing, metric, datagenerator, reduce) %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducedusing", "metric")) %>%
  filter(metric == "originalresultstime") %>%
  dcast(schema ~ dbms + datagenerator + reducedusing + variable, value.var = "value") %>%
  sortSchema()

originalresultstimesticcerTable <- originalresultstimesticcerTable[c(1,5,6,3,4,2,10,11,8,9,7)]

kable(originalresultstimesticcerTable, format = "latex", booktabs = T, label = "sticcerOriginalResultsTime",
      align=rep('r', length(originalresultstimesticcerTable[,1])),
      caption = "Median Test Suite Execution Timing In Seconds. A bold number denote significance, and $\\ast$, $\\ast\\ast$, $\\ast\\ast\\ast$ denotes small, medium, and large effect sizes respectively. The $+$ denote that \\sticcer is worse and $-$ denote \\sticcer is better.",
      col.names = c("Schemas",
                    "OTS", "\\randomR", "\\Greedy", "\\HGS", "\\sticcer",
                    "OTS", "\\randomR", "\\Greedy", "\\HGS", "\\sticcer"),
      escape = FALSE, linesep = "") %>%
  row_spec(0, align = "c") %>%
  kable_styling(latex_options = c("striped", "scale_down"), full_width = F) %>%
  add_header_above(c(" ", "\\\\AVMD" = 5, "\\\\domino" = 5), escape = FALSE) %>%
  cat(., file = "texTables/originalresultstimesticcerTable-medain.tex")

rm(originalresultstimesticcerTable)

sticcerReductionTimeTable <- comb_resultsAnalysis %>%
  select(dbms, schema, reducedusing, metric, datagenerator, reduce) %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducedusing", "metric")) %>%
  filter(metric == "reductionTime") %>%
  dcast(schema ~ dbms + datagenerator + reducedusing + variable, value.var = "value") %>%
  sortSchema()

sticcerReductionTimeTable <- sticcerReductionTimeTable[c(1,5,6,3,4,2,10,11,8,9,7)]

kable(sticcerReductionTimeTable, format = "latex", booktabs = T, label = "reductionTimeTable",
      align=rep('r', length(sticcerReductionTimeTable[,1])),
      caption = "Median Reduction Timing In Seconds. A bold number denote significance, and $\\ast$, $\\ast\\ast$, $\\ast\\ast\\ast$ denotes small, medium, and large effect sizes respectively. The $+$ denote that \\sticcer is worse and $-$ denote \\sticcer is better.",
      col.names = c("Schemas",
                    "OTS", "\\randomR", "\\Greedy", "\\HGS", "\\sticcer",
                    "OTS", "\\randomR", "\\Greedy", "\\HGS", "\\sticcer"),
      escape = FALSE, linesep = "") %>%
  row_spec(0, align = "c") %>%
  kable_styling(latex_options = c("striped", "scale_down"), full_width = F) %>%
  add_header_above(c(" ", "\\\\AVMD" = 5, "\\\\domino" = 5), escape = FALSE) %>%
  cat(., file = "texTables/sticcerReductionTimeTable-medain.tex")

rm(sticcerReductionTimeTable)
#rm(comb_resultsAnalysis)
totaltimes <- ots_vs_sticcer_total_time(mutationanalysis, func = "median")
#browser()

fullsticceranalysis <- totaltimes %>% 
  select(dbms, schema, reducedusing, metric, datagenerator, reduce) %>%
  melt(id.vars = c("dbms", "schema", "datagenerator", "reducedusing", "metric")) %>%
  dcast(schema ~ dbms + datagenerator + reducedusing + metric, value.var = "value") %>%
  sortSchema() %>%
  select(schema,
         SQLite_avsDefaults_original_totalTime,
         SQLite_avsDefaults_sticcer_reductionTime, SQLite_avsDefaults_sticcer_mutationanalysistime, SQLite_avsDefaults_sticcer_totalTime,
         SQLite_dominoRandom_original_totalTime,
         SQLite_dominoRandom_sticcer_reductionTime, SQLite_dominoRandom_sticcer_mutationanalysistime, SQLite_dominoRandom_sticcer_totalTime)

kable(fullsticceranalysis, format = "latex", booktabs = T, label = "fullTimingTable",
      align=rep('r', length(fullsticceranalysis[,1])),
      caption = "Median of Reduction Timing, Mutation Analysis Timing,  and Total time In Seconds. A bold number denote significance, and $\\ast$, $\\ast\\ast$, $\\ast\\ast\\ast$ denotes small, medium, and large effect sizes respectively. The $+$ denote that \\sticcer is worse and $-$ denote \\sticcer is better.",
      col.names = c("Schemas",
                    "Total",
                    "RT", "MT", "Total",
                    "Total",
                    "RT", "MT", "Total"),
      escape = FALSE, linesep = "") %>%
  row_spec(0, align = "c") %>%
  kable_styling(latex_options = c("striped", "scale_down"), full_width = F) %>%
  add_header_above(c(" ",
                     "OTS" = 1,
                     "\\\\sticcer" = 3,
                     "OTS" = 1,
                     "\\\\sticcer" = 3), escape = FALSE) %>%
  add_header_above(c(" ", "\\\\AVMD" = 4, "\\\\domino" = 4), escape = FALSE) %>%
  cat(., file = "texTables/OTSvsSticcerTotalTimingTable-medain.tex")

#fullsticceranalysis <- fullsticceranalysis[c(1,2,5,6,3,4,7,10,11,8,9)]

# Merges
mergesFile %>%
  group_by(datagenerator, schema) %>%
  summarise(medianMergedTests = median(merges)) %>%
  melt(id.vars = c("schema", "datagenerator"))  %>%
  mutate(value = round(value, 0)) %>%
  dcast(schema ~ datagenerator + variable, value.var = "value") %>%
  sortSchema() %>%
  kable(format = "latex", booktabs = T, label = "mergesTable",
        caption = "Median Number of Merges",
        col.names = c("Schemas",
                      "\\AVMD No. Merges", "\\domino No. Merges"),
        escape = FALSE, linesep = "") %>%
  kable_styling(latex_options = c("striped", "scale_down"), full_width = F)  %>%
  cat(., file = "texTables/merges-medain.tex")

mergesPlot <- mergesFile %>%
  group_by(datagenerator, schema) %>%
  summarise(medianMergedTests = median(merges)) %>%
  arrange(schema) %>%
  ggplot(aes(x = factor(schema, levels = rev(levels(factor(schema)))),
             y=medianMergedTests, fill=datagenerator)) +
  geom_bar(stat="identity", position=position_dodge2(reverse = TRUE)) +
  coord_flip() +
  scale_y_continuous(trans='log2') +
  xlab(element_blank()) + ylab(bquote('Number of Merges '(log[2]))) +
  scale_fill_manual(name="Test Data\nGenerators",
                    values=c('black','lightgray'),
                    breaks=c("avsDefaults", "dominoRandom"),
                    labels=c("AVM-D", "DOMINO")) + theme_bw()

mergesPlot %>% ggsave(device = "pdf", filename="plots/MergesBarChar.pdf", 
                      width = 5, height = 6, units = "in")
