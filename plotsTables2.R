

rm(list = ls())

# This script takes all the final model fits (one for each population) and
# generates and saves every plot and table used.



# libraries
library(purrr)
library(kableExtra)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)


# This section produces the summary tables in the paper and the set of 
# population level plots


# processed data sets: chronics and episodics
chronics = readRDS("chronics.Rds")
episodics = readRDS("episodics.Rds")
# transformed data used for modeling
chronics_transformed = readRDS("final results and code/chronics_transformed_data.Rds")
episodics_transformed = readRDS("final results and code/episodics_transformed_data.Rds")

# draws
chron.fit = readRDS("./final results and code/cloglog_full_chronics.Rds")
chron.draws = as.data.frame(chron.fit, pars = c("thetaf", "theta", "theta_random", "lambda"))

epis.fit = readRDS("./final results and code/cloglog_full_episodics.Rds")
epis.draws = as.data.frame(epis.fit, pars = c("thetaf", "theta", "theta_random", "lambda"))

# Note: these functions extract draws for summaries, and do not convert anything
convert_scale.epis = function(data, draws) {
  lambda.draws = draws[,grepl("lamb", colnames(draws))]
  lambda.draws = apply(lambda.draws, 2, function(nu) {
    1 - exp(-exp(nu))
  })
  thetaf.draws = draws[,grepl("f", colnames(draws))]
  theta.draws = draws[,grepl("theta\\[", colnames(draws))]
  theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
  patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
  patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                 paste0("theta_random", pattern, "2]"),
                                                                                                 paste0("theta_random", pattern, "3]"),
                                                                                                 paste0("theta_random", pattern, "4]"),
                                                                                                 paste0("theta_random", pattern, "5]")))) }
  theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
  theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}),
                                        function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
  return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
}

convert_scale.chron = function(data, draws) {
  lambda.draws = draws[,grepl("lamb", colnames(draws))]
  lambda.draws = apply(lambda.draws, 2, function(nu) {
    1 - exp(-exp(nu))
  })
  thetaf.draws = draws[,grepl("f", colnames(draws))]
  theta.draws = draws[,grepl("theta\\[", colnames(draws))]
  theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
  patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
  patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                 paste0("theta_random", pattern, "2]")))) }
  theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
  theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}),
                                        function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
  return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
  
}

# posterior summaries function
post.summary = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), sd = sd(x), 
                                  quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))
}

# exrracted draws in the needed format for summaries:
converted.draws.chron = convert_scale.chron(data = chronics, draws = chron.draws)
converted.draws.epis = convert_scale.epis(data = episodics, draws = epis.draws)

chronics.table.baseline = post.summary(converted.draws.chron$lambda.draws)
chronics.table.baseline.df = as.data.frame(chronics.table.baseline)
chronics.table.baseline.df$days = as.numeric(levels(unique(chronics$interval)))
chronics.table.baseline.df$id = "chronics"

episodics.table.baseline = post.summary(converted.draws.epis$lambda.draws)
episodics.table.baseline.df = as.data.frame(episodics.table.baseline)
episodics.table.baseline.df$days = as.numeric(levels(unique(episodics$interval)))
episodics.table.baseline.df$id = "episodics"

total.table.baseline = rbind(chronics.table.baseline.df, episodics.table.baseline.df)
total.table.baseline$id = as.factor(total.table.baseline$id)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
baseffects = ggplot(data = total.table.baseline, aes(x = days, y = mean, color = id))  +
  geom_line() + geom_point() + scale_colour_manual(values=cbbPalette) +
  geom_line(aes(x = days, y = `2.5%`, color = id), linetype = "dashed") + 
  geom_line(aes(x = days, y = `97.5%`, color = id), linetype = "dashed") + theme_bw() + 
  ylab("Probability of migraine attack") + xlab("Days between attacks") + guides(color=guide_legend(title="Group")) + 
                                               theme(plot.title = element_text(hjust = 0.5, size = 16),
                                                     axis.text.x = element_text(size = 12),
                                                     axis.text.y = element_text(size = 12),  
                                                     axis.title.x = element_text(size = 14),
                                                     axis.title.y = element_text(size = 14),
                                                     legend.text=element_text(size=12),
                                                     legend.title=element_text(size=14))

chronics.table.pop = rbind(post.summary(converted.draws.chron$thetaf),
                           post.summary(converted.draws.chron$thetar))
rownames(chronics.table.pop) = c("humidity", "pressure", "temperature", "wind_gust", "wind_speed",
                                 "humidity^2", "pressure^2", "temperature^2", "wind_gust^2", "wind_speed^2")

rname = rownames(chronics.table.pop)
chronics.table.pop = as.data.frame(chronics.table.pop)
chronics.table.pop$weather = rname

episodics.table.pop = rbind(post.summary(converted.draws.epis$thetaf),
                            post.summary(converted.draws.epis$thetar))
rownames(episodics.table.pop) = c("humidity", "pressure", "temperature", "wind_gust", "wind_speed",
                                  "humidity^2", "pressure^2", "temperature^2", "wind_gust^2", "wind_speed^2")

rname = rownames(episodics.table.pop)
episodics.table.pop = as.data.frame(episodics.table.pop)
episodics.table.pop$weather = rname

geffects = ggplot(data = chronics.table.pop) + geom_segment(aes(x=!!sym('2.5%'), y=weather, xend=!!sym('97.5%'), yend=weather), color = "#000000",
                                                            size = 1.25, alpha = 0.7) + 
  geom_segment(data = episodics.table.pop, aes(x=!!sym('2.5%'), y=weather, xend=!!sym('97.5%'), yend=weather), color = "#E69F00",
               size = 1.25, alpha = 0.6) + 
  theme_bw() + xlab("Effect size") + ylab("Weather factor") + 
                                        theme(plot.title = element_text(hjust = 0.5, size = 16),
                                              axis.text.x = element_text(size = 14),
                                              axis.text.y = element_text(size = 14, angle = 0),
                                              axis.title.x = element_text(size = 14),
                                              axis.title.y = element_text(size = 14),
                                              legend.text=element_text(size=12),
                                              legend.title=element_text(size=14),
  ) + labs(color = "Population") +
  scale_color_hue(labels = c("Chronics", "Episodics"))  + scale_y_discrete(labels = c(expression({humidity}), expression({humidity}^2), expression({pressure}), expression({pressure}^2),
                                                                                      expression({temperature}), expression({temperature^2}),expression({'wind gust'}), expression({'wind gust'}^2),
                                                                                      expression({'wind speed'}), expression({'wind speed'}^2)) )


png("./final results and code/generated_tables_plots/pop_effects_comparison.png", units="in", width=10, height=5.5, res=600)

plots_arranged = ggarrange(
  plotlist = list(baseffects, geffects),
  nrow = 1, ncol = 2,
  labels = LETTERS[1:2],
  legend = "bottom",
  common.legend = TRUE
) 
plots_arranged

dev.off()


# probability that effect of weather factor X is greater than population (increased risk 
# relative to population)
get_prop2 = function(draws_ind, draws_pop) {
  lapply(as.data.frame(as.matrix(draws_ind) - as.matrix(draws_pop)), 
         function(x) {sum(x > 0)/nrow(draws_ind)})
}

chron.p.select = as.data.frame(do.call(rbind, lapply(
  map(converted.draws.chron$thetar.random, get_prop2, draws_pop = converted.draws.chron$thetar),
  unlist))[c(6,7,27,37),])
chron.p.select$individual = as.character(c(6,7,27,37))
colnames(chron.p.select) = c("Pressure", "Wind gust", "Individual")
epis.p.select = as.data.frame(do.call(rbind, lapply(
  map(converted.draws.epis$thetar.random, get_prop2, draws_pop = converted.draws.epis$thetar),
  unlist))[c(20,22,133,146),])
epis.p.select$individual = as.character(c(20,22,133,146))
colnames(epis.p.select) = c("Humidity", "Pressure", "Temperature", "wind gust", "Wind speed", "Individual")

saveRDS(chron.p.select, "./final results and code/generated_tables_plots/chron.p.select.Rds")
saveRDS(epis.p.select, "./final results and code/generated_tables_plots/epis.p.select.Rds")

saveRDS(chronics.table.baseline, "./final results and code/generated_tables_plots/chronics.table.baseline.Rds")
saveRDS(episodics.table.baseline, "./final results and code/generated_tables_plots/episodics.table.baseline.Rds")

saveRDS(chronics.table.pop, "./final results and code/generated_tables_plots/chronics.table.pop.Rds")
saveRDS(episodics.table.pop, "./final results and code/generated_tables_plots/episodics.table.pop.Rds")


chron.p = as.data.frame(do.call(rbind, lapply(
  map(converted.draws.chron$thetar.random, get_prop2, draws_pop = converted.draws.chron$thetar),
  unlist)))
chron.p$individual = as.character(1:38)
colnames(chron.p) = c("Pressure", "Wind gust", "Individual")

epis.p = as.data.frame(do.call(rbind, lapply(
  map(converted.draws.epis$thetar.random, get_prop2, draws_pop = converted.draws.epis$thetar),
  unlist)))
epis.p$Individual = as.character(1:191)
colnames(epis.p) = c("Humidity", "Pressure", "Temperature", "wind gust", "Wind speed", "Individual")

saveRDS(chron.p, "./final results and code/generated_tables_plots/chron.p.Rds")
saveRDS(epis.p, "./final results and code/generated_tables_plots/epis.p.Rds")


rm(list = ls())


# This section generates the individual level plots



# data sets: chronics and episodics
chronics = readRDS("chronics.Rds")
episodics = readRDS("episodics.Rds")
chronics_transformed = readRDS("final results and code/chronics_transformed_data.Rds")
episodics_transformed = readRDS("final results and code/episodics_transformed_data.Rds")

# draws
fit_chron = readRDS("./final results and code/cloglog_full_chronics.Rds")
chron.draws = as.data.frame(fit_chron, pars = c("thetaf", "theta", "theta_random", "lambda"))
fit_epis = readRDS("./final results and code/cloglog_full_episodics.Rds")
epis.draws = as.data.frame(fit_epis, pars = c("thetaf", "theta", "theta_random", "lambda"))

# helper functions 
inv_cloglog = function(nu) {
  1 - exp(-exp(nu))
}
post.mean = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x))))
}

unconvert_scale.epis = function(draws, data) {
  lambda.draws = draws[,grepl("lambda", colnames(draws))]
  thetaf.draws = draws[,grepl("f", colnames(draws))]
  theta.draws = draws[,grepl("theta\\[", colnames(draws))]
  theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
  patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
  patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                 paste0("theta_random", pattern, "2]"),
                                                                                                 paste0("theta_random", pattern, "3]"),
                                                                                                 paste0("theta_random", pattern, "4]"),
                                                                                                 paste0("theta_random", pattern, "5]")))) }
  theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
  theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}), 
                                        function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
  return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
}

unconvert_scale.chron = function(draws, data) {
  lambda.draws = draws[,grepl("lambda", colnames(draws))]
  thetaf.draws = draws[,grepl("f", colnames(draws))]
  theta.draws = draws[,grepl("theta\\[", colnames(draws))]
  theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
  patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
  patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                 paste0("theta_random", pattern, "2]")))) }
  theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
  theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}), 
                                        function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
  return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
}


unconverted.draws.epis = unconvert_scale.epis(draws = epis.draws, data = episodics)
unconverted.draws.chron = unconvert_scale.chron(draws = chron.draws, data = chronics)

# posterior summaries function
post.summary = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), sd = sd(x), 
                                  quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))
}
post.summary.confint = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))))
}

get_draws_indiv = function(draws, patient2ID) {
  cbind(draws$thetaf, draws$thetar.random[[patient2ID]], draws$lambda.draws)
}

x = episodics_transformed
colnames(x) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
x = as.data.frame(x)
x = cbind(x,episodics$patient2,episodics$interval)
colnames(x)[11:12] = c("patient2", "interval")

predict_episod = x %>% group_by(patient2, interval) %>% 
  summarize_at(vars(contains("VarSt")), mean)



predict_indiv = function(patient2ID) {
  data_indiv = predict_episod %>% filter(patient2 == patient2ID)
  draws_indiv = get_draws_indiv(draws = unconverted.draws.epis, patient2ID = patient2ID)
  
  lamb_vals = 1:length(data_indiv$interval)
  lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
  
  pred_indiv = 
    do.call(rbind,
            map(lamb_vals, 
                function(i) {
                  cbind(
                    as.data.frame(map2( draws_indiv %>% select(contains("theta")),  
                                        (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,],
                                        `*`)),
                    lam_indiv[,i]
                  ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() %>% post.summary.confint()
                }))
  
  pred_indiv = as.data.frame(pred_indiv)
  pred_indiv$patient2 = patient2ID
  pred_indiv$interval = as.numeric(paste(unique(predict_episod$interval)))[1:nrow(pred_indiv)]
  pred_indiv
}


predictions = do.call(rbind, map(unique(predict_episod$patient2), predict_indiv))
predictions$patient2 = as.factor(predictions$patient2)


predict_smooths = predictions %>% filter(interval < 11) %>% group_by(patient2) %>%
  filter(all(mean <= 0.35))
predict_high_risk = predictions %>% filter(interval < 11) %>% group_by(patient2) %>%
  filter(!all(mean <= 0.35))

palette.mike.jackie <- c("burlywood", "burlywood4", 
                         "lightsteelblue", "lightsteelblue4", 
                         "lightpink2", "lightpink4", 
                         "#FFDB6D", "#C4961A", 
                         "tan2", "#D16103", 
                         "darkseagreen3", "#52854C", 
                         "#4E84C4", "#293352")
p_epis = 
  ggplot(data = predict_high_risk, aes( x = interval, y = mean, color = patient2, group = patient2 )) + 
  stat_smooth(geom = "line", data = predict_smooths, aes(x = interval, y = mean, group = patient2), se=F, color = "gray", size = 0.4,
              alpha = 0.3) +
  geom_point(size = 1.5) + geom_line(size = 1) +  
  ylab("Probability of migraine attack") + xlab("Days between attacks") + theme_bw() +
  ggtitle("Episodic migraineurs") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                          axis.text.x = element_text(size = 14),
                                          axis.text.y = element_text(size = 14),  
                                          axis.title.x = element_text(size = 14),
                                          axis.title.y = element_text(size = 14),
                                          legend.text=element_text(size=12),
                                          legend.title=element_text(size=14),
  ) + labs(color = "Individual ID") + 
  scale_x_continuous(breaks = seq(0, 10, by = 1)) + scale_color_manual(values = palette.mike.jackie)



x = chronics_transformed
colnames(x) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
x = as.data.frame(x)
x = cbind(x,chronics$patient2,chronics$interval)
colnames(x)[11:12] = c("patient2", "interval")

predict_chron = x %>% group_by(patient2, interval) %>% 
  summarize_at(vars(contains("VarSt")), mean)



predict_indiv = function(patient2ID) {
  data_indiv = predict_chron %>% filter(patient2 == patient2ID)
  data_indiv = data_indiv %>% select(patient2, interval, humidityVarSt, pressureVarSt, temperatureVarSt, windGustVarSt, windSpeedVarSt, humidityVarSt_sq, temperatureVarSt_sq,
                                     windSpeedVarSt_sq, pressureVarSt_sq, windGustVarSt_sq)
  draws_indiv = get_draws_indiv(draws = unconverted.draws.chron, patient2ID = patient2ID)
  
  lamb_vals = 1:length(data_indiv$interval)
  lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
  
  pred_indiv = 
    do.call(rbind,
            map(lamb_vals, 
                function(i) {
                  cbind(
                    as.data.frame(map2( draws_indiv %>% select(contains("theta")),  
                                        (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,],
                                        `*`)),
                    lam_indiv[,i]
                  ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() %>% post.summary.confint()
                }))
  
  pred_indiv = as.data.frame(pred_indiv)
  pred_indiv$patient2 = patient2ID
  pred_indiv$interval = as.numeric(paste(unique(predict_chron$interval)))[1:nrow(pred_indiv)]
  pred_indiv
}


predictions = do.call(rbind, map(unique(predict_chron$patient2), predict_indiv))
predictions$patient2 = as.factor(predictions$patient2)


predict_smooths = predictions %>% filter(interval < 11) %>% group_by(patient2) %>%
  filter(all(mean <= 0.58))
predict_high_risk = predictions %>% filter(interval < 11) %>% group_by(patient2) %>%
  filter(!all(mean <= 0.58))

palette.mike.jackie <- c("burlywood", "burlywood4", 
                         "lightsteelblue", "lightsteelblue4", 
                         "lightpink2", "lightpink4", 
                         "#FFDB6D", "#C4961A", 
                         "tan2", "#D16103", 
                         "darkseagreen3", "#52854C", 
                         "#4E84C4", "#293352")
p_chron = 
  ggplot(data = predict_high_risk, aes( x = interval, y = mean, color = patient2, group = patient2 )) + 
  stat_smooth(geom = "line", data = predict_smooths, aes(x = interval, y = mean, group = patient2), se=F, color = "gray", size = 0.4,
              alpha = 0.3) + 
  geom_point(size = 1.5) + geom_line(size = 1) +  
  ylab("Probability of migraine attack") + xlab("Days between attacks") + theme_bw() +
  ggtitle("Chronic migraineurs") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                          axis.text.x = element_text(size = 14),
                                          axis.text.y = element_text(size = 14),  
                                          axis.title.x = element_text(size = 14),
                                          axis.title.y = element_text(size = 14),
                                          legend.text=element_text(size=12),
                                          legend.title=element_text(size=14),
  ) + labs(color = "Individual ID") + 
  scale_x_continuous(breaks = seq(0, 10, by = 1)) + scale_color_manual(values = palette.mike.jackie)




png("./final results and code/generated_tables_plots/probability_plots.png", units="in", width=10, height=14, res=600)

ggarrange(
  p_epis, 
  p_chron, 
  nrow = 2, 
  labels = c("A", "B")
) 

plots_arranged = ggarrange(
  plotlist = list(p_epis, p_chron),
  nrow = 2, ncol = 1,
  labels = LETTERS[1:2]
) 

plots_arranged

dev.off()



rm(list = ls())



# breakdown



# data sets: chronics and episodics
episodics = readRDS("episodics.Rds")
episodics_transformed = readRDS("final results and code/episodics_transformed_data.Rds")

# draws
fit_epis = readRDS("./final results and code/cloglog_full_episodics.Rds")
epis.draws = as.data.frame(fit_epis, pars = c("thetaf", "theta", "theta_random", "lambda"))

# helper functions 
inv_cloglog = function(nu) {
  1 - exp(-exp(nu))
}
post.mean = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x))))
}

unconvert_scale.epis = function(draws, data) {
  lambda.draws = draws[,grepl("lambda", colnames(draws))]
  thetaf.draws = draws[,grepl("f", colnames(draws))]
  theta.draws = draws[,grepl("theta\\[", colnames(draws))]
  theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
  patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
  patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                 paste0("theta_random", pattern, "2]"),
                                                                                                 paste0("theta_random", pattern, "3]"),
                                                                                                 paste0("theta_random", pattern, "4]"),
                                                                                                 paste0("theta_random", pattern, "5]")))) }
  theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
  theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}), 
                                        function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
  return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
}



unconverted.draws.epis = unconvert_scale.epis(draws = epis.draws, data = episodics)

# posterior summaries function
post.summary = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), sd = sd(x), 
                                  quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))
}
post.summary.confint = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))))
}

get_draws_indiv = function(draws, patient2ID) {
  cbind(draws$thetaf, draws$thetar.random[[patient2ID]], draws$lambda.draws)
}

x = episodics_transformed
colnames(x) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
x = as.data.frame(x)
x = cbind(x,episodics$patient2,episodics$interval)
colnames(x)[11:12] = c("patient2", "interval")

predict_episod = x %>% group_by(patient2, interval) %>% 
  summarize_at(vars(contains("VarSt")), mean)



predict_indiv_mi = function(patient2ID, mi) {
  mk = c(mi, mi + 5)
  data_indiv = predict_episod %>% filter(patient2 == patient2ID)
  draws_indiv = get_draws_indiv(draws = unconverted.draws.epis, patient2ID = patient2ID)
  
  lamb_vals = 1:length(data_indiv$interval)
  lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
  
  pred_indiv_mi = 
    do.call(rbind,
            
            map(lamb_vals, 
                function(i) {
                  (
                    cbind(
                      as.data.frame(map2( (draws_indiv %>% select(contains("theta"))),  
                                          (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,],
                                          `*`)),
                      lam_indiv[,i]
                    ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() - 
                      
                      cbind(
                        as.data.frame(map2( (draws_indiv %>% select(contains("theta")))[,-mk],  
                                            (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,][,-mk],
                                            `*`)),
                        lam_indiv[,i]
                      ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() 
                    
                  ) %>% post.summary.confint()
                  
                })
            
    )
  
  pred_indiv_mi = as.data.frame(pred_indiv_mi)
  pred_indiv_mi$patient2 = patient2ID
  pred_indiv_mi$interval = as.numeric(paste(unique(predict_episod$interval)))[1:nrow(pred_indiv_mi)]
  pred_indiv_mi
}

predictions_hum = do.call(rbind, 
                          c( map(unique(predict_episod$patient2), predict_indiv_mi, mi = 1)  )
)

predictions_press = do.call(rbind, 
                            c( map(unique(predict_episod$patient2), predict_indiv_mi, mi = 2)  )
)

predictions_temp = do.call(rbind, 
                           c( map(unique(predict_episod$patient2), predict_indiv_mi, mi = 3)  )
)

predictions_windgust = do.call(rbind, 
                               c( map(unique(predict_episod$patient2), predict_indiv_mi, mi = 4)  )
)

predictions_windspeed = do.call(rbind, 
                                c( map(unique(predict_episod$patient2), predict_indiv_mi, mi = 5)  )
)

breakdown = data.frame(
  dHumidity = predictions_hum$mean,
  dPressure = predictions_press$mean,
  dTemperature = predictions_temp$mean,
  dWind_gust = predictions_windgust$mean,
  dWind_speed = predictions_windspeed$mean,
  patient2 = predictions_hum$patient2,
  interval = predictions_hum$interval
)  

breakdown = breakdown %>% select(patient2, interval, starts_with("d")) %>%
  pivot_longer(cols = starts_with("d"), names_to = "weather_factor", values_to = "delta_probability")


palette.mike.jackie <- c("burlywood", 
                         "lightsteelblue4", 
                         "lightpink2", 
                         "#FFDB6D", 
                         "tan2", 
                         "#52854C")



plot_list = vector(mode = "list", length = 8)
names(plot_list) = paste0("p", seq(1,8,by=1))

ID2 = 55
plot_list$p1 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 55") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 89
plot_list$p2 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 89") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 105
plot_list$p3 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 105") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 123
plot_list$p4 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 123") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 130
plot_list$p5 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 130") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 137
plot_list$p6 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 137") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 142
plot_list$p7 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 142") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 146
plot_list$p8 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 146") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)

png("./final results and code/generated_tables_plots/breakdown_probability_plots_episodics.png", units="in", width=10, height=14, res=600)


plots_arranged = ggarrange(
  plotlist = plot_list,
  nrow = 4, ncol = 2,
  labels = LETTERS[1:8],
  legend = "bottom",
  common.legend = TRUE
) 

plots_arranged 

dev.off()






rm(list = ls())





# data sets: chronics and episodics
chronics = readRDS("chronics.Rds")
chronics_transformed = readRDS("final results and code/chronics_transformed_data.Rds")

# draws
fit_chron = readRDS("./final results and code/cloglog_full_chronics.Rds")
chron.draws = as.data.frame(fit_chron, pars = c("thetaf", "theta", "theta_random", "lambda"))

# helper functions 
inv_cloglog = function(nu) {
  1 - exp(-exp(nu))
}
post.mean = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x))))
}

unconvert_scale.chron = function(draws, data) {
  lambda.draws = draws[,grepl("lambda", colnames(draws))]
  thetaf.draws = draws[,grepl("f", colnames(draws))]
  theta.draws = draws[,grepl("theta\\[", colnames(draws))]
  theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
  patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
  patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                 paste0("theta_random", pattern, "2]")))) }
  theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
  theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}), 
                                        function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
  return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
}


unconverted.draws.chron = unconvert_scale.chron(draws = chron.draws, data = chronics)

# posterior summaries function
post.summary = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), sd = sd(x), 
                                  quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))
}
post.summary.confint = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))))
}

get_draws_indiv = function(draws, patient2ID) {
  cbind(draws$thetaf, draws$thetar.random[[patient2ID]], draws$lambda.draws)
}

x = chronics_transformed
colnames(x) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
x = as.data.frame(x)
x = cbind(x,chronics$patient2,chronics$interval)
colnames(x)[11:12] = c("patient2", "interval")

predict_chron = x %>% group_by(patient2, interval) %>% 
  summarize_at(vars(contains("VarSt")), mean)



predict_indiv_mi = function(patient2ID, mi) {
  mk = c(mi, mi + 5)
  data_indiv = predict_chron %>% filter(patient2 == patient2ID)
  draws_indiv = get_draws_indiv(draws = unconverted.draws.chron, patient2ID = patient2ID)
  
  lamb_vals = 1:length(data_indiv$interval)
  lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
  
  pred_indiv_mi = 
    do.call(rbind,
            
            map(lamb_vals, 
                function(i) {
                  (
                    cbind(
                      as.data.frame(map2( (draws_indiv %>% select(contains("theta"))),  
                                          (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,],
                                          `*`)),
                      lam_indiv[,i]
                    ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() - 
                      
                      cbind(
                        as.data.frame(map2( (draws_indiv %>% select(contains("theta")))[,-mk],  
                                            (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,][,-mk],
                                            `*`)),
                        lam_indiv[,i]
                      ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() 
                    
                  ) %>% post.summary.confint()
                  
                })
            
    )
  
  pred_indiv_mi = as.data.frame(pred_indiv_mi)
  pred_indiv_mi$patient2 = patient2ID
  pred_indiv_mi$interval = as.numeric(paste(unique(predict_chron$interval)))[1:nrow(pred_indiv_mi)]
  pred_indiv_mi
}

predictions_hum = do.call(rbind, 
                          c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 1)  )
)

predictions_press = do.call(rbind, 
                            c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 2)  )
)

predictions_temp = do.call(rbind, 
                           c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 3)  )
)

predictions_windgust = do.call(rbind, 
                               c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 4)  )
)

predictions_windspeed = do.call(rbind, 
                                c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 5)  )
)

breakdown = data.frame(
  dHumidity = predictions_hum$mean,
  dPressure = predictions_press$mean,
  dTemperature = predictions_temp$mean,
  dWind_gust = predictions_windgust$mean,
  dWind_speed = predictions_windspeed$mean,
  patient2 = predictions_hum$patient2,
  interval = predictions_hum$interval
)  

breakdown = breakdown %>% select(patient2, interval, starts_with("d")) %>%
  pivot_longer(cols = starts_with("d"), names_to = "weather_factor", values_to = "delta_probability")


palette.mike.jackie <- c("burlywood", 
                         "lightsteelblue4", 
                         "lightpink2", 
                         "#FFDB6D", 
                         "tan2", 
                         "#52854C")



plot_list = vector(mode = "list", length = 8)
names(plot_list) = paste0("p", seq(1,8,by=1))

ID2 = 1
plot_list$p1 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 1") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 3
plot_list$p2 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 3") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                   axis.text.x = element_text(size = 14),
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   legend.text=element_text(size=12),
                                   legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 6
plot_list$p3 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 6") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14),
                                    axis.title.x = element_text(size = 14),
                                    axis.title.y = element_text(size = 14),
                                    legend.text=element_text(size=12),
                                    legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 8
plot_list$p4 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 8") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14),
                                    axis.title.x = element_text(size = 14),
                                    axis.title.y = element_text(size = 14),
                                    legend.text=element_text(size=12),
                                    legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 22
plot_list$p5 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 22") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14),
                                    axis.title.x = element_text(size = 14),
                                    axis.title.y = element_text(size = 14),
                                    legend.text=element_text(size=12),
                                    legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 24
plot_list$p6 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 24") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14),
                                    axis.title.x = element_text(size = 14),
                                    axis.title.y = element_text(size = 14),
                                    legend.text=element_text(size=12),
                                    legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 26
plot_list$p7 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 26") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14),
                                    axis.title.x = element_text(size = 14),
                                    axis.title.y = element_text(size = 14),
                                    legend.text=element_text(size=12),
                                    legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)
ID2 = 34
plot_list$p8 = 
  ggplot(data = breakdown %>% filter(patient2 == ID2), aes(x = interval, y = delta_probability, color = weather_factor)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) + 
  theme_bw() + xlab("Days between migraine attacks") + ylab(expression(Delta*"Probability")) +
  ggtitle("Individual 34") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14),
                                    axis.title.x = element_text(size = 14),
                                    axis.title.y = element_text(size = 14),
                                    legend.text=element_text(size=12),
                                    legend.title=element_text(size=14),
  ) + labs(color = "Weather factor") +
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + 
  scale_color_manual(labels = c("Humidity", "Pressure", 
                                "Temperature", "Wind gust", "Wind speed"),
                     values = palette.mike.jackie)

png("./final results and code/generated_tables_plots/breakdown_probability_plots_chronics.png", units="in", width=10, height=14, res=600)


plots_arranged = ggarrange(
  plotlist = plot_list,
  nrow = 4, ncol = 2,
  labels = LETTERS[1:8],
  legend = "bottom",
  common.legend = TRUE
) 

plots_arranged 

dev.off()



rm(list = ls())




# This section selects one individual and plots their curve over time and breakdown by factor with 
# 95% posterior probability confidence intervals on chronic Individual 21

# data sets: chronics and episodics
chronics = readRDS("chronics.Rds")
chronics_transformed = readRDS("final results and code/chronics_transformed_data.Rds")

# draws
fit_chron = readRDS("./final results and code/cloglog_full_chronics.Rds")
chron.draws = as.data.frame(fit_chron, pars = c("thetaf", "theta", "theta_random", "lambda"))

# helper functions 
inv_cloglog = function(nu) {
  1 - exp(-exp(nu))
}
post.mean = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x))))
}

unconvert_scale.chron = function(draws, data) {
  lambda.draws = draws[,grepl("lambda", colnames(draws))]
  thetaf.draws = draws[,grepl("f", colnames(draws))]
  theta.draws = draws[,grepl("theta\\[", colnames(draws))]
  theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
  patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
  patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                 paste0("theta_random", pattern, "2]")))) }
  theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
  theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}), 
                                        function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
  return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
}


unconverted.draws.chron = unconvert_scale.chron(draws = chron.draws, data = chronics)

# posterior summaries function
post.summary = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), sd = sd(x), 
                                  quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))
}
post.summary.confint = function(draws) {
  t(apply(draws, 2, function(x) c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))))
}

get_draws_indiv = function(draws, patient2ID) {
  cbind(draws$thetaf, draws$thetar.random[[patient2ID]], draws$lambda.draws)
}

x = chronics_transformed
colnames(x) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
x = as.data.frame(x)
x = cbind(x,chronics$patient2,chronics$interval)
colnames(x)[11:12] = c("patient2", "interval")

predict_chron = x %>% group_by(patient2, interval) %>% 
  summarize_at(vars(contains("VarSt")), mean)

predict_indiv = function(patient2ID) {
  data_indiv = predict_chron %>% filter(patient2 == patient2ID)
  data_indiv = data_indiv %>% select(patient2, interval, humidityVarSt, pressureVarSt, temperatureVarSt, windGustVarSt, windSpeedVarSt, humidityVarSt_sq, temperatureVarSt_sq,
                                     windSpeedVarSt_sq, pressureVarSt_sq, windGustVarSt_sq)
  draws_indiv = get_draws_indiv(draws = unconverted.draws.chron, patient2ID = patient2ID)
  
  lamb_vals = 1:length(data_indiv$interval)
  lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
  
  pred_indiv = 
    do.call(rbind,
            map(lamb_vals, 
                function(i) {
                  cbind(
                    as.data.frame(map2( draws_indiv %>% select(contains("theta")),  
                                        (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,],
                                        `*`)),
                    lam_indiv[,i]
                  ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() %>% post.summary.confint()
                }))
  
  pred_indiv = as.data.frame(pred_indiv)
  pred_indiv$patient2 = patient2ID
  pred_indiv$interval = as.numeric(paste(unique(predict_chron$interval)))[1:nrow(pred_indiv)]
  pred_indiv
}


predictions = do.call(rbind, map(unique(predict_chron$patient2), predict_indiv))
predictions$patient2 = as.factor(predictions$patient2)


predict_indiv_mi = function(patient2ID, mi) {
  mk = c(mi, mi + 5)
  data_indiv = predict_chron %>% filter(patient2 == patient2ID)
  draws_indiv = get_draws_indiv(draws = unconverted.draws.chron, patient2ID = patient2ID)
  
  lamb_vals = 1:length(data_indiv$interval)
  lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
  
  pred_indiv_mi = 
    do.call(rbind,
            
            map(lamb_vals, 
                function(i) {
                  (
                    cbind(
                      as.data.frame(map2( (draws_indiv %>% select(contains("theta"))),  
                                          (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,],
                                          `*`)),
                      lam_indiv[,i]
                    ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() - 
                      
                      cbind(
                        as.data.frame(map2( (draws_indiv %>% select(contains("theta")))[,-mk],  
                                            (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,][,-mk],
                                            `*`)),
                        lam_indiv[,i]
                      ) %>% rowSums() %>% inv_cloglog() %>% as.matrix() 
                    
                  ) %>% post.summary.confint()
                  
                })
            
    )
  
  pred_indiv_mi = as.data.frame(pred_indiv_mi)
  pred_indiv_mi$patient2 = patient2ID
  pred_indiv_mi$interval = as.numeric(paste(unique(predict_chron$interval)))[1:nrow(pred_indiv_mi)]
  pred_indiv_mi
}

predictions_hum = do.call(rbind, 
                          c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 1)  )
)

predictions_press = do.call(rbind, 
                            c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 2)  )
)

predictions_temp = do.call(rbind, 
                           c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 3)  )
)

predictions_windgust = do.call(rbind, 
                               c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 4)  )
)

predictions_windspeed = do.call(rbind, 
                                c( map(unique(predict_chron$patient2), predict_indiv_mi, mi = 5)  )
)

palette.mike.jackie <- c("burlywood", 
                         "lightsteelblue4", 
                         "lightpink2", 
                         "#FFDB6D", 
                         "tan2", 
                         "#52854C")



ID2 = 21
gp1 = ggplot() +
  geom_line(data = predictions_hum %>% filter(patient2 == ID2), aes(x = interval, y = mean),
            color = palette.mike.jackie[1], size = 1 ) + 
  geom_ribbon(data = predictions_hum %>% filter(patient2 == ID2), aes(x = interval, ymin = `2.5%`, ymax = `97.5%`),
              fill = palette.mike.jackie[1], alpha = 0.2) + 
  geom_line(linetype = "dashed", data = predictions_hum %>% filter(patient2 == ID2), aes(x = interval, y = `2.5%`),
            color = palette.mike.jackie[1], size = 0.7 ) +
  geom_line(linetype = "dashed", data = predictions_hum %>% filter(patient2 == ID2), aes(x = interval, y = `97.5%`),
            color = palette.mike.jackie[1] ) +  theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  ylab(expression(Delta*"Probability")) + xlab("Days between attacks") + 
  ggtitle(expression(Delta*"Probability with humidity")) + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                                                 axis.text.x = element_text(size = 14),
                                                                 axis.text.y = element_text(size = 14),  
                                                                 axis.title.x = element_text(size = 14),
                                                                 axis.title.y = element_text(size = 14),
                                                                 legend.text=element_text(size=12),
                                                                 legend.title=element_text(size=14),
  ) + 
  scale_x_continuous(breaks = seq(0, 8, by = 1))

gp2 = ggplot() +
  geom_line(data = predictions_press %>% filter(patient2 == ID2), aes(x = interval, y = mean),
            color = palette.mike.jackie[2], size = 1) +
  geom_ribbon(data = predictions_press %>% filter(patient2 == ID2), aes(x = interval, ymin = `2.5%`, ymax = `97.5%`),
              fill = palette.mike.jackie[2], alpha = 0.2) + 
  geom_line(linetype = "dashed", data = predictions_press %>% filter(patient2 == ID2), aes(x = interval, y = `2.5%`),
            color = palette.mike.jackie[2] ) +
  geom_line(linetype = "dashed", data = predictions_press %>% filter(patient2 == ID2), aes(x = interval, y = `97.5%`),
            color = palette.mike.jackie[2] ) +  theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ylab(expression(Delta*"Probability")) + xlab("Days between attacks") + 
  ggtitle(expression(Delta*"Probability with pressure")) + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                                                 axis.text.x = element_text(size = 14),
                                                                 axis.text.y = element_text(size = 14),  
                                                                 axis.title.x = element_text(size = 14),
                                                                 axis.title.y = element_text(size = 14),
                                                                 legend.text=element_text(size=12),
                                                                 legend.title=element_text(size=14),
  ) + 
  scale_x_continuous(breaks = seq(0, 8, by = 1))

gp3 = ggplot() +
  geom_line(data = predictions_temp %>% filter(patient2 == ID2), aes(x = interval, y = mean),
            color = palette.mike.jackie[3], size = 1) +
  geom_ribbon(data = predictions_temp %>% filter(patient2 == ID2), aes(x = interval, ymin = `2.5%`, ymax = `97.5%`),
              fill = palette.mike.jackie[3], alpha = 0.2) + 
  geom_line(linetype = "dashed", data = predictions_temp %>% filter(patient2 == ID2), aes(x = interval, y = `2.5%`),
            color = palette.mike.jackie[3] ) +
  geom_line(linetype = "dashed", data = predictions_temp %>% filter(patient2 == ID2), aes(x = interval, y = `97.5%`),
            color = palette.mike.jackie[3] ) +  theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  ylab(expression(Delta*"Probability")) + xlab("Days between attacks") +
  ggtitle(expression(Delta*"Probability with temperature")) + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                                                    axis.text.x = element_text(size = 14),
                                                                    axis.text.y = element_text(size = 14),  
                                                                    axis.title.x = element_text(size = 14),
                                                                    axis.title.y = element_text(size = 14),
                                                                    legend.text=element_text(size=12),
                                                                    legend.title=element_text(size=14),
  ) + 
  scale_x_continuous(breaks = seq(0, 8, by = 1))

gp4 = ggplot() +
  geom_line(data = predictions_windgust %>% filter(patient2 == ID2), aes(x = interval, y = mean),
            color = palette.mike.jackie[4], size = 1) +
  geom_ribbon(data = predictions_windgust %>% filter(patient2 == ID2), aes(x = interval, ymin = `2.5%`, ymax = `97.5%`),
              fill = palette.mike.jackie[4], alpha = 0.2) + 
  geom_line(linetype = "dashed", data = predictions_windgust %>% filter(patient2 == ID2), aes(x = interval, y = `2.5%`),
            color = palette.mike.jackie[4] ) +
  geom_line(linetype = "dashed", data = predictions_windgust %>% filter(patient2 == ID2), aes(x = interval, y = `97.5%`),
            color = palette.mike.jackie[4] ) +  theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  ylab(expression(Delta*"Probability")) + xlab("Days between attacks") + 
  ggtitle(expression(Delta*"Probability with wind gust")) + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                                                  axis.text.x = element_text(size = 14),
                                                                  axis.text.y = element_text(size = 14),  
                                                                  axis.title.x = element_text(size = 14),
                                                                  axis.title.y = element_text(size = 14),
                                                                  legend.text=element_text(size=12),
                                                                  legend.title=element_text(size=14),
  ) + 
  scale_x_continuous(breaks = seq(0, 8, by = 1))

gp5 = ggplot() +
  geom_line(data = predictions_windspeed %>% filter(patient2 == ID2), aes(x = interval, y = mean),
            color = palette.mike.jackie[5], size = 1) + 
  geom_ribbon(data = predictions_windspeed %>% filter(patient2 == ID2), aes(x = interval, ymin = `2.5%`, ymax = `97.5%`),
              fill = palette.mike.jackie[5], alpha = 0.2) + 
  geom_line(linetype = "dashed", data = predictions_windspeed %>% filter(patient2 == ID2), aes(x = interval, y = `2.5%`),
            color = palette.mike.jackie[5] ) +
  geom_line(linetype = "dashed", data = predictions_windspeed %>% filter(patient2 == ID2), aes(x = interval, y = `97.5%`),
            color = palette.mike.jackie[5] ) +  theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  ylab(expression(Delta*"Probability")) + xlab("Days between attacks") + theme_bw() +
  ggtitle(expression(Delta*"Probability with wind speed")) + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                                                   axis.text.x = element_text(size = 14),
                                                                   axis.text.y = element_text(size = 14),  
                                                                   axis.title.x = element_text(size = 14),
                                                                   axis.title.y = element_text(size = 14),
                                                                   legend.text=element_text(size=12),
                                                                   legend.title=element_text(size=14),
  ) + 
  scale_x_continuous(breaks = seq(0, 8, by = 1))

gp6 = ggplot() +
  geom_line(data = predictions %>% filter(patient2 == ID2), aes(x = interval, y = mean),
            color = palette.mike.jackie[6], size = 1) + 
  geom_ribbon(data = predictions %>% filter(patient2 == ID2), aes(x = interval, ymin = `2.5%`, ymax = `97.5%`),
              fill = palette.mike.jackie[6], alpha = 0.2) + 
  geom_line(linetype = "dashed", data = predictions %>% filter(patient2 == ID2), aes(x = interval, y = `2.5%`),
            color = palette.mike.jackie[6] ) +
  geom_line(linetype = "dashed", data = predictions %>% filter(patient2 == ID2), aes(x = interval, y = `97.5%`),
            color = palette.mike.jackie[6] ) +  theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  ylab("Probability of attack") + xlab("Days between attacks") + theme_bw() +
  ggtitle("Overall probability") + theme(plot.title = element_text(hjust = 0.5, size = 16),
                                         axis.text.x = element_text(size = 14),
                                         axis.text.y = element_text(size = 14),  
                                         axis.title.x = element_text(size = 14),
                                         axis.title.y = element_text(size = 14),
                                         legend.text=element_text(size=12),
                                         legend.title=element_text(size=14),
  ) + 
  scale_x_continuous(breaks = seq(0, 8, by = 1))


png("./final results and code/generated_tables_plots/change_breakdown_one_patient.png", units="in", width=10, height=14, res=600)


plots_arranged = ggarrange(
  plotlist = list(gp1, gp2, gp3, gp4, gp5, gp6),
  nrow = 3, ncol = 2,
  labels = LETTERS[1:6],
  legend = "bottom",
  common.legend = TRUE
) 

plots_arranged 

dev.off()



rm(list = ls())

# RMSE tables for 3 models X 2 data sets = 6 RMSE's
RMSE_dist_list = vector(mode = "list", length = 6)
names(RMSE_dist_list) = c("epis_full", "chron_full", "epis_FE", "chron_FE", "epis_null", "chron_null")

# processed data sets: chronics and episodics
chronics = readRDS("chronics.Rds")
episodics = readRDS("episodics.Rds")
# transformed data used for modeling
chronics_transformed = readRDS("final results and code/chronics_transformed_data.Rds")
episodics_transformed = readRDS("final results and code/episodics_transformed_data.Rds")


rmse_full = function(chronics, episodics, chronics_transformed, episodics_transformed, chron.draws, epis.draws) {
  
  # helper functions 
  inv_cloglog = function(nu) {
    1 - exp(-exp(nu))
  }
  
  get_draws_indiv = function(draws, patient2ID) {
    cbind(draws$thetaf, draws$thetar.random[[patient2ID]], draws$lambda.draws)
  }
  
  get_RMSE = function(x, unconverted.draws, which_set) {
    
    predict_ci_indiv = function(patient2ID) {
      data_indiv = x %>% filter(patient2 == patient2ID)
      if(which_set == "chron") {
        data_indiv = data_indiv %>% select(patient2, interval, event, humidityVarSt, pressureVarSt, temperatureVarSt, 
                                         windGustVarSt, windSpeedVarSt, humidityVarSt_sq, temperatureVarSt_sq,
                                         windSpeedVarSt_sq, pressureVarSt_sq, windGustVarSt_sq)
      } else if(which_set == "epis") {
        data_indiv = data_indiv %>% select(patient2, interval, event, contains("VarSt"))
      }
      
      draws_indiv = get_draws_indiv(draws = unconverted.draws, patient2ID = patient2ID)
      
      lamb_vals = 1:length(unique(data_indiv$interval))
      lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
      
      draw_ci = function(pi) {
        rbinom(n = length(pi), size = 1, prob = pi)
      }
      
      p_i_indiv = 
        do.call(cbind,
                map(1:nrow(data_indiv), 
                    function(i) {
                      cbind(
                        as.data.frame(map2( draws_indiv %>% select(contains("theta")),  
                                            (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,],
                                            `*`)),
                        lam_indiv[,data_indiv$interval[i]]
                      ) %>% rowSums() %>% inv_cloglog() %>% as.matrix()
                    })) %>% as.data.frame()
      
      
      predicted_ci_indiv = do.call(cbind, map(p_i_indiv, draw_ci))
      predicted_ci_indiv
    }
    
    predictions_ci = do.call(cbind, map(unique(x$patient2), predict_ci_indiv))
    
    bayesRMSE = function(predictions_ci, observed_ci) {
      rmse = function(row_predictions_ci) sum((row_predictions_ci - observed_ci)^2) / length(observed_ci)
      RMSE_dist = apply(predictions_ci, 1, rmse)
      return(RMSE_dist)
    }
    
    RMSE_dist_est = bayesRMSE(predictions_ci = predictions_ci, observed_ci = x$event)
    RMSE_dist_est
    
  }
  
  unconvert_scale.chron = function(draws, data) {
    lambda.draws = draws[,grepl("lambda", colnames(draws))]
    thetaf.draws = draws[,grepl("f", colnames(draws))]
    theta.draws = draws[,grepl("theta\\[", colnames(draws))]
    theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
    patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
    patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                   paste0("theta_random", pattern, "2]")))) }
    theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
    theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}), 
                                          function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
    return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
  }
  
  unconvert_scale.epis = function(draws, data) {
    lambda.draws = draws[,grepl("lambda", colnames(draws))]
    thetaf.draws = draws[,grepl("f", colnames(draws))]
    theta.draws = draws[,grepl("theta\\[", colnames(draws))]
    theta.random.draws = draws[,grepl("theta_random", colnames(draws))]
    patient.pattern = map(1:length(unique(data$patient2)), function(id) paste0("[", id, ","))
    patient.colnames = function(patient.pattern) { unlist(map(patient.pattern, function(pattern) c(paste0("theta_random", pattern, "1]"),
                                                                                                   paste0("theta_random", pattern, "2]"),
                                                                                                   paste0("theta_random", pattern, "3]"),
                                                                                                   paste0("theta_random", pattern, "4]"),
                                                                                                   paste0("theta_random", pattern, "5]")))) }
    theta.random.draws.reorder = theta.random.draws[,patient.colnames(patient.pattern)]
    theta.random.draws.reorder.list = map(map(patient.pattern, function(x) { paste0("\\", x)}), 
                                          function(x) { theta.random.draws.reorder[,grepl(x, colnames(theta.random.draws.reorder))] })
    return(list(thetaf = thetaf.draws, thetar = theta.draws, thetar.random = theta.random.draws.reorder.list, lambda.draws = lambda.draws))
  }
  
  unconverted.draws.chron = unconvert_scale.chron(draws = chron.draws, data = chronics)
  unconverted.draws.epis = unconvert_scale.epis(draws = epis.draws, data = episodics)
  
  colnames(chronics_transformed) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                  "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
  chronics_transformed = as.data.frame(chronics_transformed)
  chronics_transformed = cbind(chronics_transformed, chronics$patient2, chronics$interval, chronics$event)
  colnames(chronics_transformed)[11:13] = c("patient2", "interval", "event")
  
  colnames(episodics_transformed) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                  "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
  episodics_transformed = as.data.frame(episodics_transformed)
  episodics_transformed = cbind(episodics_transformed, episodics$patient2, episodics$interval, episodics$event)
  colnames(episodics_transformed)[11:13] = c("patient2", "interval", "event")
  
  
  list("epis" = get_RMSE(x = episodics_transformed, unconverted.draws = unconverted.draws.epis, which_set = "epis"),
       "chron" = get_RMSE(x = chronics_transformed, unconverted.draws = unconverted.draws.chron, which_set = "chrom"))
  
}


# for full models
chron.fit_full = readRDS("./final results and code/cloglog_full_chronics.Rds")
chron.draws_full = as.data.frame(chron.fit_full, pars = c("thetaf", "theta", "theta_random", "lambda"))

epis.fit_full = readRDS("./final results and code/cloglog_full_episodics.Rds")
epis.draws_full = as.data.frame(epis.fit_full, pars = c("thetaf", "theta", "theta_random", "lambda"))


RMSE_set = rmse_full(chronics = chronics, episodics = episodics, 
                chronics_transformed = chronics_transformed, episodics_transformed = episodics_transformed, 
                chron.draws = chron.draws_full, epis.draws = epis.draws_full)
RMSE_dist_list$epis_full = RMSE_set$epis
RMSE_dist_list$chron_full = RMSE_set$chron


rmse_FE = function(chronics, episodics, chronics_transformed, episodics_transformed, chron.draws, epis.draws) {
  
  # helper functions 
  inv_cloglog = function(nu) {
    1 - exp(-exp(nu))
  }
  
  get_draws_indiv = function(draws, patient2ID) {
    cbind(draws$thetaf, draws$lambda.draws)
  }
  
  get_RMSE = function(x, unconverted.draws, which_set) {
    
    predict_ci_indiv = function(patient2ID) {
      data_indiv = x %>% filter(patient2 == patient2ID)
      if(which_set == "chron") {
        data_indiv = data_indiv %>% select(patient2, interval, event, humidityVarSt, pressureVarSt, temperatureVarSt, 
                                           windGustVarSt, windSpeedVarSt, humidityVarSt_sq, temperatureVarSt_sq,
                                           windSpeedVarSt_sq, pressureVarSt_sq, windGustVarSt_sq)
      } else if(which_set == "epis") {
        data_indiv = data_indiv %>% select(patient2, interval, event, contains("VarSt"))
      }
      
      draws_indiv = get_draws_indiv(draws = unconverted.draws, patient2ID = patient2ID)
      
      lamb_vals = 1:length(unique(data_indiv$interval))
      lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
      
      draw_ci = function(pi) {
        rbinom(n = length(pi), size = 1, prob = pi)
      }
      
      p_i_indiv = 
        do.call(cbind,
                map(1:nrow(data_indiv), 
                    function(i) {
                      cbind(
                        as.data.frame(map2( draws_indiv %>% select(contains("theta")),  
                                            (data_indiv %>% ungroup() %>% select(contains("VarSt")))[i,],
                                            `*`)),
                        lam_indiv[,data_indiv$interval[i]]
                      ) %>% rowSums() %>% inv_cloglog() %>% as.matrix()
                    })) %>% as.data.frame()
      
      
      predicted_ci_indiv = do.call(cbind, map(p_i_indiv, draw_ci))
      predicted_ci_indiv
    }
    
    predictions_ci = do.call(cbind, map(unique(x$patient2), predict_ci_indiv))
    
    bayesRMSE = function(predictions_ci, observed_ci) {
      rmse = function(row_predictions_ci) sum((row_predictions_ci - observed_ci)^2) / length(observed_ci)
      RMSE_dist = apply(predictions_ci, 1, rmse)
      return(RMSE_dist)
    }
    
    RMSE_dist_est = bayesRMSE(predictions_ci = predictions_ci, observed_ci = x$event)
    RMSE_dist_est
    
  }
  
  unconvert_scale.chron = function(draws, data) {
    lambda.draws = draws[,grepl("lambda", colnames(draws))]
    thetaf.draws = draws[,grepl("f", colnames(draws))]
    return(list(thetaf = thetaf.draws, lambda.draws = lambda.draws))
  }
  
  unconvert_scale.epis = function(draws, data) {
    lambda.draws = draws[,grepl("lambda", colnames(draws))]
    thetaf.draws = draws[,grepl("f", colnames(draws))]
    return(list(thetaf = thetaf.draws, lambda.draws = lambda.draws))
  }
  
  unconverted.draws.chron = unconvert_scale.chron(draws = chron.draws, data = chronics)
  unconverted.draws.epis = unconvert_scale.epis(draws = epis.draws, data = episodics)
  
  colnames(chronics_transformed) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                                     "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
  chronics_transformed = as.data.frame(chronics_transformed)
  chronics_transformed = cbind(chronics_transformed, chronics$patient2, chronics$interval, chronics$event)
  colnames(chronics_transformed)[11:13] = c("patient2", "interval", "event")
  
  colnames(episodics_transformed) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                                      "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
  episodics_transformed = as.data.frame(episodics_transformed)
  episodics_transformed = cbind(episodics_transformed, episodics$patient2, episodics$interval, episodics$event)
  colnames(episodics_transformed)[11:13] = c("patient2", "interval", "event")
  
  
  list("epis" = get_RMSE(x = episodics_transformed, unconverted.draws = unconverted.draws.epis, which_set = "epis"),
       "chron" = get_RMSE(x = chronics_transformed, unconverted.draws = unconverted.draws.chron, which_set = "chrom"))
  
}


# for FE models
chron.fit_FE = readRDS("./final results and code/cloglog_FE_chronics.Rds")
chron.draws_FE = as.data.frame(chron.fit_FE, pars = c("thetaf", "lambda"))

epis.fit_FE = readRDS("./final results and code/cloglog_FE_episodics.Rds")
epis.draws_FE = as.data.frame(epis.fit_FE, pars = c("thetaf", "lambda"))


RMSE_set = rmse_FE(chronics = chronics, episodics = episodics, 
                chronics_transformed = chronics_transformed, episodics_transformed = episodics_transformed, 
                chron.draws = chron.draws_FE, epis.draws = epis.draws_FE)
RMSE_dist_list$epis_FE = RMSE_set$epis
RMSE_dist_list$chron_FE = RMSE_set$chron


rmse_null = function(chronics, episodics, chronics_transformed, episodics_transformed, chron.draws, epis.draws) {
  
  # helper functions 
  inv_cloglog = function(nu) {
    1 - exp(-exp(nu))
  }
  
  get_draws_indiv = function(draws, patient2ID) {
    cbind(draws$lambda.draws)
  }
  
  get_RMSE = function(x, unconverted.draws, which_set) {
    
    predict_ci_indiv = function(patient2ID) {
      data_indiv = x %>% filter(patient2 == patient2ID)
      if(which_set == "chron") {
        data_indiv = data_indiv %>% select(patient2, interval, event, humidityVarSt, pressureVarSt, temperatureVarSt, 
                                           windGustVarSt, windSpeedVarSt, humidityVarSt_sq, temperatureVarSt_sq,
                                           windSpeedVarSt_sq, pressureVarSt_sq, windGustVarSt_sq)
      } else if(which_set == "epis") {
        data_indiv = data_indiv %>% select(patient2, interval, event, contains("VarSt"))
      }
      
      draws_indiv = get_draws_indiv(draws = unconverted.draws, patient2ID = patient2ID)
      
      lamb_vals = 1:length(unique(data_indiv$interval))
      lam_indiv = (draws_indiv %>% select(contains("lambda")))[,lamb_vals, drop = FALSE]
      
      draw_ci = function(pi) {
        rbinom(n = length(pi), size = 1, prob = pi)
      }
      
      p_i_indiv = 
        do.call(cbind,
                map(1:nrow(data_indiv), 
                    function(i) {
                      (
                        lam_indiv[,data_indiv$interval[i], drop = FALSE]
                      ) %>% rowSums() %>% inv_cloglog() %>% as.matrix()
                    })) %>% as.data.frame()
      
      
      predicted_ci_indiv = do.call(cbind, map(p_i_indiv, draw_ci))
      predicted_ci_indiv
    }
    
    predictions_ci = do.call(cbind, map(unique(x$patient2), predict_ci_indiv))
    
    bayesRMSE = function(predictions_ci, observed_ci) {
      rmse = function(row_predictions_ci) sum((row_predictions_ci - observed_ci)^2) / length(observed_ci)
      RMSE_dist = apply(predictions_ci, 1, rmse)
      return(RMSE_dist)
    }
    
    RMSE_dist_est = bayesRMSE(predictions_ci = predictions_ci, observed_ci = x$event)
    RMSE_dist_est
    
  }
  
  unconvert_scale.chron = function(draws, data) {
    lambda.draws = draws[,grepl("lambda", colnames(draws))]
    return(list(lambda.draws = lambda.draws))
  }
  
  unconvert_scale.epis = function(draws, data) {
    lambda.draws = draws[,grepl("lambda", colnames(draws))]
    return(list(lambda.draws = lambda.draws))
  }
  
  unconverted.draws.chron = unconvert_scale.chron(draws = chron.draws, data = chronics)
  unconverted.draws.epis = unconvert_scale.epis(draws = epis.draws, data = episodics)
  
  colnames(chronics_transformed) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                                     "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
  chronics_transformed = as.data.frame(chronics_transformed)
  chronics_transformed = cbind(chronics_transformed, chronics$patient2, chronics$interval, chronics$event)
  colnames(chronics_transformed)[11:13] = c("patient2", "interval", "event")
  
  colnames(episodics_transformed) = c("humidityVarSt", "pressureVarSt", "temperatureVarSt", "windGustVarSt", "windSpeedVarSt",
                                      "humidityVarSt_sq", "pressureVarSt_sq", "temperatureVarSt_sq", "windGustVarSt_sq", "windSpeedVarSt_sq")
  episodics_transformed = as.data.frame(episodics_transformed)
  episodics_transformed = cbind(episodics_transformed, episodics$patient2, episodics$interval, episodics$event)
  colnames(episodics_transformed)[11:13] = c("patient2", "interval", "event")
  
  
  list("epis" = get_RMSE(x = episodics_transformed, unconverted.draws = unconverted.draws.epis, which_set = "epis"),
       "chron" = get_RMSE(x = chronics_transformed, unconverted.draws = unconverted.draws.chron, which_set = "chrom"))
  
}


# for null models
chron.fit_null = readRDS("./final results and code/cloglog_null_chronics.Rds")
chron.draws_null = as.data.frame(chron.fit_null, pars = c("lambda"))

epis.fit_null = readRDS("./final results and code/cloglog_null_episodics.Rds")
epis.draws_null = as.data.frame(epis.fit_null, pars = c("lambda"))

RMSE_set = rmse_null(chronics = chronics, episodics = episodics, 
                chronics_transformed = chronics_transformed, episodics_transformed = episodics_transformed, 
                chron.draws = chron.draws_null, epis.draws = epis.draws_null)
RMSE_dist_list$epis_null = RMSE_set$epis
RMSE_dist_list$chron_null = RMSE_set$chron

saveRDS(RMSE_dist_list, "./final results and code/generated_tables_plots/RMSE.Rds")


rm(list = ls())


