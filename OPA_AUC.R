### Trish Milletich 
#Flow Cytometry IC50 for Mandi and Brianna 

library(ggplot2)
#library(drc)
library(ggpubr)
library(DescTools)

setwd("C:/Users/PMilletich/OneDrive - University of Maryland School of Medicine/Desktop/OPA Analysis & data/")

#Create Dilution table to match Letters with concentration 
Dilution_values = data.frame("Dilution" = c("A", "B", "C", 
                                            "D", "E", "F"), 
                             "Concentration" = c(1/8, 1/16, 1/32, 1/64, 1/128, 1/256),
                             "Descriptions" = c("1:008", "1:016", "1:032", 
                                                "1:064", "1:128", "1:256"))
#Create a empty list to add dilution titer for all patients 
Dilution_list = c()
Bar_data = data.frame()
R2_list = c()
current_data = "July_2"
Model_data = data.frame()

#Load in file names and list of subjects for each file
OG_data = read.csv("Experiment 02-JULY-2024.csv")
Patient_IDlist = c("01-101 01", "01-101 02", "01-107 01", "01-107 02")

#create index cikumn 
OG_data$i = 1:nrow(OG_data)

#Calculations = =(Geometric Mean*RFP+)/10000
#Subset for the RFP + data 
OG_data_RFP = subset(OG_data, grepl("RFP +", OG_data$Name) == T)
current_pat = Patient_IDlist[1]
#Cycle through all the subjects and average the three replicates by dilution 
for (current_pat in Patient_IDlist) {
  pat_subset = subset(OG_data, grepl(current_pat, OG_data$Name))
  
  pat_averages = data.frame()
  pat_all = data.frame()
  for (Dilution in c("A", "B", "C", "D", "E", "F")) {
    search_list = paste(current_pat, Dilution, sep = "-")
    search_list <- unique(grep(paste(search_list,collapse="|"), 
                               pat_subset$Name, value=TRUE))  
    
    Dilution_subset = subset(pat_subset, pat_subset$Name %in% search_list)
    RFP_subset = subset(Dilution_subset, grepl("RFP +", Dilution_subset$Name) == T)
    
    Metric_list = c()
    for (current_i in unique(RFP_subset$i)) {
      i_subset = subset(OG_data, OG_data$i %in% current_i:(current_i+1))
      Geometric_Mean = as.numeric(i_subset[i_subset$i != current_i,"Statistic"])
      if (grepl("Geometric Mean", i_subset[i_subset$i != current_i,"Name"]) == F) {
        print("Bad News")
        break
      }
      RFP = as.numeric(i_subset[i_subset$i == current_i,"Statistic"])
      
      cell_count_name = strsplit(i_subset[i_subset$i == current_i,"Name"], "/")[[1]][1]
      Max_cells = as.numeric(OG_data[OG_data$Name == cell_count_name,4])
      
      if (Max_cells != 10000) {
        print(paste(cell_count_name, Max_cells))
      }
      
      Metric = (Geometric_Mean*RFP)/Max_cells
      Metric_list = c(Metric_list, Metric)
    }
    
    avg_Metric = ifelse(length(Metric_list) == 0, NA, mean(Metric_list, na.rm = T))
    
    current_row = data.frame("Patient" = current_pat,
                             "Dilution" = Dilution,
                             "Metric" = avg_Metric)
    pat_averages = rbind(pat_averages, current_row)
    
    if (length(Metric_list) != 0) {
      #Save raw total data to new data.frame 
      current_row = data.frame("Patient" = current_pat,
                               "Dilution" = Dilution,
                               "Metric" = Metric_list)
      pat_all = rbind(pat_all, current_row)
    }
    
  }
  
  #Combined the numeric data with Dilution values to convert letters to dilutions 
  pat_averages = merge(pat_averages, Dilution_values)
  pat_all = merge(pat_all, Dilution_values)
  
  #create a scaled metric so the highest value is 1 
  pat_averages$Metric_Scaled = pat_averages$Metric/max(pat_averages$Metric, na.rm = T)
  pat_all$Metric_Scaled = pat_all$Metric/max(pat_all$Metric, na.rm = T)
  pat_all_OG = pat_all
  #Fit linear model on "right" half of the data if it's not a smooth gradient 
  Max_Dilution = subset(pat_averages, pat_averages$Metric == max(pat_averages$Metric, na.rm = T))$Dilution
  if (Max_Dilution != "A") {
    letter_list = c("A", "B", "C", "D", "E", "F")
    letter_list = letter_list[match(Max_Dilution, letter_list):6]
    pat_all = subset(pat_all, pat_all$Dilution %in% letter_list)
  }
  
  #fit a linear regression model 
  fit_1 <- lm(log10(Concentration) ~ Metric_Scaled, data = pat_all)
  #Predict the log10 Concentration where the scaled metric = 0.5
  Predicted = predict(fit_1, data.frame(Metric_Scaled = 0.5))
  #Convert from log10 to regular number
  Predicted_Dilution = 10^(Predicted); Predicted_Dilution
  
  #Find the R2 of the linear model 
  #https://www.dataquest.io/blog/statistical-learning-for-predictive-modeling-r/ 
  fit_1_sum = summary(fit_1)
  R2 = fit_1_sum$r.squared
  
  #Plot the linear model 
  linear_plot = ggplot(data = pat_all, aes(x = log10(Concentration), y = Metric_Scaled,
                                           color = Descriptions)) +
    scale_color_manual(breaks = c("1:008", "1:016", "1:032",
                                  "1:064", "1:128", "1:256"), 
                       values = c("darkred", "brown3", "firebrick1",
                                  "coral", "pink", "mistyrose"))+ 
    stat_smooth(method = "lm", col = "dodgerblue3",formula = y ~ x) +
    geom_point() +
    theme(panel.background = element_rect(fill = "white"),
          axis.line.x=element_line(),
          axis.line.y=element_line()) +
    scale_x_reverse() + 
    annotate("text", label = paste("R2 = ", round(R2, 3), sep = ""), y = 1, x = -2) + 
    theme(legend.position = "none"); linear_plot
  
  
  Description_value = ifelse(round(round(1/Predicted_Dilution)) < 100, paste("1:0", round(1/Predicted_Dilution), sep=""),
                             paste("1:", round(1/Predicted_Dilution), sep=""))
  
  AUC = AUC(pat_averages$Concentration, pat_averages$Metric)
  
  pat_averages_1 = pat_averages
  #Add the new predicted data 
  average_data = data.frame("Dilution" = "Predicted", 
                            "Patient" = current_pat, 
                            "Metric" = (max(pat_averages$Metric, na.rm = T)/2), 
                            "Concentration" = Predicted_Dilution,
                            "Descriptions" = Description_value,
                            "Metric_Scaled" = 0.5)
  pat_averages = rbind(pat_averages, 
                       average_data)
  
  AUC_plot = ggplot(pat_averages, aes(x = Descriptions, y = Metric, 
                                      fill = Descriptions,
                                      group = Patient)) + 
    geom_area(data = pat_averages_1, 
              aes(x = Descriptions, y = Metric),
              fill = "grey", alpha = 0.5) + 
    geom_bar(stat = "identity", position = "dodge", width =0.75) + 
    geom_hline(yintercept = max(pat_averages$Metric, na.rm = T)/2) + 
    theme_bw() + 
    annotate("text", label = paste("IC50=", round(Predicted_Dilution, 3), "\n",
                                   "AUC=", round(AUC,3),
                                   sep = ""),
             x = 6, y = max(pat_averages$Metric,na.rm = T)-(max(pat_averages$Metric, na.rm= T)/4), size = 3.5, colour = "grey34")+ 
    scale_color_manual(breaks = c("1:008", "1:016", "1:032",
                                  "1:064", "1:128", "1:256", 
                                  Description_value), 
                       values = c("darkred", "brown3", "firebrick1",
                                  "coral", "pink", "mistyrose", "blue")) + 
    scale_fill_manual(breaks = c("1:008", "1:016", "1:032",
                                 "1:064", "1:128", "1:256", 
                                 Description_value), 
                      values = c("darkred", "brown3", "firebrick1",
                                 "coral", "pink", "mistyrose", "blue")) + 
    theme(legend.position = "bottom"); suppressWarnings(AUC_plot)
  
  
  Combined = ggarrange(linear_plot, AUC_plot, nrow = 2, heights = c(1,2))
  Combined = annotate_figure(Combined, current_pat)
  jpeg(paste("./Figures/",current_data, "_", gsub(" ", ".", current_pat), "_3.jpeg", sep = ""), 
       res = 400, height = 3500, width = 3500)
  plot(Combined)
  dev.off()
  
  Dilution_list = c(Dilution_list, Predicted_Dilution)
  pat_averages$Date = current_data
  Bar_data = rbind(Bar_data, pat_averages)
  R2_list = c(R2_list, R2)
  
  
  model_row = data.frame("date" = current_data, 
                         "Subject_ID" = current_pat, 
                         "IC50" = round(Predicted_Dilution, 6), 
                         "AUC" = round(AUC,6))
  
  Model_data = rbind(Model_data, model_row)
}
