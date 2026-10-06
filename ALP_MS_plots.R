

### ALP

#Fig 3/D
##Effects of ALP on egg maturation

library(readxl)
ALP_treatment_repeat <- read_excel("~/Projekt2/ALP_treatment_repeat.xlsx", sheet = "Results_sum")# col_types = c("numeric",  # "text", "text", "numeric"))
View(ALP_treatment_repeat)

library(ggplot2)

library(ggsignif)
library(dplyr)
library(ggpubr)

library(ggnewscale)


ALP_treatment_repeat$treatment <- factor(
  ALP_treatment_repeat$treatment,
  levels = c("DMSO", "ALP")
)



# Individual backround colours for facets
background_colors <- data.frame(
  tray = unique(ALP_treatment_repeat$tray),  
  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
  fill_color = c("#ffe156", "#ff6f61", "#ffcab1", "#a8d5ba"))
  
  pvals <- data.frame(
    tray= c("1", "2", "3", "4"),
    group1 = "DMSO",
    group2 = "ALP",
    p = c(0.00009584, 0.00009584, 0.289, 0.612),
    y.position = c(68, 68, 68, 68)  
  )
  
pvals$label <- ifelse(pvals$p < 0.001,
                      "p < 0.001",
                      paste0("p = ", formatC(pvals$p, format = "f", digits = 3)))
  
  # Facetelt boxplot háttérszínekkel
ggplot(ALP_treatment_repeat, aes(x = treatment, y = egg_maturation_time)) +
    geom_rect(data = background_colors, aes(
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color
    ), alpha = 0.7, inherit.aes = FALSE) +
    scale_fill_identity()+
    new_scale_fill()+
    
    # Transparency settings
    geom_boxplot(aes(fill = treatment), outlier.shape = 16, size=0.8, color="black") + 
    labs(fill="Treatment")+
    stat_summary(
      fun.min = min,
      fun.max = max,
      geom = "errorbar",
      width = 0.6,       
      color = "black",
      lwd = 0.8
    )+
    scale_fill_manual(values=c("DMSO"="#737373","ALP"="#091edd"))+
    new_scale_fill()+
    geom_jitter(aes(fill = treatment),
                shape = 21, 
                width = 0.2, 
                size = 2.5, 
                color = "black", 
                stroke = 0.8,
                show.legend = FALSE)+
    scale_fill_manual(values = c("DMSO" = "black",
                                 "ALP" = "#00008B")) + 
    facet_wrap(~ tray, ncol=4, labeller=as_labeller(c("1"="D1", "2"="D7", "3"="D14", "4"="D21"))) +  
    expand_limits(y = max(ALP_treatment_repeat$egg_maturation_time)*1.06) +
    stat_pvalue_manual(pvals, label = "label", tip.length = 0.01, size = 5) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Effects of ALP on egg maturation",
         x = "Treatment",
         y = "Egg maturation time",
         color = "Treatment") +
    
    theme(
      panel.background = element_rect(fill = "white"),  
      strip.background = element_rect(fill = "gray90", colour = "black"), 
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin=unit(c(1,1,1,1),"mm"),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(colour="black",fill=NA,linewidth = 1),
      axis.title = element_text(face = "bold"),  # Bold axis titles
      axis.text = element_text(face = "bold"),    # Bold axis text
      axis.ticks = element_line(linewidth = 1),  # Customize tick line size
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  
  
ggsave("EMT_sum_FINAL.png", width = 20, height = 10, units = "cm")


  
hist(ALP_treatment_repeat$egg_maturation_time, breaks=20)
  
W1 <- ALP_treatment_repeat%>%filter(time=="D1-D4")
wilcox.test(W1$egg_maturation_time~W1$treatment)
  
W2 <- ALP_treatment_repeat%>%filter(time=="D7-D10")
wilcox.test(W2$egg_maturation_time~W2$treatment)
  
W3 <- ALP_treatment_repeat%>%filter(time=="D14-D17")
wilcox.test(W3$egg_maturation_time~W3$treatment)
  
W4 <- ALP_treatment_repeat%>%filter(time=="D21-D24")
wilcox.test(W4$egg_maturation_time~W4$treatment)
  
p.adjust(c(3.134e-05, 4.792e-05, 0.2168, 0.6122), method="BH")

#_____________________________________________________________________________________  
#
#_____________________________________________________________________________________


# Fig 3/E
##Effects of ALP on fecundity

library(readxl)
ALP_treatment_repeat <- read_excel("~/Projekt2/ALP_treatment_repeat.xlsx", sheet = "Results_sum")
View(ALP_treatment_repeat)

library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)


library(ggnewscale)

ALP_treatment_repeat$treatment <- factor(
  ALP_treatment_repeat$treatment,
  levels = c("DMSO", "ALP")
)


# Individual colors for facet backround
background_colors <- data.frame(
  tray = unique(ALP_treatment_repeat$tray),  
  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,  
  fill_color = c("#ffe156", "#ff6f61", "#ffcab1", "#a8d5ba") 
)


ALP_treatment_repeat$log_totaleggnumber <- log(ALP_treatment_repeat$totaleggnumber)

pvals <- data.frame(
  tray= c("1", "2", "3", "4"),
  group1 = "DMSO",
  group2 = "ALP",
  p = c(0.471, 0.662, 0.013,0.459),
  y.position = c(5, 5, 5, 5)  # minden celltype-hoz külön!
)

pvals$label <- ifelse(pvals$p < 0.001,
                      "p < 0.001",
                      paste0("p = ", formatC(pvals$p, format = "f", digits = 3)))

# Facetelt boxplot háttérszínekkel
ggplot(ALP_treatment_repeat, aes(x = treatment, y = log_totaleggnumber)) +
  geom_rect(data = background_colors, aes(
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color
  ), alpha = 0.7, inherit.aes = FALSE) +
  scale_fill_identity()+
  new_scale_fill()+
  
  # Átlátszóság állítása
  geom_boxplot(aes(fill = treatment), outlier.shape = 16, size=0.8, color="black") + 
  labs(fill="Treatment")+
  stat_summary(
    fun.min = min,
    fun.max = max,
    geom = "errorbar",
    width = 0.6,       
    color = "black",
    lwd = 0.8
  )+
  scale_fill_manual(values=c("DMSO"="#737373","ALP"="#091edd"))+
  new_scale_fill() +
  geom_jitter(aes(fill = treatment),
              shape = 21, 
              width = 0.2, 
              size = 2.5, 
              color = "black", 
              stroke = 0.8,
              show.legend = FALSE)+
  scale_fill_manual(values = c("DMSO" = "black",
                               "ALP" = "#00008B")) +
  facet_wrap(~ tray, ncol=4, labeller=as_labeller(c("1"="D1", "2"="D7", "3"="D14", "4"="D21"))) + 
  #geom_blank(aes(y = ratio_lim))+
  expand_limits(y = max(ALP_treatment_repeat$log_totaleggnumber)+0.5)+
  stat_pvalue_manual(pvals, label = "label", tip.length = 0.01, size = 5) +
  theme_minimal(base_size = 14) +
  labs(title= "Effects of ALP on fecundity",
       x= "Treatment",
       y= "Fecundity",
       color="Treatment")+
  theme(
    panel.background = element_rect(fill = "white"),  
    strip.background = element_rect(fill = "gray90", colour = "black"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin=unit(c(1,1,1,1),"mm"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(colour="black",fill=NA,linewidth = 1),
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),    # Bold axis text
    axis.ticks = element_line(linewidth = 1),  # Customize tick line size
    #legend.position = "top",  # Position legend at the top
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



ggsave("TEN_BIGGER.png", width = 15, height = 9, units = "cm", dpi=300)


hist(ALP_treatment_repeat$totaleggnumber)

library(dplyr)
W1 <- ALP_treatment_repeat%>%filter(tray=="1")
wilcox.test(W1$totaleggnumber~W1$treatment)

W2 <- ALP_treatment_repeat%>%filter(tray=="2")
wilcox.test(W2$totaleggnumber~W2$treatment)

W3 <- ALP_treatment_repeat%>%filter(tray=="3")
wilcox.test(W3$totaleggnumber~W3$treatment)

W4 <- ALP_treatment_repeat%>%filter(tray=="4")
wilcox.test(W4$totaleggnumber~W4$treatment)

p.adjust(c(0.3535, 0.6629, 0.003382, 0.2295), method="BH")

#_____________________________________________________________________________________  
#
#_____________________________________________________________________________________

#Fig 3/F
##Mortality curve

library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(lubridate)
library(patchwork)



ALP_3_0 <- read_excel(
  "Projekt3/ALP/ALP_X11_14_repeat3.0.xlsx",
  sheet="Munka3"
)



#______________________________________________________________________________
# HANDLING DATES
#______________________________________________________________________________

ALP_3_0 <- ALP_3_0 %>%
  mutate(
    observation_date = as.Date(observation_date)
  )
#______________________________________________________________________________
#Treatment column as character
#______________________________________________________________________________

ALP_3_0 <- ALP_3_0 %>%
  mutate(
    treatment = as.character(treatment)
  )


ALP_3_0 <- ALP_3_0 %>%
  mutate(
    treatment = factor(treatment,
                       levels = c("DMSO1", "ALP1", "DMSO4", "ALP4")))
#______________________________________________________________________________
#SURVIVAL DATABASE
#______________________________________________________________________________

survival_ALP_3_0 <- ALP_3_0 %>%
  
  arrange(ID, observation_date) %>%
  
  group_by(ID, treatment, treated) %>%
  
  summarise(
    
    death_date = if (any(survival == "x")) {
      min(observation_date[survival == "x"])
    } else {
      as.Date(NA)
    },
    
    last_date = max(observation_date),
    
    .groups = "drop"
  )  %>%  
  
  
  mutate(
    
    event = ifelse(is.na(death_date), 0, 1),
    
    final_date = coalesce(death_date, last_date),
    
    time_days = as.numeric(
      final_date - min(ALP_3_0$observation_date))
  )


fit_1 <- survfit(
  Surv(time_days, event) ~ treatment,
  data = filter(survival_ALP_3_0, treated == 1)
)

p1<-ggsurvplot(
  fit_1,
  data = filter(survival_ALP_3_0, treated == 1),
  palette = c("#737373","#091edd"),
  size = 1.5,
  risk.table = TRUE,
  pval = TRUE,
  pval.size = 7,
  conf.int = TRUE,
  surv.scale = "percent",
  ggtheme = theme(axis.line = element_line(linewidth = 1),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks = element_line(linewidth = 1),
                  axis.text = element_text(face = "bold"),
                  axis.title = element_text(face = "bold")),
  font.x = 18,
  font.y = 18,
  font.tickslab = 18,
  font.legend = 18,
  legend.title = "Treatment:",
  legend.labs = c("ALP 1x", "DMSO 1x"),
  xlab = "Days after SI",
  ylab = "Survival (%)",
  xlim = c(0,180)
)



fit_2 <- survfit(
  Surv(time_days, event) ~ treatment,
  data = filter(survival_ALP_3_0, treated == 2)
)

p2<-ggsurvplot(
  fit_2,
  data = filter(survival_ALP_3_0, treated == 2),
  palette = c("#737373","#091edd"),
  size = 1.5,
  linetype = c("dashed"),
  risk.table = TRUE,
  pval = TRUE,
  pval.size = 7,
  conf.int = TRUE,
  surv.scale = "percent",
  ggtheme = theme(axis.line = element_line(linewidth = 1),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks = element_line(linewidth = 1),
                  axis.text = element_text(face = "bold"),
                  axis.title = element_text(face = "bold")),
  font.tickslab = 18,
  font.legend = 18,
  font.x = 18,
  legend.title = "",
  legend.labs = c("ALP 2x", "DMSO 2x"),
  xlab = "Days after SI",
  ylab = "",
  xlim = c(0,180)
)

pic<-(p1$plot | p2$plot)

pic <- ggsave("survival_curve8.png", width = 35, height = 11, units = "cm")

#-----------------------------------------------------------
#statistics
#-----------------------------------------------------------

fit_1 <- survfit(
  Surv(time_days, event) ~ treatment,
  data = filter(survival_ALP_3_0, treated == 1)
)
summary(fit_1)$table

fit_2 <- survfit(
  Surv(time_days, event) ~ treatment,
  data = filter(survival_ALP_3_0, treated == 1)
)
summary(fit_2)$table

#_________________________________________________________________________________
#
#_________________________________________________________________________________

#Fig. 4/A
## Effects of ALP on cellular composition

library(readxl)
cell_ratio <- read_excel("~/Projekt2/Cell_counting_ALP_DMSO_treated.xlsx", sheet="Munka2")
View((cell_ratio))



library(ggplot2)
library(ggsignif)
library(ggpubr)
library(dplyr)
library(ggnewscale)

cell_ratio$treatment <-factor(
  cell_ratio$treatment,
  levels = c("DMSO", "ALP")
)


cell_ratio <- cell_ratio %>%
  group_by(celltype) %>%
  mutate(ratio_max = max(ratio),
         ratio_lim = ratio_max * 1.1)  


pvals <- data.frame(
  celltype = c("IE", "NB", "NC"),
  group1 = "ALP",
  group2 = "DMSO",
  p = c(0.222, 0.008, 0.010),
  y.position = c(0.485, 0.56, 1.14)  
)

pvals$label <- ifelse(pvals$p < 0.001,
                      "p < 0.001",
                      paste0("p = ", formatC(pvals$p, format = "f", digits = 3)))

ggplot(cell_ratio, aes(x=treatment, y=ratio))+
  geom_boxplot(aes(fill = treatment), outlier.shape = NA, size=0.8, color="black") +
  stat_summary(
    fun.min = min,
    fun.max = max,
    geom = "errorbar",
    width = 0.6,      
    color = "black",
    lwd = 0.8
  )+
  labs(fill="Treatment")+
  scale_fill_manual(values = c("ALP" = "#091edd",   
                               "DMSO" = "#737373")) + 
  new_scale_fill() +   # <-- ÚJ fill skála a jitterhez
  
 
  geom_jitter(
    aes(fill = treatment),
    shape = 21, 
    width = 0.2, 
    size = 4, 
    color = "black", 
    stroke = 0.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("ALP" = "#00008B",   
                               "DMSO" = "black")) + 
  #scale_fill_manual(values=c("ALP"="#091edd", "DMSO"="#737373"))+
  facet_wrap("celltype", ncol=4, scales="free_y", labeller = labeller(celltype = c(
    "IE" = "Interstitial Cells",
    "NB" = "Nematoblasts",
    "NC" = "Nurse Cells"
  )))+
  geom_blank(aes(y = ratio_lim))+
  stat_pvalue_manual(pvals, label = "label", tip.length = 0.01, size = 5) +
  labs(title= "Effects of ALP on cellular composition",
       x= "Treatment",
       y= "Cell ratio",
       color="Treatment")+
  theme_minimal(base_size = 19) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin=unit(c(1,1,1,1),"mm"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(colour="black",fill=NA,linewidth = 1),
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),    # Bold axis text
    axis.ticks = element_line(linewidth = 1),  # Customize tick line size
    legend.position = "left",  # Position legend at the top
    plot.title = element_text(
      hjust = 0.5, face = "bold", margin = margin(t = 80, r = 20, b = 20, l = 20), vjust = 20, size = 21)
  )

ggsave(filename = "~/Projekt2/pics/ALP_cell_counting15.png", width=20, height=15, units="cm")

#statistics

NB <- cell_ratio%>%filter(celltype=="NB")
wilcox.test(NB$ratio~NB$treatment)

IE <- cell_ratio%>%filter(celltype=="IE")
wilcox.test(IE$ratio~IE$treatment)

NC <- cell_ratio%>%filter(celltype=="NC")
wilcox.test(NC$ratio~NC$treatment)




