#Dexamethasone analysis project - Brechtje van der Molen
#clear working directory 
rm(list=ls()) 

#1 Set the working directory to a prefered location
setwd("C:/Users/brech/OneDrive - Maastricht University/MSP/PRA/PRA2026")

########open data set
#read.delim opens the tsv file (different than csv file.)
#header=TRUE tells R the first row is column names 
#check.names=FALSE prevents R from 'cleaning up' column  names. 
#1 Make sure the downloaded data is stored in the working directory you previously set
Cytokine.data<- read.delim("dex_cytokine.tsv", header=TRUE, check.names = FALSE)


#######add column that specifies whether patient received dexamethasone
#identify no-treatment patients
no.dex.patients <- c ("MVIR1-HS1", "MVIR1-HS107", "MVIR1-HS2", "MVIR1-HS23", "MVIR1-HS33", "MVIR1-HS43", "MVIR1-HS50", "MVIR1-HS53", "MVIR1-HS59", "MVIR1-HS72", "MVIR1-HS77", "MVIR1-HS78", "MVIR1-HS89", "MVIR1-HS91", "MVIR1-HS94")
length(no.dex.patients)
#add new column called Treatment that puts no.dex.patients in "no dex" category.ifesle(condition, yes, no) 
Cytokine.data$Treatment <- ifelse(Cytokine.data[["Patient id"]] %in% no.dex.patients, "No Dex", "Dex")

######calculate mean per cytokine per treatment group 
cytokine.cols <- colnames(Cytokine.data)
cytokine.cols <- cytokine.cols[!(cytokine.cols %in% c("Patient id", "Treatment"))]

#Check that all the data is stored numerically 
str(Cytokine.data[ , cytokine.cols])

#==================================================================
##################Data analysis###################################
#==================================================================

####Check whether each cytokine per treatment group has normal distribution(affects the statistical tests)
#Check kurtosis and skewness
kurtosis <- function (x){
    m<- mean(x)
    s<- sd(x)
    n<- length(x)
    sum((x-m)^4)/((n-1)*s^4)
}
    #Calculate skewness of each subsample 
    #sapply applies function(x) skewness() to each subsample dataframe.
    subsample.kurtosis <- sapply(Cytokine.data[ , cytokine.cols], function(x) kurtosis(x))
subsample.kurtosis

#Establish the skewness function
skewness<- function(x){
    m<- mean(x)
    s<- sd(x)
    n<- length(x)
    sum((x-m)^3)/((n-1)* s^3)
}

subsample.skewness <- sapply(Cytokine.data[ , cytokine.cols], function(x) skewness(x))
subsample.skewness

#Create data table in which the p values will be printed
normality.results <- data.frame(   #make a data frame consisting of variables:
  Cytokine = character(),          #cytokine variable with characters inside of it
  Treatment_Group = character(),   #Group variable which will distinguish treatment or not
  p.value = numeric()              #P vlaues where #P < 0.05 --> not normal distribution, and #P > 0.05 --> normal distribution
)
for(i in cytokine.cols){
       #Split treatment groups 
       Cytokine_No_dex <- as.numeric(Cytokine.data[Cytokine.data$Treatment == "No Dex", i])
       Cytokine_dex <- as.numeric(Cytokine.data[Cytokine.data$Treatment == "Dex", i])

       #remove NA from the data
       Cytokine_No_dex <- Cytokine_No_dex[!is.na(Cytokine_No_dex)]
       Cytokine_dex<- Cytokine_dex[!is.na(Cytokine_dex)]

       # only test if at least 3 unique values
       p_no_dex <- if (length(unique(Cytokine_No_dex)) >= 3) shapiro.test(Cytokine_No_dex)$p.value else NA
       p_dex    <- if (length(unique(Cytokine_dex)) >= 3) shapiro.test(Cytokine_dex)$p.value else NA

       #create the table that ads the results each i round.
       normality.results<-rbind(
       normality.results, #Keep old rows
       data.frame(Cytokine=i, Group="No dex", p.value=p_no_dex),
       data.frame(Cytokine = i, Group = "Dex",    p.value = p_dex))
}

#view the results
normality.results

#Results show that for all cytokine groups (except for protein C) either one or both treatment groups are not normal
#not normal distributions can be analyzed with Mann-whitney test
#Create data table in which the results will be printed
mw_results <- data.frame(
       Cytokine = character(),
       p.value = numeric()
)

for(i in cytokine.cols){
        #Split treatment groups 
       Cytokine_No_dex <- as.numeric(Cytokine.data[Cytokine.data$Treatment == "No Dex", i])
       Cytokine_dex <- as.numeric(Cytokine.data[Cytokine.data$Treatment == "Dex", i])

       #remove NA from the data
       Cytokine_No_dex <- Cytokine_No_dex[!is.na(Cytokine_No_dex)]
       Cytokine_dex<- Cytokine_dex[!is.na(Cytokine_dex)]

       #run Mann-Whitney test
       p<- wilcox.test(Cytokine_No_dex, Cytokine_dex)$p.value

       #add results to mw.data table
       mw_results<- rbind(
              mw_results,
              data.frame(Cytokine=i, p.value = p)
       )
}

#with so many tests run (18 cytokines), the multiple testing problem can occur (aka random false positives) 
#False discovery rate using BH curves the p-values so that false positives are less likely 
mw_results$FDR<- p.adjust(mw_results$p.value, method = "BH")
mw_results
#The siginificant cytokines are those with an adjusted p-value < 0.05
significant<- mw_results[mw_results$FDR < 0.05, ]
significant<- mw_results$Cytokine[mw_results$FDR < 0.05]         #shows the significant cytokines are IFN-y, IL-10, and IL-6
significant

#take the cytokine columns from cytokine data, calculate the mean per treatment
# split cytokine data into Dex and No Dex groups
split_data <- split(Cytokine.data[cytokine.cols], Cytokine.data$Treatment)

# calculate means per treatment
Cytokine.means <- lapply(split_data, function(df) sapply(df, mean, na.rm = TRUE))

# calculate SDs per treatment
Cytokine.sd <- lapply(split_data, function(df) sapply(df, sd, na.rm = TRUE))

#======================================================================================
########################Calculate values for significant cytokines#####################
#======================================================================================
###How many patients were treated vs untreaated
n_dex    <- sum(Cytokine.data$Treatment == "Dex")
n_nodex  <- sum(Cytokine.data$Treatment == "No Dex")

###Mean of treated vs untreated
#IL6
mean.IL6.Tr <- Cytokine.means$Dex["IL-6"]
mean.IL6.UTr <- Cytokine.means$"No Dex"["IL-6"]

#IL10
mean.IL10.Tr <- Cytokine.means$Dex["IL-10"]
mean.IL10.Utr <- Cytokine.means$"No Dex"["IL-10"]

#IFNg
mean.IFN.Tr <- Cytokine.means$Dex["IFN-gamma"]
mean.IFN.Utr <- Cytokine.means$"No Dex"["IFN-gamma"]

###SD of treated vs untreated
#IL6
sd.IL6.Tr<- Cytokine.sd$Dex["IL-6"]
sd.IL6.Utr<- Cytokine.sd$"No Dex"["IL-6"]

#IL10
sd.IL10.Tr<- Cytokine.sd$Dex["IL-10"]
sd.IL10.Utr<- Cytokine.sd$"No Dex"["IL-10"]

#IFNy
sd.IFN.Tr<- Cytokine.sd$Dex["IFN-gamma"]
sd.IFN.Utr<- Cytokine.sd$"No Dex"["IFN-gamma"]

###Standard error of treated vs untreated
#IL6
se_IL6.Tr  <- Cytokine.sd$Dex["IL-6"]    / sqrt(n_dex)
se_IL6.Utr <- Cytokine.sd$"No Dex"["IL-6"] / sqrt(n_nodex)

#IL10
se_IL10.Tr  <- Cytokine.sd$Dex["IL-10"]    / sqrt(n_dex)
se_IL10.Utr <- Cytokine.sd$"No Dex"["IL-10"] / sqrt(n_nodex)

#IFNgamma 
se_IFN.Tr  <- Cytokine.sd$Dex["IFN-gamma"]    / sqrt(n_dex)
se_IFN.Utr <- Cytokine.sd$"No Dex"["IFN-gamma"] / sqrt(n_nodex)

#================================
#STORE THE data for phython 
#Create summary table for significant cytokine bar graphs
significant.cytokine.summary <- data.frame(
  Cytokine = rep(c("IL-6", "IL-10", "IFN-gamma"), each = 2),
  Treatment = rep(c("No Dex", "Dex"), 3),

  Mean = c(
    mean.IL6.UTr, mean.IL6.Tr,
    mean.IL10.Utr, mean.IL10.Tr,
    mean.IFN.Utr, mean.IFN.Tr
  ),

  SEM = c(
    se_IL6.Utr, se_IL6.Tr,
    se_IL10.Utr, se_IL10.Tr,
    se_IFN.Utr, se_IFN.Tr
  )
)

#Add Mann-Whitney p-values and FDR values
significant.cytokine.summary <- merge(
  significant.cytokine.summary,
  mw_results,
  by = "Cytokine",
  all.x = TRUE
)

#Export results for Python
write.csv(
  significant.cytokine.summary,
  "significant_cytokine_summary.csv",
  row.names = FALSE
)

#=======================================================================
###########################Bar Graph with IL6##########################
#=======================================================================
#########plot bar graph with segment 
#need to pick the x locations of the bars
x.base <- 1:2
x.labels<- c("Untreated", "Treated")

#Establish highest and lowest y points 
####ADD standard deviation bars using arrow
dexam.error.bar.upper <- (mean.IL6.Tr + se_IL6.Tr)
dexam.error.bar.lower <- (mean.IL6.Tr - se_IL6.Tr)

nodexam.error.upper <- (mean.IL6.UTr + se_IL6.Utr)
nodexam.error.lower <- (mean.IL6.UTr - se_IL6.Utr)

png("Sig.Cytokine.IL6.plot.png", width = 5, height = 6, units = "in", res = 120) ## saves a png in the folder "images" within the working directory
#draw empty box plot
par(mar = c(4, 2, 5, 10))    #margins (bottom, left, top, right)
par(pin = c(2, 4))          #Recangular
par(xpd = NA)
plot(NA, xlim = c(1, 2),
     ylim=c(0.5, max(dexam.error.bar.upper, nodexam.error.upper)+50), 
     xaxt="n", xlab="", ylab="Concentration (pg/mL)",
     main="IL-6 Cytokine concentration\n by Dexamethasone Treatment")

axis(1, at=c(1.3, 1.7), labels=x.labels, las=2, cex.axis=0.8)

baseline<- par("usr")[3]

#design bars for untreated cytokine means
segments(x0=1.3, y0=baseline, x1=1.3, y1=mean.IL6.UTr,
         lwd=20,
         col="#f7b257",
         lend="butt")

#design bar for treated cytokine means
segments(x0=1.7, y0=baseline, x1=1.7, y1=mean.IL6.Tr,
         lwd=20,                #thickness of bar
         col="#e36f60",           #color of bar
         lend="butt")          #makes bars blunt ended

#error bar  treated
arrows(x0=1.7, y0=dexam.error.bar.lower, x1=1.7, y1=dexam.error.bar.upper,
       code=3, #draw an arrow on both ends
       angle=90, #have the arrows at 90 degrees from one another
       length=0.03, #length of the arrow feet
       col="#ba4357", #color
       lwd=1.5) #thickness

#error bar untreated
arrows(x0=1.3, y0=nodexam.error.lower, x1=1.3, y1=nodexam.error.upper,
       code=3,
       angle=90,
       length=0.03, 
       col="#e5765f",
       lwd=1.5)

#add a legend
legend(x=2.04, y= 49,
       c("w/o Dexamethasone","Dexamethasone"),
       pch = c(15,15),
       col = c("#f7b257","#e36f60"),
       horiz = FALSE, #on top of each other or next to each other
       cex = 0.6)

dev.off()

#=======================================================================
###########################Bar Graph with IL10#########################
#=======================================================================
#########plot bar graph with segment 
#need to pick the x locations of the bars
x.base <- 1:2
x.labels<- c("Untreated", "Treated")

#Establish highest and lowest y points 
####ADD standard deviation bars using arrow
dexam.error.bar.upper <- (mean.IL10.Tr + se_IL10.Tr)
dexam.error.bar.lower <- (mean.IL10.Tr - se_IL10.Tr)

nodexam.error.upper <- (mean.IL10.Utr + se_IL10.Utr)
nodexam.error.lower <- (mean.IL10.Utr - se_IL10.Utr)

png("Sig.Cytokine.plot.IL10.png", width = 5, height = 6, units = "in", res = 120) ## saves a png in the folder "images" within the working directory
#draw empty box plot
par(mar = c(4, 2, 5, 10))    #margins (bottom, left, top, right)
par(pin = c(2, 4))          #Recangular
par(xpd = NA)             #allows things to be drawn outside of the plot (legens)
plot(NA, xlim = c(1, 2),
     ylim=c(1.21, max(dexam.error.bar.upper, nodexam.error.upper)*1.5), 
     xaxt="n", xlab="", ylab="Concentration (pg/mL)",
     main="IL-10 Cytokine concentration\n by Dexamethasone Treatment")

axis(1, at=c(1.3, 1.7), labels=x.labels, las=2, cex.axis=0.8)

baseline<- 0

#design bars for untreated cytokine means
segments(x0=1.3, y0=baseline, x1=1.3, y1=mean.IL10.Utr,
         lwd=20,
         col="#d085df",
         lend="butt")

#design bar for treated cytokine means
segments(x0=1.7, y0=baseline, x1=1.7, y1=mean.IL10.Tr,
         lwd=20,                #thickness of bar
         col="#cb4f73",           #color of bar
         lend="butt")          #makes bars blunt ended

#error bar  treated
arrows(x0=1.7, y0=dexam.error.bar.lower, x1=1.7, y1=dexam.error.bar.upper,
       code=3, #draw an arrow on both ends
       angle=90, #have the arrows at 90 degrees from one another
       length=0.03, #length of the arrow feet
       col="#b82828", #color
       lwd=1.5) #thickness

#error bar untreated
arrows(x0=1.3, y0=nodexam.error.lower, x1=1.3, y1=nodexam.error.upper,
       code=3,
       angle=90,
       length=0.03, 
       col="#bf5b8e",
       lwd=1.5)
       
#add a legend
legend(x=2.05, y= 3 ,
       c("w/o Dexamethasone","Dexamethasone"),
       pch = c(15,15),
       col = c("#d085df","#cb4f73"),
       horiz = FALSE, #on top of each other or next to each other
       cex = 0.6
       )    #size
dev.off()

#=========================================================================
###########################Bar Graph with IFN-gamma#######################
#=========================================================================
#########plot bar graph with segment 
#need to pick the x locations of the bars
x.base <- 1:2
x.labels<- c("Untreated", "Treated")

#Establish highest and lowest y points 
####ADD standard deviation bars using arrow
dexam.error.bar.upper <- (mean.IFN.Tr + se_IFN.Tr)
dexam.error.bar.lower <- (mean.IFN.Tr - se_IFN.Tr)

nodexam.error.upper <- (mean.IFN.Utr + se_IFN.Utr)
nodexam.error.lower <- (mean.IFN.Utr - se_IFN.Utr)

png("Sig.Cytokine.IFN.plot.png", width = 5, height = 6, units = "in", res = 120) ## saves a png in the folder "images" within the working directory
#draw empty box plot
par(mar = c(4, 2, 5, 10))    #margins (bottom, left, top, right)
par(pin = c(2, 4))          #Recangular
par(xpd = NA)
plot(NA, xlim = c(1, 2),
     ylim=c(2.4, max(dexam.error.bar.upper, nodexam.error.upper)+4), 
     xaxt="n", xlab="", ylab="Concentration (pg/mL)",
     main="ILF-gamma Cytokine concentration\n by Dexamethasone Treatment")

axis(1, at=c(1.3, 1.7), labels=x.labels, las=2, cex.axis=0.8)

baseline<- par("usr")[3]

#design bars for untreated cytokine means
segments(x0=1.3, y0=baseline, x1=1.3, y1=mean.IFN.Utr,
         lwd=20,
         col="#89eed0",
         lend="butt")

#design bar for treated cytokine means
segments(x0=1.7, y0=baseline, x1=1.7, y1=mean.IFN.Tr,
         lwd=20,                #thickness of bar
         col="#a0a4d1",           #color of bar
         lend="butt")          #makes bars blunt ended

#error bar  treated
arrows(x0=1.7, y0=dexam.error.bar.lower, x1=1.7, y1=dexam.error.bar.upper,
       code=3, #draw an arrow on both ends
       angle=90, #have the arrows at 90 degrees from one another
       length=0.03, #length of the arrow feet
       col="#8d77b1", #color
       lwd=1.5) #thickness

#error bar untreated
arrows(x0=1.3, y0=nodexam.error.lower, x1=1.3, y1=nodexam.error.upper,
       code=3,
       angle=90,
       length=0.03, 
       col="#82a2b2",
       lwd=1.5)

#add a legend
legend(x=2.04, y= 3.4,
       c("w/o Dexamethasone","Dexamethasone"),
       pch = c(15,15),
       col = c("#89eed0", "#a0a4d1"),
       horiz = FALSE, #on top of each other or next to each other
       cex = 0.6)
 dev.off()      

#=========================================================================
################FOLD CHANGE Graph#########################################
#=========================================================================
# fold change
# fold change
fold.change <- Cytokine.means$"Dex"[cytokine.cols] / Cytokine.means$"No Dex"[cytokine.cols]

log2.fc <- log2(fold.change)

# error propagation
n.dex   <- sum(Cytokine.data$Treatment == "Dex")
n.nodex <- sum(Cytokine.data$Treatment == "No Dex")

sd.dex   <- Cytokine.sd$"Dex"[cytokine.cols]
sd.nodex <- Cytokine.sd$"No Dex"[cytokine.cols]

sem.dex   <- sd.dex   / sqrt(n.dex)
sem.nodex <- sd.nodex / sqrt(n.nodex)

mean.dex   <- Cytokine.means$"Dex"[cytokine.cols]
mean.nodex <- Cytokine.means$"No Dex"[cytokine.cols]

sem.fc <- fold.change * sqrt( (sem.dex/mean.dex)^2 + (sem.nodex/mean.nodex)^2 )
sem.log2.fc <- sem.fc / (fold.change * log(2))

# x positions
number.cytokines <- length(cytokine.cols)
x.base.f <- 1:number.cytokines

#Group future bars by color
Proinflammatory   <- c("IFN-gamma", "IL-18", "IL-6", "IL-8", "TNF R1", "IP-10", "TREM-1", "MMP-8")
Antiinflammatory  <- c("IL-10")
Endothelial_Injury<- c("Ang-1", "Ang-2", "Thrombomodulin", "VEGF", "Protein C", "ICAM-1", "PAI-1")
Tissue_damage     <- c("RAGE", "SP-D")

cytokine.group <- ifelse(cytokine.cols %in% Proinflammatory, "Proinflammatory",
                         ifelse(cytokine.cols %in% Antiinflammatory, "Antiinflammatory",
                                ifelse(cytokine.cols %in% Endothelial_Injury, "Endothelial_Injury",
                                       "Tissue_damage")))
group.colors <- ifelse(cytokine.group == "Proinflammatory", "#EE799F",
                       ifelse(cytokine.group == "Antiinflammatory", "#FFDEAD",
                              ifelse(cytokine.group == "Endothelial_Injury", "#ffdae0",
                                     "#DDA0DD")))

# Draw the plot 
png("Fold.Change.DexGraph9.png", width = 5, height = 6, units = "in", res = 120) ## saves a png in the folder "images" within the working directory
par(mar = c(8, 5, 4, 5))
plot(NA,
     xlim = c(0.75, number.cytokines + 0.25),
     ylim = c(min(log2.fc - sem.log2.fc),
              max(log2.fc + sem.log2.fc)),
     xaxt = "n",
     xlab = "",
     ylab = "Log2 Fold Change (Dex / No Dex)",
     main = "Log2 Fold Change in Cytokine Concentration\n after Dexamethasone treatment")

usr <- par("usr")
print(usr)
segments(x0 = x.base.f, y0 = usr[3], x1 = x.base.f, y1 = 0,
         col = "lightgrey", lwd = 3,lty = 3, xpd = TRUE)

axis(1, at = x.base.f, labels = cytokine.cols, las = 2, cex.axis = 0.8)

mtext("Immune-system Cytokines", 
      side=1,
      line=6.2,
      at=9,
      cex=1)

# bars
segments(x0 = x.base.f,
         y0 = 0,
         x1 = x.base.f,
         y1 = log2.fc,
         lwd = 15,
         col = group.colors,
         lend = "butt")

abline(h = 0, lwd = 2)

# error bars
upper <- log2.fc + sem.log2.fc
lower <- log2.fc - sem.log2.fc

arrows(x0 = x.base.f,
       y0 = lower,
       x1 = x.base.f,
       y1 = upper,
       code = 3,
       angle = 90,
       length = 0.032,
       col = "#8B475D",
       lwd = 1.5)

# legend
legend("bottomleft",
       legend = c("Vascular Injury markers", "Pro-inflammatory cyto-/chemokines", "Anti-inflammatory cytokines", "Tissue damage markers"),
       fill = c("#ffdae0", "#EE799F", "#FFDEAD","#DDA0DD" ),
       cex = 0.55,
       bty = "o",
       bg ="white")

#ADD ASTERISKS TO THE SIGNIFICANT BARS 
#want the asteriks above IL 10 but below IL6 and IFNy so it has to be sorted
significant.positive<- "IL-10"
significant.negative<- significant[!(significant %in% significant.positive)] 

sig.locations.positive <- which(cytokine.cols %in% significant.positive)
sig.locations.negative <- which(cytokine.cols %in% significant.negative)

#add asteriks above positive
text(x= sig.locations.positive,
     y= upper[sig.locations.positive]+0.1,
     labels = "*",
     cex=1,
     col= "black")

#add asteriks below the negative deflections
text(x = sig.locations.negative,
     y=lower[sig.locations.negative]-0.1,
     labels = "*",
     cex = 1,
     col= "black")

dev.off()

#Assign each cytokine to its biological group
cytokine.group.name <- ifelse(
  cytokine.cols %in% Proinflammatory,
  "Pro-inflammatory cyto-/chemokines",
  ifelse(
    cytokine.cols %in% Antiinflammatory,
    "Anti-inflammatory cytokines",
    ifelse(
      cytokine.cols %in% Endothelial_Injury,
      "Vascular Injury markers",
      "Tissue damage markers"
    )
  )
)

#Create table containing all values needed by Python
fold.change.results <- data.frame(
  Cytokine = cytokine.cols,
  log2_FC = as.numeric(log2.fc),
  SEM_log2_FC = as.numeric(sem.log2.fc),
  Upper = as.numeric(upper),
  Lower = as.numeric(lower),
  Group = cytokine.group.name,
  Significant = cytokine.cols %in% significant
)

#Export results for Python
write.csv(
  fold.change.results,
  "fold_change_results.csv",
  row.names = FALSE
)

####################################################




