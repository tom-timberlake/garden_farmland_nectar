###########################################################################################################################################
######################          GARDENS FILL SEASONAL HUNGER GAPS FOR FARMLAND POLLINATORS           ######################################
###########################################################################################################################################


##########################
## Study background    
#########################

#Our study investigated the potential for gardens to have a widespread role in buffering pollinator populations against seasonal shortages 
#in farmland floral resource supply by answering the following three questions: 

#1) Do gardens fill seasonal gaps in the resource supply of farmland landscapes? 
#2) Do pollinators increase their use of gardens during gaps in farmland resource availability?
#3) Do bumblebees respond more strongly to changes in the timing than the total quantity of garden resources? 

#Using bumblebees as a model group, and Southwest UK as a study region, we test these three hypotheses using a combination of 
#empirical field data and  agent-based-modelling

#The code in this document shows how each prediction was tested and visualised.


##########################
## General housekeeping    
#########################

rm(list=ls()) #clear workspace 

#### Load packages
library(tinytex) 
library(mgcv) 
library(ggplot2) 
library(ggpubr) 
library(scales) 
library(nlme) 
library(MuMIn) 
library(emmeans) 
library(multcomp) 
library(contrast) 
library(dplyr) 
library(plotrix) 
library(broom) 
library(Rmisc) 
library(gridExtra)
library(svglite)
library(ggridges)
library(lme4)
library(marginaleffects)
library(patchwork)

# Set the input and output file paths (Edit these before running)
input.path <- "INPUT YOUR FILE PATH HERE"
output.path <- "INPUT YOUR FILE PATH HERE"


## Running code for each of the 3 questions:

##############################################################################################
### **QUESTION 1: Do gardens fill seasonal gaps in the nectar supply of farmland landscapes?**
##############################################################################################

#### Explanation:

#Generalised additive models (GAM) were used to model the smooth, non-linear trend in garden availability over time for each farmland habitat, and gardens. 
#A thin-plate regression spline was used to model day of the year, with the degree of smoothing selected using the default 
#generalised cross-validation method (Wood 2011). 


#Load in the raw garden nectar values
nectar_data_garden <- read.csv(file.path(input.path,"Garden_nectar_production.csv")) 

#We varied both the link function and the number of knots (from 2-12), to test for the best fitting GAM. The best-fitting model is shown below 
model3 <- gam(sugar_grams_km2 + 0.0001 ~s(day, fx = T, k= 4), ## K=4 was the best model, assessed by comparing AIC scores and goodness of fit to data
              family = Gamma(link = "log"),
              data = nectar_data_garden)

summary(model3) 
gam.check(model3)
AIC(model3)

##Using pred function to predict values for garden nectar for each day of the year
pdat <- expand.grid(day = seq(-28,260,1)) #predicts values for garden nectar for each day of the sampling season (Day 1 = March 01)
pred <- predict (model3, newdata = pdat, na.rm = T, type= "response", se.fit = TRUE)
predframe <- data.frame (pdat,level=0, preds = pred$fit, se = pred$se.fit)

#Save the model predictions
write.csv(predframe, file.path(output.path,"Garden_nectar_predictions.csv"))

#### Plot the nectar phenology of each habitat (Figure 1)

#The following code is used to produce time series kite graphs from abundance data (in this case nectar abundance). 
#The code is the same as used in Russo et al.2013 - Network-based phenological matching which is described in more detail here: 
#https://rpubs.com/thoughtfulbloke/kitegraph

#Kite graphs are produced for each habitat separately (to ensure equivalent scaling) and then joined later to produce a single plot (Figure 2a)

#In addition to the kite plots, there is also code for plotting the proportional contribution of gardens to total farmland nectar supply (Figure 2b) 
#as well as some code for plotting the actual nectar production of each habitat, rather than just relative production (Figure S4).

#Load in habitat phenology data in correct format
habitat_nectar_phenology <- read.csv(file.path(input.path,"habitat_nectar_phenology.csv"))

#Split up dataset by habitat
garden_phenology<-habitat_nectar_phenology[,c(2:3)]
hedge_phenology<-habitat_nectar_phenology[,c(2,11)]
pasture_phenology<-habitat_nectar_phenology[,c(2,7)]
wood_phenology<-habitat_nectar_phenology[,c(2,5)]
margin_phenology<-habitat_nectar_phenology[,c(2,9)]


##############################
##### Plotting woodland data
############################## 

###Setting up to print as an svg
svg(file.path(output.path,"Woodland phenology kite diagram.svg"), bg="transparent", width=10,height=5)
par(mfrow=c(1,1))

graphType <- "raw" #sets if the data is raw numbers or abundance estimates

###Selectively rescale the data
#because it is raw data, it divides the data (all but the first column) by the biggest entry (making everything between 0 and 1) then divides that result by 2 (making it all between 0 and 0.5)
if (graphType == "raw"){
  wood_phenology[,2:ncol(wood_phenology)] <- (wood_phenology[,2:ncol(wood_phenology)] / max(wood_phenology[,2:ncol(wood_phenology)]))/1.9
}

wood_phenology <- wood_phenology[order(wood_phenology[,1]),] # makes sure the data is in the right order - This just makes sure all the rows of the data are in numeric order based on the first column (the distance)

wood_phenology2 <- wood_phenology #make a full copy of the data (Because we need to go there and back, we are making a second copy.)

#reverse the order in the second copy
wood_phenology2 <- wood_phenology2[order(wood_phenology2[,1], decreasing = TRUE),] #Because we need to go back again, we are reversing the order.

#make the values negative
wood_phenology2[,2:ncol(wood_phenology2)] <- wood_phenology2[,2:ncol(wood_phenology2)] * -1 #We want to come back below the step level at the equivalent height so make our second set of data negative

#stick the two sets of data together to make one
wood_phenology <- rbind(wood_phenology,wood_phenology2)

###Make a graph
#Because we are making a blank graph and building it up with the kites, we will start with working out the dimensions
leftedge <- min(wood_phenology[,1])
rightedge <- max(wood_phenology[,1])
bottomedge <- 0
topedge <- ncol(wood_phenology)

#we will also store the old woods and set new ones as we are tinkering with the amount of space on the left
oldwoods <- par()$mar
par(mar=c(2,8,12,0.5))

#make the initial graph with no contents or y axis (this will appear in the plots tab in the lower right of RStudio)
plot(c(leftedge,rightedge), c(bottomedge,topedge), type= "n", xlab="Date", frame.plot=F, yaxt="n", ylab="",xaxt="n")
axis(2, labels=names(wood_phenology)[2:ncol(wood_phenology)], at=1:(ncol(wood_phenology)-1), las=2, lty=0, font=3) #cex.axis=1

#Now loop through each column, making the kite plot for each by drawing a connect the points polygon
xValues = wood_phenology[,1]
for (i in 2:ncol(wood_phenology)){
  yValues = i + wood_phenology[,i] - 1
  polygon(xValues,yValues, col=rgb(0.36,0.44,0.20))
}

#Finally, because we played with the wood setting, we will put them back to the defaults
par(mar=oldwoods)

#Add vertical dotted lines for each month
abline(v=c(1,32,62,93,123,154,185,215,246), col="lightgray", lty="dotted")

dev.off() #End plot print

###########################
### Plotting hedgerow data
###########################

svg(file.path(output.path,"Hedgerow phenology kite diagram.svg"), bg="transparent", width=10,height=5)
par(mfrow=c(1,1))
graphType <- "raw"

if (graphType == "raw"){
  hedge_phenology[,2:ncol(hedge_phenology)] <- (hedge_phenology[,2:ncol(hedge_phenology)] / max(hedge_phenology[,2:ncol(hedge_phenology)]))/1.9
}

###Setting up the data
hedge_phenology <- hedge_phenology[order(hedge_phenology[,1]),] 
hedge_phenology2 <- hedge_phenology
hedge_phenology2 <- hedge_phenology2[order(hedge_phenology2[,1], decreasing = TRUE),]
hedge_phenology2[,2:ncol(hedge_phenology2)] <- hedge_phenology2[,2:ncol(hedge_phenology2)] * -1 
hedge_phenology <- rbind(hedge_phenology,hedge_phenology2)

###Make a graph
leftedge <- min(hedge_phenology[,1])
rightedge <- max(hedge_phenology[,1])
bottomedge <- 0
topedge <- ncol(hedge_phenology)

oldMargins <- par()$mar
par(mar=c(2,8,0,0.5))

#make the initial graph with no contents or y axis (this will appear in the plots tab in the lower right of RStudio)
plot(c(leftedge,rightedge), c(bottomedge,topedge), type= "n", xlab="Date", frame.plot=F, yaxt="n", ylab="",xaxt="n")
axis(2, labels=names(hedge_phenology)[2:ncol(hedge_phenology)], at=1:(ncol(hedge_phenology)-1), las=2, lty=0, font=3) #cex.axis=1
xValues = hedge_phenology[,1]
for (i in 2:ncol(hedge_phenology)){
  yValues = i + hedge_phenology[,i] - 1
  polygon(xValues,yValues, col=rgb(0.56,0.22,0.21))
}

par(mar=oldMargins)
abline(v=c(1,32,62,93,123,154,185,215,246), col="lightgray", lty="dotted")
dev.off()


################################### 
##### Plotting Field margin data
###################################

svg(file.path(output.path,"Margin phenology kite diagram.svg"), bg="transparent", width=10,height=5)
par(mfrow=c(1,1))
graphType <- "raw" 

###Selectively rescale the data
if (graphType == "raw"){
  margin_phenology[,2:ncol(margin_phenology)] <- (margin_phenology[,2:ncol(margin_phenology)] / max(margin_phenology[,2:ncol(margin_phenology)]))/1.9
}

###Setting up the data
margin_phenology <- margin_phenology[order(margin_phenology[,1]),] 
margin_phenology2 <- margin_phenology
margin_phenology2 <- margin_phenology2[order(margin_phenology2[,1], decreasing = TRUE),]
margin_phenology2[,2:ncol(margin_phenology2)] <- margin_phenology2[,2:ncol(margin_phenology2)] * -1 
margin_phenology <- rbind(margin_phenology,margin_phenology2)

###Make a graph
leftedge <- min(margin_phenology[,1])
rightedge <- max(margin_phenology[,1])
bottomedge <- 0
topedge <- ncol(margin_phenology)
oldMargins <- par()$mar
par(mar=c(2,8,12,0.5))

plot(c(leftedge,rightedge), c(bottomedge,topedge), type= "n", xlab="Date", frame.plot=F, yaxt="n", ylab="",xaxt="n")
axis(2, labels=names(margin_phenology)[2:ncol(margin_phenology)], at=1:(ncol(margin_phenology)-1), las=2, lty=0, font=3) #cex.axis=1

xValues = margin_phenology[,1]
for (i in 2:ncol(margin_phenology)){
  yValues = i + margin_phenology[,i] - 1
  polygon(xValues,yValues, col=rgb(0.54,0.32,0.14))
}

par(mar=oldMargins)
abline(v=c(1,32,62,93,123,154,185,215,246), col="lightgray", lty="dotted")
dev.off()


##############################
##### Plotting Pasture data
##############################

svg(file.path(output.path,"Pasture phenology kite diagram.svg"), bg="transparent", width=10,height=5)
par(mfrow=c(1,1))
graphType <- "raw" #sets if the data is raw numbers or abundance estimates

###Selectively rescale the data
if (graphType == "raw"){
  pasture_phenology[,2:ncol(pasture_phenology)] <- (pasture_phenology[,2:ncol(pasture_phenology)] / max(pasture_phenology[,2:ncol(pasture_phenology)]))/1.9
}

###Setting up the data
pasture_phenology <- pasture_phenology[order(pasture_phenology[,1]),] 
pasture_phenology2 <- pasture_phenology 
pasture_phenology2 <- pasture_phenology2[order(pasture_phenology2[,1], decreasing = TRUE),] 
pasture_phenology2[,2:ncol(pasture_phenology2)] <- pasture_phenology2[,2:ncol(pasture_phenology2)] * -1 
pasture_phenology <- rbind(pasture_phenology,pasture_phenology2)

###Make a graph
leftedge <- min(pasture_phenology[,1])
rightedge <- max(pasture_phenology[,1])
bottomedge <- 0
topedge <- ncol(pasture_phenology)

oldMargins <- par()$mar
par(mar=c(2,8,0,0.5))

plot(c(leftedge,rightedge), c(bottomedge,topedge), type= "n", xlab="Date", frame.plot=F, yaxt="n", ylab="",xaxt="n")
axis(2, labels=names(pasture_phenology)[2:ncol(pasture_phenology)], at=1:(ncol(pasture_phenology)-1), las=2, lty=0, font=3) #cex.axis=1

xValues = pasture_phenology[,1]
for (i in 2:ncol(pasture_phenology)){
  yValues = i + pasture_phenology[,i] - 1
  polygon(xValues,yValues, col=rgb(0.19,0.32,0.48))
}

par(mar=oldMargins)
abline(v=c(1,32,62,93,123,154,185,215,246), col="lightgray", lty="dotted")
dev.off()

##############################
##### Plotting Garden data
##############################

svg(file.path(output.path,"Garden phenology kite diagram.svg"), bg="transparent", width=10,height=5)
par(mfrow=c(1,1))
graphType <- "raw" 

###Selectively rescale the data
if (graphType == "raw"){
  garden_phenology[,2:ncol(garden_phenology)] <- (garden_phenology[,2:ncol(garden_phenology)] / max(garden_phenology[,2:ncol(garden_phenology)]))/1.9
}

###Setting up the data
garden_phenology <- garden_phenology[order(garden_phenology[,1]),] 
garden_phenology2 <- garden_phenology 
garden_phenology2 <- garden_phenology2[order(garden_phenology2[,1], decreasing = TRUE),] 
garden_phenology2[,2:ncol(garden_phenology2)] <- garden_phenology2[,2:ncol(garden_phenology2)] * -1 
garden_phenology <- rbind(garden_phenology,garden_phenology2)

###Make a graph
leftedge <- min(garden_phenology[,1])
rightedge <- max(garden_phenology[,1])
bottomedge <- 0
topedge <- ncol(garden_phenology)
oldMargins <- par()$mar
par(mar=c(2,8,0,0.5))

plot(c(leftedge,rightedge), c(bottomedge,topedge), type= "n", xlab="Date", frame.plot=F, yaxt="n", ylab="",xaxt="n")
axis(2, labels=names(garden_phenology)[2:ncol(garden_phenology)], at=1:(ncol(garden_phenology)-1), las=2, lty=0, font=3) #cex.axis=1

xValues = garden_phenology[,1]
for (i in 2:ncol(garden_phenology)){
  yValues = i + garden_phenology[,i] - 1
  polygon(xValues,yValues, col=rgb(0.19,0.32,0.48))
}

par(mar=oldMargins)
abline(v=c(1,32,62,93,123,154,185,215,246), col="lightgray", lty="dotted")
dev.off()

###################################################################################
#Plotting proportional contribution of gardens to total landscape nectar (Fig. 1b)
###################################################################################

#Load in data
farm_garden_nectar <- read.csv(file.path(input.path,"Farm_garden_nectar_production.csv")) 
#Subset to only the garden nectar data
garden_data <- farm_garden_nectar[farm_garden_nectar$habitat == "garden",] 

#Plot proportional contribution of garden nectar to total landscape nectar (=daily garden nectar/SUM(garden+farmland nectar))
Garden_prop_nectar <- ggplot(garden_data, aes(x=day, y=prop_sugar)) + 
  geom_line(size=1.5) +
  geom_line(aes(y = prop_sugar-prop_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = prop_sugar+prop_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = 0.5), color = "red", linetype = "dashed", size=1) + #Adds a line to show the point at which garden nectar = farmland nectar
  theme_classic(base_size = 22) + 
  labs(x=bquote(""),
       y=bquote("Garden nectar contribution")) + 
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ylim(0,1)+
  geom_text(x=100, y=0.55, size=6, label="Equal garden & farmland nectar")

#Set up and print the plot with a manually added x axis
svg(file.path(output.path,"Garden_prop_nectar_km2.svg"), bg="transparent", width=11,height=5)
par(mar=c(5,8,2,1)) #sets the size of the margins around the plot to be slightly larger than the default settings
Garden_prop_nectar + scale_x_continuous(name="Date", breaks=c(1,32,62,93,123,154,185,215,246), 
                                        labels=c( "March","April","May","June","July","Aug","Sept","Oct","Nov"))
dev.off()


###################################################################################
#Plotting actual (rather than relative) values  per m2 for each habitat (Figure S4)
###################################################################################

pasture_plot <- ggplot(habitat_nectar_phenology, aes(x=day, y=pasture_grams_m2)) + 
  geom_line(size=1.5) +
  geom_line(aes(y = pasture_grams_m2-pasture_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = pasture_grams_m2+pasture_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = 0), color = "red", linetype = "dotted", size=1) +
  theme_classic(base_size = 22) + 
  labs(x=bquote(""),
       y=bquote("Sugar/m"^2*"/day (grams)")) + #sets plot theme - there are various other themes available
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ggtitle("Pasture")+
  ylim(0,0.0065)+
  scale_x_continuous(name="Date", breaks=c(1,32,62,93,123,154,185,215,246), labels=c( "March","April","May","June","July","Aug","Sept","Oct","Nov"))

hedgerow_plot <- ggplot(habitat_nectar_phenology, aes(x=day, y=hedgerow_grams_m2)) + 
  geom_line(size=1.5) +
  geom_line(aes(y = hedgerow_grams_m2-hedgerow_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = hedgerow_grams_m2+hedgerow_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = 0), color = "red", linetype = "dotted", size=1) +
  theme_classic(base_size = 22) + 
  labs(x=bquote(""),
       y=bquote("Sugar/m"^2*"/day (grams)")) + #sets plot theme - there are various other themes available
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ggtitle("Hedgerow")+
  ylim(0,0.065)+
  scale_x_continuous(name="Date", breaks=c(1,32,62,93,123,154,185,215,246), labels=c( "March","April","May","June","July","Aug","Sept","Oct","Nov"))

margin_plot <- ggplot(habitat_nectar_phenology, aes(x=day, y=margin_grams_m2)) + 
  geom_line(size=1.5) +
  geom_line(aes(y = margin_grams_m2-margin_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = margin_grams_m2+margin_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = 0), color = "red", linetype = "dotted", size=1) +
  theme_classic(base_size = 22) + 
  labs(x=bquote(""),
       y=bquote("Sugar/m"^2*"/day (grams)")) + #sets plot theme - there are various other themes available
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ggtitle("Margin")+
  ylim(0,0.065)+
  scale_x_continuous(name="Date", breaks=c(1,32,62,93,123,154,185,215,246), labels=c( "March","April","May","June","July","Aug","Sept","Oct","Nov"))

woodland_plot <- ggplot(habitat_nectar_phenology, aes(x=day, y=woodland_grams_m2)) + 
  geom_line(size=1.5) +
  geom_line(aes(y = woodland_grams_m2-woodland_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = woodland_grams_m2+woodland_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = 0), color = "red", linetype = "dotted", size=1) +
  theme_classic(base_size = 22) + 
  labs(x=bquote(""),
       y=bquote("Sugar/m"^2*"/day (grams)")) + #sets plot theme - there are various other themes available
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ggtitle("Woodland")+
  ylim(0,0.065)+
  scale_x_continuous(name="Date", breaks=c(1,32,62,93,123,154,185,215,246), labels=c( "March","April","May","June","July","Aug","Sept","Oct","Nov"))

garden_plot <- ggplot(habitat_nectar_phenology, aes(x=day, y=garden_grams_m2)) + 
  geom_line(size=1.5) +
  geom_line(aes(y = garden_grams_m2-garden_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = garden_grams_m2+garden_grams_m2_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = 0), color = "red", linetype = "dotted", size=1) +
  theme_classic(base_size = 22) + 
  labs(x=bquote(""),
       y=bquote("Sugar/m"^2*"/day (grams)")) + #sets plot theme - there are various other themes available
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  ggtitle("Garden")+
  ylim(0,0.065)+
  scale_x_continuous(name="Date", breaks=c(1,32,62,93,123,154,185,215,246), labels=c( "March","April","May","June","July","Aug","Sept","Oct","Nov"))


habitat_nectar_plots <- ggarrange(hedgerow_plot, margin_plot, woodland_plot, garden_plot,pasture_plot, nrow = 3,ncol=2, labels="auto",  common.legend = TRUE, legend="right")

ggsave(plot=habitat_nectar_plots, filename=file.path(output.path,"Habitat_nectar_m2.svg"), width=13, height=12, dpi=500)
ggsave(plot=habitat_nectar_plots, filename=file.path(output.path,"Habitat_nectar_m2.png"), width=13, height=12, dpi=500)


#### Identifying hunger gaps for bumblebees (Figure 2)

#This code plots the supply of farmland nectar per bee (i.e. nectar availability divided by bumblebee abundance at each timepoint) 
#with and without the contribution of gardens. Superimposed on the graph is the total daily estimated nectar requirements of each bee, 
#taken from Rotheray et al. 2017.

#Load in supply-demand data
supply_demand_data <- read.csv(file.path(input.path,"Farm_garden_nectar_supply_demand.csv")) 

# Plotting nectar supply vs demand + Std Err
supply_demand_plot<- ggplot(supply_demand_data, aes(x=day, y=total_sugar_per_bee, color="coral1")) + 
  geom_line(size=1.5) +
  geom_point( size=3)+
  theme_classic(base_size = 22) +  #sets plot theme - there are various other themes available
  geom_line(aes(y = total_sugar_per_bee-total_sugar_per_bee_se), color = "coral1", linetype = "dashed", size=1) +
  geom_line(aes(y = total_sugar_per_bee+total_sugar_per_bee_se), color = "coral1", linetype = "dashed", size=1) +
  geom_line(aes(y = bee_demand), color = "black",  size=1.5) +
  geom_point(aes(y = bee_demand), color = "black",  size=3)+
  geom_line(aes(y = bee_demand-bee_demand_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = bee_demand+bee_demand_se), color = "black", linetype = "dashed", size=1) +
  geom_line(aes(y = farm_sugar_per_bee), color = "cyan3",  size=1.5) +
  geom_point(aes(y = farm_sugar_per_bee), color = "cyan3",  size=3) +
  geom_line(aes(y = farm_sugar_per_bee-farm_sugar_per_bee_se), color = "cyan3",linetype = "dashed",  size=1) +
  geom_line(aes(y = farm_sugar_per_bee+farm_sugar_per_bee_se), color = "cyan3",linetype = "dashed",  size=1) +
  labs(x=bquote("Month"),
       y=bquote("Daily sugar per bee/km" ^2* "(grams)")) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        legend.position = "none")

#Print supply-demand plot with manual x-axis
svg(file.path(output.path,"Farm_garden_Supply_demand_plusSE.svg"), bg="transparent", width=8,height=4.5)
print(supply_demand_plot) + 
  scale_y_continuous(trans='log10',labels = comma, limits=c(0.01, 100))+
  scale_x_continuous(name="Date", breaks= c(1,32,62,93,123,154,185,215,246), labels= c("March","April","May","June","July","Aug","Sept","Oct","Nov"))
dev.off()


####################################################################################################################
### **Question 2: do pollinators increase their use of gardens during gaps in farmland resource availability?** ###
####################################################################################################################

#Desription
#Here we use BeeWalk data from transects in South-West UK to plot the pattern of bumblebee abundance in gardens vs farmland through the year

#Load data
beewalk_data <- read.csv(file.path(input.path,"BeeWalk data_South-West.csv")) 

##Split dataset up by species
BeeWalk_B.terr <- beewalk_data[beewalk_data$species == "B.terrestris",] 
BeeWalk_B.lap <- beewalk_data[beewalk_data$species == "B.lapidarius",] 
BeeWalk_B.pasc <- beewalk_data[beewalk_data$species == "B.pascuorum",] 
BeeWalk_B.prat <- beewalk_data[beewalk_data$species == "B.pratorum",] 

#### Plot each species' seasonal activity pattern as a 4-panel plot
svg(file.path(output.path,"Farm-garden_abundance_ratio_scaled.svg"), bg="transparent", width=13,height=9) #Start svg print

par(mfrow=c(2,2))

## Plotting B.terrestris
plot(gf_ratio_w_q_scaled ~ month, data = BeeWalk_B.terr,
     type="o", lty=1, lwd=4, font.lab=2, cex.axis=1.5,
     ylab = "Garden:Farm abundance ratio", 
     xlab = "B. terrestris", # labels x axis
     cex.lab = 1.5,  las = 1, col = "black", pch = 20,
     xaxt = "n",
     xlim=c(3,10), ylim=c(0,12),
     bty = "L") 
abline(h=1, col="darkgray", lty="dashed", lwd=3) #Add a line to show equal garden and farmland abundance

#Adding labels for months
axis(side = 1, at = c(3,4,5,6,7,8,9,10), 
     labels = c("March","April","May","June","July","Aug","Sept","Oct"),cex.axis=1.5 )

## Plotting B.lapidarius
plot(gf_ratio_w_q_scaled ~ month, data = BeeWalk_B.lap,
     type="o", lty=1, lwd=4, font.lab=2, cex.axis=1.5,
     ylab = "Garden:Farm abundance ratio", 
     xlab = "B. lapidarius", # labels x axis
     cex.lab = 1.5, las = 1, col = "black", pch = 20, 
     xaxt = "n", 
     xlim=c(3,10), ylim=c(0,12),
     bty = "L")
abline(h=1, col="darkgray", lty="dashed", lwd=3) #Add a line to show equal garden and farmland abundance

#Adding labels for months
axis(side = 1, at = c(3,4,5,6,7,8,9,10), 
     labels = c("March","April","May","June","July","Aug","Sept","Oct"),cex.axis=1.5 )


##Plotting B.pascuorum
plot(gf_ratio_w_q_scaled ~ month, data = BeeWalk_B.pasc,
     type="o", lty=1, lwd=4, font.lab=2, cex.axis=1.5,
     ylab = "Garden:Farm abundance ratio", 
     xlab = "B. pascuorum", # labels x axis
     cex.lab = 1.5,  las = 1, col = "black",  pch = 20,
     xaxt = "n", # suppresses the x axis labels
     xlim=c(3,10), ylim=c(0,12), 
     bty = "L") # changes the box type to an L rather than square
abline(h=1, col="darkgray", lty="dashed", lwd=3) #Add a line to show equal garden and farmland abundance

#Adding labels for months
axis(side = 1, at = c(3,4,5,6,7,8,9,10), 
     labels = c("March","April","May","June","July","Aug","Sept","Oct"),cex.axis=1.5 )

##Plotting B.pratorum
plot(gf_ratio_w_q_scaled ~ month,  data = BeeWalk_B.prat,
     type="o", lty=1, lwd=4,  font.lab=2, cex.axis=1.5,
     ylab = "Garden:Farm abundance ratio", 
     xlab = "B. pratorum", # labels x axis
     cex.lab = 1.5, las = 1, col = "black", pch = 20, 
     xaxt = "n", # suppresses the x axis labels
     xlim=c(3,10), ylim=c(0,12),
     bty = "L")
abline(h=1, col="darkgray", lty="dashed", lwd=3) #Add a line to show equal garden and farmland abundance

#Adding labels for months
axis(side = 1, at = c(3,4,5,6,7,8,9,10), 
     labels = c("March","April","May","June","July","Aug","Sept","Oct"),cex.axis=1.5 )

dev.off() #End svg print


#######################################################################################################################
### **Question 3: are farmland pollinator populations impacted by the removal or alteration of garden nectar supply?**#
#######################################################################################################################

#Description:
#Here we analyse the outputs of an in-silico experiment using the agent-based model ‘BEE-STEWARD’ (Twiston‐Davies et al. 2021) 
#implemented in Netlogo. The model simulates the foraging behaviour, life history and colony growth of Bombus terrestris foraging for nectar 
#and pollen in a spatially explicit, user-defined landscape. Our model landscapes were based on 12 real circular farmland landscapes 
#with a 1 km radius (3.14 km2 area) in Southwest UK which were characterised and mapped by Anonymous  et al. (2021a)

####Loading in data from BEE-STEWARD model simulations implemented in Netlogo 

experiment_outputs <- read.csv(file.path(input.path,"ALL_BEE-STEWARD_RESULTS.csv")) 

#define the different variables
experiment_outputs$farm=as.factor(experiment_outputs$farm)
experiment_outputs$treatment=as.factor(experiment_outputs$treatment)
experiment_outputs$endpoint=as.factor(experiment_outputs$endpoint)
experiment_outputs$month=as.factor(experiment_outputs$month)

##Subset the data by the different variables measured in the experiment
colony_density <- experiment_outputs[experiment_outputs$endpoint == "ColonyDensity_km2",] 
colony_number <- experiment_outputs[experiment_outputs$endpoint == "TotalColonies",] 
adult_queens <- experiment_outputs[experiment_outputs$endpoint == "TotalAdultQueens",] 
adult_workers <- experiment_outputs[experiment_outputs$endpoint == "TotalAdultWorkers",]
total_males <- experiment_outputs[experiment_outputs$endpoint == "TotalMales",] 
total_bees <- experiment_outputs[experiment_outputs$endpoint == "TotalBeesEverProduced",] 
total_colonies <- experiment_outputs[experiment_outputs$endpoint == "TotalColoniesEverProduced",] 
total_active_bees <- experiment_outputs[experiment_outputs$endpoint == "TotalActiveBees",] 
new_queens <- experiment_outputs[experiment_outputs$endpoint == "new_queens",] 
queen_survival <- experiment_outputs[experiment_outputs$endpoint == "queen_survival",] 
annual_cd <- total_active_bees[total_active_bees$month == "March",]

#Remove irrelevant part of the data (e.g. month by month values for annual metrics like queen survival)

colony_density_months <- subset(colony_density, !(month %in% c('Overall_max', 'Overall_mean'))) #Create dataset of month by month colony density values
colony_density_max <- subset(colony_density, (month %in% 'Overall_max')) #Create dataset of maximum colony density values
colony_number_max <- subset(colony_number, (month %in% 'Overall_max')) #Create dataset of maximum colony number values
new_queens <- subset(new_queens, (month %in% 'Overall_mean')) #Remove month by month values (which are all NA), leaving only the mean values 
queen_survival <- subset(queen_survival, (month %in% 'Overall_mean')) #Remove month by month values (which are all NA), leaving only the mean values 
total_bees_max <- subset(total_bees, (month %in% 'Overall_max')) #Remove month by month values, leaving only the max values 


####################################################################
##### Plotting colony density data by treatment and month ########
####################################################################

##Summarising data at the farm-level
cd_month_farm <- colony_density_months %>% dplyr::group_by(farm, treatment, month, month_number) %>%
  dplyr::summarize(farm_mean = mean(value)) # calculate the mean for each farm

##Summarising data
cd_month_summary <- cd_month_farm %>% dplyr::group_by(treatment, month) %>%
  dplyr::summarize(mean_months = mean(farm_mean), # calculate the mean for each farm
                   sd_months = sd(farm_mean)) # calculate the SD for each farm

##Reordering factors by name
cd_month_summary$month <- factor(cd_month_summary$month, levels = c("March","April","May","June","July","Aug","Sept","Oct"))
cd_month_summary$treatment <- factor(cd_month_summary$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

##Write summary to csv
#write.csv(cd_month_summary , file="BEESTEWARD_colony density_monthly summary.csv")

##Plotting the seasonality of colony density
CD<- ggplot(cd_month_summary, aes(x=month, y=mean_months, group=treatment, color=treatment)) + 
  geom_line(size=1) +
  geom_point( size=2, position=position_dodge(0.2))+
  theme_classic() +  #sets plot theme - there are various other themes available
  labs(x=bquote("Month"),
       y=bquote("Mean colony density (nests/km"^2*")")) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))+
  geom_errorbar(aes(ymin=mean_months-sd_months, ymax=mean_months+sd_months), size=1, width=0,
                position=position_dodge(0.2))

#Add correct colour scales etc
seasonal_colony_density <- CD + scale_colour_manual(values = c("#ee7621ff", "#f5dc83ff",  "#8fa33fff",  "#cdd4dcff"),
                               name="Treatment",
                               labels=c("1) Garden", "2) Garden with pasture phenology", "3) Pasture replacement", "4) No garden"))+
  theme(legend.position = c(0.8, 0.8))


####################################################################
#####              Plotting maximum colony numbers          ########
####################################################################

# Define custom colors
colours <- c("#ee7621ff", "#f5dc83ff",  "#8fa33fff",  "#cdd4dcff")

##Summarising data at the farm-level
colony_number_farm <- colony_number_max %>% dplyr::group_by(farm, treatment) %>%
  dplyr::summarize(farm_mean = mean(value)) # calculate the mean for each farm

##Summarising data at the treatment-level
colony_number_treatment <- colony_number_farm %>% dplyr::group_by(treatment) %>%
  dplyr::summarize(overall_mean = mean(farm_mean),
                   overall_sd = sd(farm_mean)) # calculate the mean and SD for each treatment

colony_number_treatment$metric <- "max_colony_number"

colony_number_farm$treatment <- factor(colony_number_farm$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

max_colonies <- ggplot(colony_number_farm, aes(x = treatment, y = farm_mean, color = treatment)) +
  geom_point(size=2, position=position_jitter(width=0.1, height=0)) +  # Plot farm_mean values for each farm
  geom_point(data = colony_number_treatment, aes(x = treatment, y = overall_mean), shape = 9, size = 2, color = "black") +  # Overlay treatment means
  theme_classic() + ylab("Maximum colony number") + xlab("") + #This sets the plot theme (black and white) and labels the x and y axis manually so you can change the labels if you want
  scale_color_manual(values = colours) +  # Specify custom colors
  labs(title = "",
       x = "Treatment") +
  scale_x_discrete(labels = c("1", "2", "3", "4")) +   # Set x-axis labels to numbers
  theme(legend.position = "none")



####################################################################
#####              Plotting new queens produced          ########
####################################################################

##Summarising data at the farm-level
new_queens_farm <- new_queens %>% dplyr::group_by(farm, treatment) %>%
  dplyr::summarize(farm_mean = mean(value)) # calculate the mean for each farm

##Summarising data at the treatment-level
new_queens_treatment <- new_queens_farm %>% dplyr::group_by(treatment) %>%
  dplyr::summarize(overall_mean = mean(farm_mean),
                   overall_sd = sd(farm_mean)) # calculate the mean and SD for each treatment

new_queens_treatment$metric <- "new_queens"

new_queens_farm$treatment <- factor(new_queens_farm$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

new_queens <- ggplot(new_queens_farm, aes(x = treatment, y = farm_mean, color = treatment)) +
  geom_point(size=2, position=position_jitter(width=0.1, height=0)) +  # Plot farm_mean values for each farm
  geom_point(data = new_queens_treatment, aes(x = treatment, y = overall_mean), shape = 9, size = 2, color = "black") +  # Overlay treatment means
  theme_classic() + ylab("New queens produced") + xlab("") + #This sets the plot theme (black and white) and labels the x and y axis manually so you can change the labels if you want
  scale_color_manual(values = colours) +  # Specify custom colors
  labs(title = "",
       x = "Treatment") +
  scale_x_discrete(labels = c("1", "2", "3", "4")) +   # Set x-axis labels to numbers
  theme(legend.position = "none")


####################################################################
#####              Plotting spring queen survival          ########
####################################################################

##Summarising data at the farm-level
queen_survival_farm <- queen_survival %>% dplyr::group_by(farm, treatment) %>%
  dplyr::summarize(farm_mean = mean(value)) # calculate the mean for each farm

##Summarising data at the treatment-level
queen_survival_treatment <- queen_survival_farm %>% dplyr::group_by(treatment) %>%
  dplyr::summarize(overall_mean = mean(farm_mean),
                   overall_sd = sd(farm_mean)) # calculate the mean and SD for each treatment

queen_survival_treatment$metric <- "queen_survival"

queen_survival_farm$treatment <- factor(queen_survival_farm$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

queen_survival <- ggplot(queen_survival_farm, aes(x = treatment, y = farm_mean, color = treatment)) +
  geom_point(size=2, position=position_jitter(width=0.1, height=0)) +  # Plot farm_mean values for each farm
  geom_point(data = queen_survival_treatment, aes(x = treatment, y = overall_mean), shape = 9, size = 2, color = "black") +  # Overlay treatment means
  theme_classic() + ylab("Spring queen survival") + xlab("") + #This sets the plot theme (black and white) and labels the x and y axis manually so you can change the labels if you want
  scale_color_manual(values = colours) +  # Specify custom colors
  labs(title = "",
       x = "Treatment") +
  scale_x_discrete(labels = c("1", "2", "3", "4")) +   # Set x-axis labels to numbers
  scale_y_continuous(labels = scales::percent) + 
  theme(legend.position = "none")


####################################################################
#####              Total  number of active   bees         ########
####################################################################

##Summarising data at the farm-level
total_bees_farm <- total_bees_max %>% dplyr::group_by(farm, treatment) %>%
  dplyr::summarize(farm_mean = mean(value)) # calculate the mean for each farm

##Summarising data at the treatment-level
total_bees_treatment <- total_bees_farm %>% dplyr::group_by(treatment) %>%
  dplyr::summarize(overall_mean = mean(farm_mean),
                   overall_sd = sd(farm_mean)) # calculate the mean and SD for each treatment

total_bees_treatment$metric <- "total_active_bees"

total_bees_farm$treatment <- factor(total_bees_farm$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

total_bees <- ggplot(total_bees_farm, aes(x = treatment, y = farm_mean, color = treatment)) +
  geom_point(size=2, position=position_jitter(width=0.1, height=0)) +  # Plot farm_mean values for each farm
  geom_point(data = total_bees_treatment, aes(x = treatment, y = overall_mean), shape = 9, size = 2, color = "black") +  # Overlay treatment means
  theme_classic() + ylab("Total active bees through season") + xlab("") + #This sets the plot theme (black and white) and labels the x and y axis manually so you can change the labels if you want
  scale_color_manual(values = colours) +  # Specify custom colors
  labs(title = "",
       x = "Treatment") +
  scale_x_discrete(labels = c("1", "2", "3", "4")) +   # Set x-axis labels to numbers
  theme(legend.position = "none")



#PRODUCE ONE PANEL OF ALL PLOTS

# Combine the plots using patchwork
combined_plot <- seasonal_colony_density / (max_colonies | total_bees) / (new_queens | queen_survival) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 14, face = "bold"))

#Save combined plot as a PNG and SVG
ggsave(plot=combined_plot, filename=file.path(output.path,"BEESTEWARD_model_outputs.png"), width=6, height=8, dpi=500, bg="white")
ggsave(plot=combined_plot, filename=file.path(output.path,"BEESTEWARD_model_outputs.svg"), width=6, height=8, dpi=500, bg="white")
 
#Bind together summary results and output as csv file

all_results <- rbind(colony_number_treatment, new_queens_treatment, queen_survival_treatment, total_bees_treatment)
all_results <- select(all_results, metric, treatment, overall_mean, overall_sd)

write.csv(all_results, file.path(output.path,"BEESTEWARD_model_result_summary.csv"))


########################################################################
####  Statistical analysis of BEE-STEWARD model results
########################################################################

###################################
##### Colony density analysis
###################################

##Reordering factors by name
colony_density_max$treatment <- factor(colony_density_max$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

GLMM_CD <- lme(value ~ treatment, random = ~1|farm, method="REML", data = colony_density_max)
summary(GLMM_CD)

#Check model residuals
plot(GLMM_CD)

#Post-hoc test
summary(glht(GLMM_CD, emm(pairwise ~ treatment)), test=adjusted(type="bonferroni"))

# Generate pairwise comparisons and plot (Figure S7)
marginal_means_CD <- emmeans(GLMM_CD, specs = "treatment")
pairwise_comparisons_CD <- pairs(marginal_means_CD)
pairwise_comparisons_CD_plot <- plot(pairwise_comparisons_CD) +
  geom_vline(xintercept = 0, linetype = "dashed") + # Add dashed vertical line at x=0
  #xlim(-3, 3)+ # Specify x axis limits (adjust these limits as needed)
  labs(title = "Colony density",
    x = "Estimated Effect Size",
    y = "Pairwise contrast") # Add title and axis labels
  
################################### 
##### Max colony number analysis
###################################

##GLMM of max yearly colony numbers
#Reordering factors by name
colony_number_max$treatment <- factor(colony_number_max$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

GLMM_colony_number <- lme(value ~ treatment, random = ~1|farm, method="REML", data = colony_number_max)
summary(GLMM_colony_number)

#Check model residuals
plot(GLMM_colony_number)

#Post-hoc test
summary(glht(GLMM_colony_number, emm(pairwise ~ treatment)), test=adjusted(type="bonferroni"))

# Generate pairwise comparisons and plot
marginal_means_colony_number <- emmeans(GLMM_colony_number, specs = "treatment")
pairwise_comparisons_colony_number <- pairs(marginal_means_colony_number)
pairwise_comparisons_colony_number_plot <- plot(pairwise_comparisons_colony_number) +
  geom_vline(xintercept = 0, linetype = "dashed") + # Add dashed vertical line at x=0
  #xlim(-4, 4)+ # Specify x axis limits (adjust these limits as needed)
  labs(title = "Maximum colony number",
       x = "Estimated Effect Size",
       y = "Pairwise contrast") # Add title and axis labels 

################################### 
##### New queen analysis
###################################

##GLMM of new queen numbers
#Reordering factors by name
new_queens$treatment <- factor(new_queens$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

GLMM_new_queens <- lme(value ~ treatment, random = ~1|farm, method="REML", data = new_queens)
summary(GLMM_new_queens)

#Check model residuals
plot(GLMM_new_queens)

#Post-hoc test
summary(glht(GLMM_new_queens, emm(pairwise ~ treatment)), test=adjusted(type="bonferroni"))

# Generate pairwise comparisons and plot
marginal_means_new_queens <- emmeans(GLMM_new_queens, specs = "treatment")
pairwise_comparisons_new_queens <- pairs(marginal_means_new_queens)
pairwise_comparisons_new_queens_plot <- plot(pairwise_comparisons_new_queens) +
  geom_vline(xintercept = 0, linetype = "dashed") + # Add dashed vertical line at x=0
  #xlim(-4, 4)+ # Specify x axis limits (adjust these limits as needed)
  labs(title = "New queens",
       x = "Estimated Effect Size",
       y = "Pairwise contrast") # Add title and axis labels 


################################### 
##### Queen survival analysis
###################################

##GLMM of queen survival
#Reordering factors by name
queen_survival$treatment <- factor(queen_survival$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

GLMM_queen_survival <- lme(value ~ treatment, random = ~1|farm, method="REML", data = queen_survival)
summary(GLMM_queen_survival)

#Check model residuals
plot(GLMM_queen_survival)

#Post-hoc test
summary(glht(GLMM_queen_survival, emm(pairwise ~ treatment)), test=adjusted(type="bonferroni"))

# Generate pairwise comparisons and plot
marginal_means_queen_survival <- emmeans(GLMM_queen_survival, specs = "treatment")
pairwise_comparisons_queen_survival <- pairs(marginal_means_queen_survival)
pairwise_comparisons_queen_survival_plot <- plot(pairwise_comparisons_queen_survival) +
  geom_vline(xintercept = 0, linetype = "dashed") + # Add dashed vertical line at x=0
  #xlim(-4, 4)+ # Specify x axis limits (adjust these limits as needed)
  labs(title = "Queen survival",
       x = "Estimated Effect Size",
       y = "Pairwise contrast") # Add title and axis labels 

################################### 
##### Total bee numbers analysis
###################################

##GLMM of max yearly active bee numbers
#Reordering factors by name
total_bees_max$treatment <- factor(total_bees_max$treatment, levels = c("gdn","alt gdn","gdn pasture","no gdn"))

GLMM_total_bees <- lme(value ~ treatment, random = ~1|farm, method="REML", data = total_bees_max)
summary(GLMM_total_bees)

#Check model residuals
plot(GLMM_total_bees)

#Post-hoc test
summary(glht(GLMM_total_bees, emm(pairwise ~ treatment)), test=adjusted(type="bonferroni"))

# Generate pairwise comparisons and plot
marginal_means_total_bees <- emmeans(GLMM_total_bees, specs = "treatment")
pairwise_comparisons_total_bees <- pairs(marginal_means_total_bees)
pairwise_comparisons_total_bees_plot <- plot(pairwise_comparisons_total_bees) +
  geom_vline(xintercept = 0, linetype = "dashed") + # Add dashed vertical line at x=0
  #xlim(-4, 4)+ # Specify x axis limits (adjust these limits as needed)
  labs(title = "Total active bees",
       x = "Estimated Effect Size",
       y = "Pairwise contrast") # Add title and axis labels 

#TO PRODUCE ONE PANEL OF ALL PAIRWISE PLOTS

# Combine the plots using patchwork
combined_plot_pairwise <- (pairwise_comparisons_colony_number_plot | pairwise_comparisons_total_bees_plot) / 
  (pairwise_comparisons_new_queens_plot | pairwise_comparisons_queen_survival_plot) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 14, face = "bold"))

#Save combined plot as a PNG and SVG
ggsave(plot=combined_plot_pairwise, filename=file.path(output.path,"BEESTEWARD_model_contrasts.png"), width=8, height=8, dpi=500, bg="white")
ggsave(plot=combined_plot_pairwise, filename=file.path(output.path,"BEESTEWARD_model_contrasts.svg"), width=8, height=8, dpi=500, bg="white")



##########################################################################################################
############# **Additional analysis - testing representativeness of study sites** ########################
##########################################################################################################

#Description:
#Here we test the representativeness of our three farmland study sites by investigating how their nectar production and habitat areas compare 
# with those of  12  farmland sites in southwest UK, sampled in Anonymous et al. 2021

####Loading in data from 12 other farms  in southwest UK to test generalisability of our study sites
willow_phenology <- read.csv(file.path(input.path,"Twelve farms_seasonal nectar.csv")) 
farm_areas <- read.csv(file.path(input.path,"Twelve farms_habitat areas.csv")) 

willow_phenology$month<-factor(willow_phenology$month, levels=c("March", "May", "July", "September"))
farm_areas$Ha<-factor(farm_areas$Habitat, levels=c("Hedgerow", "Margin", "Pasture", "Woodland"))

##Plotting both graphs together as box and whisker plots with values for the three study farms overlain

svg(file.path(output.path,"Nectar and habitat area_boxplots.svg"), bg="white", width=10,height=5) ##Setting up to print plot as svg

par(mfrow=c(1,2))

par(mar=c(5,5,2,1)) #sets the size of the margins around the plot to be slightly larger than the default settings

####Plotting seasonal nectar supply

boxplot((mean_sugar_grams) ~ month,
        data=willow_phenology,
        xlab="Month",
        cex.lab=1.5,
        font.lab=2,
        ylab="",
        ylim = c(0,8000))

title(ylab= expression("Sugar/day/km"^2* "(grams)"), line=3, cex.lab=1.5, font.lab=2)

#Adding the points for the three study farms (Birches, Eastwood and Elmtree)
points(x=0.92,y=0.45009, pch= 22, col="red", bg="red") ##Inserting points to show Birches seasonal values
points(x=2,y=3832, pch= 22, col="red", bg="red")
points(x=3,y=1104, pch= 22, col="red", bg="red")
points(x=4,y=1589, pch= 22, col="red", bg="red")

points(x=0.97,y=30.5, pch= 21, col="blue", bg="blue") ##Inserting points to show Eastwood seasonal values
points(x=2,y=1559.3, pch= 21, col="blue", bg="blue")
points(x=3,y=5408, pch= 21, col="blue", bg="blue")
points(x=4,y=2425.7, pch= 21, col="blue", bg="blue")

points(x=1.03,y=1.72, pch= 24, col="orange", bg="orange") ##Inserting points to show Elmtree seasonal values
points(x=2,y=2521, pch= 24, col="orange", bg="orange")
points(x=3,y=3992.37, pch= 24, col="orange", bg="orange")
points(x=4.05,y=763, pch= 24, col="orange", bg="orange")


####Plotting habitat areas

par(mar=c(5,6,2,1)) #sets the size of the margins around the plot to be slightly larger than the default settings

boxplot((Proportional.area) ~ Habitat,
        data=farm_areas,
        xlab="Habitat",
        cex.lab=1.5,
        font.lab=2,
        las = 1,
        ylab="",
        ylim = c(0,0.8))

title(ylab= expression("Proportional area on farm"), line=3, cex.lab=1.5, font.lab=2)

points(x=1.0,y=0.09, pch= 22, col="red", bg="red") ##Inserting points to show Birches habitat areas
points(x=1.92,y=0.0027, pch= 22, col="red", bg="red")
points(x=3,y=0.47, pch= 22, col="red", bg="red")
points(x=4,y=0.0856, pch= 22, col="red", bg="red")

points(x=1,y=0.0519, pch= 21, col="blue", bg="blue") ##Inserting points to show Eastwood seasonal values
points(x=1.97,y=0, pch= 21, col="blue", bg="blue")
points(x=3,y=0.6158, pch= 21, col="blue", bg="blue")
points(x=4,y=0.0551, pch= 21, col="blue", bg="blue")

points(x=1.05,y=0.0846, pch= 24, col="orange", bg="orange") ##Inserting points to show Elmtree seasonal values
points(x=2.03,y=0.0068, pch= 24, col="orange", bg="orange")
points(x=3,y=0.7052, pch= 24, col="orange", bg="orange")
points(x=3.95,y=0, pch= 24, col="orange", bg="orange")

dev.off()  #finishing the svg print


############################
######### References
##############################

#Rotheray, E.L., Osborne, J.L. & Goulson, D. (2017). Quantifying the food requirements and effects of food stress on bumble bee colony development. J Apicult Res, 56, 288-299.

#Russo, L., DeBarros, N., Yang, S., Shea, K. and Mortensen, D., 2013. Supporting crop pollinators with floral resources: network‐based phenological matching. Ecology and evolution, 3(9), pp.3125-3140.

#[Anonymous 2019] Details omitted for double-blind reviewing.

#[Anonymous 2021] Details omitted for double-blind reviewing.

#Twiston‐Davies, G., Becher, M.A. & Osborne, J.L. (2021). BEE‐STEWARD: A research and decision‐support software for effective land management to promote bumblebee populations. Methods Ecol Evol, 12, 1809-1815.