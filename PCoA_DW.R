library(vegan)
library(ape)
library(ggplot2)
library(grid)

data <- read.csv("P1.csv", head=TRUE,sep=",",row.names = 1)
data <- t(data)
data[is.na(data)] <- 0
data <- vegdist(data)
pcoa<- pcoa(data, correction = "none", rn = NULL)
groups <- read.csv("G1.csv",sep = ",",header = F,colClasses = c("character"))
groups <- as.list(groups)

#Define some parameters needed for the drawing
length=length(unique(as.character(groups$V1)))
times1=length%/%8
res1=length%%8
times2=length%/%5
res2=length%%5
col1=rep(1:8,times1)
col=c(col1,1:res1)
pich1=rep(c(21:24),times2)
pich=c(pich1,15:(15+res2))
pcoa<- pcoa(data, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]

cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",
                "#ADD1E5")

#Use the following code to perform PCoA analysis and create a PCoA drawing data file
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,groups$V2)
colnames(plotdata) <-c("sample","PC1","PC2","Group")
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
#plotdata$Group <- factor(plotdata$Group,levels = c("setosa","versicolor","virginica"))
plotdata$Group <- factor(plotdata$Group,levels = c("AF","HF","WA","NT"))

##Significance tests for PC1 and PC2
#The following code is used to test the difference between PC1 and PC2 in the PCoA result obtained in the previous step, so this step must be run after the PCoA analysis
library(dplyr)
yf <- plotdata
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1

fit1 <- aov(PC1~Group,data = plotdata)

library(multcomp)
tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)

fit2 <- aov(PC2~Group,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)


test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group)
#test$Group <- factor(test$Group,levels = c("setosa","versicolor","virginica"))
test$Group <- factor(test$Group,levels = c("AF","HF","WA","NT"))

##Phase diagram drawing
#In particular, be sure to draw the top and right side of the phase diagram!!!
#There is a detail here, that is, because the phase diagram is added with the difference test letter, it will cause the axis range of the phase diagram and PCoA scatter diagram to be inconsistent, if it is directly merged, it will cause the image to be distorted, and the box cannot accurately correspond to the distribution of PCoA midpoint.
#Therefore, it is necessary to draw the phase diagram first, and then call the coordinate axis range of the two phase diagrams in the subsequent drawing process of the PCoA diagram to achieve the perfect match of the four diagrams. p1 &lt; - ggplot(plotdata,aes(Group,PC1)) 
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 7,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=20,face = "bold"),
        axis.text.x=element_blank(),
        legend.position = "none")

p3 <- ggplot(plotdata,aes(Group,PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 7,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=20,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        axis.text.y=element_blank(),
        legend.position = "none")

##PCoA plotting
p2<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill=Group),size=8,pch = 21)+
  scale_fill_manual(values=cbbPalette,name = "Group")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=30))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=10),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=10,vjust = 7),
        axis.title.y=element_text(colour='black', size=10,vjust = -2),
        axis.text=element_text(colour='black',size=10),
        legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size=8),
        legend.key=element_blank(),legend.position = c(0.88,0.13),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1,"cm")) +
  guides(fill = guide_legend(ncol = 1))
p2

##PERMANOVA anlysis
otu.adonis=adonis(data~V2,data = groups,distance = "bray")

p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",otu.adonis$aov.tab$Df[1],
                                               "\nR2 = ",round(otu.adonis$aov.tab$R2[1],4),
                                               "\np-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),
            size = 7) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())


##Image Stitching
library(patchwork)
p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)

pdf("PCoA.pdf",height=12,width=15)
p5
png(filename="PCoA.png",res=600,height=7000,width=9000)
p5
dev.off()
