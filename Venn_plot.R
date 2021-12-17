library(ggplot2)
library("ggVennDiagram")
ast <- read.csv("D:/pfc/astmarker-teenvsadult.csv", sep=",")
ast <- subset(ast, ast$avg_log2FC > 0.1 & ast$p_val_adj < 0.05)
inn <- read.csv("D:/pfc/InNmarker-teenvsadult.csv", sep=",")
inn <- subset(inn, inn$avg_log2FC > 0.1 & inn$p_val_adj < 0.05)
exn <- read.csv("D:/pfc/exnmarker-teenvsadult.csv", sep=",")
exn <- subset(exn, exn$avg_log2FC > 0.1 & exn$p_val_adj < 0.05)
mic <- read.csv("D:/pfc/Micmarker-teenvsadult.csv", sep=",")
mic <- subset(mic, mic$avg_log2FC > 0.1 & mic$p_val_adj < 0.05)
olig <- read.csv("D:/pfc/oligmarker-teenvsadult.csv", sep=",")
olig <- subset(olig, olig$avg_log2FC > 0.1 & olig$p_val_adj < 0.05)
mic <- paste(mic[,1], sep = "")
ast <- paste(ast[,1], sep = "")
inn <- paste(inn[,1], sep = "")
exn <- paste(exn[,1], sep = "")
olig <- paste(olig[,1], sep = "")

x <- list(
  A = ast,
  B = mic,
  C = olig,
  D = inn,
  E = exn
)


ggVennDiagram(
  x, label_alpha = 0,label = c("count"),label_size = 12,set_size = 15,
  category.names = c("Ast","Mic","Olig", "InN", "Exn")
) +
  ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#FF0000")+
  theme(legend.title = element_text(color = "black", size =30), legend.key.size =  unit(1.3, 'cm'), legend.text = element_text(colour="black", size=22),
        legend.position = "right")
ggsave("D:/vennplot-up.png", height = 24, width = 22, dpi = 300)

ast <- read.csv("D:/pfc/astmarker-teenvsadult.csv", sep=",")
ast <- subset(ast, ast$avg_log2FC < -0.1 & ast$p_val_adj < 0.05)
inn <- read.csv("D:/pfc/InNmarker-teenvsadult.csv", sep=",")
inn <- subset(inn, inn$avg_log2FC < -0.1 & inn$p_val_adj < 0.05)
exn <- read.csv("D:/pfc/exnmarker-teenvsadult.csv", sep=",")
exn <- subset(exn, exn$avg_log2FC < -0.1 & exn$p_val_adj < 0.05)
mic <- read.csv("D:/pfc/Micmarker-teenvsadult.csv", sep=",")
mic <- subset(mic, mic$avg_log2FC < -0.1 & mic$p_val_adj < 0.05)
olig <- read.csv("D:/pfc/oligmarker-teenvsadult.csv", sep=",")
olig <- subset(olig, olig$avg_log2FC < -0.1 & olig$p_val_adj < 0.05)
mic <- paste(mic[,1], sep = "")
ast <- paste(ast[,1], sep = "")
inn <- paste(inn[,1], sep = "")
exn <- paste(exn[,1], sep = "")
olig <- paste(olig[,1], sep = "")

y <- list(
  A = ast,
  B = mic,
  C = olig,
  D = inn,
  E = exn
)


ggVennDiagram(
  y, label_alpha = 0,label = c("count"),label_size = 12,set_size = 15,
  category.names = c("Ast","Mic","Olig", "InN", "Exn")
) +
  ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#008000")+
  theme(legend.title = element_text(color = "black", size = 30), legend.key.size =  unit(1.3, 'cm'), legend.text = element_text(colour="black", size=22),
        legend.position = "right")
ggsave("D:/vennplot-down.png", height = 24, width = 22, dpi = 300)


library(cowplot)

library(patchwork)

library(BiocManager)

library(multtest)

library(metap)

library(dplyr)

library(ggplot2)

library(data.table)

library(magrittr)

library(DESeq2)

library(EnhancedVolcano)

res=read.csv("D:/pfc/InNmarker-latevsearly.csv", header = T, sep=",", row.names = 1)

keyvals <- ifelse(
  
  res$avg_log2FC < -0.1 & res$p_val_adj < 0.05, 'royalblue',
  
  
  
  ifelse(res$avg_log2FC > 0.1 & res$p_val_ad < 0.05, 'red',
         
         
         
         'black'))



keyvals[is.na(keyvals)] <- 'black'

names(keyvals)[keyvals == 'red'] <- 'Upregulated'

names(keyvals)[keyvals == 'black'] <- 'Nonsig'

names(keyvals)[keyvals == 'royalblue'] <- 'Downregulated'

EnhancedVolcano(res,
                
                lab = rownames(res),
                
                x = 'avg_log2FC',
                
                y = 'p_val_adj',
                
                selectLab = "",
                
                xlab = bquote(~Log[2]~ 'fold change'),
                
                title = 'Late vs Early Adolescence',
                
                titleLabSize = 35,
                
                subtitle = "Inhibitory neurons",
                
                subtitleLabSize = 35,
                
                pCutoff = 0.05,
                
                FCcutoff = 0.1,
                
                pointSize = 2,
                
                labSize = 10,
                
                xlim = c(-1.2, 2),
                
                ylim = c(-5, 120),
                
                colCustom = keyvals,
                
                colAlpha = 1,
                
                legendLabSize = 30,
                
                legendIconSize =10,
                
                drawConnectors = TRUE,
                
                widthConnectors = 1.0,
                
                colConnectors = 'black',
                
                arrowheads = FALSE,
                
                gridlines.major = TRUE,
                
                gridlines.minor = FALSE,
                
                border = 'partial',
                
                borderWidth = 1.5,
                
                borderColour = 'black') +theme(axis.text.y = element_text(size = 30), axis.title.y = element_text(size = rel(2), hjust=0.5), axis.title.x = element_text(size = rel(2)), axis.text.x = element_text(size = 30),aspect.ratio = 0.7)

ggsave("D:/volcano-inn.png", dpi = 600)
