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
