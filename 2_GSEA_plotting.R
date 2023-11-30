##################### 
### GSEA Plots #####
#####################

#### APPLY THE GSEA based on rank script !!!!

### Libraries and WD 
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(fgsea)
library(org.Hs.eg.db)
library(ggExtra)
library(plyr)
library(ggnewscale)
library(treemap)
library(plotly)

source("./functions/minestrone.R")
source("./functions/treemaping.R")

### Test that you have he data to do the plot
# ranked lists are already annotated
head(rnk_g1)

# .gmt file aready imported
head(names(pathways))

# main paths
head(mainPathways1)

#### Plotting the results (nicely)
## Clean path names when plotting 
# save all pathway names
mainPathways1$old_path <- mainPathways1$pathway

# read the names and add more parsing rules as new "problems or patterns could occour"
mainPathways1 <- mainPathways1 %>% 
  mutate(pathway = gsub("%.*%","",pathway)) %>% 
  mutate(pathway = gsub("GO\\:","",pathway)) %>% 
  mutate(pathway = gsub("[0-9]{5,7}$","",pathway)) %>% 
  mutate(pathway = gsub("R\\-HSA\\-[0-9]{5,7}\\.[0-9]{1}$","",pathway)) %>%
  mutate(pathway = gsub("HALLMARK\\_","",pathway)) %>% 
  mutate(pathway = gsub("\\_"," ",pathway))
mainPathways1 <- mainPathways1[order(mainPathways1$pathway),]
mainPathways1 <- mainPathways1[!duplicated(mainPathways1$pathway),]

### NES values
nessie <- ggplot(mainPathways1)+
  geom_col(aes(x = reorder(pathway,NES), y = NES, fill = NES < 0))+
  geom_line(aes(group=1, x = reorder(pathway,NES), y = size/100), col = "darkgreen")+
  theme_classic()+
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#D9DDDC"))+
  labs(title = "Enriched pathways", x = "Pathway")+
  scale_fill_manual(values = c("#0077B6","#FF0800"))+
  coord_flip()
nessie

# Save the plot
ggsave(filename = "../plots/GSEA_NES.tiff", plot = nessie,
       units = "cm", height = 20, width = 30)

### Dotplot 
# follow this when more than a group is being studied, the example can handle 4 gorups but more
# can be added
# Create the data frame
names_gseas <- c("Group 1","Group 2","Group 3","Group 4")

mainPathways1$group <- names_gseas[1]
mainPathways2$group <- names_gseas[2]
mainPathways3$group <- names_gseas[3]
mainPathways4$group <- names_gseas[4]

## At first, only P-Value, size and group will be required, create the data
cnames <- c("pathway","old_path","pval","padj","size","NES","group")
mainPathways1 <- mainPathways1 %>% 
  select(all_of(cnames))

mainPathways2 <- mainPathways2 %>% 
  select(all_of(cnames))

mainPathways3 <- mainPathways3 %>% 
  select(all_of(cnames))

mainPathways4 <- mainPathways4 %>% 
  select(all_of(cnames)) 

list_of_gseas <- list(mainPathways1,mainPathways2,mainPathways3,mainPathways4)
gseas_to_dot <- ldply(list_of_gseas,rbind)

## Get the gene ratio
gseas_to_dot$gene_ratio <- 0

for (i in 1:length(gseas_to_dot$old_path)){
  look_to_this <- which(names(pathways) == gseas_to_dot$old_path[i])
  tamany <- length(unlist(pathways[look_to_this]))
  gseas_to_dot$gene_ratio[i] <- gseas_to_dot$size[i]/tamany
}

gseas_to_dot <- gseas_to_dot %>% 
  filter(NES > 0)

# Set the order of the dotplot by group 
gseas_to_dot <- gseas_to_dot %>% 
  arrange(group)

rownames(gseas_to_dot) <- NULL
gseas_to_dot$order <- rownames(gseas_to_dot)

for (i in 1:length(gseas_to_dot$order)){
  gseas_to_dot$order[i] <- min(which(gseas_to_dot$pathway[i] == gseas_to_dot$pathway))
}

gseas_to_dot$order <- as.numeric(gseas_to_dot$order)

# Do the plot
sulley <- ggplot(gseas_to_dot)+
  geom_point(aes(x = group, y = reorder(pathway, order), size = gene_ratio, color = padj))+ ### FORCE THE ORDER AT THE Y AXIS !!!
  scale_color_gradient(low = "red", high = "blue")+
  labs(title = "GSEA enrichment",x = "Pacient group", y = "Path", 
       size = "GeneRatio", color = "Adjusted P.Val")+
  theme_bw()+
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "left",
        plot.margin = unit(c(1,6,1,1),"cm"))+
  geom_tile(aes(x = group, y = -1, fill = group), show.legend = F)+
  new_scale_colour()+
  scale_fill_manual(values = c("#0d0887","#ff1493","#00ff00","#a020f0"))
sulley

# Save the plot
ggsave(filename = "../plots/GSEA_dots.tiff", plot = sulley,
       units = "cm", height = 20, width = 30)

### Letter soup
soup1 <- minestrone_for_GSEA(collapsed_list = collapsedPathways1, original_data = fgseaRes1)

# set the size
soup1$size <- 1

# set the size (nicely)
x <- soup1$child
soup1$size <- fgseaRes1$size[match(x, fgseaRes11$Description)]

# Do the plot and save it
png("../plots/GSEA_collapse.png", height = 1400, width = 1400)
treemaping(soup = soup1, graph_title = "Group1 vs Group2 GSEA collapse")
dev.off()

