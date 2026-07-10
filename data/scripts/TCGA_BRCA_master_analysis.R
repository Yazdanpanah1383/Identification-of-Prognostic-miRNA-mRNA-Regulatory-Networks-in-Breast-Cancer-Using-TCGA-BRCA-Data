####نصب پکیج ها ####
install.packages("data.table")
install.packages("survival")
install.packages("survminer")
install.packages("ggplot2")
install.packages("igraph")

####بارگذاری پکیج ها ####
library(data.table)
library(survival)
library(survminer)
library(ggplot2)

####دانلود فایل ####
file.choose()
mirna<-fread(file.choose())

####بررسی ابعاد داده ####
dim(mirna)
####بررسی ساختار فایل ####
head(mirna)

####فهمسدن تعداده ژن های یکتا  ####
colnames(mirna)
length(unique(mirna$Target))

####هر ژن چند بار تکرار شده####

gene_count <- sort(table(mirna$Target), decreasing = TRUE)
head(gene_count, 20)


View(gene_count)
#### TCGAدانلود فایل ####
file.choose()
rna<-fread(file.choose())
dim(rna)
head(rna)
####پیدا کردن 9 ژن داخل TCGA####
hub_genes <- c(
  "CDK6","HMGA2","HMGB1",
  "HUWE1","IMPDH1","MYC",
  "PDGFRA","PER1","SMARCC1"
)

hub_genes %in% rna$sample
####استخراج ژن ها####
hub_exp <- rna[sample %in% hub_genes]
dim(hub_exp)

#### پیدا کردن نمونه های توموری و سرطانی####
samples <- colnames(hub_exp)[-1]

tumor <- grep("-01$", samples, value = TRUE)
normal <- grep("-11$", samples, value = TRUE)

length(tumor)
length(normal)

length(tumor)
length(normal)

####ژنMYC####
myc <- hub_exp[sample == "MYC"]

tumor_exp <- as.numeric(myc[, ..tumor])
normal_exp <- as.numeric(myc[, ..normal])

boxplot(
  tumor_exp,
  normal_exp,
  names = c("Tumor", "Normal"),
  main = "MYC expression in TCGA-BRCA",
  ylab = "Expression (log2 RSEM)"
)
wilcox.test(tumor_exp, normal_exp)

####محاسبه  پی ولیو هر 9 ژن ####
results <- data.frame(
  Gene = character(),
  Tumor_Median = numeric(),
  Normal_Median = numeric(),
  P_value = numeric(),
  Direction = character()
)

for (g in hub_genes) {
  
  gene_row <- hub_exp[sample == g]
  
  tumor_exp <- as.numeric(gene_row[, ..tumor])
  normal_exp <- as.numeric(gene_row[, ..normal])
  
  p <- wilcox.test(tumor_exp, normal_exp)$p.value
  
  med_tumor <- median(tumor_exp, na.rm = TRUE)
  med_normal <- median(normal_exp, na.rm = TRUE)
  
  direction <- ifelse(
    med_tumor > med_normal,
    "Up in Tumor",
    "Down in Tumor"
  )
  
  results <- rbind(
    results,
    data.frame(
      Gene = g,
      Tumor_Median = med_tumor,
      Normal_Median = med_normal,
      P_value = p,
      Direction = direction
    )
  )
}

results
####مرتب کردن بر اساس معنی داری اماری####
results <- results[order(results$P_value), ]
results
####ذخیره جدول ####
write.csv(results,
          "Hub_genes_TCGA_BRCA_results.csv",
          row.names = FALSE)
####ژن های برتر کدام میکرو هدف میگیرند####
key_genes <- c("HMGA2", "PER1", "PDGFRA", "IMPDH1")

key_network <- mirna[Target %in% key_genes]
dim(key_network)
names(key_network)
key_network[, .(ID, Target)]
####ساخت جدول و شبکه####
library(igraph)

edges <- key_network[, .(ID, Target)]
edges

net <- graph_from_data_frame(edges, directed = TRUE)

V(net)$color <- ifelse(
  V(net)$name %in% unique(edges$ID),
  "skyblue",
  "tomato"
)
V(net)$size <- ifelse(
  V(net)$name %in% unique(edges$ID),
  35,
  30
)
plot(
  net,
  vertex.label.cex = 0.9,
  vertex.label.color = "black",
  edge.arrow.size = 0.4,
  layout = layout_with_fr(net)
)


png(
  "miRNA_mRNA_network.png",
  width = 2200,
  height = 1800,
  res = 300
)

plot(
  net,
  vertex.label.cex = 1,
  vertex.label.color = "black",
  edge.arrow.size = 0.4,
  layout = layout_with_fr(net)
)

dev.off()
####ارتباط بقا ####

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")




BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

clinical <- GDCquery_clinic(
  project = "TCGA-BRCA",
  type = "clinical"
)

dim(clinical)
head(clinical)

names(clinical)

save.image("TCGA_BRCA_Project.RData")
savehistory("TCGA_BRCA_history.Rhistory")
