library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggsci)
mouseImmu <- ImmGenData()
mouseRNA <- MouseRNAseqData()
source("colors.R")

DimPlot(seurat_joint_harmony, reduction = "umap",raster = F,
        label.size = 0.5,pt.size = 0.1,
        cols = colors,raster.dpi=c(1024,1024))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = "right")+
  labs(title = "Immune_Gut")

phe=t_sce@meta.data
write_csv(phe, "t_metadata_all.csv")

seurat_joint_harmony@active.ident=as.factor(seurat_joint_harmony$mouseImmu)
seurat_joint_harmony=JoinLayers(seurat_joint_harmony)
save(seurat_joint_harmony,file="joined_seurat.Rda")

cluster_markers <- FindAllMarkers(object = b_sce, only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  thresh.use = 0.25)

# display the top 30 per cluster
cluster_markers %>%
  dplyr::group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::select(cluster, gene, everything()) %>%
  DT::datatable(cluster_markers, filter = "top")

write_tsv(cluster_markers, file.path("~/BR seq/code/", "cluster_markers_b_annoed.tsv"))

# Display a heatmap of the top 10 per cluster
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(seurat_joint_harmony, features = top10$gene, group.colors = seurat_joint_harmony@misc$colours) +
  NoLegend() +
  scale_fill_gradientn(colors = c("#2166AC", "#E5E0DC", "#B2182B"))

rm(t_sce)
t_sce=seurat_joint_harmony[,seurat_joint_harmony@meta.data$mouseRNA %in% c("T cells")]
b_sce=seurat_joint_harmony[,seurat_joint_harmony@meta.data$mouseRNA %in% c("B cells")]
myeloid=seurat_joint_harmony[,seurat_joint_harmony@meta.data$celltype %in% c("cDC","NK cells","M2","PMN MDSC","ILC","Nav Mac","M MDSC","Mig DC","pDC")]

save(t_sce,file="t.Rda")
save(b_sce,file="b.Rda")
save(myeloid,file="myeloid.Rda")

t_sce <- AddModuleScore(t_sce, features = list(neuropeptide_genes), name = "neuropeptide_score")
t_sce <- AddModuleScore(t_sce, features = list(camp_pathway_genes), name = "camp_score")
t_sce <- AddModuleScore(t_sce, features = list(nfkb_pathway_genes), name = "nfkb_score")
t_sce <- AddModuleScore(t_sce, features = list(mtor_pathway_genes), name = "mtor_score")
t_sce <- AddModuleScore(t_sce, features = list(cgrp_related_genes), name = "cgrp_score")
t_sce <- AddModuleScore(t_sce, features = list(vip_related_genes), name = "vip_score")

b_sce <- AddModuleScore(b_sce, features = list(neuropeptide_genes), name = "neuropeptide_score")
b_sce <- AddModuleScore(b_sce, features = list(camp_pathway_genes), name = "camp_score")
b_sce <- AddModuleScore(b_sce, features = list(nfkb_pathway_genes), name = "nfkb_score")
b_sce <- AddModuleScore(b_sce, features = list(mtor_pathway_genes), name = "mtor_score")
b_sce <- AddModuleScore(b_sce, features = list(cgrp_related_genes), name = "cgrp_score")
b_sce <- AddModuleScore(b_sce, features = list(vip_related_genes), name = "vip_score")

myeloid <- AddModuleScore(myeloid, features = list(neuropeptide_genes), name = "neuropeptide_score")
myeloid <- AddModuleScore(myeloid, features = list(camp_pathway_genes), name = "camp_score")
myeloid <- AddModuleScore(myeloid, features = list(nfkb_pathway_genes), name = "nfkb_score")
myeloid <- AddModuleScore(myeloid, features = list(mtor_pathway_genes), name = "mtor_score")
myeloid <- AddModuleScore(myeloid, features = list(cgrp_related_genes), name = "cgrp_score")
myeloid <- AddModuleScore(myeloid, features = list(vip_related_genes), name = "vip_score")


library(Seurat)
library(ggplot2)
library(dittoSeq)
library(viridis)

neuropeptide_related_genes <- c(
  # 主要神经肽
  "Vip",    # Vasoactive Intestinal Peptide (VIP)
  "Tac1",   # 编码 Substance P
  "Calca",  # 编码 Calcitonin Gene-Related Peptide (cGRP)
  "Avp",    # 编码 Arginine Vasopressin (AVP)
  "Oxt",    # 编码 Oxytocin
  "Npy",    # 编码 Neuropeptide Y (NPY)
  "Pdyn",   # 编码 Prodynorphin (PDYN)
  "Penk",   # 编码 Proenkephalin (PENK)
  "Som",    # 编码 Somatostatin (SST)
  "Mrgprx2",# 编码 Mas-related G-protein-coupled receptor X2
  "Cckbr",  # 编码 Cholecystokinin B receptor (CCKBR)
  
  # 神经肽受体
  "Npr1",   # 编码 Natriuretic peptide receptor 1 (NPR1), 对肽类信号的响应
  "Npr2",   # 编码 Natriuretic peptide receptor 2 (NPR2)
  "Adcy1",  # 编码 腺苷酸环化酶1，涉及cAMP信号转导
  "Adcy8",  # 编码 腺苷酸环化酶8，涉及cAMP信号转导
  "Prkar1a",# 编码蛋白激酶A的调节亚单位
  "Prkacb", # 编码蛋白激酶A的催化亚单位
  "Gnas",   # 编码G蛋白α亚单位，涉及神经肽的下游信号转导
  "Plcb1",  # 编码磷脂酶Cβ1，涉及IP3信号通路
  "Plcb4",  # 编码磷脂酶Cβ4
  "Cnr1",   # 编码 大麻素受体1，涉及痛觉调节等
  "Trpv1",  # 编码温度和痛觉感受器
  "TrpA1",  # 编码 TRPA1 痛觉受体
  "Ntrk1",  # 编码神经营养因子受体 TrkA，参与疼痛感知
  "Fos",    # 脱氧核糖核酸转录因子，涉及神经应激反应
  "Jun",    # 早期基因，调节神经应激和炎症反应
  "Egr1",   # 早期生长反应基因1，涉及神经系统的反应性变化
  "P2rx3",  # ATP受体，参与神经传导
  "P2ry1",  # ATP受体，参与神经传导
  "Tnfaip3", # 在免疫反应中的重要因子
  "Cxcl12"  # 趋化因子，调节神经细胞的迁移
)

camp_related_genes <- c(
  # cAMP合成
  "Adcy1",  # 编码腺苷酸环化酶1，负责将ATP转化为cAMP
  "Adcy2",  # 编码腺苷酸环化酶2
  "Adcy3",  # 编码腺苷酸环化酶3
  "Adcy4",  # 编码腺苷酸环化酶4
  "Adcy5",  # 编码腺苷酸环化酶5
  "Adcy6",  # 编码腺苷酸环化酶6
  "Adcy7",  # 编码腺苷酸环化酶7
  "Adcy8",  # 编码腺苷酸环化酶8
  "Adcy9",  # 编码腺苷酸环化酶9
  
  # cAMP降解
  "Pde1a",  # 编码磷酸二酯酶1a，降解cAMP
  "Pde2a",  # 编码磷酸二酯酶2a，降解cAMP
  "Pde3a",  # 编码磷酸二酯酶3a，降解cAMP
  "Pde4a",  # 编码磷酸二酯酶4a，降解cAMP
  "Pde5a",  # 编码磷酸二酯酶5a，降解cAMP
  "Pde7a",  # 编码磷酸二酯酶7a，降解cAMP
  "Pde8a",  # 编码磷酸二酯酶8a，降解cAMP
  
  # cAMP效应器（蛋白激酶A）
  "Prkar1a",  # 编码蛋白激酶A的调节亚单位
  "Prkacb",   # 编码蛋白激酶A的催化亚单位
  "Prkaca",   # 编码蛋白激酶A的催化亚单位
  "Prkacg",   # 编码蛋白激酶A的催化亚单位
  
  # cAMP相关G蛋白
  "Gnas",   # 编码G蛋白α亚单位，涉及cAMP信号通路
  "Gnb1",   # 编码G蛋白β亚单位
  "Gng2",   # 编码G蛋白γ亚单位
  
  # 受体与cAMP信号转导
  "Adora1",  # 编码腺苷受体A1，受体与cAMP信号通路相关
  "Adora2a", # 编码腺苷受体A2a，受体与cAMP信号通路相关
  "Adora2b", # 编码腺苷受体A2b，受体与cAMP信号通路相关
  "Drd1",    # 编码多巴胺受体D1，激活cAMP信号
  "Drd5",    # 编码多巴胺受体D5，激活cAMP信号
  "Cnr1",    # 编码大麻素受体1，涉及cAMP信号
  
  # 与cAMP信号相关的下游效应器
  "Creb1",   # 编码cAMP反应元件结合蛋白，受cAMP信号调控
  "Fos",     # 与cAMP调节的基因表达相关
  "Jun",     # 与cAMP信号转导的调节相关
  "Egr1",    # 早期基因1，参与cAMP信号的转录调控
  
  # 其他调节因子
  "Rgs2",    # 编码调节G蛋白的因子，参与cAMP信号的调控
  "Rap1a",   # 编码Ras相关蛋白，调节cAMP信号
  "Cacna1c", # 编码L型钙通道，涉及cAMP的下游效应
  "Mapk1",   # 编码MAPK1，参与cAMP相关信号通路
  "Mapk3",   # 编码MAPK3，参与cAMP相关信号通路
  "Stat3",   # 编码信号转导与转录激活因子3，cAMP信号的下游调控
  "Tnfaip3", # 在免疫反应中的作用，与cAMP信号相关
  "Cxcl12"   # 趋化因子，cAMP调节的免疫细胞迁移因子
)

nfkb_related_genes <- c(
  # NF-κB家族成员
  "Nfkb1",   # 编码NF-κB p50亚单位，NF-κB家族的主要成员
  "Nfkb2",   # 编码NF-κB p52亚单位
  "RelA",    # 编码NF-κB p65亚单位，是最常见的NF-κB亚单位之一
  "RelB",    # 编码NF-κB RelB亚单位，参与免疫和发育过程
  "C-Rel",   # 编码C-Rel亚单位，是NF-κB家族的成员之一
  
  # IκB抑制因子
  "Ikbkb",   # 编码IκB激酶β，关键的NF-κB激活调节因子
  "Ikbka",   # 编码IκB激酶α，调节NF-κB通路的关键因子
  "Ikbke",   # 编码IκBε，作为NF-κB通路的调节因子
  "Ikbzg",   # 编码IκBζ，涉及NF-κB抑制
  
  # NF-κB信号传导调控基因
  "Tnfa",    # 编码肿瘤坏死因子α，激活NF-κB信号通路
  "Tnfrsf1a",# 编码肿瘤坏死因子受体1A，介导NF-κB信号传导
  "Tnfrsf1b",# 编码肿瘤坏死因子受体1B，介导NF-κB信号传导
  "Map3k1",  # 编码MAP激酶激酶激酶1，参与NF-κB信号转导
  "Map3k7",  # 编码TAK1，NF-κB信号通路的关键激酶
  "Ripk1",   # 编码受体相互作用蛋白激酶1，调节NF-κB的激活
  "Ripk3",   # 编码受体相互作用蛋白激酶3，参与调控NF-κB通路
  
  # NF-κB下游基因
  "Fos",     # 编码Fos蛋白，是NF-κB调节的早期基因之一
  "Jun",     # 编码Jun蛋白，与Fos共同形成AP-1，参与NF-κB相关的基因转录
  "Egr1",    # 编码早期生长反应基因1，NF-κB信号的转录效应因子
  "Cxcl2",   # 编码趋化因子CXCL2，是NF-κB下游的炎症因子
  "Il1b",    # 编码白介素1β，是NF-κB调控的主要炎症因子
  "Il6",     # 编码白介素6，参与NF-κB调控的免疫反应
  "Ccl2",    # 编码趋化因子C-C基序配体2（MCP-1），NF-κB介导的免疫细胞趋化因子
  "Tnf",     # 编码肿瘤坏死因子，是NF-κB信号通路的重要调节因子
  "Nos2",    # 编码诱导型一氧化氮合成酶，是NF-κB调控的产物
  "Il8",     # 编码白介素8，NF-κB激活下的炎症反应因子
  "Nfkbia",  # 编码IκBα，NF-κB的抑制因子，负向调节NF-κB信号
  
  # NF-κB调节基因
  "Traf2",   # 编码TNF受体相关因子2，参与NF-κB的激活
  "Traf6",   # 编码TNF受体相关因子6，NF-κB信号通路的关键调节因子
  "Relb",    # 编码RelB，NF-κB的一个亚单位
  "Akt1",    # 编码蛋白激酶B（Akt1），参与NF-κB的正向调控
  "P53",     # 编码肿瘤抑制蛋白p53，调节NF-κB和细胞应激反应
  "Stat3",   # 编码信号转导与转录激活因子3，涉及NF-κB的信号转导
  "Pik3r1",  # 编码PI3K的调节亚单位，参与NF-κB信号的上游调节
  "Mapk14"   # 编码p38 MAP激酶，参与NF-κB信号通路的调节
)

mtor_related_genes <- c(
  # mTOR复合物成员
  "Mtor",    # 编码mTOR，mTORC1和mTORC2复合物的核心组件
  "Raptor",  # 编码Raptor，是mTORC1的必需组分，调节mTORC1的功能
  "Rictor",  # 编码Rictor，是mTORC2的必需组分，调节mTORC2的功能
  "Mlst8",   # 编码mLST8（GβL），参与mTORC1的稳定性和功能
  
  # mTOR信号传导的上游激活因子
  "Akt1",    # 编码Akt1，激活mTORC1信号通路
  "Tsc1",    # 编码Tuberin，是TSC1/TSC2复合物的组成部分，负向调节mTORC1
  "Tsc2",    # 编码Tuberin，参与TSC1/TSC2复合物，抑制mTORC1活性
  "Pten",    # 编码PTEN，负向调节PI3K/Akt/mTOR通路
  
  # 直接调节mTORC1的因子
  "S6k1",    # 编码S6激酶1，mTORC1的下游靶标，参与蛋白质合成
  "4E-BP1",  # 编码eIF4E结合蛋白1，mTORC1的下游效应分子，抑制蛋白合成
  "Rheb",    # 编码Rheb，小GTP酶，直接激活mTORC1
  "Gbp2",    # 编码GTP结合蛋白2，参与mTORC1信号传导
  
  # mTORC2的下游效应
  "Akt2",    # 编码Akt2，mTORC2通过Akt2调节细胞存活和代谢
  "Prkci",   # 编码PKCι（蛋白激酶C），是mTORC2的下游效应分子，参与细胞代谢和存活
  
  # 参与mTOR调节的转录因子
  "Hif1a",   # 编码低氧诱导因子1α，mTOR调节下的转录因子，参与细胞代谢调节
  "Myc",     # 编码Myc，转录因子，mTOR调节细胞增殖和生长
  "Pax3",    # 编码Paired box蛋白3，调节肌肉发育过程中mTOR的功能
  
  # 参与自噬的基因
  "Atg1",    # 编码自噬相关基因1，mTOR负向调节自噬
  "Atg13",   # 编码自噬相关基因13，参与mTOR调控的自噬过程
  "Atg5",    # 编码自噬相关基因5，参与mTOR调节的自噬过程
  
  # mTOR与能量代谢相关的基因
  "Cpt1a",   # 编码肉碱棕榈酰转移酶1，mTOR参与脂肪酸氧化调节
  "Sirt1",   # 编码去乙酰化酶Sirtuin1，mTOR与Sirt1调节代谢功能
  "Nrf2",    # 编码Nrf2，抗氧化反应的转录因子，mTOR调节代谢和应激反应
  "Ampk",    # 编码AMP-活化蛋白激酶，mTOR调节代谢与能量平衡
  
  # 参与mTOR相关的细胞周期基因
  "CyclinD1",# 编码细胞周期蛋白D1，mTOR在细胞周期调节中的作用
  "Cdkn1a",  # 编码p21，细胞周期的抑制因子，mTOR对细胞增殖的调节
  "CyclinE1",# 编码细胞周期蛋白E1，mTOR调节G1/S期的过渡
  "Rb1",     # 编码视网膜母细胞瘤蛋白，调控细胞周期的关键基因
  
  # 其他与mTOR调控相关的基因
  "Gsk3b",   # 编码GSK3β，mTOR调节的关键信号分子
  "Mcl1",    # 编码Mcl-1，抗凋亡蛋白，mTOR信号调节细胞生存
  "Fasl",    # 编码FasL，凋亡相关的蛋白，mTOR参与细胞凋亡调控
  "Bcl2",    # 编码Bcl-2，抗凋亡蛋白，mTOR调节细胞生存
  "Bcl2l1"   # 编码Bcl-xL，调节细胞凋亡，mTOR的下游效应
)

cgrp_vip=c("Calcrl","Calcr","Ramp1","Ramp2","Ramp3","Tacr1","Tacr2","Npy1r","Npy2r","Npy4r","Npy5r","Sstr1","Sstr2","Sstr3",
           "Sstr4","Sstr5","Vipr1","Vipr2","Nmur1","Nmur2","Galr1","Galr2","Galr3")

DotPlot(seurat_joint_harmony, features = cgrp_vip,group.by="celltype")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  scale_size_continuous(range = c(0, 10))

DotPlot(t_sce, features = mtor_related_genes,group.by="celltype")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  scale_size_continuous(range = c(0, 10))

DotPlot(b_sce, features = mtor_related_genes,group.by="celltype")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  scale_size_continuous(range = c(0, 10))

DotPlot(myeloid, features = mtor_related_genes,group.by="celltype")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  scale_size_continuous(range = c(0, 10))

DotPlot(seurat_joint_harmony, features = mtor_related_genes,group.by="Group")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  scale_size_continuous(range = c(0, 10))
