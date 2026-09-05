neuropeptide_genes <- c(
  "Calcrl", "Calcr", "Ramp1", "Ramp2", "Ramp3", "Tacr1", "Tacr2", "Mrgpra1", 
  "Mrgprb2", "Npy1r", "Npy2r", "Npy4r", "Npy5r", "Sstr1", "Sstr2", "Sstr3", 
  "Sstr4", "Sstr5", "Vipr1", "Vipr2", "Nmur1", "Nmur2", "Galr1", "Galr2", 
  "Galr3", "Nk1r", "Nk2r", "Nk3r", "Oprm1", "Oprd1", "Oprk1", "Oprl1", "Bdkrb1", "Bdkrb2", "Pac1", 
  "Vpac1", "Vpac2"
)

load("t.Rda")
levels(t_sce)
new.cluster.ids <- c("Proliferating CD4",
                     "Naive Th CD4",
                     "Gut-resident CD8",
                     "Type I IFN-stimulated CD8",
                     "Effector Memory CD8",
                     "Activated Cytotoxic CD8",
                     "Proliferating CD4",
                     "Tcm CD4",
                     "Activated Treg",
                     "Activated Th1",
                     "Effector γδT",
                     "LN-resident RORγ+ Treg",
                     "Effector-like Killer γδT",
                     "Type I IFN-stimulated CD8",
                     "Activated Memory CD8",
                     "Naive CD4",
                     "Tfh", 
                     "Tfh",
                     "γδT",
                     "ILC2"
)

new.cluster.ids <- c("Activated MZB",
                     "GC LZ B",
                     "IFN-Responsive B",
                     "Early Pro-B",
                     "Plasma cell",
                     "Activated B",
                     "Early Pro-B",
                     "Activated Memory B",
                     "Activated Naive B",
                     "Activated Naive B",
                     "Naive B",
                     "Plasmablast",
                     "GC B",
                     "IFN-Responsive B",
                     "Macrophage-like B",
                     "Activated Memory B",
                     "Activated Plasma cell",
                     "Plasma cell",
                     "GC B",
                     "GC B"
)
names(new.cluster.ids) <- levels(b_sce)
b_sce <- RenameIdents(b_sce, new.cluster.ids)
b_sce$celltype=b_sce@active.ident

FeaturePlot(t_sce,feature="Npy2r")

DimPlot(b_sce, reduction = "umap",raster = F,
        label.size = 0.5,pt.size = 0.1,
        cols = colors,raster.dpi=c(1024,1024))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = "right")+
  labs(title = "B cell type")

DotPlot(b_sce, features = neuropeptide_genes,group.by="celltype")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  scale_size_continuous(range = c(0, 10))

phe=b_sce@meta.data
write_csv(phe, "b_metadata_all.csv")
