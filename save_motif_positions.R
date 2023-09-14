library(glue)

base_dir = glue('{dirname(rstudioapi::getActiveDocumentContext()$path)}/motif_positions')


motif_positions_all <- readRDS("/gstore/project/ar_ligands/NE/reprogram_seq/multiome_arrayed/OUTPUT/doubletremoved/Annotations/Motif-Positions-In-Peaks.rds")
for (motif_name in c('Gata6.Zf_109', 'Nkx2.1.Homeobox_182')){
  motif_positions<- motif_positions_all[[motif_name]]
  write.csv(motif_positions, glue('{base_dir}/{motif_name}.csv'))
}


motif_name_AR = 'AR_689'

motif_positions_all = readRDS('/gstore/project/ar_ligands/AR/scRNAseq/nonpipeline/OUTPUT/ArchRProject/Annotations/Motif-Positions-In-Peaks.rds')
motif_positions = motif_positions_all[[motif_name_AR]]

write.csv(motif_positions, glue('{base_dir}/{motif_name_AR}.csv')) 
