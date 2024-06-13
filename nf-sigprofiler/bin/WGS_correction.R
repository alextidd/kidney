#----------Corrected for WGS, SBS---------
NanoSeq_mat <- read.table('/Users/yw2/Library/CloudStorage/OneDrive-UniversityofCambridge/PhD/kidney/SBS96/sig_matrix_full.txt',check.names = F,header=T)
mut_rate <- read.table('/Users/yw2/Documents/FARM/kidney/Nanoseq_summary/mutation_rate.txt',check.names = F,header=F)

selected = NanoSeq_mat[,-1]
for (i in 1:ncol(selected )){
  patient = substr(colnames(selected)[i],1,15)
  selected[,i] = selected[,i]*6e9*mut_rate$V5[grep(patient,mut_rate$V1)]/mut_rate$V3[grep(patient,mut_rate$V1)]
}
selected=round(selected)
NanoSeq_mat_corrected = cbind(NanoSeq_mat[,1],selected)
colnames(NanoSeq_mat_corrected)[1] = 'MutationType'

write.table(NanoSeq_mat_corrected,'/Users/yw2/Library/CloudStorage/OneDrive-UniversityofCambridge/PhD/kidney/blood/Nanoseq_blood_matrix_WGS.txt',quote = F,sep='\t',row.names = F)
