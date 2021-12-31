library(DESeq2)

# Load pregenerated count data
count_matrix <- as.matrix(read.csv('../../../Analysis/feature_counts_all/counts.csv', header = TRUE, row.names = 1))

# Generate a metadata table. type is analysis type of immunoprecipitated mRNA (IP) or whole animal input (whole). tissue is strain JAC127 (muscle), JAC382 (intestine), JAC126 (pan-neuronal), JAC433 (serotonergic neuron), or JAC379 (dopaminergic neuron).
metadata <- data.frame(row.names = colnames(count_matrix), type = factor(c('IP', 'IP', 'IP', 'IP', 'IP', 'IP', 'IP', 'IP', 'IP', 'IP', 'whole', 'whole', 'whole', 'whole', 'whole', 'whole', 'whole', 'whole', 'whole', 'whole')), tissue = factor(c('127', '127', '382', '382', '126', '126', '433', '433', '379', '379', '127', '127', '382', '382', '126', '126', '433', '433', '379', '379')))

# Run differential expresion between analysis types, all as IP / whole
deseq_data_type_all <- DESeqDataSetFromMatrix(count_matrix, metadata, ~ type)
write.csv(as.data.frame(results(DESeq(deseq_data_type_all), contrast = c('type', 'IP', 'whole'))), '../../../Analysis/diff_exp_all/type/results_all.csv')

deseq_data_type_127 <- DESeqDataSetFromMatrix(count_matrix[, c(1, 2, 11, 12)], metadata[c(1, 2, 11, 12), ], ~ type)
write.csv(as.data.frame(results(DESeq(deseq_data_type_127), contrast = c('type', 'IP', 'whole'))), '../../../Analysis/diff_exp_all/type/results_127.csv')

deseq_data_type_382 <- DESeqDataSetFromMatrix(count_matrix[, c(3, 4, 13, 14)], metadata[c(3, 4, 13, 14), ], ~ type)
write.csv(as.data.frame(results(DESeq(deseq_data_type_382), contrast = c('type', 'IP', 'whole'))), '../../../Analysis/diff_exp_all/type/results_382.csv')

deseq_data_type_126 <- DESeqDataSetFromMatrix(count_matrix[, c(5, 6, 15, 16)], metadata[c(5, 6, 15, 16), ], ~ type)
write.csv(as.data.frame(results(DESeq(deseq_data_type_126), contrast = c('type', 'IP', 'whole'))), '../../../Analysis/diff_exp_all/type/results_126.csv')

deseq_data_type_433 <- DESeqDataSetFromMatrix(count_matrix[, c(7, 8, 17, 18)], metadata[c(7, 8, 17, 18), ], ~ type)
write.csv(as.data.frame(results(DESeq(deseq_data_type_433), contrast = c('type', 'IP', 'whole'))), '../../../Analysis/diff_exp_all/type/results_433.csv')

deseq_data_type_379 <- DESeqDataSetFromMatrix(count_matrix[, c(9, 10, 19, 20)], metadata[c(9, 10, 19, 20), ], ~ type)
write.csv(as.data.frame(results(DESeq(deseq_data_type_379), contrast = c('type', 'IP', 'whole'))), '../../../Analysis/diff_exp_all/type/results_379.csv')

# Run differential expression between tissue types, with file names numerator_vs_denominator
deseq_data_tissue_127_vs_382 <- DESeqDataSetFromMatrix(count_matrix[, c(1, 2, 3, 4)], metadata[c(1, 2, 3, 4), ], ~ tissue)
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_127_vs_382), contrast = c('tissue', '127', '382'))), '../../../Analysis/diff_exp_all/tissue/results_127_vs_382.csv')
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_127_vs_382), contrast = c('tissue', '382', '127'))), '../../../Analysis/diff_exp_all/tissue/results_382_vs_127.csv')

deseq_data_tissue_127_vs_126 <- DESeqDataSetFromMatrix(count_matrix[, c(1, 2, 5, 6)], metadata[c(1, 2, 5, 6), ], ~ tissue)
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_127_vs_126), contrast = c('tissue', '127', '126'))), '../../../Analysis/diff_exp_all/tissue/results_127_vs_126.csv')
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_127_vs_126), contrast = c('tissue', '126', '127'))), '../../../Analysis/diff_exp_all/tissue/results_126_vs_127.csv')

deseq_data_tissue_382_vs_126 <- DESeqDataSetFromMatrix(count_matrix[, c(3, 4, 5, 6)], metadata[c(3, 4, 5, 6), ], ~ tissue)
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_382_vs_126), contrast = c('tissue', '382', '126'))), '../../../Analysis/diff_exp_all/tissue/results_382_vs_126.csv')
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_382_vs_126), contrast = c('tissue', '126', '382'))), '../../../Analysis/diff_exp_all/tissue/results_126_vs_382.csv')

deseq_data_tissue_126_vs_433 <- DESeqDataSetFromMatrix(count_matrix[, c(5, 6, 7, 8)], metadata[c(5, 6, 7, 8), ], ~ tissue)
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_126_vs_433), contrast = c('tissue', '126', '433'))), '../../../Analysis/diff_exp_all/tissue/results_126_vs_433.csv')
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_126_vs_433), contrast = c('tissue', '433', '126'))), '../../../Analysis/diff_exp_all/tissue/results_433_vs_126.csv')

deseq_data_tissue_126_vs_379 <- DESeqDataSetFromMatrix(count_matrix[, c(5, 6, 9, 10)], metadata[c(5, 6, 9, 10), ], ~ tissue)
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_126_vs_379), contrast = c('tissue', '126', '379'))), '../../../Analysis/diff_exp_all/tissue/results_126_vs_379.csv')
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_126_vs_379), contrast = c('tissue', '379', '126'))), '../../../Analysis/diff_exp_all/tissue/results_379_vs_126.csv')

deseq_data_tissue_433_vs_379 <- DESeqDataSetFromMatrix(count_matrix[, c(7, 8, 9, 10)], metadata[c(7, 8, 9, 10), ], ~ tissue)
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_433_vs_379), contrast = c('tissue', '433', '379'))), '../../../Analysis/diff_exp_all/tissue/results_433_vs_379.csv')
write.csv(as.data.frame(results(DESeq(deseq_data_tissue_433_vs_379), contrast = c('tissue', '379', '433'))), '../../../Analysis/diff_exp_all/tissue/results_379_vs_433.csv')

