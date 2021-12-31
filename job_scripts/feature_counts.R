library(Rsubread)

alignments <- c('../../../maps/SRR6238092_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238093_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238094_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238095_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238096_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238097_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238098_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238099_wbcel235.104_all/Aligned.out.bam','../../../maps/SRR6238100_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238101_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238102_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238103_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238104_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238105_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238106_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238107_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238108_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238109_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238110_wbcel235.104_all/Aligned.out.bam', '../../../maps/SRR6238111_wbcel235.104_all/Aligned.out.bam')

feature_counts <- featureCounts(files = alignments, annot.ext = '../../../annotations/wbcel235.104.gtf', isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, allowMultiOverlap = TRUE, isPairedEnd = TRUE, minFragLength = 10, nthreads = 40, reportReads = 'BAM', reportReadsPath = '../../../Analysis/feature_counts_all/reportReads/')

counts_matrix <- feature_counts[[1]]
counts_annot <- feature_counts[[2]]
counts_stats <- feature_counts[[4]]

write.csv(counts_matrix, '../../../Analysis/feature_counts_all/counts.csv')
write.csv(counts_annot, '../../../Analysis/feature_counts_all/annotations.csv')
write.csv(counts_stats, '../../../Analysis/feature_counts_all/stats.csv')
