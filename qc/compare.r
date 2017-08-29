library(qrqc)

# 指定fastq文件路径
fq_files = c(rawdata="../raw_data/untreated.fq",
             trimmomatic="u.trimmomatic.fq",
             soapnuke="u.soapnuke.fq",
             sickle="u.sickle.fq",
             seqtk_trimfq="u.trimfq.fq")


# fq_files = c(rawdata="raw_data/untreated.fq",
#              soapnuke_l20_q0.5="u.soapnuke.fq",
#              soapnuke_l20_q0.2="u.soapnuke_q0.2.fq")

# Load each fastq file in by using qrqc's readSeqFile
seq_read = lapply(fq_files, function(file) {
  readSeqFile(file, hash = FALSE, kmer = FALSE)
})

quals = mapply(function(sfq, name) {
  qs = getQual(sfq)
  qs$trimmer = name
  qs
}, seq_read, names(fq_files), SIMPLIFY = FALSE)

d = do.call(rbind, quals)

# 画图
library(ggplot2)
position = d$position
mean = d$mean
trimmer = d$trimmer
p1 = qplot(position, mean, color=trimmer, geom=c("line", "point"))
p1 = p1 + ylab("Mean quality (Phred33)") + theme_bw()
print(p1)

p2 = qualPlot(seq_read, quartile.color = NULL, mean.color = NULL) + theme_bw()
p2 = p2 + scale_y_continuous("Quality (Phred33)")
print(p2)








