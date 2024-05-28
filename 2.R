args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
reference <- args[2]

library(data.table)
library(Biostrings)

# 读取基因组文件
if (reference == "hg19") {
  genomefile <- "./hg19/hg19.fasta"
} else {
  stop("unsupported reference version")
}

seqDb <- readDNAStringSet(genomefile, format="fasta")

print(paste("reading genome file", reference))
print(paste("seqDb info:", length(seqDb), "sequences"))
print(paste("id info:", names(seqDb)))

outdir <- "./Fasta/"
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

name <- sub(".bed$", "", basename(input))
print(paste("name:", name))

data <- fread(input, header = FALSE)

output_file <- file.path(outdir, paste0(name, ".fasta"))
fasta_conn <- file(output_file, "w")

for (i in 1:nrow(data)) {
  row <- data[i]
  chr <- as.character(row$V1)
  start <- as.integer(row$V2)
  end <- as.integer(row$V3)
  if (chr == "chrY") next
  
  seqname <- paste(row$V1, row$V2, row$V3, row$V5, sep = "-")
  writeLines(paste0(">", seqname), con = fasta_conn)
  
  subseq <- as.character(subseq(seqDb[[chr]], start=start, end=end))
  writeLines(subseq, con = fasta_conn)
}

close(fasta_conn)

# 调用外部的 Perl 脚本
system(paste("perl 3-fimo.pl", name))
