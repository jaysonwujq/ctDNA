library(cghFLasso)
#' stats.of.ct
#'
#' get list of mean, standard deviation, median, pk
#'
#' @param ct normCT values
#' @param cl parallel cluster used in package 'parallel'
#' @param n.sigma     n sigma of normal distribution, eg. 1.96 means 95.45,
#'                     2.575829 means 99(the default)
#
getProgramName<-function(arguments){
         args <- commandArgs(trailingOnly = FALSE)
         sub("--file=", "", args[grep("--file=", args)])
}
program <- getProgramName()

args <- commandArgs(trailingOnly = TRUE)
indir <- args[3]
outdir<-args[2]
stats.of.ct <- function(ct, n.sigma = 2.575829) {
t(apply(as.matrix(ct),1, function(x) {
stats.row.of.ct <- function(x, n = 2.575829) {
u = !is.na(x)
xu = x[u]
if(sum(u)>1) {
m   = mean(xu)
std = sd(xu)
md  = median(xu)
should.trim = xu > m + n*std | xu < m - n*std
if ( sum(should.trim) > 0) {
x = xu[!should.trim]
stats.row.of.ct(x, n)
} else {
re = c(m,std,md)
re
}
} else {
m   = NA
std = NA
md  = NA
re = c(m,std,md)
re
}
}
return(stats.row.of.ct(x, n = n.sigma))
}));
}

ct = read.table(args[1], header = T)
samples <- names(ct)
samples <- substr(samples, 1 + as.numeric(grepl("^X", samples)),100)
samples <- gsub("\\.", "-", samples)
names(ct) <- samples

msdct <- stats.of.ct(ct)
zct <- (ct - msdct[,1]) / msdct[,2]       # zvalue of each row
fct <- ct / msdct[,3]                     # x/meanvalue foreach row

for (sample in samples) {
#  sample = samples[sample]
  data = data.frame(ct = ct[, sample], fct = fct[, sample], zct = zct[, sample])

  # trim NA values from tail
  data.isna = rle(is.na(data$ct))
  data.length = length(data$ct)
  data.isna.length = length(data.isna$values)
  if (data.isna$values[data.isna.length]) {
    data = data[seq(1, data.length - data.isna$lengths[data.isna.length]),]
  }

  data$chromosome = 1
  data$position = as.numeric(rownames(data))
  data.length = length(data$position)
  data$log2 <- log2(data$fct)  # log2 values
  any.na <- function(x) { any(is.na(x)) }
  data$is.na = apply(data, MARGIN = 1, any.na)
  data$not.na = !data$is.na    # not.na is called valid

  flasso.log2 = cghFLasso(data$log2, nucleotide.position=data$position, FDR = 0.01)

  data$copyn.log2 = flasso.log2$Esti.CopyN
  write.table(data$copyn.log2,paste(outdir,'/',sample,'.',args[1],sep=''))
  png(file=paste(outdir,'/',sample,'.',args[1],'_ct_cnv.png',sep=''));
  plot( data$copyn.log2)
  dev.off()
}
