library(ggplot2)
library(readr)
library(optparse)

option_list = list(
  make_option("--intable", type="character", help="Tabulated file with abundances. Header: (sampling reference position coverage)"),
  make_option("--threshold", type="integer", help="The treshold used for the sampling"),
  make_option("--outpdf", type="character", help="The path used to write the pdf results")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


REFERENCES <-
  c(
    'EU286989.1.1503',
    'AM980830.1.1445',
    'JQMC01000035.331.1802',
    'JF219847.1.1310',
    'KM356606.1.1440',
    'JQ427863.1.1441',
    'AJ412686.1.1278',
    'KC286759.1.1477',
    'KF712645.1.1504',
    'AY221072.1.1451',
    'JN860392.1.1464',
    'AJ318187.1.1488',
    'GU118209.1.1413'
  )

print("Parse input table")
coverage <- read_delim(opt$intable, "\t", escape_double = FALSE, col_types = cols(sampling = col_factor(levels = c("before", "after"))), trim_ws = TRUE)

print("Plot results")
pdf(opt$outpdf)
for ( reference in REFERENCES ) {
  subset_cov <- coverage[coverage$reference == reference, ]
  gplot <- ggplot(subset_cov, aes(x=position, y=coverage, color=sampling))  + 
    geom_line() + 
    coord_cartesian(ylim=c(0, 10000)) + 
    geom_ribbon(aes(ymin=0, ymax=coverage, fill=sampling), alpha=0.4) + 
    geom_hline(aes(yintercept=opt$threshold), colour="#AA0000", linetype="dashed") + 
    facet_grid(reference ~ .)
  print(gplot)
}
dev.off()