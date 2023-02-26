raw_vcf <- snakemake@input[['RAW']]
umi_vcf <- snakemake@input[['UMI']]
out_file <- snakemake@output[['OUTF']]

require(vcfR)
require(tidyverse)
require(plotly)
require(htmlwidgets)
raw_vcf <- read.vcfR(raw_vcf, verbose = FALSE)
umi_vcf <- read.vcfR(umi_vcf, verbose = FALSE)

raw_vcf <- vcfR2tidy(raw_vcf) 
raw_vcf <- raw_vcf$gt %>% 
  select(ChromKey, POS, gt_FREQ) %>% 
  left_join(raw_vcf$fix) %>% 
  mutate(gt_FREQ = as.numeric(gsub("\\%", "", gt_FREQ)))

umi_vcf <- vcfR2tidy(umi_vcf) 
umi_vcf <- umi_vcf$gt %>% 
  select(ChromKey, POS, gt_FREQ) %>% 
  left_join(umi_vcf$fix) %>% 
  mutate(gt_FREQ = as.numeric(gsub("\\%", "", gt_FREQ)))


raw_vcf <- raw_vcf %>% 
  filter(CHROM == unique(umi_vcf$CHROM))

vcf <- umi_vcf %>% 
  full_join(raw_vcf, by = join_by(ChromKey, POS, CHROM, ID, REF, ALT))


vcf_plt <- vcf %>% 
  plot_ly(x = ~POS, y = ~gt_FREQ.x, color = ~ALT, type = "bar") %>% 
  add_trace(y = ~gt_FREQ.y*-1) %>% 
  layout(yaxis = list(title = "Raw variants        |  Percentage  |        UMI corrected"),
         xaxis = list(title = "Position in amplicon"),
         title = unique(umi_vcf$CHROM)
         ) %>% 
  layout(barmode = 'group', yaxis = list(range = c(-100, 100)))


saveWidget(as_widget(vcf_plt), paste0(out_file, ".html"), selfcontained = T)



