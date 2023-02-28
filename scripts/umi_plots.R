raw_vcf <- snakemake@input[['RAW']]
umi_vcf <- snakemake@input[['UMI']]
out_file <- snakemake@output[['OUTF']]
# system('export HOME="/home"') #for pandoc

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
  full_join(raw_vcf, by = c("ChromKey", "POS", "CHROM", "ID", "REF", "ALT"))

vcf_plt <- vcf %>% 
  plot_ly(x = ~POS, y = ~gt_FREQ.x, color = ~ALT, type = "bar",
          hovertemplate = paste('<i>Percentage</i>: %{y}',
                                '<br><br>')) %>% 
  add_trace(y = ~gt_FREQ.y*-1) %>% 
  layout(yaxis = list(title = "Raw variants        |  Percentage  |        UMI corrected"),
         xaxis = list(title = "Position in amplicon"),
         title = unique(umi_vcf$CHROM)
         ) %>% 
  layout(barmode = 'group', yaxis = list(range = c(-100, 100))) %>% 
  layout(hovermode = "x unified")



vcf_lng <- vcf %>% 
  select(POS, gt_FREQ.x, gt_FREQ.y) %>% 
  pivot_longer(cols = c(gt_FREQ.x, gt_FREQ.y),
               values_to = "freq"
               )

vcf_lng <- vcf_lng %>% 
  group_by(name) %>% 
  mutate(gt_et_1 = sum(freq >= 1, na.rm = T)) %>% 
  mutate(gt_et_5 = sum(freq >= 5, na.rm = T)) %>% 
  mutate(gt_et_10 = sum(freq >= 10, na.rm = T))

vcf_lng <- vcf_lng %>% 
  ungroup() %>% 
  select(name, gt_et_1, gt_et_5, gt_et_10) %>% 
  unique()

vcf_lng$name[vcf_lng$name == "gt_FREQ.x"] <- "UMI"  
vcf_lng$name[vcf_lng$name == "gt_FREQ.y"] <- "RAW"


vcf_tbl <- plot_ly(
  type = 'table',
  header = list(
    values = c("≥1%","≥5%", "≥10%"),
    align = c("center", "center", "center"),
    line = list(width = 1, color = 'black'),
    fill = list(color = c("grey", "grey", "grey")),
    font = list(family = "Arial", size = 14, color = "white")
  ),
  cells = list(
    values = rbind(vcf_lng$gt_et_1, vcf_lng$gt_et_5, vcf_lng$gt_et_10),
    align = c("center", "center", "center"),
    line = list(color = "black", width = 1),
    font = list(family = "Arial", size = 12, color = c("black"))
  ))

vcf_plt <- subplot(vcf_plt, vcf_tbl,nrows = 1)

saveWidget(as_widget(vcf_plt), out_file, selfcontained = F)



