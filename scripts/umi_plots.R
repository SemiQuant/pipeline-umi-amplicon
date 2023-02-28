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



hline <- function(y = 0, color = "gray") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, dash="dot", width = 0.1)
  )
}


vcf_plt <- vcf_plt %>% 
  layout(shapes = list(hline(1),
                       hline(5),
                       hline(10),
                       hline(-1),
                       hline(-5),
                       hline(-10)
                       )) %>% 
  add_annotations(x = 0, 
                   y = c(1, 5, 10,
                         -1, -5, -10),
                   xref = "x",
                   yref = "y",
                   text = c(as.character(vcf_lng[vcf_lng$name == "UMI", "gt_et_1"]),
                            as.character(vcf_lng[vcf_lng$name == "UMI", "gt_et_5"]),
                            as.character(vcf_lng[vcf_lng$name == "UMI", "gt_et_10"]),
                            as.character(vcf_lng[vcf_lng$name == "RAW", "gt_et_1"]),
                            as.character(vcf_lng[vcf_lng$name == "RAW", "gt_et_5"]),
                            as.character(vcf_lng[vcf_lng$name == "RAW", "gt_et_10"])),
                   font = list(color = 'grey'), 
                   showarrow = F)


saveWidget(as_widget(vcf_plt), out_file, selfcontained = F)



