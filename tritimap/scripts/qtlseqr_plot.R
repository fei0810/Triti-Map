
###############################################################################
#
# Author contact:
# zhaofei920810@gmail.com
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
###############################################################################

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

if (!requireNamespace(c("QTLseqr"), quietly = TRUE)) {
  install.packages("scripts/QTLseqr-master/", repos = NULL, type = "source")
}

library(tidyr)
library(dplyr)
library(ggplot2)
library(QTLseqr)
library(readr)

qtldat <- snakemake@input[["qtlseqr"]]
snpraw <- snakemake@input[["snpindex"]]

pool1 <- snakemake@params[["bulk1"]]
pool2 <- snakemake@params[["bulk2"]]
pop <- snakemake@params[["pop_struc"]]
bulksize <- as.numeric(snakemake@params[["bulksize"]])
winsize <- as.numeric(snakemake@params[["winsize"]])

filter_percentage <- as.numeric(snakemake@params[["filterpercentage"]])
f_pvalue <- as.numeric(snakemake@params[["pvalue"]])
min_length <- as.numeric(snakemake@params[["min_length"]])

output_dir <- snakemake@params[["dirname"]]

qtl_out <- snakemake@output[["qtlout"]]
qtl_regionout <- snakemake@output[["qtlregion"]]
qtl_rawregion <- snakemake@output[["qtlrawregion"]]
qtl_snpinfo <- snakemake@output[["qtlsnpinfo"]]
qtl_countp <- snakemake@output[["qtlcount_p"]]
qtl_countl <- snakemake@output[["qtlcount_l"]]
qtl_indexp <- snakemake@output[["qtlindex_p"]]
qtl_indexl <- snakemake@output[["qtlindex_l"]]

### run qtlseqr

snpdat <- read_tsv(snpraw)

dat <- importFromTable(qtldat,
  highBulk = pool1,
  lowBulk = pool2, sep = "\t"
)
datqtl <- runQTLseqAnalysis(
  SNPset = dat,
  windowSize = winsize,
  popStruc = pop,
  bulkSize = c(bulksize, bulksize),
  replications = 10000,
  filter = 0.3,
  intervals = 99,
  maxk = 10000
)

### fisher test
snpfishertest <- function(dat) {
  dat$fisher_p <- apply(
    dat[, c(3, 4, 7, 8)], 1,
    function(x) {
      fisher.test(
        matrix(
          data = c(x[1], x[3], x[2], x[4]),
          nrow = 2, ncol = 2, byrow = TRUE,
          dimnames = NULL
        )
      )$p.value
    }
  )
  return(dat)
}

datqtl <- snpfishertest(datqtl)

### get candidate region

names(datqtl) <- gsub("LOW", pool2, names(datqtl))
names(datqtl) <- gsub("HIGH", pool1, names(datqtl))


rawout <- getQTLTable(SNPset = datqtl, method = "QTLseq", interval = 99, export = F)

# filter raw region by snp count and region length, 10 is suitable for EMS data?
out <- rawout %>% filter(nSNPs >= 10 & avgSNPs_Mb >= 10 &
  length >= min_length)

# if no region left
if (nrow(out) < 1) {
  cat("Error! No candidate region find! \n")
  stop("No candidate region find!", call. = F)
  quit(status = 1)
}

# sort by abs index
out <- out[order(abs(out$avgDeltaSNP), decreasing = T), ]

# snp conunt > xx% and snp index > xx%
filter_out <- out[abs(out$avgDeltaSNP) >= quantile(abs(out$avgDeltaSNP), probs = filter_percentage) & out$avgSNPs_Mb >= quantile(out$avgSNPs_Mb, probs = filter_percentage), ]

# if no region left
if (nrow(filter_out) < 1) {
  cat("Error! The filtration standard is too strict. Please lower the filtration standard. \n")
  stop("The filtration standard is too strict. Please lower the filtration standard.", call. = F)
  quit(status = 1)
}

print("qtlseq is down")

### get snp in candidate region

snpfilter <- data.frame(matrix(ncol = length(datqtl), nrow = 0))
names(snpfilter) <- names(datqtl)
for (i in seq_len(nrow(filter_out))) {
  snptemp <- datqtl[datqtl$CHROM == filter_out$CHROM[i] & datqtl$POS >= filter_out$start[i] & datqtl$POS <= filter_out$end[i] & datqtl$fisher_p < f_pvalue, ]
  snpfilter <- rbind(snpfilter, snptemp)
}

print("snpinfo is down")

snpfilter$CHROM <- as.character(snpfilter$CHROM)
snpfilter$POS <- as.double(snpfilter$POS)
snpdat$`#CHROM` <- as.character(snpdat$`#CHROM`)
snpdat$POS <- as.double(snpdat$POS)

snpmerge <- left_join(snpfilter, snpdat, by = c("CHROM" = "#CHROM", "POS" = "POS")) %>% select(c(1, 2, 20, 21, 3:18))

write.table(datqtl,
  file = qtl_out,
  quote = F, sep = "\t", row.names = F, col.names = T
)

write.table(out,
  file = qtl_rawregion,
  quote = F, sep = "\t", row.names = F, col.names = T
)

write.table(filter_out,
  file = qtl_regionout,
  quote = F, sep = "\t", row.names = F, col.names = T
)

write.table(snpmerge,
  file = qtl_snpinfo,
  quote = F, sep = "\t", row.names = F, col.names = T
)

### plot all figure

qtlplot <- function(rawdat, var = "deltaSNP", type = "point", chroms = NULL) {
  rawdat <-
    if (is.null(chroms)) {
      rawdat
    } else {
      rawdat %>% filter(CHROM == chroms)
    }

  rawdat <- rawdat %>%
    mutate(candidate = ifelse(
      abs(tricubeDeltaSNP) > abs(CI_99),
      "pass",
      "fail"
    ))

  p <- ggplot(
    data = rawdat
  ) +
    scale_x_continuous(
      breaks = seq(
        0, max(rawdat$POS),
        10^(floor(log10(max(rawdat$POS))))
      ),
      labels = function(x) {
        format(round(x / 1e6, 1),
          trim = T, scientific = F
        )
      },
      name = "Genomic Position (Mb)"
    )

  if (var == "deltaSNP") {
    if (type == "point") {
      p <- p + ggplot2::geom_point(aes(POS, tricubeDeltaSNP, color = candidate), size = 0.5, alpha = 0.8) +
        ggplot2::geom_line(ggplot2::aes(POS, CI_99), color = "blue", alpha = 0.7) +
        ggplot2::geom_line(ggplot2::aes(POS, -CI_99), color = "blue", alpha = 0.7) +
        labs(
          title = bquote(.(pool1) ~ .(pool2) ~ Delta * "SNPindex"),
          y = expression(Delta * "SNPindex"),
          x = "Genomic Position (Mb)"
        ) + scale_color_manual(
          values = c(
            pass = "black",
            fail = "grey"
          )
        )
    }
    if (type == "line") {
      p <- p + ggplot2::geom_line(aes(POS, tricubeDeltaSNP), alpha = 0.8) +
        ggplot2::geom_line(ggplot2::aes(POS, CI_99), color = "blue", alpha = 0.7) +
        ggplot2::geom_line(ggplot2::aes(POS, -CI_99), color = "blue", alpha = 0.7) +
        labs(
          title = bquote(.(pool1) ~ .(pool2) ~ Delta * "SNPindex"),
          y = expression(Delta * "SNPindex"),
          x = "Genomic Position (Mb)"
        )
    }
  }
  if (var == "nSNPs") {
    if (type == "line") {
      p <- p + ggplot2::geom_line(aes(POS, nSNPs), alpha = 0.8) +
        labs(
          title = paste(pool1, pool2, "SNP Counts"),
          y = "SNP Counts",
          x = "Genomic Position (Mb)"
        )
    }
    if (type == "point") {
      p <- p + ggplot2::geom_point(aes(POS, nSNPs), size = 0.5, alpha = 0.8) +
        labs(
          title = paste(pool1, pool2, "SNP Counts"),
          y = "SNP Counts",
          x = "Genomic Position (Mb)"
        )
    }
  }
  p <- p +
    theme_bw() +
    facet_wrap(. ~ CHROM, nrow = 4) +
    ggplot2::theme(
      axis.text.x = element_text(angle = 90)
    )
  p
}

countp_pointp <- qtlplot(datqtl, var = "nSNPs", type = "point")
countp_linep <- qtlplot(datqtl, var = "nSNPs", type = "line")
index_pointp <- qtlplot(datqtl, var = "deltaSNP", type = "point")
index_linep <- qtlplot(datqtl, var = "deltaSNP", type = "line")

suppressMessages(
  ggsave(qtl_countp, countp_pointp)
)

suppressMessages(
  ggsave(qtl_countl, countp_linep)
)

suppressMessages(
  ggsave(qtl_indexp, index_pointp)
)

suppressMessages(
  ggsave(qtl_indexl, index_linep)
)
### plot candidate region

if (length(filter_out$CHROM) > 0) {
  for (i in 1:length(filter_out$CHROM)) {
    chrom <- as.character(filter_out$CHROM[i])
    start <- filter_out$start[i]
    end <- filter_out$end[i]
    zoomp <- qtlplot(datqtl, var = "deltaSNP", type = "point", chroms = chrom) +
      scale_x_continuous(
        limits = c(start - 3 * 1e6, end + 3 * 1e6),
        labels = function(x) {
          format(round(x / 1e6, 1),
            trim = T, scientific = F
          )
        },
        name = "Genomic Position (Mb)"
      ) +
      geom_line(aes(POS, tricubeDeltaSNP)) +
      geom_point(aes(POS, tricubeDeltaSNP)) +
      annotate("rect",
        xmin = c(start), xmax = c(end),
        ymin = -Inf, ymax = Inf,
        alpha = 0.2, fill = c("#b4aee8")
      ) +
      labs(
        title = bquote(.(pool1) ~ "vs" ~ .(pool2) ~
        .(chrom) ~ ":" ~ .(start) ~ "-" ~ .(end)),
        y = expression(Delta * "SNPindex"),
        x = "Genomic Position (Mb)"
      ) +
      theme(legend.position = "none")
    ggsave(paste0(output_dir, "/06_regionout/", pool1, "_vs_", pool2, "_candidateregion_", i, ".pdf"), plot = zoomp)
  }
}