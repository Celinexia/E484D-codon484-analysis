# Figure 10, panels A and B -- ggplot2 versions for publication.
#
#   A  amino-acid substitutions reachable from the parental codon 484 (GAA) in
#      a single nucleotide step, against their frequency among GenBank genomes.
#   B  allele fraction of E484D in every setting where it has been observed.
#
# Panel C (brain / lung E484A, Richner group) is supplied separately.
# No methods text on the figure: description belongs in the legend and Methods.
# Uses only ggplot2 + scales + grid, so it runs on a stock R install.
#
#   Rscript figure10_panelsAB.R

suppressPackageStartupMessages({
  library(ggplot2); library(scales); library(grid)
})

FONT   <- "Arial"
INK    <- "#0b0b0b"
INK2   <- "#52514e"
MUTED  <- "#898781"
RULE   <- "#b9b8b0"
FAINT  <- "#e8e7e1"
BLUE   <- "#2a78d6"        # 23014 A>C  -- these two hues mean the nucleotide
ORANGE <- "#eb6834"        # 23014 A>T     route at 23014 and nothing else
TEAL   <- "#0F9D9A"        # E484D focus, as in the Figure 10A legend

base_theme <- theme_minimal(base_size = 7, base_family = FONT) +
  theme(
    text              = element_text(colour = INK, family = FONT),
    axis.text.x       = element_text(colour = INK, size = 6.6, family = FONT,
                                     face = "bold"),
    axis.text.y       = element_blank(),
    axis.title.x      = element_text(colour = INK, size = 7.2, family = FONT,
                                     face = "bold", margin = margin(t = 5)),
    axis.title.y      = element_blank(),
    panel.grid        = element_blank(),
    panel.background  = element_blank(),
    plot.background   = element_rect(fill = "white", colour = NA),
    axis.line         = element_blank(),
    axis.ticks.x      = element_line(colour = MUTED, linewidth = 0.4),
    axis.ticks.y      = element_blank(),
    axis.ticks.length = unit(1.6, "pt"),
    legend.position   = "none",
    plot.margin       = margin(2, 4, 1, 3)
  )

# =============================================================== panel A =====
N_GENOMES <- 9438404
PARENTAL  <- "GAA"

A <- data.frame(
  aa    = c("Parental", "E484A", "E484K", "E484Q", "E484D", "E484G", "E484V"),
  nt    = c("E484 reference", "23013 A>C", "23012 G>A", "23012 G>C",
            "23014 A>C / A>T", "23013 A>G", "23013 A>T"),
  codon = c(PARENTAL, "GCA", "AAA", "CAA", "GAC / GAT", "GGA", "GTA"),
  n     = c(NA, 4306700, 370685, 8237, 229, 170, NA),
  stringsAsFactors = FALSE)

A$parental <- A$aa == "Parental"
A$lab_aa   <- ifelse(A$parental, "Parental (E484)", A$aa)

# Lineage colour code. Blue and orange are reserved for the 23014 route, so
# they do not appear here. Focus (E484D) is teal, as in the legend.
LIN_COL <- c(
  "Omicron"          = "#7B3294",
  "Beta/Gamma/Mu"    = "#1B7837",
  "Kappa"            = "#C51B7D",
  "no lineage"       = "#898781",
  "reference"        = "#0b0b0b",
  "not observed"     = "#898781")
A$lineage <- c(Parental = "reference", E484A = "Omicron",
               E484K = "Beta/Gamma/Mu", E484Q = "Kappa", E484D = "no lineage",
               E484G = "no lineage", E484V = "not observed")[A$aa]
A$fillc   <- ifelse(A$aa == "E484D", TEAL, unname(LIN_COL[A$lineage]))
A$freq     <- A$n / N_GENOMES
A$focus    <- A$aa == "E484D"
A$ypos     <- rev(seq_len(nrow(A)))
A$lab      <- ifelse(A$parental, "",
              ifelse(is.na(A$freq), "not observed",
              ifelse(A$freq >= 0.01, sprintf("%.1f%%",  100 * A$freq),
              ifelse(A$freq >= 1e-4, sprintf("%.3f%%",  100 * A$freq),
                                     sprintf("%.4f%%",  100 * A$freq)))))
yP <- A$ypos[A$parental]

FLOOR <- 5e-6                 # left edge of the data axis
XMIN  <- 10^-7.0              # left edge of the reserved label zone
X_AA  <- 10^-7.0              # substitution name
X_COD <- 10^-5.95             # codon
X_NT  <- 10^-6.1              # nucleotide change
YLO   <- 0.42
YTOP  <- nrow(A) + 0.55
yD    <- A$ypos[A$aa == "E484D"]

pA <- ggplot(A) +
  annotate("rect", xmin = XMIN, xmax = 3, ymin = yD - 0.45, ymax = yD + 0.45,
           fill = FAINT, colour = NA) +
  annotate("segment", x = FLOOR, xend = 1, y = YLO, yend = YLO,
           colour = MUTED, linewidth = 0.4) +
  annotate("segment", x = FLOOR, xend = 1, y = yP, yend = yP,
           colour = RULE, linewidth = 0.45) +
  geom_segment(data = subset(A, !parental),
               aes(x = FLOOR, xend = freq, y = ypos, yend = ypos, colour = fillc),
               linewidth = 0.5, alpha = 0.7, na.rm = TRUE) +
  geom_point(data = subset(A, !parental),
             aes(x = freq, y = ypos, fill = fillc),
             size = 1.7, shape = 21, colour = INK, stroke = 0.45, na.rm = TRUE) +
  # lineage key, bottom right of the data region
  annotate("point", x = 0.005, y = 1.90, shape = 21, size = 1.6,
           fill = LIN_COL[["Omicron"]], colour = INK, stroke = 0.4) +
  annotate("text", x = 0.0068, y = 1.90, label = "Omicron (E484A)", hjust = 0,
           size = 2.2, family = FONT, colour = INK2, fontface = "bold") +
  annotate("point", x = 0.06, y = 1.90, shape = 21, size = 1.6,
           fill = LIN_COL[["Beta/Gamma/Mu"]], colour = INK, stroke = 0.4) +
  annotate("text", x = 0.074, y = 1.90, label = "Beta/Gamma/Mu (E484K)", hjust = 0,
           size = 2.2, family = FONT, colour = INK2, fontface = "bold") +
  annotate("point", x = 0.005, y = 1.45, shape = 21, size = 1.6,
           fill = LIN_COL[["Kappa"]], colour = INK, stroke = 0.4) +
  annotate("text", x = 0.0068, y = 1.45, label = "Kappa (E484Q)", hjust = 0,
           size = 2.2, family = FONT, colour = INK2, fontface = "bold") +
  annotate("point", x = 0.06, y = 1.45, shape = 21, size = 1.6,
           fill = LIN_COL[["no lineage"]], colour = INK, stroke = 0.4) +
  annotate("text", x = 0.074, y = 1.45, label = "no lineage (E484G)", hjust = 0,
           size = 2.2, family = FONT, colour = INK2, fontface = "bold") +
  annotate("point", x = 0.005, y = 1.00, shape = 21, size = 1.6,
           fill = INK, colour = INK, stroke = 0.4) +
  annotate("text", x = 0.0068, y = 1.00, label = "E484D (no lineage)", hjust = 0,
           size = 2.2, family = FONT, colour = INK, fontface = "bold") +
  geom_point(data = subset(A, is.na(freq) & !parental), aes(x = FLOOR, y = ypos),
             shape = 4, size = 1.4, colour = MUTED, stroke = 0.5) +
  geom_text(aes(x = ifelse(is.na(freq), FLOOR * 2.2, freq * 1.6),
                y = ypos, label = lab),
            hjust = 0, size = 2.2, family = FONT, colour = INK2) +
  geom_text(aes(x = X_AA, y = ypos, label = lab_aa),
            fontface = "bold", hjust = 0, size = 2.6, family = FONT, colour = INK) +
  geom_text(aes(x = X_COD, y = ypos, label = codon),
            fontface = "bold",
            hjust = 0, size = 2.2, family = FONT, colour = INK) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_size_manual(values = c(`FALSE` = 1.45, `TRUE` = 2.15)) +
  scale_x_log10(limits = c(XMIN, 3),
                breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
                labels = c("0.001%", "0.01%", "0.1%", "1%", "10%", "100%"),
                expand = c(0, 0)) +
  scale_y_continuous(limits = c(YLO, YTOP), expand = c(0, 0)) +
  labs(x = "Frequency among SARS-CoV-2 genomes") +
  base_theme

# =============================================================== panel B =====
d <- read.csv("panelB_data.csv", stringsAsFactors = FALSE)

ma1 <- data.frame(
  row = c("Mouse lung, MA1 adaptation", "Cell supernatant, MA1"),
  vaf = c(0.96, 0.92), stringsAsFactors = FALSE)

rows <- data.frame(
  kind = c("group","sra","sra","sra","ma1","sra",
           "group","sra","ma1","sra","group","sra"),
  row  = c("IN VIVO","Clinical surveillance","Wild white-tailed deer",
           "Experimental animal","Mouse lung, MA1 adaptation","Treatment-selection",
           "EX VIVO","Cell culture","Cell supernatant, MA1","Wastewater",
           "CONTEXT NOT RECORDED","Unspecified"),
  ctx  = c(NA,"Clinical surveillance","Wild animal","Experimental animal",NA,
           "Treatment-selection",NA,"Cell culture",NA,"Wastewater",NA,"Unspecified"),
  stringsAsFactors = FALSE)
rows$y <- rev(seq_len(nrow(rows)))

pts <- merge(d, rows[rows$kind == "sra", c("row","ctx","y")],
             by.x = "context", by.y = "ctx")
mp  <- merge(ma1, rows[rows$kind == "ma1", c("row","y")], by = "row")

set.seed(3)
pts$yj <- pts$y + runif(nrow(pts), -0.23, 0.23)

med <- aggregate(vaf ~ y, data = pts, FUN = median)
cnt <- aggregate(vaf ~ y, data = pts, FUN = length); names(cnt)[2] <- "n"
med <- merge(med, cnt, by = "y")

lab_df <- rows[rows$kind != "group", ]
grp_df <- rows[rows$kind == "group", ]

BXMIN <- -0.42
X_ROW <- -0.41
X_N   <-  1.02
BYLO  <-  0.30
LEG_Y <- nrow(rows) + 0.75

pB <- ggplot() +
  geom_segment(data = lab_df, aes(x = 0, xend = 1, y = y, yend = y),
               colour = FAINT, linewidth = 0.3) +
  geom_segment(data = grp_df, aes(x = X_ROW, xend = 1, y = y, yend = y),
               colour = RULE, linewidth = 0.4) +
  geom_text(data = grp_df, aes(x = X_ROW, y = y + 0.20, label = row),
            hjust = 0, size = 1.95, family = FONT, colour = INK, fontface = "bold") +
  annotate("segment", x = 0, xend = 1, y = BYLO, yend = BYLO,
           colour = MUTED, linewidth = 0.4) +
  geom_point(data = pts, aes(x = vaf, y = yj, colour = route),
             size = 0.72, alpha = 0.8, stroke = 0) +
  geom_point(data = mp, aes(x = vaf, y = y), shape = 23, size = 1.45,
             fill = "white", colour = INK, stroke = 0.45) +
  geom_segment(data = med, aes(x = vaf, xend = vaf, y = y - 0.3, yend = y + 0.3),
               colour = INK, linewidth = 0.55) +
  geom_text(data = lab_df, aes(x = X_ROW, y = y, label = row),
            hjust = 0, size = 2.5, family = FONT, colour = INK, fontface = "bold") +
  geom_text(data = med, aes(x = X_N, y = y, label = n),
            hjust = 0, size = 2.7, family = FONT, colour = INK, fontface = "bold") +
  geom_text(data = data.frame(y = mp$y), aes(x = X_N, y = y, label = "1"),
            hjust = 0, size = 2.7, family = FONT, colour = INK, fontface = "bold") +
  scale_colour_manual(values = c(`A>C` = BLUE, `A>T` = ORANGE)) +
  scale_x_continuous(limits = c(BXMIN, 1.09), breaks = seq(0, 1, 0.25),
                     labels = c("0", "0.25", "0.50", "0.75", "1.0"),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(BYLO, LEG_Y + 0.35), expand = c(0, 0)) +
  labs(x = "Fraction of reads carrying E484D") +
  base_theme + theme(plot.margin = margin(2, 12, 1, 3))

pB <- pB +
  annotate("point", x = 0.00, y = LEG_Y, colour = BLUE, size = 0.95) +
  annotate("text",  x = 0.02, y = LEG_Y, label = "23014 A>C",
           hjust = 0, size = 1.95, family = FONT, colour = INK2) +
  annotate("point", x = 0.17, y = LEG_Y, colour = ORANGE, size = 0.95) +
  annotate("text",  x = 0.19, y = LEG_Y, label = "23014 A>T",
           hjust = 0, size = 1.95, family = FONT, colour = INK2) +
  annotate("point", x = 0.34, y = LEG_Y, shape = 23, size = 1.35,
           fill = "white", colour = INK, stroke = 0.45) +
  annotate("text",  x = 0.36, y = LEG_Y, label = "MA1 lineage",
           hjust = 0, size = 1.95, family = FONT, colour = INK2) +
  annotate("segment", x = 0.52, xend = 0.52, y = LEG_Y - 0.12, yend = LEG_Y + 0.12,
           colour = INK, linewidth = 0.55) +
  annotate("text",  x = 0.54, y = LEG_Y, label = "median",
           hjust = 0, size = 1.95, family = FONT, colour = INK2)

# =============================================================== output ======
W  <- 7.09                      # 180 mm, two-column
HA <- 1.85
HB <- 2.60

# Arial for PDF: quartz embeds the system Arial. Falls back to pdf() if absent.
pdf_arial <- function(file, plot, width, height) {
  ok <- tryCatch({
    grDevices::quartz(file = file, type = "pdf", width = width, height = height)
    print(plot); grDevices::dev.off(); TRUE
  }, error = function(e) FALSE)
  if (!ok) ggsave(file, plot, width = width, height = height, device = "pdf")
}
pdf_arial("Figure10-panelA-color.pdf", pA, W, HA)
ggsave("Figure10-panelA-color.png", pA, width = W, height = HA, dpi = 600, bg = "white")

cat(sprintf("panel A: %d substitutions; E484D = %.4f%% of %s genomes\n",
            nrow(A), 100 * A$freq[A$aa == "E484D"], format(N_GENOMES, big.mark = ",")))
cat(sprintf("panel B: %d SRA runs, %d contexts, + 2 MA1 points\n",
            nrow(pts), length(unique(pts$context))))
