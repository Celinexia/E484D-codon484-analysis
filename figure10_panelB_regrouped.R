# Figure 10B, regrouped: sample contexts along the x axis, E484D allele fraction
# on the y axis. Five columns.
#
#   IN VIVO          clinical surveillance (incl. one patient mAb-treatment run)
#                    animal (wild deer and experimental animals; MA1 lung point)
#   EX VIVO          cell culture, incl. drug-selection passages (molnupiravir,
#                    remdesivir) and MA1 cell supernatant point
#                    wastewater
#   NOT RECORDED     unspecified
#
# Input: panelB_data_regrouped.csv (see the regrouping step in the notes).
# No methods text on the figure. Uses only ggplot2 + scales.
#
#   Rscript figure10_panelB_regrouped.R

suppressPackageStartupMessages({ library(ggplot2); library(scales) })

FONT   <- "Helvetica"      # metric-compatible with Arial
INK    <- "#000000"
INK2   <- "#52514e"
MUTED  <- "#898781"
RULE   <- "#b9b8b0"
FAINT  <- "#e8e7e1"
BLUE   <- "#2a78d6"        # 23014 A>C
ORANGE <- "#eb6834"        # 23014 A>T

base_theme <- theme_minimal(base_size = 6, base_family = FONT) +
  theme(
    text               = element_text(colour = INK, family = FONT),
    axis.text.x        = element_text(colour = "#000000", face = "bold", size = 6.0,
                                      lineheight = 0.95, margin = margin(t = 2)),
    axis.text.y        = element_text(colour = "#000000", face = "bold", size = 6.5),
    axis.title.y       = element_text(colour = INK, face = "bold", size = 7.5,
                                      margin = margin(r = 3)),
    axis.title.x       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_rect(fill = "white", colour = NA),
    axis.line          = element_blank(),
    axis.ticks         = element_line(colour = "#333333", linewidth = 0.35),
    axis.ticks.length  = unit(3, "pt"),
    legend.position    = "none",
    plot.margin        = margin(4, 6, 2, 4)
  )

d <- read.csv("panelB_data_regrouped.csv", stringsAsFactors = FALSE)

# ---- five columns; each context maps to one column
COL <- c(clinical = 1, animal = 1.85, culture = 2.7, wastewater = 3.55, unspec = 4.4)
ctx_to_col <- c("Clinical surveillance" = unname(COL["clinical"]),
                "Wild animal" = unname(COL["animal"]), "Experimental animal" = unname(COL["animal"]),
                "Cell culture" = unname(COL["culture"]),
                "Drug selection (cell culture)" = unname(COL["culture"]),
                "Wastewater" = unname(COL["wastewater"]),
                "Unspecified" = unname(COL["unspec"]))
d$x <- unname(ctx_to_col[d$context])
pts <- d[!is.na(d$x), ]
set.seed(3)
pts$xj <- pts$x + runif(nrow(pts), -0.24, 0.24)

# MA1 points, offset to the right of their column
MA1_OFF <- 0.30
ma1 <- data.frame(x = c(unname(COL["animal"]) + 0.2, unname(COL["culture"]) + MA1_OFF),
                  vaf = c(0.96, 0.92))

cols <- data.frame(
  x     = unname(COL),
  label = c("Clinical\nsurveillance", "Animal", "Cell culture\n& drug\nselection",
            "Wastewater", "Unspecified"),
  n     = c(227, 17, 39, 14, 111),   # totals include the MA1 points
  stringsAsFactors = FALSE)
cols$label2 <- paste0(cols$label, "\n(n = ", cols$n, ")")

med <- aggregate(vaf ~ x, data = pts, FUN = median)

# group brackets, spanning the columns of each group
grp <- data.frame(
  group = c("IN VIVO", "EX VIVO", "NOT RECORDED"),
  x0    = c(0.6, 2.45, 4.15),
  x1    = c(2.1, 3.85, 4.65),
  stringsAsFactors = FALSE)
GROUP_Y <- 1.06
LEG_Y1  <- 1.25

pB <- ggplot() +
  geom_segment(data = grp, aes(x = x0, xend = x1, y = GROUP_Y, yend = GROUP_Y),
               colour = RULE, linewidth = 0.4) +
  annotate("segment", x = 0.3, xend = 0.3, y = -0.02, yend = 1.0,
           colour = "#333333", linewidth = 0.35) +
  annotate("segment", x = 0.3, xend = 4.9, y = -0.02, yend = -0.02,
           colour = "#333333", linewidth = 0.35) +
  geom_text(data = grp, aes(x = (x0 + x1) / 2, y = GROUP_Y + 0.03, label = group),
            size = 2.2, family = FONT, colour = "#000000", fontface = "bold", vjust = 0) +
  geom_point(data = pts, aes(x = xj, y = vaf, colour = route),
             size = 1.8, alpha = 1, stroke = 0) +
  geom_point(data = ma1, aes(x = x, y = vaf), shape = 23, size = 2.2,
             fill = "white", colour = INK, stroke = 0.4) +
  geom_segment(data = med, aes(x = x - 0.3, xend = x + 0.3, y = vaf, yend = vaf),
               colour = INK, linewidth = 0.5) +
  scale_colour_manual(values = c(`A>C` = BLUE, `A>T` = ORANGE)) +
  scale_x_continuous(breaks = cols$x, labels = cols$label2,
                     limits = c(0.3, 5.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.02, 1.30), breaks = seq(0, 1, 0.25),
                     labels = c("0", "0.25", "0.50", "0.75", "1.0"),
                     expand = c(0, 0)) +
  labs(y = "Fraction of reads carrying E484D") +
  base_theme

# one-row key, evenly spaced across the plot width
pB <- pB +
  annotate("point", x = 0.5, y = LEG_Y1, colour = BLUE, size = 1.6) +
  annotate("text",  x = 0.7, y = LEG_Y1, label = "23014 A>C", hjust = 0,
           size = 2.2, family = FONT, colour = "#000000", fontface = "bold") +
  annotate("point", x = 1.7, y = LEG_Y1, colour = ORANGE, size = 1.6) +
  annotate("text",  x = 1.9, y = LEG_Y1, label = "23014 A>T", hjust = 0,
           size = 2.2, family = FONT, colour = "#000000", fontface = "bold") +
  annotate("point", x = 2.9, y = LEG_Y1, shape = 23, size = 2.0,
           fill = "white", colour = INK, stroke = 0.4) +
  annotate("text",  x = 3.1, y = LEG_Y1, label = "MA1 lineage", hjust = 0,
           size = 2.2, family = FONT, colour = "#000000", fontface = "bold") +
  annotate("segment", x = 4.1, xend = 4.4, y = LEG_Y1, yend = LEG_Y1,
           colour = INK, linewidth = 0.6) +
  annotate("text",  x = 4.5, y = LEG_Y1, label = "median", hjust = 0,
           size = 2.2, family = FONT, colour = "#000000", fontface = "bold")

W <- 3.9
H <- 4.2
ggsave("Figure10-panelB-regrouped.pdf", pB, width = W, height = H, device = "pdf", bg = "white")
ggsave("Figure10-panelB-regrouped.png", pB, width = W, height = H, dpi = 300, bg = "white")
cat("wrote Figure10-panelB-regrouped.pdf\n")
