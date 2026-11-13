# =====================================================================
# DTS 204 - Statistical Computing, Inference and Modelling
# Lab Assessment: Palmer Penguins Analysis
# =====================================================================
# Required packages: readxl, ggplot2  (palmerpenguins optional)
# install.packages(c("readxl", "ggplot2", "palmerpenguins"))

library(ggplot2)

# Folder for figures
dir.create("figures", showWarnings = FALSE)

# =====================================================================
# TASK 1: DATA EXPLORATION AND PREPROCESSING
# =====================================================================

# 1. Load the dataset (from the supplied Excel file; falls back to the
#    palmerpenguins package if the file is not found)
if (file.exists("DTS_204_Dataset.xlsx")) {
  library(readxl)
  penguins <- as.data.frame(read_excel("DTS_204_Dataset.xlsx"))
  penguins$species <- factor(penguins$species)
  penguins$island  <- factor(penguins$island)
  penguins$sex     <- factor(penguins$sex)
} else {
  library(palmerpenguins)
  penguins <- as.data.frame(palmerpenguins::penguins)
}

# 2. First six rows
cat("\n--- First six rows ---\n")
print(head(penguins, 6))
str(penguins)

# 3. Missing values
cat("\n--- Missing values per column ---\n")
print(colSums(is.na(penguins)))
cat("Rows with at least one missing value:", sum(!complete.cases(penguins)), "\n")

# Handling: the 2 rows missing all four measurements are unusable for
# modelling, so they are dropped. 9 further rows only lack 'sex'; they are
# kept for the measurement-based analyses and excluded only from the sex
# proportion test. (Complete-case analysis; imputation is not justified
# because <3.5% of rows are affected.)
num_vars <- c("bill_length_mm", "bill_depth_mm",
              "flipper_length_mm", "body_mass_g")
penguins_clean <- penguins[complete.cases(penguins[, num_vars]), ]
cat("Rows after removing missing measurements:", nrow(penguins_clean), "\n")

# 4. Summary statistics
summary_stats <- data.frame(
  variable = num_vars,
  n      = sapply(penguins[num_vars], function(x) sum(!is.na(x))),
  mean   = sapply(penguins[num_vars], mean,   na.rm = TRUE),
  median = sapply(penguins[num_vars], median, na.rm = TRUE),
  min    = sapply(penguins[num_vars], min,    na.rm = TRUE),
  max    = sapply(penguins[num_vars], max,    na.rm = TRUE),
  sd     = sapply(penguins[num_vars], sd,     na.rm = TRUE)
)
summary_stats$range <- summary_stats$max - summary_stats$min
cat("\n--- Summary statistics ---\n")
print(format(summary_stats, digits = 4), row.names = FALSE)

# =====================================================================
# TASK 2: HYPOTHESIS TESTING
# =====================================================================

# i. Two-sample t-test: body mass, Adelie vs Gentoo
# H0: mean mass (Adelie) = mean mass (Gentoo);  H1: they differ
adelie <- penguins$body_mass_g[penguins$species == "Adelie"  & !is.na(penguins$body_mass_g)]
gentoo <- penguins$body_mass_g[penguins$species == "Gentoo"  & !is.na(penguins$body_mass_g)]

# Assumption checks
cat("\nShapiro-Wilk Adelie:\n"); print(shapiro.test(adelie))
cat("Shapiro-Wilk Gentoo:\n");   print(shapiro.test(gentoo))
cat("Variance test (F test):\n"); print(var.test(adelie, gentoo))

t_res <- t.test(adelie, gentoo, alternative = "two.sided",
                var.equal = FALSE, conf.level = 0.95)  # Welch test
cat("\n--- Welch two-sample t-test ---\n")
print(t_res)
cat("p-value:", format.pval(t_res$p.value, digits = 4), "\n")

# ii. One-sample proportion test: proportion of males = 0.5 ?
sex_known <- penguins$sex[!is.na(penguins$sex)]
n_male  <- sum(sex_known == "male")
n_total <- length(sex_known)
cat("\nMales:", n_male, "of", n_total, "=", round(n_male / n_total, 4), "\n")

prop_res <- prop.test(x = n_male, n = n_total, p = 0.5,
                      alternative = "two.sided", conf.level = 0.95)
cat("\n--- One-sample proportion test ---\n")
print(prop_res)
cat("Exact binomial check:\n")
print(binom.test(n_male, n_total, p = 0.5))

# =====================================================================
# TASK 3: CORRELATION ANALYSIS
# =====================================================================

# 1. Correlation matrix (Pearson, complete cases)
cor_mat <- cor(penguins_clean[, num_vars], method = "pearson")
cat("\n--- Correlation matrix ---\n")
print(round(cor_mat, 3))

# Significance tests for each pair
for (i in 1:3) for (j in (i + 1):4) {
  ct <- cor.test(penguins_clean[[num_vars[i]]], penguins_clean[[num_vars[j]]])
  cat(sprintf("%s vs %s: r = %.3f, p = %s\n", num_vars[i], num_vars[j],
              ct$estimate, format.pval(ct$p.value, digits = 3)))
}

# 2. Heatmap
cor_long <- as.data.frame(as.table(cor_mat))
names(cor_long) <- c("Var1", "Var2", "r")
p_heat <- ggplot(cor_long, aes(Var1, Var2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  labs(title = "Correlation matrix of numeric variables", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_heat)
ggsave("figures/correlation_heatmap.png", p_heat, width = 6, height = 5, dpi = 200)

# 3. Discussion: see report. Strongest: flipper length ~ body mass (r ~ 0.87);
#    strongest negative: flipper length ~ bill depth (r ~ -0.58).

# =====================================================================
# TASK 4: REGRESSION MODELING
# =====================================================================

# i. Simple linear regression
m1 <- lm(body_mass_g ~ flipper_length_mm, data = penguins_clean)
cat("\n--- Simple linear regression ---\n")
print(summary(m1))
cat("Equation: body_mass_g =", round(coef(m1)[1], 2), "+",
    round(coef(m1)[2], 2), "* flipper_length_mm\n")
cat("R-squared:", round(summary(m1)$r.squared, 4),
    " Adj R-squared:", round(summary(m1)$adj.r.squared, 4), "\n")

# ii. Multiple linear regression
m2 <- lm(body_mass_g ~ flipper_length_mm + bill_length_mm + bill_depth_mm,
         data = penguins_clean)
cat("\n--- Multiple linear regression ---\n")
print(summary(m2))
print(confint(m2))

# Model comparison
cat("\n--- Model comparison ---\n")
print(anova(m1, m2))
print(data.frame(model = c("Simple", "Multiple"),
                 R2     = c(summary(m1)$r.squared,     summary(m2)$r.squared),
                 AdjR2  = c(summary(m1)$adj.r.squared, summary(m2)$adj.r.squared),
                 RSE    = c(sigma(m1), sigma(m2)),
                 AIC    = c(AIC(m1), AIC(m2))))

# Diagnostics: multicollinearity (VIF computed by hand) and residuals
vif_manual <- sapply(c("flipper_length_mm", "bill_length_mm", "bill_depth_mm"),
  function(v) {
    others <- setdiff(c("flipper_length_mm", "bill_length_mm", "bill_depth_mm"), v)
    1 / (1 - summary(lm(reformulate(others, v), data = penguins_clean))$r.squared)
  })
cat("\nVIF:\n"); print(round(vif_manual, 2))
cat("Shapiro-Wilk on MLR residuals:\n"); print(shapiro.test(residuals(m2)))

png("figures/regression_diagnostics.png", width = 1400, height = 600, res = 150)
par(mfrow = c(1, 2))
plot(fitted(m2), residuals(m2), xlab = "Fitted", ylab = "Residual",
     main = "Residuals vs fitted (MLR)", pch = 19, col = rgb(0, 0, 0, 0.4))
abline(h = 0, col = "red")
qqnorm(residuals(m2), main = "Normal Q-Q (MLR residuals)"); qqline(residuals(m2), col = "red")
dev.off()

# =====================================================================
# TASK 5: VISUALIZATION
# =====================================================================

# 1. Scatter plot with regression line
p_scatter <- ggplot(penguins_clean, aes(flipper_length_mm, body_mass_g)) +
  geom_point(alpha = 0.6, colour = "#1f4e79") +
  geom_smooth(method = "lm", se = TRUE, colour = "#c0392b") +
  labs(title = "Body mass vs flipper length with regression line",
       x = "Flipper length (mm)", y = "Body mass (g)") +
  theme_minimal()
print(p_scatter)
ggsave("figures/scatter_regression.png", p_scatter, width = 6.5, height = 4.5, dpi = 200)

# 2. Boxplot by species
p_box <- ggplot(penguins_clean, aes(species, body_mass_g, fill = species)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Body mass distribution by species",
       x = "Species", y = "Body mass (g)") +
  theme_minimal()
print(p_box)
ggsave("figures/boxplot_species.png", p_box, width = 6.5, height = 4.5, dpi = 200)

# 3. Histogram with normal curve
bm <- penguins_clean$body_mass_g
p_hist <- ggplot(penguins_clean, aes(body_mass_g)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20,
                 fill = "#7fb3d5", colour = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(bm), sd = sd(bm)),
                colour = "#c0392b", linewidth = 1) +
  labs(title = "Histogram of body mass with normal curve",
       x = "Body mass (g)", y = "Density") +
  theme_minimal()
print(p_hist)
ggsave("figures/histogram_normal.png", p_hist, width = 6.5, height = 4.5, dpi = 200)

cat("\nShapiro-Wilk on body mass:\n"); print(shapiro.test(bm))

# =====================================================================
# END OF SCRIPT
# =====================================================================
