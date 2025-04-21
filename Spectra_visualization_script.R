################################
## Load required packages
################################

library(readxl)
library(prospectr) 
library (stringr)

#######################################
## Load Data
#######################################

raw_data <- read_excel("C:/~~/~~/raw_data.xlsx")

##Average of replicas according to the origin (ID)
raw_data <- raw_data %>%
  mutate(ID = gsub("_R1|_R2|_R3", "", ID))

avg_data <- raw_data %>%
  group_by(ID) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  mutate(Origin = str_sub(ID, 1, 2))

avg_origin <- avg_data %>%
  group_by(Origin) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long_avg_origin <- avg_origin %>%
  pivot_longer(cols = -Origin, names_to = "Drift_time", values_to = "Intensity") %>%
  mutate(Wavelength = as.numeric(Wavelength))

origin_color <- c("ES" = "#F58684", "IT" = "#71BCFF", "MR" = "#F5AC3D", "PT" = "#87F951")


  spectra_plot <- ggplot(long_avg_origen, aes(x = Drift_time, y = Intensity, color = Origin, group = Origin)) +
  geom_line(size = 1.2, alpha = 0.9) +
  scale_color_manual(name = "Origin", values = origin_color) +
  labs(title = "Ion Mobility Sum Spectra",
       x = "Drift time [RIP relative] (ms)",
       y = "Intensy (V)") +
  theme_minimal()

print(spectra_plot)
ggsave("figures/spectra_raw.png", spectra_plot, width = 10, height = 6, dpi = 300)


############################################################
##  Automatically save loaded packages and their versions
############################################################

loaded_pkgs <- sessionInfo()$otherPkgs

pkg_versions <- sapply(loaded_pkgs, function(pkg) {
  paste0(pkg$Package, "==", pkg$Version)
})

req_file <- "requirements.txt"

if (file.exists(req_file)) {
  existing_lines <- readLines(req_file)
  existing_pkgs <- sub("(.*)==.*", "\\1", existing_lines)
} else {
  existing_pkgs <- character(0)
}

new_pkgs <- pkg_versions[!(names(pkg_versions) %in% existing_pkgs)]



if (length(new_pkgs) > 0) {
  cat(new_pkgs, file = req_file, sep = "\n", append = TRUE)
  message("New packaged added to requirementss.txt")
} else {
  message("No new packages to add.")
}

cat(paste0(new_pkgs, "\n"), file = req_file, sep = "", append = TRUE)

}