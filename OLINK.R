
#OLINK CODE:

#required packages for this analysis: OLinkAnalyze, dplyr, tidyverse, stringr, data.table, purrr, ggplot2, ggrepel


#Set PrimaryDirectory

dirname(rstudioapi::getActiveDocumentContext()$path)            
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# Set input directory

setwd(PrimaryDirectory)
setwd("../data/")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

# Set metadata directory

setwd(PrimaryDirectory)
setwd("../metadata/")
MetaDirectory <- getwd()
setwd(PrimaryDirectory)


# Create output directory

setwd(PrimaryDirectory)
dir.create("Output_Olink", showWarnings = FALSE)
setwd("Output_Olink")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)



#Read in Data


setwd(PrimaryDirectory)
df <- read.csv("E2528888_MM_StudyB_Extended_NPX_2025-06-11.csv", sep = ";")  
str(df)

npx <- df

#Create frame without controls
npx_samples <- npx %>%
  filter(SampleType == "SAMPLE") 

#PCAA
p_pca <- olink_pca_plot(
  df          = npx_samples,
  color_g     = "PlateID",
  label_samples = TRUE,
  quiet       = TRUE
)

p_pca[[1]]

p_pca_out <- olink_pca_plot(
  df          = npx_samples,
  color_g     = "PlateID",
  label_samples = TRUE,
  outlierDefX = 2.5,            # SD threshold on PC1
  outlierDefY = 4,              # SD threshold on PC2
  outlierLines = TRUE,
  quiet       = TRUE
)

p_pca_out[[1]]   # plot with labelled outliers

PCA object from:
  p_pca_out <- olink_pca_plot(..., outlierDefX = ..., outlierDefY = ..., outlierLines = TRUE)

pca_df <- p_pca_out[[1]]$data

pca_outlier_df <- pca_df %>%
  dplyr::filter(Outlier == 1) %>%
  dplyr::select(SampleID, PCX, PCY, Outlier)

pca_outlier_df

p_qc <- olink_qc_plot(
  df               = npx_samples,
  color_g          = "PlateID",     # or QC_Warning, etc.
  label_outliers   = TRUE,
  IQR_outlierDef   = 3,
  median_outlierDef = 3,
  outlierLines     = TRUE
)

qc_df <- p_qc$data

qc_outliers <- qc_df %>%
  dplyr::filter(Outlier == 1) %>%
  dplyr::select(SampleID, Panel, sample_median, IQR, Outlier)

qc_outliers

names(qc_df)

samples_pca <- unique(pca_outlier_df$SampleID)

samples_qc  <- unique(qc_outliers$SampleID)   # if you used olink_qc_plot


npx_samples <- npx_samples %>%
  mutate(SampleQC_status = case_when(
    SampleQC %in% c("Warning", "WARN", "Warn") ~ "Warning",
    SampleQC %in% c("Pass", "PASS")           ~ "Pass", TRUE ~ "Other")  
    
    
    p_dist <- olink_dist_plot(
      df      = npx_samples,        
      color_g = "SampleQC_status")
    
    p_dist
    
    p_dist +
      ggplot2::labs(
        title = "NPX distributions by sample",
        subtitle = "Samples coloured by Sample QC (Pass vs Warning)",
        x = "Sample",
        y = "NPX", save.to.disk = TRUE)
    
    setwd(MetaDirectory)
    
    meta.dat <- fread("OlinkMetadata.csv")
    meta.dat
    
    #Add metadata to data.table
    
    metadata <- meta.dat %>%
      filter(SampleID %in% npx_samples$SampleID)
    
    npx_merged <- npx_samples %>%
      left_join(metadata, by = "SampleID")
    
    sum(is.na(npx_merged$Group))
    
    #T-Tests for each
    long_clean <- npx_merged %>%
      filter(SampleType == "SAMPLE") %>%                      # only real samples
      filter(Time %in% c("Baseline", "Follow-Up")) %>%        # adjust if names differ
      mutate(
        Time    = factor(Time, levels = c("Baseline", "Follow-Up")),
        Group   = factor(Group),
        Patient = factor(Patient)
      ) %>%
      
      group_by(Patient, Group, Time, OlinkID, Assay) %>%
      summarise(NPX = mean(NPX, na.rm = TRUE), .groups = "drop")
    
    run_within_group_ttests <- function(dat, group_name,
                                        fdr_cutoff = 0.05,
                                        effect_cutoff = 0) {
      gdat <- dat %>%
        filter(Group == group_name) %>%
        droplevels()
      
      # wide per patient: Baseline & Follow-Up columns
      wide <- gdat %>%
        pivot_wider(
          names_from  = Time,
          values_from = NPX
        ) %>%
        # keep only complete pairs
        filter(!is.na(Baseline), !is.na(`Follow-Up`))
      
      if (nrow(wide) == 0) {
        warning("No complete Baseline/Follow-Up pairs for group: ", group_name)
        return(tibble())
      }
      
      # split by Assay (one test per protein)
      by_assay <- split(wide, wide$Assay)
      
      res_list <- map(by_assay, function(d) {
        # need at least 2 paired patients for a sensible paired t-test
        if (nrow(d) < 2) return(NULL)
        
        delta <- d$`Follow-Up` - d$Baseline
        mean_delta <- mean(delta, na.rm = TRUE)
        
        tt <- try(t.test(d$Baseline, d$`Follow-Up`, paired = TRUE), silent = TRUE)
        if (inherits(tt, "try-error")) return(NULL)
        
        data.frame(
          Group      = group_name,
          Assay      = d$Assay[1],
          OlinkID    = d$OlinkID[1],
          mean_delta = mean_delta,       # Follow-Up – Baseline
          pvalue     = tt$p.value,
          n_pairs    = nrow(d)
        )
      })
      
      res <- bind_rows(res_list)
      if (nrow(res) == 0) return(res)
      
      res <- res %>%
        mutate(
          adj_pvalue = p.adjust(pvalue, method = "BH"),
          sig = !is.na(adj_pvalue) &
            adj_pvalue < fdr_cutoff &
            abs(mean_delta) > effect_cutoff
        )
      
      res
    }
    
    
    res_ici      <- run_within_group_ttests(long_clean, "ICI",      fdr_cutoff = 0.05, effect_cutoff = 0)
    res_vegfi    <- run_within_group_ttests(long_clean, "VEGFI",    fdr_cutoff = 0.05, effect_cutoff = 0)
    res_icivegfi <- run_within_group_ttests(long_clean, "ICI/VEGFI", fdr_cutoff = 0.05, effect_cutoff = 0)
    
    
    plot_volcano <- function(res_df, title) {
      if (nrow(res_df) == 0) {
        warning("Result data frame is empty for: ", title)
        return(NULL)
      }
      
      # remove NA sig so ggplot doesn't create an 'NA' legend entry
      res_df <- res_df %>% filter(!is.na(sig))
      
      ggplot(res_df, aes(x = mean_delta, y = -log10(adj_pvalue))) +
        geom_point(aes(color = sig), alpha = 0.7) +
        scale_color_manual(
          values = c(`FALSE` = "grey70", `TRUE` = "red"),
          labels = c("Not significant", "Significant"),
          name   = "Significance"
        ) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        theme_minimal() +
        labs(
          title = title,
          x = "Mean ΔNPX (Follow-Up – Baseline)",
          y = "-log10(FDR-adjusted p-value)"
        ) +
        # label only significant proteins
        geom_text_repel(
          data = subset(res_df, sig),
          aes(label = Assay),
          size = 3,
          max.overlaps = 25
        )
    }
    
    vol_ici      <- plot_volcano(res_ici,      "ICI: Baseline vs Follow-Up")
    vol_vegfi    <- plot_volcano(res_vegfi,    "VEGFI: Baseline vs Follow-Up")
    vol_icivegfi <- plot_volcano(res_icivegfi, "ICI/VEGFI: Baseline vs Follow-Up")
    
    vol_ici
    vol_vegfi
    vol_icivegfi
    
    