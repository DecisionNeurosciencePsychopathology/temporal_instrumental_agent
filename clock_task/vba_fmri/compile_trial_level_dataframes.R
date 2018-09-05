source("parse_sceptic_outputs.R")
library(dplyr)
library(readr)

trial_out_dir <- "/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs"

datasets <- c("mmclock_meg", "mmclock_fmri")
#datasets <- c("mmclock_fmri")
models <- c("decay", "fixed", "decay_factorize")
fitting <- c("ffx", "mfx")

for (d in datasets) {
  for (m in models) {
    for (f in fitting) {
      outdir <- file.path("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out", d, f, m)
      subj_dir <- file.path("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/subjects", d)

      sceptic_stats <- parse_sceptic_outputs(outdir, subj_dir)
      write_csv(sceptic_stats, path=file.path(trial_out_dir, paste(d, m, f, "trial_statistics.csv.gz", sep="_")))
    }
  }
}
