source("parse_sceptic_outputs.R")
library(dplyr)
library(readr)

trial_out_dir <- "/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs"

#datasets <- c("mmclock_meg", "mmclock_fmri")
#datasets <- c("specc")
#datasets <- c("mmclock_fmri")
datasets <- c("mmclock_meg")

#models <- c("decay", "fixed", "decay_factorize", "decay_uniform", "decay_ps_equate", "decay_uniform_ps_equate")
models <- c("decay")
fitting <- c("ffx", "mfx")
factorize <- c("factorize", "separate")
decay <- c("uniform", "selective")
ps <- c("psnarrow", "psequate")

mdf <- do.call(expand.grid, list(models, factorize, decay, ps))
models <- apply(mdf, 1, paste, collapse="_")
#models <- c(models, "decay_factorize_selective_psequate_fixedparams_meg")
#models <- c(models, "decay_factorize_selective_psequate_fixedparams_fmri")
#models <- c("decay_factorize_selective_psequate_fixedparams")
#models <- c("specc_decay_factorize_selective_psequate_fixedparams")
#models <- c("fixed_uv_ureset", "fixed_uv_baked_ureset")
#models <- c("fixed_uv_ureset_fixedparams_fmri")
models <- c("fixed_uv_ureset_fixedparams_meg")
#mdf$model <- sub("_", "/", mdf$model) #replace just the first instance since it's organized by ffx and mfx folders

for (d in datasets) {
  for (m in models) {
    for (f in fitting) {
      outdir <- file.path("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out", d, f, m)
      if (!dir.exists(outdir)) { next } #skip non-existent outputs
      subj_dir <- file.path("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/subjects", d)
      
      file.copy(file.path(outdir, paste0(d, "_", m, "_", f, "_sceptic_trial_outputs_by_timestep.csv")), trial_out_dir, overwrite=TRUE)
      file.copy(file.path(outdir, paste0(d, "_", m, "_", f, "_sceptic_global_statistics.csv")), trial_out_dir, overwrite=TRUE)
      file.copy(file.path(outdir, paste0(d, "_", m, "_", f, "_sceptic_basis.csv")), trial_out_dir, overwrite=TRUE)
      sceptic_stats <- parse_sceptic_outputs(outdir, subj_dir)
      write_csv(sceptic_stats, path=file.path(trial_out_dir, paste(d, m, f, "trial_statistics.csv.gz", sep="_")))
    }
  }
}

print(warnings())
