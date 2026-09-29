Extended Data Fig. 3 source data (posterior expected fledglings, first vs last study year)

Computed from the ordinal nest-success model Q17 (female_fit_Q17.rds / male_fit_Q17.rds) at the first (2003)
and last (2023) study year, other predictors at their centred means (0), random effects excluded.
Expected fledglings = sum over categories (0,1,2,3) of category x posterior category probability, per draw.

EDFig3_expected_fledglings_draws.csv : all posterior draws plotted as violins. sex; year; draw (posterior draw index); expected_fledglings.
EDFig3_summary.csv : per sex x year: n_draws; posterior_mean (plotted point); posterior_median; lower_2.5, upper_97.5 (plotted error bar).
EDFig3_pct_change_2003_2023.csv : per-draw percent change 100*(2023 - 2003)/2003, summarised as median and 2.5/97.5% quantiles.
