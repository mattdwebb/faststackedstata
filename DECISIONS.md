# Decisions

- Saturated official default, not uncorrected stacking or interacted FE.
- Required full integer event window; reference -1; all nonnegative events averaged.
- Missing/+Inf never treated; zero remains an actual date.
- Strict complete panel support and time-constant panel weights within stack.
- Contemporary additive covariates labelled as an adjustment extension.
- CR1 with original clusters; cluster-t p-values rather than R helper normal CIs.
- Fail if deletion or assignment changes eligible treated-cohort support.
- Supplied RI assignments only; standalone generation is explicitly unsupported.
- Point-only fitting skips clustered inference.
- Original project code MIT, copyright 2026 Matt Webb.

