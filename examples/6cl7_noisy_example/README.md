# 6cl7 Noisy Example

This is a recorded worked example showing that the neutral-atom Mott-Bethe wrapper can run a full protein-scale `phenix.refine` job on a noisy `6cl7` electron-diffraction dataset. 

## Dataset

The example used a previously generated dataset in the main `Charges` workspace:

- model: neutral `6cl7`
- resolution: `2.0 A`
- noise: about `20%` achieved `Rmerge` proxy
- empirical thickness setting: `100 nm`

Charges is a workspace not publicly available yet, but was used to generate simulated structure factors from a perfect 100nm thick crystal with an overall R merge of 20% non-uniformly over the resolution range. 
The thickness approximation of dynamical scattering is not a full multi-slice but an approximate toy model for corrupting Is/Fs to see how the refinement changes. 

In the current workspace the source files were:

- `/path/to/runs/6cl7_stress_test/thickness_runs/neutral_rmerge_20_thickness_100nm/reference_model.pdb`
- `/path/to/runs/6cl7_stress_test/thickness_runs/neutral_rmerge_20_thickness_100nm/simulated_obs.mtz`

Replace those with your own model PDB and MTZ.

## Wrapper command

```bash
bin/phenix.refine.mb \
  --mb-xray-table wk1995 \
  --mb-electron-voltage-kv 200 \
  reference_model.pdb \
  simulated_obs.mtz \
  data_manager.model.type=electron \
  data_manager.miller_array.labels.name=FSIM,SIGFSIM \
  data_manager.miller_array.labels.name=FREE \
  refinement.main.scattering_table=electron \
  refinement.main.number_of_macro_cycles=3 \
  refinement.output.write_maps=False \
  refinement.output.write_map_coefficients=True \
  overwrite=true \
  output.prefix=phenix_mb_run
```

## Result

Reference ordinary `phenix.refine` run on the same input:

- start `Rwork/Rfree = 0.2853 / 0.2830`
- final `Rwork/Rfree = 0.2748 / 0.3154`

Neutral Mott-Bethe wrapper run:

- start `Rwork/Rfree = 0.2924 / 0.2890`
- final `Rwork/Rfree = 0.2748 / 0.3127`

## Interpretation

- the wrapper completed `6cl7` refinement successfully
- in this specific case the final `Rwork` matched the ordinary Phenix run to four decimal places
- the final `Rfree` was slightly lower with the Mott-Bethe wrapper

This is a worked example, not a general claim that the wrapper will always improve refinement statistics. In fact, our experience with SOURCE EC MB in REFMAC5 indicates it is essnetially identical to using ED SF libraries. 

See [summary.json]for a machine-readable version of the same result.
