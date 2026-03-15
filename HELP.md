# Help and Troubleshooting
Most of this is stuff that came up trying to get this up and running on a Mac-silicon device. YMMV.

## Wrapper-specific flags

- `--mb-xray-table {wk1995,it1992}`
- `--mb-electron-voltage-kv FLOAT`
- `--mb-verbose`
- `--mb-help`

All other arguments are passed through to `phenix.refine`.

## Common command pattern

```bash
bin/phenix.refine.mb \
  --mb-xray-table wk1995 \
  --mb-electron-voltage-kv 200 \
  model.pdb data.mtz \
  data_manager.model.type=electron \
  data_manager.miller_array.labels.name=F,SIGF \
  data_manager.miller_array.labels.name=FREE \
  refinement.main.scattering_table=electron \
  overwrite=true
```

## Common problems

### `Phenix Python was not found`

Set:

```bash
export PHENIX_PYTHON=/path/to/phenix/python
```

### OpenMP shared-memory errors

The launcher already sets:

- `OMP_NUM_THREADS=1`
- `KMP_AFFINITY=disabled`
- `TMPDIR=/tmp`

If your environment still overrides these, reset them before running.

### Wrong MTZ labels

This wrapper does not guess your observation columns. Pass the same label PHIL you would use for ordinary `phenix.refine`.

Typical simulated-amplitude case:

```bash
data_manager.miller_array.labels.name=F,SIGF
data_manager.miller_array.labels.name=FREE
```

### Does this support fractional charges?

No. This standalone wrapper is for neutral-atom Mott-Bethe refinement only. Our research tool for this should be out eventually. Even then, do not trust this. If you have data with charges in it, please go use the SHELXL implementation described in: https://github.com/CF-CSA/iSFAC_code  

### Does this replace all Phenix scattering behavior everywhere?

No. It patches the main reciprocal-space `F_calc` and atomic gradient path used during refinement. It is intentionally a lab hack, not a complete rewrite of Phenix internals. Ain't nobody got time for that. 
