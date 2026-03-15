# Phenix Neutral-Atom Mott-Bethe Hack

This directory is a standalone export of the neutral-atom `phenix.refine` Mott-Bethe wrapper built in the main `Charges` workspace. (Not public yet but soon - MWM)

It is intentionally small:
- `bin/phenix.refine.mb`
- `scripts/phenix_refine_mb_neutral.py`
- docs and examples

The goal is simple: run ordinary `phenix.refine`, but intercept the internal neutral X-ray `F_calc` and gradient path so the refinement uses neutral-atom electron scattering from the Mott-Bethe relation instead of the usual built-in electron table.

## What it does

For neutral atoms it replaces the internal scattering calculation with:

`F_e(h) = C / s^2 * (F_n(h) - F_x(h))`

where:

- `F_x(h)` is the usual neutral X-ray structure factor already computed by Phenix/CCTBX
- `F_n(h)` is a nuclear pseudo-structure built from atomic numbers
- `C` is the relativistic Mott-Bethe prefactor

The same chain rule is applied to the atomic gradients, so `phenix.refine` can still refine coordinates, isotropic `B` values, and occupancies under the patched model.

## Files

- `bin/phenix.refine.mb`: launcher
- `scripts/phenix_refine_mb_neutral.py`: monkeypatch implementation
- `HOWTO.md`: setup and usage
- `LIMITATIONS.md`: current caveats
- `examples/6cl7_noisy_example.md`: worked example on a noisy `6cl7` dataset

## Quick start

```bash
bin/phenix.refine.mb \
  --mb-xray-table wk1995 \
  --mb-electron-voltage-kv 200 \
  model.pdb data.mtz \
  data_manager.model.type=electron \
  data_manager.miller_array.labels.name=FSIM,SIGFSIM \
  data_manager.miller_array.labels.name=FREE \
  refinement.main.scattering_table=electron \
  overwrite=true
```

Use `--mb-help` to see the wrapper-specific options.
