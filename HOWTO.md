# How To Use

## 1. Requirements

You need:

- a working Phenix installation. This was built using 1.21.1-5286.
- the path to the Python executable inside that Phenix installation
- a model file (`PDB` or `mmCIF`)
- an MTZ with amplitudes or intensities plus `FREE` flags, just as for ordinary `phenix.refine`

The launcher assumes:

```bash
/path/to/phenix-1.21.1-5286/python.app/Contents/MacOS/python
```
- This will be different on non-mac systems.
  
If that path is wrong on your machine, set:

```bash
export PHENIX_PYTHON=/path/to/phenix/python
```

## 2. Wrapper-only options

These options belong to the Mott-Bethe wrapper itself:

- `--mb-xray-table {wk1995,it1992}`
- `--mb-electron-voltage-kv FLOAT`
- `--mb-verbose`
- `--mb-help`

Everything else is passed directly to `phenix.refine`. This only works on the command line, so no gui support. Its a hack for research, not anything condoned by the Phenix people. 

## 3. Minimal command

```bash
bin/phenix.refine.mb \
  --mb-xray-table wk1995 \
  --mb-electron-voltage-kv 200 \
  model.pdb data.mtz \
  data_manager.model.type=electron \
  data_manager.miller_array.labels.name=F,SIGF \
  data_manager.miller_array.labels.name=FREE \
  refinement.main.scattering_table=electron \
  refinement.main.number_of_macro_cycles=1 \
  overwrite=true \
  output.prefix=mb_trial
```

## 4. What is actually patched

The wrapper patches two internal Phenix/CCTBX points:

1. the `F_calc` computation
2. the gradient calculation with respect to atomic parameters

That is enough for ordinary neutral-atom reciprocal-space refinement to proceed. Anisotropic Bs have not been tested. Nor has occupancies yet. 

## 5. Recommended use

Use this hack when:

- you want try neutral-atom Mott-Bethe scattering inside an otherwise standard `phenix.refine` run
- you want something light enough to try immediately without modifying Phenix source. 
- you are testing whether neutral Mott-Bethe alone changes refinement behavior on your electron-diffraction data

## 6. Not the right tool for

This standalone hack is not the right tool when:
- you need arbitrary fractional charge refinement (yet)
- you need a maintained upstream Phenix feature rather than a research patch
- you need a formally benchmarked production method across many Phenix versions
- you want to actually publish the results
