# Limitations

## Scope
- Neutral atoms only. For now. 
- No general arbitrary fractional-charge refinement inside Phenix. For now. 
- No claim that this is an upstream-supported Phenix mode.

## Technical caveats

- This is an unofficial patch wrapper. It depends on internal Phenix/CCTBX hook points staying similar across versions. They likely will not.
- The wrapper intercepts `F_calc` and atomic gradients, but it does not rewrite every possible auxiliary code path in Phenix. Other versions will likely break. 
- It assumes standard neutral X-ray factors are a suitable starting point for the Mott-Bethe transform. This is generally close to true. 
- Low-angle reflections can still be numerically delicate because the Mott-Bethe factor scales like `1/s^2`. Have not found to be an issue on most real data. 

## Practical caveats

- Expect some version sensitivity if you move between major Phenix builds. E.g. completely breaking. 
- Treat this as a research/lab tool, not a polished official backend. It is a thing I wanted to get working for another project to be posted soon. Do not use for real science yet. 
- Always compare against ordinary `phenix.refine` or 'REFMAC5' on the same input so you can see what the patch actually changes.

## Suggested validation

For any real project:

1. run ordinary `phenix.refine`
2. run `phenix.refine.mb`
3. compare `R_work`, `R_free`, maps, and coordinate/B-factor shifts
4. inspect whether the patched result is chemically and physically more plausible, not just numerically different
