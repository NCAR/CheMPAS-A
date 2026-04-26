# CheMPAS-A Plan Index

## Document Status

- `Historical Context:` Top-level index for active and historical planning docs.
- `Current State:` TUV-x integration now extends through Phase 3 cloud
  attenuation on the idealized supercell development case. The active planning
  focus is Phase 3 hardening, result consolidation, and later photolysis
  realism work.
- `Use This As:` Navigation entry point. Detailed implementation plans live
  under `docs/chempas/plans/`.

## Current Focus

1. Phase 3 follow-up: rebuild/retest after the recent TUV-x wavelength-grid
   ownership fix, tighten remaining cloud-path guards, and finish the deferred
   chemistry-response documentation.
2. Keep new plan execution details in dated files under `docs/chempas/plans/`.

## Active Plans

- [Photolysis and Tropospheric Chemistry](docs/chempas/plans/2026-03-06-tuvx-photolysis-integration.md) — LNOx-O3 mechanism and TUV-x integration through Phase 3 cloud attenuation; active work is Phase 3 hardening and later-scope follow-up

## Historical Plans

- [Generalize Chemistry Coupling (Phases 1 and 2)](docs/chempas/plans/2026-03-06-generalize-chemistry-coupling.md)
- [Runtime Chemistry Tracer Allocation](docs/chempas/plans/2026-03-06-runtime-tracer-allocation.md)

## Plan Authoring Rules

1. Create one plan file per initiative under `docs/chempas/plans/` using
   `YYYY-MM-DD-topic.md`.
2. Include a `## Document Status` section with:
   - `Historical Context`
   - `Current State`
   - `Use This As`
3. When work is complete, update `Current State` with completion date and merge
   status.
4. Keep `PLAN.md` concise; avoid duplicating full task logs here.
