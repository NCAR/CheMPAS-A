# CheMPAS Plan Index

## Document Status

- `Historical Context:` Top-level index for active and historical planning docs.
- `Current State:` Runtime MICM tracer support is complete and merged (as of
  March 6, 2026).
- `Use This As:` Navigation entry point. Detailed implementation plans live
  under `docs/plans/`.

## Current Focus

1. Follow-on runtime tracer work tracked in [TODO.md](TODO.md).
2. Keep new plan execution details in dated files under `docs/plans/`.

## Active Plans

- [Photolysis and Tropospheric Chemistry](docs/plans/2026-03-06-tuvx-photolysis-integration.md) — LNOx-O3 mechanism for supercell domain, then TUV-x integration (Phase 0 next)

## Historical Plans

- [Generalize Chemistry Coupling (Phases 1 and 2)](docs/plans/2026-03-06-generalize-chemistry-coupling.md)
- [Runtime Chemistry Tracer Allocation](docs/plans/2026-03-06-runtime-tracer-allocation.md)

## Plan Authoring Rules

1. Create one plan file per initiative under `docs/plans/` using
   `YYYY-MM-DD-topic.md`.
2. Include a `## Document Status` section with:
   - `Historical Context`
   - `Current State`
   - `Use This As`
3. When work is complete, update `Current State` with completion date and merge
   status.
4. Keep `PLAN.md` concise; avoid duplicating full task logs here.
