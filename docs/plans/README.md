# Planning Docs Conventions

Use this directory for detailed implementation plans and execution logs.

## Canonical Structure

1. `PLAN.md` (repo root) is the high-level index and current-focus summary.
2. `docs/plans/*.md` are detailed, dated plans for individual initiatives.
3. Detailed plan content should not be duplicated in `PLAN.md`.

## File Naming

Use:

`YYYY-MM-DD-short-topic.md`

Example:

`2026-03-06-runtime-tracer-allocation.md`

## Required Status Block

Every plan file should begin with:

```markdown
## Document Status

- `Historical Context:` ...
- `Current State:` ...
- `Use This As:` ...
```

For a ready-to-copy starter, use [TEMPLATE.md](TEMPLATE.md).
