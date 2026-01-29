# CheMPAS Manifesto

## The Problem

Large-scale earth system models like MPAS are developed under institutional governance structures where human PR review is the primary quality gate. This worked when contributions came at human pace. It no longer does.

The [MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model) repository currently has 75 open pull requests, with the oldest dating back to August 2017. Approximately 40 PRs older than a year remain unmerged. The pattern is consistent with batch merges around release cycles rather than continuous integration — PRs queue up, wait months, and get flushed during release pushes.

This is not a criticism of the MPAS-Dev team. The bottleneck is structural: a small number of qualified reviewers, each with their own research and development responsibilities, gatekeeping a large and complex Fortran/MPI codebase with multiple cores (atmosphere, ocean, sea ice, land ice) and institutional stakeholders (NCAR, LANL).

## The Catalyst

Frontier coding agents have made this bottleneck dramatically worse. An agent can now produce well-formed, compilable, tested PRs at a rate that no human review team can absorb. The throughput mismatch between code generation and code review has gone from inconvenient to untenable.

The upstream repos face an impossible choice: either refuse agent-generated contributions (losing velocity) or drown in a review backlog that grows faster than it can be processed.

## The CheMPAS Response

CheMPAS was created as a clean fork from [NCAR/MPAS-Model-ACOM-dev](https://github.com/NCAR/MPAS-Model-ACOM-dev) (itself a fork of [MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model)) with no intent to sync back. This deliberate divergence allows CheMPAS to adopt an agent-driven development model that would be incompatible with the upstream governance structure.

The core insight: **if agents are generating the code, agents should be reviewing it too, with humans reserved for what humans are actually better at** — scientific judgment, architectural vision, and physical intuition.

## Agent-Driven Development

### Roles

| Role | Agent | Responsibility |
|------|-------|----------------|
| Development | Claude Opus 4.5 (Max) | Feature implementation, bug fixes, refactoring, test writing |
| Review & Planning | Codex 5.2 (Extra High) | Merge reviews, architectural planning, cross-cutting analysis |
| CI Verification | TBD | Build validation, test execution, result reporting |

### Why Different Models

The development and review agents are deliberately different models from different vendors. A single model reviewing its own output has correlated blind spots — it will be confident about the same things it got wrong. Cross-vendor review provides genuine independence of perspective.

### What Humans Do

Agents handle velocity. Humans handle judgment. Specifically, human review is required for:

- **Scientific correctness**: Physics parameterizations, chemistry solvers, numerical methods, tendency calculations. An agent can write code that compiles, runs, and produces plausible-looking output while being physically wrong. Unit errors, sign flips in tendencies, incorrect operator stencils — these require domain expertise to catch.
- **Architectural decisions**: New modules, major refactors, dependency changes, build system modifications.
- **Registry changes**: The MPAS Registry.xml is the metadata backbone of the model. Bad edits break the build through code generation in non-obvious ways.
- **External interfaces**: MUSICA-Fortran coupling, MPI communication patterns, I/O format changes.

Everything else — style, correctness of straightforward logic, test coverage, documentation, build fixes — agents can review faster and more consistently than humans.

## Principles

1. **Agents build and test before committing.** No exceptions. Every PR includes evidence of compilation and test results.
2. **No self-review.** The development agent and review agent are always different models from different vendors.
3. **Scientific correctness is non-negotiable.** Physics and chemistry changes are always flagged for human review.
4. **One PR, one purpose.** Minimal, focused changes. Resist scope creep and over-engineering.
5. **Read before writing.** The codebase is the source of truth. Agents understand existing code before modifying it.
6. **Continuous flow over batch releases.** The point of removing the review bottleneck is to enable continuous integration, not to accumulate a different kind of backlog.

## On Divergence

CheMPAS will diverge from upstream MPAS. This is intentional. The upstream project serves a broad community with conservative governance appropriate to its role. CheMPAS serves a focused research objective — coupled atmospheric chemistry on MPAS meshes — and can move at the pace that objective demands.

The two projects share a common ancestor. They do not need to share a future.
