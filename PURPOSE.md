# CheMPAS Development Workflow

CheMPAS is a fork of
[NCAR/MPAS-Model-ACOM-dev](https://github.com/NCAR/MPAS-Model-ACOM-dev),
itself a fork of [MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model).
It is maintained independently for coupled atmospheric chemistry work on MPAS
meshes and is not intended to track upstream development closely.

## Tools In Use

| Function | Tool | Current use |
|----------|------|-------------|
| Development | Claude Opus 4.6 (Max) | Feature work, bug fixes, refactoring, test scaffolding, doc edits |
| Review | Codex 5.4 (Extra High) | Code review, plan review, design checks |
| Deep analysis | Gemini 3 Pro (Deep Think) | Full-tree analysis, scientific literature review, architecture work |
| Build and run | Local shell tools | Compile MPAS, run test cases, inspect logs and outputs |
| Verification scripts | Repo `scripts/` | Initialization, diagnostics, plotting, gate checks |

The tool mix can change over time. The repo documentation should reflect the
tools currently in use rather than treat them as fixed policy.

## Workflow

1. Define the task or issue.
2. Gather local code and documentation context before editing.
3. Implement the change in a focused patch.
4. Build and run the relevant case when the change affects code behavior.
5. Review the change with a different model/vendor than the one used for
   implementation.
6. Keep plans, run notes, and architecture docs in sync with the code.
7. Escalate to human review for scientific, architectural, Registry, or
   external-interface changes.

## Human Review Gates

Human review is required for:

- physics, chemistry, and numerical-method changes
- `Registry.xml` edits
- new modules, major refactors, or build-system changes
- MUSICA, MPI, and I/O interface changes

## Working Norms

- build and test before committing
- keep changes small and specific
- treat the codebase as the source of truth
- record nontrivial validation results in the repo
- prefer explicit review evidence over informal confidence
