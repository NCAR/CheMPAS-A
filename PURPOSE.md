# CheMPAS Approach

## The CheMPAS Response

CheMPAS was created as a clean fork from [NCAR/MPAS-Model-ACOM-dev](https://github.com/NCAR/MPAS-Model-ACOM-dev) (itself a fork of [MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model)) with no intent to sync back. This deliberate divergence allows CheMPAS to adopt an agent-driven development model that would be incompatible with the upstream governance structure.

The core insight: **if agents are generating the code, agents should be reviewing it too, with humans reserved for what humans are actually better at** — scientific judgment, architectural vision, and physical intuition.

## Agent-Driven Development

### Roles

| Role | Agent | Responsibility |
|------|-------|----------------|
| Development | Claude Opus 4.6 (Max) | Feature implementation, bug fixes, refactoring, test writing |
| Review | Codex 5.4 (Extra High) | Merge reviews, planning review |
| Deep Research | Gemini 3 Pro (Deep Think) | Codebase-wide analysis, scientific review, architecture planning |
| CI Verification | TBD | Build validation, test execution, result reporting |

### Why Three Models from Three Vendors

The development, review, and research agents are deliberately different models from different vendors (Anthropic, OpenAI, Google DeepMind). A single model reviewing its own output has correlated blind spots — it will be confident about the same things it got wrong. Three-vendor independence ensures no single architectural bias, training artifact, or failure mode can pass through unchallenged.

Each model is assigned to the role that matches its strengths: Claude Opus 4.6's debugging accuracy and SWE-bench performance for development, Codex 5.4's methodical code review at Extra High reasoning for PR review, and Gemini 3 Pro's 1M+ token context and graduate-level science reasoning (GPQA Diamond 93.8%) for deep research and scientific validation. See [BENCHMARKS.md](BENCHMARKS.md) for the full comparison.

### What Humans Do

Agents handle velocity. Humans handle judgment. Specifically, human review is required for:

- **Scientific correctness**: Physics parameterizations, chemistry solvers, numerical methods, tendency calculations. An agent can write code that compiles, runs, and produces plausible-looking output while being physically wrong. Unit errors, sign flips in tendencies, incorrect operator stencils — these require domain expertise to catch.
- **Architectural decisions**: New modules, major refactors, dependency changes, build system modifications.
- **Registry changes**: The MPAS Registry.xml is the metadata backbone of the model. Bad edits break the build through code generation in non-obvious ways.
- **External interfaces**: MUSICA-Fortran coupling, MPI communication patterns, I/O format changes.

The deep research agent adds a layer that neither the development nor review agent can provide: full-codebase comprehension in a single context. Before major changes, it can ingest the entire source tree alongside scientific literature to identify coupling risks, unit inconsistencies, and deviations from published methods. Its role is to inform, not to write or approve code.

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
