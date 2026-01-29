# CheMPAS Agent-Driven Development

## Motivation

The upstream MPAS-Dev repositories rely on human PR review, which creates bottlenecks that limit development velocity. CheMPAS removes that constraint by making AI agents first-class participants in the development lifecycle — writing code, reviewing PRs, and running CI — with human oversight reserved for where it matters most.

## Agent Roles

| Role | Agent | Responsibility |
|------|-------|----------------|
| Development | Claude Opus 4.5 (Max) | Feature implementation, bug fixes, refactoring, test writing |
| Review & Planning | Codex 5.2 (Extra High) | Merge reviews, architectural planning, cross-cutting analysis |
| CI Verification | TBD | Build validation, test execution, result reporting on PRs |

### Development Agent

The development agent works from issues or direct instructions from the project lead. It operates on feature branches, opens PRs with clear summaries and test plans, and provides evidence that the code compiles and passes relevant tests.

### Review Agent

The review agent evaluates PRs for correctness, style, performance, and scientific integrity. It approves, requests changes, or flags items for human attention. The review agent must be a different model/vendor from the development agent to provide genuine independence of perspective.

### CI Agent

A dedicated CI role (agent or automation) that builds the code, runs test cases, and reports results on PRs. This decouples "does it work" from "is it good," keeping the review agent focused on design and correctness rather than mechanical verification.

## Workflow

1. **Issue or instruction** defines the work to be done.
2. **Development agent** creates a feature branch, implements the change, builds, tests, and opens a PR.
3. **CI agent** runs on the PR: builds with and without MUSICA, runs test cases, reports pass/fail.
4. **Review agent** evaluates the PR: code quality, architectural fit, correctness, and whether the change requires human review.
5. **Merge** proceeds once the review agent approves, unless the change is flagged for human review.
6. **Human oversight** is retained for the categories listed below.

## Human Review Gates

Certain categories of changes must always be flagged for human review, regardless of agent confidence:

- **Scientific correctness**: Changes to physics parameterizations, chemistry solvers, numerical methods, or tendency calculations.
- **Registry.xml modifications**: Bad Registry edits break the entire build in non-obvious ways. The Registry is the metadata backbone of MPAS and changes propagate through code generation.
- **Architectural decisions**: New modules, major refactors, dependency changes, or build system modifications.
- **External interfaces**: Changes to MUSICA-Fortran coupling, MPI communication patterns, or I/O formats.

## Principles

- **Agents build and test before committing.** Every PR must include evidence that the code compiles and passes relevant tests.
- **Agents do not merge to `main` without review.** The development agent and review agent must be different models/vendors.
- **Scientific correctness is non-negotiable.** Changes to physics, chemistry, or numerical methods are always flagged for human review.
- **Minimal, focused changes.** One PR, one purpose. Resist scope creep and over-engineering.
- **The codebase is the source of truth.** Agents read before writing, understand before modifying.
- **Context matters.** Review agents must have access to project documentation (CLAUDE.md, ARCHITECTURE.md, etc.) to provide meaningful reviews.

## Considerations

### Cross-Model Review Independence

Using different models/vendors for development and review is deliberate. A single model reviewing its own output has correlated blind spots. Cross-model review surfaces different failure modes — one model's confident mistake may be another's obvious flag.

### Fortran and Domain Expertise

Not all models handle Fortran equally well. The development and review agents should be validated on their ability to:

- Read and write Fortran 2008 with MPI
- Understand MPAS Registry conventions and code generation
- Reason about array indexing, unit conversions, and physical consistency
- Navigate the MPAS framework's pool/field data structures

### Failure Modes to Watch

- **Physically wrong but compilable code**: Incorrect units, sign errors in tendencies, wrong operator stencils. Tests help but may not catch subtle numerical issues.
- **Shallow reviews**: An agent that approves because the code "looks reasonable" without verifying logic against the existing codebase.
- **Over-engineering**: Agents tend toward abstraction and generalization. For a Fortran codebase, straightforward is better.
- **Registry drift**: Accumulating unused variables or inconsistent metadata in Registry.xml without cleanup.
