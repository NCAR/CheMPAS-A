# Agent Model Benchmarks

Comparison of the AI models used in CheMPAS agent-driven development.

## Model Roles

| Role | Model | Provider |
|------|-------|----------|
| Development | Claude Opus 4.6 (Max) | Anthropic |
| Review | Codex 5.4 (Extra High) | OpenAI |
| Deep Research | Gemini 3 Pro (Deep Think) | Google DeepMind |

## Why This Three-Model Approach Works for CheMPAS

- **Claude Opus 4.6 for development**: Strong SWE-bench performance, debugging
  accuracy, and terminal proficiency — relevant for writing and testing
  Fortran/MPI code.
- **Codex 5.4 Extra High for review**: Thorough, methodical analysis across large
  codebases. Abstract reasoning advantage may catch issues that Claude misses.
- **Gemini 3 Pro for deep research**: 1M+ token context enables full-codebase
  reasoning, scientific literature review, and architectural analysis that
  neither Claude nor Codex can perform in a single pass. Top scores on
  graduate-level science benchmarks make it well-suited for assessing
  scientific correctness of implementations against published methods.
- **Three-vendor independence**: Three models from three vendors (Anthropic,
  OpenAI, Google DeepMind) with different training data, architectures, and
  failure modes. No single correlated blind spot can pass through all three.

## Fortran and HPC

No model in this lineup has published benchmarks for Fortran 2008, MPI, or
earth system domain code. The closest available study is an [evaluation of the original
OpenAI Codex (GPT-3 era) for HPC parallel programming](https://arxiv.org/html/2306.15121)
(ICPP 2023), which found:

- Codex could generate correct Fortran HPC kernels (AXPY, GEMM, Jacobi, CG)
- Using the `subroutine` keyword in prompts was critical for Fortran quality
- OpenMP and OpenACC parallel Fortran produced better results than other
  parallel programming models

Current-generation models represent a substantial leap from GPT-3-era Codex,
but Fortran-specific validation remains an open gap for all three models.

## Gap

No model has been validated on MPAS-specific patterns: Registry.xml
conventions, pool/field data structures, MPAS operator stencils, or
coupled physics-chemistry tendencies. The human review gates defined in
[AGENTS.md](../../AGENTS.md) exist specifically to cover this gap.

## Role Suitability Summary

| Capability | Best Model | Why |
|------------|-----------|-----|
| Writing Fortran/MPI code | Claude Opus 4.6 | SWE-bench performance, debugging accuracy |
| PR code review | Codex 5.4 Extra High | Methodical, precise bug-finding |
| Full-codebase analysis | Gemini 3 Pro | 1M token context, entire repo in one pass |
| Scientific literature review | Gemini 3 Pro | Graduate-level science reasoning, Deep Think |
| Architecture planning | Gemini 3 Pro + Codex 5.4 | Large context for analysis, abstract reasoning for design |
| CI / build verification | Any (automated) | Mechanical task, model quality less critical |

## Sources

- [Evaluation of OpenAI Codex for HPC Parallel Programming (ICPP 2023)](https://arxiv.org/html/2306.15121)
- [Gemini 3 Pro — Google DeepMind](https://deepmind.google/models/gemini/pro/)
