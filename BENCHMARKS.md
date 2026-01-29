# Agent Model Benchmarks

Comparison of the AI models used in CheMPAS agent-driven development, based on
publicly available benchmarks as of January 2026.

## Model Roles

| Role | Model | Provider |
|------|-------|----------|
| Development | Claude Opus 4.5 (Max) | Anthropic |
| Review | GPT-5.2 Codex (Extra High) | OpenAI |
| Deep Research | Gemini 3 Pro (Deep Think) | Google DeepMind |

## Coding Benchmarks

| Benchmark | Claude Opus 4.5 | GPT-5.2 Codex | Gemini 3 Pro | Notes |
|-----------|-----------------|---------------|--------------|-------|
| SWE-bench Verified | **80.9%** | 80.0% | — | Real-world GitHub issue resolution |
| SWE-bench Pro | — | **56.4%** | — | Harder subset of SWE-bench |
| HumanEval | **94.2%** | 91.7% | — | Function-level code generation |
| Terminal-Bench | **59.3%** | 47.6% | — | Command-line proficiency |
| Debugging accuracy | **89%** | 84% | — | Multi-threaded race conditions, memory leaks |
| LMArena Elo | — | — | **1501** | Overall model quality (top of leaderboard) |

## Reasoning Benchmarks

| Benchmark | Claude Opus 4.5 | GPT-5.2 Codex | Gemini 3 Pro (Deep Think) | Notes |
|-----------|-----------------|---------------|---------------------------|-------|
| ARC-AGI-2 | 37.6% | **~53%** | 45.1% (w/ code exec) | Abstract reasoning |
| AIME 2025 | 92.8% | **100%** | — | Mathematical reasoning |
| GPQA Diamond | — | — | **93.8%** | Graduate-level science questions |
| Humanity's Last Exam | — | — | **41.0%** | Hardest general knowledge benchmark |

## Operational Characteristics

| Metric | Claude Opus 4.5 | GPT-5.2 Codex |
|--------|-----------------|---------------|
| Code generation speed | Slower | ~23% faster |
| Token efficiency | ~76% fewer tokens per task | Higher token volume |
| Code volume (Sonar study) | 639K LOC generated | 974K LOC generated |
| Functional pass rate (Sonar) | **83.62%** | 80.66% |

## Context Window

| Model | Context Window | Output Limit | Notes |
|-------|---------------|--------------|-------|
| Claude Opus 4.5 | 200K tokens | — | Standard coding agent context |
| GPT-5.2 Codex | — | — | Multi-context via compaction |
| Gemini 3 Pro | **1M+ tokens** | 64K tokens | 77% recall at 128K in benchmarks |

Gemini 3 Pro's context window is the key differentiator for the deep research
role. The full CheMPAS source tree (~1,815 files) can be loaded in a single
context, enabling codebase-wide analysis that the other models cannot perform.

## Agentic Coding (Qualitative)

Results from practical agentic coding tests are mixed:

- **Claude Opus 4.5** tends to ship working implementations with better
  architecture, but sometimes with rough edges requiring integration work.
- **GPT-5.2 Codex** is faster and can deliver production-ready code, but
  is more prone to API and version mismatches in some test scenarios.
- **GPT-5.2 xhigh reasoning** is described as "extremely thorough, extremely
  precise" for code review, particularly effective at finding bugs and
  inconsistencies across large codebases.

## Fortran and HPC

No model in this lineup has published benchmarks for Fortran 2008, MPI, or
earth system domain code. The closest available study is an [evaluation of the original
OpenAI Codex (GPT-3 era) for HPC parallel programming](https://arxiv.org/html/2306.15121)
(ICPP 2023), which found:

- Codex could generate correct Fortran HPC kernels (AXPY, GEMM, Jacobi, CG)
- Using the `subroutine` keyword in prompts was critical for Fortran quality
- OpenMP and OpenACC parallel Fortran produced better results than other
  parallel programming models
- Fortran's domain-specific nature and large legacy codebase helped
  compensate for lower training data volume compared to mainstream languages

The GPT-5.2 generation represents a substantial leap from GPT-3-era Codex,
but Fortran-specific validation remains an open gap for all three models.

## Why This Three-Model Approach Works for CheMPAS

- **Claude for development**: Higher accuracy on SWE-bench Verified, better
  debugging, stronger terminal proficiency — relevant for writing and testing
  Fortran/MPI code.
- **GPT-5.2 xhigh for review**: Thorough, methodical analysis across large
  codebases. Abstract reasoning advantage may catch issues that Claude misses.
- **Gemini 3 Pro for deep research**: 1M+ token context enables full-codebase
  reasoning, scientific literature review, and architectural analysis that
  neither Claude nor GPT-5.2 can perform in a single pass. Top scores on
  graduate-level science benchmarks (GPQA Diamond 93.8%) make it well-suited
  for assessing scientific correctness of implementations against published
  methods.
- **Three-vendor independence**: Three models from three vendors (Anthropic,
  OpenAI, Google DeepMind) with different training data, architectures, and
  failure modes. No single correlated blind spot can pass through all three.

## Gap

No model has been validated on MPAS-specific patterns: Registry.xml
conventions, pool/field data structures, MPAS operator stencils, or
coupled physics-chemistry tendencies. The human review gates defined in
[AGENTS.md](AGENTS.md) exist specifically to cover this gap.

## Role Suitability Summary

| Capability | Best Model | Why |
|------------|-----------|-----|
| Writing Fortran/MPI code | Claude Opus 4.5 | Highest SWE-bench, best debugging |
| PR code review | GPT-5.2 Codex xhigh | Methodical, precise bug-finding |
| Full-codebase analysis | Gemini 3 Pro | 1M token context, entire repo in one pass |
| Scientific literature review | Gemini 3 Pro | GPQA Diamond 93.8%, Deep Think reasoning |
| Architecture planning | Gemini 3 Pro + GPT-5.2 | Large context for analysis, abstract reasoning for design |
| CI / build verification | Any (automated) | Mechanical task, model quality less critical |

## Sources

- [Claude Opus 4.5 vs GPT-5.2 Codex: Best AI for Coding 2026](https://vertu.com/lifestyle/claude-opus-4-5-vs-gpt-5-2-codex-head-to-head-coding-benchmark-comparison/)
- [Claude 4.5 Opus vs. Gemini 3 Pro vs. GPT-5.2-codex-max (Composio)](https://composio.dev/blog/claude-4-5-opus-vs-gemini-3-pro-vs-gpt-5-codex-max-the-sota-coding-model)
- [New data on code quality: GPT-5.2 high, Opus 4.5, Gemini 3 (Sonar)](https://www.sonarsource.com/blog/new-data-on-code-quality-gpt-5-2-high-opus-4-5-gemini-3-and-more/)
- [Introducing GPT-5.2-Codex (OpenAI)](https://openai.com/index/introducing-gpt-5-2-codex/)
- [Building more with GPT-5.1-Codex-Max (OpenAI)](https://openai.com/index/gpt-5-1-codex-max/)
- [Evaluation of OpenAI Codex for HPC Parallel Programming (ICPP 2023)](https://arxiv.org/html/2306.15121)
- [GPT 5.2 is incredible - GitHub Discussion](https://github.com/openai/codex/discussions/8712)
- [Gemini 3 Pro — Google DeepMind](https://deepmind.google/models/gemini/pro/)
- [Gemini 3: Introducing the latest Gemini AI model (Google Blog)](https://blog.google/products-and-platforms/products/gemini/gemini-3/)
- [Testing Gemini 3.0 Pro's 1M Token Context Window](https://vertu.com/lifestyle/testing-gemini-3-0-pros-1-million-token-context-window/)
- [Gemini 3 Pro 1M Token Context Window Explained](https://www.sentisight.ai/gemini-1-million-token-context-window/)
