#!/usr/bin/env julia

using TOML

const ROOT = @__DIR__

fail(message) = error("AIMORA package check failed: $(message)")

const REQUIRED_PATHS = (
    "Project.toml",
    "README.md",
    "src/AIMORA.jl",
    "src/core",
    "src/io",
    "src/models",
    "src/solver_api/backend.jl",
    "src/solver_api/worker.jl",
    "src/studies",
    "test/runtests.jl",
    "test/study_worker_contracts.jl",
)

for relative_path in REQUIRED_PATHS
    ispath(joinpath(ROOT, relative_path)) || fail("missing $(relative_path)")
end

project = TOML.parsefile(joinpath(ROOT, "Project.toml"))
get(project, "name", nothing) == "AIMORA" || fail("Project.toml name is not AIMORA")

for forbidden_path in (
    "AGENTS.md",
    "MEMORY.md",
    "MAP.md",
    "LEDGER.md",
    "LEDGER",
    "TRANSLATION_MAP.md",
    "TRANSLATION_LEDGER.md",
    ".codex",
    "runs",
    "validation",
    "examples",
    "docs",
)
    !ispath(joinpath(ROOT, forbidden_path)) ||
        fail("internal development path is present: $(forbidden_path)")
end

tracked_solver = read(`git -C $ROOT ls-files --stage src/solvers`, String)
isempty(tracked_solver) || fail("public repository retains a nested solver Gitlink")

if isfile(joinpath(ROOT, ".gitmodules"))
    gitmodules = read(joinpath(ROOT, ".gitmodules"), String)
    !contains(gitmodules, "src/solvers") ||
        fail("public repository metadata retains the private solver path")
end

println("AIMORA package boundary check passed")
