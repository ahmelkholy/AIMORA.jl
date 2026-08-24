export FixedStepReport, run_fixed_step, write_realtime_summary

struct FixedStepReport
    steps::Int
    dt_s::Float64
    duration_s::Float64
    mean_step_s::Float64
    max_step_s::Float64
    overruns::Int
    realtime::Bool
end

function run_fixed_step(step!::F; state = nothing, dt_s::Float64, duration_s::Float64, realtime::Bool = false, disable_gc::Bool = true) where {F}
    steps = Int(round(duration_s / dt_s))
    steps > 0 || error("Fixed-step run requires at least one step")
    total_elapsed = 0.0
    max_elapsed = 0.0
    overruns = 0
    old_gc = GC.enable(!disable_gc)
    t0 = time_ns()
    try
        for n in 1:steps
            target_ns = t0 + round(Int, n * dt_s * 1.0e9)
            t = (n - 1) * dt_s
            before = time_ns()
            step!(state, t, dt_s)
            after = time_ns()
            elapsed = (after - before) / 1.0e9
            total_elapsed += elapsed
            max_elapsed = max(max_elapsed, elapsed)
            elapsed > dt_s && (overruns += 1)
            if realtime
                while time_ns() < target_ns
                    # Legacy software pacing; the admitted absolute-clock owner is separate.
                end
            end
        end
    finally
        GC.enable(old_gc)
    end
    return FixedStepReport(
        steps,
        dt_s,
        steps * dt_s,
        total_elapsed / steps,
        max_elapsed,
        overruns,
        realtime,
    )
end

function write_realtime_summary(path::AbstractString, report::FixedStepReport)
    isdir(dirname(path)) || mkpath(dirname(path))
    open(path, "w") do io
        println(io, "{")
        @printf(io, "  \"steps\": %d,\n", report.steps)
        @printf(io, "  \"dt_s\": %.9f,\n", report.dt_s)
        @printf(io, "  \"duration_s\": %.9f,\n", report.duration_s)
        @printf(io, "  \"mean_step_s\": %.9f,\n", report.mean_step_s)
        @printf(io, "  \"max_step_s\": %.9f,\n", report.max_step_s)
        @printf(io, "  \"overruns\": %d,\n", report.overruns)
        @printf(io, "  \"realtime\": %s\n", report.realtime ? "true" : "false")
        println(io, "}")
    end
end
