export RealtimeTimingSample,
       RealtimeTimingSummary,
       realtime_release_ns,
       realtime_timing_sample,
       summarize_realtime_timing

struct RealtimeTimingSample
    step::Int
    release_ns::Int64
    start_ns::Int64
    completion_ns::Int64
    deadline_ns::Int64
    jitter_ns::Int64
    computation_ns::Int64
    response_ns::Int64
    slack_ns::Int64
    overrun::Bool
end

function realtime_release_ns(epoch_ns::Integer, period_ns::Integer, step::Integer)
    step >= 0 || throw(ArgumentError("real-time step index must be nonnegative"))
    period_ns > 0 || throw(ArgumentError("real-time period must be positive"))
    return Base.Checked.checked_add(
        Int64(epoch_ns),
        Base.Checked.checked_mul(Int64(period_ns), Int64(step)),
    )
end

function realtime_timing_sample(
    step::Integer,
    release_ns::Integer,
    start_ns::Integer,
    completion_ns::Integer,
    period_ns::Integer,
)
    index = Int(step)
    index >= 0 || throw(ArgumentError("real-time step index must be nonnegative"))
    release = Int64(release_ns)
    start = Int64(start_ns)
    completion = Int64(completion_ns)
    period = Int64(period_ns)
    period > 0 || throw(ArgumentError("real-time period must be positive"))
    start >= release || throw(ArgumentError(
        "real-time start cannot precede its monotonic release",
    ))
    completion >= start || throw(ArgumentError(
        "real-time completion cannot precede its start",
    ))
    deadline = Base.Checked.checked_add(release, period)
    jitter = start - release
    computation = completion - start
    response = completion - release
    slack = deadline - completion
    return RealtimeTimingSample(
        index,
        release,
        start,
        completion,
        deadline,
        jitter,
        computation,
        response,
        slack,
        response > period,
    )
end

struct RealtimeTimingSummary
    samples::Int
    minimum_response_ns::Int64
    median_response_ns::Float64
    response_percentile_90_ns::Float64
    response_percentile_99_ns::Float64
    response_percentile_99_9_ns::Float64
    maximum_response_ns::Int64
    maximum_jitter_ns::Int64
    overruns::Int
end

function summarize_realtime_timing(samples::AbstractVector{RealtimeTimingSample})
    isempty(samples) && throw(ArgumentError(
        "real-time timing summary requires at least one sample",
    ))
    steps = getfield.(samples, :step)
    issorted(steps) && length(unique(steps)) == length(steps) || throw(ArgumentError(
        "real-time timing samples require unique increasing step indices",
    ))
    responses = getfield.(samples, :response_ns)
    jitters = getfield.(samples, :jitter_ns)
    return RealtimeTimingSummary(
        length(samples),
        minimum(responses),
        median(responses),
        quantile(responses, 0.90),
        quantile(responses, 0.99),
        quantile(responses, 0.999),
        maximum(responses),
        maximum(jitters),
        count(getfield.(samples, :overrun)),
    )
end
