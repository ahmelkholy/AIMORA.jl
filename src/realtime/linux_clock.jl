export monotonic_time_ns, wait_until_monotonic_ns

const _LINUX_CLOCK_MONOTONIC = Cint(1)
const _LINUX_TIMER_ABSTIME = Cint(1)
const _LINUX_EINTR = Cint(4)
const _NANOSECONDS_PER_SECOND = Int64(1_000_000_000)

struct _LinuxTimespec
    seconds::Clong
    nanoseconds::Clong
end

function _linux_timespec_from_ns(timestamp_ns::Int64)
    timestamp_ns >= 0 || throw(ArgumentError(
        "monotonic timestamp must be nonnegative",
    ))
    seconds, nanoseconds = divrem(timestamp_ns, _NANOSECONDS_PER_SECOND)
    return _LinuxTimespec(Clong(seconds), Clong(nanoseconds))
end

function monotonic_time_ns()
    Sys.islinux() || throw(ArgumentError(
        "the admitted absolute monotonic clock is available only on Linux",
    ))
    timestamp = Ref(_LinuxTimespec(0, 0))
    status = ccall(
        :clock_gettime,
        Cint,
        (Cint, Ref{_LinuxTimespec}),
        _LINUX_CLOCK_MONOTONIC,
        timestamp,
    )
    status == 0 || throw(SystemError("clock_gettime", true))
    value = timestamp[]
    return Base.Checked.checked_add(
        Base.Checked.checked_mul(Int64(value.seconds), _NANOSECONDS_PER_SECOND),
        Int64(value.nanoseconds),
    )
end

function wait_until_monotonic_ns(target_ns::Integer)
    Sys.islinux() || throw(ArgumentError(
        "the admitted absolute monotonic wait is available only on Linux",
    ))
    target = Int64(target_ns)
    deadline = Ref(_linux_timespec_from_ns(target))
    while true
        status = ccall(
            :clock_nanosleep,
            Cint,
            (Cint, Cint, Ref{_LinuxTimespec}, Ptr{Cvoid}),
            _LINUX_CLOCK_MONOTONIC,
            _LINUX_TIMER_ABSTIME,
            deadline,
            C_NULL,
        )
        status == 0 && return monotonic_time_ns()
        status == _LINUX_EINTR && continue
        Base.Libc.errno(status)
        throw(SystemError("clock_nanosleep", true))
    end
end
