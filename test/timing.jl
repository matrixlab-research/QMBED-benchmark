mutable struct VerificationTiming
    path::String
    scope::String
    status::String
    wall_seconds::Float64
    compile_seconds::Float64
    recompile_seconds::Float64
    gc_seconds::Float64
    allocated_bytes::Int
    lock_conflicts::Int
end

const VERIFICATION_TIMINGS = VerificationTiming[]
const VERIFICATION_SUITE_START_NS = time_ns()
const VERIFICATION_TIMING_OUTPUT =
    get(ENV, "QUSPIN_TEST_TIMING_OUTPUT", "")

function timed_include(path::AbstractString, scope::AbstractString)
    started_ns = time_ns()
    try
        measurement = @timed include(path)
        push!(
            VERIFICATION_TIMINGS,
            VerificationTiming(
                String(path),
                String(scope),
                "passed",
                measurement.time,
                measurement.compile_time,
                measurement.recompile_time,
                measurement.gctime,
                measurement.bytes,
                measurement.lock_conflicts,
            ),
        )
        return measurement.value
    catch
        push!(
            VERIFICATION_TIMINGS,
            VerificationTiming(
                String(path),
                String(scope),
                "failed",
                (time_ns() - started_ns) / 1.0e9,
                NaN,
                NaN,
                NaN,
                0,
                0,
            ),
        )
        rethrow()
    end
end

function _csv_cell(value)
    text = string(value)
    return occursin(r"[\",\n]", text) ?
        "\"" * replace(text, "\"" => "\"\"") * "\"" :
        text
end

function write_verification_timings()
    isempty(VERIFICATION_TIMING_OUTPUT) && return
    mkpath(dirname(VERIFICATION_TIMING_OUTPUT))
    total_wall = (time_ns() - VERIFICATION_SUITE_START_NS) / 1.0e9
    open(VERIFICATION_TIMING_OUTPUT, "w") do io
        println(
            io,
            "path,scope,status,wall_seconds,compile_seconds,recompile_seconds," *
            "gc_seconds,allocated_bytes,lock_conflicts",
        )
        for timing in VERIFICATION_TIMINGS
            values = (
                timing.path,
                timing.scope,
                timing.status,
                timing.wall_seconds,
                timing.compile_seconds,
                timing.recompile_seconds,
                timing.gc_seconds,
                timing.allocated_bytes,
                timing.lock_conflicts,
            )
            println(io, join(_csv_cell.(values), ","))
        end
        total_values = (
            "__suite__",
            "all",
            all(timing -> timing.status == "passed", VERIFICATION_TIMINGS) ?
                "passed" :
                "failed",
            total_wall,
            sum(timing.compile_seconds for timing in VERIFICATION_TIMINGS),
            sum(timing.recompile_seconds for timing in VERIFICATION_TIMINGS),
            sum(timing.gc_seconds for timing in VERIFICATION_TIMINGS),
            sum(timing.allocated_bytes for timing in VERIFICATION_TIMINGS),
            sum(timing.lock_conflicts for timing in VERIFICATION_TIMINGS),
        )
        println(io, join(_csv_cell.(total_values), ","))
    end
end

atexit(write_verification_timings)
