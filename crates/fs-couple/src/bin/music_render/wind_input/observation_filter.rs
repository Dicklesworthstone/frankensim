//! The admitted observation filter and its numerical latency, before radiation.

use fs_couple::pcm_wav::observation::resample::BandlimitedResampler;

pub(super) fn provenance(
    filter: &BandlimitedResampler,
    physical_intervals: u64,
    output_samples: u64,
) -> String {
    let policy = filter.policy();
    let report = filter.report();
    format!(
        concat!(
            ",\"time_projection\":{{",
            "\"physical_interval_s\":{dt:e},\"declared_physical_intervals\":{physical_intervals},",
            "\"declared_physical_duration_s\":{physical_duration:e},",
            "\"observation_rate_hz\":{rate},\"declared_observation_samples\":{output_samples},",
            "\"declared_observation_duration_s\":{output_duration:e},",
            "\"method\":\"causal windowed-sinc reconstruction of original interval means\",",
            "\"boundary_control\":\"declared hold or linear ramp sampled at each original physical interval start\",",
            "\"bandlimited_antialias_filter\":true,",
            "\"filter\":{{\"window\":\"four-term Blackman-Harris\",",
            "\"passband_hz\":{pass:e},\"stopband_hz\":{stop:e},",
            "\"channels\":{channels},\"taps\":{taps},\"half_support_input_frames\":{half},",
            "\"phase_count\":{phases},\"phase_interpolation\":\"linear\",",
            "\"phase_spacing_input_samples\":{phase_spacing:e},",
            "\"delay_s\":{delay:e},\"delay_observation_samples\":{delay_samples:e},",
            "\"history_frames\":{history_frames},\"history_values\":{history_values},",
            "\"preparation_terms\":{actual_preparation},\"output_terms\":{actual_work},",
            "\"max_preparation_terms\":{preparation},\"max_history_values\":{history},",
            "\"max_output_terms\":{work},\"initial_history\":\"zero\",",
            "\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}},",
            "\"scope\":\"pressure, DG terminal flow and selected reservoir tile flow share one filter before radiation; interval-mean samples retain their original averaging; finite FIR and phase-table approximation; source-clock aliasing is not recovered; physical source and control clocks are unchanged\"}}"
        ),
        dt = report.input_period_s,
        physical_intervals = physical_intervals,
        physical_duration = report.input_period_s * physical_intervals as f64,
        rate = report.output_rate_hz,
        output_samples = output_samples,
        output_duration = output_samples as f64 / f64::from(report.output_rate_hz),
        pass = policy.passband_hz,
        stop = policy.stopband_hz,
        channels = report.channels,
        taps = report.taps,
        half = report.half_support,
        phases = report.phase_count,
        phase_spacing = report.phase_spacing_input_samples,
        delay = report.latency_s,
        delay_samples = report.latency_output_samples,
        history_frames = report.history_frames,
        history_values = report.history_values,
        actual_preparation = report.preparation_terms,
        actual_work = report.output_terms,
        preparation = policy.max_preparation_terms,
        history = policy.max_history_values,
        work = policy.max_output_terms,
    )
}
