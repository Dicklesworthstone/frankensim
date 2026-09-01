//! Lossless distributed-line relations for resonant electrical conductors.

/// Parameters for an electrically long secondary conductor.
#[derive(Debug, Clone)]
pub struct QuarterWaveParams {
    /// Electrical-disturbance frequency in hertz.
    pub frequency_hz: f64,
    /// Propagation speed along the circuit in metres per second.
    pub propagation_speed_mps: f64,
    /// Developed conductor length in metres.
    pub conductor_length_m: f64,
}

/// Source-computable quarter-wave geometry. The normalized profile is not an
/// absolute voltage prediction: that requires excitation, impedance, loss,
/// loading, and coupling data that a topology-only source may not provide.
#[derive(Debug, Clone)]
pub struct QuarterWaveResult {
    /// Distributed wavelength in metres.
    pub wavelength_m: f64,
    /// One quarter of the distributed wavelength in metres.
    pub quarter_wave_length_m: f64,
    /// Phase accumulated along the developed conductor in radians.
    pub electrical_length_rad: f64,
    /// Signed phase error relative to exactly pi/2 radians.
    pub quarter_wave_error_rad: f64,
    /// Signed conductor-length error relative to one quarter wavelength.
    pub length_error_m: f64,
    /// Developed conductor length divided by the quarter-wave target.
    pub length_ratio: f64,
    /// Normalized open-end standing-wave profile; never an absolute voltage.
    pub remote_terminal_profile_fraction: f64,
}

/// Compute the distributed-wave coordinates without inventing a lumped coil,
/// quality factor, breakdown law, or voltage gain.
#[must_use]
pub fn step_quarter_wave(params: &QuarterWaveParams) -> QuarterWaveResult {
    let wavelength_m = params.propagation_speed_mps / params.frequency_hz;
    let quarter_wave_length_m = wavelength_m / 4.0;
    let electrical_length_rad =
        2.0 * std::f64::consts::PI * params.conductor_length_m / wavelength_m;

    QuarterWaveResult {
        wavelength_m,
        quarter_wave_length_m,
        electrical_length_rad,
        quarter_wave_error_rad: electrical_length_rad - std::f64::consts::FRAC_PI_2,
        length_error_m: params.conductor_length_m - quarter_wave_length_m,
        length_ratio: params.conductor_length_m / quarter_wave_length_m,
        remote_terminal_profile_fraction: electrical_length_rad.sin().abs(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn teslas_printed_example_is_exactly_quarter_wave() {
        let result = step_quarter_wave(&QuarterWaveParams {
            frequency_hz: 925.0,
            propagation_speed_mps: 185_000.0 * 1_609.344,
            conductor_length_m: 50.0 * 1_609.344,
        });
        assert!((result.quarter_wave_length_m - 80_467.2).abs() < 1.0e-8);
        assert!(result.length_error_m.abs() < 1.0e-8);
        assert!((result.electrical_length_rad - std::f64::consts::FRAC_PI_2).abs() < 1.0e-12);
        assert!((result.remote_terminal_profile_fraction - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn frequency_changes_required_wire_length_without_changing_the_law() {
        let low = step_quarter_wave(&QuarterWaveParams {
            frequency_hz: 500.0,
            propagation_speed_mps: 185_000.0 * 1_609.344,
            conductor_length_m: 50.0 * 1_609.344,
        });
        let high = step_quarter_wave(&QuarterWaveParams {
            frequency_hz: 1_000.0,
            propagation_speed_mps: 185_000.0 * 1_609.344,
            conductor_length_m: 50.0 * 1_609.344,
        });
        assert!((low.quarter_wave_length_m / high.quarter_wave_length_m - 2.0).abs() < 1.0e-12);
    }
}
