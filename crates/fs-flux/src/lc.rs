#[derive(Debug, Clone)]
pub struct TeslaCoilParams {
    pub resonant_freq_khz: f64,
    pub input_kv: f64,
    pub spark_gap_mm: f64,
    pub q_factor: f64,
}

#[derive(Debug, Clone)]
pub struct TeslaCoilResult {
    pub resonant_freq_khz: f64,
    pub secondary_potential_mv: f64,
    pub streamer_length_inches: f64,
    pub streamer_length_meters: f64,
}

pub fn step_tesla_coil(params: &TeslaCoilParams) -> TeslaCoilResult {
    let primary_l = 0.012; // mH
    let secondary_l = 85.0; // mH
    let transformation_ratio = (secondary_l / primary_l as f64).sqrt();
    let k = 0.18;
    let secondary_potential_mv =
        ((params.input_kv * transformation_ratio * k * params.q_factor.sqrt()) / 1000.0)
            * (params.spark_gap_mm / 15.0);
    let streamer_length_inches = secondary_potential_mv * 28.0;

    TeslaCoilResult {
        resonant_freq_khz: params.resonant_freq_khz,
        secondary_potential_mv: (secondary_potential_mv * 100.0).round() / 100.0,
        streamer_length_inches: (streamer_length_inches * 10.0).round() / 10.0,
        streamer_length_meters: (((streamer_length_inches * 2.54) / 100.0) * 100.0).round() / 100.0,
    }
}
