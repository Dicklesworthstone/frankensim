//! Analytic point-source receiver potentials drive the actual derivative,
//! physical delay, rate converter and ensemble. They are boundary data for an
//! observation test, not a replacement mechanical engine.
use super::*;
use crate::pcm_wav::decimate::Decimator;
use core::f64::consts::TAU;

const RATE: u32 = 48_000;
const SAMPLES: u64 = 1024;
const BLOCK: usize = 97;

struct HarmonicPressure {
    rate: u32,
    maximum_samples: u64,
    completed: u64,
    amplitude_pa: f64,
    omega: f64,
    derivative: Decimator,
    propagation: FractionalDelay,
}
impl HarmonicPressure {
    fn new(ratio: u32, amplitude_pa: f64, fraction: f64, physical_travel: f64) -> Self {
        let rate = RATE * ratio;
        Self { rate, maximum_samples: SAMPLES * u64::from(ratio), completed: 0,
            amplitude_pa, omega: TAU * fraction * f64::from(RATE),
            derivative: Decimator::new_differentiating(1, 1, rate).unwrap(),
            propagation: FractionalDelay::new(physical_travel * f64::from(ratio), 1024).unwrap() }
    }
}
impl PressureRenderer for HarmonicPressure {
    fn samples_rendered(&self) -> u64 { self.completed }
    fn max_block_len(&self) -> usize { BLOCK }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate == self.rate { Ok(()) } else { Err(sizing("analytic pressure source clock")) }
    }
    fn validate_sample_count(&self, count: u64) -> Result<(), RenderError> {
        if count <= self.maximum_samples - self.completed { Ok(()) }
        else { Err(sizing("analytic pressure source horizon")) }
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > BLOCK { return Err(sizing("analytic source callback capacity")); }
        self.validate_sample_count(output.len() as u64)?;
        for pressure in output {
            // A monopole's receiver potential is rho*Q/(4*pi*r) [Pa s].
            // This supplied sinusoidal potential differentiates to A*cos(wt)
            // pascals; its physical travel remains independent of filter latency.
            let time = (self.completed + 1) as f64 / f64::from(self.rate);
            let potential = self.amplitude_pa / self.omega * (self.omega * time).sin();
            let emitted = self.derivative.preview(&[potential]).map_err(sizing)?[0];
            *pressure = self.propagation.step(emitted);
            self.derivative.commit();
            self.completed += 1;
        }
        Ok(())
    }
}
fn config() -> PressureEnsembleConfig {
    PressureEnsembleConfig { sample_rate_hz: RATE, samples: SAMPLES, max_block: BLOCK, max_parts: 4 }
}
fn part(ratio: u32, amplitude: f64, frequency: f64, travel: f64) -> DecimatedRenderer<HarmonicPressure> {
    DecimatedRenderer::new(HarmonicPressure::new(ratio, amplitude, frequency, travel),
        RATE * ratio, RATE, BLOCK).unwrap()
}
fn parts(frequency: f64) -> Vec<DecimatedRenderer<HarmonicPressure>> {
    vec![part(3, 0.75, frequency, 0.125), part(6, -0.25, frequency, 0.875)]
}
fn render(mut ensemble: PressureEnsemble<HarmonicPressure>, partition: usize)
    -> (PressureEnsemble<HarmonicPressure>, Vec<f64>)
{
    let mut output = vec![0.0; SAMPLES as usize];
    for block in output.chunks_mut(partition) { ensemble.block(block).unwrap(); }
    (ensemble, output)
}

#[test]
fn three_and_six_times_clocks_align_actual_derivative_pressure_and_original_travel() {
    for frequency in [0.05, 0.2, 0.4, 0.45] {
        assert!(PressureEnsemble::new_with_source_delays(parts(frequency), config(),
            &[80.0 / 3.0, 80.0 / 6.0]).is_err());
        let mut reference = Vec::new();
        for partition in [1, 31, BLOCK] {
            let ensemble = PressureEnsemble::new_with_bandlimited_source_delays(
                parts(frequency), config(), &[80.0 / 3.0, 80.0 / 6.0],
                0.45 * f64::from(RATE)).unwrap();
            let common = ensemble.delay_output_samples();
            assert_eq!(common, 139); // ceil(max(80 + 80/3, 44 + 80/6)) + 32
            assert_eq!(ensemble.alignment_band_limit_hz(), Some(0.45 * f64::from(RATE)));
            for (info, ratio) in ensemble.part_info().iter().zip([3, 6]) {
                assert_eq!(info.observation.ratio, ratio);
                assert_eq!(info.alignment_filter_delay_samples, 32);
                let source = info.source_filter_delay_output_samples as f64
                    + info.source_filter_delay_fraction_samples;
                assert!((source - 80.0 / ratio as f64).abs() < 1e-13);
                let latency = info.total_filter_delay_output_samples as f64
                    + info.source_filter_delay_fraction_samples
                    + info.alignment_delay_samples as f64 + info.alignment_delay_fraction_samples
                    + info.alignment_filter_delay_samples as f64;
                assert!((latency - common as f64).abs() < 1e-13);
            }
            let (ensemble, output) = render(ensemble, partition);
            if reference.is_empty() { reference = output.clone(); }
            else { assert_eq!(output, reference, "frequency {frequency}, partition {partition}"); }
            let (mut error, mut signal) = (0.0, 0.0);
            for (n, &actual) in output.iter().enumerate().skip(384) {
                let exact: f64 = [(0.75, 0.125), (-0.25, 0.875)].into_iter().map(|(a, travel)|
                    a * (TAU * frequency * ((n + 1) as f64 - common as f64 - travel)).cos()).sum();
                error += (actual - exact).powi(2);
                signal += exact * exact;
            }
            let relative = (error / signal).sqrt();
            assert!(relative < 5e-4, "frequency {frequency}, relative waveform error {relative}");
            if partition == 1 { eprintln!("3x/6x physical pressure: fraction={frequency}, relative_error={relative:e}, retained_delay={common}"); }
            assert_eq!(ensemble.samples_rendered(), SAMPLES);
            assert_eq!(ensemble.remaining_samples(), 0);
            for (part, ratio) in ensemble.parts().iter().zip([3_u64, 6]) {
                assert_eq!(part.samples_rendered(), SAMPLES);
                assert_eq!(part.source().samples_rendered(), SAMPLES * ratio);
            }
        }
    }
}

#[test]
fn integral_parts_share_the_common_reconstruction_in_a_fractional_ensemble() {
    let frequency = 0.4;
    let ensemble = PressureEnsemble::new_with_bandlimited_source_delays(
        vec![part(1, 0.5, frequency, 0.125), part(3, 0.5, frequency, 0.125)],
        config(), &[80.0, 80.0 / 3.0], 0.45 * f64::from(RATE)).unwrap();
    assert_eq!(ensemble.delay_output_samples(), 139);
    let integral = ensemble.part_info()[0];
    assert_eq!(integral.source_filter_delay_fraction_samples, 0.0);
    assert_eq!(integral.alignment_delay_samples, 27);
    assert_eq!(integral.alignment_delay_fraction_samples, 0.0);
    assert_eq!(integral.alignment_filter_delay_samples, 32);
    let (_, pressure) = render(ensemble, BLOCK);
    let (mut error, mut signal) = (0.0, 0.0);
    for (n, &actual) in pressure.iter().enumerate().skip(384) {
        let exact = (TAU * frequency * ((n + 1) as f64 - 139.0 - 0.125)).cos();
        error += (actual - exact).powi(2);
        signal += exact * exact;
    }
    assert!((error / signal).sqrt() < 5e-4);
}

#[test]
fn fractional_alignment_preserves_weak_pressure_between_cancelling_physical_sources() {
    let epsilon = 2.0_f64.powi(-60);
    let amplitudes = [1.0, epsilon, -1.0];
    for order in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
        let sources = order.map(|i| part(3, amplitudes[i], 0.4, 0.125)).into_iter().collect();
        let ensemble = PressureEnsemble::new_with_bandlimited_source_delays(
            sources, config(), &[80.0 / 3.0; 3], 0.45 * f64::from(RATE)).unwrap();
        let common = ensemble.delay_output_samples();
        let (_, pressure) = render(ensemble, 31);
        let (mut error, mut signal) = (0.0, 0.0);
        for (n, &actual) in pressure.iter().enumerate().skip(384) {
            let exact = epsilon * (TAU * 0.4 * ((n + 1) as f64 - common as f64 - 0.125)).cos();
            error += (actual - exact).powi(2);
            signal += exact * exact;
        }
        assert!((error / signal).sqrt() < 5e-4, "lost weak coherent pressure in source order {order:?}");
    }
}

#[test]
fn admitted_integral_alignment_retains_the_original_samples_and_latency() {
    let make = || vec![part(1, 0.75, 0.2, 0.125), part(1, -0.25, 0.2, 0.875)];
    let original = PressureEnsemble::new_with_source_delays(make(), config(), &[80.0, 80.0]).unwrap();
    let bandlimited = PressureEnsemble::new_with_bandlimited_source_delays(
        make(), config(), &[80.0, 80.0], 0.45 * f64::from(RATE)).unwrap();
    assert_eq!(original.delay_output_samples(), bandlimited.delay_output_samples());
    assert_eq!(original.part_info(), bandlimited.part_info());
    assert!(bandlimited.part_info().iter().all(|p| p.alignment_filter_delay_samples == 0));
    let (_, before) = render(original, 1);
    let (_, after) = render(bandlimited, BLOCK);
    assert_eq!(before.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        after.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
}

#[test]
fn fractional_band_clock_and_complete_history_are_admitted_before_observation() {
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY, 0.450001 * f64::from(RATE)] {
        assert!(PressureEnsemble::new_with_bandlimited_source_delays(parts(0.2), config(),
            &[80.0 / 3.0, 80.0 / 6.0], frequency).is_err());
    }
    for delays in [vec![0.25], vec![-0.25, 0.0], vec![f64::NAN, 0.0], vec![0.0, f64::INFINITY]] {
        assert!(PressureEnsemble::new_with_bandlimited_source_delays(parts(0.2), config(),
            &delays, 0.45 * f64::from(RATE)).is_err());
    }
    // Bound final latency AND the full FIR support, not just its center delay.
    let make = || vec![part(1, 0.75, 0.2, 0.125), part(1, -0.25, 0.2, 0.875)];
    let cap = MAX_TOTAL_FILTER_DELAY_OUTPUT_SAMPLES as f64;
    for delays in [[0.0, cap - 31.5], [0.0, cap - 64.5]] {
        assert!(PressureEnsemble::new_with_bandlimited_source_delays(make(), config(),
            &delays, 0.45 * f64::from(RATE)).is_err());
    }
    let mut admitted = PressureEnsemble::new_with_bandlimited_source_delays(make(), config(),
        &[0.25, cap - 64.5], 0.45 * f64::from(RATE)).unwrap();
    assert_eq!(admitted.delay_output_samples(), MAX_TOTAL_FILTER_DELAY_OUTPUT_SAMPLES - 32);
    let mut unchanged = [17.0; BLOCK + 1];
    assert!(admitted.block(&mut unchanged).is_err());
    assert_eq!(unchanged, [17.0; BLOCK + 1]);
    assert_eq!(admitted.samples_rendered(), 0);
    assert!(admitted.parts().iter().all(|p| p.source().samples_rendered() == 0));
}
