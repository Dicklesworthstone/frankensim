use super::*;
use crate::pcm_wav::baffled::gas::GasPropagationLimits;
use fs_material::gas::GasSpec;

const RATE: u32 = 16_000;
const MICROPHONE: [f64; 3] = [0.0, 0.0, 0.2];
fn samples(count: usize) -> Vec<FlowSurfaceSample> {
    (0..count).map(|i| FlowSurfaceSample {
        position_m: [i as f64 * 1e-4, 0.0, 0.0], area_m2: 1e-8,
    }).collect()
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap()
}
fn basis(samples: &[FlowSurfaceSample]) -> Vec<SurfaceSample> {
    samples.iter().enumerate().map(|(i, tile)| {
        let mut shape = vec![0.0; samples.len()];
        shape[i] = 1.0 / tile.area_m2;
        SurfaceSample { position_m: tile.position_m, area_m2: tile.area_m2, mode_shape: shape }
    }).collect()
}
fn receiver(samples: &[FlowSurfaceSample]) -> BaffledTilePressure {
    BaffledTilePressure::new_differentiating(samples, RATE, MICROPHONE, (&gas()).into(),
        TilePressureLimits::surface_reference(), &mut || false).unwrap()
}
fn assert_accepted_same(a: &BaffledTilePressure, b: &BaffledTilePressure) {
    assert_eq!(a.started, b.started);
    for (a, b) in a.receivers.iter().zip(&b.receivers) {
        assert_eq!(a.head, b.head);
        assert_eq!(a.history, b.history);
        assert_eq!(a.taps, b.taps);
        assert_eq!(a.previous_velocity, b.previous_velocity);
    }
}

#[test]
fn one_piece_retains_original_pressure_bits_and_complete_volume_flow_basis() {
    let gas = gas();
    for count in [1, 37, MAX_RAYLEIGH_INPUTS] {
        let samples = samples(count);
        let mut bank = receiver(&samples);
        let mut original = BaffledPressure::new_differentiating(&basis(&samples),
            count, RATE, MICROPHONE, (&gas).into()).unwrap();
        assert_eq!(bank.receiver_count(), 1);
        assert_eq!(bank.delay_samples, original.delay_samples);
        assert_eq!(bank.filter_delay_samples(), original.filter_delay_samples());
        let mut live = false;
        for step in 0..320 {
            let flows: Vec<_> = (0..count).map(|i|
                1e-8 * ((step as f64 * 0.17) + i as f64 * 0.37).sin()).collect();
            let actual = bank.step_trace(&flows).unwrap();
            let expected = original.step_trace(&flows).unwrap();
            assert_eq!(actual.to_bits(), expected.to_bits());
            live |= actual.abs() > 1e-12;
        }
        assert!(live);
    }
}

#[test]
fn every_tile_beyond_the_modal_cap_keeps_its_actual_path_and_signed_flow() {
    let samples = samples(2 * MAX_RAYLEIGH_INPUTS + 1);
    let mut bank = receiver(&samples);
    let last = *samples.last().unwrap();
    let mut reference = BaffledPressure::new_differentiating(
        &[SurfaceSample { position_m: last.position_m, area_m2: last.area_m2,
            mode_shape: vec![1.0 / last.area_m2] }],
        1, RATE, MICROPHONE, (&gas()).into()).unwrap();
    reference.differentiate_pressure().unwrap();
    assert_eq!(bank.input_count(), samples.len());
    assert_eq!(bank.receiver_count(), 3);
    let mut flows = vec![0.0; samples.len()];
    let mut live = false;
    for step in 0..320 {
        *flows.last_mut().unwrap() = 1e-8 * (step as f64 * 0.19).sin();
        let actual = bank.step_trace(&flows).unwrap();
        let expected = reference.step_trace(&[*flows.last().unwrap()]).unwrap();
        assert!((actual - expected).abs() <= 1e-13 * expected.abs().max(1e-20));
        live |= actual.abs() > 1e-12;
    }
    assert!(live, "the last physical tile must not be truncated or averaged away");
}

#[test]
fn cancellation_invalid_tail_and_late_piece_overflow_retain_all_accepted_histories() {
    let mut bank = receiver(&samples(MAX_RAYLEIGH_INPUTS + 1));
    let mut flows = vec![1e-8; bank.input_count()];
    for _ in 0..128 { bank.step_trace(&flows).unwrap(); }
    let mut before = bank.clone();
    *flows.last_mut().unwrap() = f64::NAN;
    assert!(bank.step_trace(&flows).is_err());
    assert_accepted_same(&before, &bank);
    flows.fill(1e-8);
    let mut polls = 0;
    assert!(bank.step_trace_with_cancel(&flows, &mut || {
        polls += 1; polls == 2
    }).is_err());
    assert_accepted_same(&before, &bank);
    assert_eq!(bank.step_trace(&flows).unwrap().to_bits(),
        before.step_trace(&flows).unwrap().to_bits());

    let mut overflowed = false;
    *flows.last_mut().unwrap() = f64::MAX;
    for _ in 0..128 {
        let before = bank.clone();
        if bank.step_trace(&flows).is_err() {
            assert_accepted_same(&before, &bank);
            overflowed = true;
            break;
        }
    }
    assert!(overflowed, "finite extreme flow must refuse before any piece accepts");
}

#[test]
fn retained_storage_and_material_work_are_bounded_over_the_entire_surface() {
    let samples = samples(MAX_RAYLEIGH_INPUTS + 1);
    let limits = TilePressureLimits::surface_reference();
    let bank = receiver(&samples);
    let (histories, kernels) = bank.retained_values();
    for limits in [
        TilePressureLimits { maximum_tiles: samples.len() - 1, ..limits },
        TilePressureLimits { maximum_history_values: histories - 1, ..limits },
        TilePressureLimits { maximum_kernel_values: kernels - 1, ..limits },
    ] {
        assert!(BaffledTilePressure::new_differentiating(&samples, RATE, MICROPHONE,
            (&gas()).into(), limits, &mut || false).is_err());
    }
    let gas = gas();
    let propagation = GasPropagation::classical(GasPropagationLimits::surface_reference());
    let (_, first) = BaffledTilePressure::new_gas_differentiating(
        &samples[..MAX_RAYLEIGH_INPUTS], RATE, MICROPHONE,
        &gas, &propagation, limits, &mut || false).unwrap();
    for bounded in [
        GasPropagationLimits { max_unique_paths: first.filters.len(), ..propagation.limits },
        GasPropagationLimits { max_filter_taps: first.filters.iter().map(|p| p.taps).sum(),
            ..propagation.limits },
        GasPropagationLimits { max_kernel_terms: first.kernel_terms, ..propagation.limits },
    ] {
        let mut limited = propagation.clone();
        limited.limits = bounded;
        assert!(BaffledTilePressure::new_gas_differentiating(
            &samples, RATE, MICROPHONE, &gas, &limited, limits, &mut || false).is_err());
    }
    let (mut bank, mut report) = BaffledTilePressure::new_gas_differentiating(
        &samples, RATE, MICROPHONE, &gas, &propagation, limits, &mut || false).unwrap();
    assert_eq!(report.gas, gas);
    assert_eq!(report.law, propagation.law);
    assert_eq!(report.surface_filter_indices.len(), samples.len());
    for (tile, &index) in samples.iter().zip(&report.surface_filter_indices) {
        let path = &report.filters[index];
        let distance = fs_math::det::sqrt(tile.position_m[0].powi(2) + MICROPHONE[2].powi(2));
        assert_eq!(path.range_m.to_bits(), distance.to_bits());
    }
    let support = bank.delay_samples;
    let filter_delay = bank.filter_delay_samples();
    let latency = report.additional_latency_samples + 3;
    bank.align_gas_latency(&mut report, latency).unwrap();
    assert_eq!(bank.delay_samples, (support.0 + 3, support.1 + 3));
    assert_eq!(bank.filter_delay_samples(), filter_delay + 3.0);
    assert_eq!(report.additional_latency_samples, latency);
    assert!(bank.align_gas_latency(&mut report, latency - 1).is_err());
}

#[test]
fn one_piece_material_paths_and_pressure_match_the_original_owner() {
    let samples = samples(37);
    let gas = gas();
    let propagation = GasPropagation::classical(GasPropagationLimits::surface_reference());
    let (mut bank, report) = BaffledTilePressure::new_gas_differentiating(
        &samples, RATE, MICROPHONE, &gas, &propagation,
        TilePressureLimits::surface_reference(), &mut || false).unwrap();
    let (mut original, expected) = BaffledPressure::new_multirate_gas_differentiating(
        &basis(&samples), samples.len(), RATE, RATE, MICROPHONE,
        &gas, &propagation, &mut || false).unwrap();
    assert_eq!(report.surface_filter_indices, expected.surface_filter_indices);
    assert_eq!(report.additional_latency_samples, expected.additional_latency_samples);
    assert_eq!(report.kernel_terms, expected.kernel_terms);
    assert_eq!(bank.delay_samples, original.delay_samples);
    for step in 0..320 {
        let flows: Vec<_> = (0..samples.len()).map(|i|
            1e-8 * (step as f64 * 0.11 + i as f64 * 0.13).sin()).collect();
        assert_eq!(bank.step_trace(&flows).unwrap().to_bits(),
            original.step_trace(&flows).unwrap().to_bits());
    }
}
