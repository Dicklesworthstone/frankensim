//! Source law -> individual paths -> actual shared pressure convolution.
//! Transport substitutions are manufactured inputs, not measured gas data.
use super::*;
use fs_material::gas::{ConductivityModel, GasSpec};

fn gas(viscosity_factor: f64) -> GasState {
    let mut source = GasSpec::dry_air_ussa1976();
    source.sutherland_beta *= viscosity_factor;
    source.conductivity = ConductivityModel::Eucken;
    GasState::try_new(&source, 293.15, 101_325.0).unwrap()
}
fn surface() -> Vec<SurfaceSample> {
    vec![
        SurfaceSample {
            position_m: [0.0; 3],
            area_m2: 0.2,
            mode_shape: vec![1.0, 0.2],
        },
        SurfaceSample {
            position_m: [0.8, 0.4, 0.0],
            area_m2: 0.1,
            mode_shape: vec![-0.7, 0.4],
        },
    ]
}
fn controls() -> GasPropagation {
    GasPropagation::classical(GasPropagationLimits::surface_reference())
}
fn complex(taps: &[f64], omega: f64) -> [f64; 2] {
    taps.iter()
        .enumerate()
        .fold([0.0; 2], |[re, im], (n, &tap)| {
            [
                re + tap * (omega * n as f64).cos(),
                im - tap * (omega * n as f64).sin(),
            ]
        })
}
fn kernel(receiver: &BaffledPressure, input: usize, omega: f64) -> [f64; 2] {
    receiver
        .kernels
        .iter()
        .fold([0.0; 2], |[re, im], (n, weights)| {
            [
                re + weights[input] * (omega * *n as f64).cos(),
                im - weights[input] * (omega * *n as f64).sin(),
            ]
        })
}

#[test]
fn every_signed_surface_path_uses_its_own_material_transfer_and_travel_time() {
    let gas = gas(100.0);
    let surface = surface();
    let position = [0.2, 0.1, 1.0];
    let rate = 48_000;
    let (receiver, report) = BaffledPressure::new_multirate_gas_differentiating(
        &surface,
        2,
        rate,
        rate,
        position,
        &gas,
        &controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.surface_filter_indices, [0, 1]);
    assert_ne!(report.filters[0].range_m, report.filters[1].range_m);
    for frequency in [100.0, 4000.0, 12_000.0, 21_000.0] {
        let omega = TAU * frequency / f64::from(rate);
        for input in 0..2 {
            let (mut expected, mut scale) = ([0.0; 2], 0.0);
            for (sample, &index) in surface.iter().zip(&report.surface_filter_indices) {
                let receipt = report.filters[index];
                let filter =
                    ObserverFilter::for_law(report.law, &gas, None, receipt.range_m, rate).unwrap();
                let h = complex(filter.taps(), omega);
                // Independent complex Rayleigh path, before summing surfaces.
                let phase = -omega
                    * (receipt.receiver_delay_s * f64::from(rate)
                        + report.additional_latency_samples as f64
                        + path::HALF as f64);
                let gain = gas.density * sample.area_m2 * sample.mode_shape[input]
                    / (TAU * receipt.range_m);
                expected[0] += gain * (h[0] * phase.cos() - h[1] * phase.sin());
                expected[1] += gain * (h[0] * phase.sin() + h[1] * phase.cos());
                scale += gain.abs() * h[0].hypot(h[1]);
                let analytic =
                    (-gas.stokes_kirchhoff_absorption(TAU * frequency) * receipt.range_m).exp();
                assert!((h[0].hypot(h[1]) - analytic).abs() < 2e-3);
            }
            let actual = kernel(&receiver, input, omega);
            assert!((actual[0] - expected[0]).hypot(actual[1] - expected[1]) < 5e-4 * scale);
        }
    }
}

#[test]
fn gas_transport_data_changes_real_pressure_while_geometry_density_and_speed_stay_fixed() {
    let gases = [gas(1.0), gas(100.0)];
    assert_eq!(gases[0].density, gases[1].density);
    assert_eq!(gases[0].sound_speed, gases[1].sound_speed);
    let patch = [SurfaceSample {
        position_m: [0.0; 3],
        area_m2: 0.2,
        mode_shape: vec![1.0],
    }];
    let mut receivers: Vec<_> = gases
        .iter()
        .map(|gas| {
            BaffledPressure::new_multirate_gas_differentiating(
                &patch,
                1,
                48_000,
                48_000,
                [0.0, 0.0, 2.0],
                gas,
                &controls(),
                &mut || false,
            )
            .unwrap()
            .0
        })
        .collect();
    let mut energy = [0.0; 2];
    for frame in 0..4096 {
        let velocity = [(TAU * 12_000.0 * frame as f64 / 48_000.0).sin()];
        for (i, receiver) in receivers.iter_mut().enumerate() {
            let pressure = receiver.step(&velocity).unwrap();
            if frame >= 2048 {
                energy[i] += pressure * pressure;
            }
        }
    }
    let alpha = gases.map(|gas| gas.stokes_kirchhoff_absorption(TAU * 12_000.0));
    let ratio = (-2.0 * (alpha[1] - alpha[0])).exp();
    assert!(energy[0] > 1.0 && energy[1] < 0.95 * energy[0]);
    assert!(((energy[1] / energy[0]).sqrt() - ratio).abs() < 2e-3);
    let mut copy = receivers[1].clone();
    for frame in 0..128 {
        let velocity = [(0.13 * frame as f64).cos()];
        assert_eq!(
            copy.step(&velocity).unwrap(),
            receivers[1].step(&velocity).unwrap()
        );
    }
}

#[test]
fn lossless_material_path_matches_original_shared_receiver() {
    let gas = gas(1.0);
    let surface = surface();
    let position = [0.2, 0.1, 1.0];
    let mut selection = controls();
    selection.law = GasAbsorptionModel::Lossless;
    let (mut material, report) = BaffledPressure::new_multirate_gas_differentiating(
        &surface,
        2,
        48_000,
        48_000,
        position,
        &gas,
        &selection,
        &mut || false,
    )
    .unwrap();
    let mut original = BaffledPressure::new_differentiating(
        &surface,
        2,
        48_000,
        position,
        RayleighMedium::from(&gas),
    )
    .unwrap();
    assert_eq!(report.additional_latency_samples, 0);
    assert_eq!(
        material.filter_delay_samples(),
        original.filter_delay_samples()
    );
    for frame in 0..1024 {
        let velocity = [(0.17 * frame as f64).sin(), (0.31 * frame as f64).cos()];
        let (a, b) = (
            material.step(&velocity).unwrap(),
            original.step(&velocity).unwrap(),
        );
        assert!((a - b).abs() < 1e-11 * b.abs().max(1.0));
    }
}

#[test]
fn common_padding_preserves_near_receiver_physical_delay_without_future_samples() {
    let (gas, air) =
        fs_material::iso9613::TerrestrialAir::moist_air(293.15, 101_325.0, 0.5).unwrap();
    let mut selection = controls();
    selection.law = GasAbsorptionModel::Iso9613Air {
        relative_humidity: 0.5,
        air,
    };
    let patch = [SurfaceSample {
        position_m: [0.0; 3],
        area_m2: 0.01,
        mode_shape: vec![1.0],
    }];
    let (receiver, report) = BaffledPressure::new_multirate_gas_differentiating(
        &patch,
        1,
        48_000,
        48_000,
        [0.0, 0.0, 0.05],
        &gas,
        &selection,
        &mut || false,
    )
    .unwrap();
    assert!(report.additional_latency_samples > 0);
    let (gain, moment) = receiver
        .kernels
        .iter()
        .fold((0.0, 0.0), |(g, m), (n, row)| {
            (g + row[0], m + *n as f64 * row[0])
        });
    let expected = 0.05 / gas.sound_speed * 48_000.0
        + report.additional_latency_samples as f64
        + path::HALF as f64;
    assert!((moment / gain - expected).abs() < 1e-4);
    assert!(receiver.kernels.last().unwrap().0 < MAX_HISTORY);
}

#[test]
fn separately_prepared_apertures_keep_one_latency_and_their_physical_relative_phase() {
    let (gas, air) =
        fs_material::iso9613::TerrestrialAir::moist_air(293.15, 101_325.0, 0.5).unwrap();
    let mut selection = controls();
    selection.law = GasAbsorptionModel::Iso9613Air { relative_humidity: 0.5, air };
    let surfaces = [
        SurfaceSample { position_m: [0.0; 3], area_m2: 0.01, mode_shape: vec![1.0] },
        SurfaceSample { position_m: [2.0, 0.0, 0.0], area_m2: 0.02, mode_shape: vec![-0.4] },
    ];
    let mut receivers: Vec<_> = surfaces.iter().map(|surface|
        BaffledPressure::new_multirate_gas_differentiating(
            core::slice::from_ref(surface), 1, 48_000, 48_000, [0.0, 0.0, 0.05],
            &gas, &selection, &mut || false,
        ).unwrap()
    ).collect();
    assert_ne!(receivers[0].1.additional_latency_samples, receivers[1].1.additional_latency_samples);
    let latency = receivers.iter().map(|(_, report)| report.additional_latency_samples).max().unwrap();
    assert!(latency > 0);
    for (receiver, report) in &mut receivers {
        receiver.align_gas_latency(report, latency).unwrap();
        let (gain, moment) = receiver.kernels.iter().fold((0.0, 0.0), |(g,m), (n,row)|
            (g + row[0], m + *n as f64 * row[0]));
        let expected = report.filters[0].range_m / gas.sound_speed * 48_000.0
            + latency as f64 + path::HALF as f64;
        assert!((moment/gain - expected).abs() < 1e-4);
        assert_eq!(report.additional_latency_samples, latency);
    }
    assert_eq!(receivers[0].0.filter_delay_samples(), receivers[1].0.filter_delay_samples());
    // Independent physical contributions must reconstruct the single signed
    // surface integral, including frequency-dependent gas phase and flight time.
    let (together, report) = BaffledPressure::new_multirate_gas_differentiating(
        &surfaces, 1, 48_000, 48_000, [0.0, 0.0, 0.05], &gas, &selection, &mut || false,
    ).unwrap();
    assert_eq!(report.additional_latency_samples, latency);
    for frequency in [100.0, 4000.0, 12_000.0, 21_000.0] {
        let omega = TAU * frequency / 48_000.0;
        let a = kernel(&receivers[0].0, 0, omega);
        let b = kernel(&receivers[1].0, 0, omega);
        let actual = kernel(&together, 0, omega);
        assert!((actual[0]-a[0]-b[0]).hypot(actual[1]-a[1]-b[1]) < 1e-12);
    }
}

#[test]
fn exact_distance_sharing_and_all_work_source_and_cancellation_refusals_are_explicit() {
    let gas = gas(1.0);
    let surface = surface();
    let make = |selection: &GasPropagation| {
        BaffledPressure::new_multirate_gas_differentiating(
            &surface,
            2,
            48_000,
            48_000,
            [0.2, 0.1, 1.0],
            &gas,
            selection,
            &mut || false,
        )
    };
    let mut selected = controls();
    selected.limits.max_unique_paths = 1;
    assert!(make(&selected).is_err());
    selected = controls();
    selected.limits.max_filter_taps = 1;
    assert!(make(&selected).is_err());
    selected = controls();
    selected.limits.max_kernel_terms = 1;
    assert!(make(&selected).is_err());
    selected = controls();
    selected.law = GasAbsorptionModel::MolecularRelaxation;
    assert!(make(&selected).is_err());
    let (_, air) = fs_material::iso9613::TerrestrialAir::moist_air(293.15, 101_325.0, 0.5).unwrap();
    selected.law = GasAbsorptionModel::Iso9613Air {
        relative_humidity: 0.5,
        air,
    };
    assert!(
        make(&selected).is_err(),
        "an air admission cannot be applied to a different gas"
    );
    assert!(
        BaffledPressure::new_multirate_gas_differentiating(
            &surface,
            2,
            48_000,
            48_000,
            [0.2, 0.1, 1.0],
            &gas,
            &controls(),
            &mut || true
        )
        .is_err()
    );
    let repeated = [surface[0].clone(), surface[0].clone()];
    let (_, report) = BaffledPressure::new_multirate_gas_differentiating(
        &repeated,
        2,
        48_000,
        48_000,
        [0.2, 0.1, 1.0],
        &gas,
        &controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.filters.len(), 1);
    assert_eq!(report.surface_filter_indices, [0, 0]);
}

#[test]
fn opening_array_keeps_each_signed_path_and_material_transport_at_fixed_density_and_speed() {
    use super::super::{CircularOutletReceiver, openings::CircularOpening};
    let openings = [
        CircularOpening { center_m: [-0.05, 0.0], radius_m: 0.001 },
        CircularOpening { center_m: [0.05, 0.0], radius_m: 0.001 },
    ];
    let config = CircularOutletReceiver { position_m: [0.02, 0.0, 0.5],
        radial_rings: 1, angular_points: 8, maximum_frequency_hz: 12_000.0 };
    let build = |gas: &GasState, selection: &GasPropagation| {
        BaffledPressure::circular_openings_multirate_gas_differentiating(
            &openings, 48_000, 48_000, config, gas, selection, &mut || false,
        ).unwrap()
    };
    let gases = [gas(1.0), gas(100.0)];
    let (mut low, report) = build(&gases[0], &controls());
    let (mut high, altered) = build(&gases[1], &controls());
    assert_eq!(report.surface_filter_indices.len(), 16);
    assert!(report.filters.len() > 1, "actual path lengths cannot collapse to one range");
    assert_eq!(report.gas.density, altered.gas.density);
    assert_eq!(report.gas.sound_speed, altered.gas.sound_speed);
    assert_ne!(report.gas.dynamic_viscosity, altered.gas.dynamic_viscosity);
    let mut difference = 0.0;
    for sample in 0..1024 {
        let v = (TAU * 12_000.0 * sample as f64 / 48_000.0).sin() * 1e-6;
        let a = low.step_trace(&[v, 0.0]).unwrap();
        let b = high.step_trace(&[v, 0.0]).unwrap();
        difference += (a-b).powi(2);
    }
    assert!(difference > 1e-14);
    let mut lossless = controls();
    lossless.law = GasAbsorptionModel::Lossless;
    let (mut explicit, _) = build(&gases[0], &lossless);
    let mut original = BaffledPressure::circular_openings_multirate_differentiating(
        &openings, 48_000, 48_000, config, (&gases[0]).into(),
    ).unwrap();
    for sample in 0..512 {
        let v = (sample as f64 * 0.3).sin() * 1e-6;
        let flow = [v, -v];
        let (a,b) = (explicit.step_trace(&flow).unwrap(), original.step_trace(&flow).unwrap());
        assert!((a-b).abs() < 1e-12, "original signed cancellation is retained");
    }
}

#[test]
fn shared_native_gas_record_rejects_implicit_species_laws_and_unbounded_work() {
    let read = |text: &str| GasPropagation::read_record(&mut text.split_ascii_whitespace());
    assert_eq!(read("stokes-kirchhoff 4096 1000000 50000000").unwrap().law,
        GasAbsorptionModel::StokesKirchhoff);
    assert_eq!(read("lossless 1 1 1000").unwrap().law, GasAbsorptionModel::Lossless);
    for text in ["", "iso9613 1 1 1", "molecular 1 1 1",
        "stokes-kirchhoff 0 1 1", "stokes-kirchhoff 16385 1 1",
        "stokes-kirchhoff 1 4000001 1", "stokes-kirchhoff 1 1 500000001",
        "lossless 1 1", "lossless 1 1 1 extra"] {
        assert!(read(text).is_err(), "{text}");
    }
}
