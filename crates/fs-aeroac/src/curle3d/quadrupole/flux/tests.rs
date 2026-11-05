use super::*;
use fs_flux::compressible::spatial::mixture::state_into;
use fs_flux::compressible::spatial::{AxisBoundaries, State, StepControl};
use fs_gas::dynamics::CaloricDomain;
use fs_qty::{Temperature, Time};
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, cantera_h2o2};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};

fn frozen(species: usize) -> FrozenGas {
    // The existing pinned NASA source supplies all calorics; no substitute EOS.
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let model = mechanism.thermodynamics().models()[species].clone();
    let source = FrozenIdealGasMixtureModelV1::new(
        vec![model],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    FrozenGas::new(
        source,
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap()
}

fn config(gas: &FrozenGas) -> CompressibleVolumeConfig3 {
    let p0 = gas.gas_constant() * 500.0;
    let ambient = State::from_primitive(gas, 1.0, [0.0; 3], p0)
        .unwrap()
        .primitive(gas, || false)
        .unwrap();
    CompressibleVolumeConfig3 {
        acoustic: FluidAcousticState::try_new(1.0, ambient.sound_speed).unwrap(),
        reference_pressure_pa: ambient.pressure,
        origin_m: [0.2, 0.3, -0.1],
        grid_to_acoustic: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
    }
}

fn frozen_flow(gas: FrozenGas, grid: CartesianGrid) -> Euler3d {
    let state = State::from_primitive(
        &gas,
        1.2,
        [2.0, -3.0, 4.0],
        1.2 * gas.gas_constant() * 620.0,
    )
    .unwrap();
    let count: usize = grid.shape().into_iter().product();
    Euler3d::new(gas, grid, vec![state; count], PERIODIC, || false).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_frozen_native_units_rotation_and_accepted_clock() {
    let gas = frozen(6); // original N2 species
    let settings = config(&gas);
    let grid = CartesianGrid::new(2, 2, 2, 0.02, 0.03, 0.04).unwrap();
    let mut flow = frozen_flow(gas.clone(), grid);
    let selected = [0, 2, 7];
    let states_before = flow.states().to_vec();
    let bytes = PreparedCompressibleVolume3::required_workspace_bytes(selected.len()).unwrap();
    let mut sample =
        PreparedCompressibleVolume3::from_frozen(&flow, &selected, settings, bytes, || false)
            .unwrap();
    assert_eq!(flow.states(), states_before);
    assert_eq!(sample.cell_indices(), &selected);
    assert_eq!(sample.frozen_source(), Some(flow.gas()));
    assert_eq!(sample.sample_time_s(), flow.time().value());
    assert!(sample.mixture_source().is_none());
    let rho = 1.2;
    let pressure = rho * gas.gas_constant() * 620.0;
    let c = settings.acoustic.sound_speed_m_s();
    let isotropic = (pressure - settings.reference_pressure_pa) - c * c * (rho - 1.0);
    let volume = 0.02 * 0.03 * 0.04;
    let u = [3.0, 2.0, 4.0];
    let expected = [
        (rho * u[0] * u[0] + isotropic) * volume,
        (rho * u[1] * u[1] + isotropic) * volume,
        (rho * u[2] * u[2] + isotropic) * volume,
        rho * u[0] * u[1] * volume,
        rho * u[0] * u[2] * volume,
        rho * u[1] * u[2] * volume,
    ];
    for (i, &cell) in selected.iter().enumerate() {
        let [x, y, z] = grid.cell_center(cell).unwrap();
        assert_eq!(sample.source_positions_m()[i], [0.2 - y, 0.3 + x, -0.1 + z]);
        for (observed, wanted) in sample.stress_moments_n_m()[i].into_iter().zip(expected) {
            close(observed, wanted);
        }
    }
    flow.step(Time::new(1e-6), StepControl::default(), || false)
        .unwrap();
    sample.sample_frozen(&flow, || false).unwrap();
    assert!(sample.sample_time_s() > 0.0);
    assert_eq!(sample.sample_time_s(), flow.time().value());
    for row in sample.stress_moments_n_m() {
        for (&actual, wanted) in row.iter().zip(expected) {
            close(actual, wanted);
        }
    }
}

#[test]
fn g1_evolving_species_use_original_calorics_and_native_pressure() {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap();
    let grid = CartesianGrid::new(3, 1, 1, 0.04, 0.02, 0.03).unwrap();
    let n = gas.species_count();
    let mut states = vec![0.0; 3 * (n + 4)];
    for (i, row) in states.chunks_exact_mut(n + 4).enumerate() {
        let mut density = vec![0.0; n];
        density[0] = i as f64 * 0.35; // AR
        density[6] = 0.7 - density[0]; // N2
        state_into(
            &gas,
            &density,
            [1.0 + 0.1 * i as f64, -0.5, 0.25],
            100000.0,
            row,
        )
        .unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let mut flow = MixtureEuler3d::new(gas, grid, states, PERIODIC, bytes, || false).unwrap();
    let settings = config(&frozen(6));
    let mut sample = PreparedCompressibleVolume3::from_mixture(
        &flow,
        &[0, 1, 2],
        settings,
        PreparedCompressibleVolume3::required_workspace_bytes(3).unwrap(),
        || false,
    )
    .unwrap();
    assert_eq!(sample.mixture_source(), Some(flow.gas()));
    assert!(sample.frozen_source().is_none());
    for i in 0..3 {
        let velocity = [0.5, 1.0 + 0.1 * i as f64, 0.25];
        let isotropic = (100000.0 - settings.reference_pressure_pa)
            - settings.acoustic.sound_speed_m_s().powi(2) * (0.7 - 1.0);
        for a in 0..3 {
            close(
                sample.stress_moments_n_m()[i][a],
                (0.7 * velocity[a] * velocity[a] + isotropic) * grid.cell_volume(),
            );
        }
        close(
            sample.stress_moments_n_m()[i][3],
            0.7 * velocity[0] * velocity[1] * grid.cell_volume(),
        );
    }
    // Changing composition remains in the original solver; extraction follows
    // its accepted pressure and density, not a frozen gamma or a new clock.
    flow.step(Time::new(1e-7), StepControl::default(), || false)
        .unwrap();
    let before = flow.states().to_vec();
    sample.sample_mixture(&flow, || false).unwrap();
    assert_eq!(flow.states(), before);
    assert_eq!(sample.sample_time_s(), flow.time().value());
    for i in 0..3 {
        let p = flow
            .cell(i)
            .unwrap()
            .primitive(flow.gas(), || false)
            .unwrap();
        let psi = (p.pressure - settings.reference_pressure_pa)
            - settings.acoustic.sound_speed_m_s().powi(2) * (p.density - 1.0);
        let u = [-p.velocity[1], p.velocity[0], p.velocity[2]];
        for a in 0..3 {
            close(
                sample.stress_moments_n_m()[i][a],
                (p.density * u[a] * u[a] + psi) * grid.cell_volume(),
            );
        }
    }
}

#[test]
fn g4_source_layout_clock_budget_and_cancellation_preserve_accepted_sample() {
    let gas = frozen(6);
    let settings = config(&gas);
    let grid = CartesianGrid::new(2, 2, 2, 0.02, 0.03, 0.04).unwrap();
    let mut flow = frozen_flow(gas.clone(), grid);
    let cells = [0, 2, 7];
    let bytes = PreparedCompressibleVolume3::required_workspace_bytes(cells.len()).unwrap();
    assert!(matches!(
        PreparedCompressibleVolume3::from_frozen(&flow, &cells, settings, bytes - 1, || false),
        Err(CompressibleVolumeError3::Workspace { .. })
    ));
    for invalid in [&[][..], &[2, 0][..], &[1, 1][..], &[8][..]] {
        assert!(
            PreparedCompressibleVolume3::from_frozen(&flow, invalid, settings, bytes, || false)
                .is_err()
        );
    }
    let mut bad = settings;
    bad.grid_to_acoustic[2][2] = -1.0;
    assert!(PreparedCompressibleVolume3::from_frozen(&flow, &cells, bad, bytes, || false).is_err());
    let mut sample =
        PreparedCompressibleVolume3::from_frozen(&flow, &cells, settings, bytes, || false).unwrap();
    let first = sample.stress_moments_n_m().to_vec();
    assert!(matches!(
        sample.sample_frozen(&flow, || false),
        Err(CompressibleVolumeError3::SampleClock { .. })
    ));
    let other_source = frozen_flow(frozen(0), grid);
    assert_eq!(
        sample.sample_frozen(&other_source, || false),
        Err(CompressibleVolumeError3::SourceMismatch)
    );
    let other_grid = frozen_flow(gas, CartesianGrid::new(2, 2, 2, 0.021, 0.03, 0.04).unwrap());
    assert_eq!(
        sample.sample_frozen(&other_grid, || false),
        Err(CompressibleVolumeError3::GeometryMismatch)
    );
    flow.step(Time::new(1e-6), StepControl::default(), || false)
        .unwrap();
    let mut calls = 0;
    assert_eq!(
        sample.sample_frozen(&flow, || {
            calls += 1;
            calls >= 3
        }),
        Err(CompressibleVolumeError3::Cancelled)
    );
    assert_eq!(sample.sample_time_s(), 0.0);
    assert_eq!(sample.stress_moments_n_m(), first);
    sample.sample_frozen(&flow, || false).unwrap();
    assert_eq!(sample.sample_time_s(), flow.time().value());
}

#[test]
fn g1_native_flow_history_reaches_coherent_acoustic_bin_without_clock_or_phase_replacement() {
    use super::super::recorded::{
        QuadrupoleRecordingBudget, QuadrupoleSampleClock, QuadrupoleWelchStream,
    };
    use super::super::{CoherenceBudget, QuadrupolePhasor};
    use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
    use fs_math::c64::C64;

    let gas = frozen(6);
    let settings = config(&gas);
    let grid = CartesianGrid::new(2, 2, 2, 0.02, 0.03, 0.04).unwrap();
    let states = (0..8)
        .map(|cell| {
            let sign = if cell % 2 == 0 { -1.0 } else { 1.0 };
            State::from_primitive(
                &gas,
                1.0,
                [sign, 0.3, -0.2],
                gas.gas_constant() * 600.0 * (1.0 + 0.01 * sign),
            )
            .unwrap()
        })
        .collect();
    let mut flow = Euler3d::new(gas, grid, states, PERIODIC, || false).unwrap();
    let cells = [0, 2, 7];
    let mut source = PreparedCompressibleVolume3::from_frozen(
        &flow,
        &cells,
        settings,
        PreparedCompressibleVolume3::required_workspace_bytes(cells.len()).unwrap(),
        || false,
    )
    .unwrap();
    let count = 8;
    let dt = 1e-6;
    let mut recorder = QuadrupoleWelchStream::new(
        source.source_positions_m(),
        QuadrupoleSampleClock::new(source.sample_time_s(), dt, 1e-15).unwrap(),
        count,
        WelchConfig {
            segment_samples: count,
            hop_samples: count,
            window: WelchWindow::Rectangular,
            detrend: WelchDetrend::None,
        },
        QuadrupoleRecordingBudget {
            maximum_sources: cells.len(),
            welch: WelchBudget {
                maximum_channels: 6 * cells.len(),
                maximum_samples_per_channel: count,
                maximum_transforms: 6 * cells.len(),
                maximum_workspace_bytes: 1 << 20,
            },
        },
        || false,
    )
    .unwrap();
    let mut history = Vec::new();
    for n in 0..count {
        if n != 0 {
            flow.step(Time::new(dt), StepControl::default(), || false)
                .unwrap();
            source.sample_frozen(&flow, || false).unwrap();
        }
        source.record_sample(&mut recorder, || false).unwrap();
        history.push(source.stress_moments_n_m().to_vec());
    }
    let spectrum = recorder.finish().unwrap();
    let observer = spectrum
        .prepare_observer(settings.acoustic, [1.1, -0.7, 0.9], 1, || false)
        .unwrap();
    let result = spectrum
        .evaluate_bin(
            &observer,
            1,
            CoherenceBudget {
                max_modes: 1,
                max_source_mode_pairs: cells.len(),
            },
            || false,
        )
        .unwrap();
    // Independent direct DFT of the actual native stress samples. With one
    // rectangular segment, 2/N times its conjugate DFT is the peak factor.
    // This exercises the public owner chain, not a spatial/noise qualification.
    let mut peaks = vec![QuadrupolePhasor::default(); cells.len()];
    for (n, sample) in history.iter().enumerate() {
        let theta = core::f64::consts::TAU * n as f64 / count as f64;
        let phase = C64::new(theta.cos(), theta.sin()).scale(2.0 / count as f64);
        for (peak, row) in peaks.iter_mut().zip(sample) {
            for (component, &moment) in peak.stress_moment_n_m.iter_mut().zip(row) {
                *component = *component + phase.scale(moment);
            }
        }
    }
    let direct = observer.evaluate(&peaks).unwrap();
    let expected_pressure = 0.5 * (direct.pressure_pa.re.powi(2) + direct.pressure_pa.im.powi(2));
    assert!(expected_pressure > 0.0);
    assert!(
        (result.mean_square_pressure_pa2 - expected_pressure).abs() <= 2e-8 * expected_pressure
    );
    for a in 0..3 {
        let u = direct.particle_velocity_m_s[a];
        let expected = 0.5 * (u.re * u.re + u.im * u.im);
        assert!((result.velocity_covariance_m2_s2[a][a].re - expected).abs() <= 2e-8 * expected);
        assert!(
            (result.active_intensity_w_m2[a] - direct.active_intensity_w_m2[a]).abs()
                <= 2e-8 * (direct.active_intensity_w_m2[a].abs() + 1e-30)
        );
    }
}
