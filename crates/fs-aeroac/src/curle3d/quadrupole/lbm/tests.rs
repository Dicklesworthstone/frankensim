use super::*;
use fs_lbm::d3q19::{BoundarySpec3, Q3, equilibrium3};
use fs_lbm::material::plan_scaling_from_material;
use fs_material::fluid_transport::FluidTransportState;

fn material() -> FluidTransportState {
    FluidTransportState::try_new(1.2, 1000.0, 0.025, 1.8e-5).unwrap()
}

fn config() -> LbmReynoldsVolumeConfig3 {
    LbmReynoldsVolumeConfig3 {
        lower: [0; 3],
        upper: [4; 3],
        cell_zero_center_m: [0.0; 3],
        grid_to_source: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        max_cells: 64,
        max_workspace_bytes: 16_384,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-20),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g0_native_force_centered_tensor_has_original_si_units_and_frame() {
    let material = material();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 2.0).unwrap();
    let uniform_force = [0.001, 0.002, -0.0015];
    let mut grid = BoundaryGrid3::new(
        4,
        4,
        4,
        scaling.flow().tau,
        uniform_force,
        BoundarySpec3::periodic(),
    );
    let density = 1.25;
    let raw_velocity = [0.02, -0.01, 0.03];
    let populations = equilibrium3(density, raw_velocity);
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                grid.set_populations(x, y, z, &populations);
            }
        }
    }
    let local_force = |[x, y, z]: [usize; 3]| {
        [
            0.0001 * (x + 1) as f64,
            -0.0001 * (y + 1) as f64,
            0.0002 * (z + 1) as f64,
        ]
    };
    grid.set_force_field(local_force);
    let mut setup = config();
    setup.lower = [1, 1, 1];
    setup.upper = [3, 2, 3];
    setup.cell_zero_center_m = [0.7, -0.4, 0.2];
    setup.grid_to_source = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let mut volume = PreparedLbmReynoldsVolume3::new(&grid, &scaling, setup).unwrap();
    assert_eq!(volume.source_count(), 4);
    assert!(volume.stress_moments_n_m().is_none());
    assert_eq!(volume.sample_time_s(), None);
    assert_eq!(volume.timestep_s(), 0.001);
    volume.sample(&grid, 1.25).unwrap();
    let moments = volume.stress_moments_n_m().unwrap();
    for (index, [x, y, z]) in [[1, 1, 1], [2, 1, 1], [1, 1, 2], [2, 1, 2]]
        .into_iter()
        .enumerate()
    {
        // Independent equilibrium/Guo arithmetic, not a copy of grid.velocity().
        let local = local_force([x, y, z]);
        let (native_density, native_velocity) = grid.macroscopic(x, y, z).unwrap();
        assert_eq!(native_density, grid.density(x, y, z));
        assert_eq!(native_velocity, grid.velocity(x, y, z));
        let grid_velocity: [f64; 3] = core::array::from_fn(|axis| {
            (raw_velocity[axis] + 0.5 * (uniform_force[axis] + local[axis]) / density)
                * (0.02 / 0.001)
        });
        let [u, v, w] = [-grid_velocity[1], grid_velocity[0], grid_velocity[2]];
        let cell_mass = 1.2 * density * (0.02 * 0.02 * 0.02);
        for (actual, expected) in
            moments[index]
                .into_iter()
                .zip([u * u, v * v, w * w, u * v, u * w, v * w])
        {
            close(actual, cell_mass * expected);
        }
        for (actual, expected) in volume.source_positions_m()[index].into_iter().zip([
            0.7 - y as f64 * 0.02,
            -0.4 + x as f64 * 0.02,
            0.2 + z as f64 * 0.02,
        ]) {
            close(actual, expected);
        }
        assert_eq!(grid.populations(x, y, z), populations);
    }
    assert_eq!(volume.sample_time_s(), Some(1.25));
}

#[test]
fn g1_accepted_native_force_step_updates_the_original_volume_source() {
    let material = material();
    let scaling = plan_scaling_from_material(&material, 0.01, 0.002, 0.5).unwrap();
    let force = [1e-4, -2e-4, 3e-4];
    let mut grid = BoundaryGrid3::new(
        4,
        4,
        4,
        scaling.flow().tau,
        force,
        BoundarySpec3::periodic(),
    );
    let mut volume = PreparedLbmReynoldsVolume3::new(&grid, &scaling, config()).unwrap();
    let initial = volume.sample(&grid, 0.0).unwrap().to_vec();
    let original_mass = grid.total_mass();
    grid.step();
    volume.sample(&grid, scaling.timestep_s()).unwrap();
    let cell_mass = 1.2 * 0.01 * 0.01 * 0.01;
    let [u, v, w] = force.map(|component| 1.5 * component * (0.01 / 0.002));
    let expected = [u * u, v * v, w * w, u * v, u * w, v * w];
    for (after, before) in volume.stress_moments_n_m().unwrap().iter().zip(initial) {
        for ((actual, before), expected) in after.iter().zip(before).zip(expected) {
            close(*actual, cell_mass * expected);
            close(*actual, 9.0 * before);
        }
    }
    close(grid.total_mass(), original_mass);
    assert_eq!(volume.source_positions_m()[0], [0.0; 3]);
    assert_eq!(volume.source_positions_m()[63], [0.03; 3]);
}

#[test]
fn g4_late_topology_data_clock_and_cancellation_refusals_keep_accepted_sample() {
    let material = material();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 2.0).unwrap();
    let mut grid = BoundaryGrid3::new(
        4,
        4,
        4,
        scaling.flow().tau,
        [0.001, -0.002, 0.003],
        BoundarySpec3::periodic(),
    );
    let mut volume = PreparedLbmReynoldsVolume3::new(&grid, &scaling, config()).unwrap();
    let accepted = volume.sample(&grid, 0.0).unwrap().to_vec();
    for time in [f64::NAN, f64::INFINITY, 0.0, -0.001] {
        assert!(volume.sample(&grid, time).is_err());
    }
    let mut changed_topology = BoundaryGrid3::new(
        4,
        4,
        4,
        scaling.flow().tau,
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    changed_topology.voxelize_sdf(|p| if p == [3.5; 3] { -1.0 } else { 1.0 });
    assert!(volume.sample(&changed_topology, 0.001).is_err());
    let changed_size = BoundaryGrid3::new(
        8,
        4,
        4,
        scaling.flow().tau,
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    assert!(volume.sample(&changed_size, 0.001).is_err());
    let last_cell = grid.populations(3, 3, 3);
    grid.set_populations(3, 3, 3, &[0.0; Q3]);
    assert!(grid.macroscopic(3, 3, 3).is_err());
    assert!(volume.sample(&grid, 0.001).is_err());
    grid.set_populations(3, 3, 3, &last_cell);
    let mut polls = 0;
    assert!(matches!(
        volume.sample_cancellable(&grid, 0.001, || {
            polls += 1;
            polls >= 3
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(volume.sample_time_s(), Some(0.0));
    assert_eq!(volume.stress_moments_n_m().unwrap(), accepted);
    assert_eq!(volume.sample(&grid, 0.001).unwrap(), accepted);
    assert_eq!(volume.sample_time_s(), Some(0.001));
    assert_eq!(grid.populations(3, 3, 3), last_cell);
}

#[test]
fn g0_preparation_keeps_fixed_fluid_membership_and_admits_frame_and_budgets() {
    let material = material();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 2.0).unwrap();
    let mut grid = BoundaryGrid3::new(
        4,
        4,
        4,
        scaling.flow().tau,
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    grid.voxelize_sdf(|p| if p == [0.5; 3] { -1.0 } else { 1.0 });
    let mut volume = PreparedLbmReynoldsVolume3::new(&grid, &scaling, config()).unwrap();
    assert_eq!(volume.source_count(), 63);
    assert_eq!(volume.source_positions_m()[0], [0.02, 0.0, 0.0]);
    assert!(
        volume
            .sample(&grid, 0.0)
            .unwrap()
            .iter()
            .all(|q| *q == [0.0; 6])
    );

    let mut bad = config();
    bad.max_cells = 63;
    assert!(PreparedLbmReynoldsVolume3::new(&grid, &scaling, bad).is_err());
    bad = config();
    bad.max_workspace_bytes = volume.workspace_bytes() - 1;
    assert!(PreparedLbmReynoldsVolume3::new(&grid, &scaling, bad).is_err());
    bad = config();
    bad.upper = [1; 3];
    assert!(PreparedLbmReynoldsVolume3::new(&grid, &scaling, bad).is_err());
    bad = config();
    bad.lower[0] = bad.upper[0];
    assert!(PreparedLbmReynoldsVolume3::new(&grid, &scaling, bad).is_err());
    bad = config();
    bad.upper[2] = 5;
    assert!(PreparedLbmReynoldsVolume3::new(&grid, &scaling, bad).is_err());
    bad = config();
    bad.grid_to_source[1][1] = -1.0;
    assert!(PreparedLbmReynoldsVolume3::new(&grid, &scaling, bad).is_err());
    bad = config();
    bad.grid_to_source[0][1] = 0.1;
    assert!(PreparedLbmReynoldsVolume3::new(&grid, &scaling, bad).is_err());
    bad = config();
    bad.cell_zero_center_m[0] = f64::NAN;
    assert!(PreparedLbmReynoldsVolume3::new(&grid, &scaling, bad).is_err());
    assert!(matches!(
        PreparedLbmReynoldsVolume3::new_cancellable(&grid, &scaling, config(), || true),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn g1_native_periodic_flow_to_recorded_quadrupoles_retains_quadratic_frequency_and_coherence() {
    use crate::curle3d::CoherenceBudget;
    use crate::curle3d::quadrupole::QuadrupolePhasor;
    use crate::curle3d::quadrupole::recorded::{QuadrupoleRecordingBudget, QuadrupoleSampleClock};
    use fs_exec::CancelGate;
    use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
    use fs_material::fluid_acoustics::FluidAcousticState;
    use fs_math::{c64::C64, det};

    const SAMPLES: usize = 64;
    const VELOCITY_BIN: usize = 3;
    let material = material();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 2.0).unwrap();
    let clock = QuadrupoleSampleClock::new(0.25, scaling.timestep_s(), 0.0).unwrap();
    let mut grid = BoundaryGrid3::new(
        4,
        4,
        4,
        scaling.flow().tau,
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    let amplitude = [0.01, 0.004, -0.002];
    let theta = core::f64::consts::TAU * VELOCITY_BIN as f64 / SAMPLES as f64;
    // Native endpoint forcing integrates (F_n + F_(n+1))/2. This forcing
    // therefore gives the exact discrete uniform solution u_n=A*sin(n theta).
    let force_gain = 2.0 * det::sin(theta * 0.5) / det::cos(theta * 0.5);
    grid.rebind_force_field(|_| amplitude.map(|a| a * force_gain));
    let mut selection = config();
    selection.lower = [1; 3];
    selection.upper = [3, 2, 2];
    let mut volume = PreparedLbmReynoldsVolume3::new(&grid, &scaling, selection).unwrap();
    let window = WelchConfig {
        segment_samples: SAMPLES,
        hop_samples: SAMPLES,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    };
    let budget = QuadrupoleRecordingBudget {
        maximum_sources: 2,
        welch: WelchBudget {
            maximum_channels: 12,
            maximum_samples_per_channel: SAMPLES,
            maximum_transforms: 12,
            maximum_workspace_bytes: 100_000,
        },
    };
    let mut record = QuadrupoleWelchStream::new(
        volume.source_positions_m(),
        clock,
        SAMPLES,
        window,
        budget,
        || false,
    )
    .unwrap();
    assert!(volume.record_sample(&mut record, || false).is_err());
    let mut wrong_clock = QuadrupoleWelchStream::new(
        volume.source_positions_m(),
        QuadrupoleSampleClock::new(clock.epoch_s(), 2.0 * scaling.timestep_s(), 0.0).unwrap(),
        SAMPLES,
        window,
        budget,
        || false,
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    for n in 0..SAMPLES {
        volume.sample(&grid, clock.time_at(n).unwrap()).unwrap();
        if n == 0 {
            assert!(volume.record_sample(&mut wrong_clock, || false).is_err());
            assert_eq!(wrong_clock.accepted_sample_count(), 0);
            assert!(volume.record_sample(&mut record, || true).is_err());
            assert_eq!(record.accepted_sample_count(), 0);
            assert_eq!(volume.sample_time_s(), Some(clock.epoch_s()));
        }
        volume.record_sample(&mut record, || false).unwrap();
        if n + 1 < SAMPLES {
            grid.step_with_endpoint_force_scale(det::cos(theta * (n + 1) as f64), &gate)
                .unwrap();
        }
    }
    let spectrum = record.finish().unwrap();
    let acoustic = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let stress_bin = 2 * VELOCITY_BIN;
    let observer = spectrum
        .prepare_observer(acoustic, [0.4, -0.3, 0.2], stress_bin, || false)
        .unwrap();
    let work = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 2,
    };
    let observed = spectrum
        .evaluate_bin(&observer, stress_bin, work, || false)
        .unwrap();

    // sin²(theta)=1/2-cos(2 theta)/2: transform the actual tensor, not velocity.
    let [u, v, w] = amplitude.map(|a| a * (0.02 / 0.001));
    let scale = -0.5 * 1.2 * (0.02 * 0.02 * 0.02);
    let physical_source = QuadrupolePhasor {
        stress_moment_n_m: [u * u, v * v, w * w, u * v, u * w, v * w]
            .map(|entry| C64::from_re(scale * entry)),
    };
    let coherent = observer.evaluate(&[physical_source; 2]).unwrap();
    let expected_power = 0.5 * coherent.pressure_pa.norm_sq();
    assert!(expected_power > 0.0);
    assert!((observed.mean_square_pressure_pa2 / expected_power - 1.0).abs() < 2e-9);
    for i in 0..3 {
        let expected_velocity_power = 0.5 * coherent.particle_velocity_m_s[i].norm_sq();
        assert!(
            (observed.velocity_covariance_m2_s2[i][i].re / expected_velocity_power - 1.0).abs()
                < 2e-9
        );
    }
    let first_harmonic = spectrum
        .prepare_observer(acoustic, [0.4, -0.3, 0.2], VELOCITY_BIN, || false)
        .unwrap();
    let absent = spectrum
        .evaluate_bin(&first_harmonic, VELOCITY_BIN, work, || false)
        .unwrap();
    assert!(absent.mean_square_pressure_pa2 < 1e-16 * observed.mean_square_pressure_pa2);
}
