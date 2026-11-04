use super::*;
use crate::curle3d::array::CompactSourcePhasor;
use crate::slot_jet_3d::distributed::SlotJet3dPatch;
use fs_fft::welch::{WelchDetrend, WelchWindow};

fn diagnostics(samples: usize) -> SlotJet3dDiagnostics {
    SlotJet3dDiagnostics {
        mach_max_lattice: 0.0,
        flux_plate_plane: 0.0,
        flux_fringe_plane: 0.0,
        reynolds: 0.0,
        record_len: samples,
        strouhal_bin_width: 0.0,
    }
}

fn spectrum_budget() -> WelchBudget {
    WelchBudget {
        maximum_channels: 8,
        maximum_samples_per_channel: 128,
        maximum_transforms: 16,
        maximum_workspace_bytes: 1_000_000,
    }
}

#[test]
fn g1_opposite_spatial_forces_radiate_with_measured_cross_phase_and_fz() {
    let gas = GasState::try_new_moist_air(300.0, 97_000.0, 0.2).unwrap();
    let forces = [
        C64::new(0.02, 0.0),
        C64::new(0.0, -0.01),
        C64::new(0.03, 0.01),
    ];
    let dt = 1.0 / 64_000.0;
    let mut history = Vec::new();
    for j in 0..128 {
        let (s, c) = (core::f64::consts::TAU * j as f64 / 64.0).sin_cos();
        let f = forces.map(|v| v.re * c + v.im * s);
        history.push(f);
        history.push(f.map(|v| -v));
    }
    let record = MaterialSlotJet3dPatchRun {
        material: &gas,
        acoustic: gas.acoustic_state().unwrap(),
        timestep_s: dt,
        positions_m: vec![[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]],
        radii_m: vec![0.001; 2],
        force_on_fluid_n: history,
        final_plane_physical_mach_max: 0.0,
        lattice_diagnostics: diagnostics(128),
        sgs: None,
    };
    let spectrum = record
        .estimate_welch(
            WelchConfig {
                segment_samples: 64,
                hop_samples: 64,
                window: WelchWindow::Rectangular,
                detrend: WelchDetrend::None,
            },
            spectrum_budget(),
            || false,
        )
        .unwrap();
    assert_eq!(spectrum.cross_spectrum().channel_count(), 8);
    let at = [2.0, 0.7, 0.4];
    let limits = CompactPatchLimits {
        maximum_ka: 0.1,
        maximum_radius_over_range: 0.01,
    };
    let work = CoherenceBudget {
        max_modes: 2,
        max_source_mode_pairs: 4,
    };
    let actual = spectrum.observe_bin(1, at, limits, work, || false).unwrap();
    let direct = PreparedCompactSourceObserver::new(
        record.acoustic(),
        1000.0,
        at,
        record.source_positions_m(),
        2,
    )
    .unwrap();
    let expected = direct
        .evaluate(&[
            CompactSourcePhasor::force_dipole(forces),
            CompactSourcePhasor::force_dipole(forces.map(|f| f.scale(-1.0))),
        ])
        .unwrap();
    let expected_p2 = 0.5 * (expected.pressure_pa.re.powi(2) + expected.pressure_pa.im.powi(2));
    assert!(expected_p2 > 1e-8); // Summing the two force histories first gives zero.
    assert!((actual.mean_square_pressure_pa2 - expected_p2).abs() < 1e-12 * expected_p2);
    for a in 0..3 {
        let expected_flux = expected.active_intensity_w_m2[a];
        assert!(
            (actual.active_intensity_w_m2[a] - expected_flux).abs()
                < 1e-11 * expected_flux.abs().max(1e-20)
        );
    }
    assert!(actual.scope.contains("periodic 3D flow cell"));
    assert!(
        spectrum
            .observe_bin(
                1,
                at,
                CompactPatchLimits {
                    maximum_ka: 1e-6,
                    ..limits
                },
                work,
                || false
            )
            .is_err()
    );
    assert!(
        spectrum
            .observe_bin(1, [0.0; 3], limits, work, || false)
            .is_err()
    );
    assert!(spectrum.observe_bin(0, at, limits, work, || false).is_err());
    assert!(matches!(
        spectrum.observe_bin(1, at, limits, work, || true),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn g1_physical_patch_conversion_preserves_native_force_sign_and_original_material() {
    let gas = GasState::try_new_moist_air(300.0, 100_000.0, 0.0).unwrap();
    let setup = SlotJet3dSetup {
        nx: 64,
        ny: 32,
        nz: 8,
        slot_half: 2.5,
        higher_order_rate: 1.6,
        nozzle_thickness: 1,
        edge_distance: 8,
        plate_length: 2,
        fringe_width: 8,
        fringe_sigma: 0.4,
        seed_amplitude: 0.02,
        steps_settle: 1,
        steps_record: 64,
    };
    let physical = JetPhysicalState {
        cell_size_m: 1e-4,
        timestep_s: 1e-5,
        jet_speed_m_s: 0.4,
        max_speed_m_s: 1.2,
        max_cells: 64 * 32 * 8,
    };
    let source = MaterialSlotJet3d::try_new(&gas, &setup, physical).unwrap();
    let native = [-0.04, 0.012, 0.003];
    let run = source
        .physical_patch_run(SlotJet3dPatchRun {
            patches: vec![SlotJet3dPatch {
                center_lu: [10.0, 16.0, 4.0],
                radius_lu: 2.0,
                measured_links: 20,
            }],
            force_on_fluid_lu: vec![native; 64],
            diagnostics: diagnostics(64),
            sgs: None,
        })
        .unwrap();
    let unit = gas.density * physical.cell_size_m.powi(4) / physical.timestep_s.powi(2);
    for a in 0..3 {
        assert!(
            (run.force_on_fluid_n()[0][a] - native[a] * unit).abs()
                < 1e-14 * (native[a] * unit).abs()
        );
    }
    assert_eq!(run.source_positions_m()[0], [0.001, 0.0016, 0.0004]);
    assert_eq!(run.enclosing_radii_m(), [0.0002]);
    assert_eq!(run.acoustic(), gas.acoustic_state().unwrap());
    assert_eq!(run.timestep_s(), physical.timestep_s);
    assert!(core::ptr::eq(run.material(), &gas));
}
