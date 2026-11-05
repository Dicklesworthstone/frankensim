//! G0/G1/G3: native failed-plane release retains physical stress and history.

use fs_material::brittle::{
    BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings, BrittleError,
};
use fs_material::elastic::IsotropicElastic;
use fs_material::tensor::{contract, rotate, rotation};
use fs_material::{SmallStrainLaw, Voigt};

fn law(normal: [f64; 3]) -> BrittleCrackBand {
    BrittleCrackBand::new(
        IsotropicElastic::new(3e9, 0.25, 0.05).unwrap(),
        BrittleCrackBandParameters {
            tensile_strength_pa: 3e6,
            shear_strength_pa: 6e6,
            mode_i_fracture_energy_j_m2: 120.0,
            mode_ii_fracture_energy_j_m2: 360.0,
            bk_exponent: 1.6,
        },
        BrittleCrackBandSettings {
            band_width_m: 0.004,
            normal,
            maximum_strain: 0.05,
            energy_tolerance_j_m3: 1e-6,
        },
    )
    .unwrap()
}

fn near(a: f64, b: f64, floor: f64) {
    assert!(
        (a - b).abs() <= 2e-10 * a.abs().max(b.abs()).max(floor),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_open_failed_planes_retain_in_plane_storage_and_pure_mode_fracture_work() {
    let band = law([1.0, 0.0, 0.0]);
    for (rupture, fracture) in [
        ([0.022, 0.0, 0.0, 0.0, 0.0, 0.0], 120.0),
        ([0.0, 0.0, 0.0, 0.017, 0.0, 0.0], 360.0),
    ] {
        let failed = band.trial(&rupture, &band.initial_state()).unwrap().state;
        assert!(failed.failed());
        let strain = [0.025, 0.003, -0.001, 0.002, 0.001, -0.003];
        let held = band.trial(&strain, &failed).unwrap();
        let saved = held.state;
        let d = band.elastic().youngs / (1.0 - band.elastic().poisson.powi(2));
        let expected_stress = [
            0.0,
            d * (strain[1] + 0.25 * strain[2]),
            d * (strain[2] + 0.25 * strain[1]),
            0.0,
            band.elastic().youngs / 1.25 * strain[4],
            0.0,
        ];
        for normal in [[2.0, 0.0, 0.0], [-3.0, 0.0, 0.0]] {
            let release = band
                .release_failed_plane(&held.state, normal, 0.0, 0.0)
                .unwrap();
            assert_eq!(release.traction_pa, [0.0; 3]);
            for (actual, expected) in release.stress_pa.into_iter().zip(expected_stress) {
                near(actual, expected, 1.0);
            }
            near(
                release.parallel_stored_energy_j_m3,
                0.5 * contract(&strain, &expected_stress),
                1.0,
            );
            assert_eq!(
                release.parallel_stored_energy_j_m3,
                held.state.stored_energy_j_m3()
            );
            near(
                release.fracture_dissipation_j_m3 * band.band_width_m(),
                fracture,
                1.0,
            );
            assert_eq!(held.state, saved);
        }
        // Release also admits fully unloaded failed bonds without healing or
        // requiring a fictitious geometric opening to prove their history.
        let unloaded = band.trial(&[0.0; 6], &held.state).unwrap().state;
        let release = band
            .release_failed_plane(&unloaded, [1.0, 0.0, 0.0], 0.0, 0.0)
            .unwrap();
        assert_eq!(release.parallel_stored_energy_j_m3, 0.0);
        assert_eq!(
            release.fracture_dissipation_j_m3,
            held.state.dissipated_j_m3()
        );
    }
}

#[test]
fn g3_rotated_tensile_and_shear_rupture_preserve_traction_energy_and_normal_orientation() {
    let q = rotation(
        [
            1.0 / 6.0_f64.sqrt(),
            2.0 / 6.0_f64.sqrt(),
            -1.0 / 6.0_f64.sqrt(),
        ],
        0.73,
    );
    let normal = [q[0][0], q[1][0], q[2][0]];
    let a = law([1.0, 0.0, 0.0]);
    let b = law(normal);
    for rupture in [
        [0.022, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.017, 0.0, 0.0],
    ] {
        let sa = a.trial(&rupture, &a.initial_state()).unwrap().state;
        let sb = b
            .trial(&rotate(&rupture, &q), &b.initial_state())
            .unwrap()
            .state;
        let strain = [0.025, 0.003, -0.001, 0.002, 0.001, -0.003];
        let sa = a.trial(&strain, &sa).unwrap().state;
        let sb = b.trial(&rotate(&strain, &q), &sb).unwrap().state;
        let ra = a
            .release_failed_plane(&sa, [1.0, 0.0, 0.0], 0.0, 0.0)
            .unwrap();
        for scale in [1e-200, -1e200] {
            let rb = b
                .release_failed_plane(&sb, normal.map(|x| x * scale), 2e-14, 1e-6)
                .unwrap();
            for (actual, expected) in rb.stress_pa.into_iter().zip(rotate(&ra.stress_pa, &q)) {
                near(actual, expected, 1.0);
            }
            near(
                rb.parallel_stored_energy_j_m3,
                ra.parallel_stored_energy_j_m3,
                1.0,
            );
            near(
                rb.fracture_dissipation_j_m3,
                ra.fracture_dissipation_j_m3,
                1.0,
            );
            let s = rb.stress_pa;
            let [nx, ny, nz] = rb.face_normal;
            let actual_traction = [
                s[0] * nx + s[3] * ny + s[5] * nz,
                s[3] * nx + s[1] * ny + s[4] * nz,
                s[5] * nx + s[4] * ny + s[2] * nz,
            ];
            for (actual, expected) in rb.traction_pa.into_iter().zip(actual_traction) {
                assert!((actual - expected).abs() < 1e-7);
                assert!(actual.abs() < 1e-6);
            }
            assert!((rb.face_normal[0] * normal[0]).signum() == scale.signum());
        }
    }
}

#[test]
fn g0_surviving_bonds_and_any_active_compression_refuse_even_with_loose_traction_gate() {
    let band = law([1.0, 0.0, 0.0]);
    let initial = band.initial_state();
    let damaged = band
        .trial(&[0.012, 0.0, 0.0, 0.0, 0.0, 0.0], &initial)
        .unwrap()
        .state;
    assert!(damaged.damage() > 0.0 && !damaged.failed());
    for state in [initial, damaged] {
        assert!(
            band.release_failed_plane(&state, band.normal(), 0.0, 1e300)
                .is_err()
        );
    }
    let failed = band
        .trial(&[0.022, 0.0, 0.0, 0.0, 0.0, 0.0], &damaged)
        .unwrap()
        .state;
    for compression in [-1e-12, -0.01] {
        let compressed = band
            .trial(&[compression, 0.0, 0.0, 0.0, 0.0, 0.0], &failed)
            .unwrap()
            .state;
        let saved = compressed;
        assert!(compressed.failed());
        assert!(
            band.release_failed_plane(&compressed, band.normal(), 0.0, 1e300)
                .is_err()
        );
        assert_eq!(compressed, saved);
    }
    assert!(
        band.release_failed_plane(&failed, band.normal(), 0.0, 0.0)
            .is_ok()
    );
}

#[test]
fn g0_original_owner_domain_face_direction_and_actual_traction_are_independent_gates() {
    let band = law([1.0, 0.0, 0.0]);
    let state = band
        .trial(
            &[0.022, 0.003, -0.001, 0.0, 0.001, 0.0],
            &band.initial_state(),
        )
        .unwrap()
        .state;
    let saved = state;
    let clone = band.clone();
    assert!(matches!(
        clone.release_failed_plane(&state, band.normal(), 0.0, 0.0),
        Err(BrittleError::StateMismatch)
    ));
    for normal in [
        [0.0; 3],
        [f64::NAN, 0.0, 0.0],
        [f64::INFINITY, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ] {
        assert!(
            band.release_failed_plane(&state, normal, 0.0, 1e300)
                .is_err()
        );
    }
    for tolerance in [-1.0, 1.0, f64::INFINITY, f64::NAN] {
        assert!(
            band.release_failed_plane(&state, band.normal(), tolerance, 0.0)
                .is_err()
        );
    }
    for tolerance in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            band.release_failed_plane(&state, band.normal(), 0.0, tolerance)
                .is_err()
        );
    }
    let tilted = [1.0, 1e-4, 0.0];
    assert!(
        band.release_failed_plane(&state, tilted, 0.0, 1e300)
            .is_err()
    );
    // Permitting a small geometric angle cannot silently discard the actual
    // traction carried by the remaining in-plane stress through that face.
    assert!(
        band.release_failed_plane(&state, tilted, 1e-3, 1e-6)
            .is_err()
    );
    let invalid = band.small_strain_law().update_state(&[0.1; 6], &state);
    assert!(!invalid.is_valid());
    assert!(
        band.release_failed_plane(&invalid, band.normal(), 0.0, 0.0)
            .is_err()
    );
    assert_eq!(state, saved);
    assert_eq!(
        band.release_failed_plane(&state, band.normal(), 0.0, 0.0)
            .unwrap(),
        band.release_failed_plane(&saved, band.normal(), 0.0, 0.0)
            .unwrap()
    );
}

#[test]
fn g3_thermal_failure_and_cooling_keep_the_accumulated_loss_instead_of_current_toughness() {
    let cold = law([1.0, 0.0, 0.0]);
    let p = cold.parameters();
    let hot = BrittleCrackBand::new(
        IsotropicElastic::new(1.5e9, 0.25, 0.05).unwrap(),
        BrittleCrackBandParameters {
            tensile_strength_pa: 0.1 * p.tensile_strength_pa,
            shear_strength_pa: 0.1 * p.shear_strength_pa,
            mode_i_fracture_energy_j_m2: 0.03 * p.mode_i_fracture_energy_j_m2,
            mode_ii_fracture_energy_j_m2: 0.03 * p.mode_ii_fracture_energy_j_m2,
            ..p
        },
        cold.settings(),
    )
    .unwrap();
    let strain: Voigt = [0.012, 0.0002, -0.0001, 0.0003, 0.0004, -0.0002];
    let old = cold.trial(&strain, &cold.initial_state()).unwrap().state;
    assert!(old.damage() > 0.0 && !old.failed());
    let heated = old.transfer_temperature(&hot, 1e-6).unwrap();
    assert!(heated.state.failed());
    let released = hot
        .release_failed_plane(&heated.state, hot.normal(), 0.0, 0.0)
        .unwrap();
    assert_eq!(
        released.fracture_dissipation_j_m3,
        heated.state.dissipated_j_m3()
    );
    assert!(
        released.fracture_dissipation_j_m3
            > hot.parameters().mode_i_fracture_energy_j_m2 / hot.band_width_m()
    );
    let cooled = heated.state.transfer_temperature(&cold, 1e-6).unwrap();
    assert!(cooled.state.failed());
    let cooled_release = cold
        .release_failed_plane(&cooled.state, cold.normal(), 0.0, 0.0)
        .unwrap();
    assert_eq!(
        cooled_release.fracture_dissipation_j_m3,
        released.fracture_dissipation_j_m3
    );
    near(
        cooled_release.parallel_stored_energy_j_m3,
        2.0 * released.parallel_stored_energy_j_m3,
        1.0,
    );
    assert_eq!(cooled.damage_dissipated_j_m3, 0.0);
}
