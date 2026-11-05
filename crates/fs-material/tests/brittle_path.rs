//! G0/G1/G3: exact fracture-path work and derivatives of its work-conjugate stress.
use fs_material::brittle::{
    BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings, BrittleError,
    BrittleState,
};
use fs_material::elastic::IsotropicElastic;
use fs_material::tensor::{contract, rotate, rotation};
use fs_material::{SmallStrainLaw, Tangent6, Voigt};

const E: f64 = 3.0e9;
const NU: f64 = 0.25;

fn parameters() -> BrittleCrackBandParameters {
    BrittleCrackBandParameters {
        tensile_strength_pa: 3e6,
        shear_strength_pa: 6e6,
        mode_i_fracture_energy_j_m2: 120.0,
        mode_ii_fracture_energy_j_m2: 360.0,
        bk_exponent: 1.6,
    }
}
fn settings(width: f64, normal: [f64; 3]) -> BrittleCrackBandSettings {
    BrittleCrackBandSettings {
        band_width_m: width,
        normal,
        maximum_strain: 0.05,
        energy_tolerance_j_m3: 1e-6,
    }
}
fn law(width: f64, normal: [f64; 3]) -> BrittleCrackBand {
    BrittleCrackBand::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        parameters(),
        settings(width, normal),
    )
    .unwrap()
}
fn close(a: f64, b: f64, relative: f64, floor: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(floor),
        "{a:e} != {b:e}"
    );
}
fn close_tensor(a: &Voigt, b: &Voigt, relative: f64, floor: f64) {
    for i in 0..6 {
        close(a[i], b[i], relative, floor);
    }
}
fn elastic_strain(stress: Voigt, e: f64, nu: f64) -> Voigt {
    [
        (stress[0] - nu * (stress[1] + stress[2])) / e,
        (stress[1] - nu * (stress[0] + stress[2])) / e,
        (stress[2] - nu * (stress[0] + stress[1])) / e,
        (1.0 + nu) * stress[3] / e,
        (1.0 + nu) * stress[4] / e,
        (1.0 + nu) * stress[5] / e,
    ]
}
fn action(tangent: &Tangent6, direction: &Voigt) -> Voigt {
    core::array::from_fn(|i| (0..6).map(|j| tangent[i][j] * direction[j]).sum())
}

#[test]
fn g1_intact_mean_is_exact_elastic_average_and_its_tangent_is_half_stiffness() {
    let old_stress = [2e4, 1.3e4, 1.5e4, 2e3, -5e2, 7e2];
    let new_stress = [8e4, 3e4, 4e4, -4e3, 1e3, 2e3];
    for nu in [-0.8, 0.0, 0.25, 0.49] {
        for normal in [[1.0, 0.0, 0.0], [1.0, 2.0, -3.0]] {
            let elastic = IsotropicElastic::new(E, nu, 0.05).unwrap();
            let band =
                BrittleCrackBand::new(elastic.clone(), parameters(), settings(0.004, normal))
                    .unwrap();
            let a = elastic_strain(old_stress, E, nu);
            let b = elastic_strain(new_stress, E, nu);
            let old = band.trial(&a, &band.initial_state()).unwrap().state;
            let saved = old;
            let (path, tangent) = band.trial_with_mean_stress_tangent(&b, &old).unwrap();
            assert_eq!(path, band.trial_with_mean_stress(&b, &old).unwrap());
            assert_eq!(path.endpoint, band.trial(&b, &old).unwrap());
            assert_eq!(old, saved);
            close_tensor(
                &path.mean_stress_pa,
                &core::array::from_fn(|i| 0.5 * (old_stress[i] + new_stress[i])),
                3e-12,
                1.0,
            );
            let tangent = tangent.expect("strictly tensile intact path is smooth");
            let expected = elastic.stiffness();
            for i in 0..6 {
                for j in 0..6 {
                    close(tangent[i][j], 0.5 * expected[i][j], 3e-12, E);
                }
            }
            let delta_storage = 0.5 * (contract(&new_stress, &b) - contract(&old_stress, &a));
            close(path.path_work_j_m3, delta_storage, 3e-12, 1.0);
            assert_eq!(path.endpoint.energy.dissipated, 0.0);
            assert!(path.path_work_defect_j_m3.abs() <= 1e-6);
        }
    }
}

#[test]
fn g1_one_overshooting_step_dissipates_gc_per_width_despite_zero_endpoint_traction() {
    for width in [0.004, 0.008, 0.012] {
        let band = law(width, [1.0, 0.0, 0.0]);
        for shear in [false, true] {
            let (strength, energy) = if shear { (6e6, 360.0) } else { (3e6, 120.0) };
            let mode_end = 1.4 * (2.0 * energy / (strength * width));
            let mut end = [0.0; 6];
            end[if shear { 3 } else { 0 }] = if shear { 0.5 * mode_end } else { mode_end };
            let initial = band.initial_state();
            let (path, tangent) = band.trial_with_mean_stress_tangent(&end, &initial).unwrap();
            assert!(path.endpoint.state.failed());
            assert_eq!(path.endpoint.stress_pa, [0.0; 6]);
            close(path.path_work_j_m3, energy / width, 3e-12, 1.0);
            close(path.endpoint.energy.dissipated, energy / width, 3e-12, 1.0);
            close(
                path.mean_stress_pa[if shear { 3 } else { 0 }],
                energy / (width * mode_end),
                3e-12,
                1.0,
            );
            assert!(path.mean_stress_pa[if shear { 3 } else { 0 }] > 0.0);
            if shear {
                // The full classical derivative is absent at a held normal-contact
                // switch, even though this pure-shear force/work is well defined.
                assert!(tangent.is_none());
            } else {
                close(
                    tangent.unwrap()[0][0],
                    -energy / (width * mode_end * mode_end),
                    3e-12,
                    1.0,
                );
            }
            assert_eq!(initial.damage(), 0.0);
        }
    }
}

#[test]
fn g3_unloading_reversal_failure_and_reclosure_keep_exact_work_and_old_history() {
    let band = law(0.004, [1.0, 0.0, 0.0]);
    let mut state = band.initial_state();
    let mut path_work = 0.0;
    let mut dissipated = 0.0;
    for strain in [
        [0.012, 0.0, 0.0, 0.002, 0.0, 0.0],
        [0.004, 0.0, 0.0, 0.0005, 0.0, 0.0],
        [0.01, 0.0002, -0.0001, -0.002, 0.0, 0.001],
        [0.028, 0.0, 0.0, 0.002, 0.0, 0.0],
        [-0.006, 0.0015, 0.0015, 0.0, 0.0, 0.0],
        [0.03, 0.0, 0.0, 0.003, 0.0, 0.0],
    ] {
        let saved = state;
        let path = band.trial_with_mean_stress(&strain, &state).unwrap();
        assert_eq!(state, saved);
        assert_eq!(path.endpoint, band.trial(&strain, &state).unwrap());
        assert!(path.endpoint.state.damage() >= state.damage());
        let change = core::array::from_fn(|i| strain[i] - state.strain()[i]);
        close(
            contract(&path.mean_stress_pa, &change),
            path.path_work_j_m3,
            1e-13,
            1.0,
        );
        close(
            path.path_work_j_m3,
            path.endpoint.state.stored_energy_j_m3() - state.stored_energy_j_m3()
                + path.endpoint.energy.dissipated,
            2e-11,
            1.0,
        );
        path_work += path.path_work_j_m3;
        dissipated += path.endpoint.energy.dissipated;
        state = path.endpoint.state;
    }
    assert!(state.failed());
    close(state.stored_energy_j_m3(), 0.0, 1e-12, 1.0);
    close(path_work, dissipated, 2e-11, 1.0);
    close(dissipated, state.dissipated_j_m3(), 2e-11, 1.0);
}

fn check_tangent(band: &BrittleCrackBand, old: &BrittleState<'_>, end: Voigt) {
    let (path, tangent) = band.trial_with_mean_stress_tangent(&end, old).unwrap();
    let tangent = tangent.expect("selected mixed path has no held switch");
    let h = 1e-8;
    for j in 0..6 {
        let mut plus = end;
        let mut minus = end;
        plus[j] += h;
        minus[j] -= h;
        let plus = band
            .trial_with_mean_stress(&plus, old)
            .unwrap()
            .mean_stress_pa;
        let minus = band
            .trial_with_mean_stress(&minus, old)
            .unwrap()
            .mean_stress_pa;
        for i in 0..6 {
            close(tangent[i][j], (plus[i] - minus[i]) / (2.0 * h), 3e-6, 3e3);
        }
    }
    assert_eq!(path.endpoint, band.trial(&end, old).unwrap());
}

#[test]
fn g1_mean_tangent_differentiates_first_bk_initiation_and_frozen_mixed_history() {
    let band = law(0.004, [1.0, 2.0, 3.0]);
    let old = band
        .trial(
            &[0.0002, 0.0, 0.0, 0.00004, 0.0, 0.0],
            &band.initial_state(),
        )
        .unwrap()
        .state;
    assert_eq!(old.damage(), 0.0);
    let end = [0.006, -0.001, 0.0005, 0.002, 0.0008, -0.001];
    check_tangent(&band, &old, end);
    let damaged = band.trial(&end, &old).unwrap().state;
    assert!(damaged.damage() > 0.0 && damaged.damage() < 1.0);
    check_tangent(&band, &damaged, end.map(|value| 1.25 * value));
    check_tangent(&band, &damaged, end.map(|value| 0.65 * value));
}

#[test]
fn g1_mixed_mean_matches_independent_composite_gauss_integration_of_endpoint_stress() {
    let band = law(0.004, [1.0, 2.0, -1.0]);
    let start = [0.0002, 0.0001, 0.0, 0.00004, 0.0, 0.0];
    let end = [0.008, 0.001, -0.0004, 0.002, -0.0003, 0.001];
    let old = band.trial(&start, &band.initial_state()).unwrap().state;
    let path = band.trial_with_mean_stress(&end, &old).unwrap();
    assert!(path.endpoint.state.damage() > 0.0);
    // Independent two-point Gauss quadrature of the actual point response.
    // This samples native prefix histories; it never calls the path integral.
    let mut integral = [0.0; 6];
    let segments = 1024;
    let gauss = 1.0 / fs_math::det::sqrt(3.0);
    for interval in 0..segments {
        for sign in [-1.0, 1.0] {
            let fraction = (f64::from(interval) + 0.5 * (1.0 + sign * gauss)) / f64::from(segments);
            let strain = core::array::from_fn(|i| (1.0 - fraction) * start[i] + fraction * end[i]);
            let point = band.trial(&strain, &old).unwrap();
            for i in 0..6 {
                integral[i] += point.stress_pa[i] / (2.0 * f64::from(segments));
            }
        }
    }
    close_tensor(&path.mean_stress_pa, &integral, 2e-6, 1e5);
}

#[test]
fn g3_rotating_the_plane_and_strain_path_rotates_mean_stress_and_tangent_actions() {
    let axis = [1.0, 2.0, -1.0].map(|value| value / fs_math::det::sqrt(6.0));
    let q = rotation(axis, 0.73);
    let a = law(0.004, [1.0, 0.0, 0.0]);
    let b = law(0.004, [q[0][0], q[1][0], q[2][0]]);
    let mut sa = a.initial_state();
    let mut sb = b.initial_state();
    let direction = [0.0003, -0.0002, 0.0001, 0.00015, -0.00007, 0.00012];
    for end in [
        [0.004, 0.0002, -0.0001, 0.0007, 0.0002, -0.0003],
        [0.009, -0.0003, 0.0004, -0.001, 0.0002, 0.002],
        [0.002, 0.0001, 0.0, 0.0002, -0.0001, 0.0],
        [-0.004, 0.0002, -0.0001, 0.001, 0.0001, -0.0007],
    ] {
        let (pa, ja) = a.trial_with_mean_stress_tangent(&end, &sa).unwrap();
        let (pb, jb) = b
            .trial_with_mean_stress_tangent(&rotate(&end, &q), &sb)
            .unwrap();
        close_tensor(
            &pb.mean_stress_pa,
            &rotate(&pa.mean_stress_pa, &q),
            2e-10,
            1e5,
        );
        close(pa.path_work_j_m3, pb.path_work_j_m3, 2e-10, 1.0);
        close(
            pa.endpoint.state.damage(),
            pb.endpoint.state.damage(),
            2e-10,
            1.0,
        );
        let expected = rotate(&action(&ja.unwrap(), &direction), &q);
        let actual = action(&jb.unwrap(), &rotate(&direction, &q));
        close_tensor(&actual, &expected, 2e-9, 1e4);
        sa = pa.endpoint.state;
        sb = pb.endpoint.state;
    }
}

#[test]
fn g3_temperature_transferred_history_preserves_accumulated_loss_on_the_next_path() {
    let cold = law(0.004, [1.0, 0.0, 0.0]);
    let warm = BrittleCrackBand::new(
        IsotropicElastic::new(0.7 * E, NU, 0.05).unwrap(),
        BrittleCrackBandParameters {
            mode_i_fracture_energy_j_m2: 90.0,
            mode_ii_fracture_energy_j_m2: 270.0,
            ..parameters()
        },
        settings(0.004, [1.0, 0.0, 0.0]),
    )
    .unwrap();
    let loaded = cold
        .trial(
            &[0.008, 0.0002, -0.0001, 0.001, 0.0, 0.0],
            &cold.initial_state(),
        )
        .unwrap()
        .state;
    let thermal = loaded.transfer_temperature(&warm, 1e-6).unwrap();
    let old_loss = thermal.state.dissipated_j_m3();
    let (path, _) = warm
        .trial_with_mean_stress_tangent(
            &[0.014, 0.0003, -0.0002, -0.002, 0.0, 0.001],
            &thermal.state,
        )
        .unwrap();
    close(
        path.endpoint.state.dissipated_j_m3() - old_loss,
        path.endpoint.energy.dissipated,
        2e-11,
        1.0,
    );
    close(
        path.path_work_j_m3,
        path.endpoint.state.stored_energy_j_m3() - thermal.state.stored_energy_j_m3()
            + path.endpoint.energy.dissipated,
        2e-11,
        1.0,
    );
    assert_eq!(thermal.state.dissipated_j_m3(), old_loss);
}

#[test]
fn g0_owner_and_frobenius_refusals_leave_accepted_history_reusable() {
    let band = law(0.004, [1.0, 0.0, 0.0]);
    let initial = band.initial_state();
    let foreign = band.clone();
    assert_eq!(
        foreign.trial_with_mean_stress(&[0.0; 6], &initial),
        Err(BrittleError::StateMismatch)
    );
    assert_eq!(
        foreign.trial_with_mean_stress_tangent(&[0.0; 6], &initial),
        Err(BrittleError::StateMismatch)
    );
    for invalid in [[0.02; 6], [f64::NAN; 6], [f64::INFINITY; 6]] {
        assert!(band.trial_with_mean_stress(&invalid, &initial).is_err());
        assert!(
            band.trial_with_mean_stress_tangent(&invalid, &initial)
                .is_err()
        );
    }
    let marker = band.small_strain_law().update_state(&[0.02; 6], &initial);
    assert!(!marker.is_valid());
    assert!(band.trial_with_mean_stress(&[0.0; 6], &marker).is_err());
    assert!(
        band.trial_with_mean_stress_tangent(&[0.0; 6], &marker)
            .is_err()
    );
    assert_eq!(initial, band.initial_state());
    let held = band.trial_with_mean_stress(&[0.0; 6], &initial).unwrap();
    assert_eq!(held.path_work_j_m3, 0.0);
    assert_eq!(held.path_work_defect_j_m3, 0.0);
    assert_eq!(held.endpoint.state, initial);
}
