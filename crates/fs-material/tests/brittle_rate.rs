//! Physical fracture-plane viscosity retains native irreversible fracture work.
use fs_material::brittle::{
    BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings, BrittleError,
};
use fs_material::cohesive::CohesiveViscosity;
use fs_material::elastic::IsotropicElastic;
use fs_material::tensor::{contract, rotate, rotation};
use fs_material::{SmallStrainLaw, Tangent6, Voigt};

const E: f64 = 3e9;
const ETA_NORMAL: f64 = 2e8;
const ETA_TANGENTIAL: f64 = 1e8;

fn intrinsic(width: f64, normal: [f64; 3]) -> BrittleCrackBand {
    BrittleCrackBand::new(
        IsotropicElastic::new(E, 0.25, 0.05).unwrap(),
        BrittleCrackBandParameters {
            tensile_strength_pa: 3e6,
            shear_strength_pa: 6e6,
            mode_i_fracture_energy_j_m2: 120.0,
            mode_ii_fracture_energy_j_m2: 360.0,
            bk_exponent: 1.6,
        },
        BrittleCrackBandSettings {
            band_width_m: width,
            normal,
            maximum_strain: 0.05,
            energy_tolerance_j_m3: 1e-6,
        },
    )
    .unwrap()
}

fn law(width: f64, normal: [f64; 3]) -> BrittleCrackBand {
    intrinsic(width, normal)
        .with_viscosity(CohesiveViscosity::new(ETA_NORMAL, ETA_TANGENTIAL, 1.0).unwrap())
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

fn action(tangent: &Tangent6, direction: &Voigt) -> Voigt {
    core::array::from_fn(|i| (0..6).map(|j| tangent[i][j] * direction[j]).sum())
}

#[test]
fn g1_intact_total_stress_work_and_tangent_keep_physical_interface_units() {
    let end = [1e-4, 2e-5, -1e-5, 3e-5, 2e-5, -4e-5];
    let dt = 0.02;
    for width in [0.004, 0.008] {
        let band = law(width, [1.0, 0.0, 0.0]);
        let initial = band.initial_state();
        let (path, mean_tangent) = band
            .trial_over_time_with_mean_stress_tangent(&end, &initial, dt)
            .unwrap();
        let step = &path.endpoint;
        let elastic = band.elastic().stress(&end, &());
        let (lambda, mu) = band.elastic().lame();
        let ratio = lambda / (lambda + 2.0 * mu);
        let normal_mode = [1.0, ratio, ratio, 0.0, 0.0, 0.0];
        let normal_strain = contract(&normal_mode, &end);
        let normal_traction = ETA_NORMAL * width * normal_strain / dt;
        let mut viscous: Voigt = normal_mode.map(|value| value * normal_traction);
        viscous[3] = ETA_TANGENTIAL * width * 2.0 * end[3] / dt;
        viscous[5] = ETA_TANGENTIAL * width * 2.0 * end[5] / dt;
        close_tensor(
            &step.stress_pa,
            &core::array::from_fn(|i| elastic[i] + viscous[i]),
            3e-12,
            1.0,
        );
        close_tensor(
            &path.mean_stress_pa,
            &core::array::from_fn(|i| 0.5 * elastic[i] + viscous[i]),
            3e-12,
            1.0,
        );
        let expected_loss = width / dt
            * (ETA_NORMAL * normal_strain * normal_strain
                + 4.0 * ETA_TANGENTIAL * (end[3] * end[3] + end[5] * end[5]));
        close(step.viscous_dissipated_j_m3, expected_loss, 3e-12, 1.0);
        close(
            path.path_work_j_m3,
            0.5 * contract(&elastic, &end) + expected_loss,
            3e-12,
            1.0,
        );
        assert_eq!(step.equilibrium.energy.dissipated, 0.0);
        assert_eq!(step.equilibrium.state.damage(), 0.0);
        let endpoint_tangent = step.tangent_pa.expect("intact tensile endpoint is smooth");
        let mean_tangent = mean_tangent.expect("intact tensile path is smooth");
        let stiffness = band.elastic().stiffness();
        for i in 0..6 {
            for j in 0..6 {
                let mut viscous_tangent = ETA_NORMAL * width / dt * normal_mode[i] * normal_mode[j];
                if i == j && [3, 5].contains(&i) {
                    viscous_tangent += 2.0 * ETA_TANGENTIAL * width / dt;
                }
                close(
                    endpoint_tangent[i][j],
                    stiffness[i][j] + viscous_tangent,
                    3e-12,
                    E,
                );
                close(
                    mean_tangent[i][j],
                    0.5 * stiffness[i][j] + viscous_tangent,
                    3e-12,
                    E,
                );
            }
        }

        // An accepted endpoint velocity is independent of this interval's
        // secant velocity. Observing it must not advance damage or charge work.
        let accepted = step.equilibrium.state;
        let saved = accepted;
        close_tensor(
            &band.stress_at_rate(&accepted, &[0.0; 6]).unwrap(),
            &elastic,
            3e-12,
            1.0,
        );
        close_tensor(
            &band
                .stress_at_rate(&accepted, &end.map(|value| value / dt))
                .unwrap(),
            &step.stress_pa,
            3e-12,
            1.0,
        );
        close_tensor(
            &band
                .stress_at_rate(&accepted, &end.map(|value| 2.0 * value / dt))
                .unwrap(),
            &core::array::from_fn(|i| elastic[i] + 2.0 * viscous[i]),
            3e-12,
            1.0,
        );
        assert_eq!(accepted, saved);
        assert_eq!(initial, band.initial_state());
    }
}

#[test]
fn g3_rupture_keeps_intrinsic_toughness_separate_from_integrated_viscous_loss() {
    let width = 0.004;
    let dt = 0.02;
    for shear in [false, true] {
        let band = law(width, [1.0, 0.0, 0.0]);
        let mut end = [0.0; 6];
        end[if shear { 3 } else { 0 }] = if shear { 0.018 } else { 0.024 };
        let toughness = if shear { 360.0 } else { 120.0 };
        let initial = band.initial_state();
        let path = band
            .trial_over_time_with_mean_stress(&end, &initial, dt)
            .unwrap();
        assert!(path.endpoint.equilibrium.state.failed());
        assert_eq!(path.endpoint.stress_pa, [0.0; 6]);
        close(
            path.endpoint.equilibrium.energy.dissipated * width,
            toughness,
            3e-12,
            1.0,
        );
        assert!(path.endpoint.viscous_dissipated_j_m3 > 0.0);
        close(
            path.path_work_j_m3,
            toughness / width + path.endpoint.viscous_dissipated_j_m3,
            3e-12,
            1.0,
        );

        // Refining the same constant-rate history must preserve the complete
        // analytic work, including the interval that first ruptures the bond.
        let mut state = initial;
        let mut work = 0.0;
        let mut viscous_loss = 0.0;
        for i in 1..=32 {
            let strain = end.map(|value| value * f64::from(i) / 32.0);
            let step = band
                .trial_over_time_with_mean_stress(&strain, &state, dt / 32.0)
                .unwrap();
            assert!(step.endpoint.equilibrium.state.damage() >= state.damage());
            work += step.path_work_j_m3;
            viscous_loss += step.endpoint.viscous_dissipated_j_m3;
            state = step.endpoint.equilibrium.state;
        }
        assert!(state.failed());
        close(work, path.path_work_j_m3, 3e-11, 1.0);
        close(
            viscous_loss,
            path.endpoint.viscous_dissipated_j_m3,
            3e-11,
            1.0,
        );
        close(state.dissipated_j_m3() * width, toughness, 3e-12, 1.0);

        // Once failure was accepted, detached opening/shear and independent
        // compression are not constrained by the defunct adhesive rate cap.
        let reopened = [0.025, 0.0, 0.0, 0.003, 0.0, -0.001];
        let free = band.trial_over_time(&reopened, &state, 1e-9).unwrap();
        assert_eq!(free.stress_pa, [0.0; 6]);
        assert_eq!(free.viscous_dissipated_j_m3, 0.0);
        let released = band
            .release_failed_plane(&free.equilibrium.state, [1.0, 0.0, 0.0], 1e-12, 1e-6)
            .unwrap();
        assert_eq!(released.stress_pa, [0.0; 6]);
        close(
            released.fracture_dissipation_j_m3 * width,
            toughness,
            3e-12,
            1.0,
        );
        let compressed = [-0.01, 0.0025, 0.0025, 0.0, 0.0, 0.0];
        let contact = band
            .trial_over_time(&compressed, &free.equilibrium.state, 1e-9)
            .unwrap();
        close_tensor(
            &contact.stress_pa,
            &band.elastic().stress(&compressed, &()),
            3e-12,
            E,
        );
        assert_eq!(contact.viscous_dissipated_j_m3, 0.0);
        assert_eq!(contact.equilibrium.energy.dissipated, 0.0);
        assert!(
            band.release_failed_plane(&contact.equilibrium.state, [1.0, 0.0, 0.0], 1e-12, E,)
                .is_err()
        );
        close_tensor(
            &band
                .stress_at_rate(&contact.equilibrium.state, &[1e300; 6])
                .unwrap(),
            &contact.stress_pa,
            3e-12,
            E,
        );
        assert!(
            band.stress_at_rate(&contact.equilibrium.state, &[f64::INFINITY; 6])
                .is_err()
        );
        assert_eq!(initial.damage(), 0.0);
    }
}

#[test]
fn g1_mixed_endpoint_and_mean_tangents_differentiate_integrity_and_physical_rate() {
    let band = law(0.004, [1.0, 2.0, 3.0]);
    let dt = 0.03;
    let initial = band
        .trial_over_time(
            &[0.0002, 0.0, 0.0, 0.00004, 0.0, 0.0],
            &band.initial_state(),
            dt,
        )
        .unwrap()
        .equilibrium
        .state;
    let end = [0.006, -0.001, 0.0005, 0.002, 0.0008, -0.001];
    let damaged = band
        .trial_over_time(&end, &initial, dt)
        .unwrap()
        .equilibrium
        .state;
    assert!(damaged.damage() > 0.0 && damaged.damage() < 1.0);
    let endpoint_law = band.small_strain_law_over_time(dt).unwrap();
    let mean_law = band.mean_small_strain_law_over_time(dt).unwrap();
    assert_eq!(endpoint_law.duration_s(), dt);
    assert_eq!(mean_law.duration_s(), dt);
    for (old, target) in [
        (initial, end),
        (damaged, end.map(|value| 1.25 * value)),
        (damaged, end.map(|value| 0.65 * value)),
    ] {
        let saved = old;
        let (path, mean_tangent) = band
            .trial_over_time_with_mean_stress_tangent(&target, &old, dt)
            .unwrap();
        let endpoint_tangent = path
            .endpoint
            .tangent_pa
            .expect("selected endpoint is smooth");
        let mean_tangent = mean_tangent.expect("selected mixed path has no held switch");
        for j in 0..6 {
            let mut plus = target;
            let mut minus = target;
            let h = 1e-8;
            plus[j] += h;
            minus[j] -= h;
            let plus = band
                .trial_over_time_with_mean_stress(&plus, &old, dt)
                .unwrap();
            let minus = band
                .trial_over_time_with_mean_stress(&minus, &old, dt)
                .unwrap();
            for i in 0..6 {
                close(
                    endpoint_tangent[i][j],
                    (plus.endpoint.stress_pa[i] - minus.endpoint.stress_pa[i]) / (2.0 * h),
                    3e-6,
                    3e3,
                );
                close(
                    mean_tangent[i][j],
                    (plus.mean_stress_pa[i] - minus.mean_stress_pa[i]) / (2.0 * h),
                    3e-6,
                    3e3,
                );
            }
        }
        assert_eq!(endpoint_law.stress(&target, &old), path.endpoint.stress_pa);
        assert_eq!(endpoint_law.tangent(&target, &old), endpoint_tangent);
        assert_eq!(mean_law.stress(&target, &old), path.mean_stress_pa);
        assert_eq!(mean_law.tangent(&target, &old), mean_tangent);
        assert_eq!(
            endpoint_law.update_state(&target, &old),
            path.endpoint.equilibrium.state
        );
        assert_eq!(
            mean_law.update_state(&target, &old),
            path.endpoint.equilibrium.state
        );
        assert_eq!(old, saved);
    }
}

#[test]
fn g3_rotated_rate_histories_preserve_total_work_stress_and_tangent_action() {
    let axis = [1.0, 2.0, -1.0].map(|value| value / fs_math::det::sqrt(6.0));
    let q = rotation(axis, 0.73);
    let a = law(0.004, [1.0, 0.0, 0.0]);
    let b = law(0.004, [q[0][0], q[1][0], q[2][0]]);
    let mut sa = a.initial_state();
    let mut sb = b.initial_state();
    let direction = [0.0002, -0.0001, 0.00005, 0.00003, 0.00002, -0.00004];
    for end in [
        [0.0001, 0.0, 0.0, 0.0002, 0.0001, 0.0003],
        [0.006, -0.001, 0.0005, 0.002, 0.0008, -0.001],
        [0.003, 0.001, 0.0005, -0.003, 0.0005, 0.002],
    ] {
        let (pa, ta) = a
            .trial_over_time_with_mean_stress_tangent(&end, &sa, 0.02)
            .unwrap();
        let (pb, tb) = b
            .trial_over_time_with_mean_stress_tangent(&rotate(&end, &q), &sb, 0.02)
            .unwrap();
        close_tensor(
            &pb.endpoint.stress_pa,
            &rotate(&pa.endpoint.stress_pa, &q),
            3e-11,
            1e6,
        );
        close_tensor(
            &pb.mean_stress_pa,
            &rotate(&pa.mean_stress_pa, &q),
            3e-11,
            1e6,
        );
        close_tensor(
            &action(&tb.unwrap(), &rotate(&direction, &q)),
            &rotate(&action(&ta.unwrap(), &direction), &q),
            3e-10,
            1e3,
        );
        close(pa.path_work_j_m3, pb.path_work_j_m3, 3e-11, 1.0);
        close(
            pa.endpoint.viscous_dissipated_j_m3,
            pb.endpoint.viscous_dissipated_j_m3,
            3e-11,
            1.0,
        );
        sa = pa.endpoint.equilibrium.state;
        sb = pb.endpoint.equilibrium.state;
        close_tensor(
            &b.stress_at_rate(&sb, &rotate(&direction, &q)).unwrap(),
            &rotate(&a.stress_at_rate(&sa, &direction).unwrap(), &q),
            3e-11,
            1e6,
        );
    }
}

#[test]
fn g0_duration_rate_owner_and_untimed_refusals_leave_original_history_untouched() {
    let band = law(0.004, [1.0, 0.0, 0.0]);
    let initial = band.initial_state();
    let end = [0.024, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(band.trial(&end, &initial).is_err());
    assert!(band.trial_with_mean_stress(&end, &initial).is_err());
    assert!(band.trial_with_mean_stress_tangent(&end, &initial).is_err());
    let partial = band
        .trial_over_time(&[0.006, 0.0, 0.0, 0.0, 0.0, 0.0], &initial, 0.02)
        .unwrap()
        .equilibrium
        .state;
    assert!(
        band.release_failed_plane(&partial, [1.0, 0.0, 0.0], 1e-12, E)
            .is_err()
    );
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(band.trial_over_time(&end, &initial, dt).is_err());
        assert!(
            band.trial_over_time_with_mean_stress(&end, &initial, dt)
                .is_err()
        );
        assert!(band.small_strain_law_over_time(dt).is_err());
        assert!(band.mean_small_strain_law_over_time(dt).is_err());
    }
    // Newly reached failure cannot bypass the rate bound of the original
    // surviving bonds over the interval that breaks them.
    assert!(band.trial_over_time(&end, &initial, 1e-6).is_err());
    assert!(
        band.stress_at_rate(&initial, &[1e4, 0.0, 0.0, 0.0, 0.0, 0.0])
            .is_err()
    );
    assert!(band.stress_at_rate(&initial, &[f64::NAN; 6]).is_err());
    let clone = band.clone();
    let foreign = clone.initial_state();
    assert_eq!(
        band.trial_over_time(&end, &foreign, 0.02),
        Err(BrittleError::StateMismatch)
    );
    assert_eq!(
        band.stress_at_rate(&foreign, &[0.0; 6]),
        Err(BrittleError::StateMismatch)
    );
    let endpoint = band.small_strain_law_over_time(0.02).unwrap();
    assert!(
        endpoint
            .stress(&end, &foreign)
            .iter()
            .all(|value| value.is_nan())
    );
    assert!(
        endpoint
            .update_state(&end, &foreign)
            .strain()
            .iter()
            .all(|value| value.is_nan())
    );
    let mean = band.mean_small_strain_law_over_time(0.02).unwrap();
    let (held, tangent) = mean.trial(initial.strain(), &initial).unwrap();
    assert!(tangent.is_none());
    assert!(
        mean.tangent(initial.strain(), &initial)
            .iter()
            .flatten()
            .all(|value| value.is_nan())
    );
    assert_eq!(held.path_work_j_m3, 0.0);
    assert_eq!(held.endpoint.viscous_dissipated_j_m3, 0.0);
    assert_eq!(initial, band.initial_state());
}

#[test]
fn g0_absent_viscosity_retains_original_untimed_response_and_history() {
    let band = intrinsic(0.004, [1.0, 0.0, 0.0]);
    assert!(band.viscosity().is_none());
    let mut state = band.initial_state();
    for end in [
        [0.006, 0.0, 0.0, 0.002, 0.0, 0.0],
        [0.003, 0.0, 0.0, -0.001, 0.0, 0.0],
        [0.03, 0.0, 0.0, 0.003, 0.0, 0.0],
        [-0.01, 0.0025, 0.0025, 0.0, 0.0, 0.0],
    ] {
        let original = band.trial_with_mean_stress(&end, &state).unwrap();
        let timed = band
            .trial_over_time_with_mean_stress(&end, &state, 0.02)
            .unwrap();
        assert_eq!(timed.endpoint.equilibrium, original.endpoint);
        assert_eq!(timed.endpoint.stress_pa, original.endpoint.stress_pa);
        assert_eq!(timed.mean_stress_pa, original.mean_stress_pa);
        assert_eq!(timed.path_work_j_m3, original.path_work_j_m3);
        assert_eq!(timed.endpoint.viscous_dissipated_j_m3, 0.0);
        assert_eq!(
            timed.endpoint.total_work_j_m3(),
            original.endpoint.energy.work
        );
        state = timed.endpoint.equilibrium.state;
        assert_eq!(
            band.stress_at_rate(&state, &[1e300; 6]).unwrap(),
            original.endpoint.stress_pa
        );
    }
}
