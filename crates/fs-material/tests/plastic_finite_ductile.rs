//! Native multiplicative plasticity, full damage derivative and physical work.
//! Synthetic coefficients verify the law, not an empirical alloy calibration.
#![deny(unsafe_code)]

use fs_material::ductile::evolution::{
    DifferentiableDuctileInitiationLaw, DuctileJ2Conditions, DuctileJ2Settings,
    DuctileLocusGradient,
};
use fs_material::ductile::{
    DuctileDomain, DuctileError, DuctileInitiationLaw, DuctileStressState, TabulatedDuctileLocus,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::{FiniteDuctileError, FiniteDuctileJ2};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};

fn native() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: 500.0,
            bulk_modulus_pa: 1000.0,
            yield_stress_pa: 3.0,
            hardening_modulus_pa: 40.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}
fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-50.0, 50.0],
        strain_rate_s_inv: [1e-12, 1000.0],
        temperature_k: [300.0, 300.0],
    }
}
fn settings() -> DuctileJ2Settings {
    DuctileJ2Settings {
        characteristic_length_m: 0.02,
        post_initiation_plastic_work_j_m2: 0.004,
        maximum_damage_increment: 1.0,
        energy_tolerances: EnergyTolerances {
            energy_j_m3: 1e-10,
            yield_pa: 1e-8,
        },
    }
}
fn conditions() -> DuctileJ2Conditions {
    DuctileJ2Conditions {
        duration_s: 0.2,
        temperature_k: 300.0,
    }
}
fn material(onset: f64) -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    FiniteDuctileJ2::new(
        native(),
        TabulatedDuctileLocus::new(vec![(-50.0, onset), (50.0, onset)], domain()).unwrap(),
        settings(),
    )
    .unwrap()
}
fn coaxial(e: f64, dilation: f64) -> [f64; 9] {
    [
        (dilation + e).exp(),
        0.0,
        0.0,
        0.0,
        (dilation - 0.5 * e).exp(),
        0.0,
        0.0,
        0.0,
        (dilation - 0.5 * e).exp(),
    ]
}
fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}
fn transpose(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|i| a[3 * (i % 3) + i / 3])
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}; difference {:.8e}",
        (actual - expected).abs()
    );
}

#[test]
fn intact_branch_preserves_original_metric_stress_storage_and_all_tangent_entries() {
    let m = material(10.0);
    let mut old = m.initial_state();
    for e in [0.001, 0.014, 0.028, 0.027] {
        let f = coaxial(e, 0.001);
        let original = m.native().response(&f, old.native_state()).unwrap();
        let r = m.trial(&f, &old, conditions()).unwrap();
        assert_eq!(r.state.native_state(), &original.state);
        assert_eq!(r.first_piola_pa, original.first_piola_pa);
        assert_eq!(r.kirchhoff_pa, original.kirchhoff_pa);
        assert_eq!(r.tangent_pa, original.tangent_pa);
        assert_eq!(r.state.stored_energy_j_m3(), original.stored_energy_j_m3);
        assert_eq!(
            r.energy.plastic_dissipation_j_m3,
            original.plastic_dissipation_j_m3
        );
        assert_eq!(m.stress(&f, &old, conditions()).unwrap(), r.first_piola_pa);
        assert_eq!(r.state.damage(), 0.0);
        old = r.state;
    }
}

#[test]
fn linear_and_voce_attenuation_match_independent_integrated_kirchhoff_resistance() {
    for saturation in [0.0, 15.0] {
        let onset = 0.002;
        let m = FiniteDuctileJ2::new(
            native().with_voce_hardening(saturation, 12.0).unwrap(),
            TabulatedDuctileLocus::new(vec![(-50.0, onset), (50.0, onset)], domain()).unwrap(),
            settings(),
        )
        .unwrap();
        let mut state = m.initial_state();
        let mut equivalent = 0.0;
        let mut degraded = 0.0;
        for index in 1..=20 {
            let f = coaxial(0.002 + 0.04 * index as f64 / 20.0, 0.0);
            let r = m.trial(&f, &state, conditions()).unwrap();
            let alpha = r.state.native_state().alpha;
            // Independent scalar endpoint yield relation in logarithmic strain.
            near(
                1500.0 * (0.002 + 0.04 * index as f64 / 20.0 - alpha),
                3.0 + 40.0 * alpha + saturation * (1.0 - (-12.0 * alpha).exp()),
                2e-11,
                2e-11,
            );
            let delta = (alpha - onset).max(0.0);
            let work = delta * (3.0 + 40.0 * (onset + 0.5 * delta))
                + saturation * (delta + (-12.0 * onset).exp() * (-12.0 * delta).exp_m1() / 12.0);
            let expected_g = (-settings().characteristic_length_m * work
                / settings().post_initiation_plastic_work_j_m2)
                .exp();
            near(r.state.integrity(), expected_g, 3e-13, 3e-13);
            equivalent += r.energy.equivalent_plastic_work_increment_j_m2;
            degraded += r.energy.degraded_equivalent_plastic_work_increment_j_m2;
            near(
                equivalent,
                settings().characteristic_length_m * work,
                3e-15,
                3e-12,
            );
            near(
                degraded,
                settings().post_initiation_plastic_work_j_m2 * (1.0 - expected_g),
                3e-15,
                3e-12,
            );
            near(
                r.state.equivalent_plastic_work_j_m2(),
                equivalent,
                3e-15,
                3e-12,
            );
            near(
                r.state.degraded_equivalent_plastic_work_j_m2(),
                degraded,
                3e-15,
                3e-12,
            );
            assert!(
                r.energy.plastic_dissipation_j_m3 >= 0.0 && r.energy.damage_release_j_m3 >= 0.0
            );
            state = r.state;
        }
    }
}

#[test]
fn tabulated_damage_uses_integrated_resistance_across_knots_and_plateaus() {
    use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};

    let knots = [(0.0, 3.0), (0.01, 5.0), (0.03, 6.0), (0.1, 6.0), (3.0, 20.0)];
    let table = TabulatedHardening::new(knots.iter().map(|&(plastic_strain, yield_stress_pa)| {
        HardeningPoint { plastic_strain, yield_stress_pa }
    }).collect()).unwrap();
    let reference = native();
    let backbone = FiniteJ2::new(
        FiniteJ2Parameters { hardening_modulus_pa: 0.0, ..*reference.parameters() },
        *reference.limits(),
    ).unwrap().with_tabulated_hardening(table).unwrap();
    let onset = 0.002;
    let m = FiniteDuctileJ2::new(
        backbone,
        TabulatedDuctileLocus::new(vec![(-50.0, onset), (50.0, onset)], domain()).unwrap(),
        settings(),
    ).unwrap();
    let mut old = m.initial_state();
    for alpha in [0.007_f64, 0.021, 0.055, 0.08] {
        // Integrate the independently supplied line segments, including the
        // partially consumed onset segment and the later perfect-flow plateau.
        let mut resistance = 0.0;
        let mut work = 0.0;
        for pair in knots.windows(2) {
            let [(a, ya), (b, yb)] = [pair[0], pair[1]];
            let slope = (yb - ya) / (b - a);
            if alpha >= a && alpha <= b {
                resistance = ya + slope * (alpha - a);
            }
            let lower = a.max(onset);
            let upper = b.min(alpha);
            if upper > lower {
                work += (ya + slope * (0.5 * lower + 0.5 * upper - a)) * (upper - lower);
            }
        }
        let f = coaxial(alpha + resistance / 1500.0, 0.002);
        let r = m.trial(&f, &old, conditions()).unwrap();
        near(r.state.native_state().alpha, alpha, 3e-13, 3e-12);
        let equivalent = settings().characteristic_length_m * work;
        let integrity = (-equivalent / settings().post_initiation_plastic_work_j_m2).exp();
        near(r.state.equivalent_plastic_work_j_m2(), equivalent, 3e-15, 3e-12);
        near(r.state.integrity(), integrity, 3e-13, 3e-12);
        near(
            r.state.degraded_equivalent_plastic_work_j_m2(),
            settings().post_initiation_plastic_work_j_m2 * (1.0 - integrity),
            3e-15,
            3e-12,
        );
        for j in 0..9 {
            let h = 2e-7;
            let mut plus = f;
            let mut minus = f;
            plus[j] += h;
            minus[j] -= h;
            let rp = m.stress(&plus, &old, conditions()).unwrap();
            let rm = m.stress(&minus, &old, conditions()).unwrap();
            for i in 0..9 {
                near(r.tangent_pa[i][j], (rp[i] - rm[i]) / (2.0 * h), 2e-5, 2e-6);
            }
        }
        old = r.state;
    }
    let hold = m.trial(old.deformation(), &old, conditions()).unwrap();
    assert_eq!(hold.state, old);
    assert_eq!(hold.energy.equivalent_plastic_work_increment_j_m2, 0.0);
}

#[derive(Clone, Debug, PartialEq)]
struct SmoothLocus(DuctileDomain);
impl DuctileInitiationLaw for SmoothLocus {
    fn domain(&self) -> &DuctileDomain {
        &self.0
    }
    fn fracture_strain(&self, _: f64, _: f64, _: f64) -> Result<f64, DuctileError> {
        Err(DuctileError::Invalid("full stress state required"))
    }
    fn fracture_strain_at(
        &self,
        s: DuctileStressState,
        rate: f64,
        _: f64,
    ) -> Result<f64, DuctileError> {
        Ok(0.003
            * (1.0 + 0.06 * s.triaxiality + 0.15 * s.lode_parameter + 0.025 * (rate / 0.01).ln()))
    }
}
impl DifferentiableDuctileInitiationLaw for SmoothLocus {
    fn fracture_strain_gradient(
        &self,
        _: f64,
        _: f64,
        _: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        Err(DuctileError::Invalid("full stress state required"))
    }
    fn fracture_strain_gradient_at(
        &self,
        _: DuctileStressState,
        rate: f64,
        _: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        Ok(DuctileLocusGradient {
            triaxiality: 0.003 * 0.06,
            lode_parameter: 0.003 * 0.15,
            strain_rate_s_inv: 0.003 * 0.025 / rate,
        })
    }
}

#[test]
fn all_81_tangent_entries_include_moving_lode_rate_onset_voce_and_compression() {
    let m = FiniteDuctileJ2::new(
        native().with_voce_hardening(15.0, 12.0).unwrap(),
        SmoothLocus(domain()),
        settings(),
    )
    .unwrap();
    let virgin = m.initial_state();
    let before = m
        .trial(&coaxial(0.0035, 0.0), &virgin, conditions())
        .unwrap()
        .state;
    assert!(!before.initiation().initiated());
    let f = [
        1.028, 0.017, -0.006, 0.004, 0.987, 0.009, 0.001, -0.003, 0.995,
    ];
    let crossed = m.trial(&f, &before, conditions()).unwrap();
    assert!(crossed.state.damage() > 0.0);
    // Add deviatoric loading as well as compression so both perturbations
    // remain on the plastic branch, away from the elastic/yield kink.
    let compressed = multiply(&coaxial(0.012, -0.04), &f);
    for (old, endpoint) in [(&before, f), (&crossed.state, compressed)] {
        let response = m.trial(&endpoint, old, conditions()).unwrap();
        let mut asymmetry = 0.0_f64;
        for j in 0..9 {
            let h = 2e-7;
            let mut plus = endpoint;
            let mut minus = endpoint;
            plus[j] += h;
            minus[j] -= h;
            let rp = m.stress(&plus, old, conditions()).unwrap();
            let rm = m.stress(&minus, old, conditions()).unwrap();
            for i in 0..9 {
                near(
                    response.tangent_pa[i][j],
                    (rp[i] - rm[i]) / (2.0 * h),
                    2e-5,
                    2e-6,
                );
                asymmetry =
                    asymmetry.max((response.tangent_pa[i][j] - response.tangent_pa[j][i]).abs());
            }
        }
        if core::ptr::eq(old, &before) {
            assert!(asymmetry > 1e-3);
        }
    }
}

#[test]
fn rigid_rotation_preserves_damage_metric_and_storage_with_signed_work_defect() {
    let m = material(0.002);
    let old = m.initial_state();
    let f = coaxial(0.04, -0.003);
    let r = m.trial(&f, &old, conditions()).unwrap();
    let angle = 0.7_f64;
    let (s, c) = angle.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let rotated_f = multiply(&rotation, &f);
    let same_trial = m.trial(&rotated_f, &old, conditions()).unwrap();
    let expected_p = multiply(&rotation, &r.first_piola_pa);
    let expected_tau = multiply(&multiply(&rotation, &r.kirchhoff_pa), &transpose(&rotation));
    for i in 0..9 {
        near(same_trial.first_piola_pa[i], expected_p[i], 2e-11, 2e-11);
        near(same_trial.kirchhoff_pa[i], expected_tau[i], 2e-11, 2e-11);
    }
    near(same_trial.state.damage(), r.state.damage(), 2e-12, 2e-12);
    for (a, b) in same_trial
        .state
        .native_state()
        .inverse_plastic_cauchy_green
        .iter()
        .zip(r.state.native_state().inverse_plastic_cauchy_green)
    {
        near(*a, b, 2e-12, 2e-12);
    }
    // A subsequent rigid rotation is lossless but endpoint Piola quadrature
    // can have negative work. It must remain a signed numerical account.
    let rotated = m.trial(&rotated_f, &r.state, conditions()).unwrap();
    assert_eq!(rotated.state.native_state(), r.state.native_state());
    assert_eq!(rotated.state.damage(), r.state.damage());
    assert_eq!(rotated.energy.plastic_dissipation_j_m3, 0.0);
    assert_eq!(rotated.energy.damage_release_j_m3, 0.0);
    near(
        rotated.state.stored_energy_j_m3(),
        r.state.stored_energy_j_m3(),
        2e-12,
        2e-12,
    );
    assert!(rotated.energy.endpoint_work_defect_j_m3 < 0.0);
    let work: f64 = (0..9)
        .map(|i| rotated.first_piola_pa[i] * (rotated_f[i] - f[i]))
        .sum();
    near(work, rotated.energy.endpoint_work_defect_j_m3, 2e-12, 2e-12);
}

#[test]
fn compression_keeps_bulk_strength_and_exact_holds_and_unloading_cannot_heal() {
    let m = material(0.002);
    let old = m.initial_state();
    let damaged = m
        .trial(&coaxial(0.04, 0.0), &old, conditions())
        .unwrap()
        .state;
    let hold = m
        .trial(damaged.deformation(), &damaged, conditions())
        .unwrap();
    assert_eq!(hold.state, damaged);
    assert_eq!(hold.energy.mechanical_work_j_m3, 0.0);
    assert_eq!(hold.energy.endpoint_work_defect_j_m3, 0.0);
    let compressed = coaxial(0.0398, -0.02);
    let original = m
        .native()
        .response(&compressed, damaged.native_state())
        .unwrap();
    let r = m.trial(&compressed, &damaged, conditions()).unwrap();
    assert_eq!(r.state.damage(), damaged.damage());
    let native_pressure =
        (original.kirchhoff_pa[0] + original.kirchhoff_pa[4] + original.kirchhoff_pa[8]) / 3.0;
    let actual_pressure = (r.kirchhoff_pa[0] + r.kirchhoff_pa[4] + r.kirchhoff_pa[8]) / 3.0;
    near(actual_pressure, -60.0, 3e-11, 3e-12);
    near(actual_pressure, native_pressure, 2e-11, 2e-12);
    for i in [0, 4, 8] {
        near(
            r.kirchhoff_pa[i] - actual_pressure,
            damaged.integrity() * (original.kirchhoff_pa[i] - native_pressure),
            2e-11,
            2e-12,
        );
    }
}

#[test]
fn physical_losses_close_independent_work_and_refinement_reduces_absolute_defect() {
    let m = material(0.002);
    let mut defects = Vec::new();
    for steps in [20, 80] {
        let mut state = m.initial_state();
        let mut sum = 0.0;
        for index in 1..=steps {
            let f = coaxial(
                0.04 * index as f64 / steps as f64,
                0.002 * index as f64 / steps as f64,
            );
            let r = m.trial(&f, &state, conditions()).unwrap();
            let work: f64 = (0..9)
                .map(|i| r.first_piola_pa[i] * (f[i] - state.deformation()[i]))
                .sum();
            near(
                work,
                r.state.stored_energy_j_m3() - state.stored_energy_j_m3()
                    + r.energy.plastic_dissipation_j_m3
                    + r.energy.damage_release_j_m3
                    + r.energy.endpoint_work_defect_j_m3,
                2e-12,
                2e-12,
            );
            sum += r.energy.endpoint_work_defect_j_m3.abs();
            state = r.state;
        }
        assert!(state.dissipated_plastic_j_m3() > 0.0);
        assert!(state.dissipated_damage_j_m3() > 0.0);
        defects.push(sum);
    }
    assert!(defects[1] < 0.35 * defects[0], "{defects:?}");
}

#[test]
fn owner_damage_domain_and_geometry_refusals_preserve_history() {
    let m = material(0.002);
    let other = m.clone();
    let old = m.initial_state();
    let retained = old.clone();
    assert!(matches!(
        other.trial(&coaxial(0.04, 0.0), &old, conditions()),
        Err(FiniteDuctileError::StateMismatch)
    ));
    assert!(
        m.trial(
            &[-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            &old,
            conditions()
        )
        .is_err()
    );
    assert!(
        m.trial(
            old.deformation(),
            &old,
            DuctileJ2Conditions {
                temperature_k: 301.0,
                ..conditions()
            }
        )
        .is_err()
    );
    let limited = FiniteDuctileJ2::new(
        native(),
        m.initiation_law().clone(),
        DuctileJ2Settings {
            maximum_damage_increment: 0.001,
            ..settings()
        },
    )
    .unwrap();
    let virgin = limited.initial_state();
    assert!(matches!(
        limited.trial(&coaxial(0.04, 0.0), &virgin, conditions()),
        Err(FiniteDuctileError::DamageIncrement { .. })
    ));
    assert!(matches!(
        limited.stress(&coaxial(0.04, 0.0), &virgin, conditions()),
        Err(FiniteDuctileError::DamageIncrement { .. })
    ));
    assert_eq!(old, retained);
    assert_eq!(virgin, limited.initial_state());
}
