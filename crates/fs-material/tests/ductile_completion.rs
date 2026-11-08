//! Independent completion and work oracles for the opt-in finite-support law.
//! Synthetic coefficients isolate constitutive behavior; they are not an alloy
//! calibration, a fracture-energy identification, or a topology-deletion rule.
#![deny(unsafe_code)]

use fs_material::IsotropicElastic;
use fs_material::ductile::evolution::{
    DuctileEvolutionError, DuctileJ2, DuctileJ2Conditions, DuctileJ2Settings, DuctileWorkEvolution,
};
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::{FiniteDuctileError, FiniteDuctileJ2};
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::plastic::reduction::{
    PlaneStressReduction, ReducedState, ReductionConfig, ReductionError,
};

const MU: f64 = 512.0;
const FINITE_K: f64 = 1024.0;
const SMALL_K: f64 = 1280.0 / 1.5;
const Y: f64 = 3.0;
const WIDTH: f64 = 1.0 / 32.0;
const GP: f64 = 3.0 / 32768.0;
const ONSET: f64 = 1.0 / 1024.0;

fn locus() -> TabulatedDuctileLocus {
    TabulatedDuctileLocus::new(
        vec![(-50.0, ONSET), (50.0, ONSET)],
        DuctileDomain {
            triaxiality: [-50.0, 50.0],
            strain_rate_s_inv: [1e-12, 1000.0],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap()
}

fn settings() -> DuctileJ2Settings {
    DuctileJ2Settings {
        characteristic_length_m: WIDTH,
        post_initiation_plastic_work_j_m2: GP,
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

fn small_default() -> DuctileJ2<TabulatedDuctileLocus> {
    // These binary-representable E and nu give mu=512 exactly.
    DuctileJ2::new(
        J2Plasticity::new(IsotropicElastic::new(1280.0, 0.25, 0.05).unwrap(), Y, 0.0).unwrap(),
        locus(),
        settings(),
    )
    .unwrap()
}

fn finite_default() -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    FiniteDuctileJ2::new(
        FiniteJ2::new(
            FiniteJ2Parameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: FINITE_K,
                yield_stress_pa: Y,
                hardening_modulus_pa: 0.0,
            },
            FiniteJ2Limits {
                minimum_jacobian: 0.1,
                maximum_elastic_log_strain: 0.5,
                maximum_plastic_log_strain: 2.0,
                maximum_accumulated_plastic_strain: 3.0,
            },
        )
        .unwrap(),
        locus(),
        settings(),
    )
    .unwrap()
}

fn small() -> DuctileJ2<TabulatedDuctileLocus> {
    small_default().with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic)
}

fn finite() -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    finite_default().with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic)
}

fn strain(e: f64, dilation: f64) -> [f64; 6] {
    [
        dilation + e,
        dilation - 0.5 * e,
        dilation - 0.5 * e,
        0.0,
        0.0,
        0.0,
    ]
}

fn deformation(e: f64, dilation: f64) -> [f64; 9] {
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

// For monotone perfect flow, alpha=e-Y/(3*mu), hence the independently
// integrated coordinate is x=WIDTH*Y*(alpha-ONSET)/GP=1024*alpha-1.
fn e_at(x: f64) -> f64 {
    (x + 3.0) / 1024.0
}

fn integrity_at(x: f64) -> f64 {
    ((3.0 - x.max(0.0)).max(0.0) / 3.0).powi(2)
}

fn capacity_at(x: f64) -> f64 {
    // Integral of 1-2*x/3+x*x/9, evaluated without onset cancellation.
    let x = x.clamp(0.0, 3.0);
    x * (1.0 - x / 3.0 + x * x / 27.0)
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}; difference {:.8e}",
        (actual - expected).abs()
    );
}

fn tensor_work(stress: &[f64; 6], delta: &[f64; 6]) -> f64 {
    (0..6)
        .map(|i| stress[i] * delta[i] * if i < 3 { 1.0 } else { 2.0 })
        .sum()
}

#[test]
fn both_native_owners_integrate_the_chosen_work_capacity_and_continue_the_original_clock() {
    let s = small();
    let f = finite();
    assert_eq!(s.evolution(), DuctileWorkEvolution::FiniteSupportQuadratic);
    assert_eq!(f.evolution(), DuctileWorkEvolution::FiniteSupportQuadratic);
    let mut ss = s.initial_state();
    let mut fs = f.initial_state();
    let mut small_integral = 0.0;
    let mut finite_integral = 0.0;
    let dilation = 0.0002;
    for x in [-0.25_f64, 0.25, 1.0, 2.0, 2.8, 3.25, 4.0] {
        let sr = s
            .trial(&strain(e_at(x), dilation), &ss, conditions())
            .unwrap();
        let fr = f
            .trial(&deformation(e_at(x), dilation), &fs, conditions())
            .unwrap();
        let g = integrity_at(x);
        let area_work = GP * capacity_at(x);
        near(
            sr.state.native_state().alpha,
            (x + 1.0) / 1024.0,
            2e-14,
            2e-12,
        );
        near(
            fr.state.native_state().alpha,
            (x + 1.0) / 1024.0,
            2e-13,
            2e-11,
        );
        near(sr.state.integrity(), g, 2e-12, 2e-12);
        near(fr.state.integrity(), g, 2e-11, 2e-11);
        near(
            sr.state.equivalent_plastic_work_j_m2(),
            GP * x.max(0.0),
            2e-15,
            2e-11,
        );
        near(
            fr.state.equivalent_plastic_work_j_m2(),
            GP * x.max(0.0),
            2e-14,
            2e-10,
        );
        small_integral += sr.energy.degraded_equivalent_plastic_work_increment_j_m2;
        finite_integral += fr.energy.degraded_equivalent_plastic_work_increment_j_m2;
        near(small_integral, area_work, 2e-15, 2e-11);
        near(finite_integral, area_work, 2e-14, 2e-10);
        near(
            sr.state.degraded_equivalent_plastic_work_j_m2(),
            area_work,
            2e-15,
            2e-11,
        );
        near(
            fr.state.degraded_equivalent_plastic_work_j_m2(),
            area_work,
            2e-14,
            2e-10,
        );
        for i in 0..3 {
            let deviator = if i == 0 { 2.0 * Y / 3.0 } else { -Y / 3.0 };
            near(
                sr.stress_pa[i],
                g * (3.0 * SMALL_K * dilation + deviator),
                2e-10,
                2e-10,
            );
            near(
                fr.kirchhoff_pa[4 * i],
                g * (3.0 * FINITE_K * dilation + deviator),
                2e-10,
                2e-10,
            );
        }
        assert_eq!(sr.state.is_fully_failed(), x > 3.0);
        assert_eq!(fr.state.is_fully_failed(), x > 3.0);
        assert!(sr.state.integrity() <= ss.integrity());
        assert!(fr.state.integrity() <= fs.integrity());
        if x == 1.0 {
            near(sr.state.damage(), 5.0 / 9.0, 2e-12, 2e-12);
            near(fr.state.damage(), 5.0 / 9.0, 2e-11, 2e-11);
            // Capacity is Phi=19/27 here, not damage=5/9.
            near(small_integral / GP, 19.0 / 27.0, 2e-11, 2e-11);
            near(finite_integral / GP, 19.0 / 27.0, 2e-10, 2e-10);
        }
        if ss.is_fully_failed() {
            assert_eq!(
                sr.energy.degraded_equivalent_plastic_work_increment_j_m2,
                0.0
            );
            assert!(sr.energy.equivalent_plastic_work_increment_j_m2 > 0.0);
        }
        if fs.is_fully_failed() {
            assert_eq!(
                fr.energy.degraded_equivalent_plastic_work_increment_j_m2,
                0.0
            );
            assert!(fr.energy.equivalent_plastic_work_increment_j_m2 > 0.0);
        }
        ss = sr.state;
        fs = fr.state;
    }
    near(small_integral, GP, 2e-15, 2e-11);
    near(finite_integral, GP, 2e-14, 2e-10);
}

#[test]
fn tiny_remaining_capacity_is_not_lost_to_subtraction_of_rounded_cumulative_receipts() {
    let s = small();
    let f = finite();
    let x0 = 3.0 - 2e-6;
    let x1 = 3.0 - 1e-6;
    let s0 = s
        .trial(&strain(e_at(x0), 0.0002), &s.initial_state(), conditions())
        .unwrap();
    let f0 = f
        .trial(
            &deformation(e_at(x0), 0.0002),
            &f.initial_state(),
            conditions(),
        )
        .unwrap();
    let s1 = s
        .trial(&strain(e_at(x1), 0.0002), &s0.state, conditions())
        .unwrap();
    let f1 = f
        .trial(&deformation(e_at(x1), 0.0002), &f0.state, conditions())
        .unwrap();
    for (before, after, increment) in [
        (
            s0.state.equivalent_plastic_work_j_m2() / GP,
            s1.state.equivalent_plastic_work_j_m2() / GP,
            s1.energy.degraded_equivalent_plastic_work_increment_j_m2,
        ),
        (
            f0.state.equivalent_plastic_work_j_m2() / GP,
            f1.state.equivalent_plastic_work_j_m2() / GP,
            f1.energy.degraded_equivalent_plastic_work_increment_j_m2,
        ),
    ] {
        let r0 = (3.0 - before) / 3.0;
        let r1 = (3.0 - after) / 3.0;
        let expected = GP * (r0 - r1) * (r0 * r0 + r0 * r1 + r1 * r1);
        assert!(expected > 0.0 && increment > 0.0);
        near(increment, expected, 1e-30, 2e-5);
    }
    assert!(!s1.state.is_fully_failed());
    assert!(!f1.state.is_fully_failed());
}

#[test]
fn complete_small_and_finite_tangents_include_quadratic_softening_in_tension_and_compression() {
    let s = small();
    let f = finite();
    let so = s
        .trial(&strain(0.001, 0.0001), &s.initial_state(), conditions())
        .unwrap()
        .state;
    let fo = f
        .trial(
            &deformation(0.001, 0.0001),
            &f.initial_state(),
            conditions(),
        )
        .unwrap()
        .state;
    let saved_small = so.clone();
    let saved_finite = fo.clone();
    for compress in [false, true] {
        let mut eps = [0.0045, -0.0014, -0.0016, 0.0009, -0.0003, 0.0006];
        let mut endpoint = [
            1.0048, 0.0018, -0.0006, 0.0002, 0.9989, 0.0009, 0.0004, -0.0003, 0.9986,
        ];
        if compress {
            for value in &mut eps[..3] {
                *value -= 0.004;
            }
            for value in &mut endpoint {
                *value *= 0.99;
            }
        }
        let sr = s.trial(&eps, &so, conditions()).unwrap();
        let fr = f.trial(&endpoint, &fo, conditions()).unwrap();
        assert!(sr.state.integrity() > 0.0 && sr.state.integrity() < 1.0);
        assert!(fr.state.integrity() > 0.0 && fr.state.integrity() < 1.0);
        for j in 0..6 {
            let h = 1e-8;
            let mut plus = eps;
            let mut minus = eps;
            plus[j] += h;
            minus[j] -= h;
            let p = s.trial(&plus, &so, conditions()).unwrap();
            let m = s.trial(&minus, &so, conditions()).unwrap();
            for i in 0..6 {
                near(
                    sr.branch_tangent_pa[i][j],
                    (p.stress_pa[i] - m.stress_pa[i]) / (2.0 * h),
                    2e-5,
                    3e-6,
                );
            }
        }
        for j in 0..9 {
            let h = 1e-7;
            let mut plus = endpoint;
            let mut minus = endpoint;
            plus[j] += h;
            minus[j] -= h;
            let p = f.stress(&plus, &fo, conditions()).unwrap();
            let m = f.stress(&minus, &fo, conditions()).unwrap();
            for i in 0..9 {
                near(fr.tangent_pa[i][j], (p[i] - m[i]) / (2.0 * h), 3e-5, 3e-6);
            }
        }
    }
    assert_eq!(so, saved_small);
    assert_eq!(fo, saved_finite);
}

#[test]
fn terminal_release_closes_physical_work_without_counting_the_area_diagnostic_twice() {
    let s = small();
    let f = finite();
    let so = s
        .trial(&strain(e_at(1.0), 0.001), &s.initial_state(), conditions())
        .unwrap()
        .state;
    let fo = f
        .trial(
            &deformation(e_at(1.0), 0.001),
            &f.initial_state(),
            conditions(),
        )
        .unwrap()
        .state;
    assert!(so.stored_energy_j_m3() > 0.0 && fo.stored_energy_j_m3() > 0.0);
    let eps = strain(e_at(3.25), 0.001);
    let endpoint = deformation(e_at(3.25), 0.001);
    let sr = s.trial(&eps, &so, conditions()).unwrap();
    let fr = f.trial(&endpoint, &fo, conditions()).unwrap();
    assert!(sr.state.is_fully_failed() && fr.state.is_fully_failed());
    assert_eq!(sr.stress_pa, [0.0; 6]);
    assert_eq!(sr.branch_tangent_pa, [[0.0; 6]; 6]);
    assert_eq!(fr.first_piola_pa, [0.0; 9]);
    assert_eq!(fr.tangent_pa, [[0.0; 9]; 9]);
    assert_eq!(sr.state.stored_energy_j_m3(), 0.0);
    assert_eq!(fr.state.stored_energy_j_m3(), 0.0);
    assert_eq!(sr.energy.plastic_dissipation_j_m3, 0.0);
    assert_eq!(fr.energy.plastic_dissipation_j_m3, 0.0);
    assert_eq!(sr.energy.numerical_dissipation_j_m3, 0.0);
    assert_eq!(fr.energy.endpoint_work_defect_j_m3, 0.0);
    near(
        sr.energy.damage_release_j_m3,
        so.stored_energy_j_m3(),
        2e-12,
        2e-12,
    );
    near(
        fr.energy.damage_release_j_m3,
        fo.stored_energy_j_m3(),
        2e-12,
        2e-12,
    );
    let ds = core::array::from_fn(|i| eps[i] - so.strain()[i]);
    let small_work = tensor_work(&sr.stress_pa, &ds);
    let finite_work: f64 = (0..9)
        .map(|i| fr.first_piola_pa[i] * (endpoint[i] - fo.deformation()[i]))
        .sum();
    near(
        small_work,
        sr.energy.final_stored_j_m3 - sr.energy.initial_stored_j_m3
            + sr.energy.plastic_dissipation_j_m3
            + sr.energy.damage_release_j_m3
            + sr.energy.numerical_dissipation_j_m3
            + sr.energy.residual_j_m3,
        2e-12,
        2e-12,
    );
    near(
        finite_work,
        fr.energy.final_stored_j_m3 - fr.energy.initial_stored_j_m3
            + fr.energy.plastic_dissipation_j_m3
            + fr.energy.damage_release_j_m3
            + fr.energy.endpoint_work_defect_j_m3
            + fr.energy.residual_j_m3,
        2e-12,
        2e-12,
    );
    assert!(sr.energy.degraded_equivalent_plastic_work_increment_j_m2 > 0.0);
    assert!(fr.energy.degraded_equivalent_plastic_work_increment_j_m2 > 0.0);
    // A one-step jump from zero storage has no endpoint mechanical work or
    // damage release. Exhausting Phi therefore cannot itself mint a heat or
    // fracture-energy credit; physical work remains the separately audited path.
    let sj = s.trial(&eps, &s.initial_state(), conditions()).unwrap();
    let fj = f
        .trial(&endpoint, &f.initial_state(), conditions())
        .unwrap();
    assert_eq!(sj.energy.damage_release_j_m3, 0.0);
    assert_eq!(fj.energy.damage_release_j_m3, 0.0);
    assert_eq!(sj.energy.mechanical_work_j_m3, 0.0);
    assert_eq!(fj.energy.mechanical_work_j_m3, 0.0);
    near(
        sj.state.degraded_equivalent_plastic_work_j_m2(),
        GP,
        2e-15,
        2e-12,
    );
    near(
        fj.state.degraded_equivalent_plastic_work_j_m2(),
        GP,
        2e-15,
        2e-12,
    );
}

#[test]
fn failed_material_retains_exact_compressive_bulk_and_cannot_heal_on_unloading_or_holds() {
    let s = small();
    let f = finite();
    let ss = s
        .trial(&strain(e_at(3.25), 0.001), &s.initial_state(), conditions())
        .unwrap()
        .state;
    let fs = f
        .trial(
            &deformation(e_at(3.25), 0.001),
            &f.initial_state(),
            conditions(),
        )
        .unwrap()
        .state;
    let original_small = ss.clone();
    let original_finite = fs.clone();
    let eps = [-0.004, -0.004, -0.004, 0.0, 0.0, 0.0];
    let a = 0.98_f64;
    let compressed = [a, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, a];
    let sr = s.trial(&eps, &ss, conditions()).unwrap();
    let fr = f.trial(&compressed, &fs, conditions()).unwrap();
    let small_pressure = SMALL_K * -0.012;
    let log_j = 3.0 * a.ln();
    let finite_pressure = FINITE_K * log_j;
    assert!(sr.state.is_fully_failed() && fr.state.is_fully_failed());
    assert_eq!(sr.state.integrity(), 0.0);
    assert_eq!(fr.state.integrity(), 0.0);
    near(
        sr.state.stored_energy_j_m3(),
        0.5 * SMALL_K * 0.012_f64.powi(2),
        2e-12,
        2e-12,
    );
    near(
        fr.state.stored_energy_j_m3(),
        0.5 * FINITE_K * log_j * log_j,
        2e-11,
        2e-11,
    );
    near(
        fr.state.compression_stored_energy_j_m3(),
        fr.state.stored_energy_j_m3(),
        2e-12,
        2e-12,
    );
    for i in 0..6 {
        near(
            sr.stress_pa[i],
            if i < 3 { small_pressure } else { 0.0 },
            2e-12,
            2e-12,
        );
        for j in 0..6 {
            near(
                sr.branch_tangent_pa[i][j],
                if i < 3 && j < 3 { SMALL_K } else { 0.0 },
                2e-11,
                2e-12,
            );
        }
    }
    // Independent d[K*ln(J)*F^(-T)]/dF, all 81 entries. A finite
    // compression tangent includes the geometric -p term as well as K.
    for ij in 0..9 {
        let (i, j) = (ij / 3, ij % 3);
        near(
            fr.first_piola_pa[ij],
            if i == j { finite_pressure / a } else { 0.0 },
            2e-10,
            2e-11,
        );
        for kl in 0..9 {
            let (k, l) = (kl / 3, kl % 3);
            let bulk = if i == j && k == l {
                FINITE_K / (a * a)
            } else {
                0.0
            };
            let geometric = if i == l && k == j {
                -finite_pressure / (a * a)
            } else {
                0.0
            };
            near(fr.tangent_pa[ij][kl], bulk + geometric, 2e-9, 2e-11);
        }
    }
    assert!(sr.state.equivalent_plastic_work_j_m2() >= ss.equivalent_plastic_work_j_m2());
    assert!(fr.state.equivalent_plastic_work_j_m2() >= fs.equivalent_plastic_work_j_m2());
    assert_eq!(
        sr.energy.degraded_equivalent_plastic_work_increment_j_m2,
        0.0
    );
    assert_eq!(
        fr.energy.degraded_equivalent_plastic_work_increment_j_m2,
        0.0
    );
    let sh = s.trial(sr.state.strain(), &sr.state, conditions()).unwrap();
    let fh = f
        .trial(fr.state.deformation(), &fr.state, conditions())
        .unwrap();
    assert_eq!(sh.state, sr.state);
    assert_eq!(fh.state, fr.state);
    assert_eq!(sh.energy.mechanical_work_j_m3, 0.0);
    assert_eq!(fh.energy.mechanical_work_j_m3, 0.0);
    let su = s
        .trial(&strain(0.0, 0.001), &sr.state, conditions())
        .unwrap();
    let fu = f
        .trial(&deformation(0.0, 0.001), &fr.state, conditions())
        .unwrap();
    assert!(su.state.is_fully_failed() && fu.state.is_fully_failed());
    assert_eq!(su.stress_pa, [0.0; 6]);
    assert_eq!(fu.first_piola_pa, [0.0; 9]);
    assert_eq!(su.state.stored_energy_j_m3(), 0.0);
    assert_eq!(fu.state.stored_energy_j_m3(), 0.0);
    assert_eq!(ss, original_small);
    assert_eq!(fs, original_finite);
}

#[test]
fn completion_uses_profile_support_and_never_rounded_damage_or_a_tiny_integrity_cutoff() {
    let s = small();
    let f = finite();
    let x = 3.0 - 1e-8;
    let sr = s
        .trial(&strain(e_at(x), 0.0002), &s.initial_state(), conditions())
        .unwrap();
    let fr = f
        .trial(
            &deformation(e_at(x), 0.0002),
            &f.initial_state(),
            conditions(),
        )
        .unwrap();
    assert!(sr.state.integrity() > 0.0 && sr.state.integrity() < 1e-15);
    assert!(fr.state.integrity() > 0.0 && fr.state.integrity() < 1e-15);
    assert!(!sr.state.is_fully_failed());
    assert!(!fr.state.is_fully_failed());
    let es = small_default();
    let ef = finite_default();
    // x=38 has strictly positive exp(-x), but 1-exp(-x) rounds to one.
    // The small-strain Frobenius norm remains below its declared 0.05 bound.
    let se = es
        .trial(&strain(e_at(38.0), 0.0), &es.initial_state(), conditions())
        .unwrap();
    let fe = ef
        .trial(
            &deformation(e_at(38.0), 0.0),
            &ef.initial_state(),
            conditions(),
        )
        .unwrap();
    assert_eq!(se.state.damage(), 1.0);
    assert_eq!(fe.state.damage(), 1.0);
    assert!(se.state.integrity() > 0.0);
    assert!(fe.state.integrity() > 0.0);
    assert!(!se.state.is_fully_failed());
    assert!(!fe.state.is_fully_failed());
    near(se.state.integrity(), (-38.0_f64).exp(), 1e-28, 2e-10);
    near(fe.state.integrity(), (-38.0_f64).exp(), 1e-28, 2e-9);
}

#[test]
fn default_and_explicit_exponential_preserve_stress_history_tangent_and_work_receipts() {
    let s = small_default();
    let explicit_s = small_default().with_evolution(DuctileWorkEvolution::Exponential);
    let f = finite_default();
    let explicit_f = finite_default().with_evolution(DuctileWorkEvolution::Exponential);
    assert_eq!(s.evolution(), DuctileWorkEvolution::Exponential);
    assert_eq!(f.evolution(), DuctileWorkEvolution::Exponential);
    // Histories are born from their own immutable material owner.
    let mut so = s.initial_state();
    let mut se = explicit_s.initial_state();
    let mut fo = f.initial_state();
    let mut fe = explicit_f.initial_state();
    for (e, dilation) in [
        (0.001, 0.0001),
        (e_at(0.25), 0.0002),
        (e_at(1.0), 0.0002),
        (e_at(1.0), 0.0002),
        (e_at(1.0) - 0.0001, -0.001),
        (e_at(2.0), -0.001),
    ] {
        let sr = s.trial(&strain(e, dilation), &so, conditions()).unwrap();
        let er = explicit_s
            .trial(&strain(e, dilation), &se, conditions())
            .unwrap();
        assert_eq!(sr.stress_pa, er.stress_pa);
        assert_eq!(sr.branch_tangent_pa, er.branch_tangent_pa);
        assert_eq!(sr.energy, er.energy);
        assert_eq!(sr.state.native_state(), er.state.native_state());
        assert_eq!(sr.state.initiation(), er.state.initiation());
        assert_eq!(sr.state.integrity(), er.state.integrity());
        assert_eq!(sr.state.damage(), er.state.damage());
        let fr = f
            .trial(&deformation(e, dilation), &fo, conditions())
            .unwrap();
        let ef = explicit_f
            .trial(&deformation(e, dilation), &fe, conditions())
            .unwrap();
        assert_eq!(fr.first_piola_pa, ef.first_piola_pa);
        assert_eq!(fr.kirchhoff_pa, ef.kirchhoff_pa);
        assert_eq!(fr.tangent_pa, ef.tangent_pa);
        assert_eq!(fr.energy, ef.energy);
        assert_eq!(fr.state.native_state(), ef.state.native_state());
        assert_eq!(fr.state.initiation(), ef.state.initiation());
        assert_eq!(fr.state.integrity(), ef.state.integrity());
        assert_eq!(fr.state.damage(), ef.state.damage());
        let sx = sr.state.equivalent_plastic_work_j_m2() / GP;
        let fx = fr.state.equivalent_plastic_work_j_m2() / GP;
        near(sr.state.integrity(), (-sx).exp(), 2e-15, 2e-13);
        near(fr.state.integrity(), (-fx).exp(), 2e-15, 2e-13);
        near(
            sr.state.degraded_equivalent_plastic_work_j_m2(),
            GP * -(-sx).exp_m1(),
            2e-16,
            2e-13,
        );
        near(
            fr.state.degraded_equivalent_plastic_work_j_m2(),
            GP * -(-fx).exp_m1(),
            2e-16,
            2e-13,
        );
        so = sr.state;
        se = er.state;
        fo = fr.state;
        fe = ef.state;
    }
}

#[test]
fn plane_stress_refuses_terminal_loss_of_closure_and_keeps_the_accepted_history() {
    let s = small();
    let f = finite();
    let ss = s
        .trial(&strain(e_at(3.25), 0.001), &s.initial_state(), conditions())
        .unwrap()
        .state;
    let fs = f
        .trial(
            &deformation(e_at(3.25), 0.001),
            &f.initial_state(),
            conditions(),
        )
        .unwrap()
        .state;
    let saved_small = ss.clone();
    let saved_finite = fs.clone();
    let reduced = PlaneStressReduction::new(
        s.small_strain_law(conditions()).unwrap(),
        ReductionConfig {
            max_iterations: 32,
            stress_tolerance_pa: 1e-10,
            strain_limit: 0.05,
        },
    )
    .unwrap();
    let old = ReducedState {
        strain: *ss.strain(),
        history: ss,
    };
    // The stress already vanishes, but an underdetermined transverse strain
    // is not an admissible condensed response; the factorization must refuse.
    assert!(matches!(
        reduced.evaluate([old.strain[0], old.strain[1], old.strain[3]], &old),
        Err(ReductionError::Linear(
            fs_la::LuWorkspaceError::Singular { .. }
        ))
    ));
    assert!(matches!(
        f.plane_stress_trial(
            &[1.01, 0.0, 0.0, 1.01],
            &fs,
            [0.0, 0.0, 1.0],
            conditions(),
            PlaneStressSettings::default(),
        ),
        Err(FiniteDuctileError::PlaneStressTerminalFailure)
    ));
    assert_eq!(old.history, saved_small);
    assert_eq!(fs, saved_finite);
}

#[test]
fn damage_increment_gate_uses_loss_of_integrity_and_rejections_remain_transactional() {
    for limit in [0.5, 0.6] {
        let base_s = small_default();
        let s = DuctileJ2::new(
            base_s.native().clone(),
            locus(),
            DuctileJ2Settings {
                maximum_damage_increment: limit,
                ..settings()
            },
        )
        .unwrap()
        .with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic);
        let base_f = finite_default();
        let f = FiniteDuctileJ2::new(
            base_f.native().clone(),
            locus(),
            DuctileJ2Settings {
                maximum_damage_increment: limit,
                ..settings()
            },
        )
        .unwrap()
        .with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic);
        let so = s.initial_state();
        let fo = f.initial_state();
        let ss = so.clone();
        let fs = fo.clone();
        let sr = s.trial(&strain(e_at(1.0), 0.0002), &so, conditions());
        let fr = f.trial(&deformation(e_at(1.0), 0.0002), &fo, conditions());
        if limit < 5.0 / 9.0 {
            assert!(matches!(
                sr,
                Err(DuctileEvolutionError::DamageIncrement { .. })
            ));
            assert!(matches!(
                fr,
                Err(FiniteDuctileError::DamageIncrement { .. })
            ));
            assert!(matches!(
                f.stress(&deformation(e_at(1.0), 0.0002), &fo, conditions()),
                Err(FiniteDuctileError::DamageIncrement { .. })
            ));
        } else {
            // Delta damage=5/9 fits; Delta Phi=19/27 does not. These gates
            // cannot accidentally use the newly distinct capacity increment.
            near(sr.unwrap().state.damage(), 5.0 / 9.0, 2e-12, 2e-12);
            near(fr.unwrap().state.damage(), 5.0 / 9.0, 2e-11, 2e-11);
        }
        assert_eq!(so, ss);
        assert_eq!(fo, fs);
    }
}

#[test]
fn physical_overstress_hold_crosses_terminal_support_without_recalibrating_capacity_or_time() {
    let viscosity = 256.0;
    let dt = conditions().duration_s;
    let m = finite()
        .with_overstress(FiniteOverstressParameters {
            viscosity_pa_s: viscosity,
            stress_scale_pa: Y,
            exponent: 1.0,
        })
        .unwrap();
    let alpha0 = (2.5 + 1.0) / 1024.0;
    let e = Y / (3.0 * MU) + alpha0 * (1.0 + viscosity / (3.0 * MU * dt));
    let endpoint = deformation(e, 0.0002);
    let old = m.initial_state();
    let loaded = m.trial(&endpoint, &old, conditions()).unwrap();
    let q0 = Y + viscosity * alpha0 / dt;
    near(loaded.state.native_state().alpha, alpha0, 3e-13, 2e-10);
    near(loaded.state.integrity(), integrity_at(2.5), 3e-11, 2e-10);
    near(
        loaded.state.equivalent_plastic_work_j_m2(),
        WIDTH * Y * (alpha0 - ONSET),
        2e-14,
        2e-10,
    );
    near(
        loaded.energy.plastic_dissipation_j_m3,
        loaded.state.integrity() * q0 * alpha0,
        2e-11,
        2e-10,
    );
    assert!(!loaded.state.is_fully_failed());
    let saved = loaded.state.clone();
    // Independent backward-Euler relaxation at fixed total logarithmic strain.
    // One trial consumes exactly one dt from the original accepted metric.
    let delta = (3.0 * MU * (e - alpha0) - Y) / (3.0 * MU + viscosity / dt);
    let alpha1 = alpha0 + delta;
    let held = m.trial(&endpoint, &loaded.state, conditions()).unwrap();
    near(held.state.native_state().alpha, alpha1, 5e-13, 3e-10);
    assert!(1024.0 * alpha1 - 1.0 > 3.0);
    assert!(held.state.is_fully_failed());
    assert_eq!(held.first_piola_pa, [0.0; 9]);
    assert_eq!(held.energy.mechanical_work_j_m3, 0.0);
    assert_eq!(held.energy.plastic_dissipation_j_m3, 0.0);
    near(
        held.energy.damage_release_j_m3,
        loaded.state.stored_energy_j_m3(),
        2e-12,
        2e-11,
    );
    near(
        held.state.degraded_equivalent_plastic_work_j_m2(),
        GP,
        2e-15,
        2e-12,
    );
    let later = m.trial(&endpoint, &held.state, conditions()).unwrap();
    assert!(later.state.native_state().alpha > held.state.native_state().alpha);
    assert!(later.state.is_fully_failed());
    assert_eq!(
        later.energy.degraded_equivalent_plastic_work_increment_j_m2,
        0.0
    );
    assert_eq!(later.energy.plastic_dissipation_j_m3, 0.0);
    assert_eq!(later.energy.damage_release_j_m3, 0.0);
    assert_eq!(loaded.state, saved);
}
