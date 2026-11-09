//! Native mixed creep/overstress and independent differential/recovery oracles.
//! Synthetic coefficients test the implementation, not a calibrated specimen.

use fs_material::creep::StrainHardeningCreep;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::creep::recovery::{
    FiniteRecoveringNortonCreep, FiniteRecoveringNortonJ2,
};
use fs_material::plastic::finite::mixed::MixedFiniteInelasticMaterial;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const DT: f64 = 0.2;
const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const F: [f64; 9] = [1.08, 0.04, 0.01, -0.02, 0.96, 0.03, 0.01, 0.0, 1.02];

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.05,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 2.0,
        maximum_accumulated_plastic_strain: 0.8,
    }
}

fn backbone(kappa: f64) -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: kappa,
            yield_stress_pa: 2e6,
            hardening_modulus_pa: 4e6,
        },
        limits(),
    )
    .unwrap()
}

fn kinetics() -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: 40e6,
        stress_scale_pa: 2e6,
        exponent: 1.6,
    }
}

fn primary() -> StrainHardeningCreep {
    StrainHardeningCreep {
        reference_strain: 0.015,
        exponent: 0.8,
        residual_rate_fraction: 0.2,
    }
}

fn perzyna(kappa: f64) -> FiniteOverstressJ2 {
    FiniteOverstressJ2::new(
        backbone(kappa).with_voce_hardening(6e6, 12.0).unwrap(),
        kinetics(),
        DT,
    )
    .unwrap()
}

fn norton(kappa: f64) -> FiniteNortonCreep {
    FiniteNortonCreep::new(MU, kappa, limits(), kinetics(), DT).unwrap()
}

fn combined(kappa: f64) -> FiniteNortonJ2 {
    FiniteNortonJ2::new(backbone(kappa), kinetics(), DT, 0.8).unwrap()
}

fn recovering_norton(kappa: f64) -> FiniteRecoveringNortonCreep {
    FiniteRecoveringNortonCreep::new(norton(kappa).with_strain_hardening(primary()).unwrap(), 1.5)
        .unwrap()
}

fn recovering_combined(kappa: f64) -> FiniteRecoveringNortonJ2 {
    FiniteRecoveringNortonJ2::new(
        combined(kappa).with_strain_hardening(primary()).unwrap(),
        1.5,
    )
    .unwrap()
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= scale * tolerance,
        "{a:e} versus {b:e}, scale {scale:e}, tolerance {tolerance:e}"
    );
}

fn original_pressure_parity<L: MixedFiniteInelasticMaterial>(law: &L)
where
    L::State: PartialEq,
{
    let old = law.initial_state();
    let before = old.clone();
    let native = law.response(&F, &old).unwrap();
    let remainder = law.remainder_response(&F, &old).unwrap();
    let j = determinant(&F);
    let mixed = law
        .freeze_remainder(&old)
        .unwrap()
        .mixed_response(&F, K * j.ln() / j, VolumetricTreatment::OriginalMaterial)
        .unwrap();
    assert_eq!(law.duration_s(), Some(DT));
    assert_eq!(old, before);
    assert_eq!(native.state, remainder.state);
    assert_eq!(native.plastic_increment, remainder.plastic_increment);
    assert_eq!(native.creep_increment, remainder.creep_increment);
    assert_eq!(native.flow_potential_j_m3, remainder.flow_potential_j_m3);
    assert_eq!(
        native.inelastic_dissipation_j_m3(),
        remainder.inelastic_dissipation_j_m3()
    );
    assert!(native.inelastic_dissipation_j_m3() > native.flow_potential_j_m3);
    near(mixed.state.pressure_residual, 0.0, 1.0, 2e-14);
    near(
        remainder.stored_energy_j_m3 + mixed.state.bulk_stored_pa,
        native.stored_energy_j_m3,
        MU,
        2e-12,
    );
    near(
        mixed.state.potential_pa,
        native.incremental_potential_j_m3,
        MU,
        2e-12,
    );
    law.validate_remainder_committed(&F, &remainder.state, 1e-4)
        .unwrap();
    near(
        law.stored_remainder_energy(&F, &remainder.state).unwrap(),
        remainder.stored_energy_j_m3,
        MU,
        2e-12,
    );
    for row in 0..9 {
        near(
            mixed.state.first_piola_pa[row],
            native.first_piola_pa[row],
            MU,
            2e-12,
        );
        for col in 0..9 {
            let condensed = mixed.tangent_pa[row][col]
                - mixed.pressure_coupling[row] * mixed.pressure_coupling[col]
                    / mixed.pressure_tangent_per_pa;
            near(condensed, native.tangent_pa[row][col], K, 3e-11);
        }
    }
}

#[test]
fn g1_all_original_kinetic_owners_recover_full_stress_energy_history_and_tangent() {
    original_pressure_parity(&perzyna(K));
    original_pressure_parity(&norton(K));
    original_pressure_parity(&norton(K).with_strain_hardening(primary()).unwrap());
    original_pressure_parity(&combined(K));
    original_pressure_parity(&combined(K).with_strain_hardening(primary()).unwrap());
    original_pressure_parity(&recovering_norton(K));
    original_pressure_parity(&recovering_combined(K));
    let table = TabulatedHardening::new(vec![
        HardeningPoint {
            plastic_strain: 0.0,
            yield_stress_pa: 2e6,
        },
        HardeningPoint {
            plastic_strain: 0.2,
            yield_stress_pa: 5e6,
        },
        HardeningPoint {
            plastic_strain: 0.8,
            yield_stress_pa: 9e6,
        },
    ])
    .unwrap();
    let mut parameters = *backbone(K).parameters();
    parameters.hardening_modulus_pa = 0.0;
    let tabulated = FiniteJ2::new(parameters, limits())
        .unwrap()
        .with_tabulated_hardening(table)
        .unwrap();
    original_pressure_parity(&FiniteNortonJ2::new(tabulated, kinetics(), DT, 0.8).unwrap());
}

fn independent_derivatives<L: MixedFiniteInelasticMaterial>(law: &L) {
    let old = law.response(&F, &law.initial_state()).unwrap().state;
    let frozen = law.freeze_remainder(&old).unwrap();
    for treatment in [
        VolumetricTreatment::OriginalMaterial,
        VolumetricTreatment::ExactIncompressible,
    ] {
        let pressure = 3e6;
        let response = frozen.mixed_response(&F, pressure, treatment).unwrap();
        for col in 0..9 {
            let h = 2e-6;
            let (mut plus, mut minus) = (F, F);
            plus[col] += h;
            minus[col] -= h;
            let plus = frozen.mixed_state(&plus, pressure, treatment).unwrap();
            let minus = frozen.mixed_state(&minus, pressure, treatment).unwrap();
            near(
                (plus.potential_pa - minus.potential_pa) / (2.0 * h),
                response.state.first_piola_pa[col],
                MU,
                3e-8,
            );
            for row in 0..9 {
                near(
                    (plus.first_piola_pa[row] - minus.first_piola_pa[row]) / (2.0 * h),
                    response.tangent_pa[row][col],
                    MU,
                    3e-7,
                );
            }
        }
    }
}

#[test]
fn g0_kinetic_mixed_value_and_stress_match_independent_differences_at_fixed_history() {
    independent_derivatives(&perzyna(K));
    independent_derivatives(&norton(K).with_strain_hardening(primary()).unwrap());
    independent_derivatives(&recovering_combined(K));
}

fn bulk_independence<L: MixedFiniteInelasticMaterial>(ordinary: &L, enormous: &L)
where
    L::State: PartialEq,
{
    let f = [0.7, 0.02, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.7];
    let old = ordinary.initial_state();
    let expected = ordinary.remainder_response(&f, &old).unwrap();
    let actual = enormous.remainder_response(&f, &old).unwrap();
    assert_eq!(actual, expected);
    enormous
        .validate_remainder_committed(&f, &actual.state, 1e-4)
        .unwrap();
    let a = ordinary
        .freeze_remainder(&old)
        .unwrap()
        .mixed_response(&f, 3e6, VolumetricTreatment::ExactIncompressible)
        .unwrap();
    let b = enormous
        .freeze_remainder(&old)
        .unwrap()
        .mixed_response(&f, 3e6, VolumetricTreatment::ExactIncompressible)
        .unwrap();
    assert_eq!(a.state.first_piola_pa, b.state.first_piola_pa);
    assert_eq!(a.tangent_pa, b.tangent_pa);
    assert_eq!(b.state.bulk_stored_pa, 0.0);
    assert!(enormous.response(&f, &old).is_err());
}

#[test]
fn g3_exact_kinetic_constraint_never_evaluates_unused_bulk_penalty() {
    bulk_independence(&perzyna(K), &perzyna(1e308));
    bulk_independence(&norton(K), &norton(1e308));
    bulk_independence(&combined(K), &combined(1e308));
    bulk_independence(&recovering_norton(K), &recovering_norton(1e308));
    bulk_independence(&recovering_combined(K), &recovering_combined(1e308));
}

#[test]
fn g1_unloaded_recovery_retains_permanent_history_and_inspection_does_not_advance_time() {
    let law = recovering_norton(K);
    let mut old = law.initial_state();
    old.native.alpha = 0.1;
    old.resistance_strain = 0.07;
    let before = old.clone();
    law.validate_remainder_committed(&I, &old, 0.0).unwrap();
    assert_eq!(law.stored_remainder_energy(&I, &old).unwrap(), 0.0);
    let frozen = law.freeze_remainder(&old).unwrap();
    for _ in 0..3 {
        let response = frozen
            .mixed_state(&I, 0.0, VolumetricTreatment::ExactIncompressible)
            .unwrap();
        assert_eq!(response.first_piola_pa, [0.0; 9]);
    }
    let returned = law.remainder_response(&I, &old).unwrap();
    assert_eq!(returned.state.native, old.native);
    near(
        returned.state.resistance_strain,
        old.resistance_strain / (1.0 + 1.5 * DT),
        1.0,
        2e-15,
    );
    assert_eq!(returned.inelastic_dissipation_j_m3(), 0.0);
    assert_eq!(returned.flow_potential_j_m3, 0.0);
    assert_eq!(old, before);
    let mut invalid = old.clone();
    invalid.resistance_strain = 2.0 * invalid.native.alpha;
    assert!(law.freeze_remainder(&invalid).is_err());
    assert!(law.stored_remainder_energy(&I, &invalid).is_err());
    assert_eq!(old, before);
}
