//! Exact axial creep/plasticity, independent solutions and original source lineage.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::plastic::source::NortonJ2Source;
use fs_material::creep::plastic::uniaxial::UniaxialNortonJ2;
use fs_material::creep::source::*;
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::creep::{CreepTolerances, NortonCreep, NortonCreepParameters};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_material::{IsotropicElastic, Uniaxial};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const E: f64 = 1.2e9;
const NU: f64 = 0.2;
const Y: f64 = 2.0e6;
const H: f64 = 2.0e8;
const RATE: f64 = 2.0e-3;
const REFERENCE: f64 = 3.0e6;

fn law(n: f64, dt: f64, hardening: f64) -> UniaxialNortonJ2 {
    let elastic = IsotropicElastic::new(E, NU, 0.05).unwrap();
    let creep = NortonCreep::new(
        elastic.clone(),
        NortonCreepParameters {
            reference_stress_pa: REFERENCE,
            reference_rate_s: RATE,
            stress_exponent: n,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 500.0,
            temperature_range_k: [300.0, 900.0],
        },
        500.0,
        dt,
    )
    .unwrap();
    UniaxialNortonJ2::new(
        NortonJ2::new(creep, J2Plasticity::new(elastic, Y, hardening).unwrap()).unwrap(),
    )
}
fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 2e-6,
        energy_j_m3: 1e-7,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; tolerance={tolerance}"
    );
}

#[test]
fn linear_axial_solution_uses_e_with_separate_flow_and_native_energy() {
    let dt = 0.25;
    let a = dt * RATE / REFERENCE;
    for hardening in [0.0, H] {
        let law = law(1.0, dt, hardening);
        for strain in [-0.006_f64, 0.006] {
            let old = law.initial_state();
            let next = law.response(strain, &old).unwrap();
            let q = (hardening * E * strain.abs() + E * Y) / (E + hardening + E * hardening * a);
            let dc = a * q;
            let dp = strain.abs() - q / E - dc;
            close(next.stress_pa, strain.signum() * q, 2e-8);
            close(
                next.tangent_pa,
                E * hardening / (E + hardening + E * hardening * a),
                2e-6,
            );
            close(next.state.creep_strain(), strain.signum() * dc, 2e-17);
            close(next.state.plastic_strain(), strain.signum() * dp, 2e-17);
            close(next.state.equivalent_creep_strain(), dc, 2e-17);
            close(next.state.equivalent_plastic_strain(), dp, 2e-17);
            close(
                next.strain[1],
                -NU * next.stress_pa / E - 0.5 * strain.signum() * (dc + dp),
                2e-17,
            );
            let checked = law
                .advance_with_energy(0.0, &old, strain, tolerances())
                .unwrap();
            close(checked.stress_pa[0], next.stress_pa, 2e-8);
            for stress in &checked.stress_pa[1..] {
                close(*stress, 0.0, 2e-8);
            }
            close(checked.energy.creep_dissipation_j_m3, q * dc, 2e-8);
            close(checked.energy.plastic_dissipation_j_m3, Y * dp, 2e-8);
            let stored = q * q / (2.0 * E) + 0.5 * hardening * dp * dp;
            close(checked.energy.final_stored_j_m3, stored, 2e-8);
            close(checked.energy.numerical_dissipation_j_m3, stored, 2e-8);
            close(
                checked.energy.mechanical_work_j_m3,
                next.stress_pa * strain,
                2e-8,
            );
            assert_eq!(old, law.initial_state());
        }
    }
}

#[test]
fn inactive_plasticity_and_zero_stress_match_the_original_axial_norton_law() {
    let law = law(1.0, 0.25, H);
    let original = UniaxialNorton::new(law.backbone().creep().clone());
    let old = law.initial_state();
    for strain in [0.0, 0.001, -0.001] {
        let next = law.response(strain, &old).unwrap();
        let reference = original
            .response(strain, &original.initial_state())
            .unwrap();
        assert_eq!(next.stress_pa, reference.stress_pa);
        assert_eq!(next.tangent_pa, reference.tangent_pa);
        assert_eq!(next.state.creep_strain(), reference.state.creep_strain());
        assert_eq!(
            next.state.equivalent_creep_strain(),
            reference.state.equivalent_creep_strain()
        );
        assert_eq!(next.state.plastic_strain(), 0.0);
        assert_eq!(next.state.equivalent_plastic_strain(), 0.0);
    }
}

#[test]
fn nonlinear_axial_tangent_is_the_free_sided_tensor_condensation() {
    let law = law(3.0, 0.125, H);
    let first = law.response(0.004, &law.initial_state()).unwrap();
    for strain in [0.007, -0.007] {
        let next = law.response(strain, &first.state).unwrap();
        let step = 1e-9;
        let fd = (law.response(strain + step, &first.state).unwrap().stress_pa
            - law.response(strain - step, &first.state).unwrap().stress_pa)
            / (2.0 * step);
        close(fd, next.tangent_pa, 2e-6 * next.tangent_pa.abs());
        let full = law
            .backbone()
            .response(&next.strain, &first.state.tensor_history())
            .unwrap();
        close(full.stress_pa[0], next.stress_pa, 2e-8);
        for stress in &full.stress_pa[1..] {
            close(*stress, 0.0, 2e-8);
        }
        let c = full.tangent_pa;
        let determinant = c[1][1] * c[2][2] - c[1][2] * c[2][1];
        let v1 = (c[2][2] * c[1][0] - c[1][2] * c[2][0]) / determinant;
        let v2 = (c[1][1] * c[2][0] - c[2][1] * c[1][0]) / determinant;
        close(next.tangent_pa, c[0][0] - c[0][1] * v1 - c[0][2] * v2, 2e-6);
        let checked = law
            .advance_with_energy(0.004, &first.state, strain, tolerances())
            .unwrap();
        close(
            checked.state.creep.equivalent_creep_strain,
            next.state.equivalent_creep_strain(),
            2e-17,
        );
        close(
            checked.state.plastic.alpha,
            next.state.equivalent_plastic_strain(),
            2e-17,
        );
    }
}

#[test]
fn constant_stress_adds_creep_without_repeating_initial_plastic_strain() {
    let dt = 0.1;
    let law = law(3.0, dt, H);
    let sigma: f64 = 3.0e6;
    let plastic = (sigma - Y) / H;
    let rate = RATE * (sigma / REFERENCE).powi(3);
    let mut old = law.initial_state();
    let mut old_strain = 0.0;
    for step in 1..=6 {
        let creep = step as f64 * dt * rate;
        let strain = sigma / E + plastic + creep;
        let next = law.response(strain, &old).unwrap();
        let checked = law
            .advance_with_energy(old_strain, &old, strain, tolerances())
            .unwrap();
        close(next.stress_pa, sigma, 2e-8);
        close(next.state.creep_strain(), creep, 2e-17);
        close(next.state.plastic_strain(), plastic, 2e-17);
        close(next.state.equivalent_plastic_strain(), plastic, 2e-17);
        close(
            checked.energy.creep_dissipation_j_m3,
            sigma * dt * rate,
            2e-8,
        );
        close(
            checked.energy.plastic_dissipation_j_m3,
            if step == 1 { Y * plastic } else { 0.0 },
            2e-8,
        );
        old = next.state;
        old_strain = strain;
    }
    let held = law.response(old_strain, &old).unwrap();
    assert!(held.stress_pa < sigma);
    assert!(held.state.creep_strain() > old.creep_strain());
    assert_eq!(
        held.state.equivalent_plastic_strain(),
        old.equivalent_plastic_strain()
    );
    let checked = law
        .advance_with_energy(old_strain, &old, old_strain, tolerances())
        .unwrap();
    close(checked.energy.mechanical_work_j_m3, 0.0, 1e-8);
    close(
        checked.energy.initial_stored_j_m3 - checked.energy.final_stored_j_m3,
        checked.energy.physical_dissipation_j_m3 + checked.energy.numerical_dissipation_j_m3,
        2e-8,
    );
}

#[test]
fn retained_history_activation_ties_and_adjacent_strains_remain_admissible() {
    let dt = 0.1;
    let law = law(3.0, dt, H);
    let sigma = REFERENCE;
    let plastic = (sigma - Y) / H;
    let first_strain = sigma / E + plastic + dt * RATE;
    let first = law.response(first_strain, &law.initial_state()).unwrap();
    let transition = sigma / E + plastic + 2.0 * dt * RATE;
    for offset in -8_i64..=8 {
        let strain = f64::from_bits((transition.to_bits() as i64 + offset) as u64);
        let next = law.response(strain, &first.state).unwrap();
        let checked = law
            .advance_with_energy(first_strain, &first.state, strain, tolerances())
            .unwrap();
        close(next.stress_pa, sigma, 1e-7);
        close(next.state.equivalent_plastic_strain(), plastic, 2e-17);
        assert!(checked.energy.plastic_dissipation_j_m3 >= 0.0);
    }
    let below = law.response(transition - 1e-7, &first.state).unwrap();
    let above = law.response(transition + 1e-7, &first.state).unwrap();
    assert_eq!(
        below.state.equivalent_plastic_strain(),
        first.state.equivalent_plastic_strain()
    );
    assert!(
        above.state.equivalent_plastic_strain() > first.state.equivalent_plastic_strain() + 1e-8
    );
}

#[test]
fn unloading_and_reversal_keep_both_signed_histories_and_grow_both_paths() {
    let law = law(1.0, 0.25, H);
    let first = law.response(0.006, &law.initial_state()).unwrap();
    let permanent = first.state.creep_strain() + first.state.plastic_strain();
    let unloaded = law.response(permanent, &first.state).unwrap();
    close(unloaded.stress_pa, 0.0, 2e-8);
    close(
        unloaded.state.creep_strain(),
        first.state.creep_strain(),
        2e-17,
    );
    close(
        unloaded.state.plastic_strain(),
        first.state.plastic_strain(),
        2e-17,
    );
    let checked = law
        .advance_with_energy(0.006, &first.state, permanent, tolerances())
        .unwrap();
    close(checked.energy.physical_dissipation_j_m3, 0.0, 1e-12);
    let reversed = law.response(-0.008, &unloaded.state).unwrap();
    assert!(reversed.state.creep_strain() < 0.0 && reversed.state.plastic_strain() < 0.0);
    assert!(reversed.state.equivalent_creep_strain() > first.state.equivalent_creep_strain());
    assert!(reversed.state.equivalent_plastic_strain() > first.state.equivalent_plastic_strain());
    law.advance_with_energy(permanent, &unloaded.state, -0.008, tolerances())
        .unwrap();
}

fn sourced_card() -> MaterialCard {
    let coefficients = [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), NU),
        (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, 2e6),
        (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 1e-5),
        (CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 3.0),
        (
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            30_000.0,
        ),
        (
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            400.0,
        ),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, Y),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H),
    ];
    let mut claims = ClaimSet::new();
    for (name, dims, value) in coefficients {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    500.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured axial combined-law reference".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic axial solid".into(),
            phase: "solid".into(),
            process: "manufactured reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn source_bound_axial_law_shares_all_receipts_and_continues_both_histories() {
    let card = sourced_card();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
            .unwrap()
    };
    let source = NortonJ2Source::new(
        &card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.02,
        },
        HardeningPropertyBasis::PlasticStrainModulus,
        || false,
    )
    .unwrap();
    let resolved = source.at_temperature(300.0, 0.1).unwrap();
    let cold = UniaxialNortonJ2::from_source(resolved.clone());
    assert_eq!(cold.source().unwrap(), &resolved);
    assert_eq!(
        cold.material().unwrap().card_identity(),
        card.content_hash()
    );
    assert_eq!(cold.material().unwrap().properties().len(), 9);
    let cloned = cold.clone();
    assert!(std::ptr::eq(
        cold.source().unwrap(),
        cloned.source().unwrap()
    ));
    let first = cold.response(0.006, &cold.initial_state()).unwrap();
    cold.advance_with_energy(0.0, &cold.initial_state(), 0.006, tolerances())
        .unwrap();
    let hot = UniaxialNortonJ2::from_source(source.at_temperature(500.0, 0.2).unwrap());
    let next = hot.response(0.006, &first.state).unwrap();
    let checked = hot
        .advance_with_energy(0.006, &first.state, 0.006, tolerances())
        .unwrap();
    assert!(next.state.equivalent_creep_strain() > first.state.equivalent_creep_strain());
    close(
        next.state.equivalent_plastic_strain(),
        first.state.equivalent_plastic_strain(),
        2e-17,
    );
    assert!(next.stress_pa < first.stress_pa);
    assert!(checked.energy.creep_dissipation_j_m3 > 0.0);
    close(checked.energy.mechanical_work_j_m3, 0.0, 1e-8);
    for property in hot.material().unwrap().properties() {
        let original = cold
            .material()
            .unwrap()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            original.answer().receipt.selected
        );
    }
}

#[test]
fn repeated_trials_and_refusals_do_not_advance_accepted_fiber_history() {
    let law = law(3.0, 0.125, H);
    let old = law.initial_state();
    let first = law.response(0.006, &old).unwrap();
    assert_eq!(first, law.response(0.006, &old).unwrap());
    assert_eq!(
        law.tensor_strain(0.006, &first.state).unwrap(),
        first.strain
    );
    let checked = law
        .advance_with_energy(0.0, &old, 0.006, tolerances())
        .unwrap();
    assert_eq!(
        checked,
        law.advance_with_energy(0.0, &old, 0.006, tolerances())
            .unwrap()
    );
    for invalid in [f64::NAN, f64::INFINITY, 0.06] {
        assert!(law.response(invalid, &old).is_err());
        assert!(law.stress(invalid, &old).is_nan());
    }
    let bad = CreepTolerances {
        energy_j_m3: 0.0,
        ..tolerances()
    };
    assert!(law.advance_with_energy(0.0, &old, 0.006, bad).is_err());
    assert_eq!(old, law.initial_state());
    assert!(law.source().is_none() && law.material().is_none());
}
