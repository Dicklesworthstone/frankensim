//! G0/G1: free-sided creep, original tensor equivalence and sourced continuation.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::source::{NortonCreepSource, NortonSourceQuery, norton_requirements};
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::creep::{CreepTolerances, NortonCreep, NortonCreepParameters};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::{IsotropicElastic, SmallStrainLaw, Uniaxial};

const E: f64 = 1.2e9;
const NU: f64 = 0.2;
const REFERENCE_STRESS: f64 = 3e6;
const REFERENCE_RATE: f64 = 2e-3;
fn law(n: f64, dt: f64) -> UniaxialNorton {
    UniaxialNorton::new(
        NortonCreep::new(
            IsotropicElastic::new(E, NU, 0.05).unwrap(),
            NortonCreepParameters {
                reference_stress_pa: REFERENCE_STRESS,
                reference_rate_s: REFERENCE_RATE,
                stress_exponent: n,
                activation_energy_j_mol: 0.0,
                reference_temperature_k: 500.0,
                temperature_range_k: [300.0, 900.0],
            },
            500.0,
            dt,
        )
        .unwrap(),
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
fn g1_linear_axial_creep_uses_young_modulus_and_complete_lateral_strain() {
    let law = law(1.0, 0.25);
    let old = law.initial_state();
    let factor = 1.0 + E * 0.25 * REFERENCE_RATE / REFERENCE_STRESS;
    for strain in [-0.004, 0.004] {
        let result = law.response(strain, &old).unwrap();
        let sigma = E * strain / factor;
        let creep = strain - sigma / E;
        close(result.stress_pa, sigma, 2e-8);
        close(result.tangent_pa, E / factor, 2e-6);
        close(result.state.creep_strain(), creep, 2e-17);
        close(result.state.equivalent_creep_strain(), creep.abs(), 2e-17);
        close(result.strain[1], -NU * sigma / E - 0.5 * creep, 2e-17);
        let checked = law
            .advance_with_energy(0.0, &old, strain, tolerances())
            .unwrap();
        close(
            checked.energy.physical_dissipation_j_m3,
            sigma * creep,
            2e-8,
        );
        close(
            checked.energy.final_stored_j_m3,
            sigma * sigma / (2.0 * E),
            2e-8,
        );
        close(
            checked.energy.numerical_dissipation_j_m3,
            sigma * sigma / (2.0 * E),
            2e-8,
        );
        close(checked.energy.mechanical_work_j_m3, sigma * strain, 2e-8);
        for transverse in &checked.stress_pa[1..] {
            close(*transverse, 0.0, 2e-8);
        }
    }
    let zero = law.response(0.0, &old).unwrap();
    assert_eq!(zero.state, old);
    close(zero.stress_pa, 0.0, 0.0);
    close(zero.tangent_pa, E / factor, 2e-6);
}

#[test]
fn g1_constant_stress_creeps_permanently_and_reversal_retains_path_length() {
    let law = law(3.0, 2.0);
    let stress: f64 = 1e6;
    let delta = 2.0 * REFERENCE_RATE * (stress / REFERENCE_STRESS).powi(3);
    let mut old = law.initial_state();
    let mut old_strain = stress / E;
    for step in 1..=8 {
        let strain = stress / E + step as f64 * delta;
        let response = law.response(strain, &old).unwrap();
        let checked = law
            .advance_with_energy(old_strain, &old, strain, tolerances())
            .unwrap();
        close(response.stress_pa, stress, 2e-8);
        close(response.state.creep_strain(), step as f64 * delta, 2e-17);
        close(
            checked.energy.physical_dissipation_j_m3,
            stress * delta,
            1e-8,
        );
        old = response.state;
        old_strain = strain;
    }
    let unloaded = law.response(old.creep_strain(), &old).unwrap();
    assert_eq!(unloaded.state, old);
    assert_eq!(unloaded.stress_pa, 0.0);
    let checked = law
        .advance_with_energy(old_strain, &old, old.creep_strain(), tolerances())
        .unwrap();
    close(checked.energy.physical_dissipation_j_m3, 0.0, 0.0);
    let reversed = law
        .response(old.creep_strain() - delta - stress / E, &old)
        .unwrap();
    close(reversed.stress_pa, -stress, 2e-8);
    close(
        reversed.state.creep_strain(),
        old.creep_strain() - delta,
        2e-17,
    );
    close(
        reversed.state.equivalent_creep_strain(),
        old.equivalent_creep_strain() + delta,
        2e-17,
    );
}

#[test]
fn g1_held_axial_strain_relaxation_converges_to_continuous_norton_reference() {
    let strain = 0.005;
    let sigma0 = E * strain;
    let duration = 0.4;
    let exact = sigma0
        / (1.0 + 2.0 * E * REFERENCE_RATE * sigma0.powi(2) / REFERENCE_STRESS.powi(3) * duration)
            .sqrt();
    let mut errors = Vec::new();
    for steps in [4, 16, 64] {
        let law = law(3.0, duration / steps as f64);
        let mut state = law.initial_state();
        let mut stress = sigma0;
        for _ in 0..steps {
            let response = law.response(strain, &state).unwrap();
            let checked = law
                .advance_with_energy(strain, &state, strain, tolerances())
                .unwrap();
            assert!(response.stress_pa < stress);
            assert!(checked.energy.physical_dissipation_j_m3 > 0.0);
            close(checked.energy.mechanical_work_j_m3, 0.0, 1e-8);
            state = response.state;
            stress = response.stress_pa;
        }
        errors.push((stress - exact).abs());
    }
    assert!(errors[1] < 0.32 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.28 * errors[1], "{errors:?}");
}

#[test]
fn g0_axial_tangent_matches_finite_differences_and_native_schur_complement() {
    for exponent in [1.0, 3.0, 8.0] {
        let law = law(exponent, 0.2);
        let old = law.response(0.004, &law.initial_state()).unwrap().state;
        for strain in [-0.003, 0.0, 0.006] {
            let response = law.response(strain, &old).unwrap();
            let h = 1e-8;
            let numerical =
                (law.stress(strain + h, &old) - law.stress(strain - h, &old)) / (2.0 * h);
            close(response.tangent_pa, numerical, 2e-5 * E);
            let native = law
                .backbone()
                .tangent(&response.strain, &old.tensor_history());
            let determinant = native[1][1] * native[2][2] - native[1][2] * native[2][1];
            let dy = (native[2][2] * native[1][0] - native[1][2] * native[2][0]) / determinant;
            let dz = (native[1][1] * native[2][0] - native[2][1] * native[1][0]) / determinant;
            let condensed = native[0][0] - native[0][1] * dy - native[0][2] * dz;
            close(response.tangent_pa, condensed, 3e-6);
            let tensor_stress = law
                .backbone()
                .stress(&response.strain, &old.tensor_history());
            close(response.stress_pa, tensor_stress[0], 3e-8);
            for transverse in &tensor_stress[1..] {
                close(*transverse, 0.0, 3e-8);
            }
        }
    }
}

#[test]
fn g0_source_receipts_and_history_survive_temperature_and_duration_changes() {
    let mut claims = ClaimSet::new();
    for (requirement, value) in norton_requirements().unwrap().iter().zip([
        E,
        NU,
        REFERENCE_STRESS,
        REFERENCE_RATE,
        3.0,
        45_000.0,
        500.0,
    ]) {
        let dims = requirement.dims();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(requirement.name(), dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 900.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured axial creep fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic axial creep solid".into(),
            phase: "solid".into(),
            process: "manufactured reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let low = QueryPoint::new().with("T", 300.0).unwrap();
    let high = QueryPoint::new().with("T", 900.0).unwrap();
    let source = NortonCreepSource::new(
        &card,
        NortonSourceQuery {
            lower: &low,
            upper: &high,
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.05,
        },
        || false,
    )
    .unwrap();
    let cold = UniaxialNorton::from_source(source.at_temperature(500.0, 0.2).unwrap());
    let state = cold.response(0.004, &cold.initial_state()).unwrap().state;
    let hot = UniaxialNorton::from_source(source.at_temperature(600.0, 0.2).unwrap());
    let cold_next = cold.response(0.004, &state).unwrap();
    let hot_next = hot.response(0.004, &state).unwrap();
    assert!(hot_next.state.creep_strain() > cold_next.state.creep_strain());
    assert!(hot_next.stress_pa < cold_next.stress_pa);
    for law in [&cold, &hot] {
        let material = law.material().unwrap();
        assert_eq!(material.card_identity(), card.content_hash());
        for property in material.properties() {
            assert_eq!(
                property.answer().receipt.selected,
                source
                    .reference()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
    let short = UniaxialNorton::from_source(source.at_temperature(600.0, 0.1).unwrap());
    assert!(
        short.response(0.004, &state).unwrap().state.creep_strain() < hot_next.state.creep_strain()
    );
}

#[test]
fn g0_trials_do_not_age_accepted_history_and_refuse_invalid_full_tensor_kinematics() {
    let law = law(3.0, 0.2);
    let old = law.initial_state();
    let snapshot = old.clone();
    let response = law.response(0.004, &old).unwrap();
    for _ in 0..5 {
        assert_eq!(law.response(0.004, &old).unwrap(), response);
        assert_eq!(law.stress(0.004, &old), response.stress_pa);
        assert_eq!(law.tangent(0.004, &old), response.tangent_pa);
        assert_eq!(law.update_state(0.004, &old), response.state);
        assert_eq!(old, snapshot);
    }
    for bad in [f64::NAN, f64::INFINITY, 0.06, -0.06] {
        assert!(law.response(bad, &old).is_err());
        assert!(
            law.advance_with_energy(0.0, &old, bad, tolerances())
                .is_err()
        );
        assert_eq!(old, snapshot);
    }
    // Axial strain alone fits; its mandatory lateral strain exceeds the envelope.
    assert!(law.tensor_strain(0.05, &old).is_err());
    assert!(
        law.advance_with_energy(
            0.0,
            &old,
            0.004,
            CreepTolerances {
                stress_pa: 0.0,
                energy_j_m3: 1e-7
            }
        )
        .is_err()
    );
}
