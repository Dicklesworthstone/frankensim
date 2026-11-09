use super::*;
use crate::eos::IdealGasMolarVolumeQuantityV1 as MolarVolume;
use crate::eos::peng_robinson::saturation::PengRobinsonSaturationV1;
use crate::eos::peng_robinson::{
    PengRobinsonComponentV1, saturation::PengRobinsonSaturationOptionsV1,
};
use crate::{
    ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    Nasa9RegionV1, SpeciesId, StandardStateConventionV1, required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use fs_qty::Pressure;
use std::collections::BTreeMap;

fn component() -> PengRobinsonComponentV1 {
    PengRobinsonComponentV1 {
        critical_temperature: Temperature::new(190.564),
        critical_pressure: Pressure::new(4.5992e6),
        acentric_factor: 0.01142,
        molar_mass: MolarMass::new(0.016_042_46),
    }
}
fn model() -> PengRobinsonMixtureV1 {
    PengRobinsonMixtureV1::new(vec![component()], vec![1.0], vec![]).unwrap()
}
fn source(formation: f64, entropy: f64) -> Nasa9StandardStateModelV1 {
    let mut a = [0.0; 9];
    a[2] = 4.0;
    a[7] = formation / R;
    a[8] = entropy / R;
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1e5,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 100.0, 500.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "constant-Cp pure-fluid caloric oracle".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new("CH4").unwrap(),
        component().molar_mass,
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1e5),
            ElementalReferenceIdV1::new("pure-fluid-test").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn saturation(model: &PengRobinsonMixtureV1, t: f64) -> PengRobinsonSaturationV1 {
    model
        .pure_saturation_at_temperature(
            Temperature::new(t),
            [Pressure::new(1e5), Pressure::new(4.4e6)],
            PengRobinsonSaturationOptionsV1 {
                log_fugacity_tolerance: 1e-13,
                relative_pressure_tolerance: 1e-13,
                ..Default::default()
            },
        )
        .unwrap()
}
fn state(model: &PengRobinsonMixtureV1, t: f64, volume: f64) -> PengRobinsonSaturatedMixtureV1 {
    saturation(model, t)
        .at_molar_volume(MolarVolume::new(volume))
        .unwrap()
}
fn reference<'a>(
    model: &'a PengRobinsonMixtureV1,
    source: &'a Nasa9StandardStateModelV1,
    initial: PengRobinsonSaturatedMixtureV1,
) -> PengRobinsonIsochoricCaloricReferenceV1<'a> {
    let density = Density::new(initial.mass_density_kg_m3());
    PengRobinsonIsochoricCaloricReferenceV1::new(
        model,
        source,
        initial,
        density,
        [Temperature::new(140.0), Temperature::new(160.0)],
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e} (allowance {tolerance})"
    );
}
fn absolute_u(s: &PengRobinsonSaturatedMixtureV1) -> f64 {
    (3.0 * R * s.saturation().temperature().value() + s.residual_internal_energy_j_mol())
        / component().molar_mass.value()
}
fn absolute_s(s: &PengRobinsonSaturatedMixtureV1) -> f64 {
    let sat = s.saturation();
    (4.0 * R * sat.temperature().value().ln() - R * (sat.pressure().value() / 1e5).ln()
        + s.residual_entropy_j_mol_k())
        / component().molar_mass.value()
}

#[test]
fn g1_rigid_coexistence_retains_original_energy_entropy_and_volume() {
    let model = model();
    let source = source(0.0, 0.0);
    let initial = saturation(&model, 150.0).at_vapor_fraction(0.35).unwrap();
    let volume = initial.molar_volume().value();
    let path = reference(&model, &source, initial.clone());
    for t in [140.0, 145.0, 150.0, 155.0, 160.0] {
        let final_state = state(&model, t, volume);
        let change = path.change_to(&final_state, || false).unwrap();
        close(
            change.specific_internal_energy_change_j_kg,
            absolute_u(&final_state) - absolute_u(&initial),
            2e-7,
        );
        close(
            change.source_entropy_change_j_kg_k,
            absolute_s(&final_state) - absolute_s(&initial),
            2e-8,
        );
        assert!(change.relative_volume_residual.abs() <= 128.0 * f64::EPSILON);
        assert!(change.response.molar_cv_j_mol_k >= 3.0 * R);
        assert_eq!(
            change.caloric_evaluations,
            if final_state == initial { 0 } else { 4 }
        );
        if t != 150.0 {
            assert_ne!(final_state.vapor_fraction(), initial.vapor_fraction());
            assert_ne!(
                final_state.saturation().pressure(),
                initial.saturation().pressure()
            );
            assert_eq!(
                change
                    .specific_internal_energy_change_j_kg
                    .is_sign_positive(),
                t > 150.0
            );
        }
    }
}

#[test]
fn g1_moving_phase_cv_and_entropy_are_the_fixed_volume_derivatives() {
    let model = model();
    let source = source(0.0, 0.0);
    let initial = saturation(&model, 150.0).at_vapor_fraction(0.35).unwrap();
    let volume = initial.molar_volume().value();
    let path = reference(&model, &source, initial.clone());
    let step = 0.001;
    let minus = path
        .change_to(&state(&model, 150.0 - step, volume), || false)
        .unwrap();
    let plus = path
        .change_to(&state(&model, 150.0 + step, volume), || false)
        .unwrap();
    let response = path.response_at(&initial, || false).unwrap();
    let cv = response.molar_cv_j_mol_k / component().molar_mass.value();
    close(
        (plus.specific_internal_energy_change_j_kg - minus.specific_internal_energy_change_j_kg)
            / (2.0 * step),
        cv,
        0.002,
    );
    close(
        (plus.source_entropy_change_j_kg_k - minus.source_entropy_change_j_kg_k) / (2.0 * step),
        cv / 150.0,
        2e-5,
    );
    let frozen = 3.0 * R
        + (1.0 - initial.vapor_fraction()) * initial.saturation().liquid().residual_cv_j_mol_k()
        + initial.vapor_fraction() * initial.saturation().vapor().residual_cv_j_mol_k();
    assert!(response.molar_cv_j_mol_k > frozen);
}

#[test]
fn g3_large_reference_datums_preserve_weak_local_heat_and_entropy() {
    let model = model();
    let ordinary = source(0.0, 0.0);
    let shifted = source(1e14, -1e12);
    let initial = saturation(&model, 150.0).at_vapor_fraction(0.35).unwrap();
    let final_state = state(&model, 150.0 + 1e-8, initial.molar_volume().value());
    let a = reference(&model, &ordinary, initial.clone())
        .change_to(&final_state, || false)
        .unwrap();
    let b = reference(&model, &shifted, initial.clone())
        .change_to(&final_state, || false)
        .unwrap();
    assert_eq!(a, b);
    assert!(a.specific_internal_energy_change_j_kg > 0.0);
    let dt = final_state.saturation().temperature().value() - 150.0;
    close(
        a.specific_internal_energy_change_j_kg,
        a.response.molar_cv_j_mol_k / component().molar_mass.value() * dt,
        2e-6,
    );
    close(
        a.source_entropy_change_j_kg_k,
        a.specific_internal_energy_change_j_kg / 150.0,
        2e-8,
    );
    let absolute = (1e14 + 3.0 * R * 150.0 + initial.residual_internal_energy_j_mol())
        / component().molar_mass.value();
    assert_eq!(absolute + a.specific_internal_energy_change_j_kg, absolute);
}

#[test]
fn g3_direct_source_changes_compose_and_reverse_without_a_second_latent_account() {
    let model = model();
    let source = source(0.0, 0.0);
    let a = saturation(&model, 145.0).at_vapor_fraction(0.25).unwrap();
    let volume = a.molar_volume().value();
    let b = state(&model, 150.0, volume);
    let c = state(&model, 155.0, volume);
    let pa = reference(&model, &source, a.clone());
    let ab = pa.change_to(&b, || false).unwrap();
    let ac = pa.change_to(&c, || false).unwrap();
    let bc = reference(&model, &source, b)
        .change_to(&c, || false)
        .unwrap();
    let ca = reference(&model, &source, c)
        .change_to(&a, || false)
        .unwrap();
    close(
        ac.specific_internal_energy_change_j_kg,
        ab.specific_internal_energy_change_j_kg + bc.specific_internal_energy_change_j_kg,
        2e-7,
    );
    close(
        ac.source_entropy_change_j_kg_k,
        ab.source_entropy_change_j_kg_k + bc.source_entropy_change_j_kg_k,
        2e-8,
    );
    close(
        ac.specific_internal_energy_change_j_kg,
        -ca.specific_internal_energy_change_j_kg,
        2e-7,
    );
    close(
        ac.source_entropy_change_j_kg_k,
        -ca.source_entropy_change_j_kg_k,
        2e-8,
    );
    let mean = ac.source_entropy_change_j_kg_k / ac.specific_internal_energy_change_j_kg;
    assert!(mean >= 1.0 / 155.0 && mean <= 1.0 / 145.0);
}

#[test]
fn g4_source_density_phase_domain_and_late_cancellation_refuse_atomically() {
    let model = model();
    let source = source(0.0, 0.0);
    let sat = saturation(&model, 150.0);
    let initial = sat.at_vapor_fraction(0.35).unwrap();
    let density = Density::new(initial.mass_density_kg_m3());
    let domain = [Temperature::new(140.0), Temperature::new(160.0)];
    assert!(
        PengRobinsonIsochoricCaloricReferenceV1::new(
            &model,
            &source,
            initial.clone(),
            Density::new(density.value() * 1.001),
            domain,
            || false
        )
        .is_err()
    );
    assert!(
        PengRobinsonIsochoricCaloricReferenceV1::new(
            &model,
            &source,
            initial.clone(),
            density,
            [Temperature::new(140.0), component().critical_temperature],
            || false
        )
        .is_err()
    );
    for q in [0.0, 1.0] {
        let end = sat.at_vapor_fraction(q).unwrap();
        let rho = Density::new(end.mass_density_kg_m3());
        assert!(
            PengRobinsonIsochoricCaloricReferenceV1::new(&model, &source, end, rho, domain, || {
                false
            })
            .is_err()
        );
    }
    let path = reference(&model, &source, initial.clone());
    assert!(
        path.change_to(&sat.at_vapor_fraction(0.4).unwrap(), || false)
            .is_err()
    );
    assert!(
        path.change_to(
            &state(&model, 161.0, initial.molar_volume().value()),
            || false
        )
        .is_err()
    );
    let final_state = state(&model, 155.0, initial.molar_volume().value());
    let mut count = 0;
    let expected = path
        .change_to(&final_state, || {
            count += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        path.change_to(&final_state, || {
            calls += 1;
            calls == count
        }),
        Err(E::Cancelled)
    ));
    assert_eq!(count, calls);
    assert_eq!(path.change_to(&final_state, || false).unwrap(), expected);
    let zero = path.change_to(path.initial_state(), || false).unwrap();
    assert_eq!(zero.specific_internal_energy_change_j_kg, 0.0);
    assert_eq!(zero.source_entropy_change_j_kg_k, 0.0);
    assert_eq!(zero.caloric_evaluations, 0);
}

use crate::eos::peng_robinson::saturation::isochoric::sealed::{
    PengRobinsonEquilibriumDensityV1 as SealedState, PengRobinsonIsochoricPhaseV1 as SealedPhase,
    PengRobinsonSealedCaloricReferenceV1 as SealedReference,
};

fn sealed_state(model: &PengRobinsonMixtureV1, t: f64, density: Density) -> SealedState {
    model
        .pure_equilibrium_density_at_saturation(saturation(model, t), density, || false)
        .unwrap()
}
fn sealed_reference<'a>(
    model: &'a PengRobinsonMixtureV1,
    source: &'a Nasa9StandardStateModelV1,
    state: SealedState,
) -> SealedReference<'a> {
    SealedReference::new(
        model,
        source,
        state,
        [Temperature::new(140.0), Temperature::new(180.0)],
        || false,
    )
    .unwrap()
}
fn sealed_absolute(state: &SealedState) -> (f64, f64) {
    let (u, s) = match state.phase() {
        SealedPhase::Coexistence(m) => (
            m.residual_internal_energy_j_mol(),
            m.residual_entropy_j_mol_k(),
        ),
        SealedPhase::Liquid(s) | SealedPhase::Vapor(s) => (
            s.root().residual_internal_energy_j_mol(),
            s.root().residual_entropy_j_mol_k(),
        ),
    };
    let t = state.temperature().value();
    let m = component().molar_mass.value();
    (
        (3.0 * R * t + u) / m,
        (4.0 * R * t.ln() - R * (state.pressure().value() / 1e5).ln() + s) / m,
    )
}
fn boundary_density(model: &PengRobinsonMixtureV1, vapor: bool) -> Density {
    let sat = saturation(model, 160.0);
    let volume = if vapor { sat.vapor() } else { sat.liquid() }.molar_volume_m3_per_mol();
    let m = component().molar_mass.value();
    let rho = m / volume;
    for rho in [rho, rho.next_down(), rho.next_up()] {
        if m / rho == volume {
            return Density::new(rho);
        }
    }
    panic!("manufactured boundary must be representable in the conserved density coordinate");
}

#[test]
fn g1_sealed_source_crosses_both_dome_boundaries_with_exact_original_calorics() {
    let model = model();
    let source = source(0.0, 0.0);
    for vapor in [false, true] {
        let density = boundary_density(&model, vapor);
        let initial = sealed_state(&model, 140.0, density);
        assert!(matches!(initial.phase(), SealedPhase::Coexistence(_)));
        let base = sealed_absolute(&initial);
        let path = sealed_reference(&model, &source, initial);
        for t in [
            140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0,
        ] {
            let state = sealed_state(&model, t, density);
            let change = path.change_to(&state, || false).unwrap();
            let value = sealed_absolute(&state);
            close(
                change.specific_internal_energy_change_j_kg,
                value.0 - base.0,
                2e-7,
            );
            close(change.source_entropy_change_j_kg_k, value.1 - base.1, 2e-8);
            assert_eq!(state.density(), density);
            assert!(change.relative_volume_residual.abs() <= 128.0 * f64::EPSILON);
            assert!(change.response.molar_cv_j_mol_k >= 3.0 * R);
            assert!(change.response.dp_dt_v_pa_k > 0.0);
            assert_eq!(
                change.caloric_evaluations,
                path.caloric_evaluations_to(&state)
            );
            if t < 160.0 {
                assert!(matches!(state.phase(), SealedPhase::Coexistence(_)));
            }
            if t > 160.0 {
                assert_eq!(state.homogeneous_evaluations(), 1);
                assert_eq!(matches!(state.phase(), SealedPhase::Vapor(_)), vapor);
                assert_eq!(matches!(state.phase(), SealedPhase::Liquid(_)), !vapor);
            }
        }
    }
}

#[test]
fn g1_sealed_boundary_retains_distinct_one_sided_cv_and_pressure_tangents() {
    let model = model();
    let source = source(0.0, 0.0);
    for vapor in [false, true] {
        let density = boundary_density(&model, vapor);
        let boundary = sealed_state(&model, 160.0, density);
        assert!(boundary.is_phase_boundary());
        let path = sealed_reference(&model, &source, boundary.clone());
        let response = path.response_at(&boundary, || false).unwrap();
        assert!(response.phase_boundary);
        let dt = 1e-4;
        let lower = sealed_state(&model, 160.0 - dt, density);
        let upper = sealed_state(&model, 160.0 + dt, density);
        let a = path.change_to(&lower, || false).unwrap();
        let b = path.change_to(&upper, || false).unwrap();
        let m = component().molar_mass.value();
        close(
            -a.specific_internal_energy_change_j_kg / dt,
            response.molar_cv_j_mol_k / m,
            0.2,
        );
        close(
            b.specific_internal_energy_change_j_kg / dt,
            b.response.molar_cv_j_mol_k / m,
            0.2,
        );
        close(
            -a.source_entropy_change_j_kg_k / dt,
            response.molar_cv_j_mol_k / (m * 160.0),
            0.002,
        );
        close(
            b.source_entropy_change_j_kg_k / dt,
            b.response.molar_cv_j_mol_k / (m * 160.0),
            0.002,
        );
        assert!(response.molar_cv_j_mol_k > b.response.molar_cv_j_mol_k);
        assert!(!a.response.phase_boundary && !b.response.phase_boundary);
        let outside = match upper.phase() {
            SealedPhase::Liquid(s) | SealedPhase::Vapor(s) => s,
            _ => panic!("outside dome"),
        };
        assert_eq!(
            b.response.molar_cv_j_mol_k,
            3.0 * R + outside.root().residual_cv_j_mol_k()
        );
        assert_eq!(
            b.response.dp_dt_v_pa_k,
            outside.root().pressure_derivatives().0
        );
    }
}

#[test]
fn g3_sealed_cross_phase_increments_compose_reverse_and_cancel_large_datums() {
    let model = model();
    let ordinary = source(0.0, 0.0);
    let shifted = source(1e14, -1e12);
    for vapor in [false, true] {
        let density = boundary_density(&model, vapor);
        let a = sealed_state(&model, 145.0, density);
        let b = sealed_state(&model, 160.0, density);
        let c = sealed_state(&model, 175.0, density);
        let pa = sealed_reference(&model, &ordinary, a.clone());
        let ab = pa.change_to(&b, || false).unwrap();
        let ac = pa.change_to(&c, || false).unwrap();
        let bc = sealed_reference(&model, &ordinary, b)
            .change_to(&c, || false)
            .unwrap();
        let ca = sealed_reference(&model, &ordinary, c.clone())
            .change_to(&a, || false)
            .unwrap();
        close(
            ac.specific_internal_energy_change_j_kg,
            ab.specific_internal_energy_change_j_kg + bc.specific_internal_energy_change_j_kg,
            2e-7,
        );
        close(
            ac.source_entropy_change_j_kg_k,
            ab.source_entropy_change_j_kg_k + bc.source_entropy_change_j_kg_k,
            2e-8,
        );
        close(
            ac.specific_internal_energy_change_j_kg,
            -ca.specific_internal_energy_change_j_kg,
            2e-7,
        );
        close(
            ac.source_entropy_change_j_kg_k,
            -ca.source_entropy_change_j_kg_k,
            2e-8,
        );
        let shifted_path = sealed_reference(&model, &shifted, a);
        assert_eq!(ac, shifted_path.change_to(&c, || false).unwrap());
        let weak = sealed_state(&model, 175.0 + 1e-8, density);
        let x = sealed_reference(&model, &ordinary, c.clone())
            .change_to(&weak, || false)
            .unwrap();
        let y = sealed_reference(&model, &shifted, c)
            .change_to(&weak, || false)
            .unwrap();
        assert_eq!(x, y);
        assert!(x.specific_internal_energy_change_j_kg > 0.0);
        assert_eq!(x.caloric_evaluations, 2);
        close(
            x.source_entropy_change_j_kg_k,
            x.specific_internal_energy_change_j_kg / 175.0,
            1e-10,
        );
        let mean = ac.source_entropy_change_j_kg_k / ac.specific_internal_energy_change_j_kg;
        assert!(mean >= 1.0 / 175.0 && mean <= 1.0 / 145.0);
    }
}

#[test]
fn g0_sealed_phase_admission_excludes_mechanically_stable_metastability() {
    let model = model();
    let sat = saturation(&model, 150.0);
    let v = sat.liquid().molar_volume_m3_per_mol() * 1.01;
    let rho = Density::new(component().molar_mass.value() / v);
    let homogeneous = model
        .evaluate_density_temperature(Temperature::new(150.0), rho)
        .unwrap();
    assert!(homogeneous.root().pressure_derivatives().1 < 0.0);
    let actual = model
        .pure_equilibrium_density_at_saturation(sat, rho, || false)
        .unwrap();
    assert!(matches!(actual.phase(), SealedPhase::Coexistence(_)));
    assert_eq!(actual.homogeneous_evaluations(), 0);
    assert_ne!(actual.pressure(), homogeneous.pressure());
}

#[test]
fn g4_sealed_binding_domain_work_and_final_cancellation_remain_atomic() {
    let model = model();
    let source = source(0.0, 0.0);
    let rho = boundary_density(&model, true);
    let initial = sealed_state(&model, 145.0, rho);
    let path = sealed_reference(&model, &source, initial.clone());
    let final_state = sealed_state(&model, 175.0, rho);
    assert!(
        path.change_to(
            &sealed_state(&model, 175.0, Density::new(rho.value().next_up())),
            || false
        )
        .is_err()
    );
    assert!(
        path.change_to(&sealed_state(&model, 181.0, rho), || false)
            .is_err()
    );
    let mut other = component();
    other.acentric_factor += 0.01;
    let other_model = PengRobinsonMixtureV1::new(vec![other], vec![1.0], vec![]).unwrap();
    assert!(
        other_model
            .pure_equilibrium_density_at_saturation(saturation(&model, 150.0), rho, || false)
            .is_err()
    );
    let mut count = 0;
    let expected = path
        .change_to(&final_state, || {
            count += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        path.change_to(&final_state, || {
            calls += 1;
            calls == count
        }),
        Err(E::Cancelled)
    ));
    assert_eq!(path.change_to(&final_state, || false).unwrap(), expected);
    assert_eq!(
        path.change_to(&initial, || false)
            .unwrap()
            .caloric_evaluations,
        0
    );
    assert!(
        SealedReference::new(
            &model,
            &source,
            initial,
            [Temperature::new(140.0), component().critical_temperature],
            || false
        )
        .is_err()
    );
}
