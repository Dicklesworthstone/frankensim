//! Source-integrated eigenstrain, physical restraint, native work and refusal.
use super::{card, near, properties, scalar, settings};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyValue, QueryPoint,
};
use fs_material::SmallStrainLaw;
use fs_material::brittle::source::expansion::{ResolvedThermalBrittle, resolve_thermal_brittle};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, YOUNG_MODULUS_PROPERTY,
};
use fs_material::tensor::contract;
use fs_qty::{Pressure, QuantitySpec, Temperature};

fn query(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
}

fn curve_card() -> MaterialCard {
    let mut claims = properties();
    let mut alpha = scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        1e-5,
    );
    alpha.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(250.0, 0.5e-5), (300.0, 1e-5), (350.0, 3e-5), (400.0, 1e-5)],
        dims: INVERSE_TEMPERATURE_DIMS,
    };
    alpha.validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        250.0,
        400.0,
    );
    alpha.interpolation = InterpolationPolicy::LinearInside;
    claims.push(alpha);
    card(claims)
}

fn resolve(source: &MaterialCard, reference: f64, current: f64) -> ResolvedThermalBrittle {
    resolve_thermal_brittle(
        source,
        &query(reference),
        &query(current),
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
    )
    .unwrap()
}

#[test]
fn integral_crosses_every_alpha_knot_and_retains_the_original_reference() {
    let source = curve_card();
    let start = resolve(&source, 300.0, 300.0);
    let hot = start.at_temperature(&source, 400.0).unwrap();
    near(start.free_linear_strain(), 0.0);
    // Two trapezoids: endpoint alpha * delta T would incorrectly give 0.001.
    near(hot.free_linear_strain(), 0.002);
    near(hot.expansion().current().linear_coefficient_per_k(), 1e-5);
    assert_eq!(hot.band().material().properties().len(), 8);
    assert_eq!(
        hot.expansion().selected_claim(),
        start.expansion().selected_claim()
    );
    assert_eq!(
        hot.expansion().reference().resolved().identity(),
        start.expansion().reference().resolved().identity(),
    );
    let cool = hot.at_temperature(&source, 250.0).unwrap();
    near(cool.free_linear_strain(), -0.000375);
    near(
        resolve(&source, 350.0, 250.0).free_linear_strain(),
        -0.001375,
    );
    assert!(hot.at_temperature(&source, 400.1).is_err());
    assert_eq!(
        start.band().material().card_identity(),
        source.content_hash()
    );
}

#[test]
fn free_cooling_is_stress_free_while_restraint_drives_tensile_fracture_without_healing() {
    let source = curve_card();
    let start = resolve(&source, 400.0, 400.0);
    let cold = start.at_temperature(&source, 300.0).unwrap();
    let virgin = start.band().law().initial_state();
    let saved = virgin;

    // Coefficient-only preparation permits the displacement solve to expand or
    // contract freely, without accepting a fictitious restrained damage trial.
    let prepared = start.transfer_coefficients(&virgin, &cold, 1e-6).unwrap();
    near(prepared.parameter_work_j_m3, 0.0);
    near(prepared.damage_dissipated_j_m3, 0.0);
    let free = cold.free_linear_strain();
    let expanded = cold
        .trial_total(&[free, free, free, 0.0, 0.0, 0.0], &prepared.state)
        .unwrap();
    assert_eq!(expanded.stress_pa, [0.0; 6]);
    assert_eq!(expanded.state.damage(), 0.0);
    assert_eq!(expanded.state.strain(), &[0.0; 6]);
    near(expanded.energy.work, 0.0);

    // Holding total strain at zero instead creates positive mechanical opening.
    let restrained = start.transfer_temperature(&virgin, &cold, 1e-6).unwrap();
    assert!(restrained.stress_pa[0] > 0.0);
    assert!(restrained.state.damage() > 0.0);
    assert!(restrained.eigenstrain_damage_dissipated_j_m3 > 0.0);
    near(restrained.coefficient_parameter_work_j_m3, 0.0);
    near(restrained.coefficient_damage_dissipated_j_m3, 0.0);
    near(
        restrained.eigenstrain_work_j_m3,
        restrained.stored_change_j_m3 + restrained.damage_dissipated_j_m3,
    );
    near(restrained.eigenstrain_work_defect_j_m3, 0.0);
    near(restrained.defect_j_m3, 0.0);
    for component in cold.total_strain(&restrained.state).unwrap() {
        near(component, 0.0);
    }

    let reheated = cold
        .transfer_temperature(&restrained.state, &start, 1e-6)
        .unwrap();
    near(reheated.state.damage(), restrained.state.damage());
    near(
        reheated.state.dissipated_j_m3(),
        restrained.state.dissipated_j_m3(),
    );
    near(reheated.damage_dissipated_j_m3, 0.0);
    assert_eq!(virgin, saved);
}

#[test]
fn coefficient_and_eigenstrain_work_remain_separate_and_use_native_mean_stress() {
    let mut claims = properties();
    let modulus = claims
        .iter_mut()
        .find(|p| p.key.name() == YOUNG_MODULUS_PROPERTY)
        .unwrap();
    modulus.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1e9), (400.0, 0.5e9)],
        dims: Pressure::DIMS,
    };
    modulus.validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        300.0,
        400.0,
    );
    modulus.interpolation = InterpolationPolicy::LinearInside;
    claims.push(scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        1e-6,
    ));
    let source = card(claims);
    let start = resolve(&source, 300.0, 300.0);
    let hot = start.at_temperature(&source, 400.0).unwrap();
    let total = [0.0003, 0.0002, -0.0001, 0.0001, 0.0002, -0.00015];
    let loaded = start
        .trial_total(&total, &start.band().law().initial_state())
        .unwrap();
    assert_eq!(loaded.state.damage(), 0.0);
    let coefficient = start
        .transfer_coefficients(&loaded.state, &hot, 1e-6)
        .unwrap();
    near(
        coefficient.parameter_work_j_m3,
        -0.5 * loaded.state.stored_energy_j_m3(),
    );

    let mechanical = hot.mechanical_strain(&total).unwrap();
    assert_eq!(&mechanical[3..], &total[3..]);
    let native = hot
        .band()
        .law()
        .trial_with_mean_stress(&mechanical, &coefficient.state)
        .unwrap();
    let thermal = start
        .transfer_temperature(&loaded.state, &hot, 1e-6)
        .unwrap();
    assert_eq!(thermal.state, native.endpoint.state);
    assert_eq!(thermal.stress_pa, native.endpoint.stress_pa);
    assert_eq!(thermal.mean_stress_pa, native.mean_stress_pa);
    near(
        thermal.coefficient_parameter_work_j_m3,
        coefficient.parameter_work_j_m3,
    );
    let free = hot.free_linear_strain();
    near(
        thermal.eigenstrain_work_j_m3,
        contract(
            &native.mean_stress_pa,
            &[-free, -free, -free, 0.0, 0.0, 0.0],
        ),
    );
    near(thermal.eigenstrain_work_j_m3, native.path_work_j_m3);
    near(
        thermal.coefficient_parameter_work_j_m3 + thermal.eigenstrain_work_j_m3,
        thermal.stored_change_j_m3 + thermal.damage_dissipated_j_m3,
    );
    assert!(thermal.coefficient_parameter_work_j_m3 < 0.0);
    assert!(thermal.eigenstrain_work_j_m3 < 0.0);
    assert_eq!(thermal.damage_dissipated_j_m3, 0.0);
    let adapter = hot.small_strain_law();
    assert_eq!(
        adapter.stress(&total, &coefficient.state),
        native.endpoint.stress_pa
    );
    assert_eq!(
        adapter.tangent(&total, &coefficient.state),
        native.endpoint.branch_tangent_pa
    );
    assert_eq!(
        adapter.update_state(&total, &coefficient.state),
        native.endpoint.state
    );
}

#[test]
fn negative_expansion_is_admitted_and_heating_under_restraint_can_create_tension() {
    let mut claims = properties();
    claims.push(scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        -2e-5,
    ));
    let source = card(claims);
    let start = resolve(&source, 300.0, 300.0);
    let hot = start.at_temperature(&source, 400.0).unwrap();
    near(hot.free_linear_strain(), -0.002);
    let step = start
        .transfer_temperature(&start.band().law().initial_state(), &hot, 1e-6)
        .unwrap();
    assert!(step.stress_pa[0] > 0.0);
    assert!(step.state.damage() > 0.0);
    assert!(step.eigenstrain_work_j_m3 > 0.0);
}

#[test]
fn source_reference_context_and_small_strain_refusals_leave_history_retryable() {
    let source = curve_card();
    let start = resolve(&source, 300.0, 300.0);
    let hot = start.at_temperature(&source, 400.0).unwrap();
    let virgin = start.band().law().initial_state();
    let saved = virgin;
    let reset_reference = resolve(&source, 350.0, 400.0);
    assert!(
        start
            .transfer_coefficients(&virgin, &reset_reference, 1e-6)
            .is_err()
    );
    assert!(start.transfer_temperature(&virgin, &hot, f64::NAN).is_err());

    let changed_context = query(400.0).with("pressure", 1e5).unwrap();
    assert!(
        resolve_thermal_brittle(
            &source,
            &query(300.0),
            &changed_context,
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        )
        .is_err()
    );
    let changed_schema = QueryPoint::new().with("T", 400.0).unwrap();
    assert!(
        resolve_thermal_brittle(
            &source,
            &query(300.0),
            &changed_schema,
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        )
        .is_err()
    );

    let missing = card(properties());
    assert!(
        resolve_thermal_brittle(
            &missing,
            &query(300.0),
            &query(400.0),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        )
        .is_err()
    );
    let mut wrong = properties();
    wrong.push(scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        Pressure::DIMS,
        1e-5,
    ));
    assert!(
        resolve_thermal_brittle(
            &card(wrong),
            &query(300.0),
            &query(400.0),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        )
        .is_err()
    );
    assert!(start.at_temperature(&missing, 400.0).is_err());

    // Each norm is tested independently: these points violate only one bound.
    assert!(
        hot.mechanical_strain(&[0.03, 0.03, 0.03, 0.0, 0.0, 0.0])
            .is_err()
    );
    assert!(
        hot.mechanical_strain(&[-0.028, -0.028, -0.028, 0.0, 0.0, 0.0])
            .is_err()
    );
    let cloned = hot.clone();
    let foreign = cloned.band().law().initial_state();
    assert!(hot.trial_total(&[0.0; 6], &foreign).is_err());
    let marker = hot.small_strain_law().update_state(&[0.0; 6], &foreign);
    assert!(!marker.is_valid());
    assert!(hot.trial_total(&[0.0; 6], &marker).is_err());
    assert_eq!(virgin, saved);
    assert!(start.transfer_coefficients(&virgin, &hot, 1e-6).is_ok());
}

#[test]
fn numerically_identical_alpha_cannot_silently_switch_selected_evidence() {
    let mut claims = ClaimSet::new();
    let mut pins = Vec::new();
    for property in properties() {
        let name = property.key.name().to_owned();
        let selected = claims.insert_claim(property).unwrap();
        pins.push((name, selected));
    }
    let alpha = scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        1e-5,
    );
    let mut other_alpha = alpha.clone();
    other_alpha.provenance.source = "different manufactured expansion evidence".into();
    let first = claims.insert_claim(alpha).unwrap();
    let second = claims.insert_claim(other_alpha).unwrap();
    assert_ne!(first, second);
    let source = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle coupon".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let selected = |claim| {
        let mut selected = pins.clone();
        selected.push((
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(),
            claim,
        ));
        MaterialPropertySelection::PinnedByProperty(selected)
    };
    let a = resolve_thermal_brittle(
        &source,
        &query(300.0),
        &query(300.0),
        selected(first),
        settings(0.01),
    )
    .unwrap();
    let b = resolve_thermal_brittle(
        &source,
        &query(300.0),
        &query(400.0),
        selected(second),
        settings(0.01),
    )
    .unwrap();
    let state = a.band().law().initial_state();
    assert!(a.transfer_coefficients(&state, &b, 1e-6).is_err());
    assert!(a.transfer_temperature(&state, &b, 1e-6).is_err());
    assert_eq!(state.damage(), 0.0);
}

#[path = "expansion/dynamics.rs"]
mod dynamics;
