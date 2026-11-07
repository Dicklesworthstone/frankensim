//! Source receipts, reference binding and explicit instantaneous thermal rates.
use super::super::{card, near, properties, scalar, settings};
use super::query;
use fs_evidence::ValidityDomain;
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyValue};
use fs_material::brittle::source::BrittleSourceError;
use fs_material::brittle::source::expansion::rate::{
    resolve_thermal_viscous_brittle, resolve_thermal_viscous_brittle_dynamic_material,
    thermal_viscous_brittle_dynamic_requirements, thermal_viscous_brittle_requirements,
};
use fs_material::brittle::source::expansion::{ResolvedThermalBrittle, resolve_thermal_brittle};
use fs_material::brittle::source::rate::{
    BRITTLE_NORMAL_VISCOSITY_PROPERTY, BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
    BRITTLE_VISCOSITY_DIMS,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection,
};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};

fn claims() -> Vec<PropertyClaim> {
    let mut claims = properties();
    for (name, dims, values) in [
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            [1e-5, 3e-5],
        ),
        (
            BRITTLE_NORMAL_VISCOSITY_PROPERTY,
            BRITTLE_VISCOSITY_DIMS,
            [1e8, 2e8],
        ),
        (
            BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
            BRITTLE_VISCOSITY_DIMS,
            [2e8, 3e8],
        ),
        (DENSITY_PROPERTY, Density::DIMS, [1000.0, 800.0]),
    ] {
        let mut claim = scalar(name, dims, values[0]);
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, values[0]), (400.0, values[1])],
            dims,
        };
        claim.validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            400.0,
        );
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.push(claim);
    }
    claims
}

fn dynamic(
    source: &MaterialCard,
    reference: f64,
    current: f64,
    maximum_rate: f64,
) -> Result<ResolvedThermalBrittle, BrittleSourceError> {
    resolve_thermal_viscous_brittle_dynamic_material(
        source,
        &query(reference),
        &query(current),
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
        maximum_rate,
    )
}

#[test]
fn thermal_viscous_resolution_retains_every_original_claim_and_reference() {
    let source = card(claims());
    let start = dynamic(&source, 300.0, 300.0, 1.0).unwrap();
    let hot = start.at_temperature(&source, 400.0).unwrap();
    assert_eq!(hot, dynamic(&source, 300.0, 400.0, 1.0).unwrap());
    assert_eq!(thermal_viscous_brittle_requirements().unwrap().len(), 10);
    assert_eq!(
        thermal_viscous_brittle_dynamic_requirements()
            .unwrap()
            .len(),
        11
    );
    assert_eq!(start.band().material().properties().len(), 11);
    assert_eq!(hot.band().material().card_identity(), source.content_hash());
    near(hot.free_linear_strain(), 0.002);
    near(start.band().reference_density_kg_m3().unwrap(), 1000.0);
    near(hot.band().reference_density_kg_m3().unwrap(), 800.0);
    let viscosity = hot.band().law().viscosity().unwrap();
    assert_eq!(viscosity.normal_pa_s_m(), 2e8);
    assert_eq!(viscosity.tangential_pa_s_m(), 3e8);
    assert_eq!(viscosity.maximum_rate_m_s(), 1.0);
    assert_eq!(
        start.expansion().reference().resolved().identity(),
        hot.expansion().reference().resolved().identity()
    );
    assert_eq!(
        start.expansion().selected_claim(),
        hot.expansion().selected_claim()
    );
    for (before, after) in start
        .band()
        .material()
        .properties()
        .iter()
        .zip(hot.band().material().properties())
    {
        assert_eq!(before.requirement(), after.requirement());
        assert_eq!(
            before.answer().receipt.selected,
            after.answer().receipt.selected
        );
    }
    assert!(hot.at_temperature(&source, 400.1).is_err());
    // Additional claims in the card do not silently select a physical loss law.
    let intrinsic = resolve_thermal_brittle(
        &source,
        &query(300.0),
        &query(400.0),
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
    )
    .unwrap();
    assert!(intrinsic.band().law().viscosity().is_none());
    let state = hot.band().law().initial_state();
    let free = hot.free_linear_strain();
    // Removing the blanket source rejection does not enable an untimed update.
    assert!(
        hot.trial_total(&[free, free, free, 0.0, 0.0, 0.0], &state)
            .is_err()
    );
    assert_eq!(state.strain(), &[0.0; 6]);
}

#[test]
fn actual_temperature_rate_is_distinct_from_the_integrated_interval_rate() {
    let source = card(claims());
    let start = dynamic(&source, 300.0, 300.0, 1.0).unwrap();
    let hot = start.at_temperature(&source, 400.0).unwrap();
    let dt = 50.0;
    let temperature_rate = (400.0 - 300.0) / dt;
    let interval_rate = (hot.free_linear_strain() - start.free_linear_strain()) / dt;
    let endpoint_rate = hot
        .free_strain_rate_from_temperature_rate(temperature_rate)
        .unwrap();
    near(interval_rate, 4e-5);
    near(endpoint_rate, 6e-5);
    assert!(endpoint_rate > interval_rate);

    let warm_reference = dynamic(&source, 400.0, 400.0, 1.0).unwrap();
    let cold = warm_reference.at_temperature(&source, 300.0).unwrap();
    near(cold.free_linear_strain(), -0.002);
    near(
        cold.free_strain_rate_from_temperature_rate(-temperature_rate)
            .unwrap(),
        -2e-5,
    );
    near(
        cold.free_strain_rate_from_temperature_rate(0.0).unwrap(),
        0.0,
    );
    for unsupported in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::from_bits(1),
    ] {
        assert!(
            cold.free_strain_rate_from_temperature_rate(unsupported)
                .is_err()
        );
    }
    // Rate observation does not change any source state or integrate temperature.
    assert_eq!(cold, dynamic(&source, 400.0, 300.0, 1.0).unwrap());
}

#[test]
fn source_temperature_rate_preserves_zero_and_negative_expansion() {
    for alpha in [0.0, -1e-5] {
        let mut claims = claims();
        let coefficient = claims
            .iter_mut()
            .find(|claim| claim.key.name() == LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .unwrap();
        *coefficient = scalar(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            alpha,
        );
        let source = card(claims);
        let resolved = dynamic(&source, 300.0, 350.0, 1.0).unwrap();
        near(
            resolved
                .free_strain_rate_from_temperature_rate(2.0)
                .unwrap(),
            2.0 * alpha,
        );
        near(
            resolved
                .free_strain_rate_from_temperature_rate(-2.0)
                .unwrap(),
            -2.0 * alpha,
        );
    }
}

#[test]
fn explicit_thermal_rate_requirements_refuse_missing_or_misidentified_physics() {
    for name in [
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        BRITTLE_NORMAL_VISCOSITY_PROPERTY,
        BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
    ] {
        let mut incomplete = claims();
        incomplete.retain(|claim| claim.key.name() != name);
        assert!(dynamic(&card(incomplete), 300.0, 300.0, 1.0).is_err());
    }
    let mut wrong = claims();
    let eta = wrong
        .iter_mut()
        .find(|claim| claim.key.name() == BRITTLE_NORMAL_VISCOSITY_PROPERTY)
        .unwrap();
    // Bulk dynamic viscosity in Pa s cannot substitute for plane viscosity.
    *eta = scalar(
        BRITTLE_NORMAL_VISCOSITY_PROPERTY,
        Dims([-1, 1, -1, 0, 0, 0]),
        1e8,
    );
    assert!(dynamic(&card(wrong), 300.0, 300.0, 1.0).is_err());

    let mut without_density = claims();
    without_density.retain(|claim| claim.key.name() != DENSITY_PROPERTY);
    let source = card(without_density);
    assert!(dynamic(&source, 300.0, 300.0, 1.0).is_err());
    assert!(
        resolve_thermal_viscous_brittle(
            &source,
            &query(300.0),
            &query(400.0),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
            1.0
        )
        .is_ok()
    );
    let source = card(claims());
    for maximum in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(dynamic(&source, 300.0, 300.0, maximum).is_err());
    }
    let changed_context = query(400.0).with("pressure", 1e5).unwrap();
    assert!(
        resolve_thermal_viscous_brittle_dynamic_material(
            &source,
            &query(300.0),
            &changed_context,
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
            1.0
        )
        .is_err()
    );
}

#[test]
fn thermal_rate_successors_retain_irreversible_history_without_dropping_the_rate_owner() {
    let source = card(claims());
    let start = dynamic(&source, 300.0, 300.0, 1.0).unwrap();
    let hot = start.at_temperature(&source, 400.0).unwrap();
    let virgin = start.band().law().initial_state();
    let loaded = start
        .band()
        .law()
        .trial_over_time(&[0.003, 0.0, 0.0, 0.0, 0.0, 0.0], &virgin, 0.1)
        .unwrap();
    let state = loaded.equilibrium.state;
    let saved = state;
    assert!(state.damage() > 0.0 && !state.failed());
    for invalid in [
        dynamic(&source, 350.0, 400.0, 1.0).unwrap(),
        dynamic(&source, 300.0, 400.0, 0.5).unwrap(),
        resolve_thermal_brittle(
            &source,
            &query(300.0),
            &query(400.0),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        )
        .unwrap(),
    ] {
        assert!(start.transfer_coefficients(&state, &invalid, 1e-6).is_err());
        assert_eq!(state, saved);
    }
    let prepared = start.transfer_coefficients(&state, &hot, 1e-6).unwrap();
    assert_eq!(prepared.state.strain(), state.strain());
    assert_eq!(prepared.state.damage(), state.damage());
    assert_eq!(prepared.state.dissipated_j_m3(), state.dissipated_j_m3());
    // Only eta and the free-strain snapshot changed; viscosity stores no energy.
    near(prepared.parameter_work_j_m3, 0.0);
    near(prepared.damage_dissipated_j_m3, 0.0);
    assert_eq!(state, saved);
}
