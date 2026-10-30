//! Native finite-fluid convection: exact linear limit and material substitution.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::links::observer::SourceLinkConvectionObserver;
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};

fn constant_fluid_source(
    cp: f64,
    transport_lower_k: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = fluid_source(1.0);
    let mut claims = ClaimSet::new();
    for (_, value) in original.claims().claims_ordered() {
        let mut claim = value.clone();
        let value = match claim.key.name() {
            SPECIFIC_HEAT_CAPACITY_PROPERTY => Some(cp),
            THERMAL_CONDUCTIVITY_PROPERTY => Some(0.5),
            DYNAMIC_VISCOSITY_PROPERTY => Some(0.001),
            _ => None,
        };
        if let Some(value) = value {
            claim.value = PropertyValue::Scalar {
                value,
                dims: claim.key.dims(),
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        if [THERMAL_CONDUCTIVITY_PROPERTY, DYNAMIC_VISCOSITY_PROPERTY].contains(&claim.key.name()) {
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("T", temperature_quantity(), transport_lower_k, 320.0)
                .with("p", 100000.0, 100000.0)
                .with("solute_mass_fraction", 0.2, 0.2);
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &context(280.0),
            upper: &context(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        -20.0 * cp,
        1e-8,
    )
    .unwrap();
    (card, caloric)
}
fn controls(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-8,
    }
}
fn conductance(cp: f64) -> f64 {
    LINK_AREA * 0.664 * (1000.0_f64 * 0.1 / 0.001).sqrt() * (0.001 * cp / 0.5).powf(1.0 / 3.0) * 0.5
        / 0.1
}
fn exact(cp: f64, t: f64) -> [f64; 2] {
    let wall_capacity = MASS * CP;
    let fluid_capacity = MASS * cp;
    let total = wall_capacity + fluid_capacity;
    let mean = (wall_capacity * 310.0 + fluid_capacity * 290.0) / total;
    let difference =
        20.0 * (-conductance(cp) * (1.0 / wall_capacity + 1.0 / fluid_capacity) * t).exp();
    [
        mean + fluid_capacity / total * difference,
        mean - wall_capacity / total * difference,
    ]
}

#[test]
fn g1_g3_one_heat_capacity_source_drives_both_transport_and_storage_on_native_steps() {
    let wall_chart = curve();
    let mut substitutions = Vec::new();
    for cp in [4000.0, 8000.0] {
        let (card, caloric) = constant_fluid_source(cp, 280.0);
        let law = liquid_law(&card);
        let thermal = pair(&wall_chart, &caloric, false, 0.0);
        let saved = thermal.cells().clone();
        let bindings = [link()];
        let observer = SourceLinkConvectionObserver::new(
            &bindings,
            1,
            |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
                law.at_temperature(t)
            },
        );
        let reference = exact(cp, 0.2);
        let mut errors = Vec::new();
        for dt in [0.02, 0.01, 0.005] {
            let (next, accepted, audit) = thermal
                .advance_observed(&observer, 0.2, controls(dt), || Ok(()))
                .unwrap();
            let actual = IDS.map(|id| next.temperature_k(id).unwrap());
            errors.push((actual[0] - reference[0]).abs());
            near(next.energy_j(), thermal.energy_j(), 1e-8);
            near(audit.external_heat_j, 0.0, 1e-12);
            let sample = &accepted.sample().unwrap().samples()[0];
            near(sample.wall_temperature_k, actual[0], 1e-12);
            near(sample.fluid_temperature_k, actual[1], 1e-12);
            near(sample.thermal.transport().specific_heat_j_kg_k(), cp, 1e-12);
            near(sample.conductance_w_k, conductance(cp), 1e-12);
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
            let applied = &accepted.applied_sample().unwrap().samples()[0];
            assert!(applied.wall_temperature_k > sample.wall_temperature_k);
            assert!(applied.fluid_temperature_k < sample.fluid_temperature_k);
            assert_eq!(next.links()[&5].conductance_w_k, 0.0);
        }
        assert!(errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1]);
        assert_eq!(thermal.cells(), &saved);
        assert!(observer.sample().is_none() && observer.applied_sample().is_none());
        // Its hA belongs in the ORIGINAL conductive CFL. A huge requested dt
        // therefore needs multiple native substeps and cannot fit a one-step budget.
        assert!(matches!(
            thermal.advance_observed_limited(&observer, 0.3, controls(1.0), 1, || Ok(()),),
            Err(SourceCellConvectionError::Convection(
                ConvectionSourceError::Thermal(DepositionError::Budget(_))
            ))
        ));
        substitutions.push(reference);
    }
    // Same initial temperatures, geometry and solver. Changing Cp data changes
    // the thermal capacity AND the convection correlation's Prandtl number.
    assert!(substitutions[1][0] < substitutions[0][0]);
    assert!(substitutions[1][1] < substitutions[0][1]);
}

#[test]
fn g0_late_material_refusal_cancellation_and_competing_link_owners_are_atomic() {
    let chart = curve();
    let (limited_card, limited_caloric) = constant_fluid_source(4000.0, 299.0);
    let limited = liquid_law(&limited_card);
    let thermal = pair(&chart, &limited_caloric, false, 0.0);
    let saved = thermal.cells().clone();
    let bindings = [link()];
    let observer = SourceLinkConvectionObserver::new(
        &bindings,
        1,
        |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
            limited.at_temperature(t)
        },
    );
    // The caloric law supports both cells; only the new film endpoint leaves
    // transport support. Neither finite body's proposed energy may be published.
    assert!(matches!(thermal.advance_observed(
        &observer, 0.04, controls(0.04), || Ok(()),
    ), Err(SourceCellConvectionError::Source { property_temperature_k, .. })
        if property_temperature_k < 299.0));
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none() && observer.applied_sample().is_none());

    let (card, caloric) = constant_fluid_source(4000.0, 280.0);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let saved = thermal.cells().clone();
    let observer = SourceLinkConvectionObserver::new(
        &bindings,
        1,
        |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
            law.at_temperature(t)
        },
    );
    let duplicated = CoupledObservers::new(observer.clone(), observer.clone());
    assert!(matches!(
        thermal.advance_observed(&duplicated, 0.04, controls(0.02), || Ok(()),),
        Err(ObserverCouplingError::ConflictingConductance(5))
    ));
    let mut polls = 0;
    thermal
        .advance_observed(&observer, 0.04, controls(0.02), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut at = 0;
    assert!(matches!(
        thermal.advance_observed(&observer, 0.04, controls(0.02), || {
            at += 1;
            if at == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Cancelled
        ))
    ));
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none() && observer.applied_sample().is_none());
}
