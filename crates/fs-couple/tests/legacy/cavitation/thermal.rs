//! G1/G3/G4 original-source temperature continuation of retained film mass.
//! These manufactured curves test coupling, not measured fluid accuracy.
#[path = "compressible_film.rs"]
mod compressible_film;
#[path = "power_law.rs"]
mod power_law;

use super::{AMBIENT_PA, AREA_M2, CHANNEL_GAP_M, INITIAL_GAP_M, film, near, try_film};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::cavitation::{
    CavitationPropertyRequirements, ResolvedCavitationState, resolve_cavitation_state,
};
use fs_material::liquid::surface::SURFACE_TENSION_DIMS;
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_DIMS};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature, Velocity};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), AMBIENT_PA)
        .unwrap()
}

struct Source {
    card: MaterialCard,
    requirements: CavitationPropertyRequirements,
    pins: Vec<(String, ClaimId)>,
    other_viscosity: ClaimId,
}
impl Source {
    fn resolve(&self, temperature: f64) -> ResolvedCavitationState {
        self.resolve_with(temperature, self.pins.clone())
    }
    fn resolve_with(
        &self,
        temperature: f64,
        pins: Vec<(String, ClaimId)>,
    ) -> ResolvedCavitationState {
        resolve_cavitation_state(
            &self.card,
            &point(temperature),
            &self.requirements,
            "temperature",
            MaterialPropertySelection::PinnedByProperty(pins),
        )
        .unwrap()
    }
}

fn source(basis: LiquidViscosityBasis, constant: bool, revision: u32) -> Source {
    let (viscosity_name, viscosity_dims, viscosity) = match basis {
        LiquidViscosityBasis::Dynamic => {
            ("dynamic_viscosity", DYNAMIC_VISCOSITY_DIMS, [0.001, 0.0005])
        }
        LiquidViscosityBasis::Kinematic => (
            "kinematic_viscosity",
            KINEMATIC_VISCOSITY_DIMS,
            [1e-6, 0.5e-6],
        ),
    };
    let rows = [
        ("density", Density::DIMS, [1000.0, 800.0]),
        (viscosity_name, viscosity_dims, viscosity),
        ("surface_tension", SURFACE_TENSION_DIMS, [0.07, 0.06]),
        (
            "saturation_vapor_pressure",
            Pressure::DIMS,
            [2000.0, 20_000.0],
        ),
        ("sound_speed", Velocity::DIMS, [1500.0, 1400.0]),
        ("saturated_vapor_density", Density::DIMS, [0.02, 0.2]),
        (
            "latent_vaporization_enthalpy",
            Dims([2, 0, -2, 0, 0, 0]),
            [2.5e6, 2.2e6],
        ),
    ];
    let mut claims = ClaimSet::new();
    let mut pins = Vec::new();
    let mut other_viscosity = None;
    for (name, dims, values) in rows {
        let hot = if constant { values[0] } else { values[1] };
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, values[0]), (400.0, hot)],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    300.0,
                    400.0,
                )
                .with_quantity(
                    "pressure",
                    quantity(QuantityKind::Pressure),
                    90_000.0,
                    120_000.0,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured film temperature curve".into(),
                license: "test-only".into(),
                artifact: None,
            },
        };
        pins.push((name.into(), claims.insert_claim(claim.clone()).unwrap()));
        if name == viscosity_name {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                for (_, value) in knots {
                    *value *= 2.0;
                }
            }
            other_viscosity = Some(claims.insert_claim(claim).unwrap());
        }
    }
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured film liquid".into(),
            phase: "liquid with equilibrium vapor data".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    for _ in 0..revision {
        card = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    }
    Source {
        card,
        requirements: CavitationPropertyRequirements::canonical(basis, true).unwrap(),
        pins,
        other_viscosity: other_viscosity.unwrap(),
    }
}

#[test]
fn sourced_heating_and_cooling_conserve_original_mass_and_refresh_dynamic_or_kinematic_viscosity() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let source = source(basis, false, 0);
        let cold = source.resolve(300.0);
        let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
        let original = film(&cold);
        let state = original.initial_state(&[0.0]).unwrap();
        let saved = state.clone();
        let next = original
            .advance_with_temperature_source(&hot, &state, &[0.0], &[0.0], 0.5, || false)
            .unwrap();
        let volume = AREA_M2 * INITIAL_GAP_M;
        let flow = (1000.0 / 800.0 - 1.0) * volume / 0.5;
        let mu = if basis == LiquidViscosityBasis::Dynamic {
            0.0005
        } else {
            800.0 * 0.5e-6
        };
        let conductance = 1e-3 * CHANNEL_GAP_M.powi(3) / (12.0 * mu * 1e-2);
        assert!(std::ptr::eq(next.film.source(), &hot));
        assert_eq!(next.film.source().temperature_axis(), "temperature");
        assert_eq!(next.film.material(), hot.state());
        assert_eq!(next.film.film().geometry().viscosity_pa_s(), Some(mu));
        assert_eq!(next.step.state.liquid_fraction(), &[1.0]);
        near(next.step.channel_volume_flow_m3_s[0], flow);
        near(
            next.step.state.pressure_pa()[0],
            AMBIENT_PA + flow / conductance,
        );
        assert_eq!(next.step.report.liquid_mass_before_kg, 1000.0 * volume);
        near(next.step.report.liquid_mass_after_kg, 800.0 * volume);
        near(next.step.report.boundary_mass_out_kg, 200.0 * volume);
        assert!(next.step.report.mass_residual_kg.abs() < 1e-19);
        assert_eq!(state, saved);
        assert_eq!(original.material(), cold.state());
        assert!(hot.material().property("specific_heat_capacity").is_none());
        // Return to the original source temperature without resetting filling.
        let cooled = next
            .film
            .advance_with_temperature_source(&cold, &next.step.state, &[0.0], &[0.0], 0.5, || false)
            .unwrap();
        assert_eq!(
            cooled.step.report.liquid_mass_before_kg,
            next.step.report.liquid_mass_after_kg
        );
        near(cooled.step.report.liquid_mass_after_kg, 1000.0 * volume);
        near(cooled.step.report.boundary_mass_out_kg, -200.0 * volume);
        assert_eq!(cooled.step.state.liquid_fraction(), &[1.0]);
        assert!(std::ptr::eq(cooled.film.source(), &cold));
    }
}

#[test]
fn source_changed_cavity_uses_new_rupture_pressure_and_viscosity_without_changing_old_mass() {
    let source = source(LiquidViscosityBasis::Dynamic, false, 0);
    let cold = source.resolve(300.0);
    let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
    let original = film(&cold);
    let state = original.initialize(&[0.0], &[0.2]).unwrap();
    let next = original
        .advance_with_temperature_source(&hot, &state, &[0.0], &[0.0], 0.01, || false)
        .unwrap();
    let volume = AREA_M2 * INITIAL_GAP_M;
    let inflow = 1e-3 * CHANNEL_GAP_M.powi(3) * (AMBIENT_PA - 20_000.0) / (12.0 * 0.0005 * 1e-2);
    let old_mass = 1000.0 * (volume * 0.2);
    let next_mass = old_mass + 800.0 * 0.01 * inflow;
    assert_eq!(next.step.state.pressure_pa(), &[20_000.0]);
    assert_eq!(state.pressure_pa(), &[2000.0]);
    near(next.step.channel_volume_flow_m3_s[0], -inflow);
    near(
        next.step.state.liquid_fraction()[0],
        next_mass / (800.0 * volume),
    );
    assert_eq!(next.step.report.liquid_mass_before_kg, old_mass);
    near(next.step.report.liquid_mass_after_kg, next_mass);
    near(
        next.step.generalized_force[0],
        AREA_M2 * (20_000.0 - AMBIENT_PA),
    );
}

#[test]
fn changed_card_claim_or_nonthermal_context_and_dense_vapor_refuse_before_acceptance() {
    let source = source(LiquidViscosityBasis::Dynamic, false, 0);
    let cold = source.resolve(300.0);
    let original = film(&cold);
    let state = original.initialize(&[0.0], &[0.2]).unwrap();
    let saved = state.clone();
    let foreign = self::source(LiquidViscosityBasis::Dynamic, false, 1).resolve(400.0);
    let changed_context = cold
        .at_query_point(
            &source.card,
            &point(400.0).with("pressure", 110_000.0).unwrap(),
        )
        .unwrap();
    let mut pins = source.pins.clone();
    pins.iter_mut()
        .find(|(name, _)| name == "dynamic_viscosity")
        .unwrap()
        .1 = source.other_viscosity;
    let switched_claim = source.resolve_with(400.0, pins);
    for candidate in [&foreign, &changed_context, &switched_claim] {
        assert!(
            original
                .advance_with_temperature_source(candidate, &state, &[0.0], &[0.0], 0.01, || false,)
                .is_err()
        );
    }
    assert!(cold.at_query_point(&source.card, &point(401.0)).is_err());
    let dilute_only = try_film(&cold, 1e-4).unwrap();
    let dilute_state = dilute_only.initialize(&[0.0], &[0.2]).unwrap();
    let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
    assert!(
        dilute_only
            .advance_with_temperature_source(&hot, &dilute_state, &[0.0], &[0.0], 0.01, || false,)
            .is_err()
    );
    assert_eq!(state, saved);
    assert_eq!(original.material(), cold.state());
    assert!(
        original
            .advance(&state, &[0.0], &[0.0], 0.01, || false)
            .is_ok()
    );
}

#[test]
fn unchanged_source_coefficients_keep_native_bits_and_late_cancel_keeps_the_old_pair() {
    let source = source(LiquidViscosityBasis::Dynamic, true, 0);
    let cold = source.resolve(300.0);
    let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
    let original = film(&cold);
    let state = original.initial_state(&[0.0]).unwrap();
    let expected = original
        .advance(&state, &[-INITIAL_GAP_M], &[0.0], 0.01, || false)
        .unwrap();
    let mut polls = 0;
    let next = original
        .advance_with_temperature_source(&hot, &state, &[-INITIAL_GAP_M], &[0.0], 0.01, || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(next.step, expected);
    assert_eq!(next.film.material().temperature_k(), 400.0);
    let saved = state.clone();
    let mut calls = 0;
    assert!(
        original
            .advance_with_temperature_source(&hot, &state, &[-INITIAL_GAP_M], &[0.0], 0.01, || {
                calls += 1;
                calls == polls
            },)
            .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(state, saved);
    assert_eq!(original.material(), cold.state());
    let retry = original
        .advance_with_temperature_source(&hot, &state, &[-INITIAL_GAP_M], &[0.0], 0.01, || false)
        .unwrap();
    assert_eq!(retry.step, expected);
}
