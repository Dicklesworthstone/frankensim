//! Density-table substitutions through the real convection and enthalpy owners.

use super::*;
use fs_airflow::graph::thermal::medium::convection::evaluate_vertical_plate_from_thermal_material;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    UncertaintyModel,
};
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::state_point::*;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec};

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

// Every material has exactly the SAME transport at 300 K. Only the sourced
// density slope differs, isolating the expansion contribution to convection.
fn liquid_card(density_drop_per_k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, scalar) in [
        (DENSITY_PROPERTY, Density::DIMS, 0.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            4000.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.5,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 0.001),
        ("sound_speed", fs_qty::Velocity::DIMS, 1500.0),
    ] {
        let density = name == DENSITY_PROPERTY;
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: if density {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: [280.0, 320.0]
                            .map(|t| (t, 1000.0 - density_drop_per_k * (t - 300.0)))
                            .to_vec(),
                        dims,
                    }
                } else {
                    PropertyValue::Scalar {
                        value: scalar,
                        dims,
                    }
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity(),
                    280.0,
                    320.0,
                ),
                interpolation: if density {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic liquid density law; not observations".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "test liquid".into(),
            phase: "liquid".into(),
            process: "source-substitution fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn liquid_law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    let point = QueryPoint::new()
        .with_quantity("T", temperature_quantity(), 300.0)
        .unwrap();
    LiquidTemperatureLaw::try_new(
        card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn acoustic_liquid_binding_reaches_actual_convection_without_losing_its_density_law() {
    use fs_couple::vibroacoustic::{AcousticMedium, rectangular_cavity_modes};
    use fs_material::fluid_acoustics::{
        FluidAcousticParameter, FluidAcoustics, ResolvedAcousticLiquidState,
        resolve_acoustic_liquid_state,
    };

    let chart = curve();
    let heat = initial(&chart, [310.0, 290.0], 0.0, 0.0);
    let mut endpoints = Vec::new();
    for slope in [0.2, 1.0] {
        let card = liquid_card(slope);
        let point = QueryPoint::new()
            .with_quantity("T", temperature_quantity(), 300.0)
            .unwrap();
        let parameter = FluidAcousticParameter::sound_speed().unwrap();
        let acoustic = resolve_acoustic_liquid_state(
            &card,
            &point,
            &parameter,
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let liquid = acoustic
            .liquid()
            .clone()
            .with_density_expansion(&card, "T", Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let medium = ResolvedAcousticLiquidState::try_from_liquid(liquid, &parameter).unwrap();
        let convection = evaluate_vertical_plate_from_thermal_material(
            &medium,
            9.81,
            Temperature::new(310.0),
            Length::new(0.1),
        )
        .unwrap();
        near(
            convection.coefficient().value.value(),
            expected_coefficient(slope, 310.0),
            1e-8,
        );
        assert!(std::ptr::eq(convection.fluid(), &medium));
        let cavity = rectangular_cavity_modes(
            0.2,
            0.1,
            0.1,
            AcousticMedium::try_from_material(&medium).unwrap(),
            0.0,
            2,
            &[[0.0, 0.0]],
        )
        .unwrap();
        near(cavity.omegas[1], core::f64::consts::PI * 1500.0 / 0.2, 1e-9);
        assert_eq!(
            cavity.rho0,
            medium.acoustic_state().unwrap().density_kg_m3()
        );

        let dt = 0.1;
        let step = BackwardEulerStep::new(&heat, dt, || Ok(())).unwrap();
        let increments = step
            .cell_ids()
            .iter()
            .map(|&id| {
                let old = heat.temperature_k(id).unwrap();
                MASS * CP * (expected_endpoint(slope, old, dt) - old)
            })
            .collect::<Vec<_>>();
        let trial = step
            .evaluate(
                &increments,
                |endpoint| {
                    let mut sources = BTreeMap::new();
                    for id in IDS {
                        let wall = endpoint.temperature_k(id).unwrap();
                        let convection = evaluate_vertical_plate_from_thermal_material(
                            &medium,
                            9.81,
                            Temperature::new(wall),
                            Length::new(0.1),
                        )
                        .map_err(ConvectionSourceError::Material)?;
                        sources.insert(
                            id,
                            endpoint.exposed_area_m2(id).unwrap()
                                * convection.coefficient().value.value()
                                * (300.0 - wall),
                        );
                    }
                    Ok::<_, ConvectionSourceError>(StepCoefficients {
                        sources_w: sources,
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    })
                },
                || Ok(()),
            )
            .unwrap();
        let (next, (), audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
        near(
            next.energy_j() - heat.energy_j(),
            audit.external_heat_j,
            1e-8,
        );
        near(
            next.temperature_k(IDS[0]).unwrap(),
            expected_endpoint(slope, 310.0, dt),
            1e-10,
        );
        endpoints.push(next.temperature_k(IDS[0]).unwrap());
    }
    assert!(endpoints[1] < endpoints[0] && endpoints[1] > 300.0);
}

fn expected_coefficient(slope: f64, wall: f64) -> f64 {
    let rho = 1000.0;
    let fluid = FluidTransportState::try_new(rho, 4000.0, 0.5, 0.001).unwrap();
    natural_h(fluid, wall - 300.0, slope / rho)
}

fn expected_endpoint(slope: f64, old: f64, dt: f64) -> f64 {
    let (mut lower, mut upper) = (old.min(300.0), old.max(300.0));
    for _ in 0..70 {
        let t = 0.5 * lower + 0.5 * upper;
        let power = AREA * expected_coefficient(slope, t) * (300.0 - t);
        let residual = MASS * CP * (t - old) - dt * power;
        if residual > 0.0 {
            upper = t;
        } else {
            lower = t;
        }
    }
    0.5 * lower + 0.5 * upper
}

#[test]
fn g1_g3_density_law_alone_changes_convection_and_accepted_thermal_history() {
    let chart = curve();
    let heat = initial(&chart, [310.0, 290.0], 0.0, 0.0);
    let dt = 0.1;
    let mut cooled = Vec::new();
    for slope in [0.2, 1.0] {
        let card = liquid_card(slope);
        let law = liquid_law(&card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let medium = law.at_temperature(Temperature::new(300.0)).unwrap();
        let convection = evaluate_vertical_plate_from_thermal_material(
            &medium,
            9.81,
            Temperature::new(310.0),
            Length::new(0.1),
        )
        .unwrap();
        near(
            convection.coefficient().value.value(),
            expected_coefficient(slope, 310.0),
            1e-8,
        );
        assert_eq!(convection.transport().density_kg_m3(), 1000.0);
        let boundary = convection.robin_boundary(Temperature::new(300.0)).unwrap();
        assert_eq!(boundary.coefficient().value, convection.coefficient().value);

        let step = BackwardEulerStep::new(&heat, dt, || Ok(())).unwrap();
        let increments = step
            .cell_ids()
            .iter()
            .map(|&id| {
                let old = heat.temperature_k(id).unwrap();
                MASS * CP * (expected_endpoint(slope, old, dt) - old)
            })
            .collect::<Vec<_>>();
        let trial = step
            .evaluate(
                &increments,
                |endpoint| {
                    let mut sources = BTreeMap::new();
                    for id in IDS {
                        let wall = endpoint.temperature_k(id).unwrap();
                        let convection = evaluate_vertical_plate_from_thermal_material(
                            &medium,
                            9.81,
                            Temperature::new(wall),
                            Length::new(0.1),
                        )
                        .map_err(ConvectionSourceError::Material)?;
                        near(
                            convection.coefficient().value.value(),
                            expected_coefficient(slope, wall),
                            1e-8,
                        );
                        sources.insert(
                            id,
                            endpoint.exposed_area_m2(id).unwrap()
                                * convection.coefficient().value.value()
                                * (300.0 - wall),
                        );
                    }
                    Ok::<_, ConvectionSourceError>(StepCoefficients {
                        sources_w: sources,
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    })
                },
                || Ok(()),
            )
            .unwrap();
        let (next, (), audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
        for id in IDS {
            near(
                next.temperature_k(id).unwrap(),
                expected_endpoint(slope, heat.temperature_k(id).unwrap(), dt),
                1e-10,
            );
        }
        near(
            next.energy_j() - heat.energy_j(),
            audit.external_heat_j,
            1e-8,
        );
        cooled.push(next.temperature_k(IDS[0]).unwrap());
    }
    assert!(cooled[1] < cooled[0] && cooled[1] > 300.0);
    assert_eq!(heat.temperature_k(IDS[0]).unwrap(), 310.0);
}

#[test]
fn g0_unknown_expansion_refuses_and_source_temperature_reaches_the_reservoir() {
    let card = liquid_card(0.5);
    let law = liquid_law(&card);
    let unknown = law.at_temperature(Temperature::new(305.0)).unwrap();
    assert!(
        evaluate_vertical_plate_from_thermal_material(
            &unknown,
            9.81,
            Temperature::new(310.0),
            Length::new(0.1)
        )
        .is_err()
    );
    let sourced = law
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let medium = sourced.at_temperature(Temperature::new(305.0)).unwrap();
    let convection = evaluate_vertical_plate_from_thermal_material(
        &medium,
        9.81,
        Temperature::new(310.0),
        Length::new(0.1),
    )
    .unwrap();
    let thermal = medium.thermal_transport_state().unwrap();
    let ra = thermal
        .transport()
        .rayleigh(
            9.81,
            thermal.volumetric_expansion_per_k().unwrap(),
            5.0,
            0.1,
        )
        .unwrap();
    near(convection.nusselt().groups()["Ra"], ra, 1e-8);
    assert_eq!(convection.transport(), &thermal.transport());
}
