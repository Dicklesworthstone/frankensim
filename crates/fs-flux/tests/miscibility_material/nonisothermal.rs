//! Actual local material temperatures, original calorics and atomic source gates.
//! These source tables are numerical fixtures, not measured liquid properties.

use super::*;
use fs_flux::miscibility::adaptive::BinaryMiscibilityTimeControl;
use fs_flux::miscibility::material::nonisothermal::BinaryMiscibilityMaterialNonisothermalError as Error;
use fs_flux::miscibility::nonisothermal::{
    BinaryNonisothermalFace, BinaryNonisothermalSolve,
    adaptive::BinaryNonisothermalTimeControl,
};
use fs_material::liquid::solution::caloric::BinaryLiquidCaloricLaw;
use fs_material::thermal::transport::isobaric::{IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery};

fn mixture(card: &MaterialCard, interaction_k: f64) -> BinaryLiquidSolutionLaw<'_> {
    let point = point(300.0);
    let reference = resolve_liquid_state_from_kinematic(
        card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    BinaryLiquidSolutionLaw::try_new(
        pair(interaction_k),
        [MolarMass::new(0.018), MolarMass::new(0.10)],
        card,
        reference,
        "T",
        "x_a",
    )
    .unwrap()
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("x_a", composition(), 0.4)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
}

fn pure_cp(species: &str, molar_mass: f64, molar_cp: f64) -> IsobaricHeatCapacityLaw {
    let key = PropertyKey::new("pure_cp", SPECIFIC_HEAT_CAPACITY_DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (290.0, molar_cp / molar_mass),
                    (600.0, molar_cp / molar_mass),
                ],
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    quantity(QuantityKind::AbsoluteTemperature),
                    290.0,
                    600.0,
                )
                .with_quantity(
                    "pressure",
                    quantity(QuantityKind::Pressure),
                    100_000.0,
                    100_000.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: format!("synthetic pure {species}"),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: species.into(),
            phase: "liquid".into(),
            process: "constant caloric fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    IsobaricHeatCapacityLaw::try_new(
        &card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &key,
            temperature_axis: "T",
            pressure_axis: "pressure",
            lower: &point(290.0),
            upper: &point(600.0),
            selection: fs_matdb::ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
        },
        300.0,
    )
    .unwrap()
}

fn caloric<'liquid, 'card>(
    liquid: &'liquid BinaryLiquidSolutionLaw<'card>,
) -> BinaryLiquidCaloricLaw<'liquid, 'card> {
    BinaryLiquidCaloricLaw::try_new(
        liquid,
        [
            (SpeciesId::new("A").unwrap(), pure_cp("A", 0.018, 40.0)),
            (SpeciesId::new("B").unwrap(), pure_cp("B", 0.10, 80.0)),
        ],
        "pressure",
    )
    .unwrap()
}
fn faces() -> Vec<BinaryNonisothermalFace> {
    vec![BinaryNonisothermalFace {
        thermal_conductivity_w_m_k: 0.03,
        heat_of_transport_j_per_mol: 0.0,
    }]
}
fn bracket() -> [Temperature; 2] {
    [Temperature::new(290.0), Temperature::new(600.0)]
}
fn time() -> BinaryNonisothermalTimeControl {
    BinaryNonisothermalTimeControl {
        composition: BinaryMiscibilityTimeControl {
            initial_step_s: 0.05,
            minimum_step_s: 1e-7,
            maximum_step_s: 0.1,
            maximum_attempts: 1000,
            absolute_mole_fraction_tolerance: 1e-5,
            relative_mole_fraction_tolerance: 1e-4,
        },
        absolute_temperature_tolerance_k: 1e-4,
        relative_temperature_tolerance: 1e-6,
    }
}

#[test]
fn original_sources_follow_local_temperature_composition_and_adaptive_evolution() {
    let card = temperature_card(1.0, 600.0);
    let liquid = mixture(&card, 0.0);
    let caloric = caloric(&liquid);
    let model = gradient(300.0, 10_000.0, 0.0, 0.0);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model.clone(), 1e-12).unwrap();
    let net = network(model, 1e-5);
    let initial = net
        .state_from_mole_fractions(vec![[0.25, 0.75], [0.75, 0.25]], || false)
        .unwrap();
    let resolved = net
        .resolve_liquid_field(&initial, &source, 2, || false)
        .unwrap();
    let local = resolved
        .bind_nonisothermal(
            &caloric,
            vec![Temperature::new(320.0), Temperature::new(360.0)],
            faces(),
            bracket(),
            1 << 20,
            || false,
        )
        .unwrap();
    let flow = net
        .prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: 1e-5,
                reverse_m3_s: 1e-5,
            }],
            0.0,
            || false,
        )
        .unwrap();
    let advected = local
        .advect_candidate(
            &flow,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            || false,
        )
        .unwrap();
    // One percent of each original parcel trades places. Recover temperature
    // from its carried enthalpy and actual species Cp, not an average of T.
    near(advected.field().field().component_moles()[0][0], 0.255);
    near(
        advected.field().field().temperatures()[0].value(),
        300.0 + 1416.0 / 69.8,
    );
    near(
        advected.field().field().temperatures()[1].value(),
        300.0 + 2984.0 / 50.2,
    );
    let transported = local
        .transport_candidate(
            &flow,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            || false,
        )
        .unwrap();
    assert_eq!(transported.field().field().elapsed_s(), 0.1);
    assert!(std::ptr::eq(transported.field().caloric_law(), &caloric));
    for (i, cell) in transported.field().cells().iter().enumerate() {
        assert_eq!(
            cell.thermodynamics().temperature(),
            transported.field().field().temperatures()[i]
        );
    }
    let mut polls = 0;
    let candidate = local
        .step_candidate(0.1, BinaryNonisothermalSolve::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let end = candidate.field();
    assert_eq!(end.field().elapsed_s(), 0.1);
    assert!(std::ptr::eq(end.caloric_law(), &caloric));
    assert!(std::ptr::eq(end.source(), &source));
    assert_eq!(end.maximum_resolved_cells(), 2);
    assert_eq!(end.field().transport_faces(), faces());
    let mut enthalpy = 0.0;
    let mut totals = [0.0; 2];
    for (i, cell) in end.cells().iter().enumerate() {
        let t = end.field().temperatures()[i].value();
        let n = end.field().component_moles()[i];
        let x = n[0] / (n[0] + n[1]);
        near(cell.phase().volume_m3(), 1e-4);
        assert_eq!(cell.phase().component_moles(), n);
        let transport = cell.phase().transport_state().unwrap();
        near(transport.density_kg_m3(), 1000.0 - 820.0 * x);
        near(
            transport.kinematic_viscosity_m2_s(),
            6e-6 - 5e-6 * (t - 300.0) / 300.0,
        );
        near(
            transport.specific_heat_j_kg_k(),
            3500.0 + 500.0 * (t - 300.0) / 300.0,
        );
        near(
            transport.thermal_conductivity_w_m_k(),
            0.8 - 0.4 * (t - 300.0) / 300.0,
        );
        let material = cell.phase().liquid().material();
        assert_eq!(material.state_coordinate("T"), Some(t));
        assert_eq!(material.state_coordinate("pressure"), Some(100_000.0));
        assert_eq!(
            material.card_identity(),
            liquid.reference_state().material().card_identity()
        );
        enthalpy += (40.0 * n[0] + 80.0 * n[1]) * (t - 300.0);
        for s in 0..2 {
            totals[s] += n[s];
        }
    }
    near(totals[0], 1.0);
    near(totals[1], 1.0);
    assert!((enthalpy - (70.0 * 20.0 + 50.0 * 60.0)).abs() < 1e-8);
    for stop in [1, polls] {
        let mut calls = 0;
        assert!(
            local
                .step_candidate(0.1, BinaryNonisothermalSolve::default(), || {
                    calls += 1;
                    calls == stop
                })
                .is_err()
        );
        assert_eq!(local.field().component_moles(), initial.component_moles());
        assert_eq!(local.field().elapsed_s(), 0.0);
    }
    let adaptive = end
        .advance_candidate(0.2, time(), BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    assert_eq!(adaptive.field().field().elapsed_s(), 0.1 + 0.2);
    assert!(adaptive.advance().accepted_intervals > 0);
    assert!(std::ptr::eq(adaptive.field().caloric_law(), &caloric));
    for (i, cell) in adaptive.field().cells().iter().enumerate() {
        assert_eq!(
            cell.thermodynamics().temperature(),
            adaptive.field().field().temperatures()[i]
        );
    }
    let repeated = end
        .advance_candidate(0.2, time(), BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    assert_eq!(adaptive.field().field(), repeated.field().field());
    assert_eq!(adaptive.field().cells(), repeated.field().cells());
    assert_eq!(adaptive.advance(), repeated.advance());
    let coupled = end
        .advance_transport_candidate(
            &flow,
            0.2,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            time(),
            || false,
        )
        .unwrap();
    assert_eq!(coupled.field().field().elapsed_s(), 0.1 + 0.2);
    assert!(std::ptr::eq(coupled.field().caloric_law(), &caloric));
    for (i, cell) in coupled.field().cells().iter().enumerate() {
        assert_eq!(
            cell.thermodynamics().temperature(),
            coupled.field().field().temperatures()[i]
        );
    }
}

#[test]
fn exact_original_caloric_liquid_binding_is_required() {
    let card = temperature_card(1.0, 600.0);
    let liquid = mixture(&card, 0.0);
    let other = mixture(&card, 0.0);
    let other_caloric = caloric(&other);
    let model = gradient(300.0, 10_000.0, 0.0, 0.0);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model.clone(), 1e-12).unwrap();
    let net = network(model, 1e-5);
    let initial = net
        .state_from_mole_fractions(vec![[0.5, 0.5]; 2], || false)
        .unwrap();
    let resolved = net
        .resolve_liquid_field(&initial, &source, 2, || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        resolved.bind_nonisothermal(
            &other_caloric,
            vec![Temperature::new(320.0); 2],
            faces(),
            bracket(),
            1 << 20,
            || {
                polls += 1;
                false
            },
        ),
        Err(Error::Invalid(_))
    ));
    assert_eq!(polls, 1);
    assert_eq!(resolved.state(), &initial);
}

#[test]
fn numerical_mixing_heat_cannot_publish_unsupported_density_or_temperature() {
    for (density_scale, maximum_temperature) in [(0.9, 600.0), (1.0, 300.000001)] {
        let card = temperature_card(density_scale, maximum_temperature);
        let liquid = mixture(&card, -500.0);
        let caloric = caloric(&liquid);
        let model = gradient(300.0, 10_000.0, 0.0, -500.0);
        let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model.clone(), 1e-12).unwrap();
        let net = network(model, 1e-5);
        let initial = net
            .state_from_mole_fractions(vec![[0.25, 0.75], [0.75, 0.25]], || false)
            .unwrap();
        let resolved = net
            .resolve_liquid_field(&initial, &source, 2, || false)
            .unwrap();
        let local = resolved
            .bind_nonisothermal(
                &caloric,
                vec![Temperature::new(300.0); 2],
                faces(),
                bracket(),
                1 << 20,
                || false,
            )
            .unwrap();
        let flow = net
            .prepare_advection(
                vec![BinaryMiscibilityAdvectionFace {
                    forward_m3_s: 1e-5,
                    reverse_m3_s: 1e-5,
                }],
                0.0,
                || false,
            )
            .unwrap();
        let numerical = local
            .field()
            .step_candidate(0.1, BinaryNonisothermalSolve::default(), || false)
            .unwrap();
        assert!(
            numerical
                .field()
                .temperatures()
                .iter()
                .any(|t| t.value() > 300.000001)
        );
        assert!(matches!(
            local.step_candidate(0.1, BinaryNonisothermalSolve::default(), || false),
            Err(Error::Material(
                BinaryMiscibilityMaterialError::Source { .. }
            ))
        ));
        assert!(matches!(
            local.advance_candidate(0.1, time(), BinaryNonisothermalSolve::default(), || false),
            Err(Error::Material(
                BinaryMiscibilityMaterialError::Source { .. }
            ))
        ));
        assert!(matches!(
            local.advect_candidate(
                &flow,
                0.1,
                BinaryMiscibilityAdvectionControls::default(),
                BinaryNonisothermalSolve::default(),
                || false,
            ),
            Err(Error::Material(
                BinaryMiscibilityMaterialError::Source { .. }
            ))
        ));
        assert!(matches!(
            local.transport_candidate(
                &flow,
                0.1,
                BinaryMiscibilityAdvectionControls::default(),
                BinaryNonisothermalSolve::default(),
                || false,
            ),
            Err(Error::Material(
                BinaryMiscibilityMaterialError::Source { .. }
            ))
        ));
        assert!(matches!(
            local.advance_transport_candidate(
                &flow,
                0.1,
                BinaryMiscibilityAdvectionControls::default(),
                BinaryNonisothermalSolve::default(),
                time(),
                || false,
            ),
            Err(Error::Material(
                BinaryMiscibilityMaterialError::Source { .. }
            ))
        ));
        assert_eq!(local.field().elapsed_s(), 0.0);
        assert_eq!(local.field().component_moles(), initial.component_moles());
        assert_eq!(local.cells(), resolved.cells());
    }
}
