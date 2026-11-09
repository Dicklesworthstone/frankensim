//! Original caloric identities, local-temperature receipts and atomic source gates.
//! Every property and interaction below is synthetic numerical test data.

use super::*;
use fs_flux::miscibility::multicomponent::nonisothermal::{
    MulticomponentNonisothermalError, MulticomponentNonisothermalFace,
    MulticomponentNonisothermalField, MulticomponentNonisothermalSolve,
    adaptive::MulticomponentNonisothermalTimeControl,
    material::MulticomponentMiscibilityMaterialNonisothermalError as ThermalError,
};
use fs_material::liquid::solution::multicomponent::caloric::MulticomponentLiquidCaloricLaw;
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const CP: [f64; 3] = [40.0, 80.0, 120.0];
const CP_SLOPE: [f64; 3] = [0.1, 0.2, 0.3];
const B: [f64; 9] = [0.0, 120.0, -60.0, 120.0, 0.0, 80.0, -60.0, 80.0, 0.0];

fn temperature_source_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, spec, values) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            [44.0, 2.0, 28.0],
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            [3500.0, 4000.0, 0.0],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
            [0.3, 0.6, 0.0],
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            [1e-6, 3e-6, 0.0],
        ),
    ] {
        let value = if name == DENSITY_PROPERTY {
            PropertyValue::Surface2D {
                axes: ["x_a".into(), "x_b".into()],
                axis_dims: [Dimensionless::DIMS; 2],
                samples: [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
                    .into_iter()
                    .zip(values)
                    .collect(),
                triangles: vec![[0, 1, 2]],
                dims: spec.dims(),
            }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, values[0]), (600.0, values[1])],
                dims: spec.dims(),
            }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, spec),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        300.0,
                        600.0,
                    )
                    .with_quantity("x_a", composition(), 0.0, 1.0)
                    .with_quantity("x_b", composition(), 0.0, 1.0)
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
                    source: "synthetic local-temperature ternary liquid".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic A+B+C".into(),
            phase: "liquid".into(),
            process: "original thermal transport fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn pure_cp(index: usize) -> IsobaricHeatCapacityLaw {
    let key = PropertyKey::new("original_pure_cp", SPECIFIC_HEAT_CAPACITY_DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [300.0, 450.0, 600.0]
                    .into_iter()
                    .map(|t| {
                        (
                            t,
                            (CP[index] + CP_SLOPE[index] * (t - 300.0)) / MASSES[index],
                        )
                    })
                    .collect(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    300.0,
                    600.0,
                )
                .with_quantity(
                    "source_pressure",
                    quantity(QuantityKind::Pressure),
                    100_000.0,
                    100_000.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: format!("synthetic original pure component {index}"),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("pure component {index}"),
            phase: "liquid".into(),
            process: "linear molar Cp fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                t,
            )
            .unwrap()
            .with_quantity(
                "source_pressure",
                quantity(QuantityKind::Pressure),
                100_000.0,
            )
            .unwrap()
    };
    IsobaricHeatCapacityLaw::try_new(
        &card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &key,
            temperature_axis: "temperature",
            pressure_axis: "source_pressure",
            lower: &point(300.0),
            upper: &point(600.0),
            selection: fs_matdb::ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
        },
        300.0,
    )
    .unwrap()
}

fn caloric<'liquid, 'card>(
    liquid: &'liquid MulticomponentLiquidSolutionLaw<'card>,
) -> MulticomponentLiquidCaloricLaw<'liquid, 'card> {
    MulticomponentLiquidCaloricLaw::try_new(
        liquid,
        [2, 0, 1]
            .into_iter()
            .map(|i| (species(i), pure_cp(i)))
            .collect(),
        "pressure",
    )
    .unwrap()
}
fn faces() -> Vec<MulticomponentNonisothermalFace> {
    vec![MulticomponentNonisothermalFace {
        thermal_conductivity_w_m_k: 0.03,
        heat_of_transport_j_per_mol: vec![100.0, -60.0, 20.0],
    }]
}
fn bracket() -> [Temperature; 2] {
    [Temperature::new(300.0), Temperature::new(600.0)]
}
fn thermal_time() -> MulticomponentNonisothermalTimeControl {
    MulticomponentNonisothermalTimeControl {
        composition: MulticomponentMiscibilityTimeControl {
            initial_step_s: 0.001,
            minimum_step_s: 1e-8,
            maximum_step_s: 0.01,
            absolute_mole_fraction_tolerance: 1e-5,
            relative_mole_fraction_tolerance: 1e-4,
            maximum_attempts: 1000,
        },
        absolute_temperature_tolerance_k: 1e-4,
        relative_temperature_tolerance: 1e-6,
    }
}

// Independent integral of the authored pure Cp and original B/T interactions.
// No original caloric total or numerical energy receipt enters this oracle.
fn enthalpy(field: &MulticomponentNonisothermalField) -> f64 {
    field
        .temperatures()
        .iter()
        .enumerate()
        .map(|(cell, t)| {
            let n = field.cell_component_moles(cell).unwrap();
            let dt = t.value() - 300.0;
            let total: f64 = n.iter().sum();
            let sensible: f64 = (0..3)
                .map(|i| n[i] * (CP[i] * dt + 0.5 * CP_SLOPE[i] * dt * dt))
                .sum();
            let excess: f64 = (0..3)
                .flat_map(|i| ((i + 1)..3).map(move |j| R * B[3 * i + j] * n[i] * n[j] / total))
                .sum();
            sensible + excess
        })
        .sum()
}

#[test]
fn original_ternary_calorics_follow_actual_local_temperatures_and_replay_atomically() {
    let card = temperature_source_card();
    let model = MulticomponentRegularSolutionV1::new(
        (0..3).map(species).collect(),
        vec![0.0, 0.2, -0.3, 0.2, 0.0, 0.1, -0.3, 0.1, 0.0],
        B.into(),
        Temperature::new(290.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let liquid = source(&card, model.clone(), MASSES);
    let caloric = caloric(&liquid);
    let original = numerical();
    let network = MulticomponentMiscibilityNetwork::new(
        model.clone(),
        original.network().transport().clone(),
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        original.network().cells().to_vec(),
        original.network().faces().to_vec(),
        vec![vec![0.0, 1e-6, 1e-4, 1e-6, 0.0, 1e-4, 1e-4, 1e-4, 0.0]],
        16 * 1024 * 1024,
        || false,
    )
    .unwrap();
    let numerical = network
        .bind(vec![vec![0.3, 0.3, 0.4], vec![0.31, 0.1, 0.59]], || false)
        .unwrap();
    let resolved = numerical
        .resolve_liquid_field(&liquid, 2, 1e-12, || false)
        .unwrap();
    let local = resolved
        .bind_nonisothermal(
            &caloric,
            vec![Temperature::new(320.0), Temperature::new(500.0)],
            faces(),
            bracket(),
            16 * 1024 * 1024,
            || false,
        )
        .unwrap();
    // Equal values from another original liquid instance do not grant identity.
    let other_liquid = source(&card, model, MASSES);
    let other_caloric = self::caloric(&other_liquid);
    assert!(matches!(
        resolved.bind_nonisothermal(
            &other_caloric,
            local.field().temperatures().to_vec(),
            faces(),
            bracket(),
            16 * 1024 * 1024,
            || false,
        ),
        Err(ThermalError::Invalid(_))
    ));
    let before = local.field().clone();
    let before_cells = local.cells().to_vec();
    let controls = MulticomponentNonisothermalSolve::default();
    let mut polls = 0;
    let candidate = local
        .step_candidate(0.001, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let end = candidate.field();
    assert!(std::ptr::eq(end.source(), &liquid));
    assert!(std::ptr::eq(end.caloric_law(), &caloric));
    assert_eq!(end.field().caloric_model(), caloric.model());
    assert_eq!(end.network(), &network);
    assert_eq!(end.maximum_resolved_cells(), 2);
    assert_eq!(end.maximum_relative_volume_mismatch(), 1e-12);
    assert_eq!(end.field().transport_faces(), faces());
    assert_eq!(end.field().elapsed_s(), 0.001);
    assert_ne!(end.field().temperatures(), before.temperatures());
    near(enthalpy(end.field()), enthalpy(&before));
    assert!(candidate.step().entropy_change_j_per_k >= -controls.entropy_tolerance_j_per_k);
    for (cell, source) in end.cells().iter().enumerate() {
        let t = end.field().temperatures()[cell].value();
        let n = end.field().cell_component_moles(cell).unwrap();
        let x = end.field().cell_mole_fractions(cell).unwrap();
        assert_eq!(source.component_moles(), n);
        assert_eq!(source.thermodynamics().temperature(), Temperature::new(t));
        assert_eq!(source.species(), liquid.species());
        near(source.volume_m3(), 0.001);
        let transport = source.transport_state().unwrap();
        near(transport.density_kg_m3(), 44.0 - 42.0 * x[0] - 16.0 * x[1]);
        near(
            transport.specific_heat_j_kg_k(),
            3500.0 + 500.0 * (t - 300.0) / 300.0,
        );
        near(
            transport.thermal_conductivity_w_m_k(),
            0.3 + 0.3 * (t - 300.0) / 300.0,
        );
        near(
            transport.kinematic_viscosity_m2_s(),
            1e-6 + 2e-6 * (t - 300.0) / 300.0,
        );
        let receipt = source.liquid().material();
        assert_eq!(receipt.card_identity(), card.content_hash());
        assert_eq!(receipt.state_coordinate("T"), Some(t));
        assert_eq!(receipt.state_coordinate("pressure"), Some(100_000.0));
        near(receipt.state_coordinate("x_a").unwrap(), x[0]);
        near(receipt.state_coordinate("x_b").unwrap(), x[1]);
    }
    for i in 0..3 {
        near(
            end.field().cell_component_moles(0).unwrap()[i]
                + end.field().cell_component_moles(1).unwrap()[i],
            before.cell_component_moles(0).unwrap()[i] + before.cell_component_moles(1).unwrap()[i],
        );
    }
    // The final poll follows complete numerical and original-source evaluation.
    for stop in [1, polls] {
        let mut calls = 0;
        assert!(matches!(
            local.step_candidate(0.001, controls, || {
                calls += 1;
                calls == stop
            }),
            Err(ThermalError::Numerical(
                MulticomponentNonisothermalError::Cancelled
            ))
        ));
        assert_eq!(local.field(), &before);
        assert_eq!(local.cells(), before_cells);
    }
    let replay = local.step_candidate(0.001, controls, || false).unwrap();
    assert_eq!(replay.field().field(), end.field());
    assert_eq!(replay.field().cells(), end.cells());
    assert_eq!(replay.step(), candidate.step());
    let advanced = local
        .advance_candidate(0.01, thermal_time(), controls, || false)
        .unwrap();
    let raw = before
        .advance_candidate(0.01, thermal_time(), controls, || false)
        .unwrap();
    assert_eq!(advanced.field().field(), raw.field());
    assert_eq!(advanced.advance(), raw.advance());
    assert!(std::ptr::eq(advanced.field().caloric_law(), &caloric));
    near(enthalpy(advanced.field().field()), enthalpy(&before));
    for (i, cell) in advanced.field().cells().iter().enumerate() {
        assert_eq!(
            cell.thermodynamics().temperature(),
            advanced.field().field().temperatures()[i]
        );
    }
}

#[test]
fn completed_thermal_candidates_cannot_escape_original_composition_or_volume_sources() {
    for (upper, perturbation, expect_source) in [(0.31, 0.0, true), (1.0, 1.0, false)] {
        let card = source_card(upper, perturbation);
        let liquid = source(&card, model([2, 0, 1]), MASSES);
        let caloric = caloric(&liquid);
        let numerical = numerical();
        let original = numerical
            .resolve_liquid_field(&liquid, 2, 1e-12, || false)
            .unwrap();
        let local = original
            .bind_nonisothermal(
                &caloric,
                vec![Temperature::new(300.0); 2],
                vec![MulticomponentNonisothermalFace {
                    thermal_conductivity_w_m_k: 0.0,
                    heat_of_transport_j_per_mol: vec![0.0; 3],
                }],
                bracket(),
                16 * 1024 * 1024,
                || false,
            )
            .unwrap();
        let before = local.field().clone();
        let before_cells = local.cells().to_vec();
        let controls = MulticomponentNonisothermalSolve::default();
        let numerical_step = before.step_candidate(0.1, controls, || false).unwrap();
        assert!(numerical_step.field().cell_mole_fractions(1).unwrap()[0] > 0.31);
        let refusal = local.step_candidate(0.1, controls, || false).unwrap_err();
        if expect_source {
            assert!(matches!(
                refusal,
                ThermalError::Material(Error::Source { cell: 1, .. })
            ));
        } else {
            assert!(matches!(
                refusal,
                ThermalError::Material(Error::VolumeMismatch { .. })
            ));
        }
        assert_eq!(local.field(), &before);
        assert_eq!(local.cells(), before_cells);
        assert!(std::ptr::eq(local.source(), &liquid));
        assert!(std::ptr::eq(local.caloric_law(), &caloric));
    }
}
