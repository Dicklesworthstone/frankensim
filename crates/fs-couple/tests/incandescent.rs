//! Manufactured source laws exercise native electrothermal feedback.
//! These fixtures are numerical verification, not measured filament data.

#[path = "incandescent/filament.rs"]
mod filament;
#[path = "incandescent/history.rs"]
mod history;
#[path = "incandescent/phase_transport.rs"]
mod phase_transport;

use fs_blake3::ContentHash;
use fs_conduction::ConductionError;
use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_couple::incandescent::PreparedIncandescentSource;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::conductor::{
    ElectricalPropertyBasis, OhmicDrive, ResolvedConductorState, graph::MaterialOhmicNode,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn point(t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with_quantity("pressure", pressure_quantity(), pressure)
        .unwrap()
}

fn knot(h: f64, t: f64, fraction: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: fraction,
        bulk_density_kg_m3: 1000.0,
    }
}

fn fixture(basis: ElectricalPropertyBasis) -> (MaterialOhmicNode, EquilibriumEnthalpyPhaseCurve) {
    source_fixture(basis, false)
}

fn source_fixture(
    basis: ElectricalPropertyBasis,
    phase_dependent: bool,
) -> (MaterialOhmicNode, EquilibriumEnthalpyPhaseCurve) {
    thermal_source_fixture(basis, phase_dependent, 1e6)
}

fn thermal_source_fixture(
    basis: ElectricalPropertyBasis,
    phase_dependent: bool,
    liquid_conductivity: f64,
) -> (MaterialOhmicNode, EquilibriumEnthalpyPhaseCurve) {
    // q(T)=T/300 or q(f)=1+2f in the selected independent electrical basis.
    let key = PropertyKey::new("authored-filament-source", basis.dimensions());
    let mut claims = ClaimSet::new();
    let validity = ValidityDomain::unconstrained()
        .with("T", 300.0, 900.0)
        .with_quantity("pressure", pressure_quantity(), 1e5, 2e5);
    let validity = if phase_dependent {
        validity.with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
    } else {
        validity
    };
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: if phase_dependent { "f" } else { "T" }.into(),
                abscissa_dims: if phase_dependent {
                    Dims::NONE
                } else {
                    Temperature::DIMS
                },
                knots: if phase_dependent {
                    vec![(0.0, 1.0), (1.0, 3.0)]
                } else {
                    vec![(300.0, 1.0), (900.0, 3.0)]
                },
                dims: basis.dimensions(),
            },
            validity,
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured incandescent feedback fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    for (name, dims, value) in [
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, 10.0),
        ("density", Density::DIMS, 1000.0),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1e6),
        (
            fs_conduction::radiation::SURFACE_EMISSIVITY_PROPERTY,
            Dims::NONE,
            0.8,
        ),
    ] {
        let phase_curve = phase_dependent
            && (name == fs_conduction::radiation::SURFACE_EMISSIVITY_PROPERTY
                || (name == "thermal_conductivity" && liquid_conductivity != value));
        let value = if phase_curve {
            PropertyValue::Curve {
                abscissa: "f".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![
                    (0.0, value),
                    (
                        1.0,
                        if name == "thermal_conductivity" {
                            liquid_conductivity
                        } else {
                            0.2
                        },
                    ),
                ],
                dims,
            }
        } else {
            PropertyValue::Scalar { value, dims }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 900.0)
                    .with_quantity("pressure", pressure_quantity(), 1e5, 2e5),
                interpolation: if phase_curve {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured incandescent thermal fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored resistor".into(),
            phase: "declared caloric chart".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = if phase_dependent {
        point(300.0, 1e5)
            .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
            .unwrap()
    } else {
        point(300.0, 1e5)
    };
    let material = resolve_material_state_point(
        &card,
        &query,
        &[basis.requirement_for(&key).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original =
        ResolvedConductorState::try_from_material_with_property(material, basis, &key).unwrap();
    let (node, curve) = if phase_dependent {
        (
            MaterialOhmicNode::try_new_with_liquid_mass_fraction(
                &card,
                original.with_liquid_mass_fraction_axis("f").unwrap(),
                "T",
                [300.0, 900.0],
                [0.0, 1.0],
                0.01,
                0.01,
            )
            .unwrap(),
            EquilibriumEnthalpyPhaseCurve::try_new(
                card.content_hash(),
                vec![
                    knot(0.0, 300.0, 0.0),
                    knot(3000.0, 600.0, 0.0),
                    knot(5000.0, 600.0, 1.0),
                    knot(8000.0, 900.0, 1.0),
                ],
            )
            .unwrap(),
        )
    } else {
        (
            MaterialOhmicNode::try_new(&card, original, "T", [300.0, 900.0], 0.01, 0.01).unwrap(),
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                card.content_hash(),
                SolidLiquidPhase::Solid,
                vec![knot(0.0, 300.0, 0.0), knot(6000.0, 900.0, 0.0)],
            )
            .unwrap(),
        )
    };
    (node, curve)
}

fn environment() -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k: 300.0,
        radiation_temperature_k: 300.0,
        convection_w_per_m2_k: 0.0,
        transport: LumpedThermalTransport::try_declared(1e6, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-9,
        maximum_thermal_residual_j: 1e-9,
    }
}

fn input(curve: &EquilibriumEnthalpyPhaseCurve, h: f64, dt: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: 0.1,
        volume_m3: 0.0001,
        surface_area_m2: 0.01,
        internal_heat_j: 0.0,
        duration_s: dt,
    }
}

#[test]
fn both_electrical_source_bases_close_their_distinct_implicit_drive_laws() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let (node, curve) = fixture(basis);
        for drive in [
            OhmicDrive::Current(2.0),
            OhmicDrive::Current(-2.0),
            OhmicDrive::Voltage(2.0),
            OhmicDrive::Voltage(-2.0),
        ] {
            let source = PreparedIncandescentSource::try_new(&node, drive).unwrap();
            let proposal = input(&curve, 0.0, 1.0);
            let step = source.advance(&environment(), proposal, || Ok(())).unwrap();
            let direct = matches!(
                (basis, drive),
                (ElectricalPropertyBasis::Resistivity, OhmicDrive::Current(_))
                    | (
                        ElectricalPropertyBasis::Conductivity,
                        OhmicDrive::Voltage(_)
                    )
            );
            // C=1 J/K. Direct P=T/75; reciprocal P=1200/T.
            let expected = if direct {
                300.0 / (1.0 - 1.0 / 75.0)
            } else {
                0.5 * (300.0 + (300.0_f64.powi(2) + 4800.0).sqrt())
            };
            assert!((step.thermal.state.temperature_k() - expected).abs() < 1e-8);
            assert_eq!(step.electrical.conductor.source_basis(), basis);
            assert_eq!(
                step.electrical.conductor.source_key(),
                node.source().original().source_key()
            );
            assert_eq!(
                step.electrical.conductor.material().query().axes()["T"],
                step.thermal.state.temperature_k()
            );
            assert_eq!(
                step.electrical
                    .conductor
                    .source_property()
                    .answer()
                    .receipt
                    .selected,
                node.source()
                    .original()
                    .source_property()
                    .answer()
                    .receipt
                    .selected
            );
            assert_eq!(
                step.thermal.additional_internal_heat_j,
                step.electrical.joule_power_w * proposal.duration_s
            );
            assert!(step.thermal.endpoint.residual_j.abs() <= 1e-9);
            assert!(step.electrical.electrical_power_residual_w.abs() < 1e-12);
            assert_eq!(proposal.initial.temperature_k(), 300.0);
        }
    }
}

#[test]
fn prescribed_heat_and_electrical_work_are_deposited_once() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let source = PreparedIncandescentSource::try_new(&node, OhmicDrive::Current(2.0)).unwrap();
    let mut proposal = input(&curve, 0.0, 1.0);
    proposal.internal_heat_j = 10.0;
    let step = source.advance(&environment(), proposal, || Ok(())).unwrap();
    assert!((step.thermal.state.temperature_k() - 310.0 / (1.0 - 1.0 / 75.0)).abs() < 1e-8);
    assert!(
        (step.thermal.endpoint.stored_energy_change_j
            - proposal.internal_heat_j
            - step.thermal.additional_internal_heat_j)
            .abs()
            < 1e-9
    );
}

#[test]
fn turning_off_preserves_stored_heat_and_cools_through_original_boundaries() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let on = PreparedIncandescentSource::try_new(&node, OhmicDrive::Voltage(2.0)).unwrap();
    let off = PreparedIncandescentSource::try_new(&node, OhmicDrive::Voltage(0.0)).unwrap();
    let mut bath = environment();
    bath.convection_w_per_m2_k = 2.0;
    bath.transport = LumpedThermalTransport::try_declared(1e6, 0.8).unwrap();
    let mut proposal = input(&curve, 0.0, 0.1);
    let mut supplied = 0.0;
    let mut peak = 300.0;
    for index in 0..40 {
        let source = if index < 20 { &on } else { &off };
        let previous = proposal.initial;
        let step = source.advance(&bath, proposal, || Ok(())).unwrap();
        supplied += step.thermal.additional_internal_heat_j
            + step.thermal.endpoint.convection_heat_j
            + step.thermal.endpoint.radiation_heat_j;
        if index < 20 {
            assert!(step.thermal.state.temperature_k() > previous.temperature_k());
            peak = step.thermal.state.temperature_k();
        } else {
            assert_eq!(step.electrical.joule_power_w, 0.0);
            assert!(step.thermal.state.temperature_k() < previous.temperature_k());
            assert!(step.thermal.state.temperature_k() > 300.0);
        }
        proposal.initial = step.thermal.state;
    }
    assert!(proposal.initial.temperature_k() < peak);
    assert!((proposal.mass_kg * proposal.initial.specific_enthalpy_j_kg() - supplied).abs() < 4e-8);
}

#[test]
fn temperature_only_electrical_source_preserves_native_latent_storage() {
    let (node, original) = fixture(ElectricalPropertyBasis::Resistivity);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        original.material_card_identity(),
        vec![
            knot(0.0, 300.0, 0.0),
            knot(3000.0, 600.0, 0.0),
            knot(5000.0, 600.0, 1.0),
            knot(8000.0, 900.0, 1.0),
        ],
    )
    .unwrap();
    let source = PreparedIncandescentSource::try_new(&node, OhmicDrive::Current(2.0)).unwrap();
    let step = source
        .advance(&environment(), input(&curve, 4000.0, 0.01), || Ok(()))
        .unwrap();
    assert_eq!(step.thermal.state.temperature_k(), 600.0);
    assert!(step.thermal.state.liquid_mass_fraction() > 0.5);
    assert!((step.thermal.additional_internal_heat_j - 0.08).abs() < 1e-12);
    assert!((step.thermal.state.specific_enthalpy_j_kg() - 4000.8).abs() < 1e-8);
}

#[test]
fn incompatible_geometry_source_identity_and_uncovered_chart_refuse() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let source = PreparedIncandescentSource::try_new(&node, OhmicDrive::Current(0.0)).unwrap();
    let mut wrong_volume = input(&curve, 0.0, 1.0);
    wrong_volume.volume_m3 *= 2.0;
    assert!(
        source
            .advance(&environment(), wrong_volume, || Ok(()))
            .is_err()
    );
    let foreign = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([9; 32]),
        SolidLiquidPhase::Solid,
        vec![knot(0.0, 300.0, 0.0), knot(6000.0, 900.0, 0.0)],
    )
    .unwrap();
    assert!(
        source
            .advance(&environment(), input(&foreign, 0.0, 1.0), || Ok(()))
            .is_err()
    );
    let uncovered = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        curve.material_card_identity(),
        SolidLiquidPhase::Solid,
        vec![knot(0.0, 300.0, 0.0), knot(7000.0, 1000.0, 0.0)],
    )
    .unwrap();
    assert!(
        source
            .advance(&environment(), input(&uncovered, 0.0, 1.0), || Ok(()))
            .is_err()
    );
}

#[test]
fn sourced_step_limit_and_cancellation_refuse_without_mutating_inputs() {
    let (node, curve) = fixture(ElectricalPropertyBasis::Resistivity);
    let source = PreparedIncandescentSource::try_new(&node, OhmicDrive::Current(2.0)).unwrap();
    let proposal = input(&curve, 0.0, 1.0);
    let limit = source
        .maximum_step_s(&environment(), proposal, || Ok(()))
        .unwrap();
    assert!(limit < 37.5 && limit > 37.499_999);
    let overlong = UniformEnthalpyStepInput {
        duration_s: limit * 1.01,
        ..proposal
    };
    assert!(source.advance(&environment(), overlong, || Ok(())).is_err());
    assert!(
        source
            .advance(&environment(), proposal, || Err(
                ConductionError::Cancelled {
                    stage: "incandescent test",
                    at: 0
                }
            ))
            .is_err()
    );
    assert_eq!(proposal.initial.temperature_k(), 300.0);
    assert!(PreparedIncandescentSource::try_new(&node, OhmicDrive::Current(f64::NAN)).is_err());
}

#[test]
fn same_card_caloric_and_transport_sources_retain_their_fixed_conditions() {
    let (node, authored) = fixture(ElectricalPropertyBasis::Resistivity);
    let source = PreparedIncandescentSource::try_new(&node, OhmicDrive::Voltage(0.0)).unwrap();
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    for pressure in [1e5, 2e5] {
        let lower = point(300.0, pressure);
        let upper = point(900.0, pressure);
        let caloric = ResolvedHeatCapacityCurve::try_new(
            node.source().card(),
            HeatCapacitySourceQuery {
                heat_capacity: &cp,
                density: &density,
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            0.0,
            1e-8,
        )
        .unwrap();
        let thermal_result =
            source.advance(&environment(), input(caloric.curve(), 0.0, 1.0), || Ok(()));
        assert_eq!(thermal_result.is_ok(), pressure == 1e5);
        let mut bath = environment();
        bath.transport = LumpedThermalTransport::from_material_card_without_radiation(
            node.source().card(),
            &conductivity,
            "T",
            &[lower, upper],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        let transport_result = source.advance(&bath, input(&authored, 0.0, 1.0), || Ok(()));
        assert_eq!(transport_result.is_ok(), pressure == 1e5);
    }
}

#[test]
fn g1_melting_resistance_and_drive_feedback_follow_native_latent_enthalpy() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let (node, curve) = source_fixture(basis, true);
        for drive in [OhmicDrive::Current(2.0), OhmicDrive::Voltage(-2.0)] {
            assert!(PreparedIncandescentSource::try_new(&node, drive).is_err());
            let source =
                PreparedIncandescentSource::try_new_with_caloric_curve(&node, drive, &curve)
                    .unwrap();
            let proposal = input(&curve, 4000.0, 1.0);
            let direct = matches!(
                (basis, drive),
                (ElectricalPropertyBasis::Resistivity, OhmicDrive::Current(_))
                    | (
                        ElectricalPropertyBasis::Conductivity,
                        OhmicDrive::Voltage(_)
                    )
            );
            let limit = source
                .maximum_step_s(&environment(), proposal, || Ok(()))
                .unwrap();
            if direct {
                // dT/dh=0 on this plateau, but dP/dh=0.004 W/(J/kg).
                assert!(limit < 12.5 && limit > 12.499_999);
            } else {
                assert_eq!(limit, f64::INFINITY);
            }
            let step = source.advance(&environment(), proposal, || Ok(())).unwrap();
            // Independently solve m*delta_h=P(h0+delta_h)*dt.
            let delta_h = if direct {
                8.0 / (0.1 - 0.004)
            } else {
                80_000.0 / (2000.0 + 4_160_000.0_f64.sqrt())
            };
            let expected_fraction = 0.5 + delta_h / 2000.0;
            assert_eq!(step.thermal.state.temperature_k(), 600.0);
            assert!((step.thermal.state.specific_enthalpy_j_kg() - 4000.0 - delta_h).abs() < 1e-8);
            assert!((step.thermal.state.liquid_mass_fraction() - expected_fraction).abs() < 1e-11);
            assert_eq!(
                step.electrical.conductor.material().query().axes()["f"],
                step.thermal.state.liquid_mass_fraction()
            );
            let q = 1.0 + 2.0 * expected_fraction;
            let resistance = if basis == ElectricalPropertyBasis::Resistivity {
                q
            } else {
                1.0 / q
            };
            assert!((step.electrical.conductor.resistance_ohm() - resistance).abs() < 1e-10);
            assert_eq!(
                step.thermal.additional_internal_heat_j,
                step.electrical.joule_power_w * proposal.duration_s
            );
            assert!(step.thermal.endpoint.residual_j.abs() <= 1e-9);
        }
    }
}

#[test]
fn g4_phase_feedback_rejects_a_different_chart_and_frozen_phase_transport() {
    let (node, curve) = source_fixture(ElectricalPropertyBasis::Resistivity, true);
    let source = PreparedIncandescentSource::try_new_with_caloric_curve(
        &node,
        OhmicDrive::Current(2.0),
        &curve,
    )
    .unwrap();
    let different_chart = EquilibriumEnthalpyPhaseCurve::try_new(
        curve.material_card_identity(),
        vec![
            knot(0.0, 300.0, 0.0),
            knot(3000.0, 600.0, 0.0),
            knot(4000.0, 600.0, 1.0),
            knot(7000.0, 900.0, 1.0),
        ],
    )
    .unwrap();
    assert!(
        source
            .advance(&environment(), input(&different_chart, 3500.0, 1.0), || Ok(
                ()
            ))
            .is_err()
    );
    let queries = [300.0, 900.0].map(|t| {
        point(t, 1e5)
            .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
            .unwrap()
    });
    let mut bath = environment();
    bath.transport = LumpedThermalTransport::from_material_card_without_radiation(
        node.source().card(),
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        "T",
        &queries,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        source
            .advance(&bath, input(&curve, 4000.0, 1.0), || Ok(()))
            .is_err()
    );
    // A wholly liquid chart has a fixed f=1 even when the electrical node's
    // original reference query was solid. Match the actual chart, not that
    // historical reference fraction, when admitting thermal transport.
    let liquid = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        curve.material_card_identity(),
        SolidLiquidPhase::Liquid,
        vec![knot(0.0, 300.0, 1.0), knot(6000.0, 900.0, 1.0)],
    )
    .unwrap();
    let source = PreparedIncandescentSource::try_new_with_caloric_curve(
        &node,
        OhmicDrive::Current(0.0),
        &liquid,
    )
    .unwrap();
    assert!(
        source
            .advance(&bath, input(&liquid, 0.0, 1.0), || Ok(()))
            .is_err()
    );
    let queries = queries.map(|query| query.with("f", 1.0).unwrap());
    bath.transport = LumpedThermalTransport::from_material_card_without_radiation(
        node.source().card(),
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        "T",
        &queries,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        source
            .advance(&bath, input(&liquid, 0.0, 1.0), || Ok(()))
            .is_ok()
    );
}
