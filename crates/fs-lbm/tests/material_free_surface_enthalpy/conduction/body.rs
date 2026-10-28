//! G1/G3/G4: two finite material stores share heat, source feedback and rollback.

#[path = "body/network.rs"]
mod network;

use super::*;
use fs_conduction::enthalpy::CellGeometry;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceBodyContact, FreeSurfaceBodyInterface, FreeSurfaceThermalBody,
};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::thermal::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RESISTANCE_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RESISTANCE, ThermalInterfacePoint,
};
use fs_material::thermal::source::ThermalMaterialPoint;

const WALL_LIQUID: [usize; 3] = [0, 3, 3];

pub(super) fn wall_card(scale: f64, upper_k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(
                "thermal_conductivity",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                )),
            ),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, scale * 100.0), (500.0, scale * 300.0)],
                dims: Dims([1, 1, -3, -1, 0, 0]),
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, upper_k),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic finite solid thermal fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

pub(super) fn wall_chart(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 150000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2000.0,
            },
        ],
    )
    .unwrap()
}

pub(super) fn wall<'a>(
    card: &'a MaterialCard,
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    temperature_k: f64,
) -> FreeSurfaceThermalBody<'a> {
    FreeSurfaceThermalBody::try_new(
        42,
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        chart,
        CellGeometry {
            start_m: [-0.001, 0.0035, 0.0035],
            end_m: [0.0, 0.0035, 0.0035],
            volume_m3: 1e-9,
            surface_area_m2: 6e-6,
        },
        (temperature_k - 300.0) * 500.0,
        &PropertyKey::with_quantity(
            "thermal_conductivity",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
        ),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

pub(super) fn contact() -> FreeSurfaceBodyContact {
    FreeSurfaceBodyContact {
        body_id: 42,
        liquid: WALL_LIQUID,
        wall: [-1, 3, 3],
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }
}

pub(super) fn interface_card(
    liquid: &MaterialCard,
    solid: &MaterialCard,
    scale: f64,
    liquid_upper_k: f64,
    fallback: bool,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    let claim = PropertyClaim {
        key: PropertyKey::new(RESISTANCE, RESISTANCE_DIMS),
        value: PropertyValue::Curve {
            abscissa: "meanT".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, scale), (500.0, 2.0 * scale)],
            dims: RESISTANCE_DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "liquidT",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                liquid_upper_k,
            )
            .with_quantity(
                "solidT",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                500.0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic ordered thermal interface".into(),
            license: "test-only".into(),
            artifact: None,
        },
    };
    claims.insert_claim(claim.clone()).unwrap();
    if fallback {
        let mut other = claim;
        other.validity = ValidityDomain::unconstrained()
            .with_quantity(
                "liquidT",
                QuantitySpec::dimensional(Temperature::DIMS),
                400.000001,
                500.0,
            )
            .with_quantity(
                "solidT",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                500.0,
            );
        claims.insert_claim(other).unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: liquid.id().clone(),
            texture_frame: "liquid".into(),
        },
        SurfaceSpec {
            material: solid.id().clone(),
            texture_frame: "wall".into(),
        },
        SystemContext {
            medium: "declared contact".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "fixed surface".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

pub(super) fn interface_point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "liquidT",
            QuantitySpec::dimensional(Temperature::DIMS),
            400.0,
        )
        .unwrap()
        .with_quantity(
            "solidT",
            QuantitySpec::dimensional(Temperature::DIMS),
            450.0,
        )
        .unwrap()
        .with_quantity("meanT", QuantitySpec::dimensional(Temperature::DIMS), 425.0)
        .unwrap()
}

pub(super) fn interface<'a>(
    liquid: &'a MaterialCard,
    solid: &'a MaterialCard,
    source: &'a InterfaceSystemCard,
    context: &'a QueryPoint,
) -> FreeSurfaceBodyInterface<'a> {
    FreeSurfaceBodyInterface::try_new(
        0,
        ThermalInterfacePoint {
            interface: source,
            materials: [liquid, solid],
            context,
            temperature_axes: ["liquidT", "solidT"],
            liquid_fraction_axes: [None; 2],
            mean_temperature_axis: Some("meanT"),
        },
        [
            caloric(liquid)
                .state_at_specific_enthalpy(200000.0)
                .unwrap(),
            wall_chart(solid)
                .state_at_specific_enthalpy(75000.0)
                .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn g1_interface_and_two_bulk_layers_share_the_actual_thermal_substeps() {
    let liquid = conducting_card(1.0);
    let solid = wall_card(0.5, 500.0);
    let chart = wall_chart(&solid);
    let context = point(450.0);
    let source = interface_card(&liquid, &solid, 1e-5, 500.0, false);
    let source_point = interface_point();
    let binding = interface(&liquid, &solid, &source, &source_point);
    with_liquid(&liquid, false, |thermal, gate| {
        let mut bodies = [wall(&solid, &chart, &context, 450.0)];
        let before = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j();
        let report = thermal
            .step_with_body_interfaces(
                &[],
                &[],
                controls(),
                &mut bodies,
                &[contact()],
                &[binding],
                &[],
                gate,
            )
            .unwrap();
        let (mut tl, mut ts) = (400.0, 450.0);
        let mut applied_g = 0.0;
        for _ in 0..4 {
            let resistance = 1e-5 * (1.0 + (0.5 * (tl + ts) - 300.0) / 200.0);
            applied_g = 1e-6 / (0.0005 / (tl - 200.0) + 0.0005 / (0.5 * (ts - 200.0)) + resistance);
            let q = 0.00025 * applied_g * (ts - tl);
            tl += q / 0.002;
            ts -= q / 0.001;
        }
        assert_eq!(report.thermal.substeps, 4);
        assert_eq!(report.flow.external_heat_j, 0.0);
        let body_report = report.bodies.unwrap();
        near(
            body_report.last_applied[0].conductance_w_k,
            applied_g,
            2e-13,
        );
        near(
            temperature(thermal.material_at(WALL_LIQUID).unwrap()),
            tl,
            2e-13,
        );
        near(
            bodies[0].specific_enthalpy_j_kg(),
            (ts - 300.0) * 500.0,
            2e-13,
        );
        let endpoint = body_report.thermal_endpoint[0].interface.as_ref().unwrap();
        assert_eq!(endpoint.card_identity(), source.content_hash());
        near(endpoint.state_coordinate("liquidT").unwrap(), tl, 2e-13);
        near(endpoint.state_coordinate("solidT").unwrap(), ts, 2e-13);
        near(
            endpoint.state_coordinate("meanT").unwrap(),
            0.5 * (tl + ts),
            2e-13,
        );
        near(
            thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j(),
            before,
            2e-12,
        );
        assert!(body_report.combined_balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn g3_interface_data_changes_liquid_material_state_and_native_motion() {
    let run = |resistance| {
        let liquid = conducting_card(1.0);
        let solid = wall_card(1.0, 500.0);
        let chart = wall_chart(&solid);
        let context = point(450.0);
        let source = interface_card(&liquid, &solid, resistance, 500.0, false);
        let source_point = interface_point();
        let binding = interface(&liquid, &solid, &source, &source_point);
        with_liquid(&liquid, true, |thermal, gate| {
            let mut bodies = [wall(&solid, &chart, &context, 450.0)];
            let before = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j();
            for _ in 0..3 {
                thermal
                    .step_with_body_interfaces(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[contact()],
                        std::slice::from_ref(&binding),
                        &[],
                        gate,
                    )
                    .unwrap();
            }
            near(
                thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j(),
                before,
                2e-11,
            );
            (
                thermal.flow().lattice().grid().state_bits(),
                temperature(thermal.material_at(WALL_LIQUID).unwrap()),
                thermal
                    .material_at(WALL_LIQUID)
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                bodies[0].specific_enthalpy_j_kg(),
            )
        })
    };
    let conductive = run(1e-6);
    let resistive = run(1e-3);
    assert_ne!(
        conductive.0, resistive.0,
        "interface data must reach native populations"
    );
    assert!(conductive.1 > resistive.1);
    assert!(conductive.2 < resistive.2);
    assert!(conductive.3 < resistive.3);
}

#[test]
fn g4_interface_support_exit_does_not_switch_claims_or_publish_either_body() {
    let liquid = conducting_card(1.0);
    let solid = wall_card(1.0, 500.0);
    let chart = wall_chart(&solid);
    let context = point(450.0);
    // A replacement exists at warmer liquid temperatures. The source selected
    // at 400 K is pinned, so the actual first heat increment must still refuse.
    let source = interface_card(&liquid, &solid, 1e-5, 400.0, true);
    let source_point = interface_point();
    let binding = interface(&liquid, &solid, &source, &source_point);
    with_liquid(&liquid, true, |thermal, gate| {
        let mut bodies = [wall(&solid, &chart, &context, 450.0)];
        let bits = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        let body_h = bodies[0].specific_enthalpy_j_kg();
        let body_source = bodies[0].material().clone();
        for interfaces in [
            vec![binding.clone(), binding.clone()],
            vec![binding.clone()],
        ] {
            assert!(
                thermal
                    .step_with_body_interfaces(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[contact()],
                        &interfaces,
                        &[],
                        gate
                    )
                    .is_err()
            );
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), body_h);
            assert_eq!(bodies[0].material(), &body_source);
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            thermal
                .step_with_body_interfaces(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[contact()],
                    &[binding],
                    &[],
                    &cancelled
                )
                .is_err()
        );
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), body_h);
    });
}

#[test]
fn g1_finite_body_heat_matches_independent_two_material_calculation() {
    let liquid = conducting_card(1.0);
    let solid = wall_card(0.5, 500.0);
    let chart = wall_chart(&solid);
    let context = point(450.0);
    with_liquid(&liquid, false, |thermal, gate| {
        let mut bodies = [wall(&solid, &chart, &context, 450.0)];
        let energy = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j();
        let report = thermal
            .step_with_body_conduction(&[], &[], controls(), &mut bodies, &[contact()], gate)
            .unwrap();
        // Independent half-cell resistance and finite masses: liquid 1 mg,
        // cp=2000; solid 2 mg, cp=500. Four explicit substeps, no reservoir.
        let (mut tl, mut ts) = (400.0, 450.0);
        let mut last_g = 0.0;
        for _ in 0..4 {
            last_g = 1e-6 / (0.0005 / (tl - 200.0) + 0.0005 / (0.5 * (ts - 200.0)));
            let q = 0.00025 * last_g * (ts - tl);
            tl += q / 0.002;
            ts -= q / 0.001;
        }
        assert_eq!(report.thermal.substeps, 4);
        assert_eq!(report.thermal.external_heat_j, 0.0);
        assert_eq!(report.flow.external_heat_j, 0.0);
        assert_eq!(report.flow.viscous_heat_j, 0.0);
        let body_report = report.bodies.unwrap();
        near(body_report.last_applied[0].conductance_w_k, last_g, 2e-13);
        near(
            temperature(thermal.material_at(WALL_LIQUID).unwrap()),
            tl,
            2e-13,
        );
        near(
            bodies[0].specific_enthalpy_j_kg(),
            (ts - 300.0) * 500.0,
            2e-13,
        );
        near(report.flow.body_heat_j, 0.001 * (450.0 - ts), 2e-12);
        near(body_report.heat_to_liquid_j, report.flow.body_heat_j, 1e-14);
        assert_eq!(bodies[0].material().card_identity(), solid.content_hash());
        assert_eq!(
            body_report.thermal_endpoint[0]
                .liquid
                .material()
                .card_identity(),
            liquid.content_hash()
        );
        assert_eq!(body_report.thermal_endpoint[0].body, *bodies[0].material());
        assert!(body_report.combined_balance_defect_j.abs() < 1e-9);
        near(
            thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j(),
            energy,
            2e-12,
        );
    });
}

#[test]
fn g3_body_material_swap_changes_native_flow_without_replacing_any_solver() {
    let run = |scale| {
        let liquid = conducting_card(1.0);
        let solid = wall_card(scale, 500.0);
        let chart = wall_chart(&solid);
        let context = point(450.0);
        with_liquid(&liquid, true, |thermal, gate| {
            let mut bodies = [wall(&solid, &chart, &context, 450.0)];
            let initial = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j();
            for _ in 0..3 {
                let report = thermal
                    .step_with_body_conduction(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[contact()],
                        gate,
                    )
                    .unwrap();
                assert_eq!(report.flow.external_heat_j, 0.0);
                assert!(report.flow.body_heat_j > 0.0);
            }
            near(
                thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j(),
                initial,
                2e-11,
            );
            (
                thermal.flow().lattice().grid().state_bits(),
                temperature(thermal.material_at(WALL_LIQUID).unwrap()),
                thermal
                    .material_at(WALL_LIQUID)
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                bodies[0].specific_enthalpy_j_kg(),
            )
        })
    };
    let slow = run(0.1);
    let fast = run(2.0);
    assert_ne!(
        slow.0, fast.0,
        "the source data must reach native populations"
    );
    assert!(fast.1 > slow.1);
    assert!(fast.2 < slow.2);
    assert!(fast.3 < slow.3, "the finite source must actually cool");
}

#[test]
fn g3_insulating_body_reduces_exactly_to_original_liquid_conduction() {
    let liquid = conducting_card(1.0);
    let reference = with_liquid(&liquid, true, |thermal, gate| {
        thermal
            .step_with_conduction(&[], &[], controls(), gate)
            .unwrap();
        (
            thermal.flow().lattice().grid().state_bits(),
            thermal.total_enthalpy_j().unwrap(),
        )
    });
    let solid = wall_card(0.0, 500.0);
    let chart = wall_chart(&solid);
    let context = point(450.0);
    with_liquid(&liquid, true, |thermal, gate| {
        let mut bodies = [wall(&solid, &chart, &context, 450.0)];
        let before = bodies[0].specific_enthalpy_j_kg();
        let report = thermal
            .step_with_body_conduction(&[], &[], controls(), &mut bodies, &[contact()], gate)
            .unwrap();
        assert_eq!(thermal.flow().lattice().grid().state_bits(), reference.0);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), reference.1);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), before);
        assert_eq!(report.flow.body_heat_j, 0.0);
    });
}

#[test]
fn g4_late_body_source_and_geometry_refusals_preserve_both_participants() {
    let liquid = conducting_card(1.0);
    // This body's conductivity is supported at its initial 350 K but not at
    // the warmer state produced by contact. The caloric chart supports both.
    let solid = wall_card(1.0, 350.0);
    let chart = wall_chart(&solid);
    let context = point(350.0);
    with_liquid(&liquid, true, |thermal, gate| {
        let mut bodies = [wall(&solid, &chart, &context, 350.0)];
        let bits = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        let h = bodies[0].specific_enthalpy_j_kg();
        let source = bodies[0].material().clone();
        for contacts in [
            vec![contact(), contact()],
            vec![FreeSurfaceBodyContact {
                wall: [1, 3, 3],
                ..contact()
            }],
            vec![contact()],
        ] {
            assert!(
                thermal
                    .step_with_body_conduction(&[], &[], controls(), &mut bodies, &contacts, gate)
                    .is_err()
            );
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), h);
            assert_eq!(bodies[0].material(), &source);
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            thermal
                .step_with_body_conduction(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[contact()],
                    &cancelled
                )
                .is_err()
        );
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), h);
    });
}

#[test]
fn g4_native_flow_refusal_discards_completed_solid_liquid_exchange() {
    let liquid = conducting_card(1.0);
    let law = LiquidTemperatureLaw::try_new(
        &liquid,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &liquid,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&liquid);
    let solid = wall_card(1.0, 500.0);
    let solid_chart = wall_chart(&solid);
    let context = point(450.0);
    let interface_source = interface_card(&liquid, &solid, 1e-5, 500.0, false);
    let interface_context = interface_point();
    let binding = interface(&liquid, &solid, &interface_source, &interface_context);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial(true), cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls {
                    max_mass_routes: 1,
                    ..FreeSurfaceEnthalpyControls::default()
                },
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        let mut bodies = [wall(&solid, &solid_chart, &context, 450.0)];
        let bits = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        let body_h = bodies[0].specific_enthalpy_j_kg();
        let body_source = bodies[0].material().clone();
        assert!(matches!(
            thermal.step_with_body_interfaces(
                &[],
                &[],
                FreeSurfaceConductionControls {
                    maximum_faces: 1,
                    ..controls()
                },
                &mut bodies,
                &[contact()],
                &[binding],
                &[],
                &gate,
            ),
            Err(FreeSurfaceEnthalpyError::Flow(_))
        ));
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), body_h);
        assert_eq!(bodies[0].material(), &body_source);
    });
}

#[test]
fn g1_body_exchange_external_heat_and_viscous_heat_are_distinct_accounts() {
    let liquid = conducting_card(1.0);
    let solid = wall_card(1.0, 500.0);
    let chart = wall_chart(&solid);
    let context = point(450.0);
    let interface_source = interface_card(&liquid, &solid, 1e-5, 500.0, false);
    let interface_context = interface_point();
    let binding = interface(&liquid, &solid, &interface_source, &interface_context);
    with_liquid(&liquid, true, |thermal, gate| {
        let mut bodies = [wall(&solid, &chart, &context, 450.0)];
        let before = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j();
        let mut external = 0.0;
        let mut molecular = 0.0;
        for _ in 0..3 {
            let report = thermal
                .step_with_body_interfaces_and_viscous_heating(
                    &[],
                    &[FreeSurfaceCellHeat {
                        coordinates: [6, 6, 6],
                        energy_j: 0.001,
                    }],
                    controls(),
                    &mut bodies,
                    &[contact()],
                    std::slice::from_ref(&binding),
                    &[],
                    gate,
                )
                .unwrap();
            near(report.conduction.flow.external_heat_j, 0.001, 1e-14);
            molecular += report.conduction.flow.viscous_heat_j;
            external += report.conduction.flow.external_heat_j;
            assert!(report.conduction.flow.body_heat_j > 0.0);
        }
        assert!(molecular > 0.0);
        near(
            thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j() - before,
            external + molecular,
            2e-8,
        );
    });
}
