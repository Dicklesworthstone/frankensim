//! Real thermal substeps -> pinned liquid sources -> native waterline and buoyancy.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    liquid_temperature::{LiquidCaloricError, LiquidTemperatureLaw, LiquidViscosityBasis},
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
};
use fs_mbd_material::hydrostatic::{FiniteTankThermalControls, FiniteTankThermalError};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, TimedEvent,
};

fn liquid_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100_000.0)
        .unwrap()
}
fn liquid_card(short_viscosity: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Density::DIMS, [1000.0, 800.0]),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            [2000.0, 4000.0],
        ),
        (
            "thermal_conductivity",
            Dims([1, 1, -3, -1, 0, 0]),
            [0.5, 1.0],
        ),
        (
            "dynamic_viscosity",
            Dims([-1, 1, -1, 0, 0, 0]),
            [0.001, 0.002],
        ),
    ] {
        let high = if short_viscosity && name == "dynamic_viscosity" {
            350.0
        } else {
            400.0
        };
        // LiquidTemperatureLaw requires thermal conductivity as a static semantic
        // quantity, not by its dimensions alone; the other properties are dimensional.
        let quantity = if name == "thermal_conductivity" {
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
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (high, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.0,
                        high,
                    )
                    .with_quantity(
                        "p",
                        QuantitySpec::dimensional(Pressure::DIMS),
                        90_000.0,
                        110_000.0,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured thermal tank oracle, not empirical data".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "thermal gravity fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &liquid_point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
    .unwrap()
}
fn chart(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new_exact(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &liquid_point(300.0),
            upper: &liquid_point(400.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
    )
    .unwrap()
}
fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-3,
    }
}
fn thermal_controls(work: f64) -> FiniteTankThermalControls {
    FiniteTankThermalControls {
        maximum_source_density_defect_kg_m3: 1e-8,
        maximum_absolute_parameter_work_j: work,
    }
}
fn thermal(
    curve: &fs_material::phase::EquilibriumEnthalpyPhaseCurve,
    mass: f64,
    volume: f64,
    specific_power: f64,
    reverse: bool,
) -> DepositionSimulation<'_> {
    let boundary = Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: mass * specific_power,
    };
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(11),
            segment_id: 31,
            body_id: 41,
            material_id: MaterialId(21),
            mass_kg: mass,
            specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [0.0, 0.0, 1.0],
                volume_m3: volume,
                surface_area_m2: 32.0,
            },
            boundary,
        }),
    }];
    if reverse {
        events.push(TimedEvent {
            time_s: 1.0,
            event: Event::Boundary(
                CellId(11),
                Boundary {
                    source_w: -boundary.source_w,
                    ..boundary
                },
            ),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(21), curve)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, step(), || Ok(()))
    .unwrap()
    .0
}
// Independent first moment of liquid: tank prism minus submerged cube,
// with tank A=16, cube cross-section B=4, bottom z=0.5, partial immersion.
fn liquid_energy(rho: f64) -> f64 {
    let height = (16000.0 / rho - 2.0) / 12.0;
    rho * 10.0 * (6.0 * height * height + 0.5)
}

#[test]
fn accepted_heat_changes_original_density_waterline_and_buoyancy_with_separate_work() {
    let card = liquid_card(false);
    let (law, caloric) = (law(&card), chart(&card));
    let tank = FiniteHydrostaticTank::new(
        law.reference_state(),
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let mesh = cube();
    let before = thermal(
        caloric.curve(),
        tank.liquid_mass_kg(),
        16.0,
        125_000.0,
        false,
    );
    let observer = tank
        .thermal_density_observer(
            &before,
            CellId(11),
            caloric.curve(),
            &law,
            &mesh,
            pose(1.5),
            thermal_controls(1e6),
            || Ok(()),
        )
        .unwrap();
    let (after, observed, heat) = before
        .advance_observed(&observer, 1.0, step(), || Ok(()))
        .unwrap();
    let current = observed.current_tank().unwrap();
    let state = observed.endpoint().state();
    assert_eq!(
        state,
        caloric
            .curve()
            .state_at_specific_enthalpy(after.cells()[&CellId(11)].specific_h_j_kg)
            .unwrap()
    );
    assert_eq!(observed.time_s(), after.time_s());
    close(state.temperature_k(), 350.0);
    close(current.transport().density_kg_m3(), 900.0);
    close(current.source().specific_heat_j_kg_k(), 3000.0);
    close(current.source().thermal_conductivity_w_m_k(), 0.75);
    close(current.source().dynamic_viscosity_pa_s(), 0.0015);
    close(
        current.source().volumetric_expansion_per_k().unwrap(),
        2.0 / 900.0,
    );
    assert!(std::ptr::eq(
        current.source(),
        observed.endpoint().law().reference_state()
    ));
    assert!(after.uses_material_law(MaterialId(21), caloric.curve()));
    assert_eq!(current.liquid_mass_kg(), 16000.0);
    close(current.liquid_volume_m3(), 16000.0 / 900.0);
    let fill = current.fill_at_pose(&mesh, pose(1.5), || false).unwrap();
    close(fill.surface_height_m, (16000.0 / 900.0 - 2.0) / 12.0);
    let solid = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let principal = solid.principal_body(1e-12, || false).unwrap();
    let rigid = principal
        .state_from_kinematics(pose(1.5), Twist::zero())
        .unwrap();
    let load = current.load(&mesh, &principal, rigid, || false).unwrap();
    close(
        load.wrench().force_world.z,
        9000.0 * 4.0 * (fill.surface_height_m - 0.5),
    );
    let audit = observed.audit();
    let work = liquid_energy(900.0) - liquid_energy(1000.0);
    close(audit.gravitational_parameter_work_j, work);
    close(audit.absolute_gravitational_parameter_work_j, work);
    close(audit.caloric_storage_change_j, 16000.0 * 125_000.0);
    close(heat.external_heat_j, audit.caloric_storage_change_j);
    assert_eq!(audit.source_updates, heat.substeps);
    assert_eq!(audit.source_updates, 4);
    // Gravity storage is a separate external debit, not heat credited a second time.
    assert!(audit.gravitational_parameter_work_j > 10_000.0);
    assert_eq!(
        after.cells()[&CellId(11)].geometry,
        before.cells()[&CellId(11)].geometry
    );
    for receipt in caloric.receipts() {
        assert_eq!(
            current
                .source()
                .material()
                .property(&receipt.property)
                .unwrap()
                .answer()
                .receipt
                .selected,
            receipt.selected
        );
    }
    assert_eq!(observer.endpoint().state().temperature_k(), 300.0);
    assert_eq!(tank.transport().density_kg_m3(), 1000.0);
}

#[test]
fn dry_displacement_still_requires_complete_empty_bath_gravity_work() {
    let card = liquid_card(false);
    let (law, caloric) = (law(&card), chart(&card));
    let tank = FiniteHydrostaticTank::new(
        law.reference_state(),
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let mesh = cube();
    let before = thermal(caloric.curve(), 16000.0, 16.0, 125_000.0, false);
    let observer = tank
        .thermal_density_observer(
            &before,
            CellId(11),
            caloric.curve(),
            &law,
            &mesh,
            pose(8.0),
            thermal_controls(1e6),
            || Ok(()),
        )
        .unwrap();
    let (_, after, _) = before
        .advance_observed(&observer, 1.0, step(), || Ok(()))
        .unwrap();
    let expected = 16000.0 * 10.0 * ((16000.0 / 900.0) / 16.0 - 1.0) / 2.0;
    close(after.audit().gravitational_parameter_work_j, expected);
    let last = after.last_density_change().unwrap();
    close(last.displacement_potential_change_j, 0.0);
    assert!(last.empty_tank_gravitational_energy_change_j > 0.0);
    assert_eq!(
        after
            .current_tank()
            .unwrap()
            .fill_at_pose(&mesh, pose(8.0), || false)
            .unwrap()
            .displaced
            .volume_m3(),
        0.0
    );
    let zero = tank
        .thermal_density_observer(
            &before,
            CellId(11),
            caloric.curve(),
            &law,
            &mesh,
            pose(8.0),
            thermal_controls(0.0),
            || Ok(()),
        )
        .unwrap();
    assert!(matches!(
        before.advance_observed(&zero, 1.0, step(), || Ok(())),
        Err(FiniteTankThermalError::ParameterWork { .. })
    ));
}

#[test]
fn cooling_cancellation_and_work_refusal_roll_back_both_native_and_source_candidates() {
    let card = liquid_card(false);
    let (law, caloric) = (law(&card), chart(&card));
    let tank = FiniteHydrostaticTank::new(
        law.reference_state(),
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let mesh = cube();
    let before = thermal(caloric.curve(), 16000.0, 16.0, 125_000.0, true);
    let observer = tank
        .thermal_density_observer(
            &before,
            CellId(11),
            caloric.curve(),
            &law,
            &mesh,
            pose(1.5),
            thermal_controls(1e6),
            || Ok(()),
        )
        .unwrap();
    let mut polls = 0;
    let (after, observed, heat) = before
        .advance_observed(&observer, 2.0, step(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    close(observed.audit().gravitational_parameter_work_j, 0.0);
    close(
        observed.audit().absolute_gravitational_parameter_work_j,
        2.0 * (liquid_energy(900.0) - liquid_energy(1000.0)),
    );
    close(observed.audit().caloric_storage_change_j, 0.0);
    close(heat.external_heat_j, 0.0);
    assert_eq!(observed.audit().source_updates, 8);
    // Back at 300 K the tank holds a refreshed source, whose receipts record
    // pinned-claim selection, while the law's reference was resolved under
    // SingleClaimOnly; fs-material's rule is that their physical tuples, not
    // their policy receipts, agree. The card, the query and every selected
    // claim and value are the reference's.
    let source = observed.current_tank().unwrap().source().material();
    let reference = law.reference_state().material();
    assert_eq!(source.card_identity(), reference.card_identity());
    assert_eq!(source.query(), reference.query());
    assert_eq!(source.properties().len(), reference.properties().len());
    for (refreshed, original) in source.properties().iter().zip(reference.properties()) {
        assert_eq!(
            refreshed.answer().receipt.selected,
            original.answer().receipt.selected
        );
        assert_eq!(
            refreshed.value_si().to_bits(),
            original.value_si().to_bits()
        );
    }
    let limited = tank
        .thermal_density_observer(
            &before,
            CellId(11),
            caloric.curve(),
            &law,
            &mesh,
            pose(1.5),
            thermal_controls(15_000.0),
            || Ok(()),
        )
        .unwrap();
    assert!(matches!(
        before.advance_observed(&limited, 2.0, step(), || Ok(())),
        Err(FiniteTankThermalError::ParameterWork { .. })
    ));
    let mut count = 0;
    assert!(matches!(
        before.advance_observed(&observer, 2.0, step(), || {
            count += 1;
            if count == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(FiniteTankThermalError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(before.time_s(), 0.0);
    assert_eq!(before.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    assert_eq!(observer.time_s(), 0.0);
    assert_eq!(limited.endpoint().state().temperature_k(), 300.0);
    assert_eq!(observer.audit().source_updates, 0);
    assert_eq!(observer.current_tank().unwrap().liquid_volume_m3(), 16.0);
    // A source observer from another accepted clock cannot restart this field.
    assert!(matches!(
        before.advance_observed(&observed, 1.0, step(), || Ok(())),
        Err(FiniteTankThermalError::Invalid(_))
    ));
    assert_eq!(after.time_s(), 2.0);
}

#[test]
fn source_support_tank_overflow_and_foreign_bindings_publish_no_replacement() {
    for short in [false, true] {
        let card = liquid_card(short);
        let (law, caloric) = (law(&card), chart(&card));
        let height = if short { 4.0 } else { 1.25 };
        let tank = FiniteHydrostaticTank::new(
            law.reference_state(),
            gravity(),
            tank_geometry(Vec3::new(4.0, 4.0, height)),
            16.0,
            volume_controls(),
        )
        .unwrap();
        let mesh = cube();
        let power = if short { 206_250.0 } else { 125_000.0 };
        let before = thermal(caloric.curve(), 16000.0, 16.0, power, false);
        let observer = tank
            .thermal_density_observer(
                &before,
                CellId(11),
                caloric.curve(),
                &law,
                &mesh,
                pose(1.5),
                thermal_controls(1e6),
                || Ok(()),
            )
            .unwrap();
        let error = before.advance_observed(&observer, 1.0, step(), || Ok(()));
        if short {
            assert!(matches!(
                error,
                Err(FiniteTankThermalError::Caloric(
                    LiquidCaloricError::Material(_)
                ))
            ));
        } else {
            assert!(matches!(
                error,
                Err(FiniteTankThermalError::Hydrostatic(
                    HydrostaticError::TankOverflow { .. }
                ))
            ));
        }
        assert_eq!(before.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(observer.endpoint().state().temperature_k(), 300.0);
        assert_eq!(observer.current_tank().unwrap().liquid_mass_kg(), 16000.0);
        let foreign_chart = caloric.curve().clone();
        let foreign = thermal(&foreign_chart, 16000.0, 16.0, 0.0, false);
        assert!(
            tank.thermal_density_observer(
                &foreign,
                CellId(11),
                caloric.curve(),
                &law,
                &mesh,
                pose(1.5),
                thermal_controls(1e6),
                || Ok(())
            )
            .is_err()
        );
        let initial = caloric.curve().state_at_specific_enthalpy(0.0).unwrap();
        // The same liquid knots in a chart of another identity: try_new admits
        // only a curve from fully solid to fully liquid, so a liquid chart is
        // built by the single-phase ingress.
        let other = fs_material::phase::EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Liquid,
            caloric.curve().knots().to_vec(),
        )
        .unwrap();
        assert!(
            law.prepare_caloric_state_endpoint(
                caloric.curve(),
                initial,
                other.state_at_specific_enthalpy(0.0).unwrap(),
                1e-8
            )
            .is_err()
        );
        let changed_context = LiquidTemperatureLaw::try_new(
            &card,
            &liquid_point(300.0).with("p", 101_000.0).unwrap(),
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            tank.thermal_density_observer(
                &before,
                CellId(11),
                caloric.curve(),
                &changed_context,
                &mesh,
                pose(1.5),
                thermal_controls(1e6),
                || Ok(())
            )
            .is_err()
        );
        let wrong_mass = thermal(caloric.curve(), 16001.0, 16.0, 0.0, false);
        assert!(
            tank.thermal_density_observer(
                &wrong_mass,
                CellId(11),
                caloric.curve(),
                &law,
                &mesh,
                pose(1.5),
                thermal_controls(1e6),
                || Ok(())
            )
            .is_err()
        );
    }
}
