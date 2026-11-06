//! Independent finite-cell heat feedback with sourced indentation and wear.
//! These synthetic data exercise the declared fixed-inventory approximation;
//! they are not material calibration or a complete body-energy model.

#[path = "plastic_sliding_heat/area.rs"]
mod area;

#[path = "../../fs-contact/tests/source_point/common.rs"]
mod source_common;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::hash_domain;
use fs_contact::interface_binding::thermal_friction::{SourceThermalCoulomb, ThermalFrictionAxes};
use fs_contact::interface_binding::{
    ADHESION_ENERGY_DIMS, ADHESION_ENERGY_PROPERTY, bind_dry_interface_system_card,
    bind_isotropic_elastic_interface, bind_normal_interface_state,
};
use fs_contact::normal_patch::plastic_source::wear::SourcePlasticArchardWear;
use fs_contact::normal_patch::plastic_source::wear::sliding::{
    PlasticSlidingInput, SourcePlasticSlidingWear, SourcePlasticSlidingWearState,
};
use fs_contact::normal_patch::plastic_source::{
    PlasticContactControls, PlasticContactEnvelope, SourcePlasticContact,
};
use fs_contact::normal_patch::{
    FiniteGapGrid, FiniteGapPlasticDrive, FiniteGapPlasticLoadControls,
};
use fs_couple::friction_heat::plastic_wear::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    HardnessLoadStep, HardnessTestContext, InterfaceSystemCard, InterpolationPolicy, MaterialCard,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::hardness::{INDENTATION_HARDNESS_QUANTITY, IndentationHardnessStatePoint};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point, resolve_isotropic_solid_state_point,
    resolve_material_state_point,
};
use fs_qty::{Density, Dims, QtyAny, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::{InputAuthority, TangentialSolveConfig};
use std::cell::Cell;
use std::collections::BTreeMap;

const T0: f64 = 293.15;
const SIDE: f64 = 1e-4;
const AREA: f64 = SIDE * SIDE;
const LENGTH: f64 = 1e-3;
const RHO: f64 = 1000.0;
const CP: f64 = 500.0;
const CONDUCTIVITY: f64 = 10.0;
const CELLS: [CellId; 2] = [CellId(10), CellId(20)];
const LINK: u64 = 31;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 51,
                kernel_id: 205,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic hardness/sliding/finite-heat reference".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}

fn claim(key: PropertyKey, values: [f64; 2]) -> PropertyClaim {
    let dims = key.quantity().dims();
    let observations = key
        .hardness_test()
        .map(HardnessTestContext::observation)
        .into_iter()
        .collect();
    PropertyClaim {
        key,
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, values[0]), (320.0, values[1])],
            dims,
        },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations,
        provenance: provenance(),
    }
}

fn material(name: &str, side: usize) -> (MaterialCard, PropertyKey) {
    let base = source_common::solid(name, RHO, 2e7);
    let mut claims = base.claims().clone();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: name.into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: hash_domain("plastic-sliding-heat-fixture", name.as_bytes()),
            caveats: "manufactured numerical data".into(),
            provenance: provenance(),
        })
        .unwrap();
    let hardness =
        PropertyKey::with_quantity("indentation-hardness", INDENTATION_HARDNESS_QUANTITY)
            .with_hardness_test(
                HardnessTestContext::new(
                    "rigid sphere radius 1 mm",
                    vec![
                        HardnessLoadStep::new(
                            QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                            QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                        )
                        .unwrap(),
                    ],
                    "synthetic temperature-dependent indentation pressure",
                    observation,
                )
                .unwrap(),
            )
            .unwrap();
    let h = if side == 0 {
        [280.0, 320.0].map(|t| 1e6 * (1.0 - 0.01 * (t - T0)))
    } else {
        [2e6; 2]
    };
    claims.insert_claim(claim(hardness.clone(), h)).unwrap();
    claims
        .insert_claim(claim(
            PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            [CP; 2],
        ))
        .unwrap();
    (
        MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap(),
        hardness,
    )
}

struct Fixture {
    cards: [MaterialCard; 2],
    caloric: [ResolvedHeatCapacityCurve; 2],
    source: SourcePlasticSlidingWear,
}

fn fixture(cx: &Cx<'_>) -> Fixture {
    let [(a, ha), (b, hb)] = [material("heat-body-A", 0), material("heat-body-B", 1)];
    let cards = [a, b];
    let hardness_keys = [ha, hb];
    let caloric = cards.each_ref().map(|card| {
        ResolvedHeatCapacityCurve::try_new(
            card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(
                    "specific_heat_capacity",
                    SPECIFIC_HEAT_CAPACITY_DIMS,
                ),
                density: &PropertyKey::new("density", Density::DIMS),
                temperature_axis: "T",
                lower: &point(280.0),
                upper: &point(320.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            0.0,
            1e-8,
        )
        .unwrap()
    });
    let base = source_common::bind(&cards[0], &cards[1], None, 0.0);
    let old = base.normal_state().elastic().interface().card();
    let mut claims = old.claims().clone();
    let wear_keys = ["k_a", "k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
    for (key, values) in [
        (wear_keys[0].clone(), [1e-6; 2]),
        (wear_keys[1].clone(), [2e-6; 2]),
        (PropertyKey::new("mu_static", Dims::NONE), [0.9; 2]),
        (
            PropertyKey::new("mu_kinetic", Dims::NONE),
            [280.0, 320.0].map(|t| 0.4 + 0.01 * (t - T0)),
        ),
    ] {
        claims.insert_claim(claim(key, values)).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        old.surface_a().clone(),
        old.surface_b().clone(),
        old.context().clone(),
        claims,
        old.models().to_vec(),
    )
    .unwrap();
    let query = point(T0);
    let select = MaterialPropertySelection::SingleClaimOnly;
    let solids = cards
        .each_ref()
        .map(|card| resolve_isotropic_solid_state_point(card, &query, select.clone()).unwrap());
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &solids[0], &solids[1]).unwrap();
    let adhesion = resolve_interface_state_point(
        &card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        select.clone(),
    )
    .unwrap();
    let normal = bind_normal_interface_state(&elastic, &adhesion).unwrap();
    let hardness = [0, 1].map(|i| {
        let source = resolve_material_state_point(
            &cards[i],
            &query,
            &[ScalarPropertyRequirement::try_with_key(
                &hardness_keys[i],
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            select.clone(),
        )
        .unwrap();
        IndentationHardnessStatePoint::try_from_resolved(source, hardness_keys[i].clone()).unwrap()
    });
    let coefficients = resolve_interface_state_point(
        &card,
        &query,
        &wear_keys.each_ref().map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        }),
        select.clone(),
    )
    .unwrap();
    let mut gap = vec![1e-3; 25];
    gap[12] = 0.0;
    let normal = SourcePlasticContact::new(
        normal,
        hardness.clone(),
        FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: SIDE,
            cell_depth_m: SIDE,
        },
        gap,
        PlasticContactControls {
            maximum_active_set_iterations: 128,
            complementarity_tolerance_m: 1e-12,
            boundary_clearance_cells: 1,
        },
        PlasticContactEnvelope {
            maximum_approach_m: 1e-4,
            maximum_plastic_depth_m: 1e-4,
            half_space_depth_m: [LENGTH; 2],
            maximum_patch_to_depth: [0.1; 2],
        },
        cx,
    )
    .unwrap();
    let wear = SourcePlasticArchardWear::new(normal, coefficients, wear_keys, hardness, [5e-5; 2])
        .unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    let friction = SourceThermalCoulomb::resolve(
        dry,
        query,
        requirement("mu_static"),
        requirement("mu_kinetic"),
        select,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: None,
        },
    )
    .unwrap();
    Fixture {
        cards,
        caloric,
        source: SourcePlasticSlidingWear::new(wear, friction).unwrap(),
    }
}

fn network(curves: &[ResolvedHeatCapacityCurve; 2], initial_k: f64) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for side in 0..2 {
        let (start, end) = if side == 0 {
            (-LENGTH, 0.0)
        } else {
            (0.0, LENGTH)
        };
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELLS[side],
                segment_id: side as u64,
                body_id: side as u64,
                material_id: MaterialId(side as u64),
                mass_kg: RHO * AREA * LENGTH,
                specific_h_j_kg: curves[side]
                    .state_at_temperature(initial_k)
                    .unwrap()
                    .specific_enthalpy_j_kg(),
                geometry: CellGeometry {
                    start_m: [start, 0.0, 0.0],
                    end_m: [end, 0.0, 0.0],
                    volume_m3: AREA * LENGTH,
                    surface_area_m2: 2.0 * AREA + 4.0 * SIDE * LENGTH,
                },
                boundary: Boundary {
                    ambient_k: T0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: CELLS[0],
            b: CELLS[1],
            area_m2: AREA,
            conductance_w_k: series_conductance_w_k(AREA, [0.5 * LENGTH; 2], [CONDUCTIVITY; 2])
                .unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 128,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn normal_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1e-4,
        maximum_approach_m: 1e-4,
        absolute_force_tolerance_n: 1e-14,
        relative_force_tolerance: 1e-12,
        approach_tolerance_m: 1e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn controls() -> PlasticSlidingHeatControls {
    PlasticSlidingHeatControls {
        friction: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 1.0,
            force_absolute_tolerance_n: 1e-14,
            force_relative_tolerance: 1e-12,
        },
        maximum_cells: 2,
        maximum_trials: 400,
        maximum_newton_iterations: 16,
        difference_step_j: 1e-8,
        temperature_scale_j_per_k: 5e-6,
        heat: FrictionHeatTolerances {
            cell_j: 1e-11,
            thermal_j: 2e-11,
            interface_j: 1e-11,
            total_j: 3e-11,
        },
        contact_work_tolerance_j: 1e-12,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-12,
            linear_restart: 3,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn limits() -> PlasticSlidingHeatLimits {
    PlasticSlidingHeatLimits {
        maximum_removed_volume_fraction: [1e-3; 2],
        maximum_absolute_parameter_work_j: 1e-7,
    }
}

fn motion(distance: f64) -> Vec<PlasticSlidingInput> {
    let mut values = vec![
        PlasticSlidingInput {
            released_increment_m: [0.0; 2],
            compliance_m_per_n: 0.0
        };
        25
    ];
    values[12].released_increment_m[0] = distance;
    values
}

// Explicit constant conductivity is an authored fixture parameter. The payload
// records the actual caloric/interface endpoint used by this source callback.
fn transport(
    thermal: &DepositionSimulation<'_>,
    interface_k: f64,
) -> Result<PlasticSlidingHeatSource<[f64; 3]>, DepositionError> {
    Ok(PlasticSlidingHeatSource {
        conductivity_w_m_k: [CONDUCTIVITY; 2],
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: [
                thermal.temperature_k(CELLS[0])?,
                thermal.temperature_k(CELLS[1])?,
                interface_k,
            ],
        },
    })
}

fn loaded(fixture: &Fixture, cx: &Cx<'_>) -> SourcePlasticSlidingWearState {
    fixture
        .source
        .advance_normal(&fixture.source.initial_state(), 5e-5, cx)
        .unwrap()
        .candidate
}

fn bind<'a>(
    fixture: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    route: &'a FrictionHeatRoute,
    state: SourcePlasticSlidingWearState,
    limits: PlasticSlidingHeatLimits,
    cx: &Cx<'_>,
) -> PlasticSlidingHeatContact<'a> {
    PlasticSlidingHeatContact::new(
        thermal,
        fixture.source.clone(),
        state,
        [&fixture.cards[0], &fixture.cards[1]],
        [&fixture.caloric[0], &fixture.caloric[1]],
        route,
        limits,
        cx,
    )
    .unwrap()
}

#[test]
fn coupled_temperature_hardness_and_friction_match_independent_quadratic_and_wear() {
    with_cx(|cx, _| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx);
        let saved_state = contact.state().clone();
        let saved_cells = thermal.cells().clone();
        let step = contact
            .step_candidate(
                &thermal,
                0.1,
                &motion(0.01),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        // Symmetric backward Euler: theta_body = b*theta_interface, and
        // 2*C*b*theta_interface = A*H(T_body)*mu(T_interface)*distance.
        let capacity = RHO * AREA * LENGTH * CP;
        let half_g = CONDUCTIVITY * AREA / (0.5 * LENGTH);
        let b = 0.1 * half_g / (capacity + 0.1 * half_g);
        let measure = AREA * 1e6 * 0.01;
        let qa = measure * 0.01 * b * 0.01;
        let qb = 2.0 * capacity * b - measure * (0.01 - 0.4 * 0.01 * b);
        let theta = 2.0 * measure * 0.4 / (qb + (qb * qb + 4.0 * qa * measure * 0.4).sqrt());
        let body_k = T0 + b * theta;
        let pressure = 1e6 * (1.0 - 0.01 * b * theta);
        let heat = AREA * pressure * (0.4 + 0.01 * theta) * 0.01;
        let report = step.report();
        close(report.heat.interface_temperature_k, T0 + theta, 0.0, 2e-6);
        for (side, id) in CELLS.into_iter().enumerate() {
            close(step.thermal().temperature_k(id).unwrap(), body_k, 0.0, 2e-6);
            close(report.source[side], body_k, 0.0, 2e-6);
            let expected_volume =
                [1e-6, 2e-6][side] * pressure * 0.01 * AREA / [pressure, 2e6][side];
            close(report.removed_volume_m3[side], expected_volume, 2e-6, 1e-23);
            close(
                report.removed_volume_fraction[side],
                expected_volume / (AREA * LENGTH),
                2e-6,
                1e-15,
            );
        }
        close(report.frictional_work_j, heat, 2e-6, 1e-12);
        assert!(report.frictional_work_j > 1.05 * measure * 0.4);
        assert_eq!(report.heat.dissipated_j, report.frictional_work_j);
        assert!(report.normal_plastic_work_j > 0.0);
        assert!(report.temperature_parameter_work_j < 0.0);
        assert!(report.wear_parameter_work_j < 0.0);
        close(
            report.heat.thermal.stored_energy_change_j,
            report.frictional_work_j,
            0.0,
            3e-11,
        );
        assert!(
            (report.heat.thermal.stored_energy_change_j
                - report.frictional_work_j
                - report.normal_plastic_work_j)
                .abs()
                > 0.5 * report.normal_plastic_work_j
        );
        let PlasticSlidingHeatMechanics::HeldApproach(mechanics) = &report.mechanics else {
            panic!("wrong normal boundary")
        };
        close(
            mechanics.temperature.normal.response.normal.pressure_pa[12],
            pressure,
            2e-6,
            1e-5,
        );
        close(
            mechanics.sliding.cells[12].response.slip_increment_m[0],
            0.01,
            0.0,
            1e-15,
        );
        assert_eq!(
            mechanics.sliding.cells[12].response.slip_increment_m[1],
            0.0
        );
        assert!(report.trials > 1);
        assert_eq!(contact.state(), &saved_state);
        assert_eq!(contact.source(), &data.source);
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(step.contact().time_s(), 0.1);
    });
}

#[test]
fn maintained_load_and_fixed_drive_retain_accepted_shear_wear_and_clock() {
    with_cx(|cx, _| {
        let data = fixture(cx);
        for use_drive in [false, true] {
            let thermal = network(&data.caloric, T0);
            let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
            let drive = FiniteGapPlasticDrive {
                stiffness_n_per_m: 1000.0,
                position_m: 1e-5,
            };
            let boundary = if use_drive {
                PlasticSlidingHeatBoundary::FixedDrive {
                    drive,
                    controls: normal_controls(),
                }
            } else {
                PlasticSlidingHeatBoundary::MaintainedLoad {
                    force_n: 0.005,
                    controls: normal_controls(),
                }
            };
            let initial = data.source.initial_state();
            let state = if use_drive {
                data.source
                    .advance_normal_with_drive(&initial, drive, normal_controls(), cx)
                    .unwrap()
                    .candidate
            } else {
                data.source
                    .advance_normal_at_load(&initial, 0.005, normal_controls(), cx)
                    .unwrap()
                    .candidate
            };
            let contact = bind(&data, &thermal, &route, state, limits(), cx);
            let first = contact
                .step_candidate(
                    &thermal,
                    0.1,
                    &motion(0.002),
                    boundary,
                    controls(),
                    &[0.0; 2],
                    T0,
                    transport,
                    cx,
                )
                .unwrap();
            let (continued, heated, first_report) = first.into_parts();
            let saved = continued.state().clone();
            let second = continued
                .step_candidate(
                    &heated,
                    0.2,
                    &motion(-0.001),
                    boundary,
                    controls(),
                    &[0.0; 2],
                    continued.source().interface_temperature_k(),
                    transport,
                    cx,
                )
                .unwrap();
            assert_eq!(continued.state(), &saved);
            assert_eq!(continued.time_s(), 0.1);
            assert_eq!(heated.time_s(), 0.1);
            assert_eq!(second.contact().time_s(), 0.2);
            for side in 0..2 {
                assert!(
                    second.report().removed_volume_m3[side] > first_report.removed_volume_m3[side]
                );
            }
            close(
                second.report().mechanics.cells()[12]
                    .response
                    .slip_increment_m[0],
                -0.001,
                0.0,
                1e-15,
            );
            let retained_slip = second.contact().state().tangential_states()[12].slip_m();
            close(retained_slip[0], 0.001, 0.0, 1e-15);
            assert_eq!(retained_slip[1], 0.0);
            assert!(second.report().frictional_work_j > 0.0);
            match &second.report().mechanics {
                PlasticSlidingHeatMechanics::MaintainedLoad(response) => {
                    close(
                        response.sliding.boundary.achieved_force_n,
                        0.005,
                        2e-10,
                        1e-13,
                    );
                    assert_eq!(
                        response.sliding.normal.predecessor_state_identity,
                        response.temperature.normal.successor_state_identity
                    );
                }
                PlasticSlidingHeatMechanics::FixedDrive(response) => {
                    assert_eq!(response.sliding.boundary.boundary, drive);
                    assert!(
                        response.sliding.boundary.force_residual_n.abs()
                            <= response.sliding.boundary.force_tolerance_n
                    );
                    assert_eq!(response.sliding.boundary.equilibrium_evaluations, 1);
                }
                _ => panic!("normal boundary was replaced"),
            }
        }
    });
}

#[test]
fn caloric_source_clock_and_solver_refusals_do_not_publish_partial_histories() {
    with_cx(|cx, gate| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let state = loaded(&data, cx);
        let foreign_caloric = data.caloric.clone();
        assert!(
            PlasticSlidingHeatContact::new(
                &thermal,
                data.source.clone(),
                state.clone(),
                [&data.cards[0], &data.cards[1]],
                [&foreign_caloric[0], &foreign_caloric[1]],
                &route,
                limits(),
                cx
            )
            .is_err()
        );
        let warm = network(&data.caloric, T0 + 1.0);
        assert!(
            PlasticSlidingHeatContact::new(
                &warm,
                data.source.clone(),
                state.clone(),
                [&data.cards[0], &data.cards[1]],
                [&data.caloric[0], &data.caloric[1]],
                &route,
                limits(),
                cx
            )
            .is_err()
        );
        let contact = bind(&data, &thermal, &route, state, limits(), cx);
        let saved_state = contact.state().clone();
        let saved_cells = thermal.cells().clone();
        let later = thermal
            .advance_to(
                0.01,
                StepConfig {
                    max_dt_s: 0.01,
                    energy_tolerance_j: 1e-12,
                },
                || Ok(()),
            )
            .unwrap()
            .0;
        assert!(matches!(
            contact.step_candidate(
                &later,
                0.1,
                &motion(0.01),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx
            ),
            Err(PlasticSlidingHeatError::Invalid(_))
        ));
        let mut short = controls();
        short.maximum_newton_iterations = 1;
        assert!(matches!(
            contact.step_candidate(
                &thermal,
                0.1,
                &motion(0.01),
                PlasticSlidingHeatBoundary::HeldApproach,
                short,
                &[0.0; 2],
                T0,
                transport,
                cx
            ),
            Err(PlasticSlidingHeatError::NotConverged(_))
        ));
        let calls = Cell::new(0);
        let refused = contact.step_candidate(
            &thermal,
            0.1,
            &motion(0.01),
            PlasticSlidingHeatBoundary::HeldApproach,
            controls(),
            &[0.0; 2],
            T0,
            |endpoint, interface| {
                calls.set(calls.get() + 1);
                if calls.get() > 1 {
                    Err(DepositionError::InvalidInput(
                        "deliberate late source refusal",
                    ))
                } else {
                    transport(endpoint, interface)
                }
            },
            cx,
        );
        assert!(matches!(refused, Err(PlasticSlidingHeatError::Source(_))));
        assert!(calls.get() > 1);
        let zero_limits = PlasticSlidingHeatLimits {
            maximum_removed_volume_fraction: [0.0; 2],
            ..limits()
        };
        let no_wear = bind(
            &data,
            &thermal,
            &route,
            saved_state.clone(),
            zero_limits,
            cx,
        );
        assert!(matches!(
            no_wear.step_candidate(
                &thermal,
                0.1,
                &motion(0.01),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx
            ),
            Err(PlasticSlidingHeatError::Invalid(_))
        ));
        let mut sticking = motion(1e-7);
        sticking[12].compliance_m_per_n = 0.001;
        let stuck = no_wear
            .step_candidate(
                &thermal,
                0.1,
                &sticking,
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        assert_eq!(stuck.report().frictional_work_j, 0.0);
        assert_eq!(stuck.report().removed_volume_m3, [0.0; 2]);
        gate.request();
        assert!(
            contact
                .step_candidate(
                    &thermal,
                    0.1,
                    &motion(0.01),
                    PlasticSlidingHeatBoundary::HeldApproach,
                    controls(),
                    &[0.0; 2],
                    T0,
                    transport,
                    cx
                )
                .is_err()
        );
        assert_eq!(contact.state(), &saved_state);
        assert_eq!(contact.source(), &data.source);
        assert_eq!(contact.time_s(), 0.0);
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
