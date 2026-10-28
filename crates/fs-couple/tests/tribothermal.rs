//! G1/G3 finite two-body sliding references; source data is synthetic, not a
//! measured friction calibration. Native execution is required for qualification.
#[path = "circuit_enthalpy/common.rs"]
mod common;
#[path = "tribothermal/dynamics.rs"]
mod dynamics;
#[path = "tribothermal/field_motion.rs"]
mod field_motion;
#[path = "tribothermal/partition.rs"]
mod partition;

use fs_couple::tribothermal::material::{KINETIC_FRICTION_PROPERTY, STATIC_FRICTION_PROPERTY};
use fs_couple::tribothermal::*;
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use fs_tribo::{HeatPartition, InputAuthority, TangentialSolveConfig};
use std::cell::Cell;
use std::collections::BTreeMap;
use std::convert::Infallible;

const CELLS: [CellId; 2] = [CellId(30), CellId(10)];
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-6, "{a} != {b}");
}
fn other(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, Infallible> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}
fn controls() -> SlidingSolve {
    SlidingSolve {
        maximum_cells: 2,
        maximum_friction_trials: 400,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-3,
        cell_energy_tolerance_j: 1e-7,
        thermal_energy_tolerance_j: 2e-7,
        total_energy_tolerance_j: 3e-7,
        tangential: TangentialSolveConfig {
            max_iterations: 32,
            maximum_force_n: 1e5,
            force_absolute_tolerance_n: 1e-9,
            force_relative_tolerance: 1e-10,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn context() -> QueryPoint {
    // Stale authored coordinates MUST be replaced by the actual trial state.
    QueryPoint::new()
        .with_quantity("Ta", QuantitySpec::dimensional(Temperature::DIMS), 290.0)
        .unwrap()
        .with_quantity("Tb", QuantitySpec::dimensional(Temperature::DIMS), 310.0)
        .unwrap()
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 7.0)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.01)
        .unwrap()
        .with_quantity("distance", QuantitySpec::dimensional(LENGTH), 9.0)
        .unwrap()
        .with("lot", 7.0)
        .unwrap()
}
fn interface(
    cards: &[MaterialCard; 2],
    axis: &str,
    dims: Dims,
    knots: [(f64, f64); 2],
    maximum_temperature: f64,
) -> InterfaceSystemCard {
    interface_with_static(cards, 2.0, axis, dims, knots, maximum_temperature)
}
fn interface_with_static(
    cards: &[MaterialCard; 2],
    static_mu: f64,
    axis: &str,
    dims: Dims,
    knots: [(f64, f64); 2],
    maximum_temperature: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (property, value, interpolation) in [
        (
            STATIC_FRICTION_PROPERTY,
            PropertyValue::Scalar {
                value: static_mu,
                dims: Dims::NONE,
            },
            InterpolationPolicy::ConstantWithinValidity,
        ),
        (
            KINETIC_FRICTION_PROPERTY,
            PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: dims,
                knots: knots.to_vec(),
                dims: Dims::NONE,
            },
            InterpolationPolicy::LinearInside,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(property, Dims::NONE),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "Ta",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        200.0,
                        maximum_temperature,
                    )
                    .with_quantity(
                        "Tb",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        200.0,
                        maximum_temperature,
                    )
                    .with_quantity(
                        "pressure",
                        QuantitySpec::dimensional(Pressure::DIMS),
                        0.0,
                        2e5,
                    )
                    .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0, 2.0)
                    .with_quantity("distance", QuantitySpec::dimensional(LENGTH), 0.0, 10.0)
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic prescribed dry sliding fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "A finish".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "B finish".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "fixture atmosphere".into(),
            history: "original ordered sliding surfaces".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn binding<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    card: &'a InterfaceSystemCard,
    point: &'a QueryPoint,
    partition: HeatPartition,
) -> SlidingInterface<'a> {
    SlidingInterface {
        point: ThermalInterfacePoint {
            interface: card,
            materials: [&cards[0], &cards[1]],
            context: point,
            temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        caloric: [curves[0].curve(), curves[1].curve()],
        cells: CELLS,
        nominal_area_m2: 0.001,
        pressure_axis: "pressure",
        speed_axis: "speed",
        accumulated_distance_axis: Some("distance"),
        tangent_frame_id: "fixed specimen tangent axes",
        selection: MaterialPropertySelection::SingleClaimOnly,
        partition: SlidingHeatPartition::Prescribed(partition),
        outside_destination: Some("separately resolved external heat sink"),
    }
}

#[test]
fn endpoint_friction_and_finite_temperatures_match_the_independent_two_body_solution() {
    let point = context();
    let mut forces = Vec::new();
    for cp_a in [2.0, 10.0] {
        let (a, ca) = common::material(1.0, cp_a, 0.0);
        let (b, cb) = common::material(1.0, 3.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let thermal = common::initial(&curves, 1.0, None);
        let original = thermal.cells().clone();
        // mu_k(Ta) = 0.2 + 0.001*(Ta - 300 K).
        let card = interface(
            &cards,
            "Ta",
            Temperature::DIMS,
            [(200.0, 0.1), (2000.0, 1.9)],
            2000.0,
        );
        let contact = MaterialSlidingContact::new(
            &thermal,
            binding(
                &cards,
                &curves,
                &card,
                &point,
                HeatPartition::new(0.25, 0.55, 0.20).unwrap(),
            ),
            100.0,
            1.0,
            InputAuthority::SyntheticFixture,
        )
        .unwrap();
        let next = contact
            .step_candidate(
                &thermal,
                0.5,
                SlidingDrive {
                    normal_load_n: 100.0,
                    relative_increment_m: [0.3, 0.4],
                },
                controls(),
                &[0.0; 2],
                |_| {
                    Ok::<_, Infallible>(StepCoefficients {
                        sources_w: BTreeMap::from([(CELLS[0], 1.0)]),
                        conductances_w_k: BTreeMap::from([(7, 4.0)]),
                        payload: "original transport callback",
                    })
                },
                || false,
            )
            .unwrap();
        // Solve the physical 2x2 heat balance directly, independent of Newton.
        // H=N*distance=50 J, g*dt=2 J/K, alpha=0.001/K, external heat=0.5 J.
        let aa = cp_a + 2.0 - 50.0 * 0.25 * 0.001;
        let ab = -2.0;
        let ba = -2.0 - 50.0 * 0.55 * 0.001;
        let bb = 5.0;
        let rhs_a = 50.0 * 0.25 * 0.2 + 0.5;
        let rhs_b = 50.0 * 0.55 * 0.2;
        let determinant = aa * bb - ab * ba;
        let ta = 300.0 + (rhs_a * bb - ab * rhs_b) / determinant;
        let tb = 300.0 + (aa * rhs_b - ba * rhs_a) / determinant;
        let mu = 0.2 + 0.001 * (ta - 300.0);
        let work = 50.0 * mu;
        let r = next.report();
        close(next.thermal().temperature_k(CELLS[0]).unwrap(), ta);
        close(next.thermal().temperature_k(CELLS[1]).unwrap(), tb);
        close(r.friction.traction_n[0], 100.0 * mu * 0.6);
        close(r.friction.traction_n[1], 100.0 * mu * 0.8);
        close(r.drive_work_j, work);
        close(r.deposited_heat_j, 0.8 * work);
        close(r.partition.other_work_j(), 0.2 * work);
        close(r.external_heat_j, 0.5);
        close(
            r.thermal.stored_energy_change_j,
            cp_a * (ta - 300.0) + 3.0 * (tb - 300.0),
        );
        close(r.total_energy_residual_j, 0.0);
        // Omitting the outside share or crediting A twice would visibly fail.
        assert!(
            (r.thermal.stored_energy_change_j - r.drive_work_j - r.external_heat_j).abs() > 1.0
        );
        assert!(r.total_energy_residual_j + r.partition.surface_a_heat_j() > 1.0);
        assert_eq!(r.source, "original transport callback");
        assert_eq!(r.material.resolved.card_identity(), card.content_hash());
        assert_eq!(r.material.resolved.state_coordinate("lot"), Some(7.0));
        close(r.material.resolved.state_coordinate("Ta").unwrap(), ta);
        close(r.material.resolved.state_coordinate("Tb").unwrap(), tb);
        assert_eq!(r.material.pressure_pa, 1e5);
        close(r.material.speed_mps, 1.0);
        assert_eq!(
            contact.interface().texture_frame_ids(),
            ["A finish", "B finish"]
        );
        assert_eq!(
            contact.binding().selection,
            next.contact().binding().selection
        );
        assert_eq!(thermal.cells(), &original);
        assert_eq!(contact.state().slip_m(), [0.0; 2]);
        assert_eq!(contact.work().dissipated_work_j(), 0.0);
        close(next.contact().work().dissipated_work_j(), work);
        forces.push(r.friction.traction_n[0]);
    }
    assert!(forces[0] > forces[1]); // Heat-capacity data alone changes friction force.
}

#[test]
fn actual_load_area_and_slip_speed_drive_the_selected_source_curves() {
    let (a, ca) = common::material(1.0, 10.0, 0.0);
    let (b, cb) = common::material(1.0, 20.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = common::initial(&curves, 0.0, None);
    for (axis, dims, knots) in [
        ("pressure", Pressure::DIMS, [(0.0, 0.1), (2e5, 0.5)]),
        ("speed", SPEED, [(0.0, 0.1), (2.0, 0.5)]),
    ] {
        let card = interface(&cards, axis, dims, knots, 2000.0);
        let contact = MaterialSlidingContact::new(
            &thermal,
            binding(
                &cards,
                &curves,
                &card,
                &point,
                HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
            ),
            50.0,
            0.5,
            InputAuthority::SyntheticFixture,
        )
        .unwrap();
        for (load, displacement) in [(50.0, 0.25), (100.0, 0.50)] {
            let next = contact
                .step_candidate(
                    &thermal,
                    0.5,
                    SlidingDrive {
                        normal_load_n: load,
                        relative_increment_m: [displacement, 0.0],
                    },
                    controls(),
                    &[0.0; 2],
                    other,
                    || false,
                )
                .unwrap();
            let coordinate = if axis == "pressure" {
                load / 0.001
            } else {
                displacement / 0.5
            };
            let mu = 0.1 + 0.4 * coordinate / knots[1].0;
            close(next.report().friction.traction_n[0], load * mu);
            close(next.report().drive_work_j, load * mu * displacement);
            close(
                next.thermal().temperature_k(CELLS[0]).unwrap(),
                300.0 + 0.5 * load * mu * displacement / 10.0,
            );
            close(
                next.thermal().temperature_k(CELLS[1]).unwrap(),
                300.0 + 0.5 * load * mu * displacement / 20.0,
            );
            close(
                next.report()
                    .material
                    .resolved
                    .state_coordinate(axis)
                    .unwrap(),
                coordinate,
            );
        }
    }
}

#[test]
fn reversed_sliding_retains_original_source_and_advances_history_only_once() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = interface(
        &cards,
        "distance",
        LENGTH,
        [(0.0, 0.1), (10.0, 0.3)],
        2000.0,
    );
    let mut thermal = common::initial(&curves, 0.2, None);
    let mut contact = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &card,
            &point,
            HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        ),
        20.0,
        0.5,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let pins = contact.binding().selection.clone();
    for (index, increment) in [[0.3, 0.4], [-0.3, -0.4]].into_iter().enumerate() {
        let drive = SlidingDrive {
            normal_load_n: 20.0,
            relative_increment_m: increment,
        };
        let next = contact
            .step_candidate(
                &thermal,
                index as f64 + 1.0,
                drive,
                controls(),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        let replay = contact
            .step_candidate(
                &thermal,
                index as f64 + 1.0,
                drive,
                controls(),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        assert_eq!(next.thermal().cells(), replay.thermal().cells());
        assert_eq!(next.contact().state(), replay.contact().state());
        assert_eq!(next.report().friction, replay.report().friction);
        let expected_distance = 0.5 * (index + 1) as f64;
        close(
            next.report().material.accumulated_distance_m,
            expected_distance,
        );
        close(
            next.report().drive_work_j,
            10.0 * (0.1 + 0.02 * expected_distance),
        );
        assert!(next.report().friction_trials > 1);
        assert_eq!(next.contact().binding().selection, pins);
        (thermal, contact, _) = next.into_parts();
    }
    assert_eq!(contact.state().slip_m(), [0.0, 0.0]);
    close(contact.accumulated_distance_m(), 1.0);
    close(contact.work().dissipated_work_j(), 2.3);
    close(
        contact.work().surface_a_heat_j() + contact.work().surface_b_heat_j(),
        2.3,
    );
    assert_eq!(contact.time_s(), thermal.time_s());
}

#[test]
fn final_cancellation_source_domain_and_trial_budget_leave_both_histories_unchanged() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = common::initial(&curves, 0.2, None);
    let old = thermal.cells().clone();
    let card = interface(
        &cards,
        "Ta",
        Temperature::DIMS,
        [(200.0, 0.1), (2000.0, 1.9)],
        2000.0,
    );
    let contact = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &card,
            &point,
            HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        ),
        20.0,
        0.5,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let drive = SlidingDrive {
        normal_load_n: 20.0,
        relative_increment_m: [0.3, 0.4],
    };
    let polls = Cell::new(0);
    let good = contact
        .step_candidate(&thermal, 1.0, drive, controls(), &[0.0; 2], other, || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let final_poll = polls.get();
    polls.set(0);
    assert!(matches!(
        contact.step_candidate(&thermal, 1.0, drive, controls(), &[0.0; 2], other, || {
            polls.set(polls.get() + 1);
            polls.get() == final_poll
        }),
        Err(SlidingError::Cancelled)
    ));
    let mut limited = controls();
    limited.maximum_friction_trials = good.report().friction_trials - 1;
    assert!(matches!(
        contact.step_candidate(&thermal, 1.0, drive, limited, &[0.0; 2], other, || false),
        Err(SlidingError::FrictionBudget)
    ));
    let refused = contact.step_candidate(
        &thermal,
        1.0,
        drive,
        controls(),
        &[0.0; 2],
        |state| {
            if state.temperature_k(CELLS[0]).unwrap() > 300.005 {
                return Err("independent transport refused");
            }
            Ok(StepCoefficients {
                sources_w: BTreeMap::new(),
                conductances_w_k: BTreeMap::new(),
                payload: (),
            })
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(SlidingError::Source("independent transport refused"))
    ));
    let narrow = interface(
        &cards,
        "Ta",
        Temperature::DIMS,
        [(200.0, 0.1), (2000.0, 1.9)],
        300.01,
    );
    let bounded = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &narrow,
            &point,
            HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        ),
        20.0,
        0.5,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    assert!(matches!(
        bounded.step_candidate(&thermal, 1.0, drive, controls(), &[0.0; 2], other, || false),
        Err(SlidingError::Material(_))
    ));
    let mut unowned = binding(
        &cards,
        &curves,
        &card,
        &point,
        HeatPartition::new(0.4, 0.4, 0.2).unwrap(),
    );
    unowned.outside_destination = None;
    assert!(matches!(
        MaterialSlidingContact::new(
            &thermal,
            unowned,
            20.0,
            0.5,
            InputAuthority::SyntheticFixture
        ),
        Err(SlidingError::Invalid(_))
    ));
    let mut foreign = binding(
        &cards,
        &curves,
        &card,
        &point,
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
    );
    foreign.caloric.swap(0, 1);
    assert!(matches!(
        MaterialSlidingContact::new(
            &thermal,
            foreign,
            20.0,
            0.5,
            InputAuthority::SyntheticFixture
        ),
        Err(SlidingError::Invalid(_))
    ));
    assert_eq!(thermal.cells(), &old);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(contact.state().slip_m(), [0.0; 2]);
    assert_eq!(contact.time_s(), 0.0);
    assert_eq!(contact.accumulated_distance_m(), 0.0);
    assert_eq!(contact.work().dissipated_work_j(), 0.0);
}

#[path = "tribothermal/partial_slip.rs"]
mod partial_slip;

#[path = "tribothermal/normal.rs"]
mod normal;
