//! G1/G3/G4: independent motor/thermal balances, finite source feedback and rollback.
#[path = "circuit_enthalpy/common.rs"]
mod common;
#[path = "motor_enthalpy/voltage.rs"]
mod voltage;
#[path = "motor_enthalpy/heat.rs"]
mod heat;

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve,
    material::{ConductionLink, MaterialCell},
    motor::*,
};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    phase::ResolvedHeatCapacityCurve,
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_mbd::rotordynamics::{
    RotorBearingModel,
    transient::{
        RotorTransientState,
        drive::{RotorCurrentDrive, RotorDriveError, RotorDrivePowerControl},
    },
};
use fs_thermal::enthalpy::*;
use std::cell::Cell;

const WINDING: CellId = CellId(10);
const BRIDGE: CellId = CellId(20);
const R0: f64 = 0.5;
const ALPHA: f64 = 0.0039;
const DT: f64 = 0.1;

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-7,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn controls() -> CircuitEnthalpySolve {
    let mut controls = common::config();
    controls.difference_step_j = 1e-5;
    controls.newton.absolute_tolerance = 1e-11;
    controls.cell_energy_tolerance_j = 1e-9;
    controls.thermal_energy_tolerance_j = 1e-8;
    controls.total_energy_tolerance_j = 1e-8;
    controls
}
fn power(input: f64) -> RotorDrivePowerControl {
    RotorDrivePowerControl {
        maximum_source_input_w: input,
        maximum_source_return_w: 1000.0,
        absolute_power_tolerance_w: 1e-9,
        absolute_voltage_tolerance_v: 1e-10,
    }
}
fn drive() -> RotorCurrentDrive {
    RotorCurrentDrive {
        shaft_coordinate: 0,
        torque_constant_nm_per_a: 1.0,
        winding_resistance_ohm: 999.0, // Replaced by the actual material source.
        maximum_current_a: 10.0,
        maximum_terminal_voltage_v: 100.0,
        minimum_speed_rad_s: -100.0,
        maximum_speed_rad_s: 100.0,
        motoring_efficiency: 0.8,
        regenerative_efficiency: 0.6,
    }
}
fn winding() -> MotorWinding {
    MotorWinding {
        cell: WINDING,
        basis: ElectricalPropertyBasis::Resistivity,
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }
}
fn state(speed: f64) -> RotorTransientState {
    RotorTransientState {
        displacement: vec![0.0],
        velocity: vec![speed],
    }
}
fn sources<const N: usize>(
    cp: [f64; N],
    alpha: f64,
) -> ([MaterialCard; N], [ResolvedHeatCapacityCurve; N]) {
    let pairs: [_; N] = std::array::from_fn(|i| common::material(R0, cp[i], alpha));
    (
        pairs.each_ref().map(|(card, _)| card.clone()),
        pairs.map(|(_, curve)| curve),
    )
}
fn body<'a, const N: usize>(
    curves: &'a [ResolvedHeatCapacityCurve; N],
    cp: [f64; N],
    linked: bool,
    winding_cooling_w_k: f64,
) -> DepositionSimulation<'a> {
    // Ambient cooling acts on the exposed surface only: the winding's 6 m2
    // less the 1 m2 the contact covers when the pair is linked
    // (`DepositionSimulation::exposed_area_m2`). The coefficient is taken over
    // that area so the winding loses `winding_cooling_w_k` per kelvin.
    let winding_exposed_m2 = if linked { 5.0 } else { 6.0 };
    let mut events = (0..N)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: [WINDING, BRIDGE][i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: 100.0 * cp[i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: if i == 0 {
                        winding_cooling_w_k / winding_exposed_m2
                    } else {
                        0.0
                    },
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect::<Vec<_>>();
    if linked {
        assert_eq!(N, 2);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: WINDING,
                b: BRIDGE,
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
    }
    let laws = curves
        .iter()
        .enumerate()
        .map(|(i, c)| (MaterialId(i as u64), c.curve() as &dyn ThermalLaw))
        .collect::<Vec<_>>();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn bindings<'a, const N: usize>(
    cards: &'a [MaterialCard; N],
    curves: &'a [ResolvedHeatCapacityCurve; N],
    context: &'a QueryPoint,
    conductivity: Option<&'a PropertyKey>,
) -> Vec<MaterialCell<'a>> {
    (0..N)
        .map(|i| MaterialCell {
            cell: [WINDING, BRIDGE][i],
            material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: conductivity.map(|key| (key, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect()
}

#[test]
fn held_current_temperature_feedback_matches_closed_form_and_continues() {
    let cp = [1.0, 0.5];
    let (cards, curves) = sources(cp, ALPHA);
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: winding(),
        bridge_cell: BRIDGE,
    };
    let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let native = rotor
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(DT, 1e-12, || false)
        .unwrap();
    let mut thermal = body(&curves, cp, false, 0.0);
    let mut motor = network
        .bind(&thermal, &native, state(0.0), drive(), || false)
        .unwrap();
    near(motor.parameters().winding_resistance_ohm, R0);
    let (mut winding_rise, mut bridge_rise) = (0.0, 0.0);
    for index in 0..2 {
        let prior_cells = thermal.cells().clone();
        let prior_motion = motor.state().clone();
        let current = 2.0;
        let delta_t = DT * current * current * R0 * (1.0 + ALPHA * winding_rise)
            / (cp[0] - DT * current * current * R0 * ALPHA);
        winding_rise += delta_t;
        let winding_heat = cp[0] * delta_t;
        let speed = index as f64 * DT * current / 2.0;
        let expected_speed = speed + DT * current / 2.0;
        let shaft = DT * current * 0.5 * (speed + expected_speed);
        let bridge_heat = (winding_heat + shaft) * (1.0 / drive().motoring_efficiency - 1.0);
        bridge_rise += bridge_heat / cp[1];
        let result = network
            .step_candidate(
                &thermal,
                (index + 1) as f64 * DT,
                &motor,
                &[0.0],
                current,
                power(1000.0),
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        near(
            result.thermal().temperature_k(WINDING).unwrap(),
            300.0 + winding_rise,
        );
        near(
            result.thermal().temperature_k(BRIDGE).unwrap(),
            300.0 + bridge_rise,
        );
        near(result.motor().state().velocity[0], expected_speed);
        near(
            result.motor().parameters().winding_resistance_ohm,
            R0 * (1.0 + ALPHA * winding_rise),
        );
        let step = &result.report().electrical;
        near(step.drive.current_a, current);
        near(step.drive.shaft_work_j, shaft);
        near(step.drive.winding_heat_j, winding_heat);
        near(step.drive.bridge_heat_j, bridge_heat);
        near(
            step.drive.electrical_draw_j,
            (winding_heat + shaft) / drive().motoring_efficiency,
        );
        near(step.drive.electrical_return_j, 0.0);
        near(
            step.supplied_to_kinetic_and_heat_j,
            step.drive.electrical_draw_j,
        );
        near(
            result.report().deposited_joule_j,
            winding_heat + bridge_heat,
        );
        near(result.report().total_energy_residual_j, 0.0);
        assert_eq!(
            result.motor().winding().material().card_identity(),
            motor.winding().material().card_identity()
        );
        assert_eq!(thermal.cells(), &prior_cells);
        assert_eq!(motor.state(), &prior_motion);
        (thermal, motor, _) = result.into_parts();
    }
}

#[test]
fn finite_source_derates_current_inside_temperature_solve_with_shared_heat() {
    let cp = [0.02];
    let (cards, curves) = sources(cp, ALPHA);
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: winding(),
        bridge_cell: WINDING,
    };
    let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let native = rotor
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(DT, 1e-12, || false)
        .unwrap();
    let thermal = body(&curves, cp, false, 0.0);
    let motor = network
        .bind(&thermal, &native, state(0.0), drive(), || false)
        .unwrap();
    let result = network
        .step_candidate(
            &thermal,
            DT,
            &motor,
            &[0.0],
            4.0,
            power(2.0),
            controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    // For free inertia J=2, endpoint source power is I²(R+h/J)/eta.
    // Shared heat is h I² R/eta + (1/eta-1) h² I²/(2J), independently
    // reduced here before solving the scalar caloric balance by bisection.
    let inertia_term = DT / 2.0;
    let eta = drive().motoring_efficiency;
    let (mut low, mut high) = (0.0, 20.0);
    for _ in 0..80 {
        let theta = 0.5 * (low + high);
        let r = R0 * (1.0 + ALPHA * theta);
        let heat = DT * 2.0 * (r + 0.5 * (1.0 - eta) * inertia_term) / (r + inertia_term);
        if cp[0] * theta > heat {
            high = theta;
        } else {
            low = theta;
        }
    }
    let theta = 0.5 * (low + high);
    let resistance = R0 * (1.0 + ALPHA * theta);
    let current = (eta * 2.0 / (resistance + inertia_term)).sqrt();
    let cold_current = (eta * 2.0 / (R0 + inertia_term)).sqrt();
    let step = &result.report().electrical.drive;
    assert!(step.current_a < cold_current);
    near(step.current_a, current);
    near(step.current_fraction, current / 4.0);
    near(step.final_source_power_w, 2.0);
    near(result.motor().state().velocity[0], DT * current / 2.0);
    near(
        result.thermal().temperature_k(WINDING).unwrap(),
        300.0 + theta,
    );
    near(result.report().deposited_joule_j, cp[0] * theta);
    near(step.winding_heat_j, DT * current * current * resistance);
    near(
        result.report().deposited_joule_j,
        step.winding_heat_j + step.bridge_heat_j,
    );
    near(step.electrical_draw_j, step.shaft_work_j + cp[0] * theta);
    near(result.report().total_energy_residual_j, 0.0);
}

#[test]
fn regeneration_retains_source_return_and_explicit_dissipative_braking() {
    let cp = [0.2, 0.3];
    let (cards, curves) = sources(cp, ALPHA);
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: winding(),
        bridge_cell: BRIDGE,
    };
    let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let native = rotor
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(DT, 1e-12, || false)
        .unwrap();
    let thermal = body(&curves, cp, false, 0.0);
    for eta in [0.6, 0.0] {
        let motor = network
            .bind(
                &thermal,
                &native,
                state(10.0),
                RotorCurrentDrive {
                    regenerative_efficiency: eta,
                    ..drive()
                },
                || false,
            )
            .unwrap();
        let result = network
            .step_candidate(
                &thermal,
                DT,
                &motor,
                &[0.0],
                -1.0,
                power(0.0),
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let theta = DT * R0 / (cp[0] - DT * R0 * ALPHA);
        let winding_heat = cp[0] * theta;
        let shaft = -DT * 0.5 * (10.0 + 9.95);
        let return_energy = -(shaft + winding_heat) * eta;
        let bridge_heat = -(shaft + winding_heat) * (1.0 - eta);
        let step = &result.report().electrical.drive;
        near(step.motion.state_after.velocity[0], 9.95);
        near(step.winding_heat_j, winding_heat);
        near(step.bridge_heat_j, bridge_heat);
        near(step.electrical_draw_j, 0.0);
        near(step.electrical_return_j, return_energy);
        near(
            step.motion.kinetic_energy_change_j,
            9.95_f64.powi(2) - 100.0,
        );
        near(
            result.thermal().temperature_k(WINDING).unwrap(),
            300.0 + theta,
        );
        near(
            result.thermal().temperature_k(BRIDGE).unwrap(),
            300.0 + bridge_heat / cp[1],
        );
        near(
            step.motion.kinetic_energy_change_j + result.report().deposited_joule_j + return_energy,
            0.0,
        );
        near(result.report().total_energy_residual_j, 0.0);
    }
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn operator_work_stays_separate_while_native_conduction_and_cooling_remove_heat() {
    let cp = [1.0, 2.0];
    let (cards, curves) = sources(cp, 0.0);
    let context = common::point(300.0);
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, Some(&conductivity));
    let links = [ConductionLink {
        link: 7,
        distance_a_m: 0.5,
        distance_b_m: 0.5,
    }];
    let network = MotorThermalNetwork {
        cells: &cells,
        links: &links,
        winding: winding(),
        bridge_cell: BRIDGE,
    };
    let rotor = RotorBearingModel::new(1, &[2.0], &[0.7], &[3.0], &[0.0], || false).unwrap();
    let native = rotor
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(DT, 1e-12, || false)
        .unwrap();
    let thermal = body(&curves, cp, true, 0.6);
    let initial = RotorTransientState {
        displacement: vec![0.2],
        velocity: vec![4.0],
    };
    let motor = network
        .bind(&thermal, &native, initial.clone(), drive(), || false)
        .unwrap();
    let result = network
        .step_candidate(
            &thermal,
            DT,
            &motor,
            &[0.4],
            1.0,
            power(1000.0),
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let step = &result.report().electrical;
    let drive = &step.drive;
    let mean_v = 0.5 * (initial.velocity[0] + drive.motion.state_after.velocity[0]);
    let mean_q = 0.5 * (initial.displacement[0] + drive.motion.state_after.displacement[0]);
    let operator = DT * (0.7 * mean_v * mean_v + 3.0 * mean_q * mean_v);
    near(step.operator_work_j, operator);
    near(drive.external_work_j, DT * 0.4 * mean_v);
    near(drive.shaft_work_j, DT * mean_v);
    let qw = drive.winding_heat_j;
    let qb = drive.bridge_heat_j;
    // Exact two-capacity backward-Euler thermal system with G12=1 and Ga=0.6.
    let a = cp[0] + DT * 1.6;
    let b = cp[1] + DT;
    let determinant = a * b - DT * DT;
    let theta_w = (b * qw + DT * qb) / determinant;
    let theta_b = (a * qb + DT * qw) / determinant;
    let cooling = DT * 0.6 * theta_w;
    near(
        result.thermal().temperature_k(WINDING).unwrap(),
        300.0 + theta_w,
    );
    near(
        result.thermal().temperature_k(BRIDGE).unwrap(),
        300.0 + theta_b,
    );
    near(result.report().deposited_joule_j, qw + qb);
    near(result.report().non_electrical_heat_j, -cooling);
    assert!(drive.motion.velocity_force_work_j > result.report().deposited_joule_j);
    near(cp[0] * theta_w + cp[1] * theta_b + cooling, qw + qb);
    near(
        step.supplied_to_kinetic_and_heat_j,
        drive.electrical_draw_j - drive.electrical_return_j + drive.external_work_j - operator,
    );
    near(
        drive.motion.kinetic_energy_change_j + cp[0] * theta_w + cp[1] * theta_b,
        step.supplied_to_kinetic_and_heat_j - cooling,
    );
    near(result.report().total_energy_residual_j, 0.0);
}

#[test]
fn source_domain_budget_clock_receiver_and_late_cancellation_refuse_atomically() {
    let cp = [1.0, 1.0];
    let (cards, curves) = sources(cp, ALPHA);
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: winding(),
        bridge_cell: BRIDGE,
    };
    let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let native = rotor
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(DT, 1e-12, || false)
        .unwrap();
    let thermal = body(&curves, cp, false, 0.0);
    let motor = network
        .bind(&thermal, &native, state(0.0), drive(), || false)
        .unwrap();
    let prior_cells = thermal.cells().clone();
    let prior_state = motor.state().clone();
    let calls = Cell::new(0usize);
    let accepted = network
        .step_candidate(
            &thermal,
            DT,
            &motor,
            &[0.0],
            2.0,
            power(1000.0),
            controls(),
            &[0.0; 2],
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    let last = calls.get();
    calls.set(0);
    assert!(matches!(
        network.step_candidate(
            &thermal,
            DT,
            &motor,
            &[0.0],
            2.0,
            power(1000.0),
            controls(),
            &[0.0; 2],
            || {
                calls.set(calls.get() + 1);
                calls.get() >= last
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut limited = controls();
    limited.maximum_circuit_trials = 1;
    assert!(matches!(
        network.step_candidate(
            &thermal,
            DT,
            &motor,
            &[0.0],
            2.0,
            power(1000.0),
            limited,
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert!(
        network
            .step_candidate(
                &thermal,
                2.0 * DT,
                &motor,
                &[0.0],
                2.0,
                power(1000.0),
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert!(
        network
            .step_candidate(
                accepted.thermal(),
                2.0 * DT,
                &motor,
                &[0.0],
                2.0,
                power(1000.0),
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert!(
        network
            .step_candidate(
                &thermal,
                DT,
                &motor,
                &[0.0],
                2.0,
                power(1000.0),
                controls(),
                &[1e6; 2],
                || false
            )
            .is_err()
    );
    assert!(matches!(
        network.step_candidate(
            &thermal,
            DT,
            &motor,
            &[0.0],
            f64::NAN,
            power(1000.0),
            controls(),
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::Source(MotorThermalError::Drive(
            RotorDriveError::Invalid(_)
        )))
    ));
    let changed_receiver = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: winding(),
        bridge_cell: WINDING,
    };
    assert!(
        changed_receiver
            .step_candidate(
                &thermal,
                DT,
                &motor,
                &[0.0],
                2.0,
                power(1000.0),
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    let bad_geometry = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: MotorWinding {
            area_m2: 2.0,
            ..winding()
        },
        bridge_cell: BRIDGE,
    };
    assert!(
        bad_geometry
            .bind(&thermal, &native, state(0.0), drive(), || false)
            .is_err()
    );
    let missing_bridge = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: winding(),
        bridge_cell: CellId(99),
    };
    assert!(
        missing_bridge
            .bind(&thermal, &native, state(0.0), drive(), || false)
            .is_err()
    );
    let wrong_point = common::point(300.0).with("lot", 8.0).unwrap();
    let wrong_cells = bindings(&cards, &curves, &wrong_point, None);
    let wrong_source = MotorThermalNetwork {
        cells: &wrong_cells,
        links: &[],
        winding: winding(),
        bridge_cell: BRIDGE,
    };
    assert!(
        wrong_source
            .step_candidate(
                &thermal,
                DT,
                &motor,
                &[0.0],
                2.0,
                power(1000.0),
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert_eq!(thermal.cells(), &prior_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(motor.state(), &prior_state);
    assert_eq!(motor.time_s(), 0.0);
    near(motor.parameters().winding_resistance_ohm, R0);
}
