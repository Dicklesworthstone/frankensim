//! Original induction flux states and two finite material winding temperatures.
#[path = "circuit_enthalpy/common.rs"]
mod common;
#[path = "induction_motor_enthalpy/heat.rs"]
mod heat;

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    material::{ConductionLink, MaterialCell},
    motor::{MotorWinding, induction::*},
};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    phase::ResolvedHeatCapacityCurve,
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::{
    SolveBudget,
    circuit::inductance::InductorFluxCurve,
    circuit::induction::{InductionMotor, InductionMotorParameters},
};
use fs_thermal::enthalpy::*;
use std::cell::Cell;

// Deliberately opposite to sorted thermal unknown order.
const STATOR: CellId = CellId(30);
const ROTOR: CellId = CellId(10);
const CP: [f64; 2] = [0.3, 0.2];
const R0: [f64; 2] = [0.5, 0.3];
const ALPHA: [f64; 2] = [0.0039, 0.0025];
const FACTORS: [f64; 2] = [1.0 / 3.0, 2.5];

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-7,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn parameters() -> InductionMotorParameters {
    InductionMotorParameters {
        stator_leakage_inductance_h: 0.004,
        rotor_leakage_inductance_h: 0.006,
        magnetizing_inductance_h: 0.05,
        stator_resistance_ohm: R0[0],
        rotor_resistance_ohm: R0[1],
        pole_pairs: 2,
        rotor_inertia_kg_m2: 0.1,
        stator_inertia_kg_m2: 0.2,
        bearing_drag_nm_s_per_rad: 0.0,
    }
}

// Authored reversible magnetization with a twenty-six-fold reduction in
// differential inductance above the knee. This is a numerical fixture, not
// material calibration. Both polarities belong to the same isotropic source.
fn magnetizing_curve() -> InductorFluxCurve {
    InductorFluxCurve::try_new(vec![
        (-0.3, -54.0),
        (-0.1, -2.0),
        (0.0, 0.0),
        (0.1, 2.0),
        (0.3, 54.0),
    ])
    .unwrap()
}

fn core_flux(current: f64) -> f64 {
    let magnitude = current.abs();
    current.signum()
        * if magnitude <= 2.0 {
            0.05 * magnitude
        } else {
            0.1 + (magnitude - 2.0) / 260.0
        }
}

// Independent energy in the original physical leakage/core coordinates.
// It never evaluates the production eliminated flux energy or its gradient.
fn saturated_energy(motor: &MaterialInductionMotor) -> f64 {
    let p = motor.parameters();
    let point = motor.operating_point();
    let stator = point.stator_current_alpha_beta_a[0];
    let rotor = point.rotor_current_dq_a[0];
    let psi = core_flux(stator + rotor).abs();
    let core = if psi <= 0.1 {
        10.0 * psi * psi
    } else {
        let increment = psi - 0.1;
        0.1 + 2.0 * increment + 130.0 * increment * increment
    };
    0.5 * (p.stator_leakage_inductance_h * stator * stator
        + p.rotor_leakage_inductance_h * rotor * rotor)
        + core
}
fn controls() -> InductionEnthalpySolve {
    let mut thermal = common::config();
    thermal.difference_step_j = 1e-5;
    thermal.newton.absolute_tolerance = 1e-11;
    thermal.cell_energy_tolerance_j = 1e-9;
    thermal.thermal_energy_tolerance_j = 1e-8;
    thermal.total_energy_tolerance_j = 1e-8;
    InductionEnthalpySolve {
        thermal,
        machine: SolveBudget {
            relative: 2e-12,
            absolute: 2e-13,
            max_iterations: 30,
        },
        backtracking_halvings: 8,
    }
}
fn sources() -> ([MaterialCard; 2], [ResolvedHeatCapacityCurve; 2]) {
    let pairs: [_; 2] =
        std::array::from_fn(|i| common::material(R0[i] / FACTORS[i], CP[i], ALPHA[i]));
    (
        pairs.each_ref().map(|(card, _)| card.clone()),
        pairs.map(|(_, curve)| curve),
    )
}
fn winding(index: usize) -> InductionWinding {
    InductionWinding {
        conductor: MotorWinding {
            cell: [STATOR, ROTOR][index],
            basis: ElectricalPropertyBasis::Resistivity,
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 1.0,
            area_m2: 1.0,
        },
        resistance_multiplier: FACTORS[index],
    }
}
fn bindings<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    conductivity: Option<&'a PropertyKey>,
) -> Vec<MaterialCell<'a>> {
    (0..2)
        .map(|i| MaterialCell {
            cell: [STATOR, ROTOR][i],
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
fn body<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    linked: bool,
    cooling: f64,
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: [STATOR, ROTOR][i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: 100.0 * CP[i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: if i == 0 { cooling / 6.0 } else { 0.0 },
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: STATOR,
                b: ROTOR,
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
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
            max_dt_s: 0.04,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

// Independent original two-coil RL equations. Collinear current and voltage
// preserve zero torque, so all native angular states stay zero. At proposed
// endpoint temperatures, the exact native midpoint electrical equations reduce
// to this explicit 2x2 inverse; no production flux gradient or solver is used.
fn electrical_reference(
    current: [f64; 2],
    rise: [f64; 2],
    voltage: f64,
    dt: f64,
) -> ([f64; 2], [f64; 2]) {
    let p = parameters();
    let ls = p.stator_leakage_inductance_h + p.magnetizing_inductance_h;
    let lr = p.rotor_leakage_inductance_h + p.magnetizing_inductance_h;
    let lm = p.magnetizing_inductance_h;
    let resistance: [f64; 2] = std::array::from_fn(|i| R0[i] * (1.0 + ALPHA[i] * rise[i]));
    let (a, d) = (ls + 0.5 * dt * resistance[0], lr + 0.5 * dt * resistance[1]);
    let rhs = [
        ls * current[0] + lm * current[1] + 0.5 * dt * voltage,
        lm * current[0] + lr * current[1],
    ];
    let determinant = a * d - lm * lm;
    let mean = [
        (d * rhs[0] - lm * rhs[1]) / determinant,
        (a * rhs[1] - lm * rhs[0]) / determinant,
    ];
    (
        std::array::from_fn(|i| 2.0 * mean[i] - current[i]),
        std::array::from_fn(|i| dt * resistance[i] * mean[i] * mean[i]),
    )
}

// The chosen weak-feedback manufactured source is contractive. Independent
// iteration closes two explicit caloric equations, including Fourier exchange
// and the external reservoir; the production owner instead uses Newton-Krylov.
fn reference(
    current: [f64; 2],
    initial_rise: [f64; 2],
    voltage: f64,
    dt: f64,
    conductance: f64,
    cooling: f64,
) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let (a, d, off) = (
        CP[0] + dt * (conductance + cooling),
        CP[1] + dt * conductance,
        dt * conductance,
    );
    let mut rise = initial_rise;
    let mut last_change = f64::INFINITY;
    for _ in 0..200 {
        let (_, heat) = electrical_reference(current, rise, voltage, dt);
        let rhs = [
            CP[0] * initial_rise[0] + heat[0],
            CP[1] * initial_rise[1] + heat[1],
        ];
        let next = [
            (d * rhs[0] + off * rhs[1]) / (a * d - off * off),
            (a * rhs[1] + off * rhs[0]) / (a * d - off * off),
        ];
        last_change = (next[0] - rise[0]).abs().max((next[1] - rise[1]).abs());
        rise = next;
        if last_change < 1e-13 {
            break;
        }
    }
    assert!(last_change < 1e-12);
    let (current, heat) = electrical_reference(current, rise, voltage, dt);
    (current, rise, heat)
}

#[test]
fn both_winding_temperatures_feed_native_induction_with_conduction_and_cooling() {
    let (cards, curves) = sources();
    let context = common::point(300.0);
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    for linked in [false, true] {
        let cooling = if linked { 0.15 } else { 0.0 };
        let cells = bindings(&cards, &curves, &context, linked.then_some(&conductivity));
        let links = if linked {
            vec![ConductionLink {
                link: 7,
                distance_a_m: 0.5,
                distance_b_m: 0.5,
            }]
        } else {
            vec![]
        };
        let network = InductionThermalNetwork {
            cells: &cells,
            links: &links,
            stator: winding(0),
            rotor: winding(1),
            bearing_cell: STATOR,
        };
        let mut thermal = body(&curves, linked, cooling);
        let mut motor = network
            .bind(
                &thermal,
                [0.0; 7],
                InductionMotorParameters {
                    stator_resistance_ohm: 999.0,
                    rotor_resistance_ohm: 999.0,
                    ..parameters()
                },
                || false,
            )
            .unwrap();
        near(motor.parameters().stator_resistance_ohm, R0[0]);
        near(motor.parameters().rotor_resistance_ohm, R0[1]);
        let (mut current, mut rise) = ([0.0; 2], [0.0; 2]);
        let dt = 0.04;
        for index in 0..3 {
            let old_cells = thermal.cells().clone();
            let old_state = *motor.state();
            let voltage = 8.0;
            let (next_current, next_rise, heat) = reference(
                current,
                rise,
                voltage,
                dt,
                if linked { 1.0 } else { 0.0 },
                cooling,
            );
            let candidate = network
                .step_candidate(
                    &thermal,
                    (index + 1) as f64 * dt,
                    &motor,
                    &[voltage, 0.0, 0.0, 0.0],
                    controls(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let step = &candidate.report().electrical;
            near(
                candidate
                    .motor()
                    .operating_point()
                    .stator_current_alpha_beta_a[0],
                next_current[0],
            );
            near(
                candidate.motor().operating_point().rotor_current_dq_a[0],
                next_current[1],
            );
            for i in 0..2 {
                near(
                    candidate
                        .thermal()
                        .temperature_k([STATOR, ROTOR][i])
                        .unwrap(),
                    300.0 + next_rise[i],
                );
                near(
                    candidate.report().source.resistance_ohm[i],
                    R0[i] * (1.0 + ALPHA[i] * next_rise[i]),
                );
            }
            near(step.stator_copper_heat_j, heat[0]);
            near(step.rotor_copper_heat_j, heat[1]);
            assert!(step.stator_copper_heat_j > 0.0 && step.rotor_copper_heat_j > 0.0);
            assert_eq!(step.bearing_heat_j, 0.0);
            near(candidate.report().deposited_joule_j, heat[0] + heat[1]);
            near(
                candidate.report().non_electrical_heat_j,
                -dt * cooling * next_rise[0],
            );
            near(candidate.report().total_energy_residual_j, 0.0);
            assert_eq!(candidate.motor().state()[4..], [0.0; 3]);
            assert!(candidate.report().circuit_trials > 1);
            assert_eq!(
                candidate
                    .motor()
                    .stator_winding()
                    .material()
                    .card_identity(),
                motor.stator_winding().material().card_identity()
            );
            assert_eq!(
                candidate.motor().rotor_winding().material().card_identity(),
                motor.rotor_winding().material().card_identity()
            );
            assert_eq!(thermal.cells(), &old_cells);
            assert_eq!(*motor.state(), old_state);
            (thermal, motor, _) = candidate.into_parts();
            current = next_current;
            rise = next_rise;
        }
        eprintln!(
            "induction thermal linked={linked}: stator={} K, rotor={} K",
            300.0 + rise[0],
            300.0 + rise[1]
        );
    }
}

fn energy(p: InductionMotorParameters, x: &[f64; 7]) -> f64 {
    let ls = p.stator_leakage_inductance_h + p.magnetizing_inductance_h;
    let lr = p.rotor_leakage_inductance_h + p.magnetizing_inductance_h;
    let lm = p.magnetizing_inductance_h;
    let (s, c) = (f64::from(p.pole_pairs) * x[4]).sin_cos();
    let r = [c * x[2] - s * x[3], s * x[2] + c * x[3]];
    0.5 * ((lr * (x[0] * x[0] + x[1] * x[1]) + ls * (x[2] * x[2] + x[3] * x[3])
        - 2.0 * lm * (x[0] * r[0] + x[1] * r[1]))
        / (ls * lr - lm * lm)
        + x[5] * x[5] / p.rotor_inertia_kg_m2
        + x[6] * x[6] / p.stator_inertia_kg_m2)
}

#[test]
fn saturation_retains_native_flux_energy_while_both_winding_temperatures_evolve() {
    let (cards, curves) = sources();
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = InductionThermalNetwork {
        cells: &cells,
        links: &[],
        stator: winding(0),
        rotor: winding(1),
        bearing_cell: STATOR,
    };
    let p = parameters();
    let curve = magnetizing_curve();
    let mut thermal = body(&curves, false, 0.0);
    // Direct physical linkage law, already above the saturation knee.
    let initial_current = [8.0, -0.8];
    let psi = core_flux(initial_current[0] + initial_current[1]);
    let state = [
        p.stator_leakage_inductance_h * initial_current[0] + psi,
        0.0,
        p.rotor_leakage_inductance_h * initial_current[1] + psi,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    let mut motor = network
        .bind_with_magnetizing_curve(&thermal, state, p, curve.clone(), || false)
        .unwrap();
    near(motor.operating_point().stator_current_alpha_beta_a[0], initial_current[0]);
    near(motor.operating_point().rotor_current_dq_a[0], initial_current[1]);
    let dt = 0.01;
    let voltage = 4.0;
    for interval in 1..=3 {
        let before = saturated_energy(&motor) + thermal.energy_j();
        let old_state = *motor.state();
        let old_cells = thermal.cells().clone();
        let next = network
            .step_candidate(
                &thermal,
                interval as f64 * dt,
                &motor,
                &[voltage, 0.0, 0.0, 0.0],
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let machine = next.motor();
        let point = machine.operating_point();
        let currents = [point.stator_current_alpha_beta_a[0], point.rotor_current_dq_a[0]];
        let psi = core_flux(currents[0] + currents[1]);
        assert!(psi.abs() > 0.1, "trajectory must exercise saturation");
        near(machine.state()[0], p.stator_leakage_inductance_h * currents[0] + psi);
        near(machine.state()[2], p.rotor_leakage_inductance_h * currents[1] + psi);
        assert_eq!(machine.magnetizing_curve(), Some(&curve));
        assert_eq!(machine.state()[4..], [0.0; 3]);
        let step = &next.report().electrical;
        let heat = [step.stator_copper_heat_j, step.rotor_copper_heat_j];
        for i in 0..2 {
            let id = [STATOR, ROTOR][i];
            let temperature = next.thermal().temperature_k(id).unwrap();
            let resistance = R0[i] * (1.0 + ALPHA[i] * (temperature - 300.0));
            near(next.report().source.resistance_ohm[i], resistance);
            let applied = if i == 0 { voltage } else { 0.0 };
            let mean_current = (applied - (machine.state()[2 * i] - old_state[2 * i]) / dt)
                / resistance;
            near(heat[i], dt * resistance * mean_current * mean_current);
            near(
                CP[i] * (temperature - thermal.temperature_k(id).unwrap()),
                heat[i],
            );
            assert!(heat[i] > 0.0 && temperature > 300.0);
        }
        let work: f64 = step.electrical_work_j.iter().sum();
        near(saturated_energy(machine) + next.thermal().energy_j() - before, work);
        near(next.report().deposited_joule_j, heat[0] + heat[1]);
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(motor.state(), &old_state);
        assert_eq!(thermal.cells(), &old_cells);
        (thermal, motor, _) = next.into_parts();
    }

    let old_state = *motor.state();
    let old_cells = thermal.cells().clone();
    // Every nonlinear probe stays in the original source support. An extreme
    // excitation must refuse transactionally rather than extrapolating a core.
    assert!(network.step_candidate(
        &thermal, thermal.time_s() + dt, &motor, &[1e6, 0.0, 0.0, 0.0],
        controls(), &[0.0; 2], || false,
    ).is_err());
    assert_eq!(motor.state(), &old_state);
    assert_eq!(thermal.cells(), &old_cells);
}

#[test]
fn generating_shaft_work_returns_electricity_while_both_windings_and_bearing_heat() {
    let (cards, curves) = sources();
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = InductionThermalNetwork {
        cells: &cells,
        links: &[],
        stator: winding(0),
        rotor: winding(1),
        bearing_cell: STATOR,
    };
    let thermal = body(&curves, false, 0.0);
    let p = InductionMotorParameters {
        bearing_drag_nm_s_per_rad: 0.002,
        ..parameters()
    };
    // Independent balanced cage phasor above synchronous speed: shorted rotor
    // current follows 0=(Rr+j*slip*Lr)*Ir+j*slip*Lm*Is, with Is=2 A.
    let (frequency, speed, current) = (40.0, 25.0, 2.0);
    let slip = frequency - f64::from(p.pole_pairs) * speed;
    let lr = p.rotor_leakage_inductance_h + p.magnetizing_inductance_h;
    let ls = p.stator_leakage_inductance_h + p.magnetizing_inductance_h;
    let lm = p.magnetizing_inductance_h;
    let denominator = R0[1] * R0[1] + (slip * lr).powi(2);
    let rotor = [
        -slip * slip * lm * lr * current / denominator,
        -slip * lm * R0[1] * current / denominator,
    ];
    let voltage = [
        R0[0] * current - frequency * lm * rotor[1],
        frequency * (ls * current + lm * rotor[0]),
    ];
    let torque = -f64::from(p.pole_pairs) * lm * rotor[1] * current;
    let state = InductionMotor::new(p)
        .unwrap()
        .state([current, 0.0], rotor, 0.0, [speed, 0.0])
        .unwrap();
    let motor = network.bind(&thermal, state, p, || false).unwrap();
    let load = -torque + p.bearing_drag_nm_s_per_rad * speed;
    let dt = 0.001;
    let (s, c) = (0.5 * frequency * dt).sin_cos();
    let u = [
        c * voltage[0] - s * voltage[1],
        s * voltage[0] + c * voltage[1],
        load,
        -load,
    ];
    let candidate = network
        .step_candidate(&thermal, dt, &motor, &u, controls(), &[0.0; 2], || false)
        .unwrap();
    let step = &candidate.report().electrical;
    let electrical_work: f64 = step.electrical_work_j.iter().sum();
    let shaft_work: f64 = step.shaft_work_j.iter().sum();
    assert!(electrical_work < 0.0 && shaft_work > 0.0);
    assert!(
        step.stator_copper_heat_j > 0.0
            && step.rotor_copper_heat_j > 0.0
            && step.bearing_heat_j > 0.0
    );
    near(
        CP[0] * (candidate.thermal().temperature_k(STATOR).unwrap() - 300.0),
        step.stator_copper_heat_j + step.bearing_heat_j,
    );
    near(
        CP[1] * (candidate.thermal().temperature_k(ROTOR).unwrap() - 300.0),
        step.rotor_copper_heat_j,
    );
    near(
        energy(p, candidate.motor().state()) - energy(p, &state)
            + candidate.report().deposited_joule_j,
        electrical_work + shaft_work,
    );
    near(
        candidate.motor().state()[5] + candidate.motor().state()[6],
        state[5] + state[6],
    );
    assert!(candidate.motor().parameters().stator_resistance_ohm > R0[0]);
    assert!(candidate.motor().parameters().rotor_resistance_ohm > R0[1]);
}

#[test]
fn source_correspondence_budgets_and_late_cancel_preserve_both_histories_and_retry() {
    let (cards, curves) = sources();
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = InductionThermalNetwork {
        cells: &cells,
        links: &[],
        stator: winding(0),
        rotor: winding(1),
        bearing_cell: STATOR,
    };
    let thermal = body(&curves, false, 0.0);
    let motor = network
        .bind(&thermal, [0.0; 7], parameters(), || false)
        .unwrap();
    let before_cells = thermal.cells().clone();
    let before_state = *motor.state();
    let before_parameters = motor.parameters();
    let calls = Cell::new(0);
    let inputs = [6.0, 0.0, 0.0, 0.0];
    let expected = network
        .step_candidate(
            &thermal,
            0.02,
            &motor,
            &inputs,
            controls(),
            &[0.0; 2],
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    for stop in [1, calls.get()] {
        let index = Cell::new(0);
        let result = network.step_candidate(
            &thermal,
            0.02,
            &motor,
            &inputs,
            controls(),
            &[0.0; 2],
            || {
                index.set(index.get() + 1);
                index.get() == stop
            },
        );
        assert!(matches!(result, Err(CircuitEnthalpyError::Cancelled)));
    }
    let mut limited = controls();
    limited.thermal.maximum_circuit_trials = 1;
    assert!(matches!(
        network.step_candidate(&thermal, 0.02, &motor, &inputs, limited, &[0.0; 2], || {
            false
        }),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut wrong = InductionThermalNetwork {
        cells: &cells,
        links: &[],
        stator: winding(0),
        rotor: winding(1),
        bearing_cell: STATOR,
    };
    wrong.rotor.resistance_multiplier *= 2.0;
    assert!(
        wrong
            .step_candidate(
                &thermal,
                0.02,
                &motor,
                &inputs,
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    wrong.rotor = winding(1);
    wrong.rotor.conductor.basis = ElectricalPropertyBasis::Conductivity;
    assert!(
        wrong
            .step_candidate(
                &thermal,
                0.02,
                &motor,
                &inputs,
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    wrong.rotor = winding(1);
    wrong.rotor.conductor.cell = STATOR;
    assert!(
        wrong
            .bind(&thermal, [0.0; 7], parameters(), || false)
            .is_err()
    );
    // This pulse would leave the original caloric/electrical source domain;
    // there is no frozen-resistance, clipping or post-failure heat fallback.
    assert!(
        network
            .step_candidate(
                &thermal,
                0.02,
                &motor,
                &[1e6, 0.0, 0.0, 0.0],
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &before_cells);
    assert_eq!(*motor.state(), before_state);
    assert_eq!(motor.parameters(), before_parameters);
    let replay = network
        .step_candidate(
            &thermal,
            0.02,
            &motor,
            &inputs,
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(replay.motor().state(), expected.motor().state());
    assert_eq!(replay.thermal().cells(), expected.thermal().cells());
    assert_eq!(replay.report().electrical, expected.report().electrical);
}
