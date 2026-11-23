//! The equivalent-circuit battery of fs-phs as one Component: closed forms
//! for one step and on the kernel, resistances that follow their
//! temperatures, conformance and refusals.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_phs::SolveBudget;
use fs_phs::circuit::battery::BatteryParameters;
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Carrier, Component, PortRole, Proposal, Receipt};
use fs_system::SystemBuilder;
use fs_system::battery::{BatteryComponentError, BatteryLoss, EquivalentCircuitBattery};
use fs_system::circuit::ResistanceLaw;
use fs_system::elements::{FlowSource, HeatCapacity, SourceAccount, Waveform};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 43,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} vs {expected}; tolerance {tolerance}"
    );
}

/// Series and polarization resistances \[ohm\].
const R: [f64; 2] = [0.2, 0.3];
/// Polarization and bus capacitances \[F\].
const C: [f64; 2] = [2.0, 0.5];
/// A flat 4 V curve.
const FLAT: [(f64, f64); 2] = [(0.0, 4.0), (1.0, 4.0)];

fn parameters() -> BatteryParameters {
    BatteryParameters {
        series_resistance_ohm: R[0],
        polarization_resistance_ohm: R[1],
        polarization_capacitance_f: C[0],
        bus_capacitance_f: C[1],
    }
}

fn battery(capacity_c: f64, knots: &[(f64, f64)]) -> EquivalentCircuitBattery {
    let solve = SolveBudget {
        relative: 2e-12,
        absolute: 2e-13,
        max_iterations: 30,
    };
    EquivalentCircuitBattery::new(capacity_c, knots, parameters(), solve, 16, 1e-8).unwrap()
}

/// One proposal; panics on anything but a whole step.
fn step(
    battery: &EquivalentCircuitBattery,
    old: &[f64; 3],
    dt_s: f64,
    inputs: [f64; 3],
) -> ([f64; 3], [f64; 3], Receipt) {
    let mut new = [0.0; 3];
    let mut outputs = [0.0; 3];
    let receipt =
        with_cx(
            |cx| match battery.propose(old, dt_s, &inputs, &mut new, &mut outputs, cx) {
                Ok(Proposal::Step(receipt)) => receipt,
                other => panic!("{other:?}"),
            },
        );
    (new, outputs, receipt)
}

/// `r0 (1 + alpha (T - 300 K))`, with no answer below `floor_k`.
#[derive(Debug)]
struct Linear {
    r0: f64,
    alpha: f64,
    floor_k: f64,
}

impl ResistanceLaw for Linear {
    fn ohms(&self, kelvin: f64) -> Result<f64, String> {
        if kelvin < self.floor_k {
            return Err(format!("{kelvin} K is below the card's {} K", self.floor_k));
        }
        Ok(self.r0 * (1.0 + self.alpha * (kelvin - 300.0)))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            self.r0.to_bits(),
            self.alpha.to_bits(),
            self.floor_k.to_bits(),
        ]);
    }
}

#[test]
fn the_ports_are_the_terminal_then_the_series_and_polarization_heat() {
    let ports: Vec<(String, PortRole, Carrier)> = battery(10.0, &FLAT)
        .ports()
        .iter()
        .map(|port| (port.id().as_str().to_owned(), port.role(), port.carrier()))
        .collect();
    assert_eq!(
        ports,
        [
            (
                "terminal".to_owned(),
                PortRole::FlowIn,
                Carrier::Power(PortKind::ElectricalVoltageCurrent)
            ),
            ("heat_series".to_owned(), PortRole::EffortIn, Carrier::Heat),
            (
                "heat_polarization".to_owned(),
                PortRole::EffortIn,
                Carrier::Heat
            ),
        ]
    );
}

#[test]
fn a_flat_curve_holds_its_stationary_state_while_discharging_and_charging() {
    // With a flat OCV `U`, polarization voltage `-I R1` and bus voltage
    // `U + I (R0 + R1)` are stationary under a terminal current `I`: only
    // the chemical charge moves, by `I dt`, and each resistance takes
    // `R I^2`.
    let battery = battery(10.0, &FLAT);
    let dt = 0.2;
    for current in [-1.0_f64, 1.0] {
        let bus_v = 4.0 + current * (R[0] + R[1]);
        let old = battery.state(0.5, -current * R[1], bus_v).unwrap();
        let (new, outputs, receipt) = step(&battery, &old, dt, [current, 300.0, 310.0]);
        close(new[0], old[0] + dt * current, 1e-12, "chemical charge");
        close(new[1], old[1], 1e-12, "polarization charge");
        close(new[2], old[2], 1e-12, "bus charge");
        close(outputs[0], bus_v, 1e-12, "terminal voltage");
        close(
            outputs[1],
            -R[0] * current * current,
            1e-12,
            "series heat rate",
        );
        close(
            outputs[2],
            -R[1] * current * current,
            1e-12,
            "polarization heat rate",
        );
        let heat = dt * current * current;
        close(receipt.dissipated_j, heat * (R[0] + R[1]), 1e-12, "heat");
        close(
            receipt.entropy_generated_j_per_k,
            heat * (R[0] / 300.0 + R[1] / 310.0),
            1e-14,
            "entropy",
        );
        close(
            receipt.stored_change_j,
            4.0 * dt * current,
            1e-12,
            "storage",
        );
        close(receipt.numerical_loss_j, 0.0, 1e-12, "numerical loss");
    }
}

#[test]
fn resistances_that_follow_their_temperatures_meet_the_hand_eliminated_midpoint_step() {
    // Both resistances follow their ports' temperatures. At fixed
    // temperatures the step is linear: the bus and polarization midpoint
    // equations eliminate to one equation for the cell current, as in the
    // module's own oracle.
    let alpha = 0.004;
    let battery = battery(10.0, &FLAT)
        .following(
            BatteryLoss::Series,
            Arc::new(Linear {
                r0: R[0],
                alpha,
                floor_k: 1.0,
            }),
        )
        .unwrap()
        .following(
            BatteryLoss::Polarization,
            Arc::new(Linear {
                r0: R[1],
                alpha,
                floor_k: 1.0,
            }),
        )
        .unwrap();
    let (dt, input, kelvin) = (0.2, -1.0, [310.0, 320.0]);
    let (vp0, vb0) = (0.3, 3.0);
    let old = battery.state(0.8, vp0, vb0).unwrap();
    let (new, outputs, receipt) = step(&battery, &old, dt, [input, kelvin[0], kelvin[1]]);

    let r: [f64; 2] = core::array::from_fn(|i| R[i] * (1.0 + alpha * (kelvin[i] - 300.0)));
    let a = 2.0 * C[0] + dt / r[1];
    let current = (4.0 - 2.0 * C[0] * vp0 / a - vb0 - dt * input / (2.0 * C[1]))
        / (r[0] + dt / a + dt / (2.0 * C[1]));
    let vp = (2.0 * C[0] * vp0 + dt * current) / a;
    let vb = vb0 + dt * (current + input) / (2.0 * C[1]);
    let heat = [dt * r[0] * current * current, dt * vp * vp / r[1]];
    close(new[0], old[0] - dt * current, 1e-12, "chemical charge");
    close(
        new[1],
        2.0 * C[0] * vp - old[1],
        1e-12,
        "polarization charge",
    );
    close(new[2], 2.0 * C[1] * vb - old[2], 1e-12, "bus charge");
    close(outputs[0], vb, 1e-12, "terminal voltage");
    close(outputs[1], -heat[0] / dt, 1e-11, "series heat rate");
    close(outputs[2], -heat[1] / dt, 1e-11, "polarization heat rate");
    close(
        receipt.entropy_generated_j_per_k,
        heat[0] / kelvin[0] + heat[1] / kelvin[1],
        1e-14,
        "entropy",
    );
    // The laws moved the resistances by 4 % and 8 %; the fixed ones give
    // another current.
    let fixed =
        (4.0 - 2.0 * C[0] * vp0 / (2.0 * C[0] + dt / R[1]) - vb0 - dt * input / (2.0 * C[1]))
            / (R[0] + dt / (2.0 * C[0] + dt / R[1]) + dt / (2.0 * C[1]));
    assert!((current - fixed).abs() > 1e-3, "{current} vs {fixed}");
}

#[test]
fn crossing_an_ocv_knot_books_the_curves_own_chemical_energy() {
    // Knots at 3.0, 3.0 and 4.2 V over 1 C: the chemical energy is
    // `3 q + 1.2 (q - 0.5)^2` above half charge and `3 q` below it. With no
    // terminal current the bus feeds the cell, and the step crosses the knot.
    let battery = battery(1.0, &[(0.0, 3.0), (0.5, 3.0), (1.0, 4.2)]);
    let energy = |x: &[f64; 3]| {
        3.0 * x[0]
            + 1.2 * (x[0] - 0.5).max(0.0).powi(2)
            + x[1] * x[1] / (2.0 * C[0])
            + x[2] * x[2] / (2.0 * C[1])
    };
    let old = battery.state(0.6, 0.2, 0.4).unwrap();
    close(
        battery.stored_energy_j(&old).unwrap(),
        energy(&old),
        1e-14,
        "stored energy",
    );
    let dt = 0.03;
    let (new, outputs, receipt) = step(&battery, &old, dt, [0.0, 300.0, 300.0]);
    assert!(old[0] > 0.5 && new[0] < 0.5, "{old:?} -> {new:?}");
    close(
        receipt.stored_change_j,
        energy(&new) - energy(&old),
        1e-12,
        "stored change",
    );
    let current = -(new[0] - old[0]) / dt;
    let vp = f64::midpoint(old[1], new[1]) / C[0];
    let heat = [dt * R[0] * current * current, dt * vp * vp / R[1]];
    close(-outputs[1] * dt, heat[0], 1e-12, "series heat");
    close(-outputs[2] * dt, heat[1], 1e-12, "polarization heat");
    // No terminal current: the heat is what the storage gave up.
    close(
        receipt.stored_change_j + receipt.dissipated_j,
        0.0,
        1e-12,
        "balance",
    );
}

#[test]
fn a_flat_curve_discharged_on_the_kernel_heats_two_bodies_and_closes_the_ledger() {
    // A current source draws 1 A for five steps of 0.2 s from the stationary
    // state, each resistance heating its own body.
    let battery = battery(10.0, &FLAT);
    let current = -1.0_f64;
    let bus_v = 4.0 + current * (R[0] + R[1]);
    let start = battery.state(0.5, -current * R[1], bus_v).unwrap();
    let capacities = [2.0, 3.0];
    let mut builder = SystemBuilder::new();
    let cell = builder.add("battery", battery, start).unwrap();
    // Positive flow enters the source, which therefore draws the current.
    let source = builder
        .add(
            "load",
            FlowSource::in_domain(
                PortKind::ElectricalVoltageCurrent,
                Waveform::constant(-current).unwrap(),
            ),
            SourceAccount::default(),
        )
        .unwrap();
    builder
        .common_flow(
            builder.port(source, "p").unwrap(),
            &[builder.port(cell, "terminal").unwrap()],
        )
        .unwrap();
    let mut bodies = Vec::new();
    for (k, (heat, capacity)) in ["heat_series", "heat_polarization"]
        .into_iter()
        .zip(capacities)
        .enumerate()
    {
        let body = builder
            .add(
                &format!("body{k}"),
                HeatCapacity::new(capacity).unwrap(),
                300.0,
            )
            .unwrap();
        builder
            .common_effort(
                builder.port(body, "heat").unwrap(),
                &[builder.port(cell, heat).unwrap()],
            )
            .unwrap();
        bodies.push(body);
    }
    let mut system = builder.seal().unwrap();
    for heat in ["heat_series", "heat_polarization"] {
        let port = system.port(cell, heat).unwrap();
        assert!(system.guess(port, &[300.0]));
    }
    let (dt, steps) = (0.2, 5);
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).unwrap();
        }
    });
    let elapsed = dt * f64::from(steps);
    let state = *system.state::<[f64; 3]>(cell).unwrap();
    close(
        state[0],
        start[0] + elapsed * current,
        1e-11,
        "chemical charge",
    );
    close(state[1], start[1], 1e-11, "polarization charge");
    close(state[2], start[2], 1e-11, "bus charge");
    for (k, &body) in bodies.iter().enumerate() {
        close(
            *system.state::<f64>(body).unwrap(),
            300.0 + elapsed * R[k] * current * current / capacities[k],
            1e-10,
            "body temperature",
        );
    }
    let account = system.state::<SourceAccount>(source).unwrap();
    close(account.energy_in_j, elapsed * bus_v, 1e-10, "load energy");
    close(
        system.ledger().energy_residual_j(),
        0.0,
        1e-10,
        "ledger closes",
    );
}

#[test]
fn the_battery_passes_the_conformance_kit_with_both_laws() {
    let law = |r0| {
        Arc::new(Linear {
            r0,
            alpha: 0.004,
            floor_k: 1.0,
        })
    };
    let battery = battery(10.0, &[(0.0, 3.2), (0.4, 3.6), (1.0, 4.1)])
        .following(BatteryLoss::Series, law(R[0]))
        .unwrap()
        .following(BatteryLoss::Polarization, law(R[1]))
        .unwrap();
    let state = battery.state(0.7, 0.1, 3.8).unwrap();
    let stale = battery.state(0.2, -0.4, 1.0).unwrap();
    let inputs = [-0.5, 330.0, 310.0];
    with_cx(|cx| {
        let case = Case {
            state: &state,
            stale: &stale,
            dt_s: 0.01,
            inputs: &inputs,
        };
        let observed = check_component(&battery, &case, cx).expect("component conformance");
        let Observed::Stepped { outputs, .. } = observed else {
            panic!("the battery shortened the step");
        };
        assert!(
            outputs[0] > 0.0 && outputs[1] < 0.0 && outputs[2] < 0.0,
            "{outputs:?}"
        );
    });
}

#[test]
fn refusals_publish_nothing() {
    let refused = |battery: &EquivalentCircuitBattery, old: [f64; 3], dt: f64, inputs: [f64; 3]| {
        let mut new = [7.0; 3];
        let mut outputs = [7.0; 3];
        let error = with_cx(|cx| battery.propose(&old, dt, &inputs, &mut new, &mut outputs, cx))
            .expect_err("the step is refused");
        assert_eq!((new, outputs), ([7.0; 3], [7.0; 3]), "{error:?}");
        error
    };
    let law = |r0, alpha, floor_k| Arc::new(Linear { r0, alpha, floor_k });
    let plain = battery(10.0, &FLAT);
    let old = plain.state(0.5, 0.0, 4.0).unwrap();

    // A law without an answer, and one answering a negative resistance,
    // each by its channel.
    let cold_card = battery(10.0, &FLAT)
        .following(BatteryLoss::Series, law(R[0], 0.0, 320.0))
        .unwrap();
    assert!(matches!(
        refused(&cold_card, old, 0.1, [-1.0, 310.0, 330.0]),
        BatteryComponentError::Law {
            loss: BatteryLoss::Series,
            ..
        }
    ));
    let negative = battery(10.0, &FLAT)
        .following(BatteryLoss::Polarization, law(R[1], -1.0, 1.0))
        .unwrap();
    assert!(matches!(
        refused(&negative, old, 0.1, [-1.0, 310.0, 330.0]),
        BatteryComponentError::Law {
            loss: BatteryLoss::Polarization,
            ..
        }
    ));

    // A port at no temperature, and a step of no length.
    assert!(matches!(
        refused(&plain, old, 0.1, [-1.0, 0.0, 300.0]),
        BatteryComponentError::Invalid { .. }
    ));
    assert!(matches!(
        refused(&plain, old, 0.0, [-1.0, 300.0, 300.0]),
        BatteryComponentError::Invalid { .. }
    ));

    // Drawing 20 A for 0.1 s takes 2 C from a cell holding 0.1 C: the owner
    // refuses to leave its curve.
    let low = plain.state(0.01, 0.0, 4.0).unwrap();
    assert!(matches!(
        refused(&plain, low, 0.1, [-20.0, 300.0, 300.0]),
        BatteryComponentError::Step(_)
    ));

    // Construction: a law given twice, no tolerance, a curve the owner
    // refuses.
    assert!(matches!(
        battery(10.0, &FLAT)
            .following(BatteryLoss::Series, law(R[0], 0.0, 1.0))
            .unwrap()
            .following(BatteryLoss::Series, law(R[0], 0.0, 1.0)),
        Err(BatteryComponentError::Invalid { .. })
    ));
    let solve = SolveBudget::default();
    assert!(matches!(
        EquivalentCircuitBattery::new(10.0, &FLAT, parameters(), solve, 16, 0.0),
        Err(BatteryComponentError::Invalid { .. })
    ));
    assert!(matches!(
        EquivalentCircuitBattery::new(
            10.0,
            &[(0.1, 4.0), (1.0, 4.0)],
            parameters(),
            solve,
            16,
            1e-8
        ),
        Err(BatteryComponentError::Model(_))
    ));
}

#[test]
fn the_state_codec_round_trips_and_the_identity_names_curve_and_laws() {
    let law = |alpha| {
        Arc::new(Linear {
            r0: R[0],
            alpha,
            floor_k: 1.0,
        })
    };
    let one = battery(10.0, &FLAT)
        .following(BatteryLoss::Series, law(0.004))
        .unwrap();
    let state = one.state(0.4, 0.1, 3.9).unwrap();
    let mut words = Vec::new();
    assert!(one.encode_state(&state, &mut words));
    assert_eq!(one.decode_state(&words), Some(state));
    assert_eq!(one.decode_state(&words[..2]), None);
    assert_eq!(one.decode_state(&[0, 0, f64::NAN.to_bits()]), None);

    let identity = |battery: &EquivalentCircuitBattery| {
        let mut words = Vec::new();
        battery.identity_words(&mut words);
        words
    };
    let other_law = battery(10.0, &FLAT)
        .following(BatteryLoss::Series, law(0.005))
        .unwrap();
    let other_loss = battery(10.0, &FLAT)
        .following(BatteryLoss::Polarization, law(0.004))
        .unwrap();
    let other_curve = battery(10.0, &[(0.0, 4.0), (1.0, 4.1)])
        .following(BatteryLoss::Series, law(0.004))
        .unwrap();
    for other in [&other_law, &other_loss, &other_curve] {
        assert_ne!(identity(&one), identity(other));
    }
}
