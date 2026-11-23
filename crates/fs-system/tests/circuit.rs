//! A Kirchhoff circuit of fs-phs as one Component: closed forms on the
//! kernel, a resistor that follows its temperature, resistors that follow
//! their operating point, conformance and refusals.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_phs::circuit::{Branch, CircuitGraph, TransformerLink};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, PortRole, Proposal};
use fs_system::circuit::{
    CircuitComponentError, KirchhoffCircuit, OperatingPointIteration, ResistanceLaw,
};
use fs_system::elements::HeatCapacity;
use fs_system::{NodeId, System, SystemBuilder};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 41,
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

/// A charged capacitor `C` across node 1 and ground, discharged through
/// `r1` from node 1 to node 2 and `r2` from node 2 to ground.
fn series_rc(farads: f64, r1: f64, r2: f64) -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads }),
            (1, 2, Branch::Resistor { ohms: r1 }),
            (2, 0, Branch::Resistor { ohms: r2 }),
        ],
        transformers: vec![],
    }
}

/// The circuit joined to one heat capacity per resistor, each at 300 K.
fn heated(
    circuit: KirchhoffCircuit,
    state: Vec<f64>,
    capacities: &[f64],
) -> (System, NodeId, Vec<NodeId>) {
    let heats: Vec<String> = circuit
        .ports()
        .iter()
        .map(|port| port.id().as_str().to_owned())
        .filter(|id| id.starts_with("heat"))
        .collect();
    assert_eq!(heats.len(), capacities.len());
    let mut builder = SystemBuilder::new();
    let node = builder.add("circuit", circuit, state).unwrap();
    let mut bodies = Vec::new();
    for (k, (heat, &capacity)) in heats.iter().zip(capacities).enumerate() {
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
                &[builder.port(node, heat).unwrap()],
            )
            .unwrap();
        bodies.push(body);
    }
    let mut system = builder.seal().unwrap();
    for heat in &heats {
        let port = system.port(node, heat).unwrap();
        assert!(system.guess(port, &[300.0]));
    }
    (system, node, bodies)
}

#[test]
fn the_ports_are_one_per_source_then_one_heat_port_per_resistor_by_branch() {
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 1 }),
            (1, 2, Branch::Resistor { ohms: 2.0 }),
            (2, 0, Branch::CurrentSource { port: 0 }),
            (2, 0, Branch::Resistor { ohms: 5.0 }),
        ],
        transformers: vec![],
    };
    let circuit = KirchhoffCircuit::new(graph, vec![], 1.0e-12).unwrap();
    let ports: Vec<(String, PortRole)> = circuit
        .ports()
        .iter()
        .map(|port| (port.id().as_str().to_owned(), port.role()))
        .collect();
    assert_eq!(
        ports,
        [
            ("source0".to_owned(), PortRole::FlowIn),
            ("source1".to_owned(), PortRole::EffortIn),
            ("heat1".to_owned(), PortRole::EffortIn),
            ("heat3".to_owned(), PortRole::EffortIn),
        ]
    );
}

#[test]
fn a_series_rc_discharges_on_the_midpoint_rule_and_each_resistor_heats_its_own_body() {
    let (farads, r1, r2, q0) = (1.0_f64, 1.0_f64, 3.0_f64, 10.0_f64);
    let circuit = KirchhoffCircuit::new(series_rc(farads, r1, r2), vec![], 1.0e-12).unwrap();
    let state = circuit.state(&[], &[q0]).unwrap();
    let capacities = [4.0_f64, 2.0_f64];
    let (mut system, node, bodies) = heated(circuit, state, &capacities);
    let (dt, steps) = (0.25_f64, 8_i32);
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).unwrap();
        }
    });
    // Midpoint rule on q' = -q / (R C): q_{n+1} = q_n (1 - a) / (1 + a).
    let a = dt / (2.0 * (r1 + r2) * farads);
    let ratio = (1.0 - a) / (1.0 + a);
    let charge = system.state::<Vec<f64>>(node).unwrap()[0];
    close(charge, q0 * ratio.powi(steps), 1.0e-13, "charge");
    // Each step's midpoint current flows through both resistors: their heats
    // stand as r1 : r2 at every step, and sum to the lost field energy.
    let lost = (q0 * q0 - charge * charge) / (2.0 * farads);
    let heat = [lost * r1 / (r1 + r2), lost * r2 / (r1 + r2)];
    let ledger = system.ledger();
    close(
        ledger.nodes[node.index()].dissipated_j.value(),
        lost,
        1.0e-12,
        "dissipated",
    );
    close(
        ledger.nodes[node.index()].stored_change_j.value(),
        -lost,
        1.0e-12,
        "stored",
    );
    for (k, body) in bodies.iter().enumerate() {
        let kelvin = *system.state::<f64>(*body).unwrap();
        close(
            capacities[k] * (kelvin - 300.0),
            heat[k],
            1.0e-11,
            "body heat",
        );
        close(
            ledger.nodes[body.index()].stored_change_j.value(),
            heat[k],
            1.0e-11,
            "body account",
        );
    }
    close(ledger.energy_residual_j(), 0.0, 1.0e-12, "ledger closes");
    assert!(ledger.entropy_generated_j_per_k() > 0.0);
}

#[test]
fn a_voltage_source_supplies_the_inductor_field_and_the_resistor_heat() {
    // V drives node 1; R from node 1 to node 2; L from node 2 to ground.
    let (volts, ohms, henries) = (12.0_f64, 2.0_f64, 0.5_f64);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms }),
            (2, 0, Branch::Inductor { henries }),
        ],
        transformers: vec![],
    };
    let circuit = KirchhoffCircuit::new(graph, vec![], 1.0e-12).unwrap();
    let state = circuit.state(&[0.0], &[]).unwrap();
    let mut builder = SystemBuilder::new();
    let node = builder.add("circuit", circuit, state).unwrap();
    let body = builder
        .add("body", HeatCapacity::new(50.0).unwrap(), 300.0)
        .unwrap();
    builder
        .common_effort(
            builder.port(body, "heat").unwrap(),
            &[builder.port(node, "heat1").unwrap()],
        )
        .unwrap();
    builder
        .boundary(
            builder.port(node, "source0").unwrap(),
            Box::new(move |_, _, out: &mut [f64]| out[0] = volts),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let port = system.port(node, "heat1").unwrap();
    assert!(system.guess(port, &[300.0]));
    let (dt, steps) = (0.01_f64, 200_usize);
    let mut current = 0.0_f64;
    let mut supplied = 0.0_f64;
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).unwrap();
            // Midpoint rule on L i' = V - R i.
            let a = ohms * dt / (2.0 * henries);
            let next = (current * (1.0 - a) + dt * volts / henries) / (1.0 + a);
            supplied += volts * 0.5 * (current + next) * dt;
            current = next;
        }
    });
    let flux = system.state::<Vec<f64>>(node).unwrap()[0];
    close(flux / henries, current, 1.0e-12, "current");
    let field = 0.5 * henries * current * current;
    let heat = supplied - field;
    let ledger = system.ledger();
    close(
        ledger.nodes[node.index()].stored_change_j.value(),
        field,
        1.0e-10,
        "field",
    );
    close(
        ledger.nodes[node.index()].dissipated_j.value(),
        heat,
        1.0e-10,
        "heat",
    );
    close(
        50.0 * (*system.state::<f64>(body).unwrap() - 300.0),
        heat,
        1.0e-9,
        "body",
    );
    close(ledger.energy_residual_j(), 0.0, 1.0e-10, "ledger closes");
}

/// `R(T) = r0 (1 + alpha (T - 300))`, refused below `floor_k`.
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
fn a_resistor_following_its_temperature_meets_the_hand_solved_midpoint_balance() {
    // One step of an RC whose single resistor follows the step-mean
    // temperature of the body it heats. The heat capacity returns
    // T0 + Q / (2 Cb), so the step is the fixed point of
    // Q = (q0^2 - q1(R(T))^2) / 2C, T = T0 + Q / (2 Cb).
    let (farads, r0, alpha, q0, cb, dt) = (0.5_f64, 2.0_f64, 0.02_f64, 4.0_f64, 0.05_f64, 0.3_f64);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads }),
            (1, 0, Branch::Resistor { ohms: r0 }),
        ],
        transformers: vec![],
    };
    let law = Arc::new(Linear {
        r0,
        alpha,
        floor_k: 1.0,
    });
    let circuit = KirchhoffCircuit::new(graph, vec![], 1.0e-12)
        .unwrap()
        .following(1, law)
        .unwrap();
    let state = circuit.state(&[], &[q0]).unwrap();
    let (mut system, node, bodies) = heated(circuit, state, &[cb]);
    with_cx(|cx| system.step(dt, cx).unwrap());

    let heat_at = |kelvin: f64| {
        let a = dt / (2.0 * r0 * (1.0 + alpha * (kelvin - 300.0)) * farads);
        let q1 = q0 * (1.0 - a) / (1.0 + a);
        ((q0 * q0 - q1 * q1) / (2.0 * farads), q1)
    };
    // Bisection on g(T) = T - 300 - Q(T) / 2Cb, which rises with T here.
    let (mut low, mut high) = (300.0_f64, 300.0 + q0 * q0 / (2.0 * farads) / cb);
    for _ in 0..200 {
        let middle = 0.5 * (low + high);
        if middle - 300.0 - heat_at(middle).0 / (2.0 * cb) > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    let mean = 0.5 * (low + high);
    let (heat, q1) = heat_at(mean);
    // The feedback is real: the resistance moved by several percent.
    assert!(alpha * (mean - 300.0) > 0.05, "{mean}");
    let charge = system.state::<Vec<f64>>(node).unwrap()[0];
    let kelvin = *system.state::<f64>(bodies[0]).unwrap();
    close(charge, q1, 1.0e-9, "charge");
    close(kelvin, 300.0 + heat / cb, 1.0e-7, "body temperature");
    close(
        system.ledger().energy_residual_j(),
        0.0,
        1.0e-12,
        "ledger closes",
    );
}

#[test]
fn the_circuit_passes_the_conformance_kit_with_a_source_and_a_law() {
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 3.0 }),
            (2, 0, Branch::Capacitor { farads: 0.2 }),
            (2, 0, Branch::Resistor { ohms: 7.0 }),
        ],
        transformers: vec![],
    };
    let law = Arc::new(Linear {
        r0: 3.0,
        alpha: 0.004,
        floor_k: 1.0,
    });
    let circuit = KirchhoffCircuit::new(graph, vec![], 1.0e-12)
        .unwrap()
        .following(1, law)
        .unwrap();
    let state = circuit.state(&[], &[0.3]).unwrap();
    let stale = circuit.state(&[], &[-5.0]).unwrap();
    let inputs = [9.0, 350.0, 310.0];
    with_cx(|cx| {
        let case = Case {
            state: &state,
            stale: &stale,
            dt_s: 0.01,
            inputs: &inputs,
        };
        let observed = check_component(&circuit, &case, cx).expect("component conformance");
        let Observed::Stepped { outputs, .. } = observed else {
            panic!("the circuit shortened the step");
        };
        // Current into the circuit, then minus each resistor's Joule power.
        assert!(
            outputs[0] > 0.0 && outputs[1] < 0.0 && outputs[2] < 0.0,
            "{outputs:?}"
        );
    });
}

#[test]
fn a_law_without_an_answer_refuses_the_step_by_branch() {
    let law = Arc::new(Linear {
        r0: 1.0,
        alpha: 0.0,
        floor_k: 320.0,
    });
    let circuit = KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12)
        .unwrap()
        .following(2, law)
        .unwrap();
    let old = circuit.state(&[], &[1.0]).unwrap();
    let mut new = old.clone();
    let mut outputs = [0.0; 2];
    with_cx(|cx| {
        let refused = circuit.propose(&old, 0.1, &[330.0, 310.0], &mut new, &mut outputs, cx);
        assert!(
            matches!(refused, Err(CircuitComponentError::Law { branch: 2, .. })),
            "{refused:?}"
        );
        let stepped = circuit.propose(&old, 0.1, &[310.0, 330.0], &mut new, &mut outputs, cx);
        assert!(matches!(stepped, Ok(Proposal::Step(_))), "{stepped:?}");
        let cold = circuit.propose(&old, 0.1, &[0.0, 330.0], &mut new, &mut outputs, cx);
        assert!(
            matches!(cold, Err(CircuitComponentError::Invalid { .. })),
            "{cold:?}"
        );
    });
}

#[test]
fn a_netlist_whose_ports_do_not_each_drive_one_source_is_refused() {
    let twice = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (2, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 1.0 }),
        ],
        transformers: vec![],
    };
    // The owner itself admits a port shared by two sources; this boundary
    // cannot give one port two roles or two outputs.
    assert!(matches!(
        KirchhoffCircuit::new(twice, vec![], 1.0e-12),
        Err(CircuitComponentError::Invalid { .. })
    ));
    // The owner counts ports up to the highest one named, so port 0 here
    // is a port that drives nothing.
    let gap = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 1 }),
            (1, 0, Branch::Resistor { ohms: 1.0 }),
        ],
        transformers: vec![],
    };
    assert!(matches!(
        KirchhoffCircuit::new(gap, vec![], 1.0e-12),
        Err(CircuitComponentError::Invalid { .. })
    ));
    let lawless = KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12).unwrap();
    let law = Arc::new(Linear {
        r0: 1.0,
        alpha: 0.0,
        floor_k: 1.0,
    });
    assert!(matches!(
        lawless.following(0, law.clone()),
        Err(CircuitComponentError::Invalid { .. })
    ));
    let once = KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12)
        .unwrap()
        .following(1, law.clone())
        .unwrap();
    assert!(matches!(
        once.following(1, law),
        Err(CircuitComponentError::Invalid { .. })
    ));
    assert!(matches!(
        KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 0.0),
        Err(CircuitComponentError::Invalid { .. })
    ));
}

/// One proposal of a circuit with no inductor and no capacitor, every heat
/// port at 300 K. Returns the outputs and the receipt's dissipation.
fn resistive_step(graph: CircuitGraph, sources: &[f64]) -> (Vec<f64>, f64) {
    let circuit = KirchhoffCircuit::new(graph, vec![], 1.0e-12).unwrap();
    let old = circuit.state(&[], &[]).unwrap();
    assert!(old.is_empty());
    let heats = circuit.ports().len() - sources.len();
    let mut inputs = sources.to_vec();
    inputs.extend(std::iter::repeat_n(300.0, heats));
    let mut new = vec![1.0];
    let mut outputs = vec![0.0; inputs.len()];
    let receipt =
        with_cx(
            |cx| match circuit.propose(&old, 0.5, &inputs, &mut new, &mut outputs, cx) {
                Ok(Proposal::Step(receipt)) => receipt,
                other => panic!("{other:?}"),
            },
        );
    assert!(new.is_empty(), "{new:?}");
    (outputs, receipt.dissipated_j)
}

#[test]
fn a_current_source_returns_the_voltage_it_drives_its_load_at() {
    // 2 A into a 5 ohm load: the port returns 10 V and the load takes 20 W.
    let (amperes, ohms) = (2.0_f64, 5.0_f64);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::CurrentSource { port: 0 }),
            (1, 0, Branch::Resistor { ohms }),
        ],
        transformers: vec![],
    };
    let (outputs, dissipated) = resistive_step(graph, &[amperes]);
    let power = amperes * amperes * ohms;
    close(outputs[0], amperes * ohms, 1.0e-12, "port voltage");
    close(outputs[1], -power, 1.0e-12, "heat rate");
    close(dissipated, 0.5 * power, 1.0e-12, "heat");
}

#[test]
fn an_ideal_transformer_passes_the_source_power_to_the_reflected_load() {
    // fs-phs's own reflected-load case: 12 V, 100 ohm on the primary, a 3:1
    // ideal transformer, 50 ohm on the secondary. The primary current is
    // U / (R1 + N^2 R2); the secondary carries N times it.
    let (volts, n, r1, r2) = (12.0_f64, 3.0_f64, 100.0_f64, 50.0_f64);
    let graph = CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: r1 }),
            (3, 0, Branch::Resistor { ohms: r2 }),
        ],
        transformers: vec![TransformerLink {
            primary: (2, 0),
            secondary: (3, 0),
            ratio: n,
        }],
    };
    let (outputs, dissipated) = resistive_step(graph, &[volts]);
    let current = volts / (r1 + n * n * r2);
    close(outputs[0], current, 1.0e-14, "primary current");
    close(
        outputs[1],
        -current * current * r1,
        1.0e-12,
        "primary heat rate",
    );
    close(
        outputs[2],
        -(n * current) * (n * current) * r2,
        1.0e-12,
        "secondary heat rate",
    );
    close(dissipated, 0.5 * volts * current, 1.0e-12, "heat");
}

#[test]
fn the_state_codec_round_trips_and_the_identity_names_the_law() {
    let law = |alpha| {
        Arc::new(Linear {
            r0: 1.0,
            alpha,
            floor_k: 1.0,
        })
    };
    let one = KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12)
        .unwrap()
        .following(1, law(0.01))
        .unwrap();
    let two = KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12)
        .unwrap()
        .following(1, law(0.02))
        .unwrap();
    let state = one.state(&[], &[-0.75]).unwrap();
    let mut words = Vec::new();
    assert!(one.encode_state(&state, &mut words));
    assert_eq!(one.decode_state(&words), Some(state));
    assert_eq!(one.decode_state(&[]), None);
    assert_eq!(one.decode_state(&[f64::NAN.to_bits()]), None);
    let elsewhere = KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12)
        .unwrap()
        .following(2, law(0.01))
        .unwrap();
    let (mut a, mut b, mut c) = (Vec::new(), Vec::new(), Vec::new());
    one.identity_words(&mut a);
    two.identity_words(&mut b);
    elsewhere.identity_words(&mut c);
    assert_ne!(a, b);
    assert_ne!(a, c);

    // Laws whose words run together the same way still differ: [2, 9] on
    // resistor 1 against [] on resistor 1 and [9] on resistor 2.
    let identity_of = |laws: &[(usize, &[u64])]| {
        let mut circuit = KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12).unwrap();
        for &(branch, law_words) in laws {
            circuit = circuit
                .following(branch, Arc::new(Words(law_words.to_vec())))
                .unwrap();
        }
        let mut identity = Vec::new();
        circuit.identity_words(&mut identity);
        identity
    };
    assert_ne!(
        identity_of(&[(1, &[2, 9])]),
        identity_of(&[(1, &[]), (2, &[9])])
    );
}

/// A constant 1 ohm whose identity is the given words.
#[derive(Debug)]
struct Words(Vec<u64>);

impl ResistanceLaw for Words {
    fn ohms(&self, _: f64) -> Result<f64, String> {
        Ok(1.0)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend(&self.0);
    }
}

#[test]
fn the_receipt_books_each_resistors_heat_at_the_temperature_its_port_receives() {
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 3.0 }),
            (2, 0, Branch::Capacitor { farads: 0.2 }),
            (2, 0, Branch::Resistor { ohms: 7.0 }),
        ],
        transformers: vec![],
    };
    let circuit = KirchhoffCircuit::new(graph, vec![], 1.0e-12).unwrap();
    let old = circuit.state(&[], &[0.3]).unwrap();
    let mut new = old.clone();
    let mut outputs = [0.0; 3];
    let (dt, kelvin) = (0.01_f64, [350.0_f64, 310.0_f64]);
    let receipt = with_cx(|cx| {
        match circuit.propose(
            &old,
            dt,
            &[9.0, kelvin[0], kelvin[1]],
            &mut new,
            &mut outputs,
            cx,
        ) {
            Ok(Proposal::Step(receipt)) => receipt,
            other => panic!("{other:?}"),
        }
    });
    let heat = [-dt * outputs[1], -dt * outputs[2]];
    assert!(heat[0] > 0.0 && heat[1] > 0.0, "{heat:?}");
    let scale = heat[0] + heat[1];
    close(
        receipt.dissipated_j,
        heat[0] + heat[1],
        4.0 * f64::EPSILON * scale,
        "heat",
    );
    close(
        receipt.entropy_generated_j_per_k,
        heat[0] / kelvin[0] + heat[1] / kelvin[1],
        4.0 * f64::EPSILON * scale / kelvin[1],
        "entropy",
    );
}

/// `r0 + slope |x|` at any temperature, `x` the voltage across the resistor
/// or, if `current`, the current through it: a conductor of unit length and
/// section whose resistivity rises with the field or the current density.
#[derive(Debug, Clone, Copy)]
struct Operating {
    r0: f64,
    slope: f64,
    current: bool,
}

impl ResistanceLaw for Operating {
    fn ohms(&self, _: f64) -> Result<f64, String> {
        Ok(self.r0)
    }

    fn follows_operating_point(&self) -> bool {
        true
    }

    fn ohms_at(&self, _: f64, volts: f64, amperes: f64) -> Result<f64, String> {
        let x = if self.current { amperes } else { volts };
        Ok(self.r0 + self.slope * x.abs())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            self.r0.to_bits(),
            self.slope.to_bits(),
            u64::from(self.current),
        ]);
    }
}

/// The controls of the module's conductivity fixtures: 64 trials, 1e-12 ohm
/// plus 1e-12 of the resistance, 1e-8 W.
fn iteration(relaxation: f64) -> OperatingPointIteration {
    OperatingPointIteration {
        maximum_trials: 64,
        absolute_resistance_tolerance_ohm: 1e-12,
        relative_resistance_tolerance: 1e-12,
        power_tolerance_w: 1e-8,
        relaxation,
    }
}

/// One 0.2 s proposal of a 1 F capacitor at 10 C discharged through resistor
/// 1 following `law` and resistor 2 at 3 ohm, both ports at 300 K. Returns
/// the charge at the end and the two resistors' heat rates \[W\].
fn operating_step(
    law: Operating,
    iteration: OperatingPointIteration,
) -> Result<(f64, [f64; 2]), CircuitComponentError> {
    let circuit = KirchhoffCircuit::new(series_rc(1.0, 999.0, 3.0), vec![], 1.0e-12)?
        .iterating(iteration)?
        .following(1, Arc::new(law))?;
    let old = circuit.state(&[], &[10.0])?;
    let mut new = old.clone();
    let mut outputs = [0.0; 2];
    let proposal =
        with_cx(|cx| circuit.propose(&old, 0.2, &[300.0, 300.0], &mut new, &mut outputs, cx))?;
    assert!(matches!(proposal, Proposal::Step(_)), "{proposal:?}");
    Ok((circuit.charges_c(&new)[0], outputs.map(|rate| -rate)))
}

#[test]
fn a_resistor_following_its_field_meets_the_hand_solved_midpoint_step() {
    // R1 = 1 + 0.05 |v| with v = I R1 is R1 = 1 / (1 - 0.05 I). The midpoint
    // step holds the capacitor's mean voltage 10 - 0.1 I across R1 + 3 ohm;
    // the current that does so is found by bisection.
    let beta = 0.05;
    let residual =
        |current: f64| current * (1.0 / (1.0 - beta * current) + 3.0) + 0.1 * current - 10.0;
    let (mut low, mut high) = (0.0_f64, 3.0_f64);
    for _ in 0..200 {
        let middle = 0.5 * (low + high);
        if residual(middle) > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    let current = 0.5 * (low + high);
    let r1 = 1.0 / (1.0 - beta * current);
    let field = Operating {
        r0: 1.0,
        slope: beta,
        current: false,
    };
    // Each run stops once R1 is within 2e-12 of itself. On the export of
    // 1e611fdf33 with series 1-8 the charge is within 2.1e-13 C of the
    // bisection's and the heat rates within 1.5e-11 W (of 6.3 and 16.7 W),
    // the relaxed run the farther; the bounds sit about fifty times above.
    for relaxation in [1.0, 0.5] {
        let (charge, heat) = operating_step(field, iteration(relaxation)).unwrap();
        close(charge, 10.0 - 0.2 * current, 1e-11, "charge");
        close(heat[0], current * current * r1, 1e-9, "heat 1");
        close(heat[1], current * current * 3.0, 1e-9, "heat 2");
    }
    // A resistance tolerance of half the resistance leaves the power one in
    // charge: the first trial's R1 is within 0.13 ohm of the law's, but its
    // power v i is 0.73 W off the law's R1 i^2. The ninth trial meets the
    // power tolerance, 4.8e-11 C from the bisection's charge on the same
    // export; the bound sits about fifty times above that.
    let loose = OperatingPointIteration {
        absolute_resistance_tolerance_ohm: 0.0,
        relative_resistance_tolerance: 0.5,
        ..iteration(1.0)
    };
    let (charge, _) = operating_step(field, loose).unwrap();
    close(
        charge,
        10.0 - 0.2 * current,
        2.5e-9,
        "charge, power in charge",
    );
    // One trial is the circuit at zero field, where R1 is 1 ohm.
    let once = OperatingPointIteration {
        maximum_trials: 1,
        ..iteration(1.0)
    };
    let refused = operating_step(field, once);
    assert!(
        matches!(
            refused,
            Err(CircuitComponentError::OperatingPoint { trials: 1, .. })
        ),
        "{refused:?}"
    );
}

#[test]
fn a_steep_current_law_meets_its_tolerances_only_relaxed() {
    // R1 = 0.1 + 1000 |i|: the midpoint step 10 = I (R1 + 3.1) is
    // 1000 I^2 + 3.2 I = 10. Unrelaxed, each trial near the answer misses it
    // by about 0.95 of the last miss, on the other side, and 64 trials do not
    // bring the resistance within 1e-12 of itself; trials halfway to the
    // solved point do.
    let steep = Operating {
        r0: 0.1,
        slope: 1000.0,
        current: true,
    };
    let current = 20.0 / (3.2 + (3.2_f64.powi(2) + 4.0 * 1000.0 * 10.0).sqrt());
    let refused = operating_step(steep, iteration(1.0));
    assert!(
        matches!(
            refused,
            Err(CircuitComponentError::OperatingPoint { trials: 64, .. })
        ),
        "{refused:?}"
    );
    // Relaxed, the trials close in fast: on the same export the charge is
    // the closed form's to the bit and the heat rates within 2e-14 W.
    let (charge, heat) = operating_step(steep, iteration(0.5)).unwrap();
    close(charge, 10.0 - 0.2 * current, 1e-12, "charge");
    close(
        heat[0],
        current * current * (0.1 + 1000.0 * current),
        1e-12,
        "heat 1",
    );
    close(heat[1], 3.0 * current * current, 1e-12, "heat 2");
}

#[test]
fn the_iteration_comes_once_before_its_laws_with_its_controls_in_range() {
    let field = Arc::new(Operating {
        r0: 1.0,
        slope: 0.05,
        current: false,
    });
    let base = || KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12).unwrap();
    // A law that follows the operating point before the iteration, and the
    // iteration twice.
    assert!(matches!(
        base().following(1, field),
        Err(CircuitComponentError::Invalid { .. })
    ));
    assert!(matches!(
        base()
            .iterating(iteration(1.0))
            .unwrap()
            .iterating(iteration(1.0)),
        Err(CircuitComponentError::Invalid { .. })
    ));
    for bad in [
        OperatingPointIteration {
            maximum_trials: 0,
            ..iteration(1.0)
        },
        OperatingPointIteration {
            absolute_resistance_tolerance_ohm: -1e-12,
            ..iteration(1.0)
        },
        OperatingPointIteration {
            relative_resistance_tolerance: 1.0,
            ..iteration(1.0)
        },
        OperatingPointIteration {
            absolute_resistance_tolerance_ohm: 0.0,
            relative_resistance_tolerance: 0.0,
            ..iteration(1.0)
        },
        OperatingPointIteration {
            power_tolerance_w: 0.0,
            ..iteration(1.0)
        },
        OperatingPointIteration {
            relaxation: 0.0,
            ..iteration(1.0)
        },
        OperatingPointIteration {
            relaxation: 1.5,
            ..iteration(1.0)
        },
        OperatingPointIteration {
            relaxation: f64::NAN,
            ..iteration(1.0)
        },
    ] {
        assert!(
            matches!(
                base().iterating(bad),
                Err(CircuitComponentError::Invalid { .. })
            ),
            "{bad:?}"
        );
    }
}

#[test]
fn the_iteration_joins_the_identity_and_passes_the_conformance_kit() {
    let field = Arc::new(Operating {
        r0: 1.0,
        slope: 0.05,
        current: false,
    });
    let base = || KirchhoffCircuit::new(series_rc(1.0, 1.0, 3.0), vec![], 1.0e-12).unwrap();
    // A circuit whose law follows temperature alone takes one trial: the
    // iteration changes its identity, not its step.
    let linear = Arc::new(Linear {
        r0: 1.0,
        alpha: 0.01,
        floor_k: 1.0,
    });
    let plain = base().following(1, linear.clone()).unwrap();
    let iterated = base()
        .iterating(iteration(0.5))
        .unwrap()
        .following(1, linear)
        .unwrap();
    let old = plain.state(&[], &[1.0]).unwrap();
    let step = |circuit: &KirchhoffCircuit| {
        let mut new = old.clone();
        let mut outputs = [0.0; 2];
        let proposal =
            with_cx(|cx| circuit.propose(&old, 0.1, &[320.0, 310.0], &mut new, &mut outputs, cx))
                .unwrap();
        let Proposal::Step(receipt) = proposal else {
            panic!("{proposal:?}");
        };
        new.iter()
            .chain(&outputs)
            .chain([&receipt.dissipated_j, &receipt.numerical_loss_j])
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    };
    assert_eq!(step(&plain), step(&iterated));
    let words = |circuit: &KirchhoffCircuit| {
        let mut words = Vec::new();
        circuit.identity_words(&mut words);
        words
    };
    assert_ne!(words(&plain), words(&iterated));
    let relaxed = base()
        .iterating(iteration(0.5))
        .unwrap()
        .following(1, field.clone())
        .unwrap();
    let unrelaxed = base()
        .iterating(iteration(1.0))
        .unwrap()
        .following(1, field)
        .unwrap();
    assert_ne!(words(&relaxed), words(&unrelaxed));

    // The iterated circuit passes the conformance kit.
    let state = relaxed.state(&[], &[10.0]).unwrap();
    let stale = relaxed.state(&[], &[-5.0]).unwrap();
    with_cx(|cx| {
        let case = Case {
            state: &state,
            stale: &stale,
            dt_s: 0.2,
            inputs: &[300.0, 300.0],
        };
        let observed = check_component(&relaxed, &case, cx).expect("component conformance");
        assert!(
            matches!(observed, Observed::Stepped { .. }),
            "the circuit shortened the step"
        );
    });
}
