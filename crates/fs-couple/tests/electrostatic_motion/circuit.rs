//! G1/G3/G4: finite native supplies, charge/geometry work and atomic publication.

use super::*;
use fs_couple::electrostatic::motion::circuit::{CircuitMotionInput, PlanarCapacitorCircuitMotion};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve, MaterialResistor};
use fs_material::conductor::{ElectricalPropertyBasis, resolve_uniform_conductor_with_basis};
use fs_phs::circuit::{Branch, CircuitGraph};

fn resistor(branch: usize, resistance: f64) -> MaterialResistor {
    let basis = ElectricalPropertyBasis::Conductivity;
    let source = card(basis.property_name(), basis.dimensions(), 1.0 / resistance);
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor_with_basis(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            basis,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}
fn circuit_controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 1e-10,
        power_tolerance_w: 1e-10,
    }
}
fn input(voltage: &[f64]) -> CircuitMotionInput<'_> {
    CircuitMotionInput {
        electrical_ports: voltage,
        external_force_n: 0.01,
        duration_s: 0.2,
    }
}
fn isolated(charge: f64) -> MaterialCircuit {
    MaterialCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![(1, 0, Branch::Capacitor { farads: 2.0 })],
            transformers: vec![],
        },
        vec![],
        vec![],
        vec![charge],
        &[],
    )
    .unwrap()
}

#[test]
fn finite_circuit_reduces_to_the_original_isolated_field_body_motion() {
    for charge in [-0.4, 0.0, 0.4] {
        let state = PlanarCapacitorCircuitMotion::new(
            path(2.0),
            body(1.0),
            isolated(charge),
            0,
            0.0,
            0.02,
            &[],
            &mut || false,
        )
        .unwrap();
        let fixed = PlanarCapacitorMotion::new(
            path(2.0),
            body(1.0),
            ElectricalControl::IsolatedCharge { charge_c: charge },
            0.0,
            0.02,
            &mut || false,
        )
        .unwrap();
        let drive = input(&[]);
        let expected = fixed
            .step_candidate(
                MotionInput {
                    duration_s: drive.duration_s,
                    external_force_n: drive.external_force_n,
                },
                controls(),
                &mut || false,
            )
            .unwrap();
        let candidate = state
            .step_candidate(drive, controls(), circuit_controls(), &mut || false)
            .unwrap();
        near(
            candidate.state().coordinate_m(),
            expected.state().coordinate_m(),
            2e-12,
        );
        near(
            candidate.state().velocity_m_s(),
            expected.state().velocity_m_s(),
            2e-12,
        );
        near(candidate.state().charge_c(), charge, 1e-14);
        near(
            candidate.step().field_work_on_body_j,
            expected.step().field_work_on_body_j,
            1e-13,
        );
        near(candidate.step().total_energy_residual_j, 0.0, 1e-11);
        assert_eq!(candidate.step().heat_j, 0.0);
        assert_eq!(candidate.step().electrical_supply_j, 0.0);
    }
}

// The active capacitor is deliberately NOT the first electrical coordinate.
// The independent RC at the final node proves other storage and heat survive.
fn driven(inductance: Option<f64>, resistance: f64) -> (MaterialCircuit, usize) {
    let extra = if inductance.is_some() { 4 } else { 3 };
    let mut graph = CircuitGraph {
        node_count: extra + 1,
        branches: vec![
            (extra, 0, Branch::Capacitor { farads: 0.7 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 999.0 }),
        ],
        transformers: vec![],
    };
    if let Some(henries) = inductance {
        graph.branches.push((2, 3, Branch::Inductor { henries }));
    }
    let branch = graph.branches.len();
    graph.branches.push((
        if inductance.is_some() { 3 } else { 2 },
        0,
        Branch::Capacitor { farads: 2.0 },
    ));
    let extra_resistor = graph.branches.len();
    graph
        .branches
        .push((extra, 0, Branch::Resistor { ohms: 999.0 }));
    let circuit = MaterialCircuit::new(
        graph,
        vec![resistor(extra_resistor, 3.0), resistor(2, resistance)],
        inductance.map_or_else(Vec::new, |_| vec![0.08]),
        vec![-0.21, 0.4],
        &[0.9],
    )
    .unwrap();
    (circuit, branch)
}

#[test]
fn driven_rc_and_rlc_match_independent_plate_equations_and_actual_source_heat() {
    let (epsilon, mass, q0, phi0, v0, dt, external) = (2.0, 1.0, 0.4, 0.08, 0.02, 0.2, 0.01);
    for (inductance, resistance, source) in [
        (None, 1.7, 0.9),
        (None, 3.4, 0.9),
        (Some(0.5), 1.7, 0.9),
        (Some(0.5), 3.4, 0.9),
        // Reverse electrode charge within one accepted interval: force stays
        // attractive, whereas an erroneous Q0*Q1 rule would reverse it.
        (None, 0.1, -1.3),
    ] {
        let (circuit, branch) = driven(inductance, resistance);
        let original_materials = circuit.resistors().to_vec();
        let original_graph = circuit.graph().clone();
        let state = PlanarCapacitorCircuitMotion::new(
            path(epsilon),
            body(mass),
            circuit,
            branch,
            0.0,
            v0,
            &[source],
            &mut || false,
        )
        .unwrap();
        let expected_charge = |displacement: f64| {
            // Independent physical plate reduction C=epsilon*A/g, A=1.
            let inverse_c = (1.0 + 0.5 * displacement) / epsilon;
            match inductance {
                None => {
                    (q0 * (1.0 - dt * inverse_c / (2.0 * resistance)) + dt * source / resistance)
                        / (1.0 + dt * inverse_c / (2.0 * resistance))
                }
                Some(l) => {
                    q0 + (dt * source + 2.0 * phi0 - dt * q0 * inverse_c)
                        / (2.0 * l / dt + resistance + 0.5 * dt * inverse_c)
                }
            }
        };
        let force = |q1: f64| -(q0 * q0 + q1 * q1) / (4.0 * epsilon);
        // Solve the independently reduced scalar plate equation; this does
        // not invoke the production field, circuit or feedback solvers.
        let (mut lo, mut hi) = (-0.08, 0.08);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            let residual =
                mid - dt * v0 - 0.5 * dt * dt * (external + force(expected_charge(mid))) / mass;
            if residual < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let displacement = 0.5 * (lo + hi);
        let q1 = expected_charge(displacement);
        let velocity = v0 + dt * (external + force(q1)) / mass;
        let current = (q1 - q0) / dt;
        let extra_factor = (1.0 - dt / (2.0 * 3.0 * 0.7)) / (1.0 + dt / (2.0 * 3.0 * 0.7));
        let other_q = -0.21 * extra_factor;
        let other_current = (other_q + 0.21) / dt;
        let phi1 = inductance.map(|l| 2.0 * l * current - phi0);
        let candidate = state
            .step_candidate(
                input(&[source]),
                controls(),
                circuit_controls(),
                &mut || false,
            )
            .unwrap();
        let next = candidate.state();
        let step = candidate.step();
        near(next.coordinate_m(), displacement, 2e-11);
        near(next.velocity_m_s(), velocity, 2e-11);
        near(next.charge_c(), q1, 2e-11);
        near(next.circuit().charges_c()[0], other_q, 2e-11);
        if let Some(phi) = phi1 {
            near(next.circuit().fluxes_wb()[0], phi, 2e-11);
        }
        near(
            next.capacitance().capacitance_f(),
            epsilon / (1.0 + displacement),
            2e-12,
        );
        near(step.discrete_field_force_n, force(q1), 2e-10);
        near(step.electrical_supply_j, source * (q1 - q0), 2e-11);
        near(
            step.heat_j,
            dt * (resistance * current * current + 3.0 * other_current * other_current),
            2e-11,
        );
        near(step.resistors[0].current_a, current, 2e-11);
        let energy = 0.5 * q1 * q1 * (1.0 + displacement) / epsilon
            + 0.5 * other_q * other_q / 0.7
            + phi1
                .zip(inductance)
                .map_or(0.0, |(phi, l)| 0.5 * phi * phi / l);
        near(next.electrical_energy_j(), energy, 2e-11);
        near(step.field_work_on_body_j, force(q1) * displacement, 2e-12);
        near(
            step.electrical_energy_j[1] - step.electrical_energy_j[0]
                + 0.5 * mass * (velocity * velocity - v0 * v0)
                + step.heat_j
                - source * (q1 - q0)
                - external * displacement,
            0.0,
            2e-11,
        );
        near(step.total_energy_residual_j, 0.0, 1e-10);
        assert_eq!(next.circuit().resistors(), original_materials);
        for (index, (old, new)) in original_graph
            .branches
            .iter()
            .zip(&next.circuit().graph().branches)
            .enumerate()
        {
            if index != branch {
                assert_eq!(old, new);
            }
        }
        assert_eq!(state.coordinate_m(), 0.0);
        assert_eq!(state.charge_c(), q0);
        assert!(step.field_solves > 0);
        if source < 0.0 {
            assert!(next.charge_c() < 0.0 && step.discrete_field_force_n < 0.0);
        }
    }
}

#[test]
fn heterogeneous_tapered_field_retains_charge_feedback_over_the_shared_history() {
    let field_path = PlanarCapacitorPath::new(
        reference(2.0, 6.0, true),
        [vec![2, 5], vec![0, 3]],
        vec![
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
        ],
        [-0.4, 0.4],
        limits(),
        field_controls(),
        &mut || false,
    )
    .unwrap();
    let (circuit, branch) = driven(Some(0.5), 1.7);
    let mut state = PlanarCapacitorCircuitMotion::new(
        field_path,
        body(1.0),
        circuit,
        branch,
        0.0,
        -0.1,
        &[0.9],
        &mut || false,
    )
    .unwrap();
    let mut frozen = state.circuit().clone();
    let original_sources = state.path().reference().materials().to_vec();
    let original_energy = state.electrical_energy_j() + 0.5 * state.velocity_m_s().powi(2);
    let (mut supply, mut heat) = (0.0, 0.0);
    for _ in 0..8 {
        let drive = CircuitMotionInput {
            duration_s: 0.05,
            ..input(&[0.9])
        };
        let (next, step) = state
            .step_candidate(drive, controls(), circuit_controls(), &mut || false)
            .unwrap()
            .into_parts();
        let actual = next
            .path()
            .sample(next.coordinate_m(), &mut || false)
            .unwrap();
        near(
            next.capacitance().capacitance_f(),
            actual.capacitance_f(),
            1e-12,
        );
        near(
            next.voltage_v(),
            next.charge_c() / actual.capacitance_f(),
            1e-12,
        );
        assert_eq!(next.path().reference().materials(), original_sources);
        supply += step.electrical_supply_j + step.external_mechanical_work_j;
        heat += step.heat_j;
        near(
            next.electrical_energy_j() + 0.5 * next.velocity_m_s().powi(2) + heat
                - original_energy
                - supply,
            0.0,
            5e-10,
        );
        frozen = frozen
            .step_candidate(
                frozen.resistors().to_vec(),
                &[0.9],
                drive.duration_s,
                circuit_controls(),
                || false,
            )
            .unwrap()
            .into_parts()
            .0;
        state = next;
    }
    assert!((state.charge_c() - frozen.charges_c()[1]).abs() > 1e-5);
    assert!(state.coordinate_m() < -0.01);
}

#[test]
fn circuit_motion_refuses_competing_storage_clamps_and_late_cancel_atomically() {
    use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
    use fs_material::field_storage::ResolvedCapacitor;
    let source = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        2.0,
    );
    let material_c = ResolvedCapacitor::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
        1.0,
    )
    .unwrap();
    let bound = MaterialCircuit::new_with_storage(
        isolated(0.4).graph().clone(),
        vec![],
        vec![MaterialStorage {
            branch: 0,
            element: MaterialStorageElement::Capacitor(material_c),
        }],
        vec![],
        vec![0.4],
        &[],
    )
    .unwrap();
    assert!(
        PlanarCapacitorCircuitMotion::new(
            path(2.0),
            body(1.0),
            bound,
            0,
            0.0,
            0.0,
            &[],
            &mut || false
        )
        .is_err()
    );

    let clamped = MaterialCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 2.0 }),
                (1, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        vec![],
        vec![],
        vec![0.4],
        &[0.2],
    )
    .unwrap();
    let clamped = PlanarCapacitorCircuitMotion::new(
        path(2.0),
        body(1.0),
        clamped,
        0,
        0.0,
        0.0,
        &[0.2],
        &mut || false,
    )
    .unwrap();
    let saved = clamped.circuit().clone();
    assert!(
        clamped
            .step_candidate(input(&[0.2]), controls(), circuit_controls(), &mut || false)
            .is_err()
    );
    assert_eq!(clamped.circuit(), &saved);

    let (circuit, branch) = driven(Some(0.5), 1.7);
    let state = PlanarCapacitorCircuitMotion::new(
        path(2.0),
        body(1.0),
        circuit,
        branch,
        0.0,
        0.02,
        &[0.9],
        &mut || false,
    )
    .unwrap();
    let saved = state.circuit().clone();
    let mut calls = 0;
    let completed = state
        .step_candidate(input(&[0.9]), controls(), circuit_controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut actual = 0;
    assert!(matches!(
        state.step_candidate(input(&[0.9]), controls(), circuit_controls(), &mut || {
            actual += 1;
            actual == calls
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(state.circuit(), &saved);
    assert_eq!(state.coordinate_m(), 0.0);
    let retry = state
        .step_candidate(input(&[0.9]), controls(), circuit_controls(), &mut || false)
        .unwrap();
    assert_eq!(
        retry.state().coordinate_m(),
        completed.state().coordinate_m()
    );
    assert_eq!(retry.state().circuit(), completed.state().circuit());
    let narrow = MotionControls {
        displacement_bracket_m: [0.0, 0.0],
        ..controls()
    };
    assert!(
        state
            .step_candidate(input(&[0.9]), narrow, circuit_controls(), &mut || false)
            .is_err()
    );
    assert_eq!(state.circuit(), &saved);
}
