//! G1/G3/G4: the same material field drives native electrical and mechanical
//! owners. Synthetic sources test the connection, not empirical device accuracy.

use super::{config, dielectric, limits};
use fs_em::{
    electrostatic::{
        ElectrostaticControls, ElectrostaticError, PlanarElectrostatics,
        network::CapacitanceMatrixControls,
    },
    force::MaxwellForceLimits,
    transfer::PlanarFrame,
};
use fs_ga::{Mat3, Se3, Vec3, Wrench};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    free_floating_forward_dynamics,
};
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_rep_mesh::{Metric2, TriComplex2};

const DEPTH: f64 = 0.2;
const EPSILON: f64 = 2e-6;

fn controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: config(),
        charge_tolerance_c: 1e-11,
        energy_tolerance_j: 1e-10,
    }
}
fn force_limits() -> MaxwellForceLimits {
    MaxwellForceLimits {
        max_cells: 10000,
        max_edges: 20000,
    }
}
fn frame() -> PlanarFrame {
    PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-8 * b.abs().max(1e-8), "{a} != {b}");
}

fn specimen(
    n: usize,
    epsilon: f64,
    reverse: bool,
    depth: f64,
) -> (PlanarElectrostatics, Vec<usize>, Vec<usize>) {
    let points: Vec<_> = (0..=n)
        .flat_map(|y| {
            (0..=n).map(move |x| {
                [
                    -1.0 + 2.0 * x as f64 / n as f64,
                    -1.0 + 2.0 * y as f64 / n as f64,
                ]
            })
        })
        .collect();
    let mut faces = Vec::new();
    for y in 0..n {
        for x in 0..n {
            let a = (y * (n + 1) + x) as u32;
            let b = a + 1;
            let d = a + (n + 1) as u32;
            let c = d + 1;
            faces.extend(if reverse {
                [[a, c, b], [a, d, c]]
            } else {
                [[a, b, c], [a, c, d]]
            });
        }
    }
    let mut enclosed = Vec::new();
    let mut body = Vec::new();
    let ambient = dielectric(EPSILON);
    let inside = dielectric(epsilon);
    let mut sources = Vec::new();
    for (cell, face) in faces.iter().enumerate() {
        let center: [f64; 2] =
            [0, 1].map(|axis| face.iter().map(|&v| points[v as usize][axis]).sum::<f64>() / 3.0);
        let radius = center[0].abs().max(center[1].abs());
        if radius < 0.5 {
            enclosed.push(cell);
        }
        if radius < 0.25 {
            body.push(cell);
            sources.push(inside.clone());
        } else {
            sources.push(ambient.clone());
        }
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "electric-force",
        points,
        faces,
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    (
        PlanarElectrostatics::new(mesh, sources, limits(), &mut || false).unwrap(),
        enclosed,
        body,
    )
}
fn boundary(
    problem: &PlanarElectrostatics,
    potential: impl Fn([f64; 2]) -> f64,
) -> Vec<(usize, f64)> {
    problem
        .mesh()
        .boundary_trace()
        .unwrap()
        .vertices
        .into_iter()
        .map(|v| (v as usize, potential(problem.mesh().vertices()[v as usize])))
        .collect()
}

#[test]
fn exterior_stress_converges_to_charge_force_and_uniform_fields_have_no_net_load() {
    // The centered charge's self-force cancels by symmetry. Its imposed E_x=10
    // gives F_x = rho * (0.5^2 * depth) * 10, independently of stress integration.
    let exact = 1e-3 * 0.25 * DEPTH * 10.0;
    let mut previous = f64::INFINITY;
    for n in [8, 16, 32] {
        let (problem, enclosed, body) = specimen(n, EPSILON, false, DEPTH);
        let mut rho = vec![0.0; problem.mesh().faces().len()];
        for cell in body {
            rho[cell] = 1e-3;
        }
        let field = problem
            .solve(
                &rho,
                &boundary(&problem, |[x, _]| -10.0 * x),
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(field.free_charge_density_c_m3(), rho);
        let load = field
            .maxwell_wrench(
                &enclosed,
                &frame(),
                [0.0, 0.0, DEPTH / 2.0],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let error = (load.force.x - exact).abs();
        assert!(
            error < 0.6 * previous,
            "n={n}, error={error}, previous={previous}"
        );
        assert!(load.force.y.abs() < 1e-12 && load.torque.z.abs() < 1e-12);
        previous = error;
    }
    assert!(previous < 0.05 * exact);

    for reverse in [false, true] {
        let (problem, mut enclosed, _) = specimen(8, EPSILON, reverse, 3.0 * DEPTH);
        let field = problem
            .solve(
                &vec![0.0; problem.mesh().faces().len()],
                &boundary(&problem, |[x, y]| -10.0 * (2.0 * x + y)),
                controls(),
                &mut || false,
            )
            .unwrap();
        let first = field
            .maxwell_wrench(
                &enclosed,
                &frame(),
                [0.1, 0.2, 0.3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        assert!(first.to_array().iter().all(|v| v.abs() < 1e-12));
        enclosed.reverse();
        assert_eq!(
            first,
            field
                .maxwell_wrench(
                    &enclosed,
                    &frame(),
                    [0.1, 0.2, 0.3],
                    force_limits(),
                    &mut || false
                )
                .unwrap()
        );
    }
}

#[test]
fn dielectric_data_only_changes_force_and_native_body_acceleration() {
    let inertia = SpatialInertia::new(
        2.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [0.2, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.2],
        },
    )
    .unwrap();
    let mechanics = ArticulatedModel::new(vec![Link::new(
        "dielectric",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let mut acceleration = Vec::new();
    for epsilon in [EPSILON, 2.0 * EPSILON, 4.0 * EPSILON] {
        let (problem, enclosed, _) = specimen(16, epsilon, false, DEPTH);
        let field = problem
            .solve(
                &vec![0.0; problem.mesh().faces().len()],
                &boundary(&problem, |[x, y]| -100.0 * (x + 0.4 * (x * x - y * y))),
                controls(),
                &mut || false,
            )
            .unwrap();
        let placement =
            PlanarFrame::new([0.0, 0.0, -DEPTH / 2.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let load = field
            .maxwell_wrench(&enclosed, &placement, [0.0; 3], force_limits(), &mut || {
                false
            })
            .unwrap();
        let dynamics = free_floating_forward_dynamics(
            &mechanics,
            FreeFloatingBaseState::stationary(Se3::identity()),
            &[],
            &[],
            &[],
            Vec3::new(0.0, 0.0, 0.0),
            &[load],
        )
        .unwrap();
        close(
            dynamics.base_spatial_acceleration_body.linear.x,
            load.force.x / 2.0,
        );
        close(
            dynamics.base_spatial_acceleration_body.angular.z,
            load.torque.z / 0.2,
        );
        acceleration.push(dynamics.base_spatial_acceleration_body.linear.x);

        let moved = PlanarFrame::new([2.0, 3.0, 4.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
        let origin = [2.0 + DEPTH / 2.0 + 0.2, 3.3, 4.4];
        let global = field
            .maxwell_wrench(&enclosed, &moved, origin, force_limits(), &mut || false)
            .unwrap();
        let force = Vec3::new(load.force.z, load.force.x, load.force.y);
        let torque = Vec3::new(load.torque.z, load.torque.x, load.torque.y)
            + Vec3::new(-0.2, -0.3, -0.4).cross(force);
        for (actual, expected) in global
            .to_array()
            .into_iter()
            .zip(Wrench::new(torque, force).to_array())
        {
            close(actual, expected);
        }
    }
    assert!(acceleration[0].abs() < 1e-10);
    assert!(acceleration[1] > 1e-5 && acceleration[2] > 1.4 * acceleration[1]);
}

#[test]
fn native_multi_electrode_discharge_reconstructs_the_same_spatial_force() {
    let (problem, enclosed, _) = specimen(16, 4.0 * EPSILON, false, DEPTH);
    let vertices = problem.mesh().vertices();
    let left: Vec<_> = (0..vertices.len())
        .filter(|&i| vertices[i][0] == -1.0)
        .collect();
    let right: Vec<_> = (0..vertices.len())
        .filter(|&i| vertices[i][0] == 1.0)
        .collect();
    let top: Vec<_> = (0..vertices.len())
        .filter(|&i| vertices[i][1] == 1.0 && vertices[i][0].abs() <= 0.5)
        .collect();
    let terminals = [left.as_slice(), right.as_slice(), top.as_slice()];
    let matrix = problem
        .capacitance_matrix(
            &terminals,
            CapacitanceMatrixControls {
                field: controls(),
                maximum_terminals: 3,
                maximum_field_solves: 3,
                maximum_field_values: 20000,
                coefficient_tolerance_f: 1e-12,
            },
            &mut || false,
        )
        .unwrap();
    let initial_voltage = [0.0, 10.0, 30.0];
    let reduced = matrix
        .circuit_at_voltages(&[0, 1, 2], &initial_voltage, 1e-10, 1e-10, &mut || false)
        .unwrap();
    let mut branches = reduced.branches().to_vec();
    branches.extend([
        (1, 0, Branch::Resistor { ohms: 1000.0 }),
        (2, 0, Branch::Resistor { ohms: 2200.0 }),
    ]);
    let dae = assemble_circuit(&CircuitGraph {
        node_count: 3,
        branches,
        transformers: vec![],
    })
    .unwrap();
    let mut initial = vec![0.0; dae.system.state_dim()];
    for (&slot, &charge) in dae.charge_index.iter().zip(reduced.initial_charges_c()) {
        initial[slot] = charge;
    }
    let initial = dae.consistent_initial_state(&initial, &[]).unwrap();
    let (step, defect) = dae.step_audited(&initial, &[], 1e-4).unwrap();
    let initialized = dae.consistent_initial_state(&step.x, &[]).unwrap();
    // The native dt=0 descriptor initializer encodes algebraic efforts at the
    // midpoint of its input and output, while retaining physical charge.
    let voltage = [
        0.0,
        f64::midpoint(
            step.x[dae.node_potential_index[0]],
            initialized[dae.node_potential_index[0]],
        ),
        f64::midpoint(
            step.x[dae.node_potential_index[1]],
            initialized[dae.node_potential_index[1]],
        ),
    ];
    let field = matrix
        .field_at_voltages(&voltage, 1e-10, 1e-10, &mut || false)
        .unwrap();
    assert!(field.report().is_none()); // retained basis, no second Krylov solve
    close(field.energy_j(), dae.system.hamiltonian(&step.x));
    assert!(defect < 1e-10 && step.dissipated > 0.0);
    let bc: Vec<_> = terminals
        .iter()
        .zip(voltage)
        .flat_map(|(nodes, v)| nodes.iter().map(move |&node| (node, v)))
        .collect();
    let direct = problem
        .solve(
            &vec![0.0; problem.mesh().faces().len()],
            &bc,
            controls(),
            &mut || false,
        )
        .unwrap();
    let load = field
        .maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, DEPTH / 2.0],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    let original = direct
        .maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, DEPTH / 2.0],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    for (a, b) in load.to_array().into_iter().zip(original.to_array()) {
        close(a, b);
    }
    let shifted = matrix
        .field_at_voltages(&voltage.map(|v| v + 100.0), 1e-10, 1e-10, &mut || false)
        .unwrap();
    for (&a, &b) in field.voltage_v().iter().zip(shifted.voltage_v()) {
        close(a, b);
    }
    let before = matrix
        .field_at_voltages(&initial_voltage, 1e-10, 1e-10, &mut || false)
        .unwrap()
        .maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, DEPTH / 2.0],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert!((before.force.y - load.force.y).abs() > 1e-10);
    assert!(matches!(
        matrix.field_at_voltages(&voltage, 1e-10, 1e-10, &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));
}

#[test]
fn charge_material_budget_and_cancellation_refusals_leave_the_field_intact() {
    let (problem, enclosed, body) = specimen(8, 2.0 * EPSILON, false, DEPTH);
    let field = problem
        .solve(
            &vec![0.0; problem.mesh().faces().len()],
            &boundary(&problem, |[x, _]| -x),
            controls(),
            &mut || false,
        )
        .unwrap();
    let baseline = field.voltage_v().to_vec();
    for selected in [
        vec![],
        vec![enclosed[0]; 2],
        vec![usize::MAX],
        vec![0],
        body,
    ] {
        assert!(
            field
                .maxwell_wrench(&selected, &frame(), [0.0; 3], force_limits(), &mut || false)
                .is_err()
        );
    }
    assert!(
        field
            .maxwell_wrench(
                &enclosed,
                &frame(),
                [0.0; 3],
                MaxwellForceLimits {
                    max_cells: 1,
                    max_edges: 1
                },
                &mut || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        field.maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || {
            polls += 1;
            polls == 8
        }),
        Err(ElectrostaticError::Cancelled)
    ));
    assert_eq!(field.voltage_v(), baseline);
    let (homogeneous, enclosed, _) = specimen(8, EPSILON, false, DEPTH);
    // Opposite sources cancel globally but do not make a charge-free contour.
    let rho: Vec<_> = (0..homogeneous.mesh().faces().len())
        .map(|cell| if cell % 2 == 0 { 1e-6 } else { -1e-6 })
        .collect();
    let charged = homogeneous
        .solve(
            &rho,
            &boundary(&homogeneous, |[x, _]| -x),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(charged.volume_charge_c().abs() < 1e-20);
    assert!(matches!(
        charged.maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false),
        Err(ElectrostaticError::Invalid(
            "electric contour must avoid free-charge sources and electrode vertices"
        ))
    ));
}
