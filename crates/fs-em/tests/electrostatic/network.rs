//! G1/G3: source-derived multi-electrode storage reaches the native circuit.
#[path = "network/floating.rs"]
mod floating;

use super::{controls, dielectric, limits, near, slab};
use fs_em::electrostatic::{
    ElectrostaticError, PlanarElectrostatics, network::CapacitanceMatrixControls,
};
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_rep_mesh::{Metric2, TriComplex2};

const DEPTH: f64 = 0.2;
const TERMINALS: [&[usize]; 3] = [&[0, 1], &[2, 3], &[4, 5]];

fn extraction() -> CapacitanceMatrixControls {
    CapacitanceMatrixControls {
        field: controls(),
        maximum_terminals: 4,
        maximum_field_solves: 4,
        maximum_field_values: 512,
        coefficient_tolerance_f: 1e-15,
    }
}
fn hexagon(center: [f64; 2], epsilon: [f64; 6]) -> PlanarElectrostatics {
    let height = 3.0_f64.sqrt() / 2.0;
    let mesh = TriComplex2::from_indexed_triangles(
        "three-electrode-dielectric",
        vec![
            [1.0, 0.0],
            [0.5, height],
            [-0.5, height],
            [-1.0, 0.0],
            [-0.5, -height],
            [0.5, -height],
            center,
        ],
        vec![
            [0, 1, 6],
            [1, 2, 6],
            [2, 3, 6],
            [3, 4, 6],
            [4, 5, 6],
            [5, 0, 6],
        ],
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap();
    PlanarElectrostatics::new(mesh, epsilon.map(dielectric).into(), limits(), &mut || {
        false
    })
    .unwrap()
}

#[test]
fn three_electrode_matrix_matches_independent_geometry_energy_charge_and_gauge() {
    let epsilon = 4e-6;
    let problem = hexagon([0.0; 2], [epsilon; 6]);
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    assert!(std::ptr::eq(matrix.problem(), &problem));
    assert_eq!(matrix.terminals(), TERMINALS.map(<[usize]>::to_vec));
    assert_eq!(matrix.unit_fields().len(), 3);
    // Independent equilateral P1 edge weights and elimination of the one
    // center potential: each pair is 7*epsilon*depth/(6*sqrt(3)).
    let pair = 7.0 * epsilon * DEPTH / (6.0 * 3.0_f64.sqrt());
    for i in 0..3 {
        for j in 0..3 {
            near(
                matrix.maxwell_matrix_f()[3 * i + j],
                if i == j { 2.0 * pair } else { -pair },
            );
        }
    }
    let grounded = matrix.grounded_matrix_f(1).unwrap();
    for (&actual, expected) in grounded.iter().zip([2.0 * pair, -pair, -pair, 2.0 * pair]) {
        near(actual, expected);
    }
    assert!(matrix.grounded_matrix_f(3).is_err());
    let volts = [7.0, 10.0, -1.0];
    let response = matrix.response(&volts, &mut || false).unwrap();
    let boundary: Vec<_> = TERMINALS
        .iter()
        .zip(volts)
        .flat_map(|(vertices, value)| vertices.iter().map(move |&vertex| (vertex, value)))
        .collect();
    let field = problem
        .solve(&[0.0; 6], &boundary, controls(), &mut || false)
        .unwrap();
    near(response.energy_j, field.energy_j());
    for (terminal, vertices) in TERMINALS.iter().enumerate() {
        let charge: f64 = vertices
            .iter()
            .map(|&v| field.nodal_boundary_charge_c()[v].unwrap())
            .sum();
        near(response.terminal_charge_c[terminal], charge);
        near(
            charge,
            pair * (2.0 * volts[terminal] - volts[(terminal + 1) % 3] - volts[(terminal + 2) % 3]),
        );
    }
    let shifted = matrix
        .response(&volts.map(|v| v + 1e9), &mut || false)
        .unwrap();
    assert_eq!(response.terminal_charge_c, shifted.terminal_charge_c);
    assert_eq!(response.energy_j.to_bits(), shifted.energy_j.to_bits());
    let common = matrix.response(&[1e9; 3], &mut || false).unwrap();
    assert_eq!(common.terminal_charge_c, [0.0; 3]);
    assert_eq!(common.energy_j, 0.0);
    assert!(matrix.maximum_coefficient_residual_f() < 1e-18);
    assert!(matrix.maximum_gauge_residual_f() < 1e-18);

    let permutation = [2, 0, 1];
    let permuted = problem
        .capacitance_matrix(
            &permutation.map(|i| TERMINALS[i]),
            extraction(),
            &mut || false,
        )
        .unwrap();
    for i in 0..3 {
        for j in 0..3 {
            near(
                permuted.maxwell_matrix_f()[3 * i + j],
                matrix.maxwell_matrix_f()[3 * permutation[i] + permutation[j]],
            );
        }
    }

    // The original scalar reduction uses the same terminal admission.
    let slab = slab(dielectric(2e-11), dielectric(4e-11));
    let scalar = slab
        .capacitance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    let multiple = slab
        .capacitance_matrix(&[&[2, 5], &[0, 3]], extraction(), &mut || false)
        .unwrap();
    for (&actual, expected) in multiple
        .maxwell_matrix_f()
        .iter()
        .zip([1.0, -1.0, -1.0, 1.0])
    {
        near(actual, expected * scalar.capacitance_f());
    }
}

#[test]
fn regional_dielectric_swap_changes_actual_multi_capacitor_circuit_dynamics() {
    let mut outcomes = Vec::new();
    for scale in [1.0, 3.0] {
        let mut epsilon = [4e-6; 6];
        epsilon[0] *= scale;
        epsilon[1] *= scale;
        let problem = hexagon([0.0; 2], epsilon);
        let matrix = problem
            .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
            .unwrap();
        let initial_voltage = [0.0, 1.0, -0.5];
        let reduced = matrix
            .circuit_at_voltages(&[0, 1, 2], &initial_voltage, 1e-15, 1e-15, &mut || false)
            .unwrap();
        assert!(std::ptr::eq(reduced.source(), &matrix));
        assert_eq!(reduced.branches().len(), 3);
        let mut branches = reduced.branches().to_vec();
        branches.push((1, 0, Branch::Resistor { ohms: 1000.0 }));
        branches.push((2, 0, Branch::Resistor { ohms: 2200.0 }));
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
        near(
            dae.system.hamiltonian(&initial),
            reduced.field_response().energy_j,
        );
        let dt = 1e-4;
        let (step, defect) = dae.step_audited(&initial, &[], dt).unwrap();

        // Independent two-node midpoint RC equations, not the descriptor solve.
        let c = matrix.grounded_matrix_f(0).unwrap();
        let a = [
            c[0] / dt + 0.5 / 1000.0,
            c[1] / dt,
            c[2] / dt,
            c[3] / dt + 0.5 / 2200.0,
        ];
        let b = [
            (c[0] / dt - 0.5 / 1000.0) * initial_voltage[1] + c[1] / dt * initial_voltage[2],
            c[2] / dt * initial_voltage[1] + (c[3] / dt - 0.5 / 2200.0) * initial_voltage[2],
        ];
        let determinant = a[0] * a[3] - a[1] * a[2];
        let voltage = [
            0.0,
            (b[0] * a[3] - a[1] * b[1]) / determinant,
            (a[0] * b[1] - b[0] * a[2]) / determinant,
        ];
        for (index, &(i, j, branch)) in reduced.branches().iter().enumerate() {
            let Branch::Capacitor { farads } = branch else {
                panic!("only capacitors lower from the field");
            };
            near(
                step.x[dae.charge_index[index]],
                farads * (voltage[i] - voltage[j]),
            );
            near(
                reduced.initial_charges_c()[index],
                farads * (initial_voltage[i] - initial_voltage[j]),
            );
        }
        let physical = matrix.response(&voltage, &mut || false).unwrap();
        near(dae.system.hamiltonian(&step.x), physical.energy_j);
        near(
            physical.energy_j + step.dissipated,
            reduced.field_response().energy_j,
        );
        let midpoint = [(voltage[1] + 1.0) / 2.0, (voltage[2] - 0.5) / 2.0];
        near(
            step.dissipated,
            dt * (midpoint[0] * midpoint[0] / 1000.0 + midpoint[1] * midpoint[1] / 2200.0),
        );
        assert!(defect < 1e-15 && step.dissipated > 0.0);
        assert!(reduced.maximum_charge_residual_c() < 1e-15);
        assert!(reduced.energy_residual_j().abs() < 1e-15);
        for (cell, &value) in epsilon.iter().enumerate() {
            near(
                matrix.unit_fields()[0].problem().permittivity_f_m().unwrap()[cell],
                value,
            );
            assert_eq!(
                matrix.unit_fields()[0].problem().materials()[cell],
                problem.materials()[cell]
            );
        }
        outcomes.push((voltage[1], physical.energy_j));
    }
    assert!((outcomes[0].0 - outcomes[1].0).abs() > 1e-3);
    assert!((outcomes[0].1 - outcomes[1].1).abs() > 1e-8);
}

#[test]
fn nonpassive_mutual_coefficient_is_retained_for_fields_and_refuses_circuit_lowering() {
    // A deliberately poor off-center high-contrast mesh has a positive
    // discrete mutual coefficient. Its energy Gram matrix is still PSD;
    // silently dropping that positive entry would change the field problem.
    let problem = hexagon(
        [0.6, 0.0],
        [0.5, 8.0, 0.1, 0.2, 80.0, 1.0].map(|x| x * 1e-6),
    );
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    assert!(matrix.maxwell_matrix_f()[5] > 1e-8);
    let response = matrix.response(&[0.0, 1.0, -1.0], &mut || false).unwrap();
    assert!(response.energy_j > 0.0);
    let result =
        matrix.circuit_at_voltages(&[0, 1, 2], &[0.0, 1.0, -1.0], 1e-15, 1e-15, &mut || false);
    assert!(matches!(
        result,
        Err(ElectrostaticError::Invalid(
            "positive mutual capacitance coefficient cannot form a passive circuit; refine the field discretization"
        ))
    ));
    assert!(matrix.maxwell_matrix_f()[5] > 1e-8);
}

#[test]
fn incomplete_traces_work_limits_node_aliases_and_cancellation_refuse_without_publication() {
    let problem = hexagon([0.0; 2], [4e-6; 6]);
    for policy in [
        CapacitanceMatrixControls {
            maximum_terminals: 2,
            ..extraction()
        },
        CapacitanceMatrixControls {
            maximum_field_solves: 2,
            ..extraction()
        },
        CapacitanceMatrixControls {
            maximum_field_values: 131,
            ..extraction()
        },
    ] {
        assert!(
            problem
                .capacitance_matrix(&TERMINALS, policy, &mut || false)
                .is_err()
        );
    }
    assert!(
        problem
            .capacitance_matrix(&[&[0, 1], &[1, 2], &[4, 5]], extraction(), &mut || false)
            .is_err()
    );
    assert!(
        problem
            .capacitance_matrix(&[&[0], &[2, 3], &[4, 5]], extraction(), &mut || false)
            .is_err()
    );
    assert!(
        problem
            .capacitance_matrix(&[&[0, 6], &[2, 3], &[4, 5]], extraction(), &mut || false)
            .is_err()
    );
    let mut visits = 0;
    let success = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || {
            visits += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    let cancelled = problem.capacitance_matrix(&TERMINALS, extraction(), &mut || {
        count += 1;
        count == visits
    });
    assert!(matches!(cancelled, Err(ElectrostaticError::Cancelled)));
    assert!(
        success
            .circuit_at_voltages(&[0, 1, 1], &[0.0, 1.0, 0.0], 1e-15, 1e-15, &mut || false)
            .is_err()
    );
    assert!(success.response(&[0.0, 1.0], &mut || false).is_err());
    assert!(
        success
            .response(&[0.0, f64::NAN, 1.0], &mut || false)
            .is_err()
    );
    assert!(matches!(
        success.response(&[0.0; 3], &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));
    let replay = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    assert_eq!(success.maxwell_matrix_f(), replay.maxwell_matrix_f());
}
