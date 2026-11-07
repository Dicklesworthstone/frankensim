//! G0/G1/G3/G4: physical revolved conductors, trapped charge and native RC work.
use super::*;
use fs_em::electrostatic::axisymmetric::network::{
    CapacitanceMatrixControls, FloatingCapacitanceControls,
};
use fs_phs::circuit::{CircuitGraph, capacitance::assemble_circuit_with_capacitance_blocks};

fn extraction() -> CapacitanceMatrixControls {
    CapacitanceMatrixControls {
        field: controls(),
        maximum_terminals: 8,
        maximum_field_solves: 8,
        maximum_field_values: 1_000_000,
        coefficient_tolerance_f: 1e-15,
    }
}

fn reduction() -> FloatingCapacitanceControls {
    FloatingCapacitanceControls {
        maximum_floating_terminals: 4,
        maximum_dense_work: 4096,
        coefficient_tolerance_f: 1e-15,
    }
}

fn join(a: TriComplex2, b: TriComplex2, axial_offset: f64, angle: f64) -> TriComplex2 {
    let offset = a.vertices().len() as u32;
    let mut vertices = a.vertices().to_vec();
    vertices.extend(b.vertices().iter().map(|&[r, z]| [r, z + axial_offset]));
    let mut faces = a.faces().to_vec();
    faces.extend(b.faces().iter().map(|face| face.map(|v| v + offset)));
    TriComplex2::from_indexed_triangles(
        "separate dielectric pieces joined by a physical conductor",
        vertices,
        faces,
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

const DISKS: [&[usize]; 3] = [&[0, 1, 2], &[6, 7, 8, 9, 10, 11], &[15, 16, 17]];

fn disks(angle: f64, scale: f64) -> (AxisymmetricElectrostatics, [f64; 2]) {
    // Two cylindrical dielectric slabs, separated by a finite conducting disk.
    // The middle electrode joins its two disjoint physical surface traces.
    let mesh = join(
        rectangle(0.0, 2.0, 0.4, 2, 2, angle),
        rectangle(0.0, 2.0, 0.7, 2, 2, angle),
        0.6,
        angle,
    );
    let mut materials = vec![dielectric(3.0 * scale); 8];
    materials.extend(vec![dielectric(8.0 * scale); 8]);
    let p = AxisymmetricElectrostatics::new(mesh, materials, limits(), &mut || false).unwrap();
    let area = angle * 2.0;
    (p, [area * 3.0 * scale / 0.4, area * 8.0 * scale / 0.7])
}

#[test]
fn three_disks_retain_series_capacitance_trapped_charge_and_angular_span() {
    for angle in [core::f64::consts::TAU, 0.7] {
        let (p, [c1, c2]) = disks(angle, 1e-11);
        let csum = c1 + c2;
        let series = c1 * c2 / csum;
        let matrix = p
            .capacitance_matrix(&DISKS, extraction(), &mut || false)
            .unwrap();
        assert!(std::ptr::eq(matrix.problem(), &p));
        assert_eq!(matrix.terminals(), &DISKS.map(|t| t.to_vec()));
        let expected = [c1, -c1, 0.0, -c1, csum, -c2, 0.0, -c2, c2];
        for (&actual, exact) in matrix.maxwell_matrix_f().iter().zip(expected) {
            near(actual, exact, csum);
        }
        let floating = matrix
            .prepare_floating(&[1], reduction(), &mut || false)
            .unwrap();
        assert!(std::ptr::eq(floating.source(), &matrix));
        assert_eq!(floating.retained_terminals(), &[0, 2]);
        assert_eq!(floating.floating_terminals(), &[1]);
        for (&actual, sign) in floating
            .reduced_matrix_f()
            .iter()
            .zip([1.0, -1.0, -1.0, 1.0])
        {
            near(actual, sign * series, series);
        }
        for charge in [0.0, 0.7 * csum, -0.3 * csum] {
            let state = floating
                .solve(&[0.0, 6.0], &[charge], 1e-18, 1e-18, &mut || false)
                .unwrap();
            let middle = (charge + 6.0 * c2) / csum;
            near(state.component_relative_voltage_v[1], middle, 6.0);
            near(state.terminal_charge_c[0], -c1 * middle, 6.0 * csum);
            near(state.terminal_charge_c[1], charge, csum);
            near(state.terminal_charge_c[2], c2 * (6.0 - middle), 6.0 * csum);
            near(
                state.field().energy_j(),
                18.0 * series + charge * charge / (2.0 * csum),
                36.0 * csum,
            );
            assert!(std::ptr::eq(state.field().problem(), &p));

            let volts = [0.0, middle, 6.0];
            let boundary: Vec<_> = DISKS
                .iter()
                .zip(volts)
                .flat_map(|(vertices, voltage)| vertices.iter().map(move |&v| (v, voltage)))
                .collect();
            let direct = p
                .solve(&[0.0; 16], &boundary, controls(), &mut || false)
                .unwrap();
            let recovered = matrix
                .field_at_voltages(&volts, 1e-18, 1e-18, &mut || false)
                .unwrap();
            let response = matrix.response(&volts, &mut || false).unwrap();
            for (actual, exact) in state
                .field()
                .electric_field_v_m()
                .iter()
                .zip(direct.electric_field_v_m())
            {
                near(actual[0], exact[0], 20.0);
                near(actual[1], exact[1], 20.0);
            }
            near(recovered.energy_j(), direct.energy_j(), 36.0 * csum);
            near(response.energy_j, direct.energy_j(), 36.0 * csum);
            for (&a, &b) in response
                .terminal_charge_c
                .iter()
                .zip(&state.terminal_charge_c)
            {
                near(a, b, 6.0 * csum);
            }
            let biased = floating
                .solve(&[1e12, 1e12 + 6.0], &[charge], 1e-18, 1e-18, &mut || false)
                .unwrap();
            assert_eq!(biased.component_reference_voltage_v, vec![1e12; 3]);
            assert_eq!(
                biased.component_relative_voltage_v,
                state.component_relative_voltage_v
            );
            assert_eq!(biased.terminal_charge_c, state.terminal_charge_c);
            assert_eq!(biased.field().energy_j(), state.field().energy_j());
        }
        let trapped = floating
            .solve(&[0.0; 2], &[csum], 1e-18, 1e-18, &mut || false)
            .unwrap();
        near(trapped.field().energy_j(), 0.5 * csum, csum);
        let permutation = [2, 0, 1];
        let reordered = p
            .capacitance_matrix(&permutation.map(|t| DISKS[t]), extraction(), &mut || false)
            .unwrap();
        for i in 0..3 {
            for j in 0..3 {
                near(
                    reordered.maxwell_matrix_f()[3 * i + j],
                    matrix.maxwell_matrix_f()[3 * permutation[i] + permutation[j]],
                    csum,
                );
            }
        }
    }
}

#[test]
fn coaxial_floating_shell_converges_to_independent_logarithmic_capacitors() {
    let angle = 1.3;
    let height = 1.4;
    let c1 = angle * height * 3e-11 / 2.0_f64.ln();
    let c2 = angle * height * 8e-11 / (4.0_f64 / 2.5).ln();
    let expected = c1 * c2 / (c1 + c2);
    let mut previous = f64::INFINITY;
    for n in [8, 16] {
        let first = rectangle(1.0, 2.0, height, n, n, angle);
        let second = rectangle(2.5, 4.0, height, n, n, angle);
        let offset = first.vertices().len();
        let mesh = join(first, second, 0.0, angle);
        let mut materials = vec![dielectric(3e-11); 2 * n * n];
        materials.extend(vec![dielectric(8e-11); 2 * n * n]);
        let p = AxisymmetricElectrostatics::new(mesh, materials, limits(), &mut || false).unwrap();
        let inner: Vec<_> = (0..=n).map(|j| j * (n + 1)).collect();
        let middle: Vec<_> = inner
            .iter()
            .map(|&v| v + n)
            .chain(inner.iter().map(|&v| v + offset))
            .collect();
        let outer: Vec<_> = inner.iter().map(|&v| v + offset + n).collect();
        let matrix = p
            .capacitance_matrix(&[&inner, &middle, &outer], extraction(), &mut || false)
            .unwrap();
        let floating = matrix
            .prepare_floating(&[1], reduction(), &mut || false)
            .unwrap();
        let error = (floating.reduced_matrix_f()[0] / expected - 1.0).abs();
        assert!(error < 0.4 * previous, "n={n}: {error} versus {previous}");
        previous = error;
        let state = floating
            .solve(&[1.0, 0.0], &[0.0], 1e-17, 1e-17, &mut || false)
            .unwrap();
        near(
            state.terminal_charge_c[0],
            floating.reduced_matrix_f()[0],
            expected,
        );
        near(state.terminal_charge_c[1], 0.0, expected);
        near(
            state.field().energy_j(),
            0.5 * floating.reduced_matrix_f()[0],
            expected,
        );
    }
    assert!(previous < 1e-3);
}

#[test]
fn independent_biased_components_keep_small_charge_induced_fields_and_require_references() {
    let (p, [c1, c2]) = disks(0.7, 1e-11);
    let terminals: [&[usize]; 4] = [&[0, 1, 2], &[6, 7, 8], &[9, 10, 11], &[15, 16, 17]];
    let matrix = p
        .capacitance_matrix(&terminals, extraction(), &mut || false)
        .unwrap();
    let floating = matrix
        .prepare_floating(&[1, 3], reduction(), &mut || false)
        .unwrap();
    let charge = [c1 * 1e-6, -c2 * 2e-6];
    let original = floating
        .solve(&[0.0; 2], &charge, 1e-26, 1e-32, &mut || false)
        .unwrap();
    let bias = 2.0_f64.powi(50);
    let biased = floating
        .solve(&[bias, -bias], &charge, 1e-26, 1e-32, &mut || false)
        .unwrap();
    near(biased.component_relative_voltage_v[1], 1e-6, 1e-6);
    near(biased.component_relative_voltage_v[3], -2e-6, 2e-6);
    assert_eq!(
        biased.component_reference_voltage_v,
        [bias, bias, -bias, -bias]
    );
    assert_eq!(
        biased.component_relative_voltage_v,
        original.component_relative_voltage_v
    );
    assert_eq!(biased.terminal_charge_c, original.terminal_charge_c);
    assert_eq!(
        biased.field().electric_field_v_m(),
        original.field().electric_field_v_m()
    );
    near(
        biased.field().energy_j(),
        0.5 * (c1 + 4.0 * c2) * 1e-12,
        (c1 + c2) * 1e-12,
    );
    assert!(matches!(
        matrix.prepare_floating(&[2, 3], reduction(), &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
}

#[test]
fn native_circuit_keeps_floating_disk_charge_and_closes_rc_heat_against_radial_storage() {
    let (p, [c1, c2]) = disks(core::f64::consts::TAU, 1e-6);
    let csum = c1 + c2;
    let series = c1 * c2 / csum;
    let matrix = p
        .capacitance_matrix(&DISKS, extraction(), &mut || false)
        .unwrap();
    let floating = matrix
        .prepare_floating(&[1], reduction(), &mut || false)
        .unwrap();
    let charge = 0.7 * csum;
    let original = floating
        .solve(&[0.0, 6.0], &[charge], 1e-15, 1e-15, &mut || false)
        .unwrap();
    let bridge = matrix
        .coupled_circuit_at_voltages(
            0,
            &[0, 1, 2],
            &original.component_relative_voltage_v,
            1e-15,
            1e-15,
            &mut || false,
        )
        .unwrap();
    assert!(std::ptr::eq(bridge.source(), &matrix));
    assert_eq!(bridge.reference_terminal(), 0);
    assert_eq!(bridge.coordinate_terminals(), &[1, 2]);
    assert_eq!(
        bridge.capacitance_block().matrix_f(),
        matrix.grounded_matrix_f(0).unwrap()
    );
    let resistance = 1000.0;
    let mut branches = vec![(2, 0, Branch::Resistor { ohms: resistance })];
    branches.extend_from_slice(bridge.branches());
    let block = bridge.block_at_branch_offset(1).unwrap();
    assert_eq!(block.branches(), &[1, 2]);
    let dae = assemble_circuit_with_capacitance_blocks(
        &CircuitGraph {
            node_count: 3,
            branches,
            transformers: vec![],
        },
        &[block],
        &[],
        &[],
    )
    .unwrap();
    let mut state = vec![0.0; dae.system.state_dim()];
    for (&slot, &q) in dae.charge_index.iter().zip(bridge.initial_charges_c()) {
        state[slot] = q;
    }
    state = dae.consistent_initial_state(&state, &[]).unwrap();
    let dt = 0.2 * resistance * series;
    let midpoint = dt / (2.0 * resistance * series);
    let mut voltage = 6.0;
    let mut heat = 0.0;
    for _ in 0..8 {
        let (next, defect) = dae.step_audited(&state, &[], dt).unwrap();
        voltage *= (1.0 - midpoint) / (1.0 + midpoint);
        heat += next.dissipated;
        let physical = floating
            .solve(&[0.0, voltage], &[charge], 1e-15, 1e-15, &mut || false)
            .unwrap();
        for (&terminal, &slot) in bridge.coordinate_terminals().iter().zip(&dae.charge_index) {
            near(
                next.x[slot],
                physical.terminal_charge_c[terminal],
                6.0 * csum,
            );
        }
        near(next.x[dae.charge_index[0]], charge, charge);
        near(
            dae.system.hamiltonian(&next.x),
            physical.field().energy_j(),
            original.field().energy_j(),
        );
        near(
            physical.field().energy_j() + heat,
            original.field().energy_j(),
            original.field().energy_j(),
        );
        assert!(next.dissipated > 0.0);
        assert!(defect < 1e-15);
        state = next.x;
    }
}

#[test]
fn axis_only_conductors_budgets_and_cancelled_preparation_refuse_without_publication() {
    let p = problem(rectangle(0.0, 1.0, 1.0, 1, 1, 1.0), 2e-11);
    // All four nodes have prescribed potentials. A combined boundary check can
    // incorrectly lend the axis electrode physical support from top/bottom
    // edges that actually connect DIFFERENT electrodes. Check each conductor.
    assert!(matches!(
        p.capacitance(&[0, 2], &[1, 3], controls(), &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    assert!(matches!(
        p.capacitance_matrix(&[&[0, 2], &[1, 3]], extraction(), &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    let (p, _) = disks(0.7, 1e-11);
    let mut limited = extraction();
    limited.maximum_field_solves = 2;
    assert!(matches!(
        p.capacitance_matrix(&DISKS, limited, &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    let mut visits = 0;
    let original = p
        .capacitance_matrix(&DISKS, extraction(), &mut || {
            visits += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        p.capacitance_matrix(&DISKS, extraction(), &mut || {
            count += 1;
            count == visits / 2
        }),
        Err(ElectrostaticError::Cancelled)
    ));
    let repeat = p
        .capacitance_matrix(&DISKS, extraction(), &mut || false)
        .unwrap();
    assert_eq!(original.maxwell_matrix_f(), repeat.maxwell_matrix_f());
    assert!(matches!(
        repeat.prepare_floating(&[1], reduction(), &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));
}
