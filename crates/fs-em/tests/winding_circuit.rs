//! Field/material data must reach actual circuit evolution and its energy audit.
mod common;
use common::{config, limits, material, near};
use fs_em::{
    planar::{MagneticError, PlanarMagnetostatics},
    transfer::{CellSample, PlanarFrame},
    winding::{PlanarWindingBasis, WindingConfig, WindingError},
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_rep_mesh::{Metric2, TriComplex2};

fn winding_config() -> WindingConfig {
    WindingConfig {
        max_windings: 4,
        max_basis_values: 10000,
        pairing_tolerance_h: 1e-10,
        current_balance_tolerance_turns: 1e-9,
        inverse_residual_tolerance: 1e-10,
    }
}
fn square(mu_h_m: f64, depth_m: f64) -> (PlanarMagnetostatics, Vec<Vec<f64>>, Vec<usize>) {
    let n = 4_usize;
    let mut vertices = Vec::new();
    let mut boundary = Vec::new();
    for y in 0..=n {
        for x in 0..=n {
            if x == 0 || x == n || y == 0 || y == n {
                boundary.push(vertices.len());
            }
            vertices.push([x as f64 / n as f64, y as f64 / n as f64]);
        }
    }
    let mut faces = Vec::new();
    for y in 0..n {
        for x in 0..n {
            let a = (y * (n + 1) + x) as u32;
            let b = a + 1;
            let d = a + (n + 1) as u32;
            let c = d + 1;
            faces.extend([[a, b, c], [a, c, d]]);
        }
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "winding-field-circuit",
        vertices,
        faces,
        Metric2::planar(depth_m).unwrap(),
    )
    .unwrap();
    let mut profiles = vec![Vec::new(), Vec::new()];
    for face in mesh.faces() {
        let centre: [f64; 2] = [0, 1].map(|axis| {
            face.iter()
                .map(|&v| mesh.vertices()[v as usize][axis] / 3.0)
                .sum()
        });
        // Each winding has 100 outgoing and 100 return turns. The second
        // profile overlaps the first and adds an independent y-directed pair.
        let x = if centre[0] < 0.5 { 200.0 } else { -200.0 };
        let y = if centre[1] < 0.5 { 200.0 } else { -200.0 };
        profiles[0].push(x);
        profiles[1].push(x + y);
    }
    let materials = vec![material(mu_h_m); mesh.faces().len()];
    let problem = PlanarMagnetostatics::new(mesh, materials, limits(), &mut || false).unwrap();
    (problem, profiles, boundary)
}
fn basis<'a>(
    problem: &'a PlanarMagnetostatics,
    profiles: Vec<Vec<f64>>,
    boundary: &[usize],
    ends: Vec<f64>,
) -> PlanarWindingBasis<'a> {
    PlanarWindingBasis::new(
        problem,
        profiles,
        boundary,
        ends,
        winding_config(),
        config(),
        &mut || false,
    )
    .unwrap()
}
fn frame() -> PlanarFrame {
    PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap()
}
fn samples() -> [CellSample; 1] {
    [CellSample {
        triangle: 0,
        position_m: [0.2, 0.1, 0.25],
    }]
}

#[test]
fn flux_energy_gauge_and_declared_ends_share_the_original_spatial_basis() {
    let (problem, profiles, boundary) = square(2e-4, 0.5);
    let winding = basis(&problem, profiles, &boundary, vec![0.002, 0.003]);
    // Independent P1 reference on this four-by-four triangular cross-section.
    for (&actual, expected) in
        winding
            .active_inductance_h()
            .iter()
            .zip([1.0 / 28.0, 1.0 / 28.0, 1.0 / 28.0, 1.0 / 14.0])
    {
        near(actual, expected);
    }
    assert!(winding.maximum_pairing_defect_h() <= 1e-10);
    assert!(winding.inverse_residual() <= 1e-10);
    for residuals in winding.net_turns_by_component() {
        assert_eq!(residuals.len(), 1);
        near(residuals[0], 0.0);
    }

    let currents = [0.7, -0.2];
    let flux = winding.flux_linkages_wb(&currents).unwrap();
    let reconstructed = winding.currents_a(&flux).unwrap();
    for i in 0..2 {
        near(reconstructed[i], currents[i]);
    }
    let potential: Vec<_> = (0..problem.mesh().vertices().len())
        .map(|v| {
            currents[0] * winding.unit_fields()[0].potential_wb_m()[v]
                + currents[1] * winding.unit_fields()[1].potential_wb_m()[v]
        })
        .collect();
    let paired = winding
        .pair_flux_linkages_wb(&potential, &mut || false)
        .unwrap();
    let shifted: Vec<_> = potential.iter().map(|&a| a + 17.0).collect();
    let gauge = winding
        .pair_flux_linkages_wb(&shifted, &mut || false)
        .unwrap();
    for i in 0..2 {
        near(
            flux[i],
            paired[i] + winding.end_inductance_h()[i] * currents[i],
        );
        near(gauge[i], paired[i]);
    }
    let mut active_energy = 0.0;
    for cell in 0..problem.mesh().faces().len() {
        let b: [f64; 2] = [0, 1].map(|axis| {
            currents[0] * winding.unit_fields()[0].magnetic_flux_density_t()[cell][axis]
                + currents[1] * winding.unit_fields()[1].magnetic_flux_density_t()[cell][axis]
        });
        active_energy +=
            0.5 * problem.mesh().face_measure(cell).unwrap() * (b[0] * b[0] + b[1] * b[1])
                / problem.permeability_h_m()[cell];
    }
    let end_energy = 0.5 * (0.002 * currents[0] * currents[0] + 0.003 * currents[1] * currents[1]);
    near(
        0.5 * (currents[0] * flux[0] + currents[1] * flux[1]),
        active_energy + end_energy,
    );
    let sampled = winding
        .sample_for_currents(&currents, &frame(), &samples(), 1e-12, 1, &mut || false)
        .unwrap();
    for axis in 0..2 {
        near(
            sampled[0][axis],
            currents[0] * winding.unit_fields()[0].magnetic_flux_density_t()[0][axis]
                + currents[1] * winding.unit_fields()[1].magnetic_flux_density_t()[0][axis],
        );
    }
    near(sampled[0][2], 0.0);
}

#[test]
fn signed_turns_and_real_extrusion_depth_control_self_and_mutual_inductance() {
    let (problem, profiles, boundary) = square(2e-4, 0.5);
    let original = basis(&problem, profiles.clone(), &boundary, vec![0.0; 2]);
    let mut reoriented = profiles;
    for value in &mut reoriented[1] {
        *value *= -3.0;
    }
    let changed = basis(&problem, reoriented, &boundary, vec![0.0; 2]);
    let l = original.inductance_h();
    near(changed.inductance_h()[0], l[0]);
    near(changed.inductance_h()[1], -3.0 * l[1]);
    near(changed.inductance_h()[2], -3.0 * l[2]);
    near(changed.inductance_h()[3], 9.0 * l[3]);
    let (long_problem, long_profiles, long_boundary) = square(2e-4, 1.0);
    let longer = basis(&long_problem, long_profiles, &long_boundary, vec![0.0; 2]);
    for (&a, &b) in longer.inductance_h().iter().zip(l) {
        near(a, 2.0 * b);
    }
}

fn graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 4,
        // Capacitor first; secondary inductor before primary. The adapter
        // must use the original flux/charge maps, never winding position.
        branches: vec![
            (3, 0, Branch::Capacitor { farads: 0.4 }),
            (1, 0, Branch::Resistor { ohms: 1.0 }),
            (2, 0, Branch::Inductor { henries: 7.0 }),
            (3, 0, Branch::Resistor { ohms: 2.0 }),
            (1, 0, Branch::Inductor { henries: 11.0 }),
            (2, 0, Branch::Resistor { ohms: 3.0 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
        ],
        transformers: vec![],
    }
}
// Independent two-by-two midpoint equation: L di/dt = diag(R)*(source-i).
fn expected_currents(l: &[f64], step: f64, source: f64) -> [f64; 2] {
    let a = l[0] + step / 2.0;
    let b = l[1];
    let c = l[2];
    let d = l[3] + 3.0 * step / 2.0;
    let rhs0 = l[0] - step / 2.0 + step * source;
    let rhs1 = l[2];
    let determinant = a * d - b * c;
    [
        (d * rhs0 - b * rhs1) / determinant,
        (-c * rhs0 + a * rhs1) / determinant,
    ]
}

#[test]
fn material_data_changes_actual_mutual_circuit_dynamics_and_source_heat_energy_balance() {
    let mut outcomes = Vec::new();
    for mu in [2e-4, 4e-4] {
        let (problem, profiles, boundary) = square(mu, 0.5);
        let winding = basis(&problem, profiles, &boundary, vec![0.002, 0.003]);
        let circuit = winding
            .into_circuit(graph(), vec![4, 2], &mut || false)
            .unwrap();
        assert_eq!(circuit.branch_indices(), [4, 2]);
        assert_eq!(circuit.flux_coordinates(), [1, 0]);
        assert_eq!(circuit.dae().charge_index, [2]);
        assert_eq!(circuit.dae().system.state_dim(), 7);
        assert_eq!(circuit.dae().system.differential_dim(), 3);
        let total_l = circuit.basis().inductance_h();
        // Only source permeability data varies. The graph's placeholder Ls
        // cannot influence field-derived energy, response, or reconstruction.
        near(total_l[0], (mu / 2e-4) / 28.0 + 0.002);
        near(total_l[3], (mu / 2e-4) / 14.0 + 0.003);

        let source = 1.4;
        let step = 0.01;
        let q0 = 0.2;
        let initial = circuit
            .consistent_initial_state(&[1.0, 0.0], &[q0], &[source], &mut || false)
            .unwrap();
        let h0 = circuit.dae().system.hamiltonian(&initial);
        near(h0, 0.5 * total_l[0] + 0.5 * q0 * q0 / 0.4);
        let (record, defect) = circuit
            .dae()
            .step_audited(&initial, &[source], step)
            .unwrap();
        let current = circuit.currents_a(&record.x).unwrap();
        let expected = expected_currents(total_l, step, source);
        for i in 0..2 {
            near(current[i], expected[i]);
        }
        assert!(
            current[0] > 1.0 && current[1] < 0.0,
            "mutual induction must oppose the rising primary"
        );
        let q1 = q0 * (1.0 - step / (2.0 * 2.0 * 0.4)) / (1.0 + step / (2.0 * 2.0 * 0.4));
        near(record.x[circuit.dae().charge_index[0]], q1);
        let i_mid = [(1.0 + current[0]) / 2.0, current[1] / 2.0];
        let v_mid = [source - i_mid[0], -3.0 * i_mid[1]];
        let cap_voltage_mid = (q0 + q1) / (2.0 * 0.4);
        let heat = step
            * (v_mid[0] * v_mid[0]
                + v_mid[1] * v_mid[1] / 3.0
                + cap_voltage_mid * cap_voltage_mid / 2.0);
        let work = step * source * v_mid[0];
        near(record.dissipated, heat);
        near(record.supplied, work);
        near(record.delta_h, work - heat);
        near(
            circuit.dae().system.hamiltonian(&record.x) - h0,
            record.delta_h,
        );
        assert!(record.supplied > 0.0 && record.dissipated > 0.0);
        assert!(defect < 1e-9, "coupled circuit ledger: {defect}");
        assert!(record.solver_residual < 1e-8);
        let sampled = circuit
            .sample_field(&record.x, &frame(), &samples(), 1e-12, 1, &mut || false)
            .unwrap();
        let b = sampled[0];
        for axis in 0..2 {
            near(
                b[axis],
                current[0] * circuit.basis().unit_fields()[0].magnetic_flux_density_t()[0][axis]
                    + current[1]
                        * circuit.basis().unit_fields()[1].magnetic_flux_density_t()[0][axis],
            );
        }
        outcomes.push((current, b));
    }
    assert!(outcomes[1].0[0] < outcomes[0].0[0]);
    assert!(outcomes[1].0[1] > outcomes[0].0[1]);
    assert!(
        (outcomes[1].1[0] - outcomes[0].1[0]).abs() + (outcomes[1].1[1] - outcomes[0].1[1]).abs()
            > 1e-5
    );
}

#[test]
fn globally_balanced_current_cannot_cross_independent_component_gauges() {
    let mesh = TriComplex2::from_indexed_triangles(
        "separate-winding-components",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [2.0, 1.0],
        ],
        vec![[0, 1, 2], [3, 4, 5]],
        Metric2::planar(0.5).unwrap(),
    )
    .unwrap();
    let problem =
        PlanarMagnetostatics::new(mesh, vec![material(2e-4); 2], limits(), &mut || false).unwrap();
    let result = PlanarWindingBasis::new(
        &problem,
        vec![vec![1.0, -1.0]],
        &[0, 1, 2, 3, 4, 5],
        vec![0.1],
        winding_config(),
        config(),
        &mut || false,
    );
    assert!(matches!(
        result,
        Err(WindingError::Invalid(
            "unit current must close separately on every connected mesh component"
        ))
    ));
}

#[test]
fn dependent_windings_need_declared_physical_leakage_and_bad_bindings_refuse() {
    let (problem, profiles, boundary) = square(2e-4, 0.5);
    let dependent = vec![profiles[0].clone(), profiles[0].clone()];
    assert!(matches!(
        PlanarWindingBasis::new(
            &problem,
            dependent.clone(),
            &boundary,
            vec![0.0; 2],
            winding_config(),
            config(),
            &mut || false
        ),
        Err(WindingError::SingularInductance)
    ));
    // Physical end/leakage data can remove the dependency; never invented by
    // numerical regularization. The caller explicitly owns this approximation.
    let admitted = basis(&problem, dependent, &boundary, vec![0.002, 0.003]);
    assert!(admitted.inductance_h()[0] > admitted.active_inductance_h()[0]);
    assert!(matches!(
        admitted.into_circuit(graph(), vec![4, 4], &mut || false),
        Err(WindingError::Invalid("duplicate winding branch binding"))
    ));
    assert!(matches!(
        PlanarWindingBasis::new(
            &problem,
            profiles.clone(),
            &boundary,
            vec![-0.1, 0.0],
            winding_config(),
            config(),
            &mut || false
        ),
        Err(WindingError::Invalid(_))
    ));
    let mut exhausted = winding_config();
    exhausted.max_basis_values = 1;
    assert!(matches!(
        PlanarWindingBasis::new(
            &problem,
            profiles.clone(),
            &boundary,
            vec![0.0; 2],
            exhausted,
            config(),
            &mut || false
        ),
        Err(WindingError::Invalid(_))
    ));
    assert!(matches!(
        PlanarWindingBasis::new(
            &problem,
            profiles.clone(),
            &boundary,
            vec![0.0; 2],
            winding_config(),
            config(),
            &mut || true
        ),
        Err(WindingError::Magnetic(MagneticError::Cancelled))
    ));
    let recovered = basis(&problem, profiles, &boundary, vec![0.002, 0.003]);
    assert!(
        recovered
            .into_circuit(graph(), vec![4, 2], &mut || false)
            .is_ok()
    );
}
