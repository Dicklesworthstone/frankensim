//! G0/G1/G3/G4 nonzero support motion through native P1 solve/stress owners.
#[path = "tet10/support.rs"]
mod support;
use support::{block, material, with_cx};
use fs_blake3::ContentHash;
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

fn thermal(mat: &TetElasticMaterial, strain: [f64; 6]) -> TetThermalStrainState {
    TetThermalStrainState::try_new(310.0, 300.0, strain,
        mat.material_state_identity, ContentHash([17; 32]), 0.1).unwrap()
}
fn problem<'a>(mesh: &'a support::Mesh, mat: &'a TetElasticMaterial, fixed: &'a [usize])
    -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(mat), fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard() }
}
fn grips(mesh: &support::Mesh, length: f64, shift: [f64; 3], strain: f64) -> Vec<(usize, f64)> {
    let mut values = Vec::new();
    for (node, p) in mesh.nodes.iter().enumerate() {
        if p[0] == 0.0 || p[0] == length { values.push((3 * node, shift[0] + strain * p[0])); }
        if p[1] == 0.0 { values.push((3 * node + 1, shift[1])); }
        if p[2] == 0.0 { values.push((3 * node + 2, shift[2])); }
    }
    values.sort_by_key(|(dof, _)| *dof);
    values
}

#[test]
fn g1_driven_extension_combines_thermal_strain_poisson_motion_and_support_reactions() {
    with_cx(|cx| {
        for n in [1, 2, 3] {
            let mesh = block([n, 2, 2], [2.0, 0.5, 0.25]);
            let mat = material(2000.0, 0.25);
            let alpha = 0.0003;
            let field = thermal(&mat, [alpha, alpha, alpha, 0.0, 0.0, 0.0]);
            let shift = [0.001, -0.002, 0.003];
            let extension = 0.002;
            let prescribed = grips(&mesh, 2.0, shift, extension);
            let fixed: Vec<_> = prescribed.iter().map(|(d, _)| *d).collect();
            let p = problem(&mesh, &mat, &fixed);
            let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget,
                4 * mesh.tets.len(), cx).unwrap();
            let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let result = load.solve_with_prescribed_dofs(&p, TetThermalStrainField::Uniform(&field),
                &prescribed, 0.1, Default::default(), ContentHash([18; 32]), 0.1, cx).unwrap();
            let lateral = alpha - 0.25 * (extension - alpha);
            for (position, displacement) in mesh.nodes.iter().zip(result.displacement.displacement_m()) {
                let expected = [shift[0] + extension * position[0], shift[1] + lateral * position[1],
                    shift[2] + lateral * position[2]];
                for c in 0..3 { assert!((displacement[c] - expected[c]).abs() < 1e-10); }
            }
            let sigma = 2000.0 * (extension - alpha);
            for element in result.stress.elements() {
                for c in 0..6 {
                    let expected = if c == 0 { sigma } else { 0.0 };
                    assert!((element.stress_mandel_pa()[c] - expected).abs() < 1e-6);
                }
            }
            let right: f64 = mesh.nodes.iter().zip(&result.support_reactions_n)
                .filter(|(p, _)| p[0] == 2.0).map(|(_, r)| r[0]).sum();
            assert!((right - sigma * 0.5 * 0.25).abs() < 1e-7);
            assert!(result.free_force_residual_n < 1e-7);
            let expected_energy = 0.5 * 2000.0 * (extension - alpha).powi(2) * 2.0 * 0.5 * 0.25;
            assert!((result.stress.elastic_energy_j() - expected_energy).abs() < 1e-9);
        }
    });
}

#[test]
fn g1_heterogeneous_bar_preserves_series_compliance_and_nonuniform_free_strain() {
    with_cx(|cx| {
        let mesh = block([4, 1, 1], [1.0; 3]);
        let left: Vec<_> = mesh.tets.iter().map(|t|
            t.iter().map(|&i| mesh.nodes[i][0]).sum::<f64>() / 4.0 < 0.5).collect();
        let mats: Vec<_> = left.iter().map(|l| material(if *l { 1000.0 } else { 4000.0 }, 0.0)).collect();
        let states: Vec<_> = mats.iter().zip(&left).map(|(m, l)| thermal(m,
            [if *l { 0.001 } else { 0.002 }, 0.0, 0.0, 0.0, 0.0, 0.0])).collect();
        let prescribed = grips(&mesh, 1.0, [0.001, 0.0, 0.0], 0.006);
        let fixed: Vec<_> = prescribed.iter().map(|(d, _)| *d).collect();
        let p = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
            materials: TetMaterialField::PerElement(&mats), fixed_dofs: &fixed,
            budget: Default::default() };
        let field = TetThermalStrainField::PerElement(&states);
        let u = p.solve_thermal_displacement_with_prescribed_dofs(field, &[], &prescribed,
            0.1, Default::default(), cx).unwrap();
        let stress = p.recover_thermal_stress_with_prescribed_dofs(&u, field, &[], &prescribed,
            0.1, ContentHash([18; 32]), 0.1, cx).unwrap();
        let sigma = (0.006 - 0.5 * 0.001 - 0.5 * 0.002) / (0.5 / 1000.0 + 0.5 / 4000.0);
        for (point, actual) in mesh.nodes.iter().zip(u.displacement_m()) {
            let x = point[0];
            let expected = 0.001 + x.min(0.5) * (0.001 + sigma / 1000.0)
                + (x - 0.5).max(0.0) * (0.002 + sigma / 4000.0);
            assert!((actual[0] - expected).abs() < 1e-10);
        }
        for s in stress.elements() { assert!((s.stress_mandel_pa()[0] - sigma).abs() < 1e-6); }
    });
}

#[test]
fn g0_all_prescribed_affine_motion_has_exact_values_and_physical_fixed_node_loads() {
    with_cx(|cx| {
        let mesh = block([1; 3], [1.0; 3]);
        let mat = material(1200.0, 0.25);
        let field = thermal(&mat, [0.001, 0.001, 0.001, 0.0, 0.0, 0.0]);
        let prescribed: Vec<_> = mesh.nodes.iter().enumerate().flat_map(|(n, p)|
            [0.002 + 0.003 * p[0] + 0.001 * p[1], -0.004 + 0.002 * p[1], 0.001 - 0.001 * p[2]]
                .into_iter().enumerate().map(move |(c, v)| (3 * n + c, v))).collect();
        let fixed: Vec<_> = prescribed.iter().map(|(d, _)| *d).collect();
        let p = problem(&mesh, &mat, &fixed);
        let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 24, cx).unwrap();
        let no_force = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let forced = geometry.assemble(&[], VolumeLoad3::None, &[(0, 2.0)], cx).unwrap();
        let solve = |load: &fs_solid::loads3::TetLoadCase3<'_, '_>| load.solve_with_prescribed_dofs(
            &p, TetThermalStrainField::Uniform(&field), &prescribed, 0.1, Default::default(),
            ContentHash([18; 32]), 0.1, cx).unwrap();
        let a = solve(&no_force);
        let b = solve(&forced);
        assert_eq!(a.displacement.iterations(), 0);
        assert_eq!(a.displacement.true_relative_residual(), 0.0);
        assert_eq!(a.displacement.displacement_m(), b.displacement.displacement_m());
        for &(dof, value) in &prescribed {
            assert_eq!(a.displacement.displacement_m()[dof / 3][dof % 3].to_bits(), value.to_bits());
        }
        assert!((b.support_reactions_n[0][0] - a.support_reactions_n[0][0] + 2.0).abs() < 1e-12);
        assert_eq!(a.free_force_residual_n, 0.0);
        let strain = [0.002, 0.001, -0.002, 0.001 * std::f64::consts::FRAC_1_SQRT_2, 0.0, 0.0];
        let expected: [f64; 6] = std::array::from_fn(|i|
            (0..6).map(|j| mat.stiffness_mandel_pa()[i][j] * strain[j]).sum());
        for s in a.stress.elements() {
            for c in 0..6 { assert!((s.stress_mandel_pa()[c] - expected[c]).abs() < 1e-10); }
        }
        let updated = p.update_geometry_from_displacement(&a.displacement, 0.1, cx).unwrap();
        assert!(updated.nodes_m().iter().zip(&mesh.nodes).any(|(a, b)| a != b));
    });
}

#[test]
fn g3_zero_prescription_preserves_original_numerics_and_nonzero_translation_has_no_stress() {
    with_cx(|cx| {
        let mesh = block([2; 3], [1.0; 3]);
        let mat = material(1500.0, 0.2);
        let field = thermal(&mat, [0.0; 6]);
        let values = grips(&mesh, 1.0, [0.0; 3], 0.0);
        let fixed: Vec<_> = values.iter().map(|(d, _)| *d).collect();
        let p = problem(&mesh, &mat, &fixed);
        let field = TetThermalStrainField::Uniform(&field);
        let force = [(3 * (mesh.nodes.len() - 1) + 1, 0.1)];
        let old = p.solve_thermal_displacement_with_nodal_forces(field, &force, Default::default(), cx).unwrap();
        let new = p.solve_thermal_displacement_with_prescribed_dofs(field, &force, &[], 0.1,
            Default::default(), cx).unwrap();
        assert_eq!(old.displacement_m(), new.displacement_m());
        let shift = [0.001, -0.002, 0.003];
        let translated = grips(&mesh, 1.0, shift, 0.0);
        let u = p.solve_thermal_displacement_with_prescribed_dofs(field, &[], &translated,
            0.1, Default::default(), cx).unwrap();
        let s = p.recover_thermal_stress_with_prescribed_dofs(&u, field, &[], &translated,
            0.1, ContentHash([18; 32]), 0.1, cx).unwrap();
        for actual in u.displacement_m() {
            for c in 0..3 { assert!((actual[c] - shift[c]).abs() < 1e-10); }
        }
        assert!(s.elastic_energy_j() < 1e-14);
    });
}

#[test]
fn g0_stress_recovery_refuses_changed_prescriptions_or_legacy_zero_support_interpretation() {
    with_cx(|cx| {
        let mesh = block([1; 3], [1.0; 3]);
        let mat = material(1200.0, 0.25);
        let state = thermal(&mat, [0.0; 6]);
        let values = grips(&mesh, 1.0, [0.001; 3], 0.002);
        let fixed: Vec<_> = values.iter().map(|(d, _)| *d).collect();
        let p = problem(&mesh, &mat, &fixed);
        let field = TetThermalStrainField::Uniform(&state);
        let u = p.solve_thermal_displacement_with_prescribed_dofs(field, &[], &values,
            0.1, Default::default(), cx).unwrap();
        assert!(p.recover_thermal_stress_with_nodal_forces(&u, field, &[],
            ContentHash([18; 32]), 0.1, cx).is_err());
        let mut changed = values.clone(); changed[0].1 += 0.001;
        assert!(matches!(p.recover_thermal_stress_with_prescribed_dofs(&u, field, &[], &changed,
            0.1, ContentHash([18; 32]), 0.1, cx), Err(TetElasticError::StressRecoveryInputMismatch)));
        let again = p.solve_thermal_displacement_with_prescribed_dofs(field, &[], &values,
            0.1, Default::default(), cx).unwrap();
        assert_eq!(u, again);
    });
}

#[test]
fn g0_prescribed_support_admission_refuses_wrong_dofs_limits_and_unrestrained_modes() {
    with_cx(|cx| {
        let mesh = block([1; 3], [1.0; 3]);
        let mat = material(1200.0, 0.25);
        let state = thermal(&mat, [0.0; 6]);
        let values = grips(&mesh, 1.0, [0.0; 3], 0.0);
        let fixed: Vec<_> = values.iter().map(|(d, _)| *d).collect();
        let p = problem(&mesh, &mat, &fixed);
        let field = TetThermalStrainField::Uniform(&state);
        let free = (0..3 * mesh.nodes.len()).find(|d| !fixed.contains(d)).unwrap();
        for bad in [vec![(0, f64::NAN)], vec![(0, 0.2)], vec![(0, 0.0), (0, 0.0)],
            vec![(1, 0.0), (0, 0.0)], vec![(free, 0.001)], vec![(3 * mesh.nodes.len(), 0.0)]] {
            assert!(p.solve_thermal_displacement_with_prescribed_dofs(field, &[], &bad,
                0.1, Default::default(), cx).is_err());
        }
        for maximum in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(p.solve_thermal_displacement_with_prescribed_dofs(field, &[], &[],
                maximum, Default::default(), cx).is_err());
        }
        let under = problem(&mesh, &mat, &[0, 1, 2]);
        assert!(matches!(under.solve_thermal_displacement_with_prescribed_dofs(field, &[],
            &[(0, 0.001)], 0.1, Default::default(), cx), Err(TetElasticError::UnconstrainedRigidModes { .. })));
    });
}

#[test]
fn g4_prescribed_solve_observes_precancellation() {
    let mesh = block([1; 3], [1.0; 3]);
    let mat = material(1200.0, 0.25);
    let state = thermal(&mat, [0.0; 6]);
    let p = problem(&mesh, &mat, &[0, 1, 2]);
    let gate = fs_exec::CancelGate::new_clock_free(); gate.request();
    let pool = fs_alloc::ArenaPool::new(Default::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(&gate, arena, fs_exec::StreamKey {
            seed: 1, kernel_id: 1, tile: 0, iteration: 0,
        }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        assert!(matches!(p.solve_thermal_displacement_with_prescribed_dofs(
            TetThermalStrainField::Uniform(&state), &[], &[], 0.1, Default::default(), &cx),
            Err(TetElasticError::Cancelled)));
    });
}
