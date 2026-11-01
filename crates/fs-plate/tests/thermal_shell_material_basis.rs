//! G0/G1: changing source stiffness in the original retained physical basis.
//! Authored coefficients and independent work/strain oracles, not measurements.

use fs_material::sheet::SheetMatrix;
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::strain::integration_rule;
use fs_plate::shell::mitc::thermal::transient::EnrichedMitcShell;
use fs_plate::shell::mitc::thermal::{FacetThermalStrain, ShellThermalStrain};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_sparse::Csr;
use std::sync::Arc;

const AREA: f64 = 0.5;
const I0: f64 = 2.0;
const I1: f64 = 0.1;
const I2: f64 = 0.04;

fn section(changed: bool, mass: f64) -> ShearShellSection {
    let (a, b, d, shear) = if changed {
        (
            SheetMatrix::orthotropic(12000.0, 18000.0, 0.17, 4500.0).unwrap(),
            [2.0, 0.2, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0, 0.1],
            SheetMatrix::orthotropic(45.0, 7.0, 0.12, 4.5).unwrap(),
            [43.0, -4.0, 150.0],
        )
    } else {
        (
            SheetMatrix::orthotropic(20000.0, 10000.0, 0.2, 3000.0).unwrap(),
            [0.0; 9],
            SheetMatrix::orthotropic(20.0, 10.0, 0.2, 3.0).unwrap(),
            [120.0, 7.0, 80.0],
        )
    };
    ShearShellSection::new_coupled(a, b, d, shear, 1.0, mass, I1, I2, 0.0).unwrap()
}
fn model(changed: bool, mass: f64, length: f64, constrained: bool) -> Arc<MitcShellModel> {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [length, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    let constraints = if constrained {
        vec![ShellConstraint {
            node: 0,
            mask: [true, false, false, false, false, false],
            frame: ConstraintFrame::Global,
        }]
    } else {
        vec![]
    };
    Arc::new(
        assemble_mitc3p(
            &mesh,
            &directors,
            &[section(changed, mass)],
            &constraints,
            ShellMass::Consistent,
            MitcBudget {
                max_nodes: 3,
                max_triangles: 1,
            },
            if changed { 0.03 } else { 0.001 },
            || false,
        )
        .unwrap(),
    )
}
fn zero() -> [FacetThermalStrain; 1] {
    [ShellThermalStrain::try_new([0.0; 3], [0.0; 3])
        .unwrap()
        .into()]
}
fn hot() -> [FacetThermalStrain; 1] {
    [FacetThermalStrain::AtNodes([
        ShellThermalStrain::try_new([2e-4, -1e-4, 0.0], [0.0; 3]).unwrap(),
        ShellThermalStrain::try_new([2e-4, -1e-4, 0.0], [0.002, 0.0, 0.0]).unwrap(),
        ShellThermalStrain::try_new([2e-4, -1e-4, 0.0], [0.0, -0.001, 0.0]).unwrap(),
    ])]
}
fn nodal() -> [f64; 18] {
    let mut q = [0.0; 18];
    // Diagonal in-plane gradient has zero membrane spin, so theta_z=0
    // satisfies every original numerical drilling constraint.
    q[6] = 4e-4;
    q[13] = -2e-4;
    q[8] = 0.003;
    q[14] = -0.002;
    q[3] = 0.002;
    q[10] = -0.001;
    q[15] = -0.003;
    q
}
fn entry(matrix: &Csr, row: usize, column: usize) -> f64 {
    let (cols, values) = matrix.row(row);
    cols.iter()
        .zip(values)
        .find_map(|(&c, &v)| (c == column).then_some(v))
        .unwrap_or(0.0)
}
fn action(matrix: &Csr, q: &[f64]) -> Vec<f64> {
    (0..matrix.nrows())
        .map(|row| {
            let (cols, values) = matrix.row(row);
            cols.iter().zip(values).map(|(&c, &v)| v * q[c]).sum()
        })
        .collect()
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn near(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_nonproportional_material_change_preserves_strain_and_reference_momentum() {
    let reference_model = model(false, I0, 1.0, false);
    let next_model = model(true, I0, 1.0, false);
    let reference = EnrichedMitcShell::try_new(Arc::clone(&reference_model), || false).unwrap();
    let next = reference
        .with_material_model(Arc::clone(&next_model), || false)
        .unwrap();
    assert!(Arc::ptr_eq(&next.shared_model(), &next_model));
    assert!(std::ptr::eq(reference.mass(), next.mass()));
    let q = nodal();
    let v = q.map(|x| 2.3 * x);
    let initial = reference
        .state(0.2, &q, &v, &[[0.004, -0.002]], &[[0.3, -0.2]], || false)
        .unwrap();
    let same = next
        .state(
            0.2,
            &q,
            &v,
            initial.bubble_coordinates(),
            initial.bubble_velocity(),
            || false,
        )
        .unwrap();
    assert_eq!(initial.kinetic_energy_j(), same.kinetic_energy_j());
    assert_eq!(initial.coordinates(), same.coordinates());
    assert_eq!(initial.velocity(), same.velocity());

    let cold0 = reference_model.prepare_thermal(&zero(), || false).unwrap();
    let cold1 = next_model.prepare_thermal(&zero(), || false).unwrap();
    let before = reference.recover(&cold0, &initial, || false).unwrap();
    let after = next.recover(&cold1, &initial, || false).unwrap();
    for (lambda, _) in integration_rule() {
        let a = before.strain_fields[0].at_barycentric(lambda).unwrap();
        let b = after.strain_fields[0].at_barycentric(lambda).unwrap();
        assert_eq!(a, b);
    }
    assert!(before.bubble_material_force != after.bubble_material_force);
    let n = reference_model.pencil().free;
    assert!((0..n).any(|i| {
        (entry(next.stiffness(), i, n) - entry(reference.stiffness(), i, n)).abs() > 1e-4
    }));
    assert_eq!(before.numerical_gauge_force, after.numerical_gauge_force);
    assert_eq!(&after.numerical_gauge_force[n..], &[0.0, 0.0]);
    assert_eq!(
        before.drilling_gauge_energy_j,
        after.drilling_gauge_energy_j
    );
    let f0 = action(reference.stiffness(), initial.coordinates());
    let f1 = action(next.stiffness(), initial.coordinates());
    for i in 0..next.dimension() {
        near(
            after.material_generalized_force[i] + after.numerical_gauge_force[i],
            f1[i],
            2e-11,
            2e-12,
        );
    }
    let difference: Vec<_> = f1.iter().zip(f0).map(|(a, b)| a - b).collect();
    near(
        after.response.strain_energy_j - before.response.strain_energy_j,
        0.5 * dot(initial.coordinates(), &difference),
        2e-11,
        1e-13,
    );
    let back = next.with_material_model(reference_model, || false).unwrap();
    assert!(std::ptr::eq(reference.stiffness(), back.stiffness()));
    assert!(std::ptr::eq(reference.mass(), back.mass()));
}

#[test]
fn g1_original_basis_thermal_load_and_rebound_equilibrium_share_physical_work() {
    let reference_model = model(false, I0, 1.0, false);
    let next_model = model(true, I0, 1.0, false);
    let reference = EnrichedMitcShell::try_new(reference_model, || false).unwrap();
    let next = reference
        .with_material_model(Arc::clone(&next_model), || false)
        .unwrap();
    let loading = next_model.prepare_thermal(&hot(), || false).unwrap();
    let q = nodal();
    let state = next
        .state(0.0, &q, &[0.0; 18], &[[0.004, -0.002]], &[[0.0; 2]], || {
            false
        })
        .unwrap();
    let recovered = next.recover(&loading, &state, || false).unwrap();
    let load = next
        .generalized_load(&loading, &[0.0; 18], || false)
        .unwrap();
    let material = section(true, I0);
    let (a, d) = material.principal_matrices();
    let a = a.coefficients();
    let d = d.coefficients();
    let b = material.coupling_pa_m2();
    let mut work = 0.0;
    for (lambda, weight) in integration_rule() {
        let strain = recovered.strain_fields[0].at_barycentric(lambda).unwrap();
        let free = [2e-4, -1e-4, 0.0];
        let gradient = [0.002 * lambda[1], -0.001 * lambda[2], 0.0];
        for i in 0..3 {
            for j in 0..3 {
                work += AREA
                    * weight
                    * (strain.membrane[i] * (a[3 * i + j] * free[j] + b[3 * i + j] * gradient[j])
                        + strain.gradient_per_m[i]
                            * (b[3 * i + j] * free[j] + d[3 * i + j] * gradient[j]));
            }
        }
    }
    near(dot(state.coordinates(), &load), work, 2e-11, 1e-14);
    let seed = next.equilibrium_state(0.0, &loading, &q, || false).unwrap();
    let seeded = next.recover(&loading, &seed, || false).unwrap();
    for f in seeded.bubble_material_force[0] {
        near(f, 0.0, 0.0, 2e-11);
    }
    let fresh = EnrichedMitcShell::try_new(Arc::clone(&next_model), || false).unwrap();
    let current_seed = fresh
        .equilibrium_state(0.0, &loading, &q, || false)
        .unwrap();
    let current = fresh.recover(&loading, &current_seed, || false).unwrap();
    assert!(
        seed.bubble_coordinates()
            .iter()
            .flatten()
            .zip(current_seed.bubble_coordinates().iter().flatten())
            .any(|(a, b)| (a - b).abs() > 1e-6)
    );
    near(
        seeded.response.strain_energy_j,
        current.response.strain_energy_j,
        2e-11,
        1e-13,
    );
    for (lambda, _) in integration_rule() {
        let a = seeded.strain_fields[0].at_barycentric(lambda).unwrap();
        let b = current.strain_fields[0].at_barycentric(lambda).unwrap();
        for (a, b) in a
            .membrane
            .iter()
            .chain(&a.gradient_per_m)
            .chain(&a.transverse_shear)
            .zip(
                b.membrane
                    .iter()
                    .chain(&b.gradient_per_m)
                    .chain(&b.transverse_shear),
            )
        {
            near(*a, *b, 2e-11, 1e-14);
        }
    }
}

#[test]
fn g0_material_refresh_refuses_changed_reference_data_and_preserves_accepted_owners() {
    let original_model = model(false, I0, 1.0, false);
    let original = EnrichedMitcShell::try_new(Arc::clone(&original_model), || false).unwrap();
    for different in [
        model(true, 2.1, 1.0, false),
        model(true, I0, 1.01, false),
        model(true, I0, 1.0, true),
    ] {
        assert!(original.with_material_model(different, || false).is_err());
    }
    let next_model = model(true, I0, 1.0, false);
    assert!(matches!(
        original.with_material_model(Arc::clone(&next_model), || true),
        Err(PlateError::Cancelled)
    ));
    let mut checkpoints = 0;
    original
        .with_material_model(Arc::clone(&next_model), || {
            checkpoints += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        original.with_material_model(Arc::clone(&next_model), || {
            calls += 1;
            calls == checkpoints
        }),
        Err(PlateError::Cancelled)
    ));
    assert!(Arc::ptr_eq(&original.shared_model(), &original_model));
    let state = original
        .state(0.0, &nodal(), &[0.0; 18], &[[0.0; 2]], &[[0.0; 2]], || {
            false
        })
        .unwrap();
    let next = original
        .with_material_model(Arc::clone(&next_model), || false)
        .unwrap();
    let stale = original_model.prepare_thermal(&zero(), || false).unwrap();
    assert!(next.recover(&stale, &state, || false).is_err());
    assert!(next.generalized_load(&stale, &[0.0; 18], || false).is_err());
    let outsider = EnrichedMitcShell::try_new(model(false, I0, 1.0, false), || false).unwrap();
    let foreign = outsider
        .state(0.0, &nodal(), &[0.0; 18], &[[0.0; 2]], &[[0.0; 2]], || {
            false
        })
        .unwrap();
    let current = next_model.prepare_thermal(&zero(), || false).unwrap();
    assert!(next.recover(&current, &foreign, || false).is_err());
}
