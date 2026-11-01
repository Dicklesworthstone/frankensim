//! G0/G1: endpoint material work and momentum in one retained physical basis.
//! Synthetic coefficients are numerical references, not material measurements.

use fs_material::sheet::SheetMatrix;
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::transient::EnrichedMitcShell;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_sparse::Csr;
use std::sync::Arc;

const AREA: f64 = 0.6;
const I0: f64 = 2.0;
const I2: f64 = 0.04;

fn model(axial: f64, bending: f64, shear: f64, interior_only: bool) -> Arc<MitcShellModel> {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, -0.6, 0.0], [1.0, 0.6, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    // Isotropic nu=0 makes the global axial oracle independent of the native
    // facet tangent frame. The director/bending operators may still differ.
    let a = SheetMatrix::orthotropic(axial, axial, 0.0, 0.5 * axial).unwrap();
    let d = SheetMatrix::orthotropic(bending, 0.6 * bending, 0.17, 0.2 * bending).unwrap();
    let section = ShearShellSection::new_coupled(
        a,
        [0.0; 9],
        d,
        [shear, 0.12 * shear, 0.8 * shear],
        1.0,
        I0,
        0.0,
        I2,
        0.0,
    )
    .unwrap();
    let constraints: Vec<_> = (0..3)
        .map(|node| ShellConstraint {
            node,
            mask: if interior_only {
                // One algebraic drilling coordinate remains free. All physical
                // nodal motion is fixed, retaining the two interior rotations.
                [true, true, true, true, true, node != 0]
            } else {
                // Reflection exchanges free nodes 1/2, so ux=q*x is an exact
                // invariant axial mode with m=I0*A/2 and k=A11*A.
                [node == 0, true, true, true, true, true]
            },
            frame: ConstraintFrame::Global,
        })
        .collect();
    Arc::new(
        assemble_mitc3p(
            &mesh,
            &directors,
            &[section],
            &constraints,
            ShellMass::Consistent,
            MitcBudget {
                max_nodes: 3,
                max_triangles: 1,
            },
            1e-3,
            || false,
        )
        .unwrap(),
    )
}

fn near(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn entry(matrix: &Csr, row: usize, column: usize) -> f64 {
    let (cols, vals) = matrix.row(row);
    cols.iter()
        .zip(vals)
        .find_map(|(&col, &value)| (col == column).then_some(value))
        .unwrap_or(0.0)
}

#[test]
fn g1_softening_axial_mode_obeys_independent_endpoint_momentum_and_parameter_work() {
    let a0 = 20_000.0;
    let a1 = 8_000.0;
    let before = EnrichedMitcShell::try_new(model(a0, 20.0, 120.0, false), || false).unwrap();
    let after = before
        .with_material_model(model(a1, 11.0, 75.0, false), || false)
        .unwrap();
    let (q0, v0, f0, f1, dt) = (1e-4, 0.07, 0.3, -0.4, 0.003);
    let mut displacement = [0.0; 18];
    let mut velocity = [0.0; 18];
    for dof in [6, 12] {
        displacement[dof] = q0;
        velocity[dof] = v0;
    }
    let initial = before
        .state(
            0.0,
            &displacement,
            &velocity,
            &[[0.0; 2]],
            &[[0.0; 2]],
            || false,
        )
        .unwrap();
    assert_eq!(before.dimension(), 4);
    let load0 = [0.5 * f0, 0.5 * f0, 0.0, 0.0];
    let load1 = [0.5 * f1, 0.5 * f1, 0.0, 0.0];
    let step = after
        .prepare_midpoint(dt, || false)
        .unwrap()
        .advance_from(&before, &initial, &load0, &load1, || false)
        .unwrap();

    // Direct scalar trapezoidal momentum, with independently integrated mass
    // and strain stiffness. A refreshed static map cannot reproduce this
    // retained kinetic state by merely re-solving endpoint equilibrium.
    let (mass, k0, k1, h) = (I0 * AREA / 2.0, a0 * AREA, a1 * AREA, 0.5 * dt);
    let vmid = (mass * v0 + h * (0.5 * (f0 + f1) - 0.5 * (k0 + k1) * q0)) / (mass + h * h * k1);
    let q1 = q0 + dt * vmid;
    let v1 = 2.0 * vmid - v0;
    for dof in [6, 12] {
        near(step.state.nodal_displacement_6n()[dof], q1, 2e-12, 1e-15);
        near(step.state.nodal_velocity_6n()[dof], v1, 2e-12, 1e-14);
    }
    let work = 0.5 * (f0 + f1) * (q1 - q0);
    let parameter_work = 0.5 * q0 * (k1 - k0) * q1;
    let elastic_change = 0.5 * (k1 * q1 * q1 - k0 * q0 * q0);
    let kinetic_change = 0.5 * mass * (v1 * v1 - v0 * v0);
    near(step.report.load_work_j, work, 3e-12, 1e-15);
    near(
        step.report.stiffness_parameter_work_j,
        parameter_work,
        3e-12,
        1e-15,
    );
    near(
        step.report.quadratic_energy_change_j,
        elastic_change,
        3e-12,
        1e-15,
    );
    near(
        step.report.kinetic_energy_change_j,
        kinetic_change,
        3e-12,
        1e-15,
    );
    near(step.report.energy_residual_j, 0.0, 0.0, 1e-14);
    near(
        kinetic_change + elastic_change,
        work + parameter_work,
        3e-12,
        1e-15,
    );
    assert!(parameter_work < 0.0);
    assert_eq!(initial.nodal_displacement_6n(), &displacement);
    assert_eq!(initial.nodal_velocity_6n(), &velocity);
}

#[test]
fn g1_internal_modes_retain_velocity_under_independent_bending_and_shear_changes() {
    let before = EnrichedMitcShell::try_new(model(20_000.0, 20.0, 120.0, true), || false).unwrap();
    let after = before
        .with_material_model(model(15_000.0, 8.0, 180.0, true), || false)
        .unwrap();
    assert_eq!(before.dimension(), 3);
    let q0 = [0.001, -0.0007];
    let v0 = [0.013, -0.009];
    let initial = before
        .state(0.0, &[0.0; 18], &[0.0; 18], &[q0], &[v0], || false)
        .unwrap();
    let load0 = [0.0, 0.1, -0.03];
    let load1 = [0.0, -0.02, 0.07];
    let dt = 0.001;
    let h = dt * 0.5;
    let mass = AREA * I2 * 81.0 / 280.0;
    let k0: [[f64; 2]; 2] =
        core::array::from_fn(|i| core::array::from_fn(|j| entry(before.stiffness(), i + 1, j + 1)));
    let k1: [[f64; 2]; 2] =
        core::array::from_fn(|i| core::array::from_fn(|j| entry(after.stiffness(), i + 1, j + 1)));
    let rhs: [f64; 2] = core::array::from_fn(|i| {
        mass * v0[i]
            + h * (0.5 * (load0[i + 1] + load1[i + 1])
                - 0.5 * (0..2).map(|j| (k0[i][j] + k1[i][j]) * q0[j]).sum::<f64>())
    });
    // Independent closed 2-by-2 solve of the retained physical interior pair;
    // mass follows the exact triangle integral of its cubic bubble squared.
    let aa = mass + h * h * k1[0][0];
    let bb = h * h * k1[0][1];
    let cc = mass + h * h * k1[1][1];
    let determinant = aa * cc - bb * bb;
    let vmid = [
        (cc * rhs[0] - bb * rhs[1]) / determinant,
        (aa * rhs[1] - bb * rhs[0]) / determinant,
    ];
    let q1: [f64; 2] = core::array::from_fn(|i| q0[i] + dt * vmid[i]);
    let v1: [f64; 2] = core::array::from_fn(|i| 2.0 * vmid[i] - v0[i]);
    let prepared = after.prepare_midpoint(dt, || false).unwrap();
    let step = prepared
        .advance_from(&before, &initial, &load0, &load1, || false)
        .unwrap();
    for i in 0..2 {
        near(step.state.bubble_coordinates()[0][i], q1[i], 3e-12, 1e-15);
        near(step.state.bubble_velocity()[0][i], v1[i], 3e-12, 1e-14);
    }
    let parameter_work = 0.5
        * (0..2)
            .map(|i| q0[i] * (0..2).map(|j| (k1[i][j] - k0[i][j]) * q1[j]).sum::<f64>())
            .sum::<f64>();
    near(
        step.report.stiffness_parameter_work_j,
        parameter_work,
        5e-12,
        1e-15,
    );
    near(
        step.state.kinetic_energy_j(),
        0.5 * mass * v1.iter().map(|x| x * x).sum::<f64>(),
        3e-12,
        1e-15,
    );
    near(step.report.energy_residual_j, 0.0, 0.0, 1e-14);
    assert!(
        step.state.bubble_velocity()[0]
            .iter()
            .any(|v| v.abs() > 1e-5)
    );

    let mut count = 0;
    prepared
        .advance_from(&before, &initial, &load0, &load1, || {
            count += 1;
            false
        })
        .unwrap();
    let mut cancelled_count = 0;
    assert!(matches!(
        prepared.advance_from(&before, &initial, &load0, &load1, || {
            cancelled_count += 1;
            cancelled_count == count
        }),
        Err(PlateError::Cancelled)
    ));
    assert_eq!(initial.bubble_coordinates(), &[q0]);
    assert_eq!(initial.bubble_velocity(), &[v0]);
}

#[test]
fn g0_fixed_material_entry_retains_one_kernel_and_exact_zero_parameter_work() {
    let shell = EnrichedMitcShell::try_new(model(20_000.0, 20.0, 120.0, true), || false).unwrap();
    let initial = shell
        .state(
            0.0,
            &[0.0; 18],
            &[0.0; 18],
            &[[0.001, -0.002]],
            &[[0.03, 0.01]],
            || false,
        )
        .unwrap();
    let prepared = shell.prepare_midpoint(0.002, || false).unwrap();
    let load0 = [0.0, 0.1, 0.03];
    let load1 = [0.0, 0.02, -0.07];
    let old_entry = prepared
        .advance(&initial, &load0, &load1, || false)
        .unwrap();
    let new_entry = prepared
        .advance_from(&shell.clone(), &initial, &load0, &load1, || false)
        .unwrap();
    assert_eq!(old_entry.state.coordinates(), new_entry.state.coordinates());
    assert_eq!(old_entry.state.velocity(), new_entry.state.velocity());
    assert_eq!(old_entry.report, new_entry.report);
    assert_eq!(new_entry.report.stiffness_parameter_work_j, 0.0);
}

#[test]
fn g0_changed_owner_incomplete_loads_and_director_axis_effort_refuse_atomically() {
    let original_model = model(20_000.0, 20.0, 120.0, true);
    let before = EnrichedMitcShell::try_new(Arc::clone(&original_model), || false).unwrap();
    let after = before
        .with_material_model(model(18_000.0, 12.0, 150.0, true), || false)
        .unwrap();
    let initial = before
        .state(
            0.0,
            &[0.0; 18],
            &[0.0; 18],
            &[[0.0; 2]],
            &[[0.0; 2]],
            || false,
        )
        .unwrap();
    let prepared = after.prepare_midpoint(0.001, || false).unwrap();
    // Even an independently reconstructed owner of the same nodal model does
    // not share the original prepared full physical inertia/basis identity.
    let foreign = EnrichedMitcShell::try_new(original_model, || false).unwrap();
    assert!(
        prepared
            .advance_from(&foreign, &initial, &[0.0; 3], &[0.0; 3], || false)
            .is_err()
    );
    assert!(
        prepared
            .advance_from(&before, &initial, &[0.0; 2], &[0.0; 3], || false)
            .is_err()
    );
    assert!(
        prepared
            .advance_from(&before, &initial, &[0.0; 3], &[0.0, f64::NAN, 0.0], || {
                false
            })
            .is_err()
    );
    assert!(
        prepared
            .advance_from(&before, &initial, &[0.0; 3], &[1e-3, 0.0, 0.0], || false)
            .is_err()
    );
    assert!(matches!(
        prepared.advance_from(&before, &initial, &[0.0; 3], &[0.0; 3], || true),
        Err(PlateError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.coordinates(), &[0.0; 3]);
    assert_eq!(initial.velocity(), &[0.0; 3]);
}
