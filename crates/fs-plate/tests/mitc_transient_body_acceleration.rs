//! Physical graded reference mass drives every retained director coordinate.

use fs_material::sheet::SheetMatrix;
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::transient::EnrichedMitcShell;
use fs_plate::shell::mitc::thermal::{FacetThermalStrain, ShellThermalStrain};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellInertia,
    ShearShellSection, ShellConstraint, ShellMass, assemble_mitc3p_with_reference_inertia,
};
use std::sync::Arc;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn rotate(a: [f64; 3]) -> [f64; 3] {
    [[0.36, -0.48, 0.8], [0.8, 0.6, 0.0], [-0.48, 0.64, 0.6]].map(|r| dot(r, a))
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-11 * a.abs().max(b.abs()).max(1e-3),
        "observed={a:e}, expected={b:e}"
    );
}
fn model(
    supported: bool,
    rotated: bool,
    density_scale: f64,
    stiffness_scale: f64,
) -> Arc<MitcShellModel> {
    let transform = |v| if rotated { rotate(v) } else { v };
    let mesh = ShellMesh::new(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
            .map(transform)
            .to_vec(),
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![transform([0.0, 0.0, 1.0]); 3], 0.9).unwrap();
    let section = ShearShellSection::new_coupled(
        SheetMatrix::orthotropic(2e4 * stiffness_scale, 1e4, 0.2, 3e3).unwrap(),
        [0.0; 9],
        SheetMatrix::orthotropic(20.0, 10.0 * stiffness_scale, 0.2, 3.0).unwrap(),
        [120.0 * stiffness_scale, 0.0, 80.0],
        1.0,
        2.0,
        0.01,
        0.002,
        0.0,
    )
    .unwrap();
    let moments = [[1.0, 2.0, 5.0].map(|rho| {
        let rho = rho * density_scale;
        ShearShellInertia::new(rho, 0.01 * rho, 0.001 * rho).unwrap()
    })];
    let constraints: Vec<_> = if supported {
        (0..3)
            .map(|node| ShellConstraint {
                node,
                mask: [node != 0, true, true, true, true, true],
                frame: if rotated {
                    ConstraintFrame::Local([
                        transform([1.0, 0.0, 0.0]),
                        transform([0.0, 1.0, 0.0]),
                        transform([0.0, 0.0, 1.0]),
                    ])
                } else {
                    ConstraintFrame::Global
                },
            })
            .collect()
    } else {
        vec![]
    };
    Arc::new(
        assemble_mitc3p_with_reference_inertia(
            &mesh,
            &directors,
            &[section],
            &moments,
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

#[test]
fn graded_body_force_retains_constrained_columns_and_independent_director_moments() {
    let g = [3.0, -2.0, 0.5];
    for rotated in [false, true] {
        let original = model(true, rotated, 1.0, 1.0);
        let shell = EnrichedMitcShell::try_new(Arc::clone(&original), || false).unwrap();
        let world_g = if rotated { rotate(g) } else { g };
        let load = shell.body_acceleration_load(world_g, || false).unwrap();
        assert_eq!(load.len(), 3);
        // On the unit right triangle A=1/2, independently:
        // integral N0*I0 = A*(2*1+2+5)/12 = 3/8.
        // b=27*N0*N1*N2, integral b*I1 = 3*A*sum(I1_i)/20 = 0.006.
        // The retained axes are ex,ey: ex cross n=-ey; ey cross n=ex.
        near(load[0], 0.375 * g[0]);
        near(load[1], -0.006 * g[1]);
        near(load[2], 0.006 * g[0]);

        let mut restricted_rigid = original
            .project_load(
                &(0..3)
                    .flat_map(|_| [world_g[0], world_g[1], world_g[2], 0.0, 0.0, 0.0])
                    .collect::<Vec<_>>(),
            )
            .unwrap();
        restricted_rigid.extend([0.0, 0.0]);
        let mut wrong = vec![0.0; 3];
        shell.mass().spmv(&restricted_rigid, &mut wrong);
        assert!((wrong[0] - load[0]).abs() > 0.5);
        assert!((wrong[1] - load[1]).abs() + (wrong[2] - load[2]).abs() > 0.005);
    }
}

#[test]
fn complete_mass_action_and_support_projection_agree_with_original_reference_inertia() {
    let g = [0.3, -0.7, -9.81];
    for rotated in [false, true] {
        let world_g = if rotated { rotate(g) } else { g };
        let original = model(false, rotated, 1.0, 1.0);
        let full = EnrichedMitcShell::try_new(Arc::clone(&original), || false).unwrap();
        let rigid: Vec<_> = (0..3)
            .flat_map(|_| [world_g[0], world_g[1], world_g[2], 0.0, 0.0, 0.0])
            .chain([0.0, 0.0])
            .collect();
        let mut expected = vec![0.0; full.dimension()];
        full.mass().spmv(&rigid, &mut expected);
        let observed = full.body_acceleration_load(world_g, || false).unwrap();
        for (&a, &b) in observed.iter().zip(&expected) {
            near(a, b);
        }
        for axis in 0..3 {
            let total: f64 = (0..3).map(|node| observed[6 * node + axis]).sum();
            near(total, (4.0 / 3.0) * world_g[axis]);
        }
        let constrained = model(true, rotated, 1.0, 1.0);
        let projected = constrained.project_load(&expected[..18]).unwrap();
        let reduced = EnrichedMitcShell::try_new(constrained, || false).unwrap();
        let load = reduced.body_acceleration_load(world_g, || false).unwrap();
        near(load[0], projected[0]);
        near(load[1], expected[18]);
        near(load[2], expected[19]);
    }
}

#[test]
fn body_work_drives_native_midpoint_motion_and_survives_material_refresh() {
    let original = model(true, false, 1.0, 1.0);
    let shell = EnrichedMitcShell::try_new(Arc::clone(&original), || false).unwrap();
    let refreshed = shell
        .with_material_model(model(true, false, 1.0, 1.7), || false)
        .unwrap();
    let g = [3.0, -2.0, 0.5];
    let force = shell.body_acceleration_load(g, || false).unwrap();
    assert_eq!(
        force,
        refreshed.body_acceleration_load(g, || false).unwrap()
    );
    let doubled = EnrichedMitcShell::try_new(model(true, false, 2.0, 1.0), || false).unwrap();
    for (&a, b) in force
        .iter()
        .zip(doubled.body_acceleration_load(g, || false).unwrap())
    {
        near(2.0 * a, b);
    }

    let thermal = refreshed
        .model()
        .prepare_thermal(
            &[FacetThermalStrain::Uniform(
                ShellThermalStrain::try_new([0.0; 3], [0.0; 3]).unwrap(),
            )],
            || false,
        )
        .unwrap();
    let load = refreshed
        .generalized_load_with_body_acceleration(&thermal, &[0.0; 18], g, || false)
        .unwrap();
    assert_eq!(load, force);
    let initial = shell
        .state(
            0.0,
            &[0.0; 18],
            &[0.0; 18],
            &[[0.0; 2]],
            &[[0.0; 2]],
            || false,
        )
        .unwrap();
    let step = refreshed
        .prepare_midpoint(0.002, || false)
        .unwrap()
        .advance(&initial, &load, &load, || false)
        .unwrap();
    let independent_work: f64 = step
        .state
        .coordinates()
        .iter()
        .zip([1.125, 0.012, 0.018])
        .map(|(q, f)| q * f)
        .sum();
    near(step.report.load_work_j, independent_work);
    near(
        step.report.kinetic_energy_change_j + step.report.quadratic_energy_change_j,
        independent_work,
    );
    assert!(step.report.energy_residual_j.abs() < 1e-11 * independent_work.abs().max(1e-6));
    assert!(
        step.state.bubble_coordinates()[0]
            .iter()
            .any(|v| v.abs() > 1e-10)
    );
    assert_eq!(initial.coordinates(), &[0.0; 3]);

    let strain = [FacetThermalStrain::Uniform(
        ShellThermalStrain::try_new([2e-5, -1e-5, 0.0], [0.0; 3]).unwrap(),
    )];
    let thermal = refreshed
        .model()
        .prepare_thermal(&strain, || false)
        .unwrap();
    let mut external = [0.0; 18];
    external[0] = 0.125;
    let old = refreshed
        .generalized_load(&thermal, &external, || false)
        .unwrap();
    let combined = refreshed
        .generalized_load_with_body_acceleration(&thermal, &external, g, || false)
        .unwrap();
    for ((actual, old), gravity) in combined.iter().zip(old).zip(force) {
        near(*actual, old + gravity);
    }
}

#[test]
fn invalid_acceleration_overflow_and_cancellation_do_not_change_retained_state() {
    let shell = EnrichedMitcShell::try_new(model(true, false, 10.0, 1.0), || false).unwrap();
    let initial = shell
        .state(
            0.0,
            &[0.0; 18],
            &[0.0; 18],
            &[[0.0; 2]],
            &[[0.0; 2]],
            || false,
        )
        .unwrap();
    assert!(
        shell
            .body_acceleration_load([f64::NAN, 0.0, 0.0], || false)
            .is_err()
    );
    assert!(
        shell
            .body_acceleration_load([f64::INFINITY, 0.0, 0.0], || false)
            .is_err()
    );
    assert!(
        shell
            .body_acceleration_load([f64::MAX, 0.0, 0.0], || false)
            .is_err()
    );
    assert!(matches!(
        shell.body_acceleration_load([1.0; 3], || true),
        Err(PlateError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        shell.body_acceleration_load([1.0; 3], || {
            polls += 1;
            polls == 2
        }),
        Err(PlateError::Cancelled)
    ));
    assert_eq!(
        shell.body_acceleration_load([0.0; 3], || false).unwrap(),
        vec![0.0; 3]
    );
    assert_eq!(initial.coordinates(), &[0.0; 3]);
    assert_eq!(initial.kinetic_energy_j(), 0.0);
}
