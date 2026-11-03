//! G1/G3/G4: native 3D electric loads versus independent Q E and frame laws.
use super::{
    ElectrostaticVoltageFace, SpatialElectrostaticBoundary, SpatialElectrostaticError,
    SpatialElectrostatics, WeightedAssemblyLimits, cube, dielectric, source,
};
use fs_em::spatial::force::SpatialMaxwellForceLimits;
use fs_ga::{Se3, So3, So3Tangent, Vec3, Wrench};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_PROPERTY,
    spatial::{PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES, SpatialDielectricLaw},
};

const EPSILON: f64 = 2e-6;
const RHO: f64 = 1e-3;
const EXTERNAL_E: f64 = 10.0;

fn assembly_limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 30_000,
        max_dofs: 100_000,
        max_triplets: 1_100_000,
    }
}
fn force_limits() -> SpatialMaxwellForceLimits {
    SpatialMaxwellForceLimits {
        max_cells: 30_000,
        max_faces: 100_000,
    }
}
fn controls() -> super::ElectrostaticControls {
    let mut controls = super::controls();
    controls.charge_tolerance_c = 1e-11;
    controls.energy_tolerance_j = 1e-9;
    controls
}
fn region(problem: &SpatialElectrostatics, radius: f64) -> Vec<usize> {
    problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .filter_map(|(cell, tet)| {
            (0..3)
                .all(|axis| {
                    let center: f64 = tet
                        .iter()
                        .map(|&v| problem.positions_m()[v as usize][axis] / 4.0)
                        .sum();
                    center.abs() < radius
                })
                .then_some(cell)
        })
        .collect()
}
fn specimen(n: usize) -> SpatialElectrostatics {
    let xs: Vec<_> = (0..=n).map(|i| i as f64 / n as f64).collect();
    let (points, tets) = cube(&xs, n, n);
    let materials = vec![dielectric(EPSILON); tets.len()];
    SpatialElectrostatics::new(
        points
            .into_iter()
            .map(|p| p.map(|x| 4.0 * x - 2.0))
            .collect(),
        tets,
        materials,
        assembly_limits(),
        &mut || false,
    )
    .unwrap()
}
fn boundary(problem: &SpatialElectrostatics) -> SpatialElectrostaticBoundary {
    SpatialElectrostaticBoundary {
        voltage: problem
            .exterior_faces()
            .map(|face| ElectrostaticVoltageFace {
                face,
                voltage_v: problem.mesh().faces[face]
                    .map(|v| -EXTERNAL_E * problem.positions_m()[v as usize][0]),
            })
            .collect(),
        displacement: Vec::new(),
    }
}
fn error(a: Vec3, b: Vec3) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}
fn same_wrench(actual: Wrench, expected: Wrench, tolerance: f64) {
    assert!(
        error(actual.force, expected.force) < tolerance,
        "force {:?} != {:?}",
        actual.force,
        expected.force
    );
    assert!(
        error(actual.torque, expected.torque) < tolerance,
        "torque {:?} != {:?}",
        actual.torque,
        expected.torque
    );
}

#[test]
fn enclosed_charge_force_refines_toward_qe_and_transforms_as_a_native_wrench() {
    // The centered charged unit cube has no net self-force by inversion symmetry.
    // The externally imposed uniform E therefore gives F=Q E=(.01,0,0) N,
    // independently of the exterior-side P1 stress quadrature.
    let expected = Vec3::new(RHO * EXTERNAL_E, 0.0, 0.0);
    let mut previous = None;
    for n in [8, 16] {
        let problem = specimen(n);
        let selected = region(&problem, 1.0);
        let mut rho = vec![0.0; problem.mesh().tets.len()];
        for cell in region(&problem, 0.5) {
            rho[cell] = RHO;
        }
        let field = problem
            .solve(&rho, &boundary(&problem), controls(), &mut || false)
            .unwrap();
        assert!((field.volume_charge_c() - RHO).abs() < 1e-15);
        let actual = field
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let actual_error = error(actual.force, expected);
        assert!(
            actual.force.x > 0.0,
            "force must follow positive charge and imposed E"
        );
        if let Some((coarse, coarse_error)) = previous {
            assert!(
                actual_error < 0.65 * coarse_error,
                "exterior stress must refine: fine={actual_error}, coarse={coarse_error}"
            );
            let coarse: Wrench = coarse;
            let extrapolated = Vec3::new(
                2.0 * actual.force.x - coarse.force.x,
                2.0 * actual.force.y - coarse.force.y,
                2.0 * actual.force.z - coarse.force.z,
            );
            // Test-only first-order extrapolation checks the known continuum
            // limit; the public API returns the unmodified surface integral.
            assert!(error(extrapolated, expected) < 0.03 * expected.x);
            assert!(actual_error < 0.2 * expected.x);
        }
        previous = Some((actual, actual_error));

        let mut reversed = selected.clone();
        reversed.reverse();
        assert_eq!(
            actual,
            field
                .maxwell_wrench(
                    &reversed,
                    &Se3::identity(),
                    [0.0; 3],
                    force_limits(),
                    &mut || false,
                )
                .unwrap()
        );
        let placement = Se3::from_parts(
            So3::exp(So3Tangent {
                angular: Vec3::new(0.3, -0.2, 0.4),
            })
            .unwrap(),
            Vec3::new(4.0, -3.0, 2.0),
        )
        .unwrap();
        let origin = Vec3::new(0.4, -0.7, 0.2);
        let transformed = actual.transform_by(&placement);
        let shift = origin.cross(transformed.force);
        let expected_frame = Wrench::new(
            Vec3::new(
                transformed.torque.x - shift.x,
                transformed.torque.y - shift.y,
                transformed.torque.z - shift.z,
            ),
            transformed.force,
        );
        let placed = field
            .maxwell_wrench(
                &selected,
                &placement,
                [origin.x, origin.y, origin.z],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        same_wrench(placed, expected_frame, 1e-12);
        assert!(placed.torque.norm() > 1e-3);
    }
}

#[test]
fn homogeneous_uniform_field_is_load_free_and_refusals_leave_it_unchanged() {
    let problem = specimen(8);
    let selected = region(&problem, 1.0);
    let rho = vec![0.0; problem.mesh().tets.len()];
    let field = problem
        .solve(&rho, &boundary(&problem), controls(), &mut || false)
        .unwrap();
    let accepted = field
        .maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    same_wrench(accepted, Wrench::default(), 1e-12);
    let mut duplicate = selected.clone();
    duplicate.push(selected[0]);
    for invalid in [vec![], vec![0], duplicate, region(&problem, 1.5)] {
        assert!(
            field
                .maxwell_wrench(
                    &invalid,
                    &Se3::identity(),
                    [0.0; 3],
                    force_limits(),
                    &mut || false,
                )
                .is_err()
        );
    }
    assert!(
        field
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                SpatialMaxwellForceLimits {
                    max_cells: problem.mesh().tets.len() - 1,
                    ..force_limits()
                },
                &mut || false,
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        field.maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || {
                polls += 1;
                polls > selected.len() + 25
            },
        ),
        Err(SpatialElectrostaticError::Cancelled)
    ));
    assert_eq!(
        accepted,
        field
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap()
    );
}

#[test]
fn stress_collar_refuses_charge_different_complete_sources_and_orthotropic_laws() {
    let problem = specimen(8);
    let selected = region(&problem, 1.0);
    let mut rho = vec![0.0; problem.mesh().tets.len()];
    rho[selected[0]] = RHO;
    let charged = problem
        .solve(&rho, &boundary(&problem), controls(), &mut || false)
        .unwrap();
    assert!(
        charged
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .is_err()
    );

    // Same epsilon is insufficient: one collar cell has a different complete
    // source bundle, so the original homogeneous-medium assumption is absent.
    let mut materials = problem.materials().to_vec();
    materials[selected[0]] = source(&[
        (ELECTRIC_PERMITTIVITY_PROPERTY, EPSILON),
        ("auxiliary_permittivity", EPSILON),
    ]);
    let changed = SpatialElectrostatics::new(
        problem.positions_m().to_vec(),
        problem.mesh().tets.clone(),
        materials,
        assembly_limits(),
        &mut || false,
    )
    .unwrap();
    let field = changed
        .solve(
            &vec![0.0; rho.len()],
            &boundary(&changed),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(
        field
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .is_err()
    );

    let mut properties = vec![(ELECTRIC_PERMITTIVITY_PROPERTY, EPSILON)];
    properties.extend(PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES.map(|name| (name, EPSILON)));
    let mut laws = vec![SpatialDielectricLaw::Isotropic; rho.len()];
    laws[selected[0]] =
        SpatialDielectricLaw::orthotropic([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
    let oriented = SpatialElectrostatics::new_with_laws(
        problem.positions_m().to_vec(),
        problem.mesh().tets.clone(),
        vec![source(&properties); rho.len()],
        laws,
        assembly_limits(),
        &mut || false,
    )
    .unwrap();
    let field = oriented
        .solve(
            &vec![0.0; rho.len()],
            &boundary(&oriented),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(
        field
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .is_err()
    );
}
