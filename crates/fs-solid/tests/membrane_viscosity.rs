//! Synthetic membrane damping checks; these are not measured sheet parameters.
use fs_evidence::ValidityDomain;
use fs_material::plane_stress_viscosity::{
    MaterialPlane, PlaneStressViscosity, PLANE_STRESS_VISCOSITY_PROPERTIES as KEYS,
};
use fs_material::sheet::{AREAL_DENSITY, MEMBRANE_DIMS, ResolvedSheet, THICKNESS};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Dims, DynViscosity};
use fs_solid::membrane::vibration::{MembraneModes, VibrationSettings};
use fs_solid::membrane::viscosity::{MembraneViscosity, MembraneViscosityError};
use fs_solid::membrane::{MembraneBudget, MembraneDomain, MembraneError, MembranePatch};
use fs_solid::mesh2::Mesh2;

fn patch(
    values: [f64; 6],
    thickness: f64,
    grain: f64,
) -> (MembranePatch, Vec<PlaneStressViscosity>) {
    let mut properties = vec![
        (AREAL_DENSITY, 6.0, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, thickness, Dims([1, 0, 0, 0, 0, 0])),
        ("warp_stiffness", 1200.0, MEMBRANE_DIMS),
        ("weft_stiffness", 600.0, MEMBRANE_DIMS),
        ("poisson_warp_weft", 0.0, Dims::NONE),
        ("shear_stiffness", 250.0, MEMBRANE_DIMS),
    ];
    properties.extend(KEYS.into_iter().zip(values).map(|(key, value)| {
        (key, value, DynViscosity::DIMS)
    }));
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, value, dims) in properties {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Scalar { value, dims },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            validity: ValidityDomain::unconstrained()
                .with("T", 290.0, 310.0).with("omega", 0.0, 1000.0),
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic membrane viscosity; no specimen authority".into(),
                license: "test".into(),
                artifact: None,
            },
        }).unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic sheet".into(),
        phase: "solid".into(),
        process: "damping test".into(),
        revision: 0,
    }, claims, vec![]).unwrap();
    let source = resolve_material_state_point(
        &card, &QueryPoint::new().with("T", 300.0).unwrap().with("omega", 0.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let law = PlaneStressViscosity::from_material(
        source.clone(), MaterialPlane::Axes12, KEYS,
    ).unwrap();
    let sheet = ResolvedSheet::try_from_material(source).unwrap();
    let patch = MembranePatch::try_new(
        Mesh2 {
            nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            elems: vec![vec![0, 1, 2]],
            patches: vec![],
        },
        vec![sheet], vec![grain],
        MembraneDomain { minimum_area_ratio: 0.01, maximum_abs_strain: 0.5 },
        MembraneBudget { max_nodes: 3, max_triangles: 1 },
        || false,
    ).unwrap();
    (patch, vec![law])
}

fn deformed() -> Vec<[f64; 3]> {
    vec![[0.0; 3], [0.08, 0.12, 0.19], [-0.04, -0.06, 0.25]]
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-14 + 3e-11 * a.abs().max(b.abs()), "{a:e} != {b:e}");
}

fn modes<'a>(
    patch: &'a MembranePatch,
    displacement: &[[f64; 3]],
    dofs: &[usize],
) -> MembraneModes<'a> {
    let n = dofs.len();
    let mut stiffness = vec![0.0; n * n];
    let mut direction = vec![[0.0; 3]; patch.mesh().nodes.len()];
    for (column, &dof) in dofs.iter().enumerate() {
        direction[dof / 3][dof % 3] = 1.0;
        let action = patch.tangent_apply(displacement, &direction, || false).unwrap();
        direction[dof / 3][dof % 3] = 0.0;
        for (row, &index) in dofs.iter().enumerate() {
            stiffness[row * n + column] = action[index / 3][index % 3];
        }
    }
    patch.modes_from_tangent(
        displacement, dofs, &stiffness,
        VibrationSettings { max_free_dofs: 16, relative_residual_tolerance: 1e-8 },
        || false,
    ).unwrap()
}

#[test]
fn curved_sheet_matches_independent_green_rate_force_and_power() {
    let (p, laws) = patch([4.0, 3.0, 2.0, 1.0, 0.4, 0.2], 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let velocity = [[0.0; 3], [0.0, 0.0, 1.3], [0.0, 0.0, -0.7]];
    let response = damping.evaluate(&deformed(), &velocity, || false).unwrap();
    // A*h = .001 and Green-rate columns [.19,0,.25], [0,.25,.19].
    // Independent transverse C = [.0003074,.00016944;.00016944,.0002787].
    close(response.resisting_force_n[1][2], 0.000281012);
    close(response.resisting_force_n[2][2], 0.000025182);
    close(response.power_w, 0.0003476882);
    close(response.element_power_w[0], response.power_w);
    let nodal_work: f64 = velocity.iter().flatten()
        .zip(response.resisting_force_n.iter().flatten()).map(|(v, f)| v * f).sum();
    close(nodal_work, response.power_w);
    for axis in 0..3 {
        close(response.resisting_force_n.iter().map(|f| f[axis]).sum(), 0.0);
    }

    // Membrane loss scales with h, independently of areal inertia or elastic A.
    let (thicker, laws) = patch([4.0, 3.0, 2.0, 1.0, 0.4, 0.2], 0.004, 0.0);
    let twice = MembraneViscosity::new(&thicker, laws, (0.0, 1000.0), || false)
        .unwrap().evaluate(&deformed(), &velocity, || false).unwrap();
    close(twice.power_w, 2.0 * response.power_w);
    for (a, b) in twice.resisting_force_n.iter().flatten()
        .zip(response.resisting_force_n.iter().flatten())
    {
        close(*a, 2.0 * b);
    }
    assert_eq!(thicker.reference_mass_kg(), p.reference_mass_kg());
    close(thicker.evaluate(&deformed(), || false).unwrap().strain_energy_j,
        p.evaluate(&deformed(), || false).unwrap().strain_energy_j);
}

#[test]
fn rigid_translation_spin_and_rotated_observations_do_not_create_loss() {
    let (p, laws) = patch([4.0, 3.0, 2.0, 1.0, 0.4, 0.2], 0.002, 0.37);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let u = deformed();
    let translation = damping.evaluate(&u, &[[0.7, -0.4, 0.2]; 3], || false).unwrap();
    assert_eq!(translation.power_w, 0.0);
    assert!(translation.resisting_force_n.iter().flatten().all(|f| *f == 0.0));

    let positions: Vec<[f64; 3]> = p.mesh().nodes.iter().zip(&u)
        .map(|(x, u)| [x[0] + u[0], x[1] + u[1], u[2]]).collect();
    let spin: Vec<_> = positions.iter().map(|x| [
        -1.1 * x[2] - 0.4 * x[1],
        0.4 * x[0] - 0.7 * x[2],
        0.7 * x[1] + 1.1 * x[0],
    ]).collect();
    let rigid = damping.evaluate(&u, &spin, || false).unwrap();
    assert!(rigid.power_w.abs() < 1e-28);
    assert!(rigid.resisting_force_n.iter().flatten().all(|f| f.abs() < 1e-14));

    // A proper finite 3-D rotation of both current geometry and velocity.
    let rotate = |x: [f64; 3]| [x[2], x[0], x[1]];
    let velocity = [[0.2, -0.3, 0.1], [-0.4, 0.1, 0.8], [0.3, 0.6, -0.2]];
    let before = damping.evaluate(&u, &velocity, || false).unwrap();
    let rotated_u: Vec<_> = positions.iter().zip(&p.mesh().nodes).map(|(x, reference)| {
        let x = rotate(*x);
        [x[0] + 0.2 - reference[0], x[1] - 0.3 - reference[1], x[2] + 0.7]
    }).collect();
    let after = damping.evaluate(&rotated_u, &velocity.map(rotate), || false).unwrap();
    close(after.power_w, before.power_w);
    for (a, b) in after.resisting_force_n.iter().zip(&before.resisting_force_n) {
        for (a, b) in a.iter().zip(rotate(*b)) { close(*a, b); }
    }
}

#[test]
fn directional_nulls_and_prestress_are_not_replaced_by_artificial_loss() {
    let (p, laws) = patch([4.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let flat = [[0.0; 3]; 3];
    let weft_rate = [[0.0; 3], [0.0; 3], [0.0, 1.0, 0.0]];
    assert_eq!(damping.evaluate(&flat, &weft_rate, || false).unwrap().power_w, 0.0);
    let warp_rate = [[0.0; 3], [1.0, 0.0, 0.0], [0.0; 3]];
    close(damping.evaluate(&flat, &warp_rate, || false).unwrap().power_w, 0.004);

    let (p, laws) = patch([4.0, 3.0, 2.0, 1.0, 0.4, 0.2], 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let stretch = [[0.0; 3], [0.02, 0.0, 0.0], [0.0, 0.02, 0.0]];
    let transverse = [[0.0; 3], [0.0, 0.0, 1.0], [0.0; 3]];
    assert!(p.tangent_apply(&stretch, &transverse, || false).unwrap()[1][2] > 0.0);
    let response = damping.evaluate(&stretch, &transverse, || false).unwrap();
    assert_eq!(response.power_w, 0.0);
    assert!(response.resisting_force_n.iter().flatten().all(|f| *f == 0.0));
    let basis = modes(&p, &stretch, &[5]);
    let projected = damping.project_modes(&basis, 100_000, || false).unwrap();
    assert_eq!(projected.matrix_per_s(), &[0.0]);
    assert!(projected.dissipation_rows().iter().all(|x| *x == 0.0));
}

#[test]
fn full_modal_projection_keeps_original_work_cross_couplings_and_geometry() {
    let (p, laws) = patch([4.0, 3.0, 2.0, 1.0, 0.4, 0.2], 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let u = deformed();
    let basis = modes(&p, &u, &[5, 8]);
    assert_eq!(basis.linearization_displacement_m(), u.as_slice());
    let projected = damping.project_modes(&basis, 100_000, || false).unwrap();
    assert_eq!(projected.mode_count(), 2);
    assert_eq!(projected.row_count(), 6);
    let matrix = projected.matrix_per_s();
    assert_eq!(matrix[1], matrix[2]);
    assert!(matrix[1].abs() > 1e-5);
    let physical = [0.0003074, 0.00016944, 0.00016944, 0.0002787];
    for i in 0..2 {
        for j in 0..2 {
            let expected: f64 = (0..2).map(|a| (0..2).map(|b| {
                basis.modes()[i].phi[a] * physical[2 * a + b] * basis.modes()[j].phi[b]
            }).sum::<f64>()).sum();
            close(matrix[2 * i + j], expected);
            let rows: f64 = (0..projected.row_count()).map(|row| {
                projected.dissipation_rows()[2 * row + i]
                    * projected.dissipation_rows()[2 * row + j]
            }).sum();
            close(rows, expected);
        }
    }
    let qdot = [1.3, -0.7];
    let velocity = basis.reconstruct(&qdot).unwrap();
    let response = damping.evaluate(&u, &velocity, || false).unwrap();
    let modal_power = qdot[0] * (matrix[0] * qdot[0] + matrix[1] * qdot[1])
        + qdot[1] * (matrix[2] * qdot[0] + matrix[3] * qdot[1]);
    close(modal_power, response.power_w);
    let diagonal_power = matrix[0] * qdot[0].powi(2) + matrix[3] * qdot[1].powi(2);
    assert!((diagonal_power - modal_power).abs() > 0.1 * modal_power);
    let repeated = damping.project_modes(&basis, 100_000, || false).unwrap();
    assert_eq!(matrix, repeated.matrix_per_s());
    assert_eq!(projected.dissipation_rows(), repeated.dissipation_rows());
}

#[test]
fn source_axes_band_budget_foreign_geometry_and_cancelled_trials_refuse() {
    let values = [4.0, 3.0, 2.0, 1.0, 0.4, 0.2];
    let (p, laws) = patch(values, 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws.clone(), (0.0, 1000.0), || false).unwrap();
    assert!(MembraneViscosity::new(&p, vec![], (0.0, 1000.0), || false).is_err());
    assert!(MembraneViscosity::new(&p, laws.clone(), (0.0, 1001.0), || false).is_err());
    assert!(MembraneViscosity::new(&p, laws.clone(), (10.0, 1.0), || false).is_err());
    let wrong_axes = PlaneStressViscosity::from_material(
        p.sheets()[0].material().clone(), MaterialPlane::Axes21, KEYS,
    ).unwrap();
    assert!(MembraneViscosity::new(&p, vec![wrong_axes], (0.0, 1000.0), || false).is_err());
    let (_, foreign_source) = patch(values.map(|v| 2.0 * v), 0.002, 0.0);
    assert!(MembraneViscosity::new(&p, foreign_source, (0.0, 1000.0), || false).is_err());

    let u = deformed();
    let basis = modes(&p, &u, &[5, 8]);
    let narrow = MembraneViscosity::new(
        &p, laws, (0.0, 0.5 * basis.modes()[0].lambda.sqrt()), || false,
    ).unwrap();
    assert!(narrow.project_modes(&basis, 100_000, || false).is_err());
    let copy = p.clone();
    let foreign_modes = modes(&copy, &u, &[5, 8]);
    assert!(damping.project_modes(&foreign_modes, 100_000, || false).is_err());
    assert!(matches!(damping.project_modes(&basis, 0, || false),
        Err(MembraneViscosityError::Membrane(MembraneError::Budget(_)))));
    assert!(matches!(damping.project_modes(&basis, 100_000, || true),
        Err(MembraneViscosityError::Membrane(MembraneError::Cancelled))));
    assert!(damping.evaluate(&u, &[[0.0; 3]; 2], || false).is_err());
    let collapsed = [[0.0; 3], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]];
    assert!(damping.evaluate(&collapsed, &[[0.0; 3]; 3], || false).is_err());
    let mut polls = 0;
    assert!(matches!(damping.evaluate(&u, &[[0.0; 3]; 3], || {
        polls += 1;
        polls == 3
    }), Err(MembraneViscosityError::Membrane(MembraneError::Cancelled))));
    assert_eq!(basis.linearization_displacement_m(), u.as_slice());
    damping.evaluate(&u, &[[0.0; 3]; 3], || false).unwrap();
}

#[path = "membrane_viscosity/modal.rs"]
mod modal;

#[path = "membrane_viscosity/natural_reference.rs"]
mod natural_reference;
