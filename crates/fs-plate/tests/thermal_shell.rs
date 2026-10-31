//! G1 thermal patch and G3 work/rotation checks on the native MITC shell.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::sheet::SheetMatrix;
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
    integrate_isotropic_thermal_expansion, integrate_orthotropic_thermal_expansion,
};
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::{FacetThermalStrain, ShellThermalStrain};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_qty::Temperature;

fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() < 2e-8 * scale.max(a.abs()).max(b.abs()),
        "{a:e} != {b:e}"
    );
}

fn section() -> ShearShellSection {
    ShearShellSection::from_orthotropic_solid(
        2e8,
        2e8,
        0.25,
        8e7,
        8e7,
        8e7,
        0.1,
        900.0,
        0.0,
        5.0 / 6.0,
    )
    .unwrap()
}

fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.6, 0.0],
            [1.0, 0.6, 0.0],
        ],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap()
}

fn model(mesh: &ShellMesh, section: ShearShellSection, fixed: bool) -> MitcShellModel {
    let normal = mesh.facet(0).unwrap().frame[2];
    let directors = NodalDirectors::new(mesh, vec![normal; mesh.nodes.len()], 0.9).unwrap();
    let constraints = if fixed {
        vec![ShellConstraint {
            node: 0,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        }]
    } else {
        vec![]
    };
    assemble_mitc3p(
        mesh,
        &directors,
        &vec![section; mesh.tris.len()],
        &constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 4,
            max_triangles: 2,
        },
        1e-3,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_free_heating_and_cooling_change_geometry_without_stress_or_mass() {
    let mesh = mesh();
    let model = model(&mesh, section(), true);
    let mass = model.pencil().m.clone();
    for expansion in [0.0012, -0.0008] {
        let strain = ShellThermalStrain::try_new([expansion, expansion, 0.0], [0.0; 3]).unwrap();
        let prepared = model
            .prepare_thermal(&[strain.into(); 2], || false)
            .unwrap();
        let response = prepared.solve_static(&[0.0; 24], || false).unwrap();
        for (node, p) in mesh.nodes.iter().enumerate() {
            near(response.displacement_6n[6 * node], expansion * p[0], 0.001);
            near(
                response.displacement_6n[6 * node + 1],
                expansion * p[1],
                0.001,
            );
        }
        assert!(response.strain_energy_j < 1e-18);
        for facet in response.facets {
            assert!(facet.membrane_force_n_m.iter().all(|v| v.abs() < 1e-6));
        }
    }
    for row in 0..model.pencil().free {
        assert_eq!(mass.row(row), model.pencil().m.row(row));
    }
}

#[test]
fn g1_restrained_heating_recovers_plane_stress_force_and_storage() {
    let mesh = mesh();
    let model = model(&mesh, section(), false);
    let expansion = 0.0012;
    let strain = ShellThermalStrain::try_new([expansion, expansion, 0.0], [0.0; 3]).unwrap();
    let prepared = model
        .prepare_thermal(&[strain.into(); 2], || false)
        .unwrap();
    let response = prepared.recover(&[0.0; 24], || false).unwrap();
    // Independent plane-stress biaxial compression: N=-E*h*eps/(1-nu).
    let n = -2e8 * 0.1 * expansion / (1.0 - 0.25);
    for facet in response.facets {
        near(facet.membrane_force_n_m[0], n, 1.0);
        near(facet.membrane_force_n_m[1], n, 1.0);
        near(facet.membrane_force_n_m[2], 0.0, 1.0);
    }
    near(response.strain_energy_j, -0.6 * n * expansion, 1.0);
    for (force, equivalent) in response
        .material_nodal_force_6n
        .iter()
        .zip(prepared.load_6n())
    {
        near(*force, -*equivalent, 1.0);
    }
}

#[test]
fn g1_thickness_gradient_drives_compatible_curvature() {
    let mesh = mesh();
    let model = model(&mesh, section(), true);
    let curvature = 0.002;
    let strain = ShellThermalStrain::from_faces(
        [-0.05 * curvature, -0.05 * curvature, 0.0],
        [0.05 * curvature, 0.05 * curvature, 0.0],
        0.1,
    )
    .unwrap();
    let prepared = model
        .prepare_thermal(&[strain.into(); 2], || false)
        .unwrap();
    let response = prepared.solve_static(&[0.0; 24], || false).unwrap();
    for (node, p) in mesh.nodes.iter().enumerate() {
        near(
            response.displacement_6n[6 * node + 2],
            -0.5 * curvature * (p[0] * p[0] + p[1] * p[1]),
            0.001,
        );
        near(
            response.displacement_6n[6 * node + 3],
            -curvature * p[1],
            0.001,
        );
        near(
            response.displacement_6n[6 * node + 4],
            curvature * p[0],
            0.001,
        );
    }
    assert!(response.strain_energy_j < 1e-18);
}

#[test]
fn g1_coupled_section_preserves_signed_thermal_moments() {
    // Explicit independent section: A=100I, B=2I, D=3I is SPD.
    let matrix = |s| SheetMatrix::orthotropic(s, s, 0.0, s).unwrap();
    let section = ShearShellSection::new_coupled(
        matrix(100.0),
        [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
        matrix(3.0),
        [30.0, 0.0, 30.0],
        1.0,
        1.0,
        0.0,
        0.1,
        0.0,
    )
    .unwrap();
    let mesh = mesh();
    let model = model(&mesh, section, false);
    let strain = ShellThermalStrain::try_new([0.01, 0.02, 0.0], [0.03, -0.01, 0.0]).unwrap();
    let response = model
        .prepare_thermal(&[strain.into(); 2], || false)
        .unwrap()
        .recover(&[0.0; 24], || false)
        .unwrap();
    for facet in response.facets {
        near(facet.membrane_force_n_m[0], -1.06, 1.0);
        near(facet.membrane_force_n_m[1], -1.98, 1.0);
        near(facet.bending_moment_n[0], -0.11, 1.0);
        near(facet.bending_moment_n[1], -0.01, 1.0);
    }
}

#[test]
fn g3_spatial_gradient_retains_thermal_bubble_offset_and_energy_derivative() {
    let mesh = mesh();
    let model = model(&mesh, section(), false);
    let fields: Vec<_> = mesh
        .tris
        .iter()
        .map(|tri| {
            FacetThermalStrain::AtNodes(tri.map(|node| {
                let x = mesh.nodes[node][0];
                ShellThermalStrain::try_new([0.0; 3], [0.02 * x; 3]).unwrap()
            }))
        })
        .collect();
    let prepared = model.prepare_thermal(&fields, || false).unwrap();
    let zero = prepared.recover(&[0.0; 24], || false).unwrap();
    // Zero physical nodal motion still has an interior thermal director field.
    // Its whole-element observer must include the affine offset; an observer
    // reconstructed solely from the cold condensed map would return zero.
    let mut maximum_internal_strain = 0.0_f64;
    for element in 0..mesh.tris.len() {
        let field = prepared
            .strain_field(element, &[0.0; 18], || false)
            .unwrap();
        let bound = field
            .reduced_elastic_layer_bound([-0.05, 0.05], [[0.0; 3]; 2], &mut || false)
            .unwrap();
        maximum_internal_strain = maximum_internal_strain.max(bound);
    }
    assert!(maximum_internal_strain > 1e-8);
    // Condensation must shift the interior director response under thermal
    // forcing. A homogeneous-only strain expansion would return -f_n instead.
    for (force, equivalent) in zero.material_nodal_force_6n.iter().zip(prepared.load_6n()) {
        near(*force, -*equivalent, 1.0);
    }
    for dof in [0, 2, 3, 7, 10, 16, 21] {
        let mut plus = [0.0; 24];
        let mut minus = plus;
        plus[dof] = 1e-6;
        minus[dof] = -1e-6;
        let derivative = (prepared.recover(&plus, || false).unwrap().strain_energy_j
            - prepared.recover(&minus, || false).unwrap().strain_energy_j)
            / 2e-6;
        near(derivative, zero.material_nodal_force_6n[dof], 1.0);
    }
    assert!(zero.strain_energy_j > 0.0);
}

#[test]
fn g3_world_rotation_preserves_thermal_work_and_resultants() {
    let mesh = mesh();
    // Proper cyclic-axis rotation (x,y,z)->(z,x,y), exactly representable.
    let rotated = ShellMesh::new(
        mesh.nodes.iter().map(|p| [p[2], p[0], p[1]]).collect(),
        mesh.tris.clone(),
    )
    .unwrap();
    let a = model(&mesh, section(), false);
    let b = model(&rotated, section(), false);
    let strain =
        ShellThermalStrain::try_new([0.002, -0.001, 0.0005], [0.01, 0.02, -0.005]).unwrap();
    let ta = a.prepare_thermal(&[strain.into(); 2], || false).unwrap();
    let tb = b.prepare_thermal(&[strain.into(); 2], || false).unwrap();
    let ra = ta.recover(&[0.0; 24], || false).unwrap();
    let rb = tb.recover(&[0.0; 24], || false).unwrap();
    near(ra.strain_energy_j, rb.strain_energy_j, 1.0);
    for node in 0..4 {
        for block in [0, 3] {
            for c in 0..3 {
                near(
                    ta.load_6n()[6 * node + block + c],
                    tb.load_6n()[6 * node + block + (c + 1) % 3],
                    1.0,
                );
            }
        }
    }
}

#[test]
fn g0_bad_fields_and_cancellation_publish_no_partial_thermal_state() {
    let mesh = mesh();
    let model = model(&mesh, section(), false);
    let strain = ShellThermalStrain::try_new([0.0; 3], [0.0; 3]).unwrap();
    assert!(ShellThermalStrain::try_new([f64::NAN, 0.0, 0.0], [0.0; 3]).is_err());
    assert!(ShellThermalStrain::from_faces([0.0; 3], [0.0; 3], 0.0).is_err());
    assert!(matches!(
        model.prepare_thermal(&[], || false),
        Err(PlateError::SectionCount { .. })
    ));
    assert!(matches!(
        model.prepare_thermal(&[strain.into(); 2], || true),
        Err(PlateError::Cancelled)
    ));
    let prepared = model
        .prepare_thermal(&[strain.into(); 2], || false)
        .unwrap();
    assert!(prepared.recover(&[0.0; 3], || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        prepared.recover(&[0.0; 24], || {
            polls += 1;
            polls == 5
        }),
        Err(PlateError::Cancelled)
    ));
    assert_eq!(
        prepared
            .recover(&[0.0; 24], || false)
            .unwrap()
            .strain_energy_j,
        0.0
    );
}

#[test]
fn g1_temperature_dependent_source_paths_drive_the_actual_shell() {
    let mut claims = ClaimSet::new();
    for (name, scale) in [
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, 1.0),
        (ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[0], 1.0),
        (
            ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[1],
            -0.5,
        ),
        (ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[2], 3.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, INVERSE_TEMPERATURE_DIMS),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (300.0, scale * 1e-5),
                        (350.0, scale * 4e-5),
                        (400.0, scale * 2e-5),
                    ],
                    dims: INVERSE_TEMPERATURE_DIMS,
                },
                interpolation: InterpolationPolicy::LinearInside,
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic shell CTE integration oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "test shell".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let at = |t| QueryPoint::new().with("T", t).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let isotropic =
        integrate_isotropic_thermal_expansion(&card, &at(300.0), &at(400.0), selection.clone())
            .unwrap();
    // Two exact trapezoids; alpha(400)*100 would incorrectly produce 0.002.
    let expected = 50.0 * (1e-5 + 4e-5) / 2.0 + 50.0 * (4e-5 + 2e-5) / 2.0;
    let strain = ShellThermalStrain::from_isotropic(&isotropic).unwrap();
    near(strain.membrane()[0], expected, 0.001);
    let mesh = mesh();
    let model = model(&mesh, section(), true);
    let response = model
        .prepare_thermal(&[strain.into(); 2], || false)
        .unwrap()
        .solve_static(&[0.0; 24], || false)
        .unwrap();
    near(response.displacement_6n[6], expected, 0.001);
    near(response.displacement_6n[19], 0.6 * expected, 0.001);
    let orthotropic =
        integrate_orthotropic_thermal_expansion(&card, &at(300.0), &at(400.0), selection.clone())
            .unwrap();
    let oriented =
        ShellThermalStrain::from_orthotropic_in_plane(&orthotropic, core::f64::consts::FRAC_PI_4)
            .unwrap();
    near(oriented.membrane()[0], 0.25 * expected, 0.001);
    near(oriented.membrane()[1], 0.25 * expected, 0.001);
    near(oriented.membrane()[2], 1.5 * expected, 0.001);
    assert!(
        integrate_isotropic_thermal_expansion(&card, &at(300.0), &at(410.0), selection).is_err()
    );
}
