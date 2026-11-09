//! G1 bimetal bending and G3 work/objectivity in the original native shell.
//! Authored numerical material fixtures; no measured laminate qualification.

use fs_material::laminate::thermal::nodal::PreparedNodalLaminateThermal;
use fs_material::laminate::thermal::{LaminaThermalStrain, PreparedLaminateThermal};
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::prestress::MembranePrestressLimits;
use fs_plate::shell::mitc::thermal::{
    FacetThermalStrain, LaminateThermalField, NodalLaminateThermalField, ShellThermalStrain,
};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use std::sync::Arc;

const E: f64 = 2.0e8;
const H: f64 = 0.02;
const AREA: f64 = 0.6;

fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 2.0e-8 * scale.max(a.abs()).max(b.abs()),
        "{a:e} != {b:e}"
    );
}

fn layer(e: f64, h: f64, rho: f64) -> LaminaLayer {
    LaminaLayer {
        thickness_m: h,
        orientation_rad: 0.0,
        density_kg_m3: rho,
        material: None,
        law: LayerLaw::PlaneStress {
            stiffness_pa: [e, 0.0, 0.0, 0.0, e, 0.0, 0.0, 0.0, 0.5 * e],
            transverse_shear_pa: [0.5 * e; 2],
        },
    }
}

fn source() -> Arc<LaminateSection> {
    Arc::new(
        LaminateSection::try_new(
            vec![layer(E, 0.5 * H, 900.0), layer(E, 0.5 * H, 1200.0)],
            ReferenceSurface::Midsurface,
        )
        .unwrap(),
    )
}

fn thermal(source: Arc<LaminateSection>, expansion: f64, angle: f64) -> LaminateThermalField {
    let profiles = [
        LaminaThermalStrain::uniform([0.0; 3]).unwrap(),
        LaminaThermalStrain::uniform([expansion, expansion, 0.0]).unwrap(),
    ];
    let state = Arc::new(PreparedLaminateThermal::try_new(source, &profiles, || false).unwrap());
    LaminateThermalField::try_new(state, angle).unwrap()
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

fn model(mesh: &ShellMesh, source: &LaminateSection, angle: f64, fixed: bool) -> MitcShellModel {
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
    let section = ShearShellSection::from_laminate(source, angle).unwrap();
    assemble_mitc3p(
        mesh,
        &directors,
        &vec![section; mesh.tris.len()],
        &constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: mesh.nodes.len(),
            max_triangles: mesh.tris.len(),
        },
        1.0e-3,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_uniform_heating_bends_bimetal_and_retains_irreducible_layer_energy() {
    let mesh = mesh();
    let source = source();
    let model = model(&mesh, &source, 0.0, true);
    let mass = model.pencil().m.clone();
    for expansion in [0.0002, -0.0002] {
        let field = thermal(source.clone(), expansion, 0.0);
        let prepared = model
            .prepare_laminate_thermal(&[field.clone(), field], || false)
            .unwrap();
        let response = prepared.solve_static(&[0.0; 24], || false).unwrap();
        // Independent equal-stiffness bilayer solution: epsilon0=delta/2,
        // kappa=3 delta/(2h). Piecewise constant free strains cannot be made
        // stress-free by an affine compatible deformation: U/A=E h delta²/16
        // for the two uncoupled in-plane normal directions.
        let strain = 0.5 * expansion;
        let curvature = 1.5 * expansion / H;
        for (node, p) in mesh.nodes.iter().enumerate() {
            near(response.displacement_6n[6 * node], strain * p[0], 0.001);
            near(response.displacement_6n[6 * node + 1], strain * p[1], 0.001);
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
        near(
            response.strain_energy_j,
            AREA * E * H * expansion * expansion / 16.0,
            1.0e-8,
        );
        assert!(response.strain_energy_j > 0.0);
        for facet in response.facets {
            for value in facet
                .membrane_force_n_m
                .into_iter()
                .chain(facet.bending_moment_n)
            {
                near(value, 0.0, 1.0);
            }
        }
    }
    for row in 0..model.pencil().free {
        assert_eq!(mass.row(row), model.pencil().m.row(row));
    }
}

#[test]
fn g1_restrained_layers_preserve_thermal_force_moment_and_storage() {
    let source = source();
    let model = model(&mesh(), &source, 0.0, false);
    let expansion = 0.0002;
    let field = thermal(source, expansion, 0.0);
    let prepared = model
        .prepare_laminate_thermal(&[field.clone(), field], || false)
        .unwrap();
    let response = prepared.recover(&[0.0; 24], || false).unwrap();
    for facet in response.facets {
        for i in 0..2 {
            near(facet.membrane_force_n_m[i], -E * H * expansion / 2.0, 1.0);
            near(facet.bending_moment_n[i], -E * H * H * expansion / 8.0, 1.0);
        }
    }
    near(
        response.strain_energy_j,
        AREA * E * H * expansion * expansion / 2.0,
        1.0e-8,
    );
    for (force, load) in response
        .material_nodal_force_6n
        .iter()
        .zip(prepared.load_6n())
    {
        near(*force, -*load, 1.0);
    }
}

#[test]
fn g3_layered_energy_derivative_matches_native_material_forces() {
    let source = Arc::new(
        LaminateSection::try_new(
            vec![layer(E, H * 0.3, 800.0), layer(E * 3.0, H * 0.7, 1600.0)],
            ReferenceSurface::Offset(H * 0.4),
        )
        .unwrap(),
    );
    let angle = 0.37;
    let model = model(&mesh(), &source, angle, false);
    let profiles = [
        LaminaThermalStrain::from_faces([0.0001, -0.0002, 0.0003], [0.0004, 0.0001, -0.0002])
            .unwrap(),
        LaminaThermalStrain::uniform([-0.0003, 0.0005, 0.0001]).unwrap(),
    ];
    let state = Arc::new(PreparedLaminateThermal::try_new(source, &profiles, || false).unwrap());
    let field = LaminateThermalField::try_new(state, angle).unwrap();
    let prepared = model
        .prepare_laminate_thermal(&[field.clone(), field], || false)
        .unwrap();
    let u = std::array::from_fn::<_, 24, _>(|i| 0.00001 * ((i % 7) as f64 - 3.0));
    let response = prepared.recover(&u, || false).unwrap();
    for dof in [0, 2, 3, 7, 10, 16, 21] {
        let mut plus = u;
        let mut minus = u;
        plus[dof] += 1.0e-6;
        minus[dof] -= 1.0e-6;
        let derivative = (prepared.recover(&plus, || false).unwrap().strain_energy_j
            - prepared.recover(&minus, || false).unwrap().strain_energy_j)
            / 2.0e-6;
        near(derivative, response.material_nodal_force_6n[dof], 1.0);
    }
}

#[test]
fn g3_layered_loads_and_energy_follow_rigid_world_rotation() {
    let mesh = mesh();
    let rotated = ShellMesh::new(
        mesh.nodes.iter().map(|p| [p[2], p[0], p[1]]).collect(),
        mesh.tris.clone(),
    )
    .unwrap();
    let source = source();
    let field = thermal(source.clone(), 0.0002, 0.31);
    let fields = [field.clone(), field];
    let a = model(&mesh, &source, 0.31, false);
    let b = model(&rotated, &source, 0.31, false);
    let pa = a.prepare_laminate_thermal(&fields, || false).unwrap();
    let pb = b.prepare_laminate_thermal(&fields, || false).unwrap();
    near(
        pa.recover(&[0.0; 24], || false).unwrap().strain_energy_j,
        pb.recover(&[0.0; 24], || false).unwrap().strain_energy_j,
        1.0e-8,
    );
    for node in 0..4 {
        for block in [0, 3] {
            for c in 0..3 {
                near(
                    pa.load_6n()[6 * node + block + c],
                    pb.load_6n()[6 * node + block + (c + 1) % 3],
                    1.0,
                );
            }
        }
    }
}

#[test]
fn g3_one_layer_reduces_to_original_affine_thermal_shell() {
    let source = Arc::new(
        LaminateSection::try_new(vec![layer(E, H, 1000.0)], ReferenceSurface::Midsurface).unwrap(),
    );
    let bottom = [-0.0001, 0.0002, 0.0001];
    let top = [0.0003, -0.0001, -0.0002];
    let model = model(&mesh(), &source, 0.0, false);
    let state = Arc::new(
        PreparedLaminateThermal::try_new(
            source,
            &[LaminaThermalStrain::from_faces(bottom, top).unwrap()],
            || false,
        )
        .unwrap(),
    );
    let field = LaminateThermalField::try_new(state, 0.0).unwrap();
    let layered = model
        .prepare_laminate_thermal(&[field.clone(), field], || false)
        .unwrap();
    let affine = ShellThermalStrain::from_faces(bottom, top, H).unwrap();
    let affine = model
        .prepare_thermal(&[affine.into(); 2], || false)
        .unwrap();
    for (a, b) in layered.load_6n().iter().zip(affine.load_6n()) {
        near(*a, *b, 1.0);
    }
    let u = std::array::from_fn::<_, 24, _>(|i| 0.00001 * ((i % 5) as f64 - 2.0));
    let a = layered.recover(&u, || false).unwrap();
    let b = affine.recover(&u, || false).unwrap();
    near(a.strain_energy_j, b.strain_energy_j, 1.0e-8);
    for (a, b) in a
        .material_nodal_force_6n
        .iter()
        .zip(&b.material_nodal_force_6n)
    {
        near(*a, *b, 1.0);
    }
}

#[test]
fn g1_oriented_layer_face_stress_uses_the_work_conjugate_tensor_transform() {
    let mut ply = layer(E, H, 1000.0);
    ply.law = LayerLaw::PlaneStress {
        stiffness_pa: [90.0e6, 0.0, 0.0, 0.0, 30.0e6, 0.0, 0.0, 0.0, 12.0e6],
        transverse_shear_pa: [12.0e6, 15.0e6],
    };
    let source =
        Arc::new(LaminateSection::try_new(vec![ply], ReferenceSurface::Midsurface).unwrap());
    let state = Arc::new(
        PreparedLaminateThermal::try_new(
            source.clone(),
            &[LaminaThermalStrain::uniform([0.002, -0.001, 0.0008]).unwrap()],
            || false,
        )
        .unwrap(),
    );
    let field = LaminateThermalField::try_new(state, core::f64::consts::FRAC_PI_4).unwrap();
    // Principal restraint stress [-180000,30000,-9600] Pa rotated by 45°.
    // Actual stress shear is NOT doubled like engineering strain shear.
    let expected = [-65400.0, -84600.0, -105000.0];
    let stress = field.layer_stress(0, [0.0; 3], [0.0; 3]).unwrap();
    for i in 0..3 {
        near(stress.bottom_pa[i], expected[i], 1.0);
        near(stress.top_pa[i], expected[i], 1.0);
        near(field.membrane_force_n_m()[i], -H * expected[i], 1.0);
    }
    let wrong_angle = model(&mesh(), &source, 0.0, false);
    assert!(
        wrong_angle
            .prepare_laminate_thermal(&[field.clone(), field], || false)
            .is_err()
    );
}

#[test]
fn g0_foreign_section_and_cancellation_refuse_without_mutating_source() {
    let source = source();
    let field = thermal(source.clone(), 0.0002, 0.0);
    let model = model(&mesh(), &source, 0.0, false);
    assert!(matches!(
        model.prepare_laminate_thermal(&[], || false),
        Err(PlateError::SectionCount { .. })
    ));
    let other = Arc::new(
        LaminateSection::try_new(
            vec![layer(2.0 * E, H / 2.0, 900.0), layer(E, H / 2.0, 1200.0)],
            ReferenceSurface::Midsurface,
        )
        .unwrap(),
    );
    let foreign = thermal(other, 0.0002, 0.0);
    assert!(
        model
            .prepare_laminate_thermal(&[foreign.clone(), foreign], || false)
            .is_err()
    );
    let fields = [field.clone(), field];
    assert!(matches!(
        model.prepare_laminate_thermal(&fields, || true),
        Err(PlateError::Cancelled)
    ));
    let prepared = model.prepare_laminate_thermal(&fields, || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        prepared.recover(&[0.0; 24], || {
            polls += 1;
            polls == 4
        }),
        Err(PlateError::Cancelled)
    ));
    assert!(
        prepared
            .recover(&[0.0; 24], || false)
            .unwrap()
            .strain_energy_j
            > 0.0
    );
}

fn nodal_fields(
    mesh: &ShellMesh,
    source: &Arc<LaminateSection>,
    profiles: &[Vec<LaminaThermalStrain>],
    angle: f64,
) -> Vec<NodalLaminateThermalField> {
    let nodes: Vec<_> = profiles
        .iter()
        .map(|profile| {
            Arc::new(PreparedLaminateThermal::try_new(source.clone(), profile, || false).unwrap())
        })
        .collect();
    mesh.tris
        .iter()
        .map(|tri| {
            let thermal =
                PreparedNodalLaminateThermal::try_new(tri.map(|node| nodes[node].clone())).unwrap();
            NodalLaminateThermalField::try_new(Arc::new(thermal), angle).unwrap()
        })
        .collect()
}

#[test]
fn g3_nodal_single_layer_reuses_original_p1_thermal_load_bubble_and_storage() {
    let mesh = mesh();
    let source = Arc::new(
        LaminateSection::try_new(vec![layer(E, H, 1000.0)], ReferenceSurface::Midsurface).unwrap(),
    );
    let profiles: Vec<_> = mesh
        .nodes
        .iter()
        .map(|p| {
            vec![
                LaminaThermalStrain::from_faces(
                    [0.0001 + 0.0002 * p[0], -0.0002 * p[1], 0.0003 * p[0]],
                    [-0.0002 * p[0], 0.0003 + 0.0001 * p[1], -0.0001 * p[1]],
                )
                .unwrap(),
            ]
        })
        .collect();
    let fields = nodal_fields(&mesh, &source, &profiles, 0.0);
    let model = model(&mesh, &source, 0.0, false);
    let layered = model
        .prepare_nodal_laminate_thermal(&fields, || false)
        .unwrap();
    let affine: Vec<_> = mesh
        .tris
        .iter()
        .map(|tri| {
            FacetThermalStrain::AtNodes(tri.map(|node| {
                let [bottom, top] = profiles[node][0].faces();
                ShellThermalStrain::from_faces(bottom, top, H).unwrap()
            }))
        })
        .collect();
    let original = model.prepare_thermal(&affine, || false).unwrap();
    for (a, b) in layered.load_6n().iter().zip(original.load_6n()) {
        near(*a, *b, 1.0);
    }
    let u = std::array::from_fn::<_, 24, _>(|i| 0.00001 * ((i % 7) as f64 - 3.0));
    let actual = layered.recover(&u, || false).unwrap();
    let expected = original.recover(&u, || false).unwrap();
    near(actual.strain_energy_j, expected.strain_energy_j, 1e-8);
    for (a, b) in actual
        .material_nodal_force_6n
        .iter()
        .zip(&expected.material_nodal_force_6n)
    {
        near(*a, *b, 1.0);
    }
    // The energy derivative includes the SAME inhomogeneous director bubble.
    for dof in [2, 4, 9, 16] {
        let mut plus = u;
        let mut minus = u;
        plus[dof] += 1e-6;
        minus[dof] -= 1e-6;
        let derivative = (layered.recover(&plus, || false).unwrap().strain_energy_j
            - layered.recover(&minus, || false).unwrap().strain_energy_j)
            / 2e-6;
        near(derivative, actual.material_nodal_force_6n[dof], 1.0);
    }
    assert!(matches!(
        model.prepare_nodal_laminate_thermal(&[], || false),
        Err(PlateError::SectionCount { .. })
    ));
    assert!(matches!(
        model.prepare_nodal_laminate_thermal(&fields, || true),
        Err(PlateError::Cancelled)
    ));
}

#[test]
fn g1_nodal_layer_heating_preserves_spatial_storage_and_actual_constrained_prestress() {
    let mesh = mesh();
    let source = Arc::new(
        LaminateSection::try_new(
            vec![
                layer(E, 0.25 * H, 900.0),
                layer(E, 0.5 * H, 1200.0),
                layer(E, 0.25 * H, 1500.0),
            ],
            ReferenceSurface::Midsurface,
        )
        .unwrap(),
    );
    let expansion: Vec<_> = mesh
        .nodes
        .iter()
        .map(|p| 0.0001 + 0.0004 * p[0] + 0.0003 * p[1])
        .collect();
    let profiles: Vec<_> = expansion
        .iter()
        .map(|&value| {
            let hot = LaminaThermalStrain::uniform([value, value, 0.0]).unwrap();
            vec![hot, LaminaThermalStrain::uniform([0.0; 3]).unwrap(), hot]
        })
        .collect();
    let fields = nodal_fields(&mesh, &source, &profiles, 0.0);
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 4], 0.9).unwrap();
    // Restrain every in-plane motion/director; retain transverse coordinates
    // with one translation anchor. Symmetric heating has no thermal moment.
    let constraints: Vec<_> = (0..4)
        .map(|node| ShellConstraint {
            node,
            mask: [true, true, node == 0, true, true, true],
            frame: ConstraintFrame::Global,
        })
        .collect();
    let model = assemble_mitc3p(
        &mesh,
        &directors,
        &vec![ShearShellSection::from_laminate(&source, 0.0).unwrap(); 2],
        &constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 4,
            max_triangles: 2,
        },
        1e-3,
        || false,
    )
    .unwrap();
    let mass = model.pencil().m.clone();
    let prepared = model
        .prepare_nodal_laminate_thermal(&fields, || false)
        .unwrap();
    let prestressed = prepared
        .solve_membrane_prestress(
            &[0.0; 24],
            MembranePrestressLimits {
                maximum_displacement_gradient: 0.01,
                maximum_rotation_rad: 0.01,
            },
            || false,
        )
        .unwrap();
    let response = prestressed.reference_state();
    let mut energy = 0.0;
    let mut centroid_energy = 0.0;
    let mut geometric_work = 0.0;
    for (e, tri) in mesh.tris.iter().enumerate() {
        let delta = tri.map(|node| expansion[node]);
        let mean = delta.iter().sum::<f64>() / 3.0;
        // Exact simplex moments: average(lambda_i²)=1/6,
        // average(lambda_i lambda_j)=1/12 for i != j.
        let mean_square = (delta.iter().map(|v| v * v).sum::<f64>()
            + delta[0] * delta[1]
            + delta[0] * delta[2]
            + delta[1] * delta[2])
            / 6.0;
        let area = mesh.facet(e).unwrap().area_m2;
        let force = -0.5 * E * H * mean;
        for c in 0..2 {
            near(response.facets[e].membrane_force_n_m[c], force, 1.0);
            near(response.facets[e].bending_moment_n[c], 0.0, 1e-3);
        }
        energy += area * E * H * mean_square / 2.0;
        centroid_energy += area * E * H * mean * mean / 2.0;
        geometric_work += area * force * 5.0; // w=x+2y.
        let lambda = [0.2, 0.3, 0.5];
        let at = (0..3).map(|i| lambda[i] * delta[i]).sum::<f64>();
        let outer = fields[e]
            .layer_stress(0, lambda, [0.0; 3], [0.0; 3])
            .unwrap();
        let inner = fields[e]
            .layer_stress(1, lambda, [0.0; 3], [0.0; 3])
            .unwrap();
        near(outer.bottom_pa[0], -E * at, 1.0);
        near(inner.bottom_pa[0], 0.0, 1.0);
    }
    near(response.strain_energy_j, energy, 1e-8);
    assert!(
        energy > 1.05 * centroid_energy,
        "facet averaging must lose spatial storage"
    );
    let mut probe = [0.0; 24];
    for (node, p) in mesh.nodes.iter().enumerate() {
        probe[6 * node + 2] = p[0] + 2.0 * p[1];
        for c in 0..6 {
            near(response.displacement_6n[6 * node + c], 0.0, 1e-5);
        }
    }
    let reduced = prestressed.project_load(&probe).unwrap();
    let observed: f64 = (0..reduced.len())
        .map(|row| {
            let (columns, values) = prestressed.geometric_stiffness().row(row);
            reduced[row]
                * columns
                    .iter()
                    .zip(values)
                    .map(|(&col, &v)| v * reduced[col])
                    .sum::<f64>()
        })
        .sum();
    near(observed, geometric_work, 1.0);
    assert!(
        observed < 0.0,
        "thermal compression must soften the signed pencil"
    );
    for row in 0..model.pencil().free {
        assert_eq!(mass.row(row), prestressed.pencil().m.row(row));
    }
}

#[test]
fn g3_nodal_layer_loading_follows_world_rotation_and_rejects_a_foreign_frame() {
    let mesh = mesh();
    let rotated = ShellMesh::new(
        mesh.nodes.iter().map(|p| [p[2], p[0], p[1]]).collect(),
        mesh.tris.clone(),
    )
    .unwrap();
    let mut ply = layer(E, H, 1000.0);
    ply.law = LayerLaw::PlaneStress {
        stiffness_pa: [90e6, 7e6, 0.0, 7e6, 30e6, 0.0, 0.0, 0.0, 12e6],
        transverse_shear_pa: [12e6, 15e6],
    };
    let source =
        Arc::new(LaminateSection::try_new(vec![ply], ReferenceSurface::Midsurface).unwrap());
    let profiles: Vec<_> = (0..4)
        .map(|node| {
            vec![
                LaminaThermalStrain::uniform([
                    0.0001 * node as f64,
                    -0.0002,
                    0.00003 * (node + 1) as f64,
                ])
                .unwrap(),
            ]
        })
        .collect();
    let fields = nodal_fields(&mesh, &source, &profiles, 0.31);
    let a = model(&mesh, &source, 0.31, false);
    let b = model(&rotated, &source, 0.31, false);
    let pa = a.prepare_nodal_laminate_thermal(&fields, || false).unwrap();
    let pb = b.prepare_nodal_laminate_thermal(&fields, || false).unwrap();
    near(
        pa.recover(&[0.0; 24], || false).unwrap().strain_energy_j,
        pb.recover(&[0.0; 24], || false).unwrap().strain_energy_j,
        1e-8,
    );
    for node in 0..4 {
        for block in [0, 3] {
            for c in 0..3 {
                near(
                    pa.load_6n()[6 * node + block + c],
                    pb.load_6n()[6 * node + block + (c + 1) % 3],
                    1.0,
                );
            }
        }
    }
    let foreign = model(&mesh, &source, 0.0, false);
    assert!(
        foreign
            .prepare_nodal_laminate_thermal(&fields, || false)
            .is_err()
    );
}
