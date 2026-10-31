//! G1 bimetal bending and G3 work/objectivity in the original native shell.
//! Authored numerical material fixtures; no measured laminate qualification.

use fs_material::laminate::thermal::{LaminaThermalStrain, PreparedLaminateThermal};
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::{LaminateThermalField, ShellThermalStrain};
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
