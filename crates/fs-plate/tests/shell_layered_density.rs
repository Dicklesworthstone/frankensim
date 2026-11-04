//! Layered reference inventory in the existing CST/DKT mass and pressure owners.
//! Independent volume moments test inertia; laminate stiffness remains explicit.
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use fs_material::sheet::SheetMatrix;
use fs_plate::PlateSection;
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace, ShellPressurePort};
use fs_plate::shell::profile::ProfileBudget;
use fs_plate::shell::reduction::rigid::mass_normalized_modes;
use fs_plate::shell::survey::{IsotropicMaterial, MeshShell};
use fs_plate::shell::{
    ReferenceDensityMoments, ShellMesh, ShellModel, ShellSection, ShellSupport,
    assemble_shell_materials_with_reference_density_moments,
    assemble_shell_sections_with_reference_density_moments,
    assemble_shell_with_reference_density_moments,
};

const RHO: f64 = 1000.0;
const H: f64 = 0.02;
const E: f64 = 3e9;

fn budget() -> ProfileBudget {
    ProfileBudget {
        max_nodes: 3,
        max_triangles: 1,
        max_feature_evaluations: 0,
    }
}
fn mesh(z: f64, reverse: bool) -> ShellMesh {
    ShellMesh::new(
        vec![[0.0, 0.0, z], [1.0, 0.0, z], [0.0, 1.0, z]],
        vec![if reverse { [0, 2, 1] } else { [0, 1, 2] }],
    )
    .unwrap()
}
fn section() -> PlateSection {
    PlateSection::isotropic(E, 0.0, H, 2.5 * RHO).unwrap()
}
fn laminate(scale: f64, reverse: bool, reference: ReferenceSurface) -> LaminateSection {
    let layer = |rho| LaminaLayer {
        thickness_m: H / 2.0,
        orientation_rad: 0.0,
        density_kg_m3: rho,
        material: None,
        law: LayerLaw::PlaneStress {
            stiffness_pa: [E, 0.0, 0.0, 0.0, E, 0.0, 0.0, 0.0, 0.5 * E],
            transverse_shear_pa: [0.5 * E; 2],
        },
    };
    let mut layers = vec![layer(scale * RHO), layer(scale * 4.0 * RHO)];
    if reverse {
        layers.reverse();
    }
    LaminateSection::try_new(layers, reference).unwrap()
}
fn moments(scale: f64, reverse: bool, reference: ReferenceSurface) -> ReferenceDensityMoments {
    ReferenceDensityMoments::from_laminate(&laminate(scale, reverse, reference)).unwrap()
}
fn graded(scale: f64) -> [[ReferenceDensityMoments; 3]; 1] {
    [[1.0, 4.0, 1.0].map(|value| moments(scale * value, false, ReferenceSurface::Midsurface))]
}
fn survey() -> MeshShell {
    MeshShell::new(
        mesh(0.0, false),
        vec![H, 2.0 * H, H],
        &[IsotropicMaterial {
            young_pa: E,
            poisson: 0.0,
            density_kg_m3: 2.5 * RHO,
        }],
        budget(),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1e-3),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}
fn rigid(mesh: &ShellMesh, axis: usize) -> Vec<f64> {
    let mut out = vec![0.0; 6 * mesh.nodes.len()];
    for (i, p) in mesh.nodes.iter().enumerate() {
        if axis < 3 {
            out[6 * i + axis] = 1.0;
        } else {
            let mut w = [0.0; 3];
            w[axis - 3] = 1.0;
            out[6 * i] = w[1] * p[2] - w[2] * p[1];
            out[6 * i + 1] = w[2] * p[0] - w[0] * p[2];
            out[6 * i + 2] = w[0] * p[1] - w[1] * p[0];
            out[6 * i + axis] = 1.0;
        }
    }
    out
}
fn bilinear(model: &ShellModel, a: &[f64], b: &[f64]) -> f64 {
    let mut mb = vec![0.0; model.free];
    model.m.spmv(b, &mut mb);
    a.iter().zip(mb).map(|(u, v)| u * v).sum()
}
// Standard rigid-body spatial inertia from independently integrated volume moments.
fn spatial(mass: f64, first: [f64; 3], second: [f64; 6]) -> [[f64; 6]; 6] {
    let [mx, my, mz] = first;
    let [xx, yy, zz, xy, xz, yz] = second;
    [
        [mass, 0.0, 0.0, 0.0, mz, -my],
        [0.0, mass, 0.0, -mz, 0.0, mx],
        [0.0, 0.0, mass, my, -mx, 0.0],
        [0.0, -mz, my, yy + zz, -xy, -xz],
        [mz, 0.0, -mx, -xy, xx + zz, -yz],
        [-my, mx, 0.0, -xz, -yz, xx + yy],
    ]
}
fn graded_spatial() -> [[f64; 6]; 6] {
    let (q0, q1, q2) = (2.5 * RHO, 3.0 * RHO / 8.0, 5.0 * RHO / 24.0);
    spatial(
        q0 * H * 17.0 / 12.0,
        [
            q0 * H * 13.0 / 20.0,
            q0 * H * 23.0 / 60.0,
            q1 * H.powi(2) * 31.0 / 15.0,
        ],
        [
            q0 * H * 23.0 / 60.0,
            q0 * H / 6.0,
            q2 * H.powi(3) * 31.0 / 10.0,
            q0 * H * 2.0 / 15.0,
            q1 * H.powi(2) * 31.0 / 30.0,
            q1 * H.powi(2) * 31.0 / 60.0,
        ],
    )
}
fn check_spatial(mesh: &ShellMesh, model: &ShellModel, expected: [[f64; 6]; 6]) {
    let basis: Vec<_> = (0..6).map(|axis| rigid(mesh, axis)).collect();
    for i in 0..6 {
        for j in 0..6 {
            close(bilinear(model, &basis[i], &basis[j]), expected[i][j]);
        }
    }
}

#[test]
fn measured_layered_shell_retains_all_spatial_moments_and_density_thickness_correlation() {
    let shell = survey();
    let field = graded(1.0);
    close(field[0][0].zeroth_kg_m3(), 2.5 * RHO);
    close(field[0][0].first_kg_m3(), 3.0 * RHO / 8.0);
    close(field[0][0].second_kg_m3(), 5.0 * RHO / 24.0);
    let model = shell
        .assemble_with_reference_density_moments(&field, &[], ShellSupport::Free)
        .unwrap();
    // Equal-thickness light/heavy layers, all q-fields proportional to 1+3x,
    // h=H(1+x), on the unit right triangle. These rational volume integrals
    // include the displaced z center and xz/yz products of inertia.
    let expected = graded_spatial();
    check_spatial(&shell.mesh, &model, expected);
    close(
        shell
            .mass_with_reference_density_moments_kg(&field)
            .unwrap(),
        expected[0][0],
    );
    assert!(expected[0][4] / expected[0][0] > 0.0);
    let mut translation = vec![0.0; model.free];
    let mut director = translation.clone();
    for (i, p) in shell.mesh.nodes.iter().enumerate() {
        translation[6 * i] = p[0];
        director[6 * i + 4] = p[0];
    }
    // Degree-five first moment and degree-six second moment, independently integrated.
    close(
        bilinear(&model, &translation, &director),
        (3.0 * RHO / 8.0) * H.powi(2) * 67.0 / 105.0,
    );
    close(
        bilinear(&model, &director, &director),
        (5.0 * RHO / 24.0) * H.powi(3) * 43.0 / 40.0,
    );
}

#[test]
fn layer_reversal_changes_cross_inertia_and_reversing_the_mesh_normal_preserves_the_body() {
    let forward = mesh(0.0, false);
    let backward = mesh(0.0, true);
    let original = moments(1.0, false, ReferenceSurface::Midsurface);
    let mirrored = moments(1.0, true, ReferenceSurface::Midsurface);
    close(original.zeroth_kg_m3(), mirrored.zeroth_kg_m3());
    close(original.first_kg_m3(), -mirrored.first_kg_m3());
    close(original.second_kg_m3(), mirrored.second_kg_m3());
    let build = |mesh: &ShellMesh, q| {
        assemble_shell_with_reference_density_moments(
            mesh,
            &section(),
            &[[q; 3]],
            &[],
            ShellSupport::Free,
        )
        .unwrap()
    };
    let a = build(&forward, original);
    let b = build(&forward, mirrored);
    let same_body = build(&backward, mirrored);
    assert!(a.m.get(0, 4) > 0.0);
    assert!(a.m.get(1, 3) < 0.0);
    close(a.m.get(0, 4), -b.m.get(0, 4));
    close(a.m.get(1, 3), -b.m.get(1, 3));
    for i in 0..a.free {
        close(a.m.get(i, i), b.m.get(i, i));
        for j in 0..a.free {
            close(a.m.get(i, j), same_body.m.get(i, j));
        }
    }
}

#[test]
fn translating_the_reference_surface_keeps_the_same_physical_inventory() {
    // Physical layers occupy z=[0,H/2] and [H/2,H], regardless of mesh plane.
    let i0 = 2.5 * RHO * H;
    let i1 = 13.0 * RHO * H.powi(2) / 8.0;
    let i2 = 29.0 * RHO * H.powi(3) / 24.0;
    let expected = spatial(
        i0 / 2.0,
        [i0 / 6.0, i0 / 6.0, i1 / 2.0],
        [
            i0 / 12.0,
            i0 / 12.0,
            i2 / 2.0,
            i0 / 24.0,
            i1 / 6.0,
            i1 / 6.0,
        ],
    );
    for offset in [0.0, H / 2.0, H] {
        let mesh = mesh(offset, false);
        let field = [[moments(1.0, false, ReferenceSurface::Offset(offset)); 3]];
        let model = assemble_shell_sections_with_reference_density_moments(
            &mesh,
            &[section()],
            &field,
            &[],
            ShellSupport::Free,
        )
        .unwrap();
        check_spatial(&mesh, &model, expected);
    }
}

#[test]
fn homogeneous_limits_reuse_the_original_scalar_density_and_independent_section_paths() {
    let shell = survey();
    let rho = [[RHO, 4.0 * RHO, RHO]];
    let field =
        rho.map(|row| row.map(|value| ReferenceDensityMoments::from_uniform(value).unwrap()));
    let old = shell
        .assemble_with_reference_density(&rho, &[], ShellSupport::Free)
        .unwrap();
    let new = shell
        .assemble_with_reference_density_moments(&field, &[], ShellSupport::Free)
        .unwrap();
    for i in 0..old.free {
        assert_eq!(old.m.row(i), new.m.row(i));
        assert_eq!(old.k.row(i), new.k.row(i));
    }
    assert_eq!(
        shell.mass_with_reference_density_kg(&rho).unwrap(),
        shell
            .mass_with_reference_density_moments_kg(&field)
            .unwrap()
    );

    let mesh = mesh(0.0, false);
    let layered = [[moments(1.0, false, ReferenceSurface::Midsurface); 3]];
    let ordinary = assemble_shell_with_reference_density_moments(
        &mesh,
        &section(),
        &layered,
        &[],
        ShellSupport::Free,
    )
    .unwrap();
    let independent = ShellSection::new(
        section(),
        SheetMatrix::orthotropic(E * H, E * H, 0.0, 0.5 * E * H).unwrap(),
    )
    .unwrap();
    let material = assemble_shell_materials_with_reference_density_moments(
        &mesh,
        &[independent],
        &layered,
        &[],
        ShellSupport::Free,
    )
    .unwrap();
    for i in 0..ordinary.free {
        assert_eq!(ordinary.m.row(i), material.m.row(i));
    }
}

#[test]
fn layered_density_reaches_native_elastic_modes_and_mass_normalized_pressure_ports() {
    let shell = survey();
    let light = graded(1.0);
    let heavy = graded(4.0);
    let assemble = |field: &[[ReferenceDensityMoments; 3]], nodes: &[usize], support| {
        shell
            .assemble_with_reference_density_moments(field, nodes, support)
            .unwrap()
    };
    let a = assemble(&light, &[0, 1], ShellSupport::Clamped);
    let b = assemble(&heavy, &[0, 1], ShellSupport::Clamped);
    let modes = |model: &ShellModel| {
        let mut k = vec![0.0; model.free * model.free];
        let mut m = k.clone();
        for i in 0..model.free {
            for j in 0..model.free {
                k[i * model.free + j] = model.k.get(i, j);
                m[i * model.free + j] = model.m.get(i, j);
            }
        }
        fs_modal::eigh_gen_dense(&k, &m, model.free).unwrap()
    };
    for i in 0..a.free {
        assert_eq!(a.k.row(i), b.k.row(i));
        for j in 0..a.free {
            close(b.m.get(i, j), 4.0 * a.m.get(i, j));
        }
    }
    for (a, b) in modes(&a).iter().zip(modes(&b)) {
        assert!(a.lambda > 0.0 && b.lambda > 0.0);
        assert!((b.lambda / a.lambda - 0.25).abs() < 2e-9);
    }
    let port = ShellPressurePort::new(
        &shell.mesh,
        vec![OrientedShellFacet {
            triangle: 0,
            side: ShellFace::Positive,
        }],
        budget(),
        || false,
    )
    .unwrap();
    let motion = [0.3, 0.0, 1.0, 0.0, 2.0, 0.0];
    let participation = |field: &[[ReferenceDensityMoments; 3]]| {
        let model = assemble(field, &[], ShellSupport::Free);
        let mode = mass_normalized_modes(&shell.mesh, &model, &[motion]).unwrap();
        port.project_modes(&[mode[0].phi.clone()], 100, || false)
            .unwrap()[0]
    };
    let mut norm_squared = 0.0;
    let spatial = graded_spatial();
    for i in 0..6 {
        for j in 0..6 {
            norm_squared += motion[i] * spatial[i][j] * motion[j];
        }
    }
    // Integral of normal rigid velocity 1-2x over the unit triangle is 1/6.
    close(participation(&light), (1.0 / 6.0) / norm_squared.sqrt());
    close(participation(&heavy), 0.5 * participation(&light));
}

#[test]
fn impossible_mass_moments_and_incomplete_fields_refuse_without_changing_the_source() {
    for [q0, q1, q2] in [
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [2.0, 1.0, 0.5],
        [1.0, -2.0, 1.0],
        [f64::NAN, 0.0, 1.0],
        [1.0, f64::NAN, 1.0],
        [1.0, 0.0, f64::INFINITY],
    ] {
        assert!(ReferenceDensityMoments::new(q0, q1, q2).is_err());
    }
    assert!(ReferenceDensityMoments::new(2.0, -0.5, 1.0).is_ok());
    assert!(ReferenceDensityMoments::new(1e300, 0.0, 1e-300).is_ok());
    let shell = survey();
    for field in [vec![], vec![graded(1.0)[0]; 2]] {
        assert!(
            shell
                .assemble_with_reference_density_moments(&field, &[], ShellSupport::Free)
                .is_err()
        );
        assert!(
            shell
                .mass_with_reference_density_moments_kg(&field)
                .is_err()
        );
    }
    let huge = ReferenceDensityMoments::new(f64::MAX, 0.0, f64::MAX).unwrap();
    let mut large = shell.clone();
    for p in &mut large.mesh.nodes {
        p[0] *= 1000.0;
        p[1] *= 1000.0;
    }
    assert!(
        large
            .assemble_with_reference_density_moments(&[[huge; 3]], &[], ShellSupport::Free)
            .is_err()
    );
    assert!(
        large
            .mass_with_reference_density_moments_kg(&[[huge; 3]])
            .is_err()
    );
    assert_eq!(shell.nodal_thickness_m, [H, 2.0 * H, H]);
    assert_eq!(shell.sections[0].density, 2.5 * RHO);
}
