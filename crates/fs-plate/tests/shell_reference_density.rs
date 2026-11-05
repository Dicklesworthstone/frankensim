//! G1/G3: exact graded reference inventory reaches the native shell pencil,
//! mass-normalized modes and reciprocal acoustic pressure port.
use fs_material::sheet::SheetMatrix;
use fs_plate::PlateSection;
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace, ShellPressurePort};
use fs_plate::shell::profile::ProfileBudget;
use fs_plate::shell::reduction::rigid::mass_normalized_modes;
use fs_plate::shell::survey::{IsotropicMaterial, MeshShell};
use fs_plate::shell::{
    CstDktMass, ShellMesh, ShellModel, ShellSection, ShellSupport,
    assemble_shell_materials_with_reference_density, assemble_shell_sections_with_mass,
    assemble_shell_sections_with_reference_density, assemble_shell_with_mass,
    assemble_shell_with_reference_density,
};

const RHO: f64 = 1000.0;
const H: f64 = 0.02;

fn budget() -> ProfileBudget {
    ProfileBudget {
        max_nodes: 4,
        max_triangles: 2,
        max_feature_evaluations: 0,
    }
}
fn triangle() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}
fn section() -> PlateSection {
    PlateSection::isotropic(3e9, 0.25, H, RHO).unwrap()
}
fn survey() -> MeshShell {
    MeshShell::new(
        triangle(),
        vec![H, 2.0 * H, H],
        &[IsotropicMaterial {
            young_pa: 3e9,
            poisson: 0.25,
            density_kg_m3: RHO,
        }],
        budget(),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-3),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}
fn bilinear(model: &ShellModel, left: &[f64], right: &[f64]) -> f64 {
    let mut action = vec![0.0; model.free];
    model.m.spmv(right, &mut action);
    left.iter().zip(action).map(|(a, b)| a * b).sum()
}
fn rigid(mesh: &ShellMesh, axis: usize) -> Vec<f64> {
    let mut field = vec![0.0; 6 * mesh.nodes.len()];
    for (i, p) in mesh.nodes.iter().enumerate() {
        if axis < 3 {
            field[6 * i + axis] = 1.0;
        } else {
            let mut w = [0.0; 3];
            w[axis - 3] = 1.0;
            field[6 * i] = w[1] * p[2] - w[2] * p[1];
            field[6 * i + 1] = w[2] * p[0] - w[0] * p[2];
            field[6 * i + 2] = w[0] * p[1] - w[1] * p[0];
            field[6 * i + axis] = 1.0;
        }
    }
    field
}
fn identical(a: &ShellModel, b: &ShellModel) {
    assert_eq!(a.dof_map, b.dof_map);
    assert_eq!(a.free, b.free);
    for row in 0..a.free {
        assert_eq!(a.k.row(row), b.k.row(row));
        assert_eq!(a.m.row(row), b.m.row(row));
    }
}

#[test]
fn correlated_density_and_thickness_retain_the_complete_rigid_spatial_inertia() {
    let shell = survey();
    let rho = [[RHO, 4.0 * RHO, RHO]];
    let model = shell
        .assemble_with_reference_density(&rho, &[], ShellSupport::Free)
        .unwrap();
    // Unit right triangle, rho=r(1+3x), h=t(1+x), symmetric through thickness.
    // Direct polynomial volume integrals; no production quadrature is used.
    let mass = RHO * H * 17.0 / 12.0;
    let mx = RHO * H * 13.0 / 20.0;
    let my = RHO * H * 23.0 / 60.0;
    let x2 = RHO * H * 23.0 / 60.0;
    let y2 = RHO * H / 6.0;
    let xy = RHO * H * 2.0 / 15.0;
    let z2 = RHO * H.powi(3) * 31.0 / 120.0;
    let expected = [
        [mass, 0.0, 0.0, 0.0, 0.0, -my],
        [0.0, mass, 0.0, 0.0, 0.0, mx],
        [0.0, 0.0, mass, my, -mx, 0.0],
        [0.0, 0.0, my, y2 + z2, -xy, 0.0],
        [0.0, 0.0, -mx, -xy, x2 + z2, 0.0],
        [-my, mx, 0.0, 0.0, 0.0, x2 + y2],
    ];
    let motions: Vec<_> = (0..6).map(|i| rigid(&shell.mesh, i)).collect();
    for i in 0..6 {
        for j in 0..6 {
            close(bilinear(&model, &motions[i], &motions[j]), expected[i][j]);
        }
    }
    close(shell.mass_with_reference_density_kg(&rho).unwrap(), mass);
    close(mx / mass, 39.0 / 85.0);
    close(my / mass, 23.0 / 85.0);
    let product_of_means = (2.0 * RHO) * (4.0 * H / 3.0) * 0.5;
    assert!(
        (mass - product_of_means).abs() > 1.0,
        "mean fields must not erase correlation"
    );
    // An affine director rate probes degree SIX, beyond the former seven-point
    // thickness rule: integral rho*h^3*x^2/12 = r*t^3*43/480.
    let mut director = vec![0.0; model.free];
    for (node, p) in shell.mesh.nodes.iter().enumerate() {
        director[6 * node + 4] = p[0];
    }
    close(
        bilinear(&model, &director, &director),
        RHO * H.powi(3) * 43.0 / 480.0,
    );
}

#[test]
fn equal_mass_redistribution_changes_center_and_rotary_inertia() {
    let mesh = triangle();
    let build = |rho| {
        assemble_shell_with_reference_density(&mesh, &section(), &[rho], &[], ShellSupport::Free)
            .unwrap()
    };
    let right = build([RHO, 4.0 * RHO, RHO]);
    let left = build([3.0 * RHO, 0.5 * RHO, 2.5 * RHO]);
    let translation = rigid(&mesh, 0);
    let spin = rigid(&mesh, 5);
    close(bilinear(&right, &translation, &translation), RHO * H);
    close(bilinear(&left, &translation, &translation), RHO * H);
    assert!(
        (bilinear(&right, &translation, &spin) - bilinear(&left, &translation, &spin)).abs() > 0.1
    );
    assert!((bilinear(&right, &spin, &spin) - bilinear(&left, &spin, &spin)).abs() > 0.1);
}

#[test]
fn uniform_limits_and_sharp_facet_jumps_retain_existing_matrix_bits() {
    let shell = survey();
    let ordinary = shell.assemble(&[], ShellSupport::Free).unwrap();
    let explicit = shell
        .assemble_with_reference_density(&[[RHO; 3]], &[], ShellSupport::Free)
        .unwrap();
    identical(&ordinary, &explicit);
    let mesh = triangle();
    identical(
        &assemble_shell_with_mass(
            &mesh,
            &section(),
            &[],
            ShellSupport::Free,
            CstDktMass::Consistent,
        )
        .unwrap(),
        &assemble_shell_with_reference_density(
            &mesh,
            &section(),
            &[[RHO; 3]],
            &[],
            ShellSupport::Free,
        )
        .unwrap(),
    );

    let mesh = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [1, 3, 2]],
    )
    .unwrap();
    let mut dense = section();
    dense.density = 3.0 * RHO;
    let original = assemble_shell_sections_with_mass(
        &mesh,
        &[section(), dense],
        &[0],
        ShellSupport::Clamped,
        CstDktMass::Consistent,
    )
    .unwrap();
    let rho = [[RHO; 3], [3.0 * RHO; 3]]; // Different densities at shared vertices remain sharp.
    let explicit = assemble_shell_sections_with_reference_density(
        &mesh,
        &[section(); 2],
        &rho,
        &[0],
        ShellSupport::Clamped,
    )
    .unwrap();
    identical(&original, &explicit);

    let plate = section();
    let membrane = plate.d.map(|d| 12.0 * d / (H * H));
    let independent = ShellSection::new(
        plate,
        SheetMatrix::orthotropic(membrane[0], membrane[4], membrane[1], membrane[8]).unwrap(),
    )
    .unwrap();
    let model = assemble_shell_materials_with_reference_density(
        &mesh,
        &[independent; 2],
        &rho,
        &[0],
        ShellSupport::Clamped,
    )
    .unwrap();
    for i in 0..model.free {
        assert_eq!(model.m.row(i), explicit.m.row(i));
    }
}

#[test]
fn density_scaling_reaches_original_elastic_frequencies_and_acoustic_pressure_ports() {
    let shell = survey();
    let rho = [[RHO, 4.0 * RHO, RHO]];
    let heavy = rho.map(|row| row.map(|value| 4.0 * value));
    let assemble = |density: &[[f64; 3]], nodes: &[usize], support| {
        shell
            .assemble_with_reference_density(density, nodes, support)
            .unwrap()
    };
    let base = assemble(&rho, &[0, 1], ShellSupport::Clamped);
    let scaled = assemble(&heavy, &[0, 1], ShellSupport::Clamped);
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
    for i in 0..base.free {
        assert_eq!(base.k.row(i), scaled.k.row(i));
        for j in 0..base.free {
            close(scaled.m.get(i, j), 4.0 * base.m.get(i, j));
        }
    }
    for (a, b) in modes(&base).iter().zip(modes(&scaled)) {
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
    let amplitude = |density: &[[f64; 3]]| {
        let model = assemble(density, &[], ShellSupport::Free);
        let modes =
            mass_normalized_modes(&shell.mesh, &model, &[[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]]).unwrap();
        port.project_modes(&[modes[0].phi.clone()], 100, || false)
            .unwrap()[0]
    };
    let light_port = amplitude(&rho);
    close(
        light_port,
        0.5 / shell.mass_with_reference_density_kg(&rho).unwrap().sqrt(),
    );
    close(amplitude(&heavy), 0.5 * light_port);
}

#[test]
fn incomplete_or_invalid_density_never_reaches_a_shell_model() {
    let shell = survey();
    for rho in [
        vec![],
        vec![[RHO; 3]; 2],
        vec![[RHO, 0.0, RHO]],
        vec![[RHO, -1.0, RHO]],
        vec![[RHO, f64::NAN, RHO]],
        vec![[RHO, f64::INFINITY, RHO]],
        vec![[f64::from_bits(1); 3]], // Positive input, but no representable section mass.
    ] {
        assert!(
            shell
                .assemble_with_reference_density(&rho, &[], ShellSupport::Free)
                .is_err()
        );
        assert!(shell.mass_with_reference_density_kg(&rho).is_err());
    }
    let mut large = shell.clone();
    for point in &mut large.mesh.nodes {
        point[0] *= 1000.0;
        point[1] *= 1000.0;
    }
    let overflow = [[f64::MAX, 0.5 * f64::MAX, 0.25 * f64::MAX]];
    assert!(
        large
            .assemble_with_reference_density(&overflow, &[], ShellSupport::Free)
            .is_err()
    );
    assert!(large.mass_with_reference_density_kg(&overflow).is_err());
    let mut changed = shell.clone();
    changed.nodal_thickness_m = vec![H];
    assert!(
        changed
            .assemble_with_reference_density(&[[RHO; 3]], &[], ShellSupport::Free)
            .is_err()
    );
    assert!(changed.mass_with_reference_density_kg(&[[RHO; 3]]).is_err());
    assert_eq!(shell.nodal_thickness_m, [H, 2.0 * H, H]);
}
