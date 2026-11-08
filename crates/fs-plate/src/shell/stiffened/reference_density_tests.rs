//! Graded panel and bonded-beam inventory in the original sparse shell pencil.
use super::super::{
    CstDktMass, assemble_shell_sections_with_mass, assemble_shell_sections_with_reference_density,
};
use super::*;
use fs_material::sheet::SheetMatrix;

const PANEL_DENSITY: [[f64; 3]; 1] = [[200.0, 500.0, 1100.0]];
const BEAM_DENSITY: [[f64; 2]; 2] = [[700.0, 1700.0], [1200.0, 400.0]];
const THICKNESS: f64 = 0.02;

fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}
fn section() -> PlateSection {
    PlateSection::isotropic(3e9, 0.25, THICKNESS, 1000.0).unwrap()
}
fn beams() -> [ShellBeam; 2] {
    let first = ShellBeam {
        nodes: [0, 1],
        offsets_m: [[0.02, -0.03, 0.04]; 2],
        section_z: [0.0, 0.0, 1.0],
        section: BeamSection {
            young_pa: 12e9,
            shear_pa: 0.8e9,
            density_kg_m3: 450.0,
            area_m2: 0.001,
            iy_m4: 1.3e-7,
            iz_m4: 4.5e-8,
            torsion_m4: 8e-8,
        },
    };
    let mut second = first;
    second.nodes = [0, 2];
    second.offsets_m = [[-0.015, 0.025, -0.03]; 2];
    second.section.area_m2 = 0.0015;
    second.section.iy_m4 = 2.2e-7;
    second.section.iz_m4 = 8.1e-8;
    [first, second]
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(0.01),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}
fn bilinear(matrix: &Csr, left: &[f64], right: &[f64]) -> f64 {
    let mut action = vec![0.0; left.len()];
    matrix.spmv(right, &mut action);
    left.iter().zip(action).map(|(a, b)| a * b).sum()
}
fn rigid(mesh: &ShellMesh, axis: usize) -> Vec<f64> {
    let mut q = vec![0.0; 6 * mesh.nodes.len()];
    for (node, &position) in mesh.nodes.iter().enumerate() {
        if axis < 3 {
            q[6 * node + axis] = 1.0;
        } else {
            let mut omega = [0.0; 3];
            omega[axis - 3] = 1.0;
            q[6 * node..6 * node + 3].copy_from_slice(&cross(omega, position));
            q[6 * node + axis] = 1.0;
        }
    }
    q
}
fn rotate(v: [f64; 3]) -> [f64; 3] {
    [v[2], v[0], v[1]]
}

// Independent volume moments on the reference triangle and two prismatic ribs.
// These are polynomial antiderivatives, not production mass quadrature.
fn physical_moments() -> (f64, [f64; 3], [[f64; 3]; 3]) {
    // rho = 200 + 300*x + 900*y, h=0.02 on the unit right triangle.
    let mut mass = 6.0;
    let mut first = [23.0 / 12.0, 29.0 / 12.0, 0.0];
    let mut second = [
        [14.0 / 15.0, 17.0 / 30.0, 0.0],
        [17.0 / 30.0, 4.0 / 3.0, 0.0],
        [0.0, 0.0, 1.0 / 5000.0],
    ];
    let mesh = mesh();
    for (beam, rho) in beams().into_iter().zip(BEAM_DENSITY) {
        let start: [f64; 3] =
            core::array::from_fn(|a| mesh.nodes[beam.nodes[0]][a] + beam.offsets_m[0][a]);
        let end: [f64; 3] =
            core::array::from_fn(|a| mesh.nodes[beam.nodes[1]][a] + beam.offsets_m[1][a]);
        let chord: [f64; 3] = core::array::from_fn(|a| end[a] - start[a]);
        let length = dot(chord, chord).sqrt();
        let tangent = chord.map(|v| v / length);
        // Both authored section hints are perpendicular unit axes.
        let z = beam.section_z;
        let y = cross(z, tangent);
        let s = beam.section;
        let m0 = length * (rho[0] + rho[1]) / 2.0;
        let m1 = length * (rho[0] + 2.0 * rho[1]) / 6.0;
        let m2 = length * (rho[0] + 3.0 * rho[1]) / 12.0;
        mass += s.area_m2 * m0;
        for a in 0..3 {
            first[a] += s.area_m2 * (m0 * start[a] + m1 * chord[a]);
            for b in 0..3 {
                second[a][b] += s.area_m2
                    * (m0 * start[a] * start[b]
                        + m1 * (start[a] * chord[b] + chord[a] * start[b])
                        + m2 * chord[a] * chord[b])
                    + m0 * (s.iz_m4 * y[a] * y[b] + s.iy_m4 * z[a] * z[b]);
            }
        }
    }
    (mass, first, second)
}
fn spatial_inertia(mass: f64, first: [f64; 3], second: [[f64; 3]; 3]) -> [[f64; 6]; 6] {
    let mut expected = [[0.0; 6]; 6];
    let trace = second[0][0] + second[1][1] + second[2][2];
    for axis in 0..3 {
        expected[axis][axis] = mass;
        let mut omega = [0.0; 3];
        omega[axis] = 1.0;
        let moment = cross(omega, first);
        for a in 0..3 {
            expected[a][axis + 3] = moment[a];
            expected[axis + 3][a] = moment[a];
            expected[a + 3][axis + 3] = if a == axis {
                trace - second[a][axis]
            } else {
                -second[a][axis]
            };
        }
    }
    expected
}
fn identical(a: &ShellModel, b: &ShellModel) {
    assert_eq!(a.dof_map, b.dof_map);
    assert_eq!(a.free, b.free);
    assert_eq!(a.k, b.k);
    assert_eq!(a.m, b.m);
}

#[test]
fn coupled_graded_inventory_preserves_complete_spatial_inertia_and_frame_changes() {
    let (mass, first, second) = physical_moments();
    near(mass, 8.4);
    for rotated in [false, true] {
        let mut mesh = mesh();
        let mut beams = beams();
        let (first, second) = if rotated {
            for node in &mut mesh.nodes {
                *node = rotate(*node);
            }
            for beam in &mut beams {
                beam.offsets_m = beam.offsets_m.map(rotate);
                beam.section_z = rotate(beam.section_z);
            }
            let axes = [2, 0, 1];
            (
                rotate(first),
                core::array::from_fn(|a| core::array::from_fn(|b| second[axes[a]][axes[b]])),
            )
        } else {
            (first, second)
        };
        let model = assemble_stiffened_shell_with_reference_density(
            &mesh,
            &[section()],
            &PANEL_DENSITY,
            &[],
            ShellSupport::Free,
            &beams,
            &BEAM_DENSITY,
            BeamMass::Consistent,
        )
        .unwrap();
        let expected = spatial_inertia(mass, first, second);
        let motions: Vec<_> = (0..6).map(|axis| rigid(&mesh, axis)).collect();
        for a in 0..6 {
            for b in 0..6 {
                near(bilinear(&model.m, &motions[a], &motions[b]), expected[a][b]);
            }
        }
        // Translation/rotation cross terms also determine the torque of M*g.
        // A replacement of either density field by its mean loses these terms.
        near(bilinear(&model.m, &motions[1], &motions[5]), first[0]);
        near(bilinear(&model.m, &motions[0], &motions[5]), -first[1]);
    }
}

#[test]
fn offset_beam_density_velocity_correlation_reaches_the_selected_coupled_mass() {
    let mesh = mesh();
    let beam = beams()[0];
    let rho = BEAM_DENSITY[0];
    let plain = assemble_shell_sections_with_reference_density(
        &mesh,
        &[section()],
        &PANEL_DENSITY,
        &[],
        ShellSupport::Free,
    )
    .unwrap();
    // Centroidal v=(x,x^3,x^2), omega=(x,-2x,3x^2) on the unit beam.
    // Subtract omega cross offset to supply the actual shell-node velocity.
    let mut q = vec![0.0; 18];
    let omega = [1.0, -2.0, 3.0];
    let arm = cross(omega, beam.offsets_m[1]);
    for a in 0..3 {
        q[6 + a] = 1.0 - arm[a];
        q[9 + a] = omega[a];
    }
    let moment =
        |power: i32| rho[0] / f64::from(power + 1) + (rho[1] - rho[0]) / f64::from(power + 2);
    let s = beam.section;
    let consistent = s.area_m2 * (moment(2) + moment(4) + moment(6))
        + (5.0 * s.iy_m4 + s.iz_m4) * moment(2)
        + 9.0 * s.iz_m4 * moment(4);
    let lumped = (rho[0] + 2.0 * rho[1]) / 6.0 * (3.0 * s.area_m2 + 5.0 * s.iy_m4 + 10.0 * s.iz_m4);
    for (image, expected) in [
        (BeamMass::Consistent, consistent),
        (BeamMass::Lumped, lumped),
    ] {
        let model = assemble_stiffened_shell_with_reference_density(
            &mesh,
            &[section()],
            &PANEL_DENSITY,
            &[],
            ShellSupport::Free,
            &[beam],
            &[rho],
            image,
        )
        .unwrap();
        near(
            bilinear(&model.m, &q, &q) - bilinear(&plain.m, &q, &q),
            expected,
        );
    }
    assert!(lumped > consistent);
}

#[test]
fn density_redistribution_and_scaling_reach_supported_native_modes() {
    let mesh = mesh();
    let beam = beams()[0];
    let assemble = |panels: &[[f64; 3]], beam_density: &[[f64; 2]], boundary: &[usize]| {
        assemble_stiffened_shell_with_reference_density(
            &mesh,
            &[section()],
            panels,
            boundary,
            ShellSupport::Clamped,
            &[beam],
            beam_density,
            BeamMass::Consistent,
        )
        .unwrap()
    };
    let supported = assemble(&PANEL_DENSITY, &BEAM_DENSITY[..1], &[0, 2]);
    let full = assemble(&PANEL_DENSITY, &BEAM_DENSITY[..1], &[]);
    let reverse = [[BEAM_DENSITY[0][1], BEAM_DENSITY[0][0]]];
    let relocated = assemble(&PANEL_DENSITY, &reverse, &[0, 2]);
    let relocated_full = assemble(&PANEL_DENSITY, &reverse, &[]);
    let translation = rigid(&mesh, 0);
    near(
        bilinear(&full.m, &translation, &translation),
        bilinear(&relocated_full.m, &translation, &translation),
    );
    assert_eq!(supported.k, relocated.k);
    for (a, &row) in supported.dof_map.iter().enumerate() {
        if let Some(row) = row {
            for (b, &column) in supported.dof_map.iter().enumerate() {
                if let Some(column) = column {
                    near(supported.m.get(row, column), full.m.get(a, b));
                    near(supported.k.get(row, column), full.k.get(a, b));
                }
            }
        }
    }
    let spectrum = |model: &ShellModel| {
        let n = model.free;
        let k: Vec<_> = (0..n * n).map(|i| model.k.get(i / n, i % n)).collect();
        let m: Vec<_> = (0..n * n).map(|i| model.m.get(i / n, i % n)).collect();
        fs_modal::eigh_gen_dense(&k, &m, n).unwrap()
    };
    let original_modes = spectrum(&supported);
    let relocated_modes = spectrum(&relocated);
    assert!(
        original_modes
            .iter()
            .zip(&relocated_modes)
            .any(|(a, b)| (b.lambda / a.lambda - 1.0).abs() > 1e-4)
    );
    let heavy_panels = PANEL_DENSITY.map(|rho| rho.map(|v| 4.0 * v));
    let heavy_beams = [BEAM_DENSITY[0].map(|v| 4.0 * v)];
    let heavy = assemble(&heavy_panels, &heavy_beams, &[0, 2]);
    assert_eq!(supported.k, heavy.k);
    for (a, b) in original_modes.iter().zip(spectrum(&heavy)) {
        assert!(a.lambda > 0.0 && b.lambda > 0.0);
        assert!((b.lambda / a.lambda - 0.25).abs() < 2e-9);
    }
}

#[test]
fn uniform_beam_bits_material_sections_and_empty_reinforcement_keep_native_owners() {
    let mesh = mesh();
    let beams = beams();
    for image in [BeamMass::Consistent, BeamMass::Lumped] {
        for beam in beams {
            let ordinary = element_with_mass(&mesh, &beam, image).unwrap();
            let explicit = element_with_reference_density(
                &mesh,
                &beam,
                image,
                Some([beam.section.density_kg_m3; 2]),
            )
            .unwrap();
            assert_eq!(ordinary, explicit);
        }
    }
    let original =
        assemble_stiffened_shell(&mesh, &[section()], &[0], ShellSupport::Clamped, &beams).unwrap();
    let explicit_legacy = assemble_stiffened_shell_with_mass(
        &mesh,
        &[section()],
        &[0],
        ShellSupport::Clamped,
        &beams,
        BeamMass::Lumped,
    )
    .unwrap();
    identical(&original, &explicit_legacy);

    let empty = assemble_stiffened_shell_with_reference_density(
        &mesh,
        &[section()],
        &PANEL_DENSITY,
        &[0],
        ShellSupport::Clamped,
        &[],
        &[],
        BeamMass::Consistent,
    )
    .unwrap();
    let plain = assemble_shell_sections_with_reference_density(
        &mesh,
        &[section()],
        &PANEL_DENSITY,
        &[0],
        ShellSupport::Clamped,
    )
    .unwrap();
    identical(&empty, &plain);
    let uniform = assemble_stiffened_shell_with_reference_density(
        &mesh,
        &[section()],
        &[[section().density; 3]],
        &[0],
        ShellSupport::Clamped,
        &[],
        &[],
        BeamMass::Consistent,
    )
    .unwrap();
    identical(
        &uniform,
        &assemble_shell_sections_with_mass(
            &mesh,
            &[section()],
            &[0],
            ShellSupport::Clamped,
            CstDktMass::Consistent,
        )
        .unwrap(),
    );

    let plate = section();
    let membrane = plate.d.map(|d| 12.0 * d / (THICKNESS * THICKNESS));
    let material = ShellSection::new(
        plate,
        SheetMatrix::orthotropic(membrane[0], membrane[4], membrane[1], membrane[8]).unwrap(),
    )
    .unwrap();
    let isotropic = assemble_stiffened_shell_with_reference_density(
        &mesh,
        &[plate],
        &PANEL_DENSITY,
        &[0],
        ShellSupport::Clamped,
        &beams,
        &BEAM_DENSITY,
        BeamMass::Consistent,
    )
    .unwrap();
    let independent = assemble_stiffened_shell_materials_with_reference_density(
        &mesh,
        &[material],
        &PANEL_DENSITY,
        &[0],
        ShellSupport::Clamped,
        &beams,
        &BEAM_DENSITY,
        BeamMass::Consistent,
    )
    .unwrap();
    assert_eq!(isotropic.dof_map, independent.dof_map);
    assert_eq!(isotropic.m, independent.m);
}

#[test]
fn incomplete_or_invalid_fields_refuse_before_returning_a_coupled_model() {
    let mesh = mesh();
    let beams = beams();
    let assemble = |panels: &[[f64; 3]], reinforcement: &[[f64; 2]]| {
        assemble_stiffened_shell_with_reference_density(
            &mesh,
            &[section()],
            panels,
            &[],
            ShellSupport::Free,
            &beams,
            reinforcement,
            BeamMass::Consistent,
        )
    };
    for panels in [
        vec![],
        vec![[500.0; 3]; 2],
        vec![[200.0, 0.0, 1100.0]],
        vec![[200.0, -1.0, 1100.0]],
        vec![[200.0, f64::NAN, 1100.0]],
        vec![[200.0, f64::INFINITY, 1100.0]],
    ] {
        assert!(assemble(&panels, &BEAM_DENSITY).is_err());
    }
    for reinforcement in [
        vec![],
        vec![[700.0; 2]],
        vec![[700.0; 2]; 3],
        vec![[0.0, 1700.0], BEAM_DENSITY[1]],
        vec![[700.0, -1.0], BEAM_DENSITY[1]],
        vec![[700.0, f64::NAN], BEAM_DENSITY[1]],
        vec![[700.0, f64::INFINITY], BEAM_DENSITY[1]],
    ] {
        assert!(assemble(&PANEL_DENSITY, &reinforcement).is_err());
    }
    // Refuse the authored beam field before the malformed geometry can reach
    // any shell sparse allocation.
    let mut invalid_mesh = mesh.clone();
    invalid_mesh.nodes.clear();
    assert!(matches!(
        assemble_stiffened_shell_with_reference_density(
            &invalid_mesh,
            &[section()],
            &PANEL_DENSITY,
            &[],
            ShellSupport::Free,
            &beams,
            &[[0.0, 1700.0], BEAM_DENSITY[1]],
            BeamMass::Consistent,
        ),
        Err(PlateError::BadStiffener { .. }),
    ));
    // Kernel and sparse arithmetic refusals remain native.
    let mut large = mesh.clone();
    large.nodes[1][0] = 2.0;
    assert!(
        assemble_stiffened_shell_with_reference_density(
            &large,
            &[section()],
            &PANEL_DENSITY,
            &[],
            ShellSupport::Free,
            &beams[..1],
            &[[f64::MAX; 2]],
            BeamMass::Consistent,
        )
        .is_err()
    );
}
