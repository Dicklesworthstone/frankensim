use super::*;
use crate::PlateSection;
use crate::shell::stiffened::{
    BeamMass, BeamSection, ShellBeam, assemble_stiffened_shell_materials_with_reference_density,
    assemble_stiffened_shell_materials_with_reference_density_and_body_acceleration,
    assemble_stiffened_shell_with_reference_density,
    assemble_stiffened_shell_with_reference_density_and_body_acceleration,
};
use crate::shell::{ShellMesh, ShellModel, ShellSection, ShellSupport};
use fs_material::sheet::SheetMatrix;

const PANEL: [[f64; 3]; 1] = [[200.0, 500.0, 1100.0]];
const RIB: [[f64; 2]; 1] = [[700.0, 1700.0]];

fn fixture() -> (ShellMesh, PlateSection, ShellBeam) {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let plate = PlateSection::isotropic(3e9, 0.25, 0.02, 1000.0).unwrap();
    let beam = ShellBeam {
        nodes: [0, 1],
        offsets_m: [[0.04, -0.03, 0.05]; 2],
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
    (mesh, plate, beam)
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * b.abs().max(0.01),
        "{a:.16e} != {b:.16e}"
    );
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn resultant(mesh: &ShellMesh, load: &[f64]) -> ([f64; 3], [f64; 3]) {
    let mut force = [0.0; 3];
    let mut moment = [0.0; 3];
    for (node, &position) in mesh.nodes.iter().enumerate() {
        let f = core::array::from_fn(|a| load[6 * node + a]);
        let orbital = cross(position, f);
        for a in 0..3 {
            force[a] += f[a];
            moment[a] += orbital[a] + load[6 * node + 3 + a];
        }
    }
    (force, moment)
}
fn identical(a: &ShellModel, b: &ShellModel) {
    assert_eq!(a.dof_map, b.dof_map);
    assert_eq!(a.free, b.free);
    assert_eq!(a.k, b.k);
    assert_eq!(a.m, b.m);
}

#[test]
fn full_mass_acceleration_matches_independent_weight_centroid_and_offset_torque() {
    let (mesh, plate, beam) = fixture();
    let acceleration = [1.2, -3.4, -9.81];
    // Exact polynomial volume integrals, independent of the native mass rule:
    // panel rho=200+300x+900y, h=.02 on the unit right triangle.
    let mass = 6.0 + 0.001 * (700.0 + 1700.0) / 2.0;
    let first = [
        23.0 / 12.0 + 1.2 * 0.04 + 0.001 * (700.0 + 2.0 * 1700.0) / 6.0,
        29.0 / 12.0 - 1.2 * 0.03,
        1.2 * 0.05,
    ];
    let expected_moment = cross(first, acceleration);
    for image in [BeamMass::Consistent, BeamMass::Lumped] {
        let (model, response) =
            assemble_stiffened_shell_with_reference_density_and_body_acceleration(
                &mesh,
                &[plate],
                &PANEL,
                &[],
                ShellSupport::Free,
                &[beam],
                &RIB,
                image,
            )
            .unwrap();
        let load = response.full_load(acceleration).unwrap();
        let rigid: Vec<_> = (0..18)
            .map(|i| if i % 6 < 3 { acceleration[i % 6] } else { 0.0 })
            .collect();
        let mut direct = vec![0.0; 18];
        model.m.spmv(&rigid, &mut direct);
        for (&actual, expected) in load.iter().zip(direct) {
            near(actual, expected);
        }
        let (force, moment) = resultant(&mesh, &load);
        for a in 0..3 {
            near(force[a], mass * acceleration[a]);
            near(moment[a], expected_moment[a]);
        }
        // The eccentric rib supplies physical nodal moments even under uniform g.
        assert!(
            load.iter()
                .enumerate()
                .any(|(i, f)| i % 6 >= 3 && f.abs() > 0.01)
        );
    }
}

#[test]
fn support_restriction_keeps_all_mass_columns_and_original_pencil_bits() {
    let (mesh, plate, beam) = fixture();
    let acceleration = [0.8, -0.3, -9.81];
    for image in [BeamMass::Consistent, BeamMass::Lumped] {
        let (_, full) = assemble_stiffened_shell_with_reference_density_and_body_acceleration(
            &mesh,
            &[plate],
            &PANEL,
            &[],
            ShellSupport::Free,
            &[beam],
            &RIB,
            image,
        )
        .unwrap();
        let expected = full.full_load(acceleration).unwrap();
        for support in [ShellSupport::Pinned, ShellSupport::Clamped] {
            let (model, response) =
                assemble_stiffened_shell_with_reference_density_and_body_acceleration(
                    &mesh,
                    &[plate],
                    &PANEL,
                    &[0, 2],
                    support,
                    &[beam],
                    &RIB,
                    image,
                )
                .unwrap();
            let ordinary = assemble_stiffened_shell_with_reference_density(
                &mesh,
                &[plate],
                &PANEL,
                &[0, 2],
                support,
                &[beam],
                &RIB,
                image,
            )
            .unwrap();
            identical(&model, &ordinary);
            let actual = response.full_load(acceleration).unwrap();
            assert_eq!(actual, expected);
            assert!(
                model
                    .dof_map
                    .iter()
                    .zip(&actual)
                    .any(|(d, f)| d.is_none() && f.abs() > 0.1)
            );
            let mut free_acceleration = vec![0.0; model.free];
            for (full, free) in model.dof_map.iter().enumerate() {
                if let Some(i) = free {
                    if full % 6 < 3 {
                        free_acceleration[*i] = acceleration[full % 6];
                    }
                }
            }
            let mut reduced_only = vec![0.0; model.free];
            model.m.spmv(&free_acceleration, &mut reduced_only);
            assert!(model.dof_map.iter().enumerate().any(|(full, free)| {
                free.is_some_and(|i| (reduced_only[i] - actual[full]).abs() > 0.1)
            }));
        }
    }
}

#[test]
fn rigid_frame_changes_rotate_both_force_and_reference_node_moment() {
    let (mesh, plate, beam) = fixture();
    let acceleration = [1.2, -3.4, -9.81];
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let (_, response) = assemble_stiffened_shell_with_reference_density_and_body_acceleration(
        &mesh,
        &[plate],
        &PANEL,
        &[0],
        ShellSupport::Clamped,
        &[beam],
        &RIB,
        BeamMass::Consistent,
    )
    .unwrap();
    let expected = response.full_load(acceleration).unwrap();
    let mut moved_mesh = mesh.clone();
    for node in &mut moved_mesh.nodes {
        let r = rotate(*node);
        *node = core::array::from_fn(|a| r[a] + [2.0, -3.0, 1.25][a]);
    }
    let mut moved_beam = beam;
    moved_beam.offsets_m = beam.offsets_m.map(rotate);
    moved_beam.section_z = rotate(beam.section_z);
    let (_, moved) = assemble_stiffened_shell_with_reference_density_and_body_acceleration(
        &moved_mesh,
        &[plate],
        &PANEL,
        &[0],
        ShellSupport::Clamped,
        &[moved_beam],
        &RIB,
        BeamMass::Consistent,
    )
    .unwrap();
    let actual = moved.full_load(rotate(acceleration)).unwrap();
    for node in 0..3 {
        for block in 0..2 {
            let transformed = rotate(core::array::from_fn(|a| expected[6 * node + 3 * block + a]));
            for a in 0..3 {
                near(actual[6 * node + 3 * block + a], transformed[a]);
            }
        }
    }
}

#[test]
fn graded_inventory_changes_weight_moment_and_scales_every_applied_load() {
    let (mesh, plate, beam) = fixture();
    let acceleration = [0.0, 0.0, -10.0];
    let prepare = |panels: &[[f64; 3]], ribs: &[[f64; 2]]| {
        assemble_stiffened_shell_with_reference_density_and_body_acceleration(
            &mesh,
            &[plate],
            panels,
            &[0, 2],
            ShellSupport::Clamped,
            &[beam],
            ribs,
            BeamMass::Consistent,
        )
        .unwrap()
    };
    let (model, original) = prepare(&PANEL, &RIB);
    let (relocated_model, relocated) = prepare(&PANEL, &[[1700.0, 700.0]]);
    let (scaled_model, scaled) = prepare(&[[400.0, 1000.0, 2200.0]], &[[1400.0, 3400.0]]);
    assert_eq!(model.k, relocated_model.k);
    assert_eq!(model.k, scaled_model.k);
    let initial = original.full_load(acceleration).unwrap();
    let shifted = relocated.full_load(acceleration).unwrap();
    let doubled = scaled.full_load(acceleration).unwrap();
    let (force, moment) = resultant(&mesh, &initial);
    let (shifted_force, shifted_moment) = resultant(&mesh, &shifted);
    for a in 0..3 {
        near(force[a], shifted_force[a]);
    }
    near(
        moment[1] - shifted_moment[1],
        10.0 * 0.001 * (1700.0 - 700.0) / 6.0,
    );
    assert!((initial[8] - shifted[8]).abs() > 0.1);
    for (&a, b) in initial.iter().zip(doubled) {
        near(b, 2.0 * a);
    }
}

#[test]
fn material_sections_empty_reinforcement_and_refusals_share_native_owners() {
    let (mesh, plate, beam) = fixture();
    let membrane = plate
        .d
        .map(|d| 12.0 * d / (plate.thickness * plate.thickness));
    let material = ShellSection::new(
        plate,
        SheetMatrix::orthotropic(membrane[0], membrane[4], membrane[1], membrane[8]).unwrap(),
    )
    .unwrap();
    let (model, response) =
        assemble_stiffened_shell_materials_with_reference_density_and_body_acceleration(
            &mesh,
            &[material],
            &PANEL,
            &[0],
            ShellSupport::Clamped,
            &[],
            &[],
            BeamMass::Consistent,
        )
        .unwrap();
    let ordinary = assemble_stiffened_shell_materials_with_reference_density(
        &mesh,
        &[material],
        &PANEL,
        &[0],
        ShellSupport::Clamped,
        &[],
        &[],
        BeamMass::Consistent,
    )
    .unwrap();
    identical(&model, &ordinary);
    let gravity = response.full_load([0.0, 0.0, -9.81]).unwrap();
    near(resultant(&mesh, &gravity).0[2], -6.0 * 9.81);
    assert_eq!(response.full_load([0.0; 3]).unwrap(), vec![0.0; 18]);
    for invalid in [
        [f64::NAN, 0.0, 0.0],
        [0.0, f64::INFINITY, 0.0],
        [0.0, 0.0, f64::MAX],
    ] {
        assert!(response.full_load(invalid).is_err());
        assert_eq!(response.full_load([0.0, 0.0, -9.81]).unwrap(), gravity);
    }
    let prepare = |panels: &[[f64; 3]], ribs: &[[f64; 2]]| {
        assemble_stiffened_shell_with_reference_density_and_body_acceleration(
            &mesh,
            &[plate],
            panels,
            &[0],
            ShellSupport::Clamped,
            &[beam],
            ribs,
            BeamMass::Consistent,
        )
    };
    assert!(prepare(&[], &RIB).is_err());
    assert!(prepare(&PANEL, &[]).is_err());
    assert!(prepare(&[[0.0, 500.0, 1100.0]], &RIB).is_err());
    assert!(prepare(&PANEL, &[[700.0, f64::NAN]]).is_err());
}
