//! G0/G3: signed beam tangents through the original nonmatching, eccentric mortar.
use super::*;
use crate::shell::stiffened::timoshenko::ShearAreas;
use crate::shell::stiffened::timoshenko::{assemble_beam_chain, assemble_prestressed_beam_chain};

fn fixture() -> (ShellMesh, ShellModel, TieBudget) {
    let mesh = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.4, 0.0],
            [0.0, 0.4, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();
    // A declared manufactured reference pencil isolates the transfer identity.
    // No shell element or measured-material accuracy is inferred from this test.
    // Node zero is supported, so its beam endpoint follows the actual elimination.
    let mut k = Coo::new(18, 18);
    let mut m = Coo::new(18, 18);
    for i in 0..18 {
        k.push(i, i, 100.0 * (i + 1) as f64);
        m.push(i, i, 0.01 * (i + 1) as f64);
    }
    let model = ShellModel {
        k: k.assemble(),
        m: m.assemble(),
        free: 18,
        dof_map: (0..24)
            .map(|i| if i < 6 { None } else { Some(i - 6) })
            .collect(),
    };
    let budget = TieBudget {
        max_segments: 64,
        max_slave_dofs: 32,
        gap_tolerance_m: 1e-12,
        angle_tolerance_rad: 1e-6,
    };
    (mesh, model, budget)
}

fn beam(nodes: &[[f64; 3]]) -> BeamOnShell<'_> {
    BeamOnShell {
        beam_nodes_m: nodes,
        section_z: [0.0, 0.0, 1.0],
        section: BeamSection {
            young_pa: 12e9,
            shear_pa: 0.75e9,
            density_kg_m3: 450.0,
            area_m2: 6e-4,
            iy_m4: 4.5e-8,
            iz_m4: 2e-8,
            torsion_m4: 3e-8,
        },
        kind: BeamKind::Timoshenko(ShearAreas {
            y_m2: 5e-4,
            z_m2: 5e-4,
        }),
        master: TieSide {
            component: 0,
            nodes: vec![0, 1],
        },
        side: 1,
        centroid_height_m: 0.035,
    }
}

fn product(a: &Csr, x: &[f64], y: &[f64]) -> f64 {
    let mut ay = vec![0.0; y.len()];
    a.spmv(y, &mut ay);
    x.iter().zip(ay).map(|(x, y)| x * y).sum()
}

#[test]
fn prestress_uses_the_original_eccentric_map_with_signed_virtual_work() {
    let (mesh, model, budget) = fixture();
    let thickness = [0.01; 4];
    let components = [ShellComponent {
        id: 0,
        mesh: &mesh,
        model: &model,
        nodal_thickness_m: &thickness,
    }];
    let tied = tie_shells(&components, &[], budget, || false).unwrap();
    let nodes = [[0.0, 0.0, 0.04], [0.35, 0.0, 0.04], [1.0, 0.0, 0.04]];
    let beams = [beam(&nodes)];
    let original = tie_beams_to_shells(&tied, &components, &beams, budget, || false).unwrap();
    let zero = tie_prestressed_beams_to_shells(
        &tied,
        &components,
        &beams,
        &[&[0.0, -0.0]],
        budget,
        || false,
    )
    .unwrap();
    for r in 0..model.free {
        assert_eq!(original.model.k.row(r), zero.model.k.row(r));
        assert_eq!(original.model.m.row(r), zero.model.m.row(r));
    }
    let u: Vec<_> = (0..model.free)
        .map(|i| ((i * 13 + 1) as f64).sin() * 0.003)
        .collect();
    let v: Vec<_> = (0..model.free)
        .map(|i| ((i * 7 + 4) as f64).cos() * 0.004)
        .collect();
    let bu = original.expand_beam(0, &u).unwrap();
    let bv = original.expand_beam(0, &v).unwrap();
    assert!(bu[..6].iter().all(|x| x.abs() < 1e-16));
    // The eccentric beam endpoint must retain rotation-induced transverse motion.
    assert!((bu[13] - (u[1] - 0.04 * u[3])).abs() < 1e-15);
    let (unloaded, _) =
        assemble_beam_chain(&nodes, beams[0].section_z, beams[0].section, beams[0].kind).unwrap();
    for forces in [[120.0, 120.0], [-120.0, -120.0], [120.0, -70.0]] {
        let loaded =
            tie_prestressed_beams_to_shells(&tied, &components, &beams, &[&forces], budget, || {
                false
            })
            .unwrap();
        let (bk, _) = assemble_prestressed_beam_chain(
            &nodes,
            beams[0].section_z,
            beams[0].section,
            beams[0].kind,
            &forces,
        )
        .unwrap();
        assert_eq!(loaded.expand_beam(0, &u).unwrap(), bu);
        for r in 0..model.free {
            assert_eq!(loaded.model.m.row(r), original.model.m.row(r));
        }
        let transferred = product(&loaded.model.k, &u, &v) - product(&original.model.k, &u, &v);
        let physical = product(&bk, &bu, &bv) - product(&unloaded, &bu, &bv);
        assert!(physical.abs() > 1e-7);
        assert!(
            (transferred - physical).abs() < 2e-10 * physical.abs(),
            "{transferred:e} != {physical:e}"
        );
        let reverse = product(&loaded.model.k, &v, &u) - product(&original.model.k, &v, &u);
        assert!((reverse - transferred).abs() < 2e-10 * physical.abs());
        if forces[0] == forces[1] {
            let work = product(&loaded.model.k, &u, &u) - product(&original.model.k, &u, &u);
            assert_eq!(work.is_sign_positive(), forces[0].is_sign_positive());
        }
    }
}

#[test]
fn prestressed_ties_refuse_incomplete_force_fields_budgets_and_cancellation() {
    let (mesh, model, budget) = fixture();
    let thickness = [0.01; 4];
    let components = [ShellComponent {
        id: 0,
        mesh: &mesh,
        model: &model,
        nodal_thickness_m: &thickness,
    }];
    let tied = tie_shells(&components, &[], budget, || false).unwrap();
    let nodes = [[0.0, 0.0, 0.04], [0.35, 0.0, 0.04], [1.0, 0.0, 0.04]];
    let beams = [beam(&nodes)];
    assert!(
        tie_prestressed_beams_to_shells(&tied, &components, &beams, &[], budget, || false).is_err()
    );
    for forces in [
        &[1.0][..],
        &[1.0, 2.0, 3.0][..],
        &[f64::NAN, 0.0][..],
        &[0.0, f64::INFINITY][..],
    ] {
        assert!(
            tie_prestressed_beams_to_shells(&tied, &components, &beams, &[forces], budget, || {
                false
            })
            .is_err()
        );
    }
    let forces = [1.0, 2.0];
    assert!(matches!(
        tie_prestressed_beams_to_shells(
            &tied,
            &components,
            &beams,
            &[&forces],
            TieBudget {
                max_slave_dofs: 17,
                ..budget
            },
            || false
        ),
        Err(PlateError::Budget { .. })
    ));
    assert!(matches!(
        tie_prestressed_beams_to_shells(&tied, &components, &beams, &[&forces], budget, || true),
        Err(PlateError::Cancelled)
    ));
    let mut calls = 0;
    assert!(matches!(
        tie_prestressed_beams_to_shells(&tied, &components, &beams, &[&forces], budget, || {
            calls += 1;
            calls > 4
        }),
        Err(PlateError::Cancelled)
    ));
    // The immutable original body remains usable after either refusal.
    assert!(
        tie_prestressed_beams_to_shells(&tied, &components, &beams, &[&forces], budget, || false)
            .is_ok()
    );
}
