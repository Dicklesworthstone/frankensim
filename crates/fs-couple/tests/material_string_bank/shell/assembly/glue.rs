//! Adhesive data reach assembled modal loss and the original pressure history.
//! Geometry and coefficients are manufactured integration fixtures.
use super::*;
use fs_plate::shell::junction::{FaceTie, GlueLine, ShellComponent, tie_shells_with_faces};
use fs_plate::shell::mitc::assemble_mitc3p;

fn glued_body(face: bool, eta: f64, gate: &CancelGate) -> SharedMaterialShellAssemblyBody {
    let mut result = body(gate);
    let glue = GlueLine {
        normal_n_m3: 2e7,
        shear_n_m3: [3e7, 8e6],
        viscosity_n_s_m3: Some([eta, 0.3 * eta, 1.7 * eta]),
        width_m: 0.012,
    };
    if face {
        // Coincident midsurfaces are an explicit zero-offset numerical fixture.
        // Unequal inertias prevent an accidentally undamped in-phase subspace.
        let mut right = part(0, 1400.0, gate);
        right.name = "part-1".into();
        result.specimen.components[1] = right;
        result.specimen.edge_ties.clear();
        result.specimen.face_ties = vec![FaceTie {
            master_component: 0,
            master_tris: (0..18).collect(),
            slave_component: 1,
            slave_tris: (0..18).collect(),
            offset: TieOffset::Midsurface,
            rotations: RotationCoupling::Full,
            kind: TieKind::Compliant(glue),
        }];
        result.specimen.tie_budget.max_segments = 256;
        result.specimen.tie_budget.max_slave_dofs = 192;
    } else {
        result.specimen.edge_ties[0].kind = TieKind::Compliant(glue);
    }
    for (index, rho) in [700.0, if face { 1400.0 } else { 700.0 }]
        .into_iter()
        .enumerate()
    {
        result.specimen = result
            .specimen
            .with_region_material(
                &format!("part-{index}"),
                "wood",
                &shell_material(rho, 0.0),
                gate,
            )
            .unwrap();
    }
    result
}

#[test]
fn g1_edge_and_face_glue_keep_the_original_physical_dissipation_rate() {
    with_cx(|_, gate| {
        for face in [false, true] {
            let body = glued_body(face, 2000.0, gate);
            assert!(body.specimen.has_viscous_glue());
            let image = body
                .specimen
                .reduce(&body.reduction, &body.force_ports, gate)
                .unwrap();
            let sources = &body.specimen.components;
            let models: Vec<_> = sources
                .iter()
                .map(|part| {
                    assemble_mitc3p(
                        part.specimen.mesh(),
                        &part.directors,
                        part.specimen.sections(),
                        &part.constraints,
                        part.mass,
                        part.assembly_budget,
                        part.drilling_alpha,
                        || false,
                    )
                    .unwrap()
                })
                .collect();
            let thickness: Vec<_> = sources
                .iter()
                .map(|part| {
                    vec![part.specimen.regions()[0].thickness_m(); part.specimen.mesh().nodes.len()]
                })
                .collect();
            let components: Vec<_> = sources
                .iter()
                .enumerate()
                .map(|(id, part)| ShellComponent {
                    id,
                    mesh: part.specimen.mesh(),
                    model: models[id].pencil(),
                    nodal_thickness_m: &thickness[id],
                })
                .collect();
            let tied = tie_shells_with_faces(
                &components,
                &body.specimen.edge_ties,
                &body.specimen.face_ties,
                body.specimen.tie_budget,
                || false,
            )
            .unwrap();
            let n = body.reduction.n_modes;
            let mut basis = vec![vec![0.0; tied.model.free]; n];
            for (mode, reduced) in basis.iter_mut().enumerate() {
                let physical = image
                    .physical_modes()
                    .iter()
                    .flat_map(|part| part[mode].iter());
                for (value, free) in physical.zip(&tied.model.dof_map) {
                    if let Some(index) = free {
                        reduced[*index] = *value;
                    }
                }
            }
            let matrix = image.damping().matrix_per_s();
            assert!(matrix.iter().any(|entry| entry.abs() > 1e-12));
            for velocity in [[1.0, -0.7, 0.4], [0.2, 0.3, -0.8]] {
                let mut physical = vec![0.0; tied.model.free];
                for (weight, mode) in velocity.iter().zip(&basis) {
                    for (value, shape) in physical.iter_mut().zip(mode) {
                        *value += weight * shape;
                    }
                }
                let expected = tied.glue_dissipation_rate(&physical).unwrap();
                let actual: f64 = (0..n)
                    .map(|i| {
                        (0..n)
                            .map(|j| velocity[i] * matrix[n * i + j] * velocity[j])
                            .sum::<f64>()
                    })
                    .sum();
                assert!(expected >= 0.0);
                assert!(
                    (actual - expected).abs() <= 5e-11 * expected.abs().max(actual.abs()).max(1.0),
                    "face={face}: modal {actual}, physical {expected}"
                );
            }
            let twice = glued_body(face, 4000.0, gate);
            let twice = twice
                .specimen
                .reduce(&twice.reduction, &twice.force_ports, gate)
                .unwrap();
            assert_eq!(image.physical_modes(), twice.physical_modes());
            assert_eq!(image.force_ports(), twice.force_ports());
            assert_eq!(
                image.volume_shapes_m2_per_sqrt_kg(),
                twice.volume_shapes_m2_per_sqrt_kg()
            );
            for (a, b) in matrix.iter().zip(twice.damping().matrix_per_s()) {
                assert!((2.0 * a - b).abs() <= 5e-11 * a.abs().max(b.abs()).max(1.0));
            }
        }
    });
}

#[test]
fn g3_glue_viscosity_changes_shared_motion_and_pressure_with_unchanged_elastic_modes() {
    with_cx(|cx, gate| {
        let dry = glued_body(false, 0.0, gate);
        assert!(!dry.specimen.has_viscous_glue());
        let mut dry = build(cx, gate, dry);
        let mut wet = build(cx, gate, glued_body(false, 2000.0, gate));
        assert_eq!(
            dry.source()
                .material_shell_assembly_image()
                .unwrap()
                .physical_modes(),
            wet.source()
                .material_shell_assembly_image()
                .unwrap()
                .physical_modes()
        );
        assert!(
            dry.source()
                .material_shell_assembly_image()
                .unwrap()
                .damping()
                .matrix_per_s()
                .iter()
                .all(|x| *x == 0.0)
        );
        let zero = vec![0.0; dry.system().mode_count()];
        let mut dry_pressure = Vec::new();
        let mut wet_pressure = Vec::new();
        for _ in 0..128 {
            let a = dry.step(&zero, drive()).unwrap();
            let b = wet.step(&zero, drive()).unwrap();
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
            dry_pressure.push(a.observer_pressure_pa);
            wet_pressure.push(b.observer_pressure_pa);
        }
        assert_ne!(dry_pressure, wet_pressure);
        assert!(wet_pressure.iter().any(|p| p.abs() > 1e-15));
        assert_ne!(
            dry.system().components()[2].states(),
            wet.system().components()[2].states()
        );
    });
}

#[test]
fn g0_compliant_supports_remain_active_and_glue_loss_reserves_runtime_channels() {
    with_cx(|cx, gate| {
        let mut body = glued_body(false, 2000.0, gate);
        let node = body.specimen.edge_ties[0].b.nodes[0];
        body.specimen.components[1]
            .constraints
            .push(ShellConstraint {
                node,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            });
        let image = body
            .specimen
            .reduce(&body.reduction, &body.force_ports, gate)
            .unwrap();
        for mode in &image.physical_modes()[1] {
            assert!(mode[6 * node..6 * node + 6].iter().all(|v| *v == 0.0));
        }
        let mut rigid = body.clone();
        rigid.specimen.edge_ties[0].kind = TieKind::Rigid;
        assert!(
            rigid
                .specimen
                .reduce(&rigid.reduction, &rigid.force_ports, gate)
                .is_err()
        );
        let mut limited = options(0);
        limited.network.max_connections = 2;
        let error = ReciprocalStringBank::from_material_shell_assembly(
            cx,
            vec![member(7800.0).source, layered().source],
            gas(101_325.0),
            body,
            None,
            law(),
            limited,
            gate,
        )
        .err()
        .expect("glue damping channels must be admitted before assembly");
        assert!(
            error
                .to_string()
                .contains("shell assembly loss, bridges and cavity")
        );
    });
}
