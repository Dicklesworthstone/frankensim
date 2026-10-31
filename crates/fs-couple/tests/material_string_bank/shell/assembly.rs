//! Original mortar, material and reciprocal owners on a synthetic two-part body.
//! These are numerical composition checks, not instrument qualification.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::shell::SharedMaterialShellAssemblyBody;
use fs_couple::shear_shell::assembly::{
    MaterialShellAssembly, MaterialShellComponent, ShellAssemblyForcePort,
    ShellAssemblyReductionOptions,
};
use fs_couple::shear_shell::pressure::ShellPressureSurfaces;
use fs_plate::shell::junction::{
    EdgeTie, RotationCoupling, TieBudget, TieKind, TieOffset, TieSide,
};
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace, ShellPressurePort};
use fs_plate::shell::profile::ProfileBudget;

fn surfaces(triangles: usize, exterior: bool) -> ShellPressureSurfaces {
    let faces = |side| {
        (0..triangles)
            .map(|triangle| OrientedShellFacet { triangle, side })
            .collect()
    };
    ShellPressureSurfaces {
        exterior: if exterior {
            faces(ShellFace::Positive)
        } else {
            Vec::new()
        },
        cavity: faces(ShellFace::Negative),
    }
}

fn part(index: usize, rho: f64, gate: &CancelGate) -> MaterialShellComponent {
    let flat = PlateMesh::rectangle(0.2, 0.3, 3, 3);
    let nodes = flat
        .nodes
        .iter()
        .map(|&(x, y)| [x + 0.2 * index as f64, y, 0.0])
        .collect();
    let mesh = ShellMesh::new(nodes, flat.tris).unwrap();
    let directors =
        NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; mesh.nodes.len()], 0.99).unwrap();
    let constraints = if index == 0 {
        mesh.nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p[0] == 0.0)
            .map(|(node, _)| ShellConstraint {
                node,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            })
            .collect()
    } else {
        Vec::new()
    };
    let budget = MitcBudget {
        max_nodes: mesh.nodes.len(),
        max_triangles: mesh.tris.len(),
    };
    let region = ShellRegionMaterial {
        region: PlateRegionMaterial {
            region: PlateRegion {
                name: "wood".into(),
                triangle_indices: (0..mesh.tris.len()).collect(),
            },
            material: shell_material(rho, 1000.0),
            model: PlateMaterialModel::Orthotropic12,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.003),
            material_angle_rad: 0.17,
        },
        shear_correction: 5.0 / 6.0,
        axis1_global: [1.0, 0.0, 0.0],
    };
    let pressure = surfaces(mesh.tris.len(), true);
    let specimen = compile_material_shear_shell(mesh, vec![region], budget, gate)
        .unwrap()
        .with_pressure_surfaces(pressure, gate)
        .unwrap();
    MaterialShellComponent {
        name: format!("part-{index}"),
        specimen,
        directors,
        constraints,
        mass: ShellMass::Consistent,
        drilling_alpha: 1e-3,
        assembly_budget: budget,
    }
}

fn body(gate: &CancelGate) -> SharedMaterialShellAssemblyBody {
    let left = part(0, 700.0, gate);
    let right = part(1, 700.0, gate);
    let seam = |part: &MaterialShellComponent| {
        part.specimen
            .mesh()
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| (p[0] - 0.2).abs() < 1e-12)
            .map(|(node, _)| node)
            .collect()
    };
    let tie = EdgeTie {
        a: TieSide {
            component: 0,
            nodes: seam(&left),
        },
        b: TieSide {
            component: 1,
            nodes: seam(&right),
        },
        offset: TieOffset::Midsurface,
        rotations: RotationCoupling::Full,
        kind: TieKind::Rigid,
    };
    let width = 6 * left.specimen.mesh().nodes.len();
    let force_ports = (0..2)
        .map(|part| {
            let mut rows = vec![vec![0.0; width]; 2];
            rows[part][6 * 5 + 2] = 1.0;
            ShellAssemblyForcePort { components: rows }
        })
        .collect();
    SharedMaterialShellAssemblyBody {
        specimen: MaterialShellAssembly {
            components: vec![left, right],
            edge_ties: vec![tie],
            face_ties: Vec::new(),
            tie_budget: TieBudget {
                max_segments: 32,
                max_slave_dofs: 24,
                gap_tolerance_m: 1e-10,
                angle_tolerance_rad: 1e-10,
            },
        },
        reduction: ShellAssemblyReductionOptions {
            eigenvalue_window: (0.0, 1e9),
            n_modes: 3,
            projection_budget: ShellProjectionBudget {
                max_modes: 8,
                max_projection_terms: 80_000_000,
                reference_resistance_n_s_m: 1.0,
            },
            max_total_nodes: 32,
            max_total_triangles: 36,
        },
        force_ports,
        bridges: vec![
            SharedPlateBridge {
                station_fraction: 0.8,
                stiffness_n_m: 3000.0,
                damping_n_s_m: 0.1,
                rest_extension_m: 0.0
            };
            2
        ],
        modal_budget: ModalAcousticTimeBudget::audible_reference(),
        listener_m: 1.0,
    }
}

fn build(
    cx: &Cx<'_>,
    gate: &CancelGate,
    body: SharedMaterialShellAssemblyBody,
) -> ReciprocalStringBank {
    let mut options = options(0);
    options.network.max_connections = 16;
    options.network.max_setup_terms = 131_072;
    let gas = gas(101_325.0);
    let cavity = fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity {
        cavity: fs_couple::cavity_modal::HelmholtzModalCavity::new(
            fs_scenario::HelmholtzCavity {
                volume_m3: 0.01,
                neck_radius_m: 0.01,
                neck_length_m: 0.02,
                inner_end: fs_scenario::NeckEnd::Flanged,
                outer_end: fs_scenario::NeckEnd::Flanged,
            },
            gas,
            fs_couple::cavity_modal::CavityNeckLoss::Resistance(100.0),
        )
        .unwrap(),
        modal_budget: ModalAcousticTimeBudget::audible_reference(),
        listener_m: 1.0,
    };
    ReciprocalStringBank::from_material_shell_assembly(
        cx,
        vec![member(7800.0).source, layered().source],
        gas,
        body,
        Some(cavity),
        law(),
        options,
        gate,
    )
    .unwrap()
}

#[test]
fn local_component_torque_uses_the_original_mortar_assembled_rotation_basis() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::shell::ShellRotationGeometry;
    with_cx(|cx,gate| {
        let bank = build(cx,gate,body(gate));
        let source = bank.source();
        let image = source.material_shell_assembly_image().unwrap();
        for component in 0..2 {
            let geometry = ShellRotationGeometry { shell_component: component,triangle: 0,
                barycentric: [0.2,0.3,0.5],axis: [0.6,0.8,0.0] };
            let rotation = source.shell_rotation(&geometry,gate).unwrap();
            // Original parts share body zero; component is not a new modal body.
            assert_eq!(rotation.component,0);
            let tri = source.material_shell_assembly().unwrap().specimen.components[component]
                .specimen.mesh().tris[0];
            for (mode,&angle) in image.physical_modes()[component].iter().zip(&rotation.angular_shapes) {
                let expected: f64 = (0..3).map(|a|geometry.barycentric[a]*(0..3)
                    .map(|i|geometry.axis[i]*mode[6*tri[a]+3+i]).sum::<f64>()).sum();
                assert!((angle-expected).abs() < 1e-12*expected.abs().max(1.0));
            }
        }
    });
}

#[test]
fn g3_tied_material_parts_feed_one_string_body_cavity_and_pressure_basis() {
    with_cx(|cx, gate| {
        let mut bank = build(cx, gate, body(gate));
        assert_eq!(bank.source().bridge_surface_materials(1).unwrap(),
            [material(15_600.0, 100e9).material().clone(), shell_material(700.0, 1000.0).material().clone()]);
        let image = bank.source().material_shell_assembly_image().unwrap();
        assert_eq!(image.source().components.len(), 2);
        assert_eq!(image.tie_receipt().slave_dofs, 24);
        assert_eq!(
            bank.system().components().len(),
            4,
            "two strings, one assembled body, one cavity"
        );
        let tie = &image.source().edge_ties[0];
        for mode in 0..3 {
            for (&a, &b) in tie.a.nodes.iter().zip(&tie.b.nodes) {
                for dof in 0..6 {
                    let left = image.physical_modes()[0][mode][6 * a + dof];
                    let right = image.physical_modes()[1][mode][6 * b + dof];
                    assert!((left - right).abs() <= 1e-11 * left.abs().max(right.abs()).max(1.0));
                }
            }
            let exterior: f64 = image
                .source()
                .components
                .iter()
                .enumerate()
                .map(|(i, part)| {
                    let selected = part.specimen.pressure_surfaces().unwrap().exterior.clone();
                    let row = ShellPressurePort::new(
                        part.specimen.mesh(),
                        selected,
                        ProfileBudget {
                            max_nodes: 16,
                            max_triangles: 18,
                            max_feature_evaluations: 0,
                        },
                        || false,
                    )
                    .unwrap();
                    row.physical_row_m2()
                        .iter()
                        .zip(&image.physical_modes()[i][mode])
                        .map(|(a, b)| a * b)
                        .sum::<f64>()
                })
                .sum();
            assert!((image.volume_shapes_m2_per_sqrt_kg()[mode] - exterior).abs() < 1e-12);
            assert!((image.cavity_volume_shapes_m2_per_sqrt_kg()[mode] + exterior).abs() < 1e-12);
        }
        for (row, member) in image.force_ports().iter().zip(bank.source().strings()) {
            assert_eq!(row, &member.connections[0].body.shapes);
        }
        let cavity = bank.source().material_cavity().unwrap();
        let expected = cavity
            .cavity
            .connections(0, 1, image.cavity_volume_shapes_m2_per_sqrt_kg())
            .unwrap();
        assert_eq!(
            &bank.source().body_connections()[..expected.len()],
            &expected
        );
        let zero = vec![0.0; bank.system().mode_count()];
        let mut peak = 0.0_f64;
        for _ in 0..96 {
            let frame = bank.step(&zero, drive()).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            peak = peak.max(frame.observer_pressure_pa.abs());
        }
        assert!(peak > 1e-15);
        assert!(
            bank.system().components()[1]
                .states()
                .iter()
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0)
        );
        assert!(
            bank.system().components()[3]
                .states()
                .iter()
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0)
        );
        let original_modes = bank
            .source()
            .material_shell_assembly_image()
            .unwrap()
            .physical_modes()
            .to_vec();
        let original_cavity = bank
            .source()
            .material_shell_assembly_image()
            .unwrap()
            .cavity_volume_shapes_m2_per_sqrt_kg()
            .to_vec();
        let multiple = bank.with_additional_bows(vec![
            fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MaterialBowContact {
                location: StringBowLocation { string: 1, station_fraction: 0.31 },
                law: law(), friction: options(0).friction,
            }], fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig {
                max_ports: 4, max_sweeps: 128, max_setup_terms: 131_072,
            }, gate).unwrap();
        let source = multiple.source().bank();
        assert_eq!(
            source
                .material_shell_assembly()
                .unwrap()
                .specimen
                .components
                .len(),
            2
        );
        assert_eq!(
            source
                .material_shell_assembly_image()
                .unwrap()
                .physical_modes(),
            original_modes
        );
        assert_eq!(
            source
                .material_shell_assembly_image()
                .unwrap()
                .cavity_volume_shapes_m2_per_sqrt_kg(),
            original_cavity
        );
    });
}

#[test]
fn g3_internal_part_keeps_mechanics_and_full_tensor_loss_without_fluid_double_counting() {
    with_cx(|_, gate| {
        let full = body(gate);
        let mut internal = full.clone();
        internal.specimen.components[1].specimen = internal.specimen.components[1]
            .specimen
            .with_pressure_surfaces(
                ShellPressureSurfaces {
                    exterior: Vec::new(),
                    cavity: Vec::new(),
                },
                gate,
            )
            .unwrap();
        let a = full
            .specimen
            .reduce(&full.reduction, &full.force_ports, gate)
            .unwrap();
        let b = internal
            .specimen
            .reduce(&internal.reduction, &internal.force_ports, gate)
            .unwrap();
        assert_eq!(a.physical_modes(), b.physical_modes());
        assert_eq!(a.force_ports(), b.force_ports());
        assert_eq!(a.damping().matrix_per_s(), b.damping().matrix_per_s());
        assert!(b.source().components[1].specimen.mass_kg() > 0.0);
        assert!(b.physical_modes()[1].iter().flatten().any(|v| *v != 0.0));
        let part = &internal.specimen.components[0];
        let row = ShellPressurePort::new(
            part.specimen.mesh(),
            part.specimen.pressure_surfaces().unwrap().exterior.clone(),
            ProfileBudget {
                max_nodes: 16,
                max_triangles: 18,
                max_feature_evaluations: 0,
            },
            || false,
        )
        .unwrap();
        for mode in 0..3 {
            let expected: f64 = row
                .physical_row_m2()
                .iter()
                .zip(&b.physical_modes()[0][mode])
                .map(|(p, u)| p * u)
                .sum();
            assert!((b.volume_shapes_m2_per_sqrt_kg()[mode] - expected).abs() < 1e-12);
            assert!((b.cavity_volume_shapes_m2_per_sqrt_kg()[mode] + expected).abs() < 1e-12);
        }
        let mut left_only = full.clone();
        left_only.specimen = full
            .specimen
            .with_region_material("part-1", "wood", &shell_material(700.0, 0.0), gate)
            .unwrap();
        let mut right_only = full.clone();
        right_only.specimen = full
            .specimen
            .with_region_material("part-0", "wood", &shell_material(700.0, 0.0), gate)
            .unwrap();
        let left = left_only
            .specimen
            .reduce(&left_only.reduction, &left_only.force_ports, gate)
            .unwrap();
        let right = right_only
            .specimen
            .reduce(&right_only.reduction, &right_only.force_ports, gate)
            .unwrap();
        assert_eq!(a.physical_modes(), left.physical_modes());
        assert_eq!(a.physical_modes(), right.physical_modes());
        for (whole, (l, r)) in a.damping().matrix_per_s().iter().zip(
            left.damping()
                .matrix_per_s()
                .iter()
                .zip(right.damping().matrix_per_s()),
        ) {
            assert!(
                (whole - l - r).abs() <= 1e-11 * whole.abs().max(l.abs()).max(r.abs()).max(1.0)
            );
        }
        assert!((0..3).any(|i| {
            (0..3).any(|j| i != j && a.damping().matrix_per_s()[3 * i + j].abs() > 1e-10)
        }));
    });
}

#[test]
fn g3_one_part_material_swap_changes_the_common_body_without_replacing_other_sources() {
    with_cx(|_, gate| {
        let original = body(gate);
        let mut changed = original.clone();
        changed.specimen = original
            .specimen
            .with_region_material("part-1", "wood", &shell_material(1400.0, 2000.0), gate)
            .unwrap();
        let a = original
            .specimen
            .reduce(&original.reduction, &original.force_ports, gate)
            .unwrap();
        let b = changed
            .specimen
            .reduce(&changed.reduction, &changed.force_ports, gate)
            .unwrap();
        assert_eq!(
            b.source().components[0].specimen.mass_kg(),
            a.source().components[0].specimen.mass_kg()
        );
        assert!(
            (b.source().components[1].specimen.mass_kg()
                / a.source().components[1].specimen.mass_kg()
                - 2.0)
                .abs()
                < 1e-12
        );
        for i in 0..2 {
            assert_eq!(
                a.source().components[i].specimen.mesh(),
                b.source().components[i].specimen.mesh()
            );
            assert_eq!(
                a.source().components[i].specimen.pressure_surfaces(),
                b.source().components[i].specimen.pressure_surfaces()
            );
        }
        let midpoint = |x: (f64, f64)| 0.5 * (x.0 + x.1);
        assert!(
            midpoint(b.eigenvalue_intervals()[0]) < 0.999 * midpoint(a.eigenvalue_intervals()[0])
        );
        assert_ne!(a.damping().matrix_per_s(), b.damping().matrix_per_s());
        assert_ne!(a.force_ports()[0], b.force_ports()[0]);
        assert_ne!(
            a.volume_shapes_m2_per_sqrt_kg(),
            b.volume_shapes_m2_per_sqrt_kg()
        );
    });
}

#[test]
fn g3_observer_only_part_selection_preserves_the_entire_assembled_mechanics() {
    with_cx(|cx, gate| {
        let original = body(gate);
        let mut quiet = original.clone();
        for part in &mut quiet.specimen.components {
            part.specimen = part
                .specimen
                .with_pressure_surfaces(surfaces(part.specimen.mesh().tris.len(), false), gate)
                .unwrap();
        }
        let mut a = build(cx, gate, original);
        let mut b = build(cx, gate, quiet);
        let ia = a.source().material_shell_assembly_image().unwrap();
        let ib = b.source().material_shell_assembly_image().unwrap();
        assert_eq!(ia.physical_modes(), ib.physical_modes());
        assert_eq!(ia.damping().matrix_per_s(), ib.damping().matrix_per_s());
        assert_eq!(
            ia.cavity_volume_shapes_m2_per_sqrt_kg(),
            ib.cavity_volume_shapes_m2_per_sqrt_kg()
        );
        assert!(ib.volume_shapes_m2_per_sqrt_kg().iter().all(|x| *x == 0.0));
        assert_eq!(a.source().body_connections(), b.source().body_connections());
        let zero = vec![0.0; a.system().mode_count()];
        let mut pressure_difference = 0.0_f64;
        for _ in 0..96 {
            let pa = a.step(&zero, drive()).unwrap().observer_pressure_pa;
            let pb = b.step(&zero, drive()).unwrap().observer_pressure_pa;
            pressure_difference = pressure_difference.max((pa - pb).abs());
        }
        assert!(pressure_difference > 1e-15);
        for (left, right) in a.system().components().iter().zip(b.system().components()) {
            assert_eq!(left.states(), right.states());
        }
    });
}

#[test]
fn g0_assembly_refuses_ambiguous_frames_slave_supports_budgets_and_missing_fluid_roles() {
    with_cx(|_, gate| {
        for case in 0..7 {
            let mut input = body(gate);
            match case {
                0 => input.reduction.max_total_nodes -= 1,
                1 => input.reduction.projection_budget.max_projection_terms = 1,
                2 => {
                    input.specimen.components[0].constraints[0].frame =
                        ConstraintFrame::Local([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
                }
                3 => input.specimen.components[1]
                    .constraints
                    .push(ShellConstraint {
                        node: input.specimen.edge_ties[0].b.nodes[0],
                        mask: [true; 6],
                        frame: ConstraintFrame::Global,
                    }),
                4 => input.force_ports[0].components.pop().map(|_| ()).unwrap(),
                5 => {
                    let old = &input.specimen.components[1];
                    let source = compile_material_shear_shell(
                        old.specimen.mesh().clone(),
                        old.specimen
                            .regions()
                            .iter()
                            .map(|r| r.input().clone())
                            .collect(),
                        old.assembly_budget,
                        gate,
                    )
                    .unwrap();
                    input.specimen.components[1].specimen = source;
                }
                _ => input.specimen.tie_budget.max_slave_dofs -= 1,
            }
            assert!(
                input
                    .specimen
                    .reduce(&input.reduction, &input.force_ports, gate)
                    .is_err(),
                "case {case}"
            );
        }
        let input = body(gate);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            input
                .specimen
                .reduce(&input.reduction, &input.force_ports, &cancelled)
                .is_err()
        );
    });
}

#[path = "assembly/glue.rs"]
mod glue;
