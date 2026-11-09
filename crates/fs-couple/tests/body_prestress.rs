//! Native-owner parity and analytic force/stability oracles for assembled preload.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::body_assembly::{
    BodyEquilibriumTolerance, BodyPrestressBudget, ThermalShellComponent,
    equilibrate_thermal_shell_body,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::sheet::SheetMatrix;
use fs_modal::{SliceOptions, WindowPlan};
use fs_plate::PlateError;
use fs_plate::shell::arching::ReferenceLoadState;
use fs_plate::shell::junction::{
    EdgeTie, FaceTie, GlueLine, RotationCoupling, TieBudget, TieKind, TieOffset, TieSide,
};
use fs_plate::shell::mitc::thermal::prestress::MembranePrestressLimits;
use fs_plate::shell::mitc::thermal::{FacetThermalStrain, ShellThermalStrain};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_plate::shell::{ShellMesh, ShellModel};
use fs_sparse::Csr;

const A: f64 = 1.0e6;
const UNLOADED: ReferenceLoadState = ReferenceLoadState::Unloaded;

fn section() -> ShearShellSection {
    ShearShellSection::new(
        SheetMatrix::orthotropic(A, A, 0.0, 0.5 * A).unwrap(),
        SheetMatrix::orthotropic(1.0, 1.0, 0.0, 0.5).unwrap(),
        [100.0, 0.0, 100.0],
        1.0,
        2.0,
        0.001,
        0.0,
    )
    .unwrap()
}
fn model(mesh: &ShellMesh, constraints: &[ShellConstraint]) -> MitcShellModel {
    let directors =
        NodalDirectors::new(mesh, vec![[0.0, 0.0, 1.0]; mesh.nodes.len()], 0.9).unwrap();
    assemble_mitc3p(
        mesh,
        &directors,
        &vec![section(); mesh.tris.len()],
        constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 32,
            max_triangles: 32,
        },
        1e-3,
        || false,
    )
    .unwrap()
}
fn fixed(node: usize) -> ShellConstraint {
    ShellConstraint {
        node,
        mask: [true; 6],
        frame: ConstraintFrame::Global,
    }
}
fn uniform(free: f64, count: usize) -> Vec<FacetThermalStrain> {
    vec![
        ShellThermalStrain::try_new([free, free, 0.0], [0.0; 3])
            .unwrap()
            .into();
        count
    ]
}
fn limits() -> MembranePrestressLimits {
    MembranePrestressLimits {
        maximum_displacement_gradient: 0.01,
        maximum_rotation_rad: 0.01,
    }
}
fn budget() -> BodyPrestressBudget {
    BodyPrestressBudget {
        maximum_components: 8,
        maximum_total_nodes: 64,
        maximum_component_free_coordinates: 384,
        ties: TieBudget {
            max_segments: 128,
            max_slave_dofs: 128,
            gap_tolerance_m: 1e-12,
            angle_tolerance_rad: 1e-10,
        },
    }
}
fn tolerance() -> BodyEquilibriumTolerance {
    BodyEquilibriumTolerance {
        maximum_residual_force_n: 1e-8,
        moment_length_m: 1.0,
    }
}
fn panel() -> ShellMesh {
    ShellMesh::new(
        vec![
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
            [0.5, 0.5, 0.],
        ],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
    )
    .unwrap()
}
fn one_free_panel(mesh: &ShellMesh) -> MitcShellModel {
    model(
        mesh,
        &(0..5)
            .map(|node| ShellConstraint {
                node,
                mask: [true, true, node != 4, true, true, true],
                frame: ConstraintFrame::Global,
            })
            .collect::<Vec<_>>(),
    )
}
fn entry(a: &Csr, i: usize, j: usize) -> f64 {
    a.row(i)
        .0
        .iter()
        .zip(a.row(i).1)
        .find_map(|(&c, &v)| (c == j).then_some(v))
        .unwrap_or(0.0)
}
fn near(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn matrices(a: &ShellModel, b: &ShellModel) {
    assert_eq!(a.free, b.free);
    for i in 0..a.free {
        for j in 0..a.free {
            near(entry(&a.k, i, j), entry(&b.k, i, j), 2e-11, 2e-9);
            near(entry(&a.m, i, j), entry(&b.m, i, j), 2e-12, 1e-13);
        }
    }
}
fn with_cx<R>(run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 64,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn window(lo: f64, hi: f64) -> WindowPlan {
    WindowPlan {
        windows: vec![(lo, hi)],
        threads: 1,
        maximum_windows: 1,
        maximum_bytes: 64 << 20,
        maximum_relative_halfwidth: 1e-7,
        coordinates: None,
    }
}

#[test]
fn one_part_matches_original_static_thermal_recovery_and_tangent() {
    let mesh = panel();
    let model = one_free_panel(&mesh);
    let thermal = model.prepare_thermal(&uniform(-1e-5, 4), || false).unwrap();
    let mut force = vec![0.0; 30];
    force[26] = 0.001;
    let original = thermal
        .solve_membrane_prestress(&force, limits(), || false)
        .unwrap();
    let body = equilibrate_thermal_shell_body(
        &[ThermalShellComponent {
            thermal: &thermal,
            dead_load_6n: &force,
            nodal_thickness_m: &[0.01; 5],
            reference_state: &UNLOADED,
        }],
        &[],
        &[],
        budget(),
        tolerance(),
        limits(),
        || false,
    )
    .unwrap();
    matrices(&body.assembly().model, original.pencil());
    for (&a, &b) in body.displacements().unwrap()[0]
        .iter()
        .zip(&original.reference_state().displacement_6n)
    {
        near(a, b, 2e-13, 1e-15);
    }
    for (a, b) in body.components()[0]
        .reference_state()
        .facets
        .iter()
        .zip(&original.reference_state().facets)
    {
        for i in 0..3 {
            near(
                a.membrane_force_n_m[i],
                b.membrane_force_n_m[i],
                2e-13,
                1e-12,
            );
        }
    }
    assert!(body.residual_force_n() < 1e-14);
    assert_eq!(body.assembly().model.m.row(0), model.pencil().m.row(0));
}

#[test]
fn analytic_tension_compression_and_negative_modes_remain_signed() {
    let mesh = panel();
    let model = one_free_panel(&mesh);
    let k0 = entry(&model.pencil().k, 0, 0);
    let mass = entry(&model.pencil().m, 0, 0);
    let critical = k0 / (4.0 * A);
    for ratio in [-0.5, 0.0, 0.5, 2.0] {
        let thermal = model
            .prepare_thermal(&uniform(ratio * critical, 4), || false)
            .unwrap();
        let body = equilibrate_thermal_shell_body(
            &[ThermalShellComponent {
                thermal: &thermal,
                dead_load_6n: &[0.0; 30],
                nodal_thickness_m: &[0.01; 5],
                reference_state: &UNLOADED,
            }],
            &[],
            &[],
            budget(),
            tolerance(),
            limits(),
            || false,
        )
        .unwrap();
        // Four quarter-area P1 facets have |grad(N_center)|=2, hence K_G=4*N.
        // N=-A*epsilon_thermal, so lambda=(1-ratio)*K0/M exactly for this panel.
        let expected = (1.0 - ratio) * k0 / mass;
        near(
            entry(&body.assembly().model.k, 0, 0) / mass,
            expected,
            5e-13,
            1e-10,
        );
        assert_eq!(body.assembly().model.m.row(0), model.pencil().m.row(0));
        let span = 0.25 * expected.abs();
        with_cx(|cx| {
            let modes = body
                .eigen_windows(
                    &window(expected - span, expected + span),
                    &SliceOptions {
                        mass_diagonal_equilibration: true,
                        ..SliceOptions::default()
                    },
                    cx,
                )
                .unwrap();
            assert_eq!(modes.spectrum().reports[0].modes.len(), 1);
            let actual = modes.spectrum().reports[0].modes[0].lambda;
            near(actual, expected, 2e-11, 1e-9);
            assert_eq!(actual.is_sign_negative(), ratio > 1.0);
            let mut load = vec![0.0; 30];
            load[26] = 3.0;
            let participation = modes.project_loads(&[&load]).unwrap()[0][0];
            near(
                participation,
                3.0 * modes.expand_mode(0, 0).unwrap()[0][26],
                1e-14,
                1e-14,
            );
        });
    }
}

#[test]
fn tied_native_panels_match_monolithic_equilibrium_tangent_and_virtual_work() {
    let mesh = ShellMesh::new(
        vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
        vec![[0, 1, 2], [1, 3, 2]],
    )
    .unwrap();
    let left = ShellMesh::new(mesh.nodes[..3].to_vec(), vec![[0, 1, 2]]).unwrap();
    let right = ShellMesh::new(
        vec![mesh.nodes[1], mesh.nodes[3], mesh.nodes[2]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let monolithic = model(&mesh, &[fixed(0), fixed(2)]);
    let lm = model(&left, &[fixed(0), fixed(2)]);
    let rm = model(&right, &[]);
    let whole = monolithic
        .prepare_thermal(&uniform(-1e-5, 2), || false)
        .unwrap();
    let lt = lm.prepare_thermal(&uniform(-1e-5, 1), || false).unwrap();
    let rt = rm.prepare_thermal(&uniform(-1e-5, 1), || false).unwrap();
    let mut load = vec![0.; 24];
    load[18] = 1.;
    load[20] = 0.001;
    let mut right_load = vec![0.; 18];
    right_load[6] = 1.;
    right_load[8] = 0.001;
    let original = whole
        .solve_membrane_prestress(&load, limits(), || false)
        .unwrap();
    let tie = EdgeTie {
        a: TieSide {
            component: 0,
            nodes: vec![1, 2],
        },
        b: TieSide {
            component: 1,
            nodes: vec![0, 2],
        },
        offset: TieOffset::Midsurface,
        rotations: RotationCoupling::Full,
        kind: TieKind::Rigid,
    };
    let body = equilibrate_thermal_shell_body(
        &[
            ThermalShellComponent {
                thermal: &lt,
                dead_load_6n: &[0.; 18],
                nodal_thickness_m: &[0.01; 3],
                reference_state: &UNLOADED,
            },
            ThermalShellComponent {
                thermal: &rt,
                dead_load_6n: &right_load,
                nodal_thickness_m: &[0.01; 3],
                reference_state: &UNLOADED,
            },
        ],
        &[tie],
        &[],
        budget(),
        tolerance(),
        limits(),
        || false,
    )
    .unwrap();
    matrices(&body.assembly().model, original.pencil());
    let physical = body.displacements().unwrap();
    for (component, map) in [(0, vec![0, 1, 2]), (1, vec![1, 3, 2])] {
        for (node, original_node) in map.into_iter().enumerate() {
            for d in 0..6 {
                near(
                    physical[component][6 * node + d],
                    original.reference_state().displacement_6n[6 * original_node + d],
                    2e-9,
                    1e-12,
                );
            }
        }
    }
    assert!(
        body.component_free_residuals()
            .iter()
            .flatten()
            .any(|x| x.abs() > 0.1)
    );
    assert!(body.residual_force_n() < 1e-8);
    let velocity: Vec<_> = (0..body.assembly().model.free)
        .map(|i| (i as f64 - 4.) * 0.017)
        .collect();
    let world = body.assembly().expand_components(&velocity).unwrap();
    let forces = [
        (0..18).map(|i| (i as f64 - 5.) * 0.2).collect::<Vec<_>>(),
        (0..18)
            .map(|i| (i as f64 + 2.) * (-0.13))
            .collect::<Vec<_>>(),
    ];
    let generalized = body
        .assembly()
        .project_loads(&[&forces[0], &forces[1]])
        .unwrap();
    let reduced_power: f64 = velocity.iter().zip(generalized).map(|(u, f)| u * f).sum();
    let physical_power: f64 = world
        .iter()
        .zip(&forces)
        .flat_map(|(u, f)| u.iter().zip(f))
        .map(|(u, f)| u * f)
        .sum();
    near(reduced_power, physical_power, 1e-12, 1e-13);
}

#[test]
fn incompatible_reference_and_supports_refuse_without_changing_sources() {
    let mesh = panel();
    let model = one_free_panel(&mesh);
    let thermal = model.prepare_thermal(&uniform(1e-5, 4), || false).unwrap();
    let loaded = ReferenceLoadState::Loaded {
        description: "measured under string tension".into(),
    };
    let mut component = ThermalShellComponent {
        thermal: &thermal,
        dead_load_6n: &[0.; 30],
        nodal_thickness_m: &[0.01; 5],
        reference_state: &loaded,
    };
    assert!(
        equilibrate_thermal_shell_body(
            &[component],
            &[],
            &[],
            budget(),
            tolerance(),
            limits(),
            || false
        )
        .is_err()
    );
    component.reference_state = &UNLOADED;
    assert!(matches!(
        equilibrate_thermal_shell_body(
            &[component],
            &[],
            &[],
            budget(),
            tolerance(),
            limits(),
            || true
        ),
        Err(PlateError::Cancelled)
    ));
    let mut cap = budget();
    cap.maximum_total_nodes = 4;
    assert!(matches!(
        equilibrate_thermal_shell_body(&[component], &[], &[], cap, tolerance(), limits(), || {
            false
        }),
        Err(PlateError::Budget { .. })
    ));
    let before = model.pencil().k.clone();
    let _ = equilibrate_thermal_shell_body(
        &[component],
        &[],
        &[],
        budget(),
        tolerance(),
        limits(),
        || false,
    )
    .unwrap();
    for row in 0..before.nrows() {
        assert_eq!(before.row(row), model.pencil().k.row(row));
    }
}

#[test]
fn compliant_interface_is_retained_once_in_equilibrium_and_tangent() {
    let mesh = panel();
    let model = one_free_panel(&mesh);
    let thermal = model.prepare_thermal(&uniform(-1e-5, 4), || false).unwrap();
    let mut load = vec![0.0; 30];
    load[26] = 0.001;
    let parts = [
        ThermalShellComponent {
            thermal: &thermal,
            dead_load_6n: &load,
            nodal_thickness_m: &[0.01; 5],
            reference_state: &UNLOADED,
        },
        ThermalShellComponent {
            thermal: &thermal,
            dead_load_6n: &[0.0; 30],
            nodal_thickness_m: &[0.01; 5],
            reference_state: &UNLOADED,
        },
    ];
    let face = FaceTie {
        master_component: 0,
        master_tris: vec![0, 1, 2, 3],
        slave_component: 1,
        slave_tris: vec![0, 1, 2, 3],
        offset: TieOffset::Midsurface,
        rotations: RotationCoupling::Full,
        kind: TieKind::Compliant(GlueLine {
            normal_n_m3: 100.0,
            shear_n_m3: [100.0; 2],
            viscosity_n_s_m3: Some([0.1; 3]),
            width_m: 0.01,
        }),
    };
    let body = equilibrate_thermal_shell_body(
        &parts,
        &[],
        &[face],
        budget(),
        tolerance(),
        limits(),
        || false,
    )
    .unwrap();
    let assembly = body.assembly();
    assert_eq!(assembly.model.free, 2);
    let k0 = entry(&model.pencil().k, 0, 0);
    let kg = entry(&assembly.glue_stiffness, 0, 0);
    assert!(kg > 0.0);
    near(entry(&assembly.glue_stiffness, 0, 1), -kg, 1e-13, 1e-12);
    for i in 0..2 {
        for j in 0..2 {
            let bare = if i == j { k0 } else { 0.0 };
            let geometric = if i == j { 40.0 } else { 0.0 };
            near(
                entry(body.elastic_stiffness(), i, j),
                bare + entry(&assembly.glue_stiffness, i, j),
                1e-13,
                1e-12,
            );
            near(
                entry(&assembly.model.k, i, j),
                bare + geometric + entry(&assembly.glue_stiffness, i, j),
                1e-12,
                1e-11,
            );
        }
    }
    let q = body.coordinates();
    // Two equal grounded elastic cells connected by their physical glue.
    near(
        q[0],
        0.001 * (k0 + kg) / (k0 * (k0 + 2.0 * kg)),
        1e-12,
        1e-15,
    );
    near(q[1], 0.001 * kg / (k0 * (k0 + 2.0 * kg)), 1e-12, 1e-15);
    assert!(assembly.glue_damping.is_some());
    assert!(body.residual_force_n() < 1e-12);
}

#[test]
fn oblique_or_supported_slave_constraints_and_unsupported_body_refuse() {
    let mesh = panel();
    let model = one_free_panel(&mesh);
    let thermal = model.prepare_thermal(&uniform(0.0, 4), || false).unwrap();
    let part = ThermalShellComponent {
        thermal: &thermal,
        dead_load_6n: &[0.0; 30],
        nodal_thickness_m: &[0.01; 5],
        reference_state: &UNLOADED,
    };
    let tie = EdgeTie {
        a: TieSide {
            component: 0,
            nodes: vec![0, 1],
        },
        b: TieSide {
            component: 1,
            nodes: vec![0, 1],
        },
        offset: TieOffset::Midsurface,
        rotations: RotationCoupling::Full,
        kind: TieKind::Rigid,
    };
    assert!(matches!(
        equilibrate_thermal_shell_body(
            &[part, part],
            &[tie],
            &[],
            budget(),
            tolerance(),
            limits(),
            || false
        ),
        Err(PlateError::BadConstraint {
            what: "supported rigid-tie slaves require native constraint composition"
        })
    ));
    let oblique = model_with_oblique_supports(&mesh);
    let field = oblique.prepare_thermal(&uniform(0.0, 4), || false).unwrap();
    assert!(matches!(
        equilibrate_thermal_shell_body(
            &[ThermalShellComponent {
                thermal: &field,
                ..part
            }],
            &[],
            &[],
            budget(),
            tolerance(),
            limits(),
            || false
        ),
        Err(PlateError::BadConstraint {
            what: "assembled preload requires original global support coordinates"
        })
    ));
    let free_body = crate::model(&mesh, &[]);
    let field = free_body
        .prepare_thermal(&uniform(0.0, 4), || false)
        .unwrap();
    let mut unsupported = vec![0.0; 30];
    unsupported[26] = 1.0;
    assert!(
        equilibrate_thermal_shell_body(
            &[ThermalShellComponent {
                thermal: &field,
                dead_load_6n: &unsupported,
                ..part
            }],
            &[],
            &[],
            budget(),
            tolerance(),
            limits(),
            || false
        )
        .is_err()
    );
}

fn model_with_oblique_supports(mesh: &ShellMesh) -> MitcShellModel {
    model(
        mesh,
        &(0..5)
            .map(|node| ShellConstraint {
                node,
                mask: [true, true, node != 4, true, true, true],
                frame: ConstraintFrame::Local([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
            })
            .collect::<Vec<_>>(),
    )
}
