//! G1/G3: native tensor-Prony shell coordinates, two wet faces and one pressure jump.
//! The structural mesh is fixed; the small acoustic ladder is not a continuum
//! shell-accuracy or edge-singularity convergence certificate.

use fs_bem::helmholtz::Medium;
use fs_bem::helmholtz::bc::{Admittance, SurfaceBc};
use fs_bem::near_field::Options;
use fs_bem::reciprocal_radiation::open_screen::{OpenOptions, OpenPanels, ScreenPortLimits};
use fs_couple::tensor_shell_radiation::TensorShellRadiationError;
use fs_couple::tensor_shell_radiation::open::{
    OpenTensorShellRadiation, OpenTensorShellRadiationControls,
    OpenTensorShellRadiationFrequencyControls,
};
use fs_exec::{CancelGate, Cx};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicPortControls, HarmonicResponseControls};
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::viscosity::harmonic::{
    TensorPronyShell, TensorPronyShellSection, TensorShellHarmonicError,
};
use fs_plate::shell::{CstDktMass, ShellMesh, ShellSupport};

const OMEGA: f64 = 1000.0;

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 47,
                    kernel_id: 59,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn source(height: f64, cx: &Cx<'_>) -> TensorPronyShell {
    // Explicit numerical polymer-like material, not a sourced material claim.
    let c: Mandel6 = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            1e6 * if i < 3 && j < 3 {
                4.0 + if i == j { 6.0 } else { 0.0 }
            } else if i == j {
                6.0
            } else {
                0.0
            }
        })
    });
    let scale = |s: f64| c.map(|r| r.map(|x| x * s));
    let law = TensorProny::new_mandel(
        c,
        vec![
            PronyBranch::new(scale(0.4), 0.001).unwrap(),
            PronyBranch::new(scale(0.2), 0.02).unwrap(),
        ],
        0.02,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(scale(1e-5))
    .unwrap()
    .with_band((0.0, 1500.0))
    .unwrap();
    let section = TensorPronyShellSection::new(law, 0.01, 1200.0).unwrap();
    let mesh = ShellMesh::new(
        vec![
            [-0.2, -0.15, 0.0],
            [0.2, -0.15, 0.0],
            [0.2, 0.15, 0.0],
            [-0.2, 0.15, 0.0],
            [0.0, 0.0, height],
        ],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
    )
    .unwrap();
    TensorPronyShell::assemble(
        mesh,
        vec![section; 4],
        &[0, 1, 2, 3],
        ShellSupport::Clamped,
        CstDktMass::Consistent,
        0.0,
        ShellViscosityBudget {
            max_elements: 4,
            max_full_dofs: 30,
            max_dofs: 6,
            max_modes: 6,
            max_work: 20_000_000,
        },
        cx,
    )
    .unwrap()
}

fn medium() -> Medium {
    Medium {
        density: 1000.0,
        sound_speed: 1480.0,
    }
}
fn geometry_controls(level: u32) -> OpenTensorShellRadiationControls {
    OpenTensorShellRadiationControls {
        acoustic_refinements: level,
        maximum_panels: 64,
        maximum_coordinates: 6,
        maximum_setup_bytes: 1 << 20,
        maximum_setup_work: 10_000_000,
    }
}
fn boundary_options() -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 0,
        triangle_order: 8,
        angular_refinement: 0,
        maximum_unknowns: 64,
        maximum_cache_bytes: 16 << 20,
        maximum_work: 2_000_000_000,
        workers: 1,
    }
}
fn frequency_controls() -> OpenTensorShellRadiationFrequencyControls {
    OpenTensorShellRadiationFrequencyControls {
        mechanical: HarmonicResponseControls {
            maximum_coordinates: 6,
            maximum_payload_bytes: 16 << 20,
            maximum_relative_residual: 2e-11,
        },
        ports: HarmonicPortControls {
            maximum_ports: 6,
            maximum_projection_terms: 1_000_000,
            maximum_payload_bytes: 4 << 20,
        },
        boundary: ScreenPortLimits {
            maximum_ports: 6,
            maximum_matrix_entries: 36,
            maximum_jump_entries: 6 * 64,
            maximum_work: 2_000_000_000,
        },
        maximum_boundary_payload_bytes: 16 << 20,
        maximum_radiated_power_work: 512_000_000,
        power_relative_tolerance: 1e-8,
    }
}
fn force() -> [C64; 6] {
    [
        C64::new(1e-4, -2e-5),
        C64::new(-3e-5, 1e-5),
        C64::new(2e-5, 4e-5),
        C64::new(1e-6, -2e-7),
        C64::new(3e-7, 1e-7),
        C64::new(-5e-7, 1e-7),
    ]
}
fn near(a: C64, b: C64) {
    assert!((a - b).abs() <= 3e-8 * b.abs().max(1e-23), "{a:?} != {b:?}");
}
fn near_power(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-8 * b.abs().max(1e-26), "{a} != {b}");
}

// Independent closed-form center shape on each of the four original facets.
// It remains this shape on the shallow tent: height changes facet normals.
fn center_shape(point: [f64; 3]) -> f64 {
    1.0 - (point[0].abs() / 0.2).max(point[1].abs() / 0.15)
}
fn area_normal(points: [[f64; 3]; 3]) -> (f64, [f64; 3]) {
    let a: [f64; 3] = core::array::from_fn(|d| points[1][d] - points[0][d]);
    let b: [f64; 3] = core::array::from_fn(|d| points[2][d] - points[0][d]);
    let n = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let length = n.iter().map(|x| x * x).sum::<f64>().sqrt();
    (0.5 * length, n.map(|x| x / length))
}
// Degree-two quadrature integrates P1 pressure times the ORIGINAL P1 velocity
// exactly. This does not read adapter rows or the production triangle mass map.
fn projected_effort(surface: &OpenPanels, jump: &[C64]) -> [C64; 6] {
    let mut effort = [C64::ZERO; 6];
    for triangle in surface.triangles() {
        let points = triangle.map(|v| surface.vertices()[v]);
        let (area, normal) = area_normal(points);
        for q in [
            [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
            [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
            [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
        ] {
            let position = core::array::from_fn(|d| (0..3).map(|a| q[a] * points[a][d]).sum());
            let pressure = (0..3).fold(C64::ZERO, |sum, a| sum + jump[triangle[a]].scale(q[a]));
            for d in 0..3 {
                effort[d] =
                    effort[d] + pressure.scale(area / 3.0 * center_shape(position) * normal[d]);
            }
        }
    }
    effort
}

#[test]
fn full_open_shell_preserves_native_material_coordinates_and_two_sided_work() {
    with_cx(|cx, _| {
        let shell = source(0.0, cx);
        assert_eq!(shell.model().free, 6);
        assert_eq!(shell.mesh().tris.len(), 4);
        for d in 0..6 {
            assert_eq!(shell.model().dof_map[24 + d], Some(d));
        }
        let geometry =
            OpenTensorShellRadiation::new(&shell, medium(), geometry_controls(1), cx).unwrap();
        assert!(core::ptr::eq(geometry.shell(), &shell));
        assert_eq!(geometry.surface().triangles().len(), 16);
        assert_eq!(geometry.normal_velocity_rows().len(), 6);
        for d in [0, 1, 3, 4, 5] {
            assert!(
                geometry.normal_velocity_rows()[d]
                    .iter()
                    .flatten()
                    .all(|&v| v == 0.0)
            );
        }
        for (ids, values) in geometry
            .surface()
            .triangles()
            .iter()
            .zip(&geometry.normal_velocity_rows()[2])
        {
            for a in 0..3 {
                assert!(
                    (values[a] - center_shape(geometry.surface().vertices()[ids[a]])).abs() < 2e-15
                );
            }
        }
        let prepared = geometry.prepare(2000.0, boundary_options(), cx).unwrap();
        let loaded = prepared
            .prepare_frequency(OMEGA, frequency_controls(), cx)
            .unwrap();
        assert!(core::ptr::eq(loaded.shell_preparation().source(), &shell));
        assert_eq!(loaded.radiation().ports(), 6);
        assert_eq!(loaded.radiation().solutions().len(), 6);
        let result = loaded.solve(&force(), cx).unwrap();
        let motion = result.motion();
        assert_eq!(motion.displacement.len(), 6);
        assert!(motion.displacement[3..].iter().all(|x| x.abs() > 0.0));
        assert!(motion.average_viscous_power_w > 0.0);
        assert!(result.radiated_power_w() > 0.0);
        near_power(
            motion.average_input_power_w,
            motion.average_viscous_power_w + result.radiated_power_w(),
        );
        assert!(result.relative_coupled_power_defect() < 1e-8);
        assert!(result.relative_port_work_defect() < 1e-8);
        assert!(result.relative_flux_defect() < 1e-8);

        // Independently invoke the source shell's complete dispersive pencil.
        let mut identity = [0.0; 36];
        for d in 0..6 {
            identity[6 * d + d] = 1.0;
        }
        let native = shell
            .prepare_harmonic_response_with_port_impedance(
                OMEGA,
                &identity,
                6,
                loaded.impedance(),
                frequency_controls().mechanical,
                frequency_controls().ports,
                cx,
            )
            .unwrap()
            .solve(&force(), cx)
            .unwrap();
        for (&a, &b) in motion.displacement.iter().zip(&native.displacement) {
            near(a, b);
        }
        near_power(
            motion.average_viscous_power_w,
            native.average_viscous_power_w,
        );
        let bare = shell
            .prepare_harmonic_response(OMEGA, frequency_controls().mechanical, cx)
            .unwrap()
            .solve(&force(), cx)
            .unwrap();
        assert!(
            (motion.displacement[2] - bare.displacement[2]).abs()
                > 0.01 * bare.displacement[2].abs()
        );
        for d in [0, 1, 3, 4, 5] {
            near(motion.displacement[d], bare.displacement[d]);
        }

        let jump = result.pressure_jump(cx).unwrap();
        let effort = projected_effort(geometry.surface(), &jump);
        let port = motion.port_load.as_ref().unwrap();
        for (&a, &b) in effort.iter().zip(&port.effort) {
            near(a, b);
        }
        let work = 0.5
            * motion
                .velocity
                .iter()
                .zip(effort)
                .map(|(&v, f)| (v.conj() * f).re)
                .sum::<f64>();
        near_power(work, result.jump_power_w());
        near_power(work, result.radiated_power_w());
        // One direct boundary solve driven by the same physical center motion,
        // constructed without the adapter's source rows or column combination.
        let drive: Vec<_> = geometry
            .surface()
            .vertices()
            .iter()
            .map(|&p| motion.velocity[2].scale(center_shape(p)))
            .collect();
        let direct = prepared
            .boundary()
            .solve(
                OMEGA / medium().sound_speed,
                medium(),
                SurfaceBc {
                    drive_velocity: &drive,
                    admittance: Admittance::Rigid,
                },
                || false,
            )
            .unwrap();
        for (&a, &b) in jump.iter().zip(direct.pressure_jump_pa()) {
            near(a, b);
        }
        near_power(work, direct.driving_power_w);
        let receivers = [[0.04, -0.03, 0.4], [0.04, -0.03, -0.4]];
        let pressure = result
            .pressure_at(&receivers, 0.01, Options::default(), cx)
            .unwrap();
        let reference = direct
            .pressure_at(&receivers, 0.01, Options::default(), || false)
            .unwrap();
        for (&a, &b) in pressure.iter().zip(&reference) {
            near(a, b);
        }
        assert!(pressure[0].abs() > 0.0);
        near(pressure[0], -pressure[1]);
        near(
            result.far_field([0.0, 0.0, 1.0], cx).unwrap(),
            direct.far_field([0.0, 0.0, 1.0]).unwrap(),
        );
        near(
            result.far_field([0.0, 0.0, 1.0], cx).unwrap(),
            -result.far_field([0.0, 0.0, -1.0], cx).unwrap(),
        );
    });
}

#[test]
fn folded_shell_keeps_facet_normal_jumps_and_reuses_coherent_sources_for_forces() {
    with_cx(|cx, _| {
        let shell = source(0.02, cx);
        let geometry =
            OpenTensorShellRadiation::new(&shell, medium(), geometry_controls(1), cx).unwrap();
        for (panel, ids) in geometry.surface().triangles().iter().enumerate() {
            let points = ids.map(|v| geometry.surface().vertices()[v]);
            let (_, normal) = area_normal(points);
            for d in 0..3 {
                for a in 0..3 {
                    assert!(
                        (geometry.normal_velocity_rows()[d][panel][a]
                            - normal[d] * center_shape(points[a]))
                        .abs()
                            < 2e-15
                    );
                }
            }
        }
        assert!(
            geometry.normal_velocity_rows()[0]
                .iter()
                .flatten()
                .any(|&v| v > 0.0)
        );
        assert!(
            geometry.normal_velocity_rows()[0]
                .iter()
                .flatten()
                .any(|&v| v < 0.0)
        );
        assert!(
            geometry.normal_velocity_rows()[3..]
                .iter()
                .flatten()
                .flatten()
                .all(|&v| v == 0.0)
        );
        let prepared = geometry.prepare(2000.0, boundary_options(), cx).unwrap();
        let loaded = prepared
            .prepare_frequency(OMEGA, frequency_controls(), cx)
            .unwrap();
        let original = loaded.impedance().to_vec();
        let first = loaded.solve(&force(), cx).unwrap();
        let gain = C64::new(-0.7, 0.4);
        let second = loaded.solve(&force().map(|v| v * gain), cx).unwrap();
        for (&a, &b) in second
            .motion()
            .displacement
            .iter()
            .zip(&first.motion().displacement)
        {
            near(a, b * gain);
        }
        near_power(
            second.radiated_power_w(),
            first.radiated_power_w() * gain.norm_sq(),
        );
        near_power(
            second.motion().average_viscous_power_w,
            first.motion().average_viscous_power_w * gain.norm_sq(),
        );
        let jump = first.pressure_jump(cx).unwrap();
        for (&a, &b) in second.pressure_jump(cx).unwrap().iter().zip(&jump) {
            near(a, b * gain);
        }
        let effort = projected_effort(geometry.surface(), &jump);
        for (&a, &b) in effort
            .iter()
            .zip(&first.motion().port_load.as_ref().unwrap().effort)
        {
            near(a, b);
        }
        let receivers = [[0.04, -0.03, 0.4], [0.04, -0.03, -0.4]];
        let a = first
            .pressure_at(&receivers, 0.01, Options::default(), cx)
            .unwrap();
        let b = second
            .pressure_at(&receivers, 0.01, Options::default(), cx)
            .unwrap();
        for (&a, &b) in b.iter().zip(&a) {
            near(a, b * gain);
        }
        assert_eq!(loaded.impedance(), original.as_slice());
        let zero = loaded.solve(&[C64::ZERO; 6], cx).unwrap();
        assert_eq!(zero.radiated_power_w(), 0.0);
        assert_eq!(zero.jump_power_w(), 0.0);
        assert_eq!(zero.relative_coupled_power_defect(), 0.0);
        assert_eq!(zero.relative_port_work_defect(), 0.0);
        assert_eq!(zero.relative_flux_defect(), 0.0);
    });
}

#[test]
fn acoustic_refinement_changes_the_load_without_changing_original_shell_coordinates() {
    with_cx(|cx, _| {
        let shell = source(0.0, cx);
        let original = shell.mesh().clone();
        let mut impedances = Vec::new();
        for level in 0..=2 {
            let geometry =
                OpenTensorShellRadiation::new(&shell, medium(), geometry_controls(level), cx)
                    .unwrap();
            assert_eq!(geometry.surface().triangles().len(), 4 * 4usize.pow(level));
            assert_eq!(geometry.normal_velocity_rows().len(), 6);
            let prepared = geometry.prepare(2000.0, boundary_options(), cx).unwrap();
            let loaded = prepared
                .prepare_frequency(OMEGA, frequency_controls(), cx)
                .unwrap();
            let result = loaded.solve(&force(), cx).unwrap();
            near_power(
                result.motion().average_input_power_w,
                result.motion().average_viscous_power_w + result.radiated_power_w(),
            );
            let z = loaded.impedance()[2 * 6 + 2];
            assert!(
                z.re > 0.0 && z.im < 0.0,
                "radiation resistance and added mass: {z:?}"
            );
            impedances.push(z);
        }
        let coarse_change = (impedances[1] - impedances[0]).abs();
        let fine_change = (impedances[2] - impedances[1]).abs();
        assert!(
            coarse_change > 1e-4 * impedances[2].abs(),
            "acoustic refinement must be live"
        );
        assert!(
            fine_change < coarse_change,
            "three-rung load ladder: {impedances:?}"
        );
        assert_eq!(shell.mesh(), &original);
        assert_eq!(shell.model().free, 6);
    });
}

#[test]
fn parallel_screen_payload_cannot_omit_native_angular_or_pair_storage() {
    with_cx(|cx, _| {
        let shell = source(0.0, cx);
        let geometry =
            OpenTensorShellRadiation::new(&shell, medium(), geometry_controls(2), cx).unwrap();
        let panels = geometry.surface().triangles().len();
        assert_eq!(panels, 64);
        for (workers, angular_refinement) in [(2, 2), (16, 0)] {
            let prepared = geometry
                .prepare(
                    2000.0,
                    OpenOptions {
                        workers,
                        angular_refinement,
                        ..boundary_options()
                    },
                    cx,
                )
                .unwrap();
            let pressure_unknowns = prepared.boundary().unknowns();
            assert_eq!(pressure_unknowns, 25);
            // Actual allocations in the native parallel owner, independently
            // of the adapter's new conservative bound. Neither includes H,
            // factors, sources or bookkeeping, so neither suffices for the
            // whole frequency. The former 168,416-byte admission accepted both.
            let one_native_allocation = if angular_refinement == 2 {
                // One complex vector per angular direction, all retained.
                16 * pressure_unknowns * (64 * 4usize.pow(angular_refinement))
            } else {
                // With 16 workers all 64 target rows share the first block:
                // ten complex values per upper-triangular facet pair.
                160 * panels * (panels + 1) / 2
            };
            match prepared.prepare_frequency(
                OMEGA,
                OpenTensorShellRadiationFrequencyControls {
                    maximum_boundary_payload_bytes: one_native_allocation,
                    maximum_radiated_power_work: 4_000_000_000,
                    ..frequency_controls()
                },
                cx,
            ) {
                Err(TensorShellRadiationError::Budget {
                    quantity,
                    required,
                    maximum,
                }) => {
                    assert_eq!(quantity, "open shell frequency payload");
                    assert_eq!(maximum, one_native_allocation as u128);
                    assert!(required > maximum);
                }
                Err(error) => panic!("unexpected refusal: {error:?}"),
                Ok(_) => panic!("frequency admitted less than its native parallel workspace"),
            }
        }
    });
}

#[test]
fn full_source_caps_material_domains_receivers_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let shell = source(0.0, cx);
        let original = shell.mesh().clone();
        assert!(matches!(
            OpenTensorShellRadiation::new(
                &shell,
                medium(),
                OpenTensorShellRadiationControls {
                    maximum_coordinates: 5,
                    ..geometry_controls(1)
                },
                cx
            ),
            Err(TensorShellRadiationError::Budget { .. })
        ));
        assert!(
            OpenTensorShellRadiation::new(
                &shell,
                medium(),
                OpenTensorShellRadiationControls {
                    maximum_panels: 15,
                    ..geometry_controls(1)
                },
                cx
            )
            .is_err()
        );
        assert!(
            OpenTensorShellRadiation::new(
                &shell,
                medium(),
                OpenTensorShellRadiationControls {
                    maximum_setup_bytes: 1,
                    ..geometry_controls(1)
                },
                cx
            )
            .is_err()
        );
        assert!(
            OpenTensorShellRadiation::new(
                &shell,
                Medium {
                    density: f64::NAN,
                    ..medium()
                },
                geometry_controls(1),
                cx
            )
            .is_err()
        );
        let geometry =
            OpenTensorShellRadiation::new(&shell, medium(), geometry_controls(1), cx).unwrap();
        let prepared = geometry.prepare(2000.0, boundary_options(), cx).unwrap();
        assert!(
            prepared
                .prepare_frequency(2001.0, frequency_controls(), cx)
                .is_err()
        );
        assert!(
            prepared
                .prepare_frequency(
                    OMEGA,
                    OpenTensorShellRadiationFrequencyControls {
                        maximum_boundary_payload_bytes: 1,
                        ..frequency_controls()
                    },
                    cx
                )
                .is_err()
        );
        // Five dark coordinates are still original physical source columns.
        // Neither the source-count nor full 6x6 matrix cap may drop them.
        for boundary in [
            ScreenPortLimits {
                maximum_ports: 5,
                ..frequency_controls().boundary
            },
            ScreenPortLimits {
                maximum_matrix_entries: 35,
                ..frequency_controls().boundary
            },
        ] {
            assert!(
                prepared
                    .prepare_frequency(
                        OMEGA,
                        OpenTensorShellRadiationFrequencyControls {
                            boundary,
                            ..frequency_controls()
                        },
                        cx
                    )
                    .is_err()
            );
        }
        assert!(matches!(
            prepared.prepare_frequency(1600.0, frequency_controls(), cx),
            Err(TensorShellRadiationError::Shell(
                TensorShellHarmonicError::Material { .. }
            ))
        ));
        let loaded = prepared
            .prepare_frequency(OMEGA, frequency_controls(), cx)
            .unwrap();
        assert!(matches!(
            loaded.solve(&force().map(|v| v.scale(1e12)), cx),
            Err(TensorShellRadiationError::Shell(
                TensorShellHarmonicError::Material { .. }
            ))
        ));
        let response = loaded.solve(&force(), cx).unwrap();
        assert!(
            response
                .pressure_at(&[[0.0, 0.0, 0.0]], 0.001, Options::default(), cx)
                .is_err()
        );
        assert!(response.far_field([0.0; 3], cx).is_err());
        assert_eq!(shell.mesh(), &original);
        gate.request();
        assert!(matches!(
            loaded.solve(&force(), cx),
            Err(TensorShellRadiationError::Cancelled)
        ));
        assert!(matches!(
            response.pressure_jump(cx),
            Err(TensorShellRadiationError::Cancelled)
        ));
        assert!(matches!(
            geometry.prepare(2000.0, boundary_options(), cx),
            Err(TensorShellRadiationError::Cancelled)
        ));
    });
}
