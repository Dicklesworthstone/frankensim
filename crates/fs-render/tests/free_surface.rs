use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_lbm::d3q19::freesurface3::{ContactModel3, FreeSurface3, FreeSurfaceInitialCell3};
use fs_render::charts::{Ray, TriMesh};
use fs_render::dielectric::{BeerLambertAbsorption, DielectricEvent, sample_smooth_dielectric};
use fs_render::free_surface::{
    FreeSurfaceMeshError, FreeSurfaceMeshOptions, FreeSurfaceVolumeOptions, free_surface_mesh,
    volume_matched_free_surface_mesh,
};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn block(cx: &Cx<'_>, top_fill: f64, density: f64) -> FreeSurface3 {
    let mut cells = Vec::new();
    for z in 3..9 {
        for y in 3..9 {
            for x in 3..9 {
                cells.push(FreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: if z == 8 { top_fill } else { 1.0 },
                    density,
                    velocity: [0.0; 3],
                });
            }
        }
    }
    FreeSurface3::from_cells(
        [12; 3],
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        &[],
        &cells,
        cx,
    )
    .unwrap()
}

fn options() -> FreeSurfaceMeshOptions {
    FreeSurfaceMeshOptions {
        origin_m: [0.0; 3],
        cell_size_m: 0.01,
        max_cells: 12 * 12 * 12,
        max_triangles: 10_000,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

#[test]
fn g1_actual_fraction_moves_the_optical_boundary_without_changing_liquid_state() {
    with_cx(|cx, _| {
        for (fill, height) in [(0.25, 7.5 + 2.0 / 3.0), (0.75, 8.5 + 1.0 / 3.0)] {
            let state = block(cx, fill, 1.0);
            let mass = state.ledger_mass();
            let bits = state.grid().state_bits();
            let mesh = free_surface_mesh(&state, options(), cx).unwrap();
            let replay = free_surface_mesh(&state, options(), cx).unwrap();
            assert_eq!(mesh.vertices, replay.vertices);
            assert_eq!(mesh.triangles, replay.triangles);
            assert_eq!(mesh.bvh_fingerprint(), replay.bvh_fingerprint());
            assert_eq!(state.ledger_mass(), mass);
            assert_eq!(state.grid().state_bits(), bits);
            let hit = mesh
                .intersect(&Ray {
                    origin: Point3::new(0.056, 0.057, 0.11),
                    dir: Vec3::new(0.0, 0.0, -1.0),
                })
                .unwrap();
            close(hit.point.z, height * 0.01, 2e-16);
            assert_eq!(hit.normal.unwrap(), Vec3::new(0.0, 0.0, 1.0));
            // Liquid density controls hydrodynamics, not the observed fill geometry.
            let scaled = free_surface_mesh(&block(cx, fill, 2.0), options(), cx).unwrap();
            assert_eq!(mesh.vertices, scaled.vertices);
            assert_eq!(mesh.triangles, scaled.triangles);
        }
    });
}

fn check_native_liquid_transport(mesh: &TriMesh, cx: &Cx<'_>, lower: f64, upper: f64) -> f64 {
    let incoming = Ray {
        origin: Point3::new(0.045, 0.056, 0.0),
        dir: Vec3::new(0.2, 0.0, 0.96_f64.sqrt()),
    };
    let first = mesh.intersect_with_cx(cx, &incoming).unwrap().unwrap();
    close(first.point.z, lower, 2e-16);
    assert_eq!(first.normal.unwrap(), Vec3::new(0.0, 0.0, -1.0));
    let entry = sample_smooth_dielectric(
        first.normal.unwrap(),
        incoming.dir.scale(-1.0),
        1.0,
        1.5,
        0.99,
    )
    .unwrap();
    assert_eq!(entry.event, DielectricEvent::Transmission);
    close(1.5 * entry.direction.x, incoming.dir.x, 3e-16);
    let offset = 1e-8;
    let internal = Ray {
        origin: first.point.offset(entry.direction.scale(offset)),
        dir: entry.direction,
    };
    let last = mesh.intersect_with_cx(cx, &internal).unwrap().unwrap();
    close(last.point.z, upper, 2e-16);
    assert_eq!(last.normal.unwrap(), Vec3::new(0.0, 0.0, 1.0));
    let distance = last.t + offset;
    close(
        distance,
        (upper - lower) / (1.0_f64 - (0.2_f64 / 1.5).powi(2)).sqrt(),
        3e-16,
    );
    let exit = sample_smooth_dielectric(
        last.normal.unwrap().scale(-1.0),
        entry.direction.scale(-1.0),
        1.5,
        1.0,
        0.99,
    )
    .unwrap();
    close(exit.direction.x, incoming.dir.x, 3e-16);
    close(exit.direction.z, incoming.dir.z, 3e-16);
    close(entry.radiance_weight * exit.radiance_weight, 1.0, 3e-16);
    let absorption = BeerLambertAbsorption::try_constant(4.0).unwrap();
    close(
        absorption.transmittance(550.0, distance).unwrap(),
        (-4.0 * distance).exp(),
        2e-15,
    );
    distance
}

#[test]
fn g1_reconstructed_liquid_enters_and_exits_native_snell_and_beer_transport() {
    with_cx(|cx, _| {
        let state = block(cx, 1.0, 1.0);
        let mesh = free_surface_mesh(&state, options(), cx).unwrap();
        check_native_liquid_transport(&mesh, cx, 0.03, 0.09);
    });
}

#[test]
fn g0_surface_budgets_and_unresolved_physical_boundaries_refuse() {
    with_cx(|cx, gate| {
        let state = block(cx, 0.75, 1.0);
        assert!(matches!(
            free_surface_mesh(
                &state,
                FreeSurfaceMeshOptions {
                    max_cells: 1,
                    ..options()
                },
                cx
            ),
            Err(FreeSurfaceMeshError::CellBudgetExceeded { .. })
        ));
        assert!(matches!(
            free_surface_mesh(
                &state,
                FreeSurfaceMeshOptions {
                    max_triangles: 1,
                    ..options()
                },
                cx
            ),
            Err(FreeSurfaceMeshError::Isosurface(
                fs_viz::IsoSurfaceError::TriangleBudgetExceeded { .. }
            ))
        ));
        let wall_contact = FreeSurface3::new(
            12,
            12,
            12,
            0.8,
            [0.0; 3],
            0.0,
            ContactModel3::Neutral,
            |x, y, z| x < 6 && (3..9).contains(&y) && (3..9).contains(&z),
        )
        .unwrap();
        assert!(matches!(
            free_surface_mesh(&wall_contact, options(), cx),
            Err(FreeSurfaceMeshError::WetDomainBoundary { .. })
        ));
        let solid = FreeSurface3::with_solid_tiles(
            12,
            12,
            12,
            0.8,
            [0.0; 3],
            0.0,
            ContactModel3::Neutral,
            &[(0, 0, 0)],
            |x, y, z| (5..9).contains(&x) && (5..9).contains(&y) && (5..9).contains(&z),
        )
        .unwrap();
        assert!(matches!(
            free_surface_mesh(&solid, options(), cx),
            Err(FreeSurfaceMeshError::UnsupportedSolid { .. })
        ));
        gate.request();
        assert!(matches!(
            free_surface_mesh(&state, options(), cx),
            Err(FreeSurfaceMeshError::Cancelled)
        ));
    });
}

fn volume_options() -> FreeSurfaceVolumeOptions {
    FreeSurfaceVolumeOptions {
        absolute_tolerance_m3: 0.0,
        relative_tolerance: 1e-9,
        max_volume_evaluations: 64,
    }
}

fn single_cell(cx: &Cx<'_>, fill_fraction: f64, density: f64) -> FreeSurface3 {
    FreeSurface3::from_cells(
        [12; 3],
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        &[],
        &[FreeSurfaceInitialCell3 {
            coordinates: [5; 3],
            fill_fraction,
            density,
            velocity: [0.0; 3],
        }],
        cx,
    )
    .unwrap()
}

fn mesh_volume(mesh: &TriMesh, reference: [f64; 3]) -> f64 {
    mesh.triangles
        .iter()
        .map(|triangle| {
            let [a, b, c] = triangle.map(|index| {
                core::array::from_fn::<_, 3, _>(|axis| {
                    mesh.vertices[index as usize][axis] - reference[axis]
                })
            });
            (a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

#[test]
fn g1_volume_match_recovers_native_fractional_drop_volume_and_si_scaling() {
    with_cx(|cx, _| {
        for fill in [0.125, 0.25, 1.0] {
            let state = single_cell(cx, fill, 1.0);
            let before = state.grid().state_bits();
            let mass = state.ledger_mass();
            let matched =
                volume_matched_free_surface_mesh(&state, options(), volume_options(), cx).unwrap();
            let expected = fill * 0.01_f64.powi(3);
            close(matched.report.target_volume_m3, expected, 1e-21);
            close(matched.report.enclosed_volume_m3, expected, expected * 1e-9);
            close(
                mesh_volume(&matched.mesh, [0.0; 3]),
                expected,
                expected * 1e-9,
            );
            // An interior grid node belongs to 24 equal native tetrahedra.
            // Its nodal hat has volume 4*(1-level/fill)^3 in cell units.
            let expected_level = fill * (1.0 - (fill / 4.0).cbrt());
            close(matched.report.fill_isovalue, expected_level, 2e-9 * fill);
            assert!(matched.report.volume_evaluations > 1);
            assert!(matched.report.volume_evaluations <= 64);
            assert_eq!(state.grid().state_bits(), before);
            assert_eq!(state.ledger_mass(), mass);
            let replay =
                volume_matched_free_surface_mesh(&state, options(), volume_options(), cx).unwrap();
            assert_eq!(matched.mesh.vertices, replay.mesh.vertices);
            assert_eq!(matched.mesh.triangles, replay.mesh.triangles);
            assert_eq!(
                matched.mesh.bvh_fingerprint(),
                replay.mesh.bvh_fingerprint()
            );
            assert_eq!(matched.report, replay.report);
            let density_scaled = volume_matched_free_surface_mesh(
                &single_cell(cx, fill, 2.0),
                options(),
                volume_options(),
                cx,
            )
            .unwrap();
            assert_eq!(matched.mesh.vertices, density_scaled.mesh.vertices);
            assert_eq!(matched.report, density_scaled.report);
        }
        // Fixed half fill misses half of a full single-cell drop's volume.
        let state = single_cell(cx, 1.0, 1.0);
        let half = free_surface_mesh(&state, options(), cx).unwrap();
        close(mesh_volume(&half, [0.0; 3]), 0.5e-6, 1e-20);
        assert!(matches!(
            free_surface_mesh(&single_cell(cx, 0.25, 1.0), options(), cx),
            Err(FreeSurfaceMeshError::UnresolvedSurface)
        ));
        let origin = [100.0, -20.0, 30.0];
        let placed = volume_matched_free_surface_mesh(
            &single_cell(cx, 0.25, 1.0),
            FreeSurfaceMeshOptions {
                origin_m: origin,
                cell_size_m: 0.03,
                ..options()
            },
            volume_options(),
            cx,
        )
        .unwrap();
        let expected = 0.25 * 0.03_f64.powi(3);
        close(placed.report.target_volume_m3, expected, 1e-20);
        close(mesh_volume(&placed.mesh, origin), expected, expected * 1e-9);
    });
}

#[test]
fn g1_volume_matched_surface_changes_native_refracted_beer_path_consistently() {
    with_cx(|cx, _| {
        let state = block(cx, 1.0, 1.0);
        let matched =
            volume_matched_free_surface_mesh(&state, options(), volume_options(), cx).unwrap();
        let level = matched.report.fill_isovalue;
        assert!(level < 0.5);
        close(matched.report.target_volume_m3, 216e-6, 1e-19);
        close(mesh_volume(&matched.mesh, [0.0; 3]), 216e-6, 216e-15);
        let lower = (2.5 + level) * 0.01;
        let upper = (9.5 - level) * 0.01;
        let distance = check_native_liquid_transport(&matched.mesh, cx, lower, upper);
        assert!(distance > 0.06 / (1.0_f64 - (0.2_f64 / 1.5).powi(2)).sqrt());
    });
}

#[test]
fn g0_volume_matching_refuses_exhausted_passes_unattainable_plateaus_and_lost_precision() {
    with_cx(|cx, gate| {
        let state = single_cell(cx, 1.0, 1.0);
        assert!(matches!(
            volume_matched_free_surface_mesh(
                &state,
                options(),
                FreeSurfaceVolumeOptions {
                    max_volume_evaluations: 1,
                    ..volume_options()
                },
                cx,
            ),
            Err(FreeSurfaceMeshError::VolumeMatchExhausted {
                volume_evaluations: 1,
                ..
            })
        ));
        let mut cells = Vec::new();
        for z in 3..9 {
            for y in 3..9 {
                for x in 3..9 {
                    cells.push(FreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 0.25,
                        density: 1.0,
                        velocity: [0.0; 3],
                    });
                }
            }
        }
        let plateau = FreeSurface3::from_cells(
            [12; 3],
            0.8,
            [0.0; 3],
            0.0,
            ContactModel3::Neutral,
            &[],
            &cells,
            cx,
        )
        .unwrap();
        // The requested 54 cells lies strictly between zero and the 125-cell
        // constant plateau: no single global threshold has that volume.
        assert!(matches!(
            volume_matched_free_surface_mesh(&plateau, options(), volume_options(), cx),
            Err(FreeSurfaceMeshError::VolumeMatchExhausted { .. })
        ));
        assert!(matches!(
            volume_matched_free_surface_mesh(
                &state,
                FreeSurfaceMeshOptions {
                    origin_m: [1e9; 3],
                    ..options()
                },
                FreeSurfaceVolumeOptions {
                    relative_tolerance: 1e-12,
                    ..volume_options()
                },
                cx,
            ),
            Err(FreeSurfaceMeshError::VolumeMismatch { .. })
        ));
        gate.request();
        assert!(matches!(
            volume_matched_free_surface_mesh(&state, options(), volume_options(), cx),
            Err(FreeSurfaceMeshError::Cancelled)
        ));
    });
}
