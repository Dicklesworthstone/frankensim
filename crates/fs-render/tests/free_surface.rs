use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_lbm::d3q19::freesurface3::{ContactModel3, FreeSurface3, FreeSurfaceInitialCell3};
use fs_render::charts::Ray;
use fs_render::dielectric::{BeerLambertAbsorption, DielectricEvent, sample_smooth_dielectric};
use fs_render::free_surface::{FreeSurfaceMeshError, FreeSurfaceMeshOptions, free_surface_mesh};

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

#[test]
fn g1_reconstructed_liquid_enters_and_exits_native_snell_and_beer_transport() {
    with_cx(|cx, _| {
        let state = block(cx, 1.0, 1.0);
        let mesh = free_surface_mesh(&state, options(), cx).unwrap();
        let incoming = Ray {
            origin: Point3::new(0.045, 0.056, 0.0),
            dir: Vec3::new(0.2, 0.0, 0.96_f64.sqrt()),
        };
        let first = mesh.intersect_with_cx(cx, &incoming).unwrap().unwrap();
        close(first.point.z, 0.03, 2e-16);
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
        close(last.point.z, 0.09, 2e-16);
        assert_eq!(last.normal.unwrap(), Vec3::new(0.0, 0.0, 1.0));
        let distance = last.t + offset;
        close(
            distance,
            0.06 / (1.0_f64 - (0.2_f64 / 1.5).powi(2)).sqrt(),
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
