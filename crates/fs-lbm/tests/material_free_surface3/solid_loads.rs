//! Physical source, scaling and reciprocal load tests on immersed solids.
use super::*;
use fs_exec::{Budget, Cx, ExecMode, StreamKey};
use fs_lbm::d3q19::TILE;
use fs_lbm::d3q19::freesurface3::SubmergedSolid3;
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialFreeSurfaceInitialCell3, MaterialSubmergedSolid3,
    SubmergedSolidLoadControls3,
};

fn controls() -> SubmergedSolidLoadControls3 {
    SubmergedSolidLoadControls3 {
        maximum_solid_tiles: 1,
        maximum_boundary_links: 1_000,
        maximum_relative_mass_defect: 1e-11,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn hydrostatic<'a>(
    setup: MaterialFreeSurfaceScaling<'a>,
    reference: f64,
    gravity: [f64; 3],
) -> MaterialFreeSurface3<'a> {
    let solids = [(1, 1, 1)];
    let input = MaterialFreeSurfaceDomain3 {
        solid_tiles: &solids,
        gravity_m_s2: gravity,
        ..domain()
    };
    let mut initial = Vec::new();
    for z in 0..input.cells[2] {
        for y in 0..input.cells[1] {
            for x in 0..input.cells[0] {
                if (x / TILE, y / TILE, z / TILE) == (1, 1, 1) {
                    continue;
                }
                initial.push(MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    velocity_m_s: [0.0; 3],
                });
            }
        }
    }
    let density = setup.material().density_kg_m3();
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 43,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        setup
            .free_surface3_from_cells_with_pressure(
                input,
                &initial,
                0.1,
                |p| reference + density * (0..3).map(|a| gravity[a] * p[a]).sum::<f64>(),
                &cx,
            )
            .unwrap()
    })
}

#[test]
fn actual_si_buoyancy_and_torque_match_archimedes_across_physical_scales() {
    let source = liquid(0, 300.0, None);
    let gravity = [0.5, -0.25, -9.80665];
    let solid = [(1, 1, 1)];
    let gate = CancelGate::new_clock_free();
    for (dx, dt) in [(0.001, 0.0002), (0.002, 0.0004)] {
        let setup = plan_free_surface_from_material(&source, dx, dt, 0.1)
            .unwrap()
            .with_pressure_reference(100_000.0)
            .unwrap();
        let mut flow = hydrostatic(setup, 100_000.0, gravity);
        let mut native = flow.lattice().clone();
        let origin = [2.0 * dx, 3.0 * dx, dx];
        let raw = native
            .step_with_submerged_solid_loads(
                &[SubmergedSolid3 {
                    solid_tiles: &solid,
                    moment_origin: origin.map(|x| x / dx),
                }],
                controls(),
                &gate,
            )
            .unwrap();
        let result = flow
            .step_with_submerged_solid_loads_si(
                &[MaterialSubmergedSolid3 {
                    solid_tiles: &solid,
                    moment_origin_m: origin,
                }],
                controls(),
                &gate,
            )
            .unwrap();
        let body = result.bodies[0];
        let volume = (TILE as f64 * dx).powi(3);
        let expected = gravity.map(|g| -source.density_kg_m3() * volume * g);
        let arm = core::array::from_fn::<_, 3, _>(|a| 1.5 * TILE as f64 * dx - origin[a]);
        let expected_torque = [
            arm[1] * expected[2] - arm[2] * expected[1],
            arm[2] * expected[0] - arm[0] * expected[2],
            arm[0] * expected[1] - arm[1] * expected[0],
        ];
        assert_eq!(
            native.grid().state_bits(),
            flow.lattice().grid().state_bits()
        );
        assert!(std::ptr::eq(flow.scaling().material(), &source));
        assert_eq!(result.geometry_epoch, flow.lattice().geometry_epoch());
        for a in 0..3 {
            near(body.force_n[a], expected[a], 3e-12);
            near(body.torque_n_m[a], expected_torque[a], 3e-14);
            near(body.impulse_kg_m_s[a], body.force_n[a] * dt, 1e-20);
            near(
                body.angular_impulse_kg_m2_s[a],
                body.torque_n_m[a] * dt,
                1e-22,
            );
            near(
                body.impulse_kg_m_s[a],
                raw.bodies[0].wall_impulse[a]
                    * setup.mass_unit_kg()
                    * setup.velocity_si(1.0).unwrap(),
                1e-20,
            );
        }
        near(result.mass_after_kg, result.mass_before_kg, 1e-12);
    }
}

#[test]
fn closed_loads_are_pressure_datum_independent_and_si_refusals_publish_nothing() {
    let source = liquid(0, 300.0, None);
    let solids = [(1, 1, 1)];
    let objects = [MaterialSubmergedSolid3 {
        solid_tiles: &solids,
        moment_origin_m: [0.006; 3],
    }];
    let gate = CancelGate::new_clock_free();
    let mut reference = None;
    for pressure in [100_000.0, 105_000.0] {
        let setup = plan(&source).with_pressure_reference(pressure).unwrap();
        let mut flow = hydrostatic(setup, pressure, [0.0, 0.0, -9.80665]);
        let report = flow
            .step_with_submerged_solid_loads_si(&objects, controls(), &gate)
            .unwrap();
        let state = (flow.lattice().grid().state_bits(), report);
        if let Some(expected) = &reference {
            assert_eq!(&state, expected);
        }
        reference = Some(state);
    }

    // Binary ingress has physical u=g*dt/2. The first native candidate stays
    // well inside the lattice Mach range but exceeds this source's 0.1 m/s
    // scenario speed, so the material owner must reject AFTER native stepping.
    let input = MaterialFreeSurfaceDomain3 {
        solid_tiles: &solids,
        gravity_m_s2: [500.0, 0.0, 0.0],
        ..domain()
    };
    let mut flow = plan(&source)
        .free_surface3(input, |x, y, z| (x / TILE, y / TILE, z / TILE) != (1, 1, 1))
        .unwrap();
    let bits = flow.lattice().grid().state_bits();
    let mass = flow.ledger_mass_kg().unwrap();
    assert!(matches!(
        flow.step_with_submerged_solid_loads_si(&objects, controls(), &gate,),
        Err(MaterialLbmError::Invalid(_))
    ));
    assert_eq!(flow.lattice().grid().state_bits(), bits);
    assert_eq!(flow.ledger_mass_kg().unwrap(), mass);
    flow.rebind_gravity([0.0; 3], &gate).unwrap();
    flow.step_with_submerged_solid_loads_si(&objects, controls(), &gate)
        .unwrap();
}
