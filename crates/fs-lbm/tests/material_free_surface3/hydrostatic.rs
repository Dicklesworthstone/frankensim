//! Balanced rest and actual pressure/traction under the material rho*a source.

use super::*;
use fs_exec::{Budget, Cx, ExecMode, StreamKey};
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialFreeSurfaceInitialCell3, MaterialSubmergedSolid3,
    SubmergedSolidLoadControls3,
};

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 47,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn small_domain(gravity: [f64; 3]) -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [8; 3],
        gravity_m_s2: gravity,
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 8 * 8 * 8,
    }
}

fn make<'a>(
    setup: MaterialFreeSurfaceScaling<'a>,
    domain: MaterialFreeSurfaceDomain3<'_>,
    wall: PlanarWall3,
    depth: usize,
) -> MaterialFreeSurface3<'a> {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        setup
            .free_surface3_hydrostatic(domain, wall, depth, 0.1, cx)
            .unwrap()
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn assert_same_mass(left: &FreeSurface3, right: &FreeSurface3) {
    let [nx, ny, nz] = left.domain_cells();
    assert_eq!(left.domain_cells(), right.domain_cells());
    assert_eq!(left.carry().to_bits(), right.carry().to_bits());
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                assert_eq!(
                    left.sample_at(x, y, z).map(|sample| sample.cell),
                    right.sample_at(x, y, z).map(|sample| sample.cell),
                    "cell kind at ({x}, {y}, {z})"
                );
                // A one-cell region retains the actual interface mass, even
                // when its sampled fill fraction is clamped to zero or one.
                let cell = |a, b, c| (a, b, c) == (x, y, z);
                assert_eq!(
                    left.region_mass(cell).to_bits(),
                    right.region_mass(cell).to_bits(),
                    "cell mass at ({x}, {y}, {z})"
                );
            }
        }
    }
}

#[test]
fn actual_material_density_gravity_stays_at_rest_in_all_six_directions() {
    let source = liquid(0, 300.0, None);
    let setup = plan(&source).with_pressure_reference(100_000.0).unwrap();
    for (wall, axis, sign) in [
        (PlanarWall3::XMin, 0, -1.0),
        (PlanarWall3::XMax, 0, 1.0),
        (PlanarWall3::YMin, 1, -1.0),
        (PlanarWall3::YMax, 1, 1.0),
        (PlanarWall3::ZMin, 2, -1.0),
        (PlanarWall3::ZMax, 2, 1.0),
    ] {
        let mut gravity = [0.0; 3];
        gravity[axis] = sign * 9.80665;
        let acceleration = gravity.map(|g| setup.gravity_lattice(g).unwrap());
        // Five layers put the pressure plane inside a tile, not at a special
        // tile boundary. Both walls and diagonal gas links must balance it.
        let mut flow = make(setup, small_domain(gravity), wall, 5);
        assert!(std::ptr::eq(flow.scaling().material(), &source));
        assert_eq!(flow.lattice().body_acceleration(), Some(acceleration));
        let mut initial = Vec::new();
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..8 {
                    if let Some(sample) = flow.lattice().sample_at(x, y, z) {
                        initial.push(([x, y, z], sample));
                    }
                }
            }
        }
        assert_eq!(initial.len(), 8 * 8 * 5);
        let mass = flow.ledger_mass_kg().unwrap();
        for _ in 0..64 {
            flow.step().unwrap();
        }
        for (position, before) in initial {
            let [x, y, z] = position;
            let after = flow.lattice().sample_at(x, y, z).unwrap();
            assert!(
                after.velocity.iter().all(|u| u.abs() < 2e-11),
                "{wall:?} {position:?}: {after:?}"
            );
            near(after.density, before.density, 2e-11);
            near(after.fill_fraction, 1.0, 2e-11);
            let force = flow.lattice().force_density_at(x, y, z).unwrap();
            for a in 0..3 {
                near(force[a], after.density * acceleration[a], 2e-17);
            }
        }
        near(flow.ledger_mass_kg().unwrap(), mass, 2e-11 * mass);
    }
}

#[test]
fn zero_gravity_is_the_original_uniform_pressure_initial_state() {
    let source = liquid(0, 300.0, None);
    let setup = plan(&source).with_pressure_reference(100_000.0).unwrap();
    let domain = small_domain([0.0; 3]);
    let mut balanced = make(setup, domain, PlanarWall3::ZMin, 5);
    let initial: Vec<_> = (0..5)
        .flat_map(|z| {
            (0..8).flat_map(move |y| {
                (0..8).map(move |x| MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    velocity_m_s: [0.0; 3],
                })
            })
        })
        .collect();
    let mut ordinary = with_cx(&CancelGate::new_clock_free(), |cx| {
        setup
            .free_surface3_from_cells(domain, &initial, cx)
            .unwrap()
    });
    for _ in 0..8 {
        assert_eq!(
            balanced.lattice().grid().state_bits(),
            ordinary.lattice().grid().state_bits()
        );
        assert_same_mass(balanced.lattice(), ordinary.lattice());
        balanced.step().unwrap();
        ordinary.step().unwrap();
    }
    for state in initial {
        let [x, y, z] = state.coordinates;
        let sample = balanced.sample_at_si(x, y, z).unwrap().unwrap();
        near(sample.pressure_perturbation_pa, 0.0, 2e-10);
        assert!(sample.velocity_m_s.iter().all(|u| u.abs() < 2e-12));
    }
}

#[test]
fn changing_the_absolute_pressure_datum_does_not_change_rest_or_loading() {
    let source = liquid(0, 300.0, None);
    let domain = small_domain([0.0, 0.0, -9.80665]);
    let mut low = make(
        plan(&source).with_pressure_reference(100_000.0).unwrap(),
        domain,
        PlanarWall3::ZMin,
        5,
    );
    let mut high = make(
        plan(&source).with_pressure_reference(105_000.0).unwrap(),
        domain,
        PlanarWall3::ZMin,
        5,
    );
    for _ in 0..8 {
        assert_eq!(
            low.lattice().grid().state_bits(),
            high.lattice().grid().state_bits()
        );
        assert_same_mass(low.lattice(), high.lattice());
        let rho = low.lattice().sample_at(3, 3, 0).unwrap().density;
        near(
            high.scaling().absolute_pressure_si(rho).unwrap()
                - low.scaling().absolute_pressure_si(rho).unwrap(),
            5_000.0,
            2e-11,
        );
        low.step().unwrap();
        high.step().unwrap();
    }
}

#[test]
fn native_pressure_and_submerged_traction_converge_as_compressibility_decreases() {
    let source = liquid(0, 300.0, None);
    let dx = 0.001;
    let gravity = 9.80665;
    let solid = [(1, 1, 1)];
    let domain = MaterialFreeSurfaceDomain3 {
        cells: [16; 3],
        gravity_m_s2: [0.0, 0.0, -gravity],
        solid_tiles: &solid,
        max_cells: 16 * 16 * 16,
        contact: ContactModel3::Neutral,
    };
    let origin = [2.0 * dx, 3.0 * dx, dx];
    let object = [MaterialSubmergedSolid3 {
        solid_tiles: &solid,
        moment_origin_m: origin,
    }];
    let controls = SubmergedSolidLoadControls3 {
        maximum_solid_tiles: 1,
        maximum_boundary_links: 1_000,
        maximum_relative_mass_defect: 1e-11,
    };
    let gate = CancelGate::new_clock_free();
    let volume = (4.0_f64 * dx).powi(3);
    let archimedes = source.density_kg_m3() * volume * gravity;
    let continuum_pressure = source.density_kg_m3() * gravity * 11.5 * dx;
    let mut pressure_errors = Vec::new();
    let mut force_errors = Vec::new();
    // Keep physical geometry and gravity fixed. Decreasing dt/dx removes
    // the weak-compressibility error; spatial refinement alone would not.
    for dt in [0.0004, 0.0002, 0.0001] {
        let setup = plan_free_surface_from_material(&source, dx, dt, 0.1)
            .unwrap()
            .with_pressure_reference(100_000.0)
            .unwrap();
        let mut flow = make(setup, domain, PlanarWall3::ZMin, 12);
        let pressure = flow
            .sample_at_si(1, 1, 0)
            .unwrap()
            .unwrap()
            .pressure_perturbation_pa;
        pressure_errors.push((pressure / continuum_pressure - 1.0).abs());
        let mass = flow.ledger_mass_kg().unwrap();
        let mut last_force = 0.0;
        for _ in 0..8 {
            let report = flow
                .step_with_submerged_solid_loads_si(&object, controls, &gate)
                .unwrap();
            let load = report.bodies[0];
            assert!(report.maximum_speed_m_s < 2e-10);
            near(load.force_n[0], 0.0, 2e-12);
            near(load.force_n[1], 0.0, 2e-12);
            near(
                load.torque_n_m[0],
                (6.0 * dx - origin[1]) * load.force_n[2],
                2e-14,
            );
            near(
                load.torque_n_m[1],
                -(6.0 * dx - origin[0]) * load.force_n[2],
                2e-14,
            );
            near(load.torque_n_m[2], 0.0, 2e-14);
            last_force = load.force_n[2];
        }
        force_errors.push((last_force / archimedes - 1.0).abs());
        near(flow.ledger_mass_kg().unwrap(), mass, 2e-11 * mass);
    }
    for errors in [&pressure_errors, &force_errors] {
        assert!(errors[0] < 0.04, "{errors:?}");
        assert!(errors[2] < 0.003, "{errors:?}");
        for pair in errors.windows(2) {
            assert!(pair[1] < 0.28 * pair[0], "{errors:?}");
        }
    }
}

#[test]
fn incompatible_or_cancelled_columns_refuse_without_changing_accepted_state() {
    let source = liquid(0, 300.0, None);
    let setup = plan(&source).with_pressure_reference(100_000.0).unwrap();
    let domain = small_domain([0.0, 0.0, -9.80665]);
    let accepted = make(setup, domain, PlanarWall3::ZMin, 5);
    let bits = accepted.lattice().grid().state_bits();
    let mass = accepted.ledger_mass_kg().unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        assert!(
            plan(&source)
                .free_surface3_hydrostatic(domain, PlanarWall3::ZMin, 5, 0.1, cx)
                .is_err()
        );
        for depth in [0, 8, usize::MAX] {
            assert!(
                setup
                    .free_surface3_hydrostatic(domain, PlanarWall3::ZMin, depth, 0.1, cx)
                    .is_err()
            );
        }
        for bound in [0.0, 1.0, f64::NAN, 1e-6] {
            assert!(
                setup
                    .free_surface3_hydrostatic(domain, PlanarWall3::ZMin, 5, bound, cx)
                    .is_err()
            );
        }
        for gravity in [
            [0.0, 0.0, 9.80665],
            [0.1, 0.0, -9.80665],
            [0.0, 0.0, f64::NAN],
            [0.0, 0.0, -1e9],
        ] {
            assert!(
                setup
                    .free_surface3_hydrostatic(small_domain(gravity), PlanarWall3::ZMin, 5, 0.1, cx)
                    .is_err()
            );
        }
        let oversized = MaterialFreeSurfaceDomain3 {
            max_cells: 511,
            ..domain
        };
        assert!(
            setup
                .free_surface3_hydrostatic(oversized, PlanarWall3::ZMin, 5, 0.1, cx)
                .is_err()
        );
        let wetting = MaterialFreeSurfaceDomain3 {
            contact: ContactModel3::Wetting,
            ..domain
        };
        assert!(
            setup
                .free_surface3_hydrostatic(wetting, PlanarWall3::ZMin, 5, 0.1, cx)
                .is_err()
        );
        let capillary_source = liquid(0, 300.0, Some(0.05));
        let capillary = plan(&capillary_source)
            .with_pressure_reference(100_000.0)
            .unwrap();
        assert!(
            capillary
                .free_surface3_hydrostatic(domain, PlanarWall3::ZMin, 5, 0.1, cx)
                .is_err()
        );
    });
    gate.request();
    let refusal = with_cx(&gate, |cx| {
        setup.free_surface3_hydrostatic(domain, PlanarWall3::ZMin, 5, 0.1, cx)
    });
    assert!(matches!(
        refusal,
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    ));
    assert_eq!(accepted.lattice().grid().state_bits(), bits);
    assert_eq!(accepted.ledger_mass_kg().unwrap(), mass);
    let retry = make(setup, domain, PlanarWall3::ZMin, 5);
    assert_eq!(retry.lattice().grid().state_bits(), bits);
}
