//! Actual resolved surface traction on the conservative native mass routes.

use super::super::thermocapillary::HeightGraphSurface3;
use super::*;

fn graph() -> HeightGraphSurface3 {
    HeightGraphSurface3 {
        axis: 2,
        outward_sign: 1,
        maximum_slope: 0.5,
        minimum_fill: 0.001,
    }
}

fn surface(density: impl Fn([usize; 3]) -> f64) -> FreeSurface3 {
    let mut cells = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                // Keep every initial cut volume >= 0.02; integer-height
                // columns can otherwise produce unresolved vanishing films.
                let height = 3.52 + 0.12 * (x as f64 - 3.5) - 0.08 * (y as f64 - 3.5);
                let fill = (height - z as f64).clamp(0.0, 1.0);
                if fill > 0.0 {
                    cells.push(FreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: fill,
                        density: density([x, y, z]),
                        velocity: [0.0; 3],
                    });
                }
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 83,
                kernel_id: 37,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        FreeSurface3::from_cells(
            [8; 3],
            0.8,
            [0.0; 3],
            0.002,
            ContactModel3::Neutral,
            &[],
            &cells,
            &cx,
        )
        .unwrap()
        .with_conversion_redistribution(ConversionRedistribution3::ConservativeLocal)
    })
}

fn tension(surface: &FreeSurface3) -> BTreeMap<(u64, usize), f64> {
    surface
        .mass_snapshot()
        .keys()
        .filter_map(|&location| {
            let MassLocation3::Cell(key, lane) = location else {
                return None;
            };
            let (x, y, _) = FreeSurface3::global_of(key, lane);
            Some(((key, lane), 0.002 + 1e-5 * x as f64 - 5e-6 * y as f64))
        })
        .collect()
}

fn population_bits(surface: &FreeSurface3) -> Vec<u64> {
    let mut bits = Vec::new();
    for &key in surface.grid.active_keys() {
        let slot = surface.grid.slot_of(key).unwrap();
        for lane in 0..TILE_CELLS {
            bits.extend(surface.grid.populations(slot, lane).map(f64::to_bits));
        }
    }
    bits
}

fn momentum(surface: &FreeSurface3) -> [f64; 3] {
    let mut total = [0.0; 3];
    for (location, mass) in surface.mass_snapshot() {
        let Some([x, y, z]) = FreeSurface3::transport_coordinates(location) else {
            continue;
        };
        let velocity = surface.sample_at(x, y, z).unwrap().velocity;
        for d in 0..3 {
            total[d] += mass * velocity[d];
        }
    }
    total
}

fn check_routes(trace: &MassTransport3) {
    let mut balance = trace.before.clone();
    for route in &trace.routes {
        *balance.entry(route.from).or_default() -= route.mass;
        *balance.entry(route.to).or_default() += route.mass;
    }
    for (&location, &mass) in &trace.after {
        *balance.entry(location).or_default() -= mass;
    }
    for (location, residual) in balance {
        assert!(residual.abs() < 4e-13, "{location:?}: {residual}");
    }
}

fn check_traction(
    surface: &FreeSurface3,
    forces: &BTreeMap<(u64, usize), [f64; 3]>,
    gate: &CancelGate,
) {
    let patches = surface.height_graph_geometry(graph(), gate).unwrap();
    let mut integrated = BTreeMap::<[usize; 2], [f64; 3]>::new();
    for (&(key, lane), force) in forces {
        let (x, y, z) = FreeSurface3::global_of(key, lane);
        let sample = surface
            .sample_at(x as usize, y as usize, z as usize)
            .unwrap();
        let total = integrated.entry([x as usize, y as usize]).or_default();
        // Native force is per occupied volume. The density-aware mass form
        // m*(F/rho) must give the same physical patch traction.
        let mass = surface.tiles[&key].mass[lane];
        for d in 0..3 {
            total[d] += mass * force[d] / sample.density;
        }
    }
    for patch in patches {
        let gradient = [1e-5, -5e-6, 0.0];
        let normal = patch.outward_normal;
        let projected: f64 = gradient.iter().zip(normal).map(|(g, n)| g * n).sum();
        for d in 0..3 {
            let expected = (gradient[d] - projected * normal[d]) * patch.area;
            assert!((integrated[&patch.column][d] - expected).abs() < 2e-17);
        }
    }
}

#[test]
fn g1_resolved_thermocapillary_step_reuses_conservative_routes_density_and_conversion_receipt() {
    let gate = CancelGate::new_clock_free();
    // A weak numerical pressure-density field, not a thermodynamic density law.
    let mut surface = surface(|p| 1.0 + 2e-5 * (p[0] as f64 - 3.5));
    let mass = surface.ledger_mass();
    for _ in 0..3 {
        let mut sigma = tension(&surface);
        let traction = surface
            .height_graph_marangoni_forces(graph(), &mut sigma, &gate)
            .unwrap();
        check_traction(&surface, &traction, &gate);
        surface.rebind_surface_tensions(sigma, &gate).unwrap();
        let forces = surface
            .mass_snapshot()
            .keys()
            .filter_map(|location| match *location {
                MassLocation3::Cell(key, lane) => Some((
                    (key, lane),
                    traction.get(&(key, lane)).copied().unwrap_or([0.0; 3]),
                )),
                MassLocation3::Carry => None,
            })
            .collect();
        surface.rebind_force_field(forces, &gate).unwrap();
        let mut trace = MassTransport3::new(surface.mass_snapshot(), 100_000);
        surface
            .step_inner(&gate, Some(&mut trace), None, None)
            .unwrap();
        let before = trace.before.clone();
        let routes = trace.routes.len();
        let receipt = surface.conversion_momentum_report().copied().unwrap();
        let populations = population_bits(&surface);
        surface.resolve_interface_fills(&mut trace, &gate).unwrap();
        assert_eq!(surface.conversion_momentum_report(), Some(&receipt));
        assert_eq!(population_bits(&surface), populations);
        assert_eq!(trace.routes.len(), routes);
        assert_eq!(trace.before, before);
        check_routes(&trace);
        let mut endpoint_sigma = tension(&surface);
        let endpoint = surface
            .height_graph_marangoni_forces(graph(), &mut endpoint_sigma, &gate)
            .unwrap();
        check_traction(&surface, &endpoint, &gate);
        assert!((surface.ledger_mass() - mass).abs() < 2e-12 * mass);
    }
}

#[test]
fn g1_resolved_volume_traction_is_independent_of_the_mass_density_representation() {
    let gate = CancelGate::new_clock_free();
    for scale in [0.5, 1.0, 2.0] {
        let mut surface = surface(|p| scale * (0.8 + 0.1 * (p[0] % 3) as f64));
        let mut trace = MassTransport3::new(surface.mass_snapshot(), 1000);
        let bits = population_bits(&surface);
        surface.resolve_interface_fills(&mut trace, &gate).unwrap();
        assert_eq!(population_bits(&surface), bits);
        assert_eq!(trace.before, trace.after);
        assert!(trace.routes.is_empty());
        assert!(surface.conversion_momentum_report().is_none());
        let mut sigma = tension(&surface);
        let forces = surface
            .height_graph_marangoni_forces(graph(), &mut sigma, &gate)
            .unwrap();
        check_traction(&surface, &forces, &gate);
    }
}

#[test]
fn g1_unresolved_conservative_geometry_appends_real_momentum_routes_once() {
    let gate = CancelGate::new_clock_free();
    let mut surface = surface(|_| 1.0);
    let [x, y, z] = [3, 3, 3];
    let (slot, lane) = surface.grid.resolve_source(x, y, z).unwrap();
    let key = surface.grid.active_keys()[slot];
    let rho = surface.grid.populations(slot, lane).iter().sum::<f64>();
    let previous = surface
        .sample_at(x as usize, y as usize, z as usize)
        .unwrap()
        .velocity;
    let populations = super::super::forcing::retarget_velocity(
        surface.grid.populations(slot, lane),
        previous,
        [0.03, -0.02, 0.01],
        surface.force_at_lane(key, lane),
    )
    .unwrap();
    surface.grid.set_populations(slot, lane, populations);
    surface.tiles.get_mut(&key).unwrap().mass[lane] = rho + 0.001;
    let mut trace = MassTransport3::new(surface.mass_snapshot(), 1000);
    let before = trace.before.clone();
    let initial_momentum = momentum(&surface);
    surface.resolve_interface_fills(&mut trace, &gate).unwrap();
    assert_eq!(trace.before, before);
    assert!(!trace.routes.is_empty());
    check_routes(&trace);
    let final_momentum = momentum(&surface);
    for d in 0..3 {
        assert!((final_momentum[d] - initial_momentum[d]).abs() < 2e-14);
    }
    let receipt = surface.conversion_momentum_report().copied().unwrap();
    assert!(receipt.kinetic_energy_change < 0.0);
    let routes = trace.routes.len();
    let bits = population_bits(&surface);
    surface.resolve_interface_fills(&mut trace, &gate).unwrap();
    assert_eq!(trace.routes.len(), routes);
    assert_eq!(surface.conversion_momentum_report(), Some(&receipt));
    assert_eq!(population_bits(&surface), bits);
}

#[test]
fn g4_conservative_geometry_checks_actual_density_fill_carry_and_cancellation() {
    let gate = CancelGate::new_clock_free();
    for invalid in 0..4 {
        let mut surface = surface(|_| 1.0);
        let mut trace = MassTransport3::new(surface.mass_snapshot(), 1000);
        let (slot, lane) = surface.grid.resolve_source(3, 3, 3).unwrap();
        let key = surface.grid.active_keys()[slot];
        match invalid {
            0 => surface.tiles.get_mut(&key).unwrap().mass[lane] = f64::NAN,
            1 => surface.grid.set_populations(slot, lane, [0.0; Q3]),
            2 => surface.carry = 0.001,
            _ => surface.redistribution = ConversionRedistribution3::LegacyGlobal,
        }
        let before = trace.before.clone();
        let bits = population_bits(&surface);
        assert!(surface.resolve_interface_fills(&mut trace, &gate).is_err());
        assert_eq!(trace.before, before);
        assert!(trace.routes.is_empty());
        assert_eq!(population_bits(&surface), bits);
    }
    let mut surface = surface(|_| 1.0);
    let mut trace = MassTransport3::new(surface.mass_snapshot(), 1000);
    gate.request();
    assert!(matches!(
        surface.resolve_interface_fills(&mut trace, &gate),
        Err(FreeSurfaceError3::Cancelled)
    ));
    assert!(trace.routes.is_empty());
}
