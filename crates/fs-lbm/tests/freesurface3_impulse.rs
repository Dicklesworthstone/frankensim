//! G0/G1/G3: impulses act on original retained liquid mass without a VOF step.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_lbm::d3q19::freesurface3::{Cell3, ContactModel3, FreeSurface3, FreeSurfaceInitialCell3};
use fs_lbm::d3q19::{E3, Q3};
use fs_lbm::impulse::{ImpulseControls, ImpulseError};

const FORCE: [f64; 3] = [0.0002, -0.0001, 0.0003];

fn controls() -> ImpulseControls {
    ImpulseControls {
        max_cells: 256,
        maximum_speed: 0.15,
        maximum_momentum_defect: 1e-11,
        maximum_work_defect: 1e-12,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn coordinates(index: usize) -> [usize; 3] {
    [index % 4, (index / 4) % 4, index / 16]
}

fn seeded(initial: &[FreeSurfaceInitialCell3], force: [f64; 3]) -> FreeSurface3 {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        FreeSurface3::from_cells(
            [4; 3],
            0.83,
            force,
            0.0,
            ContactModel3::Neutral,
            &[],
            initial,
            &cx,
        )
        .unwrap()
    })
}

fn partial_fixture() -> FreeSurface3 {
    let initial: Vec<_> = (0..64)
        .map(|i| FreeSurfaceInitialCell3 {
            coordinates: coordinates(i),
            fill_fraction: match i {
                21 => 0.25,
                42 => 0.5,
                _ => 1.0,
            },
            density: 1.2,
            velocity: [0.04, -0.02, 0.01],
        })
        .collect();
    seeded(&initial, FORCE)
}

fn mass_at(flow: &FreeSurface3, position: [usize; 3]) -> f64 {
    flow.region_mass(|x, y, z| [x, y, z] == position)
}

fn geometry_held(before: &FreeSurface3, after: &FreeSurface3) {
    assert_eq!(after.domain_cells(), before.domain_cells());
    assert_eq!(
        after.ledger_mass().to_bits(),
        before.ledger_mass().to_bits()
    );
    assert_eq!(after.carry().to_bits(), before.carry().to_bits());
    assert_eq!(after.conversions(), before.conversions());
    assert_eq!(after.geometry_epoch(), before.geometry_epoch());
    assert_eq!(after.grid().active_tiles(), before.grid().active_tiles());
    assert_eq!(after.grid().steps(), before.grid().steps());
    let [nx, ny, nz] = before.domain_cells();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                assert_eq!(
                    after.force_density_at(x, y, z),
                    before.force_density_at(x, y, z)
                );
                match (before.sample_at(x, y, z), after.sample_at(x, y, z)) {
                    (Some(a), Some(b)) => {
                        assert_eq!(b.cell, a.cell);
                        assert_eq!(b.density.to_bits(), a.density.to_bits());
                        assert_eq!(b.fill_fraction.to_bits(), a.fill_fraction.to_bits());
                    }
                    (None, None) => {}
                    _ => panic!("impulse changed occupancy at ({x},{y},{z})"),
                }
            }
        }
    }
}

// The public sparse snapshot is tile-major, then direction, then lane. All
// moment-oracle fixtures use exactly one 4^3 tile. Integrate the original
// populations directly, independently of the impulse kernel or its observer.
fn moments(bits: &[u64], lane: usize, force: [f64; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    assert_eq!(bits.len(), Q3 * 64);
    let f: [f64; Q3] = core::array::from_fn(|q| f64::from_bits(bits[q * 64 + lane]));
    let rho = f.iter().sum::<f64>();
    let e = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    let velocity: [f64; 3] = core::array::from_fn(|a| {
        (f.iter().zip(e).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a]) / rho
    });
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            f.iter().zip(e).map(|(f, e)| f * e[a] * e[b]).sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + velocity[b] * force[a])
                - rho * velocity[a] * velocity[b]
                - if a == b { rho / 3.0 } else { 0.0 }
        })
    });
    (rho, velocity, stress)
}

#[test]
fn partial_cell_inertia_and_actual_liquid_work_use_retained_mass() {
    let mut flow = partial_fixture();
    let before = flow.clone();
    let original = flow.grid().state_bits();
    let targets = [0, 21, 42];
    let j = [0.006, -0.003, 0.0015];
    let mut impulses = vec![[0.0; 3]; 64];
    for i in targets {
        impulses[i] = j;
    }
    let report = flow
        .apply_impulse_density(&impulses, controls(), &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(report.affected_cells, 3);
    assert_eq!(report.mass_change, 0.0);
    geometry_held(&before, &flow);
    let after_bits = flow.grid().state_bits();
    let mut observed = [0.0; 3];
    let mut work = 0.0;
    let mut kinetic = 0.0;
    for i in 0..64 {
        if !targets.contains(&i) {
            for q in 0..Q3 {
                assert_eq!(after_bits[q * 64 + i], original[q * 64 + i]);
            }
            continue;
        }
        let [x, y, z] = coordinates(i);
        let old = before.sample_at(x, y, z).unwrap();
        let next = flow.sample_at(x, y, z).unwrap();
        let mass = mass_at(&before, [x, y, z]);
        assert_eq!(
            old.cell,
            if i == 0 {
                Cell3::Fluid
            } else {
                Cell3::Interface
            }
        );
        for a in 0..3 {
            let change = mass * (next.velocity[a] - old.velocity[a]);
            close(next.velocity[a] - old.velocity[a], j[a] / mass, 3e-15);
            observed[a] += change;
            work += j[a] * (next.velocity[a] + old.velocity[a]) * 0.5;
            kinetic += 0.5 * mass * (next.velocity[a].powi(2) - old.velocity[a].powi(2));
        }
    }
    for a in 0..3 {
        close(report.applied_impulse[a], 3.0 * j[a], 1e-16);
        close(report.momentum_change[a], observed[a], 3e-15);
    }
    close(report.work, work, 1e-15);
    close(report.kinetic_energy_change, kinetic, 1e-15);
    assert!(report.work > 0.0);
}

#[test]
fn signed_braking_preserves_actual_viscous_stress_and_held_guo_force() {
    let initial: Vec<_> = (0..64)
        .map(|i| {
            let [x, y, z] = coordinates(i);
            FreeSurfaceInitialCell3 {
                coordinates: [x, y, z],
                fill_fraction: 1.0,
                density: 1.0,
                velocity: [
                    0.035 + 0.001 * x as f64,
                    -0.02 + 0.001 * y as f64,
                    0.01 - 0.001 * z as f64,
                ],
            }
        })
        .collect();
    let mut flow = seeded(&initial, FORCE);
    flow.step().unwrap(); // Native streaming/bounce-back creates nonzero viscous stress.
    let before = flow.clone();
    let old_bits = flow.grid().state_bits();
    let j = [-0.0025, 0.00125, -0.0008];
    let report = flow
        .apply_impulse_density(&vec![j; 64], controls(), &CancelGate::new_clock_free())
        .unwrap();
    geometry_held(&before, &flow);
    let new_bits = flow.grid().state_bits();
    let mut maximum_stress = 0.0_f64;
    let mut kinetic = 0.0;
    for i in 0..64 {
        let old = moments(&old_bits, i, FORCE);
        let next = moments(&new_bits, i, FORCE);
        assert_eq!(next.0.to_bits(), old.0.to_bits());
        for a in 0..3 {
            close(old.0 * (next.1[a] - old.1[a]), j[a], 2e-15);
            kinetic += 0.5 * old.0 * (next.1[a].powi(2) - old.1[a].powi(2));
            for b in 0..3 {
                maximum_stress = maximum_stress.max(old.2[a][b].abs());
                close(next.2[a][b], old.2[a][b], 3e-14);
            }
        }
    }
    assert!(
        maximum_stress > 1e-6,
        "fixture needs actual nonequilibrium stress"
    );
    assert!(
        report.work < 0.0,
        "a braking impulse must report extracted work"
    );
    close(report.kinetic_energy_change, kinetic, 3e-14);
    close(report.work, kinetic, 3e-14);
}

#[test]
fn partial_cell_impulses_subdivide_reverse_and_preserve_zero_bits() {
    let initial = partial_fixture();
    let mut whole = initial.clone();
    let mut split = initial.clone();
    let mut impulses = vec![[0.0; 3]; 64];
    impulses[21] = [0.006, -0.003, 0.0015];
    let whole_report = whole
        .apply_impulse_density(&impulses, controls(), &CancelGate::new_clock_free())
        .unwrap();
    let half: Vec<_> = impulses.iter().map(|j| j.map(|v| 0.5 * v)).collect();
    let first = split
        .apply_impulse_density(&half, controls(), &CancelGate::new_clock_free())
        .unwrap();
    let second = split
        .apply_impulse_density(&half, controls(), &CancelGate::new_clock_free())
        .unwrap();
    close(first.work + second.work, whole_report.work, 2e-15);
    for (a, b) in whole
        .grid()
        .state_bits()
        .iter()
        .zip(split.grid().state_bits())
    {
        close(f64::from_bits(*a), f64::from_bits(b), 3e-15);
    }
    let reversed: Vec<_> = impulses.iter().map(|j| j.map(|v| -v)).collect();
    let reverse_report = whole
        .apply_impulse_density(&reversed, controls(), &CancelGate::new_clock_free())
        .unwrap();
    close(whole_report.work + reverse_report.work, 0.0, 2e-15);
    for (a, b) in whole
        .grid()
        .state_bits()
        .iter()
        .zip(initial.grid().state_bits())
    {
        close(f64::from_bits(*a), f64::from_bits(b), 3e-15);
    }
    let accepted = whole.clone();
    let zero = whole
        .apply_impulse_density(
            &vec![[0.0; 3]; 64],
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(zero.affected_cells, 0);
    assert_eq!(zero.work, 0.0);
    assert_eq!(whole.grid().state_bits(), accepted.grid().state_bits());
    geometry_held(&initial, &whole);
    geometry_held(&initial, &split);
}

#[test]
fn sparse_unsupported_targets_and_late_thin_cell_refusals_publish_nothing() {
    let mut sparse = FreeSurface3::with_solid_tiles(
        16,
        4,
        4,
        0.83,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        &[(3, 0, 0)],
        |x, _, _| x < 2,
    )
    .unwrap();
    assert_eq!(sparse.grid().active_tiles(), 2);
    let original = sparse.clone();
    // x=3 is active gas, x=8 is inactive, x=12 is a permanent solid tile.
    for x in [3, 8, 12] {
        let mut impulses = vec![[0.0; 3]; 256];
        impulses[0] = [0.001, 0.0, 0.0];
        impulses[x] = [0.001, 0.0, 0.0];
        assert!(matches!(
            sparse.apply_impulse_density(&impulses, controls(), &CancelGate::new_clock_free()),
            Err(ImpulseError::Invalid(_))
        ));
        assert_eq!(sparse.grid().state_bits(), original.grid().state_bits());
        geometry_held(&original, &sparse);
    }
    let mut too_small = controls();
    too_small.max_cells = 255;
    assert!(
        sparse
            .apply_impulse_density(
                &vec![[0.0; 3]; 256],
                too_small,
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    assert!(
        sparse
            .apply_impulse_density(
                &vec![[0.0; 3]; 255],
                controls(),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );

    let mut thin = seeded(
        &[
            FreeSurfaceInitialCell3 {
                coordinates: [0; 3],
                fill_fraction: 1.0,
                density: 1.0,
                velocity: [0.0; 3],
            },
            FreeSurfaceInitialCell3 {
                coordinates: [3; 3],
                fill_fraction: 2.0_f64.powi(-20),
                density: 1.0,
                velocity: [0.0; 3],
            },
        ],
        [0.0; 3],
    );
    let original = thin.clone();
    let mut impulses = vec![[0.0; 3]; 64];
    impulses[0] = [0.001, 0.0, 0.0];
    impulses[63] = [0.001, 0.0, 0.0];
    assert!(matches!(
        thin.apply_impulse_density(&impulses, controls(), &CancelGate::new_clock_free()),
        Err(ImpulseError::Cell { index: 63 })
    ));
    assert_eq!(thin.grid().state_bits(), original.grid().state_bits());
    geometry_held(&original, &thin);

    impulses[0] = [0.0; 3];
    impulses[63] = [mass_at(&thin, [3; 3]) * 0.04, 0.0, 0.0];
    thin.apply_impulse_density(&impulses, controls(), &CancelGate::new_clock_free())
        .unwrap();
    close(thin.sample_at(3, 3, 3).unwrap().velocity[0], 0.04, 2e-15);
    geometry_held(&original, &thin);
    let accepted = thin.clone();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert_eq!(
        thin.apply_impulse_density(&impulses, controls(), &gate),
        Err(ImpulseError::Cancelled)
    );
    impulses[63][1] = f64::NAN;
    assert!(matches!(
        thin.apply_impulse_density(&impulses, controls(), &CancelGate::new_clock_free()),
        Err(ImpulseError::Cell { index: 63 })
    ));
    assert_eq!(thin.grid().state_bits(), accepted.grid().state_bits());
    geometry_held(&accepted, &thin);
}
