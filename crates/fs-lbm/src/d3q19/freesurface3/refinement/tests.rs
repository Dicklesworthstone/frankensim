//! G0/G1/G3/G4: retain the actual coarse inventory, force frame and stress.

use super::*;
use crate::d3q19::{W3, collide_cell3};

fn controls() -> DyadicRefinementControls3 {
    DyadicRefinementControls3 {
        max_cells: 4096,
        max_active_tiles: 64,
        maximum_speed: 0.15,
        maximum_mass_defect: 1e-11,
        maximum_first_moment_defect: 1e-10,
        maximum_momentum_defect: 1e-11,
        maximum_angular_momentum_defect: 1e-10,
        maximum_kinetic_energy_defect: 1e-12,
        maximum_viscous_stress_defect: 2e-13,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn fixture(boost: [f64; 3]) -> FreeSurface3 {
    let mut flow = FreeSurface3::new(
        4,
        4,
        4,
        0.83,
        [0.0; 3],
        0.001,
        ContactModel3::Neutral,
        |_, _, _| true,
    )
    .unwrap();
    let key = flow.grid.active_keys()[0];
    let mut forces = BTreeMap::new();
    for lane in 0..64 {
        let (x, y, z) = FreeSurface3::global_of(key, lane);
        let r = [x as f64 - 1.5, y as f64 - 1.5, z as f64 - 1.5];
        let rotation = cross([0.003, -0.002, 0.004], r);
        let u = core::array::from_fn(|a| [0.025, -0.015, 0.01][a] + boost[a] + rotation[a]);
        let rho = 0.97 + lane as f64 * 0.001;
        let force = [0.0001 + lane as f64 * 1e-6, -0.00015, 0.0002];
        let mut f = equilibrium3(rho, u);
        for q in 0..Q3 {
            f[q] += 9.0 * W3[q] * 0.0004 * f64::from(E3[q].0 * E3[q].1);
        }
        f = forcing::reframe(f, u, [0.0; 3], force).unwrap();
        // An actual native BGK-evolved state, not equilibrium-only copying.
        // No streaming or full VOF timestep is claimed by this unit fixture.
        f = collide_cell3(f, CollisionModel3::Bgk { tau: 0.83 }, force).unwrap();
        flow.grid.set_populations(0, lane, f);
        let density = f.iter().sum::<f64>();
        let fill = match lane % 5 {
            0 => 0.2,
            1 => 0.5,
            _ => 1.0,
        };
        let tile = flow.tiles.get_mut(&key).unwrap();
        tile.cells[lane] = if fill < 1.0 {
            Cell3::Interface
        } else {
            Cell3::Fluid
        };
        tile.mass[lane] = density * fill;
        forces.insert((key, lane), force);
    }
    flow.local_force = Some(forces);
    flow
}

// Independent original public sampling, with explicitly physical cell volumes
// and positions. This does not reuse the refinement receipt or accumulator.
fn observed(flow: &FreeSurface3, width: f64, mass_unit: f64, origin: [f64; 3]) -> [f64; 11] {
    let mut result = [0.0; 11];
    for (&key, tile) in &flow.tiles {
        for lane in 0..TILE_CELLS {
            if tile.cells[lane] == Cell3::Gas {
                continue;
            }
            let (x, y, z) = FreeSurface3::global_of(key, lane);
            let sample = flow.sample_at(x as usize, y as usize, z as usize).unwrap();
            let mass = mass_unit
                * if sample.cell == Cell3::Fluid {
                    sample.density
                } else {
                    tile.mass[lane]
                };
            let position =
                core::array::from_fn(|a| width * ([x, y, z][a] as f64 + 0.5) - origin[a]);
            let momentum = sample.velocity.map(|u| mass * u);
            let angular = cross(position, momentum);
            result[0] += mass;
            for a in 0..3 {
                result[1 + a] += mass * position[a];
                result[4 + a] += momentum[a];
                result[7 + a] += angular[a];
                result[10] += 0.5 * mass * sample.velocity[a].powi(2);
            }
        }
    }
    result
}

#[test]
fn evolved_partial_liquid_retains_physical_moments_pressure_and_viscous_stress() {
    for boost in [[0.0; 3], [0.017, -0.009, 0.005]] {
        let source = fixture(boost);
        let bits = source.grid.state_bits();
        let before = observed(&source, 1.0, 1.0, [-2.25, 3.5, 0.75]);
        let (target, report) = source
            .refine_dyadic(controls(), &CancelGate::new_clock_free())
            .unwrap();
        let after = observed(&target, 0.5, 0.125, [-2.25, 3.5, 0.75]);
        assert_eq!(source.grid.state_bits(), bits);
        assert_eq!(target.domain_cells(), [8; 3]);
        assert_eq!(target.grid.active_tiles(), 8 * source.grid.active_tiles());
        assert_eq!(target.geometry_epoch(), source.geometry_epoch() + 1);
        assert_eq!(target.conversions(), source.conversions());
        assert_eq!(target.carry(), 0.0);
        assert_eq!(target.sigma, 2.0 * source.sigma);
        assert_eq!(target.collision_model(), CollisionModel3::Bgk { tau: 1.16 });
        for i in 0..11 {
            near(after[i], before[i], 2e-11);
        }
        near(report.before.mass, before[0], 2e-13);
        near(report.after.mass, after[0], 2e-12);
        near(report.kinetic_energy_change, after[10] - before[10], 3e-15);
        let mut maximum_stress = 0.0_f64;
        let source_key = source.grid.active_keys()[0];
        for lane in 0..64 {
            let (x, y, z) = FreeSurface3::global_of(source_key, lane);
            let original = source
                .sample_at(x as usize, y as usize, z as usize)
                .unwrap();
            let old_force = source.force_at_lane(source_key, lane);
            let old = source.grid.populations(0, lane);
            let (rho, old_u) = state(&old, old_force, 0.15).unwrap();
            let old_stress = force_corrected_stress3(&old, rho, old_u, old_force);
            for child in 0..8 {
                let p = [
                    2 * x + (child & 1) as i64,
                    2 * y + ((child >> 1) & 1) as i64,
                    2 * z + (child >> 2) as i64,
                ];
                let sample = target
                    .sample_at(p[0] as usize, p[1] as usize, p[2] as usize)
                    .unwrap();
                assert_eq!(sample.density.to_bits(), original.density.to_bits());
                assert_eq!(
                    sample.fill_fraction.to_bits(),
                    original.fill_fraction.to_bits()
                );
                assert_eq!(sample.cell, original.cell);
                let (slot, cl) = target.grid.resolve_source(p[0], p[1], p[2]).unwrap();
                let ck = target.grid.active_keys()[slot];
                let force = target.force_at_lane(ck, cl);
                assert_eq!(force, old_force.map(|f| f * 0.5));
                let f = target.grid.populations(slot, cl);
                let (_, u) = state(&f, force, 0.15).unwrap();
                let stress = force_corrected_stress3(&f, rho, u, force);
                for a in 0..3 {
                    near(u[a], old_u[a], 2e-15);
                }
                for a in 0..6 {
                    maximum_stress = maximum_stress.max(old_stress[a].abs());
                    near(
                        stress[a] * (1.16 - 0.5) / 1.16,
                        old_stress[a] * (0.83 - 0.5) / 0.83,
                        3e-16,
                    );
                }
            }
        }
        assert!(maximum_stress > 1e-5);
    }
}

#[test]
fn affine_acceleration_is_sampled_at_new_centers_without_an_impulse() {
    let mut source = fixture([0.0; 3]);
    let field = BodyAccelerationField3::Affine {
        origin: [1.25, 2.0, -0.75],
        acceleration: [0.0002, -0.0001, 0.00015],
        gradient: [[1e-5, 2e-5, 0.0], [0.0, -1e-5, 1e-5], [2e-5, 0.0, 1e-5]],
    };
    let key = source.grid.active_keys()[0];
    for lane in 0..64 {
        let old = source.grid.populations(0, lane);
        let force = source.force_at_lane(key, lane);
        let (rho, u) = state(&old, force, 0.15).unwrap();
        let next = field.at_cell(key, lane).unwrap().map(|a| rho * a);
        source
            .grid
            .set_populations(0, lane, forcing::reframe(old, u, force, next).unwrap());
        source
            .local_force
            .as_mut()
            .unwrap()
            .insert((key, lane), next);
    }
    source.body_acceleration = Some(field);
    let before = observed(&source, 1.0, 1.0, [0.0; 3]);
    let (target, report) = source
        .refine_dyadic(controls(), &CancelGate::new_clock_free())
        .unwrap();
    let after = observed(&target, 0.5, 0.125, [0.0; 3]);
    for i in 0..11 {
        near(after[i], before[i], 2e-11);
    }
    let mut differs_from_parent_hold = false;
    for (&key, tile) in &target.tiles {
        for lane in 0..64 {
            if tile.cells[lane] == Cell3::Gas {
                continue;
            }
            let (x, y, z) = FreeSurface3::global_of(key, lane);
            let sample = target
                .sample_at(x as usize, y as usize, z as usize)
                .unwrap();
            let position = [
                x as f64 * 0.5 + 0.25,
                y as f64 * 0.5 + 0.25,
                z as f64 * 0.5 + 0.25,
            ];
            let BodyAccelerationField3::Affine {
                origin,
                acceleration,
                gradient,
            } = field
            else {
                unreachable!()
            };
            let force = target.force_at_lane(key, lane);
            for a in 0..3 {
                let physical = acceleration[a]
                    + (0..3)
                        .map(|b| gradient[a][b] * (position[b] - origin[b]))
                        .sum::<f64>();
                near(force[a], 0.5 * sample.density * physical, 3e-19);
            }
            let parent_lane = (x as usize / 2) + 4 * (y as usize / 2) + 16 * (z as usize / 2);
            let held = source.force_at_lane(0, parent_lane).map(|f| f * 0.5);
            differs_from_parent_hold |= (force[0] - held[0]).abs() > 1e-7;
        }
    }
    assert!(differs_from_parent_hold);
    assert!(report.absolute_momentum_defect < 1e-12);
    assert!(report.maximum_viscous_stress_defect < 2e-15);
}

#[test]
fn sparse_fixed_solids_piston_and_contact_keep_their_physical_locations() {
    let mut source = FreeSurface3::with_solid_tiles(
        8,
        8,
        8,
        0.9,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        &[(1, 1, 1)],
        |x, y, z| x >= 2 && x < 4 && y < 4 && z < 4,
    )
    .unwrap();
    source.planar_piston = Some(PlanarPiston3 {
        wall: PlanarWall3::XMin,
        depth_cells: 1,
    });
    source.planar_contact = Some(PlanarContactAngle3 {
        wall: PlanarWall3::YMin,
        angle_radians: core::f64::consts::FRAC_PI_3,
        cotangent: 1.0 / 3.0_f64.sqrt(),
    });
    let (target, report) = source
        .refine_dyadic(controls(), &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(target.domain_cells(), [16; 3]);
    assert_eq!(target.grid.active_tiles(), source.grid.active_tiles() * 8);
    assert_eq!(target.solid.len(), 8 * source.solid.len());
    for x in 2..4 {
        for y in 2..4 {
            for z in 2..4 {
                assert!(target.solid.contains(&morton3(x, y, z)));
                assert!(!target.grid.is_active(x, y, z));
            }
        }
    }
    assert_eq!(target.planar_piston.unwrap().depth_cells, 2);
    assert_eq!(
        target.planar_contact.unwrap().angle_radians,
        source.planar_contact.unwrap().angle_radians
    );
    near(report.before.mass, report.after.mass, 1e-12);
    for i in 0..3 {
        near(
            report.before.first_mass_moment[i],
            report.after.first_mass_moment[i],
            1e-12,
        );
    }
}

#[test]
fn rejected_models_fields_inventory_limits_and_cancellation_keep_source_bits() {
    let base = fixture([0.0; 3]);
    let bits = base.grid.state_bits();
    let mut cases = Vec::new();
    let mut carry = base.clone();
    carry.carry = 0.01;
    cases.push(carry);
    let mut local = base.clone();
    local.local_sigma = Some(BTreeMap::new());
    cases.push(local);
    let mut bulk = base.clone();
    bulk.local_bulk_relaxation = Some(BTreeMap::new());
    cases.push(bulk);
    let mut collision = base.clone();
    collision
        .grid
        .set_collision_model(CollisionModel3::Trt {
            tau: 0.83,
            magic: 0.1875,
        })
        .unwrap();
    cases.push(collision);
    let mut fill = base.clone();
    fill.tiles.get_mut(&0).unwrap().mass[0] = 2.0;
    cases.push(fill);
    let mut capillary = base.clone();
    capillary.sigma = 0.05;
    cases.push(capillary);
    for source in cases {
        let original = source.grid.state_bits();
        assert!(
            source
                .refine_dyadic(controls(), &CancelGate::new_clock_free())
                .is_err()
        );
        assert_eq!(source.grid.state_bits(), original);
    }
    for limits in [
        DyadicRefinementControls3 {
            max_cells: 511,
            ..controls()
        },
        DyadicRefinementControls3 {
            max_active_tiles: 7,
            ..controls()
        },
        DyadicRefinementControls3 {
            maximum_viscous_stress_defect: 0.0,
            ..controls()
        },
        DyadicRefinementControls3 {
            maximum_speed: 0.001,
            ..controls()
        },
    ] {
        assert!(
            base.refine_dyadic(limits, &CancelGate::new_clock_free())
                .is_err()
        );
        assert_eq!(base.grid.state_bits(), bits);
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        base.refine_dyadic(controls(), &gate),
        Err(FreeSurfaceError3::Cancelled)
    ));
    assert_eq!(base.grid.state_bits(), bits);
}

#[test]
fn repeated_refinement_retains_the_same_physical_inventory_and_viscosity() {
    let source = fixture([0.0; 3]);
    let gate = CancelGate::new_clock_free();
    let (once, _) = source.refine_dyadic(controls(), &gate).unwrap();
    let (twice, _) = once.refine_dyadic(controls(), &gate).unwrap();
    let before = observed(&source, 1.0, 1.0, [1.25, -2.5, 0.75]);
    let after = observed(&twice, 0.25, 1.0 / 64.0, [1.25, -2.5, 0.75]);
    for i in 0..11 {
        near(after[i], before[i], 1e-10);
    }
    near(
        twice.collision_model().kinematic_viscosity() * 0.25,
        source.collision_model().kinematic_viscosity(),
        1e-16,
    );
    assert_eq!(twice.geometry_epoch(), source.geometry_epoch() + 2);
}

#[test]
fn retained_out_of_domain_solids_refuse_before_child_morton_encoding() {
    let gate = CancelGate::new_clock_free();
    let half_limit = 1u32 << (MORTON_COORD_BITS - 1);
    for axis in 0..3 {
        // An ordinary outside tile plus the first overflowing parent and the
        // highest representable source key all pass the original constructor.
        for coordinate in [1, half_limit, 2 * half_limit - 1] {
            let mut solid = [0u32; 3];
            solid[axis] = coordinate;
            let source = FreeSurface3::with_solid_tiles(
                4,
                4,
                4,
                0.83,
                [0.0; 3],
                0.0,
                ContactModel3::Neutral,
                &[(solid[0], solid[1], solid[2])],
                |_, _, _| true,
            )
            .unwrap();
            let original_bits = source.grid.state_bits();
            let original_solid = source.solid.clone();
            let original_mass = source.ledger_mass().to_bits();
            assert!(matches!(
                source.refine_dyadic(controls(), &gate),
                Err(FreeSurfaceError3::Transport(
                    "dyadic source solid is outside the physical or refinable Morton domain"
                ))
            ));
            assert_eq!(source.grid.state_bits(), original_bits);
            assert_eq!(source.solid, original_solid);
            assert_eq!(source.ledger_mass().to_bits(), original_mass);
            assert_eq!(source.geometry_epoch(), 0);
        }
    }
}
