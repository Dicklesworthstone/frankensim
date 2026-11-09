use super::*;
use fs_lbm::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;
use fs_lbm::lattice3::open_boundaries::forcing::resistance::{
    DarcyForchheimer3, LinearResistance3, ResistanceOpenGrid3,
};
use fs_lbm::thermal::boundary::{FiniteThermalWall, ThermalWall};
use fs_lbm::thermal::open::{ThermalOpenBoundaries3, ThermalOpenFace3, ThermalOpenLimits3};

const DIMS: [usize; 3] = [6, 5, 8];
const CAP: usize = 1024 * 1024;

fn make<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    upper: f64,
) -> OpenBoundaryGrid3<L, Q> {
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        DIMS,
        collision,
        [0.0; 3],
        [true, false, true],
        budget(),
        |_| false,
        |_| (1.0, [0.0; 3]),
    )
    .unwrap();
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[2] = FaceBoundary3::stationary_wall();
    faces[3] = FaceBoundary3::Wall {
        velocity: [upper, 0.0, 0.0],
    };
    OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
}

fn drive(point: [f64; 3], mean: f64) -> [f64; 3] {
    [
        mean + 0.01 * (core::f64::consts::TAU * point[2] / DIMS[2] as f64).cos(),
        0.0,
        0.0,
    ]
}

fn uniform<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let mut held = make::<L, Q>(cx, collision, 0.04);
        let mut traced = make::<L, Q>(cx, collision, 0.0);
        let gate = CancelGate::new_clock_free();
        let report = traced
            .set_wall_velocity_profile(Face3::YMax, CAP, &gate, |_| Ok([0.04, 0.0, 0.0]))
            .unwrap();
        assert!(report.sampled_links > DIMS[0] * DIMS[2]);
        assert_eq!(report.maximum_speed_lu, 0.04);
        for _ in 0..400 {
            assert_eq!(
                held.step_serial(cx, None, [0.0; 3]).unwrap(),
                traced.step_serial(cx, None, [0.0; 3]).unwrap()
            );
            assert_eq!(held.populations(), traced.populations());
            assert_eq!(held.last_wall_exchange(), traced.last_wall_exchange());
        }
        for y in 0..DIMS[1] {
            let u = traced.macroscopic(address([2, y, 3], DIMS)).unwrap().1;
            near(u[0], 0.04 * (y as f64 + 0.5) / DIMS[1] as f64, 2e-8);
            near(u[1], 0.0, 2e-13);
            near(u[2], 0.0, 2e-13);
        }
    });
}

#[test]
fn uniform_link_trace_preserves_native_bits_and_analytic_couette() {
    uniform::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    uniform::<D3q27, 27>(bgk());
}

fn orientations<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        for face in Face3::ALL {
            let axis = face as usize / 2;
            let mut target = [0.0; 3];
            target[(axis + 1) % 3] = 0.015;
            let make = |moving: bool| {
                let mut periodic = [true; 3];
                periodic[axis] = false;
                let grid = Grid3::<L, Q>::try_new_equilibrium(
                    [4, 3, 2],
                    collision,
                    [0.0; 3],
                    periodic,
                    budget(),
                    |_| false,
                    |_| (1.0, [0.0; 3]),
                )
                .unwrap();
                let mut faces = [FaceBoundary3::Periodic; 6];
                faces[2 * axis] = FaceBoundary3::stationary_wall();
                faces[2 * axis + 1] = FaceBoundary3::stationary_wall();
                if moving {
                    faces[face as usize] = FaceBoundary3::Wall { velocity: target };
                }
                OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
            };
            let mut held = make(true);
            let mut trace = make(false);
            trace
                .set_wall_velocity_profile(face, CAP, &gate, |_| Ok(target))
                .unwrap();
            for _ in 0..3 {
                held.step_serial(cx, None, [0.0; 3]).unwrap();
                trace.step_serial(cx, None, [0.0; 3]).unwrap();
                assert_eq!(held.populations(), trace.populations());
                assert_eq!(held.last_wall_exchange(), trace.last_wall_exchange());
            }
        }
    });
}

#[test]
fn all_six_wall_planes_preserve_original_uniform_motion() {
    orientations::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    orientations::<D3q27, 27>(bgk());
}

fn spatial<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let mut flow = make::<L, Q>(cx, collision, 0.0);
        let initial = flow.populations().to_vec();
        let origin = [0.25, -0.5, 1.75];
        flow.set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| Ok(drive(p, 0.03)))
            .unwrap();
        let reservoir = flow.step_serial(cx, None, origin).unwrap();
        assert_eq!(reservoir.open_links, 0);
        let receipt = flow.last_wall_exchange().unwrap();
        let mut expected_fluid = [0.0; 3];
        let mut expected_angular = [0.0; 3];
        let mut expected_wall_work = 0.0;
        for cell in 0..initial.len() {
            let p = [
                cell % DIMS[0],
                cell / DIMS[0] % DIMS[1],
                cell / (DIMS[0] * DIMS[1]),
            ];
            if p[1] + 1 != DIMS[1] {
                continue;
            }
            let mut cell_impulse = [0.0; 3];
            for q in 1..Q {
                if L::E[q][1] != -1 {
                    continue;
                }
                let e = L::E[q].map(f64::from);
                let point = core::array::from_fn(|a| p[a] as f64 + 0.5 - 0.5 * e[a]);
                let u = drive(point, 0.03);
                let outgoing = initial[cell][L::OPP[q]];
                let delta = 2.0 * L::W[q] * e[0] * u[0] / L::CS2;
                near(flow.populations()[cell][q], outgoing + delta, 4e-16);
                for a in 0..3 {
                    cell_impulse[a] += e[a] * delta;
                }
                let arm: [f64; 3] = core::array::from_fn(|a| point[a] - origin[a]);
                for a in 0..3 {
                    let b = (a + 1) % 3;
                    let c = (a + 2) % 3;
                    expected_angular[a] += (arm[b] * e[c] - arm[c] * e[b]) * delta;
                }
                expected_wall_work += (-e[0] * (2.0 * outgoing + delta) + u[0] * delta) * u[0];
            }
            let (rho, u) = flow.macroscopic(cell).unwrap();
            near(rho, 1.0, 2e-14);
            for a in 0..3 {
                near(u[a], cell_impulse[a], 2e-14);
                expected_fluid[a] += cell_impulse[a];
            }
        }
        for a in 0..3 {
            near(
                receipt.fluid_population_impulse[a],
                expected_fluid[a],
                2e-13,
            );
            near(
                receipt.fluid_population_angular_impulse[a],
                expected_angular[a],
                3e-12,
            );
            near(
                receipt.wall_impulse[a] + receipt.fluid_population_impulse[a],
                receipt.wall_velocity_mass_impulse[a],
                2e-13,
            );
            near(
                receipt.wall_angular_impulse[a] + receipt.fluid_population_angular_impulse[a],
                receipt.wall_velocity_mass_angular_impulse[a],
                3e-12,
            );
        }
        near(receipt.wall_work, expected_wall_work, 3e-14);
        near(receipt.fluid_mass_change, 0.0, 3e-14);
        assert!(receipt.wall_work < 0.0);
        // Time-varying traces use the identical original closed-grid link law.
        let mut closed = Grid3::<L, Q>::try_new_equilibrium(
            DIMS,
            collision,
            [0.0; 3],
            [true, false, true],
            budget(),
            |_| false,
            |_| (1.0, [0.0; 3]),
        )
        .unwrap();
        let probe = closed.wall_momentum_probe(cx, origin, |_| true).unwrap();
        let mut traced = make::<L, Q>(cx, collision, 0.0);
        for mean in [0.03, -0.02, 0.01, -0.03] {
            traced
                .set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| Ok(drive(p, mean)))
                .unwrap();
            let velocities: Vec<_> = probe
                .wall_links()
                .map(|link| {
                    if link.exterior_faces[3] {
                        drive(link.position, mean)
                    } else {
                        [0.0; 3]
                    }
                })
                .collect();
            let wall = closed
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
                .unwrap();
            traced.step_serial(cx, None, origin).unwrap();
            assert_eq!(traced.populations(), closed.populations());
            assert_eq!(traced.last_wall_exchange(), Some(wall));
        }
    });
}

#[test]
fn spatial_link_velocity_controls_actual_reflection_work_and_time_updates() {
    spatial::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    spatial::<D3q27, 27>(bgk());
}

fn open_seam<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let build = || {
            let grid = Grid3::<L, Q>::try_new_equilibrium(
                DIMS,
                collision,
                [0.0; 3],
                [false, false, true],
                budget(),
                |_| false,
                |_| (1.0, [0.02, 0.0, 0.0]),
            )
            .unwrap();
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[0] = FaceBoundary3::Velocity {
                velocity: [0.02, 0.0, 0.0],
            };
            faces[1] = FaceBoundary3::Pressure { density: 1.0 };
            faces[2] = FaceBoundary3::stationary_wall();
            faces[3] = FaceBoundary3::stationary_wall();
            let mut flow = OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap();
            flow.set_velocity_profile(
                Face3::XMin,
                fs_lbm::d3q19::InletFlux3::Unconstrained,
                &gate,
                |_| Ok([0.02, 0.0, 0.0]),
            )
            .unwrap();
            let inlet = flow.velocity_profile(Face3::XMin).unwrap().to_vec();
            flow.set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| Ok(drive(p, 0.03)))
                .unwrap();
            faces[1] = FaceBoundary3::Pressure { density: 1.001 };
            flow.rebind_open_boundaries(BoundarySpec3::new(faces), &gate)
                .unwrap();
            assert_eq!(flow.velocity_profile(Face3::XMin).unwrap(), inlet);
            flow
        };
        let mut serial = build();
        let mut pooled = build();
        let initial = serial.populations()[0];
        let before = totals(&serial);
        let open = serial.step_serial(cx, None, [0.0; 3]).unwrap();
        let parallel = pooled
            .step_pooled(&TilePool::for_host(2, SEED), &gate, None, [0.0; 3])
            .unwrap();
        assert_eq!(open, parallel);
        assert_eq!(serial.populations(), pooled.populations());
        assert_eq!(serial.last_wall_exchange(), pooled.last_wall_exchange());
        let wall = serial.last_wall_exchange().unwrap();
        let after = totals(&serial);
        near(after.0 - before.0, open.mass_change, 4e-12);
        near(wall.fluid_mass_change, 0.0, 3e-14);
        for a in 0..3 {
            near(
                after.1[a] - before.1[a],
                open.impulse[a] + wall.fluid_population_impulse[a],
                4e-12,
            );
        }
        for x in [0, DIMS[0] - 1] {
            let p = [x, DIMS[1] - 1, 3];
            let cell = address(p, DIMS);
            for q in 1..Q {
                if L::E[q][1] != -1 {
                    continue;
                }
                let e = L::E[q].map(f64::from);
                let point = core::array::from_fn(|a| p[a] as f64 + 0.5 - 0.5 * e[a]);
                near(
                    serial.populations()[cell][q],
                    initial[L::OPP[q]] + 2.0 * L::W[q] * e[0] * drive(point, 0.03)[0] / L::CS2,
                    5e-16,
                );
            }
        }
    });
}

#[test]
fn wall_trace_preserves_open_seam_ownership_and_serial_pooled_identity() {
    open_seam::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    open_seam::<D3q27, 27>(bgk());
}

#[test]
fn wall_trace_refuses_invalid_sources_without_partial_control_or_state() {
    with_cx(|cx| {
        let mut flow = make::<D3q27, 27>(cx, bgk(), 0.0);
        let gate = CancelGate::new_clock_free();
        let report = flow
            .set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| Ok(drive(p, 0.03)))
            .unwrap();
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        let old = flow.wall_velocity_profile(Face3::YMax).unwrap().to_vec();
        let populations = flow.populations().to_vec();
        let receipt = flow.last_wall_exchange();
        for mode in 0..5 {
            let mut sampled = 0;
            let result = flow.set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| {
                sampled += 1;
                if mode == 4 && sampled > 3 {
                    return Err(GridError3::Input("source refused"));
                }
                Ok(match mode {
                    0 => [0.03, 1e-20, 0.0],
                    1 => [f64::NAN, 0.0, 0.0],
                    2 => [0.5, 0.0, 0.0],
                    3 => [0.001 * p[0], 0.0, 0.0], // Tangential but nonzero local wall mass.
                    _ => drive(p, 0.02),
                })
            });
            assert!(result.is_err());
            assert_eq!(flow.wall_velocity_profile(Face3::YMax).unwrap(), old);
        }
        assert!(
            flow.set_wall_velocity_profile(Face3::YMax, 0, &gate, |_| Ok([0.0; 3]))
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        let mut count = 0;
        assert!(matches!(
            flow.set_wall_velocity_profile(Face3::YMax, CAP, &cancelled, |p| {
                count += 1;
                if count == report.sampled_links {
                    cancelled.request();
                }
                Ok(drive(p, -0.02))
            }),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(flow.wall_velocity_profile(Face3::YMax).unwrap(), old);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(flow.last_wall_exchange(), receipt);
        let mut faces = core::array::from_fn(|i| flow.boundaries().face(Face3::ALL[i]));
        faces[3] = FaceBoundary3::Wall {
            velocity: [0.01, 0.0, 0.0],
        };
        assert!(
            flow.rebind_open_boundaries(BoundarySpec3::new(faces), &gate)
                .is_err()
        );
        faces[3] = FaceBoundary3::stationary_wall();
        faces[2] = FaceBoundary3::Wall {
            velocity: [-0.01, 0.0, 0.0],
        };
        flow.rebind_open_boundaries(BoundarySpec3::new(faces), &gate)
            .unwrap();
        assert_eq!(flow.wall_velocity_profile(Face3::YMax).unwrap(), old);
        flow.clear_wall_velocity_profile(Face3::YMax, &gate)
            .unwrap();
        assert!(flow.wall_velocity_profile(Face3::YMax).is_none());
    });
}

fn heat_faces() -> ThermalOpenBoundaries3 {
    let mut faces = [ThermalOpenFace3::Periodic; 6];
    faces[2] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
    faces[3] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
    ThermalOpenBoundaries3 {
        faces,
        internal_wall: ThermalWall::Adiabatic,
    }
}
fn heat_limits() -> ThermalOpenLimits3 {
    ThermalOpenLimits3 {
        maximum_cells: 512,
        maximum_storage_bytes: CAP,
        temperature_range: [0.0, 2.0],
        balance_tolerance: 1e-10,
    }
}

fn porous_heat<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let cells = DIMS.iter().product::<usize>();
        let law = DarcyForchheimer3::from(LinearResistance3 {
            principal: [0.2; 3],
            ..LinearResistance3::default()
        });
        let mut flow = ResistanceOpenGrid3::try_new_preserving_velocity(
            make::<L, Q>(cx, collision, 0.0),
            cx,
            vec![law; cells],
            [0.0; 3],
        )
        .unwrap();
        flow.set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| Ok(drive(p, 0.03)))
            .unwrap();
        let mut heat = flow
            .prepare_thermal_field(cx, &vec![0.8; cells], 0.9, heat_faces(), heat_limits())
            .unwrap();
        let result = flow
            .step_serial_with_thermal(&mut heat, cx, None, &vec![0.006; cells], None)
            .unwrap();
        assert_eq!(
            result.momentum.walls,
            flow.flow().last_wall_exchange().unwrap()
        );
        assert_eq!(flow.flow().completed_steps(), heat.completed_steps());
        near(result.heat.defect, 0.0, 1e-12);
        for cell in 0..cells {
            let p = [
                cell % DIMS[0],
                cell / DIMS[0] % DIMS[1],
                cell / (DIMS[0] * DIMS[1]),
            ];
            let mut raw = [0.0; 3];
            if p[1] + 1 == DIMS[1] {
                for q in 1..Q {
                    if L::E[q][1] != -1 {
                        continue;
                    }
                    let e = L::E[q].map(f64::from);
                    let point = core::array::from_fn(|a| p[a] as f64 + 0.5 - 0.5 * e[a]);
                    let delta = 2.0 * L::W[q] * e[0] * drive(point, 0.03)[0] / L::CS2;
                    for a in 0..3 {
                        raw[a] += e[a] * delta;
                    }
                }
            }
            let (rho, u) = flow.flow().macroscopic(cell).unwrap();
            near(rho, 1.0, 3e-14);
            for a in 0..3 {
                near(u[a], raw[a] / 1.1, 3e-14);
            }
            near(heat.temperature(cell).unwrap(), 0.806, 3e-14);
            let g = heat.populations()[cell];
            for a in 0..3 {
                near(g[1 + 2 * a] - g[2 + 2 * a], 0.003 * u[a], 3e-14);
            }
        }
        // A profile cannot bypass the original moving finite-wall refusal.
        flow.clear_wall_velocity_profile(Face3::YMax, &gate)
            .unwrap();
        let mut faces = heat_faces();
        faces.faces[3] = ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
            temperature: 0.8,
            heat_capacity: 2.0,
            resistance: 0.0,
            heat_per_step: 0.0,
        }));
        let mut finite = flow
            .prepare_thermal_field(cx, &vec![0.8; cells], 0.9, faces, heat_limits())
            .unwrap();
        flow.set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| Ok(drive(p, 0.02)))
            .unwrap();
        let before = (
            flow.flow().populations().to_vec(),
            flow.flow().forces().to_vec(),
            finite.populations().to_vec(),
            finite.wall_states(3).unwrap().to_vec(),
            flow.last_exchange(),
        );
        assert!(
            flow.step_serial_with_thermal(&mut finite, cx, None, &vec![0.0; cells], None)
                .is_err()
        );
        assert_eq!(flow.flow().populations(), before.0);
        assert_eq!(flow.flow().forces(), before.1);
        assert_eq!(finite.populations(), before.2);
        assert_eq!(finite.wall_states(3).unwrap(), before.3);
        assert_eq!(flow.last_exchange(), before.4);
        assert_eq!(flow.flow().completed_steps(), 1);
        assert_eq!(finite.completed_steps(), 1);
    });
}

#[test]
fn porous_and_thermal_consumers_share_resolved_wall_velocity_and_refusals() {
    porous_heat::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    porous_heat::<D3q27, 27>(bgk());
}

#[test]
fn forced_trace_matches_uniform_force_owner_without_detaching_force_epoch() {
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let mut direct = make::<D3q27, 27>(cx, bgk(), 0.0);
        let mut forced = ForcedOpenBoundaryGrid3::try_new_preserving_velocity(
            make::<D3q27, 27>(cx, bgk(), 0.0),
            cx,
            vec![[0.0; 3]; DIMS.iter().product()],
        )
        .unwrap();
        for mean in [0.03, -0.02] {
            direct
                .set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| Ok(drive(p, mean)))
                .unwrap();
            forced
                .set_wall_velocity_profile(Face3::YMax, CAP, &gate, |p| Ok(drive(p, mean)))
                .unwrap();
            assert_eq!(
                direct.step_serial(cx, None, [0.0; 3]).unwrap(),
                forced.step_serial(cx, None, [0.0; 3]).unwrap()
            );
            assert_eq!(direct.populations(), forced.populations());
            assert_eq!(direct.last_wall_exchange(), forced.last_wall_exchange());
        }
        forced
            .clear_wall_velocity_profile(Face3::YMax, &gate)
            .unwrap();
        assert!(forced.wall_velocity_profile(Face3::YMax).is_none());
    });
}
