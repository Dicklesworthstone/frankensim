use super::*;
use crate::compressible::spatial::mixture::molecular::WallAcceleration3d;

const MU: f64 = 0.75;
const ACCELERATION: f64 = 2.0;
const END: f64 = 0.4;

fn couette() -> (Coupled, ForcingWorkspace3d) {
    let grid = CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let flow = field(
        equal_calorics(),
        grid,
        Boundaries3d {
            y: WALLS,
            ..PERIODIC
        },
        |_| ([1.0, 0.0, 0.0], [0.0; 3], 1.0),
    );
    let mut molecular = prepare(&flow, law(flow.gas(), 0.001, MU, 0.0, 0.0));
    for side in [BoundarySide::YMin, BoundarySide::YMax] {
        molecular = molecular
            .with_wall(
                side,
                MolecularWall3d::NoSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                },
            )
            .unwrap();
    }
    coupled(flow, molecular)
}

fn driven_box() -> (Coupled, ForcingWorkspace3d) {
    let grid = CartesianGrid::new(3, 2, 2, 0.4, 0.5, 0.6).unwrap();
    let flow = field(
        equal_mass_gas(),
        grid,
        Boundaries3d {
            x: WALLS,
            y: WALLS,
            z: WALLS,
        },
        |_| ([0.2, 0.5, 0.3], [0.12, -0.04, 0.08], 1.2),
    );
    let mut molecular = prepare(&flow, law(flow.gas(), 0.001, 0.02, 0.005, 0.03));
    for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
        let mut velocity = [0.04, 0.07, -0.03];
        velocity[side.axis()] = 0.0;
        molecular = molecular
            .with_wall(
                side,
                MolecularWall3d::NoSlip {
                    velocity,
                    thermal: WallThermalCondition::Isothermal(Temperature::new(
                        1.1 + 0.03 * i as f64,
                    )),
                },
            )
            .unwrap();
    }
    coupled(flow, molecular)
}

fn top_drive() -> [WallAcceleration3d; 1] {
    [WallAcceleration3d {
        side: BoundarySide::YMax,
        acceleration: [ACCELERATION, 0.0, 0.0],
    }]
}

fn walls(flow: &Coupled) -> [MolecularWall3d; 6] {
    BoundarySide::ALL.map(|side| flow.molecular().wall(side))
}

fn assert_endpoints(
    initial: [MolecularWall3d; 6],
    flow: &Coupled,
    drive: &[WallAcceleration3d],
    elapsed: f64,
) {
    let mut expected = initial;
    for item in drive {
        let index = BoundarySide::ALL
            .iter()
            .position(|s| *s == item.side)
            .unwrap();
        if let MolecularWall3d::NoSlip { velocity, .. } = &mut expected[index] {
            for d in 0..3 {
                if item.acceleration[d] != 0.0 {
                    velocity[d] = item.acceleration[d].mul_add(elapsed, velocity[d]);
                }
            }
        }
    }
    assert_eq!(walls(flow), expected);
}

fn all_drives() -> [WallAcceleration3d; 6] {
    BoundarySide::ALL.map(|side| {
        let mut acceleration = [2.0, -3.0, 4.0];
        acceleration[side.axis()] = 0.0;
        WallAcceleration3d { side, acceleration }
    })
}

/// Independent exact solution of the one-cell native spatial discretization:
/// rho=H=A=1, lower wall stationary, upper speed a*t. The two half-cell
/// stresses give u' = lambda*(a*t/2-u), lambda=4*mu. This tests temporal
/// convergence, not continuum shear-profile accuracy on a one-cell mesh.
fn exact() -> [f64; 3] {
    let lambda = 4.0 * MU;
    let a = ACCELERATION;
    let t = END;
    let exp = (-lambda * t).exp();
    let velocity = 0.5 * a * (t - (1.0 - exp) / lambda);
    // Top-wall power = (2*mu)*(a*t-u)*(a*t); the stationary wall does no work.
    let work = MU
        * a
        * a
        * (t.powi(3) / 3.0 + t * t / (2.0 * lambda)
            - (1.0 - (1.0 + lambda * t) * exp) / lambda.powi(3));
    [velocity, work, work - 0.5 * velocity * velocity]
}

fn couette_error(dt: f64, composed: bool) -> [f64; 3] {
    let (mut owner, _) = couette();
    let bytes = owner.workspace_bytes();
    let initial = owner.flow().totals().energy;
    let mut work = 0.0;
    let (velocity, energy) = if composed {
        while owner.flow().time().value() < END {
            let duration = (END - owner.flow().time().value()).min(dt);
            let report = owner
                .step_with_wall_acceleration(
                    Time::new(duration),
                    MolecularStepControl::default(),
                    &top_drive(),
                    || false,
                )
                .unwrap();
            check(&report);
            assert_eq!(report.coupling_retries, 0);
            work += report.boundary_change.energy;
        }
        assert_eq!(owner.workspace_bytes(), bytes);
        close(owner.flow().time().value(), END, 1e-14);
        let end_wall = MolecularWall3d::NoSlip {
            velocity: [ACCELERATION * END, 0.0, 0.0],
            thermal: WallThermalCondition::Adiabatic,
        };
        if let (
            MolecularWall3d::NoSlip { velocity, .. },
            MolecularWall3d::NoSlip {
                velocity: expected, ..
            },
        ) = (owner.molecular().wall(BoundarySide::YMax), end_wall)
        {
            close(velocity[0], expected[0], 1e-14);
        } else {
            panic!("no-slip wall");
        }
        (
            owner.flow().cell(0).unwrap().momentum[0],
            owner.flow().totals().energy,
        )
    } else {
        let (mut flow, mut molecular) = owner.into_parts();
        let storage = molecular.workspace_bytes() + flow.workspace_bytes();
        let report = molecular
            .advance_with_wall_acceleration(
                &mut flow,
                Time::new(END),
                DissipationControl {
                    maximum_step_s: dt,
                    ..DissipationControl::default()
                },
                &top_drive(),
                || false,
            )
            .unwrap();
        super::super::check(report);
        assert_eq!(report.rejected_steps, 0);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(
            molecular.workspace_bytes() + flow.workspace_bytes(),
            storage
        );
        assert_eq!(
            molecular.wall(BoundarySide::YMax),
            MolecularWall3d::NoSlip {
                velocity: [ACCELERATION * END, 0.0, 0.0],
                thermal: WallThermalCondition::Adiabatic,
            }
        );
        work = report.boundary_change.energy;
        (flow.cell(0).unwrap().momentum[0], flow.totals().energy)
    };
    close(energy - initial, work, 3e-12);
    let heat = energy - initial - 0.5 * velocity * velocity;
    assert!(work > 0.0 && heat > 0.0);
    let reference = exact();
    std::array::from_fn(|i| ([velocity, work, heat][i] - reference[i]).abs())
}

#[test]
fn g1_accelerating_couette_momentum_work_and_heating_have_second_order_time_ladders() {
    for composed in [false, true] {
        let errors = [0.005, 0.0025, 0.00125].map(|dt| couette_error(dt, composed));
        eprintln!("accelerating Couette composed={composed}: {errors:?}");
        for component in 0..3 {
            for pair in errors.windows(2) {
                let ratio = pair[0][component] / pair[1][component];
                assert!(
                    ratio > 3.6 && ratio < 4.8,
                    "{composed} {component}: {ratio}"
                );
            }
            assert!(errors[2][component] < 8e-5);
        }
    }
}

#[test]
fn g2_accelerating_walls_and_body_sources_keep_separate_complete_ledgers() {
    let (mut owner, mut scratch) = driven_box();
    let initial = walls(&owner);
    let bytes = owner.workspace_bytes();
    let force = [0.2, -0.1, 0.15];
    let report = owner
        .step_with_forcing_and_wall_acceleration(
            Time::new(0.002),
            MolecularStepControl::default(),
            BodyForce3d::UniformForceDensity(force),
            CAP,
            &mut scratch,
            &all_drives(),
            || false,
        )
        .unwrap();
    check(&report);
    final_admission(&owner, &report, CAP);
    assert_endpoints(initial, &owner, &all_drives(), report.elapsed.value());
    assert_eq!(owner.workspace_bytes(), bytes);
    let volume = owner.flow().grid().cell_volume() * owner.flow().cell_count() as f64;
    let mut outward = [0.0; 7];
    for load in loads(&owner) {
        for d in 0..7 {
            outward[d] += load[d];
        }
    }
    for axis in 0..3 {
        close(
            report.source_change.impulse[axis],
            volume * force[axis] * report.elapsed.value(),
            3e-12,
        );
        close(
            report.after.momentum[axis] - report.before.momentum[axis] + outward[3 + axis]
                - report.source_change.impulse[axis],
            0.0,
            3e-12,
        );
    }
    close(
        report.after.energy - report.before.energy + outward[6] - report.source_change.work,
        0.0,
        3e-12,
    );
    assert!(report.boundary_change.energy.abs() > 1e-8);
    assert!(report.source_change.work.abs() > 1e-8);
}

#[test]
fn g3_wall_acceleration_retries_replay_from_original_wall_state_and_clock() {
    let (mut owner, mut scratch) = driven_box();
    let initial = walls(&owner);
    let force = BodyForce3d::UniformAcceleration([400.0, -200.0, 100.0]);
    let cap = ForceControl3d { maximum_mach: 0.8 };
    let report = owner
        .step_with_forcing_and_wall_acceleration(
            Time::new(0.2),
            MolecularStepControl::default(),
            force,
            cap,
            &mut scratch,
            &all_drives(),
            || false,
        )
        .unwrap();
    check(&report);
    assert!(report.coupling_retries > 0);
    assert_endpoints(initial, &owner, &all_drives(), report.elapsed.value());
    let (mut replay, mut replay_scratch) = driven_box();
    let expected = replay
        .step_with_forcing_and_wall_acceleration(
            report.elapsed,
            MolecularStepControl::default(),
            force,
            cap,
            &mut replay_scratch,
            &all_drives(),
            || false,
        )
        .unwrap();
    assert_eq!(expected.coupling_retries, 0);
    assert_eq!(report.first, expected.first);
    assert_eq!(report.advection, expected.advection);
    assert_eq!(report.second, expected.second);
    assert_eq!(owner.flow().states(), replay.flow().states());
    assert_eq!(owner.flow().time(), replay.flow().time());
    assert_eq!(walls(&owner), walls(&replay));
    assert_eq!(loads(&owner), loads(&replay));
    assert_eq!(molecular_loads(&owner), molecular_loads(&replay));
}

#[test]
fn g4_driven_wall_cancellation_and_work_refusal_preserve_all_accepted_owners() {
    let (mut owner, _) = driven_box();
    let (mut reference, _) = driven_box();
    for o in [&mut owner, &mut reference] {
        o.step_with_wall_acceleration(
            Time::new(0.001),
            MolecularStepControl::default(),
            &all_drives(),
            || false,
        )
        .unwrap();
    }
    let mut polls = 0;
    let expected = reference
        .step_with_wall_acceleration(
            Time::new(0.002),
            MolecularStepControl::default(),
            &all_drives(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let states = owner.flow().states().to_vec();
    let time = owner.flow().time();
    let old = *owner.last_step().unwrap();
    let wave = *owner.flow().last_step().unwrap();
    let molecular = *owner.molecular().last_step().unwrap();
    let old_loads = loads(&owner);
    let old_molecular_loads = molecular_loads(&owner);
    let old_walls = walls(&owner);
    let bytes = owner.workspace_bytes();
    let mut count = 0;
    assert_eq!(
        owner.step_with_wall_acceleration(
            Time::new(0.002),
            MolecularStepControl::default(),
            &all_drives(),
            || {
                count += 1;
                count >= polls
            },
        ),
        Err(DissipationError::Cancelled)
    );
    assert!(
        owner
            .step_with_wall_acceleration(
                Time::new(0.001),
                MolecularStepControl {
                    molecular: DissipationControl {
                        maximum_step_s: 1e-12,
                        maximum_steps: 1,
                        ..DissipationControl::default()
                    },
                    maximum_coupling_retries: 0,
                    ..MolecularStepControl::default()
                },
                &all_drives(),
                || false,
            )
            .is_err()
    );
    assert_eq!(owner.flow().states(), states);
    assert_eq!(owner.flow().time(), time);
    assert_eq!(owner.last_step(), Some(&old));
    assert_eq!(owner.flow().last_step(), Some(&wave));
    assert_eq!(owner.molecular().last_step(), Some(&molecular));
    assert_eq!(loads(&owner), old_loads);
    assert_eq!(molecular_loads(&owner), old_molecular_loads);
    assert_eq!(walls(&owner), old_walls);
    assert_eq!(owner.workspace_bytes(), bytes);
    let replay = owner
        .step_with_wall_acceleration(
            Time::new(0.002),
            MolecularStepControl::default(),
            &all_drives(),
            || false,
        )
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(owner.flow().states(), reference.flow().states());
    assert_eq!(walls(&owner), walls(&reference));
    assert_eq!(loads(&owner), loads(&reference));

    let (mut flow, mut operator) = owner.into_parts();
    let states = flow.states().to_vec();
    let walls = BoundarySide::ALL.map(|side| operator.wall(side));
    let report = *operator.last_step().unwrap();
    assert!(
        operator
            .advance_with_wall_acceleration(
                &mut flow,
                Time::new(0.01),
                DissipationControl {
                    maximum_step_s: 0.001,
                    maximum_steps: 1,
                    ..DissipationControl::default()
                },
                &all_drives(),
                || false,
            )
            .is_err()
    );
    assert_eq!(flow.states(), states);
    assert_eq!(operator.last_step(), Some(&report));
    assert_eq!(BoundarySide::ALL.map(|side| operator.wall(side)), walls);
}

#[test]
fn g0_invalid_wall_histories_refuse_and_zero_acceleration_keeps_held_bits() {
    let (mut owner, _) = couette();
    let initial = owner.flow().states().to_vec();
    let initial_walls = walls(&owner);
    for drives in [
        vec![WallAcceleration3d {
            side: BoundarySide::YMax,
            acceleration: [0.0, 1.0, 0.0],
        }],
        vec![WallAcceleration3d {
            side: BoundarySide::YMax,
            acceleration: [f64::NAN, 0.0, 0.0],
        }],
        vec![WallAcceleration3d {
            side: BoundarySide::XMin,
            acceleration: [0.0, 1.0, 0.0],
        }],
        vec![top_drive()[0]; 2],
        vec![top_drive()[0]; 7],
    ] {
        assert!(
            owner
                .step_with_wall_acceleration(
                    Time::new(0.002),
                    MolecularStepControl::default(),
                    &drives,
                    || false
                )
                .is_err()
        );
        assert_eq!(owner.flow().states(), initial);
        assert_eq!(walls(&owner), initial_walls);
        assert!(owner.last_step().is_none());
        assert_eq!(owner.flow().time().value(), 0.0);
    }
    let (mut held, _) = driven_box();
    let (mut driven, _) = driven_box();
    let zero = BoundarySide::ALL.map(|side| WallAcceleration3d {
        side,
        acceleration: [0.0, -0.0, 0.0],
    });
    let expected = held
        .step(Time::new(0.002), MolecularStepControl::default(), || false)
        .unwrap();
    let observed = driven
        .step_with_wall_acceleration(
            Time::new(0.002),
            MolecularStepControl::default(),
            &zero,
            || false,
        )
        .unwrap();
    assert_eq!(observed, expected);
    for (a, b) in held.flow().states().iter().zip(driven.flow().states()) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    assert_eq!(loads(&held), loads(&driven));
    assert_eq!(walls(&held), walls(&driven));
    let (mut h, mut hm) = held.into_parts();
    let (mut d, mut dm) = driven.into_parts();
    let expected = hm
        .advance(
            &mut h,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let observed = dm
        .advance_with_wall_acceleration(
            &mut d,
            Time::new(0.001),
            DissipationControl::default(),
            &zero,
            || false,
        )
        .unwrap();
    assert_eq!(observed, expected);
    for (a, b) in h.states().iter().zip(d.states()) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn g3_native_source_halving_reuses_original_wall_substep_time() {
    let prepare_case = || {
        let (owner, _) = couette();
        let (flow, _) = owner.into_parts();
        let gas = flow.gas();
        let data = ConstantTransport::new(
            gas.source().clone(),
            "accelerating-wall-temperature-window",
            1,
            TransportDomain {
                minimum_temperature_k: 0.001,
                maximum_temperature_k: 1.2,
                minimum_pressure_pa: 0.001,
                maximum_pressure_pa: 100.0,
            },
            TransportCoefficients {
                dynamic_viscosity: DynViscosity::new(MU),
                bulk_viscosity: DynViscosity::new(0.0),
                thermal_conductivity: ThermalConductivityQuantityV1::new(0.0),
            },
            vec![0.0, 0.001, 0.001, 0.001, 0.0, 0.001, 0.001, 0.001, 0.0],
        )
        .unwrap();
        let mut operator = prepare(&flow, data);
        for (side, speed) in [(BoundarySide::YMin, 0.0), (BoundarySide::YMax, 10.0)] {
            operator = operator
                .with_wall(
                    side,
                    MolecularWall3d::NoSlip {
                        velocity: [speed, 0.0, 0.0],
                        thermal: WallThermalCondition::Adiabatic,
                    },
                )
                .unwrap();
        }
        (flow, operator)
    };
    let drive = [WallAcceleration3d {
        side: BoundarySide::YMax,
        acceleration: [-2000.0, 0.0, 0.0],
    }];
    let (mut flow, mut operator) = prepare_case();
    // The full Euler predictor overheats beyond the authored transport-law
    // window. Two half-size SSP updates remain inside it as the wall slows.
    let report = operator
        .advance_with_wall_acceleration(
            &mut flow,
            Time::new(0.005),
            DissipationControl {
                maximum_step_s: 0.005,
                ..DissipationControl::default()
            },
            &drive,
            || false,
        )
        .unwrap();
    super::super::check(report);
    assert!(report.rejected_steps > 0);
    assert_eq!(report.accepted_steps, 2);
    let (mut replay, mut replay_operator) = prepare_case();
    let expected = replay_operator
        .advance_with_wall_acceleration(
            &mut replay,
            Time::new(0.005),
            DissipationControl {
                maximum_step_s: 0.0025,
                ..DissipationControl::default()
            },
            &drive,
            || false,
        )
        .unwrap();
    assert_eq!(expected.rejected_steps, 0);
    assert_eq!(flow.states(), replay.states());
    assert_eq!(report.after, expected.after);
    assert_eq!(report.boundary_change, expected.boundary_change);
    for side in BoundarySide::ALL {
        assert_eq!(operator.wall(side), replay_operator.wall(side));
        for (a, b) in operator
            .boundary_impulses(side)
            .zip(replay_operator.boundary_impulses(side))
        {
            assert_eq!(a.momentum, b.momentum);
            assert_eq!(a.energy, b.energy);
        }
    }
}
