use super::*;
use crate::compressible::spatial::mixture::molecular::{
    MolecularEuler3d, MolecularStepControl, MolecularStepReport3d,
};

type Coupled = MolecularEuler3d<ConstantTransport>;
const CAP: ForceControl3d = ForceControl3d { maximum_mach: 5.0 };

fn equal_calorics() -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, 40.0),
        source("B", R, 3.5, 0.0, 40.0),
        source("absent-C", R, 4.5, 0.0, -20.0),
    ])
}

fn coupled(
    flow: MixtureEuler3d,
    molecular: MolecularTransport3d<ConstantTransport>,
) -> (Coupled, ForcingWorkspace3d) {
    let scratch = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let bytes =
        Coupled::required_workspace_bytes(&flow, &molecular).unwrap() + scratch.workspace_bytes();
    let flow = Coupled::new(flow, molecular, bytes, || false).unwrap();
    assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
    (flow, scratch)
}

fn check(report: &MolecularStepReport3d) {
    assert_eq!(report.elapsed, report.advection.elapsed);
    assert_eq!(report.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(report.second.elapsed, report.first.elapsed);
    assert_eq!(report.source_change, report.advection.source_change);
    for value in report.conservation_residual.species_masses {
        close(value, 0.0, 3e-11);
    }
    for value in report.conservation_residual.momentum {
        close(value, 0.0, 3e-11);
    }
    close(report.conservation_residual.energy, 0.0, 3e-11);
    close(report.absolute_energy_residual_j, 0.0, 3e-9);
}

fn loads(flow: &Coupled) -> Vec<[f64; 7]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            flow.boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.species_masses[2],
                    q.momentum[0],
                    q.momentum[1],
                    q.momentum[2],
                    q.energy,
                ]
            })
        })
        .collect()
}

fn molecular_loads(flow: &Coupled) -> Vec<[f64; 7]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            flow.molecular().boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.species_masses[2],
                    q.momentum[0],
                    q.momentum[1],
                    q.momentum[2],
                    q.energy,
                ]
            })
        })
        .collect()
}

fn final_admission(flow: &Coupled, report: &MolecularStepReport3d, cap: ForceControl3d) {
    let mut rate: f64 = 0.0;
    for index in 0..flow.flow().cell_count() {
        let p = flow
            .flow()
            .cell(index)
            .unwrap()
            .primitive(flow.flow().gas(), || false)
            .unwrap();
        assert!(
            p.velocity[0].hypot(p.velocity[1]).hypot(p.velocity[2])
                <= cap.maximum_mach * p.sound_speed * (1.0 + 1e-13)
        );
        rate = rate.max(
            (0..3)
                .map(|axis| {
                    (p.velocity[axis].abs() + p.sound_speed) / flow.flow().grid().spacing()[axis]
                })
                .sum(),
        );
    }
    assert!(
        rate * report.elapsed.value()
            <= MolecularStepControl::default().advection.cfl * (1.0 + 1e-13)
    );
}

fn wave(n: usize) -> f64 {
    let h = 1.0 / n as f64;
    let z = std::f64::consts::PI * h;
    let average = (z.sin() / z).powi(3);
    let initial_velocity = [0.3, -0.2, 0.15];
    let acceleration = [-0.4, 0.25, -0.3];
    let diffusion = 0.02;
    let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
    let flow = field(equal_calorics(), grid, PERIODIC, |x| {
        let a = 0.5 + 0.2 * average * (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2])).sin();
        ([a, 1.0 - a, 0.0], initial_velocity, 1.0)
    });
    let mass = flow.totals().mass();
    let molecular = prepare(&flow, law(flow.gas(), diffusion, 0.015, 0.005, 0.02));
    let (mut flow, mut scratch) = coupled(flow, molecular);
    let bytes = flow.workspace_bytes();
    let forces = vec![acceleration; grid.len()];
    let duration = 0.04;
    let mut impulse = [0.0; 3];
    let mut work = 0.0;
    while flow.flow().time().value() < duration {
        let report = flow
            .step_with_forcing(
                Time::new(duration - flow.flow().time().value()),
                MolecularStepControl::default(),
                BodyForce3d::Acceleration(&forces),
                CAP,
                &mut scratch,
                || false,
            )
            .unwrap();
        check(&report);
        for axis in 0..3 {
            impulse[axis] += report.source_change.impulse[axis];
        }
        work += report.source_change.work;
    }
    assert_eq!(flow.workspace_bytes(), bytes);
    close(flow.flow().time().value(), duration, 1e-14);
    let displacement: [f64; 3] = std::array::from_fn(|axis| {
        duration * initial_velocity[axis] + 0.5 * duration * duration * acceleration[axis]
    });
    let velocity: [f64; 3] =
        std::array::from_fn(|axis| initial_velocity[axis] + duration * acceleration[axis]);
    let decay = (-diffusion * 12.0 * std::f64::consts::PI.powi(2) * duration).exp();
    let mut error = 0.0;
    for index in 0..grid.len() {
        let x = grid.cell_center(index).unwrap();
        let phase =
            2.0 * std::f64::consts::PI * (0..3).map(|a| x[a] - displacement[a]).sum::<f64>();
        let expected = 0.5 + 0.2 * average * decay * phase.sin();
        let q = flow.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        error += (q.partial_densities[0] - expected).abs();
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.density, 1.0, 2e-11);
        close(p.temperature, 1.0, 2e-10);
        close(p.pressure, 1.0, 2e-10);
        for axis in 0..3 {
            close(p.velocity[axis], velocity[axis], 2e-11);
        }
    }
    for axis in 0..3 {
        close(impulse[axis], mass * duration * acceleration[axis], 2e-11);
    }
    close(
        work,
        0.5 * mass
            * (0..3)
                .map(|axis| velocity[axis].powi(2) - initial_velocity[axis].powi(2))
                .sum::<f64>(),
        2e-11,
    );
    assert!(work < 0.0);
    error / grid.len() as f64
}

#[test]
fn g1_accelerating_three_axis_mixture_diffusion_matches_physical_wave() {
    let coarse = wave(8);
    let fine = wave(16);
    eprintln!("coupled accelerating spatial mixture L1 8/16: {coarse} / {fine}");
    assert!(fine < 0.5 * coarse);
    assert!(fine < 2e-3);
}

fn wall_case() -> (Coupled, ForcingWorkspace3d) {
    let original = field(
        equal_mass_gas(),
        CartesianGrid::new(3, 2, 2, 0.2, 0.25, 0.3).unwrap(),
        Boundaries3d {
            x: WALLS,
            y: WALLS,
            z: WALLS,
        },
        |x| {
            let rho = 1.0 + 0.1 * x[2];
            let a = 0.4 + 0.1 * x[0] - 0.05 * x[1] + 0.05 * x[2];
            (
                [rho * a, rho * (1.0 - a), 0.0],
                [0.1 + 0.05 * x[1], -0.15 + 0.03 * x[2], 0.05],
                1.0 + 0.1 * x[2],
            )
        },
    );
    let nominal = original.workspace_bytes();
    let mut states = Vec::with_capacity(original.states.len() + 31);
    states.extend_from_slice(&original.states);
    let flow = MixtureEuler3d::new(
        original.gas,
        original.grid,
        states,
        original.boundaries,
        MAX_MIXTURE_WORKSPACE_BYTES,
        || false,
    )
    .unwrap();
    assert_eq!(flow.workspace_bytes(), nominal + 31 * size_of::<f64>());
    let molecular = prepare(&flow, law(flow.gas(), 0.02, 0.03, 0.01, 0.04))
        .with_wall(
            BoundarySide::ZMin,
            MolecularWall3d::NoSlip {
                velocity: [0.2, -0.1, 0.0],
                thermal: WallThermalCondition::Isothermal(Temperature::new(1.2)),
            },
        )
        .unwrap()
        .with_wall(
            BoundarySide::ZMax,
            MolecularWall3d::NoSlip {
                velocity: [-0.1, 0.2, 0.0],
                thermal: WallThermalCondition::Isothermal(Temperature::new(0.85)),
            },
        )
        .unwrap();
    coupled(flow, molecular)
}

#[test]
fn g2_body_wall_impulse_heat_and_work_close_in_all_three_directions() {
    let (mut flow, mut scratch) = wall_case();
    let bytes = flow.workspace_bytes();
    let walls = BoundarySide::ALL.map(|side| flow.molecular().wall(side));
    let force = [0.2, -0.1, 0.15];
    let dt = 0.002;
    let report = flow
        .step_with_forcing(
            Time::new(dt),
            MolecularStepControl::default(),
            BodyForce3d::UniformForceDensity(force),
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
    check(&report);
    final_admission(&flow, &report, CAP);
    assert_eq!(report.elapsed.value(), dt);
    assert_eq!(flow.flow().time(), report.elapsed);
    assert_eq!(flow.workspace_bytes(), bytes);
    assert_eq!(flow.molecular().last_step(), Some(&report.second));
    assert_eq!(
        BoundarySide::ALL.map(|side| flow.molecular().wall(side)),
        walls
    );
    assert!(report.first.boundary_change.energy.abs() > 1e-8);
    close(report.advection.boundary_change.energy, 0.0, 2e-12);
    let mut outward = [0.0; 7];
    for row in loads(&flow) {
        for d in 0..7 {
            outward[d] += row[d];
        }
    }
    let volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
    for axis in 0..3 {
        close(
            report.source_change.impulse[axis],
            volume * dt * force[axis],
            3e-12,
        );
        close(
            report.after.momentum[axis] - report.before.momentum[axis] + outward[3 + axis]
                - report.source_change.impulse[axis],
            0.0,
            3e-12,
        );
        close(
            report.boundary_change.momentum[axis],
            report.first.boundary_change.momentum[axis]
                + report.advection.boundary_change.momentum[axis]
                + report.second.boundary_change.momentum[axis],
            3e-12,
        );
    }
    for k in 0..3 {
        close(outward[k], 0.0, 2e-12);
        close(
            report.after.species_masses[k],
            report.before.species_masses[k],
            3e-12,
        );
    }
    close(
        report.after.absolute_energy(flow.flow().gas()).unwrap()
            - report.before.absolute_energy(flow.flow().gas()).unwrap()
            + outward[6]
            - report.source_change.work,
        0.0,
        3e-11,
    );
}

#[test]
fn g3_forced_spatial_retry_replays_complete_molecular_and_wave_interval() {
    let (mut flow, mut scratch) = wall_case();
    let acceleration = [400.0, -200.0, 100.0];
    let forces = vec![acceleration; flow.flow().cell_count()];
    let cap = ForceControl3d { maximum_mach: 0.8 };
    let bytes = flow.workspace_bytes();
    let report = flow
        .step_with_forcing(
            Time::new(0.2),
            MolecularStepControl::default(),
            BodyForce3d::Acceleration(&forces),
            cap,
            &mut scratch,
            || false,
        )
        .unwrap();
    check(&report);
    final_admission(&flow, &report, cap);
    assert!(report.coupling_retries > 0);
    assert_eq!(flow.flow().time(), report.elapsed);
    assert_eq!(flow.workspace_bytes(), bytes);
    for axis in 0..3 {
        close(
            report.source_change.impulse[axis],
            report.before.mass() * acceleration[axis] * report.elapsed.value(),
            3e-12,
        );
    }
    let (mut replay, mut replay_scratch) = wall_case();
    let expected = replay
        .step_with_forcing(
            report.elapsed,
            MolecularStepControl::default(),
            BodyForce3d::Acceleration(&forces),
            cap,
            &mut replay_scratch,
            || false,
        )
        .unwrap();
    assert_eq!(expected.coupling_retries, 0);
    assert_eq!(expected.first, report.first);
    assert_eq!(expected.advection, report.advection);
    assert_eq!(expected.second, report.second);
    assert_eq!(expected.source_change, report.source_change);
    assert_eq!(replay.flow().states(), flow.flow().states());
    assert_eq!(loads(&replay), loads(&flow));
    assert_eq!(molecular_loads(&replay), molecular_loads(&flow));
}

#[test]
fn g4_spatial_refusal_preserves_both_owners_and_complete_wall_transfers() {
    let (mut flow, mut scratch) = wall_case();
    let (mut reference, mut reference_scratch) = wall_case();
    let force = BodyForce3d::UniformForceDensity([0.2, -0.1, 0.15]);
    for (flow, scratch) in [
        (&mut flow, &mut scratch),
        (&mut reference, &mut reference_scratch),
    ] {
        flow.step_with_forcing(
            Time::new(0.001),
            MolecularStepControl::default(),
            force,
            CAP,
            scratch,
            || false,
        )
        .unwrap();
    }
    let mut polls = 0;
    let expected = reference
        .step_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            force,
            CAP,
            &mut reference_scratch,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let states = flow.flow().states().to_vec();
    let time = flow.flow().time();
    let old = *flow.last_step().unwrap();
    let wave = *flow.flow().last_step().unwrap();
    let molecular = *flow.molecular().last_step().unwrap();
    let old_loads = loads(&flow);
    let old_molecular_loads = molecular_loads(&flow);
    let bytes = flow.workspace_bytes();
    let mut count = 0;
    assert_eq!(
        flow.step_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            force,
            CAP,
            &mut scratch,
            || {
                count += 1;
                count >= polls
            },
        ),
        Err(DissipationError::Cancelled)
    );
    assert!(
        flow.step_with_forcing(
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
            force,
            CAP,
            &mut scratch,
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.flow().time(), time);
    assert_eq!(flow.last_step(), Some(&old));
    assert_eq!(flow.flow().last_step(), Some(&wave));
    assert_eq!(flow.molecular().last_step(), Some(&molecular));
    assert_eq!(loads(&flow), old_loads);
    assert_eq!(molecular_loads(&flow), old_molecular_loads);
    assert_eq!(flow.workspace_bytes(), bytes);
    let replay = flow
        .step_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            force,
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(flow.flow().states(), reference.flow().states());
    assert_eq!(loads(&flow), loads(&reference));
    assert_eq!(molecular_loads(&flow), molecular_loads(&reference));
}

#[test]
fn g0_zero_source_keeps_spatial_native_bits_and_counts_all_owned_storage() {
    let (mut native, _) = wall_case();
    let (mut flow, mut scratch) = wall_case();
    let expected = native
        .step(Time::new(0.002), MolecularStepControl::default(), || false)
        .unwrap();
    let zero = vec![[0.0; 3]; flow.flow().cell_count()];
    let report = flow
        .step_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            BodyForce3d::Acceleration(&zero),
            ForceControl3d { maximum_mach: 1e-6 },
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(report, expected);
    for (left, right) in flow.flow().states().iter().zip(native.flow().states()) {
        assert_eq!(left.to_bits(), right.to_bits());
    }
    assert_eq!(loads(&flow), loads(&native));
    assert_eq!(molecular_loads(&flow), molecular_loads(&native));
    assert_eq!(flow.molecular().last_step(), native.molecular().last_step());

    let (owner, _) = wall_case();
    let (original, molecular) = owner.into_parts();
    let bytes = Coupled::required_workspace_bytes(&original, &molecular).unwrap();
    assert!(matches!(
        Coupled::new(original, molecular, bytes - 1, || false),
        Err(DissipationError::Invalid(_))
    ));
    let (owner, _) = wall_case();
    let (original, molecular) = owner.into_parts();
    let mut tight = Coupled::new(original, molecular, bytes, || false).unwrap();
    let old = tight.flow().states().to_vec();
    assert!(matches!(
        tight.step_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            BodyForce3d::Acceleration(&zero),
            CAP,
            &mut scratch,
            || false,
        ),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(tight.flow().states(), old);
    assert_eq!(tight.flow().time().value(), 0.0);
    assert_eq!(tight.last_step(), None);
    assert_eq!(tight.molecular().last_step(), None);
    tight
        .step(Time::new(0.002), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(tight.flow().states(), native.flow().states());
}

#[test]
fn g1_cold_walls_cannot_bypass_the_final_spatial_force_mach_guard() {
    let make = || {
        let flow = field(
            equal_calorics(),
            CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
            Boundaries3d {
                x: AxisBoundaries::Periodic,
                y: AxisBoundaries::Periodic,
                z: WALLS,
            },
            |_| ([0.5, 0.5, 0.0], [0.9, 0.0, 0.0], 2.0),
        );
        let mut molecular = prepare(&flow, law(flow.gas(), 0.02, 0.0, 0.0, 100.0));
        for side in [BoundarySide::ZMin, BoundarySide::ZMax] {
            molecular = molecular
                .with_wall(
                    side,
                    MolecularWall3d::NoSlip {
                        velocity: [0.9, 0.0, 0.0],
                        thermal: WallThermalCondition::Isothermal(Temperature::new(1.0)),
                    },
                )
                .unwrap();
        }
        (flow, molecular)
    };
    let dt = 0.001;
    let cap = ForceControl3d { maximum_mach: 0.55 };
    let force = BodyForce3d::UniformAcceleration([0.01, 0.0, 0.0]);
    // The original component owners independently isolate final cooling: the
    // central forced interval is admitted, then the second half exceeds Mach.
    let (mut manual, mut operator) = make();
    let mut manual_scratch = manual.forcing_workspace(usize::MAX).unwrap();
    operator
        .advance(
            &mut manual,
            Time::new(0.5 * dt),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let middle = manual
        .step_with_forcing(
            Time::new(dt),
            StepControl::default(),
            force,
            cap,
            &mut manual_scratch,
            || false,
        )
        .unwrap();
    assert_eq!(middle.elapsed.value(), dt);
    operator
        .advance(
            &mut manual,
            Time::new(0.5 * dt),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let p = manual
        .cell(0)
        .unwrap()
        .primitive(manual.gas(), || false)
        .unwrap();
    assert!(
        p.velocity[0].hypot(p.velocity[1]).hypot(p.velocity[2]) > cap.maximum_mach * p.sound_speed
    );
    let (flow, operator) = make();
    let (mut flow, mut scratch) = coupled(flow, operator);
    let old = flow.flow().states().to_vec();
    assert!(
        flow.step_with_forcing(
            Time::new(dt),
            MolecularStepControl {
                maximum_coupling_retries: 0,
                ..MolecularStepControl::default()
            },
            force,
            cap,
            &mut scratch,
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.flow().states(), old);
    assert_eq!(flow.flow().time().value(), 0.0);
    assert_eq!(flow.last_step(), None);
    assert_eq!(flow.molecular().last_step(), None);
    let report = flow
        .step_with_forcing(
            Time::new(dt),
            MolecularStepControl::default(),
            force,
            cap,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert!(report.coupling_retries > 0);
    check(&report);
    final_admission(&flow, &report, cap);
}
