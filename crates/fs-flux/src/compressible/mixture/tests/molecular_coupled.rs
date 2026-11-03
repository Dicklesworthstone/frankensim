use super::molecular::{law, prepare};
use super::reservoir::bind;
use super::*;
use crate::compressible::mixture::molecular::{
    DissipationControl, DissipationError, MolecularEuler2d, MolecularStepControl,
    MolecularStepReport,
};
use fs_gas::multicomponent::dissipation::ConstantTransport;

type Coupled = MolecularEuler2d<ConstantTransport>;
const CAP: ForceControl2d = ForceControl2d { maximum_mach: 5.0 };

fn equal_calorics() -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, 40.0),
        source("B", R, 3.5, 0.0, 40.0),
        source("absent-C", R, 4.5, 0.0, -20.0),
    ])
}

fn coupled(
    flow: MixtureEuler2d,
    diffusion: f64,
    shear: f64,
    bulk: f64,
    conductivity: f64,
) -> (Coupled, ForcingWorkspace2d) {
    let scratch = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let molecular = prepare(&flow, law(flow.gas(), diffusion, shear, bulk, conductivity));
    let bytes =
        Coupled::required_workspace_bytes(&flow, &molecular).unwrap() + scratch.workspace_bytes();
    let flow = Coupled::new(flow, molecular, bytes, || false).unwrap();
    assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
    (flow, scratch)
}

fn check(report: &MolecularStepReport) {
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

fn loads(flow: &Coupled) -> Vec<[f64; 6]> {
    SIDES
        .into_iter()
        .flat_map(|side| {
            flow.boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.species_masses[2],
                    q.momentum[0],
                    q.momentum[1],
                    q.energy,
                ]
            })
        })
        .collect()
}

fn molecular_loads(flow: &Coupled) -> Vec<[f64; 6]> {
    SIDES
        .into_iter()
        .flat_map(|side| {
            flow.molecular().boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.species_masses[2],
                    q.momentum[0],
                    q.momentum[1],
                    q.energy,
                ]
            })
        })
        .collect()
}

fn final_admission(flow: &Coupled, report: &MolecularStepReport, cap: ForceControl2d) {
    let mut rate: f64 = 0.0;
    for index in 0..flow.flow().cell_count() {
        let p = flow
            .flow()
            .cell(index)
            .unwrap()
            .primitive(flow.flow().gas(), || false)
            .unwrap();
        assert!(
            p.velocity[0].hypot(p.velocity[1]) <= cap.maximum_mach * p.sound_speed * (1.0 + 1e-13)
        );
        rate = rate.max(
            (0..2)
                .map(|axis| {
                    (p.velocity[axis].abs() + p.sound_speed) / flow.flow().grid().spacing[axis]
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
    let average = (z.sin() / z).powi(2);
    let initial_velocity = [0.35, -0.2];
    let acceleration = [-0.8, 0.4];
    let diffusion = 0.02;
    let grid = CartesianGrid::new(n, n, h, h, 0.6).unwrap();
    let flow = field(equal_calorics(), grid, PERIODIC, |x| {
        let a = 0.5 + 0.2 * average * (2.0 * std::f64::consts::PI * (x[0] + x[1])).sin();
        ([a, 1.0 - a, 0.0], initial_velocity, 1.0)
    });
    let mass = flow.totals().mass();
    let (mut flow, mut scratch) = coupled(flow, diffusion, 0.015, 0.005, 0.02);
    let bytes = flow.workspace_bytes();
    let forces = vec![acceleration; grid.len()];
    let duration = 0.08;
    let mut impulse = [0.0; 2];
    let mut work = 0.0;
    while flow.flow().time().value() < duration {
        let remaining = duration - flow.flow().time().value();
        let report = flow
            .step_with_forcing(
                Time::new(remaining),
                MolecularStepControl::default(),
                BodyForce2d::Acceleration(&forces),
                CAP,
                &mut scratch,
                || false,
            )
            .unwrap();
        check(&report);
        for axis in 0..2 {
            impulse[axis] += report.source_change.impulse[axis];
        }
        work += report.source_change.work;
    }
    assert_eq!(flow.workspace_bytes(), bytes);
    close(flow.flow().time().value(), duration, 1e-14);
    let displacement: [f64; 2] = std::array::from_fn(|axis| {
        duration * initial_velocity[axis] + 0.5 * duration * duration * acceleration[axis]
    });
    let velocity: [f64; 2] =
        std::array::from_fn(|axis| initial_velocity[axis] + duration * acceleration[axis]);
    let decay = (-diffusion * 8.0 * std::f64::consts::PI.powi(2) * duration).exp();
    let mut error = 0.0;
    for index in 0..grid.len() {
        let x = grid.cell_center(index).unwrap();
        let phase = 2.0 * std::f64::consts::PI * (x[0] + x[1] - displacement[0] - displacement[1]);
        let expected = 0.5 + 0.2 * average * decay * phase.sin();
        let q = flow.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        error += (q.partial_densities[0] - expected).abs();
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.density, 1.0, 2e-11);
        close(p.temperature, 1.0, 2e-10);
        close(p.pressure, 1.0, 2e-10);
        for axis in 0..2 {
            close(p.velocity[axis], velocity[axis], 2e-11);
        }
    }
    for axis in 0..2 {
        close(impulse[axis], mass * duration * acceleration[axis], 2e-11);
    }
    close(
        work,
        0.5 * mass
            * (0..2)
                .map(|axis| velocity[axis].powi(2) - initial_velocity[axis].powi(2))
                .sum::<f64>(),
        2e-11,
    );
    assert!(work < 0.0, "the held source brakes this moving mixture");
    error / grid.len() as f64
}

#[test]
fn g1_accelerating_multicomponent_advection_diffusion_matches_physical_wave() {
    let coarse = wave(12);
    let fine = wave(24);
    eprintln!("coupled accelerating mixture L1 12/24: {coarse} / {fine}");
    assert!(fine < 0.4 * coarse);
    assert!(fine < 1e-3);
}

fn reservoir_field(excess: usize) -> MixtureEuler2d {
    let gas = gas(vec![
        source("A", R, 3.5, 0.0, 40.0),
        source("B", R, 4.5, 0.0, -50.0),
        source("absent-C", R, 3.5, 0.0, 2.0),
    ]);
    let inputs = [
        ReservoirBoundary {
            side: BoundarySide::XMin,
            partial_densities: &[0.8, 0.2, 0.0],
            velocity: [0.3, 0.4],
            pressure_pa: 1.1,
        },
        ReservoirBoundary {
            side: BoundarySide::XMax,
            partial_densities: &[0.2, 0.8, 0.0],
            velocity: [0.3, -0.2],
            pressure_pa: 0.9,
        },
    ];
    let original = field(
        gas,
        CartesianGrid::new(3, 2, 0.2, 0.3, 0.4).unwrap(),
        Boundaries2d {
            x: OPEN,
            y: AxisBoundaries::Periodic,
        },
        |x| {
            let a = 0.5 + 0.15 * x[0];
            (
                [a, 1.0 - a, 0.0],
                [0.25, 0.1 + 0.05 * x[1]],
                1.0 + 0.1 * x[1],
            )
        },
    );
    let mut states = Vec::with_capacity(original.states.len() + excess);
    states.extend_from_slice(&original.states);
    MixtureEuler2d::new_with_reservoirs(
        original.gas,
        original.grid,
        states,
        original.boundaries,
        &inputs,
        MAX_MIXTURE_WORKSPACE_BYTES,
        || false,
    )
    .unwrap()
}

#[test]
fn g2_reservoir_and_body_transfers_close_over_the_complete_molecular_interval() {
    let (mut flow, mut scratch) = coupled(reservoir_field(31), 0.02, 0.03, 0.01, 0.04);
    let initial_bytes = flow.workspace_bytes();
    let force = [0.2, -0.1];
    let dt = 0.002;
    let report = flow
        .step_with_forcing(
            Time::new(dt),
            MolecularStepControl::default(),
            BodyForce2d::UniformForceDensity(force),
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
    check(&report);
    final_admission(&flow, &report, CAP);
    assert_eq!(report.elapsed.value(), dt);
    assert_eq!(flow.flow().time(), report.elapsed);
    assert_eq!(flow.molecular().last_step(), Some(&report.second));
    assert_eq!(flow.workspace_bytes(), initial_bytes);
    assert!(report.first.boundary_change.species_masses[0].abs() > 1e-8);
    assert!(report.first.boundary_change.energy.abs() > 1e-8);
    let mut outward = [0.0; 6];
    for row in loads(&flow) {
        for d in 0..6 {
            outward[d] += row[d];
        }
    }
    let volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
    for axis in 0..2 {
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
    let references = flow.flow().gas().energy_references();
    let boundary_energy = outward[5] + (0..3).map(|k| outward[k] * references[k]).sum::<f64>();
    close(
        report.after.absolute_energy(flow.flow().gas()).unwrap()
            - report.before.absolute_energy(flow.flow().gas()).unwrap()
            + boundary_energy
            - report.source_change.work,
        0.0,
        3e-11,
    );
    for k in 0..3 {
        close(
            report.after.species_masses[k] - report.before.species_masses[k] + outward[k],
            0.0,
            3e-12,
        );
    }
    assert_eq!(
        flow.flow().reservoir(BoundarySide::XMin).unwrap().velocity,
        [0.3, 0.4]
    );
    assert_eq!(
        flow.flow()
            .reservoir(BoundarySide::XMax)
            .unwrap()
            .partial_densities,
        [0.2, 0.8, 0.0]
    );
}

#[test]
fn g3_forced_retry_replays_both_molecular_halves_at_one_duration() {
    let make = || coupled(reservoir_field(31), 0.02, 0.03, 0.01, 0.04);
    let (mut flow, mut scratch) = make();
    let forces = vec![[400.0, -200.0]; flow.flow().cell_count()];
    let cap = ForceControl2d { maximum_mach: 0.8 };
    let bytes = flow.workspace_bytes();
    let report = flow
        .step_with_forcing(
            Time::new(0.2),
            MolecularStepControl::default(),
            BodyForce2d::Acceleration(&forces),
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
    let (mut replay, mut replay_scratch) = make();
    let expected = replay
        .step_with_forcing(
            report.elapsed,
            MolecularStepControl::default(),
            BodyForce2d::Acceleration(&forces),
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
fn g4_final_cancellation_and_work_refusal_preserve_both_published_owners() {
    let make = || coupled(reservoir_field(31), 0.02, 0.03, 0.01, 0.04);
    let (mut flow, mut scratch) = make();
    let (mut reference, mut reference_scratch) = make();
    let force = BodyForce2d::UniformForceDensity([0.2, -0.1]);
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
fn g0_zero_source_keeps_native_bits_but_still_admits_complete_storage() {
    let (mut native, _) = coupled(reservoir_field(31), 0.02, 0.03, 0.01, 0.04);
    let (mut flow, mut scratch) = coupled(reservoir_field(31), 0.02, 0.03, 0.01, 0.04);
    let expected = native
        .step(Time::new(0.002), MolecularStepControl::default(), || false)
        .unwrap();
    let zero = vec![[0.0; 2]; flow.flow().cell_count()];
    let report = flow
        .step_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            BodyForce2d::Acceleration(&zero),
            ForceControl2d { maximum_mach: 1e-6 },
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

    let original = reservoir_field(31);
    let molecular = prepare(&original, law(original.gas(), 0.02, 0.03, 0.01, 0.04));
    let bytes = Coupled::required_workspace_bytes(&original, &molecular).unwrap();
    assert!(matches!(
        Coupled::new(original, molecular, bytes - 1, || false),
        Err(DissipationError::Invalid(_))
    ));
    let original = reservoir_field(31);
    let molecular = prepare(&original, law(original.gas(), 0.02, 0.03, 0.01, 0.04));
    let mut tight = Coupled::new(original, molecular, bytes, || false).unwrap();
    let old = tight.flow().states().to_vec();
    assert!(matches!(
        tight.step_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            BodyForce2d::Acceleration(&zero),
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
fn g1_second_molecular_cooling_cannot_escape_final_force_mach_admission() {
    let make = || {
        let field = field(
            equal_calorics(),
            CartesianGrid::new(1, 1, 1.0, 1.0, 0.6).unwrap(),
            Boundaries2d {
                x: OPEN,
                y: AxisBoundaries::Periodic,
            },
            |_| ([0.5, 0.5, 0.0], [0.9, 0.0], 2.0),
        );
        let inputs = [
            ReservoirBoundary {
                side: BoundarySide::XMin,
                partial_densities: &[1.0, 1.0, 0.0],
                velocity: [0.9, 0.0],
                pressure_pa: 2.0,
            },
            ReservoirBoundary {
                side: BoundarySide::XMax,
                partial_densities: &[1.0, 1.0, 0.0],
                velocity: [0.9, 0.0],
                pressure_pa: 2.0,
            },
        ];
        bind(field, &inputs)
    };
    let dt = 0.001;
    let cap = ForceControl2d { maximum_mach: 0.55 };
    let force = BodyForce2d::UniformAcceleration([0.01, 0.0]);
    // Establish the failure independently with the original component owners:
    // central transport is admitted, but the second cold-bath half raises Mach.
    let mut manual = make();
    let mut operator = prepare(&manual, law(manual.gas(), 0.02, 0.0, 0.0, 100.0));
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
    assert!(p.velocity[0].hypot(p.velocity[1]) > cap.maximum_mach * p.sound_speed);
    let (mut flow, mut scratch) = coupled(make(), 0.02, 0.0, 0.0, 100.0);
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
