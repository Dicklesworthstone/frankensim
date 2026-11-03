use super::*;
use crate::compressible::spatial::{BodyForce3d, ForceControl3d};

#[cfg(feature = "experimental-sgs")]
#[test]
fn g1_forced_channel_sgs_uses_fresh_halves_and_closes_source_wall_energy() {
    use crate::compressible::spatial::molecular::sgs::{
        ChannelDynamicScalarConfig, ChannelDynamicSgsConfig, ChannelDynamicThermalSubgridConfig,
        SgsUse,
    };
    fn noise(i: usize, salt: u64) -> f64 {
        let mut h = ((i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ salt;
        h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        h ^= h >> 31;
        (h >> 11) as f64 / 9_007_199_254_740_992.0 * 2.0 - 1.0
    }
    let settings = ChannelDynamicThermalSubgridConfig {
        momentum: ChannelDynamicSgsConfig {
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: ChannelDynamicScalarConfig {
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    };
    let walls = [3.0, 2.0].map(|temperature| MolecularWall3d::NoSlip {
        velocity: [0.01, 0.02, 0.0],
        thermal: WallThermalCondition::Isothermal(Temperature::new(temperature)),
    });
    let make_field = || {
        field(
            gas(0.0),
            CartesianGrid::new(5, 4, 5, 0.2, 0.2, 0.2).unwrap(),
            Boundaries3d {
                z: WALL,
                ..PERIODIC
            },
            |point| {
                let c = point.map(|x| (x / 0.2) as usize);
                let cell = c[0] + 5 * (c[1] + 4 * c[2]);
                let eta = (c[2] as f64 + 0.5) / 5.0;
                let envelope = 4.0 * eta * (1.0 - eta);
                let u: [f64; 3] = core::array::from_fn(|axis| {
                    [0.01, 0.02, 0.0][axis] + 0.02 * envelope * noise(cell, 31 * axis as u64 + 7)
                });
                [
                    1.0,
                    u[0],
                    u[1],
                    u[2],
                    3.0 - eta + 0.05 * envelope * noise(cell, 101),
                ]
            },
        )
    };
    let make_operator = |flow: &Euler3d| {
        MolecularTransport3d::with_channel_dynamic_thermal_sgs(
            flow,
            law(flow.gas(), 0.002, 0.003, 0.004),
            settings,
            walls,
            256 * 100,
            || false,
        )
        .unwrap()
    };
    let initial = make_field();
    let source = law(initial.gas(), 0.002, 0.003, 0.004);
    let mut coupled = MolecularEuler3d::with_channel_dynamic_thermal_sgs(
        initial,
        source,
        settings,
        walls,
        256 * 100,
        || false,
    )
    .unwrap();
    let force = BodyForce3d::UniformForceDensity([0.2, 0.0, 0.0]);
    let cap = ForceControl3d { maximum_mach: 1.0 };
    let mut work = coupled.flow().forcing_workspace(40 * 100).unwrap();
    let report = coupled
        .step_with_forcing(
            Time::new(0.001),
            MolecularStepControl::default(),
            force,
            cap,
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(report.coupling_retries, 0);
    assert_eq!(report.source_change, report.advection.source_change);
    assert_eq!(report.first.before, report.before);
    assert_eq!(report.first.after, report.advection.before);
    assert_eq!(report.second.before, report.advection.after);
    assert_eq!(report.second.after, report.after);
    assert_eq!(report.first.elapsed, report.second.elapsed);
    assert_eq!(report.first.elapsed.value() * 2.0, report.elapsed.value());

    // Reconstruct each molecular owner from the actual incoming stage. A
    // cached pre-kick closure in the coupled owner cannot match this trajectory.
    let mut reference = make_field();
    let first = make_operator(&reference)
        .advance(
            &mut reference,
            report.first.elapsed,
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    let volume = reference.grid().cell_volume();
    let half_impulse = 0.5 * 0.2 * report.elapsed.value();
    let first_work: f64 = reference
        .states()
        .iter()
        .map(|q| half_impulse * (q.momentum[0] + 0.5 * half_impulse) / q.density * volume)
        .sum();
    let mut reference_work = reference.forcing_workspace(40 * 100).unwrap();
    let advection = reference
        .step_with_forcing(
            report.elapsed,
            StepControl::default(),
            force,
            cap,
            &mut reference_work,
            || false,
        )
        .unwrap();
    let second_work: f64 = reference
        .states()
        .iter()
        .map(|q| half_impulse * (q.momentum[0] - 0.5 * half_impulse) / q.density * volume)
        .sum();
    let second = make_operator(&reference)
        .advance(
            &mut reference,
            report.second.elapsed,
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.first, first);
    assert_eq!(report.advection, advection);
    assert_eq!(report.second, second);
    for (a, b) in coupled.flow().states().iter().zip(reference.states()) {
        assert_eq!(a.array().map(f64::to_bits), b.array().map(f64::to_bits));
    }
    close(
        report.source_change.impulse[0],
        0.2 * 100.0 * volume * report.elapsed.value(),
        3e-13,
    );
    assert_eq!(report.source_change.impulse[1..], [0.0; 2]);
    close(report.source_change.work, first_work + second_work, 3e-13);
    assert!(report.source_change.work > 0.0);
    let mut outward = [0.0; 5];
    let mut wall_energy_exchange = 0.0;
    for side in BoundarySide::ALL {
        for impulse in coupled.boundary_impulses(side) {
            for (sum, component) in outward.iter_mut().zip(impulse.array()) {
                *sum += component;
            }
            wall_energy_exchange += impulse.energy.abs();
        }
    }
    assert!(wall_energy_exchange > 0.0);
    for (outward, inward) in outward.into_iter().zip(report.boundary_change.array()) {
        close(outward + inward, 0.0, 3e-13);
    }
    close(report.after.mass, report.before.mass, 3e-13);
    for axis in 0..3 {
        close(
            report.after.momentum[axis] - report.before.momentum[axis],
            report.boundary_change.momentum[axis] + report.source_change.impulse[axis],
            3e-13,
        );
    }
    close(
        report.after.energy - report.before.energy,
        report.boundary_change.energy + report.source_change.work,
        3e-12,
    );
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 3e-12);
    }

    // The selected seeded channel has active dynamic stress/heat, not merely
    // a configured closure whose coefficients happen to vanish.
    let base = make_field();
    let source = law(base.gas(), 0.002, 0.003, 0.004);
    let mut base = MolecularEuler3d::new(base, source, || false).unwrap();
    for (side, wall) in [BoundarySide::ZMin, BoundarySide::ZMax]
        .into_iter()
        .zip(walls)
    {
        base = base.with_wall(side, wall).unwrap();
    }
    let mut base_work = base.flow().forcing_workspace(40 * 100).unwrap();
    base.step_with_forcing(
        report.elapsed,
        MolecularStepControl::default(),
        force,
        cap,
        &mut base_work,
        || false,
    )
    .unwrap();
    assert_ne!(base.flow().states(), coupled.flow().states());
}

#[test]
fn g4_forced_molecular_whole_interval_retry_and_late_cancel_replay_atomically() {
    let force = BodyForce3d::UniformAcceleration([20.0, 0.0, 0.0]);
    let cap = ForceControl3d { maximum_mach: 10.0 };
    let make = || {
        let flow = cube(3, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
        let source = law(flow.gas(), 0.0, 0.0, 0.0);
        let mut flow = MolecularEuler3d::new(flow, source, || false).unwrap();
        let mut work = flow.flow().forcing_workspace(40 * 27).unwrap();
        flow.step_with_forcing(
            Time::new(0.001),
            MolecularStepControl::default(),
            force,
            cap,
            &mut work,
            || false,
        )
        .unwrap();
        (flow, work)
    };
    let (mut probe, mut work) = make();
    let before = probe.flow().totals();
    let mut polls = 0;
    let report = probe
        .step_with_forcing(
            Time::new(1.0),
            MolecularStepControl::default(),
            force,
            cap,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(report.coupling_retries > 0);
    assert!(report.elapsed.value() < 1.0);
    assert_eq!(report.first.elapsed.value() * 2.0, report.elapsed.value());
    assert_eq!(report.second.elapsed, report.first.elapsed);
    assert_eq!(report.advection.elapsed, report.elapsed);
    let speed0 = before.momentum[0] / before.mass;
    let speed1 = speed0 + 20.0 * report.elapsed.value();
    close(
        report.source_change.impulse[0],
        before.mass * (speed1 - speed0),
        3e-13,
    );
    close(
        report.source_change.work,
        0.5 * before.mass * (speed1 * speed1 - speed0 * speed0),
        3e-13,
    );
    for q in probe.flow().states() {
        close(
            q.primitive(probe.flow().gas(), || false)
                .unwrap()
                .temperature,
            2.0,
            3e-13,
        );
    }
    let (mut direct, mut direct_work) = make();
    let short = direct
        .step_with_forcing(
            report.elapsed,
            MolecularStepControl::default(),
            force,
            cap,
            &mut direct_work,
            || false,
        )
        .unwrap();
    assert_eq!(short.coupling_retries, 0);
    assert_eq!(short.first, report.first);
    assert_eq!(short.advection, report.advection);
    assert_eq!(short.second, report.second);
    assert_eq!(direct.flow().states(), probe.flow().states());

    for cutoff in [polls / 2, polls] {
        let (mut flow, mut work) = make();
        let old = flow.flow().states().to_vec();
        let time = flow.flow().time();
        let previous = flow.last_step().copied();
        let walls = BoundarySide::ALL.map(|s| flow.boundary_impulses(s).to_vec());
        let mut calls = 0;
        assert_eq!(
            flow.step_with_forcing(
                Time::new(1.0),
                MolecularStepControl::default(),
                force,
                cap,
                &mut work,
                || {
                    calls += 1;
                    calls == cutoff
                },
            ),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(flow.flow().states(), old);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step().copied(), previous);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(flow.boundary_impulses(side), walls[i]);
        }
        let retried = flow
            .step_with_forcing(
                Time::new(1.0),
                MolecularStepControl::default(),
                force,
                cap,
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(retried, report);
        for (actual, expected) in flow.flow().states().iter().zip(probe.flow().states()) {
            assert_eq!(
                actual.array().map(f64::to_bits),
                expected.array().map(f64::to_bits)
            );
        }
    }
}
