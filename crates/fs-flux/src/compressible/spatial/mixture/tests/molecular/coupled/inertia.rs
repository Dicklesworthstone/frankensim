use super::*;
use crate::compressible::spatial::mixture::molecular::{
    InertialWall3d, WallAcceleration3d, WallInertiaControl3d, WallInertiaExchange3d,
};
use fs_qty::Mass;

const MU: f64 = 0.75;
const INITIAL: [f64; 3] = [1.0, 0.0, -0.4];

fn native(
    mass: f64,
    thermal: bool,
    controls: WallInertiaControl3d,
) -> (MixtureEuler3d, MolecularTransport3d<ConstantTransport>) {
    let flow = field(
        equal_calorics(),
        CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
        Boundaries3d {
            y: WALLS,
            ..PERIODIC
        },
        |_| ([0.7, 0.3, 0.0], [0.0; 3], 1.0),
    );
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.001, MU, 0.0, if thermal { 0.1 } else { 0.0 }),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(&flow, &source).unwrap()
        + WallInertiaControl3d::solver_workspace_bytes();
    let operator = MolecularTransport3d::new(&flow, source, bytes, || false)
        .unwrap()
        .with_wall(
            BoundarySide::YMax,
            MolecularWall3d::NoSlip {
                velocity: INITIAL,
                thermal: if thermal {
                    WallThermalCondition::Isothermal(Temperature::new(0.8))
                } else {
                    WallThermalCondition::Adiabatic
                },
            },
        )
        .unwrap()
        .with_wall_inertia(
            &[InertialWall3d {
                side: BoundarySide::YMax,
                mass: Mass::new(mass),
            }],
            controls,
        )
        .unwrap();
    assert_eq!(operator.workspace_bytes() + flow.workspace_bytes(), bytes);
    (flow, operator)
}
fn case(mass: f64, thermal: bool, controls: WallInertiaControl3d) -> (Coupled, ForcingWorkspace3d) {
    let (flow, operator) = native(mass, thermal, controls);
    coupled(flow, operator)
}
fn speed(operator: &MolecularTransport3d<ConstantTransport>) -> [f64; 3] {
    match operator.wall(BoundarySide::YMax) {
        MolecularWall3d::NoSlip { velocity, .. } => velocity,
        _ => panic!("finite no-slip wall"),
    }
}
fn exchange(report: &MolecularReport) -> WallInertiaExchange3d {
    *report
        .wall_inertia
        .as_ref()
        .unwrap()
        .wall(BoundarySide::YMax)
        .unwrap()
}
fn energy(mass: f64, u: [f64; 3]) -> f64 {
    0.5 * mass * u.iter().map(|v| v * v).sum::<f64>()
}
fn finite_balance(report: &MolecularStepReport3d) {
    super::check(report);
    let a = exchange(&report.first);
    let b = exchange(&report.second);
    assert_eq!(a.velocity_after, b.velocity_before);
    for axis in [0, 2] {
        close(
            report.after.momentum[axis] - report.before.momentum[axis]
                + a.mass.value() * (b.velocity_after[axis] - a.velocity_before[axis])
                - report.source_change.impulse[axis],
            0.0,
            4e-11,
        );
    }
    close(
        report.after.energy - report.before.energy
            + a.kinetic_energy_change_j
            + b.kinetic_energy_change_j
            + a.thermal_exchange_j
            + b.thermal_exchange_j
            - report.source_change.work,
        0.0,
        4e-11,
    );
    for q in [a, b] {
        close(q.work_residual_j, 0.0, 1e-10);
        for residual in q.momentum_residual {
            close(residual, 0.0, 1e-10);
        }
    }
    assert!(
        report.second.wall_inertia.unwrap().field_evaluations
            > report.first.wall_inertia.unwrap().field_evaluations
    );
}

#[test]
fn g1_finite_wall_recoil_refines_to_independent_heavy_and_light_mass_exponentials() {
    // One unit-volume cell, rho=A=H=1, top no-slip and bottom slip:
    // K=2*mu*A/H; relative speed decays at K*(1/mf+1/mw).
    // The light carrier demands a smaller physical step than the fluid CFL.
    for (mass, end, steps) in [
        (0.5, 0.16, [0.004, 0.002, 0.001]),
        (0.03, 0.06, [0.002, 0.001, 0.0005]),
    ] {
        let decay = (-2.0 * MU * (1.0 + 1.0 / mass) * end).exp();
        let expected_fluid = INITIAL.map(|v| mass / (1.0 + mass) * v * (1.0 - decay));
        let expected_wall = INITIAL.map(|v| v * (mass + decay) / (1.0 + mass));
        let released =
            energy(mass, INITIAL) - energy(1.0, expected_fluid) - energy(mass, expected_wall);
        let mut errors = Vec::new();
        let mut counts = Vec::new();
        for step in steps {
            let (mut flow, mut operator) = native(mass, false, WallInertiaControl3d::default());
            let initial_energy = flow.totals().energy;
            let report = operator
                .advance(
                    &mut flow,
                    Time::new(end),
                    DissipationControl {
                        maximum_step_s: step,
                        ..DissipationControl::default()
                    },
                    || false,
                )
                .unwrap();
            super::super::check(report);
            let wall = speed(&operator);
            let fluid = flow
                .cell(0)
                .unwrap()
                .primitive(flow.gas(), || false)
                .unwrap();
            let heat = flow.totals().energy - initial_energy - energy(1.0, fluid.velocity);
            let ledger = exchange(&report);
            close(
                flow.totals().energy + energy(mass, wall),
                initial_energy + energy(mass, INITIAL),
                4e-11,
            );
            for axis in [0, 2] {
                close(
                    flow.totals().momentum[axis] + mass * wall[axis],
                    mass * INITIAL[axis],
                    4e-11,
                );
            }
            close(
                heat,
                energy(mass, INITIAL) - energy(mass, wall) - energy(1.0, fluid.velocity),
                4e-11,
            );
            close(ledger.thermal_exchange_j, 0.0, 4e-11);
            assert!(heat > 0.0 && wall[0] < INITIAL[0] && fluid.velocity[0] > 0.0);
            assert_eq!(flow.cell(0).unwrap().partial_densities, [0.7, 0.3, 0.0]);
            assert_eq!(flow.time().value(), 0.0);
            errors.push([
                (fluid.velocity[0] - expected_fluid[0]).abs(),
                (wall[0] - expected_wall[0]).abs(),
                (heat - released).abs(),
            ]);
            counts.push(report.accepted_steps);
        }
        eprintln!("finite wall mass {mass}, SSP counts {counts:?}, errors {errors:?}");
        for (count, step) in counts.iter().zip(steps) {
            assert!((*count as f64 - (end / step).round()).abs() <= 1.0);
        }
        assert!(counts[1] >= 2 * counts[0] - 2 && counts[2] >= 2 * counts[1] - 2);
        for component in 0..3 {
            for pair in errors.windows(2) {
                let ratio = pair[0][component] / pair[1][component];
                assert!(
                    ratio > 3.3 && ratio < 4.8,
                    "mass {mass} component {component}: {ratio}"
                );
            }
            assert!(errors[2][component] < 4e-5);
        }
    }
}

#[test]
fn g2_finite_wall_body_source_and_thermal_reservoir_close_independent_storage() {
    let (mut flow, mut scratch) = case(0.4, true, WallInertiaControl3d::default());
    let bytes = flow.workspace_bytes();
    let report = flow
        .step_with_forcing(
            Time::new(0.02),
            MolecularStepControl::default(),
            BodyForce3d::UniformAcceleration([0.2, 0.0, -0.1]),
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
    finite_balance(&report);
    let first = exchange(&report.first);
    let second = exchange(&report.second);
    assert!(first.thermal_exchange_j + second.thermal_exchange_j > 1e-5);
    assert!(first.kinetic_energy_change_j + second.kinetic_energy_change_j < 0.0);
    assert!(
        flow.boundary_impulses(BoundarySide::YMax)
            .map(|q| q.momentum[1])
            .sum::<f64>()
            > 0.0
    );
    close(first.guide_impulse + second.guide_impulse, 0.0, 1e-12);
    assert_eq!(speed(flow.molecular()), second.velocity_after);
    assert_eq!(flow.molecular().last_step(), Some(&report.second));
    assert_eq!(flow.workspace_bytes(), bytes);
    close(
        report.source_change.impulse[0],
        0.2 * report.elapsed.value(),
        2e-13,
    );
    close(
        report.source_change.impulse[2],
        -0.1 * report.elapsed.value(),
        2e-13,
    );
}

#[test]
fn g3_finite_wall_shortening_replays_original_carrier_and_late_cancel_preserves_metadata() {
    let (mut flow, mut scratch) = case(0.4, false, WallInertiaControl3d::default());
    let force = BodyForce3d::UniformAcceleration([400.0, 0.0, -200.0]);
    let cap = ForceControl3d { maximum_mach: 0.8 };
    let report = flow
        .step_with_forcing(
            Time::new(0.2),
            MolecularStepControl::default(),
            force,
            cap,
            &mut scratch,
            || false,
        )
        .unwrap();
    finite_balance(&report);
    assert!(report.coupling_retries > 0);
    let (mut reference, mut reference_scratch) = case(0.4, false, WallInertiaControl3d::default());
    let replay = reference
        .step_with_forcing(
            report.elapsed,
            MolecularStepControl::default(),
            force,
            cap,
            &mut reference_scratch,
            || false,
        )
        .unwrap();
    assert_eq!(replay.coupling_retries, 0);
    assert_eq!(flow.flow().states(), reference.flow().states());
    assert_eq!(speed(flow.molecular()), speed(reference.molecular()));
    assert_eq!(loads(&flow), loads(&reference));
    assert_eq!(molecular_loads(&flow), molecular_loads(&reference));
    for (actual, expected) in [(report.first, replay.first), (report.second, replay.second)] {
        assert_eq!(exchange(&actual), exchange(&expected));
        assert!(
            actual.wall_inertia.unwrap().field_evaluations
                > expected.wall_inertia.unwrap().field_evaluations
        );
    }
    let before = (
        flow.flow().states().to_vec(),
        flow.flow().time(),
        flow.last_step().copied(),
        flow.flow().last_step().copied(),
        flow.molecular().last_step().copied(),
        speed(flow.molecular()),
        loads(&flow),
        molecular_loads(&flow),
        flow.workspace_bytes(),
    );
    let mut polls = 0;
    let expected = reference
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        flow.step(Time::new(0.001), MolecularStepControl::default(), || {
            count += 1;
            count >= polls
        }),
        Err(DissipationError::Cancelled)
    ));
    let after = (
        flow.flow().states().to_vec(),
        flow.flow().time(),
        flow.last_step().copied(),
        flow.flow().last_step().copied(),
        flow.molecular().last_step().copied(),
        speed(flow.molecular()),
        loads(&flow),
        molecular_loads(&flow),
        flow.workspace_bytes(),
    );
    assert_eq!(before, after);
    let recovered = flow
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(recovered, expected);
    assert_eq!(flow.flow().states(), reference.flow().states());
}

#[test]
fn g4_finite_wall_cumulative_probe_budget_and_conflicting_drive_are_atomic() {
    let (mut baseline, _) = case(0.4, false, WallInertiaControl3d::default());
    let report = baseline
        .step(Time::new(0.01), MolecularStepControl::default(), || false)
        .unwrap();
    let first_work = report.first.wall_inertia.unwrap().field_evaluations;
    let (mut limited, _) = case(
        0.4,
        false,
        WallInertiaControl3d {
            maximum_field_evaluations: first_work + 1,
            ..WallInertiaControl3d::default()
        },
    );
    let before = limited.flow().states().to_vec();
    let result = limited.step(
        Time::new(0.01),
        MolecularStepControl {
            maximum_coupling_retries: 2,
            ..MolecularStepControl::default()
        },
        || false,
    );
    assert!(format!("{result:?}").contains("cumulative field-evaluation budget"));
    assert_eq!(limited.flow().states(), before);
    assert_eq!(speed(limited.molecular()), INITIAL);
    assert!(limited.last_step().is_none() && limited.molecular().last_step().is_none());
    assert_eq!(limited.flow().time().value(), 0.0);
    let drive = [WallAcceleration3d {
        side: BoundarySide::YMax,
        acceleration: [1.0, 0.0, 0.0],
    }];
    assert!(
        limited
            .step_with_wall_acceleration(
                Time::new(0.01),
                MolecularStepControl::default(),
                &drive,
                || false
            )
            .is_err()
    );
    assert_eq!(limited.flow().states(), before);
    let (mut zero, _) = case(0.4, false, WallInertiaControl3d::default());
    let accepted = zero
        .step_with_wall_acceleration(
            Time::new(0.01),
            MolecularStepControl::default(),
            &[WallAcceleration3d {
                side: BoundarySide::YMax,
                acceleration: [0.0; 3],
            }],
            || false,
        )
        .unwrap();
    assert_eq!(accepted, report);
    assert_eq!(zero.flow().states(), baseline.flow().states());
}

#[test]
fn g2_two_finite_walls_share_one_native_shear_solve_and_conserve_total_storage() {
    let flow = field(
        equal_calorics(),
        CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
        Boundaries3d {
            y: WALLS,
            ..PERIODIC
        },
        |_| ([0.7, 0.3, 0.0], [0.0; 3], 1.0),
    );
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.001, MU, 0.0, 0.0),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(&flow, &source).unwrap()
        + WallInertiaControl3d::solver_workspace_bytes();
    let models = [
        (BoundarySide::YMin, 0.2, [-0.5, 0.0, 0.2]),
        (BoundarySide::YMax, 0.6, [0.9, 0.0, -0.4]),
    ];
    let mut op = MolecularTransport3d::new(&flow, source, bytes, || false).unwrap();
    for (side, _, velocity) in models {
        op = op
            .with_wall(
                side,
                MolecularWall3d::NoSlip {
                    velocity,
                    thermal: WallThermalCondition::Adiabatic,
                },
            )
            .unwrap();
    }
    let mut op = op
        .with_wall_inertia(
            &models.map(|(side, mass, _)| InertialWall3d {
                side,
                mass: Mass::new(mass),
            }),
            WallInertiaControl3d::default(),
        )
        .unwrap();
    let mut flow = flow;
    let initial = flow.totals();
    let report = op
        .advance(
            &mut flow,
            Time::new(0.12),
            DissipationControl {
                maximum_step_s: 0.004,
                ..DissipationControl::default()
            },
            || false,
        )
        .unwrap();
    super::super::check(report);
    let mut total_energy = flow.totals().energy - initial.energy;
    let mut momentum =
        std::array::from_fn::<_, 3, _>(|a| flow.totals().momentum[a] - initial.momentum[a]);
    for (side, mass, velocity) in models {
        let q = *report.wall_inertia.unwrap().wall(side).unwrap();
        assert_eq!(q.velocity_before, velocity);
        assert_ne!(q.velocity_after, velocity);
        total_energy += energy(mass, q.velocity_after) - energy(mass, velocity);
        for axis in [0, 2] {
            momentum[axis] += mass * (q.velocity_after[axis] - velocity[axis]);
        }
        close(q.thermal_exchange_j, 0.0, 4e-11);
    }
    close(total_energy, 0.0, 4e-11);
    for axis in [0, 2] {
        close(momentum[axis], 0.0, 4e-11);
    }
    assert!(
        flow.cell(0)
            .unwrap()
            .primitive(flow.gas(), || false)
            .unwrap()
            .temperature
            > 1.0
    );
}
