use super::*;
use crate::compressible::mixture::molecular::{DissipationControl, MolecularTransport2d};
use crate::compressible::mixture::tests::molecular::{
    impulses as molecular_impulses, law, prepare,
};
use crate::compressible::mixture::{BodyForce2d, ForceControl2d, ForcingWorkspace2d};
use fs_gas::multicomponent::dissipation::ConstantTransport;

const FORCE_CONTROL: ForceControl2d = ForceControl2d { maximum_mach: 8.0 };

fn workspace(flow: &ReactingEuler2d) -> ForcingWorkspace2d {
    let bytes = flow.flow().required_forcing_workspace_bytes().unwrap();
    let work = flow.flow().forcing_workspace(bytes).unwrap();
    assert_eq!(work.workspace_bytes(), bytes);
    work
}

fn total_bytes(
    flow: &ReactingEuler2d,
    operator: &MolecularTransport2d<ConstantTransport>,
    work: &ForcingWorkspace2d,
) -> usize {
    flow.workspace_bytes() + operator.workspace_bytes() + work.workspace_bytes()
}

#[test]
fn g1_reacting_force_retains_original_chemistry_and_signed_kinetic_work() {
    let mechanism = make_mechanism(6.0, 8.0, 4.5);
    let initial_velocity = [0.7, -0.3];
    let make = || {
        reactive(
            raw_field(
                &mechanism,
                CartesianGrid::new(4, 3, 0.25, 1.0 / 3.0, 0.7).unwrap(),
                PERIODIC,
                |_| ([1.2, 0.5, 0.0], initial_velocity, 1.7),
            ),
            &mechanism,
        )
    };
    for mode in 0..4 {
        let mut flow = make();
        let mut work = workspace(&flow);
        let bytes = flow.workspace_bytes() + work.workspace_bytes();
        let values = vec![[-3.0, 1.0]; flow.flow().cell_count()];
        let force = match mode {
            0 => BodyForce2d::UniformAcceleration([-3.0, 1.0]),
            1 => BodyForce2d::Acceleration(&values),
            2 => BodyForce2d::UniformForceDensity([-3.0, 1.0]),
            _ => BodyForce2d::ForceDensity(&values),
        };
        let report = flow
            .step_with_forcing(
                Time::new(0.005),
                ReactiveStepControl::default(),
                force,
                FORCE_CONTROL,
                &mut work,
                || false,
            )
            .unwrap();
        coupled_balance(&report);
        assert_eq!(report.source_change, report.transport.source_change);
        let scale = if mode < 2 { 1.0 } else { 1.0 / 1.7 };
        let dv = [
            -3.0 * report.elapsed.value() * scale,
            report.elapsed.value() * scale,
        ];
        let mass = report.before.mass();
        let expected_work: f64 = (0..2)
            .map(|axis| mass * dv[axis] * (initial_velocity[axis] + 0.5 * dv[axis]))
            .sum();
        assert!(expected_work < 0.0, "braking must remove mechanical energy");
        close(report.source_change.work, expected_work, 3e-12);
        for axis in 0..2 {
            close(report.source_change.impulse[axis], mass * dv[axis], 3e-12);
        }
        close(
            report.after.energy - report.before.energy,
            expected_work + report.chemistry.referenced_energy_change_j,
            3e-10,
        );
        let mut reference = make();
        reference
            .step(report.elapsed, ReactiveStepControl::default(), || false)
            .unwrap();
        for index in 0..flow.flow().cell_count() {
            let q = flow.flow().cell(index).unwrap();
            let p = q.primitive(flow.flow().gas(), || false).unwrap();
            let reference_q = reference.flow().cell(index).unwrap();
            let reference_p = reference_q
                .primitive(reference.flow().gas(), || false)
                .unwrap();
            for species in 0..3 {
                close(
                    q.partial_densities[species],
                    reference_q.partial_densities[species],
                    3e-10,
                );
            }
            assert_eq!(q.partial_densities[2], 0.0);
            close(p.temperature, reference_p.temperature, 3e-10);
            close(p.pressure, reference_p.pressure, 3e-10);
            for axis in 0..2 {
                close(p.velocity[axis], initial_velocity[axis] + dv[axis], 3e-12);
            }
        }
        assert_eq!(flow.workspace_bytes() + work.workspace_bytes(), bytes);
    }
}

#[test]
fn g4_forced_heat_release_retries_one_interval_and_zero_force_preserves_the_path() {
    let mechanism = make_mechanism(40.0, 500.0, 3.5);
    let mut flow = uniform(&mechanism);
    let mut work = workspace(&flow);
    let acceleration = [6.0, -2.0];
    let report = flow
        .step_with_forcing(
            Time::new(0.2),
            ReactiveStepControl::default(),
            BodyForce2d::UniformAcceleration(acceleration),
            FORCE_CONTROL,
            &mut work,
            || false,
        )
        .unwrap();
    coupled_balance(&report);
    assert!(report.coupling_retries > 0);
    assert_eq!(report.chemistry.cell_solves, 2 * flow.flow().cell_count());
    assert_eq!(flow.flow().time(), report.elapsed);
    for axis in 0..2 {
        close(
            report.source_change.impulse[axis],
            report.before.mass() * acceleration[axis] * report.elapsed.value(),
            3e-11,
        );
    }
    let spacing = flow.flow().grid().spacing();
    for index in 0..flow.flow().cell_count() {
        let p = flow
            .flow()
            .cell(index)
            .unwrap()
            .primitive(flow.flow().gas(), || false)
            .unwrap();
        let rate: f64 = (0..2)
            .map(|axis| (p.velocity[axis].abs() + p.sound_speed) / spacing[axis])
            .sum();
        assert!(rate * report.elapsed.value() <= 0.4 * (1.0 + 1e-13));
    }

    let mut plain = uniform(&mechanism);
    let expected = plain
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    let mut zero = uniform(&mechanism);
    let mut work = workspace(&zero);
    let values = vec![[0.0; 2]; zero.flow().cell_count()];
    let actual = zero
        .step_with_forcing(
            Time::new(0.2),
            ReactiveStepControl::default(),
            BodyForce2d::ForceDensity(&values),
            ForceControl2d {
                maximum_mach: 1e-12,
            },
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(zero.flow().states(), plain.flow().states());
    assert_eq!(zero.flow().time(), plain.flow().time());
    assert_eq!(impulses(&zero), impulses(&plain));
}

fn accelerated_reacting_contact(n: usize) -> f64 {
    let mechanism = make_mechanism(0.0, 2.0, 3.5);
    let dx = 1.0 / n as f64;
    let z = std::f64::consts::PI * dx;
    let average = (z.sin() / z).powi(2);
    let velocity = [0.4, 0.25];
    let acceleration = [12.0, -5.0];
    let exact = |point: [f64; 2], time: f64| {
        let displacement = (velocity[0] + velocity[1]) * time
            + 0.5 * (acceleration[0] + acceleration[1]) * time * time;
        0.5 + 0.2
            * average
            * (-4.0 * time).exp()
            * (2.0 * std::f64::consts::PI * (point[0] + point[1] - displacement)).sin()
    };
    let grid = CartesianGrid::new(n, n, dx, dx, 0.6).unwrap();
    let mut flow = reactive(
        raw_field(&mechanism, grid, PERIODIC, |point| {
            let a = exact(point, 0.0);
            ([a, 1.0 - a, 0.0], velocity, 1.0)
        }),
        &mechanism,
    );
    let mut work = workspace(&flow);
    let end = 0.05;
    let mut steps = 0;
    while flow.flow().time().value() < end {
        let report = flow
            .step_with_forcing(
                Time::new(end - flow.flow().time().value()),
                ReactiveStepControl::default(),
                BodyForce2d::UniformAcceleration(acceleration),
                FORCE_CONTROL,
                &mut work,
                || false,
            )
            .unwrap();
        coupled_balance(&report);
        steps += 1;
        assert!(steps < 100);
    }
    let mut error = 0.0;
    for index in 0..grid.len() {
        let q = flow.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        close(p.temperature, 1.0, 3e-10);
        close(p.pressure, 1.0, 3e-10);
        for axis in 0..2 {
            close(
                p.velocity[axis],
                velocity[axis] + acceleration[axis] * end,
                3e-11,
            );
        }
        assert_eq!(q.partial_densities[2], 0.0);
        error += (q.partial_densities[0] - exact(grid.cell_center(index).unwrap(), end)).abs();
    }
    error / grid.len() as f64
}

#[test]
fn g1_accelerating_oblique_reacting_contact_matches_the_physical_phase_and_refines() {
    let errors = [6, 12, 24].map(accelerated_reacting_contact);
    eprintln!("accelerating reacting contact L1, 6/12/24 grids: {errors:?}");
    assert!(errors[0] < 8e-3);
    assert!(errors[1] < 0.7 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.65 * errors[1], "{errors:?}");
}

#[test]
fn g4_forced_molecular_wall_impulses_and_all_histories_commit_atomically() {
    let mechanism = make_mechanism(5.0, 8.0, 3.5);
    let create = || {
        let flow = wall_field(&mechanism);
        let work = workspace(&flow);
        let operator = prepare(flow.flow(), law(flow.flow().gas(), 0.02, 0.03, 0.01, 0.02));
        (flow, operator, work)
    };
    let (mut reference, mut reference_operator, mut reference_work) = create();
    let bytes = total_bytes(&reference, &reference_operator, &reference_work);
    let molecular = DissipationControl {
        maximum_step_s: 0.0005,
        ..DissipationControl::default()
    };
    let values: Vec<_> = (0..reference.flow().cell_count())
        .map(|index| [0.3 + 0.07 * index as f64, -0.5 + 0.02 * index as f64])
        .collect();
    let force = BodyForce2d::ForceDensity(&values);
    reference
        .step_with_dissipation_and_forcing(
            &mut reference_operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            molecular,
            bytes,
            force,
            FORCE_CONTROL,
            &mut reference_work,
            || false,
        )
        .unwrap();
    let mut polls = 0;
    let expected = reference
        .step_with_dissipation_and_forcing(
            &mut reference_operator,
            Time::new(0.002),
            ReactiveStepControl::default(),
            molecular,
            bytes,
            force,
            FORCE_CONTROL,
            &mut reference_work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    coupled_balance(&expected);
    let split = expected.dissipation.unwrap();
    let mut outward = [0.0; 2];
    for (species, momentum, energy) in impulses(&reference) {
        assert_eq!(species, [0.0; 3]);
        assert_eq!(energy, 0.0);
        for axis in 0..2 {
            outward[axis] += momentum[axis];
        }
    }
    for axis in 0..2 {
        close(
            expected.source_change.impulse[axis],
            values.iter().map(|v| v[axis]).sum::<f64>()
                * reference.flow().grid().cell_volume()
                * expected.elapsed.value(),
            3e-12,
        );
        close(
            expected.after.momentum[axis] - expected.before.momentum[axis] + outward[axis]
                - expected.source_change.impulse[axis],
            0.0,
            3e-11,
        );
        close(
            expected.boundary_change.momentum[axis],
            expected.transport.boundary_change.momentum[axis]
                + split.first.boundary_change.momentum[axis]
                + split.second.boundary_change.momentum[axis],
            3e-12,
        );
    }
    assert_ne!(
        expected.boundary_change.momentum,
        expected.transport.boundary_change.momentum
    );

    let (mut flow, mut operator, mut work) = create();
    flow.step_with_dissipation_and_forcing(
        &mut operator,
        Time::new(0.001),
        ReactiveStepControl::default(),
        molecular,
        bytes,
        force,
        FORCE_CONTROL,
        &mut work,
        || false,
    )
    .unwrap();
    let cells = flow.flow().states().to_vec();
    let clock = flow.flow().time();
    let report = *flow.last_step().unwrap();
    let advection_report = *flow.flow().last_step().unwrap();
    let molecular_report = *operator.last_step().unwrap();
    let loads = impulses(&flow);
    let molecular_loads = molecular_impulses(&operator);
    let mut attempted = 0;
    let error = flow
        .step_with_dissipation_and_forcing(
            &mut operator,
            Time::new(0.002),
            ReactiveStepControl::default(),
            molecular,
            bytes,
            force,
            FORCE_CONTROL,
            &mut work,
            || {
                attempted += 1;
                attempted >= polls
            },
        )
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(flow.flow().states(), cells);
    assert_eq!(flow.flow().time(), clock);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.flow().last_step(), Some(&advection_report));
    assert_eq!(operator.last_step(), Some(&molecular_report));
    assert_eq!(impulses(&flow), loads);
    assert_eq!(molecular_impulses(&operator), molecular_loads);
    assert_eq!(total_bytes(&flow, &operator, &work), bytes);

    assert!(matches!(
        flow.step_with_dissipation_and_forcing(
            &mut operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            molecular,
            bytes - 1,
            force,
            FORCE_CONTROL,
            &mut work,
            || false,
        ),
        Err(ReactiveError::Invalid(_))
    ));
    assert_eq!(flow.flow().states(), cells);
    assert_eq!(operator.last_step(), Some(&molecular_report));
    let replay = flow
        .step_with_dissipation_and_forcing(
            &mut operator,
            Time::new(0.002),
            ReactiveStepControl::default(),
            molecular,
            bytes,
            force,
            FORCE_CONTROL,
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(flow.flow().states(), reference.flow().states());
    assert_eq!(impulses(&flow), impulses(&reference));
    assert_eq!(
        molecular_impulses(&operator),
        molecular_impulses(&reference_operator)
    );
}
