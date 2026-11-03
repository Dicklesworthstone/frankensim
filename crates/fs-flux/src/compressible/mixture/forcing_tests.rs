use super::*;
use crate::compressible::AxisBoundaries;
use crate::compressible::mixture::tests::{close, gas, source};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const PERIODIC: Boundaries2d = Boundaries2d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
};
const CONTROL: ForceControl2d = ForceControl2d { maximum_mach: 8.0 };

fn material(reference: f64, linear: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, reference),
        source("B", 2.0 * R, 3.5, linear, -2.0 * reference),
        source("absent-C", 0.5 * R, 3.5, 0.0, 3.0 * reference),
    ])
}

fn field(
    gas: MixtureGas,
    boundaries: Boundaries2d,
    initial: impl Fn(usize) -> ([f64; 3], [f64; 2], f64),
) -> MixtureEuler2d {
    let grid = CartesianGrid::new(6, 4, 0.25, 0.3, 0.4).unwrap();
    let mut states = vec![0.0; grid.len() * 6];
    for (i, row) in states.chunks_exact_mut(6).enumerate() {
        let (species, velocity, pressure) = initial(i);
        state_into(&gas, &species, velocity, pressure, row).unwrap();
    }
    let bytes = MixtureEuler2d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler2d::new(gas, grid, states, boundaries, bytes, || false).unwrap()
}

fn work(flow: &MixtureEuler2d) -> ForcingWorkspace2d {
    let bytes = flow.required_forcing_workspace_bytes().unwrap();
    let result = flow.forcing_workspace(bytes).unwrap();
    assert_eq!(result.workspace_bytes(), bytes);
    result
}

fn balance(report: &StepReport) {
    assert!(report.maximum_cfl <= StepControl::default().cfl * (1.0 + 1e-13));
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 3e-12);
    }
    for residual in report.conservation_residual.momentum {
        close(residual, 0.0, 3e-12);
    }
    close(report.conservation_residual.energy, 0.0, 3e-12);
    close(report.absolute_energy_residual_j, 0.0, 3e-10);
}

#[test]
fn g1_forced_mixture_uses_composition_mass_signed_work_and_original_caloric_gauge() {
    let velocity = [0.8, -0.3];
    for densities in [[0.2, 0.5, 0.0], [2.4, 1.0, 0.0]] {
        let mut reference: [Option<Vec<f64>>; 2] = [None, None];
        for formation in [0.0, 1e12] {
            for mode in 0..4 {
                let mut flow = field(material(formation, 0.4), PERIODIC, |_| {
                    (densities, velocity, 1.7)
                });
                let mut workspace = work(&flow);
                let initial = flow
                    .cell(0)
                    .unwrap()
                    .primitive(flow.gas(), || false)
                    .unwrap();
                let source_field = vec![[-3.0, 1.0]; flow.cell_count()];
                let force = match mode {
                    0 => BodyForce2d::UniformAcceleration([-3.0, 1.0]),
                    1 => BodyForce2d::Acceleration(&source_field),
                    2 => BodyForce2d::UniformForceDensity([-3.0, 1.0]),
                    _ => BodyForce2d::ForceDensity(&source_field),
                };
                let report = flow
                    .step_with_forcing(
                        Time::new(0.005),
                        StepControl::default(),
                        force,
                        CONTROL,
                        &mut workspace,
                        || false,
                    )
                    .unwrap();
                let scale = if mode < 2 {
                    1.0
                } else {
                    initial.density.recip()
                };
                let dv = [
                    -3.0 * report.elapsed.value() * scale,
                    report.elapsed.value() * scale,
                ];
                let mass = report.before.mass();
                let expected_work: f64 = (0..2)
                    .map(|a| mass * dv[a] * (velocity[a] + 0.5 * dv[a]))
                    .sum();
                assert!(expected_work < 0.0);
                close(report.source_change.work, expected_work, 3e-13);
                close(
                    report.after.energy - report.before.energy,
                    expected_work,
                    3e-13,
                );
                for axis in 0..2 {
                    close(report.source_change.impulse[axis], mass * dv[axis], 3e-13);
                }
                assert_eq!(report.after.species_masses, report.before.species_masses);
                for i in 0..flow.cell_count() {
                    let q = flow.cell(i).unwrap();
                    let p = q.primitive(flow.gas(), || false).unwrap();
                    assert_eq!(q.partial_densities, densities);
                    close(p.temperature, initial.temperature, 3e-11);
                    close(p.pressure, initial.pressure, 3e-11);
                    for axis in 0..2 {
                        close(p.velocity[axis], velocity[axis] + dv[axis], 3e-13);
                    }
                }
                if let Some(expected) = &reference[mode / 2] {
                    assert_eq!(flow.states(), expected);
                } else {
                    reference[mode / 2] = Some(flow.states().to_vec());
                }
                balance(&report);
            }
        }
    }
}

#[test]
fn g1_forced_mixture_spatial_sources_use_each_current_species_density() {
    for force_density in [false, true] {
        let mut flow = field(material(30.0, 0.0), PERIODIC, |i| {
            let row = (i / 6) as f64;
            let rho = 0.6 + 0.3 * row;
            let a = 0.15 + 0.15 * row;
            ([a * rho, (1.0 - a) * rho, 0.0], [0.0; 2], 1.0)
        });
        let before = flow.states().to_vec();
        let forces: Vec<_> = (0..flow.cell_count())
            .map(|i| [0.2 + 0.1 * (i / 6) as f64, 0.0])
            .collect();
        let source = if force_density {
            BodyForce2d::ForceDensity(&forces)
        } else {
            BodyForce2d::Acceleration(&forces)
        };
        let mut workspace = work(&flow);
        let report = flow
            .step_with_forcing(
                Time::new(0.002),
                StepControl::default(),
                source,
                CONTROL,
                &mut workspace,
                || false,
            )
            .unwrap();
        for (i, row) in before.chunks_exact(6).enumerate() {
            let initial = cell(row, 3).primitive(flow.gas(), || false).unwrap();
            let current = flow.cell(i).unwrap();
            let p = current.primitive(flow.gas(), || false).unwrap();
            let scale = if force_density {
                initial.density.recip()
            } else {
                1.0
            };
            close(
                p.velocity[0],
                forces[i][0] * report.elapsed.value() * scale,
                3e-13,
            );
            close(p.velocity[1], 0.0, 3e-13);
            close(p.temperature, initial.temperature, 3e-11);
            close(p.pressure, initial.pressure, 3e-11);
            assert_eq!(current.partial_densities, &row[..3]);
        }
        balance(&report);
    }

    // All active species share constant Cp/R=3.5: this material contact has a
    // common gamma despite different molecular masses and composition. Its
    // pressure equilibrium is compatible with conservative mixture advection.
    let mut moving = field(material(30.0, 0.0), PERIODIC, |i| {
        let species = if i % 6 < 3 {
            [0.15, 0.35, 0.0]
        } else {
            [1.2, 0.8, 0.0]
        };
        (species, [0.2, 0.0], 1.0)
    });
    let original = moving.states().to_vec();
    let forces = vec![[2.0, 0.0]; moving.cell_count()];
    let mut workspace = work(&moving);
    let report = moving
        .step_with_forcing(
            Time::new(0.02),
            StepControl::default(),
            BodyForce2d::Acceleration(&forces),
            CONTROL,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert!(
        moving
            .states()
            .chunks_exact(6)
            .zip(original.chunks_exact(6))
            .any(
                |(q, old)| (sum(q[..3].iter().copied()) - sum(old[..3].iter().copied())).abs()
                    > 1e-6
            )
    );
    for i in 0..moving.cell_count() {
        let q = moving.cell(i).unwrap();
        let p = q.primitive(moving.gas(), || false).unwrap();
        close(p.velocity[0], 0.2 + 2.0 * report.elapsed.value(), 3e-12);
        close(p.velocity[1], 0.0, 3e-12);
        close(p.pressure, 1.0, 3e-12);
        assert_eq!(q.partial_densities[2], 0.0);
    }
    close(
        report.source_change.impulse[0],
        report.before.mass() * 2.0 * report.elapsed.value(),
        3e-12,
    );
    balance(&report);
}

#[test]
fn g1_forced_mixture_retries_admit_rk_and_final_mach_cfl_without_reusing_impulse() {
    let make = || {
        field(material(30.0, 0.4), PERIODIC, |_| {
            ([0.4, 0.8, 0.0], [0.0; 2], 2.3)
        })
    };
    for maximum_mach in [8.0, 0.1] {
        let mut flow = make();
        let mut workspace = work(&flow);
        let force = BodyForce2d::UniformAcceleration([80.0, -40.0]);
        let cap = ForceControl2d { maximum_mach };
        let initial = flow
            .cell(0)
            .unwrap()
            .primitive(flow.gas(), || false)
            .unwrap();
        let report = flow
            .step_with_forcing(
                Time::new(0.1),
                StepControl::default(),
                force,
                cap,
                &mut workspace,
                || false,
            )
            .unwrap();
        assert!(report.halvings > 0);
        for i in 0..flow.cell_count() {
            let p = flow
                .cell(i)
                .unwrap()
                .primitive(flow.gas(), || false)
                .unwrap();
            for (axis, a) in [80.0, -40.0].into_iter().enumerate() {
                close(p.velocity[axis], a * report.elapsed.value(), 3e-13);
                close(
                    report.source_change.impulse[axis],
                    report.before.mass() * a * report.elapsed.value(),
                    3e-13,
                );
            }
            close(p.temperature, initial.temperature, 3e-11);
            assert!(
                p.velocity[0].hypot(p.velocity[1]) <= maximum_mach * p.sound_speed * (1.0 + 1e-13)
            );
        }
        balance(&report);
        let mut replay = make();
        let mut scratch = work(&replay);
        let accepted = replay
            .step_with_forcing(
                report.elapsed,
                StepControl::default(),
                force,
                cap,
                &mut scratch,
                || false,
            )
            .unwrap();
        assert_eq!(accepted.halvings, 0);
        assert_eq!(accepted.source_change, report.source_change);
        assert_eq!(replay.states(), flow.states());
    }
}

fn wall_field() -> MixtureEuler2d {
    let wall = AxisBoundaries::Ends {
        lower: BoundaryCondition::Reflecting,
        upper: BoundaryCondition::Reflecting,
    };
    field(
        material(30.0, 0.4),
        Boundaries2d { x: wall, y: wall },
        |i| {
            let a = 0.2 + 0.03 * (i % 6) as f64;
            ([a, 0.8, 0.0], [0.2, -0.1], 1.3)
        },
    )
}

#[test]
fn g0_forced_mixture_zero_source_preserves_native_bits_and_workspace_budget() {
    let mut native = wall_field();
    let mut flow = wall_field();
    let bytes = flow.required_forcing_workspace_bytes().unwrap();
    assert_eq!(
        bytes,
        flow.cell_count() * (flow.gas().species_count() + 3) * size_of::<f64>()
    );
    assert!(flow.forcing_workspace(bytes - 1).is_err());
    let mut workspace = work(&flow);
    let zero = vec![[0.0; 2]; flow.cell_count()];
    let expected = native
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let actual = flow
        .step_with_forcing(
            Time::new(0.001),
            StepControl::default(),
            BodyForce2d::Acceleration(&zero),
            CONTROL,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(actual.source_change, ForceChange2d::default());
    for (left, right) in flow.states().iter().zip(native.states()) {
        assert_eq!(left.to_bits(), right.to_bits());
    }
    assert_eq!(flow.boundary_impulses, native.boundary_impulses);
    assert_eq!(flow.time(), native.time());
    assert_eq!(workspace.workspace_bytes(), bytes);
}

#[test]
fn g4_forced_mixture_bad_sources_failure_and_last_checkpoint_preserve_prior_state() {
    let make = || {
        let mut flow = wall_field();
        flow.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        flow
    };
    let mut flow = make();
    let mut workspace = work(&flow);
    let before = flow.states().to_vec();
    let old_time = flow.time();
    let old_report = flow.last_step().copied();
    let old_loads = flow.boundary_impulses.clone();
    let malformed = vec![[0.0; 2]; flow.cell_count() - 1];
    for force in [
        BodyForce2d::UniformAcceleration([f64::NAN, 0.0]),
        BodyForce2d::UniformForceDensity([0.0, f64::INFINITY]),
        BodyForce2d::ForceDensity(&malformed),
    ] {
        assert!(
            flow.step_with_forcing(
                Time::new(0.001),
                StepControl::default(),
                force,
                CONTROL,
                &mut workspace,
                || false,
            )
            .is_err()
        );
    }
    assert!(
        flow.step_with_forcing(
            Time::new(0.01),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            BodyForce2d::UniformAcceleration([1e9, 0.0]),
            CONTROL,
            &mut workspace,
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), old_time);
    assert_eq!(flow.last_step().copied(), old_report);
    assert_eq!(flow.boundary_impulses, old_loads);

    let force = BodyForce2d::UniformAcceleration([-3.0, 1.0]);
    let mut reference = make();
    let mut reference_work = work(&reference);
    let mut calls = 0;
    let expected = reference
        .step_with_forcing(
            Time::new(0.001),
            StepControl::default(),
            force,
            CONTROL,
            &mut reference_work,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    assert_eq!(
        flow.step_with_forcing(
            Time::new(0.001),
            StepControl::default(),
            force,
            CONTROL,
            &mut workspace,
            || {
                polls += 1;
                polls == calls
            },
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), old_time);
    assert_eq!(flow.last_step().copied(), old_report);
    assert_eq!(flow.boundary_impulses, old_loads);
    let retried = flow
        .step_with_forcing(
            Time::new(0.001),
            StepControl::default(),
            force,
            CONTROL,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(retried, expected);
    assert_eq!(flow.states(), reference.states());
    assert_eq!(flow.boundary_impulses, reference.boundary_impulses);
}
