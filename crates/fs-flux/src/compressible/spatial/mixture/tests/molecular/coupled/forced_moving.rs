use super::*;
use crate::compressible::spatial::MeshMotion;
use crate::compressible::spatial::mixture::molecular::MovingMolecularStepReport3d;

fn moving_check(flow: &Coupled, report: &MovingMolecularStepReport3d, motion: MeshMotion) {
    check(&report.transport);
    close(
        report.transport.first.boundary_change.energy
            + report.transport.second.boundary_change.energy,
        report.wall_heat_j.iter().sum::<f64>() + report.molecular_wall_work_j.iter().sum::<f64>(),
        4e-11,
    );
    let mut rate: f64 = 0.0;
    let grid = flow.flow().grid();
    for index in 0..flow.flow().cell_count() {
        let p = flow
            .flow()
            .cell(index)
            .unwrap()
            .primitive(flow.flow().gas(), || false)
            .unwrap();
        assert!(
            p.velocity[0].hypot(p.velocity[1]).hypot(p.velocity[2])
                <= CAP.maximum_mach * p.sound_speed
        );
        let c = grid.coordinates(index);
        rate = rate.max(
            (0..3)
                .map(|a| {
                    ((p.velocity[a] - motion.velocity(grid, a, c[a]))
                        .abs()
                        .max((p.velocity[a] - motion.velocity(grid, a, c[a] + 1)).abs())
                        + p.sound_speed)
                        / grid.spacing()[a]
                })
                .sum(),
        );
    }
    assert!(
        rate * report.transport.elapsed.value()
            <= MolecularStepControl::default().advection.cfl * (1.0 + 1e-13)
    );
}

#[test]
fn g0_native_stationary_forcing_and_zero_force_moving_parity_include_wall_heat() {
    let force = BodyForce3d::UniformAcceleration([0.2, -0.1, 0.15]);
    let (mut expected, mut a_work) = wall_case();
    let (mut actual, mut b_work) = wall_case();
    let a = expected
        .step_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            force,
            CAP,
            &mut a_work,
            || false,
        )
        .unwrap();
    let b = actual
        .step_moving_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            MeshMotion::default(),
            force,
            CAP,
            &mut b_work,
            || false,
        )
        .unwrap();
    assert_eq!(a, b.transport);
    assert_eq!(expected.flow().states(), actual.flow().states());
    assert_eq!(loads(&expected), loads(&actual));
    assert_eq!(molecular_loads(&expected), molecular_loads(&actual));
    moving_check(&actual, &b, MeshMotion::default());
    assert!(b.wall_heat_j[4] > 0.0 && b.wall_heat_j[5] < 0.0);

    let motion = MeshMotion {
        lower_velocity: [-0.02, 0.0, 0.0],
        upper_velocity: [0.04, 0.0, 0.0],
    };
    let (mut expected, _) = wall_case();
    let (mut actual, mut scratch) = wall_case();
    let a = expected
        .step_moving(
            Time::new(0.002),
            MolecularStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
    let b = actual
        .step_moving_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            motion,
            BodyForce3d::UniformForceDensity([0.0; 3]),
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(expected.flow().states(), actual.flow().states());
    assert_eq!(loads(&expected), loads(&actual));
    assert_eq!(molecular_loads(&expected), molecular_loads(&actual));
}

#[test]
fn g1_forced_viscous_deformation_uses_actual_volumes_and_closes_all_three_momenta() {
    let value = [0.2, -0.1, 0.15];
    for axis in 0..3 {
        for density_force in [false, true] {
            let (mut flow, mut scratch) = wall_case();
            let bytes = flow.workspace_bytes() + scratch.workspace_bytes();
            let old_grid = *flow.flow().grid();
            let force = if density_force {
                BodyForce3d::UniformForceDensity(value)
            } else {
                BodyForce3d::UniformAcceleration(value)
            };
            let mut motion = MeshMotion::default();
            motion.lower_velocity[axis] = -0.02;
            motion.upper_velocity[axis] = 0.04;
            let r = flow
                .step_moving_with_forcing(
                    Time::new(0.002),
                    MolecularStepControl::default(),
                    motion,
                    force,
                    CAP,
                    &mut scratch,
                    || false,
                )
                .unwrap();
            moving_check(&flow, &r, motion);
            let dt = r.transport.elapsed.value();
            let measure = if density_force {
                0.5 * (old_grid.cell_volume() + flow.flow().grid().cell_volume())
                    * flow.flow().cell_count() as f64
            } else {
                r.transport.before.mass()
            };
            for a in 0..3 {
                close(
                    r.transport.source_change.impulse[a],
                    value[a] * measure * dt,
                    3e-12,
                );
                close(
                    r.transport.after.momentum[a]
                        - r.transport.before.momentum[a]
                        - r.transport.boundary_change.momentum[a],
                    r.transport.source_change.impulse[a],
                    3e-11,
                );
            }
            assert!(r.wall_heat_j[4] > 0.0 && r.wall_heat_j[5] < 0.0);
            assert!(r.molecular_wall_work_j.iter().any(|x| x.abs() > 1e-7));
            assert_eq!(flow.molecular().last_step(), Some(&r.transport.second));
            assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
            for side in BoundarySide::ALL {
                for q in flow.boundary_impulses(side) {
                    assert_eq!(q.species_masses, &[0.0; 3]);
                }
            }
            let next = flow
                .step_with_forcing(
                    Time::new(0.0005),
                    MolecularStepControl::default(),
                    force,
                    CAP,
                    &mut scratch,
                    || false,
                )
                .unwrap();
            check(&next);
        }
    }
}

fn accelerated_wave(n: usize) -> f64 {
    let h = 1.0 / n as f64;
    let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
    let initial_velocity = [0.3, -0.2, 0.15];
    let acceleration = [-0.4, 0.25, -0.3];
    let wave = |x: [f64; 3]| (2.0 * std::f64::consts::PI * x.iter().sum::<f64>()).sin();
    let z = std::f64::consts::PI * h;
    let amplitude = 0.15 * (z.sin() / z).powi(3);
    let original = field(equal_calorics(), grid, PERIODIC, |x| {
        let a = 0.5 + amplitude * wave(x);
        ([a, 1.0 - a, 0.0], initial_velocity, 1.0)
    });
    let molecular = prepare(&original, law(original.gas(), 0.02, 0.03, 0.01, 0.04));
    let (mut flow, mut scratch) = coupled(original, molecular);
    let end = 0.006;
    while flow.flow().time().value() < end {
        let time = flow.flow().time().value();
        let dt = (end - time).min(0.001);
        let midpoint =
            std::array::from_fn(|a| initial_velocity[a] + acceleration[a] * (time + 0.5 * dt));
        let motion = MeshMotion {
            lower_velocity: midpoint,
            upper_velocity: midpoint,
        };
        let r = flow
            .step_moving_with_forcing(
                Time::new(dt),
                MolecularStepControl::default(),
                motion,
                BodyForce3d::UniformAcceleration(acceleration),
                CAP,
                &mut scratch,
                || false,
            )
            .unwrap();
        moving_check(&flow, &r, motion);
        assert_eq!(r.transport.elapsed, Time::new(dt));
        assert_eq!(r.wall_heat_j, [0.0; 6]);
        assert_eq!(r.molecular_wall_work_j, [0.0; 6]);
        let exact_work: f64 = (0..3)
            .map(|a| r.transport.before.mass() * acceleration[a] * dt * midpoint[a])
            .sum();
        close(r.transport.source_change.work, exact_work, 2e-11);
    }
    let decay = (-0.02 * 12.0 * std::f64::consts::PI.powi(2) * end).exp();
    let mut error = 0.0;
    for index in 0..grid.len() {
        let q = flow.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        let exact = 0.5 + amplitude * decay * wave(grid.cell_center(index).unwrap());
        error += (q.partial_densities[0] - exact).abs();
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.pressure, 1.0, 3e-9);
        close(p.temperature, 1.0, 3e-9);
        for a in 0..3 {
            close(
                p.velocity[a],
                initial_velocity[a] + acceleration[a] * end,
                3e-10,
            );
        }
    }
    for a in 0..3 {
        close(
            flow.flow().grid().origin()[a],
            initial_velocity[a] * end + 0.5 * acceleration[a] * end * end,
            3e-13,
        );
    }
    error / grid.len() as f64
}

#[test]
fn g1_accelerating_mesh_species_diffusion_refines_without_spurious_thermal_work() {
    let coarse = accelerated_wave(4);
    let fine = accelerated_wave(8);
    assert!(
        fine < 0.4 * coarse,
        "accelerated moving molecular L1 {coarse} / {fine}"
    );
    assert!(fine < 2e-4, "accelerated moving molecular L1 {fine}");
}

#[test]
fn g4_forced_molecular_geometry_cancellation_preserves_binding_and_replays_all_phases() {
    let force = BodyForce3d::UniformAcceleration([0.2, -0.1, 0.15]);
    let motion = MeshMotion {
        lower_velocity: [-0.02, 0.0, 0.0],
        upper_velocity: [0.04, 0.0, 0.0],
    };
    let make = || {
        let (mut flow, mut work) = wall_case();
        flow.step_moving_with_forcing(
            Time::new(0.0005),
            MolecularStepControl::default(),
            motion,
            force,
            CAP,
            &mut work,
            || false,
        )
        .unwrap();
        (flow, work)
    };
    let (mut expected, mut expected_work) = make();
    let mut polls = 0;
    let accepted = expected
        .step_moving_with_forcing(
            Time::new(0.001),
            MolecularStepControl::default(),
            motion,
            force,
            CAP,
            &mut expected_work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let (mut flow, mut work) = make();
        let states = flow.flow().states().to_vec();
        let grid = *flow.flow().grid();
        let time = flow.flow().time();
        let report = *flow.last_step().unwrap();
        let molecular = *flow.molecular().last_step().unwrap();
        let accepted_loads = loads(&flow);
        let molecular_boundary = molecular_loads(&flow);
        let bytes = flow.workspace_bytes() + work.workspace_bytes();
        let mut seen = 0;
        assert_eq!(
            flow.step_moving_with_forcing(
                Time::new(0.001),
                MolecularStepControl::default(),
                motion,
                force,
                CAP,
                &mut work,
                || {
                    seen += 1;
                    seen == cutoff
                },
            ),
            Err(DissipationError::Cancelled)
        );
        assert_eq!(flow.flow().states(), states);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step(), Some(&report));
        assert_eq!(flow.molecular().last_step(), Some(&molecular));
        assert_eq!(loads(&flow), accepted_loads);
        assert_eq!(molecular_loads(&flow), molecular_boundary);
        assert_eq!(flow.workspace_bytes() + work.workspace_bytes(), bytes);
        let replay = flow
            .step_moving_with_forcing(
                Time::new(0.001),
                MolecularStepControl::default(),
                motion,
                force,
                CAP,
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(replay, accepted);
        assert_eq!(flow.flow().states(), expected.flow().states());
        assert_eq!(loads(&flow), loads(&expected));
        assert_eq!(molecular_loads(&flow), molecular_loads(&expected));
    }
}
