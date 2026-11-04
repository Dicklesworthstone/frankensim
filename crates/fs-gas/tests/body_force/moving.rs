//! Mechanical source work and moving-boundary work share one native interval.
use super::*;

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn mesh(count: usize, graded: bool) -> Mesh1d {
    if !graded {
        return Mesh1d::uniform(count, Length::new(1.0 / count as f64), 0.7).unwrap();
    }
    let faces: Vec<_> = (0..=count)
        .map(|i| (i as f64 / count as f64).powi(2))
        .collect();
    let volumes = faces.windows(2).map(|p| 0.7 * (p[1] - p[0])).collect();
    Mesh1d::new(faces, vec![0.7; count + 1], volumes).unwrap()
}

fn closed_frozen(count: usize, velocity: f64, graded: bool) -> Euler1d {
    let gas = frozen_gas();
    let states = (0..count)
        .map(|i| {
            gas.conserved(0.8 + 0.4 * (i as f64 + 0.5) / count as f64, velocity, 1.7)
                .unwrap()
        })
        .collect();
    Euler1d::new_on_mesh(gas, mesh(count, graded), states, CLOSED, || false).unwrap()
}

fn closed_mixture(count: usize, velocity: f64, graded: bool) -> MixtureEuler1d {
    let gas = mixture_gas();
    let states = (0..count)
        .map(|i| {
            let density = 0.8 + 0.4 * (i as f64 + 0.5) / count as f64;
            let first = if i < count / 2 { 0.8 } else { 0.3 };
            gas.conserved(
                &[first * density, (1.0 - first) * density, 0.0],
                velocity,
                2.7,
            )
            .unwrap()
        })
        .collect();
    MixtureEuler1d::new_on_mesh(gas, mesh(count, graded), states, CLOSED, || false).unwrap()
}

fn mixture_cells(flow: &MixtureEuler1d) -> Vec<(Vec<f64>, f64, f64)> {
    (0..flow.cell_count())
        .map(|i| {
            let q = flow.cell(i).unwrap();
            (q.partial_densities.to_vec(), q.momentum, q.energy)
        })
        .collect()
}

fn same_frozen(a: &Euler1d, b: &Euler1d) {
    assert_eq!(a.states(), b.states());
    assert_eq!(a.mesh(), b.mesh());
    assert_eq!(a.time(), b.time());
}

fn same_mixture(a: &MixtureEuler1d, b: &MixtureEuler1d) {
    assert_eq!(mixture_cells(a), mixture_cells(b));
    assert_eq!(a.mesh(), b.mesh());
    assert_eq!(a.time(), b.time());
}

#[test]
fn g1_accelerating_frozen_interval_preserves_internal_state_and_exact_source_work() {
    // Independent exact solution: uniform acceleration and midpoint wall speeds
    // translate every material point together. Density, pressure and temperature
    // are constant on each particle, even on this graded nonuniform mesh.
    let mut flow = closed_frozen(12, -0.5, true);
    let initial = flow.states().to_vec();
    let positions = flow.mesh().positions().to_vec();
    let volumes = flow.mesh().volumes().to_vec();
    let temperatures: Vec<_> = initial
        .iter()
        .map(|&q| flow.gas().primitive(q, || false).unwrap().temperature)
        .collect();
    let before = flow.totals();
    let (h, acceleration) = (0.00025, 100.0);
    let mut work = 0.0;
    let mut braking = false;
    for step in 0..80 {
        let velocity = -0.5 + acceleration * (step as f64 * h);
        let report = flow
            .step_with_moving_walls_and_body_force(
                Time::new(h),
                StepControl::default(),
                [velocity + 0.5 * acceleration * h; 2],
                BodyForce::UniformAcceleration(acceleration),
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(report.elapsed, Time::new(h));
        close(
            report.source_change.impulse,
            before.mass * acceleration * h,
            2e-12,
        );
        close(
            report.source_change.work,
            before.mass * acceleration * h * (velocity + 0.5 * acceleration * h),
            2e-11,
        );
        close(report.boundary_change.energy, 0.0, 2e-11);
        close(report.conservation_residual.energy, 0.0, 2e-11);
        assert_eq!(report.wall_heat, 0.0);
        braking |= report.source_change.work < 0.0;
        work += report.source_change.work;
    }
    assert!(
        braking,
        "an accelerating source must recover work while braking"
    );
    let time = flow.time().value();
    let velocity = -0.5 + acceleration * time;
    close(work, 0.5 * before.mass * (velocity * velocity - 0.25), 2e-9);
    close(flow.totals().energy - before.energy, work, 2e-9);
    assert_eq!(flow.mesh().volumes(), volumes);
    for (&position, &initial) in flow.mesh().positions().iter().zip(&positions) {
        close(
            position,
            initial - 0.5 * time + 0.5 * acceleration * time * time,
            2e-13,
        );
    }
    for (i, &q) in flow.states().iter().enumerate() {
        let p = flow.gas().primitive(q, || false).unwrap();
        close(q.density, initial[i].density, 2e-11);
        close(p.velocity, velocity, 2e-10);
        close(p.pressure, 1.7, 2e-10);
        close(p.temperature, temperatures[i], 2e-10);
    }
}

#[test]
fn g1_accelerating_mixture_retains_each_species_and_original_caloric_reference() {
    let mut flow = closed_mixture(12, -0.5, true);
    let initial = mixture_cells(&flow);
    let positions = flow.mesh().positions().to_vec();
    let volumes = flow.mesh().volumes().to_vec();
    let temperatures: Vec<_> = (0..flow.cell_count())
        .map(|i| {
            flow.gas()
                .primitive(flow.cell(i).unwrap(), || false)
                .unwrap()
                .temperature
        })
        .collect();
    let before = flow.totals();
    let (h, acceleration) = (0.00025, 100.0);
    let mut work = 0.0;
    let mut braking = false;
    for step in 0..80 {
        let velocity = -0.5 + acceleration * (step as f64 * h);
        let report = flow
            .step_with_moving_walls_and_body_force(
                Time::new(h),
                StepControl::default(),
                [velocity + 0.5 * acceleration * h; 2],
                BodyForce::UniformAcceleration(acceleration),
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(report.elapsed, Time::new(h));
        close(
            report.source_change.impulse,
            before.mass() * acceleration * h,
            2e-12,
        );
        close(
            report.source_change.work,
            before.mass() * acceleration * h * (velocity + 0.5 * acceleration * h),
            2e-11,
        );
        close(report.boundary_change.energy, 0.0, 2e-11);
        close(report.conservation_residual.energy, 0.0, 2e-11);
        assert!(
            report
                .boundary_change
                .species_masses
                .iter()
                .all(|m| *m == 0.0)
        );
        braking |= report.source_change.work < 0.0;
        work += report.source_change.work;
    }
    assert!(braking);
    let time = flow.time().value();
    let velocity = -0.5 + acceleration * time;
    close(
        work,
        0.5 * before.mass() * (velocity * velocity - 0.25),
        2e-9,
    );
    close(flow.totals().energy - before.energy, work, 2e-9);
    assert_eq!(flow.mesh().volumes(), volumes);
    for (&position, &initial) in flow.mesh().positions().iter().zip(&positions) {
        close(
            position,
            initial - 0.5 * time + 0.5 * acceleration * time * time,
            2e-13,
        );
    }
    for i in 0..flow.cell_count() {
        let q = flow.cell(i).unwrap();
        let p = flow.gas().primitive(q, || false).unwrap();
        for (actual, expected) in q.partial_densities.iter().zip(&initial[i].0) {
            close(*actual, *expected, 2e-11);
        }
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.velocity, velocity, 2e-10);
        close(p.pressure, 2.7, 2e-10);
        close(p.temperature, temperatures[i], 2e-10);
    }
}

#[test]
fn g1_compression_uses_both_source_volumes_and_independent_galilean_work() {
    let (h, density_force, boost) = (0.005, 4.0, 3.0);
    let walls = [0.2, -0.1];
    let shifted = [walls[0] + boost, walls[1] + boost];
    let mut frozen = closed_frozen(8, 0.4, false);
    let mut frozen_boost = closed_frozen(8, 0.4 + boost, false);
    let mut mixture = closed_mixture(8, 0.4, false);
    let mut mixture_boost = closed_mixture(8, 0.4 + boost, false);
    let a = frozen
        .step_with_moving_walls_and_body_force(
            Time::new(h),
            StepControl::default(),
            walls,
            BodyForce::UniformForceDensity(density_force),
            FORCE,
            || false,
        )
        .unwrap();
    let b = frozen_boost
        .step_with_moving_walls_and_body_force(
            Time::new(h),
            StepControl::default(),
            shifted,
            BodyForce::UniformForceDensity(density_force),
            FORCE,
            || false,
        )
        .unwrap();
    let c = mixture
        .step_with_moving_walls_and_body_force(
            Time::new(h),
            StepControl::default(),
            walls,
            BodyForce::UniformForceDensity(density_force),
            FORCE,
            || false,
        )
        .unwrap();
    let d = mixture_boost
        .step_with_moving_walls_and_body_force(
            Time::new(h),
            StepControl::default(),
            shifted,
            BodyForce::UniformForceDensity(density_force),
            FORCE,
            || false,
        )
        .unwrap();
    for elapsed in [a.elapsed, b.elapsed, c.elapsed, d.elapsed] {
        assert_eq!(elapsed, Time::new(h));
    }
    // V(t) = A*(1+(w_R-w_L)*t). A held force density's impulse uses
    // its exact time-integrated physical volume, not the old volume twice.
    let impulse = density_force * h * 0.7 * (1.0 + 0.5 * (walls[1] - walls[0]) * h);
    for source in [
        a.source_change,
        b.source_change,
        c.source_change,
        d.source_change,
    ] {
        close(source.impulse, impulse, 2e-14);
    }
    close(
        b.source_change.work - a.source_change.work,
        boost * a.source_change.impulse,
        2e-12,
    );
    close(
        b.boundary_change.energy - a.boundary_change.energy,
        boost * a.boundary_change.momentum,
        2e-11,
    );
    close(
        d.source_change.work - c.source_change.work,
        boost * c.source_change.impulse,
        2e-12,
    );
    close(
        d.boundary_change.energy - c.boundary_change.energy,
        boost * c.boundary_change.momentum,
        2e-11,
    );
    for (r, w) in [(&a, walls), (&b, shifted)] {
        close(
            r.left_boundary_change.energy,
            w[0] * r.left_boundary_change.momentum,
            2e-12,
        );
        close(
            r.right_boundary_change.energy,
            w[1] * r.right_boundary_change.momentum,
            2e-12,
        );
        close(
            r.after.energy - r.before.energy,
            r.boundary_change.energy + r.source_change.work,
            2e-12,
        );
        close(
            r.after.momentum - r.before.momentum,
            r.boundary_change.momentum + r.source_change.impulse,
            2e-12,
        );
        assert_eq!(r.boundary_change.mass, 0.0);
        assert_eq!(r.wall_heat, 0.0);
        assert_eq!(r.wall_work, 0.0);
    }
    for (r, w) in [(&c, walls), (&d, shifted)] {
        close(
            r.left_boundary_change.energy,
            w[0] * r.left_boundary_change.momentum,
            2e-12,
        );
        close(
            r.right_boundary_change.energy,
            w[1] * r.right_boundary_change.momentum,
            2e-12,
        );
        close(
            r.after.energy - r.before.energy,
            r.boundary_change.energy + r.source_change.work,
            2e-12,
        );
        close(
            r.after.momentum - r.before.momentum,
            r.boundary_change.momentum + r.source_change.impulse,
            2e-12,
        );
        for k in 0..3 {
            close(r.after.species_masses[k], r.before.species_masses[k], 2e-13);
            assert_eq!(r.boundary_change.species_masses[k], 0.0);
        }
    }
    assert!(a.boundary_change.energy > 0.0 && c.boundary_change.energy > 0.0);
    for i in 0..8 {
        let pa = frozen
            .gas()
            .primitive(frozen.states()[i], || false)
            .unwrap();
        let pb = frozen_boost
            .gas()
            .primitive(frozen_boost.states()[i], || false)
            .unwrap();
        close(pb.pressure, pa.pressure, 2e-10);
        close(pb.velocity, pa.velocity + boost, 2e-10);
        let pc = mixture
            .gas()
            .primitive(mixture.cell(i).unwrap(), || false)
            .unwrap();
        let pd = mixture_boost
            .gas()
            .primitive(mixture_boost.cell(i).unwrap(), || false)
            .unwrap();
        close(pd.pressure, pc.pressure, 2e-10);
        close(pd.velocity, pc.velocity + boost, 2e-10);
        assert_eq!(mixture.cell(i).unwrap().partial_densities[2], 0.0);
    }
    for (x, shifted_x) in frozen
        .mesh()
        .positions()
        .iter()
        .zip(frozen_boost.mesh().positions())
    {
        close(*shifted_x - *x, boost * h, 3e-14);
    }
}

#[test]
fn g0_zero_source_and_zero_motion_keep_original_native_paths_exactly() {
    let walls = [0.15, -0.05];
    let mut frozen = closed_frozen(8, 0.2, false);
    let mut frozen_zero = closed_frozen(8, 0.2, false);
    let mut mixture = closed_mixture(8, 0.2, false);
    let mut mixture_zero = closed_mixture(8, 0.2, false);
    for _ in 0..4 {
        let a = frozen
            .step_with_moving_walls(Time::new(0.001), StepControl::default(), walls, || false)
            .unwrap();
        let b = frozen_zero
            .step_with_moving_walls_and_body_force(
                Time::new(0.001),
                StepControl::default(),
                walls,
                BodyForce::ForceDensity(&[0.0; 8]),
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(a, b);
        same_frozen(&frozen, &frozen_zero);
        let a = mixture
            .step_with_moving_walls(Time::new(0.001), StepControl::default(), walls, || false)
            .unwrap();
        let b = mixture_zero
            .step_with_moving_walls_and_body_force(
                Time::new(0.001),
                StepControl::default(),
                walls,
                BodyForce::Acceleration(&[0.0; 8]),
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(a, b);
        same_mixture(&mixture, &mixture_zero);
    }
    let mut frozen = closed_frozen(8, 0.2, false);
    let mut frozen_still = closed_frozen(8, 0.2, false);
    let mut mixture = closed_mixture(8, 0.2, false);
    let mut mixture_still = closed_mixture(8, 0.2, false);
    let force = BodyForce::UniformAcceleration(0.5);
    for _ in 0..4 {
        let a = frozen
            .step_with_body_force(
                Time::new(0.001),
                StepControl::default(),
                force,
                FORCE,
                || false,
            )
            .unwrap();
        let b = frozen_still
            .step_with_moving_walls_and_body_force(
                Time::new(0.001),
                StepControl::default(),
                [0.0; 2],
                force,
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(a, b);
        same_frozen(&frozen, &frozen_still);
        let a = mixture
            .step_with_body_force(
                Time::new(0.001),
                StepControl::default(),
                force,
                FORCE,
                || false,
            )
            .unwrap();
        let b = mixture_still
            .step_with_moving_walls_and_body_force(
                Time::new(0.001),
                StepControl::default(),
                [0.0; 2],
                force,
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(a, b);
        same_mixture(&mixture, &mixture_still);
    }
}

#[test]
fn g4_forced_ale_retries_publish_only_the_accepted_interval_and_work() {
    let force = BodyForce::UniformAcceleration(500.0);
    let walls = [0.05, -0.02];
    let control = StepControl {
        reconstruction: fs_gas::dynamics::Reconstruction::FirstOrder,
        ..StepControl::default()
    };
    let mut frozen = closed_frozen(4, 0.0, false);
    let mass = frozen.totals().mass;
    let a = frozen
        .step_with_moving_walls_and_body_force(Time::new(0.2), control, walls, force, FORCE, || {
            false
        })
        .unwrap();
    assert!(a.halvings > 0);
    assert!(a.maximum_cfl <= control.cfl * (1.0 + 64.0 * f64::EPSILON));
    close(
        a.source_change.impulse,
        500.0 * mass * a.elapsed.value(),
        2e-11,
    );
    close(a.conservation_residual.energy, 0.0, 2e-10);
    let mut reference = closed_frozen(4, 0.0, false);
    let b = reference
        .step_with_moving_walls_and_body_force(a.elapsed, control, walls, force, FORCE, || false)
        .unwrap();
    assert_eq!(b.halvings, 0);
    assert_eq!(frozen.mesh(), reference.mesh());
    assert_eq!(frozen.time(), reference.time());
    for (q, expected) in frozen.states().iter().zip(reference.states()) {
        close(q.density, expected.density, 2e-11);
        close(q.momentum, expected.momentum, 2e-11);
        close(q.energy, expected.energy, 2e-10);
    }
    close(a.source_change.work, b.source_change.work, 2e-10);
    let mut mixture = closed_mixture(4, 0.0, false);
    let mass = mixture.totals().mass();
    let a = mixture
        .step_with_moving_walls_and_body_force(Time::new(0.2), control, walls, force, FORCE, || {
            false
        })
        .unwrap();
    assert!(a.halvings > 0);
    assert!(a.maximum_cfl <= control.cfl * (1.0 + 64.0 * f64::EPSILON));
    close(
        a.source_change.impulse,
        500.0 * mass * a.elapsed.value(),
        2e-11,
    );
    close(a.conservation_residual.energy, 0.0, 2e-10);
    let mut reference = closed_mixture(4, 0.0, false);
    let b = reference
        .step_with_moving_walls_and_body_force(a.elapsed, control, walls, force, FORCE, || false)
        .unwrap();
    assert_eq!(b.halvings, 0);
    assert_eq!(mixture.mesh(), reference.mesh());
    assert_eq!(mixture.time(), reference.time());
    for i in 0..4 {
        let q = mixture.cell(i).unwrap();
        let expected = reference.cell(i).unwrap();
        for k in 0..3 {
            close(q.partial_densities[k], expected.partial_densities[k], 2e-11);
        }
        close(q.momentum, expected.momentum, 2e-11);
        close(q.energy, expected.energy, 2e-10);
    }
    close(a.source_change.work, b.source_change.work, 2e-10);
}

#[test]
fn g4_cancelled_or_refused_ale_sources_leave_cells_geometry_and_time_unchanged() {
    let walls = [0.15, -0.05];
    let force = BodyForce::ForceDensity(&[0.5, -0.8, 0.2, 0.7, -0.1, 0.9, -0.4, 0.6]);
    let mut frozen_reference = closed_frozen(8, 0.2, false);
    let mut frozen_polls = 0;
    let expected = frozen_reference
        .step_with_moving_walls_and_body_force(
            Time::new(0.001),
            StepControl::default(),
            walls,
            force,
            FORCE,
            || {
                frozen_polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, frozen_polls / 2, frozen_polls] {
        let mut field = closed_frozen(8, 0.2, false);
        let initial = field.states().to_vec();
        let geometry = field.mesh().clone();
        let mut polls = 0;
        let error = field
            .step_with_moving_walls_and_body_force(
                Time::new(0.001),
                StepControl::default(),
                walls,
                force,
                FORCE,
                || {
                    polls += 1;
                    polls >= cutoff
                },
            )
            .unwrap_err();
        assert_eq!(error, fs_gas::dynamics::Error::Cancelled);
        assert_eq!(field.states(), initial);
        assert_eq!(field.mesh(), &geometry);
        assert_eq!(field.time(), Time::new(0.0));
        let retried = field
            .step_with_moving_walls_and_body_force(
                Time::new(0.001),
                StepControl::default(),
                walls,
                force,
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(retried.elapsed, expected.elapsed);
        close(
            retried.source_change.work,
            expected.source_change.work,
            2e-12,
        );
        for (q, reference) in field.states().iter().zip(frozen_reference.states()) {
            close(q.density, reference.density, 2e-12);
            close(q.momentum, reference.momentum, 2e-12);
            close(q.energy, reference.energy, 2e-12);
        }
        assert_eq!(field.mesh(), frozen_reference.mesh());
    }
    let mut mixture_reference = closed_mixture(8, 0.2, false);
    let mut mixture_polls = 0;
    let expected = mixture_reference
        .step_with_moving_walls_and_body_force(
            Time::new(0.001),
            StepControl::default(),
            walls,
            force,
            FORCE,
            || {
                mixture_polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, mixture_polls / 2, mixture_polls] {
        let mut field = closed_mixture(8, 0.2, false);
        let initial = mixture_cells(&field);
        let geometry = field.mesh().clone();
        let mut polls = 0;
        let error = field
            .step_with_moving_walls_and_body_force(
                Time::new(0.001),
                StepControl::default(),
                walls,
                force,
                FORCE,
                || {
                    polls += 1;
                    polls >= cutoff
                },
            )
            .unwrap_err();
        assert_eq!(error, fs_gas::multicomponent::Error::Cancelled);
        assert_eq!(mixture_cells(&field), initial);
        assert_eq!(field.mesh(), &geometry);
        assert_eq!(field.time(), Time::new(0.0));
        let retried = field
            .step_with_moving_walls_and_body_force(
                Time::new(0.001),
                StepControl::default(),
                walls,
                force,
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(retried.elapsed, expected.elapsed);
        close(
            retried.source_change.work,
            expected.source_change.work,
            2e-12,
        );
        for i in 0..8 {
            let q = field.cell(i).unwrap();
            let reference = mixture_reference.cell(i).unwrap();
            for k in 0..3 {
                close(
                    q.partial_densities[k],
                    reference.partial_densities[k],
                    2e-12,
                );
            }
            close(q.momentum, reference.momentum, 2e-12);
            close(q.energy, reference.energy, 2e-12);
        }
        assert_eq!(field.mesh(), mixture_reference.mesh());
    }
    let bad_controls = StepControl {
        maximum_halvings: 0,
        ..StepControl::default()
    };
    for force in [
        BodyForce::ForceDensity(&[1.0; 7]),
        BodyForce::UniformAcceleration(f64::NAN),
        BodyForce::UniformAcceleration(1e8),
    ] {
        let mut field = closed_frozen(8, 0.2, false);
        let original = closed_frozen(8, 0.2, false);
        assert!(
            field
                .step_with_moving_walls_and_body_force(
                    Time::new(0.01),
                    bad_controls,
                    walls,
                    force,
                    FORCE,
                    || false,
                )
                .is_err()
        );
        same_frozen(&field, &original);
        let mut field = closed_mixture(8, 0.2, false);
        let original = closed_mixture(8, 0.2, false);
        assert!(
            field
                .step_with_moving_walls_and_body_force(
                    Time::new(0.01),
                    bad_controls,
                    walls,
                    force,
                    FORCE,
                    || false,
                )
                .is_err()
        );
        same_mixture(&field, &original);
    }
}
