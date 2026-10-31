//! G1/G3/G4: analytic Basset impedance/settling, source changes and atomic memory.
use super::*;
use fs_mbd_material::creeping_flow::unsteady::{BassetResolution, UnsteadyStepControls};

fn resolution(lo: f64, hi: f64, intervals: usize) -> BassetResolution {
    BassetResolution {
        minimum_wavenumber_radius: lo,
        maximum_wavenumber_radius: hi,
        intervals,
    }
}

fn step_controls(dt: f64) -> UnsteadyStepControls {
    UnsteadyStepControls {
        duration_s: dt,
        maximum_energy_defect_j: 1e-18,
        maximum_momentum_defect_ns: 1e-16,
    }
}

#[test]
fn positive_memory_ladder_converges_to_the_analytic_unsteady_sphere_impedance() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let c = 6.0 * core::f64::consts::PI * 1.0 * RADIUS;
    let diffusion_time = RADIUS.powi(2) / 0.001;
    let added_mass = 0.5 * 1000.0 * 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    for omega in [1.0, 10.0, 100.0] {
        let history = c * (0.5_f64 * omega * diffusion_time).sqrt();
        let exact = [c + history, added_mass * omega + history];
        let mut previous = f64::INFINITY;
        for (lo, hi, n) in [
            (0.02, 40.0, 32),
            (0.001, 1000.0, 128),
            (0.0001, 10000.0, 256),
        ] {
            let motion = medium
                .unsteady_motion(&body, rest, resolution(lo, hi, n), || false)
                .unwrap();
            close(motion.added_mass_kg(), added_mass);
            assert!(
                motion
                    .branches()
                    .iter()
                    .all(|b| b.mass_kg > 0.0 && b.resistance_kg_s > 0.0)
            );
            let actual = motion.impedance_kg_s(omega).unwrap();
            let error = (actual[0] - exact[0]).hypot(actual[1] - exact[1]);
            assert!(
                error < previous * 0.15,
                "spectral ladder failed: {error} >= {previous}"
            );
            // This fixture's fine quadrature error is negligible compared with
            // the independently bounded two omitted continuum tails.
            assert!(error <= 1.001 * motion.omitted_tail_bound_kg_s(omega).unwrap());
            previous = error;
        }
        assert!(previous / exact[0].hypot(exact[1]) < 5e-4);
    }
}

#[test]
fn source_density_sets_buoyant_direction_and_actual_added_mass_acceleration() {
    for (solid_density, fluid_density) in [
        (500.0, 1000.0),
        (1000.0, 1000.0),
        (2000.0, 1000.0),
        (2000.0, 500.0),
    ] {
        let material = fluid(fluid_density, 1.0);
        let source = sphere(solid_density);
        let body = source.principal_body(1e-12, || false).unwrap();
        let medium =
            CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits())
                .unwrap();
        let rest = body
            .state_from_kinematics(Se3::identity(), Twist::zero())
            .unwrap();
        // The old quasi-steady model correctly refuses these density ratios.
        assert!(matches!(
            medium.step(&body, rest, controls(1e-10), || false),
            Err(CreepingFlowError::Unsteady { .. })
        ));
        let motion = medium
            .unsteady_motion(&body, rest, resolution(1e-4, 1000.0, 128), || false)
            .unwrap();
        let dt = 1e-10;
        let (next, audit) = motion.step(step_controls(dt), || false).unwrap();
        let velocity = next
            .state()
            .center_of_mass_velocity_world(body.properties())
            .unwrap()
            .z;
        let acceleration =
            -0.01 * (solid_density - fluid_density) / (solid_density + 0.5 * fluid_density);
        if acceleration == 0.0 {
            assert_eq!(next.state(), rest);
            assert!(
                next.fluid_velocities_world_m_s()
                    .iter()
                    .all(|u| *u == fs_mbd::Vec3::ZERO)
            );
            assert_eq!(audit.viscous_heat_j, 0.0);
        } else {
            assert!((velocity / (dt * acceleration) - 1.0).abs() < 3e-6);
            assert_eq!(velocity.is_sign_positive(), solid_density < fluid_density);
            assert!(audit.added_mass_energy_change_j > 0.0);
            assert!(audit.memory_energy_change_j > 0.0);
            assert!(audit.viscous_heat_j > 0.0);
        }
        let native = RigidBodyIntegrator::new(gravity())
            .step(rest, body.properties(), audit.fluid_wrench, dt)
            .unwrap();
        assert_eq!(native, audit.dynamics);
        assert_eq!(next.state(), native.state_after);
    }
}

/// Branch-cut inversion of V(s)=F/[s*(M*s+b*sqrt(s)+c)] from a quiescent start.
/// With alpha=b/sqrt(M*c), the normalized velocity is 1 minus this positive
/// decaying integral. Simpson quadrature is independent of the production
/// rational branch approximation and native midpoint stepping.
fn settling_fraction(time_s: f64, moving_mass: f64, c: f64, b: f64) -> f64 {
    let alpha = b / (moving_mass * c).sqrt();
    let scaled_time = time_s * c / moving_mass;
    let n = 8192;
    let dx = 12.0 / n as f64;
    let integrand =
        |x: f64| (-scaled_time * x * x).exp() / ((1.0 - x * x).powi(2) + alpha * alpha * x * x);
    let mut integral = integrand(0.0) + integrand(12.0);
    for i in 1..n {
        integral += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(i as f64 * dx);
    }
    1.0 - (2.0 * alpha / core::f64::consts::PI) * integral * dx / 3.0
}

#[test]
fn native_settling_converges_to_basset_history_with_closed_storage_and_heat() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let volume = 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    let moving_mass = 2500.0 * volume;
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let b = c * (RADIUS.powi(2) / 0.001).sqrt();
    let terminal = -1000.0 * volume * 0.01 / c;
    let target = 0.1;
    let expected = terminal * settling_fraction(target, moving_mass, c, b);
    let mut previous = f64::INFINITY;
    for (lo, hi, branches, steps) in [
        (0.01, 100.0, 32, 50),
        (0.001, 1000.0, 128, 200),
        (0.0001, 10000.0, 256, 800),
    ] {
        let mut motion = medium
            .unsteady_motion(&body, rest, resolution(lo, hi, branches), || false)
            .unwrap();
        let mut heat = 0.0;
        let mut stored = 0.0;
        let mut fluid_storage = 0.0;
        for _ in 0..steps {
            let (next, audit) = motion
                .step(step_controls(target / steps as f64), || false)
                .unwrap();
            heat += audit.viscous_heat_j;
            stored += audit.stored_energy_change_j;
            fluid_storage += audit.memory_energy_change_j + audit.added_mass_energy_change_j;
            assert!(audit.momentum_residual_ns.norm_squared().sqrt() < 1e-18);
            motion = next;
        }
        assert!(fluid_storage > 0.0 && heat > 0.0);
        assert!((stored + heat).abs() < 1e-11 * heat);
        let velocity = motion
            .state()
            .center_of_mass_velocity_world(body.properties())
            .unwrap()
            .z;
        let relative_error = (velocity / expected - 1.0).abs();
        assert!(
            relative_error < previous * 0.2,
            "history ladder: {relative_error} >= {previous}"
        );
        previous = relative_error;
    }
    assert!(previous < 3e-4);
    // Fluid history changes the trajectory, rather than just its diagnostics.
    let no_history = terminal * (1.0 - (-c * target / moving_mass).exp());
    assert!((expected / no_history - 1.0).abs() > 0.1);
}

#[test]
fn restart_and_late_cancellation_preserve_the_same_body_and_fluid_history() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let original = medium
        .unsteady_motion(&body, rest, resolution(0.001, 1000.0, 64), || false)
        .unwrap();
    let (first, _) = original.step(step_controls(0.001), || false).unwrap();
    let (whole, audit) = first.step(step_controls(0.001), || false).unwrap();
    let (resumed, repeated) = first.clone().step(step_controls(0.001), || false).unwrap();
    assert_eq!(whole.state(), resumed.state());
    assert_eq!(
        whole.fluid_velocities_world_m_s(),
        resumed.fluid_velocities_world_m_s()
    );
    assert_eq!(audit, repeated);
    let mut polls = 0;
    first
        .step(step_controls(0.001), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        first.step(step_controls(0.001), || {
            visited += 1;
            visited == polls
        }),
        Err(CreepingFlowError::Cancelled)
    ));
    let (retried, retried_audit) = first.step(step_controls(0.001), || false).unwrap();
    assert_eq!(retried.state(), whole.state());
    assert_eq!(
        retried.fluid_velocities_world_m_s(),
        whole.fluid_velocities_world_m_s()
    );
    assert_eq!(retried_audit, audit);
    assert_eq!(original.state(), rest);
    assert_eq!(original.time_s(), 0.0);
    assert!(
        medium
            .unsteady_motion(&body, whole.state(), resolution(0.001, 1000.0, 64), || {
                false
            })
            .is_err()
    );
}
