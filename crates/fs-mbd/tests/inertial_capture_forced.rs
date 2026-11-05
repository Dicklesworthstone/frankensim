//! G1/G3/G4: inertial intake and external loads on one physical trajectory.

use fs_mbd::Vec3;
use fs_mbd::rocket::exhaust::{ExhaustEnvironment, ParticleCaptureForcing, ParticleCaptureStep};
use fs_mbd::rocket::{RocketError, RocketState};

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 4.0e-13 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn close_vec(actual: Vec3, expected: Vec3, scale: f64) {
    close(actual.x, expected.x, scale);
    close(actual.y, expected.y, scale);
    close(actual.z, expected.z, scale);
}

fn state() -> RocketState {
    RocketState {
        position_m: Vec3::new(1.0, 2.0, -3.0),
        velocity_m_s: Vec3::new(12.0, -6.0, 3.0),
        dry_mass_kg: 1.0,
        propellant_mass_kg: 1.0,
    }
}

fn stream() -> ParticleCaptureForcing {
    ParticleCaptureForcing {
        mass_flow_kg_s: 2.0,
        velocity_world_m_s: Vec3::new(-2.0, 4.0, 1.0),
        internal_energy_rate_w: -7.0,
    }
}

fn environment() -> ExhaustEnvironment {
    ExhaustEnvironment {
        gravity_m_s2: Vec3::new(-0.3, -9.8, 1.2),
        external_force_n: Vec3::new(3.0, -2.0, 4.0),
    }
}

fn balances(step: ParticleCaptureStep, scale: f64) {
    close(step.mass_residual_kg, 0.0, scale);
    close_vec(step.momentum_residual_kg_m_s, Vec3::ZERO, scale);
    close_vec(step.angular_momentum_residual_kg_m2_s, Vec3::ZERO, scale);
    close(step.kinetic_energy_residual_j, 0.0, scale);
}

// An independent conservation-law integration, with no analytic log kernels.
#[derive(Clone, Copy)]
struct Reference {
    mass: f64,
    position: Vec3,
    momentum: Vec3,
    heat: f64,
    incoming_angular: Vec3,
    external_angular: Vec3,
    force_work: f64,
    gravity_work: f64,
}

impl Reference {
    fn add_scaled(self, rhs: Self, h: f64) -> Self {
        Self {
            mass: self.mass + h * rhs.mass,
            position: self.position.add(rhs.position.scale(h)),
            momentum: self.momentum.add(rhs.momentum.scale(h)),
            heat: self.heat + h * rhs.heat,
            incoming_angular: self.incoming_angular.add(rhs.incoming_angular.scale(h)),
            external_angular: self.external_angular.add(rhs.external_angular.scale(h)),
            force_work: self.force_work + h * rhs.force_work,
            gravity_work: self.gravity_work + h * rhs.gravity_work,
        }
    }

    fn rhs(self, source: ParticleCaptureForcing, loads: ExhaustEnvironment) -> Self {
        let velocity = self.momentum.scale(self.mass.recip());
        let influx = source.velocity_world_m_s.scale(source.mass_flow_kg_s);
        let gravity = loads.gravity_m_s2.scale(self.mass);
        let force = loads.external_force_n.add(gravity);
        Self {
            mass: source.mass_flow_kg_s,
            position: velocity,
            momentum: influx.add(force),
            heat: 0.5
                * source.mass_flow_kg_s
                * source.velocity_world_m_s.sub(velocity).norm_squared(),
            incoming_angular: self.position.cross(influx),
            external_angular: self.position.cross(force),
            force_work: loads.external_force_n.dot(velocity),
            gravity_work: gravity.dot(velocity),
        }
    }
}

#[test]
fn forced_capture_matches_conservation_law_integration_and_exact_momentum() {
    let initial = state();
    let source = stream();
    let loads = environment();
    let h = 1.5;
    let result = initial
        .advance_with_inertial_capture_forced(source, loads, h, || false)
        .unwrap();
    let m0 = initial.dry_mass_kg + initial.propellant_mass_kg;
    let dm = source.mass_flow_kg_s * h;
    // This endpoint follows directly from P'=q*u+F+(M0+q*t)*g.
    let final_momentum = initial
        .velocity_m_s
        .scale(m0)
        .add(source.velocity_world_m_s.scale(dm))
        .add(loads.external_force_n.scale(h))
        .add(loads.gravity_m_s2.scale(h * (m0 + 0.5 * dm)));
    close_vec(
        result.state_after.velocity_m_s,
        final_momentum.scale(1.0 / (m0 + dm)),
        20.0,
    );
    let mut reference = Reference {
        mass: m0,
        position: initial.position_m,
        momentum: initial.velocity_m_s.scale(m0),
        heat: 0.0,
        incoming_angular: Vec3::ZERO,
        external_angular: Vec3::ZERO,
        force_work: 0.0,
        gravity_work: 0.0,
    };
    let dt = h / 2048.0;
    for _ in 0..2048 {
        let k1 = reference.rhs(source, loads);
        let k2 = reference.add_scaled(k1, 0.5 * dt).rhs(source, loads);
        let k3 = reference.add_scaled(k2, 0.5 * dt).rhs(source, loads);
        let k4 = reference.add_scaled(k3, dt).rhs(source, loads);
        reference = reference
            .add_scaled(k1, dt / 6.0)
            .add_scaled(k2, dt / 3.0)
            .add_scaled(k3, dt / 3.0)
            .add_scaled(k4, dt / 6.0);
    }
    close_vec(result.state_after.position_m, reference.position, 50.0);
    close(result.thermalization_energy_j, reference.heat, 1000.0);
    close_vec(
        result.incoming_angular_momentum_world_kg_m2_s,
        reference.incoming_angular,
        1000.0,
    );
    close_vec(
        result.external_angular_momentum_world_kg_m2_s,
        reference.external_angular,
        1000.0,
    );
    close(result.external_work_j, reference.force_work, 1000.0);
    close(result.gravity_work_j, reference.gravity_work, 1000.0);
    close(
        result.retained_internal_energy_j,
        reference.heat - 7.0 * h,
        1000.0,
    );
    close(
        result.incoming_energy_world_j,
        -7.0 * h + 0.5 * dm * source.velocity_world_m_s.norm_squared(),
        1000.0,
    );
    balances(result, 1000.0);
    assert_eq!(initial, state());
}

#[test]
fn force_free_zero_rate_and_zero_duration_have_their_exact_limits() {
    let initial = state();
    assert_eq!(
        initial
            .advance_with_inertial_capture(stream(), 1.5, || false)
            .unwrap(),
        initial
            .advance_with_inertial_capture_forced(stream(), ExhaustEnvironment::ZERO, 1.5, || false)
            .unwrap(),
    );
    let loads = environment();
    let source = ParticleCaptureForcing {
        mass_flow_kg_s: 0.0,
        internal_energy_rate_w: 0.0,
        ..stream()
    };
    let h = 0.75;
    let m = initial.dry_mass_kg + initial.propellant_mass_kg;
    let acceleration = loads
        .gravity_m_s2
        .add(loads.external_force_n.scale(m.recip()));
    let displacement = initial
        .velocity_m_s
        .scale(h)
        .add(acceleration.scale(0.5 * h * h));
    let result = initial
        .advance_with_inertial_capture_forced(source, loads, h, || false)
        .unwrap();
    close_vec(
        result.state_after.velocity_m_s,
        initial.velocity_m_s.add(acceleration.scale(h)),
        20.0,
    );
    close_vec(
        result.state_after.position_m,
        initial.position_m.add(displacement),
        20.0,
    );
    close_vec(
        result.external_momentum_world_kg_m_s,
        acceleration.scale(m * h),
        50.0,
    );
    close_vec(
        result.external_angular_momentum_world_kg_m2_s,
        initial.position_m.cross(acceleration.scale(m * h)).add(
            initial
                .velocity_m_s
                .cross(acceleration)
                .scale(0.5 * m * h * h),
        ),
        200.0,
    );
    close(
        result.external_work_j,
        loads.external_force_n.dot(displacement),
        200.0,
    );
    close(
        result.gravity_work_j,
        loads.gravity_m_s2.dot(displacement) * m,
        200.0,
    );
    assert_eq!(result.thermalization_energy_j, 0.0);
    assert_eq!(result.incoming_energy_world_j, 0.0);
    balances(result, 200.0);
    let zero = initial
        .advance_with_inertial_capture_forced(stream(), loads, 0.0, || false)
        .unwrap();
    assert_eq!(zero.state_after, initial);
    assert_eq!(zero.external_momentum_world_kg_m_s, Vec3::ZERO);
    assert_eq!(zero.external_work_j, 0.0);
    assert_eq!(zero.gravity_work_j, 0.0);
    assert_eq!(zero.thermalization_energy_j, 0.0);
    balances(zero, 1.0);
}

#[test]
fn one_held_interval_composes_with_its_physical_fluxes_and_material_credit() {
    let whole = state()
        .advance_with_inertial_capture_forced(stream(), environment(), 3.0, || false)
        .unwrap();
    let mut current = state();
    let (mut heat, mut force_work, mut gravity_work) = (0.0, 0.0, 0.0);
    let (mut incoming_angular, mut external_angular) = (Vec3::ZERO, Vec3::ZERO);
    for _ in 0..12 {
        let part = current
            .advance_with_inertial_capture_forced(stream(), environment(), 0.25, || false)
            .unwrap();
        current = part.state_after;
        heat += part.thermalization_energy_j;
        force_work += part.external_work_j;
        gravity_work += part.gravity_work_j;
        incoming_angular = incoming_angular.add(part.incoming_angular_momentum_world_kg_m2_s);
        external_angular = external_angular.add(part.external_angular_momentum_world_kg_m2_s);
        balances(part, 1000.0);
    }
    close_vec(current.position_m, whole.state_after.position_m, 100.0);
    close_vec(current.velocity_m_s, whole.state_after.velocity_m_s, 40.0);
    assert_eq!(
        current.propellant_mass_kg,
        whole.state_after.propellant_mass_kg
    );
    close(heat, whole.thermalization_energy_j, 2000.0);
    close(force_work, whole.external_work_j, 2000.0);
    close(gravity_work, whole.gravity_work_j, 2000.0);
    close_vec(
        incoming_angular,
        whole.incoming_angular_momentum_world_kg_m2_s,
        2000.0,
    );
    close_vec(
        external_angular,
        whole.external_angular_momentum_world_kg_m2_s,
        2000.0,
    );
    balances(whole, 2000.0);
}

#[test]
fn forced_intake_heat_and_work_obey_galilean_and_spatial_covariance() {
    let initial = state();
    let source = stream();
    let loads = environment();
    let h = 1.5;
    let original = initial
        .advance_with_inertial_capture_forced(source, loads, h, || false)
        .unwrap();
    let rotate = |v: Vec3| Vec3::new(v.z, v.x, v.y);
    let boost = Vec3::new(7.0, -3.0, 2.0);
    let origin = Vec3::new(-4.0, 1.0, 5.0);
    let transformed = RocketState {
        position_m: rotate(initial.position_m).add(origin),
        velocity_m_s: rotate(initial.velocity_m_s).add(boost),
        ..initial
    }
    .advance_with_inertial_capture_forced(
        ParticleCaptureForcing {
            velocity_world_m_s: rotate(source.velocity_world_m_s).add(boost),
            ..source
        },
        ExhaustEnvironment {
            gravity_m_s2: rotate(loads.gravity_m_s2),
            external_force_n: rotate(loads.external_force_n),
        },
        h,
        || false,
    )
    .unwrap();
    close_vec(
        transformed.state_after.position_m,
        rotate(original.state_after.position_m)
            .add(origin)
            .add(boost.scale(h)),
        100.0,
    );
    close_vec(
        transformed.state_after.velocity_m_s,
        rotate(original.state_after.velocity_m_s).add(boost),
        40.0,
    );
    close(
        transformed.thermalization_energy_j,
        original.thermalization_energy_j,
        1000.0,
    );
    close(
        transformed.external_work_j,
        original.external_work_j + boost.dot(rotate(loads.external_force_n).scale(h)),
        1000.0,
    );
    let gravity_impulse = loads.gravity_m_s2.scale(
        h * (initial.dry_mass_kg + initial.propellant_mass_kg + 0.5 * source.mass_flow_kg_s * h),
    );
    close(
        transformed.gravity_work_j,
        original.gravity_work_j + boost.dot(rotate(gravity_impulse)),
        1000.0,
    );
    close_vec(
        transformed.external_momentum_world_kg_m_s,
        rotate(original.external_momentum_world_kg_m_s),
        100.0,
    );
    balances(transformed, 2000.0);
}

#[test]
fn tiny_intake_keeps_heat_when_force_and_gravity_almost_balance() {
    let q = 2.0_f64.powi(-48);
    let initial = RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::ZERO,
        dry_mass_kg: 1.0,
        propellant_mass_kg: 0.0,
    };
    let result = initial
        .advance_with_inertial_capture_forced(
            ParticleCaptureForcing {
                mass_flow_kg_s: q,
                velocity_world_m_s: Vec3::ZERO,
                internal_energy_rate_w: 0.0,
            },
            ExhaustEnvironment {
                gravity_m_s2: Vec3::new(1.0, 0.0, 0.0),
                external_force_n: Vec3::new(-(1.0 + q), 0.0, 0.0),
            },
            1.0,
            || false,
        )
        .unwrap();
    close(result.state_after.velocity_m_s.x, -0.5 * q / (1.0 + q), q);
    close(result.state_after.position_m.x, -q / 3.0, q);
    // v(t)=q*(t²/2-t)/(1+q*t); Q=q³/15+O(q⁴).
    close(result.thermalization_energy_j, q * q * q / 15.0, q * q * q);
    assert!(result.thermalization_energy_j > 0.0);
    close(result.kinetic_energy_residual_j, 0.0, q * q);
}

#[test]
fn dominant_intake_retains_heat_when_the_initial_mass_ratio_underflows() {
    let initial = RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::new(2.0, -1.0, 0.0),
        dry_mass_kg: 1.0e-100,
        propellant_mass_kg: 0.0,
    };
    let result = initial
        .advance_with_inertial_capture_forced(
            ParticleCaptureForcing {
                mass_flow_kg_s: 1.0e300,
                velocity_world_m_s: Vec3::ZERO,
                internal_energy_rate_w: 0.0,
            },
            ExhaustEnvironment {
                gravity_m_s2: Vec3::ZERO,
                external_force_n: Vec3::new(1.0, 0.0, 0.0),
            },
            1.0,
            || false,
        )
        .unwrap();
    // The applied-force work and final kinetic energy are O(1e-300),
    // while capture thermalizes the original 2.5e-100 J of bulk motion.
    assert_eq!(
        initial.dry_mass_kg / result.state_after.propellant_mass_kg,
        0.0
    );
    close(result.thermalization_energy_j, 2.5e-100, 2.5e-100);
    close(result.state_after.velocity_m_s.x, 1.0e-300, 1.0e-300);
    close(result.external_work_j, 1.0e-300, 1.0e-300);
}

#[test]
fn invalid_loads_and_cancellation_leave_the_borrowed_state_unchanged() {
    let initial = state();
    let mut polls = 0;
    initial
        .advance_with_inertial_capture_forced(stream(), environment(), 1.0, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        assert_eq!(
            initial.advance_with_inertial_capture_forced(stream(), environment(), 1.0, || {
                count += 1;
                count == stop
            }),
            Err(RocketError::Cancelled),
        );
        assert_eq!(initial, state());
    }
    for gravity in [f64::NAN, f64::INFINITY] {
        assert!(matches!(
            initial.advance_with_inertial_capture_forced(
                stream(),
                ExhaustEnvironment {
                    gravity_m_s2: Vec3::new(gravity, 0.0, 0.0),
                    ..environment()
                },
                1.0,
                || false
            ),
            Err(RocketError::NonFiniteInput(_)),
        ));
    }
    assert!(matches!(
        initial.advance_with_inertial_capture_forced(stream(), environment(), -1.0, || false),
        Err(RocketError::InvalidInput(_)),
    ));
    assert!(matches!(
        initial.advance_with_inertial_capture_forced(
            stream(),
            ExhaustEnvironment {
                external_force_n: Vec3::new(f64::MAX, 0.0, 0.0),
                ..environment()
            },
            2.0,
            || false
        ),
        Err(RocketError::Unrepresentable(_)),
    ));
    assert_eq!(initial, state());
}
