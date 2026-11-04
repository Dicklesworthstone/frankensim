use super::*;
use crate::UnitQuaternion;
use crate::rocket::RocketState;
use crate::rocket::exhaust::{ExhaustEnvironment, ExhaustForcing};
use crate::rocket::rigid::nozzles::BalancedNozzleGroup;
use crate::rocket::rigid::tensor::{SpecificInertia, TensorRocketInertia};
use fs_ga::Mat3;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn controls(h: f64) -> RigidRocketSolve {
    RigidRocketSolve {
        maximum_substep_s: h,
        maximum_substeps: 10_000,
    }
}
fn body(propellant: f64, specific: SpecificInertia, omega: Vec3) -> TensorRocketState {
    let translation = RocketState {
        position_m: Vec3::new(2.0, -3.0, 1.0),
        velocity_m_s: Vec3::new(0.4, -0.2, 0.6),
        dry_mass_kg: 3.0,
        propellant_mass_kg: propellant,
    };
    let inertia = TensorRocketInertia::new(
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0],
        },
        specific,
    )
    .unwrap();
    let spin = inertia
        .mass_properties(translation)
        .unwrap()
        .angular_momentum_body_checked(omega)
        .unwrap();
    TensorRocketState::new(translation, UnitQuaternion::IDENTITY, spin, inertia).unwrap()
}
fn forcing(nozzles: PreparedNozzleArray) -> NozzleTransportForcing {
    NozzleTransportForcing {
        nozzles,
        local_relative_energy_rate_w: nozzles.relative_kinetic_power_w() + 17.0,
        gravity_world_m_s2: Vec3::ZERO,
    }
}
fn central() -> PreparedNozzleArray {
    PreparedNozzleArray::new(
        &[BalancedNozzleGroup::Central {
            mass_flow_kg_s: 0.7,
            exhaust_velocity_body_m_s: Vec3::new(-3.0, 2.0, 1.0),
        }],
        || false,
    )
    .unwrap()
}

#[test]
fn nonrotating_transport_matches_original_forced_exhaust_through_burnout() {
    let state = body(2.0, SpecificInertia::ZERO, Vec3::ZERO);
    let mut input = forcing(central());
    input.gravity_world_m_s2 = Vec3::new(0.1, -0.3, 0.2);
    let force = Vec3::new(0.3, -0.2, 0.1);
    let loads = PreparedSpatialLoads::new(
        crate::Wrench {
            force_world: force,
            torque_body: Vec3::ZERO,
        },
        Vec3::ZERO,
        &[],
        || false,
    )
    .unwrap();
    let result = state
        .advance_nozzles_transport(input, &loads, Vec3::ZERO, 3.7, controls(0.2), || false)
        .unwrap();
    let native = state
        .advance_spatial(
            input.nozzles.forcing(input.gravity_world_m_s2),
            &loads,
            3.7,
            controls(0.2),
            || false,
        )
        .unwrap();
    assert_eq!(result.mechanics, native);
    let exact = state
        .translation()
        .advance_with_exhaust_forced(
            ExhaustForcing {
                mass_flow_kg_s: input.nozzles.mass_flow_kg_s(),
                relative_momentum_rate_n: input.nozzles.reaction_force_body_n().scale(-1.0),
                relative_energy_rate_w: input.local_relative_energy_rate_w,
            },
            ExhaustEnvironment {
                gravity_m_s2: input.gravity_world_m_s2,
                external_force_n: force,
            },
            3.7,
            || false,
        )
        .unwrap();
    vector(
        result.exhaust_linear_momentum_world_kg_m_s,
        exact.exhaust_momentum_world_kg_m_s,
        2e-13,
    );
    near(
        result.exhaust_energy_world_j,
        exact.exhaust_energy_world_j,
        4e-13,
    );
    near(
        result.translational_transport_energy_j,
        exact.relative_kinetic_transport_energy_j,
        4e-13,
    );
    vector(
        result.external_linear_impulse_world_n_s,
        exact.external_impulse_n_s,
        2e-13,
    );
    vector(
        result.linear_momentum_residual_world_kg_m_s,
        Vec3::ZERO,
        2e-13,
    );
    near(result.co_rotating_transport_energy_j, 0.0, 0.0);
    near(result.mechanics.burned_propellant_mass_kg, 2.0, 2e-15);
    assert!(result.mechanics.coast_duration_s > 0.0);
}

#[derive(Clone, Copy, Default)]
struct Reference {
    x: Vec3,
    v: Vec3,
    p: Vec3,
    l: Vec3,
    e: f64,
}
impl Reference {
    fn plus(self, other: Self, scale: f64) -> Self {
        Self {
            x: self.x.add(other.x.scale(scale)),
            v: self.v.add(other.v.scale(scale)),
            p: self.p.add(other.p.scale(scale)),
            l: self.l.add(other.l.scale(scale)),
            e: self.e + scale * other.e,
        }
    }
}

// Independent continuous equations: all ports share c=(-3,0,0), their
// first moments/couples cancel, K=diag(1/2,1/2,1)=K_propellant. Hence omega_z
// stays exactly 1.3 while inertia decreases. Rotation is analytic sin/cos;
// RK4 integrates COM motion and the physical particle definitions directly.
fn reference_rate(t: f64, value: Reference, local_energy: f64) -> Reference {
    let (sin, cos) = (1.3 * t).sin_cos();
    let relative = Vec3::new(-3.0 * cos, -3.0 * sin, 0.0);
    let p = value.v.add(relative).scale(0.7);
    Reference {
        x: value.v,
        v: relative.scale(-0.7 / (5.0 - 0.7 * t)),
        p,
        l: value.x.cross(p).add(Vec3::new(0.0, 0.0, 0.7 * 1.3)),
        e: local_energy
            + 0.35 * value.v.norm_squared()
            + 0.7 * value.v.dot(relative)
            + 0.35 * 1.3_f64.powi(2),
    }
}
fn reference(state: &TensorRocketState, local_energy: f64) -> Reference {
    let mut value = Reference {
        x: state.translation().position_m,
        v: state.translation().velocity_m_s,
        ..Reference::default()
    };
    let h = 1.0 / 8192.0;
    for i in 0..8192 {
        let t = f64::from(i) * h;
        let k1 = reference_rate(t, value, local_energy);
        let k2 = reference_rate(t + 0.5 * h, value.plus(k1, 0.5 * h), local_energy);
        let k3 = reference_rate(t + 0.5 * h, value.plus(k2, 0.5 * h), local_energy);
        let k4 = reference_rate(t + h, value.plus(k3, h), local_energy);
        value = value
            .plus(k1, h / 6.0)
            .plus(k2, h / 3.0)
            .plus(k3, h / 3.0)
            .plus(k4, h / 6.0);
    }
    value
}
fn rotating() -> (TensorRocketState, NozzleTransportForcing) {
    let c = Vec3::new(-3.0, 0.0, 0.0);
    let nozzles = PreparedNozzleArray::new(
        &[
            BalancedNozzleGroup::OpposedPair {
                mass_flow_per_nozzle_kg_s: 0.175,
                arm_body_m: Vec3::new(1.0, 0.0, 0.0),
                positive_exhaust_velocity_body_m_s: c,
                negative_exhaust_velocity_body_m_s: c,
            },
            BalancedNozzleGroup::OpposedPair {
                mass_flow_per_nozzle_kg_s: 0.175,
                arm_body_m: Vec3::new(0.0, 1.0, 0.0),
                positive_exhaust_velocity_body_m_s: c,
                negative_exhaust_velocity_body_m_s: c,
            },
        ],
        || false,
    )
    .unwrap();
    (
        body(
            2.0,
            nozzles.co_rotating_specific_inertia_m2(),
            Vec3::new(0.0, 0.0, 1.3),
        ),
        forcing(nozzles),
    )
}

#[test]
fn rotating_mass_loss_flux_refines_against_independent_continuous_particle_equations() {
    let (state, input) = rotating();
    let exact = reference(&state, input.local_relative_energy_rate_w);
    let run = |h| {
        state
            .advance_nozzles_transport(
                input,
                &PreparedSpatialLoads::ZERO,
                Vec3::ZERO,
                1.0,
                controls(h),
                || false,
            )
            .unwrap()
    };
    let coarse = run(0.1);
    let fine = run(0.05);
    let p_error = |s: NozzleTransportStep| {
        s.exhaust_linear_momentum_world_kg_m_s
            .sub(exact.p)
            .norm_squared()
            .sqrt()
    };
    let l_error = |s: NozzleTransportStep| {
        s.exhaust_angular_momentum_world_kg_m2_s
            .sub(exact.l)
            .norm_squared()
            .sqrt()
    };
    let e_error = |s: NozzleTransportStep| (s.exhaust_energy_world_j - exact.e).abs();
    eprintln!(
        "P/L/E transport errors h=.1: {:.9e} {:.9e} {:.9e}; h=.05: {:.9e} {:.9e} {:.9e}",
        p_error(coarse),
        l_error(coarse),
        e_error(coarse),
        p_error(fine),
        l_error(fine),
        e_error(fine)
    );
    assert!(p_error(fine) < p_error(coarse) / 3.5);
    assert!(l_error(fine) < l_error(coarse) / 3.5);
    assert!(e_error(fine) < e_error(coarse) / 3.5);
    assert!(l_error(fine) < 0.003);
    vector(
        fine.linear_momentum_residual_world_kg_m_s,
        Vec3::ZERO,
        1e-13,
    );
    assert!(
        fine.angular_momentum_residual_world_kg_m2_s
            .norm_squared()
            .sqrt()
            < 2e-8
    );
    near(
        fine.co_rotating_transport_energy_j,
        0.35 * 1.3_f64.powi(2),
        1e-13,
    );
}

#[test]
fn physical_origin_shift_and_galilean_energy_transport_keep_separate_datums() {
    let (state, mut input) = rotating();
    input.local_relative_energy_rate_w = -1e18;
    let origin = Vec3::new(-1.0, 2.0, 0.5);
    let run = |state: TensorRocketState, origin| {
        state
            .advance_nozzles_transport(
                input,
                &PreparedSpatialLoads::ZERO,
                origin,
                0.8,
                controls(0.04),
                || false,
            )
            .unwrap()
    };
    let base = run(state, Vec3::ZERO);
    let shifted = run(state, origin);
    vector(
        shifted.exhaust_angular_momentum_world_kg_m2_s,
        base.exhaust_angular_momentum_world_kg_m2_s
            .sub(origin.cross(base.exhaust_linear_momentum_world_kg_m_s)),
        2e-14,
    );
    let boost = Vec3::new(100_000.0, -200_000.0, 300_000.0);
    let boosted = TensorRocketState::new(
        RocketState {
            velocity_m_s: state.translation().velocity_m_s.add(boost),
            ..*state.translation()
        },
        state.orientation(),
        state.angular_momentum_body_kg_m2_s(),
        state.inertia(),
    )
    .unwrap();
    let result = run(boosted, Vec3::ZERO);
    let dm = base.mechanics.burned_propellant_mass_kg;
    vector(
        result.exhaust_linear_momentum_world_kg_m_s,
        base.exhaust_linear_momentum_world_kg_m_s
            .add(boost.scale(dm)),
        2e-10,
    );
    near(
        result.translational_transport_energy_j,
        base.translational_transport_energy_j
            + boost.dot(base.exhaust_linear_momentum_world_kg_m_s)
            + 0.5 * dm * boost.norm_squared(),
        1e-4,
    );
    near(
        result.co_rotating_transport_energy_j,
        base.co_rotating_transport_energy_j,
        0.0,
    );
    assert!(base.translational_transport_energy_j.abs() > 0.01);
    assert!(base.co_rotating_transport_energy_j > 0.0);
}

#[test]
fn closed_exhaust_and_zero_duration_preserve_native_motion_and_emit_nothing() {
    let (loaded, input) = rotating();
    let dry = TensorRocketState::new(
        RocketState {
            propellant_mass_kg: 0.0,
            ..*loaded.translation()
        },
        loaded.orientation(),
        loaded.angular_momentum_body_kg_m2_s(),
        loaded.inertia(),
    )
    .unwrap();
    for (state, duration) in [(loaded, 0.0), (dry, 1.0)] {
        let result = state
            .advance_nozzles_transport(
                input,
                &PreparedSpatialLoads::ZERO,
                Vec3::ZERO,
                duration,
                controls(0.1),
                || false,
            )
            .unwrap();
        assert_eq!(
            result.mechanics,
            state
                .advance_nozzles(&input.nozzles, Vec3::ZERO, duration, controls(0.1), || {
                    false
                })
                .unwrap()
        );
        assert_eq!(result.exhaust_linear_momentum_world_kg_m_s, Vec3::ZERO);
        assert_eq!(result.exhaust_angular_momentum_world_kg_m2_s, Vec3::ZERO);
        near(result.exhaust_energy_world_j, 0.0, 0.0);
    }
}

#[test]
fn invalid_energy_work_limit_and_late_cancellation_publish_no_partial_exchange() {
    let (state, input) = rotating();
    let run = |input, controls, cancelled: fn() -> bool| {
        state.advance_nozzles_transport(
            input,
            &PreparedSpatialLoads::ZERO,
            Vec3::ZERO,
            1.0,
            controls,
            cancelled,
        )
    };
    let mut invalid = input;
    invalid.local_relative_energy_rate_w = f64::NAN;
    assert!(run(invalid, controls(0.1), || false).is_err());
    assert!(matches!(
        run(
            input,
            RigidRocketSolve {
                maximum_substeps: 2,
                ..controls(0.1)
            },
            || false
        ),
        Err(RigidRocketError::SubstepBudget {
            completed: 2,
            maximum: 2
        })
    ));
    let mut polls = 0;
    let result = state.advance_nozzles_transport(
        input,
        &PreparedSpatialLoads::ZERO,
        Vec3::ZERO,
        1.0,
        controls(0.1),
        || {
            polls += 1;
            polls == 35
        },
    );
    assert!(matches!(
        result,
        Err(RigidRocketError::Cancelled)
            | Err(RigidRocketError::Translation(
                crate::rocket::RocketError::Cancelled
            ))
    ));
    assert_eq!(polls, 35);
}
