//! Causal support memory, exact interval composition and separate wave/heat work.

use super::{close, solid};
use fs_blake3::ContentHash;
use fs_couple::elastic_reflection::{
    OrientedOscillatorSupport, PreparedTransientElasticReflection, TransientElasticReflectionError,
    TransientElasticReflectionLimits,
};
use fs_ga::{So3, So3Tangent, Vec3};
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_phs::WallPin;

fn aligned(resistance: f64) -> OrientedOscillatorSupport {
    OrientedOscillatorSupport::new(
        [WallPin {
            surface_density: 1e5,
            stiffness_per_area: 1e11,
            resistance,
        }; 3],
        So3::identity(),
    )
}

fn oriented() -> OrientedOscillatorSupport {
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
    OrientedOscillatorSupport::new(
        [
            WallPin {
                surface_density: 2e4,
                stiffness_per_area: 2e10,
                resistance: 1.5e6,
            },
            WallPin {
                surface_density: 3e4,
                stiffness_per_area: 6e10,
                resistance: 2e6,
            },
            WallPin {
                surface_density: 1e5,
                stiffness_per_area: 1e11,
                resistance: 3e6,
            },
        ],
        rotation,
    )
}

#[test]
fn g1_held_pulse_matches_exact_original_oscillator_and_trajectory_work() {
    let material = solid(1.0);
    let dt = 0.001;
    let owner = PreparedTransientElasticReflection::prepare(&material, aligned(3e6), dt).unwrap();
    let mut state = owner.state(0.0, [0.0; 3], [0.0; 3]).unwrap();
    let mut workspace = owner.workspace();
    let a = 0.001;
    let frame = *owner
        .step_into(&mut state, [a, 0.0, 0.0], &mut workspace, &mut || false)
        .unwrap();
    // m*x''+(r+rho*cp)*x'+K*x=2*rho*cp*a, starting from rest.
    let mass: f64 = 1e5;
    let stiffness: f64 = 1e11;
    let impedance: f64 = 1e7;
    let resistance: f64 = 3e6;
    let force = 2.0 * impedance * a;
    let gamma = (resistance + impedance) / (2.0 * mass);
    let wd = (stiffness / mass - gamma * gamma).sqrt();
    let decay = (-gamma * dt).exp();
    let x = force / stiffness * (1.0 - decay * ((wd * dt).cos() + gamma / wd * (wd * dt).sin()));
    let velocity_scale = force / (mass * wd);
    let v = velocity_scale * decay * (wd * dt).sin();
    let alpha = 2.0 * gamma;
    let beta = 2.0 * wd;
    let e = (-alpha * dt).exp();
    let integral_one = -(-alpha * dt).exp_m1() / alpha;
    let integral_cos = (alpha * (1.0 - e * (beta * dt).cos()) + beta * e * (beta * dt).sin())
        / (alpha * alpha + beta * beta);
    let integral_velocity_squared =
        0.5 * velocity_scale * velocity_scale * (integral_one - integral_cos);
    close(state.principal_displacement_m()[2], x, 5e-17);
    close(state.principal_velocity_m_s()[2], v, 5e-14);
    assert_eq!(state.principal_displacement_m()[..2], [0.0, 0.0]);
    assert_eq!(state.principal_velocity_m_s()[..2], [0.0, 0.0]);
    close(frame.initial_reflected_velocity_m_s[0], a, 1e-16);
    close(frame.final_reflected_velocity_m_s[0], a - v, 5e-14);
    close(frame.mean_reflected_velocity_m_s[0], a - x / dt, 5e-14);
    close(
        frame.support_traction_pa[2],
        impedance * (v - 2.0 * a),
        5e-7,
    );
    // Real instantaneous incident work is h*Z*a², with no phasor factor 1/2.
    let incoming = dt * impedance * a * a;
    let outgoing = incoming - 2.0 * impedance * a * x + impedance * integral_velocity_squared;
    close(frame.work.incoming_wave_energy_j_m2, incoming, 3e-15);
    close(frame.work.outgoing_wave_energy_j_m2, outgoing, 3e-12);
    close(
        frame.work.support_dissipation_j_m2,
        resistance * integral_velocity_squared,
        3e-12,
    );
    close(
        frame.work.final_energy.total_j_m2,
        0.5 * (mass * v * v + stiffness * x * x),
        3e-12,
    );
    assert!(frame.work.balance_residual_j_m2().abs() <= frame.work.balance_tolerance_j_m2);
    // After the incident pulse ends, stored mechanical motion still emits waves.
    let tail = owner
        .step_into(&mut state, [0.0; 3], &mut workspace, &mut || false)
        .unwrap();
    assert_eq!(tail.work.incoming_wave_energy_j_m2, 0.0);
    assert!(tail.work.outgoing_wave_energy_j_m2 > 0.0);
    assert!(tail.work.support_dissipation_j_m2 > 0.0);
    assert!(tail.final_reflected_velocity_m_s[0].abs() > 1e-6);
}

#[test]
fn g1_rotated_causal_response_composes_exactly_across_held_intervals() {
    let material = solid(1.0);
    let support = oriented();
    let whole = PreparedTransientElasticReflection::prepare(&material, support, 0.0004).unwrap();
    let half = PreparedTransientElasticReflection::prepare(&material, support, 0.0002).unwrap();
    let mut one = whole
        .state(0.2, [1e-8, -2e-8, 0.5e-8], [3e-5, -2e-5, 1e-5])
        .unwrap();
    let mut two = one.clone();
    let mut scratch_one = whole.workspace();
    let mut scratch_two = half.workspace();
    let input = [0.001, -0.0003, 0.0002];
    let full = *whole
        .step_into(&mut one, input, &mut scratch_one, &mut || false)
        .unwrap();
    let first = *half
        .step_into(&mut two, input, &mut scratch_two, &mut || false)
        .unwrap();
    let second = *half
        .step_into(&mut two, input, &mut scratch_two, &mut || false)
        .unwrap();
    for axis in 0..3 {
        close(
            one.principal_displacement_m()[axis],
            two.principal_displacement_m()[axis],
            5e-17,
        );
        close(
            one.principal_velocity_m_s()[axis],
            two.principal_velocity_m_s()[axis],
            5e-14,
        );
        close(
            full.mean_reflected_velocity_m_s[axis],
            0.5 * (first.mean_reflected_velocity_m_s[axis]
                + second.mean_reflected_velocity_m_s[axis]),
            5e-14,
        );
    }
    close(one.time_s(), two.time_s(), 1e-15);
    close(
        full.work.outgoing_wave_energy_j_m2,
        first.work.outgoing_wave_energy_j_m2 + second.work.outgoing_wave_energy_j_m2,
        3e-12,
    );
    close(
        full.work.support_dissipation_j_m2,
        first.work.support_dissipation_j_m2 + second.work.support_dissipation_j_m2,
        3e-12,
    );
    close(
        full.work.stored_energy_change_j_m2,
        first.work.stored_energy_change_j_m2 + second.work.stored_energy_change_j_m2,
        3e-12,
    );
    // A normal P pulse also converts into both shear polarizations through the
    // nonaligned original support; there is no incident shear source here.
    let mut rest = whole.state(0.0, [0.0; 3], [0.0; 3]).unwrap();
    let frame = whole
        .step_into(&mut rest, [0.001, 0.0, 0.0], &mut scratch_one, &mut || {
            false
        })
        .unwrap();
    assert!(frame.final_reflected_velocity_m_s[1].abs() > 1e-8);
    assert!(frame.final_reflected_velocity_m_s[2].abs() > 1e-8);
}

#[test]
fn g1_lossless_pin_ringdown_radiates_its_storage_without_fictitious_heat() {
    let material = solid(1.0);
    let owner =
        PreparedTransientElasticReflection::prepare(&material, aligned(0.0), 0.001).unwrap();
    let mut state = owner.state(0.0, [0.0, 0.0, 1e-7], [0.0; 3]).unwrap();
    let initial = owner.energy(&state).unwrap().total_j_m2;
    let mut outgoing = 0.0;
    let mut workspace = owner.workspace();
    for _ in 0..30 {
        let frame = owner
            .step_into(&mut state, [0.0; 3], &mut workspace, &mut || false)
            .unwrap();
        assert_eq!(frame.work.support_dissipation_j_m2, 0.0);
        assert_eq!(frame.work.incoming_wave_energy_j_m2, 0.0);
        assert!(frame.work.outgoing_wave_energy_j_m2 >= -frame.work.balance_tolerance_j_m2);
        outgoing += frame.work.outgoing_wave_energy_j_m2;
    }
    let final_energy = owner.energy(&state).unwrap().total_j_m2;
    assert!(final_energy < 0.1 * initial);
    close(initial, final_energy + outgoing, 3e-12);
}

#[test]
fn g4_refusal_preserves_original_motion_and_source_binding() {
    let material = solid(1.0);
    let support = oriented();
    let owner = PreparedTransientElasticReflection::prepare(&material, support, 0.001).unwrap();
    let mut state = owner.state(0.3, [1e-8; 3], [2e-5; 3]).unwrap();
    let original = state.clone();
    let mut workspace = owner.workspace();
    let mut polls = 0;
    let result = owner.step_into(&mut state, [0.001, 0.0, 0.0], &mut workspace, &mut || {
        polls += 1;
        polls >= 3 // reject at publication after the native flow and work checks
    });
    assert!(matches!(
        result,
        Err(TransientElasticReflectionError::Cancelled)
    ));
    assert_eq!(state.time_s(), original.time_s());
    assert_eq!(
        state.principal_displacement_m(),
        original.principal_displacement_m()
    );
    assert_eq!(
        state.principal_velocity_m_s(),
        original.principal_velocity_m_s()
    );
    assert!(
        owner
            .step_into(
                &mut state,
                [f64::NAN, 0.0, 0.0],
                &mut workspace,
                &mut || false
            )
            .is_err()
    );
    let other_material = solid(1.0);
    let other =
        PreparedTransientElasticReflection::prepare(&other_material, support, 0.001).unwrap();
    assert!(
        other
            .step_into(&mut state, [0.0; 3], &mut workspace, &mut || false)
            .is_err()
    );
    let mut limits = TransientElasticReflectionLimits::default();
    limits.preparation.max_dimension = 9;
    assert!(
        PreparedTransientElasticReflection::prepare_controlled(
            &material,
            support,
            0.001,
            limits,
            &mut || false
        )
        .is_err()
    );
    for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(PreparedTransientElasticReflection::prepare(&material, support, dt).is_err());
    }
}

#[test]
fn g0_original_strain_domain_covers_the_whole_emission_interval() {
    let material = ElasticAcousticState::from_bulk_and_shear(
        2500.0,
        2500.0 * (4000.0_f64.powi(2) - 4.0 / 3.0 * 2000.0_f64.powi(2)),
        2500.0 * 2000.0_f64.powi(2),
        ContentHash([31; 32]),
        1e-12,
    )
    .unwrap();
    let owner =
        PreparedTransientElasticReflection::prepare(&material, aligned(0.0), 0.001).unwrap();
    // A zero endpoint velocity does not hide a spring reservoir that can emit
    // an out-of-domain pulse during the next interval.
    assert!(owner.state(0.0, [0.0, 0.0, 1e-7], [0.0; 3]).is_err());
    let mut state = owner.state(0.0, [0.0; 3], [0.0; 3]).unwrap();
    let mut workspace = owner.workspace();
    assert!(
        owner
            .step_into(&mut state, [0.001, 0.0, 0.0], &mut workspace, &mut || false)
            .is_err()
    );
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.principal_displacement_m(), [0.0; 3]);
    let zero = owner
        .step_into(&mut state, [0.0; 3], &mut workspace, &mut || false)
        .unwrap();
    assert_eq!(zero.work.incoming_wave_energy_j_m2, 0.0);
    assert_eq!(zero.work.outgoing_wave_energy_j_m2, 0.0);
    assert_eq!(zero.work.support_dissipation_j_m2, 0.0);
    assert_eq!(zero.work.balance_residual_j_m2(), 0.0);
}
