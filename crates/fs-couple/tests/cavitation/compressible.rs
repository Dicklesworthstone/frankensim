//! Original material, compressible bubble and observer composed without a
//! second equation owner. G1 analytic initial source; G4 atomic failure/retry.

use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits};
use fs_couple::cavitation::bubble::{MaterialBubble, MaterialBubbleError};
use fs_flux::cavitation::{
    BubbleGas, BubbleLimits, BubbleNumerics, BubblePressureSample, BubbleState, CavitationError,
    KellerMiksisBubble, KellerMiksisLimits,
};
use fs_material::cavitation::CavitationState;
use fs_material::fluid_acoustics::FluidAcousticState;

const RADIUS: f64 = 100e-6;
const AMBIENT: f64 = 101_325.0;

fn liquid() -> CavitationState {
    CavitationState::try_new(293.15, 998.2, 0.001, 0.072, 2339.0, 1482.0).unwrap()
}

fn drive(time_s: f64) -> BubblePressureSample {
    let omega = 2.0 * core::f64::consts::PI * 2000.0;
    BubblePressureSample {
        absolute_pressure_pa: AMBIENT + 200.0 * (omega * time_s).sin(),
        time_derivative_pa_s: 200.0 * omega * (omega * time_s).cos(),
    }
}

fn model_limits() -> KellerMiksisLimits {
    KellerMiksisLimits {
        bubble: BubbleLimits {
            minimum_radius_m: 0.5 * RADIUS,
            maximum_radius_m: 2.0 * RADIUS,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
        maximum_pressure_ratio: 0.01,
    }
}

fn numerics() -> BubbleNumerics {
    BubbleNumerics {
        relative_tolerance: 1e-9,
        radius_absolute_tolerance_m: 1e-14,
        velocity_absolute_tolerance_m_s: 1e-9,
        initial_step_s: 1e-8,
        minimum_step_s: 1e-12,
        maximum_step_s: 1e-7,
        maximum_attempts: 20_000,
    }
}

fn acoustic_limits() -> MonopoleLimits {
    MonopoleLimits {
        maximum_ka: 0.1,
        minimum_distance_radii: 10.0,
        maximum_pressure_ratio: 0.01,
    }
}

fn bubble(source: &CavitationState) -> MaterialBubble<'_, CavitationState, KellerMiksisBubble> {
    MaterialBubble::new_compressible(
        source,
        BubbleGas::equilibrium(*source, RADIUS, AMBIENT, 1.4).unwrap(),
        model_limits(),
        BubbleState::try_new(0.0, RADIUS, 0.0).unwrap(),
        drive(0.0),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-20),
        "{actual} != {expected}"
    );
}

#[test]
fn original_compressible_equation_and_material_feed_the_causal_observer() {
    let source = liquid();
    let mut coupled = bubble(&source);
    assert!(core::ptr::eq(coupled.source(), &source));
    assert_eq!(coupled.material(), source);
    let initial = coupled.observe(drive(0.0)).unwrap();
    // The exact initial equation includes the viscous pressure derivative in
    // the effective inertia, even though the instantaneous wall speed is zero.
    let acceleration = -drive(0.0).time_derivative_pa_s
        / (source.liquid_density_kg_m3() * source.liquid_sound_speed_m_s()
            + 4.0 * source.dynamic_viscosity_pa_s() / RADIUS);
    near(initial.wall_acceleration_m_s2, acceleration);

    let mut original_state = coupled.state();
    let original = coupled
        .model()
        .advance_controlled(&mut original_state, 50e-6, &numerics(), &drive, &mut || {
            false
        })
        .unwrap();
    let observer = coupled
        .monopole_observer(0.1, 100_000.0, acoustic_limits())
        .unwrap();
    let step = coupled
        .advance_observed(50e-6, &numerics(), &drive, &observer, &mut || false)
        .unwrap();
    assert_eq!(step.bubble, original);
    assert_eq!(coupled.state(), original_state);
    near(
        step.acoustic.pressure_pa,
        source.liquid_density_kg_m3() * original.observation.volume_acceleration_m3_s2
            / (4.0 * core::f64::consts::PI * 0.1),
    );
    assert_eq!(step.acoustic.source_time_s, original_state.time_s());
    near(
        step.acoustic.arrival_time_s,
        original_state.time_s() + 0.1 / source.liquid_sound_speed_m_s(),
    );
}

#[test]
fn compressible_observer_failures_leave_the_accepted_history_unchanged() {
    let source = liquid();
    let mut coupled = bubble(&source);
    let saved = coupled.state();
    let foreign = CompactMonopoleObserver::new(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        0.1,
        100_000.0,
        acoustic_limits(),
    )
    .unwrap();
    assert_eq!(
        coupled.advance_observed(50e-6, &numerics(), &drive, &foreign, &mut || false),
        Err(MaterialBubbleError::ObserverMediumMismatch)
    );
    assert_eq!(coupled.state(), saved);
    let mut tight = acoustic_limits();
    tight.maximum_ka = 1e-5;
    let observer = coupled.monopole_observer(0.1, 100_000.0, tight).unwrap();
    assert!(matches!(
        coupled.advance_observed(50e-6, &numerics(), &drive, &observer, &mut || false),
        Err(MaterialBubbleError::Acoustic(_))
    ));
    assert_eq!(coupled.state(), saved);
}

#[test]
fn compressible_final_cancellation_and_retry_preserve_the_same_physical_result() {
    let source = liquid();
    let mut baseline = bubble(&source);
    let observer = baseline
        .monopole_observer(0.1, 100_000.0, acoustic_limits())
        .unwrap();
    let mut polls = 0usize;
    let successful = baseline
        .advance_observed(50e-6, &numerics(), &drive, &observer, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut cancelled = bubble(&source);
    let saved = cancelled.state();
    let mut calls = 0usize;
    let error = cancelled.advance_observed(50e-6, &numerics(), &drive, &observer, &mut || {
        calls += 1;
        calls == polls
    });
    assert_eq!(
        error,
        Err(MaterialBubbleError::Dynamics(CavitationError::Cancelled))
    );
    assert_eq!(cancelled.state(), saved);
    let retry = cancelled
        .advance_observed(50e-6, &numerics(), &drive, &observer, &mut || false)
        .unwrap();
    assert_eq!(retry, successful);
    assert_eq!(cancelled.state(), baseline.state());
}
