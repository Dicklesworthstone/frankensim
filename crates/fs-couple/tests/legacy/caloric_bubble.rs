//! G1/G3/G4 material-driven bubble heating, physical volume sources, causal
//! observation and whole-interval rollback. Synthetic gas data test the model,
//! not experimental cavitation or thermal-transfer accuracy.

#![cfg(feature = "caloric-cavitation")]

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits};
use fs_couple::cavitation::bubble::caloric::{
    MaterialCaloricBubble, MaterialCaloricBubbleError,
};
use fs_flux::cavitation::caloric::{
    CaloricBubbleEquation, CaloricBubbleError, CaloricBubbleLimits, CaloricBubbleNumerics,
    CaloricGasObservation,
};
use fs_flux::cavitation::{
    BubbleLimits, BubbleNumerics, BubblePressureSample, BubbleState, CavitationError,
    KellerMiksisLimits,
};
use fs_material::cavitation::CavitationState;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

const RADIUS: f64 = 100e-6;
const AMBIENT: f64 = 100_000.0;
const TEMPERATURE: f64 = 300.0;
const END_TIME: f64 = 8e-6;

fn liquid() -> CavitationState {
    CavitationState::try_new(TEMPERATURE, 1000.0, 0.001, 0.072, 2000.0, 1500.0).unwrap()
}

fn gas() -> ClosedGasParcel {
    let model = source::load(
        "synthetic-bubble 1 0.028 3.5 0.0005 1.7e-5 0.025",
        "caloric bubble acoustic regression",
    )
    .unwrap();
    let pressure = AMBIENT - liquid().saturation_pressure_pa()
        + 2.0 * liquid().surface_tension_n_m() / RADIUS;
    let point = source::point(TEMPERATURE, pressure)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedGasParcel::new(
        model,
        point,
        (4.0 * core::f64::consts::PI / 3.0) * RADIUS.powi(3),
    )
    .unwrap()
}

fn equation(compressible: bool) -> CaloricBubbleEquation {
    let bubble = BubbleLimits {
        minimum_radius_m: 0.6 * RADIUS,
        maximum_radius_m: 1.5 * RADIUS,
        maximum_wall_mach_number: 0.01,
        maximum_vapor_density_ratio: 0.01,
    };
    if compressible {
        CaloricBubbleEquation::KellerMiksis(KellerMiksisLimits {
            bubble,
            maximum_pressure_ratio: 0.01,
        })
    } else {
        CaloricBubbleEquation::RayleighPlesset(bubble)
    }
}

fn bubble(
    source: &CavitationState,
    compressible: bool,
) -> MaterialCaloricBubble<'_, CavitationState> {
    MaterialCaloricBubble::new(
        source,
        gas(),
        equation(compressible),
        BubbleState::try_new(0.0, RADIUS, 0.0).unwrap(),
        drive(0.0),
        CaloricBubbleLimits {
            maximum_gas_liquid_density_ratio: 0.01,
            maximum_vapor_pressure_fraction: 0.1,
            maximum_gas_wall_mach_number: 0.05,
        },
        GasEnergySolve {
            minimum_temperature_k: 230.0,
            maximum_temperature_k: 1200.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 128,
        },
    )
    .unwrap()
}

fn numerics() -> CaloricBubbleNumerics {
    CaloricBubbleNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-9,
            radius_absolute_tolerance_m: 1e-14,
            velocity_absolute_tolerance_m_s: 1e-9,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-13,
            maximum_step_s: 2e-8,
            maximum_attempts: 20_000,
        },
        energy_absolute_tolerance_j: 1e-16,
        energy_balance_tolerance_j: 1e-13,
    }
}

fn drive(time_s: f64) -> BubblePressureSample {
    let omega = core::f64::consts::TAU * 10_000.0;
    BubblePressureSample {
        absolute_pressure_pa: AMBIENT + 10_000.0 * (omega * time_s).sin(),
        time_derivative_pa_s: 10_000.0 * omega * (omega * time_s).cos(),
    }
}

fn no_heat(_: f64, _: &CaloricGasObservation) -> std::result::Result<f64, CaloricBubbleError> {
    Ok(0.0)
}

fn heat_to_bath(
    _: f64,
    gas: &CaloricGasObservation,
) -> std::result::Result<f64, CaloricBubbleError> {
    // Authored thermal conductance [W/K] to a prescribed 300 K reservoir.
    // The finite liquid is not heated implicitly by this external heat port.
    Ok(1e-4 * (TEMPERATURE - gas.temperature_k))
}

fn observer_limits() -> MonopoleLimits {
    MonopoleLimits {
        maximum_ka: 0.1,
        minimum_distance_radii: 10.0,
        maximum_pressure_ratio: 0.01,
    }
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn original_caloric_trajectory_and_volume_source_are_identical_through_the_adapter() {
    let liquid = liquid();
    for compressible in [false, true] {
        let mut bubble = bubble(&liquid, compressible);
        let original_gas = bubble.gas().clone();
        let mut native = bubble.model().clone();
        let expected = native
            .advance_controlled(END_TIME, &numerics(), &drive, &heat_to_bath, &mut || false)
            .unwrap();
        let actual = bubble
            .advance_compact_source(END_TIME, &numerics(), &drive, &heat_to_bath, &mut || false)
            .unwrap();
        assert_eq!(actual.bubble, expected);
        assert_eq!(bubble.state(), native.state());
        assert_eq!(bubble.gas(), native.gas());
        assert_eq!(bubble.gas().model(), original_gas.model());
        assert_eq!(bubble.gas().mass_kg(), original_gas.mass_kg());
        assert_eq!(bubble.gas().point().axes()["specimen-lot"], 42.0);
        assert!(core::ptr::eq(bubble.source(), &liquid));
        assert_eq!(bubble.material(), liquid);
        let physical = actual.bubble.observation.bubble;
        assert_eq!(actual.source.source_time_s, END_TIME);
        assert_eq!(
            actual.source.strength.volume_flow_m3_s,
            physical.volume_rate_m3_s
        );
        assert_eq!(
            actual.source.strength.volume_flow_rate_m3_s2,
            physical.volume_acceleration_m3_s2
        );
        let radius = bubble.state().radius_m();
        let velocity = bubble.state().wall_velocity_m_s();
        near(
            actual.source.strength.volume_flow_rate_m3_s2,
            4.0 * core::f64::consts::PI * radius
                * (2.0 * velocity * velocity + radius * physical.wall_acceleration_m_s2),
            1e-14,
        );
        assert_eq!(actual.source.strength.force_on_fluid_n, [0.0; 3]);
        assert_eq!(actual.source.strength.force_rate_on_fluid_n_s, [0.0; 3]);
        assert_eq!(
            actual.source,
            bubble
                .compact_source_sample(drive(END_TIME), actual.bubble.observation.heat_power_w)
                .unwrap()
        );
    }
}

#[test]
fn finite_gas_heat_changes_the_motion_and_causal_microphone_pressure() {
    let liquid = liquid();
    let mut adiabatic = bubble(&liquid, true);
    let mut cooled = bubble(&liquid, true);
    let mass = cooled.gas().mass_kg();
    let observer = cooled
        .monopole_observer(0.1, 100_000.0, observer_limits())
        .unwrap();
    let adiabatic_step = adiabatic
        .advance_observed(
            END_TIME,
            &numerics(),
            &drive,
            &no_heat,
            &observer,
            &mut || false,
        )
        .unwrap();
    let cooled_step = cooled
        .advance_observed(
            END_TIME,
            &numerics(),
            &drive,
            &heat_to_bath,
            &observer,
            &mut || false,
        )
        .unwrap();
    assert!(adiabatic.gas().state().state().temperature > TEMPERATURE);
    assert!(cooled.gas().state().state().temperature < adiabatic.gas().state().state().temperature);
    assert!(cooled_step.bubble.supplied_heat_j < 0.0);
    assert_eq!(adiabatic_step.bubble.supplied_heat_j, 0.0);
    assert_ne!(cooled.state(), adiabatic.state());
    assert_ne!(
        cooled_step.acoustic.pressure_pa,
        adiabatic_step.acoustic.pressure_pa
    );
    for step in [adiabatic_step, cooled_step] {
        // Independently integrate this fixture's Cv/R = 2.5 + 0.0005 T.
        // Use original temperature and mass, rather than a reported residual.
        let temperature = step.bubble.observation.gas.temperature_k;
        let stored = mass * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028
            * (temperature - TEMPERATURE)
            * (2.5 + 0.00025 * (temperature + TEMPERATURE));
        assert!(
            (stored - (step.bubble.supplied_heat_j - step.bubble.work_by_gas_j)).abs()
                < numerics().energy_balance_tolerance_j
        );
        near(step.acoustic.arrival_time_s, END_TIME + 0.1 / 1500.0, 1e-14);
        near(
            step.acoustic.pressure_pa,
            1000.0 * step.bubble.observation.bubble.volume_acceleration_m3_s2
                / (4.0 * core::f64::consts::PI * 0.1),
            1e-14,
        );
        assert!(step.bubble.gas_energy_residual_j.abs() < numerics().energy_balance_tolerance_j);
        assert!(step.bubble.source_energy_residual_j.abs() < numerics().energy_balance_tolerance_j);
    }
}

#[test]
fn heat_pressure_rate_is_present_in_the_compressible_acoustic_source() {
    let liquid = liquid();
    let bubble = bubble(&liquid, true);
    let saved_gas = bubble.gas().clone();
    let saved_state = bubble.state();
    let pressure = BubblePressureSample::try_new(AMBIENT, 0.0).unwrap();
    let heat = 0.002;
    let cold = bubble.compact_source_sample(pressure, 0.0).unwrap();
    let hot = bubble.compact_source_sample(pressure, heat).unwrap();
    let gas = saved_gas.state().state();
    let cv = gas.specific_heat_cp / gas.gamma;
    let pressure_rate = gas.specific_gas_constant / cv / saved_gas.volume_m3() * heat;
    let acceleration = pressure_rate
        / (liquid.liquid_density_kg_m3() * liquid.liquid_sound_speed_m_s()
            + 4.0 * liquid.dynamic_viscosity_pa_s() / RADIUS);
    near(
        hot.strength.volume_flow_rate_m3_s2 - cold.strength.volume_flow_rate_m3_s2,
        4.0 * core::f64::consts::PI * RADIUS * RADIUS * acceleration,
        2e-12,
    );
    assert_eq!(hot.strength.volume_flow_m3_s, 0.0);
    assert_eq!(bubble.state(), saved_state);
    assert_eq!(bubble.gas(), &saved_gas);
}

#[test]
fn late_observer_refusal_preserves_the_entire_gas_and_radial_history() {
    let liquid = liquid();
    let mut bubble = bubble(&liquid, true);
    let saved_state = bubble.state();
    let saved_gas = bubble.gas().clone();
    let mut limits = observer_limits();
    limits.maximum_ka = 1e-5;
    let observer = bubble.monopole_observer(0.1, 100_000.0, limits).unwrap();
    let result = bubble.advance_observed(
        END_TIME,
        &numerics(),
        &drive,
        &heat_to_bath,
        &observer,
        &mut || false,
    );
    assert!(matches!(result, Err(MaterialCaloricBubbleError::Acoustic(_))));
    assert_eq!(bubble.state(), saved_state);
    assert_eq!(bubble.gas(), &saved_gas);
    let wrong_medium = CompactMonopoleObserver::new(
        FluidAcousticState::try_new(999.0, 1500.0).unwrap(),
        0.1,
        100_000.0,
        observer_limits(),
    )
    .unwrap();
    assert_eq!(
        bubble.advance_observed(
            END_TIME,
            &numerics(),
            &drive,
            &heat_to_bath,
            &wrong_medium,
            &mut || false,
        ),
        Err(MaterialCaloricBubbleError::ObserverMediumMismatch)
    );
    assert_eq!(bubble.state(), saved_state);
    assert_eq!(bubble.gas(), &saved_gas);
    let mut native = bubble.model().clone();
    let expected = native
        .advance_controlled(END_TIME, &numerics(), &drive, &heat_to_bath, &mut || false)
        .unwrap();
    let actual = bubble
        .advance(END_TIME, &numerics(), &drive, &heat_to_bath, &mut || false)
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(bubble.gas(), native.gas());
}

#[test]
fn final_acoustic_cancellation_is_atomic_and_the_retry_replays_both_owners() {
    let liquid = liquid();
    let mut baseline = bubble(&liquid, true);
    let observer = baseline
        .monopole_observer(0.1, 100_000.0, observer_limits())
        .unwrap();
    let mut polls = 0;
    let expected = baseline
        .advance_observed(
            END_TIME,
            &numerics(),
            &drive,
            &heat_to_bath,
            &observer,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut cancelled_bubble = bubble(&liquid, true);
    let saved_state = cancelled_bubble.state();
    let saved_gas = cancelled_bubble.gas().clone();
    let mut calls = 0;
    let result = cancelled_bubble.advance_observed(
        END_TIME,
        &numerics(),
        &drive,
        &heat_to_bath,
        &observer,
        &mut || {
            calls += 1;
            calls == polls
        },
    );
    assert_eq!(
        result,
        Err(MaterialCaloricBubbleError::Dynamics(
            CaloricBubbleError::from(CavitationError::Cancelled)
        ))
    );
    assert_eq!(cancelled_bubble.state(), saved_state);
    assert_eq!(cancelled_bubble.gas(), &saved_gas);
    let retried = cancelled_bubble
        .advance_observed(
            END_TIME,
            &numerics(),
            &drive,
            &heat_to_bath,
            &observer,
            &mut || false,
        )
        .unwrap();
    assert_eq!(retried, expected);
    assert_eq!(cancelled_bubble.gas(), baseline.gas());
    assert_eq!(cancelled_bubble.state(), baseline.state());
}
