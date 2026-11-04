use super::*;
use crate::bernoulli_aperture::BernoulliAperture;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicApertureSpec};
use crate::bernoulli_aperture::paired::BladePropagator;
use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use crate::dg_port::DgPortLoad;
use fs_dcontact::Obstacle;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_exec::CancelGate;
use fs_material::gas_species::input::{GasRecord, read_gas_record};
use fs_scenario::gesture::{GestureTarget, GestureTrack, GestureValue};

const RATE: u32 = 96_000;
const SAMPLES: usize = 512;

fn gas(text: &str) -> GasRecord {
    read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap()
}

fn bore(outlet: DgBoundary) -> DgDuct {
    let gas = gas("mixture 293.15 101325 mole 2 N2 0.8 O2 0.2");
    DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }] },
        gas.state(), DgMesh { elements: 1, order: 2 },
        DgBoundary::Port(PortId(0)), outlet,
    ).unwrap().with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
}

fn computed_load() -> DgPortLoad {
    let exterior = gas("mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2");
    let source = ComputedBaffledInput::parse(
        "20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257",
    ).unwrap();
    DgPortLoad::new_with_computed_baffled(
        bore(DgBoundary::Port(PortId(1))), RATE, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, &exterior, source, &CancelGate::new(),
    ).unwrap()
}

fn model(load: DgPortLoad) -> ApertureDg {
    let rest = 0.0002;
    let stiffness = 500.0;
    let lay = Obstacle::new(
        vec![-1.0], 1, 1, vec![0.0], vec![1.0], 0.0, 2.0,
        "DG observer test: explicit zero contact".into(),
    ).unwrap();
    let aperture = DynamicAperture::new_with_propagator(
        DynamicApertureSpec {
            aperture: BernoulliAperture {
                rest_opening_m: rest,
                width_m: 0.001,
                closing_pressure_pa: stiffness * rest / 1e-5,
            },
            mass_kg: 1e-5,
            stiffness_n_m: stiffness,
            damping_ratio: 0.02,
            density_kg_m3: load.duct().inlet_gas().density,
            impedance_pa_s_m3: load.impedance_pa_s_m3(),
            time_step_s: load.propagator().dt(),
            max_steps: SAMPLES as u64,
        },
        ApertureState { opening_m: rest, opening_velocity_m_s: 0.0 },
        lay, BladePropagator::ExactZoh,
    ).unwrap();
    ApertureDg::new(aperture, load).unwrap()
}

fn schedule() -> GestureSchedule {
    GestureSchedule::try_new(1000, vec![GestureTrack {
        id: "blow".into(),
        target: GestureTarget::BlowingPressure,
        initial: GestureValue::PressurePa(30.0),
        events: vec![],
    }]).unwrap()
}

fn config() -> AperturePerformanceConfig {
    AperturePerformanceConfig {
        sample_rate_hz: RATE, samples: SAMPLES as u64, max_block: 113,
        max_compile_work: 1 << 20, max_controls: 1024,
    }
}

fn location() -> CircularOutletReceiver {
    CircularOutletReceiver {
        position_m: [0.0, 0.0, 0.1],
        radial_rings: 2, angular_points: 16, maximum_frequency_hz: 2000.0,
    }
}

fn performance(load: DgPortLoad, observation: ApertureObservation) -> AperturePerformance {
    AperturePerformance::new(
        CoupledAperture::Dg(model(load)), observation, schedule(), config(),
    ).unwrap()
}

fn dg(performance: &AperturePerformance) -> &ApertureDg {
    let CoupledAperture::Dg(model) = performance.system() else {
        panic!("fixture retains its DG owner");
    };
    model
}

fn assert_same_physics(left: &ApertureDg, right: &ApertureDg) {
    assert_eq!(left.aperture().state(), right.aperture().state());
    assert_eq!(left.aperture().accepted_steps(), right.aperture().accepted_steps());
    assert_eq!(left.duct_state(), right.duct_state());
    assert_eq!(
        left.load().terminal().unwrap().pole_states(),
        right.load().terminal().unwrap().pole_states(),
    );
}

#[test]
fn dg_baffled_receiver_uses_accepted_outward_flow_and_retains_original_physics() {
    // One cold computed source is shared by clones; observing never refits it.
    let load = computed_load();
    let record = load.radiation_record().unwrap();
    let exterior = *record.gas_record().state();
    let medium = RayleighMedium::from(&exterior);
    assert_ne!(exterior.density.to_bits(), load.duct().inlet_gas().density.to_bits());
    assert_ne!(exterior.sound_speed.to_bits(), load.duct().inlet_gas().sound_speed.to_bits());
    let receiver = location();
    let observation = ApertureObservation::TubeBaffled(receiver);
    let mut independent = BaffledPressure::circular_openings_multirate_differentiating(
        record.apertures(), RATE, RATE, receiver, medium,
    ).unwrap();
    let mut direct = model(load.clone());
    let mut inlet = performance(load.clone(), ApertureObservation::Inlet);
    let mut observed = performance(load.clone(), observation);
    let mut partitioned = performance(load.clone(), observation);
    let mut explicit = AperturePerformance::new_with_exterior_gas(
        CoupledAperture::Dg(model(load.clone())), observation, schedule(), config(), exterior,
    ).unwrap();

    assert_eq!(observed.exterior_gas(), Some(&exterior));
    assert_eq!(observed.baffled_medium(), Some(medium));
    assert!(matches!(observed.observation_policy(), ObservationPolicy::Diagnostic));
    assert!(inlet.baffled_receiver().is_none());
    assert!(inlet.baffled_medium().is_none());
    assert!(inlet.exterior_gas().is_none());
    assert_eq!(observed.baffled_receiver().unwrap().delay_samples, independent.delay_samples);
    assert_eq!(observed.baffled_receiver().unwrap().filter_delay_samples(),
        independent.filter_delay_samples());
    assert!(independent.filter_delay_samples() > 0.0);

    // These refusals cannot alter the original source or admit an unrelated
    // observer medium/band for its already installed radiation reaction.
    assert!(AperturePerformance::new_with_exterior_gas(
        CoupledAperture::Dg(model(load.clone())), observation, schedule(), config(),
        *load.duct().inlet_gas(),
    ).is_err());
    let excessive_band = CircularOutletReceiver {
        maximum_frequency_hz: 2500.0, ..receiver
    };
    assert!(AperturePerformance::new(
        CoupledAperture::Dg(model(load.clone())),
        ApertureObservation::TubeBaffled(excessive_band), schedule(), config(),
    ).is_err());

    let mut expected = vec![0.0; SAMPLES];
    let mut distinct_flow = false;
    for slot in &mut expected {
        let frame = direct.step(TubeDrive {
            upstream_pressure_pa: 30.0, body_flow_m3_s: 0.0,
        }).unwrap();
        let terminal = frame.terminal.as_ref().unwrap();
        *slot = independent.step(&[terminal.outward_flow_m3_s]).unwrap();
        distinct_flow |= terminal.outward_flow_m3_s.to_bits()
            != frame.aperture.bore_flow_m3_s.to_bits();
        let (mut internal, mut external) = ([0.0], [0.0]);
        inlet.block(&mut internal).unwrap();
        observed.block(&mut external).unwrap();
        assert_eq!(internal[0].to_bits(), frame.aperture.bore_pressure_pa.to_bits());
        assert_eq!(external[0].to_bits(), slot.to_bits());
        assert_same_physics(&direct, dg(&inlet));
        assert_same_physics(&direct, dg(&observed));
    }
    assert!(distinct_flow, "receiver must observe propagated terminal flow");
    assert!(expected.iter().any(|sample| sample.abs() > 1e-14));

    let mut split = vec![0.0; SAMPLES];
    let mut matched = vec![0.0; SAMPLES];
    for block in split.chunks_mut(37) { partitioned.block(block).unwrap(); }
    for block in matched.chunks_mut(113) { explicit.block(block).unwrap(); }
    assert_eq!(split, expected);
    assert_eq!(matched, expected);
    assert_same_physics(&direct, dg(&partitioned));
    assert_same_physics(&direct, dg(&explicit));
    assert_eq!(inlet.applied_controls(), observed.applied_controls());
    assert_eq!(observed.applied_controls(), partitioned.applied_controls());
}

#[test]
fn dg_exterior_observation_refuses_an_unbound_ideal_terminal() {
    let load = DgPortLoad::new(
        bore(DgBoundary::IdealOpen), 1.0 / f64::from(RATE),
        1e-11, MAX_DG_PROPAGATOR_DOFS,
    ).unwrap();
    assert!(load.radiation_record().is_none());
    assert!(AperturePerformance::new(
        CoupledAperture::Dg(model(load.clone())),
        ApertureObservation::TubeBaffled(location()), schedule(), config(),
    ).is_err());
    let mut inlet = performance(load, ApertureObservation::Inlet);
    let mut output = [0.0];
    inlet.block(&mut output).unwrap();
    assert!(output[0].is_finite());
    assert!(inlet.baffled_receiver().is_none());
}
