//! Both original acoustic assemblies remain live while one physical receiver is selected.
use super::*;
use super::radiation::{ACOUSTIC_RATE, computed_primary, dx, exterior, placement};
use super::upstream::coupled;
use crate::nsf_dg::performance::ReservoirObservation;
use crate::nsf_dg::thermoviscous::performance::NsfThermoviscousDgObservationPort;
use crate::pcm_wav::baffled::CircularOutletReceiver;
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationLimits};
use crate::render::RenderError;

use NsfThermoviscousDgObservationPort::{Downstream, Upstream};

fn controls(count: usize) -> Vec<NsfDgDrive> {
    (0..count).map(|step| NsfDgDrive::FiniteReservoirActuation {
        aperture: 0, volume_rate_m3_s: if step < count / 2 { -1e-8 } else { 0.0 },
        heat_power_w: 0.0, maximum_relative_volume_change: 0.01,
    }).collect()
}

fn same_owners(a: &NsfThermoviscousDg<'_, '_, D3q27>, b: &NsfThermoviscousDg<'_, '_, D3q27>) {
    assert_eq!(a.accepted_steps(), b.accepted_steps());
    assert_eq!(a.zone().state_hash(), b.zone().state_hash());
    assert_eq!(a.zone().cumulative_ledger(), b.zone().cumulative_ledger());
    assert_eq!(a.zone().last_reservoir_flows(), b.zone().last_reservoir_flows());
    assert_eq!(a.finite_reservoir_inventory(), b.finite_reservoir_inventory());
    for (a, b) in [(a.load(), b.load()), (a.upstream_load().unwrap(), b.upstream_load().unwrap())] {
        let a = a.system(); let b = b.system();
        assert_eq!(a.accepted_steps(), b.accepted_steps());
        assert_eq!(a.state(), b.state());
        assert_eq!(a.load().state(), b.load().state());
        assert_eq!(a.load_at_port(0).unwrap().state(), b.load_at_port(0).unwrap().state());
        assert_eq!(a.load().terminal().unwrap().pole_states(), b.load().terminal().unwrap().pole_states());
        assert_eq!(a.load().last_independent_mouths_frame(), b.load().last_independent_mouths_frame());
    }
}

#[test]
fn selected_upstream_and_downstream_ports_replay_their_original_field_pressures() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let downstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let upstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&downstream, &exterior(), &gate);
    let upstream_primary = computed_primary(&upstream, &exterior(), &gate);
    let make = || coupled(&downstream, &upstream, medium.clone(), &primary,
        &upstream_primary, ThermoviscousDgLimits::default(), &gate);
    let mut original = NsfThermoviscousDgPerformance::new_with_drives(make(), 1,
        ApertureObservation::Inlet, ACOUSTIC_RATE, controls(8), 8, Arc::clone(&gate)).unwrap();
    let render = |target| NsfThermoviscousDgPerformance::new_with_drives_at(make(), target,
        ApertureObservation::Inlet, ACOUSTIC_RATE, controls(8), 8, Arc::clone(&gate)).unwrap();
    let mut explicit = render(Downstream(1));
    let mut upstream_primary = render(Upstream(1));
    let mut upstream_cork = render(Upstream(0));
    assert_eq!(original.observation_target(), Downstream(1));
    assert_eq!(upstream_primary.observation_target(), Upstream(1));
    assert_eq!(upstream_cork.observation_port(), 0);
    let mut expected = [0.0; 8];
    let mut distinct = false;
    for value in &mut expected {
        let mut actual = [0.0];
        original.block(&mut actual).unwrap();
        *value = actual[0];
        let downstream = original.system().load().system().last_frame().unwrap();
        let upstream = original.system().upstream_load().unwrap().system().last_frame().unwrap();
        assert_eq!(actual[0].to_bits(), downstream.mean_pressure_pa.to_bits());
        upstream_primary.block(&mut actual).unwrap();
        assert_eq!(actual[0].to_bits(), upstream.mean_pressure_pa.to_bits());
        upstream_cork.block(&mut actual).unwrap();
        assert_eq!(upstream.additional_ducts[0].port, 0);
        assert_eq!(actual[0].to_bits(), upstream.additional_ducts[0].mean_pressure_pa.to_bits());
        distinct |= upstream.mean_pressure_pa != downstream.mean_pressure_pa;
    }
    let mut actual = [0.0; 8];
    for chunk in actual.chunks_mut(3) { explicit.block(chunk).unwrap(); }
    assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    assert!(distinct, "the upstream target must select its own acoustic history");
    for renderer in [&explicit, &upstream_primary, &upstream_cork] {
        same_owners(original.system(), renderer.system());
        assert_eq!(renderer.samples_rendered(), 8);
    }
}

#[test]
fn upstream_material_microphone_replays_original_flows_on_an_independent_output_clock() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let downstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let upstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&downstream, &exterior(), &gate);
    let upstream_primary = computed_primary(&upstream, &exterior(), &gate);
    let make = || coupled(&downstream, &upstream, medium.clone(), &primary,
        &upstream_primary, ThermoviscousDgLimits::default(), &gate);
    let propagation = GasPropagation::classical(GasPropagationLimits::surface_reference());
    let observation = ApertureObservation::TubeBaffled(CircularOutletReceiver {
        position_m: [0.01, 0.0, 0.02], radial_rings: 4, angular_points: 16,
        maximum_frequency_hz: 1000.0,
    });
    // Upstream has its original circular exterior without downstream mouth
    // quadrature. Selecting the wrong load would refuse this distinct surface.
    let render = |block| NsfThermoviscousDgPerformance::new_resampled_with_drives_at(
        make(), Upstream(1), observation, ACOUSTIC_RATE * 3 / 2,
        controls(64), 96, block, Arc::clone(&gate)).unwrap()
        .with_observer_gas(&propagation).unwrap();
    let mut whole = render(96);
    let mut parts = render(7);
    let mut direct = make();
    let mut microphone = whole.baffled_receiver().unwrap().clone();
    assert_eq!(whole.observer_gas().unwrap().gas,
        *upstream_primary.radiation_record().unwrap().gas_record().state());
    assert_eq!(whole.reservoir_observation_count(), 0);
    let mut expected = Vec::new();
    for pair in whole.drive_history().chunks_exact(2) {
        let mut physical = Vec::new();
        for &drive in pair {
            direct.step_with_drive(drive, &gate).unwrap();
            physical.push(direct.upstream_load().unwrap().system().load().last_terminal_frames()
                .iter().map(|terminal| terminal.outward_flow_m3_s).collect::<Vec<_>>());
        }
        // Two original 8 kHz means cover three 12 kHz observer intervals.
        // The middle output straddles both source intervals equally.
        for weights in [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]] {
            let flows: Vec<_> = physical[0].iter().zip(&physical[1])
                .map(|(a, b)| weights[0] * a + weights[1] * b).collect();
            expected.push(microphone.step_trace(&flows).unwrap());
        }
    }
    let mut a = [0.0; 96]; let mut b = [0.0; 96];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(7) { parts.block(chunk).unwrap(); }
    assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    let scale = a.iter().chain(&expected).map(|x| x.abs()).fold(0.0_f64, f64::max);
    assert!(scale > 0.0, "the original upstream field must reach its microphone");
    for (actual, expected) in a.into_iter().zip(expected) {
        assert!((actual - expected).abs() <= 1024.0 * f64::EPSILON * scale);
    }
    same_owners(whole.system(), parts.system());
    same_owners(whole.system(), &direct);
    assert_eq!(whole.system().accepted_steps(), 64);
    assert_eq!(whole.system().upstream_load().unwrap().system().accepted_steps(), 64);
    assert_eq!(whole.samples_rendered(), 96);
}

#[test]
fn observation_admits_both_horizons_and_keeps_escape_geometry_downstream() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let downstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let upstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&downstream, &exterior(), &gate);
    let upstream_primary = computed_primary(&upstream, &exterior(), &gate);
    let make = |horizon| coupled(&downstream, &upstream, medium.clone(), &primary,
        &upstream_primary, ThermoviscousDgLimits {
            maximum_steps: horizon, ..ThermoviscousDgLimits::default()
        }, &gate);
    // An unobserved upstream field still participates in every physical step.
    for target in [Downstream(1), Upstream(1)] {
        assert!(NsfThermoviscousDgPerformance::new_with_drives_at(make(1), target,
            ApertureObservation::Inlet, ACOUSTIC_RATE, controls(2), 2, Arc::clone(&gate)).is_err());
    }
    let unattached = radiation::attach(radiation::system(&downstream, medium.clone(),
        &primary, 256, &gate), &gate);
    assert!(NsfThermoviscousDgPerformance::new_with_drives_at(unattached, Upstream(1),
        ApertureObservation::Inlet, ACOUSTIC_RATE, controls(1), 1, Arc::clone(&gate)).is_err());
    assert!(NsfThermoviscousDgPerformance::new_with_drives_at(make(2), Upstream(2),
        ApertureObservation::Inlet, ACOUSTIC_RATE, controls(1), 1, Arc::clone(&gate)).is_err());
    let renderer = NsfThermoviscousDgPerformance::new_with_drives_at(make(1), Upstream(1),
        ApertureObservation::Inlet, ACOUSTIC_RATE, controls(1), 1, Arc::clone(&gate)).unwrap();
    assert_eq!(renderer.system().accepted_steps(), 0);
    assert!(matches!(renderer.with_reservoir_observation_at(1, ReservoirObservation {
        origin_m: placement().origin_m, angle_rad: 0.0,
    }), Err(RenderError::Sizing { what }) if what == "resolved escape observation requires the downstream exterior"));
}
