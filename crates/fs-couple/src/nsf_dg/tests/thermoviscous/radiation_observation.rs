//! Original signed source paths and clocks; manufactured geometry, not timbre validation.
use super::*;
use super::radiation::{ACOUSTIC_RATE, attach, computed_primary, dx, exterior, placement, system};
use crate::nsf_dg::performance::ReservoirObservation;
use crate::pcm_wav::baffled::CircularOutletReceiver;
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationLimits};
use crate::render::RenderError;

/// 0.05 m is the nearest height the shared Rayleigh receiver admits
/// (`BaffledPressure::from_loaded`). It is also at least one output sample of
/// travel from every tile at the 8 kHz clock floor: 0.05 m * 8000 / c = 1.15
/// samples at c = 347 m/s, against the 1-sample causal minimum.
fn observation() -> ApertureObservation {
    ApertureObservation::TubeBaffled(CircularOutletReceiver {
        position_m: [0.0, 0.0, 0.05],
        radial_rings: placement().terminal_radial_rings,
        angular_points: placement().terminal_angular_points,
        maximum_frequency_hz: 1000.0,
    })
}

fn escape() -> ReservoirObservation {
    ReservoirObservation { origin_m: placement().origin_m, angle_rad: 0.0 }
}

fn controls(count: usize) -> Vec<NsfDgDrive> {
    (0..count).map(|step| NsfDgDrive::FiniteReservoirActuation {
        aperture: 0, volume_rate_m3_s: if step < count / 2 { -1e-8 } else { 0.0 },
        heat_power_w: 0.0, maximum_relative_volume_change: 0.01,
    }).collect()
}

fn playback<'a, 'c>(owner: NsfThermoviscousDg<'a, 'c, D3q27>, rate: u32,
    samples: u64, block: usize, gate: Arc<CancelGate>)
    -> NsfThermoviscousDgPerformance<'a, 'c, D3q27>
{
    let count = source_intervals_for(owner.load().duration_s(), rate, samples).unwrap();
    NsfThermoviscousDgPerformance::new_resampled_with_drives(owner, 1, observation(),
        rate, controls(count), samples, block, gate).unwrap()
}

fn observed<'a, 'c>(owner: NsfThermoviscousDg<'a, 'c, D3q27>, rate: u32,
    samples: u64, block: usize, gate: Arc<CancelGate>)
    -> NsfThermoviscousDgPerformance<'a, 'c, D3q27>
{
    // Observer order intentionally differs from reaction order [2, 1]. Both
    // use the original aperture identity, never a position in that other list.
    playback(owner, rate, samples, block, gate)
        .with_reservoir_observation_at(1, escape()).unwrap()
        .with_reservoir_observation_at(2, escape()).unwrap()
}

fn same_physical_owners(a: &NsfThermoviscousDg<'_, '_, D3q27>,
    b: &NsfThermoviscousDg<'_, '_, D3q27>)
{
    assert_eq!(a.zone().state_hash(), b.zone().state_hash());
    assert_eq!(a.zone().cumulative_ledger(), b.zone().cumulative_ledger());
    assert_eq!(a.zone().last_reservoir_flows(), b.zone().last_reservoir_flows());
    assert_eq!(a.finite_reservoir_inventory(), b.finite_reservoir_inventory());
    assert_eq!(a.accepted_steps(), b.accepted_steps());
    let a = a.load().system(); let b = b.load().system();
    assert_eq!(a.state(), b.state());
    assert_eq!(a.load().state(), b.load().state());
    assert_eq!(a.load_at_port(0).unwrap().state(), b.load_at_port(0).unwrap().state());
    assert_eq!(a.load().terminal().unwrap().pole_states(), b.load().terminal().unwrap().pole_states());
    assert_eq!(a.load().last_independent_mouths_frame(), b.load().last_independent_mouths_frame());
}

#[test]
fn complete_field_microphone_sums_original_signed_escape_and_bore_receivers() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let make = || attach(system(&field, medium.clone(), &primary, 256, &gate), &gate);
    // Cover the original receiver's 80-frame numerical delay and the actual
    // physical paths so the assertion reaches emitted source pressure.
    let mut renderer = observed(make(), ACOUSTIC_RATE, 96, 1, Arc::clone(&gate));
    let mut direct = make();
    let mut bore = renderer.baffled_receiver().unwrap().clone();
    let mut mouths: Vec<_> = [1, 2].map(|aperture|
        (aperture, renderer.reservoir_receiver_at(aperture).unwrap().clone())).to_vec();
    assert_eq!(renderer.reservoir_observation_count(), 2);
    assert!(renderer.reservoir_receiver().is_none());
    assert_eq!(renderer.reservoir_observation_at(1), Some(escape()));
    assert!(renderer.reservoir_receiver_at(0).is_none());
    let mut escape_peak = 0.0_f64;
    let mut distinct_escape_flow = false;
    for drive in renderer.drive_history().to_vec() {
        let frame = direct.step_with_drive(drive, &gate).unwrap();
        let terminal_flows: Vec<_> = direct.load().system().load().last_terminal_frames()
            .iter().map(|terminal| terminal.outward_flow_m3_s).collect();
        let mut expected = bore.step_trace(&terminal_flows).unwrap();
        let mut independent_flows = Vec::new();
        for (aperture, receiver) in &mut mouths {
            let flows: Vec<_> = frame.reservoir_flows.iter()
                .filter(|flow| flow.aperture == *aperture)
                .map(|flow| flow.outward_volume_flow_m3_s).collect();
            assert_eq!(flows.len(), receiver.input_count());
            independent_flows.push(flows[0]);
            let pressure = receiver.step_trace(&flows).unwrap();
            escape_peak = escape_peak.max(pressure.abs());
            expected += pressure;
        }
        distinct_escape_flow |= independent_flows[0] != independent_flows[1];
        let mut actual = [0.0];
        renderer.block(&mut actual).unwrap();
        assert_eq!(actual[0].to_bits(), expected.to_bits());
        assert_eq!(renderer.last_frame().unwrap().reservoir_flows, frame.reservoir_flows);
    }
    assert!(escape_peak > 0.0, "original escape sources must reach microphone pressure");
    assert!(distinct_escape_flow, "separate escape identities must retain their different flows");
    same_physical_owners(renderer.system(), &direct);
}

/// The repeating pattern of interval means between two regular clocks, from the
/// two rates alone. Returns the original intervals per group and, for each
/// output interval of the group, the fraction of that output interval that lies
/// inside each original interval. In units of 1 / (source * output) seconds,
/// original interval j spans [j * output, (j + 1) * output] and output interval
/// k spans [k * source, (k + 1) * source], so the fraction is their overlap over
/// `source`. No hand-typed table: the clock ratio is the only input.
fn interval_means(source_rate: u32, output_rate: u32) -> (usize, Vec<Vec<f64>>) {
    let (source, output) = (u64::from(source_rate), u64::from(output_rate));
    let (mut a, mut b) = (source, output);
    while b != 0 { (a, b) = (b, a % b); }
    let originals = (source / a) as usize;
    let outputs = (output / a) as usize;
    let weights = (0..outputs as u64).map(|k| (0..originals as u64).map(|j| {
        let low = (j * output).max(k * source);
        let high = ((j + 1) * output).min((k + 1) * source);
        high.saturating_sub(low) as f64 / source as f64
    }).collect()).collect();
    (originals, weights)
}

#[test]
fn projected_microphone_uses_the_same_interval_means_for_every_original_surface() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let make = || attach(system(&field, medium.clone(), &primary, 256, &gate), &gate);
    let propagation = GasPropagation::classical(GasPropagationLimits::surface_reference());
    // The Rayleigh receiver admits output clocks from 8 kHz, so the output is
    // faster than the 8 kHz original clock: two original intervals cover three
    // 12 kHz output intervals, as in the upstream observation test.
    let output_rate = ACOUSTIC_RATE * 3 / 2;
    let (originals_per_group, interval_weights) = interval_means(ACOUSTIC_RATE, output_rate);
    assert_eq!((originals_per_group, interval_weights.len()), (2, 3));
    let mut whole = observed(make(), output_rate, 96, 96, Arc::clone(&gate))
        .with_observer_gas(&propagation).unwrap();
    let mut parts = observed(make(), output_rate, 96, 7, Arc::clone(&gate))
        .with_observer_gas(&propagation).unwrap();
    let mut direct = make();
    let mut bore = whole.baffled_receiver().unwrap().clone();
    let mut mouths: Vec<_> = [1, 2].map(|aperture|
        (aperture, whole.reservoir_receiver_at(aperture).unwrap().clone())).to_vec();
    let mut expected = Vec::new();
    for group in whole.drive_history().chunks_exact(originals_per_group) {
        let mut physical = Vec::new();
        for &drive in group {
            let frame = direct.step_with_drive(drive, &gate).unwrap();
            let terminal_flows: Vec<_> = direct.load().system().load().last_terminal_frames()
                .iter().map(|terminal| terminal.outward_flow_m3_s).collect();
            let escape_flows: Vec<Vec<_>> = mouths.iter().map(|(aperture, receiver)| {
                let flows: Vec<_> = frame.reservoir_flows.iter()
                    .filter(|flow| flow.aperture == *aperture)
                    .map(|flow| flow.outward_volume_flow_m3_s).collect();
                assert_eq!(flows.len(), receiver.input_count());
                flows
            }).collect();
            physical.push((terminal_flows, escape_flows));
        }
        // These overlap fractions come from the two clocks, independently of
        // the production interval clock.
        for weights in &interval_weights {
            let mut terminal_flows = vec![0.0; primary.terminal_ports().len()];
            let mut escape_flows: Vec<_> = mouths.iter()
                .map(|(_, receiver)| vec![0.0; receiver.input_count()]).collect();
            for (&weight, (terminals, escapes)) in weights.iter().zip(&physical) {
                if weight == 0.0 { continue; }
                for (sum, flow) in terminal_flows.iter_mut().zip(terminals) {
                    *sum += weight * flow;
                }
                for (sums, flows) in escape_flows.iter_mut().zip(escapes) {
                    for (sum, flow) in sums.iter_mut().zip(flows) { *sum += weight * flow; }
                }
            }
            let mut pressure = bore.step_trace(&terminal_flows).unwrap();
            for ((_, receiver), flows) in mouths.iter_mut().zip(escape_flows) {
                pressure += receiver.step_trace(&flows).unwrap();
            }
            expected.push(pressure);
        }
    }
    let mut a = [0.0; 96]; let mut b = [0.0; 96];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(7) { parts.block(chunk).unwrap(); }
    assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    let scale = expected.iter().chain(&a).map(|x| x.abs()).fold(0.0_f64, f64::max);
    assert!(scale > 0.0);
    for (&actual, expected) in a.iter().zip(expected) {
        assert!((actual - expected).abs() <= 1024.0 * f64::EPSILON * scale);
    }
    // 96 outputs at 12 kHz span 96 * 8000 / 12000 original 8 kHz intervals.
    assert_eq!(whole.system().accepted_steps(), 64);
    assert_eq!(whole.samples_rendered(), 96);
    same_physical_owners(whole.system(), &direct);
    same_physical_owners(whole.system(), parts.system());
}

#[test]
fn all_original_surfaces_share_the_material_budget_and_numerical_latency() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let make = || attach(system(&field, medium.clone(), &primary, 256, &gate), &gate);
    let propagation = GasPropagation::classical(GasPropagationLimits::surface_reference());
    let renderer = observed(make(), ACOUSTIC_RATE, 1, 1, Arc::clone(&gate))
        .with_observer_gas(&propagation).unwrap();
    let bore = renderer.observer_gas().unwrap();
    for aperture in [1, 2] {
        let report = renderer.reservoir_observer_gas_at(aperture).unwrap();
        assert_eq!(report.gas, bore.gas);
        assert_eq!(report.law, bore.law);
        assert_eq!(report.additional_latency_samples, bore.additional_latency_samples);
        assert_eq!(report.surface_filter_indices.len(), 1);
        assert_eq!(report.filters.len(), 1);
    }
    assert_ne!(renderer.reservoir_observer_gas_at(1).unwrap().filters[0].range_m,
        renderer.reservoir_observer_gas_at(2).unwrap().filters[0].range_m);
    // A budget that admits exactly the bore must not be multiplied when escape
    // sources are added. The original complete path set charges it only once.
    let bounded = GasPropagation::classical(GasPropagationLimits {
        max_unique_paths: bore.filters.len(), ..propagation.limits
    });
    assert!(playback(make(), ACOUSTIC_RATE, 1, 1, Arc::clone(&gate))
        .with_observer_gas(&bounded).is_ok());
    assert!(observed(make(), ACOUSTIC_RATE, 1, 1, Arc::clone(&gate))
        .with_observer_gas(&bounded).is_err());
    assert_eq!(renderer.samples_rendered(), 0);
    assert_eq!(renderer.system().accepted_steps(), 0);
    assert!(renderer.with_reservoir_observation_at(0, escape()).is_err());
}

#[test]
fn observation_preserves_bound_geometry_and_explicit_aperture_identity() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let make = || attach(system(&field, medium.clone(), &primary, 256, &gate), &gate);
    let render = || playback(make(), ACOUSTIC_RATE, 1, 1, Arc::clone(&gate));
    assert!(render().with_reservoir_observation(escape()).is_err());
    assert!(render().with_reservoir_observation_at(0, escape()).is_err());
    for input in [
        ReservoirObservation { origin_m: [-0.07, 0.0], ..escape() },
        ReservoirObservation { angle_rad: 0.1, ..escape() },
    ] {
        assert!(render().with_reservoir_observation_at(1, input).is_err());
    }
    let duplicate = render().with_reservoir_observation_at(1, escape()).unwrap();
    assert!(duplicate.with_reservoir_observation_at(1, escape()).is_err());
    let additional = NsfThermoviscousDgPerformance::new_with_drives(make(), 0,
        ApertureObservation::Inlet, ACOUSTIC_RATE, controls(1), 1, Arc::clone(&gate)).unwrap();
    assert!(matches!(additional.with_reservoir_observation_at(1, escape()),
        Err(RenderError::Sizing { what })
            if what == "joint reservoir observation requires its primary duct's original shared exterior"));
}
