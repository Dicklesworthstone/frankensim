use super::*;
use file::FlueDgFilePerformance;
use fs_exec::CancelGate;

mod mouth;
mod channels;

fn input() -> String {
    let example = include_str!("../../../../examples/flue-dg.performance");
    let header = example.split_once("\ngesture\n").unwrap().0
        .replace("output 2048", "output 512");
    format!("{header}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t2e1\nevents\t0\n")
}

#[test]
fn g0_file_pressure_is_the_accepted_original_dg_port_mean() {
    let data = input();
    let source = FlueDgFilePerformance::from_bytes(data.as_bytes(), 37,
        &CancelGate::new_clock_free()).unwrap();
    assert_eq!(source.info().jet_card_hash, source.card().content_hash());
    let mut renderer = source.into_renderer();
    let mut squared = 0.0;
    for _ in 0..512 {
        let mut sample = [0.0];
        renderer.block(&mut sample).unwrap();
        assert_eq!(sample[0], renderer.voice().dg().unwrap().last_frame().unwrap().inlet_pressure_pa);
        squared += sample[0] * sample[0];
    }
    assert!(squared > 1e-8);
    assert_eq!(renderer.samples_rendered(), 512);
}

#[test]
fn g3_file_gas_and_callback_choices_reach_the_existing_physical_source() {
    let data = input();
    let mut a = FlueDgFilePerformance::from_bytes(data.as_bytes(), 512,
        &CancelGate::new_clock_free()).unwrap().into_renderer();
    let mut b = FlueDgFilePerformance::from_bytes(data.as_bytes(), 37,
        &CancelGate::new_clock_free()).unwrap().into_renderer();
    let changed = data.replace("293.15 101325", "293.15 60000");
    let mut c = FlueDgFilePerformance::from_bytes(changed.as_bytes(), 512,
        &CancelGate::new_clock_free()).unwrap().into_renderer();
    let mut pa = vec![0.0; 512];
    let mut pb = vec![0.0; 512];
    let mut pc = vec![0.0; 512];
    a.block(&mut pa).unwrap();
    for chunk in pb.chunks_mut(37) { b.block(chunk).unwrap(); }
    c.block(&mut pc).unwrap();
    assert_eq!(pa, pb);
    assert!(pa.iter().zip(&pc).any(|(x, y)| (x-y).abs() > 1e-8));
    assert_eq!(a.voice().dg().unwrap().load().state(), b.voice().dg().unwrap().load().state());
    assert_ne!(a.voice().dg().unwrap().load().state(), c.voice().dg().unwrap().load().state());
}

#[test]
fn g3_authored_taper_and_mixture_reach_the_same_source_and_wall_memory() {
    let data = input().replace(
        "dg_tube 0.05 0.001 open 1 2 8 50 2000 0.01 1e-11 256",
        "dg_bore 2 open 2 2 12 50 2000 0.01 1e-11 224 4\ndg_segment cylinder 0.001 0.025\ndg_segment cone 0.001 0.002 0.025",
    );
    let render = |text: &str, block| {
        let source = FlueDgFilePerformance::from_bytes(text.as_bytes(), block,
            &CancelGate::new_clock_free()).unwrap();
        let duct = source.renderer().voice().dg().unwrap().load().duct();
        assert!(duct.spatial_viscothermal_admission().is_some());
        assert_eq!(duct.total_dofs(), 204);
        assert_eq!(source.geometry().bore_area_m2, duct.inlet_area_m2());
        assert!(duct.outlet_area_m2() > duct.inlet_area_m2());
        let mut renderer = source.into_renderer();
        let mut pressure = vec![0.0; 512];
        for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        pressure
    };
    let reference = render(&data, 37);
    assert_eq!(render(&data, 512), reference);
    for changed in [
        data.replace("O2 0.2", "CO2 0.2"),
        data.replace("cone 0.001 0.002", "cone 0.001 0.0015"),
    ] {
        assert!(render(&changed, 37).iter().zip(&reference)
            .any(|(a, b)| (a - b).abs() > 1e-8));
    }
}

#[test]
fn g0_unbound_exterior_and_nonphysical_source_declarations_refuse() {
    for data in [
        input().replace("observation inlet-pressure", "observation exterior"),
        input().replace("jet_limits 1 0.3", "jet_limits 1 0.4"),
        input().replace("1e-11 256", "1e-11 1"),
        input().replace("controls 128", "controls 0"),
    ] {
        assert!(FlueDgFilePerformance::from_bytes(data.as_bytes(), 37,
            &CancelGate::new_clock_free()).is_err());
    }
}

fn baffled_input() -> String {
    input().replace(
        "dg_tube 0.05 0.001 open 1 2 8 50 2000 0.01 1e-11 256",
        concat!(
            "dg_tube 0.05 0.001 baffled 1 2 8 50 2000 0.01 1e-11 256\n",
            "dg_radiation computed-baffled 20 65 4 8 12 20 0.1 67108864 ",
            "2 100 1e-8 64 1048576 257\n",
            "dg_exterior mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2",
        ),
    ).replace(
        "observation inlet-pressure",
        "observation baffled 0 0 0.2 8 32 2000",
    ).replace(
        "events\t0\n",
        "events\t1\nevent\t1e-3\t0e0\tpressure-pa\t0e0\n",
    )
}

// Reuse the already computed cold load: rebuilding the source voice must not
// fit a second reaction operator merely to choose another observation/block.
fn replica_voice(source: &FlueDgFilePerformance) -> FlueVoice {
    use crate::flue_loop::dg::FlueDgSpec;
    FlueVoice::from_dg(
        Some(source.card()), 8.0, 1e-6, source.geometry(),
        FlueDgSpec {
            load: source.renderer().voice().dg().unwrap().load().clone(),
            minimum_jet_speed_m_s: 1.0,
            maximum_jet_mach: 0.3,
            max_jet_history_samples: 8192,
        },
        source.info().sample_rate_hz, 3, 0, &CancelGate::new_clock_free(),
    ).unwrap().with_viscous_windway(0.005).unwrap()
}

#[test]
fn baffled_pressure_uses_each_accepted_outlet_flow_and_bound_exterior_once() {
    let source = FlueDgFilePerformance::from_bytes(
        baffled_input().as_bytes(), 37, &CancelGate::new_clock_free(),
    ).unwrap();
    let original = source.renderer();
    let load = original.voice().dg().unwrap().load();
    let record = load.radiation_record().unwrap();
    let receiver = original.baffled_receiver_config().unwrap();
    assert_ne!(
        record.gas_record().state().sound_speed.to_bits(),
        load.duct().inlet_gas().sound_speed.to_bits(),
    );
    assert_eq!(record.apertures().len(), 1);
    assert_eq!(record.apertures()[0].center_m, [0.0; 2]);
    assert_eq!(record.band_hz()[1], 0.5 * f64::from(source.info().sample_rate_hz));
    assert!(record.report().maximum_fit_relative_error <= 0.1);
    assert!(record.report().maximum_held_out_relative_error <= 0.1);
    let mut expected = BaffledPressure::circular_outlet_multirate_differentiating(
        record.apertures()[0].radius_m, source.info().sample_rate_hz,
        source.info().sample_rate_hz, receiver, record.gas_record().state().into(),
    ).unwrap();
    let first_arrival = expected.delay_samples.0;
    assert_eq!(
        original.baffled_receiver().unwrap().filter_delay_samples(),
        expected.filter_delay_samples(),
    );
    let mut config = *original.config();
    config.max_block = 512;
    let mut whole = FlueDgPerformance::with_baffled_outlet(
        replica_voice(&source), original.schedule().clone(), config, receiver,
    ).unwrap();
    let mut internal = FlueDgPerformance::new(
        replica_voice(&source), original.schedule().clone(), config,
    ).unwrap();
    config.max_block = 37;
    let mut partitioned = FlueDgPerformance::with_baffled_outlet(
        replica_voice(&source), original.schedule().clone(), config, receiver,
    ).unwrap();
    let mut single = source.into_renderer();
    let mut reference = Vec::with_capacity(512);
    for _ in 0..512 {
        let mut pressure = [0.0];
        single.block(&mut pressure).unwrap();
        let flow = single.voice().dg().unwrap().load()
            .last_terminal_frame().unwrap().outward_flow_m3_s;
        assert_eq!(pressure[0], expected.step_trace(&[flow]).unwrap());
        reference.push(pressure[0]);
    }
    assert!(reference[..first_arrival].iter().all(|p| *p == 0.0));
    assert!(reference.iter().any(|p| p.abs() > 1e-14));
    let mut all_at_once = vec![0.0; 512];
    let mut chunked = vec![0.0; 512];
    whole.block(&mut all_at_once).unwrap();
    for block in chunked.chunks_mut(37) { partitioned.block(block).unwrap(); }
    internal.block(&mut vec![0.0; 512]).unwrap();
    assert_eq!(reference, all_at_once);
    assert_eq!(reference, chunked);
    for other in [&whole, &partitioned, &internal] {
        assert_eq!(single.voice().dg().unwrap().load().state(),
            other.voice().dg().unwrap().load().state());
        assert_eq!(
            single.voice().dg().unwrap().load().last_terminal_frame()
                .unwrap().outward_flow_m3_s.to_bits(),
            other.voice().dg().unwrap().load().last_terminal_frame()
                .unwrap().outward_flow_m3_s.to_bits(),
        );
    }
}

#[test]
fn baffled_source_receiver_and_finite_pressure_limits_are_enforced() {
    let unbound = input().replace(
        "observation inlet-pressure", "observation baffled 0 0 0.2 8 32 2000",
    );
    assert!(FlueDgFilePerformance::from_bytes(
        unbound.as_bytes(), 37, &CancelGate::new_clock_free(),
    ).is_err());
    let source = FlueDgFilePerformance::from_bytes(
        baffled_input().as_bytes(), 512, &CancelGate::new_clock_free(),
    ).unwrap();
    let original = source.renderer();
    let receiver = original.baffled_receiver_config().unwrap();
    for refused in [
        CircularOutletReceiver { position_m: [0.0, 0.0, 0.01], ..receiver },
        CircularOutletReceiver { radial_rings: 0, ..receiver },
        CircularOutletReceiver { maximum_frequency_hz: 47000.0, ..receiver },
    ] {
        assert!(FlueDgPerformance::with_baffled_outlet(
            replica_voice(&source), original.schedule().clone(), *original.config(), refused,
        ).is_err());
    }
    assert_eq!(original.samples_rendered(), 0);
    assert_eq!(original.voice().dg().unwrap().load().state().step_count(), 0);

    let mut config = *original.config();
    config.maximum_abs_pressure_pa = 1e-20;
    let mut limited = FlueDgPerformance::with_baffled_outlet(
        replica_voice(&source), original.schedule().clone(), config, receiver,
    ).unwrap();
    assert!(limited.block(&mut vec![0.0; 512]).is_err());
    assert!(matches!(limited.block(&mut [0.0]), Err(RenderError::Poisoned)));
}

