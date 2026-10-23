use super::*;
use file::FlueDgFilePerformance;
use fs_exec::CancelGate;

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
