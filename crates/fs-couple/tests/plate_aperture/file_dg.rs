use super::*;
use fs_couple::bernoulli_aperture::paired::BladePropagator;
use fs_couple::bernoulli_aperture::performance::{AperturePerformance, CoupledAperture};

const DG: &str = include_str!("../../examples/plate-valve-dg.performance");

fn short_source() -> String {
    let header = DG
        .split_once("\nschedule\n")
        .unwrap()
        .0
        .replace("audio 96000 9602", "audio 96000 512");
    format!("{header}\nschedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\ntrack\tmouth\tblowing-pressure\ninitial\tpressure-pa\t5e0\nevents\t0\n")
}

fn trace(text: &str, block: usize) -> (Vec<f64>, AperturePerformance) {
    let source =
        PlateValvePerformance::from_bytes(text.as_bytes(), block, &CancelGate::new()).unwrap();
    let mut result = vec![0.0; source.info().samples as usize];
    let mut renderer = source.into_renderer();
    for chunk in result.chunks_mut(block) {
        renderer.block(chunk).unwrap();
    }
    (result, renderer)
}

#[test]
fn data_selects_original_exact_dg_plate_contact_and_material_histories() {
    let input = short_source();
    let (expected, source) = trace(&input, 37);
    let CoupledAperture::Dg(dg) = source.system() else {
        panic!("DG selector replaced its physical owner")
    };
    assert_eq!(dg.duct_state().step_count(), 512);
    assert_eq!(dg.aperture().accepted_steps(), 512);
    assert_eq!(dg.aperture().blade_propagator(), BladePropagator::ExactZoh);
    assert_eq!(dg.load().duct().total_dofs(), 108);
    assert_eq!(dg.load().duct().mesh().elements, 2);
    assert_eq!(dg.load().duct().total_length(), 0.05);
    assert!(dg.load().duct().viscothermal_model().is_some());
    assert!(dg.aperture().plate_reduction().unwrap().material_chart().is_none());
    assert!(dg.aperture().plate_closure().is_some());
    assert_eq!(dg.aperture().relaxation().unwrap().branches().len(), 2);
    assert_eq!(
        dg.aperture().viscous_slit().unwrap().dynamic_viscosity_pa_s(),
        dg.load().duct().inlet_gas().dynamic_viscosity
    );
    assert!(expected.iter().any(|p| p.abs() > 1e-6));
    assert_eq!(trace(&input, 1).0, expected);
    assert_eq!(trace(&input, 512).0, expected);
}

#[test]
fn gas_solid_and_physical_force_data_change_the_live_dg_trace() {
    let input = short_source();
    let reference = trace(&input, 37).0;
    let forced = input.replace(
        "compile_limits",
        "force_ports 1\nforce_port node 13\nforce_events 2\nforce_event 64 0 -0.00001\nforce_event 128 0 0\ncompile_limits",
    );
    for changed in [
        input.replace("293.15 101325", "293.15 60000"),
        input.replace("O2 0.3", "CO2 0.3"),
        input.replace("0.0003 900 4000000000", "0.0003 1800 4000000000"),
        forced,
    ] {
        let result = trace(&changed, 37).0;
        assert!(
            reference.iter().zip(result).any(|(a, b)| (a - b).abs() > 1e-6)
        );
    }
}

#[test]
fn dg_selector_refuses_unbound_observers_budget_and_unresolved_diffusion() {
    let input = short_source();
    for changed in [
        input.replace("observation inlet", "observation terminal"),
        input.replace("observation inlet", "observation baffled-outlet 0 0 0.2 8 32 1000"),
        input.replace("1e-11 128", "1e-11 107"),
        input.replace("open 2 2 8", "open 2 2 1"),
        input.replace("open 2 2 8", "unknown 2 2 8"),
        input.replace("50 2000", "50 50000"),
        input.replace("dg_tube ", "exterior 293.15 101325 0\ndg_tube "),
    ] {
        assert!(
            PlateValvePerformance::from_bytes(changed.as_bytes(), 37, &CancelGate::new()).is_err()
        );
    }
}
