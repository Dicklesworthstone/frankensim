//! G0/G3: the same source physics reaches the original spatial material bore.
use super::*;
use crate::flue_loop::nolle_jet::{ChannelClass, NolleSpec};

const NOLLE: &str = "jet_transfer nolle-price-2015 long-channel 0\n";
const MOLECULAR: &str =
    "dg_molecular declared 2 N2 1 vibration 1 1000 1 power 10 293.15 0 273.15 323.15 O2 0\n";

fn nolle_input() -> String {
    input()
        .replace("jet_limits", &format!("{NOLLE}jet_limits"))
        .replace("\nflue ", &format!("\n{MOLECULAR}flue "))
}

fn render(text: &str, block: usize) -> (Vec<f64>, FlueDgPerformance) {
    let source =
        FlueDgFilePerformance::from_bytes(text.as_bytes(), block, &CancelGate::new_clock_free())
            .unwrap();
    let count = source.info().samples as usize;
    let mut renderer = source.into_renderer();
    let mut output = vec![0.0; count];
    for chunk in output.chunks_mut(block) {
        renderer.block(chunk).unwrap();
    }
    (output, renderer)
}

#[test]
fn nolle_file_uses_the_original_transfer_and_material_bore_without_another_gain() {
    let data = nolle_input();
    let source =
        FlueDgFilePerformance::from_bytes(data.as_bytes(), 512, &CancelGate::new_clock_free())
            .unwrap();
    let original = source.renderer();
    assert!(original.voice().nolle_jet().is_some());
    assert!(
        original
            .voice()
            .dg()
            .unwrap()
            .load()
            .duct()
            .molecular_relaxation()
            .is_some()
    );
    // The expected performer attaches the published owner directly to the
    // same cold bore. It does not use the new source-selection adapter.
    let voice = replica_voice(&source)
        .with_nolle_jet(NolleSpec::price_2015(ChannelClass::Long, 0.0).unwrap())
        .unwrap();
    let mut direct =
        FlueDgPerformance::new(voice, original.schedule().clone(), *original.config()).unwrap();
    let mut expected = vec![0.0; 512];
    direct.block(&mut expected).unwrap();
    let mut from_file = source.into_renderer();
    let mut actual = vec![0.0; 512];
    from_file.block(&mut actual).unwrap();
    assert_eq!(actual, expected);
    assert!(actual.iter().any(|p| p.abs() > 1e-12));
    assert_eq!(
        from_file.voice().dg().unwrap().load().state(),
        direct.voice().dg().unwrap().load().state()
    );
    let (partitioned, state) = render(&data, 37);
    assert_eq!(actual, partitioned);
    assert_eq!(
        from_file.voice().dg().unwrap().load().state(),
        state.voice().dg().unwrap().load().state()
    );
    let changed_gain = data.replace("jet 8 1e-6", "jet 800 1e-6");
    assert_eq!(
        render(&changed_gain, 37).0,
        actual,
        "Nolle owns receptivity; the card's authored gain must not reenter it"
    );
    let (island, _) = render(&data.replace(NOLLE, ""), 37);
    assert_ne!(island, actual);
}

#[test]
fn nolle_material_gas_changes_source_motion_and_the_original_molecular_pressure() {
    let data = nolle_input();
    let (pressure, original) = render(&data, 37);
    for changed in [
        data.replace("293.15 101325", "303.15 101325"),
        data.replace("293.15 101325", "293.15 70000"),
        data.replace("N2 0.8 O2 0.2", "N2 0.6 O2 0.4"),
    ] {
        let (altered, renderer) = render(&changed, 37);
        assert_ne!(renderer.voice().gas(), original.voice().gas());
        assert_eq!(
            renderer.voice().gas(),
            renderer.voice().dg().unwrap().load().duct().inlet_gas()
        );
        assert_ne!(
            renderer.voice().jet_speed_m_s(),
            original.voice().jet_speed_m_s()
        );
        assert_ne!(
            renderer.voice().dg().unwrap().load().state(),
            original.voice().dg().unwrap().load().state()
        );
        assert!(
            pressure
                .iter()
                .zip(altered)
                .any(|(a, b)| (a - b).abs() > 1e-12)
        );
    }
}

#[test]
fn developing_source_centreline_and_mouth_loss_share_the_actual_inlet_material() {
    let data = input()
        .replace(
            "windway slit 0.005",
            "windway slit-developing 0.015\nmouth_loss contraction 0.8 0.3 manufactured-fixture",
        )
        .replace("jet_limits", "jet_speed card-centreline\njet_limits")
        .replace("pressure-pa\t2e1", "pressure-pa\t1e2");
    let (pressure, renderer) = render(&data, 37);
    let voice = renderer.voice();
    let (law, held) = voice.developing_windway().unwrap();
    held.validity().unwrap();
    assert_eq!(law.gas(), voice.gas());
    let independent = crate::bernoulli_aperture::developing::DevelopingSlit::from_gas(
        voice.dg().unwrap().load().duct().inlet_gas(),
        0.015,
    )
    .unwrap()
    .held(0.010, 0.001, 100.0)
    .unwrap();
    assert_eq!(*held, independent);
    assert_eq!(voice.jet_speed_m_s(), held.bulk_speed_m_s);
    assert_eq!(
        voice.card_speed_m_s(),
        held.core_speed_m_s.min(1.5 * held.bulk_speed_m_s)
    );
    assert_eq!(voice.mouth_contraction().unwrap().gas(), voice.gas());
    assert_eq!(voice.mouth_contraction().unwrap().opening_m(), 0.004);
    assert!(
        voice
            .diagnostics()
            .iter()
            .any(|d| d.separation_loss_j > 0.0)
    );
    for frame in voice.diagnostics() {
        let scale = frame.source_work_j.abs()
            + frame.pipe_work_j.abs()
            + frame.stored_delta_j.abs()
            + frame.separation_loss_j;
        assert!(frame.ledger_defect_j.abs() < 1e-8 * scale.max(1e-20));
    }
    assert_eq!(render(&data, 512).0, pressure);
    let changed = data.replace("293.15 101325", "313.15 101325");
    let (altered, other) = render(&changed, 37);
    assert_ne!(held, other.voice().developing_windway().unwrap().1);
    assert_ne!(pressure, altered);
}

#[test]
fn incompatible_source_profiles_and_invalid_new_records_refuse() {
    let data = nolle_input();
    for refused in [
        data.replace("windway slit 0.005", "windway unsteady 0.005 8 8"),
        data.replace("windway slit 0.005", "windway slit-developing 0.015"),
        data.replace(NOLLE, &format!("jet_speed card-centreline\n{NOLLE}")),
        data.replace(
            NOLLE,
            &format!("{NOLLE}jet_response forced-bickley 0.1 90\n"),
        ),
        data.replace("long-channel 0", "unknown-channel 0"),
        data.replace("long-channel 0", "long-channel -0.001"),
        data.replace("long-channel 0", "long-channel 0 extra"),
        data.replace(
            "windway slit 0.005",
            "windway slit 0.005\nmouth_loss contraction 1.1 0.3 fixture",
        ),
        data.replace(
            "windway slit 0.005",
            "windway slit 0.005\nmouth_loss contraction 0.8 0.3",
        ),
    ] {
        assert!(
            FlueDgFilePerformance::from_bytes(
                refused.as_bytes(),
                37,
                &CancelGate::new_clock_free()
            )
            .is_err(),
            "admitted {refused}"
        );
    }
}

#[test]
fn omitted_flow_selection_keeps_the_same_explicit_ideal_source() {
    let explicit = input().replace("windway slit 0.005", "windway none");
    let omitted = explicit.replace("windway none\n", "");
    assert_eq!(render(&explicit, 37).0, render(&omitted, 37).0);
}

#[test]
fn canonical_nolle_source_retains_molecular_bore_and_original_exterior_microphone() {
    let data = include_str!("../../../../../examples/flue-dg-nolle-molecular.performance");
    let (pressure, renderer) = render(data, 113);
    assert!(renderer.voice().nolle_jet().is_some());
    let load = renderer.voice().dg().unwrap().load();
    assert!(load.duct().molecular_relaxation().is_some());
    assert!(load.radiation_record().is_some());
    assert!(renderer.baffled_receiver().is_some());
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-14));
    assert_eq!(load.state().step_count(), 2048);
}
