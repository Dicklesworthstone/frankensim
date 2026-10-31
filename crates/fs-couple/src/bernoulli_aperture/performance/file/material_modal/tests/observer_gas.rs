//! Material propagation through the original two-sided source and microphone.
use super::*;
use crate::air_path::GasAbsorptionModel;
use crate::pcm_wav::baffled::gas::GasPropagation;
use crate::pcm_wav::baffled::openings::duct::{DuctPressure, DuctView};

const CLASSICAL: &str = "stokes-kirchhoff 4096 1000000 50000000";

fn selected(text: &str, law: &str) -> String {
    text.replacen(
        "\nopenings ",
        &format!("\nobserver_gas {law}\nopenings "),
        1,
    )
}

fn two_openings() -> String {
    input()
        .replace(
            "receiver 0 0 0.25 4 16 1000",
            "receiver 0.02 0.01 0.25 1 8 1000",
        )
        .replace(
            "openings 1\nopening downstream 2 0 0",
            "openings 2\nopening upstream 0 -0.05 0\nopening downstream 2 0.03 0",
        )
}

#[test]
fn g1_g3_original_exterior_paths_keep_their_material_and_delay_accounts() {
    let text = two_openings();
    let plain = build(&text, 37);
    let selected = build(&selected(&text, CLASSICAL), 37);
    let source = selected.renderer().source();
    let report = source.observer_gas().unwrap();
    let gas = source.radiation().exterior.unwrap();
    assert_eq!(report.gas, gas);
    assert_eq!(report.law, GasAbsorptionModel::StokesKirchhoff);
    assert_eq!(report.surface_filter_indices.len(), 2 * 8);
    assert!(
        report.filters.len() > 1,
        "distinct physical paths must not collapse to one range"
    );
    assert_ne!(
        report.gas.temperature,
        source.system().source(ValveSide::Downstream).gases[1].temperature
    );
    for path in &report.filters {
        assert!(
            (path.receiver_delay_s + path.low_frequency_delay_s - path.range_m / gas.sound_speed)
                .abs()
                < 1e-14
        );
    }
    assert_eq!(
        source.receiver().filter_delay_samples(),
        plain.renderer().source().receiver().filter_delay_samples()
            + report.additional_latency_samples as f64
    );

    // Independently bind the actual graph/gas owners and drive signed flows
    // through every complete mechanical packet. The file path must retain the
    // original ordering, geometry and multirate kernel, not a center-range filter.
    let system = source.system();
    let views = [ValveSide::Upstream, ValveSide::Downstream].map(|side| DuctView {
        graph: system.graph(side),
        gases: system.section_gases(side),
    });
    let propagation = GasPropagation::read_record(&mut CLASSICAL.split_ascii_whitespace()).unwrap();
    let mut direct = DuctPressure::new(&views, source.radiation().clone())
        .unwrap()
        .with_observer_gas(&propagation, &mut || false)
        .unwrap();
    assert_eq!(direct.observer_gas().unwrap().filters, report.filters);
    let mut bound = source.receiver().clone();
    let ratio = source.mechanics_per_output();
    let mut trace = vec![0.0; ratio * 2];
    for frame in 0..192 {
        let mut expected = None;
        for (substep, flow) in trace.chunks_exact_mut(2).enumerate() {
            let phase = (frame * ratio + substep) as f64 * 0.03;
            flow.copy_from_slice(&[phase.sin() * 1e-6, -phase.cos() * 0.7e-6]);
            expected = direct.push_frame(|network, _| Some(flow[network])).unwrap();
        }
        assert_eq!(bound.step_trace(&trace).unwrap(), expected.unwrap());
    }
    assert_eq!(source.mechanical_samples_rendered(), 0);
}

#[test]
fn g3_observer_selection_changes_pressure_without_changing_loaded_material_mechanics() {
    let text = input();
    let (plain, original) = render(build(&text, 37));
    let (propagated, material) = render(build(&selected(&text, CLASSICAL), 37));
    assert!(propagated.iter().all(|p| p.is_finite()));
    assert!(propagated.iter().any(|p| p.abs() > 1e-15));
    assert_ne!(plain, propagated);
    assert_eq!(
        original.source().last_frame(),
        material.source().last_frame()
    );
    assert_eq!(original.applied_controls(), material.applied_controls());
    let states = |source: &ScheduledTwoSidedPressure| {
        let LoadedAperture::Modal(valve) = source.source().system().valve() else {
            unreachable!()
        };
        valve.modal_states().to_vec()
    };
    assert_eq!(states(&original), states(&material));
    assert_eq!(
        original.source().system().stored_energy_j(),
        material.source().system().stored_energy_j()
    );
    let (partitioned, _) = render(build(&selected(&text, CLASSICAL), 1));
    assert_eq!(propagated, partitioned);
}

#[test]
fn g0_g3_molecular_declaration_uses_the_exterior_species_and_refuses_incomplete_input() {
    // Explicit zero-reservoir He/O2 is a manufactured molecular-law selection,
    // not measured relaxation data. It cannot be evaluated on the N2/O2 inlet.
    let text = input()
        .replace(
            "ambient 293.15 101325 0",
            "ambient mixture 293.15 101325 mole 2 N2 0.79 O2 0.21",
        )
        .replace(
            "exterior 293.15 101325 0",
            "exterior mixture 293.15 101325 mole 2 He 0.79 O2 0.21",
        );
    let law = "molecular 4096 1000000 50000000 declared 2 He 0 O2 0";
    let source = build(&selected(&text, law), 37);
    let report = source.renderer().source().observer_gas().unwrap();
    assert_eq!(report.law, GasAbsorptionModel::MolecularRelaxation);
    assert_eq!(
        report.gas,
        source.renderer().source().radiation().exterior.unwrap()
    );
    assert_ne!(
        report.gas,
        source
            .renderer()
            .source()
            .system()
            .source(ValveSide::Upstream)
            .gases[0]
    );

    for refused in [
        selected(
            &text,
            "molecular 4096 1000000 50000000 declared 2 N2 0 O2 0",
        ),
        selected(&input(), law),
        selected(&text, "stokes-kirchhoff 0 1000000 50000000"),
        selected(&text, "lossless 4096 1000000 50000000 extra"),
        selected(&selected(&text, CLASSICAL), CLASSICAL),
    ] {
        let mut loads = 0;
        assert!(
            MaterialModalValvePerformance::from_bytes_with_sources(
                refused.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| {
                    loads += 1;
                    Ok(pack(900.0, false))
                }
            )
            .is_err()
        );
        assert_eq!(
            loads, 0,
            "malformed gas declarations must refuse before modal material work"
        );
    }
    // A valid but insufficient physical-path budget refuses at receiver admission.
    assert!(
        MaterialModalValvePerformance::from_bytes_with_sources(
            selected(&two_openings(), "lossless 1 1000000 50000000").as_bytes(),
            37,
            &CancelGate::new_clock_free(),
            None,
            |_| Ok(pack(900.0, false))
        )
        .is_err()
    );
}
