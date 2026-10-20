//! Material card -> physical slit/contact/history -> shared jet -> real WAV.
use super::*;
use fs_couple::flue_loop::network::performance::file::{
    boundary_closing_force_path, flue_source_control_paths,
};
use fs_scenario::gesture::{
    GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue,
};

fn input(force: bool) -> String {
    let original = short();
    let header = original.split_once("\ngesture\n").unwrap().0
        .replace("tube 0.3 0.006 1e-10 -0.9", concat!(
            "network 4 3\nduct_node inlet\nduct_node junction\nduct_node reflection -0.9\nduct_node driven\n",
            "duct_section 0 1 0.15 0.006 1e-10\nduct_section 1 2 0.15 0.006 1e-10\n",
            "duct_section 1 3 0.017 0.003 1e-10"))
        .replace("\nflue ", concat!(
            "\nboundaries 1\nboundary 3 instrument key-1 0 0\n",
            "boundary_plate wall.fsmcdpk 293.15 1000 0.025 0.01 0.0003 6 2 0.0002 0.0002 0\n",
            "boundary_mode 0 600 0 0.25 0.1\nboundary_contact 10000000000 2 0.5\nboundary_slit 0.0003\n",
            "boundary_loss maxwell 1 0 1000 0.1 0.1\nboundary_term boundary_modulus boundary_time\nflue "))
        .replace("viscothermal 8 300 500", "viscothermal sections 1\nloss_section 0 4 300 500")
        .replace("\ncontrols ", "\noutlets 2 0\ncontrol_source instrument flue\ncontrols ");
    let path = boundary_closing_force_path("instrument", "key-1");
    let schedule = GestureSchedule::try_new_v2(
        1000,
        vec![
            GestureTrack {
                id: "blow".into(),
                target: GestureTarget::Parameter(
                    flue_source_control_paths("instrument", "flue").pressure,
                ),
                initial: GestureValue::PressurePa(300.0),
                events: vec![],
            },
            GestureTrack {
                id: "key".into(),
                target: GestureTarget::Parameter(path),
                initial: GestureValue::ForceN(0.0),
                events: if force {
                    vec![
                        GestureEvent {
                            time_s: 0.003,
                            transition_s: 0.002,
                            value: GestureValue::ForceN(0.002),
                        },
                        GestureEvent {
                            time_s: 0.010,
                            transition_s: 0.003,
                            value: GestureValue::ForceN(0.0),
                        },
                    ]
                } else {
                    vec![]
                },
            },
        ],
    )
    .unwrap();
    format!(
        "{header}\ngesture\n{}",
        String::from_utf8(schedule.to_canonical_bytes()).unwrap()
    )
}
fn rendered(text: &str, material: &NormalizedMaterialCardPack, block: usize) -> Vec<f64> {
    let performance = loaded(text, material, block);
    assert_eq!(performance.boundary_records().len(), 1);
    let record = &performance.boundary_records()[0];
    assert_eq!(record.node(), 3);
    assert_eq!(record.pack(), material);
    let network = performance.renderer().voice().network().unwrap();
    assert_eq!(network.boundary_node(record.closing_force()), Some(3));
    let aperture = network.boundary_apertures().unwrap().aperture(3).unwrap();
    assert!(
        aperture
            .plate_reduction()
            .unwrap()
            .material_chart()
            .is_some()
    );
    assert_eq!(aperture.relaxation().unwrap().branches().len(), 1);
    assert_eq!(
        aperture.spec().density_kg_m3,
        network
            .boundary_apertures()
            .unwrap()
            .gas(3)
            .unwrap()
            .density
    );
    assert_eq!(
        aperture.viscous_slit().unwrap().gas(),
        network.boundary_apertures().unwrap().gas(3)
    );
    assert!(network.graph().sections.len() > network.source().graph.sections.len());
    let mut output = vec![0.0; performance.info().samples as usize];
    let mut renderer = performance.into_renderer();
    for samples in output.chunks_mut(block) {
        renderer.block(samples).unwrap();
        let n = renderer.voice().network().unwrap();
        let frame = n.boundary_frame().unwrap();
        let scale = frame.stored_energy_j
            + frame.storage_change_j.abs()
            + frame.dissipated_energy_j
            + frame.network.inlet_work_j.abs()
            + frame.mechanical_work_j.abs();
        assert!(frame.balance_residual_j().abs() <= 3e-8 * scale.max(1e-20));
        assert_eq!(
            n.boundary_apertures()
                .unwrap()
                .aperture(3)
                .unwrap()
                .accepted_steps(),
            renderer.voice().samples_rendered()
        );
    }
    assert!(output.iter().all(|p| p.is_finite()));
    assert!(output.iter().any(|p| p.abs() > 1e-10));
    assert!(
        renderer
            .voice()
            .network()
            .unwrap()
            .boundary_apertures()
            .unwrap()
            .aperture(3)
            .unwrap()
            .relaxation()
            .unwrap()
            .dissipated_energy_j()
            > 0.0
    );
    output
}

#[test]
fn sourced_mass_elasticity_and_force_ramps_reach_the_same_actual_microphone() {
    let source = input(true);
    let material = pack(900.0, 4e9, None, 100_000.0);
    let reference = rendered(&source, &material, 37);
    assert_eq!(reference, rendered(&source, &material, 1));
    assert_eq!(reference, rendered(&source, &material, 512));
    assert_ne!(reference, rendered(&input(false), &material, 37));
    for changed in [
        pack(1800.0, 4e9, None, 100_000.0),
        pack(900.0, 8e9, None, 100_000.0),
    ] {
        assert_ne!(
            reference,
            rendered(&source, &changed, 37),
            "a data-only material replacement must reach pressure after coupled mechanical motion"
        );
    }
}

#[test]
fn the_existing_cli_renders_material_boundary_forces_through_the_sole_pcm_owner() {
    let dir = directory();
    let text = input(true);
    let material = pack(900.0, 4e9, None, 100_000.0);
    let expected = encode_pcm16_wav(&rendered(&text, &material, 37), 48_000, 1.0)
        .unwrap()
        .0;
    let source = dir.join("moving-boundary.performance");
    std::fs::write(&source, text).unwrap();
    std::fs::write(dir.join("wall.fsmcdpk"), material.to_bytes()).unwrap();
    for block in [1, 37, 512] {
        let output = dir.join(format!("moving-boundary-{block}.wav"));
        success(run(&source, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let provenance = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in [
            "moving_boundaries",
            "exact-zoh",
            "closing_force",
            "maxwell_terms",
        ] {
            assert!(
                provenance.contains(token),
                "missing physical boundary source {token}"
            );
        }
    }
    let changed = pack(1800.0, 8e9, None, 100_000.0);
    let text = input(true);
    let changed_pcm = encode_pcm16_wav(&rendered(&text, &changed, 37), 48_000, 1.0)
        .unwrap()
        .0;
    std::fs::write(dir.join("wall.fsmcdpk"), changed.to_bytes()).unwrap();
    let output = dir.join("moving-boundary-replaced-material.wav");
    success(run(&source, &output, 37, false));
    assert_eq!(std::fs::read(output).unwrap(), changed_pcm);
    assert_ne!(
        changed_pcm, expected,
        "material substitution must survive the sole PCM quantizer"
    );
}

#[test]
fn missing_mechanics_material_terms_or_bound_force_path_refuse_during_admission() {
    let text = input(true);
    let material = pack(900.0, 4e9, None, 100_000.0);
    for invalid in [
        text.replace("duct_node driven", "duct_node reflection 1"),
        text.replace("boundary 3 instrument key-1", "boundary 2 instrument key-1"),
        text.replace(
            "boundary 3 instrument key-1",
            "boundary 3 instrument another-key",
        ),
        text.replace("boundary_time", "missing_relaxation_time"),
        text.replace("boundary_contact 10000000000", "boundary_contact -1"),
        text.replace(
            "control_source instrument flue",
            "control_source instrument another-source",
        ),
        text.replace("boundaries 1", "boundaries 65"),
    ] {
        assert!(
            FlueFilePerformance::from_bytes_with_sources(
                invalid.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(material.clone())
            )
            .is_err()
        );
    }
    let missing = pack(900.0, 4e9, Some("boundary_modulus"), 100_000.0);
    assert!(
        FlueFilePerformance::from_bytes_with_sources(
            text.as_bytes(),
            37,
            &CancelGate::new_clock_free(),
            None,
            |_| Ok(missing.clone())
        )
        .is_err()
    );
}
