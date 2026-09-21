use super::*;
use crate::{Mechanics, Stroke, drum_spec, snare::SnareSet};

#[test]
fn cavity_option_reaches_every_existing_drum_and_snare_command() {
    for family in ["drum", "drum-stretch", "drum-modal", "snare", "snare-off"] {
        for suffix in ["", "-wav", "-mic"] {
            let command = format!("{family}{suffix}");
            assert!(admit_command(true, &command).is_ok());
            let mut args = vec![command.clone(), "128".into(), "--cavity-modes".into(),
                "--strike-speed-m-s".into(), "4".into()];
            assert!(option(&mut args).unwrap());
            let (positionals, stroke) = crate::playing::parse(args).unwrap();
            assert_eq!(positionals, [command, "128".to_owned()]);
            assert_eq!(stroke.speed_m_s, 4.0);
        }
    }
    for command in ["splash", "splash-wav", "splash-mic", "unknown"] {
        assert!(admit_command(true, command).is_err());
    }
    assert!(admit_command(false, "splash").is_ok());
}

#[test]
fn twenty_strand_snare_keeps_all_contacts_while_spatial_air_changes_head_motion() {
    let stroke = Stroke { speed_m_s: 4.0, position_m: Some([0.06, 0.01]) };
    let wires = SnareSet::reference(false);
    let mut compact = crate::drum_with_air(256, 2e-6, false, true,
        Some(wires), false, stroke, false, None).unwrap();
    let mut distributed = crate::drum_with_air(256, 2e-6, false, true,
        Some(wires), false, stroke, true, None).unwrap();
    let Mechanics::Prepared(prepared) = &distributed.system else { panic!("must use the prepared contact owner"); };
    assert_eq!(prepared.contact_count(), 1+wires.strands*wires.contact_cells);
    assert_eq!(prepared.contact_count(), 241);
    let solids = compact.force.len();
    assert_eq!(distributed.air.as_ref().unwrap().coupling.structural_modes(), solids);
    assert!(solids > 160 && distributed.force.len() > solids);
    assert_eq!(compact.system.state(), &distributed.system.state()[..2*solids]);
    assert_eq!(compact.observer_a, distributed.observer_a[..solids]);
    assert_eq!(compact.observer_b, distributed.observer_b[..solids]);
    assert!(distributed.observer_a[solids-160..].iter().all(|v| *v == 0.0));
    assert!(distributed.observer_b[solids-160..].iter().all(|v| *v == 0.0));
    assert!(distributed.force.iter().all(|v| *v == 0.0));
    let initial = prepared.frame().stored_energy_j;
    let gate = CancelGate::new_clock_free();
    let mut loss = 0.0; let mut changed = 0.0_f64; let mut spatial_pressure = 0.0_f64;
    for _ in 0..256 {
        compact.system.step(&compact.force, &gate).unwrap();
        let frame = distributed.system.step(&distributed.force, &gate).unwrap();
        loss += frame.dissipated_energy_j;
        let state = distributed.system.state();
        for (a, b) in compact.system.state().iter().zip(state) { changed = changed.max((a-b).abs()); }
        let probe = distributed.air.as_ref().unwrap();
        let (a, b) = probe.points(state).unwrap(); spatial_pressure = spatial_pressure.max((a-b).abs());
        let uniform = crate::cavity_pressure(distributed.pressure.as_ref().unwrap(), state);
        assert!((probe.uniform_pressure(state).unwrap()-uniform).abs() < 1e-7*(1.0+uniform.abs()));
        assert!((frame.stored_energy_j+loss-initial).abs() < 1e-7);
    }
    assert!(changed > 1e-12 && spatial_pressure > 1e-5);
    let accepted = distributed.system.state().to_vec();
    let stopped = CancelGate::new_clock_free(); stopped.request();
    assert!(matches!(distributed.system.step(&distributed.force, &stopped),
        Err(fs_couple::render::plate::impact::ImpactError::Cancelled)));
    assert_eq!(distributed.system.state(), accepted);
}

#[test]
fn prepared_and_reference_drums_share_the_same_distributed_air_and_initial_state() {
    let stroke = Stroke { speed_m_s: 4.0, position_m: Some([0.06, 0.01]) };
    let mut reference = crate::drum_with_air(128, 2e-6, false, false,
        None, false, stroke, true, None).unwrap();
    let mut prepared = crate::drum_with_air(128, 2e-6, false, true,
        None, false, stroke, true, None).unwrap();
    assert_eq!(reference.system.state(), prepared.system.state());
    assert_eq!(reference.observer_a, prepared.observer_a);
    assert_eq!(reference.observer_b, prepared.observer_b);
    assert_eq!(reference.force, prepared.force);
    reference.system = reference.system.into_prepared_nonlinear().unwrap();
    let gate = CancelGate::new_clock_free(); let mut discrepancy = 0.0_f64;
    for _ in 0..128 {
        reference.system.step(&reference.force, &gate).unwrap();
        prepared.system.step(&prepared.force, &gate).unwrap();
        for (i, (a, b)) in reference.system.state().iter().zip(prepared.system.state()).enumerate() {
            discrepancy = discrepancy.max((a-b).abs()*if i%2 == 0 {4000.0} else {1.0});
        }
    }
    assert!(discrepancy < 1e-3, "prepared/reference onset discrepancy {discrepancy}");
}

#[test]
fn supplied_snare_off_audio_construction_preserves_head_projections_and_appends_air_last() {
    let mut spec = drum_spec::Spec::reference();
    spec.radius_m *= 1.03; spec.outer_radius_m *= 1.03; spec.depth_m *= 0.98;
    let wires = SnareSet::reference(true);
    let stroke = Stroke { speed_m_s: 1.6, position_m: Some([0.06, 0.01]) };
    let dt = crate::acoustics::MECHANICAL_DT;
    let compact = crate::drum_with_spec(32, dt, true, true, Some(wires), false,
        stroke, false, None, Some(spec)).unwrap();
    let mut distributed = crate::drum_with_spec(32, dt, true, true, Some(wires), false,
        stroke, true, None, Some(spec)).unwrap();
    assert!(compact.acoustics.is_some() && distributed.acoustics.is_some());
    let original = compact.force.len();
    assert_eq!(compact.system.state(), &distributed.system.state()[..2*original]);
    assert_eq!(compact.observer_a, distributed.observer_a[..original]);
    assert_eq!(compact.observer_b, distributed.observer_b[..original]);
    assert!(distributed.force[original..].iter().all(|f| *f == 0.0));
    assert!(distributed.observer_a[original..].iter().all(|f| *f == 0.0));
    let Mechanics::Prepared(system) = &distributed.system else { panic!(); };
    assert_eq!(system.sample_period_s().to_bits(), dt.to_bits());
    assert_eq!(system.contact_count(), 241);
    // Construction and mechanical stepping only: this does not pretend to
    // execute the separate expensive BEM bake or validate a rendered waveform.
    let gate = CancelGate::new_clock_free();
    for _ in 0..32 { distributed.system.step(&distributed.force, &gate).unwrap(); }
}

#[test]
fn unsupported_stretching_or_vented_snare_is_not_silently_downgraded() {
    let wires = Some(SnareSet::reference(false)); let stroke = Stroke::default();
    assert!(crate::drum_with_air(1, 2e-6, false, true, wires, true, stroke, true, None).is_err());
    assert!(crate::drum_with_air(1, 2e-6, false, false, wires, false, stroke, true, None).is_err());
    let neck = NeckOptions { radius_m: 0.005, effective_length_m: 0.012,
        resistance_pa_s_m3: 1000.0, azimuth_rad: 0.4, axial_position_m: 0.08 };
    for command in ["drum-modal", "snare", "snare-off", "snare-wav", "snare-mic"] {
        assert!(admit_neck_command(Some(neck), true, command).is_err());
    }
    assert!(crate::drum_with_air(1, 2e-6, false, true, wires, false, stroke, true, Some(neck)).is_err());
}
