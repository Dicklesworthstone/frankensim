//! Original source decoding -> stateful shared mechanics -> the existing PCM path.
//! Manufactured contacts test composition; they identify no real bow or rosin.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::FrictionDrive;

fn static_claims(source: &ClaimSet) -> ClaimSet {
    // The fixed Mindlin owner requires genuinely static source coefficients;
    // the general normal-stop fixture's frequency applicability is not allowed.
    // Only the native contact's selected static coefficients have no frequency
    // restriction. Original string viscosity/torsional source bands stay intact.
    let mut claims = ClaimSet::new();
    for id in source.observation_ids() {
        claims.register_observation(source.observation(id).unwrap().clone()).unwrap();
    }
    for (_, original) in source.claims_ordered() {
        let mut claim = original.clone();
        if matches!(claim.key.name(), "density" | "young_modulus" | "poisson_ratio"
            | "yield_stress" | ADHESION_ENERGY_PROPERTY | "mu_local") {
            claim.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
        }
        claims.insert_claim(claim).unwrap();
    }
    claims
}
fn contact_bulk(path: &str, young: f64) -> NormalizedPack {
    let original = bulk(path, young);
    if path == "wire.fsmatpk" || path == "receiver.fsmatpk" {
        pack_claims(path, static_claims(original.claims()))
    } else { original }
}
fn contact_interface(outer: &str) -> NormalizedInterfacePack {
    let original = interface(outer, "receiver.fsmatpk");
    let mut claims = original.claims_pack().claims().clone();
    claims.insert_claim(claim("mu_local", 0.6, Dims::NONE)).unwrap();
    NormalizedInterfacePack::new(original.card().surface_a().clone(),
        original.card().surface_b().clone(), original.card().context().clone(),
        pack_claims("mindlin-interface", static_claims(&claims))).unwrap()
}
fn specification(gate: &CancelGate, outer: &str) -> input::Spec {
    let text = String::from_utf8(spec().bytes).unwrap().replace(
        "friction stribeck 0.8 0.4 0.1 0",
        "friction contacts 1\nbow-mindlin bow receiver.fsmatpk contact.fsintpk normal.fsmodelpk mu_local 0.0002 4e-8 8e-8 0.1 0.00005 0.00005 0.01 0.02 0.03 512");
    let mut input = input::Spec::parse(text.as_bytes()).unwrap();
    input.forces = None;
    input.schedule = GestureSchedule::try_new(700, vec![GestureTrack {
        id: "bow".into(), target: fs_scenario::gesture::GestureTarget::BowStroke { string: 0 },
        initial: GestureValue::Bow { velocity_m_per_s: 0.0003, normal_force_n: 0.0, station: 0.2 },
        events: vec![GestureEvent { time_s: 0.001, transition_s: 0.0,
            value: GestureValue::Bow { velocity_m_per_s: -0.0003, normal_force_n: 0.0, station: 0.2 } }],
    }]).unwrap();
    input.friction.load_mindlin_sources(gate, |path| Ok(match path {
        "contact.fsintpk" => contact_interface(outer).to_bytes(),
        "normal.fsmodelpk" => models_for(None).to_bytes(),
        _ => panic!("unexpected original source {path}"),
    })).unwrap();
    input
}
fn ready(cx: &Cx<'_>, gate: &CancelGate, young: f64) -> Ready {
    specification(gate, "wire.fsmatpk").compile(cx, gate, |path| Ok(contact_bulk(path, young))).unwrap()
}

#[test]
fn material_swap_changes_shared_pressure_and_history_without_replacing_pcm() {
    scoped(|cx, gate| {
        let mut original = ready(cx, gate, 2e11);
        let mut changed = ready(cx, gate, 4e11);
        let a = original.performance.bank().system().material_mindlin_contact(0).unwrap();
        let b = changed.performance.bank().system().material_mindlin_contact(0).unwrap();
        assert_ne!(a.normal_source().model().identity(), b.normal_source().model().identity());
        assert_ne!(a.law().parameters().shear_moduli_pa, b.law().parameters().shear_moduli_pa);
        assert_ne!(a.normal_response().normal_force_n, b.normal_response().normal_force_n);
        assert_eq!(original.spec.bytes, changed.spec.bytes);
        for (a,b) in original.performance.bank().system().components().iter()
            .zip(changed.performance.bank().system().components()) { assert_eq!(a.modes(), b.modes()); }
        let pressure_original = pressure(&mut original, 1);
        assert!(pressure_original.iter().any(|p| *p != 0.0));
        assert_ne!(pressure_original, pressure(&mut changed, 37));
        let mut replay = ready(cx, gate, 2e11);
        assert_eq!(pressure_original, pressure(&mut replay, 128));
        assert_eq!(original.performance.bank().system().material_mindlin_state(0),
            replay.performance.bank().system().material_mindlin_state(0));
        let frame = original.performance.bank().system().last_frame().unwrap();
        let contact = frame.ports[0].mindlin.unwrap();
        assert_eq!(frame.mindlin_energy_j, contact.stored_energy_j);
        assert_eq!(frame.ports[0].dissipation_j, contact.dissipated_j);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.ports[0].normal_force_n > 0.0, "source preload was replaced by the gesture's zero additional load");
        let (expected, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&pressure_original, 48_000, 0.01).unwrap();
        let mut encoded = ready(cx, gate, 2e11);
        let (wav, summary) = encoded.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn reversal_series_drive_and_nonlinear_storage_share_one_loss_and_rollback_account() {
    scoped(|cx, gate| {
        let build = || {
            let mut input = specification(gate, "wire.fsmatpk");
            input.geometric_stretching = Some((0.2, NonlinearModalConfig {
                max_iterations: 64, maximum_generalized_force: 100000.0,
                force_absolute_tolerance: 1e-11, force_relative_tolerance: 1e-10,
            }));
            let bank = input.compile_bank_with_profiles(cx, gate, None, None, false,
                |path| Ok(contact_bulk(path, 2e11))).unwrap().bank;
            bank.into_parts().1.with_series_compliant_drives(vec![(0, 1e5, 0.02)], gate).unwrap()
        };
        let (mut original, mut replay) = (build(), build());
        let zero = vec![0.0; original.mode_count()];
        let (mut heat, mut returned) = (0.0, false);
        for sample in 0..64 {
            let drives = [FrictionDrive { speed_m_s: if sample < 16 { 0.0003 } else { -0.0003 }, normal_force_n: 0.0 }];
            let frame = original.step_under_gate(&zero, &drives, gate).unwrap().clone();
            assert_eq!(&frame, replay.step_under_gate(&zero, &drives, gate).unwrap());
            let port = frame.ports[0]; let contact = port.mindlin.unwrap();
            heat += contact.dissipated_j;
            returned |= contact.stored_change_j < 0.0;
            assert!(contact.dissipated_j >= 0.0 && contact.stored_energy_j >= 0.0);
            assert_eq!(port.dissipation_j, contact.dissipated_j);
            assert!((contact.stored_change_j + contact.dissipated_j
                - port.traction_n * port.slip_distance_m).abs() <= frame.energy_tolerance_j);
            assert_eq!(frame.friction_dissipation_j, contact.dissipated_j + port.series_drive_dissipation_j);
            assert!(frame.nonlinear.is_some());
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(heat > 0.0 && returned, "reversal must return stored shear while irreversible heat remains separate");
        let state = original.material_mindlin_state(0).unwrap().clone();
        let frame = original.last_frame().unwrap().clone();
        let energy = original.total_energy_j().unwrap();
        let bad = [FrictionDrive { speed_m_s: 0.0003, normal_force_n: 1e-3 }];
        assert!(original.step(&zero, &bad).is_err());
        let cancelled = CancelGate::new_clock_free(); cancelled.request();
        let good = [FrictionDrive { speed_m_s: 0.0003, normal_force_n: 0.0 }];
        assert!(original.step_under_gate(&zero, &good, &cancelled).is_err());
        assert_eq!(original.material_mindlin_state(0), Some(&state));
        assert_eq!(original.last_frame(), Some(&frame));
        assert_eq!(original.total_energy_j().unwrap(), energy);
        assert_eq!(original.step(&zero, &good).unwrap(), replay.step(&zero, &good).unwrap());
    });
}

#[test]
fn mismatched_source_and_changing_normal_thermal_or_geometry_inputs_refuse_cold() {
    scoped(|cx, gate| {
        assert!(specification(gate, "wrong-outer.fsmatpk")
            .compile(cx, gate, |path| Ok(contact_bulk(path, 2e11))).is_err());
        for selection in 0..4 {
            let mut input = specification(gate, "wire.fsmatpk");
            match selection {
                0 => {
                    let mut tracks = input.schedule.tracks().to_vec();
                    tracks[0].events[0].value = GestureValue::Bow {
                        velocity_m_per_s: 0.0003, normal_force_n: 0.01, station: 0.2 };
                    input.schedule = GestureSchedule::try_new(700, tracks).unwrap();
                },
                1 => input.contact_thermal = true,
                2 => input.bow_ribbon = Some((0.001, 3)),
                _ => input.bow_residual_flexibility = true,
            }
            let mut reads = 0;
            assert!(input.compile(cx, gate, |path| { reads += 1; Ok(contact_bulk(path, 2e11)) }).is_err());
            assert_eq!(reads, 0, "unsupported composition must refuse before constructing the physical bank");
        }
    });
}
