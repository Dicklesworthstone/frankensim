//! Source data -> full modal radiation -> original performer and PCM.
use super::*;

fn selected(pressure_pa: f64) -> input::Spec {
    let text = include_str!("../rayleigh-pluck.material-strings")
        .replace("directional-body.fsmatpk", "body.fsmatpk")
        .replace("293.15 101325 mole", &format!("293.15 {pressure_pa} mole"));
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    spec.samples = 48;
    spec
}

fn ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, profile: Profile) -> Ready {
    spec.compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 1.0, None))).unwrap()
}

#[test]
fn computed_air_load_changes_motion_and_streams_the_same_accepted_pressure() {
    scoped(|cx, gate| {
        let spec = selected(101_325.0);
        let mut loaded = ready(cx, gate, spec.clone(), profile());
        let report = loaded.body_radiation.as_ref().unwrap();
        assert_eq!(report.body_modes, 1);
        assert_eq!(report.surface_cells, 96);
        assert!(report.acoustic_poles > 0);
        assert!(report.peak_error <= 0.05 && report.rms_error <= 0.02);
        assert_eq!(report.gas, spec.gas);
        let mut bare_spec = spec.clone();
        bare_spec.body_radiation = None;
        let mut bare = ready(cx, gate, bare_spec, profile());
        assert_eq!(loaded.performance.bank().system().mode_count(), bare.performance.bank().system().mode_count());
        let p = pressure(&mut loaded, 7);
        let q = pressure(&mut bare, 48);
        assert!(p.iter().any(|v| v.abs() > 1e-12));
        assert!(p.iter().zip(&q).any(|(a,b)| (a-b).abs() > 1e-15), "air never reacted on the instrument");
        assert_ne!(loaded.performance.bank().system().components()[0].states(),
            bare.performance.bank().system().components()[0].states(), "radiation never reached the string");
        let mut streamed = ready(cx, gate, spec, profile());
        let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&p, streamed.spec.rate, streamed.spec.full_scale).unwrap().0;
        let (wav, summary) = streamed.write(std::io::Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), encoded);
        assert_eq!(summary.samples, 48);
        for image in [&loaded, &bare, &streamed] {
            let frame = image.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn gas_and_grain_data_rebuild_the_actual_same_basis_radiation_matrix() {
    scoped(|cx, gate| {
        let a = ready(cx, gate, selected(101_325.0), profile());
        let b = ready(cx, gate, selected(202_650.0), profile());
        let turned = Profile::parse(include_str!("../orthotropic.body")
            .replace("grain 0.17", "grain 0.87").as_bytes()).unwrap();
        let c = ready(cx, gate, selected(101_325.0), turned);
        let ar = a.body_radiation.as_ref().unwrap();
        let br = b.body_radiation.as_ref().unwrap();
        let cr = c.body_radiation.as_ref().unwrap();
        assert_eq!(ar.specimen_identity, br.specimen_identity);
        assert_ne!(ar.specimen_identity, cr.specimen_identity);
        assert_eq!(ar.gas.sound_speed, br.gas.sound_speed);
        assert_eq!(2.0*ar.gas.density, br.gas.density);
        let impedance = |ready: &Ready| ready.performance.bank().source().bank()
            .body_radiation().unwrap().model().impedance(core::f64::consts::TAU*100.0).unwrap();
        let za = impedance(&a); let zb = impedance(&b); let zc = impedance(&c);
        assert!((zb[0]-za[0].scale(2.0)).abs() <= 1e-8 * za[0].abs().max(1e-12));
        assert!((zc[0]-za[0]).abs() > 1e-10 * za[0].abs().max(1e-12));
    });
}

#[test]
fn incompatible_or_unbounded_radiation_is_refused_before_pack_io() {
    scoped(|cx, gate| {
        let mut calls = 0;
        assert!(selected(101_325.0).compile_with_body(cx, gate, None, |name| {
            calls += 1; Ok(pack(name, 1.0, None))
        }).is_err());
        assert_eq!(calls, 0);
        let text = include_str!("../rayleigh-pluck.material-strings");
        for selected in ["rayleigh 1 5 300 64 256 10000000 48",
            "rayleigh 1 5 3000 65 256 10000000 48",
            "rayleigh 1 5 300 65 999999 10000000 48",
            "rayleigh 3 5 300 65 256 10000000 48"] {
            assert!(input::Spec::parse(text.replace("rayleigh 1 5 300 65 256 10000000 48", selected).as_bytes()).is_err());
        }
    });
}
