//! Public V3 data -> original material resolution -> two-plane physics -> PCM.
use super::*;

fn v3() -> input::Spec {
    let bytes = include_str!("../polarized-pluck.material-strings")
        .replace("48000 4096 0.01 128", "48000 256 0.01 128");
    input::Spec::parse(bytes.as_bytes()).unwrap()
}

#[test]
fn two_plane_material_input_reaches_physical_pressure_and_streamed_pcm() {
    scoped(|cx, gate| {
        let mut direct = compile(cx, gate, v3(), 7800.0);
        let source = &direct.performance.bank().source().bank().strings()[0].source;
        let MaterialStringSource::WithPolarizations(p) = source else {
            panic!("lost two-plane input")
        };
        assert_eq!(p.directions().bow(), [0.6, 0.8]);
        assert_eq!(p.directions().bridge(), [1.0, 0.0]);
        let MaterialStringSource::WithTorsion(t) = p.source() else {
            panic!("lost torsion input")
        };
        assert!(matches!(
            t.bending(),
            MaterialStringSource::Concentric { .. }
        ));
        assert_eq!(direct.catalog.len(), 3);
        let component = &direct.performance.bank().system().components()[0];
        assert_eq!(component.modes().len(), 11);
        assert_eq!(component.modes()[..4], component.modes()[4..8]);
        assert!(
            component.states()[4..8]
                .iter()
                .any(|s| s.displacement_m_sqrt_kg != 0.0)
        );
        let expected = pressure(&mut direct, 37);
        assert!(expected.iter().any(|p| p.abs() > 1e-10));
        let (bytes, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48000, 0.01).unwrap();
        let mut output = compile(cx, gate, v3(), 7800.0);
        let (wav, summary) = output.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), bytes);
        assert_eq!(summary.clipped_samples, clips as u64);
        assert_eq!(summary.samples, 256);
    });
}

#[test]
fn rotating_source_directions_preserves_sound_but_material_and_physical_angle_changes_do_not() {
    scoped(|cx, gate| {
        let mut original = compile(cx, gate, v3(), 7800.0);
        let text = String::from_utf8(v3().bytes).unwrap();
        let rotated = text.replace("two 0.6 0.8 1 0", "two -0.28 0.96 0.6 0.8");
        let mut rotated = compile(
            cx,
            gate,
            input::Spec::parse(rotated.as_bytes()).unwrap(),
            7800.0,
        );
        let changed = text.replace("two 0.6 0.8 1 0", "two 1 0 1 0");
        let mut angle = compile(
            cx,
            gate,
            input::Spec::parse(changed.as_bytes()).unwrap(),
            7800.0,
        );
        let mut material = compile(cx, gate, v3(), 15_600.0);
        assert!(
            material.performance.bank().cards()[0].fundamental_hz()
                < original.performance.bank().cards()[0].fundamental_hz()
        );
        let modes = material.performance.bank().system().components()[0].modes();
        assert_eq!(modes[..4], modes[4..8]);
        let x = pressure(&mut original, 128);
        let y = pressure(&mut rotated, 37);
        let z = pressure(&mut angle, 1);
        let w = pressure(&mut material, 128);
        let norm: f64 = x.iter().map(|x| x * x).sum();
        let rotation_error: f64 = x.iter().zip(&y).map(|(x, y)| (x - y) * (x - y)).sum();
        assert!(rotation_error < 1e-16 * norm.max(1e-30));
        assert!(x.iter().zip(&z).any(|(x, z)| (x - z).abs() > 1e-12));
        assert!(x.iter().zip(&w).any(|(x, w)| (x - w).abs() > 1e-12));
    });
}

#[test]
fn v3_requires_explicit_directions_and_counts_both_planes_against_the_budget() {
    let text = include_str!("../polarized-pluck.material-strings");
    for (from, to) in [
        ("polarizations two 0.6 0.8 1 0", "polarizations two 0 0 1 0"),
        (
            "polarizations two 0.6 0.8 1 0",
            "polarizations two 0.6 0.8 2 0",
        ),
        (
            "polarizations two 0.6 0.8 1 0",
            "polarizations two NaN 0 1 0",
        ),
        ("polarizations two 0.6 0.8 1 0\n", ""),
    ] {
        assert!(
            input::Spec::parse(text.replace(from, to).as_bytes()).is_err(),
            "accepted {to}"
        );
    }
    // Both planes count. The first string count whose two-plane strings at the
    // per-field bound exceed the complete mode ceiling is derived from the
    // crate's own constants, not typed in; its one-plane twin is admitted.
    let two = "polarizations two 0.6 0.8 1 0";
    let block = "string wire.fsmatpk 0.5 0.0003 extension 0.4998 0.001 4 3 0.8 3000 0.1\n\
        section bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003\npolarizations two 0.6 0.8 1 0\n";
    let full = block.replace(
        "0.001 4 3 0.8",
        &format!("0.001 {} 3 0.8", input::MAX_MODES),
    );
    let strings = |k: usize, planes: &str| {
        text.replace(
            &format!("strings 1\n{block}"),
            &format!("strings {k}\n{}", full.repeat(k)),
        )
        .replace(two, planes)
    };
    let modes = |k: usize, planes: &str| {
        input::Spec::parse(strings(k, planes).as_bytes())
            .unwrap()
            .network_dimensions(false)
            .unwrap()
            .0
    };
    let per_string = modes(2, two) - modes(1, two);
    let base = modes(1, two) - per_string;
    let k = (2..=8)
        .find(|&k| base + k * per_string > input::MAX_TOTAL_MODES)
        .expect("the 8-string record must reach the complete mode ceiling on two planes");
    assert_eq!(modes(k - 1, two), base + (k - 1) * per_string);
    match input::Spec::parse(strings(k, two).as_bytes()) {
        Ok(_) => panic!(
            "accepted {k} two-plane strings at {} modes",
            base + k * per_string
        ),
        Err(error) => assert!(
            error.to_string().contains("including gas memory, exceeds"),
            "{error}"
        ),
    }
    assert_eq!(
        modes(k, "polarizations single"),
        base + k * (per_string - input::MAX_MODES)
    );
    assert!(
        input::Spec::parse(
            text.replace("polarizations two 0.6 0.8 1 0", "polarizations single")
                .as_bytes()
        )
        .is_ok()
    );
}
