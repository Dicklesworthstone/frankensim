//! Exterior composition reaches the actual wind and ensemble executables.
use super::*;

const REGIONAL: &str = include_str!("../../../examples/plate-valve-regional-gas.performance");
const EXTERIOR: &str = include_str!("../../../examples/plate-valve-exterior.performance");

#[test]
fn exterior_and_common_pressure_reach_cli_pcm_with_original_internal_gases() {
    // The shipped comparison changes only the exterior declaration.
    assert_eq!(EXTERIOR.lines().filter(|line| !line.starts_with("exterior "))
        .collect::<Vec<_>>(), REGIONAL.lines().collect::<Vec<_>>());
    // Hold mole fractions, not RH, fixed for the pressure-only comparison.
    let fixed = EXTERIOR.replace("ambient 293.15 101325 0",
        "ambient mixture 293.15 101325 mole 2 N2 0.7 O2 0.3")
        .replace("gas 313.15 0.5", "gas mixture 313.15 101325 mole 3 N2 0.79 O2 0.2 H2O 0.01")
        .replace("gas 283.15 0", "gas mixture 283.15 101325 mole 2 N2 0.8 O2 0.2");
    let low = fixed.replace("101325", "50662.5");
    let high_source = load(&fixed); let low_source = load(&low);
    let (CoupledAperture::Network(high), CoupledAperture::Network(low_network)) =
        (high_source.renderer().system(), low_source.renderer().system()) else { panic!("regional sources") };
    for (a, b) in high.section_gases().unwrap().iter().zip(low_network.section_gases().unwrap()) {
        assert_eq!(a.sound_speed, b.sound_speed);
        assert_eq!(a.density, 2.0 * b.density);
    }
    let dir = directory(); let mut waves = Vec::new();
    for (index, text) in [REGIONAL, EXTERIOR, &fixed, &low].into_iter().enumerate() {
        let input = dir.join(format!("exterior-{index}.performance"));
        std::fs::write(&input, text).unwrap();
        let p = load(text);
        let (expected, _) = encode_pcm16_wav(&reference(text), 48000, p.info().full_scale_pa).unwrap();
        assert!(expected[44..].iter().any(|b| *b != 0));
        for block in [37, 113] {
            let output = dir.join(format!("exterior-{index}-{block}.wav"));
            success(run(&input, &output, block, true));
            assert_eq!(std::fs::read(&output).unwrap(), expected);
            let meta = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            let receiver = meta.split_once("\"outlet_receiver\":{").unwrap().1.split_once('}').unwrap().0;
            let medium = p.renderer().baffled_medium().unwrap();
            assert!(receiver.contains(&format!("\"density_kg_m3\":{:e}", medium.density)));
            assert!(receiver.contains(&format!("\"sound_speed_m_s\":{:e}", medium.sound_speed)));
            assert_eq!(meta.contains("\"exterior_gas\":{") , p.exterior_gas().is_some());
        }
        waves.push(expected);
    }
    assert_ne!(waves[0], waves[1], "exterior composition must reach PCM, not only metadata");
    assert_ne!(waves[2], waves[3], "common static pressure must reach sound without retuning");
}

#[test]
fn exterior_ensemble_sums_the_same_pressure_without_source_scale_gains() {
    let dir = directory();
    let (a, b) = (dir.join("inside.performance"), dir.join("outside.performance"));
    std::fs::write(&a, REGIONAL).unwrap();
    let mut combined = reference(REGIONAL);
    for (p, q) in combined.iter_mut().zip(reference(EXTERIOR)) { *p += q; }
    let (expected, _) = encode_pcm16_wav(&combined, 48000, 0.02).unwrap();
    for (index, text) in [EXTERIOR.to_string(), EXTERIOR.replace("9602 0.01", "9602 100")].iter().enumerate() {
        std::fs::write(&b, text).unwrap();
        let output = dir.join(format!("exterior-ensemble-{index}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--valve").arg(&a).arg("--valve").arg(&b)
            .args(["--decimate", "--block", "37", "--full-scale-pa", "0.02"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let meta = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert_eq!(meta.matches("\"exterior_gas\":{").count(), 1);
        assert!(meta.contains("\"source_pcm_scale_applied\":false"));
    }
}

#[test]
fn exterior_refusals_create_neither_audio_nor_sidecar() {
    let dir = directory(); let input = dir.join("refused.performance");
    let declaration = "exterior mixture 293.15 101325 mole 3 N2 0.4 O2 0.2 CO2 0.4";
    for (index, text) in [
        EXTERIOR.replace(declaration, "exterior mixture 300 101325 mole 1 He 1"),
        EXTERIOR.replace(declaration, "exterior mixture 300 101325 mole 1 H2O 1"),
        EXTERIOR.replace(declaration, "exterior 293.15 50000 0"),
        EXTERIOR.replace(declaration, &format!("{declaration}\n{declaration}")),
        EXTERIOR.replace("duct_node baffled-low-ka 1000", "duct_node reflection -0.9"),
    ].iter().enumerate() {
        std::fs::write(&input, text).unwrap();
        let output = dir.join(format!("refused-{index}.wav"));
        assert!(!run(&input, &output, 37, true).status.success());
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}
