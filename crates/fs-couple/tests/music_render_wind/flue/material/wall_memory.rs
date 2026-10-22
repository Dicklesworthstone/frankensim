//! Original material spectra through wall motion, returning pressure and PCM.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::WallMaterialLoss;

fn input() -> String {
    source().replace("stress-free kelvin-voigt",
        "stress-free maxwell 1 0.1 0.2 boundary_modulus boundary_time")
}

#[test]
fn wall_relaxation_data_reaches_jet_feedback_and_the_existing_cli() {
    let text = input();
    let original = pack_with_time(1000.0, 1e6, None, 100_000.0, 0.001);
    let p = loaded(&text, &original, 37);
    let wall = &p.renderer().voice().network().unwrap().material_walls()[0];
    let WallMaterialLoss::Maxwell(spectrum) = &wall.source().loss else {
        panic!("the original material spectrum was not retained");
    };
    assert_eq!(spectrum.material().e_inf, 1e6);
    assert_eq!(spectrum.material().terms, vec![(5e5, 0.001)]);
    assert_eq!(wall.patch().wall.resistance, 0.0,
        "the source's Kelvin-Voigt property must not be counted a second time");
    assert_eq!(wall.load().compliance_terms().len(), 1);
    assert_eq!(p.wall_records()[0].pack(), &original);
    let mut renderer = p.into_renderer();
    let mut pressure = Vec::new();
    let mut maximum_memory = 0.0_f64;
    while renderer.remaining_samples() != 0 {
        let mut sample = [0.0];
        renderer.block(&mut sample).unwrap();
        pressure.push(sample[0]);
        let network = renderer.voice().network().unwrap();
        let memory = network.wall_relaxation_energy_j(1).unwrap();
        assert!(memory.is_finite() && memory >= 0.0);
        maximum_memory = maximum_memory.max(memory);
    }
    assert!(maximum_memory > 0.0);
    let reference = trace(loaded(&text, &original, 37));
    assert_eq!(pressure, reference.pressure);
    for block in [1, 37, 512] {
        let mut renderer = loaded(&text, &original, block).into_renderer();
        let mut actual = vec![0.0; pressure.len()];
        for chunk in actual.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        assert_eq!(actual, pressure);
    }
    // Keep the complete performance file, geometry, gas, equilibrium modulus
    // and branch modulus fixed. Change only a supplied relaxation-time claim.
    let substituted = pack_with_time(1000.0, 1e6, None, 100_000.0, 0.005);
    let changed = trace(loaded(&text, &substituted, 37));
    assert_ne!(reference.incoming, changed.incoming);
    assert_ne!(reference.pressure, changed.pressure);
    let dir = directory();
    let input = dir.join("wall-memory.performance");
    std::fs::write(&input, &text).unwrap();
    for (name, material, expected) in [
        ("original", original, reference.pressure),
        ("substituted", substituted, changed.pressure),
    ] {
        std::fs::write(dir.join("wall.fsmcdpk"), material.to_bytes()).unwrap();
        let output = dir.join(format!("wall-memory-{name}.wav"));
        success(run(&input, &output, 37, false));
        assert_eq!(std::fs::read(output).unwrap(),
            encode_pcm16_wav(&expected, 48_000, 1.0).unwrap().0);
    }
}

#[test]
fn incomplete_or_unresolved_wall_spectra_refuse_before_rendering() {
    let text = input();
    let gate = CancelGate::new_clock_free();
    for material in [
        pack(1000.0, 1e6, Some("boundary_modulus"), 100_000.0),
        pack(1000.0, 1e6, Some("boundary_time"), 100_000.0),
        pack_with_time(1000.0, 1e6, None, 100_000.0, 0.0),
        pack_with_time(1000.0, 1e6, None, 100_000.0, 1e-7),
    ] {
        assert!(FlueFilePerformance::from_bytes_with_sources(text.as_bytes(), 37, &gate,
            None, |_| Ok(material.clone())).is_err());
    }
    let material = pack(1000.0, 1e6, None, 100_000.0);
    for changed in [text.replace("maxwell 1", "maxwell 0"),
        text.replace("maxwell 1", "maxwell 9"),
        text.replace("boundary_modulus boundary_time", "boundary_modulus young_modulus"),
        text.replace("0.1 0.2", "0.1 0.0001"),
        text.replace("maxwell 1", "maxwell 2")
            .replace("boundary_modulus boundary_time", "boundary_modulus boundary_time boundary_modulus boundary_time")] {
        assert!(FlueFilePerformance::from_bytes_with_sources(changed.as_bytes(), 37, &gate,
            None, |_| Ok(material.clone())).is_err());
    }
}
