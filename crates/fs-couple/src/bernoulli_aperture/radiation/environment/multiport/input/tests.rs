use super::*;

mod thermochemical;

const SOURCE: &str = include_str!("../../../../../../examples/material_modal_valve/mutual.radiation");

fn decode(text: &str) -> Result<MimoRadiationArtifact, MimoRadiationInputError> {
    MimoRadiationArtifact::from_bytes_with_sources(text.as_bytes(), None, &CancelGate::new_clock_free())
}

#[test]
fn original_coefficients_gas_aperture_order_and_distinct_source_identities_survive_import() {
    let artifact = decode(SOURCE).unwrap();
    let model = artifact.source();
    assert_eq!(model.declaration_identity, [0x11; 32]);
    assert_eq!(model.realization_identity, [0x22; 32]);
    assert_ne!(artifact.content_hash().0, model.declaration_identity);
    assert_ne!(artifact.content_hash().0, model.realization_identity);
    assert_eq!(model.exterior, *artifact.gas_record().state());
    assert_eq!(model.exterior.pressure, 101325.0);
    assert_eq!(model.apertures[1].center_m, [0.03, 0.0]);
    assert_eq!(model.apertures[1].radius_m.to_bits(), 0.007_f64.to_bits());
    assert_eq!(model.impedance.ports, 2);
    assert_eq!(model.impedance.d, vec![2000.0, 300.0, 300.0, 2000.0]);
    assert_eq!(model.impedance.residues[0][1], C64::new(1_000_000.0, 0.0));
    assert_eq!(model.impedance.poles, vec![PoleTerm::Real { pole: -7000.0, residue: 0.0 }]);
    let differently_encoded = decode(&SOURCE.replace("direct 2000 ", "direct 2000.0 ")).unwrap();
    assert_eq!(model.impedance, differently_encoded.source().impedance);
    assert_eq!(model.realization_identity, differently_encoded.source().realization_identity);
    assert_ne!(artifact.content_hash(), differently_encoded.content_hash());
}

#[test]
fn shared_gas_decoder_retains_an_independently_authored_mixture() {
    let text = SOURCE.replace("exterior 293.15 101325 0",
        "exterior mixture 303.15 65000 mole 3 N2 0.70 O2 0.20 He 0.10");
    let artifact = decode(&text).unwrap();
    assert!(matches!(artifact.gas_record(), GasRecord::Mixture(_)));
    assert_eq!(artifact.source().exterior.pressure, 65000.0);
    assert_eq!(artifact.source().exterior.temperature, 303.15);
    assert_ne!(artifact.source().exterior, decode(SOURCE).unwrap().source().exterior);
    let missing = SOURCE.replace("exterior 293.15 101325 0",
        "exterior thermochemical 293.15 101325 mole 1 nitrogen 1");
    assert!(matches!(decode(&missing), Err(MimoRadiationInputError::Gas { .. })));
}

#[test]
fn complex_pair_matrix_is_retained_without_conjugating_or_refitting() {
    let text = SOURCE.replace("pole real -7000", "pole pair -7000 9000")
        .replace("residue 4000000 0 1000000 0 1000000 0 4000000 0",
            "residue 4000000 -200 1000000 20 1000000 20 4000000 -200");
    let artifact = decode(&text).unwrap();
    assert_eq!(artifact.source().impedance.poles, vec![PoleTerm::Pair {
        pole: C64::new(-7000.0, 9000.0), residue: C64::ZERO,
    }]);
    assert_eq!(artifact.source().impedance.residues[0][1], C64::new(1_000_000.0, 20.0));
}

#[test]
fn malformed_nonfinite_and_oversized_shapes_refuse_before_any_runtime_is_built() {
    for text in [
        SOURCE.replace("apertures 2", "apertures 129"),
        SOURCE.replace("apertures 2", "apertures 0"),
        SOURCE.replace("poles 1", "poles 257"),
        SOURCE.replace("0.007", "NaN"),
        SOURCE.replace("-7000", "-inf"),
        SOURCE.replace("direct 2000 300 300 2000", "direct 2000 300 2000"),
        SOURCE.replace("improper 0 0 0 0", "improper 0 0 0 0 0"),
        SOURCE.replace("reciprocal true", "reciprocal maybe"),
        SOURCE.replace("band 1 1000000", "band 1000000 1"),
        SOURCE.replace("\n", "\r\n"),
        SOURCE.trim_end().into(),
        format!("{SOURCE}direct 1 2 2 1\n"),
    ] {
        assert!(decode(&text).is_err(), "unexpected admission: {text}");
    }
    assert!(MimoRadiationArtifact::from_bytes_with_sources(
        &vec![b'\n'; MAX_MIMO_RADIATION_BYTES + 1], None, &CancelGate::new_clock_free()).is_err());
    // Repeated matrix terms may fit the byte budget but exceed the actual
    // retained real-state dimension; refuse before allocating a 129th matrix.
    let mut text = SOURCE.split("poles 1\n").next().unwrap().to_string();
    text.push_str("poles 129\n");
    for _ in 0..129 { text.push_str("pole real -7000\nresidue 4 0 1 0 1 0 4 0\n"); }
    text.push_str("direct 2000 300 300 2000\nimproper 0 0 0 0\n");
    assert!(decode(&text).unwrap_err().to_string().contains("state dimension"));
}

#[test]
fn source_selector_is_bounded_and_loader_cancellation_publishes_nothing() {
    let input = MimoRadiationInput::parse("radiation/source.txt 256 8388608 65").unwrap();
    let mut calls = 0;
    let (model, record) = input.load(None, &mut |path| {
        assert_eq!(path, "radiation/source.txt"); calls += 1; Ok(SOURCE.as_bytes().to_vec())
    }, &CancelGate::new_clock_free()).unwrap();
    assert_eq!(calls, 1);
    assert_eq!(record.path(), input.path());
    assert_eq!(record.apertures(), model.apertures);
    assert_eq!(*record.gas_record().state(), model.exterior);
    assert_eq!(record.artifact_hash(), decode(SOURCE).unwrap().content_hash());
    assert_eq!(record.limits().sampled_passivity_points, 65);
    for selection in ["../a 256 8388608 65", "/a 256 8388608 65",
        "a//b 256 8388608 65", "a 257 8388608 65", "a 256 0 65",
        "a 256 8388608 16", "a 256 8388608 65538", "a 256 8388608 65 extra"] {
        assert!(MimoRadiationInput::parse(selection).is_err());
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(input.load(None, &mut |_| panic!("cancelled before source I/O"), &gate),
        Err(MimoRadiationInputError::Cancelled)));
    let gate = CancelGate::new_clock_free();
    assert!(matches!(input.load(None, &mut |_| {
        gate.request(); Ok(SOURCE.as_bytes().to_vec())
    }, &gate), Err(MimoRadiationInputError::Cancelled)));
}
