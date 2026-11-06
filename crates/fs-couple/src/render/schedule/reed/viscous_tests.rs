use super::*;

const INPUT: &str = include_str!("../../../../examples/reed-duct.performance");
const VISCOUS_INPUT: &str = include_str!("../../../../examples/reed-viscous-slit.performance");

fn source(massive: bool) -> String {
    if massive {
        VISCOUS_INPUT.to_string()
    } else {
        VISCOUS_INPUT.replace("reed 0.0004 0.013 6000 0.00001 100 0.35",
            "reed 0.0004 0.013 6000 0 0 0.35")
    }
}

fn pressure(text: &str, block: usize) -> Vec<f64> {
    let mut source = ReedPerformance::from_bytes(text.as_bytes(), block).unwrap();
    let mut output = vec![0.0; source.info().samples as usize];
    for part in output.chunks_mut(block) { source.block(part).unwrap(); }
    assert!(output.iter().all(|p| p.is_finite()));
    output
}

// G3/G5: source-file geometry chooses the constitutive owner, alters physical
// pressure, and preserves the same exact trajectory across callback boundaries.
#[test]
fn selected_slit_reaches_massless_and_massive_performances() {
    for massive in [false, true] {
        let source = source(massive);
        let loaded = ReedPerformance::from_bytes(source.as_bytes(), 37).unwrap();
        let slit = loaded.viscous_slit().unwrap();
        assert_eq!(slit.length_m(), 0.001);
        assert_eq!(slit.gas(), Some(loaded.ambient_gas()));
        let viscous = pressure(&source, 37);
        assert_eq!(viscous, pressure(&source, 113));
        let inviscid = pressure(&source.replace("slit_length 0.001\n", ""), 37);
        let energy: f64 = inviscid.iter().map(|p| p * p).sum();
        let changed: f64 = viscous.iter().zip(inviscid).map(|(a, b)| (a-b).powi(2)).sum();
        assert!(energy > 0.0 && changed > 1e-8 * energy);
    }
}

#[test]
fn regional_inlet_supplies_channel_viscosity_instead_of_default_or_exterior_air() {
    let source = source(true)
        .replace("cylinder 0.0022 0.5", "cylinder 0.0015 0.5 gas mixture 310 101325 mole 2 N2 0.5 CO2 0.5");
    let p = ReedPerformance::from_bytes(source.as_bytes(), 37).unwrap();
    let inlet = p.segment_gases()[0];
    let slit = p.viscous_slit().unwrap();
    assert_eq!(slit.gas(), Some(&inlet));
    assert_eq!(slit.dynamic_viscosity_pa_s(), inlet.dynamic_viscosity);
    assert_ne!(inlet.dynamic_viscosity, p.ambient_gas().dynamic_viscosity);
    let actual = pressure(&source, 37);
    let different = pressure(&source.replace("gas mixture 310 101325 mole 2 N2 0.5 CO2 0.5",
        "gas mixture 310 101325 mole 2 N2 0.8 CO2 0.2"), 37);
    assert_ne!(actual, different);
}

#[test]
fn axial_field_channel_uses_the_exact_inlet_instead_of_its_first_midpoint_cell() {
    use crate::driving_point::regional::field::tests::field_with_viscosity;
    // The shared fixture's viscosity claim is flat in temperature, so it can
    // never tell the inlet from the first cell by viscosity. Give it an
    // air-like slope; the slit law consumes exactly this property.
    let field = field_with_viscosity(280.0, [1.33e-5, 4.15e-5]);
    let inlet = *field.state_at(0.0).unwrap().state();
    let source = source(true)
        .replace("ambient 293.15 101325 0", "ambient 300 101325 0")
        .replace("cylinder 0.0022 0.5", "cylinder 0.0015 0.57");
    let p = ReedPerformance::from_bytes_with_gas_field(source.as_bytes(), 37,
        ObservationPolicy::Diagnostic, None, field, 8).unwrap();
    assert_eq!(p.viscous_slit().unwrap().gas(), Some(&inlet));
    assert_eq!(p.viscous_slit().unwrap().dynamic_viscosity_pa_s(), inlet.dynamic_viscosity);
    assert_ne!(p.segment_gases()[0].temperature, inlet.temperature);
    assert_ne!(p.segment_gases()[0].dynamic_viscosity, inlet.dynamic_viscosity);
}

#[test]
fn malformed_slit_geometry_refuses_before_rendering() {
    for row in ["slit_length 0", "slit_length -0.001", "slit_length NaN",
        "slit_length inf", "slit_length", "slit_length 0.001 extra",
        "slit_length 0.001\nslit_length 0.001"] {
        let source = source(false).replace("slit_length 0.001", row);
        assert!(ReedPerformance::from_bytes(source.as_bytes(), 37).is_err(), "{row}");
    }
    assert!(ReedPerformance::from_bytes(INPUT.as_bytes(), 37).unwrap().viscous_slit().is_none());
}

// ---- staple pressure recovery (Almeida, Vergez and Caussé 2007, eqs. 14-17) ----------------------------------------

/// The oboe description's reed (rest gap 357 um on 5.6 mm: 2.0 mm^2, inside the law's domain) on the slit fixture, with
/// or without the `staple_recovery` row right after `slit_length`.
fn staple_source(row: Option<&str>) -> String {
    let base = VISCOUS_INPUT.replace("reed 0.0004 0.013 6000 0.00001 100 0.35",
        "reed 0.000357 0.0056 17720 2.2e-5 3464 0.114");
    match row {
        Some(row) => base.replace("slit_length 0.001\n", &format!("slit_length 0.001\n{row}\n")),
        None => base,
    }
}

#[test]
fn staple_recovery_row_selects_the_law_changes_the_sound_and_zero_changes_nothing() {
    let plain = pressure(&staple_source(None), 37);
    let staple_text = staple_source(Some("staple_recovery 0.8 4e-6"));
    let loaded = ReedPerformance::from_bytes(staple_text.as_bytes(), 37).unwrap();
    let recovery = loaded.staple_recovery().unwrap();
    assert_eq!((recovery.c_p(), recovery.constriction_area_m2()), (0.8, 4e-6));
    let staple = pressure(&staple_text, 37);
    assert_eq!(staple, pressure(&staple_text, 113), "callback partition");
    let energy: f64 = plain.iter().map(|p| p * p).sum();
    let changed: f64 = plain.iter().zip(&staple).map(|(a, b)| (a - b).powi(2)).sum();
    assert!(energy > 0.0 && changed > 1e-8 * energy, "the recovery must change the sound");
    // C_P = 0 is the original law: nothing installed, the same bits.
    let zero_text = staple_source(Some("staple_recovery 0 4e-6"));
    assert!(ReedPerformance::from_bytes(zero_text.as_bytes(), 37).unwrap().staple_recovery().is_none());
    let zero = pressure(&zero_text, 37);
    assert_eq!(plain.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
        zero.iter().map(|p| p.to_bits()).collect::<Vec<_>>());
    assert!(ReedPerformance::from_bytes(staple_source(None).as_bytes(), 37).unwrap().staple_recovery().is_none());
}

#[test]
fn staple_recovery_row_is_refused_for_a_massless_reed_a_malformed_row_or_an_undefined_opening() {
    use crate::acoustic_realize::AcousticRealizeError;
    use crate::bernoulli_aperture::viscous::STAPLE_RECOVERY_DOMAIN;
    let refusal = |text: &str| match ReedPerformance::from_bytes(text.as_bytes(), 37) {
        Err(ReedPerformanceError::Input { what, .. }) => what.to_string(),
        Err(ReedPerformanceError::Model(AcousticRealizeError::Reed { what })) => what.to_string(),
        other => panic!("expected a refusal, got {:?}", other.err()),
    };
    let row = Some("staple_recovery 0.8 4e-6");
    assert!(refusal(&staple_source(row).replace("reed 0.000357 0.0056 17720 2.2e-5 3464 0.114",
        "reed 0.000357 0.0056 17720 0 0 0.114")).contains("massive reed"));
    for bad in ["staple_recovery", "staple_recovery 0.8", "staple_recovery 0.8 4e-6 extra",
        "staple_recovery 1.1 4e-6", "staple_recovery -0.1 4e-6", "staple_recovery 0.8 0",
        "staple_recovery 0.8 -1e-6", "staple_recovery NaN 4e-6",
        "staple_recovery 0.8 4e-6\nstaple_recovery 0.8 4e-6"] {
        assert!(ReedPerformance::from_bytes(staple_source(Some(bad)).as_bytes(), 37).is_err(), "{bad}");
    }
    // The fixture's own clarinet-sized reed (0.4 mm x 13 mm = 5.2 mm^2) is above S_c / sqrt(C_P) = 4.47 mm^2.
    let wide = VISCOUS_INPUT.replace("slit_length 0.001\n", "slit_length 0.001\nstaple_recovery 0.8 4e-6\n");
    assert_eq!(refusal(&wide), STAPLE_RECOVERY_DOMAIN);
    // The same on the inviscid law (no slit_length): the row loads there and the same opening is refused by name.
    let wide_inviscid = INPUT.replace("reed 0.0004 0.013 6000 0 0 0.35",
        "reed 0.0004 0.013 6000 0.00001 100 0.35\nstaple_recovery 0.8 4e-6");
    assert_eq!(refusal(&wide_inviscid), STAPLE_RECOVERY_DOMAIN);
}

/// The oboe description's reed on the inviscid fixture card (no `slit_length`, as the live oboe cards), with or without the
/// `staple_recovery` row right after the reed row.
fn inviscid_staple_source(row: Option<&str>) -> String {
    let reed = "reed 0.000357 0.0056 17720 2.2e-5 3464 0.114\n";
    let base = INPUT.replace("reed 0.0004 0.013 6000 0 0 0.35\n", reed);
    match row {
        Some(row) => base.replace(reed, &format!("{reed}{row}\n")),
        None => base,
    }
}

#[test]
fn staple_recovery_row_works_on_the_inviscid_law_of_a_card_without_slit_length() {
    let plain_text = inviscid_staple_source(None);
    assert!(!plain_text.contains("slit_length"));
    let plain = pressure(&plain_text, 37);
    let staple_text = inviscid_staple_source(Some("staple_recovery 0.8 4e-6"));
    let loaded = ReedPerformance::from_bytes(staple_text.as_bytes(), 37).unwrap();
    assert!(loaded.viscous_slit().is_none(), "this card has no slit_length");
    let recovery = loaded.staple_recovery().unwrap();
    assert_eq!((recovery.c_p(), recovery.constriction_area_m2()), (0.8, 4e-6));
    let staple = pressure(&staple_text, 37);
    assert_eq!(staple, pressure(&staple_text, 113), "callback partition");
    let energy: f64 = plain.iter().map(|p| p * p).sum();
    let changed: f64 = plain.iter().zip(&staple).map(|(a, b)| (a - b).powi(2)).sum();
    assert!(energy > 0.0 && changed > 1e-8 * energy, "the recovery must change the sound");
    // C_P = 0 is the original law: nothing installed, the same bits.
    let zero_text = inviscid_staple_source(Some("staple_recovery 0 4e-6"));
    assert!(ReedPerformance::from_bytes(zero_text.as_bytes(), 37).unwrap().staple_recovery().is_none());
    assert_eq!(plain.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
        pressure(&zero_text, 37).iter().map(|p| p.to_bits()).collect::<Vec<_>>());
    // With the viscous slit the same row still loads, and then selects the viscous junction.
    let slit_text = staple_source(Some("staple_recovery 0.8 4e-6"));
    assert!(ReedPerformance::from_bytes(slit_text.as_bytes(), 37).unwrap().viscous_slit().is_some());
}

#[test]
fn contraction_row_selects_a_narrower_jet_changes_the_sound_and_one_changes_nothing() {
    let plain = pressure(&inviscid_staple_source(None), 37);
    let text = inviscid_staple_source(Some("contraction 0.611"));
    let loaded = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
    assert_eq!(loaded.jet_contraction(), Some(0.611));
    assert!(loaded.staple_recovery().is_none() && loaded.viscous_slit().is_none());
    let narrow = pressure(&text, 37);
    assert_eq!(narrow, pressure(&text, 113), "callback partition");
    let energy: f64 = plain.iter().map(|p| p * p).sum();
    let changed: f64 = plain.iter().zip(&narrow).map(|(a, b)| (a - b).powi(2)).sum();
    assert!(energy > 0.0 && changed > 1e-8 * energy, "the contraction must change the sound");
    // C_c = 1 is the original law: nothing installed, the same bits.
    let one_text = inviscid_staple_source(Some("contraction 1"));
    assert!(ReedPerformance::from_bytes(one_text.as_bytes(), 37).unwrap().jet_contraction().is_none());
    assert_eq!(plain.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
        pressure(&one_text, 37).iter().map(|p| p.to_bits()).collect::<Vec<_>>());
    // After the staple row: both load, and the sound differs from either alone.
    let both_text = inviscid_staple_source(Some("staple_recovery 0.8 4e-6\ncontraction 0.611"));
    let both = ReedPerformance::from_bytes(both_text.as_bytes(), 37).unwrap();
    assert_eq!((both.jet_contraction(), both.staple_recovery().map(|r| r.c_p())), (Some(0.611), Some(0.8)));
    let staple = pressure(&inviscid_staple_source(Some("staple_recovery 0.8 4e-6")), 37);
    let together = pressure(&both_text, 37);
    for other in [&narrow, &staple] {
        let moved: f64 = together.iter().zip(other.iter()).map(|(a, b)| (a - b).powi(2)).sum();
        assert!(moved > 1e-8 * energy, "the two rows together must differ from either alone");
    }
}

#[test]
fn contraction_row_is_refused_for_a_massless_reed_a_slit_a_malformed_row_or_a_staple_outside_its_domain() {
    use crate::acoustic_realize::AcousticRealizeError;
    use crate::bernoulli_aperture::viscous::STAPLE_RECOVERY_DOMAIN;
    let refusal = |text: &str| match ReedPerformance::from_bytes(text.as_bytes(), 37) {
        Err(ReedPerformanceError::Input { what, .. }) => what.to_string(),
        Err(ReedPerformanceError::Model(AcousticRealizeError::Reed { what })) => what.to_string(),
        other => panic!("expected a refusal, got {:?}", other.err()),
    };
    let row = Some("contraction 0.611");
    assert!(refusal(&inviscid_staple_source(row).replace("reed 0.000357 0.0056 17720 2.2e-5 3464 0.114",
        "reed 0.000357 0.0056 17720 0 0 0.114")).contains("massive reed"));
    // The viscous channel and a contraction exclude each other.
    assert!(refusal(&staple_source(row)).contains("inviscid reed law only"));
    for bad in ["contraction", "contraction 0", "contraction -0.5", "contraction 1.0000001", "contraction 2",
        "contraction NaN", "contraction 0.6 extra", "contraction 0.6\ncontraction 0.6"] {
        assert!(ReedPerformance::from_bytes(inviscid_staple_source(Some(bad)).as_bytes(), 37).is_err(), "{bad}");
    }
    // The staple's domain is taken on the jet area: the fixture's 5.2 mm^2 opening is refused by name at C_c = 1, loads at
    // 0.8 (4.16 mm^2 < 4.47 mm^2) and is refused again at 0.9 (4.68 mm^2).
    let wide = |tail: &str| INPUT.replace("reed 0.0004 0.013 6000 0 0 0.35",
        &format!("reed 0.0004 0.013 6000 0.00001 100 0.35\nstaple_recovery 0.8 4e-6{tail}"));
    assert_eq!(refusal(&wide("")), STAPLE_RECOVERY_DOMAIN);
    assert!(ReedPerformance::from_bytes(wide("\ncontraction 0.8").as_bytes(), 37).is_ok());
    assert_eq!(refusal(&wide("\ncontraction 0.9")), STAPLE_RECOVERY_DOMAIN);
}
