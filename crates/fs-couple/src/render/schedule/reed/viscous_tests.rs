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
