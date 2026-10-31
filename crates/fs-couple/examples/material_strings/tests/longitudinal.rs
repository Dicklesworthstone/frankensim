//! The ordinary material-string file selects axial physics through to PCM.
use super::*;
use fs_material::section::ConcentricKelvinVoigtSection;
#[path = "longitudinal/loss.rs"]
mod loss;
#[path = "longitudinal/bending_memory.rs"]
mod bending_memory;
#[path = "longitudinal/axial_memory.rs"]
mod axial_memory;
const RECORD: &str = "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2";
fn text() -> String {
    include_str!("../longitudinal-cylinder-bow.material-strings")
        .replace("48000 4096 0.01 128", "48000 129 0.01 128")
}
#[test]
fn data_selected_axial_modes_reach_bonded_strings_cavity_gas_and_bit_stable_pcm() {
    scoped(|cx, gate| {
        let build = || {
            compile(
                cx,
                gate,
                input::Spec::parse(text().as_bytes()).unwrap(),
                7800.0,
            )
        };
        let mut direct = build();
        let source = direct.performance.bank().source().bank();
        assert!(source.geometric_stretching().is_none());
        let longitudinal = source.longitudinal_waves().unwrap();
        assert_eq!(longitudinal.strings()[0].model().axial_modes(), &[1, 2]);
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.material_cavity().is_some());
        assert_eq!(
            direct
                .performance
                .bank()
                .system()
                .components()
                .iter()
                .map(|c| c.modes().len())
                .collect::<Vec<_>>(),
            vec![3, 2, 1, 12, 2]
        );
        let expected = pressure(&mut direct, 37);
        assert!(expected.iter().any(|p| p.abs() > 1e-15));
        assert!(
            direct.performance.bank().system().components()[4].states()[1].displacement_m_sqrt_kg
                != 0.0
        );
        for block in [1, 128] {
            assert_eq!(pressure(&mut build(), block), expected);
        }
        let frame = direct.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let (output, summary) = build().write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let wav = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48000, 0.01).unwrap();
        assert_eq!(output.into_inner(), wav.0);
        assert_eq!(summary.samples, 129);
    });
}
#[test]
fn axial_material_substitution_changes_real_output_and_preserves_source_density_scaling() {
    scoped(|cx, gate| {
        let mut a = compile(
            cx,
            gate,
            input::Spec::parse(text().as_bytes()).unwrap(),
            7800.0,
        );
        let mut b = compile(
            cx,
            gate,
            input::Spec::parse(text().as_bytes()).unwrap(),
            15600.0,
        );
        let law = |ready: &Ready| {
            ready
                .performance
                .bank()
                .source()
                .bank()
                .longitudinal_waves()
                .unwrap()
                .strings()[0]
                .clone()
        };
        assert!(law(&b).model().axial_omegas()[0] < law(&a).model().axial_omegas()[0]);
        assert_eq!(law(&b).axial_rigidity_n(), law(&a).axial_rigidity_n());
        assert!(
            pressure(&mut a, 37)
                .iter()
                .zip(pressure(&mut b, 128))
                .any(|(x, y)| (x - y).abs() > 1e-18)
        );
    });
}
// Manufactured material response, not a measured alloy or moisture-uptake law.
// Both runs keep the same source artifact; only their explicit held states vary.
fn environment_pack(name: &str, band: f64, typed: bool, frequency_response: bool) -> NormalizedPack {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, MoistureBasis, QuantityKind, SemanticType, ValueForm}};
    let original = pack(name, 7800.0, band, None);
    let mut claims = ClaimSet::new();
    for id in original.claims().observation_ids() {
        assert_eq!(claims.register_observation(original.claims().observation(id).unwrap().clone()).unwrap(), id);
    }
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 280.0, 320.0)
            .with_quantity("omega", quantity(QuantityKind::Frequency(FrequencyConvention::Angular)), 0.0, band);
        if matches!(claim.key.name(), "density" | "young_modulus") {
            claim.validity = if typed {
                claim.validity.with_quantity("moisture_content",
                    quantity(QuantityKind::Moisture(MoistureBasis::DryMass)), 0.0, 0.2)
            } else { claim.validity.with("moisture_content", 0.0, 0.2) };
            let PropertyValue::Scalar { value, dims } = claim.value else { unreachable!() };
            let scale = if claim.key.name() == "density" { 1.2 } else { 0.8 };
            let (axis, axis_dims, end) = if frequency_response && claim.key.name() == "young_modulus" {
                ("omega", Dims([0, 0, -1, 0, 0, 0]), band)
            } else { ("moisture_content", Dims::NONE, 0.2) };
            claim.value = PropertyValue::Curve { abscissa: axis.into(), abscissa_dims: axis_dims,
                knots: vec![(0.0, value), (end, scale * value)], dims };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    NormalizedPack::new(name, "manufactured-axial-environment-v1", original.source_artifact(),
        "CC0-1.0", claims, vec![], vec![]).unwrap()
}

fn environment_text(core: f64, sleeve: f64) -> String {
    text().replace("state 293.15 1000", &format!(
        "state 293.15 1000\nmaterial-states 2\nmaterial-state wire.fsmatpk 1 moisture_content {core}\nmaterial-state sleeve.fsmatpk 1 moisture_content {sleeve}"
    ))
}

fn environment_section(source: &MaterialStringSource) -> &ConcentricKelvinVoigtSection {
    match source {
        MaterialStringSource::WithPolarizations(p) => environment_section(p.source()),
        MaterialStringSource::WithTorsion(t) => environment_section(t.bending()),
        MaterialStringSource::Concentric { section, .. } => section,
        _ => panic!("retain the original bonded core and sleeve"),
    }
}

#[test]
fn held_layer_states_reach_longitudinal_mass_rigidity_and_coupled_pressure() {
    scoped(|cx, gate| {
        let prepare = |raw: String, band, typed, frequency_response| {
            input::Spec::parse(raw.as_bytes()).unwrap().compile(cx, gate, |name|
                Ok(if name == "body.fsmatpk" { pack(name, 7800.0, 200_000.0, None) }
                    else { environment_pack(name, band, typed, frequency_response) }))
        };
        let mut low = prepare(environment_text(0.0, 0.0), 200_000.0, true, false).unwrap();
        let mut high = prepare(environment_text(0.2, 0.2), 200_000.0, true, false).unwrap();
        let law = |ready: &Ready| ready.performance.bank().source().bank()
            .longitudinal_waves().unwrap().strings()[0].clone();
        let a = law(&low);
        let b = law(&high);
        assert_eq!(a.resolution(), b.resolution());
        assert_eq!(a.model().axial_modes(), &[1, 2]);
        assert_eq!(a.elastic_frequency_band_rad_s(), b.elastic_frequency_band_rad_s());
        assert!((b.solid_linear_density_kg_m() / a.solid_linear_density_kg_m() - 1.2).abs() < 1e-12);
        assert!((b.axial_rigidity_n() / a.axial_rigidity_n() - 0.8).abs() < 1e-12);
        for (low_omega, high_omega) in a.model().axial_omegas().iter().zip(b.model().axial_omegas()) {
            assert!((high_omega.powi(2) / low_omega.powi(2) - 0.8 / 1.2).abs() < 1e-12);
        }
        for name in ["wire.fsmatpk", "sleeve.fsmatpk", "body.fsmatpk"] {
            assert_eq!(low.catalog[name].pack.content_hash(), high.catalog[name].pack.content_hash());
        }
        let low_layers = environment_section(&low.performance.bank().source().bank().strings()[0].source).elastic().layers();
        let high_layers = environment_section(&high.performance.bank().source().bank().strings()[0].source).elastic().layers();
        for (a, b) in low_layers.iter().zip(high_layers) {
            assert_eq!(a.material.card_identity(), b.material.card_identity());
            assert_ne!(a.material.identity(), b.material.identity());
            assert_eq!(a.material.state_coordinate("moisture_content"), Some(0.0));
            assert_eq!(b.material.state_coordinate("moisture_content"), Some(0.2));
        }
        // Each original layer owns its held state, even within one bonded string.
        let mixed = prepare(environment_text(0.0, 0.2), 200_000.0, true, false).unwrap();
        let mixed_layers = environment_section(&mixed.performance.bank().source().bank().strings()[0].source).elastic().layers();
        assert_eq!(mixed_layers.len(), 2);
        assert_eq!(mixed_layers[0].material.state_coordinate("moisture_content"), Some(0.0));
        assert_eq!(mixed_layers[1].material.state_coordinate("moisture_content"), Some(0.2));
        let mixed_law = law(&mixed);
        assert!(a.solid_linear_density_kg_m() < mixed_law.solid_linear_density_kg_m());
        assert!(mixed_law.solid_linear_density_kg_m() < b.solid_linear_density_kg_m());
        assert!(b.axial_rigidity_n() < mixed_law.axial_rigidity_n());
        assert!(mixed_law.axial_rigidity_n() < a.axial_rigidity_n());
        // The unchanged body and atmosphere are still in the original solve.
        assert_eq!(low.performance.bank().system().components()[1].modes(),
            high.performance.bank().system().components()[1].modes());
        assert_eq!(low.performance.bank().cards()[0].tension_n, high.performance.bank().cards()[0].tension_n);
        assert!(high.performance.bank().source().bank().cylinder_gas_loading().is_some());
        assert!(high.performance.bank().source().bank().material_cavity().is_some());
        let low_pressure = pressure(&mut low, 37);
        let high_pressure = pressure(&mut high, 128);
        assert!(high_pressure.iter().all(|p| p.is_finite()));
        assert!(high_pressure.iter().any(|p| p.abs() > 1e-15));
        assert!(low_pressure.iter().zip(&high_pressure).any(|(a, b)| (a - b).abs() > 1e-18));
        let partitioned = pressure(&mut prepare(environment_text(0.2, 0.2), 200_000.0, true, false).unwrap(), 1);
        assert_eq!(high_pressure, partitioned);
        assert_eq!(fs_couple::pcm_wav::encode_pcm16_wav(&high_pressure, 48000, 0.01).unwrap(),
            fs_couple::pcm_wav::encode_pcm16_wav(&partitioned, 48000, 0.01).unwrap());
        for ready in [&low, &high] {
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(prepare(text(), 200_000.0, true, false).is_err(), "no inferred layer moisture");
        assert!(prepare(environment_text(0.3, 0.0), 200_000.0, true, false).is_err(), "no extrapolation");
        assert!(prepare(environment_text(0.0, 0.0), 200_000.0, false, false).is_err(), "exact source semantics required");
        assert!(prepare(environment_text(0.0, 0.0), 200_000.0, true, true).is_err(), "frequency response is not a held environment");
        let error = prepare(environment_text(0.0, 0.0), 50_000.0, true, false).err().unwrap();
        assert!(error.to_string().contains("longitudinal frequency"), "{error}");
    });
}

#[test]
fn parser_counts_longitudinal_coordinates_and_refuses_conflicts_before_source_io() {
    let original = text();
    for replacement in [
        "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 2 1",
        "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 0",
        "longitudinal-waves 0.31 0.01 64 100000 1e-11 1e-10 1 2",
        "longitudinal-waves 0.2 0.11 64 100000 1e-11 1e-10 1 2",
    ] {
        assert!(input::Spec::parse(original.replace(RECORD, replacement).as_bytes()).is_err());
    }
    assert!(
        input::Spec::parse(
            original
                .replace(
                    RECORD,
                    &format!("geometric-stretching 0.2 64 100000 1e-11 1e-10\n{RECORD}")
                )
                .as_bytes()
        )
        .is_err()
    );
    assert!(
        input::Spec::parse(
            original
                .replace("prepared-linear 64", "prepared-linear 18")
                .as_bytes()
        )
        .is_err()
    );
    assert!(
        input::Spec::parse(
            original
                .replace(
                    "forces retain 0 0",
                    "forces preload 1 0\nport 0 0.2 surface 0.001"
                )
                .as_bytes()
        )
        .is_ok()
    );
}
#[test]
fn axial_elasticity_band_is_checked_instead_of_reusing_a_transverse_validity_check() {
    scoped(|cx, gate| {
        // Bending/torsion are admitted by this band, axial mode two is not.
        let spec = input::Spec::parse(text().as_bytes()).unwrap();
        let result = spec.compile(cx, gate, |name| Ok(pack(name, 7800.0, 50_000.0, None)));
        let error = match result {
            Ok(_) => panic!("axial source band must refuse"),
            Err(error) => error.to_string(),
        };
        assert!(error.contains("longitudinal frequency"), "{error}");
    });
}
