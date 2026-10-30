//! The ordinary material-string file selects axial physics through to PCM.
use super::*;
#[path = "longitudinal/loss.rs"]
mod loss;
#[path = "longitudinal/bending_memory.rs"]
mod bending_memory;
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
