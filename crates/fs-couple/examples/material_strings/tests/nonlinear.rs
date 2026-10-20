//! Performance data reaches shared geometric storage and the original writer.
use super::*;
fn text() -> String {
    include_str!("../nonlinear-cylinder-bow.material-strings")
        .replace("48000 4096 0.01 128", "48000 129 0.01 128")
}
fn input() -> input::Spec {
    input::Spec::parse(text().as_bytes()).unwrap()
}

#[test]
fn source_selected_stretching_reaches_causal_gas_body_bow_and_existing_pcm_writer() {
    scoped(|cx, gate| {
        let mut direct = compile(cx, gate, input(), 7800.0);
        let source = direct.performance.bank().source().bank();
        let nonlinear = source.geometric_stretching().unwrap();
        assert_eq!(nonlinear.potentials().len(), 1);
        assert_eq!(nonlinear.potentials()[0].transverse_planes(), 2);
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.material_cavity().is_some());
        let layout = direct
            .performance
            .bank()
            .system()
            .components()
            .iter()
            .map(|c| c.modes().len())
            .collect::<Vec<_>>();
        assert_eq!(layout, vec![3, 2, 1, 8]);
        let expected = pressure(&mut direct, 37);
        assert!(expected.iter().any(|p| p.abs() > 1e-15));
        for block in [1, 128] {
            let mut replay = compile(cx, gate, input(), 7800.0);
            assert_eq!(pressure(&mut replay, block), expected);
        }
        let frame = direct.performance.bank().system().last_frame().unwrap();
        assert!(frame.network_energy_j > 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let mut encoded = compile(cx, gate, input(), 7800.0);
        let (output, summary) = encoded.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let expected_wav = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48000, 0.01).unwrap();
        assert_eq!(output.into_inner(), expected_wav.0);
        assert_eq!(summary.samples, 129);
        assert_eq!(summary.clipped_samples, expected_wav.1 as u64);
    });
}
#[test]
fn data_only_material_gas_and_geometric_selection_change_real_pressure() {
    scoped(|cx, gate| {
        let mut nonlinear = compile(cx, gate, input(), 7800.0);
        let mut linear = compile(
            cx,
            gate,
            input::Spec::parse(
                text()
                    .replace("geometric-stretching 0.2 64 100000 1e-11 1e-10\n", "")
                    .as_bytes(),
            )
            .unwrap(),
            7800.0,
        );
        let mut dense = compile(cx, gate, input(), 15600.0);
        let mut gas = compile(
            cx,
            gate,
            input::Spec::parse(
                text()
                    .replace("mole 2 N2 0.8 O2 0.2", "mole 3 N2 0.7 O2 0.2 CO2 0.1")
                    .as_bytes(),
            )
            .unwrap(),
            7800.0,
        );
        assert_eq!(
            nonlinear.performance.bank().cards()[0].fundamental_hz(),
            linear.performance.bank().cards()[0].fundamental_hz()
        );
        let reference = pressure(&mut nonlinear, 128);
        for changed in [&mut linear, &mut dense, &mut gas] {
            let actual = pressure(changed, 37);
            assert!(
                reference
                    .iter()
                    .zip(actual)
                    .any(|(a, b)| (a - b).abs() > 1e-18)
            );
            assert_ne!(
                nonlinear.performance.bank().system().components()[1].states(),
                changed.performance.bank().system().components()[1].states()
            );
        }
    });
}
#[test]
fn impossible_geometric_inputs_and_linear_preload_refuse_before_loading_materials() {
    let original = text();
    for replacement in [
        "geometric-stretching 0.31 64 100000 1e-11 1e-10",
        "geometric-stretching 0.2 0 100000 1e-11 1e-10",
        "geometric-stretching 0.2 129 100000 1e-11 1e-10",
        "geometric-stretching 0.2 64 0 1e-11 1e-10",
        "geometric-stretching 0.2 64 100000 1e-11 1",
    ] {
        assert!(
            input::Spec::parse(
                original
                    .replace(
                        "geometric-stretching 0.2 64 100000 1e-11 1e-10",
                        replacement
                    )
                    .as_bytes()
            )
            .is_err()
        );
    }
    let preload = original.replace(
        "forces retain 0 0",
        "forces preload 1 0\nport 0 0.2 surface 0.001",
    );
    assert!(input::Spec::parse(preload.as_bytes()).is_err());
}

#[test]
fn canonical_double_bow_stretching_reaches_shared_cavity_fluid_and_pcm_without_a_new_format() {
    scoped(|cx, gate| {
        let text = include_str!("../nonlinear-double-bow.material-strings")
            .replace("48000 4096 0.01 128", "48000 129 0.01 128");
        let build = || compile(cx, gate, input::Spec::parse(text.as_bytes()).unwrap(), 7800.0);
        let mut direct = build();
        let source = direct.performance.bank().source();
        assert_eq!(source.bows().len(), 2);
        assert_eq!(source.bank().geometric_stretching().unwrap().potentials().len(), 2);
        assert!(source.bank().material_cavity().is_some());
        assert_eq!(source.bank().cylinder_gas_loading().unwrap().strings().len(), 2);
        assert_eq!(direct.performance.bank().system().components().iter()
            .map(|c| c.modes().len()).collect::<Vec<_>>(), vec![3, 3, 2, 1, 8, 8]);
        let expected = pressure(&mut direct, 37);
        assert!(expected.iter().any(|p| p.abs() > 1e-15));
        for block in [1, 128] { assert_eq!(pressure(&mut build(), block), expected); }
        let frame = direct.performance.bank().system().last_frame().unwrap();
        assert!(frame.nonlinear.unwrap().nonlinear_energy_j > 0.0);
        assert!(frame.ports.iter().all(|p| p.drive_work_j != 0.0));
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let mut encoded = build();
        let (output, summary) = encoded.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let wav = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48_000, 0.01).unwrap();
        assert_eq!(output.into_inner(), wav.0);
        assert_eq!(summary.samples, 129);
        assert_eq!(summary.clipped_samples, wav.1 as u64);
        for changed in [
            text.replace("geometric-stretching 0.2 64 100000 1e-11 1e-10\n", ""),
            text.replace("mole 2 N2 0.8 O2 0.2", "mole 3 N2 0.7 O2 0.2 CO2 0.1"),
        ] {
            let mut run = compile(cx, gate, input::Spec::parse(changed.as_bytes()).unwrap(), 7800.0);
            assert!(expected.iter().zip(pressure(&mut run, 37)).any(|(a, b)| (a - b).abs() > 1e-18));
        }
    });
}
