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

#[test]
fn moving_normal_contact_and_material_stretching_reach_the_original_pressure_writer() {
    scoped(|cx, gate| {
        let text = include_str!("../nonlinear-contact-bow.material-strings")
            .replace("48000 4096", "48000 192");
        let stops = include_str!("../contact-bow.stops").replace("knot 480", "knot 32")
            .replace("knot 960", "knot 64").replace("knot 1440", "knot 128");
        let build = |source: &str, density| input::Spec::parse(source.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None, Some(crate::stops::Profile::parse(stops.as_bytes()).unwrap()),
            |name| Ok(pack(name, if name == "wire.fsmatpk" { density } else { 7800.0 }, 200_000.0, None)),
        ).unwrap();
        let mut direct = build(&text, 7800.0);
        let bank = direct.performance.bank();
        assert!(bank.source().bank().geometric_stretching().is_some());
        assert!(bank.source().bank().material_cavity().is_some());
        assert!(bank.source().bank().cylinder_gas_loading().is_some());
        assert_eq!(bank.source().bow_contact_loads().len(), 1);
        assert_eq!(bank.system().components().iter().map(|c| c.modes().len()).collect::<Vec<_>>(), vec![3, 2, 1, 8]);
        let (mut work, mut loss, mut old_offset) = (0.0, 0.0, 0.0);
        let (mut active, mut separated, mut stretch) = (0, 0, 0.0_f64);
        let mut expected = Vec::new();
        for _ in 0..192 {
            let mut output = [0.0]; direct.performance.block(&mut output).unwrap(); expected.push(output[0]);
            let system = direct.performance.bank().system();
            let frame = system.last_frame().unwrap();
            let normal = system.last_obstacle_frames().unwrap()[0];
            assert_eq!(frame.ports[0].normal_force_n, normal.normal_force_n);
            assert_eq!(normal.actuator_work_j, -normal.normal_force_n * (normal.support_offset_m-old_offset));
            old_offset = normal.support_offset_m;
            if normal.normal_force_n == 0.0 { separated += 1; assert_eq!(frame.ports[0].traction_n, 0.0); }
            else { active += 1; }
            work += frame.external_work_j + frame.drive_work_j + frame.obstacle_work_j;
            loss += frame.network_dissipation_j + frame.friction_dissipation_j + frame.obstacle_dissipation_j;
            assert!((system.total_energy_j().unwrap() + loss - work).abs() < 1e-10);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            stretch = stretch.max(frame.nonlinear.unwrap().nonlinear_energy_j);
        }
        assert!(active > 0 && separated > 0 && stretch > 0.0);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        assert_eq!(pressure(&mut build(&text, 7800.0), 37), expected);
        let mut encoded = build(&text, 7800.0);
        let (output, summary) = encoded.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let wav = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48_000, 0.01).unwrap();
        assert_eq!(output.into_inner(), wav.0);
        assert_eq!(summary.samples, 192);
        assert_eq!(summary.clipped_samples, wav.1 as u64);
        for (source, density) in [
            (text.clone(), 15_600.0),
            (text.replace("mole 2 N2 0.8 O2 0.2", "mole 3 N2 0.7 O2 0.2 CO2 0.1"), 7800.0),
        ] {
            let actual = pressure(&mut build(&source, density), 128);
            assert!(expected.iter().zip(actual).any(|(a, b)| (a-b).abs() > 1e-14));
        }
    });
}
