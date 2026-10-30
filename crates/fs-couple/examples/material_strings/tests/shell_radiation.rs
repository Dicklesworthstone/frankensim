//! The shipped closed geometry selects the original shell/air coupling.
use super::*;

#[test]
fn closed_shell_source_retains_geometry_and_refuses_a_planar_substitute_before_pack_io() {
    let spec =
        input::Spec::parse(include_bytes!("../closed-shell-pluck.material-strings")).unwrap();
    let Some(input::BodyRadiation::ShellBem(options)) = spec.body_radiation else {
        panic!("closed shell source selected a different acoustic model");
    };
    assert_eq!(options.maximum_surface_panels, 80);
    let profile = body::Profile::parse(include_bytes!("../closed-shell.body")).unwrap();
    assert!(profile.has_shell_material_basis());
    profile.preflight(&spec).unwrap();
    scoped(|cx, gate| {
        for profile in [
            None,
            Some(body::Profile::parse(include_bytes!("../orthotropic.body")).unwrap()),
        ] {
            let mut reads = 0;
            assert!(
                spec.clone()
                    .compile_with_body::<NormalizedPack>(cx, gate, profile, |_| {
                        reads += 1;
                        Err("unexpected material read".into())
                    })
                    .is_err()
            );
            assert_eq!(reads, 0, "incompatible body must fail before material I/O");
        }
    });
}

#[test]
fn loaded_analysis_retains_shell_bem_while_prepared_playback_still_refuses() {
    use fs_math::c64::C64;
    let text = include_str!("../closed-shell-pluck.material-strings").replace(
        "body-radiation bem",
        "prepared-linear 64 67108864 200000000000\nbody-radiation bem",
    );
    let spec = input::Spec::parse(text.as_bytes()).unwrap();
    let budget = spec.prepared_linear.unwrap();
    scoped(|cx, gate| {
        let mut reads = 0;
        let refused = spec.clone().compile_with_body::<NormalizedPack>(
            cx, gate,
            Some(body::Profile::parse(include_bytes!("../closed-shell.body")).unwrap()),
            |_| { reads += 1; Err("unexpected material read".into()) },
        );
        assert_eq!(refused.err().unwrap().to_string(),
            "prepared playback admits computed body radiation only through a skin load");
        assert_eq!(reads, 0, "playback compatibility must refuse before material I/O");

        let bound = spec.compile_bank_with_profiles(
            cx, gate,
            Some(body::Profile::parse(include_bytes!("../closed-shell.body")).unwrap()),
            None, true,
            |name| Ok(if name == "closed-shell.fsmatpk" {
                super::shell_observer::shell_pack(name)
            } else { pack(name, 7800.0, 200_000.0, None) }),
        ).unwrap();
        let source = bound.bank.source().bank();
        assert_eq!(source.preparation(), ModalPreparation::HarmonicOnly(budget));
        let air = source.body_radiation().unwrap();
        let report = air.material_shell_bem().unwrap();
        assert_eq!(report.surface_panels, 80);
        assert_eq!(report.gas, bound.spec.gas);
        assert_eq!(report.body_modes, bound.spec.body.modes);
        assert!(!air.model().poles.is_empty());
        let before: Vec<_> = bound.bank.system().components().iter()
            .map(|model| model.states().to_vec()).collect();
        let response = bound.bank.material_shell_bridge_harmonic_response(
            core::f64::consts::TAU * 60.0, 0, C64::from_re(0.001),
            100_000_000, gate,
        ).unwrap();
        assert_eq!(response.released.dynamics.radiation_state.len(),
            2 * air.model().poles.len());
        assert!(response.released.dynamics.radiation_dissipation_w.is_finite());
        assert!(response.released.dynamics.radiation_dissipation_w >= 0.0);
        assert!(response.released.observer_pressure_pa.re.is_finite());
        assert!(response.released.observer_pressure_pa.im.is_finite());
        assert!(bound.bank.material_shell_bridge_harmonic_response(
            core::f64::consts::TAU * 121.0, 0, C64::from_re(0.001),
            100_000_000, gate,
        ).is_err(), "analysis must retain the fitted air source band");
        assert_eq!(before, bound.bank.system().components().iter()
            .map(|model| model.states().to_vec()).collect::<Vec<_>>());
        assert_eq!(bound.bank.system().samples_rendered(), 0);
    });
}
