//! Explicit prepared resources enable the actual >42-coordinate fluid model.
use super::*;

const PREPARED: &str = "prepared-linear 96 134217728 1000000000000";

fn text() -> String {
    include_str!("../prepared-cylinder-flow.material-strings")
        .replace("48000 4096 0.01 128", "48000 32 0.01 128")
}

#[test]
fn prepared_file_preserves_every_fluid_coordinate_and_explicit_budget_through_preload_and_pressure() {
    scoped(|cx, gate| {
        let bytes = text();
        let input = input::Spec::parse(bytes.as_bytes()).unwrap();
        let budget = input.prepared_linear.unwrap();
        assert_eq!(budget, PreparedModalBudget {
            max_modes: 96, max_dense_bytes: 134_217_728, max_multiply_terms: 1_000_000_000_000,
        });
        assert_eq!(input.network_dimensions(false).unwrap(), (78, 75, 78*76*76));

        // Both transverse planes participate in one physical body projection.
        // Removing only the explicit backend selection retains the original
        // connected-island refusal, rather than silently dropping fluid modes.
        let omitted = bytes.replace(&format!("{PREPARED}\n"), "");
        let error = input::Spec::parse(omitted.as_bytes()).unwrap()
            .compile(cx, gate, |name| Ok(pack(name, 7800.0, 200_000.0, None)))
            .err().expect("the ordinary exact block is too small").to_string();
        assert!(error.contains("dense exponential block capacity"), "{error}");

        // Invalid resource ceilings refuse during cold setup, before a sample
        // or WAV is published. They do not fall back to sampled air resistance.
        for (from, to, reason) in [
            (PREPARED, "prepared-linear 96 1 1000000000000", "max_dense_bytes"),
            (PREPARED, "prepared-linear 96 134217728 1", "max_multiply_terms"),
        ] {
            let error = input::Spec::parse(bytes.replace(from, to).as_bytes()).unwrap()
                .compile(cx, gate, |name| Ok(pack(name, 7800.0, 200_000.0, None)))
                .err().expect("insufficient prepared resources must refuse").to_string();
            assert!(error.contains(reason), "{error}");
        }

        // One prepared construction and a short actual release are enough to
        // exercise the path; its costly cold matrix work is not repeated here.
        let mut ready = compile(cx, gate, input, 7800.0);
        let source = ready.performance.bank().source().bank();
        assert_eq!(source.prepared_linear_budget(), Some(budget));
        assert_eq!(source.cylinder_gas_loading().unwrap().strings()[0].flow().branches().len(), 24);
        assert_eq!(source.cylinder_gas_loading().unwrap().strings()[0].flow().gas(), source.gas());
        assert_eq!(source.material_cavity().unwrap().cavity.gas(), source.gas());
        let MaterialStringSource::WithPolarizations(planes) = &source.strings()[0].source else {
            panic!("original physical directions disappeared");
        };
        assert_eq!(planes.directions().bridge(), [0.6, 0.8]);
        assert!(matches!(planes.source(), MaterialStringSource::WithTorsion(_)));
        assert_eq!(ready.performance.bank().system().components().iter()
            .map(|c| c.modes().len()).collect::<Vec<_>>(), vec![3, 2, 1, 72]);
        assert!(ready.performance.preload_energy_j().unwrap() > 0.0);
        assert!(ready.performance.bank().system().components()[3].states().iter()
            .all(|s| s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0));
        let p = pressure(&mut ready, 7);
        assert!(p.iter().any(|p| p.abs() > 1e-14));
        let components = ready.performance.bank().system().components();
        assert!(components[1].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        for plane in components[3].states().chunks_exact(24) {
            assert!(plane.iter().any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        }
        let frame = ready.performance.bank().system().last_frame().unwrap();
        assert_eq!(frame.external_work_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert_eq!(ready.performance.samples_rendered(), 32);
    });
}

#[test]
fn prepared_row_is_optional_strict_and_cannot_override_complete_file_mode_caps() {
    let bytes = text();
    let ceiling = fs_couple::render::schedule::force::coupled::MAX_PREPARED_MODAL_MODES;
    let over_ceiling = format!("prepared-linear {} 134217728 1000000000000", ceiling + 1);
    for invalid in [
        "prepared-linear 0 134217728 1000000000000",
        over_ceiling.as_str(),
        "prepared-linear 77 134217728 1000000000000",
        "prepared-linear 96 0 1000000000000",
        "prepared-linear 96 134217728 0",
        "prepared-linear 96 134217728 1000000000000 ignored",
        "prepared-linear 96 134217728 -1",
    ] {
        assert!(input::Spec::parse(bytes.replace(PREPARED, invalid).as_bytes()).is_err(),
            "accepted {invalid}");
    }
    assert!(input::Spec::parse(bytes.replace(PREPARED,
        &format!("{PREPARED}\n{PREPARED}")).as_bytes()).is_err());
    let reordered = bytes.replace(&format!("{PREPARED}\n"), "")
        .replace("cylinder-gas-loading", &format!("{PREPARED}\ncylinder-gas-loading"));
    assert!(input::Spec::parse(reordered.as_bytes()).is_err());
    // 78 coordinates at 24 intervals and three more per interval: one interval past the file's mode cap.
    let intervals = 24 + (input::MAX_TOTAL_MODES - 78) / 3 + 1;
    let gas = format!("cylinder-gas-loading 0.001 30 {intervals}");
    let above_file_limit = bytes
        .replace("prepared-linear 96", &format!("prepared-linear {ceiling}"))
        .replace("cylinder-gas-loading 0.001 30 24", &gas);
    assert!(input::Spec::parse(above_file_limit.as_bytes()).is_err());
    let no_fluid = bytes.replace("cylinder-gas-loading 0.001 30 24\n", "");
    assert_eq!(input::Spec::parse(no_fluid.as_bytes()).unwrap().prepared_linear,
        input::Spec::parse(bytes.as_bytes()).unwrap().prepared_linear);
}
