use super::*;

const THERMAL: &str = "air-thermal isothermal-walls 1.0 43 256 50000000 0.201\n";
const THERMAL_BODY: &str = include_str!("../../closed-shell-thermal.body");

fn profile() -> String {
    THERMAL_BODY.to_owned()
}

#[test]
fn authored_heat_field_joins_the_original_material_state_budget() {
    let body = body::Profile::parse(profile().as_bytes()).unwrap();
    assert!(body.has_thermal_cavity());
    assert_eq!(body.material_memory_states(), 43);
    assert_eq!(body.additional_body_modes(), 0);
    assert_eq!(body.additional_body_connections(), 1);
    assert!(!body::Profile::parse(BODY.as_bytes()).unwrap().has_thermal_cavity());
    let base = profile();
    for (from, to) in [
        ("isothermal-walls", "unspecified-wall"),
        ("1.0 43", "0 43"),
        ("1.0 43", "1.0 0"),
        ("1.0 43", "1.0 257"),
        ("43 256", "43 204"),
        ("50000000 0.201", "0 0.201"),
        ("50000000 0.201", "50000000 0"),
        ("50000000 0.201", "50000000 1"),
        ("50000000 0.201", "50000000 NaN"),
    ] {
        assert!(body::Profile::parse(base.replace(from, to).as_bytes()).is_err(), "admitted {to}");
    }
    assert!(body::Profile::parse(format!("{base}{THERMAL}").as_bytes()).is_err());
}

#[test]
fn thermal_source_refuses_a_backend_without_memory_before_material_reads() {
    scoped(|cx, gate| {
        let mut spec = source(PERFORMANCE);
        spec.prepared_linear = Some(PreparedModalBudget {
            max_modes: 128, max_dense_bytes: 67_108_864,
            max_multiply_terms: 200_000_000_000,
        });
        let mut reads = 0;
        let result = spec.compile_with_body::<NormalizedPack>(cx, gate,
            Some(body::Profile::parse(profile().as_bytes()).unwrap()), |_| {
                reads += 1;
                Err("material I/O after unsupported thermal propagation".into())
            });
        assert!(result.is_err());
        assert_eq!(reads, 0);
    });
}

#[test]
fn authored_thermal_air_changes_shared_motion_and_survives_source_gas_loading() {
    scoped(|cx, gate| {
        let mut loaded = ready(cx, gate, source(PERFORMANCE), THERMAL_BODY, 700.0);
        let mut replay = ready(cx, gate, source(PERFORMANCE), THERMAL_BODY, 700.0);
        let uniform = BODY.replace(RESOLUTION,
            "interior-fe-v1 2 1 1 256 80 10000000 1e-8");
        let mut adiabatic = ready(cx, gate, source(PERFORMANCE), &uniform, 700.0);
        let image = loaded.performance.bank().source().bank().spatial_cavity().unwrap();
        let thermal = image.thermal_loss().unwrap();
        assert_eq!(thermal.diffusion_modes().free_dofs, 43);
        assert_eq!(thermal.diffusion_modes().shapes.len(), 43);
        assert!((thermal.unresolved_fraction_bound() - 0.2).abs() < 1e-7);
        assert!(thermal.memory_states() > 0 && thermal.memory_states() <= 43);
        assert_eq!(loaded.performance.bank().system().material_memory().len(), thermal.memory_states());
        let diffusivity = thermal.diffusivity_m2_s();
        let eigenvalues = thermal.diffusion_modes().eigenvalues_per_m2.clone();
        let ordinary_audio = pressure(&mut loaded, 1);
        assert_eq!(ordinary_audio, pressure(&mut replay, 37));
        assert_ne!(ordinary_audio, pressure(&mut adiabatic, 128));
        let frame = loaded.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);

        let doubled = ready(cx, gate, source(&PERFORMANCE.replace("293.15 101325",
            "293.15 202650")), THERMAL_BODY, 700.0);
        let doubled_thermal = doubled.performance.bank().source().bank()
            .spatial_cavity().unwrap().thermal_loss().unwrap();
        assert!((doubled_thermal.diffusivity_m2_s() / diffusivity - 0.5).abs() < 1e-10);
        assert_eq!(doubled_thermal.diffusion_modes().eigenvalues_per_m2, eigenvalues);
        let flow = PERFORMANCE.replace("cavity none", "cylinder-gas-loading 0.01 4 2\ncavity none");
        let mut fluid = ready(cx, gate, source(&flow), THERMAL_BODY, 700.0);
        let bank = fluid.performance.bank().source().bank();
        assert!(bank.cylinder_gas_loading().is_some());
        assert!(bank.spatial_cavity().unwrap().thermal_loss().is_some());
        pressure(&mut fluid, 37);
        let frame = fluid.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}
