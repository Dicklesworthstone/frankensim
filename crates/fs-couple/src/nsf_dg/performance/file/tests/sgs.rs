//! Explicit SGS stays in the original open NSF owner and its live DG solve.
use super::*;

fn declare(text: &str, model: &str, coefficient: f64) -> String {
    text.replace(
        "wall adiabatic\n",
        &format!("wall adiabatic\nsgs exploration {model} {coefficient} 1 100 0.9\n"),
    )
}

fn open_source() -> String {
    let gas = reference(100000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let full_radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
    let patch_radius = (dx * dx / core::f64::consts::PI).sqrt();
    fixture(gas, "dry-air 300 100000")
        .replace(
            "characteristic 0 high 0.05\n",
            "characteristic 0 high 0.05\ncoupled_patch 1 2 1 2 pressure 0\n",
        )
        .replace(
            "initial rest",
            &format!(
                "reservoir_aperture 1 low {:.17e} 0 0 0 300 1\nreservoir_cell 6 0 1\ninitial rest",
                gas.density
            ),
        )
        .replace(
            "states 0",
            &format!(
                "states 1\nstate 4 1 1 {:.17e} 0.02 0.01 -0.015 300.01",
                gas.density
            ),
        )
        .replace(
            &format!("dg_tube 0.1 {full_radius:.17e}"),
            &format!("dg_tube 0.1 {patch_radius:.17e}"),
        )
}

#[test]
fn explicit_zero_sgs_preserves_open_molecular_playback_exactly() {
    let source = open_source();
    let (baseline, molecular) = render(read(&source, 8), 8);
    for model in ["smagorinsky", "wale", "vreman"] {
        let file = read(&declare(&source, model, 0.0), 3);
        assert_eq!(file.info().sgs.unwrap().closure.model.coefficient(), 0.0);
        assert_eq!(file.renderer().system().zone().sgs_config(), None);
        let (pressure, rendered) = render(file, 3);
        assert_eq!(
            pressure, baseline,
            "zero {model} must retain molecular arithmetic"
        );
        assert_eq!(
            rendered.system().zone().state_hash(),
            molecular.system().zone().state_hash()
        );
        assert_eq!(
            rendered.system().load().state(),
            molecular.system().load().state()
        );
    }
}

#[test]
fn native_sgs_reaches_live_patch_reservoir_and_conservative_coupling() {
    let source = open_source();
    let declared = declare(&source, "smagorinsky", 0.1);
    let file = read(&declared, 8);
    let config = file.info().sgs.unwrap();
    assert_eq!(file.renderer().system().zone().sgs_config(), Some(config));
    assert!(
        file.renderer()
            .system()
            .zone()
            .sgs_transport([4, 1, 1])
            .unwrap()
            .unwrap()
            .eddy_viscosity_lu
            > 0.0
    );
    let (pressure, whole) = render(file, 8);
    let (parts, partitioned) = render(read(&declared, 3), 3);
    assert_eq!(pressure, parts);
    assert_eq!(
        whole.system().zone().state_hash(),
        partitioned.system().zone().state_hash()
    );
    assert_eq!(
        whole.system().load().state(),
        partitioned.system().load().state()
    );
    let (_, molecular) = render(read(&source, 8), 8);
    assert_ne!(
        whole.system().zone().cell([4, 1, 1]),
        molecular.system().zone().cell([4, 1, 1])
    );

    let mut renderer = read(&declared, 1).into_renderer();
    let mass = renderer.system().zone().total_mass_kg();
    let energy = renderer.system().zone().total_energy_j();
    for _ in 0..SAMPLES {
        renderer.block(&mut [0.0]).unwrap();
        let frame = renderer.last_frame().unwrap();
        assert!(frame.coupled_patch.is_some());
        assert_eq!(renderer.system().zone().sgs_config(), Some(config));
        let mass_error = frame.flow_mass_change_kg
            + frame.reservoir_mass_out_kg
            + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
        let energy_error = frame.flow_energy_change_j
            - frame.flow.open_energy_j
            - frame.flow.open_boundary_energy_j
            - frame.flow.wall_heat_j
            - frame.flow.source_energy_j;
        assert!(mass_error.abs() < 5e-11 * mass, "{frame:?}");
        assert!(energy_error.abs() < 5e-11 * energy, "{frame:?}");
    }
}

#[test]
fn native_sgs_requires_explicit_valid_exploration_controls() {
    let source = fixture(reference(100000.0), "dry-air 300 100000");
    for record in [
        "sgs acceptance smagorinsky 0.1 1 100 0.9",
        "sgs exploration unknown 0.1 1 100 0.9",
        "sgs exploration wale -0.1 1 100 0.9",
        "sgs exploration vreman 0.1 0 100 0.9",
        "sgs exploration smagorinsky 0.1 1 -1 0.9",
        "sgs exploration smagorinsky 0.1 1 100 0",
        "sgs exploration smagorinsky NaN 1 100 0.9",
        "sgs exploration smagorinsky 0.1 1 100",
        "sgs exploration smagorinsky 0.1 1 100 0.9 extra",
    ] {
        let text = source.replace("wall adiabatic\n", &format!("wall adiabatic\n{record}\n"));
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                text.as_bytes(),
                8,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err(),
            "admitted {record}"
        );
    }
}
