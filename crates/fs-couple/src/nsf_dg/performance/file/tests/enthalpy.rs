//! G0/G3: explicit thermal contact reaches the original resolved/native owners.
use super::*;

fn contact_input() -> String {
    fixture(reference(100_000.0), "dry-air 300 100000")
        .replace(
            "characteristic 0 high 0.05\n",
            "characteristic 0 high 0.05\ncharacteristic_thermal enthalpy-flux upwind-entropy\n",
        )
        .replace(
            &format!("boundary 1\nhold {SAMPLES} -0.1"),
            &format!("boundary pressure-temperature 1\nhold {SAMPLES} 20 300.5"),
        )
}

#[test]
fn g3_explicit_legacy_thermal_policy_preserves_native_trajectory() {
    let old = fixture(reference(100_000.0), "dry-air 300 100000");
    let selected = old.replace(
        "characteristic 0 high 0.05\n",
        "characteristic 0 high 0.05\ncharacteristic_thermal upwind-entropy upwind-entropy\n",
    );
    let (a, original) = render(read(&old, 8), 8);
    let (b, explicit) = render(read(&selected, 3), 3);
    assert_eq!(a, b);
    assert_eq!(
        original.system().zone().state_hash(),
        explicit.system().zone().state_hash()
    );
    assert_eq!(
        original.system().load().state(),
        explicit.system().load().state()
    );
}

#[test]
fn g3_incoming_heat_uses_the_original_face_and_conservative_energy_account() {
    let mut renderer = read(&contact_input(), 8).into_renderer();
    assert_eq!(
        renderer.system().zone().characteristic_thermal_policy(),
        [
            NsfCharacteristicThermalPolicy::EnthalpyFlux,
            NsfCharacteristicThermalPolicy::UpwindEntropy
        ]
    );
    let energy = renderer.system().zone().total_energy_j();
    renderer.block(&mut [0.0]).unwrap();
    let frame = renderer.last_frame().unwrap();
    let faces = renderer
        .system()
        .zone()
        .open_boundary_traces(false)
        .unwrap();
    assert!(
        faces
            .iter()
            .all(|face| face.outward_normal_velocity_m_s < 0.0)
    );
    assert!(
        faces
            .iter()
            .all(|face| face.temperature_k > 300.0 && face.temperature_k < 300.5)
    );
    let supply = frame.flow.open_boundary_energy_j
        + frame.flow.open_energy_j
        + frame.flow.wall_heat_j
        + frame.flow.source_energy_j;
    assert!((frame.flow_energy_change_j - supply).abs() <= 2e-11 * energy);
    assert_eq!(frame.flow.source_energy_j, 0.0);
    assert_eq!(renderer.system().zone().steps(), 1);
    assert_eq!(renderer.system().load().state().step_count(), 1);
}

#[test]
fn g0_native_thermal_contact_refuses_ambiguous_fields_and_unresolved_sgs() {
    for row in [
        "characteristic_thermal enthalpy-flux",
        "characteristic_thermal upwind-entropy unknown",
        "characteristic_thermal enthalpy-flux upwind-entropy extra",
    ] {
        let text =
            contact_input().replace("characteristic_thermal enthalpy-flux upwind-entropy", row);
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                text.as_bytes(),
                8,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err(),
            "accepted {row}"
        );
    }
    let text = contact_input().replace(
        "wall adiabatic\n",
        "wall adiabatic\nsgs exploration wale 0.1 1 1000 0.7\n",
    );
    assert!(
        ResolvedDgFilePerformance::from_bytes(
            text.as_bytes(),
            8,
            Arc::new(CancelGate::new_clock_free())
        )
        .is_err()
    );
}
