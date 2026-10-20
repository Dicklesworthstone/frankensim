//! Host-only admission of a GPU plan: no CUDA feature, driver or device needed.

#[path = "../src/plan.rs"]
#[allow(dead_code)]
mod plan;

use fs_vkplate::{GasState, Instrument, kit, standard_gas};

fn drum(gas: GasState) -> Instrument {
    let mut spec = kit::kick();
    spec.grid *= 8.0;
    kit::drum_with_gas(&spec, &[], gas).unwrap().instrument
}

#[test]
fn shared_plan_refuses_a_batch_member_with_another_atmosphere() {
    let gas = standard_gas();
    let mut first = drum(gas);
    let image = first.image();
    let plan = plan::Plan::new(&image, &first, 1, 64, true).unwrap();
    assert!(plan.state(&image).is_ok());
    assert!(plan.state(&first.clone().image()).is_ok());

    let compressed = GasState::try_new_moist_air(293.15, 2.0 * gas.pressure, 0.0).unwrap();
    // A coherent state can also share rho and c while carrying different
    // transport data. Admission retains the whole material snapshot.
    let mut other_transport = gas;
    other_transport.dynamic_viscosity *= 2.0;
    other_transport.thermal_conductivity *= 2.0;
    other_transport.validate_coherence().unwrap();
    assert_eq!(other_transport.density, gas.density);
    assert_eq!(other_transport.sound_speed, gas.sound_speed);

    for other in [compressed, other_transport] {
        let mut member = drum(other);
        let other_image = member.image();
        assert!(matches!(plan.state(&other_image), Err(error) if error.contains("atmosphere")));
        // Export also refuses a mismatched image/owner pair before flattening.
        assert!(plan::Plan::new(&other_image, &first, 1, 64, true).is_err());
    }
    assert!(
        plan.state(&image).is_ok(),
        "refusal must not alter the shared plan"
    );
}
