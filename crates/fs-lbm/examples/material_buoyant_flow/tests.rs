//! G1/G3: actual source resolution and the runnable coupled scenario.
use super::*;
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Temperature;

const EXPANDING: &str = include_str!("expanding.txt");
const WEAKER: &str = include_str!("weaker_expansion.txt");

#[test]
fn table_data_reaches_the_original_source_evaluator() {
    let inputs = [MaterialInput::parse(EXPANDING).unwrap(), MaterialInput::parse(WEAKER).unwrap()];
    assert_ne!(inputs[0].card.content_hash(), inputs[1].card.content_hash());
    for (input, density) in inputs.iter().zip([0.99, 0.998]) {
        let law = LiquidTemperatureLaw::try_new(&input.card, &point(300.0), "T",
            LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap();
        let state = law.at_temperature(Temperature::new(400.0)).unwrap();
        let transport = state.transport_state().unwrap();
        assert!((transport.density_kg_m3() - density).abs() < 1e-14);
        assert!((transport.dynamic_viscosity_pa_s() - 1.5).abs() < 1e-14);
        assert!((transport.thermal_conductivity_w_m_k() - 0.15).abs() < 1e-14);
        assert_eq!(transport.specific_heat_j_kg_k(), 100.0);
        assert_eq!(state.material().card_identity(), input.card.content_hash());
        assert!(law.at_temperature(Temperature::new(801.0)).is_err());
    }
    for bad in ["100\n300 1 2 .1\n300 .99 1 .2", "100\n300 1 2 .1\n500 .99 0 .2",
        "100\n300 1 2 .1\n500 .99 1", "NaN\n300 1 2 .1\n500 .99 1 .2"] {
        assert!(MaterialInput::parse(bad).is_err());
    }
}

#[test]
fn data_only_density_swap_changes_motion_and_heat_through_the_real_entry() {
    let mut results = Vec::new();
    for text in [EXPANDING, WEAKER] {
        let material = MaterialInput::parse(text).unwrap();
        let mut emitted = Vec::new();
        let last = simulate(&material, 40, 400.0, |row| emitted.push(row)).unwrap();
        assert_eq!(emitted.len(), 40);
        assert_eq!(last.time_s, 2.0);
        assert!(last.maximum_speed_m_s > 0.0);
        assert!(last.minimum_temperature_k > 300.0);
        assert!(last.maximum_temperature_k < 400.0);
        assert!(last.physical_heat_j > 0.0);
        assert!(last.kinetic_energy_j > 0.0);
        assert!(last.numerical_loss_j > 0.0);
        assert!(last.buoyancy_work_j > last.physical_heat_j);
        assert!((last.thermal_change_j - last.physical_heat_j).abs() < 1e-6);
        for row in &emitted { assert!(row.combined_defect_j.abs() < 1e-6); }
        results.push(last);
    }
    assert!(results[0].maximum_speed_m_s > 2.0 * results[1].maximum_speed_m_s);
    assert!(results[0].physical_heat_j > 2.0 * results[1].physical_heat_j);
    assert!(results[0].maximum_temperature_k < results[1].maximum_temperature_k);
}
