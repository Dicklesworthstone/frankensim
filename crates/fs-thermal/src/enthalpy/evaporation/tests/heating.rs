//! G0/G1: the simultaneous heat/phase receipt updates original native storage.
use super::*;

#[test]
fn supplied_heat_evaporates_native_liquid_without_a_forbidden_intermediate_state() {
    let mut thermal = initial(7.0);
    let original = thermal.cells()[&CellId(1)].clone();
    let duration = 1000.0;
    let transfer = 0.01 * 100_000.0 * 0.018
        / (crate::evaporation::MOLAR_GAS_CONSTANT * 300.0)
        * (100_000.0_f64 / 96_400.0).ln()
        * original.geometry.surface_area_m2 * duration;
    let heat = transfer * LATENT;
    assert!(300.0 + heat / (original.mass_kg * CP) > 330.0);
    assert!(transfer < original.mass_kg);
    let step = thermal.evaporate_cell_with_heat(
        CellId(1), &Liquid, film(0.0), duration, heat, phase(), || Ok(()),
    ).unwrap();
    assert_eq!(thermal.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(step.supplied_heat_j, heat);
    close(step.transferred_mass_kg, transfer, 1e-16);
    close(step.latent_energy_j, heat, 1e-8);
    close(thermal.mass_kg() + step.transferred_mass_kg, original.mass_kg, 1e-17);
    close(thermal.energy_j() + step.vapor_energy_j, heat, 1e-8);
    assert_eq!(thermal.cells()[&CellId(1)].geometry, original.geometry);
    assert_eq!(thermal.cells()[&CellId(1)].boundary, original.boundary);
    assert_eq!(thermal.time_s(), 0.0);
    // The authored 7 W source is still evaluated by its ordinary heat owner.
    // This direct operation applies exactly the supplied energy receipt.
}

#[test]
fn heated_native_cell_refuses_invalid_heat_domain_and_final_cancellation_atomically() {
    let mut successful = initial(3.0);
    let mut polls = 0;
    successful.evaporate_cell_with_heat(
        CellId(1), &Liquid, film(0.0), 1.0, 200.0, phase(),
        || { polls += 1; Ok(()) },
    ).unwrap();
    for cancel_at in [1, polls] {
        let mut thermal = initial(3.0);
        let before = thermal.cells().clone();
        let mut count = 0;
        let result = thermal.evaporate_cell_with_heat(
            CellId(1), &Liquid, film(0.0), 1.0, 200.0, phase(),
            || {
                count += 1;
                if count == cancel_at { Err(DepositionError::Cancelled) } else { Ok(()) }
            },
        );
        assert_eq!(result, Err(ThermalEvaporationError::Thermal(DepositionError::Cancelled)));
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
    }
    let mut thermal = initial(3.0);
    let before = thermal.cells().clone();
    for heat in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            thermal.evaporate_cell_with_heat(
                CellId(1), &Liquid, film(0.0), 1.0, heat, phase(), || Ok(()),
            ),
            Err(ThermalEvaporationError::Surface(EvaporationError::Invalid(_)))
        ));
        assert_eq!(thermal.cells(), &before);
    }
    assert!(matches!(
        thermal.evaporate_cell_with_heat(
            CellId(1), &Liquid, film(0.0), 1.0, 200.0,
            SolverConfig {
                specific_enthalpy_bracket_j_kg: [0.0, 100.0],
                ..phase()
            },
            || Ok(()),
        ),
        Err(ThermalEvaporationError::Surface(EvaporationError::Unbracketed { .. }))
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn heated_native_cell_retains_original_inventory_admission() {
    struct FixedInventory;
    impl ThermalLaw for FixedInventory {
        fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
            if cell.mass_kg != 0.02 {
                return Err(DepositionError::InvalidInput("fixed reference inventory changed"));
            }
            Ok(())
        }
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            CALORIC.temperature_k(h)
        }
        fn maximum_temperature_slope(&self) -> f64 { 1.0 / CP }
    }
    let law = FixedInventory;
    let mut thermal = DepositionSimulation::new(
        vec![birth(0.0, cell(1, 0.0, 0.0))],
        &[(MaterialId(1), &law)],
        DepositionBudget { max_cells: 1, max_links: 0, max_events: 1, max_substeps: 4 },
    ).unwrap().advance_to(0.0, config(1.0), || Ok(())).unwrap().0;
    let before = thermal.cells().clone();
    assert_eq!(
        thermal.evaporate_cell_with_heat(
            CellId(1), &Liquid, film(0.0), 1.0, 200.0, phase(), || Ok(()),
        ),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("fixed reference inventory changed")
        )),
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}
