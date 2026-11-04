//! Mechanical work reaches the native velocity without changing caloric history.

use super::*;
use fs_lbm::impulse::ImpulseError;
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceImpulseControls3, MaterialFreeSurfaceImpulseError,
};

#[test]
fn impulse_preserves_thermal_rheology_sources_enthalpy_and_clock() {
    let card = rheological_card(0.0002, 0.003, 0.0, 1000.0);
    with_liquid(
        &card,
        true,
        true,
        [0.03, -0.02, 0.01],
        solve(),
        |p| 200000.0 + 1000.0 * p[0] as f64,
        |liquid, gate| {
            liquid.step(gate).unwrap();
            let before = liquid.clone();
            let [x, y, z] = CELL;
            let index = (z * 8 + y) * 8 + x;
            let velocity = before.flow().sample_at_si(x, y, z).unwrap().unwrap().velocity_m_s;
            let mass = before.flow().lattice().region_mass(|a, b, c| [a, b, c] == CELL)
                * before.flow().scaling().mass_unit_kg();
            let volume = before.flow().scaling().cell_size_m().powi(3);
            let controls = MaterialFreeSurfaceImpulseControls3 {
                max_cells: 512,
                maximum_speed_m_s: 0.1,
                maximum_momentum_defect_kg_m_s: 1e-16,
                maximum_work_defect_j: 1e-17,
            };
            let mut field = vec![[0.0; 3]; 512];
            field[index] = [2.0, -0.5, 0.25];
            let report = liquid.apply_impulse_density(&field, controls, gate).unwrap();
            let after = liquid.flow().sample_at_si(x, y, z).unwrap().unwrap().velocity_m_s;
            let mut work = 0.0;
            for a in 0..3 {
                let impulse = field[index][a] * volume;
                assert!((after[a] - velocity[a] - impulse / mass).abs() < 2e-13);
                assert!((report.applied_impulse_kg_m_s[a] - impulse).abs() < 1e-22);
                work += impulse * (velocity[a] + after[a]) * 0.5;
            }
            assert!((report.work_j - work).abs() < 1e-22);
            assert!((report.kinetic_energy_change_j - work).abs() < 1e-20);
            assert_eq!(liquid.accepted_steps(), before.accepted_steps());
            assert_eq!(liquid.physical_time_s(), before.physical_time_s());
            assert_eq!(liquid.total_enthalpy_j().unwrap().to_bits(), before.total_enthalpy_j().unwrap().to_bits());
            assert_eq!(liquid.flow().ledger_mass_kg().unwrap(), before.flow().ledger_mass_kg().unwrap());
            for z in 0..8 {
                for y in 0..8 {
                    for x in 0..8 {
                        let p = [x, y, z];
                        assert_eq!(liquid.specific_enthalpy_at(p), before.specific_enthalpy_at(p));
                        assert_eq!(liquid.material_at(p), before.material_at(p));
                        assert_eq!(liquid.rheology_at(p), before.rheology_at(p));
                    }
                }
            }
            let accepted = liquid.flow().lattice().grid().state_bits();
            field[0] = [1.0, 0.0, 0.0];
            field[index] = [1e9, 0.0, 0.0];
            assert!(liquid.apply_impulse_density(&field, controls, gate).is_err());
            assert_eq!(liquid.flow().lattice().grid().state_bits(), accepted);
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(matches!(
                liquid.apply_impulse_density(&field, controls, &cancelled),
                Err(MaterialFreeSurfaceImpulseError::Native(ImpulseError::Cancelled))
            ));
            assert_eq!(liquid.flow().lattice().grid().state_bits(), accepted);
            assert_eq!(liquid.accepted_steps(), before.accepted_steps());
            assert_eq!(liquid.total_enthalpy_j().unwrap(), before.total_enthalpy_j().unwrap());
        },
    );
}
