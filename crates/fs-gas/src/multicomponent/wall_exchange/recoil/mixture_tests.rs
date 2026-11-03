use super::*;
use crate::multicomponent::tests::{gas, source};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

#[test]
fn original_mixture_owner_preserves_species_and_closes_recoil_storage() {
    let gas = gas(vec![
        source("A", R, 3.5, 0.1, 1e10),
        source("B", 2.0 * R, 4.3, 0.05, -2e10),
    ]);
    let state = gas.conserved(&[0.6, 0.4], 0.9, 1.6).unwrap();
    let support = SlidingSupport {
        mass_kg: 0.6,
        velocity_m_s: -0.2,
    };
    let update = exchange_cell_with_recoil(
        &gas,
        state.as_cell(),
        0.4,
        Time::new(0.2),
        support,
        RecoilLaw {
            linear_rate_per_s: 0.8,
            quadratic_per_m: 0.7,
            gas_heat_fraction: 0.5,
        },
        WallThermal::Lumped {
            temperature_k: 3.0,
            heat_capacity_j_per_k: 0.8,
            conductance_w_per_k: 0.3,
        },
        RecoilControl {
            wall: WallControl::default(),
            momentum_tolerance_ns: 1e-12,
        },
        || false,
    )
    .unwrap();
    let actual = Cell {
        partial_densities: &state.partial_densities,
        momentum: update.gas.momentum,
        energy: update.gas.energy,
    };
    gas.primitive(actual, || false).unwrap();
    assert_eq!(actual.partial_densities, [0.6, 0.4]);
    let WallThermal::Lumped { temperature_k, .. } = update.gas.contact.thermal else {
        panic!()
    };
    let energy = 0.4 * (actual.energy - state.energy)
        + 0.5
            * support.mass_kg
            * (update.support.velocity_m_s - support.velocity_m_s)
            * (update.support.velocity_m_s + support.velocity_m_s)
        + 0.8 * (temperature_k - 3.0);
    assert!(
        energy.abs() < 1e-9,
        "actual coupled energy defect {energy:e}"
    );
    assert!(update.momentum_residual_ns.abs() < 1e-12);
}
