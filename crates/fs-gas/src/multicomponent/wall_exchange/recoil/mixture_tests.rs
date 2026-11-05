use super::*;
use crate::multicomponent::Conserved;
use crate::multicomponent::tests::{gas, source};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};
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

// Manufactured specific capacity c(T)=a+bT, with an explicitly declared work
// convention. This is an exact constitutive specimen, not measured authority.
struct LinearCapacity {
    a: f64,
    b: f64,
    basis: ThermalStorageBasis,
}
impl LinearCapacity {
    fn energy(&self, temperature: f64) -> f64 {
        self.a * temperature + 0.5 * self.b * temperature * temperature
    }
}
impl ThermalLaw for LinearCapacity {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !energy.is_finite() || energy <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(2.0 * energy / (self.a + (self.a * self.a + 2.0 * self.b * energy).sqrt()))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.a
    }
}

fn storage_cell(energy: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(17),
        segment_id: 3,
        body_id: 8,
        material_id: MaterialId(29),
        mass_kg: 0.6,
        specific_h_j_kg: energy,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: 0.1,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 2.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn material_wall(
    a: f64,
    b: f64,
    range_k: [f64; 2],
    basis: ThermalStorageBasis,
) -> MaterialWallCaloric<'static> {
    let law = LinearCapacity { a, b, basis };
    let cell = storage_cell(law.energy(3.0));
    let range = range_k.map(|temperature| law.energy(temperature));
    MaterialWallCaloric::from_thermal_law(law, cell, range).unwrap()
}

fn mixture() -> (MixtureGas, Conserved) {
    let g = gas(vec![
        source("A", R, 3.5, 0.1, 1e10),
        source("B", 2.0 * R, 4.3, 0.05, -2e10),
    ]);
    // rho=1, R_specific=0.8, so pressure 1.6 gives T=2 K. The original
    // species have cv_density(T)=2.16+0.07*T despite their large energy gauges.
    let state = g.conserved(&[0.6, 0.4], 0.9, 1.6).unwrap();
    (g, state)
}

fn material_recoil(
    gas: &MixtureGas,
    state: &Conserved,
    wall: &MaterialWallCaloric<'_>,
    conductance: f64,
    fraction: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<RecoilUpdate, WallExchangeError> {
    exchange_cell_with_recoil_and_heat_storage(
        gas,
        state.as_cell(),
        0.4,
        Time::new(0.2),
        SlidingSupport {
            mass_kg: 0.6,
            velocity_m_s: -0.2,
        },
        RecoilLaw {
            linear_rate_per_s: 0.8,
            quadratic_per_m: 0.7,
            gas_heat_fraction: fraction,
        },
        wall.initial_contact(conductance, WallDrag::None).unwrap().thermal,
        wall,
        RecoilControl {
            wall: WallControl {
                absolute_energy_tolerance_j: 1e-11,
                relative_energy_tolerance: 1e-12,
                maximum_iterations: 96,
            },
            momentum_tolerance_ns: 1e-12,
        },
        cancelled,
    )
}

fn analytic_slip_loss() -> f64 {
    // Independent reciprocal-slip solution for the fixture's frozen a+b|s|.
    let slip = 1.1;
    let reduced_mass = 0.4 * 0.6 / (0.4 + 0.6);
    let exponent: f64 = 0.8 * 0.2 * (1.0 + 0.4 / 0.6);
    let final_slip = slip / (exponent.exp() + 0.7 * slip / 0.8 * exponent.exp_m1());
    0.5 * reduced_mass * (slip * slip - final_slip * final_slip)
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 3e-10, "{actual:e} != {expected:e}");
}

#[test]
fn material_slip_heating_inverts_original_capacity_integral() {
    let (g, state) = mixture();
    let wall = material_wall(0.7, 0.3, [0.1, 10.0], ThermalStorageBasis::InternalEnergy);
    let result = material_recoil(&g, &state, &wall, 0.0, 0.25, || false).unwrap();
    let heat = 0.75 * analytic_slip_loss();
    let initial_capacity = 0.6 * (0.7 + 0.3 * 3.0);
    // Exact positive root of C0*dT + m*b*dT²/2 = Q; no production inversion.
    let increment = 2.0 * heat
        / (initial_capacity
            + (initial_capacity * initial_capacity + 2.0 * 0.6 * 0.3 * heat).sqrt());
    let actual_temperature = wall.temperature_k(result.gas.contact.thermal).unwrap();
    close(actual_temperature, 3.0 + increment);
    assert!((actual_temperature - (3.0 + heat / initial_capacity)).abs() > 1e-4);
    close(result.slip_dissipation_j, analytic_slip_loss());
    close(result.wall_slip_heating_j, heat);
    close(result.gas.report.wall_energy_change_j, heat);
    close(result.gas.report.heat_to_gas_j, -heat);
    close(result.energy_residual_j, 0.0);
}

#[test]
fn nonlinear_mixture_and_wall_close_original_storage_after_recoil_and_heat_exchange() {
    let (g, state) = mixture();
    let wall = material_wall(0.7, 0.3, [0.1, 10.0], ThermalStorageBasis::InternalEnergy);
    let result = material_recoil(&g, &state, &wall, 0.3, 0.25, || false).unwrap();
    let actual = Cell {
        partial_densities: &state.partial_densities,
        momentum: result.gas.momentum,
        energy: result.gas.energy,
    };
    let gas_temperature = g.primitive(actual, || false).unwrap().temperature;
    let wall_temperature = wall.temperature_k(result.gas.contact.thermal).unwrap();
    let wall_change =
        0.6 * (0.7 + 0.15 * (wall_temperature + 3.0)) * (wall_temperature - 3.0);
    let gas_heat = 0.4 * (2.16 + 0.035 * (gas_temperature + 2.0)) * (gas_temperature - 2.0);
    let transferred = 0.2 * 0.3 * (wall_temperature - gas_temperature);
    let loss = analytic_slip_loss();
    close(wall_change, 0.75 * loss - transferred);
    close(gas_heat, 0.25 * loss + transferred);
    close(result.gas.report.wall_energy_change_j, wall_change);
    close(result.gas.report.heat_to_gas_j, transferred - 0.75 * loss);
    let support_kinetic_change = 0.5
        * 0.6
        * (result.support.velocity_m_s + 0.2)
        * (result.support.velocity_m_s - 0.2);
    close(
        0.4 * (actual.energy - state.energy) + support_kinetic_change + wall_change,
        0.0,
    );
    close(
        0.4 * (actual.momentum - state.momentum)
            + 0.6 * (result.support.velocity_m_s + 0.2),
        0.0,
    );
    close(result.energy_residual_j, 0.0);
    assert_eq!(actual.partial_densities, [0.6, 0.4]);
    assert_eq!(result.gas.report.external_heat_to_gas_j, 0.0);
}

#[test]
fn original_internal_energy_plateau_retains_frictional_heat_without_temperature_rise() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn storage_basis(&self) -> ThermalStorageBasis {
            ThermalStorageBasis::InternalEnergy
        }
        fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
            if !(0.0..=10.0).contains(&energy) {
                return Err(DepositionError::MaterialDomain);
            }
            Ok(if energy < 3.0 {
                1.5 + 0.5 * energy
            } else if energy <= 8.0 {
                3.0
            } else {
                3.0 + energy - 8.0
            })
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0
        }
    }
    let (g, state) = mixture();
    let wall =
        MaterialWallCaloric::from_thermal_law(Plateau, storage_cell(5.0), [0.0, 10.0]).unwrap();
    let result = material_recoil(&g, &state, &wall, 0.0, 0.0, || false).unwrap();
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = result.gas.contact.thermal
    else {
        panic!("original plateau coordinate was lost")
    };
    assert_eq!(wall.temperature_k(result.gas.contact.thermal).unwrap(), 3.0);
    close(mass_kg * (specific_energy_j_kg - 5.0), analytic_slip_loss());
    close(result.gas.report.wall_energy_change_j, analytic_slip_loss());
    close(result.energy_residual_j, 0.0);
}

#[test]
fn constant_capacity_source_recovers_existing_recoil_and_thermal_limit() {
    let (g, state) = mixture();
    let wall = material_wall(0.8 / 0.6, 0.0, [0.1, 10.0], ThermalStorageBasis::InternalEnergy);
    let bound = material_recoil(&g, &state, &wall, 0.3, 0.5, || false).unwrap();
    let legacy = exchange_cell_with_recoil(
        &g,
        state.as_cell(),
        0.4,
        Time::new(0.2),
        SlidingSupport {
            mass_kg: 0.6,
            velocity_m_s: -0.2,
        },
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
    let WallThermal::Lumped { temperature_k, .. } = legacy.gas.contact.thermal else {
        panic!()
    };
    close(wall.temperature_k(bound.gas.contact.thermal).unwrap(), temperature_k);
    close(bound.gas.energy, legacy.gas.energy);
    close(bound.gas.report.wall_energy_change_j, legacy.gas.report.wall_energy_change_j);
    close(bound.gas.report.heat_to_gas_j, legacy.gas.report.heat_to_gas_j);
    assert_eq!(bound.gas.momentum, legacy.gas.momentum);
    assert_eq!(bound.support, legacy.support);
}

#[test]
fn original_wall_domain_work_constraint_and_late_cancellation_refuse_atomically() {
    let (g, state) = mixture();
    let original = state.clone();
    let wall = material_wall(0.7, 0.3, [0.1, 10.0], ThermalStorageBasis::InternalEnergy);
    let narrow =
        material_wall(0.7, 0.3, [2.9, 3.000001], ThermalStorageBasis::InternalEnergy);
    assert!(material_recoil(&g, &state, &narrow, 0.3, 0.5, || false).is_err());
    for basis in [ThermalStorageBasis::Enthalpy, ThermalStorageBasis::DeclaredSensibleStorage] {
        let unspecified = material_wall(0.7, 0.3, [0.1, 10.0], basis);
        assert_eq!(unspecified.storage_basis(), basis);
        assert!(material_recoil(&g, &state, &unspecified, 0.3, 0.5, || false).is_err());
    }
    let mut polls = 0;
    let accepted = material_recoil(&g, &state, &wall, 0.3, 0.5, || {
        polls += 1;
        false
    })
    .unwrap();
    let mut visits = 0;
    assert_eq!(
        material_recoil(&g, &state, &wall, 0.3, 0.5, || {
            visits += 1;
            visits == polls
        }),
        Err(WallExchangeError::Cancelled)
    );
    assert_eq!(state, original);
    assert_eq!(
        material_recoil(&g, &state, &wall, 0.3, 0.5, || false).unwrap(),
        accepted
    );
}
