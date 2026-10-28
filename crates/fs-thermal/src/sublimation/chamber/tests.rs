use super::super::{
    MOLAR_GAS_CONSTANT, SolidSurfaceState, SolverConfig, SublimationError, SublimationLaw,
    SublimationProperties, SurfaceKinetics, VaporReservoir,
};
use super::{ChamberSolverConfig, ChamberSurface, advance_chamber};

const MOLAR_MASS: f64 = 0.020;
const LATENT_HEAT: f64 = 500_000.0;

struct ReferenceLaw;

impl SublimationLaw for ReferenceLaw {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 300.0]
    }

    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        if !t.is_finite() || !(100.0..=300.0).contains(&t) {
            return Err(SublimationError::TemperatureOutOfRange {
                temperature_k: t,
                minimum_k: 100.0,
                maximum_k: 300.0,
            });
        }
        let scale = LATENT_HEAT * MOLAR_MASS / MOLAR_GAS_CONSTANT;
        let pressure = 100.0 * (scale * (1.0 / 250.0 - 1.0 / t)).exp();
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure * scale / t / t,
            latent_heat_j_kg: LATENT_HEAT,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: MOLAR_MASS,
        })
    }
}

fn controls() -> ChamberSolverConfig {
    ChamberSolverConfig {
        surface_solver: SolverConfig {
            relative_tolerance: 1.0e-12,
            absolute_energy_tolerance_j: 1.0e-10,
            max_iterations: 256,
        },
        relative_mass_tolerance: 1.0e-11,
        absolute_mass_tolerance_kg: 1.0e-14,
        max_iterations: 128,
        max_surfaces: 16,
    }
}

fn surface(mass: f64, area: f64) -> ChamberSurface {
    ChamberSurface {
        solid: SolidSurfaceState {
            mass_kg: mass,
            specific_enthalpy_j_kg: 0.0,
            area_m2: area,
        },
        kinetics: SurfaceKinetics {
            accommodation_coefficient: 1.0,
            mass_transfer_coefficient_m_s: None,
        },
    }
}

fn reservoir(mass: f64, t: f64) -> VaporReservoir {
    VaporReservoir::Finite {
        vapor_mass_kg: mass,
        volume_m3: 0.1,
        temperature_k: t,
    }
}

fn vapor_mass(gas: VaporReservoir) -> f64 {
    match gas {
        VaporReservoir::Finite { vapor_mass_kg, .. } => vapor_mass_kg,
        VaporReservoir::Prescribed { .. } => panic!("finite chamber required"),
    }
}

fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "actual={actual:e}, expected={expected:e}, tolerance={absolute:e}",
    );
}

// Independent simultaneous backward-Euler algebra at prescribed temperatures.
// No production flux, pressure, coefficient, transfer or nonlinear solver is
// called to compute this oracle.
fn isothermal_oracle(
    surfaces: &[ChamberSurface],
    temperatures: &[f64],
    initial_mass: f64,
    gas_temperature: f64,
    dt: f64,
) -> (f64, Vec<f64>, f64, f64) {
    let beta = MOLAR_GAS_CONSTANT * gas_temperature / MOLAR_MASS / 0.1;
    let mut source = Vec::new();
    let mut conductance = Vec::new();
    let scale = LATENT_HEAT * MOLAR_MASS / MOLAR_GAS_CONSTANT;
    for (surface, &t) in surfaces.iter().zip(temperatures) {
        assert!(surface.kinetics.mass_transfer_coefficient_m_s.is_none());
        let prefactor = surface.solid.area_m2
            * surface.kinetics.accommodation_coefficient
            * (MOLAR_MASS / (2.0 * std::f64::consts::PI * MOLAR_GAS_CONSTANT)).sqrt();
        source.push(prefactor / t.sqrt() * 100.0 * (scale * (1.0 / 250.0 - 1.0 / t)).exp());
        conductance.push(prefactor / gas_temperature.sqrt());
    }
    let sum_source: f64 = source.iter().sum();
    let sum_conductance: f64 = conductance.iter().sum();
    let pressure =
        (beta * initial_mass + beta * dt * sum_source) / (1.0 + beta * dt * sum_conductance);
    let transfers = source
        .iter()
        .zip(&conductance)
        .map(|(s, g)| dt * (s - g * pressure))
        .collect();
    (
        pressure,
        transfers,
        sum_source / sum_conductance,
        beta * sum_conductance,
    )
}

#[test]
fn g1_joint_isothermal_pressure_matches_backward_euler_algebra() {
    let surfaces = [surface(1.0, 0.01), surface(2.0, 0.025)];
    let temperatures = [260.0, 225.0];
    let initial_mass = 2.0e-5;
    let gas_temperature = 240.0;
    let dt = 0.08;
    let (pressure, transfers, _, _) =
        isothermal_oracle(&surfaces, &temperatures, initial_mass, gas_temperature, dt);
    // An imposed-T callback isolates the exact algebraic mass closure.
    let step = advance_chamber(
        &ReferenceLaw,
        &surfaces,
        reservoir(initial_mass, gas_temperature),
        dt,
        controls(),
        |i, _| Ok(temperatures[i]),
        || Ok(()),
    )
    .unwrap();
    close(
        step.reservoir.partial_pressure_pa(MOLAR_MASS).unwrap(),
        pressure,
        2.0e-6,
    );
    for (actual, expected) in step.surfaces.iter().zip(&transfers) {
        close(actual.transferred_mass_kg, *expected, 2.0e-12);
        assert_eq!(actual.reservoir, step.reservoir);
    }
    close(
        vapor_mass(step.reservoir) - initial_mass,
        transfers.iter().sum(),
        2.0e-12,
    );
    assert!(step.mass_residual_kg.abs() <= 2.0e-12);
}

#[test]
fn g1_time_refinement_converges_to_the_shared_exponential_transient() {
    let total_time = 0.03;
    let temperature = 250.0;
    let initial = [surface(1.0, 0.01)];
    let (_, _, equilibrium_pressure, rate) =
        isothermal_oracle(&initial, &[temperature], 0.0, temperature, total_time);
    let beta = MOLAR_GAS_CONSTANT * temperature / MOLAR_MASS / 0.1;
    let exact_mass = equilibrium_pressure / beta * (1.0 - (-rate * total_time).exp());
    let evolve = |count: usize| {
        let mut solid = surface(1.0, 0.01);
        let mut gas = reservoir(0.0, temperature);
        for _ in 0..count {
            let step = advance_chamber(
                &ReferenceLaw,
                &[solid],
                gas,
                total_time / count as f64,
                controls(),
                |_, _| Ok(temperature),
                || Ok(()),
            )
            .unwrap();
            gas = step.reservoir;
            solid = ChamberSurface {
                solid: step.surfaces[0].solid,
                kinetics: SurfaceKinetics {
                    accommodation_coefficient: 1.0,
                    mass_transfer_coefficient_m_s: None,
                },
            };
        }
        (vapor_mass(gas) - exact_mass).abs()
    };
    let errors = [evolve(2), evolve(4), evolve(8)];
    assert!(errors[2] < errors[1] && errors[1] < errors[0], "{errors:?}");
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((1.7..2.2).contains(&ratio), "first-order ratio={ratio}");
    }
}

#[test]
fn g1_warm_to_cold_transfer_continues_at_stationary_chamber_pressure() {
    let surfaces = [surface(1.0, 0.01), surface(1.0, 0.01)];
    let temperatures = [270.0, 200.0];
    let gas_temperature = 240.0;
    let dt = 0.2;
    let (_, _, equilibrium_pressure, _) =
        isothermal_oracle(&surfaces, &temperatures, 0.0, gas_temperature, dt);
    let beta = MOLAR_GAS_CONSTANT * gas_temperature / MOLAR_MASS / 0.1;
    let initial_mass = equilibrium_pressure / beta;
    let step = advance_chamber(
        &ReferenceLaw,
        &surfaces,
        reservoir(initial_mass, gas_temperature),
        dt,
        controls(),
        |i, _| Ok(temperatures[i]),
        || Ok(()),
    )
    .unwrap();
    assert!(step.surfaces[0].transferred_mass_kg > 0.0);
    assert!(step.surfaces[1].transferred_mass_kg < 0.0);
    close(
        step.surfaces[0].transferred_mass_kg,
        -step.surfaces[1].transferred_mass_kg,
        2.0e-12,
    );
    close(vapor_mass(step.reservoir), initial_mass, 2.0e-12);
    close(step.transferred_mass_kg, 0.0, 2.0e-12);
}

#[test]
fn g1_cold_surface_can_receive_vapor_from_a_warm_surface_in_one_step() {
    let surfaces = [surface(1.0, 0.01), surface(1.0, 0.01)];
    let temperatures = [270.0, 200.0];
    let step = advance_chamber(
        &ReferenceLaw,
        &surfaces,
        reservoir(0.0, 240.0),
        0.5,
        controls(),
        |i, _| Ok(temperatures[i]),
        || Ok(()),
    )
    .unwrap();
    // There was no initial gas inventory to copy into an independent cold bath.
    assert!(step.surfaces[0].transferred_mass_kg > 0.0);
    assert!(step.surfaces[1].transferred_mass_kg < 0.0);
    assert!(step.surfaces[1].solid.mass_kg > surfaces[1].solid.mass_kg);
    assert!(vapor_mass(step.reservoir) > 0.0);
    close(
        vapor_mass(step.reservoir),
        step.transferred_mass_kg,
        2.0e-12,
    );
}

#[test]
fn g1_multiple_depositing_surfaces_share_one_nonnegative_vapor_inventory() {
    let surfaces = [
        surface(1.0, 0.01),
        surface(1.0, 0.02),
        surface(1.0, 0.03),
        surface(1.0, 0.04),
    ];
    let initial_mass = 0.001;
    let step = advance_chamber(
        &ReferenceLaw,
        &surfaces,
        reservoir(initial_mass, 240.0),
        10.0,
        controls(),
        |_, _| Ok(200.0),
        || Ok(()),
    )
    .unwrap();
    assert!(step.surfaces.iter().all(|s| s.transferred_mass_kg < 0.0));
    let gas_mass = vapor_mass(step.reservoir);
    assert!(gas_mass >= 0.0 && gas_mass < initial_mass);
    let deposited: f64 = step.surfaces.iter().map(|s| -s.transferred_mass_kg).sum();
    assert!(deposited <= initial_mass + 2.0e-12);
    close(initial_mass - gas_mass, deposited, 2.0e-12);
    let solid_gain: f64 = step
        .surfaces
        .iter()
        .zip(&surfaces)
        .map(|(new, old)| new.solid.mass_kg - old.solid.mass_kg)
        .sum();
    close(solid_gain, deposited, 2.0e-12);
}

#[test]
fn g1_partitioning_area_and_permuting_surfaces_preserves_the_chamber() {
    let gas = reservoir(1.0e-5, 250.0);
    let whole = advance_chamber(
        &ReferenceLaw,
        &[surface(2.0, 0.02)],
        gas,
        0.04,
        controls(),
        |_, _| Ok(250.0),
        || Ok(()),
    )
    .unwrap();
    let pieces = [
        surface(0.5, 0.005),
        surface(0.5, 0.005),
        surface(0.5, 0.005),
        surface(0.5, 0.005),
    ];
    let partitioned = advance_chamber(
        &ReferenceLaw,
        &pieces,
        gas,
        0.04,
        controls(),
        |_, _| Ok(250.0),
        || Ok(()),
    )
    .unwrap();
    close(
        vapor_mass(whole.reservoir),
        vapor_mass(partitioned.reservoir),
        2.0e-12,
    );
    close(
        whole.transferred_mass_kg,
        partitioned.transferred_mass_kg,
        2.0e-12,
    );
    for part in &partitioned.surfaces {
        close(
            part.solid.specific_enthalpy_j_kg,
            whole.surfaces[0].solid.specific_enthalpy_j_kg,
            2.0e-5,
        );
    }

    let ordered = [surface(2.0, 0.003), surface(0.7, 0.013)];
    let reversed = [surface(0.7, 0.013), surface(2.0, 0.003)];
    let first = advance_chamber(
        &ReferenceLaw,
        &ordered,
        gas,
        0.12,
        controls(),
        |i, _| Ok([260.0, 210.0][i]),
        || Ok(()),
    )
    .unwrap();
    let second = advance_chamber(
        &ReferenceLaw,
        &reversed,
        gas,
        0.12,
        controls(),
        |i, _| Ok([210.0, 260.0][i]),
        || Ok(()),
    )
    .unwrap();
    close(
        vapor_mass(first.reservoir),
        vapor_mass(second.reservoir),
        2.0e-12,
    );
    for i in 0..2 {
        close(
            first.surfaces[i].transferred_mass_kg,
            second.surfaces[1 - i].transferred_mass_kg,
            2.0e-12,
        );
    }
}

#[test]
fn g0_cancellation_domain_and_surface_capacity_refuse_atomically() {
    let surfaces = [surface(1.0, 0.01), surface(1.0, 0.01)];
    let gas = reservoir(0.0, 240.0);
    let mut polls = 0;
    let cancelled = advance_chamber(
        &ReferenceLaw,
        &surfaces,
        gas,
        0.1,
        controls(),
        |_, _| Ok(250.0),
        || {
            polls += 1;
            if polls >= 3 {
                Err(SublimationError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(cancelled, Err(SublimationError::Cancelled)));
    assert!(
        advance_chamber(
            &ReferenceLaw,
            &surfaces,
            gas,
            0.1,
            controls(),
            |_, _| Ok(350.0),
            || Ok(()),
        )
        .is_err()
    );
    let mut limited = controls();
    limited.max_surfaces = 1;
    assert!(matches!(
        advance_chamber(
            &ReferenceLaw,
            &surfaces,
            gas,
            0.1,
            limited,
            |_, _| Ok(250.0),
            || Ok(()),
        ),
        Err(SublimationError::Invalid(_))
    ));
    assert!(matches!(
        advance_chamber(
            &ReferenceLaw,
            &surfaces,
            VaporReservoir::Prescribed {
                partial_pressure_pa: 0.0,
                temperature_k: 240.0
            },
            0.1,
            controls(),
            |_, _| Ok(250.0),
            || Ok(()),
        ),
        Err(SublimationError::Invalid(_))
    ));
    assert_eq!(surfaces[0].solid.mass_kg, 1.0);
    assert_eq!(surfaces[1].solid.specific_enthalpy_j_kg, 0.0);
    assert_eq!(vapor_mass(gas), 0.0);
}

#[test]
fn g1_finite_caloric_capacities_close_mass_sensible_and_latent_energy() {
    let surfaces = [surface(0.1, 0.02), surface(0.2, 0.03)];
    let initial_temperatures = [270.0, 220.0];
    let heat_capacities = [1000.0, 1800.0];
    let initial_mass = 2.0e-5;
    let step = advance_chamber(
        &ReferenceLaw,
        &surfaces,
        reservoir(initial_mass, 240.0),
        0.05,
        controls(),
        |i, h| Ok(initial_temperatures[i] + h / heat_capacities[i]),
        || Ok(()),
    )
    .unwrap();
    let mut stored_change = 0.0;
    let mut carried = 0.0;
    let mut latent = 0.0;
    for (i, (old, new)) in surfaces.iter().zip(&step.surfaces).enumerate() {
        let dm = old.solid.mass_kg - new.solid.mass_kg;
        close(dm, new.transferred_mass_kg, 2.0e-12);
        close(
            new.surface_temperature_k,
            initial_temperatures[i] + new.solid.specific_enthalpy_j_kg / heat_capacities[i],
            2.0e-10,
        );
        assert!((new.surface_temperature_k - initial_temperatures[i]) * dm <= 0.0);
        let actual_stored = new.solid.mass_kg * new.solid.specific_enthalpy_j_kg
            - old.solid.mass_kg * old.solid.specific_enthalpy_j_kg;
        close(actual_stored + new.vapor_energy_j, 0.0, 1.0e-7);
        close(new.latent_energy_j, dm * LATENT_HEAT, 1.0e-6);
        close(new.recession_m, dm / old.solid.area_m2 / 1000.0, 1.0e-12);
        stored_change += actual_stored;
        carried += new.vapor_energy_j;
        latent += new.latent_energy_j;
    }
    close(stored_change + step.vapor_energy_j, 0.0, 2.0e-7);
    close(carried, step.vapor_energy_j, 1.0e-9);
    close(latent, step.latent_energy_j, 1.0e-9);
    close(
        vapor_mass(step.reservoir) - initial_mass,
        step.transferred_mass_kg,
        2.0e-12,
    );
    assert!(step.energy_residual_j.abs() <= 2.0e-7);
}
