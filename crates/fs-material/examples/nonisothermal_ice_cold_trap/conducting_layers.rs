//! One-dimensional dense-ice layers with evolving thermal resistance.
//!
//! Reuses the parent example's NASA-9/IAPWS owners and approximate ice calorics.
//! Each film is one thermal cell: a lower face bath, an upper molecular face,
//! and a normal conduction distance derived from conserved mass. Interface
//! temperature is represented by the cell temperature. This does not establish
//! spatial convergence, porous-frost behavior, nucleation or experimental validity.

use super::*;
use fs_material::sublimation::gas::timeline::{
    GasTimelineSurface, observed::advance_deposition_with_gas,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, TimedEvent,
    conduction::{ConductionObserver, ConductiveReservoir, ConductivityLaw, PlanarPhaseLayer},
};

const AREA_M2: f64 = 1e-3;
const THICKNESS_M: f64 = 1e-5;
const DENSITY_KG_M3: f64 = 917.0;
const CONDUCTIVITY_W_M_K: f64 = 2.2;
const DURATION_S: f64 = 0.002;
const BATHS_K: [f64; 2] = [260.0, 252.0];

/// Conservative admitted envelope for this fixed 250--270 K derived chart.
/// Actual cp is checked whenever the original caloric inverse is queried.
const MINIMUM_HEAT_CAPACITY_J_KG_K: f64 = 1000.0;
struct ThermalIce<'a, 'source>(&'a IceCaloric<'source>);
impl ThermalLaw for ThermalIce<'_, '_> {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = self
            .0
            .temperature(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let cp = self
            .0
            .heat_capacity(temperature)
            .map_err(|_| DepositionError::MaterialDomain)?;
        if cp < MINIMUM_HEAT_CAPACITY_J_KG_K {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(temperature)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / MINIMUM_HEAT_CAPACITY_J_KG_K
    }
}

// Authored dense-ice scenario parameter; not a new measured source fit.
struct Conductivity;
impl ConductivityLaw for Conductivity {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(CONDUCTIVITY_W_M_K)
    }
}

#[derive(Clone, Copy)]
struct Report {
    dt_s: f64,
    temperatures_k: [f64; 3],
    pressure_pa: f64,
    warm_loss_kg: f64,
    cold_gain_kg: f64,
    thicknesses_m: [f64; 2],
    sampled_half_distances_m: [f64; 2],
    bath_heat_j: f64,
    mass_residual_kg: f64,
    energy_residual_j: f64,
    thermal_steps: usize,
    phase_steps: usize,
}

#[allow(clippy::too_many_arguments)]
fn simulate(
    thermal: &DepositionSimulation<'_>,
    surfaces: &[GasTimelineSurface<'_>],
    observer: &ConductionObserver<'_>,
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    water: &SpeciesId,
    caloric: &IceCaloric<'_>,
    phase: GasChamberSolve,
    dt_s: f64,
) -> Result<Report, Box<dyn Error>> {
    let (candidate, conduction) = advance_deposition_with_gas(
        thermal,
        surfaces,
        None,
        observer,
        gas,
        catalog,
        DURATION_S,
        StepConfig {
            max_dt_s: dt_s,
            energy_tolerance_j: 1e-8,
        },
        phase,
        || Ok(()),
    )?;
    let final_thermal = candidate.thermal();
    let final_gas = candidate.gas();
    let initial = [
        thermal.cells()[&CellId(10)].clone(),
        thermal.cells()[&CellId(20)].clone(),
    ];
    let endpoint = [
        &final_thermal.cells()[&CellId(10)],
        &final_thermal.cells()[&CellId(20)],
    ];
    let gas_mass_change = final_gas.mass_kg() - gas.mass_kg();
    let mut mass_residual = gas_mass_change;
    let mut solid_energy_change = 0.0;
    for (old, new) in initial.iter().zip(endpoint) {
        let dm = new.mass_kg - old.mass_kg;
        mass_residual += dm;
        solid_energy_change += new.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg)
            + dm * (old.specific_h_j_kg + caloric.datum);
    }
    let initial_u = gas.state().specific_internal_energy_j_kg();
    let gas_energy_change = final_gas.mass_kg()
        * (final_gas.state().specific_internal_energy_j_kg() - initial_u)
        + gas_mass_change * initial_u;
    let bath_heat_j = candidate.audit().thermal.external_heat_j;
    let energy_residual = solid_energy_change + gas_energy_change - bath_heat_j;
    let warm_loss = initial[0].mass_kg - endpoint[0].mass_kg;
    let cold_gain = endpoint[1].mass_kg - initial[1].mass_kg;
    assert!(warm_loss > 0.0 && cold_gain > 0.0);
    assert!(
        mass_residual.abs() <= 1e-14,
        "coupled mass residual: {mass_residual:e}"
    );
    assert!(
        energy_residual.abs() <= 1e-8,
        "coupled energy residual: {energy_residual:e}"
    );
    assert!(candidate.audit().mass_residual_kg.abs() <= 1e-14);
    assert!(candidate.audit().energy_residual_j.abs() <= 1e-8);
    let thicknesses = endpoint.map(|cell| cell.mass_kg / DENSITY_KG_M3 / AREA_M2);
    let sampled = [
        conduction
            .last_planar_half_distance_m(CellId(10))
            .ok_or("missing warm-layer distance")?,
        conduction
            .last_planar_half_distance_m(CellId(20))
            .ok_or("missing cold-layer distance")?,
    ];
    assert!(sampled[0] < 0.5 * THICKNESS_M && sampled[1] > 0.5 * THICKNESS_M);
    assert!(thicknesses[0] < THICKNESS_M && thicknesses[1] > THICKNESS_M);
    assert_eq!(endpoint[0].geometry, initial[0].geometry);
    assert_eq!(endpoint[1].geometry, initial[1].geometry);
    Ok(Report {
        dt_s,
        temperatures_k: [
            caloric.temperature(endpoint[0].specific_h_j_kg)?,
            caloric.temperature(endpoint[1].specific_h_j_kg)?,
            final_gas.state().state().temperature,
        ],
        pressure_pa: final_gas.species_partial_pressure_pa(water)?,
        warm_loss_kg: warm_loss,
        cold_gain_kg: cold_gain,
        thicknesses_m: thicknesses,
        sampled_half_distances_m: sampled,
        bath_heat_j,
        mass_residual_kg: mass_residual,
        energy_residual_j: energy_residual,
        thermal_steps: candidate.audit().thermal.substeps,
        phase_steps: candidate.audit().phase_steps,
    })
}

pub(super) fn run() -> Result<(), Box<dyn Error>> {
    let scale = std::env::args()
        .nth(2)
        .map(|value| value.parse::<f64>())
        .transpose()?
        .unwrap_or(1.0);
    if !scale.is_finite() || scale <= 0.0 {
        return Err("step scale must be finite and positive".into());
    }
    let water = SpeciesId::new("H2O")?;
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![water_vapor(&water)?],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0])?,
        )?,
        vec![transport(&water)?],
        Some(&water),
    )?;
    let law = WaterIceSublimation::try_new(DENSITY_KG_M3)?;
    let datum = catalog.species_specific_enthalpy_j_kg(&water, DATUM_T)?
        - law.properties(DATUM_T)?.latent_heat_j_kg;
    let caloric = IceCaloric {
        law: &law,
        vapor: &catalog,
        species: &water,
        datum,
    };
    for i in 0..=200 {
        let t = MIN_T + (MAX_T - MIN_T) * f64::from(i) / 200.0;
        assert!(caloric.heat_capacity(t)? >= MINIMUM_HEAT_CAPACITY_J_KG_K);
        assert!((caloric.temperature(caloric.enthalpy(t)?)? - t).abs() < 1e-8);
    }
    let thermal_law = ThermalIce(&caloric);
    let binding = GasSurfaceBinding {
        species: water.clone(),
        accommodation_coefficient: 1.0,
        solid_enthalpy_offset_j_kg: datum,
    };
    let mut events = Vec::new();
    for (id, temperature) in [10_u64, 20].into_iter().zip(BATHS_K) {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(1),
                mass_kg: DENSITY_KG_M3 * AREA_M2 * THICKNESS_M,
                specific_h_j_kg: caloric.enthalpy(temperature)?,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64, THICKNESS_M, 0.0],
                    volume_m3: AREA_M2 * THICKNESS_M,
                    surface_area_m2: 2.0 * AREA_M2,
                },
                boundary: Boundary {
                    ambient_k: DATUM_T,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: Some(Reservoir {
                        area_m2: AREA_M2,
                        conductance_w_k: 0.0,
                        temperature_k: temperature,
                    }),
                },
            }),
        });
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &thermal_law)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 100_000,
        },
    )?
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DURATION_S,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )?
    .0;
    let conductivity = Conductivity;
    let faces = [10_u64, 20].map(|id| ConductiveReservoir {
        cell_id: CellId(id),
        half_distance_m: 0.5 * THICKNESS_M,
        surface_resistance_m2_k_w: 0.0,
    });
    let layers = [10_u64, 20].map(|id| PlanarPhaseLayer {
        cell_id: CellId(id),
        footprint_area_m2: AREA_M2,
        density_kg_m3: DENSITY_KG_M3,
    });
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(1), &conductivity)],
        &[],
        &faces,
        || Ok(()),
    )?
    .with_planar_phase_layers(&thermal, &layers, || Ok(()))?;
    let range = [caloric.enthalpy(MIN_T)?, caloric.enthalpy(MAX_T)?];
    let surfaces = [10_u64, 20].map(|id| GasTimelineSurface {
        cell_id: CellId(id),
        law: &law,
        binding: &binding,
        pure_solid_enthalpy_range_j_kg: range,
    });
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            DATUM_T,
        )?
        .with_quantity("pressure", typed(QuantityKind::Pressure), 130.0)?;
    let gas = ClosedGasParcel::new(catalog.clone(), point, 100e-6)?;
    let phase = GasChamberSolve {
        surface: GasSurfaceSolve {
            gas: GasSpeciesBoundarySolve {
                gas: GasEnergySolve {
                    minimum_temperature_k: MIN_T,
                    maximum_temperature_k: MAX_T,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
                mass_tolerance_kg: 1e-14,
                energy_tolerance_j: 1e-8,
            },
            maximum_gross_inventory_fraction: 0.1,
            latent_compatibility_tolerance_j_kg: 1e-5,
        },
        maximum_surfaces: 2,
    };
    let mut reports = Vec::new();
    for factor in [1.0, 0.5, 0.25] {
        reports.push(simulate(
            &thermal,
            &surfaces,
            &conduction,
            &gas,
            &catalog,
            &water,
            &caloric,
            phase,
            20e-6 * scale * factor,
        )?);
    }
    println!("Conducting dense-ice layers: {DURATION_S:.6} s; fixed gas volume 100 mL");
    println!(
        "Each layer: footprint {AREA_M2:e} m2, initial thickness {:.3} micrometers, density {DENSITY_KG_M3} kg/m3",
        THICKNESS_M * 1e6
    );
    println!(
        "Explicit scenario k={CONDUCTIVITY_W_M_K} W/(m K); bottom face baths 260/252 K; one vapor face each"
    );
    println!(
        "NASA-9/IAPWS derived caloric approximation; one thermal cell per layer, no frost porosity or gas-volume work"
    );
    for row in &reports {
        println!(
            "dt={:.3} us: warm loss={:.6} ug, cold gain={:.6} ug; Tw/Tc/Tg={:.6}/{:.6}/{:.6} K, pv={:.6} Pa",
            row.dt_s * 1e6,
            row.warm_loss_kg * 1e9,
            row.cold_gain_kg * 1e9,
            row.temperatures_k[0],
            row.temperatures_k[1],
            row.temperatures_k[2],
            row.pressure_pa
        );
        println!(
            "  final thicknesses={:.6}/{:.6} um; last heat-path half-thicknesses={:.6}/{:.6} um",
            row.thicknesses_m[0] * 1e6,
            row.thicknesses_m[1] * 1e6,
            row.sampled_half_distances_m[0] * 1e6,
            row.sampled_half_distances_m[1] * 1e6
        );
        println!(
            "  bath heat={:.9} J; mass defect={:.3e} kg; independent energy defect={:.3e} J; thermal/phase steps={}/{}",
            row.bath_heat_j,
            row.mass_residual_kg,
            row.energy_residual_j,
            row.thermal_steps,
            row.phase_steps
        );
    }
    for pair in reports.windows(2) {
        let temperature_change = pair[0]
            .temperatures_k
            .into_iter()
            .zip(pair[1].temperatures_k)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        println!(
            "Refinement {:.3}->{:.3} us: max endpoint temperature change={:.6} K; cold gain change={:.6} ug",
            pair[0].dt_s * 1e6,
            pair[1].dt_s * 1e6,
            temperature_change,
            (pair[0].cold_gain_kg - pair[1].cold_gain_kg).abs() * 1e9
        );
    }
    println!(
        "Refinement differences estimate temporal sensitivity only, not spatial or experimental validation."
    );
    Ok(())
}
