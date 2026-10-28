//! Two ice films exchange vapor through one initially empty cold-trap chamber.
//!
//! Run: cargo run -p fs-material --example ice_cold_trap
//!
//! WaterIceSublimation supplies vapor pressure and latent heat. Constant
//! cp=2000 J/(kg K) and density=917 kg/m^3 are illustrative caller approximations.
//! Both caloric charts describe ONLY pure solid ice between 200 and 260 K.
//! Each 10-micrometre film has one exposed face and an insulated back; edges
//! are neglected. Alpha=1 assumes ideal accommodation with no diffusion film.
//! The chamber is well mixed and thermostatted at 230 K. Gas temperature and
//! gas caloric energy are not evolved; vapor-carried energy is accounted at
//! the solid boundaries. Geometry stays fixed; signed recession is a diagnostic.
//! Heat/phase splitting and backward Euler are first order: halve MAX_DT_S to
//! assess time error; conservation alone does not establish convergence.

use fs_material::sublimation::WaterIceSublimation;
use fs_thermal::enthalpy::sublimation::chamber::ChamberSurface;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::sublimation::chamber::ChamberSolverConfig;
use fs_thermal::sublimation::{SolverConfig, SublimationLaw, SurfaceKinetics, VaporReservoir};

const DATUM_K: f64 = 230.0;
const CP_J_KG_K: f64 = 2000.0;
const DENSITY_KG_M3: f64 = 917.0;
const FRONT_AREA_M2: f64 = 0.001;
const THICKNESS_M: f64 = 10.0e-6;
const CHAMBER_VOLUME_M3: f64 = 100.0e-6;
const MAX_DT_S: f64 = 25.0e-6;
const DURATION_S: f64 = 0.005;

struct SolidIceCaloric;
impl ThermalLaw for SolidIceCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = DATUM_K + h / CP_J_KG_K;
        if (200.0..=260.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP_J_KG_K
    }
}
static CALORIC: SolidIceCaloric = SolidIceCaloric;

fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: MAX_DT_S,
        energy_tolerance_j: 1.0e-8,
    }
}

fn initial_films() -> Result<DepositionSimulation<'static>, DepositionError> {
    let volume = FRONT_AREA_M2 * THICKNESS_M;
    let events = [(1, 240.0), (2, 220.0)].map(|(id, t)| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: id,
            material_id: MaterialId(1),
            mass_kg: DENSITY_KG_M3 * volume,
            specific_h_j_kg: CP_J_KG_K * (t - DATUM_K),
            geometry: CellGeometry {
                start_m: [id as f64 * 0.1, 0.0, 0.0],
                end_m: [id as f64 * 0.1, FRONT_AREA_M2.sqrt(), 0.0],
                volume_m3: volume,
                surface_area_m2: 2.0 * FRONT_AREA_M2,
            },
            boundary: Boundary {
                ambient_k: DATUM_K,
                htc_w_m2_k: 0.0,
                source_w: 0.0,
                reservoir: Some(Reservoir {
                    area_m2: FRONT_AREA_M2,
                    conductance_w_k: 0.0,
                    temperature_k: DATUM_K,
                }),
            },
        }),
    });
    DepositionSimulation::new(
        events.to_vec(),
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 10_000,
        },
    )?
    .advance_to(0.0, thermal_config(), || Ok(()))
    .map(|(state, _)| state)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let law = WaterIceSublimation::try_new(DENSITY_KG_M3)?;
    let initial = initial_films()?;
    let surfaces = [1, 2].map(|id| ChamberSurface {
        cell_id: CellId(id),
        kinetics: SurfaceKinetics {
            accommodation_coefficient: 1.0,
            mass_transfer_coefficient_m_s: None,
        },
    });
    let phase = ChamberSolverConfig {
        surface_solver: SolverConfig {
            relative_tolerance: 1.0e-10,
            absolute_energy_tolerance_j: 1.0e-12,
            max_iterations: 96,
        },
        relative_mass_tolerance: 1.0e-10,
        absolute_mass_tolerance_kg: 1.0e-17,
        max_iterations: 96,
        max_surfaces: 2,
    };
    let (final_state, reservoir, audit) = initial.advance_chamber_sublimating_to(
        &surfaces,
        &law,
        VaporReservoir::Finite {
            vapor_mass_kg: 0.0,
            volume_m3: CHAMBER_VOLUME_M3,
            temperature_k: DATUM_K,
        },
        DURATION_S,
        thermal_config(),
        phase,
        || Ok(()),
    )?;
    let VaporReservoir::Finite { vapor_mass_kg, .. } = reservoir else {
        unreachable!("the shared chamber remains finite");
    };
    let pressure_pa = reservoir.partial_pressure_pa(law.properties(DATUM_K)?.molar_mass_kg_mol)?;
    let warm_loss = initial.cells()[&CellId(1)].mass_kg - final_state.cells()[&CellId(1)].mass_kg;
    let cold_gain = final_state.cells()[&CellId(2)].mass_kg - initial.cells()[&CellId(2)].mass_kg;
    let mass_residual = final_state.mass_kg() + vapor_mass_kg - initial.mass_kg();
    // Vapor energy already contains latent AND transported sensible energy.
    let energy_residual =
        final_state.energy_j() - initial.energy_j() - audit.thermal.external_heat_j
            + audit.vapor_energy_j;

    assert!(
        warm_loss > 0.0 && cold_gain > 0.0,
        "cold ice must capture newly emitted vapor"
    );
    assert!(vapor_mass_kg > 0.0);
    assert!(
        mass_residual.abs() < 1.0e-15,
        "combined solid/gas mass did not close"
    );
    assert!(energy_residual.abs() < 1.0e-8, "energy did not close");

    println!("Two 10-micrometre ice films, 5 ms, initially empty 100-mL chamber at 230 K.");
    println!("Illustrative solid cp={CP_J_KG_K} J/(kg K), density={DENSITY_KG_M3} kg/m^3.");
    println!("surface  initial_K  final_K   mass_mg   signed_recession_um");
    for (id, label, start_k) in [(CellId(1), "warm", 240.0), (CellId(2), "cold", 220.0)] {
        let cell = &final_state.cells()[&id];
        let recession = audit.cells[&id].recession_m;
        assert!(cell.mass_kg > 0.0 && recession.abs() < 0.05 * THICKNESS_M);
        println!(
            "{label:7}  {start_k:9.3}  {:7.3}  {:8.5}  {:19.6}",
            final_state.temperature_k(id)?,
            cell.mass_kg * 1.0e6,
            recession * 1.0e6
        );
    }
    println!(
        "Warm loss: {:.6} micrograms; cold frost gain: {:.6} micrograms.",
        warm_loss * 1.0e9,
        cold_gain * 1.0e9
    );
    println!(
        "Shared vapor: {:.6} micrograms; pressure: {pressure_pa:.6} Pa.",
        vapor_mass_kg * 1.0e9
    );
    println!("Solid + vapor mass residual: {mass_residual:.3e} kg.");
    println!(
        "Energy residual: {energy_residual:.3e} J; vapor-carried energy: {:.9} J.",
        audit.vapor_energy_j
    );
    println!("Negative recession denotes frost growth; gas temperature is held fixed.");
    Ok(())
}
