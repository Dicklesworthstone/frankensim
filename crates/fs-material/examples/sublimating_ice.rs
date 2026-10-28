//! A cold ice film in a maintained vacuum and in a finite vapor chamber.
//!
//! Run: cargo run -p fs-material --example sublimating_ice
//!
//! Pressure and latent heat use WaterIceSublimation. Density 917 kg/m^3 and
//! constant solid cp=2000 J/(kg K) are explicit illustrative approximations,
//! not a sourced calorimetry dataset. The chart contains ONLY solid ice on
//! 200..260 K: a melting plateau is not an entirely sublimating solid inventory.
//!
//! The 10-micrometre film has uniform temperature, an insulated back face and
//! negligible edge area. Accommodation alpha=1 is the ideal kinetic limit;
//! there is no carrier-gas diffusion film. The finite chamber is well mixed and
//! held at 230 K. Its vapor mass/backpressure evolve; gas temperature, transport
//! gradients and geometry remapping are not solved. Recession is reported on
//! the fixed reference area. Halve MAX_DT_S to assess first-order coupling
//! error; conservation residuals alone do not establish physical accuracy.

use fs_material::sublimation::WaterIceSublimation;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::sublimation::{SolverConfig, SublimationLaw, SurfaceKinetics, VaporReservoir};

const INITIAL_T_K: f64 = 230.0;
const SOLID_CP_J_KG_K: f64 = 2000.0;
const SOLID_DENSITY_KG_M3: f64 = 917.0;
const FRONT_AREA_M2: f64 = 0.001;
const FILM_THICKNESS_M: f64 = 10.0e-6;
const CHAMBER_VOLUME_M3: f64 = 100.0e-6; // 100 mL
const MAX_DT_S: f64 = 25.0e-6;

struct SolidIceCaloric;
impl ThermalLaw for SolidIceCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = INITIAL_T_K + h / SOLID_CP_J_KG_K;
        if (200.0..=260.0).contains(&temperature) {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / SOLID_CP_J_KG_K
    }
}
static CALORIC: SolidIceCaloric = SolidIceCaloric;

fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: MAX_DT_S,
        energy_tolerance_j: 1.0e-8,
    }
}

fn initial_film() -> Result<DepositionSimulation<'static>, DepositionError> {
    let volume = FRONT_AREA_M2 * FILM_THICKNESS_M;
    let cell = DepositedCell {
        id: CellId(1),
        segment_id: 1,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: SOLID_DENSITY_KG_M3 * volume,
        specific_h_j_kg: 0.0, // Common energy datum: solid ice at 230 K.
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [FRONT_AREA_M2.sqrt(), 0.0, 0.0],
            volume_m3: volume,
            surface_area_m2: 2.0 * FRONT_AREA_M2,
        },
        boundary: Boundary {
            ambient_k: INITIAL_T_K,
            htc_w_m2_k: 0.0,
            source_w: 0.0,
            reservoir: Some(Reservoir {
                area_m2: FRONT_AREA_M2,
                conductance_w_k: 0.0,
                temperature_k: INITIAL_T_K,
            }),
        },
    };
    let timeline = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 10_000,
        },
    )?;
    timeline
        .advance_to(0.0, thermal_config(), || Ok(()))
        .map(|(state, _)| state)
}

struct Outcome {
    final_temperature_k: f64,
    transferred_mass_kg: f64,
}

fn run_case(
    label: &str,
    law: &WaterIceSublimation,
    mut vapor: VaporReservoir,
) -> Result<Outcome, Box<dyn std::error::Error>> {
    let mut thermal = initial_film()?;
    let initial_mass = thermal.mass_kg();
    let initial_energy = thermal.energy_j();
    let molar_mass = law.properties(INITIAL_T_K)?.molar_mass_kg_mol;
    let kinetics = SurfaceKinetics {
        accommodation_coefficient: 1.0,
        mass_transfer_coefficient_m_s: None,
    };
    let phase = SolverConfig {
        relative_tolerance: 1.0e-10,
        absolute_energy_tolerance_j: 1.0e-12,
        max_iterations: 96,
    };
    let mut transferred_mass_kg = 0.0;
    let mut vapor_energy_j = 0.0;
    let mut external_heat_j = 0.0;
    let mut recession_m = 0.0;

    println!("\n{label}");
    println!("time_ms  ice_K     ice_mg    vapor_Pa  loss_ug  recession_um  energy_residual_J");
    for sample in 1..=5 {
        let target_s = f64::from(sample) * 0.001;
        let (next, next_vapor, audit) = thermal.advance_sublimating_to(
            CellId(1),
            law,
            kinetics,
            vapor,
            target_s,
            thermal_config(),
            phase,
            || Ok(()),
        )?;
        thermal = next;
        vapor = next_vapor;
        transferred_mass_kg += audit.transferred_mass_kg;
        vapor_energy_j += audit.vapor_energy_j;
        external_heat_j += audit.thermal.external_heat_j;
        recession_m += audit.recession_m;
        let partial_pressure_pa = vapor.partial_pressure_pa(molar_mass)?;
        let energy_residual_j =
            thermal.energy_j() - initial_energy - external_heat_j + vapor_energy_j;
        let mass_residual_kg = thermal.mass_kg() - initial_mass + transferred_mass_kg;
        assert!(energy_residual_j.abs() < 1.0e-8, "energy did not close");
        assert!(mass_residual_kg.abs() < 1.0e-15, "mass did not close");
        println!(
            "{:7.3}  {:8.4}  {:8.5}  {:8.4}  {:7.3}  {:12.6}  {:18.3e}",
            target_s * 1.0e3,
            thermal.temperature_k(CellId(1))?,
            thermal.mass_kg() * 1.0e6,
            partial_pressure_pa,
            transferred_mass_kg * 1.0e9,
            recession_m * 1.0e6,
            energy_residual_j,
        );
    }
    println!("Vapor-carried energy (sensible plus latent): {vapor_energy_j:.8} J");
    println!("External sensible heat: {external_heat_j:.8} J");
    Ok(Outcome {
        final_temperature_k: thermal.temperature_k(CellId(1))?,
        transferred_mass_kg,
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let law = WaterIceSublimation::try_new(SOLID_DENSITY_KG_M3)?;
    println!("Sublimation of a 10-micrometre ice film initially at 230 K.");
    println!("Approximate solid cp={SOLID_CP_J_KG_K} J/(kg K); alpha=1; fixed reference area.");
    println!("The finite chamber stores vapor mass and is held at 230 K.");
    let vacuum = run_case(
        "Maintained vacuum (vapor partial pressure fixed at zero)",
        &law,
        VaporReservoir::Prescribed {
            partial_pressure_pa: 0.0,
            temperature_k: INITIAL_T_K,
        },
    )?;
    let chamber = run_case(
        "Finite 100-mL chamber (initially no water vapor)",
        &law,
        VaporReservoir::Finite {
            vapor_mass_kg: 0.0,
            volume_m3: CHAMBER_VOLUME_M3,
            temperature_k: INITIAL_T_K,
        },
    )?;
    assert!(vacuum.transferred_mass_kg > chamber.transferred_mass_kg);
    assert!(vacuum.final_temperature_k < chamber.final_temperature_k);
    println!(
        "\nAfter 5 ms: vacuum loses {:.3} micrograms versus {:.3} micrograms in the chamber.",
        vacuum.transferred_mass_kg * 1.0e9,
        chamber.transferred_mass_kg * 1.0e9,
    );
    println!("Vapor backpressure limits loss; sublimation cools the same finite solid inventory.");
    Ok(())
}
