//! Coherent crossing waves -> original porous loss -> native cell temperatures.
//!
//! Run `cargo run -p fs-couple --example coherent_porous_heating -- [phase_rad
//! [peak_pressure_pa [quadrature_panels]]]`. Defaults: 0, 30, 32. Compare phase
//! zero with pi, or refine 32 -> 64 -> 128 panels. CSV goes to stdout. The finite
//! window retains coherent cross-angle heat, including transverse viscous heat
//! near pressure nodes; it is not an incoherent sum of the two incident powers.
//!
//! Basotect supplies ONLY the original acoustic card. Density, Cp, conductivity
//! and surface heat transfer below are explicitly authored demonstration inputs,
//! NOT a Basotect thermal characterization. All cells use the original thermal
//! storage/contact/boundary owner. Their finite area and volume enter heat once.
//! The acoustic coefficients stay frozen at 293.15 K: this is one-way cycle-mean
//! heating, not temperature feedback or a finite-aperture scattering solution.

use fs_couple::porous_heating::PorousHeatCell;
use fs_couple::porous_heating::coherent::{
    AcousticHeatNode, CoherentHeatQuadrature, CoherentHeatingBudget, CoherentLayerHeating,
    CoherentThermalColumn,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::layers::{
    AngularSpectrumBudget, IncidentPlaneWave, Layer, TransmittingLayers,
};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::{c64::C64, dd::Dd, det};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use std::io::{BufWriter, Write};

// Explicit illustrative thermal inputs, independent of the porous acoustic card.
const DENSITY_KG_M3: f64 = 20.0;
const CP_J_KG_K: f64 = 1000.0;
const CONDUCTIVITY_W_M_K: f64 = 0.03;
const AMBIENT_HTC_W_M2_K: f64 = 5.0;
const REFERENCE_K: f64 = 293.15;
const LENGTH_M: f64 = 0.06;
const CELLS: usize = 4;

struct AuthoredCaloric;
impl ThermalLaw for AuthoredCaloric {
    fn temperature_k(&self, energy_j_kg: f64) -> Result<f64, DepositionError> {
        let temperature = REFERENCE_K + energy_j_kg / CP_J_KG_K;
        if temperature.is_finite() && temperature > 0.0 {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP_J_KG_K
    }
}

fn midpoint(panels: usize) -> Vec<AcousticHeatNode> {
    (0..panels)
        .map(|i| AcousticHeatNode {
            coordinate: (i as f64 + 0.5) / panels as f64,
            weight: 1.0 / panels as f64,
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let phase: f64 = args.next().map_or(Ok(0.0), |s| s.parse())?;
    let peak: f64 = args.next().map_or(Ok(30.0), |s| s.parse())?;
    let panels: usize = args.next().map_or(Ok(32), |s| s.parse())?;
    if args.next().is_some()
        || !phase.is_finite()
        || !peak.is_finite()
        || peak < 0.0
        || !(1..=256).contains(&panels)
    {
        return Err(
            "expected finite phase, nonnegative peak pressure and 1..=256 quadrature panels".into(),
        );
    }
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), REFERENCE_K, 101_325.0)?;
    let porous = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), gas.dynamic_viscosity)?;
    let stack = TransmittingLayers::new(
        &gas,
        &gas,
        core::f64::consts::TAU * 1000.0,
        &[Layer::Porous {
            law: &porous,
            thickness_m: LENGTH_M,
        }],
    )?;
    let (phase_sin, phase_cos) = det::sin_cos(phase);
    let beam = stack.prepare_angular_spectrum(
        &[
            IncidentPlaneWave {
                direction: [0.6, 0.0, 0.8],
                pressure_pa: C64::from_re(peak),
            },
            IncidentPlaneWave {
                direction: [-0.6, 0.0, 0.8],
                pressure_pa: C64::new(phase_cos, phase_sin).scale(peak),
            },
        ],
        AngularSpectrumBudget { max_components: 2 },
        &mut || false,
    )?;
    let field = beam.prepare_layer_field(0, &mut || false)?;
    let column = CoherentThermalColumn {
        origin_m: [0.0; 3],
        transverse_axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        lower_m: [-0.02, -0.01],
        extent_m: [0.04, 0.02],
    };
    let area = column.extent_m[0] * column.extent_m[1];
    let width = LENGTH_M / CELLS as f64;
    let mut events = Vec::new();
    let mut routes = Vec::new();
    for i in 0..CELLS {
        let start = i as f64 * width;
        let end = if i + 1 == CELLS {
            LENGTH_M
        } else {
            (i + 1) as f64 * width
        };
        let cell = CellId(i as u64 + 1);
        routes.push(PorousHeatCell {
            cell,
            end_depth_m: end,
        });
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cell,
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: DENSITY_KG_M3 * area * (end - start),
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, start],
                    end_m: [0.0, 0.0, end],
                    volume_m3: area * (end - start),
                    surface_area_m2: 2.0 * area
                        + 2.0 * (column.extent_m[0] + column.extent_m[1]) * (end - start),
                },
                boundary: Boundary {
                    ambient_k: REFERENCE_K,
                    htc_w_m2_k: AMBIENT_HTC_W_M2_K,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    for i in 1..CELLS {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                area_m2: area,
                conductance_w_k: CONDUCTIVITY_W_M_K * area / width,
            }),
        });
    }
    let caloric = AuthoredCaloric;
    let config = StepConfig {
        max_dt_s: 0.02,
        energy_tolerance_j: 1e-10,
    };
    let mut thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &caloric)],
        DepositionBudget {
            max_cells: CELLS,
            max_links: CELLS - 1,
            max_events: 2 * CELLS,
            max_substeps: 10_000,
        },
    )?
    .advance_to(0.0, config, || Ok(()))?
    .0;
    let xy_depth = midpoint(panels);
    // Both waves are uniform along y; one midpoint integrates that axis exactly.
    let y = midpoint(1);
    let heating = CoherentLayerHeating::try_new(
        &field,
        &thermal,
        column,
        &routes,
        CoherentHeatQuadrature {
            axes: [&xy_depth, &y, &xy_depth],
        },
        CoherentHeatingBudget {
            max_cells: CELLS,
            max_component_evaluations: 1_000_000,
        },
        || Ok(()),
    )?;
    let mut cell_j = [Dd::ZERO; CELLS];
    let mut acoustic_j = Dd::ZERO;
    let mut external_j = Dd::ZERO;
    let stdout = std::io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    writeln!(
        output,
        "time_s,cell_id,near_depth_m,far_depth_m,temperature_k,viscous_power_w,thermal_power_w,cell_acoustic_j,total_acoustic_j,net_external_heat_j,stored_energy_j"
    )?;
    for second in 0..=20 {
        if second > 0 {
            let step = heating.advance_to(&thermal, f64::from(second), config, || Ok(()))?;
            for ((sum, route), (id, energy)) in
                cell_j.iter_mut().zip(&routes).zip(step.cell_energy_j)
            {
                debug_assert_eq!(route.cell, id);
                *sum = *sum + Dd::from_f64(energy);
            }
            acoustic_j = acoustic_j + Dd::from_f64(step.acoustic.deposited_j);
            external_j = external_j + Dd::from_f64(step.heat_report.external_heat_j);
            thermal = step.thermal;
        }
        for (cell, energy) in heating.cells().iter().zip(cell_j) {
            writeln!(
                output,
                "{},{},{:.9e},{:.9e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                second,
                cell.cell.0,
                cell.start_depth_m,
                cell.end_depth_m,
                thermal.temperature_k(cell.cell)?,
                cell.dissipation.viscous_power_w,
                cell.dissipation.thermal_power_w,
                energy.to_f64(),
                acoustic_j.to_f64(),
                external_j.to_f64(),
                thermal.energy_j()
            )?;
        }
    }
    output.flush()?;
    Ok(())
}
