//! Equal intercepted optical pulses heat bare, coated and transparent surfaces.
//!
//! Run: cargo run -p fs-render --example surface_optical_heating \
//!   --no-default-features --features optical-heating,material-state
//!
//! Three independent synthetic 1 g cells share the same authored caloric law.
//! Their optical responses are a polished opaque metal, that same metal under
//! a coherent transparent film, and a clear dielectric interface. Each receives
//! a 192 J blue box and a 160 J red triangle. The metals melt through the native
//! 350 K latent plateau; native convection refreezes them after switch-off.
//! The clear interface reflects/transmits its power without surface heating.
//!
//! These are explicit synthetic constitutive laws, not measured material data.
//! Fixed optical constants are declared valid from 300 to 600 K; no phase/index
//! feedback is inferred. Coating and opaque-substrate loss heat one receiving
//! cell, with no claim about individual layer temperatures. Incident watts are
//! already intercepted by each surface: spatial/angular quadrature, coherence
//! and sufficient substrate optical depth remain caller assumptions. Radiation,
//! expansion and volume absorption behind the clear interface are not modeled.

use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_render::conductor::{
    CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorOptics,
    ConductorSource,
};
use fs_render::dielectric::CauchyIor;
use fs_render::thermal_surface::heating::{
    SurfaceHeatingBudget, SurfaceHeatingTolerance, SurfaceIlluminationSource,
    SurfaceOpticalHeating, SurfacePowerResponse,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
};
use std::error::Error;

const MATERIAL: MaterialId = MaterialId(1);
const CELLS: [CellId; 3] = [CellId(1), CellId(2), CellId(3)];
const WAVELENGTHS_NM: [f64; 2] = [480.0, 680.0];
const SWITCH_OFF_S: f64 = 8.5;

fn boundary(htc_w_m2_k: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k,
        reservoir: None,
        source_w: 0.0,
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let caloric = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([181; 32]),
        [
            (0.0, 300.0, 0.0),
            (50_000.0, 350.0, 0.0),
            (150_000.0, 350.0, 1.0),
            (400_000.0, 600.0, 1.0),
        ]
        .map(
            |(specific_enthalpy_j_kg, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1000.0,
            },
        )
        .to_vec(),
    )?;
    let mut events = Vec::new();
    for id in CELLS {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: id.0,
                body_id: id.0,
                material_id: MATERIAL,
                mass_kg: 0.001,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.005, 0.005],
                    end_m: [0.01, 0.005, 0.005],
                    volume_m3: 1e-6,
                    surface_area_m2: 6e-4,
                },
                boundary: boundary(10.0),
            }),
        });
    }
    for id in CELLS {
        events.push(TimedEvent {
            time_s: SWITCH_OFF_S,
            event: Event::Boundary(id, boundary(2000.0)),
        });
    }
    let step = StepConfig {
        max_dt_s: 0.05,
        energy_tolerance_j: 1e-8,
    };
    let mut thermal = DepositionSimulation::new(
        events,
        &[(MATERIAL, &caloric)],
        DepositionBudget {
            max_cells: 3,
            max_links: 0,
            max_events: 6,
            max_substeps: 1000,
        },
    )?
    .advance_to(0.0, step, || Ok(()))?
    .0;

    let substrate = ConductorOptics::try_new(
        [
            380.0, 430.0, 480.0, 530.0, 580.0, 630.0, 680.0, 730.0, 780.0,
        ]
        .map(|wavelength| ConductorIorSample::try_new(wavelength, 1.0, 2.0))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "wrong number of authored conductor knots")?,
        ConductorSource::try_new(ContentHash([182; 32]), ConductorDataStatus::Representative)?,
    )?;
    let coated = substrate.with_coating(CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(1.4)?,
        100e-9,
        ConductorSource::try_new(ContentHash([183; 32]), ConductorDataStatus::Representative)?,
    )?);
    let responses = [
        SurfacePowerResponse::from_conductor(WAVELENGTHS_NM, &substrate, [1.0; 2], 1.0, [0.5; 2])?,
        SurfacePowerResponse::from_conductor(WAVELENGTHS_NM, &coated, [1.0; 2], 1.0, [0.5; 2])?,
        SurfacePowerResponse::from_dielectric(
            WAVELENGTHS_NM,
            [1.0; 2],
            [1.5; 2],
            1.0,
            [0.5; 2],
            ContentHash([184; 32]),
        )?,
    ];
    let power = SpectralPowerSchedule::try_new(
        &[
            SpectralPowerSegment {
                start_s: 0.5,
                end_s: 4.5,
                start_power_w: [48.0, 0.0],
                end_power_w: [48.0, 0.0],
            },
            SpectralPowerSegment {
                start_s: 4.5,
                end_s: 6.5,
                start_power_w: [0.0; 2],
                end_power_w: [0.0, 80.0],
            },
            SpectralPowerSegment {
                start_s: 6.5,
                end_s: SWITCH_OFF_S,
                start_power_w: [0.0, 80.0],
                end_power_w: [0.0; 2],
            },
        ],
        3,
    )?;
    let sources: Vec<_> = CELLS
        .into_iter()
        .zip(responses)
        .map(|(id, response)| SurfaceIlluminationSource {
            source_id: id.0,
            thermal_cell: id,
            material_id: MATERIAL,
            temperature_bounds_k: (300.0, 600.0),
            response,
            power: power.clone(),
        })
        .collect();
    let mut illumination = SurfaceOpticalHeating::try_new(
        &sources,
        &thermal,
        SurfaceHeatingBudget {
            max_sources: 3,
            max_total_power_segments: 9,
            max_source_evaluations: 3000,
        },
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )?;
    eprintln!("Synthetic surface comparison: each 1 g cell receives 192 J blue + 160 J red.");
    eprintln!("Normal incidence, unpolarized light; n+ik=1+2i; coating n=1.4, d=100 nm.");
    eprintln!(
        "At 8.5 s illumination ends and native convection increases; the thermal clock continues."
    );
    for (name, response) in ["bare", "coated", "clear"].into_iter().zip(responses) {
        eprintln!(
            "{name}: R={:?}, T={:?}, A={:?}",
            response.reflected_fraction(),
            response.transmitted_fraction(),
            response.absorbed_fraction()
        );
    }
    println!(
        "time_s,bare_temperature_k,coated_temperature_k,clear_temperature_k,bare_liquid_fraction,coated_liquid_fraction,incident_j,reflected_j,transmitted_j,absorbed_j,convection_to_material_j,stored_energy_j,energy_residual_j"
    );
    let initial_energy_j = thermal.energy_j();
    let mut convection_j = 0.0;
    let mut melting = [false; 2];
    let mut melted = [false; 2];
    let mut freezing = [false; 2];
    let mut peak_temperature_k = [300.0_f64; 2];
    let mut switch_off_absorbed_j = None;
    for output in 1..=40 {
        let target_s = 0.5 * f64::from(output);
        let prior_optical_j = illumination.energy().thermal_input_j;
        let (next, accepted, audit) =
            thermal.advance_observed(&illumination, target_s, step, || Ok(()))?;
        thermal = next;
        illumination = accepted;
        let energy = illumination.energy();
        convection_j += audit.external_heat_j - (energy.thermal_input_j - prior_optical_j);
        let stored_j = thermal.energy_j() - initial_energy_j;
        let residual_j = stored_j - energy.thermal_input_j - convection_j;
        if residual_j.abs() > 1e-8 || audit.residual_j.abs() > 1e-8 {
            return Err("surface optical/thermal balance exceeded its declared tolerance".into());
        }
        let mut temperatures = [0.0; 3];
        let mut liquid = [0.0; 2];
        for (index, id) in CELLS.into_iter().enumerate() {
            let phase = caloric.state_at_specific_enthalpy(thermal.cells()[&id].specific_h_j_kg)?;
            temperatures[index] = phase.temperature_k();
            if index < 2 {
                liquid[index] = phase.liquid_mass_fraction();
                peak_temperature_k[index] = peak_temperature_k[index].max(temperatures[index]);
                melted[index] |= liquid[index] == 1.0;
                if liquid[index] > 0.0 && liquid[index] < 1.0 {
                    melting[index] |= target_s < SWITCH_OFF_S;
                    freezing[index] |= target_s > SWITCH_OFF_S;
                }
            }
        }
        if temperatures[2] != 300.0 {
            return Err("lossless transparent interface incorrectly supplied surface heat".into());
        }
        if target_s == SWITCH_OFF_S {
            switch_off_absorbed_j = Some(energy.thermal_input_j);
        } else if target_s > SWITCH_OFF_S && switch_off_absorbed_j != Some(energy.thermal_input_j) {
            return Err("incident surface heat changed during the dark cooling interval".into());
        }
        println!(
            "{target_s:.2},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{convection_j:.9},{stored_j:.9},{residual_j:.9e}",
            temperatures[0],
            temperatures[1],
            temperatures[2],
            liquid[0],
            liquid[1],
            energy.incident_j.iter().sum::<f64>(),
            energy.reflected_j.iter().sum::<f64>(),
            energy.transmitted_j.iter().sum::<f64>(),
            energy.thermal_input_j
        );
    }
    let energy = illumination.energy();
    for (lane, per_cell_fluence) in [192.0, 160.0].into_iter().enumerate() {
        let expected_absorbed = per_cell_fluence
            * responses
                .iter()
                .map(|response| response.absorbed_fraction()[lane])
                .sum::<f64>();
        if (energy.incident_j[lane] - 3.0 * per_cell_fluence).abs() > 1e-8
            || (energy.absorbed_j[lane] - expected_absorbed).abs() > 1e-8
        {
            return Err(
                "surface pulse fluence disagreed with the independent box/triangle areas".into(),
            );
        }
    }
    for index in 0..2 {
        let phase =
            caloric.state_at_specific_enthalpy(thermal.cells()[&CELLS[index]].specific_h_j_kg)?;
        if !melting[index]
            || !melted[index]
            || !freezing[index]
            || phase.liquid_mass_fraction() != 0.0
        {
            return Err("a metal receiver did not complete melting and dark refreezing".into());
        }
    }
    if peak_temperature_k[1] <= peak_temperature_k[0] {
        return Err(
            "the authored coating did not produce its expected increased absorption".into(),
        );
    }
    eprintln!(
        "Peak temperature: bare={:.6} K, coated={:.6} K; both refroze. Clear interface stayed at 300 K.",
        peak_temperature_k[0], peak_temperature_k[1]
    );
    Ok(())
}
