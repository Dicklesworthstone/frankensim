//! Phase-dependent surface absorption versus a frozen-optics control.
//!
//! Run: cargo run -p fs-render --example phase_surface_optical_heating \
//!   --no-default-features --features optical-heating,material-state
//!
//! Two synthetic 1 g cells share one native enthalpy chart and thermal clock.
//! Both receive 320 J blue and 200 J red pulses. One uses explicitly authored
//! complex indices versus the SAME specific-enthalpy chart; the other retains
//! the initial optical constants as a control. During the 350 K fusion plateau
//! n stays 1 while k changes from 2 to 1 at 480 nm and from 3 to 1.5 at 680 nm.
//! These changing indices pass through the original Fresnel owner, so optical
//! absorption can evolve while temperature remains constant.
//!
//! This is an authored constitutive example, not measured material data or an
//! inferred phase mixing rule. The h reference is bound to the caloric identity.
//! Native convection increases at switch-off and refreezes both cells. Geometry,
//! incidence and exterior index stay fixed. The model describes opaque surfaces
//! with intercepted input watts; radiation and changing geometry are not modeled.

use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_render::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use fs_render::thermal_surface::heating::feedback::SurfaceConductorLaw;
use fs_render::thermal_surface::heating::{
    StatefulSurfaceIlluminationSource, SurfaceHeatingBudget, SurfaceHeatingTolerance,
    SurfaceIlluminationSource, SurfaceIncidence, SurfaceOpticalHeating, SurfacePowerResponse,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
};
use fs_thermal::observers::CoupledObservers;
use std::error::Error;
use std::sync::Arc;

const MATERIAL: MaterialId = MaterialId(1);
const EVOLVING: CellId = CellId(1);
const FROZEN: CellId = CellId(2);
const WAVELENGTHS: [f64; 2] = [480.0, 680.0];
const SWITCH_OFF_S: f64 = 8.5;

fn boundary(htc_w_m2_k: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k,
        reservoir: None,
        source_w: 0.0,
    }
}

fn optical_source(identity: u8) -> Result<ConductorSource, Box<dyn Error>> {
    Ok(ConductorSource::try_new(
        ContentHash([identity; 32]),
        ConductorDataStatus::Representative,
    )?)
}

fn optics(extinction_scale: f64) -> Result<ConductorOptics, Box<dyn Error>> {
    let samples = [
        380.0, 430.0, 480.0, 530.0, 580.0, 630.0, 680.0, 730.0, 780.0,
    ]
    .map(|wavelength| {
        let k = extinction_scale * (2.0 + (wavelength - 480.0) / 200.0);
        ConductorIorSample::try_new(wavelength, 1.0, k)
    })
    .into_iter()
    .collect::<Result<Vec<_>, _>>()?
    .try_into()
    .map_err(|_| "wrong number of authored optical knots")?;
    Ok(ConductorOptics::try_new(samples, optical_source(221)?)?)
}

fn budget() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 1,
        max_total_power_segments: 3,
        max_source_evaluations: 1000,
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let caloric = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([222; 32]),
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
    for id in [EVOLVING, FROZEN] {
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
    for id in [EVOLVING, FROZEN] {
        events.push(TimedEvent {
            time_s: SWITCH_OFF_S,
            event: Event::Boundary(id, boundary(2000.0)),
        });
    }
    let step = StepConfig {
        max_dt_s: 0.025,
        energy_tolerance_j: 1e-8,
    };
    let mut thermal = DepositionSimulation::new(
        events,
        &[(MATERIAL, &caloric)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 4,
            max_substeps: 2000,
        },
    )?
    .advance_to(0.0, step, || Ok(()))?
    .0;
    let solid = optics(1.0)?;
    let liquid = optics(0.5)?;
    let law = Arc::new(SurfaceConductorLaw::try_enthalpy(
        optical_source(223)?,
        caloric.identity(),
        (300.0, 600.0),
        &[
            (0.0, solid),
            (50_000.0, solid),
            (150_000.0, liquid),
            (400_000.0, liquid),
        ],
        4,
    )?);
    let power = SpectralPowerSchedule::try_new(
        &[
            SpectralPowerSegment {
                start_s: 0.5,
                end_s: 4.5,
                start_power_w: [80.0, 0.0],
                end_power_w: [80.0, 0.0],
            },
            SpectralPowerSegment {
                start_s: 4.5,
                end_s: 6.5,
                start_power_w: [0.0; 2],
                end_power_w: [0.0, 100.0],
            },
            SpectralPowerSegment {
                start_s: 6.5,
                end_s: SWITCH_OFF_S,
                start_power_w: [0.0, 100.0],
                end_power_w: [0.0; 2],
            },
        ],
        3,
    )?;
    let tolerance = SurfaceHeatingTolerance {
        power_roundoff_w: 1e-10,
    };
    let evolving = SurfaceOpticalHeating::try_new_with_stateful(
        &[],
        &[StatefulSurfaceIlluminationSource {
            source_id: 1,
            thermal_cell: EVOLVING,
            material_id: MATERIAL,
            temperature_bounds_k: (300.0, 600.0),
            incidence: SurfaceIncidence {
                wavelengths_nm: WAVELENGTHS,
                incident_eta: [1.0; 2],
                incidence_cosine: 1.0,
                s_power_fraction: [0.5; 2],
            },
            optics: law.clone(),
            power: power.clone(),
        }],
        &thermal,
        budget(),
        4,
        tolerance,
        || Ok(()),
    )?;
    let frozen_response =
        SurfacePowerResponse::from_conductor(WAVELENGTHS, &solid, [1.0; 2], 1.0, [0.5; 2])?;
    let frozen = SurfaceOpticalHeating::try_new(
        &[SurfaceIlluminationSource {
            source_id: 2,
            thermal_cell: FROZEN,
            material_id: MATERIAL,
            temperature_bounds_k: (300.0, 600.0),
            response: frozen_response,
            power,
        }],
        &thermal,
        budget(),
        tolerance,
        || Ok(()),
    )?;
    // Independent receipts for the two alternatives, with one native advance,
    // one event schedule, and the original shared heat-source composition.
    let mut joint = CoupledObservers::new(evolving, frozen);
    let mut convection_j = 0.0;
    let mut peak_temperature = [300.0_f64; 2];
    let mut melting = [false; 2];
    let mut melted = [false; 2];
    let mut freezing = [false; 2];
    let mut fresnel_builds = 0usize;
    let mut plateau_extinction = Vec::new();
    let initial_energy_j = thermal.energy_j();
    eprintln!(
        "Synthetic 1 g specimens: each receives 320 J blue and 200 J red; native fusion plateau at 350 K."
    );
    eprintln!(
        "Evolving surface: authored k(h), tied to the native caloric identity. Control: initial k remains fixed."
    );
    eprintln!(
        "Both use one native clock and convection; source off/quench at 8.5 s. No measured-material claim."
    );
    println!(
        "time_s,evolving_temperature_k,frozen_temperature_k,evolving_liquid_fraction,frozen_liquid_fraction,evolving_blue_k,evolving_red_k,evolving_absorbed_j,frozen_absorbed_j,evolving_reflected_j,frozen_reflected_j,convection_to_material_j,stored_energy_j,energy_residual_j,fresnel_builds"
    );
    for output in 1..=40 {
        let target_s = 0.5 * f64::from(output);
        let prior_absorbed =
            joint.first.energy().thermal_input_j + joint.second.energy().thermal_input_j;
        let (next, accepted, audit) =
            thermal.advance_observed(&joint, target_s, step, || Ok(()))?;
        thermal = next;
        joint = accepted;
        let dynamic_energy = joint.first.energy();
        let frozen_energy = joint.second.energy();
        let absorbed = dynamic_energy.thermal_input_j + frozen_energy.thermal_input_j;
        convection_j += audit.external_heat_j - (absorbed - prior_absorbed);
        let stored_j = thermal.energy_j() - initial_energy_j;
        let residual = stored_j - absorbed - convection_j;
        if residual.abs() > 1e-8 || audit.residual_j.abs() > 1e-8 {
            return Err("phase-dependent surface/thermal energy account failed".into());
        }
        fresnel_builds += joint.first.response_evaluations();
        if target_s > SWITCH_OFF_S && joint.first.response_evaluations() != 0 {
            return Err("dark cooling rebuilt an unused Fresnel packet".into());
        }
        let mut temperatures = [0.0; 2];
        let mut fractions = [0.0; 2];
        for (index, id) in [EVOLVING, FROZEN].into_iter().enumerate() {
            let state = caloric.state_at_specific_enthalpy(thermal.cells()[&id].specific_h_j_kg)?;
            temperatures[index] = state.temperature_k();
            fractions[index] = state.liquid_mass_fraction();
            peak_temperature[index] = peak_temperature[index].max(temperatures[index]);
            melted[index] |= fractions[index] == 1.0;
            if fractions[index] > 0.0 && fractions[index] < 1.0 {
                melting[index] |= target_s < SWITCH_OFF_S;
                freezing[index] |= target_s > SWITCH_OFF_S;
            }
        }
        let current_optics = law.resolve(&thermal, EVOLVING, MATERIAL)?;
        let k_blue = current_optics.complex_ior(WAVELENGTHS[0])?.k();
        let k_red = current_optics.complex_ior(WAVELENGTHS[1])?.k();
        if target_s < SWITCH_OFF_S && fractions[0] > 0.0 && fractions[0] < 1.0 {
            plateau_extinction.push(k_blue);
        }
        println!(
            "{target_s:.2},{:.9},{:.9},{:.9},{:.9},{k_blue:.9},{k_red:.9},{:.9},{:.9},{:.9},{:.9},{convection_j:.9},{stored_j:.9},{residual:.9e},{}",
            temperatures[0],
            temperatures[1],
            fractions[0],
            fractions[1],
            dynamic_energy.thermal_input_j,
            frozen_energy.thermal_input_j,
            dynamic_energy.reflected_j.iter().sum::<f64>(),
            frozen_energy.reflected_j.iter().sum::<f64>(),
            joint.first.response_evaluations()
        );
    }
    for energy in [joint.first.energy(), joint.second.energy()] {
        for (actual, expected) in energy.incident_j.into_iter().zip([320.0, 200.0]) {
            if (actual - expected).abs() > 1e-8 {
                return Err("spectral fluence did not match the authored pulse areas".into());
            }
        }
    }
    let expected_frozen = 320.0 * 0.5 + 200.0 * (4.0 / 13.0);
    if (joint.second.energy().thermal_input_j - expected_frozen).abs() > 1e-8 {
        return Err("frozen control disagrees with analytic normal-incidence absorption".into());
    }
    if joint.first.energy().thermal_input_j <= joint.second.energy().thermal_input_j
        || peak_temperature[0] <= peak_temperature[1]
    {
        return Err(
            "authored state feedback did not increase absorption relative to the frozen control"
                .into(),
        );
    }
    if plateau_extinction.len() < 2
        || plateau_extinction.windows(2).all(|pair| pair[0] == pair[1])
        || fresnel_builds == 0
    {
        return Err("absorption did not evolve through the isothermal latent plateau".into());
    }
    for (index, id) in [EVOLVING, FROZEN].into_iter().enumerate() {
        let state = caloric.state_at_specific_enthalpy(thermal.cells()[&id].specific_h_j_kg)?;
        if !melting[index]
            || !melted[index]
            || !freezing[index]
            || state.liquid_mass_fraction() != 0.0
        {
            return Err("a specimen did not finish melting and dark refreezing".into());
        }
    }
    eprintln!(
        "Absorbed: evolving={:.9} J, frozen={:.9} J. Peak temperatures={:.6}/{:.6} K. Fresnel rebuilds={fresnel_builds}; both refroze.",
        joint.first.energy().thermal_input_j,
        joint.second.energy().thermal_input_j,
        peak_temperature[0],
        peak_temperature[1]
    );
    Ok(())
}
