//! Finite hot gas cooling into finite walls with implicit radiation feedback.
//!
//! cargo run --release -p fs-airflow --features compartment-fire \
//!   --example radiating_compartment -- [gray|spectral] [end_s=10] [max_dt_s=0.1]
//!
//! The original NASA air, TwoZone, ZoneNetwork, Wall1d, material-card emissivity
//! and participating enclosure owners perform the calculation. Gas temperatures
//! are frozen during each source candidate; wall storage and radiation are
//! solved together. Accepted gas and wall clocks advance together. No external
//! heat enters this closed example, so their combined physical energy is fixed.
//!
//! Optical coefficients are EXPLICITLY AUTHORED coupling inputs, not measured
//! nitrogen/oxygen spectroscopy. Spectral mode uses those same constant bands
//! over a finite wavelength window, providing a gray-limit reference. The
//! geometry is a reciprocal normal-chord, parallel-plate approximation, not an
//! automatically generated room enclosure. This is a numerical integration
//! example, not a validated fire or cooling experiment, soot-yield closure,
//! fully implicit gas integration or second-order time method.

#[path = "compartment_fire/air.rs"]
mod air;

use fs_airflow::compartment::radiation::spectral::{SpectralBand, SpectralOptics, SpectralSolve};
use fs_airflow::compartment::radiation::walls::{EnclosureWallSolve, WallOptics};
use fs_airflow::compartment::radiation::{
    AbsorptionCoefficient, LayerAbsorption, OpticalPath, OpticalSegment, ParticipatingEnclosure,
    RadiationError, RadiationSolve,
};
use fs_airflow::compartment::wall::{
    GasContact, HeatTransferCoefficient, SpecificEnthalpy, ThermalConductivity, Wall1d,
    WallBoundary, WallGeometry,
};
use fs_airflow::compartment::{Ambient, NetworkControl, Room, ZoneNetwork};
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ViewFactorMatrix,
};
use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::compartment::{Layer, TwoZone};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Acceleration, Area, Density, Length, Power, Pressure, Temperature, Time};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use std::error::Error;
use std::io::{self, BufWriter, Write};

struct WallCaloric;
impl ThermalLaw for WallCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / 1000.0;
        if (200.0..=1200.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}

fn surface_card(value: f64) -> Result<MaterialCard, Box<dyn Error>> {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar {
            value,
            dims: EMISSIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1200.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "authored finite-wall coupling fixture, constant diffuse-gray emissivity"
                .into(),
            license: "authored-example-input".into(),
            artifact: None,
        },
    })?;
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored wall".into(),
            phase: "solid".into(),
            process: "gray radiating face".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?)
}

fn boundary(room: &TwoZone, contacted: Layer) -> WallBoundary {
    WallBoundary {
        gas: [Layer::Lower, Layer::Upper].map(|layer| GasContact {
            area: Area::new(if layer == contacted { 2.0 } else { 0.0 }),
            temperature: room.layer(layer).temperature(),
            htc: HeatTransferCoefficient::new(8.0),
        }),
        ambient_temperature: Temperature::new(300.0),
        ambient_htc: HeatTransferCoefficient::new(0.0),
        external_radiant_power: Power::new(0.0),
    }
}

fn make_enclosure(room: &TwoZone) -> Result<ParticipatingEnclosure, RadiationError> {
    ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(2.0)
            .expect("fixed positive authored plate area"),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(2.0),
            segments: vec![
                OpticalSegment {
                    layer: Layer::Lower,
                    length: room.interface_height(),
                },
                OpticalSegment {
                    layer: Layer::Upper,
                    length: Length::new(room.height().value() - room.interface_height().value()),
                },
            ],
        }],
        "authored reciprocal normal chord through the actual accepted layer depths",
        1e-12,
        || false,
    )
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args
        .first()
        .is_some_and(|a| matches!(a.as_str(), "-h" | "--help"))
    {
        println!(
            "radiating_compartment [gray|spectral] [end_s=10] [max_dt_s=0.1]\nOriginal two-layer NASA air and two finite walls; explicit authored opacity and normal-chord geometry.\nSpectral mode resolves the same constant opacity over 1e-7..1e-3 m as a finite-window gray-limit reference.\nGas/source composition is first order; radiation and wall surface temperatures are solved together.\nNo external heat, combustion, measured air spectroscopy or experimental-validation claim."
        );
        return Ok(());
    }
    if args.len() > 3 {
        return Err("expected at most mode, end time and maximum timestep".into());
    }
    let mode = args.first().map(String::as_str).unwrap_or("spectral");
    if !matches!(mode, "gray" | "spectral") {
        return Err("mode must be gray or spectral".into());
    }
    let end = args
        .get(1)
        .map(|v| v.parse::<f64>())
        .transpose()?
        .unwrap_or(10.0);
    let maximum_dt = args
        .get(2)
        .map(|v| v.parse::<f64>())
        .transpose()?
        .unwrap_or(0.1);
    if !end.is_finite()
        || !(0.001..=300.0).contains(&end)
        || !maximum_dt.is_finite()
        || !(0.001..=0.5).contains(&maximum_dt)
        || end / maximum_dt > 50_000.0
    {
        return Err("require 0.001 <= end <= 300 s, 0.001 <= max_dt <= 0.5 s and at most 50000 requested steps".into());
    }

    let (gas, fractions) = air::source()?;
    let room = TwoZone::new(
        gas,
        Area::new(2.0),
        Length::new(2.0),
        Length::new(1.0),
        Pressure::new(101325.0),
        [Temperature::new(550.0), Temperature::new(900.0)],
        [fractions.clone(), fractions.clone()],
        || false,
    )?;
    let mut network = ZoneNetwork::new(
        vec![Room {
            gas: room,
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(300.0),
            mass_fractions: fractions,
        },
        vec![],
        vec![],
        Acceleration::new(9.81),
        || false,
    )?;
    let law = WallCaloric;
    let geometry = WallGeometry {
        area: Area::new(2.0),
        thickness: Length::new(0.02),
        density: Density::new(500.0),
        conductivity: ThermalConductivity::new(1.0),
        cells: 4,
    };
    let mut walls = [
        Wall1d::new(&law, geometry, SpecificEnthalpy::new(0.0), || false)?,
        Wall1d::new(&law, geometry, SpecificEnthalpy::new(50_000.0), || false)?,
    ];
    let cards = [surface_card(0.7)?, surface_card(0.9)?];
    let names = ["wall-0", "wall-1"];
    let enclosure = make_enclosure(&network.rooms()[0].gas)?;
    let absorption = [0.2, 0.7].map(|v| {
        LayerAbsorption::new(
            AbsorptionCoefficient::new(v),
            "authored coupling-example absorption, not measured air spectroscopy",
        )
        .unwrap()
    });
    let spectral = SpectralOptics::new(
        vec![SpectralBand {
            wavelengths: [Length::new(1e-7), Length::new(1e-3)],
            absorption: absorption.clone(),
        }],
        [None, None],
        "authored total surface emissivities held gray over the supplied window",
    )?;
    let optics = if mode == "gray" {
        WallOptics::Gray {
            absorption: &absorption,
            solve: RadiationSolve::default(),
        }
    } else {
        WallOptics::Spectral {
            optics: &spectral,
            solve: SpectralSolve {
                heat_tolerance_w: 1e-6,
                radiosity_tolerance_w_m2: 1e-6,
                ..SpectralSolve::default()
            },
        }
    };
    let mut work = enclosure.wall_workspace(&walls)?;
    let mut previous_power = [Power::new(0.0); 2];
    let initial_energy = network.rooms()[0].gas.physical_internal_energy().value()
        + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
    let initial_species: Vec<_> = [Layer::Lower, Layer::Upper]
        .into_iter()
        .map(|l| network.rooms()[0].gas.layer(l).species_masses().to_vec())
        .collect();
    let network_controls = NetworkControl {
        minimum_step_s: 1e-6,
        ..NetworkControl::default()
    };
    let mut steps = 0;
    let mut halvings = 0;
    let mut next_output = 0.0;
    let mut maximum_energy_residual = 0.0_f64;
    let mut maximum_power_residual = 0.0_f64;
    let mut maximum_temperature_residual = 0.0_f64;
    let mut total_wall_evaluations = 0;
    let mut total_radiation_evaluations = 0;
    let mut total_wavelength_evaluations = 0;
    let mut output = BufWriter::new(io::stdout().lock());
    eprintln!(
        "model={mode}; source=original NASA9 low-temperature N2/O2; opacity_m_inv=[0.2,0.7] is authored, not air spectroscopy; external_heat=0; wall_cp_J_kg_K=1000; wall_k_W_m_K=1; wall_density_kg_m3=500; optical_geometry=normal_chord_approximation"
    );
    writeln!(
        output,
        "mode,time_s,lower_temperature_k,upper_temperature_k,pressure_pa,interface_height_m,wall0_inner_face_k,wall1_inner_face_k,gas_lower_heat_w,gas_upper_heat_w,combined_energy_residual_j,maximum_endpoint_power_residual_w,maximum_endpoint_temperature_residual_k,outer_iterations,wall_evaluations,radiation_evaluations,wavelength_evaluations"
    )?;
    while network.rooms()[0].gas.time().value() < end {
        if steps >= 50_000 {
            return Err("accepted-step budget exhausted".into());
        }
        let remaining = end - network.rooms()[0].gas.time().value();
        let dt =
            if remaining > maximum_dt && remaining - maximum_dt < network_controls.minimum_step_s {
                0.5 * remaining
            } else {
                maximum_dt.min(remaining)
            };
        // Regenerate path lengths from the actual accepted interface; they
        // remain frozen only during this candidate interval and its retries.
        let enclosure = make_enclosure(&network.rooms()[0].gas)?;
        let (network_report, candidate) = network.step_with_sources(
            Time::new(dt),
            network_controls,
            |rooms, attempted, rates, cancel| {
                let candidate = enclosure.candidate_room_walls(
                    &rooms[0].gas,
                    &walls,
                    attempted,
                    &[
                        boundary(&rooms[0].gas, Layer::Lower),
                        boundary(&rooms[0].gas, Layer::Upper),
                    ],
                    optics,
                    |i, t| {
                        SurfaceEmissivity::from_card(
                            names[i],
                            &cards[i],
                            t.value(),
                            SelectionPolicy::SingleClaimOnly,
                        )
                    },
                    Some(&previous_power),
                    EnclosureWallSolve::default(),
                    &mut work,
                    &mut *cancel,
                )?;
                candidate.report.add_gas_heat(&mut rates[0])?;
                Ok(candidate)
            },
            || false,
        )?;
        let faces: [_; 2] =
            std::array::from_fn(|i| candidate.walls[i].report.surface_temperatures[0]);
        previous_power.copy_from_slice(candidate.report.radiation.surface_heat_power());
        // Install the returned original owners without another fallible operation.
        for (wall, next) in walls.iter_mut().zip(candidate.walls) {
            *wall = next.wall;
        }
        let report = candidate.report;
        let room = &network.rooms()[0].gas;
        if walls.iter().any(|w| w.time() != room.time()) {
            return Err("accepted gas/wall clocks differ".into());
        }
        let energy = room.physical_internal_energy().value()
            + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
        let defect = energy - initial_energy;
        maximum_energy_residual = maximum_energy_residual.max(defect.abs());
        maximum_power_residual = maximum_power_residual.max(report.maximum_power_residual_w);
        maximum_temperature_residual =
            maximum_temperature_residual.max(report.maximum_temperature_residual_k);
        if !defect.is_finite() || defect.abs() > 0.05 {
            return Err(
                format!("combined gas/wall physical-energy balance refused: {defect:e} J").into(),
            );
        }
        for (i, layer) in [Layer::Lower, Layer::Upper].into_iter().enumerate() {
            if room
                .layer(layer)
                .species_masses()
                .iter()
                .zip(&initial_species[i])
                .any(|(a, b)| (a - b).abs() > 1e-10)
            {
                return Err("closed example changed species inventories".into());
            }
        }
        steps += 1;
        halvings += network_report.halvings;
        total_wall_evaluations += report.work.wall_evaluations;
        total_radiation_evaluations += report.work.radiation_evaluations;
        total_wavelength_evaluations += report.work.wavelength_evaluations;
        if room.time().value() >= next_output || room.time().value() >= end {
            writeln!(
                output,
                "{mode},{:.9},{:.9},{:.9},{:.9},{:.12},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{},{},{},{}",
                room.time().value(),
                room.layer(Layer::Lower).temperature().value(),
                room.layer(Layer::Upper).temperature().value(),
                room.pressure().value(),
                room.interface_height().value(),
                faces[0].value(),
                faces[1].value(),
                report.gas_heat_power[0].value(),
                report.gas_heat_power[1].value(),
                defect,
                report.maximum_power_residual_w,
                report.maximum_temperature_residual_k,
                report.iterations,
                report.work.wall_evaluations,
                report.work.radiation_evaluations,
                report.work.wavelength_evaluations
            )?;
            output.flush()?;
            next_output = room.time().value().floor() + 1.0;
        }
    }
    eprintln!(
        "completed mode={mode} end_s={end} max_dt_s={maximum_dt} steps={steps} halvings={halvings} max_combined_energy_residual_j={maximum_energy_residual:.12e} max_endpoint_power_residual_w={maximum_power_residual:.12e} max_endpoint_temperature_residual_k={maximum_temperature_residual:.12e} wall_evaluations={total_wall_evaluations} radiation_evaluations={total_radiation_evaluations} wavelength_evaluations={total_wavelength_evaluations}"
    );
    Ok(())
}
