//! Finite hot gas cooling into finite walls with implicit radiation feedback.
//!
//! cargo run --release -p fs-airflow --features compartment-fire \
//!   --example radiating_compartment -- [gray|spectral|implicit|ventilated] [end_s=10] [max_dt_s=0.1]
//!
//! The original NASA air, TwoZone, ZoneNetwork, Wall1d, material-card emissivity
//! and participating enclosure owners perform the calculation. In `gray` and
//! `spectral` modes gas temperatures are frozen per source candidate; wall storage and radiation are
//! solved together. Accepted gas and wall clocks advance together. No external
//! heat enters the closed modes, so their combined physical energy is fixed.
//!
//! Optical coefficients are EXPLICITLY AUTHORED coupling inputs, not measured
//! nitrogen/oxygen spectroscopy. Spectral mode uses those same constant bands
//! over a finite wavelength window, providing a gray-limit reference. The
//! geometry is a reciprocal normal-chord, parallel-plate approximation, not an
//! automatically generated room enclosure. This is a numerical integration
//! example, not a validated fire or cooling experiment or soot-yield closure.
//! `implicit` closes gas and wall endpoint feedback with gray optics, rebuilding
//! paths and convective contacts from each gas trial without shortening a step.
//! `ventilated` applies that heat operator inside the original vent-pressure
//! solve; ambient mass and physical enthalpy have their own external accounts.
//! All modes retain first-order thermal integration.

#[path = "compartment_fire/air.rs"]
mod air;

use fs_airflow::compartment::radiation::spectral::{SpectralBand, SpectralOptics, SpectralSolve};
use fs_airflow::compartment::radiation::walls::closed_room::{
    ClosedRoomWallSolve, ClosedRoomWallWorkspace, candidate_closed_room_walls,
    candidate_split_room_walls,
};
use fs_airflow::compartment::radiation::walls::{
    EnclosureWallSolve, EnclosureWallWork, WallOptics,
};
use fs_airflow::compartment::radiation::{
    AbsorptionCoefficient, LayerAbsorption, OpticalPath, OpticalSegment, ParticipatingEnclosure,
    RadiationError, RadiationSolve,
};
use fs_airflow::compartment::wall::{
    GasContact, HeatTransferCoefficient, SpecificEnthalpy, ThermalConductivity, Wall1d,
    WallBoundary, WallGeometry,
};
use fs_airflow::compartment::{Ambient, Endpoint, NetworkControl, Room, WallVent, ZoneNetwork};
use fs_airflow::stratified::{DischargeCoefficient, WallOpening};
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
            "radiating_compartment [gray|spectral|implicit|ventilated] [end_s=10] [max_dt_s=0.1]\nOriginal two-layer NASA air and two finite walls; explicit authored opacity and normal-chord geometry.\nSpectral mode resolves the same constant opacity over 1e-7..1e-3 m as a finite-window gray-limit reference.\nImplicit mode solves gas/wall endpoint heat feedback with gray optics and supports requested timesteps up to 300 s.\nVentilated mode places endpoint heat feedback inside the vent-pressure solve with a 0.01 m wide, 2 m tall opening to 300 K air; mass and physical enthalpy exchange are counted explicitly.\nNo external heater, combustion, measured air spectroscopy or experimental-validation claim."
        );
        return Ok(());
    }
    if args.len() > 3 {
        return Err("expected at most mode, end time and maximum timestep".into());
    }
    let mode = args.first().map(String::as_str).unwrap_or("spectral");
    if !matches!(mode, "gray" | "spectral" | "implicit" | "ventilated") {
        return Err("mode must be gray, spectral, implicit or ventilated".into());
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
        || !(0.001..=if mode == "implicit" { 300.0 } else { 0.5 }).contains(&maximum_dt)
        || end / maximum_dt > 50_000.0
    {
        return Err("require 0.001 <= end <= 300 s, max_dt >= 0.001 s (up to 0.5 s, or 300 s for implicit) and at most 50000 requested steps".into());
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
        if mode == "ventilated" {
            vec![WallVent {
                left: Endpoint::Room(0),
                right: Endpoint::Ambient,
                opening: WallOpening::new(
                    Length::new(0.0),
                    Length::new(2.0),
                    Length::new(0.01),
                    DischargeCoefficient::sharp_edged()?,
                )?,
            }]
        } else {
            vec![]
        },
        vec![],
        Acceleration::new(9.81),
        || false,
    )?;
    let mut accepted_room = network.rooms()[0].gas.clone();
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
    let optics = if mode != "spectral" {
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
    let mut gas_work = ClosedRoomWallWorkspace::new(&accepted_room)?;
    let mut previous_gas_power = [Power::new(0.0); 2];
    let mut previous_power = [Power::new(0.0); 2];
    let initial_energy = network.rooms()[0].gas.physical_internal_energy().value()
        + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
    let initial_species: Vec<_> = [Layer::Lower, Layer::Upper]
        .into_iter()
        .map(|l| network.rooms()[0].gas.layer(l).species_masses().to_vec())
        .collect();
    let mut ambient_species = vec![0.0; initial_species[0].len()];
    let mut ambient_energy = 0.0;
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
    while accepted_room.time().value() < end {
        if steps >= 50_000 {
            return Err("accepted-step budget exhausted".into());
        }
        let remaining = end - accepted_room.time().value();
        let dt = if mode != "implicit"
            && remaining > maximum_dt
            && remaining - maximum_dt < network_controls.minimum_step_s
        {
            0.5 * remaining
        } else {
            maximum_dt.min(remaining)
        };
        // The implicit path rebuilds optical geometry at every gas endpoint;
        // the network source path freezes it only within one candidate/retry.
        let (candidate, step_halvings, step_work) = if mode == "implicit" {
            let result = candidate_closed_room_walls(
                &accepted_room,
                Time::new(dt),
                Some(previous_gas_power),
                ClosedRoomWallSolve::default(),
                &mut gas_work,
                |endpoint, cancel| {
                    make_enclosure(endpoint)?.candidate_walls(
                        &walls,
                        Time::new(dt),
                        &[
                            boundary(endpoint, Layer::Lower),
                            boundary(endpoint, Layer::Upper),
                        ],
                        [
                            endpoint.layer(Layer::Lower).temperature(),
                            endpoint.layer(Layer::Upper).temperature(),
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
                        cancel,
                    )
                },
                || false,
            )?;
            previous_gas_power = result.walls.report.gas_heat_power;
            maximum_power_residual = result
                .report
                .endpoint_power_residual_w
                .into_iter()
                .map(f64::abs)
                .fold(maximum_power_residual, f64::max);
            maximum_temperature_residual = maximum_temperature_residual
                .max(result.report.maximum_gas_temperature_residual_k)
                .max(result.report.maximum_wall_temperature_residual_k);
            let counts = result.report.successful_trial_work;
            accepted_room = result.room;
            (result.walls, 0, counts)
        } else if mode == "ventilated" {
            let mut pressure_trial_work = EnclosureWallWork::default();
            let (network_report, (), (candidate, endpoint_report)) = network
                .step_with_pressure_sources(
                    Time::new(dt),
                    network_controls,
                    |_, _, _, _| Ok(()),
                    |rooms, attempted, _, additional, cancel| {
                        let result = candidate_split_room_walls(
                            &rooms[0].gas,
                            attempted,
                            Some(previous_gas_power),
                            ClosedRoomWallSolve::default(),
                            &mut ClosedRoomWallWorkspace::new(&rooms[0].gas)?,
                            |endpoint, cancel| {
                                make_enclosure(endpoint)?.candidate_walls(
                                    &walls,
                                    attempted,
                                    &[
                                        boundary(endpoint, Layer::Lower),
                                        boundary(endpoint, Layer::Upper),
                                    ],
                                    [
                                        endpoint.layer(Layer::Lower).temperature(),
                                        endpoint.layer(Layer::Upper).temperature(),
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
                                    cancel,
                                )
                            },
                            &mut *cancel,
                        )?;
                        for power in result.walls.report.gas_heat_power {
                            additional.add_heat(power)?;
                        }
                        let counts = result.report.successful_trial_work;
                        pressure_trial_work.wall_evaluations += counts.wall_evaluations;
                        pressure_trial_work.radiation_evaluations += counts.radiation_evaluations;
                        pressure_trial_work.wavelength_evaluations += counts.wavelength_evaluations;
                        pressure_trial_work.spectral_factorizations +=
                            counts.spectral_factorizations;
                        rooms[0].gas = result.room;
                        Ok((result.walls, result.report))
                    },
                    || false,
                )?;
            accepted_room = network.rooms()[0].gas.clone();
            previous_gas_power = candidate.report.gas_heat_power;
            maximum_power_residual = endpoint_report
                .endpoint_power_residual_w
                .into_iter()
                .map(f64::abs)
                .fold(maximum_power_residual, f64::max);
            maximum_temperature_residual = maximum_temperature_residual
                .max(endpoint_report.maximum_gas_temperature_residual_k)
                .max(endpoint_report.maximum_wall_temperature_residual_k);
            // The network external account includes internal gas/wall heat.
            // Remove it once to isolate actual ambient enthalpy for the combined balance.
            ambient_energy += network_report.external_physical_energy_j
                - network_report.elapsed.value()
                    * candidate
                        .report
                        .gas_heat_power
                        .iter()
                        .map(|p| p.value())
                        .sum::<f64>();
            for (total, supplied) in ambient_species
                .iter_mut()
                .zip(network_report.external_species_change_kg)
            {
                *total += supplied;
            }
            (candidate, network_report.halvings, pressure_trial_work)
        } else {
            let enclosure = make_enclosure(&accepted_room)?;
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
            accepted_room = network.rooms()[0].gas.clone();
            let counts = candidate.report.work;
            (candidate, network_report.halvings, counts)
        };
        let faces: [_; 2] =
            std::array::from_fn(|i| candidate.walls[i].report.surface_temperatures[0]);
        previous_power.copy_from_slice(candidate.report.radiation.surface_heat_power());
        // Install the returned original owners without another fallible operation.
        for (wall, next) in walls.iter_mut().zip(candidate.walls) {
            *wall = next.wall;
        }
        let report = candidate.report;
        let room = &accepted_room;
        if walls.iter().any(|w| w.time() != room.time()) {
            return Err("accepted gas/wall clocks differ".into());
        }
        let energy = room.physical_internal_energy().value()
            + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
        let defect = energy - initial_energy - ambient_energy;
        maximum_energy_residual = maximum_energy_residual.max(defect.abs());
        maximum_power_residual = maximum_power_residual.max(report.maximum_power_residual_w);
        maximum_temperature_residual =
            maximum_temperature_residual.max(report.maximum_temperature_residual_k);
        if !defect.is_finite() || defect.abs() > 0.05 {
            return Err(
                format!("combined gas/wall physical-energy balance refused: {defect:e} J").into(),
            );
        }
        for (i, layer) in [Layer::Lower, Layer::Upper]
            .into_iter()
            .enumerate()
            .filter(|_| mode != "ventilated")
        {
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
        for (i, supplied) in ambient_species.iter().enumerate() {
            let retained = room.layer(Layer::Lower).species_masses()[i]
                + room.layer(Layer::Upper).species_masses()[i];
            if (retained - initial_species[0][i] - initial_species[1][i] - supplied).abs() > 1e-9 {
                return Err("original gas species and ambient transfer do not balance".into());
            }
        }
        steps += 1;
        halvings += step_halvings;
        total_wall_evaluations += step_work.wall_evaluations;
        total_radiation_evaluations += step_work.radiation_evaluations;
        total_wavelength_evaluations += step_work.wavelength_evaluations;
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
                step_work.wall_evaluations,
                step_work.radiation_evaluations,
                step_work.wavelength_evaluations
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
