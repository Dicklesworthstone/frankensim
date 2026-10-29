//! Heated solid -> original NASA volatile -> oxygen-limited compartment burning.
//!
//! cargo run -p fs-airflow --features compartment-fire --example pyrolyzing_compartment_fire --
//!   [all|closed|restricted|open|fan-exhaust] [end_s=30] [max_dt_s=0.02] [cells=24] [irradiance_W_m2=50000]
//!
//! This is a manufactured coupling case, NOT a PMMA/wood material or experimental
//! fire prediction. The authored condensed fuel releases pure methane, with an
//! explicitly compatible formation-energy reference. Original NASA methane and
//! product calorics are reused unchanged. Heat drives the original spatial solid
//! reaction; no mass-release history or combustion heat is prescribed.
//!
//! A 20-litre, well-mixed compartment exercises the one-active-layer TwoZone limit.
//! The front solid face exchanges convection with that actual gas, absorbs an
//! external heater, and radiates to an external 300 K blackbody environment.
//! Convection is an internal equal/opposite transfer. All other walls are
//! adiabatic. This deliberately simple optical boundary is not a smoke-radiation
//! or flame-feedback closure. The solid is fixed geometry, non-charring, and
//! non-oxidative. Burning is instantaneous complete oxidation limited by the
//! actually transported oxygen/fuel; it is not an ignition or finite-rate model.
//!
//! Solid/boundary feedback is first-order explicit at the accepted gas state.
//! Chemistry is evaluated inside every vent-pressure trial. Retried attempts
//! recompute solid release, and only the accepted solid/gas pair is installed.
//! A fourth case adds an explicitly authored pressure-dependent exhaust fan to
//! the restricted doorway. Fan work is delivered downstream to ambient; it is
//! reported but is not another heat source inside the finite room.
//!
//! This fixed-geometry example refuses when any cell retains less than 1e-4 of
//! its initial mass. Near surface depletion the unchanged cell geometry and
//! conductance require a recession/remeshing model; no seed mass, heat-capacity
//! floor or release clipping is introduced to extend the simulation.

#[path = "ventilated_compartment_fire/nasa9.rs"]
mod nasa9;

use fs_airflow::compartment::mechanical::{
    FanControl, FanDomain, FanEnergy, FanModel, FanPort, FanVent,
};
use fs_airflow::compartment::pyrolysis::{
    VolatileGasError, VolatileMixture, VolatileTransferControl,
};
use fs_airflow::compartment::{
    Ambient, Endpoint, NetworkControl, NetworkError, Room, WallVent, ZoneNetwork,
};
use fs_airflow::stratified::{DischargeCoefficient, WallOpening};
use fs_airflow::{
    FanArrangement, FanBank, FanCurve, FanPoint, LossElement, LossNetwork, LossResistance,
    SourceProvenance, ToleranceBasis,
};
use fs_gas::multicomponent::MixtureGas;
use fs_gas::multicomponent::compartment::burning::InventoryLimitedReaction;
use fs_gas::multicomponent::compartment::{Layer, TwoZone};
use fs_qty::{
    Acceleration, Amount, Area, Density, Dimensionless, Length, Pressure, Temperature, Time,
    VolumetricFlowRate,
};
use fs_thermal::enthalpy::pyrolysis::{
    Arrhenius, CondensedCaloric, PyrolysisControls, PyrolysisMaterial, PyrolysisSimulation,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
    series_conductance_w_k,
};
use fs_thermochem::ReactionId;
use std::error::Error;
use std::io::{self, BufWriter, Write};

const AMBIENT: f64 = 300.0;
const SOLID_AREA: f64 = 0.01;
const SOLID_THICKNESS: f64 = 0.006;
const SOLID_DENSITY: f64 = 850.0;
const SOLID_CP: f64 = 1500.0;
const SOLID_K: f64 = 0.2;
const SOLID_DOMAIN: [f64; 2] = [250.0, 1000.0];
const REACTION_ENTHALPY_AT_300: f64 = 1_000_000.0;
const EMISSIVITY: f64 = 0.9;
const HTC: f64 = 10.0;
const SIGMA: f64 = 5.670_374_419e-8;
const MINIMUM_STEP_S: f64 = 1e-5;
const MINIMUM_RETAINED_CELL_FRACTION: f64 = 1e-4;

struct SolidCaloric {
    h300: f64,
}
impl ThermalLaw for SolidCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = AMBIENT + (h - self.h300) / SOLID_CP;
        if t.is_finite() && (SOLID_DOMAIN[0]..=SOLID_DOMAIN[1]).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / SOLID_CP
    }
}
impl CondensedCaloric for SolidCaloric {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, DepositionError> {
        if t.is_finite() && (SOLID_DOMAIN[0]..=SOLID_DOMAIN[1]).contains(&t) {
            Ok(self.h300 + SOLID_CP * (t - AMBIENT))
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}

fn insulated() -> Boundary {
    Boundary {
        ambient_k: AMBIENT,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}

fn solid<'a>(
    law: &'a SolidCaloric,
    volatile: &'a VolatileMixture,
    cells: usize,
) -> Result<PyrolysisSimulation<'a>, DepositionError> {
    let dx = SOLID_THICKNESS / cells as f64;
    let conductance = series_conductance_w_k(SOLID_AREA, [0.5 * dx; 2], [SOLID_K; 2])?;
    let mut events = Vec::with_capacity(2 * cells - 1);
    for i in 0..cells {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: SOLID_DENSITY * SOLID_AREA * dx,
                specific_h_j_kg: law.specific_enthalpy_j_kg(AMBIENT)?,
                geometry: CellGeometry {
                    start_m: [i as f64 * dx, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * dx, 0.0, 0.0],
                    volume_m3: SOLID_AREA * dx,
                    surface_area_m2: 2.0 * SOLID_AREA,
                },
                boundary: insulated(),
            }),
        });
        if i > 0 {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: i as u64,
                    a: CellId(i as u64 - 1),
                    b: CellId(i as u64),
                    area_m2: SOLID_AREA,
                    conductance_w_k: conductance,
                }),
            });
        }
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 96,
            max_links: 95,
            max_events: 191,
            max_substeps: 20_000,
        },
    )?
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.02,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )?
    .0;
    PyrolysisSimulation::new(
        thermal,
        &[PyrolysisMaterial {
            material_id: MaterialId(1),
            virgin: law,
            char_caloric: None,
            volatile,
            char_mass_yield: 0.0,
            kinetics: Arrhenius {
                prefactor_s_inv: 1e8,
                activation_energy_j_mol: 100_000.0,
            },
            temperature_range_k: SOLID_DOMAIN,
            source: "authored non-charring methane-releasing solid; formation reference aligned to original NASA methane at 300 K with 1 MJ/kg decomposition enthalpy; not PMMA or wood",
        }],
        || Ok(()),
    )
}

#[derive(Clone, Copy)]
struct Face {
    temperature: f64,
    from_gas_w: f64,
    from_external_w: f64,
}
fn face(center: f64, gas_t: f64, cells: usize, irradiance: f64) -> Result<Face, DepositionError> {
    let dx = SOLID_THICKNESS / cells as f64;
    let k = series_conductance_w_k(SOLID_AREA, [0.25 * dx; 2], [SOLID_K; 2])?;
    let powers = |t: f64| {
        let from_gas = HTC * SOLID_AREA * (gas_t - t);
        let from_external = EMISSIVITY
            * SOLID_AREA
            * (irradiance - SIGMA * (t - AMBIENT) * (t + AMBIENT) * (t * t + AMBIENT * AMBIENT));
        (from_gas, from_external)
    };
    let residual = |t: f64| {
        let (a, b) = powers(t);
        k * (t - center) - a - b
    };
    let mut lo = center.min(gas_t).min(AMBIENT);
    let mut hi = center.max(gas_t).max(AMBIENT) + EMISSIVITY * SOLID_AREA * irradiance / k;
    if !hi.is_finite() || residual(lo) > 0.0 || residual(hi) < 0.0 {
        return Err(DepositionError::MaterialDomain);
    }
    for _ in 0..96 {
        let t = lo + 0.5 * (hi - lo);
        let r = residual(t);
        if r.abs() <= 1e-9 || t == lo || t == hi {
            if !(SOLID_DOMAIN[0]..=SOLID_DOMAIN[1]).contains(&t) {
                return Err(DepositionError::MaterialDomain);
            }
            let (from_gas_w, from_external_w) = powers(t);
            return Ok(Face {
                temperature: t,
                from_gas_w,
                from_external_w,
            });
        }
        if r < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    Err(DepositionError::MaterialDomain)
}

fn network(gas: &MixtureGas, air: &[f64], width: f64) -> Result<ZoneNetwork, Box<dyn Error>> {
    let room = TwoZone::new(
        gas.clone(),
        Area::new(0.1),
        Length::new(0.2),
        Length::new(0.0),
        Pressure::new(101325.0),
        [Temperature::new(AMBIENT); 2],
        [air.to_vec(), air.to_vec()],
        || false,
    )?;
    let discharge = DischargeCoefficient::new(
        Dimensionless::new(0.7),
        SourceProvenance::new(
            "NIST TN1889v1 section 4.1 nominal wall-opening coefficient; authored small box",
            "NIST-TN1889v1:4.1",
        ),
        0.0,
        ToleranceBasis::EngineeringAllowance,
    )?;
    Ok(ZoneNetwork::new(
        vec![Room {
            gas: room,
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(AMBIENT),
            mass_fractions: air.to_vec(),
        },
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Ambient,
            opening: WallOpening::new(
                Length::new(0.0),
                Length::new(0.2),
                Length::new(width),
                discharge,
            )?,
        }],
        vec![],
        Acceleration::new(9.81),
        || false,
    )?)
}

fn inventory(room: &TwoZone) -> [f64; 5] {
    std::array::from_fn(|s| {
        room.layer(Layer::Lower).species_masses()[s] + room.layer(Layer::Upper).species_masses()[s]
    })
}

fn install_exhaust_fan(network: &mut ZoneNetwork) -> Result<(), Box<dyn Error>> {
    let source = || {
        SourceProvenance::new(
            "authored small-box fan/duct coupling fixture; not a measured product rating",
            "pyrolyzing-compartment-authored-fan-v1",
        )
    };
    let curve = FanCurve::new(
        "authored-small-box-exhaust",
        vec![
            FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(100.0)),
            FanPoint::new(VolumetricFlowRate::new(0.01), Pressure::new(0.0)),
        ],
        source(),
        0.0,
        ToleranceBasis::Analytic,
        VolumetricFlowRate::new(0.0),
        (0.5, 2.0),
    )?;
    let losses = LossNetwork::Element(LossElement::new(
        "authored-complete-exhaust-path",
        LossResistance::new(1e8),
        0.0,
        source(),
        ToleranceBasis::Analytic,
    )?);
    let density = network.rooms()[0]
        .gas
        .layer(Layer::Upper)
        .density_kg_m3()
        .ok_or("fan inlet layer is absent")?;
    let model = FanModel::new(
        FanBank::new(curve, 1, FanArrangement::Series, 1.0)?,
        losses,
        Density::new(density),
        Area::new(0.002),
        FanEnergy {
            total_efficiency: 0.5,
            rejected_loss_fraction: 0.0,
            source: source(),
        },
        FanDomain {
            density_ratio: [0.05, 2.0],
            maximum_pressure_fraction: 0.02,
        },
        || false,
    )?;
    network.set_fan_vents(
        vec![FanVent {
            from: FanPort {
                endpoint: Endpoint::Room(0),
                height: Length::new(0.2),
            },
            to: FanPort {
                endpoint: Endpoint::Ambient,
                height: Length::new(0.2),
            },
            model,
            controls: FanControl::default(),
            closed: false,
        }],
        || false,
    )?;
    Ok(())
}

fn minimum_cell_fraction(solid: &PyrolysisSimulation<'_>, initial_cell_mass: f64) -> f64 {
    solid
        .cells()
        .values()
        .map(|cell| cell.mass_kg / initial_cell_mass)
        .fold(f64::INFINITY, f64::min)
}

fn run(
    name: &str,
    width: f64,
    end: f64,
    maximum_dt: f64,
    cells: usize,
    irradiance: f64,
    output: &mut impl Write,
) -> Result<(), Box<dyn Error>> {
    let (gas, mechanism, air) = nasa9::source()?;
    let volatile = VolatileMixture::new(&gas, &[1.0, 0.0, 0.0, 0.0, 0.0], || false)?;
    let law = SolidCaloric {
        h300: volatile
            .specific_enthalpy(Temperature::new(AMBIENT), || false)?
            .value()
            - REACTION_ENTHALPY_AT_300,
    };
    let mut solid = solid(&law, &volatile, cells)?;
    let mut network = network(&gas, &air, width)?;
    if name == "fan-exhaust" {
        install_exhaust_fan(&mut network)?;
    }
    let reaction =
        InventoryLimitedReaction::new(&gas, mechanism, &ReactionId::new(nasa9::REACTION)?, || {
            false
        })?;
    let controls = NetworkControl {
        minimum_step_s: MINIMUM_STEP_S,
        ..NetworkControl::default()
    };
    let initial_species = inventory(&network.rooms()[0].gas);
    let initial_solid_mass = solid.mass_kg();
    let initial_cell_mass = initial_solid_mass / cells as f64;
    let initial_energy =
        solid.energy_j() + network.rooms()[0].gas.physical_internal_energy().value();
    let mut emitted = 0.0;
    let mut extent = 0.0;
    let mut external_energy = 0.0;
    let mut vent_species = [0.0; 5];
    let mut max_mass_defect = 0.0_f64;
    let mut max_pressure_defect = 0.0_f64;
    let mut max_energy_defect = 0.0_f64;
    let mut max_element_defect = 0.0_f64;
    let mut next_output = 0.0;
    let mut first_oxygen_limited = None;
    let mut steps = 0usize;
    let mut halvings = 0usize;
    let mut fan_exhaust_mass = 0.0;
    let mut fan_shaft_energy = 0.0;
    while solid.time_s() < end {
        if steps >= 200_000 {
            return Err("example accepted-step budget exhausted".into());
        }
        if minimum_cell_fraction(&solid, initial_cell_mass) < MINIMUM_RETAINED_CELL_FRACTION {
            return Err("fixed-geometry example scope exhausted: retained cell mass below 1e-4 of initial; surface recession/remeshing required".into());
        }
        let remaining = end - solid.time_s();
        let dt = if remaining > maximum_dt && remaining - maximum_dt < controls.minimum_step_s {
            0.5 * remaining
        } else {
            maximum_dt.min(remaining)
        };
        let (report, (candidate, audit, face), chemistry) = network.step_with_pressure_sources(
            Time::new(dt),
            controls,
            |rooms, attempted, rates, cancel| {
                if solid.time_s() != rooms[0].gas.time().value() {
                    return Err(NetworkError::Invalid("solid/gas clocks differ"));
                }
                let face = face(
                    solid
                        .temperature_k(CellId(0))
                        .map_err(VolatileGasError::from)?,
                    rooms[0].gas.layer(Layer::Upper).temperature().value(),
                    cells,
                    irradiance,
                )
                .map_err(VolatileGasError::from)?;
                let mut prepared = solid.clone();
                let mut boundary = insulated();
                boundary.source_w = face.from_gas_w + face.from_external_w;
                prepared
                    .set_boundary(CellId(0), boundary, || {
                        if cancel() {
                            Err(DepositionError::Cancelled)
                        } else {
                            Ok(())
                        }
                    })
                    .map_err(VolatileGasError::from)?;
                let (candidate, audit) = prepared
                    .advance_to(
                        solid.time_s() + attempted.value(),
                        PyrolysisControls {
                            maximum_dt_s: attempted.value(),
                            ..PyrolysisControls::default()
                        },
                        || {
                            if cancel() {
                                Err(DepositionError::Cancelled)
                            } else {
                                Ok(())
                            }
                        },
                    )
                    .map_err(VolatileGasError::from)?;
                if minimum_cell_fraction(&candidate, initial_cell_mass) < MINIMUM_RETAINED_CELL_FRACTION {
                    return Err(NetworkError::Invalid(
                        "fixed-geometry example scope: retained cell mass below 1e-4 of initial; surface recession/remeshing required",
                    ));
                }
                for release in audit.releases.values() {
                    volatile
                        .release_rates(
                            &rooms[0].gas,
                            *release,
                            attempted,
                            VolatileTransferControl::default(),
                            &mut *cancel,
                        )?
                        .add_to(Layer::Upper, &mut rates[0])?;
                }
                rates[0][Layer::Upper.index()].referenced_energy_w -= face.from_gas_w;
                Ok((candidate, audit, face))
            },
            |rooms, attempted, _, additional, cancel| {
                let room = &mut rooms[0].gas;
                let requested = room.layer(Layer::Upper).species_masses()[0] / nasa9::MASSES[0];
                let burned = reaction.apply_split(
                    room,
                    Layer::Upper,
                    Amount::new(requested),
                    attempted,
                    Temperature::new(AMBIENT),
                    controls.zone,
                    cancel,
                )?;
                additional.add(&burned.burning.chemical_rates)?;
                Ok(burned.burning)
            },
            || false,
        )?;
        // Accepted payload is installed without another fallible/cancellation point.
        solid = candidate;
        steps += 1;
        halvings += report.halvings;
        let actual_dt = report.elapsed.value();
        let fan_point = network
            .last_fan_transfers()
            .first()
            .copied()
            .flatten()
            .and_then(|transfer| transfer.operating_point);
        let [fan_q, fan_mass, fan_shaft, fan_gas] = fan_point.map_or([0.0; 4], |point| {
            [
                point.volume_flow.value(),
                point.mass_rate.value(),
                point.shaft_power.value(),
                point.gas_power.value(),
            ]
        });
        fan_exhaust_mass += actual_dt * fan_mass;
        fan_shaft_energy += actual_dt * fan_shaft;
        emitted += audit.volatile_mass_kg;
        extent += chemistry.actual_extent.value();
        let chemical_physical = actual_dt
            * (chemistry.chemical_rates.referenced_energy_w
                + chemistry.chemical_rates.species_kg_s[..5]
                    .iter()
                    .zip(gas.energy_references())
                    .map(|(m, e)| m * e)
                    .sum::<f64>());
        // Remove internal solid/gas transfer and numerical chemical-reference
        // arithmetic from the network's aggregate authored-source account.
        let vent_energy = report.external_physical_energy_j - audit.volatile_enthalpy_j
            + actual_dt * face.from_gas_w
            - chemical_physical;
        external_energy += vent_energy + actual_dt * face.from_external_w;
        for (s, vent) in vent_species.iter_mut().enumerate() {
            *vent += report.external_species_change_kg[s]
                - if s == 0 { audit.volatile_mass_kg } else { 0.0 }
                - actual_dt * chemistry.chemical_rates.species_kg_s[s];
        }
        let room = &network.rooms()[0].gas;
        if solid.time_s() != room.time().value() {
            return Err("accepted clock mismatch".into());
        }
        if chemistry.limiting_species == Some(4) && first_oxygen_limited.is_none() {
            first_oxygen_limited = Some(solid.time_s());
        }
        let masses = inventory(room);
        let species_defect = (0..5)
            .map(|s| {
                (masses[s]
                    - initial_species[s]
                    - vent_species[s]
                    - if s == 0 { emitted } else { 0.0 }
                    - nasa9::STOICHIOMETRY[s] * nasa9::MASSES[s] * extent)
                    .abs()
            })
            .fold(0.0, f64::max);
        let mass_defect = solid.mass_kg() + masses.iter().sum::<f64>()
            - initial_solid_mass
            - initial_species.iter().sum::<f64>()
            - vent_species.iter().sum::<f64>();
        let energy_defect = solid.energy_j() + room.physical_internal_energy().value()
            - initial_energy
            - external_energy;
        let elements = gas.source().elemental_matrix();
        for e in 0..4 {
            let defect = (0..5)
                .map(|s| {
                    (masses[s] - initial_species[s] - vent_species[s]) / nasa9::MASSES[s]
                        * elements.get(e, s).unwrap() as f64
                })
                .sum::<f64>()
                + (solid.mass_kg() - initial_solid_mass) / nasa9::MASSES[0]
                    * elements.get(e, 0).unwrap() as f64;
            max_element_defect = max_element_defect.max(defect.abs());
        }
        max_mass_defect = max_mass_defect.max(mass_defect.abs());
        max_energy_defect = max_energy_defect.max(energy_defect.abs());
        max_pressure_defect = max_pressure_defect.max(report.maximum_pressure_residual_pa);
        if masses.iter().any(|m| !m.is_finite() || *m < 0.0)
            || room.layer(Layer::Lower).mass_kg() != 0.0
            || species_defect > 1e-9
            || mass_defect.abs() > 1e-9
            || energy_defect.abs() > 0.02
            || max_element_defect > 1e-7
        {
            return Err(format!("combined balance/domain refused: mass={mass_defect:e}kg energy={energy_defect:e}J elements={max_element_defect:e}mol").into());
        }
        if solid.time_s() >= next_output || solid.time_s() >= end {
            writeln!(
                output,
                "{name},{:.9},{:.9},{:.9},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                solid.time_s(),
                room.layer(Layer::Upper).temperature().value(),
                room.pressure().value(),
                face.temperature,
                solid.temperature_k(CellId(0))?,
                solid.mass_kg(),
                emitted,
                extent * nasa9::MASSES[0],
                masses[0],
                masses[4],
                chemistry.released_power.value(),
                mass_defect,
                energy_defect,
                minimum_cell_fraction(&solid, initial_cell_mass),
                fan_q,
                fan_mass,
                fan_shaft,
                fan_gas
            )?;
            next_output = solid.time_s().floor() + 1.0;
            output.flush()?;
        }
    }
    if width == 0.0 {
        let expected =
            emitted.min(initial_species[4] * nasa9::MASSES[0] / (2.0 * nasa9::MASSES[4]));
        if (extent * nasa9::MASSES[0] - expected).abs() > 1e-10 {
            return Err("closed compartment failed exact oxygen/fuel inventory limit".into());
        }
    }
    eprintln!(
        "case={name} width_m={width} steps={steps} halvings={halvings} emitted_kg={emitted:.12e} burned_kg={:.12e} first_oxygen_limited_s={first_oxygen_limited:?} minimum_retained_cell_fraction={:.12e} fan_exhaust_mass_kg={fan_exhaust_mass:.12e} fan_shaft_energy_to_ambient_j={fan_shaft_energy:.12e} max_combined_mass_residual_kg={max_mass_defect:.12e} max_combined_energy_residual_j={max_energy_defect:.12e} max_element_residual_mol={max_element_defect:.12e} max_pressure_residual_pa={max_pressure_defect:.12e}",
        extent * nasa9::MASSES[0],
        minimum_cell_fraction(&solid, initial_cell_mass)
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 5 || args.first().is_some_and(|a| a == "--help") {
        eprintln!(
            "usage: pyrolyzing_compartment_fire [all|closed|restricted|open|fan-exhaust] [end_s=30] [max_dt_s=0.02] [cells=24] [irradiance_W_m2=50000]\nAuthored methane-releasing solid coupling case, not PMMA/wood calibration or experimental validation.\nFixed geometry refuses below 1e-4 initial cell mass; minimum network step 1e-5 s.\nAuthored fan: static curve 100 Pa at 0 m3/s to 0 Pa at 0.01 m3/s; total duct resistance 1e8 Pa/(m3/s)^2; outlet area 0.002 m2; total efficiency 0.5; reference density from initial room. All shaft heat goes downstream into ambient."
        );
        return if args.len() <= 1 {
            Ok(())
        } else {
            Err("invalid argument count".into())
        };
    }
    let selected = args.first().map(String::as_str).unwrap_or("all");
    if !["all", "closed", "restricted", "open", "fan-exhaust"].contains(&selected) {
        return Err("unknown case".into());
    }
    let number = |i: usize, default: f64| -> Result<f64, Box<dyn Error>> {
        let x = args
            .get(i)
            .map(|s| s.parse::<f64>())
            .transpose()?
            .unwrap_or(default);
        if x.is_finite() && x > 0.0 {
            Ok(x)
        } else {
            Err("arguments must be finite and positive".into())
        }
    };
    let end = number(1, 30.0)?;
    let dt = number(2, 0.02)?;
    let cells = args
        .get(3)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(24);
    let irradiance = number(4, 50_000.0)?;
    if end > 120.0
        || !(MINIMUM_STEP_S..=0.1).contains(&dt)
        || end / dt > 200_000.0
        || !(2..=96).contains(&cells)
        || irradiance > 100_000.0
    {
        return Err("bounds: end <= 120 s, 1e-5 <= dt <= 0.1 s, <= 200000 requested steps, cells 2..96, irradiance <= 100000 W/m2; depletion/domain may require escalation".into());
    }
    eprintln!(
        "authored methane solid: area_m2={SOLID_AREA} thickness_m={SOLID_THICKNESS} density_kg_m3={SOLID_DENSITY} cp_j_kg_k={SOLID_CP} conductivity_w_m_k={SOLID_K} end_s={end} max_dt_s={dt} cells={cells} irradiance_w_m2={irradiance} gas_source={}",
        nasa9::SOURCE_COMMIT
    );
    let mut out = BufWriter::new(io::stdout().lock());
    writeln!(
        out,
        "case,time_s,gas_temperature_k,pressure_pa,lagged_front_face_k,front_cell_k,solid_mass_kg,emitted_fuel_kg,burned_fuel_kg,gas_fuel_kg,gas_oxygen_kg,reference300K_hrr_w,combined_mass_residual_kg,combined_energy_residual_j,minimum_retained_cell_fraction,fan_volume_rate_m3_s,fan_mass_rate_kg_s,fan_shaft_power_w,fan_downstream_gas_power_w"
    )?;
    for (name, width) in [
        ("closed", 0.0),
        ("restricted", 0.001),
        ("open", 0.02),
        ("fan-exhaust", 0.001),
    ] {
        if selected == "all" || selected == name {
            run(name, width, end, dt, cells, irradiance, &mut out)?;
        }
    }
    out.flush()?;
    Ok(())
}
