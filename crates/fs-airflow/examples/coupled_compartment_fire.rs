//! Sourced finite-rate fuel chemistry with a receding solid and finite walls.
//!
//! cargo run --release -p fs-airflow --features compartment-fire \
//!   --example coupled_compartment_fire -- [all|closed|vented] [end_s=0.1] \
//!   [max_dt_s=0.001] [solid_cells=12] [kinetic|frozen] [be|sdirk2] \
//!   [hydrogen|methane|charring-methane]
//!
//! This numerical composition example uses the original complete Cantera h2o2
//! or GRI-Mech 3.0 mechanism, physical NASA formation energies, implicit pyrolysis,
//! finite-wall conduction, molecular Planck-mean radiation and pressure-coupled
//! ventilation. The manufactured solid releases pure H2 or CH4. The charring
//! case additionally retains condensed carbon with separate density, conductivity
//! and enthalpy. These are NOT material models for PMMA, wood or a metal hydride,
//! nor experimental fire validations.
//! The hot initial gas avoids supplying an invented ignition or radical seed.
//!
//! Gas is homogeneous in one active TwoZone layer. A reciprocal normal chord
//! represents two radiating patches, not a complete geometric room enclosure.
//! Molecular emission means are explicitly adopted as gray absorption; the
//! source's 300..2500 K domain is enforced. Gas/opacity are frozen per interval,
//! while radiation, solid reaction/conduction and surface temperatures are
//! coupled implicitly. Backward Euler remains the chemistry default; optional
//! SDIRK2 changes the reactor method, while overall gas/source composition
//! remains first order.
//! Recession does not change the fixed vessel volume or supply pressure work.

use fs_airflow::compartment::pyrolysis::{VolatileMixture, VolatileTransferControl};
use fs_airflow::compartment::radiation::molecular::MolecularPlanck;
use fs_airflow::compartment::radiation::reactive::{
    SlabBoundary, SlabSurfaceSolve, ThermalSurface, ThermalSurfaceCandidate,
};
use fs_airflow::compartment::radiation::walls::{EnclosureWallSolve, WallOptics};
use fs_airflow::compartment::radiation::{
    AbsorptionCoefficient, LayerAbsorption, OpticalPath, OpticalSegment, ParticipatingEnclosure,
    RadiationSolve,
};
use fs_airflow::compartment::wall::{
    GasContact, HeatTransferCoefficient, SpecificEnthalpy, ThermalConductivity, Wall1d,
    WallBoundary, WallGeometry,
};
use fs_airflow::compartment::{Ambient, Endpoint, NetworkControl, Room, WallVent, ZoneNetwork};
use fs_airflow::stratified::{DischargeCoefficient, WallOpening};
use fs_airflow::{SourceProvenance, ToleranceBasis};
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ViewFactorMatrix,
};
use fs_evidence::ValidityDomain;
use fs_gas::chemistry::{IsochoricKineticsMethod, IsochoricKineticsSolve};
use fs_gas::multicomponent::compartment::kinetics::ZoneKinetics;
use fs_gas::multicomponent::compartment::{Layer, TwoZone};
use fs_gas::multicomponent::{CaloricDomain, MixtureGas};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Acceleration, Area, Density, Dimensionless, Length, Power, Pressure, Temperature, Time,
};
use fs_thermal::enthalpy::pyrolysis::slab::{
    CharProperties, ConductivityMixing, RecedingSlab, SlabGeometry,
};
use fs_thermal::enthalpy::pyrolysis::{
    Arrhenius, CondensedCaloric, PyrolysisControls, PyrolysisMaterial, PyrolysisSimulation,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
    series_conductance_w_k,
};
use std::error::Error;
use std::io::{self, BufWriter, Write};

#[path = "coupled_compartment_fire/fuel.rs"]
mod fuel;
use fuel::{Fuel, SpeciesAxis};

const AREA: f64 = 0.01;
const THICKNESS: f64 = 0.006;
const DENSITY: f64 = 850.0;
const CONDUCTIVITY: f64 = 0.2;
const EXTERNAL_POWER: f64 = 0.9 * AREA * 50_000.0;

struct Caloric {
    h300: f64,
    cp: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.h300) / self.cp;
        if (300.0..=1200.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}
impl CondensedCaloric for Caloric {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, DepositionError> {
        if (300.0..=1200.0).contains(&t) {
            Ok(self.h300 + self.cp * (t - 300.0))
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}

fn solid<'m>(
    law: &'m Caloric,
    char_law: &'m Caloric,
    volatile: &'m VolatileMixture,
    cells: usize,
    fuel: Fuel,
) -> Result<RecedingSlab<'m>, Box<dyn Error>> {
    let dx = THICKNESS / cells as f64;
    let conductance = series_conductance_w_k(AREA, [0.5 * dx; 2], [CONDUCTIVITY; 2])?;
    let mut events = Vec::with_capacity(2 * cells - 1);
    for i in 0..cells {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: DENSITY * AREA * dx,
                specific_h_j_kg: law.specific_enthalpy_j_kg(500.0)?,
                geometry: CellGeometry {
                    start_m: [i as f64 * dx, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * dx, 0.0, 0.0],
                    volume_m3: AREA * dx,
                    surface_area_m2: 2.0 * AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
        if i > 0 {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: i as u64,
                    a: CellId(i as u64 - 1),
                    b: CellId(i as u64),
                    area_m2: AREA,
                    conductance_w_k: conductance,
                }),
            });
        }
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: cells,
            max_links: cells,
            max_events: 2 * cells,
            max_substeps: 10_000,
        },
    )?
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.001,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )?
    .0;
    let original = PyrolysisSimulation::new(
        thermal,
        &[PyrolysisMaterial {
            material_id: MaterialId(1),
            virgin: law,
            char_caloric: (fuel.char_yield() > 0.0).then_some(char_law as &dyn CondensedCaloric),
            volatile,
            char_mass_yield: fuel.char_yield(),
            kinetics: Arrhenius {
                prefactor_s_inv: 1e8,
                activation_energy_j_mol: 100_000.0,
            },
            temperature_range_k: [300.0, 1200.0],
            source: "manufactured H2/CH4-releasing solid; optional 20% carbon residue; original volatile and char enthalpy minus 1 MJ/kg virgin at 300 K; not PMMA, wood or hydride data",
        }],
        || Ok(()),
    )?;
    let geometry = SlabGeometry {
        area_m2: AREA,
        density_kg_m3: DENSITY,
        conductivity_w_m_k: CONDUCTIVITY,
        source: "authored incompressible virgin solid; fixed back, insulated sides; no gas-volume displacement",
    };
    Ok(if fuel.char_yield() > 0.0 {
        RecedingSlab::with_char(
            original,
            geometry,
            CharProperties {
                density_kg_m3: 300.0,
                conductivity_w_m_k: 0.08,
                conductivity_mixing: ConductivityMixing::Parallel,
                source: "authored condensed carbon sensitivity; additive phase volumes and volume-fraction parallel conduction; no pore gas or oxidation",
            },
            || Ok(()),
        )?
    } else {
        RecedingSlab::new(original, geometry, || Ok(()))?
    })
}

fn surface_card(emissivity: f64) -> Result<MaterialCard, Box<dyn Error>> {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar {
            value: emissivity,
            dims: EMISSIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 2500.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "authored constant diffuse-gray example surface".into(),
            license: "authored-example-input".into(),
            artifact: None,
        },
    })?;
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured thermal surface".into(),
            phase: "solid".into(),
            process: "gray radiating face".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?)
}

fn boundary(room: &TwoZone, htc: f64) -> WallBoundary {
    WallBoundary {
        gas: [Layer::Lower, Layer::Upper].map(|layer| GasContact {
            area: Area::new(if layer == Layer::Upper { AREA } else { 0.0 }),
            temperature: room.layer(layer).temperature(),
            htc: HeatTransferCoefficient::new(htc),
        }),
        ambient_temperature: Temperature::new(300.0),
        ambient_htc: HeatTransferCoefficient::new(0.0),
        external_radiant_power: Power::new(0.0),
    }
}

fn network(
    gas: &MixtureGas,
    axis: &SpeciesAxis,
    fuel: Fuel,
    vented: bool,
) -> Result<ZoneNetwork, Box<dyn Error>> {
    let mut oxidizer = vec![0.0; gas.species_count()];
    let mean = 0.1 * axis.molar_masses[axis.oxygen] + 0.9 * axis.molar_masses[axis.nitrogen];
    oxidizer[axis.oxygen] = 0.1 * axis.molar_masses[axis.oxygen] / mean;
    oxidizer[axis.nitrogen] = 0.9 * axis.molar_masses[axis.nitrogen] / mean;
    let room = TwoZone::new(
        gas.clone(),
        Area::new(0.1),
        Length::new(0.2),
        Length::new(0.0),
        Pressure::new(101325.0),
        [Temperature::new(fuel.initial_temperature_k()); 2],
        [oxidizer.clone(), oxidizer.clone()],
        || false,
    )?;
    let vents = if vented {
        let discharge = DischargeCoefficient::new(
            Dimensionless::new(0.7),
            SourceProvenance::new(
                "authored coupling fixture nominal opening coefficient",
                "coupled-compartment-opening-v1",
            ),
            0.0,
            ToleranceBasis::EngineeringAllowance,
        )?;
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Ambient,
            opening: WallOpening::new(
                Length::new(0.0),
                Length::new(0.2),
                Length::new(0.001),
                discharge,
            )?,
        }]
    } else {
        vec![]
    };
    Ok(ZoneNetwork::new(
        vec![Room {
            gas: room,
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(300.0),
            mass_fractions: oxidizer,
        },
        vents,
        vec![],
        Acceleration::new(9.81),
        || false,
    )?)
}

fn inventory(room: &TwoZone) -> Vec<f64> {
    room.layer(Layer::Lower)
        .species_masses()
        .iter()
        .zip(room.layer(Layer::Upper).species_masses())
        .map(|(lower, upper)| lower + upper)
        .collect()
}

fn run(
    vented: bool,
    end: f64,
    maximum_dt: f64,
    cells: usize,
    kinetic: bool,
    method: IsochoricKineticsMethod,
    fuel: Fuel,
    out: &mut impl Write,
) -> Result<(), Box<dyn Error>> {
    let mechanism = fuel.mechanism()?;
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(
            Temperature::new(300.0),
            Temperature::new(fuel.maximum_temperature_k()),
        ),
        || false,
    )?;
    let axis = SpeciesAxis::new(&gas, fuel)?;
    let ns = gas.species_count();
    let mut pure_fuel = vec![0.0; ns];
    pure_fuel[axis.fuel] = 1.0;
    let volatile = VolatileMixture::new(&gas, &pure_fuel, || false)?;
    let char_law = Caloric {
        h300: 0.0,
        cp: 1000.0,
    };
    let law = Caloric {
        h300: (1.0 - fuel.char_yield())
            * volatile
                .specific_enthalpy(Temperature::new(300.0), || false)?
                .value()
            + fuel.char_yield() * char_law.h300
            - 1e6,
        cp: 1500.0,
    };
    let wall_law = Caloric {
        h300: 0.0,
        cp: 1000.0,
    };
    let mut slab = solid(&law, &char_law, &volatile, cells, fuel)?;
    let mut wall = Wall1d::new(
        &wall_law,
        WallGeometry {
            area: Area::new(AREA),
            thickness: Length::new(0.02),
            density: Density::new(500.0),
            conductivity: ThermalConductivity::new(1.0),
            cells: 4,
        },
        SpecificEnthalpy::new(0.0),
        || false,
    )?;
    let cards = [surface_card(0.9)?, surface_card(0.7)?];
    let mut network = network(&gas, &axis, fuel, vented)?;
    let controls = NetworkControl {
        minimum_step_s: 1e-9,
        ..Default::default()
    };
    let mut chemistry = ZoneKinetics::new(
        &network.rooms()[0].gas,
        mechanism,
        IsochoricKineticsSolve {
            method,
            minimum_temperature_k: 300.0,
            maximum_temperature_k: fuel.maximum_temperature_k(),
            maximum_step_s: maximum_dt,
            initial_step_s: maximum_dt.min(1e-6),
            relative_tolerance: 1e-5,
            ..Default::default()
        },
        || false,
    )?;
    let molecular = MolecularPlanck::new(
        &gas,
        "only H2O/CO2/CO/CH4 in the original mechanism participate; all other source species are treated as transparent; no soot-production closure",
    )?;
    let enclosure = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(AREA)?,
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(AREA),
            segments: vec![OpticalSegment {
                layer: Layer::Upper,
                length: Length::new(0.2),
            }],
        }],
        "authored reciprocal normal chord between two 0.01 m2 patches; 0.2 m effective gas path, not a complete room enclosure",
        1e-12,
        || false,
    )?;
    let slab_control = |dt| SlabSurfaceSolve {
        pyrolysis: PyrolysisControls {
            maximum_dt_s: dt,
            ..Default::default()
        },
        ..Default::default()
    };
    let initial_boundary = boundary(&network.rooms()[0].gas, 10.0);
    let mut work = enclosure.surface_workspace(&[
        ThermalSurface::Slab {
            slab: &slab,
            boundary: SlabBoundary {
                gas: initial_boundary.gas,
                external_radiant_power: Power::new(EXTERNAL_POWER),
            },
            solve: slab_control(maximum_dt),
        },
        ThermalSurface::Wall {
            wall: &wall,
            boundary: boundary(&network.rooms()[0].gas, 8.0),
        },
    ])?;
    let mut input_powers = vec![Power::new(EXTERNAL_POWER), Power::new(0.0)];
    let initial_species = inventory(&network.rooms()[0].gas);
    let initial_solid_mass = slab.mass_kg();
    let initial_energy = slab.energy_j()
        + wall.simulation().energy_j()
        + network.rooms()[0].gas.physical_internal_energy().value();
    let mut vent_species = vec![0.0; ns];
    let mut chemical_species = vec![0.0; ns];
    let mut emitted = 0.0;
    let mut external_energy = 0.0;
    let mut steps = 0;
    let mut halvings = 0;
    let mut chemical_steps = 0;
    let mut chemical_rejections = 0;
    let mut chemical_direct_fallbacks = 0;
    let mut remaps = 0;
    let mut max_mass = 0.0_f64;
    let mut max_energy = 0.0_f64;
    let mut max_element = 0.0_f64;
    let mut max_pressure = 0.0_f64;
    let mut max_surface_power = 0.0_f64;
    let mut next_output = 0.0;
    let name = if vented { "vented" } else { "closed" };
    let mode = if kinetic { "kinetic" } else { "frozen" };
    let fuel_name = fuel.name();
    let method_name = match method {
        IsochoricKineticsMethod::BackwardEuler => "be",
        IsochoricKineticsMethod::Sdirk2 => "sdirk2",
    };
    while slab.time_s() < end {
        if steps >= 50_000 {
            return Err("coupled example accepted-step budget exhausted".into());
        }
        let remaining = end - slab.time_s();
        let dt = if remaining > maximum_dt && remaining - maximum_dt < controls.minimum_step_s {
            0.5 * remaining
        } else {
            maximum_dt.min(remaining)
        };
        let (report, candidate, reacted) = network.step_with_pressure_sources(Time::new(dt), controls,
            |rooms, attempted, rates, cancel| {
                let room = &rooms[0].gas;
                let opacity = molecular.evaluate_room(room, &mut *cancel)?;
                let absorption = [
                    LayerAbsorption::new(AbsorptionCoefficient::new(0.0), "absent lower gas layer")?,
                    opacity[1].ok_or(fs_airflow::compartment::NetworkError::Invalid("upper layer absent"))?
                        .gray_absorption("local emission Planck mean used as gray absorption over both patch spectra and the finite 0.2 m path; approximate optical closure")?,
                ];
                let surfaces = [
                    ThermalSurface::Slab { slab: &slab, boundary: SlabBoundary { gas: boundary(room, 10.0).gas,
                        external_radiant_power: Power::new(EXTERNAL_POWER) }, solve: slab_control(attempted.value()) },
                    ThermalSurface::Wall { wall: &wall, boundary: boundary(room, 8.0) },
                ];
                let candidate = enclosure.candidate_room_surfaces(room, &surfaces, attempted,
                    WallOptics::Gray { absorption: &absorption, solve: RadiationSolve::default() },
                    |i, t| SurfaceEmissivity::from_card(if i == 0 { "solid" } else { "wall" }, &cards[i], t.value(), SelectionPolicy::SingleClaimOnly),
                    Some(&input_powers), EnclosureWallSolve::default(), &mut work, &mut *cancel)?;
                candidate.report.add_gas_heat(&mut rates[0])?;
                let ThermalSurfaceCandidate::Slab(solid) = &candidate.surfaces[0] else { unreachable!("fixed surface ordering") };
                for release in solid.report.pyrolysis.releases.values() {
                    volatile.release_rates(room, *release, attempted, VolatileTransferControl::default(), &mut *cancel)?
                        .add_to(Layer::Upper, &mut rates[0])?;
                }
                Ok(candidate)
            },
            |rooms, attempted, _, additional, cancel| {
                if !kinetic { return Ok(None); }
                let result = chemistry.apply_split(&mut rooms[0].gas, attempted, [true; 2], controls.zone, &mut *cancel)?;
                for rate in &result.chemistry.chemical_rates { additional.add(rate)?; }
                Ok(Some(result))
            }, || false)?;
        let dt = report.elapsed.value();
        let mut surfaces = candidate.surfaces.into_iter();
        let ThermalSurfaceCandidate::Slab(next_slab) = surfaces.next().expect("solid candidate")
        else {
            unreachable!("fixed surface ordering")
        };
        let ThermalSurfaceCandidate::Wall(next_wall) = surfaces.next().expect("wall candidate")
        else {
            unreachable!("fixed surface ordering")
        };
        let emitted_step = next_slab.report.pyrolysis.volatile_mass_kg;
        let volatile_energy = next_slab.report.pyrolysis.volatile_enthalpy_j;
        let slab_face = next_slab.surface_temperature.value();
        let wall_face = next_wall.report.surface_temperatures[0].value();
        remaps += next_slab.report.remaps;
        // Publish only the matching accepted payload; no thermal clock advances twice.
        slab = next_slab.slab;
        wall = next_wall.wall;
        input_powers = candidate.input_powers;
        steps += 1;
        halvings += report.halvings;
        emitted += emitted_step;
        let mut chemistry_change = vec![0.0; ns];
        let mut chemistry_energy = 0.0;
        if let Some(result) = reacted {
            for rate in result.chemistry.chemical_rates {
                chemistry_energy += dt
                    * (rate.referenced_energy_w
                        + rate.species_kg_s[..ns]
                            .iter()
                            .zip(gas.energy_references())
                            .map(|(m, e)| m * e)
                            .sum::<f64>());
                for i in 0..ns {
                    chemistry_change[i] += dt * rate.species_kg_s[i];
                }
            }
            for layer in result.chemistry.layers.into_iter().flatten() {
                chemical_steps += layer.reactor.accepted_steps;
                chemical_rejections += layer.reactor.rejected_steps;
                chemical_direct_fallbacks += layer.reactor.direct_fallbacks;
            }
        }
        // Remove INTERNAL transfers from the network's aggregate source account.
        // The remaining external gas stream is vent enthalpy. Add the heater and
        // any wall outer-boundary input once to the combined storage account.
        let vent_energy = report.external_physical_energy_j
            - volatile_energy
            - chemistry_energy
            - candidate
                .report
                .gas_heat_energy
                .iter()
                .map(|e| e.value())
                .sum::<f64>();
        external_energy += vent_energy
            + candidate.report.external_radiant_energy.value()
            + candidate.report.ambient_heat_energy.value();
        for i in 0..ns {
            chemical_species[i] += chemistry_change[i];
            vent_species[i] += report.external_species_change_kg[i]
                - chemistry_change[i]
                - if i == axis.fuel { emitted_step } else { 0.0 };
        }
        let room = &network.rooms()[0].gas;
        if slab.time_s() != room.time().value() || wall.time() != room.time() {
            return Err("accepted gas/solid/wall clocks differ".into());
        }
        let masses = inventory(room);
        let mass_defect = slab.mass_kg() - initial_solid_mass
            + (0..ns)
                .map(|i| masses[i] - initial_species[i] - vent_species[i])
                .sum::<f64>();
        let energy_defect = slab.energy_j()
            + wall.simulation().energy_j()
            + room.physical_internal_energy().value()
            - initial_energy
            - external_energy;
        let elements = gas.source().elemental_matrix();
        let virgin_mass: f64 = slab
            .state()
            .inventories()
            .values()
            .map(|s| s.virgin_mass_kg)
            .sum();
        let char_mass: f64 = slab
            .state()
            .inventories()
            .values()
            .map(|s| s.char_mass_kg)
            .sum();
        for e in 0..elements.row_count() {
            let defect = (virgin_mass - initial_solid_mass) * axis.virgin_atoms_per_kg[e]
                + char_mass * axis.char_atoms_per_kg[e]
                + (0..ns)
                    .map(|i| {
                        elements.get(e, i).unwrap() as f64
                            * (masses[i] - initial_species[i] - vent_species[i])
                            / axis.molar_masses[i]
                    })
                    .sum::<f64>();
            max_element = max_element.max(defect.abs());
        }
        let species_defect = (0..ns)
            .map(|i| {
                (masses[i]
                    - initial_species[i]
                    - vent_species[i]
                    - chemical_species[i]
                    - if i == axis.fuel { emitted } else { 0.0 })
                .abs()
            })
            .fold(0.0, f64::max);
        max_mass = max_mass.max(mass_defect.abs());
        max_energy = max_energy.max(energy_defect.abs());
        max_pressure = max_pressure.max(report.maximum_pressure_residual_pa);
        max_surface_power = max_surface_power.max(candidate.report.maximum_power_residual_w);
        if masses.iter().any(|m| !m.is_finite() || *m < 0.0)
            || room.layer(Layer::Lower).mass_kg() != 0.0
            || species_defect > 1e-10
            || max_mass > 1e-10
            || max_element > 1e-7
            || max_energy > 0.02
        {
            return Err(format!("combined balance refused: species={species_defect:e}kg mass={mass_defect:e}kg energy={energy_defect:e}J elements={max_element:e}mol").into());
        }
        if slab.time_s() >= next_output || slab.time_s() >= end {
            let opacity = molecular.evaluate_room(room, || false)?[1]
                .expect("present upper layer")
                .absorption
                .value();
            let time = slab.time_s();
            let gas_temperature = room.layer(Layer::Upper).temperature().value();
            let pressure = room.pressure().value();
            let solid_mass = slab.mass_kg();
            let gas_fuel = masses[axis.fuel];
            let hydrogen = masses[axis.h2];
            let water = masses[axis.water];
            let oxygen = masses[axis.oxygen];
            let hydrogen_atom = masses[axis.hydrogen_atom];
            let hydroxyl = masses[axis.hydroxyl];
            let carbon_monoxide = SpeciesAxis::mass(&masses, axis.co);
            let carbon_dioxide = SpeciesAxis::mass(&masses, axis.co2);
            let nitric_oxide = SpeciesAxis::mass(&masses, axis.no);
            let recession = slab.recession_m();
            let gas_radiation = candidate.report.radiation.gas_heat_power()[1].value();
            writeln!(
                out,
                "{name},{mode},{method_name},{fuel_name},{time:.9},{gas_temperature:.9},{pressure:.9},{slab_face:.9},{wall_face:.9},{solid_mass:.12e},{char_mass:.12e},{emitted:.12e},{gas_fuel:.12e},{hydrogen:.12e},{water:.12e},{oxygen:.12e},{hydrogen_atom:.12e},{hydroxyl:.12e},{carbon_monoxide:.12e},{carbon_dioxide:.12e},{nitric_oxide:.12e},{opacity:.12e},{recession:.12e},{gas_radiation:.12e},{mass_defect:.12e},{energy_defect:.12e},{max_element:.12e}"
            )?;
            next_output = slab.time_s() + end / 10.0;
            out.flush()?;
        }
    }
    eprintln!(
        "case={name} chemistry={mode} method={method_name} fuel={fuel_name} end_s={end} max_dt_s={maximum_dt} solid_cells={cells} accepted_steps={steps} network_halvings={halvings} accepted_chemical_steps={chemical_steps} chemical_rejections={chemical_rejections} chemical_direct_fallbacks={chemical_direct_fallbacks} remaps={remaps} emitted_fuel_kg={emitted:.12e} net_chemical_fuel_consumption_kg={:.12e} net_vent_mass_kg={:.12e} max_mass_residual_kg={max_mass:.12e} max_energy_residual_j={max_energy:.12e} max_element_residual_mol={max_element:.12e} max_pressure_residual_pa={max_pressure:.12e} max_surface_power_residual_w={max_surface_power:.12e}",
        -chemical_species[axis.fuel],
        vent_species.iter().sum::<f64>()
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args
        .first()
        .is_some_and(|a| matches!(a.as_str(), "-h" | "--help"))
    {
        println!(
            "coupled_compartment_fire [all|closed|vented] [end_s=0.1] [max_dt_s=0.001] [solid_cells=12] [kinetic|frozen] [be|sdirk2] [hydrogen|methane|charring-methane]\nManufactured fuel-releasing solid, finite wall and approximate molecular gray radiation. Hydrogen uses complete h2o2 (10 species/29 reactions) at 1000 K; methane uses complete GRI-Mech 3.0 (53 species/325 reactions) at 1500 K. Initial oxidizer is 10% O2 / 90% N2 by mole; no radical seed. Charring methane retains 20% condensed-carbon residue with separately declared density/conductivity. Backward Euler is default; SDIRK2 changes only the reactor method. No experimental fire, wood calibration, soot production, global entropy theorem or globally second-order coupling claim."
        );
        return Ok(());
    }
    if args.len() > 7 {
        return Err("expected at most seven arguments".into());
    }
    let selected = args.first().map(String::as_str).unwrap_or("all");
    let end = args
        .get(1)
        .map(|s| s.parse::<f64>())
        .transpose()?
        .unwrap_or(0.1);
    let dt = args
        .get(2)
        .map(|s| s.parse::<f64>())
        .transpose()?
        .unwrap_or(0.001);
    let cells = args
        .get(3)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(12);
    let mode = args.get(4).map(String::as_str).unwrap_or("kinetic");
    let method = match args.get(5).map(String::as_str).unwrap_or("be") {
        "be" => IsochoricKineticsMethod::BackwardEuler,
        "sdirk2" => IsochoricKineticsMethod::Sdirk2,
        _ => return Err("chemistry method must be be or sdirk2".into()),
    };
    let fuel = Fuel::parse(args.get(6).map(String::as_str).unwrap_or("hydrogen"))?;
    if !["all", "closed", "vented"].contains(&selected)
        || !["kinetic", "frozen"].contains(&mode)
        || !end.is_finite()
        || !(1e-5..=5.0).contains(&end)
        || !dt.is_finite()
        || !(1e-6..=0.01).contains(&dt)
        || !(2..=96).contains(&cells)
        || end / dt > 50_000.0
    {
        return Err("require an admitted case/mode, 1e-5 <= end <= 5 s, 1e-6 <= max_dt <= 0.01 s, 2..96 cells and <=50000 requested steps; source domains may refuse longer cases".into());
    }
    eprintln!(
        "source={},calorics=explicit-bounded-300K-continuity,solid=manufactured-{},char_yield={},initial_gas_K={},heater_W={EXTERNAL_POWER},optics=TNF-Planck-means-explicit-gray-approximation,geometry=two-patch-normal-chord,global_time_order=first",
        fuel.source(),
        fuel.name(),
        fuel.char_yield(),
        fuel.initial_temperature_k()
    );
    let mut out = BufWriter::new(io::stdout().lock());
    writeln!(
        out,
        "case,chemistry,method,fuel,time_s,gas_temperature_k,pressure_pa,solid_face_k,wall_face_k,solid_mass_kg,char_mass_kg,emitted_fuel_kg,gas_fuel_kg,gas_H2_kg,gas_H2O_kg,gas_O2_kg,gas_H_kg,gas_OH_kg,gas_CO_kg,gas_CO2_kg,gas_NO_kg,endpoint_kappa_m_inv,recession_m,gas_radiation_heat_w,combined_mass_residual_kg,combined_energy_residual_j,max_element_residual_mol"
    )?;
    for (name, vented) in [("closed", false), ("vented", true)] {
        if selected == "all" || selected == name {
            run(
                vented,
                end,
                dt,
                cells,
                mode == "kinetic",
                method,
                fuel,
                &mut out,
            )?;
        }
    }
    Ok(())
}
