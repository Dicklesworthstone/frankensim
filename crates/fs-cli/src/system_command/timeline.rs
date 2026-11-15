//! `thermal_timeline`: enthalpy cells, their links and their boundaries as
//! one scene element, stepped by the fs-thermal timeline owner.
//!
//! The element is [`ThermalTimeline`] over a [`DepositionSimulation`] whose
//! cells are all born at the start of the run. The owner's explicit rule
//! conducts heat through the links, exchanges it at the boundaries and
//! stores it on each cell's chart; nothing here conducts, stores or
//! integrates heat. A chart is the scene's authored sensible and latent
//! chart or a declared source caloric, the same `EquilibriumEnthalpyPhaseCurve`
//! an `enthalpy_body` uses, which `fs-material` makes a `ThermalLaw`.
//!
//! The kernel keeps component states as `'static` values and the timeline
//! borrows its laws, so the laws of a timeline are leaked: they live as long
//! as the process that runs the scene.

use std::collections::BTreeMap;

use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_ports::component::{Component, PortSpec};
use fs_system::{NodeId, SystemBuilder};
use fs_thermal::component::{ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, Reservoir, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

use super::{J, Refusal, authored_enthalpy_curve, fields, materials, number, number_or, positive};

/// What the kernel holds for a thermal timeline.
pub(super) type State = TimelineState<'static>;

/// Laws, cells, links and later events one timeline may declare.
const MAX_LAWS: usize = 64;
const MAX_CELLS: usize = 4096;
const MAX_LINKS: usize = 16_384;
const MAX_EVENTS: usize = 16_384;
/// Owner substeps one kernel step may take. The kernel cuts its steps to the
/// owner's stability bound, so a step normally takes one.
const MAX_SUBSTEPS: usize = 1024;
/// The owner's conservation gate on each step when the scene gives none.
const ENERGY_TOLERANCE_J: f64 = 1.0e-9;

/// A lowered timeline: its node and ports, the temperature each heat port
/// and published temperature starts at, and the source calorics it uses.
pub(super) struct Parts {
    pub(super) node: NodeId,
    pub(super) ports: Vec<PortSpec>,
    pub(super) port_temperatures_k: Vec<(String, f64)>,
    pub(super) material_sources: Vec<String>,
}

type Law = &'static EquilibriumEnthalpyPhaseCurve;

/// A `thermal_timeline` element: its cells born at the start with their
/// links, the later events, and the timeline with its heat ports, its
/// environment port and its published temperatures.
pub(super) fn timeline(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    materials: &materials::Materials,
    what: &str,
) -> Result<Parts, Refusal> {
    fields(
        definition,
        &[
            "kind",
            "laws",
            "cells",
            "links",
            "events",
            "ports",
            "environment",
            "publish",
            "max_step_s",
            "energy_tolerance_j",
        ],
        what,
    )?;
    let mut material_sources = Vec::new();
    let laws = laws(definition, materials, &mut material_sources, what)?;
    let cell_rows = rows(definition, "cells", 1, MAX_CELLS, what)?;
    let mut schedule = Vec::with_capacity(cell_rows.len());
    let mut start_k = BTreeMap::new();
    for (index, row) in cell_rows.iter().enumerate() {
        let (cell, kelvin) = cell(row, &laws, &format!("{what}.cells[{index}]"))?;
        if start_k.insert(cell.id, kelvin).is_some() {
            return Err(format!(
                "{what}.cells[{index}]: cell {} is declared twice",
                cell.id.0
            )
            .into());
        }
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        });
    }
    for (index, row) in rows(definition, "links", 0, MAX_LINKS, what)?
        .iter()
        .enumerate()
    {
        let at = format!("{what}.links[{index}]");
        let index = u64::try_from(index).map_err(|_| format!("{at}: too many links"))?;
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(link(row, Some(index), &at)?),
        });
    }
    for (index, row) in rows(definition, "events", 0, MAX_EVENTS, what)?
        .iter()
        .enumerate()
    {
        schedule.push(event(row, &format!("{what}.events[{index}]"))?);
    }
    let ports = cell_list(definition, "ports", what)?;
    let published = cell_list(definition, "publish", what)?;
    let environment = match definition.get("environment") {
        None => false,
        Some(J::Bool(environment)) => *environment,
        Some(_) => return Err(format!("{what}: environment is true or false").into()),
    };
    let config = StepConfig {
        max_dt_s: positive(definition, "max_step_s", what)?,
        energy_tolerance_j: number_or(definition, "energy_tolerance_j", ENERGY_TOLERANCE_J, what)?,
    };
    let mut port_temperatures_k = Vec::new();
    for (prefix, cells) in [("cell", &ports), ("temperature", &published)] {
        for cell in cells {
            // A cell the timeline does not hold is the owner's refusal at
            // its start, below; it starts at no temperature here.
            if let Some(&kelvin) = start_k.get(cell) {
                port_temperatures_k.push((format!("{prefix}{}", cell.0), kelvin));
            }
        }
    }
    let component = ThermalTimeline::new(ports, environment, config)
        .and_then(|component| component.publishing(published))
        .map_err(|error| format!("{what}: {error}"))?;
    let charts: Vec<(MaterialId, &'static dyn ThermalLaw)> = laws
        .iter()
        .zip(0_u64..)
        .map(|(&(_, law), index)| (MaterialId(index), law as &'static dyn ThermalLaw))
        .collect();
    let budget = DepositionBudget {
        max_cells: MAX_CELLS,
        max_links: MAX_LINKS + MAX_EVENTS,
        max_events: MAX_CELLS + MAX_LINKS + MAX_EVENTS,
        max_substeps: MAX_SUBSTEPS,
    };
    let simulation = DepositionSimulation::new(schedule, &charts, budget)
        .map_err(|error| format!("{what}: {error}"))?;
    let state = component
        .start(&simulation)
        .map_err(|error| format!("{what}: {error}"))?;
    let ports = component.ports().to_vec();
    let node = builder
        .add(name, component, state)
        .map_err(|error| Refusal::Admission(error.to_string()))?;
    Ok(Parts {
        node,
        ports,
        port_temperatures_k,
        material_sources,
    })
}

/// Every cell's temperature and specific enthalpy, in the order of their
/// identifiers, and the mean rate heat entered through the boundaries over
/// the last step when the timeline has an environment port.
pub(super) fn values(state: Option<&State>, environment_w: Option<f64>) -> Vec<(String, f64)> {
    let mut values = Vec::new();
    if let Some(state) = state {
        let simulation = state.simulation();
        for (id, cell) in simulation.cells() {
            let kelvin = simulation.temperature_k(*id).unwrap_or(f64::NAN);
            values.push((format!("cell{}.temperature_k", id.0), kelvin));
            values.push((
                format!("cell{}.specific_enthalpy_j_kg", id.0),
                cell.specific_h_j_kg,
            ));
        }
    }
    if let Some(watts) = environment_w {
        values.push(("environment_heat_w".to_owned(), watts));
    }
    values
}

/// The named laws, each an authored chart or a declared source caloric.
fn laws(
    definition: &J,
    materials: &materials::Materials,
    sources: &mut Vec<String>,
    what: &str,
) -> Result<Vec<(String, Law)>, String> {
    let inner = format!("{what}.laws");
    let named = definition
        .get("laws")
        .and_then(J::as_object)
        .filter(|named| (1..=MAX_LAWS).contains(&named.len()))
        .ok_or_else(|| format!("{inner}: an object of 1 to {MAX_LAWS} named charts"))?;
    let mut laws = Vec::with_capacity(named.len());
    for (name, law) in named {
        let at = format!("{inner}.{name}");
        let curve = if law.get("material").is_some() {
            fields(law, &["material"], &at)?;
            let source = law
                .str_field("material")
                .ok_or_else(|| format!("{at}.material: expected a material name"))?;
            let prepared = materials
                .get(source)
                .ok_or_else(|| format!("{at}.material: unknown material {source:?}"))?;
            if !sources.iter().any(|known| known == source) {
                sources.push(source.to_owned());
            }
            prepared.caloric.clone()
        } else {
            fields(law, &["phase", "enthalpy_curve"], &at)?;
            authored_enthalpy_curve(law, &at)?
        };
        let law: Law = Box::leak(Box::new(curve));
        laws.push((name.clone(), law));
    }
    Ok(laws)
}

/// The array under `key`, of `minimum` to `maximum` items; absent is empty.
fn rows<'j>(
    definition: &'j J,
    key: &str,
    minimum: usize,
    maximum: usize,
    what: &str,
) -> Result<&'j [J], String> {
    match definition.get(key) {
        None if minimum == 0 => Ok(&[]),
        value => value
            .and_then(J::as_array)
            .filter(|rows| (minimum..=maximum).contains(&rows.len()))
            .ok_or_else(|| format!("{what}.{key}: an array of {minimum} to {maximum} items")),
    }
}

/// A cell born at the start, and the temperature its chart gives it.
fn cell(value: &J, laws: &[(String, Law)], what: &str) -> Result<(DepositedCell, f64), String> {
    fields(
        value,
        &[
            "id",
            "law",
            "mass_kg",
            "specific_enthalpy_j_kg",
            "volume_m3",
            "surface_area_m2",
            "start_m",
            "end_m",
            "boundary",
        ],
        what,
    )?;
    let id = whole(value, "id", what)?;
    let name = value
        .str_field("law")
        .ok_or_else(|| format!("{what}.law: expected the name of a law"))?;
    let (index, law) = laws
        .iter()
        .zip(0_u64..)
        .find_map(|((known, law), index)| (known == name).then_some((index, *law)))
        .ok_or_else(|| format!("{what}.law: no law named {name:?}"))?;
    let specific_h_j_kg = number(value, "specific_enthalpy_j_kg", what)?;
    let kelvin = ThermalLaw::temperature_k(law, specific_h_j_kg)
        .map_err(|error| format!("{what}.specific_enthalpy_j_kg: {error}"))?;
    let cell = DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(index),
        mass_kg: number(value, "mass_kg", what)?,
        specific_h_j_kg,
        geometry: CellGeometry {
            start_m: point(value, "start_m", what)?,
            end_m: point(value, "end_m", what)?,
            volume_m3: number(value, "volume_m3", what)?,
            surface_area_m2: number(value, "surface_area_m2", what)?,
        },
        boundary: boundary(
            value.get("boundary"),
            kelvin,
            &[],
            &format!("{what}.boundary"),
        )?,
    };
    Ok((cell, kelvin))
}

/// A cell's boundary. Absent, the cell is insulated: no film, reservoir or
/// source, at an ambient temperature that is then never used (the owner
/// needs a positive one), its starting temperature.
fn boundary(
    value: Option<&J>,
    insulated_ambient_k: f64,
    also: &[&str],
    what: &str,
) -> Result<Boundary, String> {
    let Some(value) = value else {
        return Ok(Boundary {
            ambient_k: insulated_ambient_k,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        });
    };
    let mut allowed = vec!["ambient_k", "htc_w_m2_k", "source_w", "reservoir"];
    allowed.extend_from_slice(also);
    fields(value, &allowed, what)?;
    let reservoir = match value.get("reservoir") {
        None => None,
        Some(reservoir) => {
            let at = format!("{what}.reservoir");
            fields(
                reservoir,
                &["area_m2", "conductance_w_k", "temperature_k"],
                &at,
            )?;
            Some(Reservoir {
                area_m2: number(reservoir, "area_m2", &at)?,
                conductance_w_k: number(reservoir, "conductance_w_k", &at)?,
                temperature_k: number(reservoir, "temperature_k", &at)?,
            })
        }
    };
    Ok(Boundary {
        ambient_k: number(value, "ambient_k", what)?,
        htc_w_m2_k: number_or(value, "htc_w_m2_k", 0.0, what)?,
        reservoir,
        source_w: number_or(value, "source_w", 0.0, what)?,
    })
}

/// A link between two cells, covering `area_m2` of each (the owner needs a
/// positive area); its id is `default_id` when it gives none.
fn link(value: &J, default_id: Option<u64>, what: &str) -> Result<ThermalLink, String> {
    fields(value, &["id", "a", "b", "conductance_w_k", "area_m2"], what)?;
    let id = match (value.get("id"), default_id) {
        (None, Some(id)) => id,
        (None, None) => return Err(format!("{what}: a link made by an event names its id")),
        (Some(_), _) => whole(value, "id", what)?,
    };
    Ok(ThermalLink {
        id,
        a: CellId(whole(value, "a", what)?),
        b: CellId(whole(value, "b", what)?),
        area_m2: number(value, "area_m2", what)?,
        conductance_w_k: number(value, "conductance_w_k", what)?,
    })
}

/// A later event: a link made, a link broken, or a cell's boundary replaced.
fn event(value: &J, what: &str) -> Result<TimedEvent, String> {
    let shape = || format!("{what}: an event is {{time_s, contact | separation | boundary}}");
    let members = value.as_object().ok_or_else(shape)?;
    let time_s = number(value, "time_s", what)?;
    let mut kinds = members.iter().filter(|(key, _)| key != "time_s");
    let (Some((kind, body)), None) = (kinds.next(), kinds.next()) else {
        return Err(shape());
    };
    let at = format!("{what}.{kind}");
    let event = match kind.as_str() {
        "contact" => Event::Contact(link(body, None, &at)?),
        "separation" => Event::Separation(whole(value, "separation", what)?),
        "boundary" => Event::Boundary(
            CellId(whole(body, "cell", &at)?),
            boundary(Some(body), 0.0, &["cell"], &at)?,
        ),
        _ => return Err(shape()),
    };
    Ok(TimedEvent { time_s, event })
}

/// Cell identifiers listed under `key`; absent is none.
fn cell_list(definition: &J, key: &str, what: &str) -> Result<Vec<CellId>, String> {
    rows(definition, key, 0, MAX_CELLS, what)?
        .iter()
        .map(|row| whole_number(row).map(CellId))
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| format!("{what}.{key}: an array of cell identifiers"))
}

/// A cell or link identifier under `key`: a whole number from 0 to 2^53.
fn whole(value: &J, key: &str, what: &str) -> Result<u64, String> {
    value
        .get(key)
        .and_then(whole_number)
        .ok_or_else(|| format!("{what}: {key} must be a whole number from 0 to 2^53"))
}

fn whole_number(value: &J) -> Option<u64> {
    const LARGEST: f64 = 9_007_199_254_740_992.0;
    value
        .as_f64()
        .filter(|number| (0.0..=LARGEST).contains(number) && number.fract() == 0.0)
        // In range and whole, so the conversion is exact.
        .map(|number| number as u64)
}

/// A point of the cell's centreline [m]: three finite numbers.
fn point(value: &J, key: &str, what: &str) -> Result<[f64; 3], String> {
    match value.get(key).and_then(J::as_array) {
        Some([x, y, z]) => match (x.as_f64(), y.as_f64(), z.as_f64()) {
            (Some(x), Some(y), Some(z)) if [x, y, z].iter().all(|v| v.is_finite()) => Ok([x, y, z]),
            _ => Err(format!("{what}: {key} is three finite numbers [m]")),
        },
        _ => Err(format!("{what}: {key} is three finite numbers [m]")),
    }
}
