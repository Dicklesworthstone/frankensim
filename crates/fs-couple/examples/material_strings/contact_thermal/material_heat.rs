//! Source fiber/string-body heat in the existing finite thermal renderer.
//! Geometry and loss channels come from the original mechanical source. Profile
//! cells supply a caloric discretization, never substitute damping or a heat split.
use super::*;
use fs_couple::render::schedule::force::coupled::coefficients::heat::MaterialHeatBudget;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::{
    material::MaterialRelaxationHeatDestination, rates::MaterialRelaxationTemperature,
};

#[derive(Clone, Debug)]
pub(super) struct FiberRoute {
    strip: usize,
    layer: usize,
    cell: CellId,
    pub(super) source_temperature: bool,
}

#[derive(Clone, Debug)]
pub(super) struct StringRoute {
    string: usize,
    layer: usize,
    cell: CellId,
    source_temperature: bool,
}

#[derive(Clone, Debug)]
pub(super) struct Routes {
    budget: MaterialHeatBudget,
    pub(super) fibers: Vec<FiberRoute>,
    strings: Vec<StringRoute>,
}

impl Routes {
    pub(super) fn len(&self) -> usize { self.fibers.len() + self.strings.len() }
    pub(super) fn has_temperature_rates(&self) -> bool {
        self.fibers.iter().any(|route| route.source_temperature)
            || self.strings.iter().any(|route| route.source_temperature)
    }
}

pub(super) fn parse(
    rows: &mut Rows<'_>, cells: &[Cell], films: &[(usize, CellId)],
) -> Result<Routes> {
    let mut row = rows.row("material-heat-budget")?;
    let max_dimension: usize = row.parse()?;
    let max_dense_bytes: usize = row.parse()?;
    let max_multiply_terms: u64 = row.parse()?;
    row.end()?;
    if !(1..=256).contains(&max_dimension) || max_dense_bytes == 0 || max_multiply_terms == 0 {
        return Err("material heat requires explicit positive dimension, byte and work budgets".into());
    }
    let count = rows.count("fiber-heat", MAX_CELLS)?;
    let mut fibers = Vec::<FiberRoute>::with_capacity(count);
    for _ in 0..count {
        let mut row = rows.row("fiber-body")?;
        let route = FiberRoute {
            strip: row.parse()?, layer: row.parse()?, cell: CellId(row.parse()?),
            source_temperature: match row.words.next() {
                None => false,
                Some("source-temperature") => true,
                _ => return Err("fiber-body accepts only the explicit source-temperature rate option".into()),
            },
        };
        row.end()?;
        if !cells.iter().any(|cell| cell.id == route.cell)
            || films.iter().any(|&(_, cell)| cell == route.cell)
            || fibers.iter().any(|other| other.cell == route.cell
                || (other.strip, other.layer) == (route.strip, route.layer))
        {
            return Err("each fiber material region needs one distinct existing body cell and cannot alias a film".into());
        }
        fibers.push(route);
    }
    let mut strings = Vec::<StringRoute>::new();
    if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
        == Some("string-heat") {
        for _ in 0..rows.count("string-heat", MAX_CELLS)? {
            let mut row = rows.row("string-body")?;
            let route = StringRoute {
                string: row.parse()?, layer: row.parse()?, cell: CellId(row.parse()?),
                source_temperature: match row.words.next() {
                    None => false,
                    Some("source-temperature") => true,
                    _ => return Err("string-body accepts only the explicit source-temperature rate option".into()),
                },
            };
            row.end()?;
            if !cells.iter().any(|cell| cell.id == route.cell)
                || films.iter().any(|&(_, cell)| cell == route.cell)
                || fibers.iter().any(|fiber| fiber.cell == route.cell)
                || strings.iter().any(|other| other.cell == route.cell
                    || (other.string, other.layer) == (route.string, route.layer)) {
                return Err("each string material layer needs one distinct existing body cell and cannot alias a film or fiber".into());
            }
            strings.push(route);
        }
    }
    let count = fibers.len() + strings.len();
    // The original integral checks each complete operator. Reserve a share of
    // the declared aggregate resources for EVERY retained receiving body.
    if count > 0 && (max_dense_bytes / count == 0 || max_multiply_terms / count as u64 == 0) {
        return Err("material heat aggregate resources cannot cover all declared bodies".into());
    }
    Ok(Routes {
        // One bending block, one torsional block, and the same region's native
        // viscosity rows all consume this one original aggregate budget.
        budget: MaterialHeatBudget { max_groups: 3, max_dimension,
            max_dense_bytes, max_multiply_terms },
        fibers, strings,
    })
}

impl Prepared<'_> {
    pub(super) fn material_routes<'a>(
        &'a self, ready: &Ready, network: &DepositionSimulation<'_>, gate: &CancelGate,
    ) -> Result<(Vec<MaterialRelaxationHeatDestination<'a>>, Vec<MaterialRelaxationTemperature>)> {
        let Some(declaration) = &self.profile.material_heat else { return Ok((Vec::new(), Vec::new())); };
        let count = declaration.len();
        if count == 0 { return Ok((Vec::new(), Vec::new())); }
        let budget = MaterialHeatBudget {
            max_dense_bytes: declaration.budget.max_dense_bytes / count,
            max_multiply_terms: declaration.budget.max_multiply_terms / count as u64,
            ..declaration.budget
        };
        let mut routes = Vec::with_capacity(count);
        let mut rates = Vec::new();
        let bodies = declaration.fibers.iter().map(|route|
            (false, route.strip, route.layer, route.cell, route.source_temperature))
            .chain(declaration.strings.iter().map(|route|
                (true, route.string, route.layer, route.cell, route.source_temperature)));
        for (string, body, layer, id, temperature_rates) in bodies {
            poll(gate)?;
            let index = self.profile.cells.iter().position(|cell| cell.id == id)
                .ok_or("material heat lost its original receiving cell")?;
            let cell = &self.profile.cells[index];
            let source = &self.sources[index];
            let bank = ready.performance.bank();
            let (destination, rate) = if temperature_rates {
                let (heat, rate) = if string {
                    bank.string_material_rate_heat_destination(body, layer, id,
                        &source.card, &source.curve, "T", budget, gate)?
                } else {
                    bank.fiber_material_rate_heat_destination(body, layer, id,
                        &source.card, &source.curve, "T", budget, gate)?
                };
                (heat, Some(rate))
            } else {
                let heat = if string {
                    bank.string_material_heat_destination(body, layer, id,
                        &source.card, &source.curve, "T", budget, gate)?
                } else {
                    bank.fiber_material_heat_destination(body, layer, id,
                        &source.card, &source.curve, "T", budget, gate)?
                };
                (heat, None)
            };
            if destination.source.state_coordinate("T")
                .is_none_or(|temperature| temperature.to_bits() != cell.temperatures[1].to_bits())
            {
                return Err("material body and its finite thermal cell require one initial source temperature".into());
            }
            // Fixed routes keep their existing admission. Explicit rate routes
            // use the original source/volume/operator admission together below,
            // including exact branch-to-property roles from the fiber producer.
            if let Some(rate) = rate {
                rates.push(rate);
            } else {
                destination.admit::<std::convert::Infallible>(network)?;
            }
            routes.push(destination);
        }
        if !rates.is_empty() {
            // Read-only preflight: no advance and no output creation precedes
            // the actual native source, schedule and fixed-map admission.
            ready.performance.admit_material_relaxation_temperature(network, &routes, &rates, gate)?;
        }
        poll(gate)?;
        Ok((routes, rates))
    }
}
