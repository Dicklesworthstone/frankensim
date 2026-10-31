//! Source fiber-body heat declarations for the existing finite thermal renderer.
//! Geometry and loss channels come from the original mechanical source. Profile
//! cells supply a caloric discretization, never substitute damping or a heat split.
use super::*;
use fs_couple::render::schedule::force::coupled::coefficients::heat::MaterialHeatBudget;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::material::MaterialRelaxationHeatDestination;

#[derive(Clone, Debug)]
pub(super) struct FiberRoute {
    strip: usize,
    layer: usize,
    cell: CellId,
}

#[derive(Clone, Debug)]
pub(super) struct Routes {
    budget: MaterialHeatBudget,
    pub(super) fibers: Vec<FiberRoute>,
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
    // The original integral checks each complete operator. Reserve a share of
    // the declared aggregate resources for EVERY retained receiving body.
    if count > 0 && (max_dense_bytes / count == 0 || max_multiply_terms / count as u64 == 0) {
        return Err("material heat aggregate resources cannot cover all declared bodies".into());
    }
    Ok(Routes {
        budget: MaterialHeatBudget { max_groups: 1, max_dimension,
            max_dense_bytes, max_multiply_terms },
        fibers,
    })
}

impl Prepared {
    pub(super) fn material_routes<'a>(
        &'a self, ready: &Ready, network: &DepositionSimulation<'_>, gate: &CancelGate,
    ) -> Result<Vec<MaterialRelaxationHeatDestination<'a>>> {
        let Some(declaration) = &self.profile.material_heat else { return Ok(Vec::new()); };
        if declaration.fibers.is_empty() { return Ok(Vec::new()); }
        let count = declaration.fibers.len();
        let budget = MaterialHeatBudget {
            max_dense_bytes: declaration.budget.max_dense_bytes / count,
            max_multiply_terms: declaration.budget.max_multiply_terms / count as u64,
            ..declaration.budget
        };
        let mut routes = Vec::with_capacity(count);
        for route in &declaration.fibers {
            poll(gate)?;
            let index = self.profile.cells.iter().position(|cell| cell.id == route.cell)
                .ok_or("fiber heat lost its original receiving cell")?;
            let cell = &self.profile.cells[index];
            let source = &self.sources[index];
            let destination = ready.performance.bank().fiber_material_heat_destination(
                route.strip, route.layer, route.cell, &source.card, &source.curve,
                "T", budget, gate,
            )?;
            if destination.source.state_coordinate("T")
                .is_none_or(|temperature| temperature.to_bits() != cell.temperatures[1].to_bits())
            {
                return Err("fiber body and its finite thermal cell require one initial source temperature".into());
            }
            // Reuse the original admission, including source identity, whole
            // annular-layer volume/mass, caloric law, temperature domain and
            // complete fixed-coefficient hereditary heat. This runs BEFORE
            // the output file is opened and repeats in every physical trial.
            destination.admit::<std::convert::Infallible>(network)?;
            routes.push(destination);
        }
        poll(gate)?;
        Ok(routes)
    }
}
