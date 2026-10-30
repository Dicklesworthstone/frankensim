//! Explicit material-coordinate surfaces for the existing moving dry-heat owner.
//! Geometry and ordered participant identities are input data. No thermal mass,
//! conductivity, footprint, receiver motion or partition fraction is inferred.

use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    BowDryContactHeat, DryReceiverHeatPosition,
};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::dry::{
    DryContactFootprint, DryContactHeat, DryThermalCell,
};
use fs_material::thermal::source::ThermalMaterialPoint;

#[derive(Clone, Debug)]
struct Strip {
    cell: CellId,
    interval: [f64; 2],
    width: f64,
    half_depth: f64,
}

#[derive(Clone, Debug)]
pub(super) struct Route {
    port: usize,
    string_origin: f64,
    receiver: DryReceiverHeatPosition,
    footprint: DryContactFootprint,
    surfaces: [Vec<Strip>; 2],
}

pub(super) fn parse(rows: &mut Rows<'_>, cells: &[Cell]) -> Result<Vec<Route>> {
    let count = rows.count("dry-contacts", 64)?;
    let mut routes: Vec<Route> = Vec::with_capacity(count);
    for _ in 0..count {
        let mut row = rows.row("dry-contact")?;
        let port = row.parse()?;
        let string_origin = row.scalar()?;
        let receiver = match row.word()? {
            "fixed" => DryReceiverHeatPosition::Fixed,
            "fiber" => DryReceiverHeatPosition::Fiber { coordinate_origin_m: row.scalar()? },
            _ => return Err("dry receiver must explicitly be fixed or the original fiber".into()),
        };
        row.end()?;
        if routes.iter().any(|route| route.port == port) {
            return Err("each original dry contact has exactly one thermal route".into());
        }
        let mut row = rows.row("footprint")?;
        let a = [row.scalar()?, row.scalar()?];
        let wa = row.positive()?;
        let b = [row.scalar()?, row.scalar()?];
        let wb = row.positive()?;
        row.end()?;
        if a[0] >= a[1] || b[0] >= b[1] {
            return Err("dry footprints require positive material-coordinate lengths".into());
        }
        let mut surfaces: [Vec<Strip>; 2] = [Vec::new(), Vec::new()];
        let mut occupied = std::collections::BTreeSet::new();
        for surface in &mut surfaces {
            let count = rows.count("surface", MAX_CELLS)?;
            if count == 0 {
                return Err("each dry participant needs a finite material surface".into());
            }
            for _ in 0..count {
                let mut row = rows.row("strip")?;
                let strip = Strip {
                    cell: CellId(row.parse()?),
                    interval: [row.scalar()?, row.scalar()?],
                    width: row.positive()?,
                    half_depth: row.positive()?,
                };
                row.end()?;
                if !cells.iter().any(|cell| cell.id == strip.cell)
                    || !occupied.insert(strip.cell)
                    || strip.interval[0] >= strip.interval[1]
                    || surface.last().is_some_and(|previous| previous.interval[1] != strip.interval[0])
                {
                    return Err("dry strips require existing distinct cells and contiguous positive intervals".into());
                }
                surface.push(strip);
            }
        }
        routes.push(Route {
            port, string_origin, receiver, surfaces,
            footprint: DryContactFootprint { intervals_m: [a, b], widths_m: [wa, wb] },
        });
    }
    Ok(routes)
}

fn centre_matches(interval: [f64; 2], origin: f64, position: f64) -> bool {
    let expected = origin + position;
    let actual = 0.5 * interval[0] + 0.5 * interval[1];
    let scale = origin.abs().max(position.abs()).max(interval[0].abs()).max(interval[1].abs());
    expected.is_finite() && (actual - expected).abs() <= 64.0 * f64::EPSILON * scale
}

impl Prepared {
    pub(super) fn dry_routes(&self, ready: &Ready) -> Result<Vec<BowDryContactHeat<'_>>> {
        let mut result = Vec::with_capacity(self.profile.dry.len());
        for route in &self.profile.dry {
            let source = ready.performance.bank().source();
            let bow = source.bows().get(route.port).ok_or("dry heat names an absent bow port")?;
            let track_index = if ready.spec.bow_ribbon.is_some() { 0 } else { route.port };
            let track = ready.spec.schedule.tracks().get(track_index)
                .ok_or("dry heat names an absent original bow track")?;
            let participants = ready.spec.friction.surface_materials(&track.id)?;
            let string_position = bow.location.station_fraction
                * ready.performance.bank().cards()[bow.location.string].length_m;
            if !centre_matches(route.footprint.intervals_m[0], route.string_origin, string_position) {
                return Err("dry string footprint must begin at the actual original bow position".into());
            }
            match (route.receiver, source.fiber_ribbon()) {
                (DryReceiverHeatPosition::Fixed, None) => {}
                (DryReceiverHeatPosition::Fiber { coordinate_origin_m }, Some(fiber)) => {
                    let position = fiber.station_fraction() * fiber.model().geometry().length_m;
                    if !centre_matches(route.footprint.intervals_m[1], coordinate_origin_m, position) {
                        return Err("dry fiber footprint must begin at the actual original hair crossing".into());
                    }
                    let pack = &ready.spec.bow_fibers.as_ref()
                        .ok_or("dry fiber receiver lost its original material description")?.pack;
                    if ready.catalog.get(pack).ok_or("dry fiber material is absent")?.card.id()
                        != &participants[1]
                    {
                        return Err("dry receiver participant must be the actual fiber material".into());
                    }
                }
                _ => return Err("dry receiver motion must match the actual fixed or fiber surface".into()),
            }
            let mut surfaces: [Vec<DryThermalCell<'_>>; 2] = [Vec::new(), Vec::new()];
            for (side, strips) in route.surfaces.iter().enumerate() {
                for strip in strips {
                    let index = self.profile.cells.iter().position(|cell| cell.id == strip.cell)
                        .ok_or("dry heat cell is absent")?;
                    let source = &self.sources[index];
                    if source.card.id() != &participants[side] {
                        return Err("dry thermal cells must retain the original ordered interface material identities".into());
                    }
                    surfaces[side].push(DryThermalCell {
                        cell: strip.cell,
                        material: ThermalMaterialPoint {
                            material: &source.card,
                            context: &source.query,
                            temperature_axis: "T",
                            liquid_fraction_axis: None,
                        },
                        caloric: &source.curve,
                        transport: &source.state,
                        conductivity_property: "thermal_conductivity",
                        density_property: "density",
                        interval_m: strip.interval,
                        width_m: strip.width,
                        half_depth_m: strip.half_depth,
                    });
                }
            }
            result.push(BowDryContactHeat {
                heat: DryContactHeat { port: route.port, surfaces, footprint: route.footprint },
                string_coordinate_origin_m: route.string_origin,
                receiver: route.receiver,
            });
        }
        Ok(result)
    }
}
