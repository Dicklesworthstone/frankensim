//! Shared resolved column geometry for capillary traction and heat patches.
use super::*;

/// One physical surface patch in lattice length/area units.
#[derive(Clone, Copy, Debug)]
pub(crate) struct HeightGraphPatch3 {
    pub column: [usize; 2],
    pub centroid: [f64; 3],
    pub outward_normal: [f64; 3],
    pub area: f64,
    pub volume: f64,
    pub liquid: [usize; 3],
    pub gas: [usize; 3],
}

impl FreeSurface3 {
    /// Reconstruct the same complete graph used by stress and thermal exchange.
    /// The centroid is the height-column quadrature point, not a cut-cell
    /// liquid centroid. Thermal normal distances remain explicit model data.
    pub(crate) fn height_graph_geometry(
        &self,
        graph: HeightGraphSurface3,
        gate: &CancelGate,
    ) -> Result<Vec<HeightGraphPatch3>, FreeSurfaceError3> {
        if graph.axis >= 3
            || ![-1, 1].contains(&graph.outward_sign)
            || !graph.maximum_slope.is_finite()
            || !(0.0..=1.0).contains(&graph.maximum_slope)
            || !graph.minimum_fill.is_finite()
            || graph.minimum_fill <= 0.0
            || graph.minimum_fill > 1.0
            || !self.solid.is_empty()
            || self.planar_contact.is_some()
            || self.planar_contact_field.is_some()
            || self.contact != ContactModel3::Neutral
            || self.carry != 0.0
        {
            return Err(invalid(
                "height graph needs resolved Cartesian columns, neutral walls and bounded slopes",
            ));
        }
        let (nx, ny, nz) = self.grid.tile_dims();
        let size = [nx * TILE, ny * TILE, nz * TILE];
        let tangent = match graph.axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let shape = [size[tangent[0]], size[tangent[1]]];
        let count = shape[0]
            .checked_mul(shape[1])
            .ok_or(invalid("height-graph column count overflows"))?;
        if shape.iter().any(|&n| n < 3) || count > self.tiles.len().saturating_mul(TILE_CELLS) {
            return Err(invalid("height graph has unresolved tank columns"));
        }
        const TOLERANCE: f64 = 64.0 * f64::EPSILON;
        let mut patches = Vec::with_capacity(count);
        for b in 0..shape[1] {
            for a in 0..shape[0] {
                check_cancel(gate)?;
                let mut height = 0.0;
                let mut previous = 1.0;
                let mut volume = 0.0;
                let mut liquid = None;
                let mut gas = None;
                for depth in 0..size[graph.axis] {
                    if depth % TILE == 0 {
                        check_cancel(gate)?;
                    }
                    let mut p = [0; 3];
                    p[tangent[0]] = a;
                    p[tangent[1]] = b;
                    p[graph.axis] = if graph.outward_sign > 0 {
                        depth
                    } else {
                        size[graph.axis] - 1 - depth
                    };
                    let (fill, interface) =
                        match self
                            .grid
                            .resolve_source(p[0] as i64, p[1] as i64, p[2] as i64)
                        {
                            Some((slot, lane)) => {
                                let key = self.grid.active_keys()[slot];
                                match self.tiles[&key].cells[lane] {
                                    Cell3::Gas => (0.0, false),
                                    Cell3::Fluid => (1.0, false),
                                    Cell3::Interface => {
                                        let rho =
                                            self.grid.populations(slot, lane).iter().sum::<f64>();
                                        if !rho.is_finite() || rho <= 0.0 {
                                            return Err(invalid("invalid density in height graph"));
                                        }
                                        (self.tiles[&key].mass[lane] / rho, true)
                                    }
                                }
                            }
                            None => (0.0, false),
                        };
                    if !fill.is_finite()
                        || fill < 0.0
                        || fill > 1.0 + TOLERANCE
                        || fill > previous + TOLERANCE
                        || (depth == 0 && fill < 1.0 - TOLERANCE)
                        || (depth + 1 == size[graph.axis] && fill > 0.0)
                    {
                        return Err(invalid(
                            "height graph requires complete monotone liquid-to-gas columns",
                        ));
                    }
                    let fill = fill.min(1.0);
                    height += fill;
                    previous = fill;
                    if interface && fill > 0.0 {
                        volume += fill;
                        liquid = Some(p);
                        gas = None;
                    } else if fill == 0.0 && liquid.is_some() && gas.is_none() {
                        gas = Some(p);
                    }
                }
                if height <= 1.0
                    || height >= (size[graph.axis] - 1) as f64
                    || volume < graph.minimum_fill
                {
                    return Err(invalid(
                        "height graph needs liquid depth, gas clearance and a resolved interface volume",
                    ));
                }
                let liquid = liquid.ok_or(invalid("height graph has no liquid interface"))?;
                let gas = gas.ok_or(invalid("height graph has no exterior gas cell"))?;
                let mut centroid = [0.0; 3];
                centroid[tangent[0]] = a as f64 + 0.5;
                centroid[tangent[1]] = b as f64 + 0.5;
                centroid[graph.axis] = if graph.outward_sign > 0 {
                    height
                } else {
                    size[graph.axis] as f64 - height
                };
                patches.push(HeightGraphPatch3 {
                    column: [a, b],
                    centroid,
                    outward_normal: [0.0; 3],
                    area: 0.0,
                    volume,
                    liquid,
                    gas,
                });
            }
        }
        for b in 0..shape[1] {
            for a in 0..shape[0] {
                check_cancel(gate)?;
                let height = |p: [usize; 2]| patches[p[0] + shape[0] * p[1]].centroid[graph.axis];
                let hu = graph_derivative([a, b], 0, shape, height);
                let hv = graph_derivative([a, b], 1, shape, height);
                let square = hu * hu + hv * hv;
                if !square.is_finite() || square.sqrt() > graph.maximum_slope + TOLERANCE {
                    return Err(invalid(
                        "free surface exceeded its height-graph slope envelope",
                    ));
                }
                let patch = &mut patches[a + shape[0] * b];
                patch.area = (1.0 + square).sqrt();
                let scale = f64::from(graph.outward_sign) / patch.area;
                patch.outward_normal[tangent[0]] = -hu * scale;
                patch.outward_normal[tangent[1]] = -hv * scale;
                patch.outward_normal[graph.axis] = scale;
            }
        }
        Ok(patches)
    }
}
