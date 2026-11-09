//! Finite, lumped thermal bodies on the original immersed solid mask.
//!
//! A body has one bulk temperature and heat capacity, however many voxels or
//! fluid links form its surface. Its fixed per-link bulk-to-surface resistance
//! uses the same units and reflected-pair law as `FiniteThermalWall`. Summing
//! that law over the body's actual post-collision links gives
//! `(C + sum K) (T_new - T_old) = Q + sum K (T_trace - T_old)`,
//! where `K = 1 / (R + 1/(2w))` and `T_trace = g_out/w`.
//! One scalar backward-Euler solve advances each body; every link then uses
//! the SAME endpoint. Capacity and the external source are counted once.
//!
//! Each node has an explicitly isothermal interior. The optional reciprocal
//! conduction graph resolves temperature differences between nodes, including
//! fully buried solid nodes without fluid interfaces. The independent-body
//! path remains appropriate when internal gradients can be neglected. Moving
//! geometry is not represented. Independent bodies can bind an original
//! material enthalpy chart through `caloric`, retaining sensible and latent
//! storage in the same native transaction. There is no
//! empirical convection coefficient: surrounding fluid heat transport remains
//! the original D3Q7 field driven by native flow.

/// Reciprocal solid conduction on the original immersed thermal nodes.
pub mod conduction;
/// Source-owned sensible and latent storage on independent immersed bodies.
pub mod caloric;
/// Cartesian voxel geometry and material preparation for immersed conduction.
pub mod voxel;

use super::{
    D3Q7, FiniteThermalWallExchange, Sum, ThermalOpenError3, ThermalOpenFace3, ThermalOpenField3,
    ThermalOpenLimits3, WallReflection, check_boundary_temperature, check_wall_response, invalid,
};
use crate::thermal::boundary::FiniteThermalWall;

/// Authored grouping of original stationary solid voxels into one thermal body.
#[derive(Debug, Clone, Copy)]
pub struct ImmersedThermalBody3<'a> {
    /// Distinct x-fast native cell indices. All must already be solid; a cell
    /// cannot belong to two bodies. Interior voxels may have no fluid links.
    pub solid_cells: &'a [usize],
    /// Capacity and source are for the WHOLE body, not per voxel or link.
    /// Resistance is the common bulk-to-surface resistance of ONE exposed
    /// normal link, excluding the fluid half-cell resistance.
    pub solid: FiniteThermalWall,
}

/// Accepted body state and its independently measured last heat account.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImmersedThermalBodyState3 {
    /// Original held source and common resistance declaration. Heat capacity
    /// owns storage only when `caloric` is None; a bound enthalpy chart replaces
    /// C*dT, and this original capacity is then inactive declaration data.
    /// The resolved interface law uses `surface_resistance` below, which may
    /// additionally include a direction-dependent solid half-cell path.
    pub solid: FiniteThermalWall,
    /// Original material state and fixed mass when source enthalpy owns storage.
    /// `solid.temperature` is its derived native-temperature coordinate.
    pub caloric: Option<caloric::ImmersedCaloricState3>,
    /// Total bulk-to-surface resistance for x, y and z normal links, excluding
    /// the fluid half cell. Uniform body declarations repeat `solid.resistance`.
    pub surface_resistance: [f64; 3],
    pub interface_links: usize,
    /// `maximum_link_defect` holds this shared body's absolute storage defect;
    /// storage cannot be assigned independently to its individual links.
    pub last_exchange: Option<FiniteThermalWallExchange>,
}

struct Body {
    accepted: ImmersedThermalBodyState3,
    // Sum of the actual exposed-link conductances, including their fluid
    // half cells. No averaging of the solid conductivity enters this value.
    conductance: f64,
    trial_temperature: f64,
    reference_population: usize,
    trace_reference: f64,
    // Weighted trace differences from one ORIGINAL exposed trace. A common
    // fluid temperature is represented exactly, independent of rounded weights.
    mean_trace: Sum,
    heat: Sum,
    conducted: Sum,
    trial_exchange: FiniteThermalWallExchange,
    caloric: Option<caloric::CaloricBody3>,
}

struct Link {
    cell: usize,
    incoming: usize,
    body: usize,
    conductance: f64,
    trace_weight: f64,
}

pub(super) struct ImmersedBodies3 {
    bodies: Vec<Body>,
    // Sorted in original streaming order: a linear cursor replaces a dense
    // per-cell map or a geometry search in the hot boundary loop.
    links: Vec<Link>,
    conduction: Option<conduction::ImmersedConduction3>,
}

fn visit_links(field: &ThermalOpenField3, solid_cell: usize, mut visit: impl FnMut(usize, usize)) {
    let strides = [1, field.dims[0], field.dims[0] * field.dims[1]];
    for axis in 0..3 {
        let coordinate = solid_cell / strides[axis] % field.dims[axis];
        let (positive, negative) = D3Q7.direction_pair(axis);
        for (incoming, upper) in [(positive, true), (negative, false)] {
            let outside = if upper {
                coordinate + 1 == field.dims[axis]
            } else {
                coordinate == 0
            };
            if outside
                && !matches!(
                    field.boundaries.faces[2 * axis + usize::from(upper)],
                    ThermalOpenFace3::Periodic
                )
            {
                continue;
            }
            let fluid = if outside {
                if upper {
                    solid_cell - (field.dims[axis] - 1) * strides[axis]
                } else {
                    solid_cell + (field.dims[axis] - 1) * strides[axis]
                }
            } else if upper {
                solid_cell + strides[axis]
            } else {
                solid_cell - strides[axis]
            };
            if !field.solid[fluid] {
                visit(fluid, incoming);
            }
        }
    }
}

impl ThermalOpenField3 {
    /// Attach shared finite storage to selected original immersed solid cells.
    ///
    /// The explicit additional payload byte budget bounds retained bodies, boundary
    /// links and preparation-only cell-identity scratch. No allocation or
    /// geometry search occurs during a step. Replacement is transactional and
    /// briefly retains the previous plan; fluid fields and time are unchanged.
    /// This explicitly installs the supplied body temperatures and inventories.
    /// Unassigned solid cells retain `boundaries().internal_wall`.
    pub fn set_immersed_thermal_bodies(
        &mut self,
        declarations: &[ImmersedThermalBody3<'_>],
        maximum_storage_bytes: usize,
    ) -> Result<(), ThermalOpenError3> {
        let prepared = self.prepare_immersed_bodies(declarations, maximum_storage_bytes, false)?;
        self.immersed = Some(prepared);
        Ok(())
    }

    fn prepare_immersed_bodies(
        &self,
        declarations: &[ImmersedThermalBody3<'_>],
        maximum_storage_bytes: usize,
        allow_interior: bool,
    ) -> Result<ImmersedBodies3, ThermalOpenError3> {
        self.prepare_immersed_bodies_with_surface_resistance(
            declarations,
            maximum_storage_bytes,
            allow_interior,
            |body| [declarations[body].solid.resistance; 3],
        )
    }

    fn prepare_immersed_bodies_with_surface_resistance(
        &self,
        declarations: &[ImmersedThermalBody3<'_>],
        maximum_storage_bytes: usize,
        allow_interior: bool,
        surface_resistance: impl Fn(usize) -> [f64; 3],
    ) -> Result<ImmersedBodies3, ThermalOpenError3> {
        if declarations.is_empty() || declarations.len() > self.solid.len() {
            return Err(invalid(
                "immersed bodies require a nonempty bounded declaration",
            ));
        }
        let mut cells = 0usize;
        let mut links = 0usize;
        for (body, declaration) in declarations.iter().enumerate() {
            declaration
                .solid
                .validate()
                .map_err(ThermalOpenError3::Wall)?;
            check_boundary_temperature(6, declaration.solid.temperature, self.limits)?;
            if surface_resistance(body)
                .iter()
                .any(|r| !r.is_finite() || *r < 0.0)
            {
                return Err(invalid(
                    "immersed surface resistance must be finite and nonnegative",
                ));
            }
            cells = cells
                .checked_add(declaration.solid_cells.len())
                .ok_or_else(|| invalid("immersed body cell count overflow"))?;
            if declaration.solid_cells.is_empty() || cells > self.solid.len() {
                return Err(invalid("immersed body cells must be nonempty and disjoint"));
            }
            let mut exposed = 0usize;
            for &cell in declaration.solid_cells {
                if self.solid.get(cell) != Some(&true) {
                    return Err(invalid(
                        "immersed body cell must belong to the native solid mask",
                    ));
                }
                visit_links(self, cell, |_, _| exposed += 1);
            }
            if exposed == 0 && !allow_interior {
                return Err(invalid(
                    "each immersed body needs an exposed fluid interface",
                ));
            }
            links = links
                .checked_add(exposed)
                .ok_or_else(|| invalid("immersed body interface count overflow"))?;
        }
        let bytes = declarations
            .len()
            .checked_mul(core::mem::size_of::<Body>())
            .and_then(|n| {
                links
                    .checked_mul(core::mem::size_of::<Link>())
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| {
                cells
                    .checked_mul(core::mem::size_of::<usize>())
                    .and_then(|m| n.checked_add(m))
            })
            .ok_or_else(|| invalid("immersed body storage overflow"))?;
        if bytes > maximum_storage_bytes {
            return Err(invalid("immersed body storage budget exceeded"));
        }
        let mut claimed = Vec::new();
        claimed
            .try_reserve_exact(cells)
            .map_err(|_| ThermalOpenError3::Allocation)?;
        for declaration in declarations {
            claimed.extend_from_slice(declaration.solid_cells);
        }
        claimed.sort_unstable();
        if claimed.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid(
                "an immersed solid cell cannot have two thermal owners",
            ));
        }
        let mut prepared = ImmersedBodies3 {
            bodies: Vec::new(),
            links: Vec::new(),
            conduction: None,
        };
        prepared
            .bodies
            .try_reserve_exact(declarations.len())
            .map_err(|_| ThermalOpenError3::Allocation)?;
        prepared
            .links
            .try_reserve_exact(links)
            .map_err(|_| ThermalOpenError3::Allocation)?;
        for (body, declaration) in declarations.iter().enumerate() {
            let start = prepared.links.len();
            let resistance = surface_resistance(body);
            let normal_conductance = resistance.map(|r| 1.0 / (r + 0.5 / D3Q7.weights[1]));
            for &cell in declaration.solid_cells {
                visit_links(self, cell, |cell, incoming| {
                    prepared.links.push(Link {
                        cell,
                        incoming,
                        body,
                        conductance: normal_conductance[(incoming - 1) / 2],
                        trace_weight: 0.0,
                    });
                });
            }
            let exposed = &mut prepared.links[start..];
            let conductance = if normal_conductance
                .iter()
                .all(|k| *k == normal_conductance[0])
            {
                normal_conductance[0] * exposed.len() as f64
            } else {
                let mut total = Sum::default();
                for link in exposed.iter() {
                    total.add(link.conductance);
                }
                total.total()
            };
            if !conductance.is_finite() || (!exposed.is_empty() && conductance <= 0.0) {
                return Err(invalid("unrepresentable immersed interface conductance"));
            }
            for link in exposed.iter_mut() {
                link.trace_weight = link.conductance / conductance;
            }
            let reference_population = exposed
                .iter()
                .min_by_key(|link| (link.cell, link.incoming))
                .map_or(0, |link| 7 * link.cell + opposite(link.incoming));
            prepared.bodies.push(Body {
                accepted: ImmersedThermalBodyState3 {
                    solid: declaration.solid,
                    caloric: None,
                    surface_resistance: resistance,
                    interface_links: prepared.links.len() - start,
                    last_exchange: None,
                },
                conductance,
                trial_temperature: declaration.solid.temperature,
                reference_population,
                trace_reference: 0.0,
                mean_trace: Sum::default(),
                heat: Sum::default(),
                conducted: Sum::default(),
                trial_exchange: FiniteThermalWallExchange::default(),
                caloric: None,
            });
        }
        prepared
            .links
            .sort_unstable_by_key(|link| (link.cell, link.incoming));
        Ok(prepared)
    }

    /// Accepted state in declaration order; refused steps never alter it.
    pub fn immersed_body(&self, body: usize) -> Option<ImmersedThermalBodyState3> {
        self.immersed
            .as_ref()?
            .bodies
            .get(body)
            .map(|body| body.accepted)
    }

    /// Change the held heat source of one whole body without advancing time.
    /// Negative heat represents an explicitly authored sink. The previous
    /// accepted receipt remains an account of its original interval.
    pub fn set_immersed_body_heat(
        &mut self,
        body: usize,
        heat_per_step: f64,
    ) -> Result<(), ThermalOpenError3> {
        if !heat_per_step.is_finite() {
            return Err(invalid("immersed body source must be finite"));
        }
        let body = self
            .immersed
            .as_mut()
            .and_then(|plan| plan.bodies.get_mut(body))
            .ok_or_else(|| invalid("unknown immersed thermal body"))?;
        body.accepted.solid.heat_per_step = heat_per_step;
        Ok(())
    }
}

fn opposite(q: usize) -> usize {
    if q % 2 == 0 { q - 1 } else { q + 1 }
}

impl Body {
    fn mean_fluid_trace(&self) -> f64 {
        self.trace_reference + self.mean_trace.total()
    }

    fn reflect(&self, outgoing: f64, conductance: f64) -> WallReflection {
        let trace = outgoing / D3Q7.weights[1];
        let difference = self.trial_temperature - trace;
        let heat = if difference.is_finite() {
            conductance * difference
        } else {
            conductance * self.trial_temperature - conductance * trace
        };
        let incoming = outgoing + heat;
        WallReflection {
            incoming,
            temperature: (0.5 * incoming + 0.5 * outgoing) / D3Q7.weights[1],
            // Shared storage and source are accounted ONCE after all links.
            solid: None,
        }
    }
}

impl ImmersedBodies3 {
    pub(super) fn prepare<E: From<ThermalOpenError3>>(
        &mut self,
        post: &[f64],
        limits: ThermalOpenLimits3,
        mut poll: impl FnMut() -> Result<(), E>,
    ) -> Result<(), E> {
        for (i, body) in self.bodies.iter_mut().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            body.mean_trace = Sum::default();
            body.trace_reference = if body.accepted.interface_links == 0 {
                body.accepted.solid.temperature
            } else {
                post[body.reference_population] / D3Q7.weights[1]
            };
            body.heat = Sum::default();
            body.conducted = Sum::default();
        }
        for (i, link) in self.links.iter().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            let body = &mut self.bodies[link.body];
            let trace = post[7 * link.cell + opposite(link.incoming)] / D3Q7.weights[1];
            if !trace.is_finite() {
                return Err(invalid("nonfinite immersed interface population").into());
            }
            // Center before weighting to preserve an exactly uniform trace.
            // Scale first only when the unweighted difference overflows.
            let difference = trace - body.trace_reference;
            body.mean_trace.add(if difference.is_finite() {
                link.trace_weight * difference
            } else {
                link.trace_weight * trace - link.trace_weight * body.trace_reference
            });
        }
        if let Some(conduction) = &mut self.conduction {
            conduction.solve(&mut self.bodies, limits, &mut poll)?;
        } else {
            for (i, body) in self.bodies.iter_mut().enumerate() {
                if i % 64 == 0 {
                    poll()?;
                }
                let solid = body.accepted.solid;
                let conductance = body.conductance;
                if !conductance.is_finite() || conductance <= 0.0 {
                    return Err(invalid("unrepresentable immersed interface conductance").into());
                }
                if let Some(caloric) = &mut body.caloric {
                    let old = body.accepted.caloric.expect("bound caloric state");
                    body.trial_temperature = caloric
                        .solve(
                            i,
                            old,
                            solid.heat_per_step,
                            conductance,
                            body.trace_reference + body.mean_trace.total(),
                            limits,
                            &mut poll,
                        )?;
                    continue;
                }
                // Normalize by the larger coefficient. Never form C + sum K;
                // the bounded ratio also preserves a tiny capacity or exchange
                // coefficient without constructing an overflowing reciprocal.
                let (old_weight, fluid_weight, scale, dominant_weight) =
                    if solid.heat_capacity <= conductance {
                        let ratio = solid.heat_capacity / conductance;
                        let fluid_weight = 1.0 / (1.0 + ratio);
                        (
                            ratio * fluid_weight,
                            fluid_weight,
                            conductance,
                            fluid_weight,
                        )
                    } else {
                        let ratio = conductance / solid.heat_capacity;
                        let old_weight = 1.0 / (1.0 + ratio);
                        (
                            old_weight,
                            ratio * old_weight,
                            solid.heat_capacity,
                            old_weight,
                        )
                    };
                let source_ratio = solid.heat_per_step / scale;
                let source_temperature = if source_ratio.is_finite() {
                    source_ratio * dominant_weight
                } else {
                    // Scaling first recovers finite quotients where Q/scale
                    // overflows although Q/(C+sum K) remains representable.
                    (solid.heat_per_step * dominant_weight) / scale
                };
                let mean = body.mean_fluid_trace();
                let difference = mean - solid.temperature;
                let retained = if difference.is_finite() {
                    if fluid_weight <= old_weight {
                        fluid_weight.mul_add(difference, solid.temperature)
                    } else {
                        (-old_weight).mul_add(difference, mean)
                    }
                } else {
                    old_weight.mul_add(solid.temperature, fluid_weight * mean)
                };
                body.trial_temperature = retained + source_temperature;
                check_boundary_temperature(6, body.trial_temperature, limits)?;
            }
        }
        for (i, link) in self.links.iter().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            let body = &mut self.bodies[link.body];
            let outgoing = post[7 * link.cell + opposite(link.incoming)];
            let response = body.reflect(outgoing, link.conductance);
            check_wall_response(6, &response, limits)?;
            body.heat.add(response.incoming - outgoing);
        }
        let tolerance = self
            .conduction
            .as_ref()
            .map_or(limits.balance_tolerance, |graph| {
                limits.balance_tolerance.min(graph.balance_tolerance())
            });
        let mut total_defect = Sum::default();
        let mut total_conducted = Sum::default();
        for (i, body) in self.bodies.iter_mut().enumerate() {
            if i % 64 == 0 {
                poll()?;
            }
            let solid = body.accepted.solid;
            let stored_change = body.caloric.as_ref().map_or_else(
                || solid.heat_capacity * (body.trial_temperature - solid.temperature),
                |caloric| caloric.stored_change(body.accepted.caloric.expect("bound caloric state")),
            );
            let mut exchange = FiniteThermalWallExchange {
                heat_in: body.heat.total(),
                stored_change,
                source: solid.heat_per_step,
                conducted_in: body.conducted.total(),
                ..FiniteThermalWallExchange::default()
            };
            let defect = exchange.defect();
            let body_tolerance = body.caloric.as_ref().map_or(tolerance, |caloric| {
                tolerance.min(caloric.balance_tolerance())
            });
            if !defect.is_finite() || defect.abs() > body_tolerance {
                return Err(ThermalOpenError3::Balance { defect }.into());
            }
            total_defect.add(defect);
            total_conducted.add(exchange.conducted_in);
            exchange.maximum_link_defect = defect.abs();
            body.trial_exchange = exchange;
        }
        if self.conduction.is_some() {
            for defect in [total_defect.total(), total_conducted.total()] {
                if !defect.is_finite() || defect.abs() > tolerance {
                    return Err(ThermalOpenError3::Balance { defect }.into());
                }
            }
        }
        Ok(())
    }

    pub(super) fn reflection(
        &self,
        cell: usize,
        incoming: usize,
        outgoing: f64,
        cursor: &mut usize,
    ) -> Option<WallReflection> {
        let link = self.links.get(*cursor)?;
        if (link.cell, link.incoming) != (cell, incoming) {
            return None;
        }
        *cursor += 1;
        Some(self.bodies[link.body].reflect(outgoing, link.conductance))
    }
    pub(super) fn link_count(&self) -> usize {
        self.links.len()
    }
    pub(super) fn trial_exchanges(&self) -> impl Iterator<Item = FiniteThermalWallExchange> + '_ {
        self.bodies.iter().map(|body| body.trial_exchange)
    }
    pub(super) fn publish(&mut self) {
        for body in &mut self.bodies {
            if let Some(caloric) = &body.caloric {
                body.accepted.caloric = Some(caloric.trial_state());
            }
            body.accepted.solid.temperature = body.trial_temperature;
            body.accepted.last_exchange = Some(body.trial_exchange);
        }
    }
}

#[cfg(test)]
mod tests;
