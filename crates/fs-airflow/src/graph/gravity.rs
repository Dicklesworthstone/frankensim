//! Gravity on the original pressure graph, with one conserved reference volume.
//!
//! Static pressure obeys `p_from - p_to + rho g·dx = loss(Q)`. Solve in
//! `p - rho0 g·(x-x0)` so the uniform-density hydrostatic state is balanced
//! without iteration or subtraction of large heads inside each edge solve.
//! Only the integrated `(rho-rho0) g·dx` remains in the numerical edge law. Different edge
//! densities are a declared Boussinesq approximation: mass and heat transport
//! retain rho0, and the frozen density samples affect gravity alone.
//!
//! No density/temperature evolution, compressible hydrostatics, pressure work,
//! or self-gravity is implied. Each material is queried once at binding; all
//! graph iterations use cached coefficients and the existing loss solver.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_qty::Pressure;

use super::{FixedPressure, GraphError, LossGraph, Sum, checkpoint, finite, periodic_checkpoint};

/// Hard allocation bound for all explicit gravity sections in one graph.
pub const MAX_GRAVITY_SECTIONS: usize = 65_536;

/// One oriented straight section of an original branch, with its frozen material.
/// A bent branch retains every section instead of inventing an effective density.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphGravitySection {
    /// Section start in the node positions' Cartesian frame [m].
    pub start_m: [f64; 3],
    /// Section end in the same frame [m].
    pub end_m: [f64; 3],
    /// Original local material snapshot; density alone supplies gravity.
    pub transport: FluidTransportState,
}

#[derive(Debug, Clone, PartialEq)]
enum GravitySamples {
    Branches(Vec<FluidTransportState>),
    Sections {
        branches: Vec<Vec<GraphGravitySection>>,
        maximum_sections: usize,
    },
}

/// Frozen gravity geometry and the original material-transport snapshots.
/// A numerical snapshot does not create source authority beyond its material.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphGravity {
    acceleration_m_s2: [f64; 3],
    positions_m: Vec<[f64; 3]>,
    reference: FluidTransportState,
    samples: GravitySamples,
    maximum_relative_density_departure: f64,
    reference_pressures_pa: Vec<f64>,
    departure_heads_pa: Vec<f64>,
    absolute_departure_heads_pa: Vec<f64>,
    pressure_heads_pa: Vec<f64>,
}

impl GraphGravity {
    /// Signed acceleration in the same orthonormal frame as node positions.
    pub const fn acceleration_m_s2(&self) -> [f64; 3] {
        self.acceleration_m_s2
    }
    /// Original finite node positions, in metres.
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }
    /// One density and Cp for conserved reference mass and sensible heat.
    pub const fn reference_transport(&self) -> &FluidTransportState {
        &self.reference
    }
    /// Original single samples, when every branch has one uniform material.
    /// Distributed paths have no single effective density and return None.
    pub fn branch_materials(&self) -> Option<&[FluidTransportState]> {
        match &self.samples {
            GravitySamples::Branches(materials) => Some(materials),
            GravitySamples::Sections { .. } => None,
        }
    }
    /// Every original section and its local material, in branch/path order.
    pub fn branch_sections(&self) -> Option<&[Vec<GraphGravitySection>]> {
        match &self.samples {
            GravitySamples::Branches(_) => None,
            GravitySamples::Sections { branches, .. } => Some(branches),
        }
    }
    /// Explicit total section allocation bound, absent for uniform branch samples.
    pub const fn maximum_sections(&self) -> Option<usize> {
        match &self.samples {
            GravitySamples::Branches(_) => None,
            GravitySamples::Sections {
                maximum_sections, ..
            } => Some(*maximum_sections),
        }
    }
    /// Same gravity frame and original path partition, independent of local density.
    /// A retained momentum owner must not change its geometry during source refresh.
    pub fn same_geometry_as(&self, other: &Self) -> bool {
        self.acceleration_m_s2 == other.acceleration_m_s2
            && self.positions_m == other.positions_m
            && match (&self.samples, &other.samples) {
                (GravitySamples::Branches(a), GravitySamples::Branches(b)) => a.len() == b.len(),
                (
                    GravitySamples::Sections {
                        branches: a,
                        maximum_sections: a_max,
                    },
                    GravitySamples::Sections {
                        branches: b,
                        maximum_sections: b_max,
                    },
                ) => {
                    a_max == b_max
                        && a.len() == b.len()
                        && a.iter().zip(b).all(|(a, b)| {
                            a.len() == b.len()
                                && a.iter()
                                    .zip(b)
                                    .all(|(a, b)| a.start_m == b.start_m && a.end_m == b.end_m)
                        })
                }
                _ => false,
            }
    }
    /// Admitted `abs(rho-rho0)/rho0`; zero is the uniform-density specialization.
    pub const fn maximum_relative_density_departure(&self) -> f64 {
        self.maximum_relative_density_departure
    }
    /// Full gravitational pressure supply, signed in each edge's direction [Pa].
    pub fn pressure_heads_pa(&self) -> &[f64] {
        &self.pressure_heads_pa
    }
    /// Signed integral of `(rho-rho0) g·dx` in each original edge direction [Pa].
    pub fn departure_heads_pa(&self) -> &[f64] {
        &self.departure_heads_pa
    }
    /// Sum of absolute local anomalous heads, before any section cancellation [Pa].
    /// Multiplying by the absolute branch flow bounds omitted local buoyancy power.
    pub fn absolute_departure_heads_pa(&self) -> &[f64] {
        &self.absolute_departure_heads_pa
    }
    pub(super) fn reference_pressure(&self, node: usize) -> f64 {
        self.reference_pressures_pa[node]
    }
}

impl LossGraph {
    /// Attach uniform gravity using the SAME material state as flow/heat.
    /// Node positions and acceleration use one caller-declared Cartesian frame;
    /// zero or reversed gravity is valid. Gravity never becomes a loss term.
    ///
    /// This replaces any prior binding and leaves resistance coefficients
    /// unchanged. Authored coefficients must describe this material; material-
    /// bound Darcy coefficients are checked against the reference density.
    pub fn with_gravity<F: FluidTransport + ?Sized>(
        self,
        material: &F,
        node_positions_m: &[[f64; 3]],
        acceleration_m_s2: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<Self, GraphError> {
        checkpoint(cx)?;
        let reference = material.transport_state().map_err(GraphError::Material)?;
        let materials = vec![reference; self.branches.len()];
        self.bind_gravity(
            reference,
            GravitySamples::Branches(materials),
            node_positions_m,
            acceleration_m_s2,
            0.0,
            cx,
        )
    }

    /// Attach bounded material-density departures for natural circulation.
    ///
    /// The caller supplies one frozen material per edge and an explicit bound
    /// in `[0,1)` on `abs(rho-rho0)/rho0`. Reference volume/mass conservation
    /// and existing drag coefficients stay fixed; only gravitational head uses
    /// each actual density. This is not variable-density mass transport or a
    /// thermal-feedback solve. Updating temperatures requires rebinding with
    /// new source states and solving again. No hidden expansion coefficient is used.
    #[allow(clippy::too_many_arguments)]
    pub fn with_boussinesq_gravity<F: FluidTransport + ?Sized, G: FluidTransport>(
        self,
        reference: &F,
        branch_materials: &[G],
        node_positions_m: &[[f64; 3]],
        acceleration_m_s2: [f64; 3],
        maximum_relative_density_departure: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, GraphError> {
        checkpoint(cx)?;
        if branch_materials.len() != self.branches.len() {
            return Err(GraphError::InvalidInput(
                "gravity needs one material state per branch",
            ));
        }
        let reference = reference.transport_state().map_err(GraphError::Material)?;
        let mut materials = Vec::with_capacity(branch_materials.len());
        for (edge, material) in branch_materials.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            materials.push(material.transport_state().map_err(GraphError::Material)?);
        }
        self.bind_gravity(
            reference,
            GravitySamples::Branches(materials),
            node_positions_m,
            acceleration_m_s2,
            maximum_relative_density_departure,
            cx,
        )
    }

    /// Attach a bounded line integral of local material density along each branch.
    ///
    /// Paths are ordered from the original from-node to the to-node. Adjacent
    /// section endpoints must agree exactly; every section has finite positive
    /// length. The declared total section budget is retained with the geometry.
    /// The reference material still owns conserved mass, heat and Darcy drag.
    /// All local density departures, including horizontal sections, must obey
    /// the stated Boussinesq allowance. Binding caches both the signed head for
    /// native momentum and the absolute local sum for omitted-work envelopes.
    #[allow(clippy::too_many_arguments)]
    pub fn with_distributed_boussinesq_gravity<F: FluidTransport + ?Sized>(
        self,
        reference: &F,
        branch_sections: &[Vec<GraphGravitySection>],
        node_positions_m: &[[f64; 3]],
        acceleration_m_s2: [f64; 3],
        maximum_relative_density_departure: f64,
        maximum_sections: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, GraphError> {
        checkpoint(cx)?;
        if branch_sections.len() != self.branches.len()
            || maximum_sections == 0
            || maximum_sections > MAX_GRAVITY_SECTIONS
        {
            return Err(GraphError::InvalidInput(
                "distributed gravity needs one path per branch and a bounded section budget",
            ));
        }
        let mut count = 0usize;
        for (edge, sections) in branch_sections.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            if sections.is_empty() || sections.len() > maximum_sections - count {
                return Err(GraphError::InvalidBranch {
                    branch: edge,
                    reason: "gravity path is empty or exceeds the total section budget",
                });
            }
            count += sections.len();
        }
        let reference = reference.transport_state().map_err(GraphError::Material)?;
        self.bind_gravity(
            reference,
            GravitySamples::Sections {
                branches: branch_sections.to_vec(),
                maximum_sections,
            },
            node_positions_m,
            acceleration_m_s2,
            maximum_relative_density_departure,
            cx,
        )
    }

    /// Original gravity binding, or None for the original pressure-only graph.
    pub fn gravity(&self) -> Option<&GraphGravity> {
        self.gravity.as_deref()
    }

    #[allow(clippy::too_many_arguments)]
    fn bind_gravity(
        mut self,
        reference: FluidTransportState,
        samples: GravitySamples,
        positions: &[[f64; 3]],
        acceleration: [f64; 3],
        maximum_departure: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, GraphError> {
        if positions.len() != self.node_count()
            || positions.iter().flatten().any(|v| !v.is_finite())
            || acceleration.iter().any(|v| !v.is_finite())
            || !maximum_departure.is_finite()
            || !(0.0..1.0).contains(&maximum_departure)
        {
            return Err(GraphError::InvalidInput(
                "gravity needs finite acceleration/positions, one position per node, and a density bound in [0,1)",
            ));
        }
        let rho0 = reference.density_kg_m3();
        let mut reference_pressures_pa = Vec::with_capacity(positions.len());
        for (node, &position) in positions.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            reference_pressures_pa.push(head(acceleration, positions[0], position, rho0)?);
        }
        let mut departure_heads_pa = Vec::with_capacity(self.branches.len());
        let mut absolute_departure_heads_pa = Vec::with_capacity(self.branches.len());
        let mut pressure_heads_pa = Vec::with_capacity(self.branches.len());
        for (edge, branch) in self.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            if self.linear_losses[edge]
                .as_ref()
                .and_then(|loss| loss.material())
                .is_some_and(|binding| binding.gas.density != rho0)
            {
                return Err(GraphError::InvalidBranch {
                    branch: edge,
                    reason: "gravity reference density differs from the bound Darcy material",
                });
            }
            let (perturbation, absolute) = match &samples {
                GravitySamples::Branches(materials) => {
                    let departure =
                        density_departure(materials[edge], rho0, maximum_departure, edge)?;
                    let value = head(
                        acceleration,
                        positions[branch.from],
                        positions[branch.to],
                        departure,
                    )?;
                    (value, value.abs())
                }
                GravitySamples::Sections { branches, .. } => {
                    let mut signed = Sum::default();
                    let mut absolute = Sum::default();
                    let mut previous = positions[branch.from];
                    for (index, section) in branches[edge].iter().enumerate() {
                        periodic_checkpoint(index, cx)?;
                        if section.start_m != previous
                            || section
                                .start_m
                                .iter()
                                .chain(&section.end_m)
                                .any(|v| !v.is_finite())
                        {
                            return Err(GraphError::InvalidBranch {
                                branch: edge,
                                reason: "gravity sections must form a finite continuous original path",
                            });
                        }
                        let mut displacement = [0.0; 3];
                        for axis in 0..3 {
                            displacement[axis] = finite(
                                section.end_m[axis] - section.start_m[axis],
                                "gravity section displacement",
                            )?;
                        }
                        let scale = displacement.iter().fold(0.0_f64, |a, d| a.max(d.abs()));
                        let length = if scale == 0.0 {
                            0.0
                        } else {
                            scale
                                * fs_math::det::sqrt(
                                    displacement
                                        .iter()
                                        .map(|d| {
                                            let normalized = d / scale;
                                            normalized * normalized
                                        })
                                        .sum(),
                                )
                        };
                        if !length.is_finite() || length <= 0.0 {
                            return Err(GraphError::InvalidBranch {
                                branch: edge,
                                reason: "gravity sections need finite positive length",
                            });
                        }
                        let departure =
                            density_departure(section.transport, rho0, maximum_departure, edge)?;
                        let value = head(acceleration, section.start_m, section.end_m, departure)?;
                        signed.add(value)?;
                        absolute.add(value.abs())?;
                        previous = section.end_m;
                    }
                    if previous != positions[branch.to] {
                        return Err(GraphError::InvalidBranch {
                            branch: edge,
                            reason: "gravity path must end at its original to-node",
                        });
                    }
                    (signed.value()?, absolute.value()?)
                }
            };
            let full = finite(
                (reference_pressures_pa[branch.to] - reference_pressures_pa[branch.from])
                    + perturbation,
                "gravitational pressure head",
            )?;
            departure_heads_pa.push(perturbation);
            absolute_departure_heads_pa.push(absolute);
            pressure_heads_pa.push(full);
        }
        checkpoint(cx)?;
        self.gravity = Some(Arc::new(GraphGravity {
            acceleration_m_s2: acceleration,
            positions_m: positions.to_vec(),
            reference,
            samples,
            maximum_relative_density_departure: maximum_departure,
            reference_pressures_pa,
            departure_heads_pa,
            absolute_departure_heads_pa,
            pressure_heads_pa,
        }));
        Ok(self)
    }

    pub(super) fn reduce_boundaries(
        &self,
        boundaries: &[FixedPressure],
    ) -> Result<Option<Vec<FixedPressure>>, GraphError> {
        let Some(gravity) = &self.gravity else {
            return Ok(None);
        };
        boundaries
            .iter()
            .map(|boundary| {
                Ok(FixedPressure {
                    node: boundary.node,
                    pressure: Pressure::new(finite(
                        boundary.pressure.value() - gravity.reference_pressure(boundary.node),
                        "hydrostatic reference subtraction",
                    )?),
                })
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some)
    }

    pub(super) fn outward_head(&self, edge: usize, node: usize) -> f64 {
        let head = self
            .gravity
            .as_ref()
            .map_or(0.0, |gravity| gravity.departure_heads_pa[edge]);
        if self.branches[edge].from == node {
            head
        } else {
            -head
        }
    }

    pub(super) fn driving_pressure(
        &self,
        edge: usize,
        from: f64,
        to: f64,
    ) -> Result<f64, GraphError> {
        let difference = finite(from - to, "pressure difference")?;
        match &self.gravity {
            Some(gravity) => finite(
                difference + gravity.departure_heads_pa[edge],
                "gravity-adjusted pressure difference",
            ),
            None => Ok(difference),
        }
    }

    pub(super) fn outward_flow(
        &self,
        edge: usize,
        node: usize,
        pressure: f64,
        other: f64,
    ) -> Result<f64, GraphError> {
        let difference = finite(pressure - other, "pressure difference")?;
        let drive = if self.gravity.is_some() {
            finite(
                difference + self.outward_head(edge, node),
                "gravity-adjusted node pressure",
            )?
        } else {
            difference
        };
        self.flow_from_difference(edge, drive)
    }
}

fn density_departure(
    material: FluidTransportState,
    reference_density: f64,
    maximum_departure: f64,
    edge: usize,
) -> Result<f64, GraphError> {
    let departure = material.density_kg_m3() - reference_density;
    if !departure.is_finite() || (departure / reference_density).abs() > maximum_departure {
        Err(GraphError::InvalidBranch {
            branch: edge,
            reason: "gravity density exceeds the declared Boussinesq allowance",
        })
    } else {
        Ok(departure)
    }
}

// A local displacement, rather than two absolute potentials, preserves
// translation invariance. A zero acceleration component never forms 0*overflow.
fn head(
    acceleration: [f64; 3],
    from: [f64; 3],
    to: [f64; 3],
    density: f64,
) -> Result<f64, GraphError> {
    if density == 0.0 {
        return Ok(0.0);
    }
    let mut sum = Sum::default();
    for axis in 0..3 {
        if acceleration[axis] != 0.0 {
            sum.add(finite(
                acceleration[axis] * finite(to[axis] - from[axis], "gravity displacement")?,
                "specific gravitational head",
            )?)?;
        }
    }
    finite(density * sum.value()?, "material gravitational head")
}

#[cfg(test)]
mod tests;
