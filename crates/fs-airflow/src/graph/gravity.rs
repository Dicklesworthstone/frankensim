//! Gravity on the original pressure graph, with one conserved reference volume.
//!
//! Static pressure obeys `p_from - p_to + rho g·dx = loss(Q)`. Solve in
//! `p - rho0 g·(x-x0)` so the uniform-density hydrostatic state is balanced
//! without iteration or subtraction of large heads inside each edge solve.
//! Only `(rho-rho0) g·dx` remains in the numerical edge law. Different edge
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

/// Frozen gravity geometry and the original material-transport snapshots.
/// A numerical snapshot does not create source authority beyond its material.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphGravity {
    acceleration_m_s2: [f64; 3],
    positions_m: Vec<[f64; 3]>,
    reference: FluidTransportState,
    branch_materials: Vec<FluidTransportState>,
    maximum_relative_density_departure: f64,
    reference_pressures_pa: Vec<f64>,
    departure_heads_pa: Vec<f64>,
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
    /// Original frozen source samples in branch order. Density alone drives head.
    pub fn branch_materials(&self) -> &[FluidTransportState] {
        &self.branch_materials
    }
    /// Admitted `abs(rho-rho0)/rho0`; zero is the uniform-density specialization.
    pub const fn maximum_relative_density_departure(&self) -> f64 {
        self.maximum_relative_density_departure
    }
    /// Full gravitational pressure supply, signed in each edge's direction [Pa].
    pub fn pressure_heads_pa(&self) -> &[f64] {
        &self.pressure_heads_pa
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
            materials,
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
            materials,
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
        materials: Vec<FluidTransportState>,
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
        let mut pressure_heads_pa = Vec::with_capacity(self.branches.len());
        for (edge, (branch, material)) in self.branches.iter().zip(&materials).enumerate() {
            periodic_checkpoint(edge, cx)?;
            let departure = material.density_kg_m3() - rho0;
            if !departure.is_finite() || (departure / rho0).abs() > maximum_departure {
                return Err(GraphError::InvalidBranch {
                    branch: edge,
                    reason: "gravity density exceeds the declared Boussinesq allowance",
                });
            }
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
            let perturbation = head(
                acceleration,
                positions[branch.from],
                positions[branch.to],
                departure,
            )?;
            let full = finite(
                (reference_pressures_pa[branch.to] - reference_pressures_pa[branch.from])
                    + perturbation,
                "gravitational pressure head",
            )?;
            departure_heads_pa.push(perturbation);
            pressure_heads_pa.push(full);
        }
        checkpoint(cx)?;
        self.gravity = Some(Arc::new(GraphGravity {
            acceleration_m_s2: acceleration,
            positions_m: positions.to_vec(),
            reference,
            branch_materials: materials,
            maximum_relative_density_departure: maximum_departure,
            reference_pressures_pa,
            departure_heads_pa,
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
