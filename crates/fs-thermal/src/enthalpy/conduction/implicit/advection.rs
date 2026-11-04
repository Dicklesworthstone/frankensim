//! A prepared directed tangent for the original endpoint upwind flux owner.
use super::*;
use crate::enthalpy::advection::Advected;

pub(super) struct ImplicitAdvection {
    pub(super) flow: Advected<()>,
    // d(heat into row)/d(cell energy in column), in 1/s. Fixed mass and
    // prescribed flow make this constant even for nonlinear/latent T(h).
    pub(super) tangent: Vec<(usize, usize, f64)>,
}

impl ImplicitConductionChain<'_> {
    /// Solve prescribed bulk upwind enthalpy and nonlinear conduction together.
    ///
    /// This binds an immutable snapshot of the original `Advected` graph. Its
    /// endpoint query retains native continuity, material, caloric-datum and
    /// inlet/outlet admission. Every Newton residual transports the proposed
    /// endpoint h, including latent energy at constant temperature. The same
    /// flow must not also be supplied as an explicit source observer.
    ///
    /// The exact directed tangent is compiled once in joule coordinates. A
    /// direct chain needs a tridiagonal combined operator; other flow graphs
    /// require `ImplicitConductionGraph`. A reference-coordinate derivative is
    /// retained even for admitted continuity roundoff. No symmetric diffusion
    /// surrogate or numerical enthalpy differences are used.
    ///
    /// Only ordinary fixed-mass bulk enthalpy is admitted. Liquid carriers,
    /// stationary species, selective outlets and an attached flow entropy
    /// history require separate closures and refuse. The native thermal entropy
    /// account marks mass-carried transfers unclassified; it does not infer
    /// entropy generation or pressure work. The supplied observer's history is
    /// unchanged, including on success.
    ///
    /// No explicit flow CFL is imposed on this backward-Euler term. First-order
    /// time error and upwind spatial diffusion still require refinement; the
    /// existing adaptive h/T/phase controller can refine these same endpoints.
    /// Other observed sources retain all their original timestep restrictions.
    pub fn with_implicit_advection(
        mut self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<()>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if self.advection.is_some() {
            return Err(ImplicitConductionError::Invalid(
                "implicit advection already bound",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let rates = flow.implicit_enthalpy_tangent(&mut checkpoint)?;
        // The original physical owner performs complete mass, source, storage
        // and law admission before any compiled tangent becomes usable.
        flow.endpoint_transport(thermal, &mut checkpoint)?;
        let mut tangent = Vec::with_capacity(rates.len());
        for (row, col, rate) in rates {
            checkpoint()?;
            let row = self.index[&row];
            let col = self.index[&col];
            if matches!(self.topology, LinearTopology::Chain) && row.abs_diff(col) > 1 {
                return Err(ImplicitConductionError::Invalid(
                    "nontridiagonal enthalpy transport requires the graph solver",
                ));
            }
            let derivative = finite(rate / thermal.cells[&self.ids[col]].mass_kg)?;
            if derivative == 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            tangent.push((row, col, derivative));
        }
        self.advection = Some(ImplicitAdvection {
            flow: flow.clone(),
            tangent,
        });
        checkpoint()?;
        Ok(self)
    }
}
