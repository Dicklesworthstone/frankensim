//! A prepared directed tangent for the original endpoint upwind flux owner.
use super::*;
use crate::enthalpy::advection::Advected;
use crate::enthalpy::advection::reconstruction::LimitedEnthalpyReconstruction;

pub(super) struct ImplicitAdvection {
    pub(super) flow: Advected<()>,
    // d(heat into row)/d(cell energy in column), in 1/s. Fixed mass and
    // prescribed flow make this constant even for nonlinear/latent T(h).
    pub(super) tangent: Vec<(usize, usize, f64)>,
    // The same compiled connectivity in carried-enthalpy coordinates, divided
    // by original cell masses. Liquid mobility changes its column slopes at
    // each endpoint, so these entries must never enter the static Jacobian.
    liquid_tangent: Vec<(usize, CellId, f64)>,
    reconstruction: Option<LimitedEnthalpyReconstruction>,
}

impl ImplicitAdvection {
    pub(super) fn endpoint(
        &self,
        thermal: &DepositionSimulation<'_>,
        index: &BTreeMap<CellId, usize>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(TransportStepCoefficients<()>, Vec<(usize, usize, f64)>), ImplicitConductionError>
    {
        let Some(reconstruction) = &self.reconstruction else {
            let mut tangent = Vec::new();
            tangent
                .try_reserve_exact(self.liquid_tangent.len())
                .map_err(|_| TridiagonalError::Capacity)?;
            for &(row, cell, coefficient) in &self.liquid_tangent {
                checkpoint()?;
                let slope = self.flow.implicit_carried_enthalpy_slope(&thermal.cells[&cell])?;
                if slope != 0.0 {
                    tangent.push((row, index[&cell], finite(coefficient * slope)?));
                }
            }
            return Ok((
                self.flow.endpoint_transport(thermal, checkpoint)?,
                tangent,
            ));
        };
        let (transport, rates) =
            self.flow
                .reconstructed_endpoint(thermal, reconstruction, true, &mut checkpoint)?;
        let mut tangent = Vec::new();
        tangent
            .try_reserve_exact(rates.len())
            .map_err(|_| TridiagonalError::Capacity)?;
        for (row, col, rate) in rates {
            checkpoint()?;
            let derivative = finite(rate / thermal.cells[&col].mass_kg)?;
            if derivative == 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            tangent.push((index[&row], index[&col], derivative));
        }
        Ok((transport, tangent))
    }
}

impl ImplicitConductionChain<'_> {
    /// Solve prescribed upwind enthalpy and nonlinear conduction together.
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
    /// Ordinary fixed-mass bulk enthalpy and the original stationary-skeleton
    /// liquid carrier are admitted. For liquid flow, the carrier tangent is zero
    /// on its latent plateau and one on its fully liquid branch, reevaluated at
    /// every trial. Thus liquid circulation does not advect stationary-solid
    /// latent energy. Fully solid flow still refuses; no mobility is inferred.
    /// Stationary species, selective outlets and an attached flow entropy
    /// history require separate closures and refuse. The native thermal entropy
    /// account marks mass-carried transfers unclassified; it does not infer
    /// entropy generation or pressure work. The supplied observer's history is
    /// unchanged, including on success.
    ///
    /// No explicit flow CFL is imposed on this implicit term. The selected
    /// backward-Euler or SDIRK2 time error and upwind spatial diffusion still
    /// require refinement; the existing adaptive h/T/phase controllers can
    /// refine these same endpoints. SDIRK2 re-evaluates transport at both stages.
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
        let liquid_tangent = if flow.liquid_carrier().is_some() {
            std::mem::take(&mut tangent)
                .into_iter()
                .map(|(row, col, coefficient)| (row, self.ids[col], coefficient))
                .collect()
        } else {
            Vec::new()
        };
        self.advection = Some(ImplicitAdvection {
            flow: flow.clone(),
            tangent,
            liquid_tangent,
            reconstruction: None,
        });
        checkpoint()?;
        Ok(self)
    }

    pub(super) fn with_implicit_advection_reconstruction(
        self,
        thermal: &DepositionSimulation<'_>,
        flow: &Advected<()>,
        reconstruction: &LimitedEnthalpyReconstruction,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !matches!(self.topology, LinearTopology::Graph(_)) {
            return Err(ImplicitConductionError::Invalid(
                "limited enthalpy reconstruction requires the graph solver",
            ));
        }
        let mut result = self.with_implicit_advection(thermal, flow, &mut checkpoint)?;
        if !reconstruction.is_empty() {
            // Existing source, geometry, topology and face-domain admission is
            // authoritative. No reconstruction is frozen at the initial state.
            flow.endpoint_transport_with_reconstruction(thermal, reconstruction, &mut checkpoint)?;
            result
                .advection
                .as_mut()
                .ok_or(DepositionError::Arithmetic)?
                .reconstruction = Some(reconstruction.clone());
        }
        checkpoint()?;
        Ok(result)
    }
}
