//! Original axial enthalpy reconstruction on each live momentum graph.

use super::*;
use fs_conduction::transient::deposition::advection::{MassBoundary, MassFlow};

impl<'a, F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'a, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Retain bounded axial enthalpy reconstruction through live material flow.
    ///
    /// Attach before the first advance. Every source epoch rebuilds its actual
    /// mass graph, binds the original liquid carrier when present, and attaches
    /// this SAME shared thermal plan. Source viscosity, buoyancy, phase drag and
    /// retained momentum continue to select the mass rates; the thermal owner
    /// alone supplies limited face h, conservative heat transfers and the
    /// combined conduction/reconstructed-advection timestep restriction.
    /// Both implicit feedback entries use these same reconstructed endpoint
    /// powers inside their original joule residual and nonlinear solve.
    ///
    /// Zero-rate interior faces required by the selected stencils retain
    /// geometric adjacency at rest and during reversal and count against
    /// `FlowEnthalpyInput::max_connections`. Only actual adjacent mobile voxels
    /// are connected; dormant transverse faces add no junction. The plan describes
    /// nonbranching axial cells, using their original contiguous centreline
    /// geometry, mass and caloric material. Junctions, a lost neighbor after
    /// exact-solid locking, and changed bindings refuse the candidate. External
    /// streams keep the original donor/inlet h; no external gradient is inferred.
    ///
    /// With stationary-skeleton liquid carriage, slopes reconstruct liquid h,
    /// retaining stationary solid latent inventory in mixture storage. This
    /// does not relax the original equal-density, pure-substance, source-support
    /// or Stokes-work restrictions. The optional bulk thermal-entropy limiter
    /// remains incompatible with liquid-only carriage. Explicit integration
    /// stays first order in time and still needs temporal refinement.
    /// Empty plans preserve the original donor-cell graph, fluxes and CFL.
    pub fn with_reconstruction(
        mut self,
        reconstruction: impl Into<Arc<LimitedEnthalpyReconstruction>>,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if self.applied.is_some() {
            return Err(
                invalid("live reconstruction must be attached before the first advance").into(),
            );
        }
        let reconstruction = reconstruction.into();
        self.reconstruction = (!reconstruction.is_empty()).then_some(reconstruction);
        // No initial momentum graph exists yet. Its actual topology and the
        // plan's bound cells are admitted at the first source epoch, without
        // an extra momentum solve or a fabricated placeholder flux.
        Ok(self)
    }
}

impl<F, R, O, T> MaterialFlowFeedback<'_, F, R, O, T> {
    /// Shared numerical plan, retained unchanged across actual source epochs.
    pub fn reconstruction(&self) -> Option<&LimitedEnthalpyReconstruction> {
        self.reconstruction.as_deref()
    }

    pub(super) fn mass_graph(
        &self,
        domain: &VoxelDomain,
        field: &crate::conjugate::FlowField,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Vec<MassFlow>, Vec<MassBoundary>), DepositionError> {
        mass_graph_with_zero_interior(
            domain,
            field,
            self.density_kg_m3,
            self.input,
            self.reconstruction.as_deref(),
            check,
        )
    }

    pub(super) fn bind_advection<P>(
        &self,
        graph: Advected<P>,
    ) -> Result<ReconstructedAdvected<P>, DepositionError> {
        self.bind_liquid_carrier(graph)?
            .with_optional_reconstruction(self.reconstruction.clone())
    }

    #[cfg(feature = "implicit-enthalpy-coupling")]
    pub(super) fn endpoint_transport(
        &self,
        graph: &Advected<()>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        fs_conduction::transient::deposition::backward_euler::TransportStepCoefficients<()>,
        DepositionError,
    > {
        match &self.reconstruction {
            Some(plan) => graph.endpoint_transport_with_reconstruction(thermal, plan, checkpoint),
            None => graph.endpoint_transport(thermal, checkpoint),
        }
    }
}
