//! Generalized-Newtonian source response on the retained native MAC timeline.
//!
//! The same symmetric-stress quadrature owns the physical rate invariant and
//! molecular heat. Source viscosity is evaluated at accepted-start T/shear;
//! native inertia and the existing thermal split remain unchanged. No BGK tau,
//! cell speed, LES loss, Darcy loss or numerical time loss is a shear rate.

use super::*;

/// Immutable constitutive source on actual caloric and kinematic coordinates.
/// Implementations must retain original source claims and refuse unsupported
/// states. They must not mutate shared external physical history. Exact solids
/// retain their phase-first source semantics and receive zero native rate.
pub trait FlowRheologySource<F>: Sync {
    fn resolve(
        &self,
        voxel: usize,
        cell: &DepositedCell,
        temperature_k: f64,
        shear_rate_per_s: f64,
    ) -> Result<F, MaterialError>;
}
impl<F, Q> FlowRheologySource<F> for Q
where
    Q: Fn(usize, &DepositedCell, f64, f64) -> Result<F, MaterialError> + Sync,
{
    fn resolve(
        &self,
        voxel: usize,
        cell: &DepositedCell,
        temperature_k: f64,
        shear_rate_per_s: f64,
    ) -> Result<F, MaterialError> {
        self(voxel, cell, temperature_k, shear_rate_per_s)
    }
}

#[derive(Clone)]
pub(super) struct RheologyBinding<'a, F> {
    source: &'a dyn FlowRheologySource<F>,
    pub(super) rates: Vec<f64>,
}

impl<'a, F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'a, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Bind T/shear-dependent material response to the existing retained flow.
    ///
    /// Native symmetric-strain quadrature yields the cell-volume RMS
    /// `sqrt(2 D:D)` in physical 1/s, including graded edges and no-slip wall
    /// half cells. It observes the exact retained face velocities; no
    /// cell-speed reconstruction or constitutive relaxation parameter enters.
    /// Source coefficients are frozen over each accepted interval. The SAME
    /// viscosity drives full-stress momentum and its optional molecular heat;
    /// source k and admitted Darcy closure remain on the same material epoch.
    ///
    /// This explicitly time-lagged generalized-Newtonian constitutive law needs
    /// timestep/mesh refinement. It adds no elastic/yield history, nonlinear
    /// implicit solve, phase kinetics, moving solid skeleton or free surface.
    /// Attach retained momentum first; quasisteady flow has no accepted native
    /// starting velocity and refuses. Both pre-advance and actual endpoint T/rate
    /// are source-admitted. Any late failure discards the whole original thermal,
    /// mechanical and wrapped-observer candidate, including finite gas exchange.
    pub fn with_native_rheology(
        mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'a dyn FlowRheologySource<F>,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        self.check_map(thermal, &mut || Ok(()))?;
        self.check_unsteady_start(thermal, &mut || Ok(()))?;
        if self.unsteady.is_none() || self.rheology.is_some() {
            return Err(
                invalid("native rheology needs one source and retained unsteady momentum").into(),
            );
        }
        self.rheology = Some(RheologyBinding {
            source,
            rates: Vec::new(),
        });
        self.endpoint = self.sample_materials(thermal, &mut || Ok(()))?;
        admit_density(&self.endpoint, self.density_kg_m3)?;
        self.mushy_drag(&self.endpoint, &mut || Ok(()))?;
        self.check_phase_conduction(&self.endpoint, &mut || Ok(()))?;
        Ok(self)
    }

    /// Accepted endpoint rates in native voxel order. Source-start rates are
    /// retained separately on `applied()`, matching its physical coefficients.
    pub fn endpoint_shear_rate_per_s(&self) -> Option<&[f64]> {
        self.rheology
            .as_ref()
            .map(|binding| binding.rates.as_slice())
    }

    pub(super) fn sample_materials(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FlowMaterialSample<F, T>>, MaterialFlowError<O::Error>> {
        let Some(binding) = &self.rheology else {
            return sample(
                thermal,
                self.domain,
                self.input,
                self.gate,
                &mut self.resolve,
                self.exact_phase.then_some(self.mushy).flatten(),
                self.snapshot,
                check,
            );
        };
        let source = binding.source;
        let owner = self
            .unsteady
            .as_ref()
            .ok_or(invalid("native rheology lost its momentum owner"))?;
        let active = self.momentum_domain(&self.endpoint, check)?;
        let domain = active.as_ref().unwrap_or(self.domain);
        let mut rates = if domain.fluid_count() == 0 {
            vec![0.0; domain.cell_count()]
        } else {
            let fluid = self.reference.map(|_| self.reference_fluid()).transpose()?;
            owner.shear_rates(domain, fluid.as_ref(), self.config, self.gate)?
        };
        let phase_law = self.exact_phase.then_some(self.mushy).flatten();
        let mut resolve = |voxel: usize, cell: &DepositedCell, temperature: f64| {
            // Native strain still uses the previously admitted topology until
            // the phase handoff below. A newly exact-solid cell can inherit a
            // neighbouring edge contribution even when its own faces are at
            // rest. Its CURRENT caloric phase owns this source query: solids
            // receive zero, without clipping any live momentum or its energy.
            if let Some(law) = phase_law {
                if matches!(
                    law.classify_at_specific_enthalpy(cell.specific_h_j_kg)?,
                    PhaseMobility::ImmobileSolid { .. }
                ) {
                    rates[voxel] = 0.0;
                }
            }
            source.resolve(voxel, cell, temperature, rates[voxel])
        };
        let samples = sample(
            thermal,
            self.domain,
            self.input,
            self.gate,
            &mut resolve,
            phase_law,
            self.snapshot,
            check,
        )?;
        self.rheology.as_mut().expect("binding retained").rates = rates;
        Ok(samples)
    }
}
