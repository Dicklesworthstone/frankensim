//! Source-integrated layer coefficients on the original explicit thermal clock.
use super::*;
use fs_thermal::enthalpy::{DepositionObserver, Event};

/// Adapter from the same prepared layer network to actual finite-cell thermal
/// evolution. No heat history is owned here: one native internal-link transfer
/// debits one cell and credits the other, and the native CFL uses that SAME G.
/// Existing heat sources, events and entropy accounting stay with their owners.
///
/// Conductances are frozen from each actual starting state. Their physical
/// endpoint derivatives are not substituted into the native secant CFL. This
/// positive finite-cell reduction retains the original temperature-range bound;
/// coefficient-freezing error and unresolved internal spatial storage require
/// independent time/spatial refinement. It supplies no implicit convergence or
/// exact transient spherical/cylindrical solution claim.
#[derive(Clone)]
pub struct LayerConductionObserver<'p, 'a, 'm> {
    model: &'p PreparedLayerConduction<'a, 'm>,
    config: LayerConductionConfig,
    samples: Vec<LayerLinkSample>,
}

impl<'p, 'a, 'm> LayerConductionObserver<'p, 'a, 'm> {
    pub fn new(
        model: &'p PreparedLayerConduction<'a, 'm>,
        thermal: &DepositionSimulation<'_>,
        config: LayerConductionConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, LayerSourceError> {
        let resolved = model.resolve(thermal, config, checkpoint)?;
        Ok(Self {
            model,
            config,
            samples: resolved.payload,
        })
    }

    /// Most recent actual starting-state flux evaluation, not a claim that the
    /// same watts hold at the returned endpoint. `model.resolve` admits any
    /// endpoint diagnostic through the original constitutive owner.
    pub fn last_samples(&self) -> &[LayerLinkSample] {
        &self.samples
    }

    fn admit_endpoint(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), LayerSourceError> {
        for link in &self.model.links {
            checkpoint()?;
            let temperatures = link.validate(thermal)?;
            // The physical interface solve requires both charts across the
            // complete interval. Reject a last-step source exit immediately,
            // even when no later substep would query that accepted endpoint.
            for layer in link.binding.layers {
                for temperature in temperatures {
                    layer.conductivity_curve().conductivity_w_m_k(temperature)?;
                }
            }
        }
        checkpoint()?;
        Ok(())
    }
}

impl DepositionObserver for LayerConductionObserver<'_, '_, '_> {
    type Error = LayerSourceError;

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        let resolved = self.model.resolve(thermal, self.config, checkpoint)?;
        self.samples = resolved.payload;
        Ok(resolved.conductances_w_k)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.admit_endpoint(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.admit_endpoint(thermal, checkpoint)
    }
}
