//! One constitutive query/cache owner for material-bound FEEC flow.

use super::{MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxRheologySample};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::RefCell;

pub(super) struct RheologyField<'a, G, R> {
    centroids: &'a [[f64; 2]],
    reference: &'a FluidTransportState,
    resolve: RefCell<R>,
    samples: RefCell<Vec<Option<MaterialFluxRheologySample<G>>>>,
    source_error: RefCell<Option<MaterialFluxError>>,
}

impl<'a, G: FluidTransport, R> RheologyField<'a, G, R>
where
    R: FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
{
    pub(super) fn new(
        centroids: &'a [[f64; 2]],
        reference: &'a FluidTransportState,
        resolve: R,
    ) -> Self {
        Self {
            centroids,
            reference,
            resolve: RefCell::new(resolve),
            samples: RefCell::new((0..centroids.len()).map(|_| None).collect()),
            source_error: RefCell::new(None),
        }
    }

    pub(super) fn is_valid(&self) -> bool {
        self.source_error.borrow().is_none()
    }

    pub(super) fn viscosity(&self, cell: usize, rate: f64) -> f64 {
        let fallback = self.reference.kinematic_viscosity_m2_s();
        if !self.is_valid() {
            // Only the already-started assembly may finish. Both callers
            // inspect is_valid before solving or publishing its candidate.
            return fallback;
        }
        let mut samples = self.samples.borrow_mut();
        if let Some(sample) = &samples[cell]
            && sample.point.shear_rate_per_s.to_bits() == rate.to_bits()
        {
            return sample.transport.kinematic_viscosity_m2_s();
        }
        let response = (|| {
            if !rate.is_finite() || rate < 0.0 {
                return Err(MaterialFluxError::Invalid(
                    "physical shear rate must be finite and nonnegative",
                ));
            }
            let point = MaterialFluxRheologyPoint {
                cell,
                position_m: self.centroids[cell],
                shear_rate_per_s: rate,
            };
            let material = self.resolve.borrow_mut()(point, self.reference)?;
            let transport = material.transport_state()?;
            if transport.density_kg_m3() != self.reference.density_kg_m3() {
                return Err(MaterialFluxError::Invalid(
                    "rheological flow must retain its constant reference density",
                ));
            }
            Ok(MaterialFluxRheologySample {
                point,
                material,
                transport,
            })
        })();
        match response {
            Ok(sample) => {
                let nu = sample.transport.kinematic_viscosity_m2_s();
                samples[cell] = Some(sample);
                nu
            }
            Err(error) => {
                *self.source_error.borrow_mut() = Some(error);
                fallback
            }
        }
    }

    // Preserve the original source error first, then the solve's own refusal.
    // Only an accepted solve requires a complete final set of material owners.
    pub(super) fn finish<T>(
        self,
        result: Result<T, MaterialFluxError>,
    ) -> Result<(T, Vec<MaterialFluxRheologySample<G>>), MaterialFluxError> {
        if let Some(error) = self.source_error.into_inner() {
            return Err(error);
        }
        let result = result?;
        let samples = self
            .samples
            .into_inner()
            .into_iter()
            .map(|sample| {
                sample.ok_or(MaterialFluxError::Invalid(
                    "accepted flow is missing a cell material response",
                ))
            })
            .collect::<Result<_, _>>()?;
        Ok((result, samples))
    }
}
