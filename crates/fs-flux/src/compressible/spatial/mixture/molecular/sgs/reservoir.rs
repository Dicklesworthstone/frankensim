//! Explicit local modeled exchange with complete, source-owned mixture baths.

use super::*;
use crate::compressible::spatial::molecular::sgs::reservoir_face_gradient;

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Add local SGS stress, heat and common-species diffusion at complete baths.
    ///
    /// The original source still owns molecular transport and the admitted
    /// exterior trace. Full physical face gradients supply the selected local
    /// model; midpoint density, composition and Cp supply its transport. Species
    /// enthalpy and stress work use the prescribed bath temperature and velocity,
    /// exactly as in the molecular reservoir law. All exchange enters the same
    /// conservative boundary inventory and stage diffusion bound.
    ///
    /// This is opt-in. Dynamic and Lagrangian closures are refused because no
    /// exterior test-filter field or pathline history has been authored. A
    /// pressure-only outlet retains its declared zero diffusive flux.
    pub fn with_local_sgs_reservoir_exchange(mut self) -> Result<Self, DissipationError> {
        let Some(SubgridTransport {
            closure: ViscosityClosure::Local {
                reservoir_exchange, ..
            },
            ..
        }) = self.sgs.as_mut()
        else {
            return Err(DissipationError::Invalid(
                "mixture reservoir SGS exchange requires an authored local closure",
            ));
        };
        *reservoir_exchange = true;
        Ok(self)
    }

    /// Whether complete reservoirs use the explicitly selected local SGS law.
    pub fn local_sgs_reservoir_exchange_enabled(&self) -> bool {
        self.sgs
            .as_ref()
            .is_some_and(SubgridTransport::local_reservoir_exchange_enabled)
    }
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::mixture::molecular) fn local_reservoir_exchange_enabled(
        &self,
    ) -> bool {
        matches!(
            &self.closure,
            ViscosityClosure::Local {
                reservoir_exchange: true,
                ..
            }
        )
    }

    /// Use the already admitted trace and the same stencil as molecular stress.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::spatial::mixture::molecular) fn add_reservoir_face(
        &self,
        gas: &MixtureGas,
        cell: usize,
        primitive: [Primitive; 2],
        partial_densities: [&[f64]; 2],
        axis: usize,
        distance: f64,
        tangential_gradient: [[f64; 3]; 2],
        boundary_on_left: bool,
        receipt: TransportFaceReport,
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let ViscosityClosure::Local {
            config,
            widths,
            filter_width,
            reservoir_exchange: true,
        } = &self.closure
        else {
            return Ok(receipt);
        };
        // Preserve the exact molecular path, including its original source
        // calls and report arithmetic, when the authored model is disabled.
        if config.transport.model.coefficient() == 0.0 {
            return Ok(receipt);
        }
        let gradient = reservoir_face_gradient(
            primitive[0].velocity,
            primitive[1].velocity,
            tangential_gradient,
            axis,
            distance,
        );
        let nu = match config.transport.model {
            SgsModel::Vreman { cv } => vreman_eddy_viscosity(cv, *widths, gradient),
            model => eddy_viscosity(model, *filter_width, gradient),
        }
        .map_err(|_| {
            DissipationError::Invalid("mixture reservoir subgrid gradient or viscosity")
        })?;
        let mu = (0.5 * primitive[0].density + 0.5 * primitive[1].density) * nu;
        if mu == 0.0 {
            return Ok(receipt);
        }
        let ratio = mu / receipt.dynamic_viscosity_pa_s;
        if !ratio.is_finite() || ratio > config.transport.maximum_eddy_to_molecular {
            return Err(DissipationError::Invalid(
                "mixture reservoir subgrid eddy/molecular viscosity ratio exceeded",
            ));
        }
        poll(cancelled)?;
        let boundary = usize::from(!boundary_on_left);
        let bath = primitive[boundary];
        let mut cv = [0.0; MAX_SPECIES];
        for (k, model) in gas.source().models().iter().enumerate() {
            if k % 16 == 0 {
                poll(cancelled)?;
            }
            if partial_densities[boundary][k] != 0.0 {
                cv[k] = (partial_densities[boundary][k] / model.molar_mass().value())
                    * (model
                        .evaluate_properties(Temperature::new(bath.temperature))
                        .map_err(FlowError::from)?
                        .cp()
                        .value()
                        - R);
            }
        }
        let bath_cv = sum(cv[..gas.species_count()].iter().copied());
        if !positive(bath_cv) {
            return Err(DissipationError::Invalid(
                "unrepresentable mixture reservoir subgrid heat capacity",
            ));
        }
        let mut cv_density = [self.heat_capacity_v_density[cell]; 2];
        cv_density[boundary] = bath_cv;
        flux::add(
            gas,
            flux::Face {
                partial_densities,
                primitive,
                cv_density,
                gradient,
                axis,
                distance,
                work_velocity: bath.velocity,
            },
            mu,
            None,
            None,
            Some(config.transport.turbulent_prandtl),
            Some(config.turbulent_schmidt),
            receipt,
            controls,
            out,
            cancelled,
        )
    }
}
