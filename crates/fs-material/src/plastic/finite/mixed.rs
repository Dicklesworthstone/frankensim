//! Bulk-free multiplicative J2 returns for native mixed displacement/pressure.
//!
//! Isochoric flow leaves the original Hencky energy additively split:
//! `W = mu ||dev(log V_e)||² + H(alpha) + kappa/2 (log J)²`.
//! The SAME radial return, metric update and logarithm derivative therefore
//! supply a reduced incremental remainder independent of pressure. Its bulk
//! stress, storage and tangent are omitted at their original assembly sites;
//! no large penalty is evaluated and subtracted, no modulus is changed.
//!
//! The original log-squared pressure dual is admitted only on its locally
//! convex `log J < 1` branch. Exact incompressibility is an explicit alternative
//! and never evaluates the unused bulk penalty. Actual deformation retains all
//! original Jacobian and elastic/plastic history limits. The reconstructed bulk
//! volume must also satisfy the original Jacobian and combined elastic-strain
//! limits. Continuous spatial pressure enforces volume weakly; pointwise `J`
//! need not equal pressure-reconstructed `J(p)`.
//!
//! Linear, Voce and tabulated energetic hardening retain the same original
//! owner. This module does not convert a creep/Perzyna inspection backbone into
//! a rate-independent law or supply thermal, porous or ductile evolution.

use super::{FiniteJ2, FiniteJ2Error, FiniteJ2Response, FiniteJ2State, ReturnMap};
use crate::MaterialError;
use crate::hyper::Volumetric;
use crate::hyper::mixed::{
    MixedFiniteStrainPotential, MixedPointResponse, MixedPointState, PressurePotential,
    RemainderResponse, VolumetricTreatment, combine_response, combine_state, pressure_potential,
};

/// A borrowed original law and accepted history. Every Newton trial starts
/// here; observing this value cannot publish a plastic update.
#[derive(Clone, Copy, Debug)]
pub struct FrozenFiniteJ2Remainder<'a> {
    law: &'a FiniteJ2,
    accepted: &'a FiniteJ2State,
}

fn material_error(error: FiniteJ2Error) -> MaterialError {
    MaterialError::State {
        what: format!("mixed finite J2: {error}"),
    }
}

impl FiniteJ2 {
    fn check_mixed_selection(&self) -> Result<(), FiniteJ2Error> {
        if self.parameters.yield_stress_pa <= 0.0 || self.maximum_creep_alpha.is_some() {
            return Err(FiniteJ2Error::Parameters(
                "mixed J2 requires the original positive-yield rate-independent law, not a kinetic inspection backbone",
            ));
        }
        Ok(())
    }

    /// Freeze the original rate-independent law and accepted plastic metric.
    /// Linear, Voce and tabulated energetic hardening are preserved exactly.
    pub fn freeze_remainder<'a>(
        &'a self,
        accepted: &'a FiniteJ2State,
    ) -> Result<FrozenFiniteJ2Remainder<'a>, FiniteJ2Error> {
        self.check_mixed_selection()?;
        self.validate_history(accepted)?;
        Ok(FrozenFiniteJ2Remainder {
            law: self,
            accepted,
        })
    }

    /// The original return, with bulk-free stress, tangent and energy fields.
    ///
    /// `stored_energy_j_m3` is elastic deviatoric plus energetic-hardening
    /// storage. `incremental_potential_j_m3` adds the unchanged flow potential.
    /// The complete native history, Kirchhoff yield radius and physical yield
    /// heat are identical to the full response. A mixed consumer adds physical
    /// bulk storage from pressure and must not add `U(det F)` a second time.
    pub fn remainder_response(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<FiniteJ2Response, FiniteJ2Error> {
        self.check_mixed_selection()?;
        let returned = self.return_map_selected(f, old, None, false)?;
        self.response_from_return(f, old, returned)
    }

    /// Inspect accepted deviatoric/hardening storage without returning again
    /// or evaluating the unused volumetric penalty.
    pub fn stored_remainder_energy(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
    ) -> Result<f64, FiniteJ2Error> {
        self.check_mixed_selection()?;
        let endpoint = self.elastic_endpoint(f, state)?;
        self.check_elastic_domain(endpoint.log_j, &endpoint.deviator)?;
        self.remainder_energy(&endpoint.deviator, state.alpha)
    }

    /// Inspect native history, actual-deformation limits and equilibrium yield
    /// without evaluating bulk storage. Pressure and its reconstructed physical
    /// domain are additionally admitted by the frozen mixed response.
    pub fn validate_remainder_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        self.check_mixed_selection()?;
        self.validate_committed_selected(f, state, yield_tolerance_pa, false)
    }
}

impl FrozenFiniteJ2Remainder<'_> {
    fn prepare(&self, f: &[f64; 9]) -> Result<ReturnMap, MaterialError> {
        self.law.check_mixed_selection().map_err(material_error)?;
        self.law
            .return_map_selected(f, self.accepted, None, false)
            .map_err(material_error)
    }

    fn pressure_from_return(
        &self,
        returned: &ReturnMap,
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<PressurePotential, MaterialError> {
        if !pressure_pa.is_finite() {
            return Err(material_error(FiniteJ2Error::State(
                "mixed pressure must be finite",
            )));
        }
        let bulk = match treatment {
            VolumetricTreatment::ExactIncompressible => PressurePotential {
                complementary_pa: 0.0,
                bulk_stored_pa: 0.0,
                volume_ratio: 1.0,
                volume_increment: 0.0,
                compliance_per_pa: 0.0,
            },
            VolumetricTreatment::OriginalMaterial => {
                if returned.trial.log_j >= 1.0 {
                    return Err(material_error(FiniteJ2Error::State(
                        "actual log-squared volume requires log(det F) < 1 on the original convex branch",
                    )));
                }
                pressure_potential(
                    Volumetric::LogSquared {
                        kappa: self.law.parameters.bulk_modulus_pa,
                    },
                    pressure_pa,
                )?
            }
        };
        if bulk.volume_ratio <= self.law.limits.minimum_jacobian {
            return Err(material_error(FiniteJ2Error::Domain {
                quantity: "pressure-reconstructed J, strict lower bound",
                value: bulk.volume_ratio,
                limit: self.law.limits.minimum_jacobian,
            }));
        }
        let returned_deviator = returned.trial.deviator.map(|x| returned.scale * x);
        self.law
            .check_elastic_domain(
                super::log_one_plus(bulk.volume_increment),
                &returned_deviator,
            )
            .map_err(material_error)?;
        Ok(bulk)
    }

    fn remainder_from_return(
        &self,
        f: &[f64; 9],
        returned: &ReturnMap,
    ) -> Result<RemainderResponse, MaterialError> {
        Ok(RemainderResponse {
            potential_pa: returned.potential,
            first_piola_pa: returned.first_piola,
            tangent_pa: self
                .law
                .return_tangent(f, self.accepted, returned)
                .map_err(material_error)?,
        })
    }
}

impl MixedFiniteStrainPotential for FrozenFiniteJ2Remainder<'_> {
    fn validate_mixed(&self) -> Result<(), MaterialError> {
        self.law.check_mixed_selection().map_err(material_error)?;
        self.law
            .validate_history(self.accepted)
            .map_err(material_error)
    }

    fn remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        let returned = self.prepare(f)?;
        Ok((returned.first_piola, returned.potential))
    }

    fn remainder_response(&self, f: &[f64; 9]) -> Result<RemainderResponse, MaterialError> {
        let returned = self.prepare(f)?;
        self.remainder_from_return(f, &returned)
    }

    fn mixed_pressure_potential(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<PressurePotential, MaterialError> {
        let returned = self.prepare(f)?;
        self.pressure_from_return(&returned, pressure_pa, treatment)
    }

    fn mixed_state(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedPointState, MaterialError> {
        let returned = self.prepare(f)?;
        let bulk = self.pressure_from_return(&returned, pressure_pa, treatment)?;
        combine_state(
            f,
            pressure_pa,
            (returned.first_piola, returned.potential),
            bulk,
        )
    }

    fn mixed_response(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedPointResponse, MaterialError> {
        let returned = self.prepare(f)?;
        let bulk = self.pressure_from_return(&returned, pressure_pa, treatment)?;
        let remainder = self.remainder_from_return(f, &returned)?;
        combine_response(f, pressure_pa, remainder, bulk)
    }
}
