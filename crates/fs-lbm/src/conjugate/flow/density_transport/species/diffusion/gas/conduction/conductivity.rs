//! Original source conductivity on the existing Kirchhoff tensor interface.
//!
//! Point and interval queries retain the same original mixture. The integral
//! error is a numerical estimate in W/(m K), separate from the physical energy
//! ledger. All cells and native trials share one source point-query allowance.

use super::source::GasConductionCaloric;
use super::{CancelGate, ChtError, ConductivityTensor, TensorConductivityLaw, invalid, poll};
use fs_material::thermochemical_gas::conductivity::{
    GasConductivityError, GasConductivityQuadrature, PreparedThermochemicalConductivity,
};
use std::cell::Cell;

/// Explicit treatment of gas conductivity during the native heat interval.
#[derive(Clone, Copy, Debug)]
pub enum ThermochemicalConductivityModel {
    /// Opening original k is held per cell. A changing physical k retains
    /// first-order freezing error, including with an SDIRK2 storage formula.
    HeldAtOpening,
    /// Integrate the original composition-dependent WMS conductivity across
    /// every crossed original viscosity/conductivity source segment. Original
    /// endpoint k supplies the native Kirchhoff tangent; no fitted mixture
    /// table, midpoint substitution or extra dk/dT term is introduced.
    SourceIntegral {
        quadrature: GasConductivityQuadrature,
        /// Shared point evaluations across all gas cells, face integrals,
        /// Newton/backtracking trials and BE/SDIRK stages. Each original WMS
        /// point retains the original O(N^2) species-interaction calculation.
        /// Preparation and the opening/endpoint complete source images retain
        /// their separate bounded owners and cancellation checks.
        maximum_evaluations: usize,
    },
}

pub(super) struct GasConductivityBudget<'g> {
    model: ThermochemicalConductivityModel,
    gate: &'g CancelGate,
    evaluations: Cell<usize>,
    maximum_estimated_error_w_m_k: Cell<f64>,
}

impl<'g> GasConductivityBudget<'g> {
    pub(super) fn new(
        model: ThermochemicalConductivityModel,
        gate: &'g CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        if let ThermochemicalConductivityModel::SourceIntegral {
            quadrature,
            maximum_evaluations,
        } = model
            && (!quadrature.validate() || maximum_evaluations == 0)
        {
            return Err(invalid(
                "original gas conductivity needs valid quadrature controls and a positive point-query allowance",
            ));
        }
        Ok(Self {
            model,
            gate,
            evaluations: Cell::new(0),
            maximum_estimated_error_w_m_k: Cell::new(0.0),
        })
    }

    pub(super) fn evaluations(&self) -> usize {
        self.evaluations.get()
    }

    pub(super) fn maximum_estimated_error_w_m_k(&self) -> f64 {
        self.maximum_estimated_error_w_m_k.get()
    }

    fn checkpoint(&self) -> Result<(), GasConductivityError> {
        poll(self.gate).map_err(|_| GasConductivityError::Cancelled)
    }

    fn point(&self) -> Result<(), GasConductivityError> {
        self.checkpoint()?;
        let ThermochemicalConductivityModel::SourceIntegral {
            maximum_evaluations,
            ..
        } = self.model
        else {
            return Err(GasConductivityError::WorkLimit);
        };
        let next = self
            .evaluations
            .get()
            .checked_add(1)
            .filter(|&next| next <= maximum_evaluations)
            .ok_or(GasConductivityError::WorkLimit)?;
        self.evaluations.set(next);
        Ok(())
    }

    fn record_mean_error(&self, error_w_m_k: f64) -> Result<(), ChtError> {
        if !error_w_m_k.is_finite() || error_w_m_k < 0.0 {
            return Err(invalid(
                "original gas conductivity returned an unrepresentable integral error estimate",
            ));
        }
        self.maximum_estimated_error_w_m_k
            .set(self.maximum_estimated_error_w_m_k.get().max(error_w_m_k));
        Ok(())
    }
}

fn source_error(error: GasConductivityError) -> ChtError {
    match error {
        GasConductivityError::Cancelled => ChtError::Cancelled,
        other => invalid(format!("original gas conductivity refused: {other}")),
    }
}

fn isotropic(k: f64) -> Result<ConductivityTensor, ChtError> {
    ConductivityTensor::new([[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]])
}

pub(super) enum GasConductionConductivity<'q, 'g> {
    Held(ConductivityTensor),
    Source {
        law: PreparedThermochemicalConductivity,
        budget: &'q GasConductivityBudget<'g>,
        quadrature: GasConductivityQuadrature,
    },
}

impl<'q, 'g> GasConductionConductivity<'q, 'g> {
    pub(super) fn new(
        caloric: &GasConductionCaloric<'_, '_>,
        temperature_range_k: [f64; 2],
        opening_k_w_m_k: f64,
        budget: &'q GasConductivityBudget<'g>,
    ) -> Result<Self, ChtError> {
        poll(budget.gate)?;
        match budget.model {
            ThermochemicalConductivityModel::HeldAtOpening => {
                Ok(Self::Held(isotropic(opening_k_w_m_k)?))
            }
            ThermochemicalConductivityModel::SourceIntegral { quadrature, .. } => {
                let original = caloric.original_law();
                let law = PreparedThermochemicalConductivity::new(
                    original.model(),
                    original.context(),
                    temperature_range_k,
                    || budget.checkpoint(),
                )
                .map_err(source_error)?;
                Ok(Self::Source {
                    law,
                    budget,
                    quadrature,
                })
            }
        }
    }
}

impl TensorConductivityLaw for GasConductionConductivity<'_, '_> {
    fn conductivity_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ConductivityTensor, ChtError> {
        match self {
            Self::Held(tensor) => tensor.conductivity_at_temperature(temperature_k),
            Self::Source { law, budget, .. } => {
                let k = law
                    .conductivity_at_temperature(temperature_k, || budget.point())
                    .map_err(source_error)?;
                isotropic(k)
            }
        }
    }

    fn mean_conductivity(&self, from_k: f64, to_k: f64) -> Result<ConductivityTensor, ChtError> {
        match self {
            Self::Held(tensor) => tensor.mean_conductivity(from_k, to_k),
            Self::Source {
                law,
                budget,
                quadrature,
            } => {
                let mean = law
                    .mean_conductivity(from_k, to_k, *quadrature, || budget.point())
                    .map_err(source_error)?;
                budget.record_mean_error(mean.estimated_absolute_error_w_m_k)?;
                isotropic(mean.mean_w_m_k)
            }
        }
    }
}
