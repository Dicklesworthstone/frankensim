//! Radiative neutron-capture activation on the original irradiation history.
//!
//! A private numerical parent has the actual removal rate k_f+k_c and selects
//! capture production or a fission-removal sink. The existing positive decay
//! propagator therefore retains the exponentially declining capture source,
//! including coincident parent-removal and daughter-decay rates. The numerical
//! coordinates are not additional physical nuclides or fitted half-lives.
//!
//! Irradiation remains the sole owner of parent depletion, prompt capture
//! energy, prompt fission energy and effective delayed fission reservoirs.
//! Its capture coefficient MUST describe prompt energy only; subsequent
//! activation-product heat belongs to the original bound decay inventory.
//! The binding declares one radiative-capture product (same Z, A+1, explicit
//! isomer). It supplies no cross-section evaluation, irradiation of products,
//! fission yields, neutron balance, emitted-particle inventory or material mass
//! transport. Spontaneous decay of the irradiation parent is not inferred.

use std::sync::Arc;

use fs_la::decay::{DecayAdvance, DecayEdge, DecayLimits, DecayPlan};

use super::production::DecayProductionStep;
use super::{
    Accumulator, DecayBudget, DecayError, DecayInventory, EnergyChannels, add, channel_array,
    checked, defect_bound, export_channels, product, sum,
};
use crate::FissionError;
use crate::irradiation::{
    GroupFlux, IrradiationBudget, IrradiationState, IrradiationStep, NuclideId,
    ReactionRateConstants,
};

/// Explicit dataset association between a typed capture product and one original
/// decay-model node. The string identity of that node is never parsed as a nuclide.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CaptureProduct {
    pub node: usize,
    pub nuclide: NuclideId,
    pub data_id: String,
}

/// Original independent budgets plus joint absolute accounting limits. Numerical
/// node/edge/work limits include the two private parent/removal coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaptureActivationBudget {
    pub irradiation: IrradiationBudget,
    /// Its extensive tolerances also bound the combined noncancelling defects.
    pub decay: DecayBudget,
}

/// Accepted coupled histories; clones share only immutable data and preparation.
/// At most 62 target nodes leave room for the two numerical coordinates.
#[derive(Clone, Debug)]
pub struct CaptureActivationState {
    parent: IrradiationState,
    products: DecayInventory,
    capture_product: CaptureProduct,
    prepared: Option<Arc<PreparedCapture>>,
}

#[derive(Debug)]
struct PreparedCapture {
    rates: ReactionRateConstants,
    duration_s: f64,
    limits: DecayLimits,
    capture_fraction: f64,
    fission_fraction: f64,
    plan: DecayPlan,
}

impl CaptureActivationState {
    /// Bind separately admitted histories to an explicit radiative-capture
    /// identity. Initial products may represent an authored irradiation history;
    /// neither their populations nor previous captures are inferred.
    pub fn new(
        parent: IrradiationState,
        products: DecayInventory,
        product: CaptureProduct,
    ) -> Result<Self, CaptureActivationError> {
        if products.populations().len() > fs_la::decay::MAX_DECAY_NODES - 2 {
            return Err(CaptureActivationError::Budget("activation target nodes"));
        }
        let original = parent.model().nuclide();
        if product.data_id.trim().is_empty()
            || product.node >= products.populations().len()
            || product.nuclide.atomic_number() != original.atomic_number()
            || original.mass_number().checked_add(1) != Some(product.nuclide.mass_number())
        {
            return Err(CaptureActivationError::Invalid(
                "an identified same-Z, A+1 capture product and existing target node are required",
            ));
        }
        Ok(Self {
            parent,
            products,
            capture_product: product,
            prepared: None,
        })
    }

    #[must_use]
    pub fn parent(&self) -> &IrradiationState {
        &self.parent
    }
    #[must_use]
    pub fn products(&self) -> &DecayInventory {
        &self.products
    }
    #[must_use]
    pub fn capture_product(&self) -> &CaptureProduct {
        &self.capture_product
    }

    /// Exact frozen-flux/rate trial, subject to the original owners' numerical
    /// range and accounting admission. Changing flux or temperature is admitted
    /// through the original irradiation data on every call. Preparation is reused
    /// only for the exact rate pair, duration and numerical limits.
    pub fn trial(
        &self,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: CaptureActivationBudget,
        mut checkpoint: impl FnMut() -> Result<(), CaptureActivationError>,
    ) -> Result<(Self, CaptureActivationStep), CaptureActivationError> {
        checkpoint()?;
        budget.decay.validate()?;
        let (parent, irradiation) = bridged(&mut checkpoint, FissionError::Cancelled, |check| {
            self.parent
                .trial(duration_s, flux, temperature_k, budget.irradiation, check)
        })?;
        let rates = bridged(&mut checkpoint, FissionError::Cancelled, |check| {
            self.parent
                .model()
                .rate_constants(flux, temperature_k, budget.irradiation, check)
        })?;
        let prepared =
            self.prepare(rates, duration_s, budget.decay.propagation, &mut checkpoint)?;
        let n = self.products.populations().len();
        let mut initial = self.products.populations().to_vec();
        initial.push(self.parent.remaining_parent_nuclei());
        initial.push(0.0);
        let evolved = bridged(
            &mut checkpoint,
            fs_la::decay::DecayError::Cancelled,
            |check| {
                prepared
                    .plan
                    .apply_with_cancel(&initial, || check().is_err())
            },
        )?;

        // Compare independent extensive reaction counts and the parent endpoint.
        // A numerical sink is only a count check, never bulk material loss.
        let numerical_captures = product(evolved.decays[n], prepared.capture_fraction)?;
        let numerical_fissions = product(evolved.decays[n], prepared.fission_fraction)?;
        let mut coupling = Accumulator::new(evolved.population_balance_residual.abs());
        for residual in [
            evolved.populations[n] - parent.remaining_parent_nuclei(),
            numerical_captures - irradiation.captures,
            numerical_fissions - irradiation.fissions,
            evolved.populations[n + 1] - irradiation.fissions,
            sum(&[
                evolved.decays[n],
                -irradiation.fissions,
                -irradiation.captures,
            ]),
            evolved.decays[n].mul_add(prepared.capture_fraction, -numerical_captures),
            evolved.decays[n].mul_add(prepared.fission_fraction, -numerical_fissions),
        ] {
            add(&mut coupling, checked(residual)?.abs())?;
        }
        let coupling_roundoff_nuclei = defect_bound(coupling)?;
        let mut produced_nuclei = vec![0.0; n];
        produced_nuclei[self.capture_product.node] = irradiation.captures;
        let mut population_balance = Accumulator::new(0.0);
        for &population in &evolved.populations[..n] {
            add(&mut population_balance, population)?;
        }
        for &population in self.products.populations() {
            add(&mut population_balance, -population)?;
        }
        add(&mut population_balance, -irradiation.captures)?;
        let target = DecayAdvance {
            populations: evolved.populations[..n].to_vec(),
            decays: evolved.decays[..n].to_vec(),
            population_balance_residual: checked(population_balance.value())?,
        };
        let (products, decay) = bridged(&mut checkpoint, DecayError::Cancelled, |mut check| {
            self.products.account_evolution(
                target,
                duration_s,
                Some(&produced_nuclei),
                0.0,
                budget.decay,
                &mut check,
            )
        })?;
        let mut population_error = Accumulator::new(0.0);
        for defect in [
            irradiation.parent_roundoff_nuclei,
            decay.population_roundoff_nuclei,
            coupling_roundoff_nuclei,
        ] {
            add(&mut population_error, defect)?;
        }
        let population_roundoff_nuclei = defect_bound(population_error)?;
        if population_roundoff_nuclei > budget.decay.nuclei_tolerance {
            return Err(CaptureActivationError::Budget(
                "combined activation population accounting",
            ));
        }
        let mut energy_error = Accumulator::new(0.0);
        add(&mut energy_error, irradiation.energy_roundoff_j)?;
        add(&mut energy_error, decay.energy_roundoff_j)?;
        let released = join_released(irradiation.released, decay.released, &mut energy_error)?;
        let energy_roundoff_j = defect_bound(energy_error)?;
        if energy_roundoff_j > budget.decay.energy_tolerance_j {
            return Err(CaptureActivationError::Budget(
                "combined activation energy accounting",
            ));
        }
        let next = Self {
            parent,
            products,
            capture_product: self.capture_product.clone(),
            prepared: Some(prepared),
        };
        checkpoint()?;
        Ok((
            next,
            CaptureActivationStep {
                irradiation,
                products: DecayProductionStep {
                    produced_nuclei,
                    decay,
                },
                released,
                energy_roundoff_j,
                population_roundoff_nuclei,
                coupling_roundoff_nuclei,
            },
        ))
    }

    /// Publish both histories and their cache only after the entire trial succeeds.
    pub fn advance(
        &mut self,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: CaptureActivationBudget,
        checkpoint: impl FnMut() -> Result<(), CaptureActivationError>,
    ) -> Result<CaptureActivationStep, CaptureActivationError> {
        let (next, report) = self.trial(duration_s, flux, temperature_k, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    fn prepare(
        &self,
        rates: ReactionRateConstants,
        duration_s: f64,
        limits: DecayLimits,
        checkpoint: &mut impl FnMut() -> Result<(), CaptureActivationError>,
    ) -> Result<Arc<PreparedCapture>, CaptureActivationError> {
        if let Some(prepared) = &self.prepared {
            if prepared.rates.fission_s_inv.to_bits() == rates.fission_s_inv.to_bits()
                && prepared.rates.capture_s_inv.to_bits() == rates.capture_s_inv.to_bits()
                && prepared.duration_s.to_bits() == duration_s.to_bits()
                && prepared.limits == limits
            {
                checkpoint()?;
                return Ok(Arc::clone(prepared));
            }
        }
        let n = self.products.populations().len();
        let total_rate = rates.total_s_inv();
        let fraction = |rate: f64| -> Result<f64, CaptureActivationError> {
            if rate == 0.0 {
                return Ok(0.0);
            }
            let value = rate / total_rate;
            if value == 0.0 || !value.is_finite() {
                return Err(CaptureActivationError::Arithmetic);
            }
            Ok(value)
        };
        let capture_fraction = fraction(rates.capture_s_inv)?;
        let fission_fraction = fraction(rates.fission_s_inv)?;
        let mut numerical_rates: Vec<_> = self
            .products
            .model()
            .nuclides()
            .iter()
            .map(|node| node.law.decay_constant().value())
            .collect();
        numerical_rates.extend([total_rate, 0.0]);
        let mut edges: Vec<_> = self
            .products
            .model()
            .branches()
            .iter()
            .map(|branch| DecayEdge {
                parent: branch.parent,
                daughter: branch.daughter,
                branch_fraction: branch.probability,
            })
            .collect();
        for (daughter, branch_fraction) in [
            (self.capture_product.node, capture_fraction),
            (n + 1, fission_fraction),
        ] {
            if branch_fraction > 0.0 {
                edges.push(DecayEdge {
                    parent: n,
                    daughter,
                    branch_fraction,
                });
            }
        }
        let plan = bridged(checkpoint, fs_la::decay::DecayError::Cancelled, |check| {
            DecayPlan::prepare(&numerical_rates, &edges, duration_s, limits, || {
                check().is_err()
            })
        })?;
        Ok(Arc::new(PreparedCapture {
            rates,
            duration_s,
            limits,
            capture_fraction,
            fission_fraction,
            plan,
        }))
    }
}

/// Accepted interval accounts. Combined release includes irradiation exactly
/// once and actual product decays; it is integrated energy, not endpoint power.
#[derive(Clone, Debug, PartialEq)]
pub struct CaptureActivationStep {
    pub irradiation: IrradiationStep,
    pub products: DecayProductionStep,
    pub released: EnergyChannels,
    pub energy_roundoff_j: f64,
    pub population_roundoff_nuclei: f64,
    pub coupling_roundoff_nuclei: f64,
}

/// Refusal without accepting either physical history or a candidate cache.
#[derive(Clone, Debug, PartialEq)]
pub enum CaptureActivationError {
    Invalid(&'static str),
    Budget(&'static str),
    Arithmetic,
    Irradiation(FissionError),
    Decay(DecayError),
    Propagation(fs_la::decay::DecayError),
    Cancelled,
}
impl From<FissionError> for CaptureActivationError {
    fn from(error: FissionError) -> Self {
        if error == FissionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Irradiation(error)
        }
    }
}
impl From<DecayError> for CaptureActivationError {
    fn from(error: DecayError) -> Self {
        if error == DecayError::Cancelled {
            Self::Cancelled
        } else {
            Self::Decay(error)
        }
    }
}
impl From<fs_la::decay::DecayError> for CaptureActivationError {
    fn from(error: fs_la::decay::DecayError) -> Self {
        if error == fs_la::decay::DecayError::Cancelled {
            Self::Cancelled
        } else {
            Self::Propagation(error)
        }
    }
}
impl core::fmt::Display for CaptureActivationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "capture activation refused: {self:?}")
    }
}
impl std::error::Error for CaptureActivationError {}

fn join_released(
    first: EnergyChannels,
    second: EnergyChannels,
    defects: &mut Accumulator,
) -> Result<EnergyChannels, CaptureActivationError> {
    let mut channels = [Accumulator::new(0.0); 3];
    for energy in [first, second] {
        for (total, value) in channels.iter_mut().zip(channel_array(energy)) {
            let error = add(total, value)?.abs();
            add(defects, error)?;
        }
    }
    // Export each retained low part separately: losses and gains in different
    // energy destinations must not cancel in the total-energy residual.
    let released = export_channels(channels, defects)?;
    let residual = checked(sum(&[
        released.total_j(),
        -first.total_j(),
        -second.total_j(),
    ]))?;
    add(defects, residual.abs())?;
    Ok(released)
}

// Preserve the caller's original checkpoint error through each owner's native
// cancellation type without exposing a partially accepted inner candidate.
fn bridged<T, E: Clone + Into<CaptureActivationError>>(
    checkpoint: &mut impl FnMut() -> Result<(), CaptureActivationError>,
    cancelled: E,
    operation: impl FnOnce(&mut dyn FnMut() -> Result<(), E>) -> Result<T, E>,
) -> Result<T, CaptureActivationError> {
    let mut checkpoint_error = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => Ok(()),
        Err(error) => {
            checkpoint_error = Some(error);
            Err(cancelled.clone())
        }
    });
    if let Some(error) = checkpoint_error {
        Err(error)
    } else {
        result.map_err(Into::into)
    }
}

#[cfg(test)]
#[path = "activation/tests.rs"]
mod tests;
