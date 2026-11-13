//! Joint physical response and carrier-sensitivity resolution.

use super::super::{
    ObservationBudget, admit_storage_with_retained, input_work, minimum_observation_work,
};
use super::sensitivity::{admit_port_tolerances, observe_pair, ode_carrier_derivative};
use super::{
    C64, PhsPortSidebandCarrierDerivativeResponse, PhsPortSidebandControls, PhsPortSidebandError,
    PhsPortSidebandSource, Result, SidebandScales, admit_carrier, admit_forcing, admit_source,
    input_forcing, poll,
};
use core::cell::Cell;
use fs_orbit::sideband::sensitivity::{
    SidebandCarrierDerivativeControls,
    refinement::{
        SidebandCarrierRefinementCheck, SidebandCarrierRefinementControls,
        solve_sideband_carrier_derivative_adaptive,
    },
};

mod resonance;
pub use resonance::{
    AdaptivePhsPortResonanceResponse, PhsPortAdaptiveResonanceRequest,
    PhsPortResonanceRefinementCheck,
};

/// A complete physical pair whose original state equations and carrier
/// derivatives passed independent sampled and actual inter-order checks.
#[derive(Clone, Debug)]
pub struct AdaptivePhsPortSidebandCarrierDerivativeResponse {
    pub response: PhsPortSidebandCarrierDerivativeResponse,
    pub refinement_trace: Vec<SidebandCarrierRefinementCheck>,
}

impl PhsPortSidebandSource<'_> {
    /// Resolve the full original state and exact carrier derivative together.
    ///
    /// Both probe arrays use this source's INITIAL signed support. G maps them
    /// once; the generic owner zero-pads that same physical family and actually
    /// solves successive orders. None holds the authored probe phasors fixed.
    /// Full original-equation, derivative-equation, omitted-band, independent
    /// grid and separate state/derivative changes must all pass.
    ///
    /// The original base, phase and clock stay fixed. After state acceptance,
    /// observe BOTH complete fields through the original G^T Hxx action. Their
    /// independent physical output order remains authored; state refinement
    /// never clips the accepted state or silently expands the output bandwidth.
    /// Observation grids grow within the original maximum. The adapter's work
    /// and OUTPUT observation action caps are shared by both observations;
    /// generic base/tangent/assessment callbacks retain their separate budgets.
    /// The core and added derivative work caps are independently cumulative
    /// across attempted orders, and their payload caps add to the adapter cap.
    /// Cloned observer controls and retained joint history are charged before
    /// allocation. These finite sampled checks do not enclose an infinite tail
    /// or prove Floquet stability or sensitivity along a changing base branch.
    ///
    /// # Errors
    /// Invalid original input/slope or tolerances; generic joint resolution or
    /// numerical refusal; final output/derivative refusal; resource exhaustion
    /// or cancellation. No partial pair escapes and the source remains reusable.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_adaptive_with_carrier_derivative(
        &mut self,
        carrier: f64,
        forcing: &[C64],
        forcing_angular_derivative: Option<&[C64]>,
        scales: &SidebandScales,
        derivative_controls: SidebandCarrierDerivativeControls,
        derivative_envelope_tolerances: &[f64],
        refinement: &SidebandCarrierRefinementControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptivePhsPortSidebandCarrierDerivativeResponse> {
        poll(&mut cancelled)?;
        let system = self.problem.original().system;
        let (n, m, initial) = (
            system.state_dim(),
            system.port_dim(),
            self.controls.sideband.harmonics,
        );
        admit_carrier(self.orbit, carrier, &self.controls)?;
        admit_forcing(system, forcing, &self.controls)?;
        if let Some(slope) = forcing_angular_derivative {
            admit_forcing(system, slope, &self.controls)?;
        }
        admit_port_tolerances(m, derivative_envelope_tolerances)?;
        let retained =
            8 * m as u128 + 8 * core::mem::size_of::<SidebandCarrierRefinementCheck>() as u128;
        admit_storage_with_retained(
            n,
            m,
            &self.controls,
            self.controls.initial_samples,
            retained,
        )?;
        let mut budget = ObservationBudget {
            actions: 0,
            work: input_work(n, m, initial)
                * if forcing_angular_derivative.is_some() {
                    2
                } else {
                    1
                },
            retained_bytes: retained,
        };
        admit_observations(n, m, &self.controls, &budget, 2, 0)?;
        let derivative_work = Cell::new(0_u128);
        let exhausted = Cell::new(false);
        let resolved = {
            let rhs = input_forcing(system, forcing, initial, &mut cancelled)?;
            let slope = match forcing_angular_derivative {
                Some(slope) => input_forcing(system, slope, initial, &mut cancelled)?,
                None => vec![C64::ZERO; rhs.len()],
            };
            let solved = solve_sideband_carrier_derivative_adaptive(
                &self.problem,
                self.orbit,
                carrier,
                &rhs,
                &slope,
                self.controls.sideband,
                scales,
                derivative_controls,
                refinement,
                |_, out| {
                    let work = derivative_work.get() + 16 * n as u128 * n as u128;
                    if budget.work + work > self.controls.max_observation_work as u128 {
                        exhausted.set(true);
                        return false;
                    }
                    derivative_work.set(work);
                    ode_carrier_derivative(n, out)
                },
                &mut cancelled,
            );
            if exhausted.get() {
                return Err(PhsPortSidebandError::Budget(
                    "joint adaptive physical port derivative exceeds the arithmetic cap",
                ));
            }
            solved?
        };
        budget.work += derivative_work.get();
        let observer = observer_controls(self, resolved.harmonics, retained)?;
        admit_observations(n, m, &observer, &budget, 2, 0)?;
        let response = observe_pair(
            self,
            carrier,
            resolved.response,
            &observer,
            derivative_envelope_tolerances,
            resolved.base_residual,
            resolved.tangent_sampling_defect,
            &mut budget,
            &mut cancelled,
        )?;
        Ok(AdaptivePhsPortSidebandCarrierDerivativeResponse {
            response,
            refinement_trace: resolved.refinement_trace,
        })
    }
}

fn observer_controls(
    source: &PhsPortSidebandSource<'_>,
    harmonics: usize,
    retained: u128,
) -> Result<PhsPortSidebandControls> {
    let c = &source.controls;
    if harmonics < c.sideband.harmonics || harmonics > 64 {
        return Err(PhsPortSidebandError::Invalid(
            "bounded nondecreasing physical state order required",
        ));
    }
    let largest = (source.orbit.coeffs[0].len() - 1)
        .max(harmonics)
        .max(c.output_harmonics);
    let initial = c.initial_samples.max((2 * largest + 1).next_power_of_two());
    if initial > c.max_samples / 4 {
        return Err(PhsPortSidebandError::Budget(
            "joint refined physical fields need three resolved observation grids within the maximum",
        ));
    }
    let system = source.problem.original().system;
    let (n, m) = (system.state_dim(), system.port_dim());
    // Admit the enlarged state arrays and owned observer tolerance clone before
    // allocating that clone. The caller includes its 8*m bytes in retained.
    let enlarged = 64 * n as u128 * (harmonics - c.sideband.harmonics) as u128;
    admit_storage_with_retained(n, m, c, initial, retained + enlarged)?;
    let mut observer = c.clone();
    observer.sideband.harmonics = harmonics;
    observer.initial_samples = initial;
    admit_source(system, source.orbit, &observer)?;
    Ok(observer)
}

fn admit_observations(
    n: usize,
    m: usize,
    controls: &PhsPortSidebandControls,
    budget: &ObservationBudget,
    fields: usize,
    extra_work: u128,
) -> Result<()> {
    if budget.actions as u128 + fields as u128 * 14 * controls.initial_samples as u128
        > controls.max_hessian_actions as u128
    {
        return Err(PhsPortSidebandError::Budget(
            "joint physical observations exceed the cumulative Hessian-action cap",
        ));
    }
    if budget.work
        + extra_work
        + fields as u128 * minimum_observation_work(n, m, controls.initial_samples)
        > controls.max_observation_work as u128
    {
        return Err(PhsPortSidebandError::Budget(
            "joint physical observations exceed the cumulative arithmetic cap",
        ));
    }
    admit_storage_with_retained(
        n,
        m,
        controls,
        controls.initial_samples,
        budget.retained_bytes,
    )
}
