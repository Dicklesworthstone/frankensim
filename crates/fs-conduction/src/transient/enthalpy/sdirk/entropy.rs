//! Entropy on the actual SDIRK2 tableau, including prescribed-stage reactions.
//!
//! With tau=gamma*dt and r=(1-gamma)*dt/tau, the native equations are
//! deltaE1=tau*Q1+rho1 and deltaE2=carry+tau*Q2+rho2. The carry is the
//! original first-stage physical power times its tableau weight; it is neither
//! an extra source nor an extrapolated caloric history.
//!
//! Physical entropy is r times stage 1's tau-weighted account plus stage 2's.
//! Storage is S(h2)-S(h0) only. The signed temporal remainder is
//! deltaS - r*sum(deltaE1/T1) - sum((deltaE2-carry)/T2), while algebraic
//! residual is r*sum(rho1/T1)+sum(rho2/T2). This is a tableau-stage residual
//! decomposition, not the final nodal residual divided by final temperature.
//! Neither raw stage is an independently accepted physical interval.

use super::super::entropy::{EnthalpyEntropyBalance, EnthalpyEntropyError};
use super::*;

impl EnthalpyBackwardEuler<'_, '_> {
    /// Advance the original SDIRK2 solve and admit its integrated thermal
    /// entropy before publishing either state or account.
    ///
    /// Both authentic stage operators, material charts, films and physical
    /// contacts are retained. Source and boundary data are sampled at their
    /// original tableau times. Earlier-stage heat never becomes an external
    /// source or an invalid effective-old enthalpy. The account adds one
    /// entropy pass per converged stage and no entropy work during Newton.
    ///
    /// Each positively weighted stage must have nonnegative bulk/contact,
    /// Robin and radiative closed-reservoir production within the caller's
    /// absolute J/K allowance. The whole
    /// interval must also pass storage-minus-input and closure admission.
    /// A temporal defect can have either sign; this is an actual-state check,
    /// not an unconditional entropy theorem or a time-accuracy certificate.
    /// Pointwise frozen/material ambient radiation retains each original
    /// black reservoir and the actual stage nodal heat. This closed account
    /// does not infer photon entropy. Mean-patch and reflected-enclosure
    /// entropy are refused before either stage is solved.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2_with_entropy(
        &self,
        cx: &Cx<'_>,
        inputs: [EnthalpySdirk2Input<'_>; 2],
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        entropy_tolerance_j_k: f64,
    ) -> Result<(EnthalpySdirk2Step, EnthalpyEntropyBalance), EnthalpyEntropyError> {
        if !entropy_tolerance_j_k.is_finite() || entropy_tolerance_j_k < 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "entropy tolerance must be finite and nonnegative",
            )
            .into());
        }
        let (step, entropy) =
            self.advance_sdirk2_accounted(cx, inputs, old_h, dt_s, config, true)?;
        let entropy = entropy.expect("requested native SDIRK2 entropy");
        entropy.admit(entropy_tolerance_j_k)?;
        poll(cx, old_h.len())?;
        Ok((step, entropy.balance))
    }
}

/// Accumulated in the native finish loop from its original double-double
/// storage, carry and residual. No duplicate matrix product or caloric query.
pub(super) struct StageEntropyTerms {
    storage_over_temperature: Dd,
    dirichlet: Dd,
    algebraic: Dd,
}

impl Default for StageEntropyTerms {
    fn default() -> Self {
        Self {
            storage_over_temperature: Dd::ZERO,
            dirichlet: Dd::ZERO,
            algebraic: Dd::ZERO,
        }
    }
}

impl StageEntropyTerms {
    pub(super) fn add(
        &mut self,
        delta: Dd,
        carry: Dd,
        residual: Dd,
        temperature_k: f64,
        prescribed: bool,
    ) -> Result<(), EnthalpyError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "entropy requires positive absolute temperatures",
            ));
        }
        let t = Dd::from_f64(temperature_k);
        self.storage_over_temperature = self.storage_over_temperature + (delta - carry) / t;
        let account = if prescribed {
            &mut self.dirichlet
        } else {
            &mut self.algebraic
        };
        *account = *account + residual / t;
        checked(self.storage_over_temperature)?;
        checked(*account)?;
        Ok(())
    }
}

pub(super) struct StageEntropy {
    // Only exact state entropy and physical contributions are consumed. The
    // raw backward-Euler temporal/closure fields do not describe an SDIRK stage.
    account: EnthalpyEntropyBalance,
    storage_over_temperature: Dd,
}

pub(super) fn finish_stage(
    context: &StepContext<'_, '_, '_, '_>,
    h: &[f64],
    system: &AssembledSystem,
    temperature: &[f64],
    dofs: &DofMap,
    radiation: Option<&(Vec<f64>, pointwise::PointwiseBoundaryReport)>,
    terms: StageEntropyTerms,
) -> Result<StageEntropy, EnthalpyError> {
    let mut account = super::super::entropy::evaluate(
        context, h, system, temperature, dofs, radiation,
    )?;
    account.dirichlet_entropy_input_j_k = checked(terms.dirichlet)?;
    account.algebraic_residual_j_k = checked(terms.algebraic)?;
    Ok(StageEntropy {
        account,
        storage_over_temperature: terms.storage_over_temperature,
    })
}

pub(super) struct SdirkEntropyEvaluation {
    pub(super) balance: EnthalpyEntropyBalance,
    // The first stage has already been multiplied by the actual w1/tau.
    generation: [[f64; 3]; 2],
}

impl SdirkEntropyEvaluation {
    fn admit(&self, tolerance_j_k: f64) -> Result<(), EnthalpyEntropyError> {
        for (stage, labels) in self.generation.iter().zip([
            [
                "SDIRK2 first-stage bulk/surface/contact",
                "SDIRK2 first-stage Robin",
                "SDIRK2 first-stage radiation material/field/bath",
            ],
            [
                "SDIRK2 second-stage bulk/surface/contact",
                "SDIRK2 second-stage Robin",
                "SDIRK2 second-stage radiation material/field/bath",
            ],
        ]) {
            for (&value_j_k, account) in stage.iter().zip(labels) {
                finite(value_j_k)?;
                if value_j_k < -tolerance_j_k {
                    return Err(EnthalpyEntropyError::NegativeProduction {
                        account,
                        value_j_k,
                        tolerance_j_k,
                    });
                }
            }
        }
        self.balance.admit(tolerance_j_k)
    }
}

pub(super) fn combine(
    first: StageEntropy,
    second: StageEntropy,
    first_scale: f64,
) -> Result<SdirkEntropyEvaluation, EnthalpyError> {
    let scale = Dd::from_f64(first_scale);
    let weight = |a, b| checked(scale * Dd::from_f64(a) + Dd::from_f64(b));
    let a = first.account;
    let b = second.account;
    let mut balance = EnthalpyEntropyBalance {
        thermal_storage_change_j_k: b.thermal_storage_change_j_k,
        dirichlet_entropy_input_j_k: weight(
            a.dirichlet_entropy_input_j_k,
            b.dirichlet_entropy_input_j_k,
        )?,
        robin_entropy_input_j_k: weight(a.robin_entropy_input_j_k, b.robin_entropy_input_j_k)?,
        radiation_entropy_input_j_k: weight(
            a.radiation_entropy_input_j_k, b.radiation_entropy_input_j_k,
        )?,
        conduction_generation_j_k: weight(
            a.conduction_generation_j_k,
            b.conduction_generation_j_k,
        )?,
        robin_generation_j_k: weight(a.robin_generation_j_k, b.robin_generation_j_k)?,
        radiation_generation_j_k: weight(
            a.radiation_generation_j_k, b.radiation_generation_j_k,
        )?,
        source_at_node_j_k: weight(a.source_at_node_j_k, b.source_at_node_j_k)?,
        neumann_at_node_j_k: weight(a.neumann_at_node_j_k, b.neumann_at_node_j_k)?,
        has_unclassified_inputs: a.has_unclassified_inputs || b.has_unclassified_inputs,
        temporal_defect_j_k: checked(
            Dd::from_f64(b.thermal_storage_change_j_k)
                - scale * first.storage_over_temperature
                - second.storage_over_temperature,
        )?,
        algebraic_residual_j_k: weight(a.algebraic_residual_j_k, b.algebraic_residual_j_k)?,
        balance_residual_j_k: 0.0,
    };
    let mut closure = Dd::from_f64(balance.thermal_storage_change_j_k);
    for term in [
        balance.dirichlet_entropy_input_j_k,
        balance.robin_entropy_input_j_k,
        balance.radiation_entropy_input_j_k,
        balance.conduction_generation_j_k,
        balance.robin_generation_j_k,
        balance.radiation_generation_j_k,
        balance.source_at_node_j_k,
        balance.neumann_at_node_j_k,
        balance.temporal_defect_j_k,
        balance.algebraic_residual_j_k,
    ] {
        closure = closure - Dd::from_f64(term);
    }
    balance.balance_residual_j_k = checked(closure)?;
    Ok(SdirkEntropyEvaluation {
        balance,
        generation: [
            [
                weight(a.conduction_generation_j_k, 0.0)?,
                weight(a.robin_generation_j_k, 0.0)?,
                weight(a.radiation_generation_j_k, 0.0)?,
            ],
            [
                b.conduction_generation_j_k,
                b.robin_generation_j_k,
                b.radiation_generation_j_k,
            ],
        ],
    })
}
