//! Local isobaric calorics of one accepted multicomponent phase candidate.
//!
//! Original pure-liquid Cp sources supply the frozen-inventory response. The
//! SAME constrained Gibbs Hessian supplies phase redistribution: no perturbed
//! temperature flash, fixed-composition shortcut or independent mixture Cp fit.
//! This differentiates the retained seeded LOCAL branch at its fixed phase
//! count/support. It does not establish a global minimum, continue through a
//! disappearing phase or resolve a critical singularity. The candidate's
//! original phase-root accuracy still limits the response.

use super::{
    MulticomponentEquilibriumCandidateV1, MulticomponentEquilibriumControlV1,
    MulticomponentEquilibriumErrorV1, PhaseEvaluation, R, Variable, assemble_newton,
    largest_allocation, newton_direction, sqrt, sum,
};
use crate::solution::caloric::{BinaryCaloricErrorV1, PureLiquidHeatCapacityV1};
use fs_qty::Temperature;

/// Bounded scratch and residual admission for the original unshifted Hessian.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentCaloricResponseControlV1 {
    /// Includes the two dense matrices, vectors, copied phase responses and
    /// returned derivatives. No phase or species search is performed.
    pub max_workspace_bytes: usize,
    /// Componentwise backward error of H*direction + gradient. This is a
    /// numerical solve allowance, not a source or equilibrium error bound.
    pub relative_residual_tolerance: f64,
}
impl Default for MulticomponentCaloricResponseControlV1 {
    fn default() -> Self {
        Self {
            max_workspace_bytes: 16 * 1024 * 1024,
            relative_residual_tolerance: 1e-10,
        }
    }
}

/// Unsupported sources or an unresolved physical branch publish no response.
#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentCaloricResponseErrorV1 {
    Invalid(&'static str),
    WorkspaceBudget {
        required_bytes: usize,
    },
    Source(BinaryCaloricErrorV1),
    Equilibrium(MulticomponentEquilibriumErrorV1),
    /// The original joint Hessian is not positively resolved. A shifted or
    /// regularized Newton matrix cannot substitute for this physical derivative.
    SingularResponse,
    LinearResidual {
        maximum_scaled_residual: f64,
    },
    Unrepresentable(&'static str),
    Cancelled,
}
impl core::fmt::Display for MulticomponentCaloricResponseErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multicomponent local caloric response: {self:?}")
    }
}
impl std::error::Error for MulticomponentCaloricResponseErrorV1 {}
impl From<BinaryCaloricErrorV1> for MulticomponentCaloricResponseErrorV1 {
    fn from(error: BinaryCaloricErrorV1) -> Self {
        Self::Source(error)
    }
}
impl From<MulticomponentEquilibriumErrorV1> for MulticomponentCaloricResponseErrorV1 {
    fn from(error: MulticomponentEquilibriumErrorV1) -> Self {
        match error {
            MulticomponentEquilibriumErrorV1::Cancelled => Self::Cancelled,
            MulticomponentEquilibriumErrorV1::RegularizationLimit { .. } => Self::SingularResponse,
            other => Self::Equilibrium(other),
        }
    }
}
type Result<T> = core::result::Result<T, MulticomponentCaloricResponseErrorV1>;

/// A smooth LOCAL fixed-phase-count response, borrowing its exact candidate
/// and ordered original pure-liquid Cp curves. No global-equilibrium or
/// phase-onset qualification is inferred from this derivative.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentCaloricResponseV1<'a> {
    candidate: &'a MulticomponentEquilibriumCandidateV1,
    pure_liquid_cp: &'a [PureLiquidHeatCapacityV1],
    temperature_domain: [Temperature; 2],
    frozen_capacity_j_per_k: f64,
    redistribution_capacity_j_per_k: f64,
    capacity_j_per_k: f64,
    phase_component_derivatives: Box<[f64]>,
    phase_total_derivatives: Box<[f64]>,
    maximum_scaled_residual: f64,
}
impl<'a> MulticomponentCaloricResponseV1<'a> {
    pub const fn candidate(&self) -> &'a MulticomponentEquilibriumCandidateV1 {
        self.candidate
    }
    pub const fn pure_liquid_heat_capacities(&self) -> &'a [PureLiquidHeatCapacityV1] {
        self.pure_liquid_cp
    }
    /// Common original pair/pure-source support, without extrapolation.
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.temperature_domain
    }
    /// Original sum(n_i Cp_i(T)) [J/K]. With chi_ij=A_ij+B_ij/T, homogeneous
    /// excess enthalpy is independent of T, so its frozen contribution is zero.
    pub const fn frozen_heat_capacity_j_per_k(&self) -> f64 {
        self.frozen_capacity_j_per_k
    }
    /// Nonnegative additional dH/dT caused by phase redistribution [J/K].
    pub const fn redistribution_heat_capacity_j_per_k(&self) -> f64 {
        self.redistribution_capacity_j_per_k
    }
    /// Complete positive local isobaric dH/dT of the conserved inventory [J/K].
    /// This is not a constant-volume capacity or a phase-weighted frozen Cp.
    pub const fn isobaric_heat_capacity_j_per_k(&self) -> f64 {
        self.capacity_j_per_k
    }
    /// dn_(phase,i)/dT [mol/K], in original species order. Each species column
    /// sums to zero to roundoff; absent feed components remain exactly zero.
    pub fn phase_component_moles_temperature_derivatives_mol_per_k(
        &self,
        phase: usize,
    ) -> Option<&[f64]> {
        if phase >= self.candidate.phases().len() {
            return None;
        }
        let n = self.candidate.component_moles().len();
        Some(&self.phase_component_derivatives[phase * n..(phase + 1) * n])
    }
    /// dn_phase/dT [mol/K] in the retained candidate's phase order.
    pub fn phase_total_moles_temperature_derivatives_mol_per_k(&self) -> &[f64] {
        &self.phase_total_derivatives
    }
    /// Backward error of the unshifted linear response solve, not an accuracy
    /// certificate for the source data or the accepted equilibrium candidate.
    pub const fn maximum_scaled_residual(&self) -> f64 {
        self.maximum_scaled_residual
    }
}

impl MulticomponentEquilibriumCandidateV1 {
    /// Differentiate this accepted local phase branch at its original pressure
    /// convention and conserved component inventory. The caller must qualify
    /// the pure-liquid source evidence at that SAME pressure; identity and
    /// units alone cannot establish chemical/pressure compatibility.
    ///
    /// All curves must match canonical species order, share a datum, and
    /// support the actual temperature. The original entropy-scaled Hessian is
    /// assembled once and factored WITHOUT regularization. Singular or poorly
    /// resolved joint phase responses refuse, even if individual phases are
    /// locally stable. A one-phase candidate retains its declared homogeneous
    /// branch and does not acquire a global stability claim.
    ///
    /// For d=(phases-1)*present_species, work is O(d^3), scratch O(d^2), and
    /// cancellation is polled during source, assembly, solve and output work.
    pub fn isobaric_caloric_response<'a>(
        &'a self,
        pure_liquid_cp: &'a [PureLiquidHeatCapacityV1],
        controls: MulticomponentCaloricResponseControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentCaloricResponseV1<'a>> {
        poll(&mut cancelled)?;
        if !controls.relative_residual_tolerance.is_finite()
            || !(64.0 * f64::EPSILON..=1e-3).contains(&controls.relative_residual_tolerance)
        {
            return Err(MulticomponentCaloricResponseErrorV1::Invalid(
                "linear response controls",
            ));
        }
        let n = self.component_moles().len();
        if pure_liquid_cp.len() != n {
            return Err(MulticomponentCaloricResponseErrorV1::Invalid(
                "one original pure-liquid Cp curve per species",
            ));
        }
        let temperature = self.feed().temperature();
        let mut domain = self.feed().model().temperature_domain();
        let datum = pure_liquid_cp[0].reference_temperature();
        let mut frozen = super::super::Sum::default();
        for (i, source) in pure_liquid_cp.iter().enumerate() {
            poll(&mut cancelled)?;
            if source.species() != &self.feed().model().species()[i]
                || source.reference_temperature() != datum
            {
                return Err(MulticomponentCaloricResponseErrorV1::Invalid(
                    "ordered pure-liquid species or common datum",
                ));
            }
            let support = source.temperature_domain();
            domain[0] = Temperature::new(domain[0].value().max(support[0].value()));
            domain[1] = Temperature::new(domain[1].value().min(support[1].value()));
            let cp = source.heat_capacity_at(temperature)?.value();
            if self.component_moles()[i] > 0.0 {
                frozen.add(positive(
                    self.component_moles()[i] * cp,
                    "pure-liquid heat-capacity contribution",
                )?);
            }
        }
        let frozen_capacity_j_per_k = positive(frozen.value(), "frozen liquid heat capacity")?;
        let p = self.phases().len();
        let k = self
            .component_moles()
            .iter()
            .filter(|&&amount| amount > 0.0)
            .count();
        let d = (p - 1) * k;
        // Counts are already bounded by original model/candidate admission.
        // Account for all retained arrays before any numerical allocation.
        let required_bytes = (2 * d * d + 4 * d + p * k + 3 * k + p * n + p) * 8
            + p * core::mem::size_of::<PhaseEvaluation>()
            + d * core::mem::size_of::<Variable>()
            + 2 * k * core::mem::size_of::<usize>();
        if required_bytes > controls.max_workspace_bytes {
            return Err(MulticomponentCaloricResponseErrorV1::WorkspaceBudget { required_bytes });
        }
        let mut component_derivatives = vec![0.0; p * n];
        let mut phase_derivatives = vec![0.0; p];
        let mut redistribution = 0.0;
        let mut maximum_scaled_residual = 0.0;
        if d > 0 {
            let active: Vec<_> = self
                .component_moles()
                .iter()
                .enumerate()
                .filter_map(|(i, &amount)| (amount > 0.0).then_some(i))
                .collect();
            let roots: Vec<_> = active
                .iter()
                .map(|&i| sqrt(self.feed().mole_fractions()[i]))
                .collect();
            let allocations: Vec<_> = self
                .phases()
                .iter()
                .flat_map(|phase| {
                    active
                        .iter()
                        .map(move |&i| phase.component_moles()[i] / self.component_moles()[i])
                })
                .collect();
            let anchors: Vec<_> = (0..k)
                .map(|i| largest_allocation(&allocations, p, k, i))
                .collect();
            let variables: Vec<_> = (0..k)
                .flat_map(|i| {
                    let anchor = anchors[i];
                    (0..p)
                        .filter(move |&phase| phase != anchor)
                        .map(move |phase| Variable {
                            component: i,
                            phase,
                        })
                })
                .collect();
            let t = temperature.value();
            let rt = positive(R * t, "RT")?;
            let phases: Vec<_> = self
                .phases()
                .iter()
                .map(|phase| PhaseEvaluation {
                    state: phase.thermodynamics().clone(),
                    fraction: phase.mole_fraction_of_feed(),
                    reduced_g: phase.thermodynamics().mixing_gibbs_j_per_mol() / rt,
                })
                .collect();
            let mut matrix = vec![0.0; d * d];
            let mut factor = vec![0.0; d * d];
            let mut gradient = vec![0.0; d];
            let mut direction = vec![0.0; d];
            assemble_newton(
                &phases,
                &active,
                &roots,
                &allocations,
                &anchors,
                &variables,
                &mut matrix,
                &mut gradient,
                &mut cancelled,
            )?;
            // q_i=sqrt(z_i)*allocation_i. At the stationary point,
            // H q_T = -g_T, with T*g_T=-sqrt(z_i)*Delta(hbar_i)/(R*T).
            // Solve for T*q_T to avoid forming a potentially overflowing T².
            // Pure standard enthalpies cancel between phases of each species.
            for (row, v) in variables.iter().enumerate() {
                poll(&mut cancelled)?;
                let i = active[v.component];
                let h = phases[v.phase]
                    .state
                    .partial_molar_mixing_enthalpies_j_per_mol()[i];
                let reference_h = phases[anchors[v.component]]
                    .state
                    .partial_molar_mixing_enthalpies_j_per_mol()[i];
                gradient[row] = finite(
                    -roots[v.component] * ((h - reference_h) / rt),
                    "thermal chemical-potential derivative",
                )?;
            }
            let mut factorizations = 0;
            newton_direction(
                &matrix,
                &gradient,
                &mut factor,
                &mut direction,
                MulticomponentEquilibriumControlV1 {
                    max_regularization_steps: 1,
                    ..self.controls()
                },
                0,
                &mut factorizations,
                &mut cancelled,
            )?;
            maximum_scaled_residual =
                scaled_residual(&matrix, &gradient, &direction, &mut cancelled)?;
            if maximum_scaled_residual > controls.relative_residual_tolerance {
                return Err(MulticomponentCaloricResponseErrorV1::LinearResidual {
                    maximum_scaled_residual,
                });
            }
            // Cp_redistribution=n*R*(T*g_T)^T H^-1 (T*g_T).
            // The Cholesky norm makes its nonnegative sign explicit without
            // cancelling positive/negative species heat contributions.
            redistribution = capacity_from_factor(
                self.total_moles(),
                &factor,
                &direction,
                &mut gradient,
                &mut cancelled,
            )?;
            for (v, &value) in variables.iter().zip(&direction) {
                poll(&mut cancelled)?;
                let slope = finite(
                    (self.total_moles() * roots[v.component]) * (value / t),
                    "phase component-mole temperature derivative",
                )?;
                if slope == 0.0 && value != 0.0 {
                    return Err(MulticomponentCaloricResponseErrorV1::Unrepresentable(
                        "nonzero phase derivative underflows",
                    ));
                }
                component_derivatives[v.phase * n + active[v.component]] = slope;
            }
            // Reconstruct each anchor derivative from its actual paired
            // transfers. Conservation does not depend on subtracting two
            // independently rounded extensive thermal states.
            for (i, &species) in active.iter().enumerate() {
                poll(&mut cancelled)?;
                component_derivatives[anchors[i] * n + species] = finite(
                    -sum((0..p)
                        .filter(|&phase| phase != anchors[i])
                        .map(|phase| component_derivatives[phase * n + species])),
                    "conserved component derivative",
                )?;
            }
            for (phase, total) in phase_derivatives.iter_mut().enumerate() {
                *total = finite(
                    sum(component_derivatives[phase * n..(phase + 1) * n]
                        .iter()
                        .copied()),
                    "phase total-mole derivative",
                )?;
            }
        }
        let capacity_j_per_k = positive(
            frozen_capacity_j_per_k + redistribution,
            "total local heat capacity",
        )?;
        poll(&mut cancelled)?;
        Ok(MulticomponentCaloricResponseV1 {
            candidate: self,
            pure_liquid_cp,
            temperature_domain: domain,
            frozen_capacity_j_per_k,
            redistribution_capacity_j_per_k: redistribution,
            capacity_j_per_k,
            phase_component_derivatives: component_derivatives.into_boxed_slice(),
            phase_total_derivatives: phase_derivatives.into_boxed_slice(),
            maximum_scaled_residual,
        })
    }
}

fn scaled_residual(
    matrix: &[f64],
    gradient: &[f64],
    direction: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let d = gradient.len();
    let mut maximum = 0.0_f64;
    for row in 0..d {
        poll(cancelled)?;
        let residual = finite(
            sum(core::iter::once(gradient[row])
                .chain((0..d).map(|j| matrix[row * d + j] * direction[j]))),
            "linear response residual",
        )?;
        let scale = finite(
            gradient[row].abs() + sum((0..d).map(|j| (matrix[row * d + j] * direction[j]).abs())),
            "linear response scale",
        )?;
        if scale > 0.0 {
            maximum = maximum.max(residual.abs() / scale);
        } else if residual != 0.0 {
            return Err(MulticomponentCaloricResponseErrorV1::Unrepresentable(
                "unscaled linear residual",
            ));
        }
    }
    Ok(maximum)
}

fn capacity_from_factor(
    total_moles: f64,
    factor: &[f64],
    direction: &[f64],
    scratch: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let d = direction.len();
    let mut scale = 0.0_f64;
    for i in 0..d {
        poll(cancelled)?;
        scratch[i] = finite(
            sum((i..d).map(|j| factor[j * d + i] * direction[j])),
            "Cholesky thermal norm",
        )?;
        scale = scale.max(scratch[i].abs());
    }
    if scale == 0.0 {
        if direction.iter().any(|&v| v != 0.0) {
            return Err(MulticomponentCaloricResponseErrorV1::Unrepresentable(
                "nonzero thermal norm underflows",
            ));
        }
        return Ok(0.0);
    }
    let squares = sum(scratch.iter().map(|v| {
        let scaled = v / scale;
        scaled * scaled
    }));
    let amplitude = (sqrt(total_moles) * scale) * sqrt(R * squares);
    positive(
        amplitude * amplitude,
        "positive phase-redistribution heat capacity",
    )
}

fn finite(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MulticomponentCaloricResponseErrorV1::Unrepresentable(what))
    }
}
fn positive(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MulticomponentCaloricResponseErrorV1::Unrepresentable(what))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentCaloricResponseErrorV1::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
