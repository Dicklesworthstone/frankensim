//! Holomorphic free-mode queries on the existing complete physical polynomial.
//! T(s)=A0+s*A1+s²*A2; no local wall or material-memory denominator is formed.

use super::{CoupledPencil, C64, Dd, finite, original_memory};
use fs_la::nonlinear_eig::{NonlinearEigError, NonlinearPole, PoleEvaluation, newton_pole_with_cancel};
use fs_la::nonlinear_eig::beyn::{BeynError, BeynReport, beyn_with_cancel};
pub use fs_la::nonlinear_eig::PoleOptions as MechanicalPoleOptions;
pub use fs_la::nonlinear_eig::beyn::{
    BeynOptions as MechanicalSpectrumOptions, Ellipse as MechanicalResonanceRegion,
};

/// Additional resources for original sparse callback evaluation and observation.
/// The shared solver separately owns its dimension, work and numerical payload.
#[derive(Clone, Copy, Debug)]
pub struct MechanicalPoleLimits {
    /// Total polynomial evaluations AND final original-row audits per query.
    pub max_evaluations: usize,
    /// Callback matrices, maps and all possible returned mode coordinates;
    /// excludes the borrowed prepared owner and shared solver's own budget.
    pub max_payload_bytes: usize,
    /// Maximum coefficient-scaled backward residual in every original row.
    pub physical_residual_tolerance: f64,
}

/// Original wall motion at the returned arbitrary common modal amplitude.
#[derive(Clone, Debug)]
pub struct MechanicalPoleWall {
    pub displacement_m: Vec<C64>,
    pub normal_velocity_m_s: Vec<C64>,
    /// Every original z_j, satisfying z'_j=(sqrt(k_j)*w-z_j)/tau_j.
    pub relaxation_memory: Vec<Vec<C64>>,
    /// Every original elastic extension sqrt(k_j)*w-z_j. Retaining it avoids
    /// cancellation of two large nearly equal memories in stiff/high-Q modes.
    pub relaxation_extension: Vec<Vec<C64>>,
}

/// One numerical free mode, without source/receiver observability assumptions.
#[derive(Clone, Debug)]
pub struct MechanicalResonancePole {
    /// Signed Laplace rate in exp(s*t) [1/s]; growth is never clamped away.
    pub rate_per_s: C64,
    /// Full free-pressure/wall/extension seed order, largest magnitude one.
    pub coordinates: Vec<C64>,
    /// Original hierarchical FE pressure coefficients, including exact zeros
    /// at pressure-release constraints. This is not a cycle-averaged field.
    pub pressure_pa: Vec<C64>,
    pub walls: Vec<MechanicalPoleWall>,
    /// Maximum independently checked original row backward residual.
    pub physical_residual: f64,
    pub relative_step: f64,
    pub iterations: usize,
}

/// Numerical modes in the requested ellipse; no continuum-completeness claim.
#[derive(Clone, Debug)]
pub struct MechanicalResonanceSpectrum {
    pub poles: Vec<MechanicalResonancePole>,
    pub unrefined_rates_per_s: Vec<C64>,
    pub argument_integrals: [C64; 2],
    pub argument_count: usize,
    pub rank: usize,
    pub nodes: usize,
    pub moment_relative_change: f64,
    pub minimum_reciprocal_condition: f64,
}

/// No partial physical mode or spectrum is published on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MechanicalPoleError {
    Invalid(&'static str),
    Budget { what: &'static str, requested: u128, maximum: u128 },
    Pole(NonlinearEigError),
    Contour(BeynError),
    Residual { row: usize, observed: f64, maximum: f64 },
    Cancelled,
}
impl core::fmt::Display for MechanicalPoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete acoustic mechanical pole refused: {self:?}")
    }
}
impl std::error::Error for MechanicalPoleError {}

fn poll(cancelled: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), MechanicalPoleError> {
    if cancelled() { Err(MechanicalPoleError::Cancelled) } else { Ok(()) }
}
fn cap(what: &'static str, requested: u128, maximum: u128) -> Result<(), MechanicalPoleError> {
    if requested > maximum {
        Err(MechanicalPoleError::Budget { what, requested, maximum })
    } else { Ok(()) }
}

struct Query<'a> {
    source: &'a CoupledPencil,
    original_to_free: Vec<usize>,
    limits: MechanicalPoleLimits,
    evaluations: usize,
}
impl<'a> Query<'a> {
    fn new(source: &'a CoupledPencil, limits: MechanicalPoleLimits) -> Result<Self, MechanicalPoleError> {
        if limits.max_evaluations == 0 || !limits.physical_residual_tolerance.is_finite()
            || limits.physical_residual_tolerance <= 0.0 {
            return Err(MechanicalPoleError::Invalid("positive pole callback limits required"));
        }
        // This entry is only for the complete REAL causal mechanical owner.
        // A generic complex constant admittance has no inferred continuation.
        if (0..source.blocks[3].nrows()).any(|i| source.blocks[3].row(i).1.iter().any(|&v| v != 0.0)) {
            return Err(MechanicalPoleError::Invalid("real causal boundary polynomial required"));
        }
        let n = source.free_dofs.len() as u128;
        let full = source.blocks[0].nrows() as u128;
        cap("pole callback and result bytes", 64*n*n + 96*n*full + 128*full,
            limits.max_payload_bytes as u128)?;
        let mut original_to_free = vec![usize::MAX; full as usize];
        for (i, &original) in source.free_dofs.iter().enumerate() {
            original_to_free[original] = i;
        }
        Ok(Self { source, original_to_free, limits, evaluations: 0 })
    }

    fn admit(&mut self, cancelled: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), MechanicalPoleError> {
        poll(cancelled)?;
        cap("pole callback evaluations", self.evaluations as u128 + 1,
            self.limits.max_evaluations as u128)?;
        self.evaluations += 1;
        Ok(())
    }

    fn evaluate(&mut self, s: C64, cancelled: &mut dyn FnMut() -> bool) -> Result<PoleEvaluation, MechanicalPoleError> {
        self.admit(cancelled)?;
        if !finite(s) || !finite(s*s) {
            return Err(MechanicalPoleError::Invalid("finite representable pole argument required"));
        }
        let source = self.source;
        let n = source.free_dofs.len();
        let mut matrix = vec![C64::ZERO; n*n];
        let mut derivative = vec![C64::ZERO; n*n];
        let mut residual_scale = 0.0_f64;
        for (i, &row) in source.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let mut envelope = 0.0;
            for (block, multiplier, slope) in [
                (&source.blocks[0], C64::ONE, C64::ZERO),
                (&source.blocks[2], s, C64::ONE),
                (&source.blocks[1], s*s, s.scale(2.0)),
            ] {
                let (columns, values) = block.row(row);
                for (&col, &value) in columns.iter().zip(values) {
                    let j = self.original_to_free[col];
                    if j != usize::MAX {
                        let scale = value * source.scales[i] * source.scales[j];
                        matrix[i*n+j] = matrix[i*n+j] + multiplier.scale(scale);
                        derivative[i*n+j] = derivative[i*n+j] + slope.scale(scale);
                        envelope += multiplier.abs() * scale.abs();
                    }
                }
            }
            residual_scale = residual_scale.max(envelope);
        }
        if !residual_scale.is_finite() || residual_scale <= 0.0
            || matrix.iter().chain(&derivative).any(|&z| !finite(z)) {
            return Err(MechanicalPoleError::Invalid("unrepresentable complete pole evaluation"));
        }
        Ok(PoleEvaluation { matrix, derivative, residual_scale })
    }

    fn observe(&mut self, pole: NonlinearPole, cancelled: &mut impl FnMut() -> bool) -> Result<MechanicalResonancePole, MechanicalPoleError> {
        self.admit(cancelled)?;
        let s = pole.value;
        if s == C64::ZERO || !finite(s) {
            return Err(MechanicalPoleError::Invalid("zero-rate pressure/gauge classification is not a finite resonance query"));
        }
        let source = self.source;
        let mut coordinates: Vec<_> = pole.vector.iter().zip(&source.scales)
            .map(|(&v, &scale)| v.scale(scale)).collect();
        let amplitude = coordinates.iter().map(|v| v.abs()).fold(0.0, f64::max);
        if !amplitude.is_finite() || amplitude == 0.0 {
            return Err(MechanicalPoleError::Invalid("original pole normalization is unrepresentable"));
        }
        for v in &mut coordinates { *v = *v / C64::from_re(amplitude); }
        let mut full = vec![C64::ZERO; source.blocks[0].nrows()];
        for (&original, &value) in source.free_dofs.iter().zip(&coordinates) { full[original] = value; }
        let mut physical_residual = 0.0_f64;
        for &row in &source.free_dofs {
            poll(cancelled)?;
            let mut sum = [Dd::ZERO; 2];
            let mut envelope = 0.0;
            for (block, multiplier) in [(&source.blocks[0], C64::ONE), (&source.blocks[2], s), (&source.blocks[1], s*s)] {
                let (columns, values) = block.row(row);
                for (&col, &value) in columns.iter().zip(values) {
                    if self.original_to_free[col] == usize::MAX { continue; }
                    let z = full[col];
                    let a = Dd::from_f64(multiplier.re) * Dd::from_f64(value);
                    let b = Dd::from_f64(multiplier.im) * Dd::from_f64(value);
                    sum[0] = sum[0] + a*Dd::from_f64(z.re) - b*Dd::from_f64(z.im);
                    sum[1] = sum[1] + a*Dd::from_f64(z.im) + b*Dd::from_f64(z.re);
                    // Original coefficient row norm times ||x||_infinity=1.
                    // This stays meaningful for structurally zero coordinates.
                    envelope += value.abs() * multiplier.abs();
                }
            }
            let magnitude = C64::new(sum[0].to_f64(), sum[1].to_f64()).abs();
            let residual = if envelope == 0.0 && magnitude == 0.0 { 0.0 } else { magnitude/envelope };
            if !envelope.is_finite() || !residual.is_finite() || residual > self.limits.physical_residual_tolerance {
                return Err(MechanicalPoleError::Residual { row, observed: residual, maximum: self.limits.physical_residual_tolerance });
            }
            physical_residual = physical_residual.max(residual);
        }
        let mut walls = Vec::with_capacity(source.walls.len());
        for (wall, &offset) in source.walls.iter().zip(&source.offsets) {
            poll(cancelled)?;
            let n = wall.dofs.len();
            let displacement = &full[offset..offset+n];
            let mut relaxation_memory = Vec::with_capacity(wall.arms.len());
            let mut relaxation_extension = Vec::with_capacity(wall.arms.len());
            for (arm, &(coupling, rate)) in wall.arms.iter().enumerate() {
                let start = offset+(2+arm)*n;
                let extension = &full[start..start+n];
                let memory: Vec<C64> = if s.abs() > rate {
                    extension.iter().map(|&e| (e / s).scale(rate)).collect()
                } else { original_memory(displacement, extension, coupling) };
                if memory.iter().any(|&z| !finite(z)) {
                    return Err(MechanicalPoleError::Invalid("unrepresentable original wall memory"));
                }
                relaxation_memory.push(memory);
                relaxation_extension.push(extension.to_vec());
            }
            walls.push(MechanicalPoleWall { displacement_m: displacement.to_vec(),
                normal_velocity_m_s: full[offset+n..offset+2*n].to_vec(),
                relaxation_memory, relaxation_extension });
        }
        poll(cancelled)?;
        Ok(MechanicalResonancePole { rate_per_s: s, coordinates,
            pressure_pa: full[..source.fluid_dofs].to_vec(), walls, physical_residual,
            relative_step: pole.relative_step, iterations: pole.iterations })
    }
}

impl CoupledPencil {
    pub(crate) fn resonance_coordinate_count(&self) -> usize { self.free_dofs.len() }

    pub(crate) fn resonance_pole(&self, rate: C64, seed: &[C64], options: &MechanicalPoleOptions,
        limits: MechanicalPoleLimits, cancelled: &mut impl FnMut() -> bool) -> Result<MechanicalResonancePole, MechanicalPoleError> {
        poll(cancelled)?;
        if seed.len() != self.free_dofs.len() || seed.iter().any(|&v| !finite(v)) || rate == C64::ZERO {
            return Err(MechanicalPoleError::Invalid("nonzero rate and original complete pole seed required"));
        }
        let mut query = Query::new(self, limits)?;
        let seed: Vec<_> = seed.iter().zip(&self.scales).map(|(&v, &scale)| v/C64::from_re(scale)).collect();
        let mut refused = None;
        let result = newton_pole_with_cancel(self.free_dofs.len(), rate, &seed, options,
            |s, cancel| query.evaluate(s, cancel).map_err(|error| {
                refused = Some(error); NonlinearEigError::FunctionRefused { at: s }
            }), &mut *cancelled);
        let pole = result.map_err(|error| refused.unwrap_or_else(|| match error {
            NonlinearEigError::Cancelled => MechanicalPoleError::Cancelled,
            error => MechanicalPoleError::Pole(error),
        }))?;
        query.observe(pole, cancelled)
    }

    pub(crate) fn resonance_spectrum(&self, region: MechanicalResonanceRegion, options: &MechanicalSpectrumOptions,
        limits: MechanicalPoleLimits, cancelled: &mut impl FnMut() -> bool) -> Result<MechanicalResonanceSpectrum, MechanicalPoleError> {
        poll(cancelled)?;
        if region.contains(C64::ZERO) {
            return Err(MechanicalPoleError::Invalid("pole contour must exclude the zero-rate pressure/gauge subspace"));
        }
        let mut query = Query::new(self, limits)?;
        let mut refused = None;
        let result = beyn_with_cancel(self.free_dofs.len(), region, options,
            |s, cancel| query.evaluate(s, cancel).map_err(|error| {
                refused = Some(error); NonlinearEigError::FunctionRefused { at: s }
            }), &mut *cancelled);
        let BeynReport { poles, unrefined_values, argument_integral, coarse_argument_integral,
            argument_count, rank, nodes, moment_relative_change, minimum_reciprocal_condition } =
            result.map_err(|error| refused.unwrap_or_else(|| match error {
                BeynError::Operator(NonlinearEigError::Cancelled) => MechanicalPoleError::Cancelled,
                error => MechanicalPoleError::Contour(error),
            }))?;
        let poles = poles.into_iter().map(|pole| query.observe(pole, cancelled)).collect::<Result<_, _>>()?;
        Ok(MechanicalResonanceSpectrum { poles, unrefined_rates_per_s: unrefined_values,
            argument_integrals: [argument_integral, coarse_argument_integral], argument_count,
            rank, nodes, moment_relative_change, minimum_reciprocal_condition })
    }
}

#[cfg(test)]
mod tests;
