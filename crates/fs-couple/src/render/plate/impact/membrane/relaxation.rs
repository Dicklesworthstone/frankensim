//! Uniform proportional film memory on the original condensed strain range.
//! QR is an exact coordinate change, not mode selection: all triangular rows
//! survive, including zero rows. No inverse, eigenvalue cutoff or new time law.
use super::super::{invalid, relaxation::InitialMemory};
use super::{ImpactError, MembranePotential};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::det;
use fs_phs::{MappedRelaxationBranch, RelaxationStrainMap};
use std::{rc::Rc, sync::Arc};

/// One uniform proportional plane-stress spectrum on the actual equilibrium
/// film. A ratio multiplies the complete tensor, including every directional
/// and Poisson cross term. Installed tension and flexural loss are separate.
#[derive(Clone, Debug)]
pub struct MembraneRelaxationSpectrum {
    /// (excess membrane modulus / equilibrium modulus, relaxation time [s]).
    pub branches: Vec<(f64, f64)>,
    /// Original law's declared angular-frequency band, including DC.
    pub omega_band_rad_s: (f64, f64),
}

/// Cold work and complete constitutive state limits. No partial spectrum or
/// modal-pair field is substituted when these budgets are exhausted.
#[derive(Clone, Copy, Debug)]
pub struct MembraneRelaxationBudget {
    pub maximum_source_arms: usize,
    /// Positive arms times ALL symmetric modal pairs; native ceiling is 256.
    pub maximum_memory_states: usize,
    pub maximum_setup_terms: usize,
}

#[derive(Clone, Debug)]
struct EnergyStrains {
    pairs: Vec<(usize, usize)>,
    /// Full upper-triangular R, row major. Energy is |R*(q_i*q_j)|^2/2.
    rows: Vec<f64>,
    modal_limits: Vec<f64>,
    strain_limits: Vec<f64>,
}

#[derive(Clone, Debug)]
pub(in crate::render::plate::impact) struct FilmMemory {
    spectrum: MembraneRelaxationSpectrum,
    strain: Arc<EnergyStrains>,
    pub(in crate::render::plate::impact) initial: InitialMemory,
    states: usize,
    /// Includes original nonlinear equilibrium storage AND every excess arm.
    stiffness_bound: f64,
}

fn finite(x: f64) -> Result<f64, ImpactError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("film material projection is unrepresentable"))
    }
}

impl MembranePotential {
    /// Bind uniform proportional Maxwell stretching before constructing a host.
    /// Every subsequent ImpactSystem automatically retains the same law and
    /// every history; cloning the film cannot silently discard its damping.
    ///
    /// The original condensed strain columns satisfy in-plane equilibrium. A
    /// COMMON scalar spectrum over all facets/tensor directions preserves that
    /// subspace. Histories start either at zero or at the current reduced strain;
    /// there is no heterogeneous law or independently prescribed in-plane strain.
    /// Thin-film geometry, density, equilibrium elasticity and prestress remain
    /// frozen. No relaxing equilibrium spring or modal decay fit is added.
    ///
    /// The supplied absolute modal displacement bounds define a complete box.
    /// Its exact quadratic maps bound each strain, gradient and stress Hessian.
    /// Native history coordinates must stay inside the corresponding strain
    /// envelope; the host checks this at every accepted endpoint. Its original
    /// slope, contact, energy, cancellation and nonlinear clock gates also apply.
    /// Relaxation rates obey the existing midpoint non-ringing limit dt/tau<=2.
    /// This is not a claim that arbitrary time steps resolve material loss.
    ///
    /// Factor rows are cached once using original fs-la Cholesky/Householder QR.
    /// Per positive source arm only n(n+1)/2 histories are needed, independent of
    /// mesh refinement. All rows survive; there is no small-rank truncation.
    pub fn with_proportional_relaxation(
        mut self,
        spectrum: MembraneRelaxationSpectrum,
        maximum_modal_displacements_m_sqrt_kg: Vec<f64>,
        initial: InitialMemory,
        budget: MembraneRelaxationBudget,
    ) -> Result<Self, ImpactError> {
        let n = self.reduction.mode_count();
        let pairs = self.reduction.strain_pairs();
        let count = pairs.len();
        let facets = self.reduction.strain_facets().len();
        let physical_rows = facets
            .checked_mul(3)
            .ok_or_else(|| invalid("film strain row count overflow"))?;
        let rows = physical_rows.max(count);
        let square = count
            .checked_mul(count)
            .ok_or_else(|| invalid("film strain pair count overflow"))?;
        let work = square
            .checked_mul(facets)
            .and_then(|v| v.checked_mul(12))
            .and_then(|v| v.checked_add(square.checked_mul(rows)?.checked_mul(4)?))
            .and_then(|v| v.checked_add(square.checked_mul(count)?))
            .and_then(|v| v.checked_add(facets.checked_mul(count)?.checked_mul(9)?))
            .and_then(|v| v.checked_add(facets.checked_mul(27)?));
        if self.memory.is_some()
            || spectrum.branches.len() > budget.maximum_source_arms
            || spectrum.branches.len() > 256
            || budget.maximum_memory_states > 256
            || maximum_modal_displacements_m_sqrt_kg.len() != n
            || maximum_modal_displacements_m_sqrt_kg
                .iter()
                .any(|x| !x.is_finite() || *x <= 0.0)
            || spectrum.omega_band_rad_s.0 != 0.0
            || !spectrum.omega_band_rad_s.1.is_finite()
            || spectrum.omega_band_rad_s.1 <= 0.0
            || work.is_none_or(|v| v > budget.maximum_setup_terms)
        {
            return Err(invalid(
                "film memory needs one bounded uniform source, positive modal limits and a DC frequency band",
            ));
        }
        let law = GeneralizedMaxwell::new(1.0, spectrum.branches.clone())
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        let positive = law.terms.iter().filter(|(ratio, _)| *ratio > 0.0).count();
        let states = count
            .checked_mul(positive)
            .filter(|v| *v <= budget.maximum_memory_states)
            .ok_or_else(|| {
                invalid("complete film strain history exceeds its declared state budget")
            })?;
        let gram = self
            .reduction
            .stretching_strain_gram(budget.maximum_setup_terms)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        let mut weighted = vec![
            0.0;
            rows.checked_mul(count)
                .ok_or_else(|| invalid("film strain factor size overflow"))?
        ];
        for (facet, (metric, columns)) in self.reduction.strain_facets().enumerate() {
            let local = fs_la::factor::cholesky(metric, 3)
                .map_err(|_| invalid("original film strain metric must be positive definite"))?;
            for channel in 0..3 {
                for (pair, column) in columns.iter().enumerate() {
                    weighted[(3 * facet + channel) * count + pair] =
                        finite((channel..3).map(|i| local.l(i, channel) * column[i]).sum())?;
                }
            }
        }
        let factor = fs_la::factor::qr(&weighted, rows, count);
        let mut coefficients = vec![0.0; square];
        for i in 0..count {
            for j in i..count {
                coefficients[i * count + j] = finite(factor.r(i, j))?;
            }
        }
        // Independent original facet metric: a finite QR result is not enough.
        // This gate never changes the factor or discards a weak direction.
        for i in 0..count {
            for j in 0..=i {
                let actual = finite(
                    (0..count)
                        .map(|row| coefficients[row * count + i] * coefficients[row * count + j])
                        .sum(),
                )?;
                let scale =
                    det::sqrt(gram[i * count + i].abs()) * det::sqrt(gram[j * count + j].abs());
                if !scale.is_finite()
                    || (actual - gram[i * count + j]).abs()
                        > 4096.0
                            * f64::EPSILON
                            * (rows + count) as f64
                            * scale.max(f64::MIN_POSITIVE)
                {
                    return Err(invalid(
                        "complete film factor does not reconstruct original strain energy",
                    ));
                }
            }
        }
        let mut strain_limits = Vec::with_capacity(count);
        let mut equilibrium_bound = 0.0;
        let mut excess_bound = 0.0;
        for row in coefficients.chunks_exact(count) {
            let mut gradient = vec![0.0; n];
            let mut hessian = vec![0.0; n];
            let mut bound = 0.0;
            for (&(i, j), &coefficient) in pairs.iter().zip(row) {
                let c = coefficient.abs();
                let qi = maximum_modal_displacements_m_sqrt_kg[i];
                let qj = maximum_modal_displacements_m_sqrt_kg[j];
                bound += c * qi * qj;
                gradient[i] += c * qj;
                gradient[j] += c * qi;
                hessian[i] += c;
                hessian[j] += c;
            }
            // Outward allowance covers the bounded positive products/sums.
            let rounding = 1.0 + 64.0 * f64::EPSILON * (count + n + 1) as f64;
            let bound = finite(bound * rounding)?;
            let g2 = finite(gradient.iter().map(|g| g * g).sum::<f64>() * rounding)?;
            let h = finite(hessian.into_iter().fold(0.0_f64, f64::max) * rounding)?;
            strain_limits.push(bound);
            equilibrium_bound = finite(equilibrium_bound + g2 + bound * h)?;
            excess_bound = finite(excess_bound + g2 + 2.0 * bound * h)?;
        }
        let ratio = finite(law.terms.iter().map(|(ratio, _)| ratio).sum())?;
        let stiffness_bound = finite(equilibrium_bound + ratio * excess_bound)?;
        self.memory = Some(Arc::new(FilmMemory {
            spectrum,
            initial,
            states,
            stiffness_bound,
            strain: Arc::new(EnergyStrains {
                pairs: pairs.to_vec(),
                rows: coefficients,
                modal_limits: maximum_modal_displacements_m_sqrt_kg,
                strain_limits,
            }),
        }));
        Ok(self)
    }

    /// Original scalar ratios, individual times and applicability band.
    pub fn relaxation_spectrum(&self) -> Option<&MembraneRelaxationSpectrum> {
        self.memory.as_ref().map(|m| &m.spectrum)
    }
    pub fn relaxation_state_count(&self) -> usize {
        self.memory.as_ref().map_or(0, |m| m.states)
    }
    pub(in crate::render::plate::impact) fn memory_configuration(&self) -> Option<&FilmMemory> {
        self.memory.as_deref()
    }
    pub(in crate::render::plate::impact) fn memory_branches(
        &self,
        offset: usize,
        base_dim: usize,
    ) -> Result<Vec<MappedRelaxationBranch>, ImpactError> {
        let Some(memory) = &self.memory else {
            return Ok(vec![]);
        };
        if offset
            .checked_add(self.reduction.mode_count())
            .and_then(|n| n.checked_mul(2))
            .is_none_or(|n| n > base_dim)
        {
            return Err(invalid("film strain map lost its mechanical prefix"));
        }
        let mut branches = Vec::with_capacity(memory.states);
        for &(ratio, tau) in &memory.spectrum.branches {
            if ratio == 0.0 {
                continue;
            }
            for row in 0..memory.strain.pairs.len() {
                branches.push(MappedRelaxationBranch {
                    strain: Rc::new(FilmStrain {
                        data: Arc::clone(&memory.strain),
                        row,
                        offset,
                        base_dim,
                    }),
                    stiffness: ratio,
                    relaxation_time_s: tau,
                });
            }
        }
        Ok(branches)
    }
    pub(in crate::render::plate::impact) fn admit_memory(
        &self,
        offset: usize,
        state: &[f64],
        memory: &[f64],
    ) -> Result<(), ImpactError> {
        let Some(config) = &self.memory else {
            return Ok(());
        };
        if memory.len() != config.states || memory.iter().any(|z| !z.is_finite()) {
            return Err(invalid("complete finite film material history required"));
        }
        for (i, &limit) in config.strain.modal_limits.iter().enumerate() {
            let q = state
                .get(2 * (offset + i))
                .ok_or_else(|| invalid("film material state prefix is incomplete"))?;
            if !q.is_finite() || q.abs() > limit {
                return Err(invalid(
                    "film motion exceeds its material displacement domain",
                ));
            }
        }
        let mut at = 0;
        for &(ratio, _) in &config.spectrum.branches {
            if ratio == 0.0 {
                continue;
            }
            let root = det::sqrt(ratio);
            for &limit in &config.strain.strain_limits {
                let bound = finite(root * limit)?;
                if memory[at].abs() > bound {
                    return Err(invalid(
                        "film material history exceeds its admitted strain domain",
                    ));
                }
                at += 1;
            }
        }
        Ok(())
    }
}

impl FilmMemory {
    pub(in crate::render::plate::impact) fn stiffness_bound(&self) -> f64 {
        self.stiffness_bound
    }
    pub(in crate::render::plate::impact) fn admit_domain(
        &self,
        omega: f64,
        dt: f64,
    ) -> Result<(), ImpactError> {
        if !omega.is_finite()
            || omega > self.spectrum.omega_band_rad_s.1
            || self
                .spectrum
                .branches
                .iter()
                .any(|&(ratio, tau)| ratio > 0.0 && (!(1.0 / tau).is_finite() || dt / tau > 2.0))
        {
            return Err(invalid(
                "film strain spectrum or relaxation is unresolved by the complete mechanical clock",
            ));
        }
        Ok(())
    }
}

#[derive(Debug)]
struct FilmStrain {
    data: Arc<EnergyStrains>,
    row: usize,
    offset: usize,
    base_dim: usize,
}
impl FilmStrain {
    fn valid(&self, x: &[f64], out: usize) -> bool {
        x.len() == self.base_dim && out == self.base_dim && x.iter().all(|v| v.is_finite())
    }
    fn row(&self) -> &[f64] {
        let count = self.data.pairs.len();
        &self.data.rows[self.row * count..(self.row + 1) * count]
    }
}
impl RelaxationStrainMap for FilmStrain {
    fn state_dimension(&self) -> usize {
        self.base_dim
    }
    fn is_quadratic_strain(&self) -> bool {
        true
    }
    fn strain(&self, x: &[f64]) -> f64 {
        if !self.valid(x, x.len()) {
            return f64::NAN;
        }
        self.data
            .pairs
            .iter()
            .zip(self.row())
            .map(|(&(i, j), c)| c * x[2 * (self.offset + i)] * x[2 * (self.offset + j)])
            .sum()
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if !self.valid(x, out.len()) {
            out.fill(f64::NAN);
            return;
        }
        out.fill(0.0);
        for (&(i, j), &c) in self.data.pairs.iter().zip(self.row()) {
            out[2 * (self.offset + i)] += c * x[2 * (self.offset + j)];
            out[2 * (self.offset + j)] += c * x[2 * (self.offset + i)];
        }
    }
    fn hessian_vector(&self, x: &[f64], direction: &[f64], out: &mut [f64]) {
        if !self.valid(x, out.len()) || !self.valid(direction, out.len()) {
            out.fill(f64::NAN);
            return;
        }
        self.gradient(direction, out);
    }
}
