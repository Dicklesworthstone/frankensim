//! Multicomponent liquid regular solutions with one excess Gibbs-energy owner.
//!
//! `g_mix/(RT) = sum_i x_i ln(x_i) + sum_{i<j} chi_ij(T) x_i x_j`,
//! where every pair is explicit and `chi_ij = A_ij + B_ij/T`. Activities,
//! calorics, simplex curvature and Maxwell-Stefan forces derive from this SAME
//! potential. This is the symmetric-pair specialization of the multicomponent
//! Margules free energy documented by Cantera:
//! <https://cantera.org/stable/cxx/d7/dfe/classCantera_1_1MargulesVPSSTP.html>.
//!
//! The caller qualifies pair fits, pure-liquid references and common pressure.
//! This is a homogeneous constitutive law, including unstable states. The
//! separate equilibrium module searches explicitly seeded phase candidates;
//! it does not infer phase count or a pair interaction from material names.
//! The law excludes electrolytes, unequal-size/polymer entropy,
//! pressure/thermal diffusion, interfacial energy and evolution. Negative
//! curvature is a demixing direction, never a negative molecular diffusivity.

use super::ln_complement;
use crate::transport::maxwell_stefan::{
    MAX_MAXWELL_STEFAN_SPECIES_V1, MaxwellStefanErrorV1, MaxwellStefanEvaluationV1,
    MaxwellStefanMolarEvaluationV1, MaxwellStefanMolarSolveV1, MaxwellStefanSolveV1,
    MaxwellStefanV1, MaxwellStefanWorkspaceV1,
};
use crate::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_math::det::{ln, sqrt};
use fs_qty::{AmountConcentration, Temperature};
use std::sync::Arc;

/// Bounded deterministic phase-count and seed discovery.
pub mod discovery;
/// Conserved heat and entropy of explicitly homogeneous liquid inventories.
pub mod caloric;
/// Conserved finite-inventory phase splits at an explicit phase count and seed.
pub mod equilibrium;
/// Full-simplex local stability and finite tangent-plane comparisons.
pub mod stability;

/// Shared finite species bound; evaluation and forces use bounded stack arrays.
pub const MAX_REGULAR_SOLUTION_SPECIES_V1: usize = MAX_MAXWELL_STEFAN_SPECIES_V1;
/// Only common floating-point normalization is admitted, never a trace floor.
pub const MULTICOMPONENT_COMPOSITION_SUM_TOLERANCE_V1: f64 = 64.0 * f64::EPSILON;

#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentSolutionErrorV1 {
    Invalid(&'static str),
    OutsideTemperatureDomain,
    Unrepresentable(&'static str),
    Transport(MaxwellStefanErrorV1),
    Cancelled,
}
impl core::fmt::Display for MulticomponentSolutionErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multicomponent liquid solution: {self:?}")
    }
}
impl std::error::Error for MulticomponentSolutionErrorV1 {}
impl From<MaxwellStefanErrorV1> for MulticomponentSolutionErrorV1 {
    fn from(error: MaxwellStefanErrorV1) -> Self {
        Self::Transport(error)
    }
}
type Result<T> = core::result::Result<T, MulticomponentSolutionErrorV1>;

#[derive(Clone, Debug, PartialEq)]
struct Owner {
    species: Vec<SpeciesId>,
    constant: Vec<f64>,
    inverse_temperature_k: Vec<f64>,
    temperature_domain: [Temperature; 2],
}

/// Immutable canonical species and all original symmetric pair parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentRegularSolutionV1 {
    owner: Arc<Owner>,
}
impl MulticomponentRegularSolutionV1 {
    /// Admit 1..=128 distinct species and two complete row-major n*n matrices.
    /// Matrices refer to the OFFERED species order; both are canonicalized with
    /// the identities. Subsequent fractions/gradients follow `species()`.
    /// Finite exact symmetry and zero diagonals are required. Explicit zero
    /// means an ideal pair; missing interactions are never assumed to be zero.
    pub fn new(
        species: Vec<SpeciesId>,
        chi_constant: Vec<f64>,
        chi_inverse_temperature_k: Vec<f64>,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self> {
        let n = species.len();
        if !(1..=MAX_REGULAR_SOLUTION_SPECIES_V1).contains(&n)
            || !minimum_temperature.value().is_finite()
            || minimum_temperature.value() <= 0.0
            || !maximum_temperature.value().is_finite()
            || maximum_temperature.value() < minimum_temperature.value()
        {
            return Err(MulticomponentSolutionErrorV1::Invalid(
                "species count or positive temperature domain",
            ));
        }
        for matrix in [&chi_constant, &chi_inverse_temperature_k] {
            if matrix.len() != n * n || matrix.iter().any(|value| !value.is_finite()) {
                return Err(MulticomponentSolutionErrorV1::Invalid(
                    "complete finite pair matrix required",
                ));
            }
            for i in 0..n {
                if matrix[i * n + i] != 0.0
                    || (0..i).any(|j| matrix[i * n + j] != matrix[j * n + i])
                {
                    return Err(MulticomponentSolutionErrorV1::Invalid(
                        "symmetric pair matrix with zero diagonal required",
                    ));
                }
            }
        }
        let mut order: Vec<_> = (0..n).collect();
        order.sort_by(|&i, &j| species[i].cmp(&species[j]));
        if order
            .windows(2)
            .any(|pair| species[pair[0]] == species[pair[1]])
        {
            return Err(MulticomponentSolutionErrorV1::Invalid("duplicate species"));
        }
        let reorder = |matrix: &[f64]| {
            order
                .iter()
                .flat_map(|&i| order.iter().map(move |&j| matrix[i * n + j]))
                .collect()
        };
        Ok(Self {
            owner: Arc::new(Owner {
                species: order.iter().map(|&i| species[i].clone()).collect(),
                constant: reorder(&chi_constant),
                inverse_temperature_k: reorder(&chi_inverse_temperature_k),
                temperature_domain: [minimum_temperature, maximum_temperature],
            }),
        })
    }

    pub fn species(&self) -> &[SpeciesId] {
        &self.owner.species
    }
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        self.owner.temperature_domain
    }
    /// Original dimensionless pair constants in canonical row-major order.
    pub fn chi_constants(&self) -> &[f64] {
        &self.owner.constant
    }
    /// Original inverse-temperature coefficients [K] in the same order.
    pub fn chi_inverse_temperature_k(&self) -> &[f64] {
        &self.owner.inverse_temperature_k
    }

    /// Evaluate one homogeneous composition in O(n²), without heap allocation
    /// after construction (the immutable owner is retained by Arc). Each pair
    /// row polls cancellation. Exact absences stay absent and positive traces
    /// survive common normalization, even when the carrier rounds to one.
    pub fn evaluate(
        &self,
        temperature: Temperature,
        mole_fractions: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentSolutionStateV1> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        let n = self.species().len();
        if !t.is_finite()
            || t < self.owner.temperature_domain[0].value()
            || t > self.owner.temperature_domain[1].value()
        {
            return Err(MulticomponentSolutionErrorV1::OutsideTemperatureDomain);
        }
        if mole_fractions.len() != n
            || mole_fractions
                .iter()
                .any(|&x| !x.is_finite() || !(0.0..=1.0).contains(&x))
        {
            return Err(MulticomponentSolutionErrorV1::Invalid("mole fractions"));
        }
        let total = sum(mole_fractions.iter().copied());
        if (total - 1.0).abs() > MULTICOMPONENT_COMPOSITION_SUM_TOLERANCE_V1 {
            return Err(MulticomponentSolutionErrorV1::Invalid("composition sum"));
        }
        let mut fractions = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for (out, &x) in fractions.iter_mut().zip(mole_fractions) {
            *out = x / total;
            if x > 0.0 && *out == 0.0 {
                return Err(MulticomponentSolutionErrorV1::Unrepresentable(
                    "normalized trace",
                ));
            }
        }
        let x = &fractions[..n];
        let mut ideal = Sum::default();
        let mut logs = [None; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for i in 0..n {
            if x[i] > 0.0 {
                let logarithm = if x[i] > 0.5 {
                    ln_complement(sum((0..n).filter(|&j| j != i).map(|j| x[j])))
                } else {
                    ln(x[i])
                };
                logs[i] = Some(logarithm);
                ideal.add(x[i] * logarithm);
            }
        }
        let mut excess = Sum::default();
        let mut enthalpy_over_r = Sum::default();
        let mut entropy_over_minus_r = Sum::default();
        let mut rows = [Sum::default(); MAX_REGULAR_SOLUTION_SPECIES_V1];
        let mut enthalpy_rows = [Sum::default(); MAX_REGULAR_SOLUTION_SPECIES_V1];
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in i + 1..n {
                let index = i * n + j;
                let a = self.owner.constant[index];
                let b = self.owner.inverse_temperature_k[index];
                let chi = finite(a + b / t, "pair chi(T)")?;
                // Scale the coefficient by one fraction before the other,
                // preserving representable dilute products when x_i*x_j underflows.
                excess.add((chi * x[i]) * x[j]);
                enthalpy_over_r.add((b * x[i]) * x[j]);
                entropy_over_minus_r.add((a * x[i]) * x[j]);
                rows[i].add(chi * x[j]);
                rows[j].add(chi * x[i]);
                enthalpy_rows[i].add(b * x[j]);
                enthalpy_rows[j].add(b * x[i]);
            }
        }
        let excess = finite(excess.value(), "excess Gibbs energy")?;
        let enthalpy_over_r = finite(enthalpy_over_r.value(), "excess enthalpy")?;
        let mut activity = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        let mut partial_h = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut gamma = rows[i].value() - excess;
            let mut h = enthalpy_rows[i].value() - enthalpy_over_r;
            if x[i] > 0.5 {
                // s_i - g_ex nearly cancels for a dilute carrier. Its equivalent
                // (sum_other x)*s_i - g_ex(other,other) retains O(trace²) terms.
                let others = sum((0..n).filter(|&j| j != i).map(|j| x[j]));
                let mut other_excess = Sum::default();
                let mut other_h = Sum::default();
                for j in 0..n {
                    if j == i {
                        continue;
                    }
                    for k in j + 1..n {
                        if k == i {
                            continue;
                        }
                        let index = j * n + k;
                        let b = self.owner.inverse_temperature_k[index];
                        other_excess.add(((self.owner.constant[index] + b / t) * x[j]) * x[k]);
                        other_h.add((b * x[j]) * x[k]);
                    }
                }
                gamma = others * rows[i].value() - other_excess.value();
                h = others * enthalpy_rows[i].value() - other_h.value();
            }
            activity[i] = finite(gamma, "log activity coefficient")?;
            partial_h[i] = finite(R * h, "partial molar mixing enthalpy")?;
            if let Some(value) = logs[i] {
                logs[i] = Some(finite(value + gamma, "chemical potential")?);
            }
        }
        let state = MulticomponentSolutionStateV1 {
            model: self.clone(),
            temperature,
            fractions,
            log_activity: activity,
            reduced_chemical: logs,
            partial_mixing_enthalpy: partial_h,
            mixing_gibbs_j_per_mol: finite(
                (ideal.value() + excess) * R * t,
                "mixing Gibbs energy",
            )?,
            mixing_enthalpy_j_per_mol: finite(enthalpy_over_r * R, "mixing enthalpy")?,
            mixing_entropy_j_per_mol_k: finite(
                -R * (ideal.value() + entropy_over_minus_r.value()),
                "mixing entropy",
            )?,
        };
        poll(&mut cancelled)?;
        Ok(state)
    }
}

/// Complete homogeneous response retaining the exact original pair-law owner.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentSolutionStateV1 {
    model: MulticomponentRegularSolutionV1,
    temperature: Temperature,
    fractions: [f64; MAX_REGULAR_SOLUTION_SPECIES_V1],
    log_activity: [f64; MAX_REGULAR_SOLUTION_SPECIES_V1],
    reduced_chemical: [Option<f64>; MAX_REGULAR_SOLUTION_SPECIES_V1],
    partial_mixing_enthalpy: [f64; MAX_REGULAR_SOLUTION_SPECIES_V1],
    mixing_gibbs_j_per_mol: f64,
    mixing_enthalpy_j_per_mol: f64,
    mixing_entropy_j_per_mol_k: f64,
}
impl MulticomponentSolutionStateV1 {
    pub const fn model(&self) -> &MulticomponentRegularSolutionV1 {
        &self.model
    }
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub fn mole_fractions(&self) -> &[f64] {
        &self.fractions[..self.model.species().len()]
    }
    pub fn log_activity_coefficients(&self) -> &[f64] {
        &self.log_activity[..self.model.species().len()]
    }
    /// Chemical potentials relative to each pure-liquid standard, divided by RT.
    /// An absent component has no finite chemical potential and remains `None`.
    pub fn reduced_chemical_potentials(&self) -> &[Option<f64>] {
        &self.reduced_chemical[..self.model.species().len()]
    }
    /// Partial molar EXCESS enthalpies [J/mol], including finite dilution limits.
    /// Add original pure-liquid sensible enthalpies for a full diffusive heat flux.
    pub fn partial_molar_mixing_enthalpies_j_per_mol(&self) -> &[f64] {
        &self.partial_mixing_enthalpy[..self.model.species().len()]
    }
    pub const fn mixing_gibbs_j_per_mol(&self) -> f64 {
        self.mixing_gibbs_j_per_mol
    }
    pub const fn mixing_enthalpy_j_per_mol(&self) -> f64 {
        self.mixing_enthalpy_j_per_mol
    }
    pub const fn mixing_entropy_j_per_mol_k(&self) -> f64 {
        self.mixing_entropy_j_per_mol_k
    }

    fn chi(&self, i: usize, j: usize) -> f64 {
        let index = i * self.model.species().len() + j;
        self.model.owner.constant[index]
            + self.model.owner.inverse_temperature_k[index] / self.temperature.value()
    }

    /// Directional Hessian `v^T d²(g_mix/RT) v` on the composition simplex.
    /// `sum(v)=0` is required within a relative roundoff allowance. A negative
    /// result is a local spinodal direction; a positive result tests ONLY this
    /// direction and does not prove global miscibility or rule out metastability.
    /// Directions entering an absent component have a divergent ideal curvature
    /// and refuse instead of assigning a concentration floor. O(n²), no allocation.
    pub fn reduced_curvature(&self, direction: &[f64]) -> Result<f64> {
        validate_tangent(direction, self.mole_fractions().len())?;
        let scale = direction.iter().map(|v| v.abs()).fold(0.0, f64::max);
        if scale == 0.0 {
            return Ok(0.0);
        }
        let n = direction.len();
        let mut ideal = Sum::default();
        let mut excess = Sum::default();
        for i in 0..n {
            if self.fractions[i] == 0.0 {
                if direction[i] != 0.0 {
                    return Err(MulticomponentSolutionErrorV1::Invalid(
                        "curvature through an absent species",
                    ));
                }
                continue;
            }
            let vi = direction[i] / scale;
            // Form the physical entropy norm before scaling the interaction
            // term. v_i≈sqrt(x_i) has finite ideal curvature even for a
            // subnormal x_i, where (v_i/max(v))/x_i would overflow.
            let entropy_direction = direction[i] / sqrt(self.fractions[i]);
            ideal.add(entropy_direction * entropy_direction);
            for j in i + 1..n {
                excess.add((2.0 * vi * self.chi(i, j)) * (direction[j] / scale));
            }
        }
        finite(
            ideal.value() + (excess.value() * scale) * scale,
            "simplex curvature",
        )
    }

    /// Isothermal/isobaric `d_i = x_i grad(mu_i)/(RT)` [1/m]. The full
    /// cross-component thermodynamic factor is applied in O(n²); independent
    /// binary Fick laws cannot represent it. Ideal zero-fraction limits retain
    /// the offered gradient without evaluating a divergent chemical potential.
    /// Antisymmetric pair accumulation preserves Gibbs-Duhem balance to roundoff.
    pub fn driving_forces(
        &self,
        gradients_per_m: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentDrivingForcesV1> {
        poll(&mut cancelled)?;
        let n = self.model.species().len();
        validate_tangent(gradients_per_m, n)?;
        let mut response = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        let mut forces = [Sum::default(); MAX_REGULAR_SOLUTION_SPECIES_V1];
        for i in 0..n {
            poll(&mut cancelled)?;
            response[i] = finite(
                sum((0..n).map(|j| self.chi(i, j) * gradients_per_m[j])),
                "activity gradient",
            )?;
            forces[i].add(gradients_per_m[i]);
        }
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in i + 1..n {
                let transfer = finite(
                    ((response[i] - response[j]) * self.fractions[i]) * self.fractions[j],
                    "thermodynamic force pair",
                )?;
                forces[i].add(transfer);
                forces[j].add(-transfer);
            }
        }
        let mut values = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for i in 0..n {
            values[i] = finite(forces[i].value(), "thermodynamic force")?;
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentDrivingForcesV1 { values, length: n })
    }

    /// Original nonideal liquid forces consumed by the existing mass-frame
    /// Maxwell-Stefan owner. Species identities and canonical order must match;
    /// molecular diffusivities and liquid concentration stay explicit sources.
    #[allow(clippy::too_many_arguments)]
    pub fn mass_diffusion<'a>(
        &self,
        transport: &MaxwellStefanV1,
        gradients_per_m: &[f64],
        concentration: AmountConcentration,
        binary_diffusivities_m2_s: &[f64],
        controls: MaxwellStefanSolveV1,
        workspace: &'a mut MaxwellStefanWorkspaceV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanEvaluationV1<'a>> {
        self.check_transport(transport)?;
        let driving = self.driving_forces(gradients_per_m, &mut cancelled)?;
        Ok(transport.evaluate_driving_forces(
            self.mole_fractions(),
            driving.values(),
            concentration,
            binary_diffusivities_m2_s,
            controls,
            workspace,
            cancelled,
        )?)
    }

    /// Same physical constitutive force with the distinct zero-net-molar frame.
    /// This supplies neither spatial stability nor a phase-equilibrium result.
    #[allow(clippy::too_many_arguments)]
    pub fn molar_diffusion<'a>(
        &self,
        transport: &MaxwellStefanV1,
        gradients_per_m: &[f64],
        concentration: AmountConcentration,
        binary_diffusivities_m2_s: &[f64],
        controls: MaxwellStefanMolarSolveV1,
        workspace: &'a mut MaxwellStefanWorkspaceV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanMolarEvaluationV1<'a>> {
        self.check_transport(transport)?;
        let driving = self.driving_forces(gradients_per_m, &mut cancelled)?;
        Ok(transport.evaluate_molar_driving_forces(
            self.mole_fractions(),
            driving.values(),
            concentration,
            binary_diffusivities_m2_s,
            controls,
            workspace,
            cancelled,
        )?)
    }

    fn check_transport(&self, transport: &MaxwellStefanV1) -> Result<()> {
        if transport.species() != self.model.species() {
            Err(MulticomponentSolutionErrorV1::Invalid(
                "transport species differ from liquid owner",
            ))
        } else {
            Ok(())
        }
    }
}

/// Complete balanced force vector in the original canonical species order.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentDrivingForcesV1 {
    values: [f64; MAX_REGULAR_SOLUTION_SPECIES_V1],
    length: usize,
}
impl MulticomponentDrivingForcesV1 {
    pub fn values(&self) -> &[f64] {
        &self.values[..self.length]
    }
}

fn validate_tangent(values: &[f64], n: usize) -> Result<()> {
    if values.len() != n || values.iter().any(|v| !v.is_finite()) {
        return Err(MulticomponentSolutionErrorV1::Invalid(
            "simplex tangent dimensions or values",
        ));
    }
    let scale = values.iter().map(|v| v.abs()).fold(0.0, f64::max);
    if scale > 0.0 {
        let residual = sum(values.iter().map(|v| v / scale));
        let absolute_sum = sum(values.iter().map(|v| (v / scale).abs()));
        if residual.abs() > 64.0 * f64::EPSILON * absolute_sum {
            return Err(MulticomponentSolutionErrorV1::Invalid(
                "simplex tangent must sum to zero",
            ));
        }
    }
    Ok(())
}
fn finite(value: f64, label: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MulticomponentSolutionErrorV1::Unrepresentable(label))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentSolutionErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
#[derive(Clone, Copy, Default)]
struct Sum {
    main: f64,
    tail: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.main + value;
        self.tail += if self.main.abs() >= value.abs() {
            (self.main - next) + value
        } else {
            (value - next) + self.main
        };
        self.main = next;
    }
    fn value(self) -> f64 {
        self.main + self.tail
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut value = Sum::default();
    for term in values {
        value.add(term);
    }
    value.value()
}

#[cfg(test)]
mod tests;
