//! Neutral ideal-gas Maxwell–Stefan diffusion in the mass-averaged frame.
//!
//! With `z_i = j_i/M_i`, solve
//! `sum_j (x_j z_i - x_i z_j)/D_ij = -c [grad x_i + (x_i-Y_i) grad ln p]`
//! together with `sum_i M_i z_i = 0`. This is a coupled friction solve, not
//! independent Fick laws. It remains regular at exact zero mole fractions;
//! no division by `x_i` enters the linear system. Temperature specifies total
//! molar concentration, not a Soret force. Thermal diffusion, external
//! species forces, ions and nonideal activities are outside this closure.
//!
//! Primary formulation: COMSOL, *Multicomponent Gas Diffusion: Maxwell–Stefan
//! Description*, equations 3-75–3-80:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.chem/chem_ug_chemsptrans.08.257.html>.
//! This L1 module owns constitutive algebra and scratch storage, never evolution.

use crate::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_qty::{MolarMass, Pressure, Temperature};

/// Hard species/work bound for the dense constitutive solve.
pub const MAX_MAXWELL_STEFAN_SPECIES_V1: usize = 128;
/// Fixed operation-tree version of the barycentric friction solve.
pub const MAXWELL_STEFAN_EVALUATOR_VERSION_V1: u32 = 1;

/// Explicit checks on original constitutive equations and the mass frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaxwellStefanSolveV1 {
    /// Relative tolerance on each row's own driving/friction scale and net mass.
    pub relative_tolerance: f64,
    /// Absolute original-equation and gradient-compatibility tolerance [1/m].
    pub absolute_gradient_tolerance_per_m: f64,
    /// Absolute net diffusive mass-flux tolerance [kg/(m² s)].
    pub absolute_mass_flux_tolerance_kg_m2_s: f64,
}
impl Default for MaxwellStefanSolveV1 {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            absolute_gradient_tolerance_per_m: 1e-12,
            absolute_mass_flux_tolerance_kg_m2_s: 1e-16,
        }
    }
}

/// Refusal before a complete borrowed evaluation can be published.
#[derive(Clone, Debug, PartialEq)]
pub enum MaxwellStefanErrorV1 {
    /// Invalid species, dimensions, fractions, coefficients or controls.
    Invalid(&'static str),
    /// Nonzero physical intermediate cannot be represented by binary64.
    Arithmetic(&'static str),
    /// Scaled elimination found a zero or non-finite pivot.
    SingularSystem,
    /// An original equation or the mass-frame constraint missed its tolerance.
    Residual {
        /// Largest original Maxwell–Stefan row discrepancy [1/m].
        maximum_gradient_residual_per_m: f64,
        /// Signed sum of published species mass fluxes [kg/(m² s)].
        mass_frame_residual_kg_m2_s: f64,
    },
    /// Bounded workspace reservation failed.
    AllocationRefused,
    /// Cancellation before complete publication.
    Cancelled,
}
impl core::fmt::Display for MaxwellStefanErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "barycentric Maxwell-Stefan diffusion: {self:?}")
    }
}
impl std::error::Error for MaxwellStefanErrorV1 {}

/// Immutable canonical species and molecular masses for a diffusion closure.
/// Diffusivities remain explicit caller-supplied pair data at each actual state.
#[derive(Clone, Debug, PartialEq)]
pub struct MaxwellStefanV1 {
    species: Vec<SpeciesId>,
    molar_masses: Vec<MolarMass>,
    scaled_masses: Vec<f64>,
}
impl MaxwellStefanV1 {
    /// Canonicalize original identities and admit positive finite molar masses.
    /// Evaluation arrays follow the returned canonical `species()` order.
    pub fn new(mut species: Vec<(SpeciesId, MolarMass)>) -> Result<Self, MaxwellStefanErrorV1> {
        if !(1..=MAX_MAXWELL_STEFAN_SPECIES_V1).contains(&species.len())
            || species.iter().any(|(_, mass)| !positive(mass.value()))
        {
            return Err(MaxwellStefanErrorV1::Invalid("species count or molar mass"));
        }
        species.sort_by(|a, b| a.0.cmp(&b.0));
        if species.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(MaxwellStefanErrorV1::Invalid("duplicate species"));
        }
        let n = species.len();
        let maximum = species.iter().map(|(_, m)| m.value()).fold(0.0, f64::max);
        let mut ids = reserved(n)?;
        let mut masses = reserved(n)?;
        let mut scaled = reserved(n)?;
        for (id, mass) in species {
            let value = mass.value() / maximum;
            if !positive(value) {
                return Err(MaxwellStefanErrorV1::Arithmetic("molar-mass scaling"));
            }
            ids.push(id);
            masses.push(mass);
            scaled.push(value);
        }
        Ok(Self {
            species: ids,
            molar_masses: masses,
            scaled_masses: scaled,
        })
    }
    /// Original species identities in deterministic canonical order.
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
    /// Original molar masses [kg/mol] in canonical species order.
    pub fn molar_masses(&self) -> &[MolarMass] {
        &self.molar_masses
    }
    /// Allocate the complete reusable matrix, solve and output workspace once.
    pub fn workspace(&self) -> Result<MaxwellStefanWorkspaceV1, MaxwellStefanErrorV1> {
        let n = self.species.len();
        Ok(MaxwellStefanWorkspaceV1 {
            matrix: filled(n * n, 0.0)?,
            ratios: filled(n * n, 0.0)?,
            right: filled(n, 0.0)?,
            row_scales: filled(n, 0.0)?,
            driving: filled(n, 0.0)?,
            molar_flux: filled(n, 0.0)?,
            mass_flux: filled(n, 0.0)?,
        })
    }

    /// Evaluate one spatial direction without successful-path allocation.
    ///
    /// `binary_diffusivities_m2_s` is an exactly symmetric row-major `n*n`
    /// matrix with positive finite off-diagonals and zero diagonal. Its source
    /// temperature, pressure and validity are the caller's responsibility.
    /// Fractions may contain exact zeros; neither fractions nor fluxes are
    /// clipped or projected. Every original friction row, including the one
    /// replaced by the barycentric constraint, is checked before return.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate<'a>(
        &self,
        mole_fractions: &[f64],
        gradients_per_m: &[f64],
        temperature: Temperature,
        pressure: Pressure,
        log_pressure_gradient_per_m: f64,
        binary_diffusivities_m2_s: &[f64],
        controls: MaxwellStefanSolveV1,
        workspace: &'a mut MaxwellStefanWorkspaceV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanEvaluationV1<'a>, MaxwellStefanErrorV1> {
        poll(&mut cancelled)?;
        let n = self.species.len();
        if mole_fractions.len() != n
            || gradients_per_m.len() != n
            || binary_diffusivities_m2_s.len() != n * n
            || workspace.right.len() != n
            || !positive(temperature.value())
            || !positive(pressure.value())
            || !log_pressure_gradient_per_m.is_finite()
            || !positive(controls.relative_tolerance)
            || controls.relative_tolerance > 1e-3
            || !positive(controls.absolute_gradient_tolerance_per_m)
            || !positive(controls.absolute_mass_flux_tolerance_kg_m2_s)
            || mole_fractions
                .iter()
                .any(|x| !x.is_finite() || !(0.0..=1.0).contains(x))
            || gradients_per_m.iter().any(|x| !x.is_finite())
        {
            return Err(MaxwellStefanErrorV1::Invalid(
                "state, array dimensions or numerical controls",
            ));
        }
        let fraction_sum = sum(mole_fractions.iter().copied());
        if (fraction_sum - 1.0).abs() > 64.0 * f64::EPSILON * n as f64 {
            return Err(MaxwellStefanErrorV1::Invalid(
                "mole fractions do not close the simplex",
            ));
        }
        let gradient_sum = sum(gradients_per_m.iter().copied());
        let gradient_scale = sum(gradients_per_m.iter().map(|g| g.abs()));
        if !gradient_scale.is_finite()
            || gradient_sum.abs()
                > controls.absolute_gradient_tolerance_per_m
                    + controls.relative_tolerance * gradient_scale
        {
            return Err(MaxwellStefanErrorV1::Invalid(
                "mole-fraction gradients do not sum to zero",
            ));
        }
        let mean_mass = sum(mole_fractions
            .iter()
            .zip(&self.molar_masses)
            .map(|(x, m)| x * m.value()));
        let concentration = (pressure.value() / R) / temperature.value();
        if !positive(mean_mass) || !positive(concentration) {
            return Err(MaxwellStefanErrorV1::Arithmetic(
                "mean molar mass or total concentration",
            ));
        }
        let mut minimum_diffusion = f64::INFINITY;
        for i in 0..n {
            poll(&mut cancelled)?;
            if binary_diffusivities_m2_s[i * n + i] != 0.0 {
                return Err(MaxwellStefanErrorV1::Invalid(
                    "nonzero binary-diffusion diagonal",
                ));
            }
            for j in i + 1..n {
                let value = binary_diffusivities_m2_s[i * n + j];
                if !positive(value) || value != binary_diffusivities_m2_s[j * n + i] {
                    return Err(MaxwellStefanErrorV1::Invalid(
                        "binary diffusion is not positive symmetric",
                    ));
                }
                minimum_diffusion = minimum_diffusion.min(value);
            }
        }
        let mut drive_scale = 0.0_f64;
        let mut gauge = 0;
        for i in 0..n {
            let mass_fraction = (mole_fractions[i] * self.molar_masses[i].value()) / mean_mass;
            let drive = gradients_per_m[i]
                + (mole_fractions[i] - mass_fraction) * log_pressure_gradient_per_m;
            if !drive.is_finite() {
                return Err(MaxwellStefanErrorV1::Arithmetic(
                    "pressure/composition driving force",
                ));
            }
            workspace.driving[i] = drive;
            drive_scale = drive_scale.max(drive.abs());
            if mole_fractions[i] > mole_fractions[gauge] {
                gauge = i;
            }
        }
        workspace.mass_flux.fill(0.0);
        workspace.molar_flux.fill(0.0);
        if n == 1 || drive_scale == 0.0 {
            if drive_scale > controls.absolute_gradient_tolerance_per_m {
                return Err(MaxwellStefanErrorV1::Invalid(
                    "single-species nonzero diffusion force",
                ));
            }
            poll(&mut cancelled)?;
            return Ok(MaxwellStefanEvaluationV1 {
                mass_fluxes: &workspace.mass_flux,
                molar_fluxes: &workspace.molar_flux,
                maximum_gradient_residual_per_m: drive_scale,
                mass_frame_residual_kg_m2_s: 0.0,
                drag_entropy_production_w_m3_k: Some(0.0),
            });
        }
        workspace.matrix.fill(0.0);
        workspace.ratios.fill(0.0);
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let ratio = minimum_diffusion / binary_diffusivities_m2_s[i * n + j];
                if !positive(ratio) {
                    return Err(MaxwellStefanErrorV1::Arithmetic("binary-diffusivity ratio"));
                }
                workspace.ratios[i * n + j] = ratio;
                workspace.matrix[i * n + i] += mole_fractions[j] * ratio;
                workspace.matrix[i * n + j] = -mole_fractions[i] * ratio;
            }
            workspace.right[i] = -workspace.driving[i] / drive_scale;
        }
        workspace.matrix[gauge * n..(gauge + 1) * n].copy_from_slice(&self.scaled_masses);
        workspace.right[gauge] = 0.0;
        solve(
            &mut workspace.matrix,
            &mut workspace.right,
            &mut workspace.row_scales,
            n,
            &mut cancelled,
        )?;
        let flux_scale = product3(concentration, minimum_diffusion, drive_scale)?;
        for i in 0..n {
            let molar = workspace.right[i] * flux_scale;
            let mass = molar * self.molar_masses[i].value();
            if !molar.is_finite()
                || !mass.is_finite()
                || (workspace.right[i] != 0.0 && (molar == 0.0 || mass == 0.0))
            {
                return Err(MaxwellStefanErrorV1::Arithmetic("species diffusion flux"));
            }
            workspace.mass_flux[i] = mass;
            // All residuals concern the actually published mass flux, including
            // the final molar/mass conversion's rounding.
            workspace.molar_flux[i] = mass / self.molar_masses[i].value();
        }
        let mass_residual = sum(workspace.mass_flux.iter().copied());
        let mass_scale = sum(workspace.mass_flux.iter().map(|j| j.abs()));
        let mut maximum_residual = 0.0_f64;
        let mut rows_close = true;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut terms = [0.0; MAX_MAXWELL_STEFAN_SPECIES_V1];
            let wi = workspace.molar_flux[i] / flux_scale;
            for j in 0..n {
                if i != j {
                    let wj = workspace.molar_flux[j] / flux_scale;
                    terms[j] = (mole_fractions[j] * wi - mole_fractions[i] * wj)
                        * workspace.ratios[i * n + j]
                        * drive_scale;
                }
            }
            let residual = workspace.driving[i] + sum(terms[..n].iter().copied());
            let row_scale = workspace.driving[i].abs() + sum(terms[..n].iter().map(|t| t.abs()));
            if !residual.is_finite() || !row_scale.is_finite() {
                return Err(MaxwellStefanErrorV1::Arithmetic(
                    "original diffusion residual",
                ));
            }
            maximum_residual = maximum_residual.max(residual.abs());
            rows_close &= residual.abs()
                <= controls.absolute_gradient_tolerance_per_m
                    + controls.relative_tolerance * row_scale;
        }
        if !mass_residual.is_finite() || !mass_scale.is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic("mass-frame residual"));
        }
        if !rows_close
            || mass_residual.abs()
                > controls.absolute_mass_flux_tolerance_kg_m2_s
                    + controls.relative_tolerance * mass_scale
        {
            return Err(MaxwellStefanErrorV1::Residual {
                maximum_gradient_residual_per_m: maximum_residual,
                mass_frame_residual_kg_m2_s: mass_residual,
            });
        }
        let entropy = drag_entropy(
            mole_fractions,
            &workspace.molar_flux,
            binary_diffusivities_m2_s,
            concentration,
        );
        poll(&mut cancelled)?;
        Ok(MaxwellStefanEvaluationV1 {
            mass_fluxes: &workspace.mass_flux,
            molar_fluxes: &workspace.molar_flux,
            maximum_gradient_residual_per_m: maximum_residual,
            mass_frame_residual_kg_m2_s: mass_residual,
            drag_entropy_production_w_m3_k: entropy,
        })
    }
}

/// Reusable bounded constitutive scratch, with no evolving physical state.
#[derive(Clone, Debug)]
pub struct MaxwellStefanWorkspaceV1 {
    matrix: Vec<f64>,
    ratios: Vec<f64>,
    right: Vec<f64>,
    row_scales: Vec<f64>,
    driving: Vec<f64>,
    molar_flux: Vec<f64>,
    mass_flux: Vec<f64>,
}

/// Complete borrowed diffusive fluxes and actual equation residuals.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellStefanEvaluationV1<'a> {
    mass_fluxes: &'a [f64],
    molar_fluxes: &'a [f64],
    maximum_gradient_residual_per_m: f64,
    mass_frame_residual_kg_m2_s: f64,
    drag_entropy_production_w_m3_k: Option<f64>,
}
impl MaxwellStefanEvaluationV1<'_> {
    /// Canonical diffusive mass flux relative to the mass-averaged velocity.
    pub const fn mass_fluxes(&self) -> &[f64] {
        self.mass_fluxes
    }
    /// Canonical molar flux in that same mass-averaged frame [mol/(m² s)].
    /// Its sum is generally nonzero when molecular masses differ.
    pub const fn molar_fluxes(&self) -> &[f64] {
        self.molar_fluxes
    }
    /// Largest absolute residual among all original friction equations [1/m].
    pub const fn maximum_gradient_residual_per_m(&self) -> f64 {
        self.maximum_gradient_residual_per_m
    }
    /// Signed sum of the published mass fluxes [kg/(m² s)].
    pub const fn mass_frame_residual_kg_m2_s(&self) -> f64 {
        self.mass_frame_residual_kg_m2_s
    }
    /// Local nonnegative pair-friction entropy production [W/(m³ K)].
    /// None denotes a nonzero flux at an exact-zero composition boundary or an
    /// unrepresentable diagnostic. This is not a discrete transport inequality.
    pub const fn drag_entropy_production_w_m3_k(&self) -> Option<f64> {
        self.drag_entropy_production_w_m3_k
    }
}

fn solve(
    a: &mut [f64],
    b: &mut [f64],
    scales: &mut [f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MaxwellStefanErrorV1> {
    for i in 0..n {
        scales[i] = a[i * n..(i + 1) * n]
            .iter()
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        if !positive(scales[i]) {
            return Err(MaxwellStefanErrorV1::SingularSystem);
        }
    }
    for k in 0..n {
        poll(cancelled)?;
        let mut pivot = k;
        for i in k + 1..n {
            if a[i * n + k].abs() / scales[i] > a[pivot * n + k].abs() / scales[pivot] {
                pivot = i;
            }
        }
        if a[pivot * n + k] == 0.0 || !a[pivot * n + k].is_finite() {
            return Err(MaxwellStefanErrorV1::SingularSystem);
        }
        if pivot != k {
            for j in 0..n {
                a.swap(k * n + j, pivot * n + j);
            }
            b.swap(k, pivot);
            scales.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.0;
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    for i in (0..n).rev() {
        poll(cancelled)?;
        let rhs = b[i] - sum((i + 1..n).map(|j| a[i * n + j] * b[j]));
        b[i] = rhs / a[i * n + i];
        if !b[i].is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic("diffusion linear solve"));
        }
    }
    Ok(())
}

fn drag_entropy(x: &[f64], z: &[f64], diffusion: &[f64], concentration: f64) -> Option<f64> {
    let n = x.len();
    if x.iter().zip(z).any(|(x, z)| *x == 0.0 && *z != 0.0) {
        return None;
    }
    let mut total = 0.0;
    for i in 0..n {
        if x[i] == 0.0 {
            continue;
        }
        for j in i + 1..n {
            if x[j] == 0.0 {
                continue;
            }
            let velocity = (z[i] / concentration) / x[i] - (z[j] / concentration) / x[j];
            let term =
                (R * concentration * x[i] * x[j] / diffusion[i * n + j]) * velocity * velocity;
            if !term.is_finite() {
                return None;
            }
            total += term;
        }
    }
    total.is_finite().then_some(total)
}
fn product3(a: f64, b: f64, c: f64) -> Result<f64, MaxwellStefanErrorV1> {
    for value in [(a * b) * c, (a * c) * b, (b * c) * a] {
        if positive(value) {
            return Ok(value);
        }
    }
    let value = fs_math::det::exp(fs_math::det::ln(a) + fs_math::det::ln(b) + fs_math::det::ln(c));
    if positive(value) {
        Ok(value)
    } else {
        Err(MaxwellStefanErrorV1::Arithmetic("diffusive flux scale"))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaxwellStefanErrorV1> {
    if cancelled() {
        Err(MaxwellStefanErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn reserved<T>(count: usize) -> Result<Vec<T>, MaxwellStefanErrorV1> {
    let mut v = Vec::new();
    v.try_reserve_exact(count)
        .map_err(|_| MaxwellStefanErrorV1::AllocationRefused)?;
    Ok(v)
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, MaxwellStefanErrorV1> {
    let mut v = reserved(count)?;
    v.resize(count, value);
    Ok(v)
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        let y = value - correction;
        let next = total + y;
        correction = (next - total) - y;
        total = next;
    }
    total
}

#[cfg(test)]
mod tests;
