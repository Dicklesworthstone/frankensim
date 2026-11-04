//! Nonideal aqueous Maxwell–Stefan diffusion in the stationary-solvent frame.
//!
//! The original virial free energy supplies the solute chemical potentials.
//! Solvent concentration is `c0=1/Vbar0`, because the aqueous owner neglects
//! solute partial volumes. Logarithmic-mean face concentrations define x and c.
//! For each solute `d_i=x_i Delta(mu_i)/(RT L)`; `d_0=-sum_i d_i` is the
//! incompressibility/solvent-pressure reaction, not a gas pressure force. The
//! shared molar-frame friction owner solves these balanced forces. Finally
//! `N_i=J_i-(x_i/x_0)J_0` enforces stationary solvent, hence zero volume flux.
//! The discrete ideal chain rule `c_log Delta(ln c)=Delta(c)` keeps the force
//! finite at an absent endpoint without evaluating `0 * infinity`. Finite
//! conjugate dissipation is `-sum_solutes N_i grad(mu_i)/T`; a singular endpoint
//! diagnostic remains explicitly absent, not a fabricated finite value.
//!
//! Logarithmic-mean edge composition: Cancès, Ehrlacher and Monasse,
//! *Finite volumes for the Stefan–Maxwell cross-diffusion system*,
//! <https://doi.org/10.1093/imanum/drad032>. This constitutive face does not
//! establish the positivity or convergence theorem of their time discretization.
//!
//! This is neutral, isothermal, fixed-solvent transport on the original convex
//! concentration domain. It is not unrestricted liquid flow or an ionic current
//! law. Nonnegative endpoint concentrations are admitted; no log floor or seeded
//! trace supplies an absent species. All binary Maxwell–Stefan coefficients,
//! including solute/solute pairs, are explicit fixed numerical data in m²/s.
//! Their measurement, temperature/composition applicability and material source
//! authority belong to the caller. [`source`] binds state-dependent pairs to
//! original material claims without inferring kinetics from the free energy.

use super::{
    AqueousSolution, MAX_OSMOTIC_SOLUTES, OsmosisError, Sum,
    nonideal::{VirialSolutionModel, log_ratio},
};
use fs_qty::{AmountConcentration, MolarMass};
use fs_thermochem::transport::maxwell_stefan::{
    MaxwellStefanErrorV1, MaxwellStefanMolarSolveV1, MaxwellStefanV1, MaxwellStefanWorkspaceV1,
};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::sync::Arc;

pub mod source;

/// The original shared friction workspace, exposed through its material owner.
pub use fs_thermochem::transport::maxwell_stefan::MaxwellStefanWorkspaceV1 as SolutionDiffusionWorkspace;

/// Original source/algebra refusal; no failure changes the caller's output.
#[derive(Clone, Debug, PartialEq)]
pub enum SolutionDiffusionError {
    Invalid(&'static str),
    Thermodynamics(OsmosisError),
    Friction(MaxwellStefanErrorV1),
    Material(crate::state_point::MaterialStatePointError),
    FrameResidual {
        maximum_gradient_residual_per_m: f64,
    },
    Entropy {
        chemical_w_m3_k: f64,
        friction_w_m3_k: f64,
    },
    Cancelled,
}
impl core::fmt::Display for SolutionDiffusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "solvent-frame solution diffusion: {self:?}")
    }
}
impl std::error::Error for SolutionDiffusionError {}
impl From<OsmosisError> for SolutionDiffusionError {
    fn from(value: OsmosisError) -> Self {
        Self::Thermodynamics(value)
    }
}
impl From<MaxwellStefanErrorV1> for SolutionDiffusionError {
    fn from(value: MaxwellStefanErrorV1) -> Self {
        if value == MaxwellStefanErrorV1::Cancelled {
            Self::Cancelled
        } else {
            Self::Friction(value)
        }
    }
}
impl From<crate::state_point::MaterialStatePointError> for SolutionDiffusionError {
    fn from(value: crate::state_point::MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}
type Result<T> = core::result::Result<T, SolutionDiffusionError>;
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolutionDiffusionError::Invalid(
            "unrepresentable aqueous diffusion arithmetic",
        ))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

// Symmetric logarithmic mean with its continuous equal-state and zero-endpoint
// limits. The original relative-log helper preserves adjacent/subnormal inputs
// and avoids overflowing a ratio between widely separated concentrations.
fn logarithmic_mean_concentration(left: f64, right: f64) -> Result<f64> {
    if !left.is_finite() || !right.is_finite() || left < 0.0 || right < 0.0 {
        return Err(SolutionDiffusionError::Invalid(
            "finite nonnegative face concentrations required",
        ));
    }
    if left == 0.0 || right == 0.0 {
        return Ok(0.0);
    }
    if left == right {
        return Ok(left);
    }
    let (low, high) = (left.min(right), left.max(right));
    let delta = high - low;
    let logarithm = log_ratio(high, low);
    let mean = if delta <= low {
        // Form the O(1) ratio before restoring the concentration scale. Direct
        // delta/logarithm can round just below low for adjacent positive inputs.
        low * ((delta / low) / logarithm)
    } else {
        delta / logarithm
    };
    if !positive(mean) || mean < low || mean > high {
        return Err(SolutionDiffusionError::Invalid(
            "positive logarithmic-mean concentration is unrepresentable",
        ));
    }
    Ok(mean)
}

#[cfg(test)]
mod tests {
    use super::logarithmic_mean_concentration;

    #[test]
    fn g0_logarithmic_faces_keep_equal_adjacent_and_subnormal_limits() {
        for value in [1e-300_f64, 1.0, 100.0, 1.1351679650259807e80, 1e300] {
            assert_eq!(logarithmic_mean_concentration(value, value).unwrap(), value);
            let high = value.next_up();
            let mean = logarithmic_mean_concentration(value, high).unwrap();
            assert!((value..=high).contains(&mean));
            assert_eq!(mean, logarithmic_mean_concentration(high, value).unwrap());
        }
        let low = f64::from_bits(3);
        let high = f64::from_bits(6);
        // 3/ln(2) rounds to four least-subnormal units, not five.
        assert_eq!(
            logarithmic_mean_concentration(low, high).unwrap().to_bits(),
            4
        );
        assert_eq!(logarithmic_mean_concentration(0.0, high).unwrap(), 0.0);
        assert_eq!(logarithmic_mean_concentration(0.0, 0.0).unwrap(), 0.0);
    }

    #[test]
    fn g0_logarithmic_faces_span_extreme_finite_scales_without_ratio_overflow() {
        for (low, high) in [
            (1e-300_f64, 1e300_f64),
            (f64::from_bits(1), f64::MAX),
            (f64::MAX / 2.0, f64::MAX),
        ] {
            let mean = logarithmic_mean_concentration(low, high).unwrap();
            let expected = (high - low) / (high.ln() - low.ln());
            assert!((mean / expected - 1.0).abs() < 2e-13);
            assert!((low..=high).contains(&mean));
            assert_eq!(mean, logarithmic_mean_concentration(high, low).unwrap());
        }
        for invalid in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(logarithmic_mean_concentration(invalid, 1.0).is_err());
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SolutionDiffusionError::Cancelled)
    } else {
        Ok(())
    }
}

/// Original friction residuals and conjugate entropy agreement, all explicit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolutionDiffusionSolve {
    pub friction: MaxwellStefanMolarSolveV1,
    /// Absolute nonnegativity and force/friction agreement allowance [W/(m³ K)].
    pub entropy_tolerance_w_m3_k: f64,
}

/// Frame conversion and physical dissipation at one actual concentration face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolutionDiffusionReport {
    pub maximum_gradient_residual_per_m: f64,
    pub molar_frame_residual_mol_m2_s: f64,
    /// Stationary solvent is an exact constitutive constraint, not net molar flux.
    pub solvent_molar_flux_mol_m2_s: f64,
    /// None when a nonzero endpoint jump reaches an absent species. An
    /// identically absent spectator contributes exact zero. Finite-step
    /// free-energy checks remain the evolution owner's job.
    pub chemical_entropy_production_w_m3_k: Option<f64>,
    /// None at a singular zero-endpoint jump. Otherwise the original shared
    /// pair-friction diagnostic must be representable and agree with the
    /// chemical diagnostic within the explicitly declared budgets.
    pub friction_entropy_production_w_m3_k: Option<f64>,
}

/// Immutable solution thermodynamics and a complete constant binary-friction law.
#[derive(Clone, Debug)]
pub struct SolutionMaxwellStefan {
    thermodynamics: Arc<VirialSolutionModel>,
    species: Box<[SpeciesId]>,
    owner: MaxwellStefanV1,
    /// Canonical indices corresponding to original [solvent, solute...] order.
    original_to_canonical: Box<[usize]>,
    binary_diffusivities_m2_s: Box<[f64]>,
}
impl SolutionMaxwellStefan {
    /// Binary matrix is row-major in [solvent, original thermodynamic solutes]
    /// order. Its diagonal is zero; every other entry is positive and exactly
    /// symmetric. Species names and masses are explicit caller assignments.
    pub fn new(
        thermodynamics: Arc<VirialSolutionModel>,
        solvent: (SpeciesId, MolarMass),
        solutes: Vec<(SpeciesId, MolarMass)>,
        binary_diffusivities_m2_s: Vec<f64>,
    ) -> Result<Self> {
        let count = solutes.len();
        let n = count + 1;
        if count != thermodynamics.species_count()
            || count > MAX_OSMOTIC_SOLUTES
            || binary_diffusivities_m2_s.len() != n * n
        {
            return Err(SolutionDiffusionError::Invalid(
                "complete solvent/solute binary layout required",
            ));
        }
        for i in 0..n {
            for j in 0..n {
                let value = binary_diffusivities_m2_s[i * n + j];
                if (i == j && value != 0.0)
                    || (i != j
                        && (!positive(value) || value != binary_diffusivities_m2_s[j * n + i]))
                {
                    return Err(SolutionDiffusionError::Invalid(
                        "positive symmetric binary pairs with zero diagonal required",
                    ));
                }
            }
        }
        let species: Box<[_]> = solutes.iter().map(|(id, _)| id.clone()).collect();
        let originals: Vec<_> = core::iter::once(solvent).chain(solutes).collect();
        let owner = MaxwellStefanV1::new(originals.clone())?;
        let original_to_canonical: Box<[_]> = originals
            .iter()
            .map(|(id, _)| {
                owner
                    .species()
                    .iter()
                    .position(|candidate| candidate == id)
                    .expect("shared owner retains every admitted original species")
            })
            .collect();
        let mut canonical = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                canonical[original_to_canonical[i] * n + original_to_canonical[j]] =
                    binary_diffusivities_m2_s[i * n + j];
            }
        }
        Ok(Self {
            thermodynamics,
            species,
            owner,
            original_to_canonical,
            binary_diffusivities_m2_s: canonical.into_boxed_slice(),
        })
    }
    /// Original solute order, unchanged by the shared solver's canonicalization.
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
    pub fn thermodynamics(&self) -> &VirialSolutionModel {
        &self.thermodynamics
    }
    pub fn workspace(&self) -> Result<MaxwellStefanWorkspaceV1> {
        Ok(self.owner.workspace()?)
    }
    /// Conservative retained array/model payload, excluding this fixed-size owner.
    pub fn heap_bytes(&self) -> usize {
        let n = self.original_to_canonical.len();
        core::mem::size_of_val(self.binary_diffusivities_m2_s.as_ref())
            + core::mem::size_of_val(self.original_to_canonical.as_ref())
            + 3 * n * core::mem::size_of::<f64>()
            + 2 * n * core::mem::size_of::<SpeciesId>()
            + self
                .species
                .iter()
                .map(|id| id.as_str().len())
                .sum::<usize>()
            + self
                .owner
                .species()
                .iter()
                .map(|id| id.as_str().len())
                .sum::<usize>()
            + core::mem::size_of::<VirialSolutionModel>()
            + self.thermodynamics.heap_bytes()
            + 2 * core::mem::size_of::<usize>()
            + core::mem::align_of::<VirialSolutionModel>()
    }
    /// Check the actual original thermodynamic model, T, solvent and species layout.
    pub fn admit_solution(&self, solution: &AqueousSolution) -> Result<()> {
        if solution.nonideal_model() != Some(self.thermodynamics())
            || solution.solute_moles().len() != self.species.len()
            || solution.temperature_k().to_bits() != self.thermodynamics.temperature_k().to_bits()
            || solution.solvent_molar_volume_m3_mol().to_bits()
                != self.thermodynamics.solvent_molar_volume_m3_mol().to_bits()
        {
            return Err(SolutionDiffusionError::Invalid(
                "diffusion requires the same original solution model and ordering",
            ));
        }
        if solution
            .solute_moles()
            .iter()
            .any(|&n| !n.is_finite() || n < 0.0)
        {
            return Err(SolutionDiffusionError::Invalid(
                "diffusion requires finite nonnegative solute amounts",
            ));
        }
        Ok(())
    }

    /// Solute molar fluxes [mol/(m² s)], positive left to right, in ORIGINAL
    /// solute order. No successful-call allocations; output publishes atomically.
    #[allow(clippy::too_many_arguments)]
    pub fn flux_into(
        &self,
        left: &AqueousSolution,
        right: &AqueousSolution,
        distance_m: f64,
        controls: SolutionDiffusionSolve,
        workspace: &mut MaxwellStefanWorkspaceV1,
        output: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SolutionDiffusionReport> {
        self.flux_with_binary_into(
            left,
            right,
            distance_m,
            controls,
            &self.binary_diffusivities_m2_s,
            workspace,
            output,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn flux_with_binary_into(
        &self,
        left: &AqueousSolution,
        right: &AqueousSolution,
        distance_m: f64,
        controls: SolutionDiffusionSolve,
        binary_diffusivities_m2_s: &[f64],
        workspace: &mut MaxwellStefanWorkspaceV1,
        output: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolutionDiffusionReport> {
        poll(&mut cancelled)?;
        let count = self.species.len();
        let n = count + 1;
        if output.len() != count
            || !positive(distance_m)
            || !positive(controls.entropy_tolerance_w_m3_k)
        {
            return Err(SolutionDiffusionError::Invalid(
                "positive distance/entropy allowance and exact output layout required",
            ));
        }
        self.admit_solution(left)?;
        self.admit_solution(right)?;
        let mut a = [0.0; MAX_OSMOTIC_SOLUTES];
        let mut b = [0.0; MAX_OSMOTIC_SOLUTES];
        let mut face = [0.0; MAX_OSMOTIC_SOLUTES];
        let mut chemical_force = [0.0; MAX_OSMOTIC_SOLUTES];
        left.concentrations_into(&mut a[..count])?;
        right.concentrations_into(&mut b[..count])?;
        for (solution, concentrations) in [(left, &a[..count]), (right, &b[..count])] {
            if solution
                .solute_moles()
                .iter()
                .zip(concentrations)
                .any(|(&amount, &value)| amount > 0.0 && value == 0.0)
            {
                return Err(SolutionDiffusionError::Invalid(
                    "positive inventory concentration is unrepresentable",
                ));
            }
        }
        self.thermodynamics
            .evaluate(&a[..count])
            .map_err(OsmosisError::from)?;
        self.thermodynamics
            .evaluate(&b[..count])
            .map_err(OsmosisError::from)?;
        let mut x = [0.0; MAX_OSMOTIC_SOLUTES + 1];
        let mut drive = [0.0; MAX_OSMOTIC_SOLUTES + 1];
        let solvent = self.original_to_canonical[0];
        let c0 = finite(1.0 / self.thermodynamics.solvent_molar_volume_m3_mol())?;
        if !positive(c0) {
            return Err(SolutionDiffusionError::Invalid(
                "solvent concentration is unrepresentable",
            ));
        }
        x[solvent] = c0;
        let mut total = Sum::default();
        total.add(c0);
        for i in 0..count {
            let c = logarithmic_mean_concentration(a[i], b[i])?;
            face[i] = c;
            x[self.original_to_canonical[i + 1]] = c;
            total.add(c);
        }
        let concentration = total.finish()?;
        for value in &mut x[..n] {
            *value /= concentration;
        }
        let mut balance = Sum::default();
        let mut singular_entropy = false;
        let interactions = self.thermodynamics.interaction_matrix_m3_mol();
        for i in 0..count {
            poll(&mut cancelled)?;
            let mut excess = Sum::default();
            for j in 0..count {
                excess.add(finite(interactions[i * count + j] * (b[j] - a[j]))?);
            }
            let excess = excess.finish()?;
            let slot = self.original_to_canonical[i + 1];
            // Form c_log * Delta(mu/RT) directly. The ideal part is exactly
            // Delta(c), including an absent endpoint; the finite excess part
            // still comes from the original symmetric virial interaction law.
            let mut force = Sum::default();
            force.add(b[i] - a[i]);
            force.add(finite(face[i] * excess)?);
            drive[slot] = finite((force.finish()? / concentration) / distance_m)?;
            balance.add(drive[slot]);
            if a[i] > 0.0 && b[i] > 0.0 {
                chemical_force[i] = finite(log_ratio(b[i], a[i]) + excess)?;
            } else if a[i] != b[i] {
                singular_entropy = true;
            }
        }
        drive[solvent] = -balance.finish()?;
        let evaluated = self.owner.evaluate_molar_driving_forces(
            &x[..n],
            &drive[..n],
            AmountConcentration::new(concentration),
            binary_diffusivities_m2_s,
            controls.friction,
            workspace,
            &mut cancelled,
        )?;
        let friction_entropy = if singular_entropy {
            None
        } else {
            Some(evaluated.drag_entropy_production_w_m3_k().ok_or(
                SolutionDiffusionError::Invalid("pair-friction dissipation is unrepresentable"),
            )?)
        };
        let mut converted = [0.0; MAX_OSMOTIC_SOLUTES + 1];
        let flux = evaluated.molar_fluxes();
        for i in 0..n {
            if i != solvent {
                converted[i] = finite(flux[i] - (x[i] / x[solvent]) * flux[solvent])?;
            }
        }
        let mut maximum = evaluated.maximum_gradient_residual_per_m();
        // The shared solve checks its molar frame. Also check every original row
        // after the actual published solvent-frame conversion, without repair.
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut row = Sum::default();
            let mut scale = drive[i].abs();
            for j in 0..n {
                if i == j {
                    continue;
                }
                let term = finite(
                    ((x[j] * converted[i] - x[i] * converted[j])
                        / binary_diffusivities_m2_s[i * n + j])
                        / concentration,
                )?;
                row.add(term);
                scale += term.abs();
            }
            let residual = finite(row.finish()? + drive[i])?.abs();
            maximum = maximum.max(residual);
            if residual
                > controls.friction.absolute_gradient_tolerance_per_m
                    + controls.friction.relative_tolerance * scale
            {
                return Err(SolutionDiffusionError::FrameResidual {
                    maximum_gradient_residual_per_m: maximum,
                });
            }
        }
        let chemical_entropy = if singular_entropy {
            None
        } else {
            let mut dissipation = Sum::default();
            for i in 0..count {
                dissipation.add(
                    -converted[self.original_to_canonical[i + 1]] * (R * chemical_force[i])
                        / distance_m,
                );
            }
            Some(dissipation.finish()?)
        };
        if let (Some(chemical), Some(friction)) = (chemical_entropy, friction_entropy) {
            let allowance = controls.entropy_tolerance_w_m3_k
                + controls.friction.relative_tolerance * chemical.abs().max(friction.abs());
            if chemical < -controls.entropy_tolerance_w_m3_k
                || (chemical - friction).abs() > allowance
            {
                return Err(SolutionDiffusionError::Entropy {
                    chemical_w_m3_k: chemical,
                    friction_w_m3_k: friction,
                });
            }
        }
        poll(&mut cancelled)?;
        for (i, value) in output.iter_mut().enumerate() {
            *value = converted[self.original_to_canonical[i + 1]];
        }
        Ok(SolutionDiffusionReport {
            maximum_gradient_residual_per_m: maximum,
            molar_frame_residual_mol_m2_s: evaluated.molar_frame_residual_mol_m2_s(),
            solvent_molar_flux_mol_m2_s: 0.0,
            chemical_entropy_production_w_m3_k: chemical_entropy,
            friction_entropy_production_w_m3_k: friction_entropy,
        })
    }
}
