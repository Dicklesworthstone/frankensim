//! Consistent isothermal nonideal solution thermodynamics from one free energy.
//!
//! The finite solvent volume is the concentration basis, with negligible solute
//! partial volumes. For the caller's symmetric interaction matrix `B` [m³/mol],
//! `F_mix/(RT V) = sum(c_i*(ln(c_i/c_ref)-1)) + c^T B c/2`.
//! Consequently `Pi/(RT) = sum(c_i)+c^T B c/2`, `ln(a_w)=-Vbar_w*Pi/(RT)`,
//! `ln(gamma_i)=(B c)_i`, and `d(mu_i/RT)/dc_j=delta_ij/c_i+B_ij`.
//! The conventional osmotic second-virial matrix is **B/2**. No cross term,
//! dissociation factor, interaction coefficient, or source authority is inferred.
//!
//! An explicit concentration box and a preparation-time positive-definiteness
//! check keep this single-phase model on a convex branch. The check admits
//! attractive as well as repulsive interactions; it is not an experimental
//! phase boundary or a guarantee that a truncated virial model fits a solution.
//! The caller must supply coefficients and concentration support for the fixed
//! temperature and concentration basis. This is not a universal electrolyte,
//! Pitzer, association, precipitation, or concentrated-solution volume model.
//!
//! Water-only evaluations allocate nothing and avoid solute logarithms. Chemical
//! potentials, their derivatives, and mixing free energy are separate APIs from
//! the SAME potential. They do not retrofit a nonideal constitutive force into
//! the existing ideal solute/film transport laws. A first consumer can safely
//! use this water activity with fully retained solutes or a vapor-equilibrium
//! boundary; leaky transport needs matching chemical-potential forces.
//!
//! References: Elmoazzen, Elliott and McGann (2009), equations 1--3 and the
//! discussion of nondilute transport,
//! <https://doi.org/10.1016/j.bpj.2008.12.3929>; the free-energy/chemical-potential
//! construction of second-virial mixtures in *Phase Behavior of Ternary Polymer
//! Mixtures in a Common Solvent*, <https://doi.org/10.1021/acsomega.3c02604>.

use super::{MAX_OSMOTIC_SOLUTES, Sum};
use fs_la::factor::{FactorError, cholesky};
use fs_math::det;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Malformed model, unsupported state, unstable declared domain, or arithmetic refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum NonidealSolutionError {
    Invalid(&'static str),
    Shape {
        expected: usize,
        actual: usize,
    },
    /// A concentration exceeds the original model's declared support [mol/m³].
    OutsideDomain {
        species: usize,
        concentration_mol_m3: f64,
        maximum_mol_m3: f64,
    },
    /// The concentration box cannot be certified as a numerically stable convex branch.
    UnstableDomain {
        pivot: usize,
    },
    /// An ideal logarithm or reciprocal concentration is singular at this endpoint.
    ZeroConcentration {
        species: usize,
    },
    Arithmetic(&'static str),
}

impl core::fmt::Display for NonidealSolutionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonideal aqueous thermodynamics: {self:?}")
    }
}
impl std::error::Error for NonidealSolutionError {}

type Result<T> = core::result::Result<T, NonidealSolutionError>;

/// Water thermodynamics and the stiffness against proportional concentration change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonidealSolutionState {
    pub osmotic_pressure_pa: f64,
    pub log_water_activity: f64,
    pub water_activity: f64,
    /// `d Pi / d ln(lambda)` for `c_i -> lambda*c_i` [Pa]. At fixed moles,
    /// `d Pi / d V = -osmotic_stiffness_pa / V`.
    pub osmotic_stiffness_pa: f64,
}

/// Prepared symmetric virial interaction model at fixed temperature and solvent basis.
#[derive(Clone, Debug, PartialEq)]
pub struct VirialSolutionModel {
    temperature_k: f64,
    solvent_molar_volume_m3_mol: f64,
    reference_concentration_mol_m3: f64,
    rt: f64,
    log_reference_concentration: f64,
    maximum_concentrations_mol_m3: Box<[f64]>,
    interaction_matrix_m3_mol: Box<[f64]>,
}

impl VirialSolutionModel {
    /// Prepare an explicitly bounded model. The matrix is row-major, exactly
    /// symmetric, and in the SAME species order as the concentration caps.
    /// Its entries are twice the conventional osmotic second-virial coefficients.
    /// `c_ref` specifies dimensionless solute activities and their reference;
    /// it does not change pressure or water activity.
    pub fn new(
        temperature_k: f64,
        solvent_molar_volume_m3_mol: f64,
        reference_concentration_mol_m3: f64,
        maximum_concentrations_mol_m3: Vec<f64>,
        interaction_matrix_m3_mol: Vec<f64>,
    ) -> Result<Self> {
        let n = maximum_concentrations_mol_m3.len();
        if n > MAX_OSMOTIC_SOLUTES {
            return Err(NonidealSolutionError::Invalid(
                "independent-solute count exceeds 64",
            ));
        }
        shape(n * n, interaction_matrix_m3_mol.len())?;
        if !positive(temperature_k)
            || !positive(solvent_molar_volume_m3_mol)
            || !positive(reference_concentration_mol_m3)
            || maximum_concentrations_mol_m3
                .iter()
                .any(|&value| !positive(value))
            || interaction_matrix_m3_mol
                .iter()
                .any(|value| !value.is_finite())
        {
            return Err(NonidealSolutionError::Invalid(
                "temperature, solvent volume, reference concentration and caps must be positive; interactions finite",
            ));
        }
        let rt = finite(R * temperature_k, "RT is not representable")?;
        for i in 0..n {
            for j in 0..i {
                if interaction_matrix_m3_mol[i * n + j] != interaction_matrix_m3_mol[j * n + i] {
                    return Err(NonidealSolutionError::Invalid(
                        "interaction matrix must be exactly symmetric",
                    ));
                }
            }
        }
        admit_convex_box(&maximum_concentrations_mol_m3, &interaction_matrix_m3_mol)?;
        Ok(Self {
            temperature_k,
            solvent_molar_volume_m3_mol,
            reference_concentration_mol_m3,
            rt,
            log_reference_concentration: det::ln(reference_concentration_mol_m3),
            maximum_concentrations_mol_m3: maximum_concentrations_mol_m3.into_boxed_slice(),
            interaction_matrix_m3_mol: interaction_matrix_m3_mol.into_boxed_slice(),
        })
    }

    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn solvent_molar_volume_m3_mol(&self) -> f64 {
        self.solvent_molar_volume_m3_mol
    }
    pub const fn reference_concentration_mol_m3(&self) -> f64 {
        self.reference_concentration_mol_m3
    }
    pub fn species_count(&self) -> usize {
        self.maximum_concentrations_mol_m3.len()
    }
    pub fn maximum_concentrations_mol_m3(&self) -> &[f64] {
        &self.maximum_concentrations_mol_m3
    }
    pub fn interaction_matrix_m3_mol(&self) -> &[f64] {
        &self.interaction_matrix_m3_mol
    }

    /// Owned boxed array payload only, excluding `size_of::<Self>()`.
    pub fn heap_bytes(&self) -> usize {
        core::mem::size_of_val(self.maximum_concentrations_mol_m3.as_ref())
            + core::mem::size_of_val(self.interaction_matrix_m3_mol.as_ref())
    }

    /// Evaluate solvent pressure/activity and the proportional-dilution tangent.
    /// Repeated evaluations are O(N²), allocate nothing, and use no solute logs.
    /// Exact zero concentrations are supported, including a pure-solvent state.
    pub fn evaluate(&self, concentrations_mol_m3: &[f64]) -> Result<NonidealSolutionState> {
        let mut excess = [0.0; MAX_OSMOTIC_SOLUTES];
        self.interactions(concentrations_mol_m3, &mut excess)?;
        let (total, quadratic) = concentration_sums(concentrations_mol_m3, &excess)?;
        let osmolarity = finite(
            total + 0.5 * quadratic,
            "osmotic concentration is not representable",
        )?;
        let stiffness = finite(total + quadratic, "osmotic stiffness is not representable")?;
        if osmolarity < 0.0 || stiffness < 0.0 || (total > 0.0 && stiffness == 0.0) {
            return Err(NonidealSolutionError::Arithmetic(
                "convex-domain osmotic positivity lost to arithmetic",
            ));
        }
        let osmotic_pressure_pa = finite(
            self.rt * osmolarity,
            "osmotic pressure is not representable",
        )?;
        let osmotic_stiffness_pa = finite(
            self.rt * stiffness,
            "osmotic stiffness is not representable",
        )?;
        let log_water_activity = finite(
            -self.solvent_molar_volume_m3_mol * osmolarity,
            "log water activity is not representable",
        )?;
        let water_activity = det::exp(log_water_activity);
        if !positive(water_activity) || water_activity > 1.0 {
            return Err(NonidealSolutionError::Arithmetic(
                "positive water activity is not representable",
            ));
        }
        Ok(NonidealSolutionState {
            osmotic_pressure_pa,
            log_water_activity,
            water_activity,
            osmotic_stiffness_pa,
        })
    }

    /// Evaluate from finite solvent volume [m³] and independent-solute amounts [mol].
    /// Inventories are borrowed and unchanged; concentration staging uses the
    /// fixed 64-species stack workspace, never an allocation.
    pub fn evaluate_inventories(
        &self,
        solvent_volume_m3: f64,
        solute_moles: &[f64],
    ) -> Result<NonidealSolutionState> {
        let concentrations = self.inventory_concentrations(solvent_volume_m3, solute_moles)?;
        self.evaluate(&concentrations[..self.species_count()])
    }

    /// `ln(gamma_i)=(B c)_i` in the declared molar concentration convention.
    /// This finite excess part is defined even when a solute is absent.
    /// Failure leaves the complete output unchanged.
    pub fn log_activity_coefficients_into(
        &self,
        concentrations_mol_m3: &[f64],
        output: &mut [f64],
    ) -> Result<()> {
        shape(self.species_count(), output.len())?;
        let mut excess = [0.0; MAX_OSMOTIC_SOLUTES];
        self.interactions(concentrations_mol_m3, &mut excess)?;
        output.copy_from_slice(&excess[..self.species_count()]);
        Ok(())
    }

    /// Mixing chemical potentials `RT*(ln(c_i/c_ref)+(B c)_i)` [J/mol].
    /// This excludes caller-dependent standard-state constants. An exact absent
    /// solute has an infinite ideal chemical potential and returns a named
    /// refusal; it is never replaced by a concentration floor.
    pub fn chemical_potentials_into(
        &self,
        concentrations_mol_m3: &[f64],
        output: &mut [f64],
    ) -> Result<()> {
        shape(self.species_count(), output.len())?;
        let mut candidate = [0.0; MAX_OSMOTIC_SOLUTES];
        self.interactions(concentrations_mol_m3, &mut candidate)?;
        for (i, &c) in concentrations_mol_m3.iter().enumerate() {
            if c == 0.0 {
                return Err(NonidealSolutionError::ZeroConcentration { species: i });
            }
            candidate[i] = finite(
                self.rt * (det::ln(c) - self.log_reference_concentration + candidate[i]),
                "solute chemical potential is not representable",
            )?;
        }
        output.copy_from_slice(&candidate[..self.species_count()]);
        Ok(())
    }

    /// Right-minus-left solute chemical potentials from THIS same free energy.
    /// The common concentration/standard-state reference cancels analytically,
    /// before floating-point subtraction. Exact absent species still refuse.
    /// Output remains unchanged on any domain or arithmetic failure.
    pub fn chemical_potential_differences_into(
        &self,
        left_mol_m3: &[f64],
        right_mol_m3: &[f64],
        output: &mut [f64],
    ) -> Result<()> {
        let n = self.species_count();
        shape(n, output.len())?;
        self.admit_concentrations(left_mol_m3)?;
        self.admit_concentrations(right_mol_m3)?;
        let mut candidate = [0.0; MAX_OSMOTIC_SOLUTES];
        for i in 0..n {
            if left_mol_m3[i] == 0.0 || right_mol_m3[i] == 0.0 {
                return Err(NonidealSolutionError::ZeroConcentration { species: i });
            }
            let mut force = Sum::default();
            force.add(log_ratio(right_mol_m3[i], left_mol_m3[i]));
            for j in 0..n {
                force.add(finite(
                    self.interaction_matrix_m3_mol[i * n + j] * (right_mol_m3[j] - left_mol_m3[j]),
                    "chemical-potential difference is not representable",
                )?);
            }
            candidate[i] = finite(
                self.rt * sum(force)?,
                "chemical-potential difference is not representable",
            )?;
        }
        output.copy_from_slice(&candidate[..n]);
        Ok(())
    }

    /// Fixed-volume mixing free-energy change, with the common linear
    /// `-RT ln(c_ref) sum(c_i)` reference term cancelled BEFORE arithmetic.
    /// This is a closed, species-conserving accounting primitive; individual
    /// open compartments must retain their chemical reference transfer too.
    /// It supplies no caloric energy or mechanical pressure/volume work.
    pub fn reference_cancelled_free_energy_change_j_m3(
        &self,
        before_mol_m3: &[f64],
        after_mol_m3: &[f64],
    ) -> Result<f64> {
        self.admit_concentrations(before_mol_m3)?;
        self.admit_concentrations(after_mol_m3)?;
        let n = self.species_count();
        let mut change = Sum::default();
        for i in 0..n {
            let (a, b) = (before_mol_m3[i], after_mol_m3[i]);
            if a == b {
                continue;
            }
            let ideal = if a == 0.0 {
                b * (det::ln(b) - 1.0)
            } else if b == 0.0 {
                -a * (det::ln(a) - 1.0)
            } else {
                let delta = b - a;
                delta * det::ln(a) + b * log_ratio(b, a) - delta
            };
            change.add(finite(
                ideal,
                "mixing free-energy difference is not representable",
            )?);
            // Symmetry makes (b-a)^T B (b+a)/2 exactly the quadratic change.
            for j in 0..n {
                change.add(finite(
                    (b - a)
                        * self.interaction_matrix_m3_mol[i * n + j]
                        * (0.5 * before_mol_m3[j] + 0.5 * after_mol_m3[j]),
                    "excess free-energy difference is not representable",
                )?);
            }
        }
        finite(
            self.rt * sum(change)?,
            "mixing free-energy difference is not representable",
        )
    }

    /// Row-major `d mu_i/d c_j = RT*(delta_ij/c_i+B_ij)` [J m³/mol²].
    /// Exact zero concentrations refuse the singular derivative without changing
    /// the output. The returned matrix is symmetric on the admitted convex box.
    pub fn chemical_potential_jacobian_into(
        &self,
        concentrations_mol_m3: &[f64],
        output: &mut [f64],
    ) -> Result<()> {
        let n = self.species_count();
        shape(n * n, output.len())?;
        self.admit_concentrations(concentrations_mol_m3)?;
        let mut candidate = [0.0; MAX_OSMOTIC_SOLUTES * MAX_OSMOTIC_SOLUTES];
        for (i, &c) in concentrations_mol_m3.iter().enumerate() {
            if c == 0.0 {
                return Err(NonidealSolutionError::ZeroConcentration { species: i });
            }
            for j in 0..n {
                let value =
                    self.interaction_matrix_m3_mol[i * n + j] + if i == j { 1.0 / c } else { 0.0 };
                candidate[i * n + j] = finite(
                    self.rt * value,
                    "chemical-potential derivative is not representable",
                )?;
            }
        }
        output.copy_from_slice(&candidate[..n * n]);
        Ok(())
    }

    /// Mixing Helmholtz free-energy density [J/m³] in the declared reference.
    /// Zero-concentration `c*ln(c)` terms take their exact continuous value zero.
    /// This supplies no sensible heat, latent heat, or absolute caloric reference.
    pub fn mixing_free_energy_density_j_m3(&self, concentrations_mol_m3: &[f64]) -> Result<f64> {
        let mut excess = [0.0; MAX_OSMOTIC_SOLUTES];
        self.interactions(concentrations_mol_m3, &mut excess)?;
        let (_, quadratic) = concentration_sums(concentrations_mol_m3, &excess)?;
        let mut ideal = Sum::default();
        for &c in concentrations_mol_m3 {
            if c > 0.0 {
                ideal.add(finite(
                    c * (det::ln(c) - self.log_reference_concentration - 1.0),
                    "mixing free energy is not representable",
                )?);
            }
        }
        finite(
            self.rt * (sum(ideal)? + 0.5 * quadratic),
            "mixing free energy is not representable",
        )
    }

    fn admit_concentrations(&self, concentrations: &[f64]) -> Result<()> {
        shape(self.species_count(), concentrations.len())?;
        for (species, (&c, &maximum)) in concentrations
            .iter()
            .zip(&self.maximum_concentrations_mol_m3)
            .enumerate()
        {
            if !c.is_finite() || c < 0.0 {
                return Err(NonidealSolutionError::Invalid(
                    "concentrations must be finite and nonnegative mol/m³",
                ));
            }
            if c > maximum {
                return Err(NonidealSolutionError::OutsideDomain {
                    species,
                    concentration_mol_m3: c,
                    maximum_mol_m3: maximum,
                });
            }
        }
        Ok(())
    }

    fn interactions(
        &self,
        concentrations: &[f64],
        output: &mut [f64; MAX_OSMOTIC_SOLUTES],
    ) -> Result<()> {
        self.admit_concentrations(concentrations)?;
        let n = self.species_count();
        for (i, row) in self
            .interaction_matrix_m3_mol
            .chunks_exact(n.max(1))
            .enumerate()
        {
            let mut product = Sum::default();
            for (&b, &c) in row.iter().zip(concentrations) {
                product.add(finite(
                    b * c,
                    "excess chemical potential is not representable",
                )?);
            }
            output[i] = sum(product)?;
        }
        Ok(())
    }

    fn inventory_concentrations(
        &self,
        volume: f64,
        amounts: &[f64],
    ) -> Result<[f64; MAX_OSMOTIC_SOLUTES]> {
        shape(self.species_count(), amounts.len())?;
        if !positive(volume) {
            return Err(NonidealSolutionError::Invalid(
                "solvent volume must be finite and strictly positive m³",
            ));
        }
        let mut concentrations = [0.0; MAX_OSMOTIC_SOLUTES];
        for (i, &amount) in amounts.iter().enumerate() {
            if !amount.is_finite() || amount < 0.0 {
                return Err(NonidealSolutionError::Invalid(
                    "solute amounts must be finite and nonnegative mol",
                ));
            }
            concentrations[i] =
                finite(amount / volume, "molar concentration is not representable")?;
        }
        Ok(concentrations)
    }
}

// Relative logarithm avoids cancellation of common concentration references and
// preserves adjacent positive concentrations. The convergent atanh identity is
// used only for ratios in [1/2,2]; otherwise the original deterministic logs
// have no near-equal subtraction. No concentration floor or platform libm.
fn log_ratio(numerator: f64, denominator: f64) -> f64 {
    if numerator == denominator {
        return 0.0;
    }
    let relative = (numerator - denominator) / denominator;
    if (-0.5..=1.0).contains(&relative) {
        // Form the dimensionless relative change before any half scaling:
        // averaging subnormal concentrations would round away their ratio.
        let s = relative / (2.0 + relative);
        let square = s * s;
        let mut term = s;
        let mut result = s;
        for k in 1..24 {
            term *= square;
            result += term / f64::from(2 * k + 1);
        }
        2.0 * result
    } else {
        det::ln(numerator) - det::ln(denominator)
    }
}

// At c_max, H=diag(1/c_max)+B. For every smaller positive concentration,
// H(c)=H(c_max)+diag(1/c-1/c_max), so positive definiteness at the corner
// suffices for the complete box. Congruence by sqrt(c_max), followed by unit
// diagonal scaling, avoids using dimensional pivots as a stability threshold.
fn admit_convex_box(maximum: &[f64], interaction: &[f64]) -> Result<()> {
    let n = maximum.len();
    if n == 0 {
        return Ok(());
    }
    let mut roots = [0.0; MAX_OSMOTIC_SOLUTES];
    let mut diagonal_roots = [0.0; MAX_OSMOTIC_SOLUTES];
    let minimum_squared_pivot = 64.0 * f64::EPSILON * n as f64;
    for i in 0..n {
        roots[i] = det::sqrt(maximum[i]);
        let interaction_diagonal = maximum[i] * interaction[i * n + i];
        let diagonal = finite(
            1.0 + interaction_diagonal,
            "convex-domain diagonal is not representable",
        )?;
        if diagonal <= minimum_squared_pivot * interaction_diagonal.abs().max(1.0) {
            return Err(NonidealSolutionError::UnstableDomain { pivot: i });
        }
        diagonal_roots[i] = det::sqrt(diagonal);
    }
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        matrix[i * n + i] = 1.0;
        for j in 0..i {
            let entry = finite(
                ((roots[i] * interaction[i * n + j]) * roots[j] / diagonal_roots[i])
                    / diagonal_roots[j],
                "convex-domain cross interaction is not representable",
            )?;
            matrix[i * n + j] = entry;
            matrix[j * n + i] = entry;
        }
    }
    let factor = cholesky(&matrix, n).map_err(|error| NonidealSolutionError::UnstableDomain {
        pivot: match error {
            FactorError::NotSpd { index } | FactorError::Singular { index } => index,
        },
    })?;
    // A strictly positive rounded pivot alone cannot distinguish an unresolved
    // near-spinodal model from a stable one. This is a numerical admission
    // margin, not a change to any supplied coefficient or concentration cap.
    for i in 0..n {
        let pivot = factor.l(i, i);
        if !pivot.is_finite() || pivot * pivot <= minimum_squared_pivot {
            return Err(NonidealSolutionError::UnstableDomain { pivot: i });
        }
    }
    Ok(())
}

fn concentration_sums(concentrations: &[f64], excess: &[f64]) -> Result<(f64, f64)> {
    let mut total = Sum::default();
    let mut quadratic = Sum::default();
    for (&c, &mu_excess) in concentrations.iter().zip(excess) {
        total.add(c);
        quadratic.add(finite(
            c * mu_excess,
            "quadratic interaction is not representable",
        )?);
    }
    Ok((sum(total)?, sum(quadratic)?))
}
fn sum(value: Sum) -> Result<f64> {
    value
        .finish()
        .map_err(|_| NonidealSolutionError::Arithmetic("thermodynamic sum is not representable"))
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: f64, reason: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(NonidealSolutionError::Arithmetic(reason))
    }
}
fn shape(expected: usize, actual: usize) -> Result<()> {
    if expected == actual {
        Ok(())
    } else {
        Err(NonidealSolutionError::Shape { expected, actual })
    }
}

#[cfg(test)]
mod tests;
