//! Uniform tensor membrane elasticity for sampled von Karman plates.

use fs_material::sheet::SheetMatrix;
use fs_math::det;

use crate::{NlModalError, SampledPlateMode, SineMode, VkModel, sampled_plate_couplings};

/// Physical inputs for a sampled plate with a uniform elastic membrane.
/// Bending stiffness, supports and prestress are already represented by the
/// supplied displacement shapes and their linear frequencies.
#[derive(Debug, Clone, Copy)]
pub struct VkSampledPlateParams {
    /// Rectangle length along x [m].
    pub lx: f64,
    /// Rectangle length along y [m].
    pub ly: f64,
    /// Structural mass per reference area [kg/m²].
    pub areal_density: f64,
    /// Membrane stiffness [N/m] mapping `[εxx, εyy, 2 εxy]` to
    /// `[Nxx, Nyy, Nxy]`, in the same axes as the sampled grid.
    /// For a homogeneous plate this is `12 D / h²`, where `D` is its
    /// bending stiffness matrix. Use the actual oriented section tensor.
    pub membrane: SheetMatrix,
}

/// Build sampled plate storage with the complete plane-stress membrane law.
///
/// The retained unit-normalized Airy fields produce stress resultants
/// `b_j = [Ψ_j,yy, Ψ_j,xx, -Ψ_j,xy]`. Their complementary matrix is
/// `C_jk = ∫ b_jᵀ A⁻¹ b_k dA`, where `A` is the supplied membrane stiffness.
/// Geometric compatibility then gives `Vnl = 1/8 eᵀ C⁻¹ e`, with
/// `e_j = qᵀ E_j q` from the sampled displacement brackets. Rotated
/// anisotropy couples distinct Airy modes through normal/shear compliance;
/// dropping those off-diagonal terms changes the elastic energy.
///
/// Exact sine/cosine integrals assemble `C` without another spatial grid.
/// Its Cholesky factor transforms the coupling matrices once at construction,
/// preserving positive sum-of-squares storage and the existing allocation-free
/// effort/tangent evaluation. No dense Airy solve occurs during time stepping.
/// Isotropic elasticity recovers `C_jj = ξ_j⁴/(E h)` independently of Poisson
/// ratio. The sampled displacement Hessians and the retained Airy basis still
/// require their own spatial and modal convergence checks.
///
/// # Errors
/// Refuses invalid sampled data, nonfinite arithmetic, or an Airy compliance
/// matrix that cannot be resolved as positive definite in floating point.
pub fn von_karman_sampled_plate_with_membrane(
    params: &VkSampledPlateParams,
    disp: &[SampledPlateMode],
    stress_modes: &[SineMode],
) -> Result<VkModel, NlModalError> {
    let mut model = sampled_plate_couplings(
        params.lx,
        params.ly,
        params.areal_density,
        disp,
        stress_modes,
    )?;
    let compliance = membrane_compliance(params.membrane)?;
    let matrix = airy_compliance(params.lx, params.ly, &compliance, stress_modes);
    let count = stress_modes.len();
    let lower = positive_cholesky(&matrix, count)?;
    // If C=L Lᵀ, eᵀ C⁻¹ e = |L⁻¹ e|². Apply the same linear
    // transformation to every modal coupling matrix, retaining symmetry.
    for j in 0..count {
        let (earlier, rest) = model.storage.channels.split_at_mut(j);
        let current = &mut rest[0];
        for (k, previous) in earlier.iter().enumerate() {
            for (entry, previous_entry) in current.coupling.iter_mut().zip(&previous.coupling) {
                *entry -= lower[j * count + k] * previous_entry;
            }
        }
        for entry in &mut current.coupling {
            *entry /= lower[j * count + j];
        }
        current.coefficient = 0.5;
    }
    model.storage.validate()?;
    Ok(model)
}

fn positive_cholesky(matrix: &[f64], n: usize) -> Result<Vec<f64>, NlModalError> {
    let invalid = || NlModalError::Parameter {
        what: "finite resolvable positive membrane compliance",
    };
    if matrix.iter().any(|v| !v.is_finite()) || (0..n).any(|i| matrix[i * n + i] <= 0.0) {
        return Err(invalid());
    }
    // Diagonal congruence keeps the admission independent of units and
    // avoids products of dimensional entries with very different magnitudes.
    let scales: Vec<_> = (0..n).map(|i| det::sqrt(matrix[i * n + i])).collect();
    let mut lower = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut value = matrix[i * n + j] / scales[i] / scales[j];
            for k in 0..j {
                value -= lower[i * n + k] * lower[j * n + k];
            }
            if i == j {
                if !(value > 0.0 && value.is_finite()) {
                    return Err(invalid());
                }
                lower[i * n + j] = det::sqrt(value);
            } else {
                lower[i * n + j] = value / lower[j * n + j];
            }
        }
    }
    for i in 0..n {
        for j in 0..=i {
            lower[i * n + j] *= scales[i];
        }
    }
    if lower.iter().any(|v| !v.is_finite()) {
        return Err(invalid());
    }
    Ok(lower)
}

fn membrane_compliance(membrane: SheetMatrix) -> Result<[f64; 9], NlModalError> {
    let lower = positive_cholesky(membrane.coefficients(), 3)?;
    let mut inverse_lower = [0.0; 9];
    for column in 0..3 {
        for i in column..3 {
            let mut value = if i == column { 1.0 } else { 0.0 };
            for k in column..i {
                value -= lower[i * 3 + k] * inverse_lower[k * 3 + column];
            }
            inverse_lower[i * 3 + column] = value / lower[i * 3 + i];
        }
    }
    let mut compliance = [0.0; 9];
    for i in 0..3 {
        for j in i..3 {
            let value = (j..3)
                .map(|k| inverse_lower[k * 3 + i] * inverse_lower[k * 3 + j])
                .sum();
            compliance[i * 3 + j] = value;
            compliance[j * 3 + i] = value;
        }
    }
    Ok(compliance)
}

// Integral of the unit-normalized sine m against cosine k. Physical
// side length cancels. Integer arithmetic retains parity zeros and weak
// overlaps for widely separated modes without subtracting two fractions.
fn sine_cosine_overlap(m: usize, k: usize) -> f64 {
    let (m, k) = (m as i128, k as i128);
    if (m + k) % 2 == 0 {
        0.0
    } else {
        (4 * m) as f64 / (core::f64::consts::PI * ((m - k) * (m + k)) as f64)
    }
}

fn airy_compliance(lx: f64, ly: f64, s: &[f64; 9], modes: &[SineMode]) -> Vec<f64> {
    let count = modes.len();
    let waves: Vec<_> = modes
        .iter()
        .map(|mode| {
            let kx = mode.m as f64 * core::f64::consts::PI / lx;
            let ky = mode.n as f64 * core::f64::consts::PI / ly;
            [ky * ky, kx * kx, kx * ky]
        })
        .collect();
    let mut matrix = vec![0.0; count * count];
    for i in 0..count {
        let a = waves[i];
        matrix[i * count + i] =
            a[0] * a[0] * s[0] + a[1] * a[1] * s[4] + a[2] * a[2] * (2.0 * s[1] + s[8]);
        for j in 0..i {
            let b = waves[j];
            let ij = sine_cosine_overlap(modes[i].m, modes[j].m)
                * sine_cosine_overlap(modes[i].n, modes[j].n);
            let ji = sine_cosine_overlap(modes[j].m, modes[i].m)
                * sine_cosine_overlap(modes[j].n, modes[i].n);
            let value =
                (a[0] * s[2] + a[1] * s[5]) * b[2] * ij + (b[0] * s[2] + b[1] * s[5]) * a[2] * ji;
            matrix[i * count + j] = value;
            matrix[j * count + i] = value;
        }
    }
    matrix
}
