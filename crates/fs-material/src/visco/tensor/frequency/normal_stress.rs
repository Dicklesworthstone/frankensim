//! Normal-stress-free harmonic material retaining physical transverse shear.

use super::super::{Mandel6, MaterialPlane, TensorProny, ViscoError, solve_complex};
use super::plane_stress::projection::{self, Sum};
use fs_math::c64::C64;

/// Prepared five-component material with only the transverse NORMAL stress zero.
///
/// Ordered engineering strains are `[epsilon_u, epsilon_v, gamma_uv,
/// gamma_un, gamma_vn]`, where `u,v` are the ordered material plane and
/// `n` is its remaining axis. Work-conjugate stresses retain both transverse
/// shear tractions. Every anisotropic coupling among these five coordinates
/// survives the one-normal complex constraint.
///
/// The immutable original law is borrowed, retaining its source frequency band,
/// complete Maxwell/Kelvin tensors and full six-component strain ceiling.
/// Peak phasors use `exp(-i omega t)`. This local constitutive reduction does
/// not select a shell shear profile, thickness integration or spatial basis,
/// and does not assert traction-free transverse shear at an exterior face.
#[derive(Clone, Debug, PartialEq)]
pub struct NormalStressHarmonic<'a> {
    law: &'a TensorProny,
    storage_pa: [f64; 25],
    loss_pa: [f64; 25],
    strain_map: [[C64; 5]; 6],
    omega_rad_s: f64,
    plane: MaterialPlane,
}

fn indices(plane: MaterialPlane) -> ([usize; 5], usize) {
    match plane {
        MaterialPlane::Axes12 => ([0, 1, 3, 5, 4], 2),
        MaterialPlane::Axes21 => ([1, 0, 3, 4, 5], 2),
        MaterialPlane::Axes13 => ([0, 2, 5, 3, 4], 1),
        MaterialPlane::Axes23 => ([1, 2, 4, 3, 5], 0),
    }
}

fn invalid(what: &'static str) -> ViscoError {
    ViscoError::Parameters { what }
}

impl TensorProny {
    /// Eliminate only `sigma_nn` from the complete original complex tensor.
    ///
    /// One complex scalar solve gives the eliminated normal strain for all
    /// five retained engineering coordinates. Storage and loss are never
    /// eliminated independently, and individual relaxation branches are not
    /// separately reduced. The original scaled complex division and factored
    /// work contractions retain representable physical coefficients.
    ///
    /// The storage/loss pair is the real constrained work pullback
    /// `Re(T* C' T), Re(T* C'' T)`, algebraically equal to the real/imaginary
    /// parts of the complex Schur complement. This uses original positive
    /// work instead of subtracting large, nearly equal complex stiffnesses.
    /// The shared passive-work factor admits both resulting operators.
    ///
    /// # Errors
    /// Refuses frequencies outside the original source band, an unresolved
    /// normal constraint, nonrepresentable coefficients, nonpositive storage
    /// or nonpassive loss. No clipping, damping floor or shear modulus is added.
    pub fn prepare_normal_stress_harmonic(
        &self,
        plane: MaterialPlane,
        omega: f64,
    ) -> Result<NormalStressHarmonic<'_>, ViscoError> {
        let (storage, loss) = self.complex_modulus(omega)?;
        let (inside, normal) = indices(plane);
        let normal_block = [[C64::new(storage[normal][normal], loss[normal][normal])]];
        let normal_inside = [inside.map(|i| C64::new(storage[normal][i], loss[normal][i]))];
        let constraint = solve_complex(normal_block, normal_inside)?[0];
        let shear = core::f64::consts::FRAC_1_SQRT_2;
        let engineering = [1.0, 1.0, shear, shear, shear];
        let mut strain_map = [[C64::ZERO; 5]; 6];
        for i in 0..5 {
            strain_map[inside[i]][i] = C64::from_re(engineering[i]);
            // The original complex modulus uses C' + i C''. Physical
            // negative-time phasors require its conjugate constraint.
            strain_map[normal][i] = -constraint[i].conj().scale(engineering[i]);
        }
        let storage_pa = constrained_work(&storage, &strain_map)?;
        let loss_pa = constrained_work(&loss, &strain_map)?;
        admit_work(&storage_pa, true)?;
        admit_work(&loss_pa, false)?;
        Ok(NormalStressHarmonic {
            law: self,
            storage_pa,
            loss_pa,
            strain_map,
            omega_rad_s: omega,
            plane,
        })
    }
}

impl<'a> NormalStressHarmonic<'a> {
    /// Original immutable full tensor law; no source or branch is replaced.
    pub const fn source(&self) -> &'a TensorProny {
        self.law
    }

    /// Row-major 5x5 engineering storage modulus [Pa], including all cross terms.
    pub const fn storage_pa(&self) -> &[f64; 25] {
        &self.storage_pa
    }

    /// Row-major 5x5 engineering loss modulus [Pa], with the passive sign.
    pub const fn loss_pa(&self) -> &[f64; 25] {
        &self.loss_pa
    }

    /// Admitted physical angular frequency [rad/s].
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Ordered in-plane material axes; the remaining axis is the normal.
    pub const fn plane(&self) -> MaterialPlane {
        self.plane
    }

    /// Recover all six Mandel strain phasors and admit their full cycle.
    ///
    /// The normal strain follows the SAME complex constraint as the prepared
    /// modulus. The original source ceiling applies to the maximum Frobenius
    /// norm over the complete physical cycle, including both shear strains
    /// and the eliminated normal. No phase sampling or heap allocation occurs.
    /// This checks one physical point; spatial-domain bounds belong to the
    /// consuming shell/solid discretization. At DC all phasors must be real.
    pub fn recover_strain(&self, retained: &[C64; 5]) -> Result<[C64; 6], ViscoError> {
        if retained.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            || (self.omega_rad_s == 0.0 && retained.iter().any(|z| z.im != 0.0))
        {
            return Err(invalid(
                "normal-stress harmonic strain must be finite and real at DC",
            ));
        }
        let strain = projection::apply(&self.strain_map, retained);
        self.law.validate_harmonic_mandel_strain(&strain)?;
        Ok(strain)
    }

    /// Exact storage/loss slopes with respect to angular frequency [Pa s].
    ///
    /// Differentiates the original complete constrained tensor:
    /// `Q_omega = T^T C_omega T`, using the already prepared normal constraint.
    /// Complex symmetry removes the derivative of that constraint; the
    /// transpose is NOT a conjugate transpose. The positive-time constitutive
    /// map is recovered by conjugating the stored physical strain map.
    ///
    /// Source derivatives retain all Maxwell times and the original Kelvin
    /// viscosity. There is no perturbation, new solve, independent branch
    /// elimination or effective-viscosity subtraction. Signed slopes need
    /// not be passive. The algebraic one-sided slope is available at admitted
    /// band endpoints, including DC; a nonrepresentable tangent refuses
    /// independently of an already valid primal material.
    pub fn frequency_derivative(&self) -> Result<([f64; 25], [f64; 25]), ViscoError> {
        let (storage, loss) = super::modulus_frequency_derivative(self.law, self.omega_rad_s)?;
        let derivative: [[C64; 6]; 6] = core::array::from_fn(|i| {
            core::array::from_fn(|j| C64::new(storage[i][j], loss[i][j]))
        });
        let map = self.strain_map.map(|row| row.map(C64::conj));
        let (mut storage_rate, mut loss_rate) = ([0.0; 25], [0.0; 25]);
        for i in 0..5 {
            for j in 0..=i {
                let value = projection::project(&derivative, &map, i, j, 1.0);
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(invalid(
                        "normal-stress physical frequency derivative is not representable",
                    ));
                }
                storage_rate[5 * i + j] = value.re;
                storage_rate[5 * j + i] = value.re;
                loss_rate[5 * i + j] = value.im;
                loss_rate[5 * j + i] = value.im;
            }
        }
        Ok((storage_rate, loss_rate))
    }
}

// A real symmetric energy operator is pulled back by the complete complex
// constraint. Factored summation restores each physical coefficient only
// after its signed contractions, preserving weak independent directions.
fn constrained_work(
    tensor: &Mandel6,
    map: &[[C64; 5]; 6],
) -> Result<[f64; 25], ViscoError> {
    let mut out = [0.0; 25];
    for i in 0..5 {
        for j in 0..=i {
            let mut sum = Sum::default();
            for row in 0..6 {
                for column in 0..6 {
                    let (left, right) = (map[row][i], map[column][j]);
                    sum.add([left.re, tensor[row][column], right.re, 1.0]);
                    sum.add([left.im, tensor[row][column], right.im, 1.0]);
                }
            }
            let value = sum.value();
            if !value.is_finite() {
                return Err(invalid(
                    "normal-stress constrained work is not representable",
                ));
            }
            out[5 * i + j] = value;
            out[5 * j + i] = value;
        }
    }
    Ok(out)
}

// Reuse the original fixed-six-component passive-work owner. The final
// padded direction is exactly zero; strict storage requires all five
// physical channels. No eigenvalue tolerance tied to a global scale enters.
fn admit_work(matrix: &[f64; 25], strict: bool) -> Result<(), ViscoError> {
    let mut block = [[0.0; 6]; 6];
    for i in 0..5 {
        for j in 0..5 {
            block[i][j] = matrix[5 * i + j];
        }
    }
    let factor = crate::plane_stress_viscosity::factor_passive_work_matrix(&block).map_err(|_| {
        invalid("normal-stress reduced work is not finite, representable and passive")
    })?;
    if strict && (0..6).filter(|&j| factor.iter().any(|row| row[j] != 0.0)).count() != 5 {
        return Err(invalid(
            "normal-stress reduced storage is not positive definite",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_constraint_keeps_a_large_finite_ratio_without_rhs_rescaling() {
        let solution = solve_complex(
            [[C64::new(1e-250, 1e-250)]],
            [[C64::new(1e50, -1e50), C64::ZERO, C64::from_re(1e-250)]],
        )
        .unwrap();
        assert_eq!(solution[0][0].re, 0.0);
        assert!((solution[0][0].im / -1e300 - 1.0).abs() < 3e-15);
        assert_eq!(solution[0][1], C64::ZERO);
        assert_eq!(solution[0][2], C64::new(0.5, -0.5));
    }
}
