//! Reciprocal tensor supports on the original isotropic elastic half-space.
//!
//! Resistance and reactance need not share principal axes. A passive complex-
//! symmetric `Z=R+iX` satisfies `R >= 0`; the shared original-unit passive-work
//! factor admits R and retains its weak and null directions. The material
//! supplies every P/SV/SH trace, including outgoing evanescence.
//!
//! One native 6x6 system solves the three reflected amplitudes and the actual
//! wall velocity together: `V_out b-v=-V_in a`, `T_out b+Z v=-T_in a`.
//! Keeping wall motion as an unknown avoids cancellation at stiff supports
//! and avoids squared conditioning from normal-equation reconstruction.

use super::*;
use crate::plane_stress_viscosity::factor_passive_resultant;

/// Passive reciprocal specific mechanical impedance [Pa s/m] in the surface
/// frame `(x, y, z)`: x follows the tangential wavevector, y is SH and z points
/// outward. Peak phasors use `exp(-i omega t)` and traction is `-Z velocity`.
///
/// Both triangles must agree exactly. The real symmetric part is admitted by
/// the shared positive-semidefinite work factor; the imaginary symmetric part
/// may have either sign and need not commute with resistance. This frozen
/// frequency response supplies no causal realization, active/nonreciprocal
/// law, bulk anisotropy or experimental material authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElasticSurfaceImpedanceTensor {
    matrix_pa_s_m: [[C64; 3]; 3],
    resistance_factor: [[f64; 6]; 3],
}

impl ElasticSurfaceImpedanceTensor {
    /// Admit one finite reciprocal impedance with passive resistance.
    /// Zero and rank-deficient resistance are allowed without a loss floor.
    pub fn new(matrix_pa_s_m: [[C64; 3]; 3]) -> Result<Self, MaterialError> {
        if matrix_pa_s_m.iter().flatten().any(|&z| !finite(z))
            || (0..3).any(|i| (0..3).any(|j| matrix_pa_s_m[i][j] != matrix_pa_s_m[j][i]))
        {
            return Err(invalid(
                "elastic support tensor must be finite and reciprocal",
            ));
        }
        let resistance = core::array::from_fn(|i| matrix_pa_s_m[i / 3][i % 3].re);
        let resistance_factor = factor_passive_resultant(&resistance)?;
        Ok(Self {
            matrix_pa_s_m,
            resistance_factor,
        })
    }

    /// Retain three authored dissipative channels `G` and reciprocal reactance
    /// `X`, with `Z = G G^T + i X`. Entries of `G` have units sqrt(Pa s/m).
    ///
    /// This preserves physical null directions and weak channels that a rounded
    /// Cartesian resistance matrix alone cannot identify. The supplied channels
    /// are used directly for work; only the traction matrix is formed, with
    /// compensated products. No inferred factor, fitted loss or floor is added.
    pub fn from_resistance_gram(
        gram: [[f64; 3]; 3],
        reactance: [[f64; 3]; 3],
    ) -> Result<Self, MaterialError> {
        if gram
            .iter()
            .chain(&reactance)
            .flatten()
            .any(|v| !v.is_finite())
            || (0..3).any(|i| (0..3).any(|j| reactance[i][j] != reactance[j][i]))
        {
            return Err(invalid(
                "elastic support channels must be finite with reciprocal reactance",
            ));
        }
        let mut matrix_pa_s_m = [[C64::ZERO; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let mut resistance = Dd::ZERO;
                for channel in 0..3 {
                    resistance = resistance
                        + Dd::from_f64(gram[i][channel]) * Dd::from_f64(gram[j][channel]);
                }
                let value = resistance.to_f64();
                if !value.is_finite()
                    || (i == j && value == 0.0 && gram[i].iter().any(|&g| g != 0.0))
                {
                    return Err(invalid("elastic support resistance is unrepresentable"));
                }
                matrix_pa_s_m[i][j] = C64::new(value, reactance[i][j]);
                matrix_pa_s_m[j][i] = matrix_pa_s_m[i][j];
            }
        }
        let resistance_factor = core::array::from_fn(|i| {
            core::array::from_fn(|channel| if channel < 3 { gram[i][channel] } else { 0.0 })
        });
        Ok(Self {
            matrix_pa_s_m,
            resistance_factor,
        })
    }

    /// Complex coefficients in the surface's Cartesian frame [Pa s/m].
    /// A supplied matrix is retained exactly; authored Gram channels are reduced
    /// once for traction while their original values remain the work owner.
    pub const fn matrix_pa_s_m(&self) -> &[[C64; 3]; 3] {
        &self.matrix_pa_s_m
    }

    pub(in crate::elastic::acoustic::interface) fn traction(&self, velocity: [C64; 3]) -> [C64; 3] {
        core::array::from_fn(|axis| {
            -complex_dot((0..3).map(|j| (self.matrix_pa_s_m[axis][j], velocity[j])))
        })
    }

    pub(in crate::elastic::acoustic::interface) fn absorbed_power(
        &self,
        velocity: [C64; 3],
    ) -> Result<f64, MaterialError> {
        // R=L L^T. Sum |L^T v|^2/2 rather than cancelling Cartesian cross
        // terms or subtracting outgoing power from incoming power.
        let mut power = Dd::ZERO;
        for channel in 0..6 {
            let amplitude = complex_dot((0..3).map(|axis| {
                (
                    C64::from_re(self.resistance_factor[axis][channel]),
                    velocity[axis],
                )
            }));
            let root = amplitude.abs() * core::f64::consts::FRAC_1_SQRT_2;
            let value = root * root;
            if !finite(amplitude) || !value.is_finite() || (amplitude != C64::ZERO && value == 0.0)
            {
                return Err(invalid("tensor support absorption is unrepresentable"));
            }
            power = power + Dd::from_f64(value);
        }
        let value = power.to_f64();
        if !value.is_finite() {
            return Err(invalid("tensor support absorption is unrepresentable"));
        }
        Ok(value)
    }
}

// Preserve products before coherent cancellation, including a nearly undamped
// direction of the original Gram factor. This uses the existing DD owner.
pub(in crate::elastic::acoustic::interface) fn complex_dot(terms: impl IntoIterator<Item = (C64, C64)>) -> C64 {
    let mut re = Dd::ZERO;
    let mut im = Dd::ZERO;
    for (a, b) in terms {
        let (ar, ai) = (Dd::from_f64(a.re), Dd::from_f64(a.im));
        let (br, bi) = (Dd::from_f64(b.re), Dd::from_f64(b.im));
        re = re + ar * br - ai * bi;
        im = im + ar * bi + ai * br;
    }
    C64::new(re.to_f64(), im.to_f64())
}

impl PreparedElasticSurface<'_> {
    pub(super) fn prepare_tensor(
        &mut self,
        support: ElasticSurfaceImpedanceTensor,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MaterialError> {
        let mut matrix = [C64::ZERO; 36];
        let mut scales = [0.0_f64; 6];
        for axis in 0..3 {
            poll(cancelled)?;
            for mode in 0..3 {
                matrix[6 * axis + mode] = self.reflected[mode].velocity[axis];
                matrix[6 * (3 + axis) + mode] = self.reflected[mode].traction_pa_s_m[axis];
                matrix[6 * (3 + axis) + 3 + mode] = support.matrix_pa_s_m[axis][mode];
            }
            matrix[6 * axis + 3 + axis] = -C64::ONE;
        }
        for row in 0..6 {
            scales[row] = matrix[6 * row..6 * (row + 1)]
                .iter()
                .fold(0.0_f64, |scale, z| scale.max(z.re.abs()).max(z.im.abs()));
            if scales[row] == 0.0 || !scales[row].is_finite() {
                return Err(invalid("tensor support boundary row is unrepresentable"));
            }
            for column in 0..6 {
                matrix[6 * row + column] = divided(matrix[6 * row + column], scales[row]);
            }
        }
        let factor = lu_complex_with_cancel(&matrix, 6, &mut *cancelled).map_err(solve_error)?;
        let mut velocity = [[C64::ZERO; 3]; 3];
        for incident in 0..3 {
            if self.incoming[incident].normal_power_factor_pa_s_m == 0.0 {
                continue;
            }
            let wave = self.incoming[incident];
            let rhs: [C64; 6] = core::array::from_fn(|row| {
                -divided(
                    if row < 3 {
                        wave.velocity[row]
                    } else {
                        wave.traction_pa_s_m[row - 3]
                    },
                    scales[row],
                )
            });
            let mut solution = rhs;
            factor
                .solve_with_cancel(&mut solution, &mut *cancelled)
                .map_err(solve_error)?;
            if solution.iter().any(|&z| !finite(z)) {
                return Err(invalid("tensor support reflection is unrepresentable"));
            }
            for row in 0..6 {
                let terms: [C64; 6] = core::array::from_fn(|j| matrix[6 * row + j] * solution[j]);
                let residual = (complex_sum(terms) - rhs[row]).abs();
                let magnitude = real_sum(terms.map(C64::abs)) + rhs[row].abs();
                if !residual.is_finite() || (magnitude > 0.0 && residual / magnitude > 1e-11) {
                    return Err(invalid("tensor support boundary solve did not converge"));
                }
            }
            for axis in 0..3 {
                self.scattering[axis][incident] = solution[axis];
                velocity[axis][incident] = solution[3 + axis];
            }
        }
        self.tensor_boundary_velocity = Some(velocity);
        Ok(())
    }
}
