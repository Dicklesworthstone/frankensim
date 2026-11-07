//! Frequency tangent of the original complex plane-stress elimination.

use super::super::{MaterialPlane, TensorProny, ViscoError};
use fs_math::c64::C64;

mod projection;
mod response;
pub use response::PlaneStressHarmonic;

impl TensorProny {
    /// Differentiate the full complex plane-stress modulus with respect to omega.
    ///
    /// The returned storage/loss slopes [Pa s] are row-major 3x3 matrices in
    /// the SAME ordered engineering-shear coordinates as `plane_stress_complex`.
    /// The original primal query admits the frequency band, out-of-plane solve
    /// and reduced passivity first. Derivatives are signed; they need not be
    /// positive definite or passive operators themselves.
    ///
    /// For X=C_oo^-1 C_oi, this differentiates the actual complex constraint:
    /// X_omega=C_oo^-1(C_oi,omega-C_oo,omega X), followed by
    /// Q_omega=C_ii,omega-C_io,omega X-C_io X_omega. The complex symmetric
    /// products use no conjugation. Complex symmetry reduces this exact chain
    /// to Q_omega=T^T C_omega T with T=[I;-X]. The original small primal solve
    /// is reused; no derivative solve, inverse or perturbed material is needed.
    /// Independently eliminating storage and loss would discard their
    /// frequency-dependent transverse-strain coupling.
    ///
    /// Only storage/loss slopes are requested. An unrepresentable derivative of
    /// C''/omega does not invalidate a finite plane-stress modulus derivative.
    /// This is a local material query, not plate/shell assembly, a hereditary
    /// time-step reduction or a bound on spectral truncation.
    ///
    /// # Errors
    /// Retains `plane_stress_complex` admission errors and refuses a
    /// nonrepresentable original or reduced physical frequency derivative.
    pub fn plane_stress_complex_frequency_derivative(
        &self,
        plane: MaterialPlane,
        omega: f64,
    ) -> Result<([f64; 9], [f64; 9]), ViscoError> {
        let (_, _, x) = self.plane_stress_complex_with_constraint(plane, omega)?;
        let (storage_rate, loss_rate) = super::modulus_frequency_derivative(self, omega)?;
        if storage_rate
            .iter()
            .flatten()
            .chain(loss_rate.iter().flatten())
            .all(|v| *v == 0.0)
        {
            return Ok(([0.0; 9], [0.0; 9]));
        }
        let (inside, outside) = match plane {
            MaterialPlane::Axes12 => ([0usize, 1, 3], [2usize, 4, 5]),
            MaterialPlane::Axes21 => ([1usize, 0, 3], [2usize, 4, 5]),
            MaterialPlane::Axes13 => ([0usize, 2, 5], [1usize, 3, 4]),
            MaterialPlane::Axes23 => ([1usize, 2, 4], [0usize, 3, 5]),
        };
        let dc: [[C64; 6]; 6] = core::array::from_fn(|i| {
            core::array::from_fn(|j| C64::new(storage_rate[i][j], loss_rate[i][j]))
        });
        let mut map = [[C64::ZERO; 3]; 6];
        for i in 0..3 {
            map[inside[i]][i] = C64::from_re(1.0);
            for j in 0..3 {
                map[outside[i]][j] = -x[i][j];
            }
        }
        let engineering = [1.0, 1.0, core::f64::consts::FRAC_1_SQRT_2];
        let (mut storage_rate, mut loss_rate) = ([0.0; 9], [0.0; 9]);
        for i in 0..3 {
            for j in 0..=i {
                let scale = engineering[i] * engineering[j];
                let value = projection::project(&dc, &map, i, j, scale);
                let (re, im) = (value.re, value.im);
                if !re.is_finite() || !im.is_finite() {
                    return Err(ViscoError::Parameters {
                        what: "plane-stress physical frequency derivative is not representable",
                    });
                }
                storage_rate[3 * i + j] = re;
                storage_rate[3 * j + i] = re;
                loss_rate[3 * i + j] = im;
                loss_rate[3 * j + i] = im;
            }
        }
        Ok((storage_rate, loss_rate))
    }
}
