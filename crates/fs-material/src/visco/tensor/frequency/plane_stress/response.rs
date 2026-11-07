//! Prepared physical plane-stress response and original strain-domain admission.

use super::{MaterialPlane, TensorProny, ViscoError};
use fs_math::c64::C64;

/// Immutable frequency-specific plane-stress material for `exp(-i omega t)`.
///
/// The full original tensor is eliminated once. Storage/loss and the complete
/// transverse-strain map therefore describe the same anisotropic constraint.
/// No relaxation branch is condensed independently. The original material's
/// frequency band and reduced passivity are admitted during preparation, and
/// its strain ceiling is retained for every recovered physical response.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneStressHarmonic {
    storage_pa: [f64; 9],
    loss_pa: [f64; 9],
    strain_map: [[C64; 3]; 6],
    strain_limit: f64,
    omega_rad_s: f64,
    plane: MaterialPlane,
}

impl TensorProny {
    /// Prepare full plane-stress modulus and strain recovery at one frequency.
    ///
    /// The engineering input is `[epsilon_x, epsilon_y, gamma_xy]` in the
    /// ordered `plane`. Returned full strains use the original six-component
    /// material Mandel frame. Physical peak phasors use `exp(-i omega t)`, so
    /// stress is `(storage - i loss) * strain`. At DC only real strains admit.
    /// This is a harmonic constitutive query, not a hereditary time-step law.
    pub fn prepare_plane_stress_harmonic(
        &self,
        plane: MaterialPlane,
        omega: f64,
    ) -> Result<PlaneStressHarmonic, ViscoError> {
        let (storage_pa, loss_pa, constraint) =
            self.plane_stress_complex_with_constraint(plane, omega)?;
        let (inside, outside) = match plane {
            MaterialPlane::Axes12 => ([0usize, 1, 3], [2usize, 4, 5]),
            MaterialPlane::Axes21 => ([1usize, 0, 3], [2usize, 4, 5]),
            MaterialPlane::Axes13 => ([0usize, 2, 5], [1usize, 3, 4]),
            MaterialPlane::Axes23 => ([1usize, 2, 4], [0usize, 3, 5]),
        };
        let engineering = [1.0, 1.0, core::f64::consts::FRAC_1_SQRT_2];
        let mut strain_map = [[C64::ZERO; 3]; 6];
        for i in 0..3 {
            strain_map[inside[i]][i] = C64::from_re(engineering[i]);
            for j in 0..3 {
                // The shared constitutive elimination uses C' + i C''.
                // Conjugate its constraint for the physical negative-time
                // phasors, including complex anisotropic transverse coupling.
                strain_map[outside[i]][j] = -constraint[i][j].conj().scale(engineering[j]);
            }
        }
        Ok(PlaneStressHarmonic {
            storage_pa,
            loss_pa,
            strain_map,
            strain_limit: self.strain_limit,
            omega_rad_s: omega,
            plane,
        })
    }
}

impl PlaneStressHarmonic {
    /// Row-major 3x3 engineering-shear storage modulus [Pa].
    pub const fn storage_pa(&self) -> &[f64; 9] {
        &self.storage_pa
    }

    /// Row-major 3x3 engineering-shear loss modulus [Pa], with passive sign.
    pub const fn loss_pa(&self) -> &[f64; 9] {
        &self.loss_pa
    }

    /// Admitted angular frequency [rad/s].
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Ordered engineering-shear plane in the original material frame.
    pub const fn plane(&self) -> MaterialPlane {
        self.plane
    }

    /// Recover all six physical Mandel strain peak phasors at one point.
    ///
    /// The full-cycle maximum Frobenius norm must satisfy the original source
    /// strain bound, including every eliminated transverse component. Testing
    /// only the real and imaginary snapshots can miss the actual maximum;
    /// summing their norms can instead reject an admissible rotating strain.
    /// The scaled Gram eigenvalue admits the actual ellipse without either
    /// approximation. There is no allocation, material solve or phase sampling
    /// in this per-point query. At DC an imaginary component is refused.
    pub fn recover_strain(&self, in_plane: &[C64; 3]) -> Result<[C64; 6], ViscoError> {
        if in_plane
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
            || (self.omega_rad_s == 0.0 && in_plane.iter().any(|z| z.im != 0.0))
        {
            return Err(ViscoError::Parameters {
                what: "plane-stress harmonic strain must be finite and real at DC",
            });
        }
        let strain = super::projection::apply(&self.strain_map, in_plane);
        super::super::admit_harmonic_strain(&strain, self.strain_limit)?;
        Ok(strain)
    }
}
