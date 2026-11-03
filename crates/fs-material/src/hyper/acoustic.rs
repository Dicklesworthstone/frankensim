//! Infinitesimal elastic waves in a finitely deformed hyperelastic solid.
//!
//! For a current unit wave normal `n`, let `q=F^T n` and `A=dP/dF`.
//! The spatial acoustic tensor is `Q_ik=A_iJkL q_J q_L / J`, and mass
//! conservation gives `rho=rho0/J`. Its eigenproblem is `Q a=rho c² a`.
//! The ORIGINAL material Hessian includes the initial-stress contribution;
//! using a small-strain modulus or adding another geometric term is incorrect.
//! See Abiza, Destrade & Ogden (2012), doi:10.1016/j.wavemoti.2011.12.002.
//!
//! Preparation validates the law, deformation and density once and caches
//! the native Hessian. Directional queries only contract it and use the shared
//! `fs-la` symmetric eigensolver. Signed squared speeds retain directional
//! instability; they are never clipped or replaced by their absolute values.
//! Positive eigenvalues in one direction do not certify global ellipticity.
//! Positive modes also carry their spatial energy velocity, evaluated from
//! that SAME Hessian and deformation. This is the group velocity on a smooth
//! separated branch; a repeated root retains its chosen polarization's flux.
//! This uses the unconstrained three-dimensional compressible law, not a
//! plane-stress/incompressibility projection or a Maxwell algorithmic tangent.
//! It is a frozen constitutive prediction, not attenuation,
//! dispersion, a propagation solver, a reduced membrane wave law or measured
//! adiabatic/ultrasonic validation of isothermal material parameters.

use super::cards::HyperelasticMaterial;
use super::{Hyperelastic, det3};
use crate::MaterialError;
use fs_la::eigen::jacobi_eigh;

/// Original tangent and mass at one admitted finite-deformation state.
#[derive(Clone, Debug, PartialEq)]
pub struct HyperelasticAcousticState {
    deformation_gradient: [f64; 9],
    tangent: [[f64; 9]; 9],
    reference_density_kg_m3: f64,
    current_density_kg_m3: f64,
    jacobian: f64,
}

/// One signed acoustic eigenmode. A nonpositive squared speed does not
/// describe a propagating real phase wave in this direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HyperelasticWaveMode {
    /// Unit displacement polarization in the current spatial frame.
    pub polarization: [f64; 3],
    /// Signed acoustic eigenvalue `rho_current*c²`, in Pa.
    pub directional_modulus_pa: f64,
    /// Signed squared phase speed, in m²/s²; negative values retain instability.
    pub squared_phase_speed_m2_s2: f64,
    /// Positive phase speed, in m/s, only when the squared speed is positive.
    pub phase_speed_m_s: Option<f64>,
    /// Spatial incremental-energy flux divided by mean wave energy [m/s].
    /// Equals `grad_k omega` for a smooth separated branch, and obeys
    /// `n.dot(energy_velocity)=phase_speed`. At a repeated root it describes
    /// the selected polarization. Nonpositive squared speeds return `None`.
    pub energy_velocity_m_s: Option<[f64; 3]>,
}

/// Three modes in ascending signed squared-speed order for one direction.
/// At repeated eigenvalues only the polarization subspace is meaningful.
#[derive(Clone, Debug, PartialEq)]
pub struct DirectionalHyperelasticWaves {
    /// Normalized wave-normal direction in the current spatial frame.
    pub direction: [f64; 3],
    /// Spatial acoustic tensor `A_iJkL (F^T n)_J (F^T n)_L / J`, in Pa.
    pub acoustic_tensor_pa: [[f64; 3]; 3],
    /// Reference density divided by the original deformation's determinant.
    pub current_density_kg_m3: f64,
    /// Ascending signed modes without a stability projection.
    pub modes: [HyperelasticWaveMode; 3],
    /// Largest componentwise eigen residual divided by the acoustic matrix's
    /// largest absolute entry. Zero when the complete matrix is zero.
    pub relative_eigen_residual: f64,
}

fn state_error(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("hyperelastic acoustics: {what}"),
    }
}

impl Hyperelastic {
    /// Freeze the native 3-D compressible tangent and mass for repeated queries.
    /// `reference_density_kg_m3` is mass per ORIGINAL reference volume.
    /// Publicly mutable material parameters are readmitted before preparation;
    /// subsequent changes to the source law do not mutate this prepared state.
    pub fn acoustic_state(
        &self,
        deformation_gradient: &[f64; 9],
        reference_density_kg_m3: f64,
    ) -> Result<HyperelasticAcousticState, MaterialError> {
        HyperelasticAcousticState::from_material(
            self,
            deformation_gradient,
            reference_density_kg_m3,
        )
    }
}

impl HyperelasticAcousticState {
    /// Freeze a material's consistent tangent with respect to TOTAL F and
    /// its mass per ORIGINAL reference volume. Reference-deformation adapters
    /// retain their complete energy/stress/tangent pullback through this trait.
    ///
    /// Underlying law parameters are readmitted, while the adapter owns the
    /// elastic deformation and stretch envelope. In particular, a thermal
    /// `F=Fe Fth` need not put total F inside the bare law's elastic envelope.
    /// Total F still needs finite entries and a finite positive determinant;
    /// reference/current density and the returned tangent must be representable.
    /// Signed directional instabilities remain results of subsequent queries.
    pub fn from_material<M: HyperelasticMaterial + ?Sized>(
        material: &M,
        deformation_gradient: &[f64; 9],
        reference_density_kg_m3: f64,
    ) -> Result<Self, MaterialError> {
        if !(reference_density_kg_m3.is_finite() && reference_density_kg_m3 > 0.0) {
            return Err(MaterialError::Parameters {
                what: "hyperelastic acoustic reference density must be finite and positive".into(),
            });
        }
        let law = material.hyperelastic();
        Hyperelastic::new(law.model.clone(), law.stretch_limit)?;
        if deformation_gradient.iter().any(|value| !value.is_finite()) {
            return Err(state_error("total deformation gradient must be finite"));
        }
        let jacobian = det3(deformation_gradient);
        if !(jacobian.is_finite() && jacobian > 0.0) {
            return Err(state_error(
                "total deformation determinant must be finite and positive",
            ));
        }
        let current_density_kg_m3 = reference_density_kg_m3 / jacobian;
        if !(current_density_kg_m3.is_finite() && current_density_kg_m3 > 0.0) {
            return Err(state_error("current density is unrepresentable"));
        }
        let tangent = material.consistent_tangent(deformation_gradient)?;
        if tangent.iter().flatten().any(|value| !value.is_finite()) {
            return Err(state_error(
                "material consistent tangent is unrepresentable",
            ));
        }
        Ok(Self {
            deformation_gradient: *deformation_gradient,
            tangent,
            reference_density_kg_m3,
            current_density_kg_m3,
            jacobian,
        })
    }
    /// Original admitted deformation, in row-major order.
    #[must_use]
    pub const fn deformation_gradient(&self) -> &[f64; 9] {
        &self.deformation_gradient
    }

    /// Conserved mass per reference volume, in kg/m³.
    #[must_use]
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density_kg_m3
    }

    /// Mass per current volume, in kg/m³.
    #[must_use]
    pub const fn current_density_kg_m3(&self) -> f64 {
        self.current_density_kg_m3
    }

    /// Solve the local three-dimensional wave modes for a finite nonzero
    /// current direction. Negative and zero acoustic eigenvalues are results,
    /// not material-admission errors. Nonfinite/unresolved arithmetic refuses.
    pub fn directional_waves(
        &self,
        direction: [f64; 3],
    ) -> Result<DirectionalHyperelasticWaves, MaterialError> {
        let direction_scale = direction.iter().map(|v| v.abs()).fold(0.0, f64::max);
        if direction.iter().any(|v| !v.is_finite()) || direction_scale <= 0.0 {
            return Err(MaterialError::Parameters {
                what: "hyperelastic wave direction must be finite and nonzero".into(),
            });
        }
        let scaled = direction.map(|v| v / direction_scale);
        let length = fs_math::det::sqrt(scaled.iter().map(|v| v * v).sum());
        let direction = scaled.map(|v| v / length);
        let q: [f64; 3] = core::array::from_fn(|j| {
            (0..3)
                .map(|i| self.deformation_gradient[3 * i + j] * direction[i])
                .sum()
        });
        // K = J Q. Its eigenvalues divided by rho0 equal c² directly;
        // the cancellation of J avoids dividing two separately scaled values.
        let mut k = [0.0; 9];
        for i in 0..3 {
            for m in 0..3 {
                for j in 0..3 {
                    for l in 0..3 {
                        k[3 * i + m] += self.tangent[3 * i + j][3 * m + l] * q[j] * q[l];
                    }
                }
            }
        }
        if k.iter().any(|v| !v.is_finite()) {
            return Err(state_error("acoustic tensor is unrepresentable"));
        }
        let scale = k.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let normalizer = if scale > 0.0 { scale } else { 1.0 };
        let normalized = k.map(|v| v / normalizer);
        for i in 0..3 {
            for j in 0..i {
                if (normalized[3 * i + j] - normalized[3 * j + i]).abs() > 512.0 * f64::EPSILON {
                    return Err(state_error("original acoustic tensor is not symmetric"));
                }
            }
        }
        let (values, vectors) = jacobi_eigh(&normalized, 3);
        let mut modes = [HyperelasticWaveMode {
            polarization: [0.0; 3],
            directional_modulus_pa: 0.0,
            squared_phase_speed_m2_s2: 0.0,
            phase_speed_m_s: None,
            energy_velocity_m_s: None,
        }; 3];
        let mut relative_eigen_residual = 0.0_f64;
        for mode in 0..3 {
            let mut polarization = core::array::from_fn(|i| vectors[3 * i + mode]);
            let norm: f64 = polarization.iter().map(|v| v * v).sum();
            if !norm.is_finite() || (norm - 1.0).abs() > 1e-11 {
                return Err(state_error("wave polarization is not normalized"));
            }
            for i in 0..3 {
                let residual = (0..3)
                    .map(|j| normalized[3 * i + j] * polarization[j])
                    .sum::<f64>()
                    - values[mode] * polarization[i];
                if !residual.is_finite() || residual.abs() > 1e-11 {
                    return Err(state_error("acoustic eigenproblem did not converge"));
                }
                relative_eigen_residual = relative_eigen_residual.max(residual.abs());
            }
            let eigenvalue = values[mode] * normalizer;
            let directional_modulus_pa = eigenvalue / self.jacobian;
            let squared_phase_speed_m2_s2 = eigenvalue / self.reference_density_kg_m3;
            if !directional_modulus_pa.is_finite()
                || !squared_phase_speed_m2_s2.is_finite()
                || (values[mode] != 0.0 && eigenvalue == 0.0)
                || (eigenvalue != 0.0
                    && (directional_modulus_pa == 0.0 || squared_phase_speed_m2_s2 == 0.0))
            {
                return Err(state_error("signed wave response is unrepresentable"));
            }
            let phase_speed_m_s = (squared_phase_speed_m2_s2 > 0.0)
                .then(|| fs_math::det::sqrt(squared_phase_speed_m2_s2));
            let pivot = (1..3).fold(0, |p, i| {
                if polarization[i].abs() > polarization[p].abs() {
                    i
                } else {
                    p
                }
            });
            if polarization[pivot] < 0.0 {
                polarization = polarization.map(|v| -v);
            }
            let energy_velocity_m_s = if let Some(speed) = phase_speed_m_s {
                // In the fixed current frame, q=F^T n and the dispersion is
                // rho0*omega(k)^2 = a_i A_iJmL (F^T k)_J (F^T k)_L a_m.
                // Major symmetry of the original material Hessian gives
                // v_r = a_i A_iJmL a_m q_L F_rJ/(rho0*c). J cancels exactly
                // as in the squared speed; no extra prestress term belongs here.
                let stress: [f64; 9] = core::array::from_fn(|ij| {
                    (0..9)
                        .map(|ml| {
                            (self.tangent[ij][ml] / normalizer) * polarization[ml / 3] * q[ml % 3]
                        })
                        .sum()
                });
                let velocity = core::array::from_fn(|r| {
                    let flux: f64 = (0..9)
                        .map(|ij| {
                            polarization[ij / 3]
                                * stress[ij]
                                * self.deformation_gradient[3 * r + ij % 3]
                        })
                        .sum();
                    (flux / values[mode]) * speed
                });
                if velocity.iter().any(|v| !v.is_finite()) {
                    return Err(state_error("wave energy velocity is unrepresentable"));
                }
                Some(velocity)
            } else {
                None
            };
            modes[mode] = HyperelasticWaveMode {
                polarization,
                directional_modulus_pa,
                squared_phase_speed_m2_s2,
                phase_speed_m_s,
                energy_velocity_m_s,
            };
        }
        for i in 0..3 {
            for j in 0..i {
                let dot: f64 = (0..3)
                    .map(|k| modes[i].polarization[k] * modes[j].polarization[k])
                    .sum();
                if !dot.is_finite() || dot.abs() > 1e-11 {
                    return Err(state_error("wave polarizations are not orthogonal"));
                }
            }
        }
        let acoustic_tensor_pa =
            core::array::from_fn(|i| core::array::from_fn(|j| k[3 * i + j] / self.jacobian));
        if acoustic_tensor_pa.iter().flatten().any(|v| !v.is_finite()) {
            return Err(state_error("spatial acoustic tensor is unrepresentable"));
        }
        Ok(DirectionalHyperelasticWaves {
            direction,
            acoustic_tensor_pa,
            current_density_kg_m3: self.current_density_kg_m3,
            modes,
            relative_eigen_residual,
        })
    }
}
