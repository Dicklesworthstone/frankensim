//! Elastic waves and hydrostatic response from the admitted elastic tensor.
//!
//! In the orthonormal Mandel basis, hydrostatic stress `-p I` gives strain
//! `-p C^-1 I` and compressibility `I:C^-1:I`. A plane wave satisfies the
//! Christoffel equation `rho c^2 a_i = C_ijkl n_j n_l a_k`. Both calculations
//! use the SAME stiffness and density as the mechanical operator, including
//! normal/shear coupling and material orientation. There is no material table.
//!
//! These are predictions of a frozen, unstressed, small-strain elastic law.
//! Static/isothermal source moduli do not acquire measured adiabatic or
//! ultrasonic authority. Attenuation, dispersion, prestress, finite strain and
//! thermoelastic relaxation require their respective constitutive owners.

use crate::{MaterialError, Tangent6, tensor::mandel_basis};
use fs_blake3::ContentHash;
use fs_la::{eigen::jacobi_eigh, factor::cholesky};
use fs_math::det;

/// Source-selected density, bulk modulus and shear modulus for an explicitly
/// isotropic elastic law, retaining the original property-use receipts.
pub mod sourced;

fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters {
        what: what.to_owned(),
    }
}

/// Frozen elastic state shared by wave and pressure queries.
#[derive(Clone, Debug, PartialEq)]
pub struct ElasticAcousticState {
    stiffness_mandel_pa: Tangent6,
    density_kg_m3: f64,
    material_state_identity: ContentHash,
    maximum_absolute_mandel_strain: f64,
    hydrostatic_compliance_pa_inv: [f64; 6],
    compressibility_pa_inv: f64,
    stiffness_scale_pa: f64,
}

/// One Christoffel eigenmode; phase speed is not the group velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElasticWaveMode {
    /// Unit displacement polarization in the stiffness tensor's frame.
    pub polarization: [f64; 3],
    /// Positive directional stiffness `rho*c^2` [Pa].
    pub directional_modulus_pa: f64,
    /// Plane-wave phase speed [m/s].
    pub phase_speed_m_s: f64,
    /// Traction/velocity impedance `rho*c` on the wave-normal plane [Pa s/m].
    pub impedance_pa_s_m: f64,
}

/// Three elastic modes in ascending speed order for one propagation direction.
/// In a degenerate eigenspace only the polarization subspace is meaningful.
#[derive(Clone, Debug, PartialEq)]
pub struct DirectionalElasticWaves {
    /// Original mechanical material-state identity.
    pub material_state_identity: ContentHash,
    /// Normalized wave-normal direction in the stiffness tensor's frame.
    pub direction: [f64; 3],
    /// Symmetric Christoffel matrix `C_ijkl*n_j*n_l` [Pa].
    pub christoffel_pa: [[f64; 3]; 3],
    /// Ascending phase speeds and their displacement polarizations.
    pub modes: [ElasticWaveMode; 3],
}

/// Reversible hydrostatic response of the same linear elastic law.
#[derive(Clone, Debug, PartialEq)]
pub struct HydrostaticElasticResponse {
    /// Original mechanical material-state identity.
    pub material_state_identity: ContentHash,
    /// Pressure [Pa], positive in compression and negative in tension.
    pub pressure_pa: f64,
    /// Strain in `[xx, yy, zz, sqrt(2)xy, sqrt(2)yz, sqrt(2)zx]` ordering.
    pub strain_mandel: [f64; 6],
    /// Linearized relative volume change, the trace of strain.
    pub volumetric_strain: f64,
    /// Stored elastic energy per reference volume [J/m^3].
    pub stored_energy_density_j_m3: f64,
}

impl ElasticAcousticState {
    /// Cubic crystal stiffness in its crystallographic [100]/[010]/[001]
    /// frame. C11, C12 and C44 are the three independent tensor coefficients,
    /// not isotropic Young/Poisson parameters. The shared full-tensor admission
    /// checks positive definiteness; stable negative C12 remains admissible.
    /// The source's adiabatic/isothermal convention is retained upstream.
    pub fn from_cubic(
        density_kg_m3: f64,
        c11_pa: f64,
        c12_pa: f64,
        c44_pa: f64,
        material_state_identity: ContentHash,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, MaterialError> {
        let mut stiffness = [[0.0; 6]; 6];
        for (i, row) in stiffness.iter_mut().enumerate().take(3) {
            for (j, value) in row.iter_mut().enumerate().take(3) {
                *value = if i == j { c11_pa } else { c12_pa };
            }
        }
        for (i, row) in stiffness.iter_mut().enumerate().skip(3) {
            row[i] = 2.0 * c44_pa;
        }
        Self::try_new_mandel(
            stiffness,
            density_kg_m3,
            material_state_identity,
            maximum_absolute_mandel_strain,
        )
    }

    /// Build the isotropic tensor directly from independently selected bulk
    /// and shear moduli. Neither modulus is replaced by a rounded E/nu pair.
    /// The caller explicitly chooses isotropy and the source's response regime;
    /// static moduli are not promoted to measured ultrasonic data.
    pub fn from_bulk_and_shear(
        density_kg_m3: f64,
        bulk_modulus_pa: f64,
        shear_modulus_pa: f64,
        material_state_identity: ContentHash,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, MaterialError> {
        if [bulk_modulus_pa, shear_modulus_pa]
            .into_iter()
            .any(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(invalid("isotropic bulk and shear moduli must be finite and positive"));
        }
        let lambda = bulk_modulus_pa - (2.0 / 3.0) * shear_modulus_pa;
        let mut stiffness = [[0.0; 6]; 6];
        for (i, row) in stiffness.iter_mut().enumerate().take(3) {
            row[..3].fill(lambda);
            row[i] = bulk_modulus_pa + (4.0 / 3.0) * shear_modulus_pa;
        }
        for (i, row) in stiffness.iter_mut().enumerate().skip(3) {
            row[i] = 2.0 * shear_modulus_pa;
        }
        Self::try_new_mandel(
            stiffness,
            density_kg_m3,
            material_state_identity,
            maximum_absolute_mandel_strain,
        )
    }

    /// Use an existing isotropic law and its declared strain ceiling.
    pub fn from_isotropic(
        law: &super::IsotropicElastic,
        density_kg_m3: f64,
        material_state_identity: ContentHash,
    ) -> Result<Self, MaterialError> {
        Self::try_new_mandel(
            law.stiffness(),
            density_kg_m3,
            material_state_identity,
            law.strain_limit,
        )
    }

    /// Use an existing orthotropic law in its principal frame.
    /// Oriented consumers pass their already transformed Mandel tensor instead.
    pub fn from_orthotropic(
        law: &super::OrthotropicElastic,
        density_kg_m3: f64,
        material_state_identity: ContentHash,
    ) -> Result<Self, MaterialError> {
        Self::try_new_mandel(
            law.stiffness(),
            density_kg_m3,
            material_state_identity,
            law.strain_limit,
        )
    }

    /// Admit an exact symmetric positive-definite Mandel tensor and density.
    ///
    /// The caller retains responsibility for the source domain and declares a
    /// finite positive ceiling on each absolute Mandel strain component. This
    /// does not extend a source law's validity or silently repair its tensor.
    pub fn try_new_mandel(
        stiffness_mandel_pa: Tangent6,
        density_kg_m3: f64,
        material_state_identity: ContentHash,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, MaterialError> {
        if !density_kg_m3.is_finite() || density_kg_m3 <= 0.0 {
            return Err(invalid(
                "elastic acoustic density must be finite and positive",
            ));
        }
        if !maximum_absolute_mandel_strain.is_finite() || maximum_absolute_mandel_strain <= 0.0 {
            return Err(invalid(
                "elastic acoustic strain ceiling must be finite and positive",
            ));
        }
        if material_state_identity == ContentHash([0; 32]) {
            return Err(invalid(
                "elastic acoustic material identity must not be zero",
            ));
        }
        let mut stiffness_scale_pa = 0.0_f64;
        for i in 0..6 {
            for j in 0..6 {
                let value = stiffness_mandel_pa[i][j];
                if !value.is_finite() || value.to_bits() != stiffness_mandel_pa[j][i].to_bits() {
                    return Err(invalid(
                        "elastic acoustic stiffness must be finite and exactly symmetric",
                    ));
                }
                stiffness_scale_pa = stiffness_scale_pa.max(value.abs());
            }
        }
        if stiffness_scale_pa <= 0.0 {
            return Err(invalid(
                "elastic acoustic stiffness must be positive definite",
            ));
        }
        // Scale before factorization, preserving the original coefficient bits.
        let normalized: Vec<_> = stiffness_mandel_pa
            .iter()
            .flatten()
            .map(|value| value / stiffness_scale_pa)
            .collect();
        let factor = cholesky(&normalized, 6)
            .map_err(|_| invalid("elastic acoustic stiffness must be positive definite"))?;
        let mut hydrostatic_compliance_pa_inv = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
        factor.solve(&mut hydrostatic_compliance_pa_inv);
        for value in &mut hydrostatic_compliance_pa_inv {
            *value /= stiffness_scale_pa;
        }
        let compressibility_pa_inv = hydrostatic_compliance_pa_inv[..3].iter().sum::<f64>();
        if hydrostatic_compliance_pa_inv
            .iter()
            .any(|value| !value.is_finite())
            || !compressibility_pa_inv.is_finite()
            || compressibility_pa_inv <= 0.0
            || !compressibility_pa_inv.recip().is_finite()
        {
            return Err(invalid("elastic hydrostatic compliance is unrepresentable"));
        }
        for (i, row) in stiffness_mandel_pa.iter().enumerate() {
            let recovered: f64 = row
                .iter()
                .zip(hydrostatic_compliance_pa_inv)
                .map(|(coefficient, compliance)| coefficient * compliance)
                .sum();
            let expected = if i < 3 { 1.0 } else { 0.0 };
            if !recovered.is_finite() || (recovered - expected).abs() > 1e-10 {
                return Err(invalid(
                    "elastic hydrostatic compliance solve did not converge",
                ));
            }
        }
        Ok(Self {
            stiffness_mandel_pa,
            density_kg_m3,
            material_state_identity,
            maximum_absolute_mandel_strain,
            hydrostatic_compliance_pa_inv,
            compressibility_pa_inv,
            stiffness_scale_pa,
        })
    }

    /// Original Mandel coefficients, without symmetry projection or rotation.
    #[must_use]
    pub const fn stiffness_mandel_pa(&self) -> &Tangent6 {
        &self.stiffness_mandel_pa
    }

    /// Original inertial density [kg/m^3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Original material-state identity, unchanged by derived queries.
    #[must_use]
    pub const fn material_state_identity(&self) -> ContentHash {
        self.material_state_identity
    }

    /// Admitted absolute ceiling for each Mandel strain component.
    #[must_use]
    pub const fn maximum_absolute_mandel_strain(&self) -> f64 {
        self.maximum_absolute_mandel_strain
    }

    /// Hydrostatic volume compressibility `I:C^-1:I` [1/Pa].
    #[must_use]
    pub const fn compressibility_pa_inv(&self) -> f64 {
        self.compressibility_pa_inv
    }

    /// Hydrostatic bulk modulus, reciprocal of volume compressibility [Pa].
    #[must_use]
    pub fn bulk_modulus_pa(&self) -> f64 {
        self.compressibility_pa_inv.recip()
    }

    /// Apply hydrostatic stress `-pressure_pa*I` within the admitted strain domain.
    pub fn hydrostatic_response(
        &self,
        pressure_pa: f64,
    ) -> Result<HydrostaticElasticResponse, MaterialError> {
        let strain_mandel = self.hydrostatic_compliance_pa_inv.map(|s| -pressure_pa * s);
        let volumetric_strain = strain_mandel[..3].iter().sum::<f64>();
        let stored_energy_density_j_m3 = -0.5 * pressure_pa * volumetric_strain;
        if !pressure_pa.is_finite()
            || strain_mandel.iter().any(|value| {
                !value.is_finite() || value.abs() > self.maximum_absolute_mandel_strain
            })
            || !volumetric_strain.is_finite()
            || !stored_energy_density_j_m3.is_finite()
            || stored_energy_density_j_m3 < 0.0
        {
            return Err(invalid(
                "hydrostatic response exceeds the elastic strain domain or is unrepresentable",
            ));
        }
        Ok(HydrostaticElasticResponse {
            material_state_identity: self.material_state_identity,
            pressure_pa,
            strain_mandel,
            volumetric_strain,
            stored_energy_density_j_m3,
        })
    }

    /// Solve the Christoffel modes for any finite nonzero direction.
    ///
    /// The existing tensor owner supplies the orthonormal Mandel basis; the
    /// existing `fs-la` spectral owner solves the fixed 3x3 symmetric problem.
    /// This is a local infinitesimal phase-wave query, not a propagation solver.
    pub fn directional_waves(
        &self,
        direction: [f64; 3],
    ) -> Result<DirectionalElasticWaves, MaterialError> {
        let scale = direction.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        if direction.iter().any(|value| !value.is_finite()) || scale <= 0.0 {
            return Err(invalid("elastic wave direction must be finite and nonzero"));
        }
        let scaled = direction.map(|value| value / scale);
        let length = det::sqrt(scaled.iter().map(|value| value * value).sum());
        let direction = scaled.map(|value| value / length);
        let basis = mandel_basis();
        let projection: [[f64; 3]; 6] = core::array::from_fn(|a| {
            core::array::from_fn(|i| (0..3).map(|j| basis[a][i][j] * direction[j]).sum())
        });
        let mut normalized = [0.0; 9];
        for i in 0..3 {
            for k in i..3 {
                let mut entry = 0.0;
                for a in 0..6 {
                    for b in 0..6 {
                        entry += projection[a][i]
                            * (self.stiffness_mandel_pa[a][b] / self.stiffness_scale_pa)
                            * projection[b][k];
                    }
                }
                normalized[3 * i + k] = entry;
                normalized[3 * k + i] = entry;
            }
        }
        let (values, vectors) = jacobi_eigh(&normalized, 3);
        let mut modes = [ElasticWaveMode {
            polarization: [0.0; 3],
            directional_modulus_pa: 0.0,
            phase_speed_m_s: 0.0,
            impedance_pa_s_m: 0.0,
        }; 3];
        for mode in 0..3 {
            let mut polarization = core::array::from_fn(|i| vectors[3 * i + mode]);
            let norm: f64 = polarization.iter().map(|value| value * value).sum();
            if !norm.is_finite() || (norm - 1.0).abs() > 1e-11 {
                return Err(invalid("elastic wave polarization is not normalized"));
            }
            for i in 0..3 {
                let residual = (0..3)
                    .map(|k| normalized[3 * i + k] * polarization[k])
                    .sum::<f64>()
                    - values[mode] * polarization[i];
                if !residual.is_finite() || residual.abs() > 1e-11 {
                    return Err(invalid("elastic wave eigenproblem did not converge"));
                }
            }
            let directional_modulus_pa = values[mode] * self.stiffness_scale_pa;
            let phase_speed_m_s = det::sqrt(directional_modulus_pa) / det::sqrt(self.density_kg_m3);
            let impedance_pa_s_m = self.density_kg_m3 * phase_speed_m_s;
            if [directional_modulus_pa, phase_speed_m_s, impedance_pa_s_m]
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            {
                return Err(invalid(
                    "elastic wave response is nonpositive or unrepresentable",
                ));
            }
            let mut pivot = 0;
            for i in 1..3 {
                if polarization[i].abs() > polarization[pivot].abs() {
                    pivot = i;
                }
            }
            if polarization[pivot] < 0.0 {
                polarization = polarization.map(|value| -value);
            }
            modes[mode] = ElasticWaveMode {
                polarization,
                directional_modulus_pa,
                phase_speed_m_s,
                impedance_pa_s_m,
            };
        }
        for i in 0..3 {
            for j in 0..i {
                let dot: f64 = (0..3)
                    .map(|k| modes[i].polarization[k] * modes[j].polarization[k])
                    .sum();
                if !dot.is_finite() || dot.abs() > 1e-11 {
                    return Err(invalid("elastic wave polarizations are not orthogonal"));
                }
            }
        }
        let christoffel_pa: [[f64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| normalized[3 * i + j] * self.stiffness_scale_pa)
        });
        if christoffel_pa
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
        {
            return Err(invalid("elastic Christoffel matrix is unrepresentable"));
        }
        Ok(DirectionalElasticWaves {
            material_state_identity: self.material_state_identity,
            direction,
            christoffel_pa,
            modes,
        })
    }
}

#[cfg(test)]
mod tests;
