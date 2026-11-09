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

use crate::{MaterialError, Tangent6, Voigt, tensor::mandel_basis};
use fs_blake3::ContentHash;
use fs_la::{eigen::jacobi_eigh, factor::cholesky};
use fs_math::det;

/// Source-selected density, bulk modulus and shear modulus for an explicitly
/// isotropic elastic law, retaining the original property-use receipts.
pub mod sourced;

/// Coherent P/SV/SH reflection at a perfectly bonded isotropic interface.
pub mod interface;

/// Intrinsic absorption and coherent normal incidence on a causal solid slab.
pub mod viscoelastic;

mod constitutive;
pub use constitutive::ElasticStrainResponse;

fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters {
        what: what.to_owned(),
    }
}

// A native law's objective strain ball must not become a larger component
// box merely because the acoustic adapter uses Mandel coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ElasticStrainDomain {
    MandelComponents,
    TensorFrobenius,
}

/// Frozen elastic state shared by wave and pressure queries.
#[derive(Clone, Debug, PartialEq)]
pub struct ElasticAcousticState {
    stiffness_mandel_pa: Tangent6,
    density_kg_m3: f64,
    // A rounded dense tensor cannot retain both modes when K/G is extreme.
    // These are the original independent moduli, never inferred from a tensor.
    isotropic_moduli_pa: Option<(f64, f64)>,
    material_state_identity: ContentHash,
    maximum_absolute_mandel_strain: f64,
    strain_domain: ElasticStrainDomain,
    hydrostatic_compliance_pa_inv: [f64; 6],
    compressibility_pa_inv: f64,
    stiffness_scale_pa: f64,
}

/// One Christoffel eigenmode with phase and energy-transport velocities.
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
    /// Energy flux divided by total time-averaged wave energy density [m/s].
    /// This equals `grad_k omega` on a smooth, separated acoustic branch and
    /// satisfies `n.dot(energy_velocity) = phase_speed`. It need not be parallel
    /// to the wave normal. At a repeated eigenvalue it belongs to the returned
    /// polarization; a unique differentiable branch is not implied.
    pub energy_velocity_m_s: [f64; 3],
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
            return Err(invalid(
                "isotropic bulk and shear moduli must be finite and positive",
            ));
        }
        Self::from_isotropic_parts(
            super::bulk_shear_stiffness(bulk_modulus_pa, shear_modulus_pa),
            bulk_modulus_pa,
            shear_modulus_pa,
            density_kg_m3,
            material_state_identity,
            maximum_absolute_mandel_strain,
        )
    }

    /// Use an existing isotropic law and preserve its tensor Frobenius strain ball.
    pub fn from_isotropic(
        law: &super::IsotropicElastic,
        density_kg_m3: f64,
        material_state_identity: ContentHash,
    ) -> Result<Self, MaterialError> {
        // The public law may have been changed after construction.
        super::IsotropicElastic::new(law.youngs, law.poisson, law.strain_limit)?;
        let mut state = Self::from_isotropic_parts(
            law.stiffness(),
            law.youngs / (3.0 * (1.0 - 2.0 * law.poisson)),
            law.youngs / (2.0 * (1.0 + law.poisson)),
            density_kg_m3,
            material_state_identity,
            law.strain_limit,
        )?;
        state.strain_domain = ElasticStrainDomain::TensorFrobenius;
        Ok(state)
    }

    pub(in crate::elastic) fn from_isotropic_parts(
        stiffness_mandel_pa: Tangent6,
        bulk: f64,
        shear: f64,
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
        // K>0 and G>0 are the exact isotropic stability conditions. Factoring
        // the dense representation instead would reject stable soft shear or
        // volume modes that are below its normal-block rounding resolution.
        let compressibility_pa_inv = bulk.recip();
        if [bulk, shear, compressibility_pa_inv]
            .into_iter()
            .any(|value| !value.is_finite() || value <= 0.0)
            || stiffness_mandel_pa
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid("isotropic elastic response is unrepresentable"));
        }
        let stiffness_scale_pa = stiffness_mandel_pa
            .iter()
            .flatten()
            .fold(0.0_f64, |scale, value| scale.max(value.abs()));
        let normal_compliance = compressibility_pa_inv / 3.0;
        if normal_compliance == 0.0 {
            return Err(invalid(
                "isotropic hydrostatic compliance is unrepresentable",
            ));
        }
        Ok(Self {
            stiffness_mandel_pa,
            density_kg_m3,
            isotropic_moduli_pa: Some((bulk, shear)),
            material_state_identity,
            maximum_absolute_mandel_strain,
            strain_domain: ElasticStrainDomain::MandelComponents,
            hydrostatic_compliance_pa_inv: [
                normal_compliance,
                normal_compliance,
                normal_compliance,
                0.0,
                0.0,
                0.0,
            ],
            compressibility_pa_inv,
            stiffness_scale_pa,
        })
    }

    /// Use an existing orthotropic law in its principal frame, preserving its
    /// tensor Frobenius strain ball. Oriented consumers passing a transformed
    /// tensor to `try_new_mandel` must declare that tensor's strain domain.
    pub fn from_orthotropic(
        law: &super::OrthotropicElastic,
        density_kg_m3: f64,
        material_state_identity: ContentHash,
    ) -> Result<Self, MaterialError> {
        let mut state = Self::try_new_mandel(
            law.stiffness(),
            density_kg_m3,
            material_state_identity,
            law.strain_limit,
        )?;
        state.strain_domain = ElasticStrainDomain::TensorFrobenius;
        Ok(state)
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
            isotropic_moduli_pa: None,
            material_state_identity,
            maximum_absolute_mandel_strain,
            strain_domain: ElasticStrainDomain::MandelComponents,
            hydrostatic_compliance_pa_inv,
            compressibility_pa_inv,
            stiffness_scale_pa,
        })
    }

    /// Original Mandel coefficients, without symmetry projection or rotation.
    ///
    /// At extreme isotropic stiffness contrasts, rounding of this dense tensor
    /// can hide the weaker bulk or shear mode. The original moduli accessor and
    /// the constitutive, pressure and wave queries retain both modes; arbitrary
    /// dense tensor consumers still need to account for this conditioning.
    #[must_use]
    pub const fn stiffness_mandel_pa(&self) -> &Tangent6 {
        &self.stiffness_mandel_pa
    }

    /// Original inertial density [kg/m^3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Whether the original constructor explicitly selected an isotropic law.
    /// Generic Mandel, cubic and orthotropic admission never infer isotropy
    /// from approximately equal tensor coefficients or discard anisotropy.
    #[must_use]
    pub const fn is_explicitly_isotropic(&self) -> bool {
        self.isotropic_moduli_pa.is_some()
    }

    /// Original independent `(bulk, shear)` moduli [Pa] for explicit isotropy.
    /// Generic tensor, cubic and orthotropic constructors return `None`, even
    /// if their rounded coefficients happen to resemble an isotropic tensor.
    /// Consumers must not reconstruct a weak bulk mode by summing tensor rows.
    #[must_use]
    pub const fn isotropic_bulk_shear_pa(&self) -> Option<(f64, f64)> {
        self.isotropic_moduli_pa
    }

    /// Original material-state identity, unchanged by derived queries.
    #[must_use]
    pub const fn material_state_identity(&self) -> ContentHash {
        self.material_state_identity
    }

    /// Absolute ceiling for each Mandel strain component.
    ///
    /// Native isotropic/orthotropic adapters also retain their Frobenius norm
    /// bound: satisfying each component ceiling alone does not admit a strain.
    /// A Mandel norm no greater than this value is sufficient for either shape.
    /// Use `admits_mandel_strain` for the complete original domain.
    #[must_use]
    pub const fn maximum_absolute_mandel_strain(&self) -> f64 {
        self.maximum_absolute_mandel_strain
    }

    /// Original tensor Frobenius strain limit for a native elastic-law adapter.
    ///
    /// The Mandel Euclidean norm equals the physical tensor Frobenius norm.
    /// Explicit bulk/shear, cubic and generic Mandel constructors declare
    /// componentwise boxes and return `None`; no norm bound is inferred.
    #[must_use]
    pub const fn frobenius_strain_limit(&self) -> Option<f64> {
        match self.strain_domain {
            ElasticStrainDomain::TensorFrobenius => Some(self.maximum_absolute_mandel_strain),
            ElasticStrainDomain::MandelComponents => None,
        }
    }

    /// Admit a finite Mandel strain against the original domain, inclusively.
    /// Native elastic laws retain their objective Frobenius ball; explicitly
    /// declared componentwise domains retain their box. No allocation occurs.
    #[must_use]
    pub fn admits_mandel_strain(&self, strain: &Voigt) -> bool {
        if strain.iter().any(|value| !value.is_finite()) {
            return false;
        }
        let scale = strain.iter().fold(0.0_f64, |s, value| s.max(value.abs()));
        if scale > self.maximum_absolute_mandel_strain {
            return false;
        }
        match self.strain_domain {
            ElasticStrainDomain::MandelComponents => true,
            ElasticStrainDomain::TensorFrobenius => {
                if scale == 0.0 {
                    return true;
                }
                // Normalize before squaring: tiny/large finite calibrated
                // bounds must not underflow/overflow into a different domain.
                let squared_norm = strain.iter().map(|value| (value / scale).powi(2)).sum();
                scale * det::sqrt(squared_norm) <= self.maximum_absolute_mandel_strain
            }
        }
    }

    /// Hydrostatic volume compressibility `I:C^-1:I` [1/Pa].
    #[must_use]
    pub const fn compressibility_pa_inv(&self) -> f64 {
        self.compressibility_pa_inv
    }

    /// Hydrostatic bulk modulus, reciprocal of volume compressibility [Pa].
    #[must_use]
    pub fn bulk_modulus_pa(&self) -> f64 {
        self.isotropic_moduli_pa
            .map_or_else(|| self.compressibility_pa_inv.recip(), |(bulk, _)| bulk)
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
            || !self.admits_mandel_strain(&strain_mandel)
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
    /// Explicitly isotropic states use the analytical longitudinal and two
    /// transverse modes from their original bulk and shear moduli. For general
    /// tensors, the orthonormal Mandel basis and existing `fs-la` spectral owner
    /// solve the fixed 3x3 symmetric problem. Energy transport is analytical for
    /// isotropy and uses the tensor contraction otherwise, without angular
    /// differencing or additional eigenproblems. This is a local infinitesimal
    /// wave query, not a propagation solver.
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
        if let Some((bulk, shear)) = self.isotropic_moduli_pa {
            return self.isotropic_waves(direction, bulk, shear);
        }
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
            energy_velocity_m_s: [0.0; 3],
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
            // For unit displacement-gradient amplitude, s = C:sym(a tensor n).
            // The mean mechanical flux divided by total mean wave energy is
            // v_j = a_i s_ij/(rho*c). Normalize C first, and use rho*c^2/Cscale
            // from the SAME eigenpair, avoiding a large dimensional C*a product.
            let strain: [f64; 6] =
                core::array::from_fn(|a| (0..3).map(|i| projection[a][i] * polarization[i]).sum());
            let stress: [f64; 6] = core::array::from_fn(|a| {
                (0..6)
                    .map(|b| (self.stiffness_mandel_pa[a][b] / self.stiffness_scale_pa) * strain[b])
                    .sum()
            });
            let energy_velocity_m_s = core::array::from_fn(|j| {
                let flux: f64 = (0..6)
                    .map(|a| {
                        stress[a]
                            * (0..3)
                                .map(|i| basis[a][i][j] * polarization[i])
                                .sum::<f64>()
                    })
                    .sum();
                (flux / values[mode]) * phase_speed_m_s
            });
            if energy_velocity_m_s.iter().any(|value| !value.is_finite()) {
                return Err(invalid("elastic wave energy velocity is unrepresentable"));
            }
            modes[mode] = ElasticWaveMode {
                polarization,
                directional_modulus_pa,
                phase_speed_m_s,
                impedance_pa_s_m,
                energy_velocity_m_s,
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
    // Gamma = G I + (K + G/3) n n^T. Its two transverse eigenvalues
    // are G and its longitudinal eigenvalue is K + 4G/3. Keep these modes
    // directly: no dense eigensolve can recover G after it is rounded out of
    // Gamma. The deterministic transverse basis has no unique physical axis.
    fn isotropic_waves(
        &self,
        direction: [f64; 3],
        bulk: f64,
        shear: f64,
    ) -> Result<DirectionalElasticWaves, MaterialError> {
        let axis = (1..3).fold(0, |smallest, i| {
            if direction[i].abs() < direction[smallest].abs() {
                i
            } else {
                smallest
            }
        });
        let mut transverse = [0.0; 3];
        transverse[(axis + 1) % 3] = direction[(axis + 2) % 3];
        transverse[(axis + 2) % 3] = -direction[(axis + 1) % 3];
        let norm = det::sqrt(transverse.iter().map(|value| value * value).sum());
        transverse = transverse.map(|value| value / norm);
        let second = core::array::from_fn(|i| {
            direction[(i + 1) % 3] * transverse[(i + 2) % 3]
                - direction[(i + 2) % 3] * transverse[(i + 1) % 3]
        });
        let longitudinal = bulk + (4.0 / 3.0) * shear;
        let mut modes = [ElasticWaveMode {
            polarization: [0.0; 3],
            directional_modulus_pa: 0.0,
            phase_speed_m_s: 0.0,
            impedance_pa_s_m: 0.0,
            energy_velocity_m_s: [0.0; 3],
        }; 3];
        for (index, mut polarization) in [transverse, second, direction].into_iter().enumerate() {
            let directional_modulus_pa = if index < 2 { shear } else { longitudinal };
            let phase_speed_m_s = det::sqrt(directional_modulus_pa) / det::sqrt(self.density_kg_m3);
            let impedance_pa_s_m = self.density_kg_m3 * phase_speed_m_s;
            if [directional_modulus_pa, phase_speed_m_s, impedance_pa_s_m]
                .into_iter()
                .any(|value| !value.is_finite() || value <= 0.0)
            {
                return Err(invalid(
                    "elastic wave response is nonpositive or unrepresentable",
                ));
            }
            let pivot = (1..3).fold(0, |largest, i| {
                if polarization[i].abs() > polarization[largest].abs() {
                    i
                } else {
                    largest
                }
            });
            if polarization[pivot] < 0.0 {
                polarization = polarization.map(|value| -value);
            }
            modes[index] = ElasticWaveMode {
                polarization,
                directional_modulus_pa,
                phase_speed_m_s,
                impedance_pa_s_m,
                energy_velocity_m_s: direction.map(|value| phase_speed_m_s * value),
            };
        }
        let christoffel_pa = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (bulk + shear / 3.0).mul_add(
                    direction[i] * direction[j],
                    if i == j { shear } else { 0.0 },
                )
            })
        });
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
