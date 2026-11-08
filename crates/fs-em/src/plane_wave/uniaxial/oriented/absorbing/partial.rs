//! Explicit partial spatial coherence on the original absorbing crystal.
//!
//! Each borrowed factor is a coherent angular illumination, including its
//! square-root statistical weight in the original peak fields. Distinct
//! factors are declared mutually uncorrelated, even when their frequencies
//! coincide. Their mean work and signed energy flux add; their phasors do not.
//! This admits a positive factorization of spatial/polarization coherence
//! without a random-phase draw, an inferred bandwidth or an implicit PSD.
//! Deterministic inter-frequency beats need an explicit temporal source model.

use super::{
    AbsorbingCrystalBox, AbsorbingCrystalBoxLimits, AbsorbingCrystalCoherentBoxAbsorption,
    AbsorbingCrystalCoherentIllumination, CrystalIlluminationError, OrientedUniaxialInterface,
    ReflectionError, kernel,
};
use fs_math::dd::Dd;

type Result<T> = core::result::Result<T, CrystalIlluminationError>;

/// An explicit bank of independent coherent angular factors.
/// The original interface must be identical for every factor; each factor
/// retains its own admitted monochromatic frequency, spectrum and loss tensor.
/// Construction borrows everything and takes O(factors) work without allocation.
pub struct AbsorbingCrystalPartiallyCoherentIllumination<'f, 'p, 's, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    factors: &'f [AbsorbingCrystalCoherentIllumination<'p, 's, 'm>],
    mode_pair_count: usize,
    factor_workspace_bytes: usize,
    workspace_bytes: usize,
}

/// Ensemble-mean dielectric work and six independent signed face powers.
/// It retains the complete factorization and original physical box, not a
/// fictitious ensemble phasor. A thermal consumer owns the frozen hold interval.
pub struct AbsorbingCrystalPartiallyCoherentBoxAbsorption<'i, 'f, 'p, 's, 'm> {
    illumination: &'i AbsorbingCrystalPartiallyCoherentIllumination<'f, 'p, 's, 'm>,
    region: AbsorbingCrystalBox,
    absorbed_power_w: f64,
    outward_face_power_w: [f64; 6],
    residual_w: f64,
}

impl<'f, 'p, 's, 'm> AbsorbingCrystalPartiallyCoherentIllumination<'f, 'p, 's, 'm> {
    /// Bind explicitly uncorrelated factors to their exact original interface.
    /// Frequencies may differ BETWEEN factors. Within each factor the coherent
    /// owner has already required exactly one frequency and retained all
    /// inter-beam and ordinary/extraordinary interference.
    ///
    /// Amplitudes already include square-root statistical weights: no second
    /// weighting, bandwidth factor or peak/RMS conversion is applied here.
    /// The source supplies no implicit illumination. Empty factors and an
    /// empty bank are dark; the factor budget also bounds iteration over them.
    pub fn try_new(
        source: &'s OrientedUniaxialInterface<'m>,
        factors: &'f [AbsorbingCrystalCoherentIllumination<'p, 's, 'm>],
        maximum_factors: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        kernel::poll_illumination(&mut cancelled)?;
        if factors.len() > maximum_factors {
            return Err(ReflectionError::Invalid(
                "partially coherent crystal illumination exceeds its factor budget",
            )
            .into());
        }
        let mut mode_pair_count = 0usize;
        let mut factor_workspace_bytes = 0usize;
        for factor in factors {
            kernel::poll_illumination(&mut cancelled)?;
            if !core::ptr::eq(source, factor.source()) {
                return Err(ReflectionError::Invalid(
                    "partially coherent crystal factors require the exact original interface",
                )
                .into());
            }
            let (pairs, workspace) = kernel::coherent_box_requirements(
                factor.component_count(),
                core::mem::size_of::<AbsorbingCrystalCoherentBoxAbsorption<'_, 'p, 's, 'm>>(),
            )?;
            mode_pair_count = mode_pair_count.checked_add(pairs).ok_or(
                ReflectionError::Invalid("partially coherent crystal pair count overflows"),
            )?;
            factor_workspace_bytes = factor_workspace_bytes.max(workspace);
        }
        let workspace_bytes = factor_workspace_bytes
            .checked_add(core::mem::size_of::<
                AbsorbingCrystalPartiallyCoherentBoxAbsorption<'_, 'f, 'p, 's, 'm>,
            >())
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<[Dd; 7]>()))
            .ok_or(ReflectionError::Invalid(
                "partially coherent crystal workspace overflows",
            ))?;
        kernel::poll_illumination(&mut cancelled)?;
        Ok(Self {
            source,
            factors,
            mode_pair_count,
            factor_workspace_bytes,
            workspace_bytes,
        })
    }

    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.source
    }

    pub const fn factors(&self) -> &'f [AbsorbingCrystalCoherentIllumination<'p, 's, 'm>] {
        self.factors
    }

    pub const fn factor_count(&self) -> usize {
        self.factors.len()
    }

    /// Checked sum of each factor's (2 * component_count)^2, before any
    /// equal-wavevector coalescing. Independent factors have no cross pairs.
    pub const fn mode_pair_count(&self) -> usize {
        self.mode_pair_count
    }

    /// Integrate every original coherent factor over this SAME physical box
    /// and add its mean dielectric work and six signed Poynting face powers.
    /// Each factor uses its own original frequency-dependent loss tensor.
    ///
    /// The pair limit admits the SUM across factors before any integration.
    /// Workspace admits the aggregate result and seven DD accumulators plus
    /// the LARGEST coherent-factor workspace, including its result headers.
    /// Factors are evaluated sequentially; no array of factor results is kept.
    /// Borrowed source storage, allocator overhead and ordinary stack locals
    /// are excluded, as in the coherent query. Work is bounded by factor_count
    /// plus the sum of admitted pairs and each factor's fixed short-box rule.
    ///
    /// Cancellation or any factor refusal publishes no partial result.
    /// Equal-frequency factors remain independent. There is no additional
    /// optical boundary, inferred coherence time or temporal beat averaging.
    pub fn box_absorption(
        &self,
        region: AbsorbingCrystalBox,
        limits: AbsorbingCrystalBoxLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingCrystalPartiallyCoherentBoxAbsorption<'_, 'f, 'p, 's, 'm>> {
        kernel::poll_illumination(&mut cancelled)?;
        if self.mode_pair_count > limits.maximum_mode_pairs
            || self.workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(ReflectionError::Invalid(
                "partially coherent crystal box exceeds its total pair work or peak workspace budget",
            )
            .into());
        }
        let factor_limits = AbsorbingCrystalBoxLimits {
            maximum_mode_pairs: self.mode_pair_count,
            maximum_workspace_bytes: self.factor_workspace_bytes,
        };
        let mut sums = [Dd::ZERO; 7];
        for factor in self.factors {
            kernel::poll_illumination(&mut cancelled)?;
            let result = factor.box_absorption(region, factor_limits, &mut cancelled)?;
            sums[0] = sums[0] + Dd::from_f64(result.absorbed_power_w());
            for (sum, value) in sums[1..].iter_mut().zip(result.outward_face_power_w()) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let [absorbed_power_w, x_min, x_max, y_min, y_max, z_min, z_max] =
            sums.map(Dd::to_f64);
        let outward_face_power_w = [x_min, x_max, y_min, y_max, z_min, z_max];
        let residual_w = (-outward_face_power_w
            .iter()
            .fold(Dd::ZERO, |sum, &value| sum + Dd::from_f64(value))
            - Dd::from_f64(absorbed_power_w))
        .to_f64();
        if !absorbed_power_w.is_finite()
            || absorbed_power_w < 0.0
            || outward_face_power_w.iter().any(|value| !value.is_finite())
            || !residual_w.is_finite()
        {
            return Err(ReflectionError::Invalid(
                "partially coherent crystal box power is unrepresentable",
            )
            .into());
        }
        kernel::poll_illumination(&mut cancelled)?;
        Ok(AbsorbingCrystalPartiallyCoherentBoxAbsorption {
            illumination: self,
            region,
            absorbed_power_w,
            outward_face_power_w,
            residual_w,
        })
    }
}

impl<'i, 'f, 'p, 's, 'm> AbsorbingCrystalPartiallyCoherentBoxAbsorption<'i, 'f, 'p, 's, 'm> {
    pub const fn illumination(
        &self,
    ) -> &'i AbsorbingCrystalPartiallyCoherentIllumination<'f, 'p, 's, 'm> {
        self.illumination
    }
    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.illumination.source()
    }
    pub const fn region(&self) -> AbsorbingCrystalBox {
        self.region
    }
    pub const fn volume_m3(&self) -> f64 {
        self.region.volume_m3()
    }
    pub const fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
    /// Signed outward watts in [x_min, x_max, y_min, y_max, z_min, z_max]
    /// order. Negative power means net entry through that physical face.
    pub const fn outward_face_power_w(&self) -> [f64; 6] {
        self.outward_face_power_w
    }
    /// -sum(returned outward face watts) - returned absorbed watts. This
    /// independent numerical residual is never reassigned to material heat.
    pub const fn residual_w(&self) -> f64 {
        self.residual_w
    }
    pub const fn mode_pair_count(&self) -> usize {
        self.illumination.mode_pair_count
    }
    /// Aggregate header/accumulators plus peak coherent-factor scratch and
    /// headers. This exact admission value can be reused on later boxes.
    pub const fn workspace_bytes(&self) -> usize {
        self.illumination.workspace_bytes
    }
}
