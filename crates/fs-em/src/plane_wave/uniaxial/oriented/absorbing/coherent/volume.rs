//! Three-dimensional coherent dielectric work on an authored rectangular box.
//!
//! This is a control volume inside the original half-space, with no added
//! optical interface. Each face reports signed outward cycle-mean watts.
//! In particular, lateral energy transport is not folded into normal power.

use super::{
    AbsorbingCrystalCoherentIllumination, CrystalIlluminationError, OrientedUniaxialInterface,
    ReflectionError, Result, kernel,
};
use fs_math::dd::Dd;

/// Axis-aligned physical support in the original interface coordinates.
/// Origin and extent remain separate: a narrow box at a large coordinate
/// must not disappear when its upper face is rounded to binary64.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingCrystalBox {
    origin_m: [f64; 3],
    extent_m: [f64; 3],
    volume_m3: f64,
}

impl AbsorbingCrystalBox {
    pub fn try_new(origin_m: [f64; 3], extent_m: [f64; 3]) -> Result<Self> {
        if origin_m.iter().any(|v| !v.is_finite())
            || extent_m.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || origin_m[2] < 0.0
        {
            return Err(ReflectionError::Invalid(
                "crystal box needs finite origins, positive extents and nonnegative depth",
            )
            .into());
        }
        for axis in 0..3 {
            let upper = Dd::from_f64(origin_m[axis]) + Dd::from_f64(extent_m[axis]);
            let face = Dd::from_f64(extent_m[(axis + 1) % 3])
                * Dd::from_f64(extent_m[(axis + 2) % 3]);
            if !upper.to_f64().is_finite()
                || !face.to_f64().is_finite()
                || face.to_f64() <= 0.0
            {
                return Err(ReflectionError::Invalid(
                    "crystal box face or coordinate is unrepresentable",
                )
                .into());
            }
        }
        let volume_m3 = extent_m
            .iter()
            .fold(Dd::ONE, |volume, &width| volume * Dd::from_f64(width))
            .to_f64();
        if !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(ReflectionError::Invalid(
                "crystal box volume is unrepresentable",
            )
            .into());
        }
        Ok(Self {
            origin_m,
            extent_m,
            volume_m3,
        })
    }

    pub const fn origin_m(self) -> [f64; 3] {
        self.origin_m
    }
    pub const fn extent_m(self) -> [f64; 3] {
        self.extent_m
    }
    pub const fn volume_m3(self) -> f64 {
        self.volume_m3
    }
}

/// Admission before allocating or evaluating modal kernels.
#[derive(Clone, Copy, Debug)]
pub struct AbsorbingCrystalBoxLimits {
    /// Bounds (2 * component_count)^2 before equal-wavevector coalescing.
    /// Each pair has at most seven fixed analytic contractions. Short boxes
    /// instead use at most 8^3 volume and 6 * 8^2 surface field samples, each
    /// linear in that same admitted mode count. No adaptive loop is hidden.
    pub maximum_mode_pairs: usize,
    /// Peak logical modal buffers, explicit integration scratch and result
    /// headers; includes this query's returned result header. Borrowed source
    /// storage, allocator overhead and ordinary call-stack locals are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Original source-bound integrated work and all six independent face fluxes.
/// Illumination, optic axis and material remain frozen; a thermal consumer
/// supplies its actual hold interval and deposits the absorbed watts once.
pub struct AbsorbingCrystalCoherentBoxAbsorption<'i, 'p, 's, 'm> {
    illumination: &'i AbsorbingCrystalCoherentIllumination<'p, 's, 'm>,
    region: AbsorbingCrystalBox,
    absorbed_power_w: f64,
    outward_face_power_w: [f64; 6],
    residual_w: f64,
    mode_pair_count: usize,
    workspace_bytes: usize,
}

impl<'i, 'p, 's, 'm> AbsorbingCrystalCoherentBoxAbsorption<'i, 'p, 's, 'm> {
    pub const fn illumination(&self) -> &'i AbsorbingCrystalCoherentIllumination<'p, 's, 'm> {
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
    /// order. A negative value is net entry through that physical face.
    pub const fn outward_face_power_w(&self) -> [f64; 6] {
        self.outward_face_power_w
    }
    /// -sum(outward face watts) - absorbed watts. This numerical residual is
    /// observed independently and is never redistributed as material heat.
    pub const fn residual_w(&self) -> f64 {
        self.residual_w
    }
    /// Conservative admitted pair count before exact wavevector coalescing.
    pub const fn mode_pair_count(&self) -> usize {
        self.mode_pair_count
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
}

impl<'p, 's, 'm> AbsorbingCrystalCoherentIllumination<'p, 's, 'm> {
    /// Integrate original anisotropic dielectric work and signed Poynting
    /// transport over an actual 3D box. Every inter-beam/eigenmode cross term
    /// is retained; heat is not inferred by subtracting boundary powers.
    ///
    /// Ordinary boxes use analytic transverse harmonic and complex depth
    /// integrals. Short relative-phase boxes use a fixed positive eight-point
    /// Gauss rule per varying axis, summing fields before squaring to retain
    /// dark-fringe heat. Common normal propagation is integrated analytically
    /// even in opaque boxes. A conservative analytic truncation screen applies
    /// to the short-box loss integral; it is not a floating-point certificate.
    /// Unresolved cancellation refuses rather than clipping a negative value.
    ///
    /// This supplies spatially integrated cycle-mean work, not a new optical
    /// boundary, temporal exposure law or temperature feedback model.
    pub fn box_absorption(
        &self,
        region: AbsorbingCrystalBox,
        limits: AbsorbingCrystalBoxLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingCrystalCoherentBoxAbsorption<'_, 'p, 's, 'm>> {
        let integral = kernel::coherent_box_absorption(
            self.components.iter().map(|field| field.boundary_response()),
            region,
            limits,
            core::mem::size_of::<AbsorbingCrystalCoherentBoxAbsorption<'_, 'p, 's, 'm>>(),
            cancelled,
        )
        .map_err(CrystalIlluminationError::from)?;
        Ok(AbsorbingCrystalCoherentBoxAbsorption {
            illumination: self,
            region,
            absorbed_power_w: integral.absorbed_power_w,
            outward_face_power_w: integral.outward_face_power_w,
            residual_w: integral.residual_w,
            mode_pair_count: integral.mode_pair_count,
            workspace_bytes: integral.workspace_bytes,
        })
    }
}
