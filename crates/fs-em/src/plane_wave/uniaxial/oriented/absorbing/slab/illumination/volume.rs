//! Finite-film heat and signed Poynting flow on an actual rectangular support.
//!
//! All original forward/backward anchors and coherent angular cross terms
//! enter the shared native box integrator. An independent group bank then
//! sums those quadratic observations, retaining every original frequency.

use super::{
    AbsorbingSlabCoherentIllumination, AbsorbingSlabIllumination, CrystalIlluminationError,
    OrientedUniaxialSlab, ReflectionError, Result, crystal,
};
use super::super::super::{AbsorbingCrystalBox, AbsorbingCrystalBoxLimits};
use fs_math::dd::Dd;

/// Coherent four-wave and inter-beam work in the original finite film.
pub struct AbsorbingSlabCoherentBoxAbsorption<'i, 'p, 'f, 's, 'm> {
    illumination: &'i AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm>,
    region: AbsorbingCrystalBox,
    absorbed_power_w: f64,
    outward_face_power_w: [f64; 6],
    residual_w: f64,
}

impl<'i, 'p, 'f, 's, 'm> AbsorbingSlabCoherentBoxAbsorption<'i, 'p, 'f, 's, 'm> {
    pub const fn illumination(&self) -> &'i AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm> {
        self.illumination
    }
    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
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
    /// Signed outward watts at [x_min, x_max, y_min, y_max, z_min, z_max].
    /// Negative means entry. Lateral transport is not deposited as heat.
    pub const fn outward_face_power_w(&self) -> [f64; 6] {
        self.outward_face_power_w
    }
    /// -sum(outward_face_power_w) - absorbed_power_w; never redistributed.
    pub const fn residual_w(&self) -> f64 {
        self.residual_w
    }
    pub const fn mode_pair_count(&self) -> usize {
        self.illumination.mode_pair_count
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.illumination.box_workspace_bytes
    }
}

impl<'p, 'f, 's, 'm> AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm> {
    /// Integrate the complete original coherent field over this physical box.
    /// Its DD upper depth must not exceed the original film thickness; even
    /// an empty dark group cannot claim a volume outside its source.
    ///
    /// All four native modes per component retain their original incoming
    /// face anchors. Analytic pairs select a decaying endpoint; short boxes
    /// use the existing bounded positive coherent-field rule. No exponentially
    /// growing transfer factor, repeated material query or Maxwell solve is
    /// introduced. Both original incident ports remain part of each component.
    ///
    /// Work admits (4 * component_count)^2 before coalescing; workspace includes
    /// modal/observation buffers, explicit scratch and result headers. The
    /// numerical screens and no-certificate boundary of the native box kernel
    /// remain. A thermal consumer supplies its physical frozen hold interval.
    pub fn box_absorption(
        &self,
        region: AbsorbingCrystalBox,
        limits: AbsorbingCrystalBoxLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingSlabCoherentBoxAbsorption<'_, 'p, 'f, 's, 'm>> {
        let integral = crystal::coherent_slab_box_absorption(
            self.components.iter().map(|field| field.response()),
            region,
            self.source.thickness_m(),
            limits,
            core::mem::size_of::<AbsorbingSlabCoherentBoxAbsorption<'_, 'p, 'f, 's, 'm>>(),
            cancelled,
        )
        .map_err(CrystalIlluminationError::from)?;
        Ok(AbsorbingSlabCoherentBoxAbsorption {
            illumination: self,
            region,
            absorbed_power_w: integral.absorbed_power_w,
            outward_face_power_w: integral.outward_face_power_w,
            residual_w: integral.residual_w,
        })
    }
}

/// Expected finite-volume heat and six signed powers from independent groups.
/// The original material/frequencies and complete coherent factors stay borrowed.
pub struct AbsorbingSlabIlluminationBoxAbsorption<'i, 'g, 'p, 'f, 's, 'm> {
    illumination: &'i AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm>,
    region: AbsorbingCrystalBox,
    absorbed_power_w: f64,
    outward_face_power_w: [f64; 6],
    residual_w: f64,
}

impl<'i, 'g, 'p, 'f, 's, 'm> AbsorbingSlabIlluminationBoxAbsorption<'i, 'g, 'p, 'f, 's, 'm> {
    pub const fn illumination(&self) -> &'i AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm> {
        self.illumination
    }
    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
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
    /// Signed outward watts at [x_min, x_max, y_min, y_max, z_min, z_max].
    pub const fn outward_face_power_w(&self) -> [f64; 6] {
        self.outward_face_power_w
    }
    /// -sum(returned face watts) - returned heat; never deposited as heat.
    pub const fn residual_w(&self) -> f64 {
        self.residual_w
    }
    pub const fn mode_pair_count(&self) -> usize {
        self.illumination.mode_pair_count
    }
    /// Own header and seven DD accumulators plus the largest coherent-group
    /// workspace. Borrowed sources, allocator overhead and stack locals are
    /// excluded; the same exact admission value applies to later boxes.
    pub const fn workspace_bytes(&self) -> usize {
        self.illumination.box_workspace_bytes
    }
}

impl<'g, 'p, 'f, 's, 'm> AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm> {
    /// Observe the same physical film box in each independent coherent group.
    /// Each group retains its own original frequency-dependent loss tensor;
    /// fields interfere only within that group, even for equal frequencies.
    /// Statistical square-root weights already belong to the original fields.
    ///
    /// Limits admit the SUM of per-group four-mode pair counts and PEAK
    /// workspace before integration. Groups are evaluated sequentially with
    /// O(group_count + sum(mode_pairs)) bounded work, including each group's
    /// fixed short-box rule. The original constructor bounds empty-group
    /// traversal. Any refusal or cancellation publishes no partial result.
    pub fn box_absorption(
        &self,
        region: AbsorbingCrystalBox,
        limits: AbsorbingCrystalBoxLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingSlabIlluminationBoxAbsorption<'_, 'g, 'p, 'f, 's, 'm>> {
        crystal::poll_illumination(&mut cancelled)?;
        crystal::admit_slab_box(region, self.source.thickness_m())?;
        if self.mode_pair_count > limits.maximum_mode_pairs
            || self.box_workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(ReflectionError::Invalid(
                "film illumination box exceeds its total pair work or peak workspace budget",
            ).into());
        }
        let group_limits = AbsorbingCrystalBoxLimits {
            maximum_mode_pairs: self.mode_pair_count,
            maximum_workspace_bytes: self.group_workspace_bytes,
        };
        let mut sums = [Dd::ZERO; 7];
        for group in self.groups {
            crystal::poll_illumination(&mut cancelled)?;
            let result = group.box_absorption(region, group_limits, &mut cancelled)?;
            sums[0] = sums[0] + Dd::from_f64(result.absorbed_power_w());
            for (sum, value) in sums[1..].iter_mut().zip(result.outward_face_power_w()) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let [absorbed_power_w, x_min, x_max, y_min, y_max, z_min, z_max] =
            sums.map(Dd::to_f64);
        let outward_face_power_w = [x_min, x_max, y_min, y_max, z_min, z_max];
        let residual_w = (-outward_face_power_w.iter()
            .fold(Dd::ZERO, |sum, &value| sum + Dd::from_f64(value))
            - Dd::from_f64(absorbed_power_w)).to_f64();
        if !absorbed_power_w.is_finite() || absorbed_power_w < 0.0
            || outward_face_power_w.iter().any(|value| !value.is_finite())
            || !residual_w.is_finite()
        {
            return Err(ReflectionError::Invalid(
                "film illumination box power is unrepresentable",
            ).into());
        }
        crystal::poll_illumination(&mut cancelled)?;
        Ok(AbsorbingSlabIlluminationBoxAbsorption {
            illumination: self,
            region,
            absorbed_power_w,
            outward_face_power_w,
            residual_w,
        })
    }
}
