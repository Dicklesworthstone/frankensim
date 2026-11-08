//! Coherent mode groups and independent positive factors in a finite crystal.
//!
//! Original prepared fields carry their authored global peak amplitudes,
//! including square-root statistical weights. A group explicitly declares
//! mutual coherence at exactly one frequency. Independent groups add second
//! moments and tensor heat, never phasors; different groups may have different
//! frequencies. These are finite spectral-line powers, not an inferred PSD,
//! bandwidth, temporal beat envelope or degree of polarization.
//!
//! Both owners borrow the existing solved fields and allocate nothing. All
//! observations reuse the original faces, material tensor and anchored waves.

use super::super::{AbsorbingCrystalStatistics, CrystalIlluminationError};
use super::{
    AbsorbingCrystalField, FieldSum, OrientedCrystalField, OrientedUniaxialSlab,
    PreparedAbsorbingSlabField, ReflectionError, SlabRegion, admit_point, crystal,
};
use fs_math::dd::Dd;

mod volume;
pub use volume::{AbsorbingSlabCoherentBoxAbsorption, AbsorbingSlabIlluminationBoxAbsorption};

type Result<T> = core::result::Result<T, CrystalIlluminationError>;

/// Explicitly mutually coherent prepared beams from the exact original film.
/// Incidence directions and polarizations may differ. Relative phases share
/// the original global origin; all four internal eigenwaves of every beam
/// interfere before Poynting flux or dielectric work is evaluated.
pub struct AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm> {
    source: &'f OrientedUniaxialSlab<'s, 'm>,
    components: &'p [PreparedAbsorbingSlabField<'f, 's, 'm>],
    mode_pair_count: usize,
    box_workspace_bytes: usize,
}

impl<'p, 'f, 's, 'm> AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm> {
    /// Bind a borrowed monochromatic group without repeating material or
    /// boundary solves. Source identity and frequency must agree exactly;
    /// even one frequency ULP of difference refuses coherent composition.
    /// Empty components mean darkness and carry no implicit carrier.
    pub fn try_new(
        source: &'f OrientedUniaxialSlab<'s, 'm>,
        components: &'p [PreparedAbsorbingSlabField<'f, 's, 'm>],
        maximum_components: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        crystal::poll_illumination(&mut cancelled)?;
        if components.len() > maximum_components {
            return Err(ReflectionError::Invalid(
                "coherent film illumination exceeds its component budget",
            )
            .into());
        }
        let (mode_pair_count, box_workspace_bytes) = crystal::slab_box_requirements(
            components.len(),
            core::mem::size_of::<AbsorbingSlabCoherentBoxAbsorption<'_, 'p, 'f, 's, 'm>>(),
        )?;
        let omega = components
            .first()
            .map(|field| field.response().omega_rad_s());
        for component in components {
            crystal::poll_illumination(&mut cancelled)?;
            if !core::ptr::eq(source, component.source()) {
                return Err(ReflectionError::Invalid(
                    "coherent film components require the exact original slab",
                )
                .into());
            }
            if Some(component.response().omega_rad_s()) != omega {
                return Err(ReflectionError::Invalid(
                    "coherent film components require exactly the same frequency",
                )
                .into());
            }
        }
        crystal::poll_illumination(&mut cancelled)?;
        Ok(Self { source, components, mode_pair_count, box_workspace_bytes })
    }

    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
        self.source
    }
    pub const fn components(&self) -> &'p [PreparedAbsorbingSlabField<'f, 's, 'm>] {
        self.components
    }
    pub const fn component_count(&self) -> usize {
        self.components.len()
    }
    /// Checked four-wave box work before exact anchored-mode coalescing.
    pub const fn mode_pair_count(&self) -> usize {
        self.mode_pair_count
    }
    pub fn omega_rad_s(&self) -> Option<f64> {
        self.components
            .first()
            .map(|field| field.response().omega_rad_s())
    }

    fn electromagnetic_field(
        &self,
        region: SlabRegion,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<OrientedCrystalField> {
        crystal::poll_illumination(cancelled)?;
        admit_point(region, position_m, self.source.thickness_m())?;
        let mut sum = FieldSum::at_point(position_m).map_err(ReflectionError::Invalid)?;
        for component in self.components {
            crystal::poll_illumination(cancelled)?;
            component.response().add_waves(region, &mut sum)?;
        }
        let field = sum.finish().map_err(ReflectionError::Invalid)?;
        crystal::poll_illumination(cancelled)?;
        Ok(field)
    }

    /// Coherent global peak E/H, cycle-mean Poynting vector and local W/m³.
    /// Crystal heat contracts the common original loss tensor with the TOTAL
    /// field. Spatial fringes are retained, including dark interference nodes.
    /// A volume consumer must integrate its physical transverse extent;
    /// multiplying one point by an area is invalid for crossed-beam fringes.
    pub fn field_at(
        &self,
        region: SlabRegion,
        position_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingCrystalField> {
        let field = self.electromagnetic_field(region, position_m, &mut cancelled)?;
        let result = match self.components.first() {
            Some(reference) => reference.response().with_heating(region, field)?,
            None => AbsorbingCrystalField {
                field,
                absorbed_power_density_w_m3: 0.0,
            },
        };
        crystal::poll_illumination(&mut cancelled)?;
        Ok(result)
    }

    /// Signed energy flux without requiring representable electric squared
    /// moments or recomputing dielectric work. No per-beam fractions are added.
    pub fn poynting_w_m2(
        &self,
        region: SlabRegion,
        position_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3]> {
        self.electromagnetic_field(region, position_m, &mut cancelled)
            .map(|field| field.poynting_w_m2)
    }
}

/// Mutually uncorrelated coherent groups on the same original finite film.
/// A single prepared plane wave per group reproduces the half-space owner's
/// independent-factor convention. Multiple beams within a group preserve
/// partial spatial coherence as well as arbitrary elliptical polarization.
/// Preparation is O(groups); observation is O(groups + total components),
/// using fixed-size compensated accumulators and no allocation.
pub struct AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm> {
    source: &'f OrientedUniaxialSlab<'s, 'm>,
    groups: &'g [AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm>],
    component_count: usize,
    mode_pair_count: usize,
    group_workspace_bytes: usize,
    box_workspace_bytes: usize,
}

impl<'g, 'p, 'f, 's, 'm> AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm> {
    /// Declare independent positive factors, including any authored statistical
    /// weights in their prepared peak fields. No coherence is inferred from
    /// equal frequency, scalar intensity, direction or source proximity.
    pub fn try_new(
        source: &'f OrientedUniaxialSlab<'s, 'm>,
        groups: &'g [AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm>],
        maximum_groups: usize,
        maximum_components: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        crystal::poll_illumination(&mut cancelled)?;
        if groups.len() > maximum_groups {
            return Err(ReflectionError::Invalid(
                "film illumination exceeds its independent-group budget",
            )
            .into());
        }
        let mut component_count = 0usize;
        let mut mode_pair_count = 0usize;
        let mut group_workspace_bytes = 0usize;
        for group in groups {
            crystal::poll_illumination(&mut cancelled)?;
            if !core::ptr::eq(source, group.source()) {
                return Err(ReflectionError::Invalid(
                    "film illumination groups require the exact original slab",
                )
                .into());
            }
            component_count = component_count
                .checked_add(group.component_count())
                .filter(|&count| count <= maximum_components)
                .ok_or(ReflectionError::Invalid(
                    "film illumination exceeds its total component budget",
                ))?;
            mode_pair_count = mode_pair_count.checked_add(group.mode_pair_count()).ok_or(
                ReflectionError::Invalid("film illumination box pair count overflows"),
            )?;
            group_workspace_bytes = group_workspace_bytes.max(group.box_workspace_bytes);
        }
        let box_workspace_bytes = group_workspace_bytes
            .checked_add(core::mem::size_of::<
                AbsorbingSlabIlluminationBoxAbsorption<'_, 'g, 'p, 'f, 's, 'm>,
            >())
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<[Dd; 7]>()))
            .ok_or(ReflectionError::Invalid("film illumination box workspace overflows"))?;
        crystal::poll_illumination(&mut cancelled)?;
        Ok(Self {
            source,
            groups,
            component_count,
            mode_pair_count,
            group_workspace_bytes,
            box_workspace_bytes,
        })
    }

    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
        self.source
    }
    pub const fn groups(&self) -> &'g [AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm>] {
        self.groups
    }
    pub const fn group_count(&self) -> usize {
        self.groups.len()
    }
    pub const fn component_count(&self) -> usize {
        self.component_count
    }
    /// Sum of each independent group's checked four-wave box pair count.
    pub const fn mode_pair_count(&self) -> usize {
        self.mode_pair_count
    }

    /// Cartesian electric coherency, magnetic RMS², energy flux and tensor
    /// heat. Peak-to-cycle-mean factors are applied once by the existing
    /// statistics owner. No mean complex field is invented for the ensemble.
    /// Region and coordinate validation still applies to empty illumination.
    pub fn field_statistics(
        &self,
        region: SlabRegion,
        position_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingCrystalStatistics> {
        crystal::poll_illumination(&mut cancelled)?;
        admit_point(region, position_m, self.source.thickness_m())?;
        let mut sum = crystal::StatisticsSum::new();
        for group in self.groups {
            crystal::poll_illumination(&mut cancelled)?;
            sum.add(group.field_at(region, position_m, &mut cancelled)?);
        }
        let result = sum.finish()?;
        crystal::poll_illumination(&mut cancelled)?;
        Ok(result)
    }

    /// Ensemble energy flux, avoiding both coherency-magnitude and heat gates
    /// when a consumer only needs the signed vector power through a boundary.
    pub fn poynting_w_m2(
        &self,
        region: SlabRegion,
        position_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3]> {
        crystal::poll_illumination(&mut cancelled)?;
        admit_point(region, position_m, self.source.thickness_m())?;
        let mut sum = [Dd::ZERO; 3];
        for group in self.groups {
            crystal::poll_illumination(&mut cancelled)?;
            for (sum, value) in
                sum.iter_mut()
                    .zip(group.poynting_w_m2(region, position_m, &mut cancelled)?)
            {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let result = sum.map(Dd::to_f64);
        if result.iter().any(|value| !value.is_finite()) {
            return Err(ReflectionError::Invalid(
                "film illumination Poynting vector is unrepresentable",
            )
            .into());
        }
        crystal::poll_illumination(&mut cancelled)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
