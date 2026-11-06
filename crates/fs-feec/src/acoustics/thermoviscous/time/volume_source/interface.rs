//! Distributed mechanical and thermal drives in two original resolved regions.
//!
//! These bindings reuse the existing common-flow trajectory and its original
//! pressure, viscous traction, thermal and molecular coordinates. A passive
//! region may be joined to a volume-driven region without dummy source profiles.
//! Every original external column is retained in pressure, then per-region
//! structural, body-force and heat order. No partitioned source kick or extra
//! propagation is applied at the interface.

use super::*;

/// Borrow one original fixed-clock region with its complete source layout.
///
/// The volume cache is never exposed through the legacy boundary-only API.
/// Matching geometry, equilibrium, traces and clock are admitted by the same
/// normal or conforming interface owner used for passive regions.
#[derive(Clone, Copy)]
pub enum ThermoviscousTimeRegion<'p, 'a, 'c> {
    Boundary(&'p PreparedThermoviscousTime<'a, 'c>),
    Volume(&'p PreparedThermoviscousVolumeTime<'a, 'c>),
}

impl<'p, 'a, 'c> From<&'p PreparedThermoviscousTime<'a, 'c>>
    for ThermoviscousTimeRegion<'p, 'a, 'c>
{
    fn from(region: &'p PreparedThermoviscousTime<'a, 'c>) -> Self {
        Self::Boundary(region)
    }
}
impl<'p, 'a, 'c> From<&'p PreparedThermoviscousVolumeTime<'a, 'c>>
    for ThermoviscousTimeRegion<'p, 'a, 'c>
{
    fn from(region: &'p PreparedThermoviscousVolumeTime<'a, 'c>) -> Self {
        Self::Volume(region)
    }
}
impl<'p, 'a, 'c> ThermoviscousTimeRegion<'p, 'a, 'c> {
    fn fields(self) -> &'p PreparedThermoviscousTime<'a, 'c> {
        match self {
            Self::Boundary(region) => region,
            Self::Volume(region) => &region.fields,
        }
    }
    fn counts(self) -> [usize; 2] {
        match self {
            Self::Boundary(_) => [0, 0],
            Self::Volume(region) => [region.force_count, region.heat_count],
        }
    }
}

/// Normal-traction, insulated interface with original distributed sources.
/// The full original field on each side shares one exact held-input trajectory.
pub struct PreparedThermoviscousVolumeTimePair<'p, 'a, 'c> {
    inner: PreparedThermoviscousTimePair<'p, 'a, 'c>,
    counts: [[usize; 2]; 2],
}

impl<'p, 'a, 'c> PreparedThermoviscousVolumeTimePair<'p, 'a, 'c> {
    /// Connect selected original pressure modes, preserving their full trace
    /// basis map. Either side may be passive or volume driven. Sources are
    /// not sampled again, projected to the boundary, or merged between regions.
    pub fn connect(
        left: impl Into<ThermoviscousTimeRegion<'p, 'a, 'c>>,
        right: impl Into<ThermoviscousTimeRegion<'p, 'a, 'c>>,
        left_modes: &[usize],
        right_modes: &[usize],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousTimeError> {
        let regions = [left.into(), right.into()];
        let inner = regions[0].fields().connect(
            regions[1].fields(),
            left_modes,
            right_modes,
            limits,
            cancelled,
        )?;
        Ok(Self {
            inner,
            counts: regions.map(|region| region.counts()),
        })
    }

    pub fn system(&self) -> &PortHamiltonian {
        self.inner.system()
    }
    pub fn duration_s(&self) -> f64 {
        self.inner.duration_s()
    }
    pub fn traction_map(&self) -> &[f64] {
        self.inner.traction_map()
    }
    pub fn interface_modes(&self) -> [&[usize]; 2] {
        self.inner.interface_modes()
    }
    pub fn external_pressure_modes(&self) -> &[(usize, usize)] {
        self.inner.external_pressure_modes()
    }
    pub fn zero_states(&self) -> [ThermoviscousTimeState; 2] {
        self.inner.zero_states()
    }
    /// Force-profile counts in left, right order.
    pub fn body_force_source_counts(&self) -> [usize; 2] {
        self.counts.map(|count| count[0])
    }
    /// Heat-profile counts in left, right order.
    pub fn heat_source_counts(&self) -> [usize; 2] {
        self.counts.map(|count| count[1])
    }

    /// Evolve both complete fields with held source amplitudes [N] and [W].
    /// All inputs use original region/source ordering. Per-region work includes
    /// signed body-force and thermal available-energy work on this SAME path;
    /// internal mechanical work cancels in the enclosing pair. No factor 1/2
    /// enters causal work. Both input states remain untouched on every refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &self,
        initial: [&ThermoviscousTimeState; 2],
        external_traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: [&[f64]; 2],
        force_drive_n: [&[f64]; 2],
        heat_drive_w: [&[f64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousTimePairStep, ThermoviscousTimeError> {
        poll(cancelled)?;
        admit_sources(self.counts, force_drive_n, heat_drive_w)?;
        self.inner.step_with_volume_sources(
            initial,
            external_traction_pa,
            wall_forces_n_per_sqrt_kg,
            force_drive_n,
            heat_drive_w,
            cancelled,
        )
    }
}

/// Complete matching P2 velocity/temperature interface with distributed drives.
/// Mechanical traction and conductive flux are solved simultaneously with all
/// original fields, including molecular temperatures and conducting solids.
pub struct PreparedConformingThermoviscousVolumeTimePair<'p, 'a, 'c> {
    inner: PreparedConformingThermoviscousTimePair<'p, 'a, 'c>,
    counts: [[usize; 2]; 2],
}

impl<'p, 'a, 'c> PreparedConformingThermoviscousVolumeTimePair<'p, 'a, 'c> {
    /// Connect matching original face triangulations with the native complete
    /// P2 trace constraint. There is no mean-flow or species mixing. Active
    /// temperature and velocity constraints retain their original ownership.
    pub fn connect(
        left: impl Into<ThermoviscousTimeRegion<'p, 'a, 'c>>,
        right: impl Into<ThermoviscousTimeRegion<'p, 'a, 'c>>,
        left_faces: &[[u32; 3]],
        right_faces: &[[u32; 3]],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousTimeError> {
        let regions = [left.into(), right.into()];
        let inner = regions[0].fields().connect_conforming(
            regions[1].fields(),
            left_faces,
            right_faces,
            limits,
            cancelled,
        )?;
        Ok(Self {
            inner,
            counts: regions.map(|region| region.counts()),
        })
    }

    pub fn interface(&self) -> &ConformingThermoviscousInterface {
        self.inner.interface()
    }
    pub fn system(&self) -> &PortHamiltonian {
        self.inner.system()
    }
    pub fn duration_s(&self) -> f64 {
        self.inner.duration_s()
    }
    pub fn external_pressure_modes(&self) -> &[(usize, usize)] {
        self.inner.external_pressure_modes()
    }
    pub fn zero_states(&self) -> [ThermoviscousTimeState; 2] {
        self.inner.zero_states()
    }
    pub fn body_force_source_counts(&self) -> [usize; 2] {
        self.counts.map(|count| count[0])
    }
    pub fn heat_source_counts(&self) -> [usize; 2] {
        self.counts.map(|count| count[1])
    }

    /// The same joint constrained trajectory receives every original drive.
    /// Body-force work uses the integral of its spatial velocity dual; heat
    /// work uses the integral of its translational T'/T0 dual. These are kept
    /// separate from the two signed internal interface transfers and the
    /// identified first-order caloric flux. The entire original state is kept.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &self,
        initial: [&ThermoviscousTimeState; 2],
        external_traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: [&[f64]; 2],
        force_drive_n: [&[f64]; 2],
        heat_drive_w: [&[f64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConformingThermoviscousTimePairStep, ThermoviscousTimeError> {
        poll(cancelled)?;
        admit_sources(self.counts, force_drive_n, heat_drive_w)?;
        self.inner.step_with_volume_sources(
            initial,
            external_traction_pa,
            wall_forces_n_per_sqrt_kg,
            force_drive_n,
            heat_drive_w,
            cancelled,
        )
    }
}

fn admit_sources(
    counts: [[usize; 2]; 2],
    force_drive_n: [&[f64]; 2],
    heat_drive_w: [&[f64]; 2],
) -> Result<(), ThermoviscousTimeError> {
    if (0..2).any(|side| {
        force_drive_n[side].len() != counts[side][0] || heat_drive_w[side].len() != counts[side][1]
    }) {
        return Err(invalid(
            "complete original per-region force and heat drives required",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
