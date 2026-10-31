//! Source-normalized current constraints shared by harmonic and transient coils.
//!
//! One connected invariant-z conductor section carries `n I`. With
//! `S = integral(sigma dA)`, its uniform electric drive is `n I/S`, reciprocal
//! unit-current profile is `n sigma/S`, and DC resistance is `depth*n^2/S`.
//! Geometry/source admission and the actual projected field/heat solve remain
//! with the original planar owners; this module introduces no field solver.

use super::{EddyError, finite, poll, positive};

/// One complete connected conducting section carries `signed_turns * I`.
/// Every return section has an explicit sign. The complete winding profile
/// must remain current-balanced on every magnetic mesh component.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidWindingRegion {
    pub region: u32,
    pub winding: usize,
    pub signed_turns: f64,
}

/// Original sorted conductivity-region coordinates, shared by both time lanes.
pub(crate) struct SolidRegionLayout {
    pub bindings: Vec<SolidWindingRegion>,
    pub region_bindings: Vec<Option<usize>>,
    pub electric_per_current_v_m_a: Vec<f64>,
    pub resistance_ohm: Vec<f64>,
}
impl SolidRegionLayout {
    pub(crate) fn new(
        regions: &[u32],
        sigma_volumes: &[f64],
        windings: usize,
        depth: f64,
        mut bindings: Vec<SolidWindingRegion>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, EddyError> {
        if regions.len() != sigma_volumes.len()
            || bindings.len() > regions.len()
            || regions.windows(2).any(|pair| pair[0] >= pair[1])
            || !positive(depth)
        {
            return Err(EddyError::Invalid(
                "complete sorted solid-region geometry required",
            ));
        }
        bindings.sort_by_key(|binding| binding.region);
        let mut region_bindings = vec![None; regions.len()];
        let mut electric_per_current_v_m_a = vec![0.0; regions.len()];
        let mut resistance_ohm = vec![0.0; windings];
        for (index, binding) in bindings.iter().enumerate() {
            poll(cancelled)?;
            if binding.winding >= windings
                || !binding.signed_turns.is_finite()
                || binding.signed_turns == 0.0
                || index > 0 && bindings[index - 1].region == binding.region
            {
                return Err(EddyError::Invalid(
                    "invalid or duplicate solid winding region binding",
                ));
            }
            let group = regions.binary_search(&binding.region).map_err(|_| {
                EddyError::Invalid("solid winding region has no admitted conductor section")
            })?;
            if !positive(sigma_volumes[group]) {
                return Err(EddyError::Invalid(
                    "solid conductivity volume must be positive",
                ));
            }
            let drive = finite((binding.signed_turns * depth) / sigma_volumes[group])?;
            let resistance = finite((binding.signed_turns * depth) * drive)?;
            if drive == 0.0 || !positive(resistance) {
                return Err(EddyError::Invalid(
                    "solid winding drive or DC resistance is unrepresentable",
                ));
            }
            region_bindings[group] = Some(index);
            electric_per_current_v_m_a[group] = drive;
            resistance_ohm[binding.winding] = finite(resistance_ohm[binding.winding] + resistance)?;
        }
        Ok(Self {
            bindings,
            region_bindings,
            electric_per_current_v_m_a,
            resistance_ohm,
        })
    }

    pub(crate) fn drives(
        &self,
        currents: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>), EddyError> {
        if currents.len() != self.resistance_ohm.len() || currents.iter().any(|v| !v.is_finite()) {
            return Err(EddyError::Invalid(
                "complete finite solid winding currents required",
            ));
        }
        let mut drives = vec![0.0; self.region_bindings.len()];
        let mut targets = vec![0.0; drives.len()];
        for (group, index) in self.region_bindings.iter().enumerate() {
            poll(cancelled)?;
            if let Some(index) = index {
                let binding = self.bindings[*index];
                drives[group] =
                    finite(self.electric_per_current_v_m_a[group] * currents[binding.winding])?;
                targets[group] = finite(binding.signed_turns * currents[binding.winding])?;
            }
        }
        Ok((drives, targets))
    }
}
