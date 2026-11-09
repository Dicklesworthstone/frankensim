//! Held spatial reservoir temperatures on the original native heat carrier.

use super::{
    ThermalOpenError3, ThermalOpenFace3, ThermalOpenField3, check_boundary_temperature, invalid,
};

impl ThermalOpenField3 {
    /// Install or update a spatial reservoir temperature field.
    ///
    /// Face order is XMin/XMax/YMin/YMax/ZMin/ZMax. Entries follow x-fast plane
    /// order with the normal coordinate omitted, matching native momentum
    /// profiles: X faces use y+ny*z, Y faces x+nx*z, Z faces x+nx*y.
    /// The slice covers the whole original face, including inactive solid
    /// slots. Every value must be finite and inside the original scalar range,
    /// even if its link is currently outflow or inactive. Theta may be signed;
    /// material adapters own its conversion to physical Kelvin.
    ///
    /// Each inflowing link uses its own held temperature in the existing
    /// Danckwerts total-flux law. Outflow retains zero normal conductive flux,
    /// independent of this prescribed value; a later reversal uses this same
    /// link's profile entry. No callback or allocation occurs during a step.
    ///
    /// The profile shares the existing maximum_storage_bytes cap with the
    /// original thermal arrays and all other reservoir profiles. First install
    /// allocates one face vector; later updates reuse it. Validation/allocation
    /// failure leaves the old complete input and accepted physical state intact.
    ///
    /// The uniform boundary declaration remains the fallback. Clear the profile
    /// before updating that fallback. Face kinds are fixed for this carrier:
    /// non-reservoir targets refuse, and wall setters cannot change this face.
    /// Updating inputs neither heats the fluid nor resets its accepted history;
    /// a refused or cancelled step retains the already authored profile.
    pub fn set_reservoir_temperature_profile(
        &mut self,
        face: usize,
        temperatures: &[f64],
    ) -> Result<(), ThermalOpenError3> {
        if !matches!(
            self.boundaries.faces.get(face),
            Some(ThermalOpenFace3::Reservoir { .. })
        ) {
            return Err(invalid("thermal profile target is not a reservoir face"));
        }
        let entries = self.populations.len() / self.dims[face / 2];
        if temperatures.len() != entries {
            return Err(invalid("reservoir temperature profile shape mismatch"));
        }
        for &temperature in temperatures {
            check_boundary_temperature(face, temperature, self.limits)?;
        }
        // A retained profile's capacity remains owned even when replacing its
        // values. Other face capacities and original arrays share this same cap.
        let other_bytes = self
            .reservoir_profiles
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != face)
            .try_fold(self.base_storage_bytes, |bytes, (_, profile)| {
                profile.as_ref().map_or(Some(bytes), |values| {
                    values
                        .capacity()
                        .checked_mul(core::mem::size_of::<f64>())
                        .and_then(|size| bytes.checked_add(size))
                })
            })
            .ok_or_else(|| invalid("reservoir profile storage overflow"))?;
        let capacity = self.reservoir_profiles[face]
            .as_ref()
            .map_or(entries, Vec::capacity);
        let budget = self.limits.maximum_storage_bytes;
        let admit_capacity = |capacity: usize| {
            capacity
                .checked_mul(core::mem::size_of::<f64>())
                .and_then(|size| other_bytes.checked_add(size))
                .filter(|&bytes| bytes <= budget)
                .ok_or_else(|| invalid("thermal reservoir profile storage budget exceeded"))
        };
        admit_capacity(capacity)?;
        if let Some(profile) = &mut self.reservoir_profiles[face] {
            profile.copy_from_slice(temperatures);
        } else {
            let mut profile = Vec::new();
            profile
                .try_reserve_exact(entries)
                .map_err(|_| ThermalOpenError3::Allocation)?;
            admit_capacity(profile.capacity())?;
            profile.extend_from_slice(temperatures);
            self.reservoir_profiles[face] = Some(profile);
        }
        Ok(())
    }

    /// Remove only the held profile and restore the unchanged uniform fallback.
    /// Invalid/non-reservoir faces refuse; accepted heat, history and time stay.
    pub fn clear_reservoir_temperature_profile(
        &mut self,
        face: usize,
    ) -> Result<(), ThermalOpenError3> {
        if !matches!(
            self.boundaries.faces.get(face),
            Some(ThermalOpenFace3::Reservoir { .. })
        ) {
            return Err(invalid("thermal profile target is not a reservoir face"));
        }
        self.reservoir_profiles[face] = None;
        Ok(())
    }

    /// Held x-fast plane values; None means no profile (or an invalid face).
    /// Readback never evaluates a source or advances the thermal state.
    pub fn reservoir_temperature_profile(&self, face: usize) -> Option<&[f64]> {
        self.reservoir_profiles
            .get(face)
            .and_then(|profile| profile.as_deref())
    }
}
