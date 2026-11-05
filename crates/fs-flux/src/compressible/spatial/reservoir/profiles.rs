//! Face-resolved prescribed exterior gas on the original SSP time owner.
//!
//! These are complete Riemann traces of the original frozen gas. Sampling and
//! caloric admission occur once per explicit update by default. The opt-in
//! timed step instead samples its two physical stage times. In
//! particular, an averaged velocity is not relabeled as an average nonlinear
//! mass, momentum or energy flux. Time evolution and any turbulent field remain
//! with the caller; this owner never advances a random stream. A timed source
//! must be pure and replayable because rejected trials resample stage data.

use super::*;

mod timed;
pub(in crate::compressible::spatial) use timed::TimedSource;

/// Complete exterior primitive state on one physical face, in coherent SI.
/// Density and pressure jointly select temperature through the original gas
/// source; there is no independent temperature or caloric approximation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReservoirFaceState {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Cartesian velocity [m/s], not outward-normal components.
    pub velocity: [f64; 3],
    /// Positive absolute pressure [Pa].
    pub pressure_pa: f64,
}

/// Fixed side ranges and two reusable arrays of source-admitted face states.
#[derive(Default)]
pub(in crate::compressible::spatial) struct Profiles {
    ranges: [Option<(usize, usize)>; 6],
    values: Vec<AdmittedReservoir>,
    trial: Vec<AdmittedReservoir>,
    maximum_bytes: usize,
}

fn side_slot(side: BoundarySide) -> usize {
    2 * side.axis() + usize::from(side.sign() > 0.0)
}

fn layout(
    grid: &CartesianGrid,
    sides: &[BoundarySide],
) -> Result<([Option<(usize, usize)>; 6], usize), Error> {
    if sides.is_empty() || sides.len() > 6 {
        return Err(Error::Invalid(
            "one to six reservoir profile sides required",
        ));
    }
    let mut selected = [false; 6];
    for &side in sides {
        let slot = side_slot(side);
        if selected[slot] {
            return Err(Error::Invalid("duplicate reservoir profile side"));
        }
        selected[slot] = true;
    }
    let mut ranges = [None; 6];
    let mut count = 0_usize;
    for (slot, side) in BoundarySide::ALL.into_iter().enumerate() {
        if selected[slot] {
            let length = grid.face_count(side);
            ranges[slot] = Some((count, length));
            count = count
                .checked_add(length)
                .ok_or(Error::Invalid("reservoir profile face count overflow"))?;
        }
    }
    Ok((ranges, count))
}

fn required_bytes(count: usize) -> Result<usize, Error> {
    count
        .checked_mul(2 * core::mem::size_of::<AdmittedReservoir>())
        .ok_or(Error::Invalid("reservoir profile workspace size overflow"))
}

impl Profiles {
    pub(in crate::compressible::spatial) fn workspace_bytes(&self) -> usize {
        (self.values.capacity() + self.trial.capacity()) * core::mem::size_of::<AdmittedReservoir>()
    }

    fn admit_capacity(&self) -> Result<(), Error> {
        if self.workspace_bytes() > self.maximum_bytes {
            Err(Error::Invalid(
                "reservoir profile capacity exceeds its byte budget",
            ))
        } else {
            Ok(())
        }
    }

    pub(in crate::compressible::spatial) fn get(
        &self,
        slot: usize,
        face: usize,
    ) -> Option<&AdmittedReservoir> {
        let (start, length) = self.ranges[slot]?;
        (face < length).then(|| &self.values[start + face])
    }

    pub(in crate::compressible::spatial) fn get_staged(
        &self,
        slot: usize,
        face: usize,
    ) -> Option<&AdmittedReservoir> {
        let (start, length) = self.ranges[slot]?;
        (face < length).then(|| &self.trial[start + face])
    }

    pub(in crate::compressible::spatial) fn sides(&self) -> [bool; 6] {
        self.ranges.map(|range| range.is_some())
    }

    pub(in crate::compressible::spatial) fn try_clone(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        let mut copy = Self {
            ranges: self.ranges,
            values: allocated(self.values.len(), AdmittedReservoir::default())?,
            trial: allocated(self.values.len(), AdmittedReservoir::default())?,
            maximum_bytes: self.maximum_bytes,
        };
        copy.admit_capacity()?;
        for ((values, trial), source) in copy
            .values
            .chunks_mut(64)
            .zip(copy.trial.chunks_mut(64))
            .zip(self.values.chunks(64))
        {
            poll(cancelled)?;
            values.copy_from_slice(source);
            trial.copy_from_slice(source);
        }
        poll(cancelled)?;
        Ok(copy)
    }

    /// Fill only uncommitted storage. The original EOS supplies conserved
    /// energy and its decoded pressure/sound speed exactly as for native cells.
    pub(in crate::compressible::spatial) fn stage(
        &mut self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        mut sample: impl FnMut(BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        if self.values.is_empty() {
            return Err(Error::Invalid(
                "reservoir profiles must be prepared before updating",
            ));
        }
        for (slot, side) in BoundarySide::ALL.into_iter().enumerate() {
            let Some((start, length)) = self.ranges[slot] else {
                continue;
            };
            for face in 0..length {
                poll(cancelled)?;
                let point = grid.face_center(side, face).ok_or(Error::Invalid(
                    "reservoir profile no longer matches its grid",
                ))?;
                let input = sample(side, face, point)?;
                let state =
                    State::from_primitive(gas, input.density, input.velocity, input.pressure_pa)?;
                let primitive = state.primitive(gas, &mut *cancelled)?;
                self.trial[start + face] = AdmittedReservoir {
                    input,
                    state,
                    primitive,
                };
            }
        }
        poll(cancelled)
    }

    /// Synchronize a coupled owner's staged controls before either owner can
    /// publish them. Partial copies after cancellation remain private scratch.
    pub(in crate::compressible::spatial) fn copy_staged_from(
        &mut self,
        source: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        if self.ranges != source.ranges || self.trial.len() != source.trial.len() {
            return Err(Error::Invalid("coupled reservoir profile layouts differ"));
        }
        for (target, values) in self.trial.chunks_mut(64).zip(source.trial.chunks(64)) {
            poll(cancelled)?;
            target.copy_from_slice(values);
        }
        poll(cancelled)
    }

    /// Restore a composing owner's private bath after an earlier timed trial.
    /// Only scratch changes until the caller publishes this complete copy.
    pub(in crate::compressible::spatial) fn copy_accepted_from(
        &mut self,
        source: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        if self.ranges != source.ranges || self.trial.len() != source.values.len() {
            return Err(Error::Invalid("coupled reservoir profile layouts differ"));
        }
        for (target, values) in self.trial.chunks_mut(64).zip(source.values.chunks(64)) {
            poll(cancelled)?;
            target.copy_from_slice(values);
        }
        poll(cancelled)
    }

    pub(in crate::compressible::spatial) fn publish(&mut self) {
        std::mem::swap(&mut self.values, &mut self.trial);
    }

    pub(in crate::compressible::spatial) fn maximum_rate(
        &self,
        grid: &CartesianGrid,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        Self::rate(&self.values, grid, cancelled)
    }

    pub(in crate::compressible::spatial) fn maximum_staged_rate(
        &self,
        grid: &CartesianGrid,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        Self::rate(&self.trial, grid, cancelled)
    }

    fn rate(
        values: &[AdmittedReservoir],
        grid: &CartesianGrid,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let mut rate = 0.0_f64;
        for chunk in values.chunks(64) {
            poll(cancelled)?;
            for state in chunk {
                let p = state.primitive;
                let local: f64 = (0..3)
                    .map(|axis| (p.velocity[axis].abs() + p.sound_speed) / grid.spacing[axis])
                    .sum();
                if !local.is_finite() {
                    return Err(Error::Invalid("reservoir profile acoustic rate overflow"));
                }
                rate = rate.max(local);
            }
        }
        Ok(rate)
    }
}

impl Reservoirs {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::spatial) fn prepare_profiles(
        &self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        topology: Boundaries3d,
        sides: &[BoundarySide],
        maximum_profile_bytes: usize,
        sample: impl FnMut(BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Profiles, Error> {
        poll(cancelled)?;
        if !self.profiles.values.is_empty() {
            return Err(Error::Invalid(
                "reservoir profile geometry is already prepared",
            ));
        }
        let (ranges, count) = layout(grid, sides)?;
        for (slot, side) in BoundarySide::ALL.into_iter().enumerate() {
            if ranges[slot].is_none() {
                continue;
            }
            let condition = match topology.axis(side.axis()) {
                AxisBoundaries::Periodic => None,
                AxisBoundaries::Ends { lower, upper } => {
                    Some(if side.sign() < 0.0 { lower } else { upper })
                }
            };
            if condition != Some(BoundaryCondition::Transmissive)
                || self.states[slot].is_some()
                || self.outlets[slot].is_some()
            {
                return Err(Error::Invalid(
                    "reservoir profile requires an unbound transmissive side",
                ));
            }
        }
        if required_bytes(count)? > maximum_profile_bytes {
            return Err(Error::Invalid("reservoir profile exceeds its byte budget"));
        }
        let mut result = Profiles {
            ranges,
            values: allocated(count, AdmittedReservoir::default())?,
            trial: allocated(count, AdmittedReservoir::default())?,
            maximum_bytes: maximum_profile_bytes,
        };
        result.admit_capacity()?;
        result.stage(gas, grid, sample, cancelled)?;
        result.publish();
        Ok(result)
    }
}

impl Euler3d {
    /// Minimum additional numerical bytes for active and trial exterior states.
    /// Fixed side metadata, allocator overhead and the existing field are
    /// excluded. Molecular composition owns two full Euler fields and thus two
    /// copies of this workspace; its corresponding method admits both copies.
    pub fn required_reservoir_profile_workspace_bytes(
        &self,
        sides: &[BoundarySide],
    ) -> Result<usize, Error> {
        required_bytes(layout(&self.grid, sides)?.1)
    }

    /// Prepare a fixed set of spatially varying, complete-state exterior sides.
    ///
    /// The callback is evaluated once per face in [`BoundarySide::ALL`] order
    /// and then [`CartesianGrid::face_center`] order, independent of `sides`
    /// ordering. Positions are metres in the original grid frame. Each result
    /// is admitted by this field's full source EOS before any side is bound.
    /// Both arrays and their actual capacities must fit `maximum_profile_bytes`.
    /// Duplicate, periodic, wall, uniform-reservoir and pressure-outlet sides
    /// refuse. All requested sides bind together or none bind.
    ///
    /// Profiles remain fixed throughout `step` and its retries. The explicit
    /// [`Self::step_with_reservoir_profiles`] instead samples both stage times.
    /// They describe a prescribed infinite exterior, not finite storage,
    /// nonreflecting boundaries, molecular wall data or temporal interpolation.
    pub fn prepare_reservoir_profiles(
        &mut self,
        sides: &[BoundarySide],
        maximum_profile_bytes: usize,
        sample: impl FnMut(BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let profiles = self.reservoirs.prepare_profiles(
            &self.gas,
            &self.grid,
            self.boundaries,
            sides,
            maximum_profile_bytes,
            sample,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        self.reservoirs.profiles = profiles;
        self.last_step = None;
        Ok(())
    }

    /// Atomically replace ALL prepared face states, without allocating arrays.
    /// Callback/source/cancellation refusal preserves active traces, cells,
    /// clock, report and last accepted boundary impulses. On success the prior
    /// report is invalidated; cells, time and historical impulses are unchanged.
    /// The caller owns the sample time; advancing the external field is not
    /// part of this transaction, and no callback runs during `step`.
    pub fn update_reservoir_profiles(
        &mut self,
        sample: impl FnMut(BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.reservoirs
            .profiles
            .stage(&self.gas, &self.grid, sample, &mut cancelled)?;
        poll(&mut cancelled)?;
        self.reservoirs.profiles.publish();
        self.last_step = None;
        Ok(())
    }

    /// Prepared side mask in [`BoundarySide::ALL`] order.
    pub fn prepared_reservoir_profile_sides(&self) -> [bool; 6] {
        self.reservoirs.profiles.sides()
    }

    /// Whether the given side owns a prepared face profile.
    pub fn has_reservoir_profile(&self, side: BoundarySide) -> bool {
        self.prepared_reservoir_profile_sides()[side_slot(side)]
    }

    /// Original authored face trace. Uniform sides are also available at every
    /// valid physical face; unbound sides and out-of-range indices return None.
    pub fn reservoir_face(&self, side: BoundarySide, face: usize) -> Option<ReservoirFaceState> {
        if face >= self.grid.face_count(side) {
            return None;
        }
        self.reservoirs
            .get(side.axis(), side.sign())
            .or_else(|| self.reservoirs.profiles.get(side_slot(side), face))
            .map(|state| state.input)
    }

    /// Actual additional numerical capacity of both retained profile arrays.
    pub fn reservoir_profile_workspace_bytes(&self) -> usize {
        self.reservoirs.profiles.workspace_bytes()
    }
}

#[cfg(test)]
mod tests;
