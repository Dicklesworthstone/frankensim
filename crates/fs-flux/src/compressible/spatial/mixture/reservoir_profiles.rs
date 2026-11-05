//! Spatially varying complete exterior traces on the original mixture faces.
//!
//! A profile supplies one physical primitive trace per face, in the same order
//! as `CartesianGrid::face_center`. The original `MixtureGas` admits every
//! composition, pressure and velocity; no caloric table or species gauge is
//! duplicated here. The native HLLC flux determines inflow, outflow and reversal.
//! Profiles are held through all stages and retries until explicitly replaced.
//! They are not finite baths, characteristic outlets, or nonreflecting data.

use super::*;
use crate::compressible::AxisBoundaries;

/// Velocity and pressure returned after filling a face's canonical species row.
///
/// The callback also receives `&mut [f64]` with exactly `gas.species_count()`
/// entries and must fill every partial mass density [kg/m³], including absent
/// species. No normalization, clipping or remainder species is introduced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirFaceState {
    /// Cartesian velocity [m/s], not outward-normal components.
    pub velocity: [f64; 3],
    /// Absolute pressure [Pa].
    pub pressure_pa: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Span {
    start: usize,
    count: usize,
}

/// Two fixed numerical payloads make replacement a single publication point.
pub(super) struct Profiles {
    slots: [Option<Span>; 6],
    rows: Vec<f64>,
    primitive: Vec<Primitive>,
    trial_rows: Vec<f64>,
    trial_primitive: Vec<Primitive>,
    rate: f64,
    trial_rate: f64,
}

impl MixtureEuler3d {
    /// Additional bytes for active and trial face rows and primitive caches.
    /// The original field's numerical storage is separate and remains counted.
    pub fn required_reservoir_profile_workspace_bytes(
        &self,
        sides: &[BoundarySide],
    ) -> Result<usize, Error> {
        let (_, count) = selected_faces(self, sides)?;
        Profiles::required_bytes(self.gas(), count)
    }

    /// Prepare complete prescribed states on selected transmissive sides once.
    ///
    /// The callback sees canonical side order, then increasing physical face
    /// index and its centre [m] relative to the grid origin. Its species slice
    /// starts as NaNs so an omitted species cannot reuse a previous face's gas.
    /// Every supplied state must pass the original source admission. Uniform
    /// reservoirs, walls, periodic sides and pressure-only outlets cannot also
    /// own a profile. The entire preparation is transactional.
    ///
    /// `maximum_profile_bytes` bounds only additional active/trial numerical
    /// arrays. Their actual capacities also enter this field's stored byte
    /// count and the common total numerical ceiling. No per-face allocation is
    /// performed during subsequent replacement or native timestepping.
    pub fn prepare_reservoir_profiles(
        &mut self,
        sides: &[BoundarySide],
        maximum_profile_bytes: usize,
        mut sample: impl FnMut(
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        let (slots, count) = selected_faces(self, sides)?;
        let required = Profiles::required_bytes(&self.gas, count)?;
        admit_additional(self.workspace_bytes, required, maximum_profile_bytes)?;
        let width = self.gas.species_count() + 4;
        let mut profiles = Profiles {
            slots,
            rows: allocated(count * width, 0.0)?,
            primitive: allocated(count, Primitive::default())?,
            trial_rows: allocated(count * width, 0.0)?,
            trial_primitive: allocated(count, Primitive::default())?,
            rate: 0.0,
            trial_rate: 0.0,
        };
        let bytes = profiles.workspace_bytes();
        let total = admit_additional(self.workspace_bytes, bytes, maximum_profile_bytes)?;
        profiles.stage(&self.gas, &self.grid, &mut sample, &mut cancelled)?;
        stepping::poll::<Error>(&mut cancelled)?;
        profiles.commit();
        self.reservoirs.profiles = Some(profiles);
        self.workspace_bytes = total;
        Ok(())
    }

    /// Replace every prepared face using existing scratch, or publish nothing.
    ///
    /// Profiles remain held across the complete next native step, including
    /// shortening retries. This method neither advances a stochastic source
    /// nor claims stage-time interpolation. A provider should use one coherent
    /// immutable realization; provider side effects are outside the rollback.
    /// A refused replacement preserves cells, clock, loads and all active faces.
    pub fn update_reservoir_profiles(
        &mut self,
        mut sample: impl FnMut(
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        let profiles = self.reservoirs.profiles.as_mut().ok_or(Error::Invalid(
            "mixture reservoir profiles are not prepared",
        ))?;
        profiles.stage(&self.gas, &self.grid, &mut sample, &mut cancelled)?;
        stepping::poll::<Error>(&mut cancelled)?;
        profiles.commit();
        Ok(())
    }

    /// Exterior state at one physical face, including uniform reservoirs.
    /// Returns `None` for unbound sides or an out-of-range physical face.
    pub fn reservoir_face(&self, side: BoundarySide, face: usize) -> Option<ReservoirBoundary<'_>> {
        if face >= self.grid.face_count(side) {
            return None;
        }
        self.reservoirs
            .boundary_face(side, face, self.gas.species_count())
    }

    /// Whether a nonuniform-capable profile owns this side.
    pub fn has_reservoir_profile(&self, side: BoundarySide) -> bool {
        self.prepared_reservoir_profile_sides()[slot(side)]
    }

    /// Profile ownership in `BoundarySide::ALL` order.
    pub fn prepared_reservoir_profile_sides(&self) -> [bool; 6] {
        self.reservoirs
            .profiles
            .as_ref()
            .map_or([false; 6], |p| p.slots.map(|span| span.is_some()))
    }

    /// Actual additional active/trial numerical capacity [bytes].
    pub fn reservoir_profile_workspace_bytes(&self) -> usize {
        self.reservoirs
            .profiles
            .as_ref()
            .map_or(0, Profiles::workspace_bytes)
    }
}

fn selected_faces(
    flow: &MixtureEuler3d,
    sides: &[BoundarySide],
) -> Result<([Option<Span>; 6], usize), Error> {
    if sides.is_empty() || sides.len() > 6 || flow.reservoirs.profiles.is_some() {
        return Err(Error::Invalid(
            "prepare one nonempty set of mixture reservoir profiles",
        ));
    }
    let mut selected = [false; 6];
    for &side in sides {
        let condition = match flow.boundaries.axis(side.axis()) {
            AxisBoundaries::Periodic => None,
            AxisBoundaries::Ends { lower, upper } => {
                Some(if side.sign() < 0.0 { lower } else { upper })
            }
        };
        if selected[slot(side)]
            || condition != Some(BoundaryCondition::Transmissive)
            || flow.reservoir(side).is_some()
        {
            return Err(Error::Invalid(
                "mixture profile needs a unique unbound open side",
            ));
        }
        selected[slot(side)] = true;
    }
    let mut slots = [None; 6];
    let mut count = 0usize;
    for side in BoundarySide::ALL {
        if selected[slot(side)] {
            let faces = flow.grid.face_count(side);
            slots[slot(side)] = Some(Span {
                start: count,
                count: faces,
            });
            count = count.checked_add(faces).ok_or(Error::Invalid(
                "mixture reservoir profile face count overflow",
            ))?;
        }
    }
    Ok((slots, count))
}

fn admit_additional(base: usize, extra: usize, maximum: usize) -> Result<usize, Error> {
    let total = base.checked_add(extra).ok_or(Error::Invalid(
        "mixture reservoir profile workspace overflow",
    ))?;
    if extra > maximum || total > MAX_MIXTURE_WORKSPACE_BYTES {
        return Err(Error::Invalid("mixture reservoir profile workspace budget"));
    }
    Ok(total)
}

fn slot(side: BoundarySide) -> usize {
    2 * side.axis() + usize::from(side.sign() > 0.0)
}

impl Profiles {
    fn required_bytes(gas: &MixtureGas, count: usize) -> Result<usize, Error> {
        let row = (gas.species_count() + 4)
            .checked_mul(size_of::<f64>())
            .and_then(|bytes| bytes.checked_add(size_of::<Primitive>()))
            .ok_or(Error::Invalid(
                "mixture reservoir profile row size overflow",
            ))?;
        count
            .checked_mul(row)
            .and_then(|bytes| bytes.checked_mul(2))
            .ok_or(Error::Invalid(
                "mixture reservoir profile workspace overflow",
            ))
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        (self.rows.capacity() + self.trial_rows.capacity()) * size_of::<f64>()
            + (self.primitive.capacity() + self.trial_primitive.capacity()) * size_of::<Primitive>()
    }

    fn index(&self, side: BoundarySide, face: usize) -> Option<usize> {
        self.slots[slot(side)]
            .filter(|span| face < span.count)
            .map(|span| span.start + face)
    }

    pub(super) fn row(&self, side: BoundarySide, face: usize, species: usize) -> Option<&[f64]> {
        let width = species + 4;
        self.index(side, face)
            .map(|index| &self.rows[index * width..(index + 1) * width])
    }

    pub(super) fn primitive(&self, side: BoundarySide, face: usize) -> Option<Primitive> {
        self.index(side, face).map(|index| self.primitive[index])
    }

    pub(super) fn maximum_rate(&self) -> f64 {
        self.rate
    }

    fn stage(
        &mut self,
        gas: &MixtureGas,
        grid: &CartesianGrid,
        sample: &mut impl FnMut(
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<ReservoirFaceState, Error>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let n = gas.species_count();
        let width = n + 4;
        let mut conservative = [0.0; WIDTH_LIMIT];
        let mut rate = 0.0_f64;
        for side in BoundarySide::ALL {
            let Some(span) = self.slots[slot(side)] else {
                continue;
            };
            for face in 0..span.count {
                stepping::poll::<Error>(cancelled)?;
                let index = span.start + face;
                let row = &mut self.trial_rows[index * width..(index + 1) * width];
                row[..n].fill(f64::NAN);
                let input = sample(
                    side,
                    face,
                    grid.face_center(side, face)
                        .expect("prepared physical face"),
                    &mut row[..n],
                )?;
                row[n..n + 3].copy_from_slice(&input.velocity);
                row[n + 3] = input.pressure_pa;
                let p = gas.spatial_conserved_into(
                    normal_primitive(row, n, 0),
                    &mut conservative[..width],
                )?;
                self.trial_primitive[index] = Primitive {
                    density: p.density,
                    velocity: input.velocity,
                    pressure: p.pressure,
                    temperature: p.temperature,
                    sound_speed: p.sound_speed,
                };
                let local = (0..3)
                    .map(|axis| (input.velocity[axis].abs() + p.sound_speed) / grid.spacing[axis])
                    .sum::<f64>();
                if !local.is_finite() {
                    return Err(Error::Invalid(
                        "mixture reservoir profile wave rate overflow",
                    ));
                }
                rate = rate.max(local);
            }
        }
        stepping::poll::<Error>(cancelled)?;
        self.trial_rate = rate;
        Ok(())
    }

    fn commit(&mut self) {
        std::mem::swap(&mut self.rows, &mut self.trial_rows);
        std::mem::swap(&mut self.primitive, &mut self.trial_primitive);
        std::mem::swap(&mut self.rate, &mut self.trial_rate);
    }

    pub(super) fn try_clone(&self, cancelled: &mut impl FnMut() -> bool) -> Result<Self, Error> {
        let mut result = Self {
            slots: self.slots,
            rows: allocated(self.rows.len(), 0.0)?,
            primitive: allocated(self.primitive.len(), Primitive::default())?,
            trial_rows: allocated(self.trial_rows.len(), 0.0)?,
            trial_primitive: allocated(self.trial_primitive.len(), Primitive::default())?,
            rate: self.rate,
            trial_rate: self.trial_rate,
        };
        copy_values(&self.rows, &mut result.rows, cancelled)?;
        copy_values(&self.trial_rows, &mut result.trial_rows, cancelled)?;
        copy_primitives(&self.primitive, &mut result.primitive, cancelled)?;
        copy_primitives(
            &self.trial_primitive,
            &mut result.trial_primitive,
            cancelled,
        )?;
        stepping::poll::<Error>(cancelled)?;
        Ok(result)
    }

    fn copy_candidate_from(
        &mut self,
        source: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        if self.slots != source.slots
            || self.trial_rows.len() != source.trial_rows.len()
            || self.trial_primitive.len() != source.trial_primitive.len()
        {
            return Err(Error::Invalid(
                "coupled mixture reservoir profile layout differs",
            ));
        }
        copy_values(&source.trial_rows, &mut self.trial_rows, cancelled)?;
        copy_primitives(
            &source.trial_primitive,
            &mut self.trial_primitive,
            cancelled,
        )?;
        self.trial_rate = source.trial_rate;
        Ok(())
    }
}

fn copy_primitives(
    source: &[Primitive],
    target: &mut [Primitive],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    for (source, target) in source.chunks(1024).zip(target.chunks_mut(1024)) {
        stepping::poll::<Error>(cancelled)?;
        target.copy_from_slice(source);
    }
    Ok(())
}

/// Shared control transaction for composed accepted/trial field owners.
/// Providers are evaluated once, then their admitted candidate is copied before
/// either owner publishes. Neither flow clock or conservation report advances.
pub(super) fn update_pair(
    first: &mut MixtureEuler3d,
    second: &mut MixtureEuler3d,
    mut sample: impl FnMut(
        BoundarySide,
        usize,
        [f64; 3],
        &mut [f64],
    ) -> Result<ReservoirFaceState, Error>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(), Error> {
    stepping::poll::<Error>(&mut cancelled)?;
    if first.gas != second.gas || first.grid != second.grid || first.boundaries != second.boundaries
    {
        return Err(Error::Invalid(
            "coupled mixture reservoir profile binding differs",
        ));
    }
    let a = first.reservoirs.profiles.as_mut().ok_or(Error::Invalid(
        "mixture reservoir profiles are not prepared",
    ))?;
    let b = second.reservoirs.profiles.as_mut().ok_or(Error::Invalid(
        "coupled mixture reservoir profiles are not prepared",
    ))?;
    a.stage(&first.gas, &first.grid, &mut sample, &mut cancelled)?;
    b.copy_candidate_from(a, &mut cancelled)?;
    stepping::poll::<Error>(&mut cancelled)?;
    a.commit();
    b.commit();
    Ok(())
}
