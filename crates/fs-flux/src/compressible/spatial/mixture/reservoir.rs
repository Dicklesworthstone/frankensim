//! Source-admitted exterior composition held on physical open sides.

use super::*;
use crate::compressible::AxisBoundaries;
use crate::compressible::mixture::finite::FiniteReservoirs;

mod finite;
pub use finite::{FiniteReservoir, FiniteReservoirReport};

/// Uniform externally maintained gas on one original physical side.
///
/// Every density uses the field's canonical source order, including absent
/// species. The source alone supplies temperature, calorics and energy gauges.
/// The original HLLC solve determines inflow or outflow from the two traces.
/// These reservoirs have prescribed infinite storage, not a second time owner,
/// a nonreflecting outlet, or a molecular species/thermal boundary condition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirBoundary<'a> {
    pub side: BoundarySide,
    /// Original canonical species mass densities [kg/m³].
    pub partial_densities: &'a [f64],
    /// Cartesian velocity [m/s], not outward-normal components.
    pub velocity: [f64; 3],
    /// Absolute pressure [Pa].
    pub pressure_pa: f64,
}

/// One flat primitive row per side and its original-source caloric cache.
pub(super) struct Reservoirs {
    slots: [Option<usize>; 6],
    stagnation: [bool; 6],
    reservoir_energy: [f64; 6],
    states: Vec<f64>,
    primitive: Vec<Primitive>,
    pub(super) profiles: Option<reservoir_profiles::Profiles>,
    pub(super) finite: FiniteReservoirs<3, 6>,
}

impl Reservoirs {
    pub(super) fn required_bytes(gas: &MixtureGas, count: usize) -> Result<usize, Error> {
        if count > 6 {
            return Err(Error::Invalid(
                "at most six volumetric mixture reservoir sides",
            ));
        }
        count
            .checked_mul((gas.species_count() + 4) * size_of::<f64>() + size_of::<Primitive>())
            .ok_or(Error::Invalid(
                "volumetric reservoir workspace size overflow",
            ))
    }

    pub(super) fn new(
        gas: &MixtureGas,
        topology: Boundaries3d,
        inputs: &[ReservoirBoundary<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::required_bytes(gas, inputs.len())?;
        let n = gas.species_count();
        let width = n + 4;
        let mut result = Self {
            slots: [None; 6],
            stagnation: [false; 6],
            reservoir_energy: [0.0; 6],
            states: allocated(inputs.len() * width, 0.0)?,
            primitive: allocated(inputs.len(), Primitive::default())?,
            profiles: None,
            finite: FiniteReservoirs::empty(),
        };
        let mut conservative = [0.0; WIDTH_LIMIT];
        for (index, input) in inputs.iter().enumerate() {
            stepping::poll::<Error>(cancelled)?;
            let axis = input.side.axis();
            let sign = input.side.sign();
            let slot = side_index(axis, sign);
            let condition = match topology.axis(axis) {
                AxisBoundaries::Periodic => None,
                AxisBoundaries::Ends { lower, upper } => {
                    Some(if sign < 0.0 { lower } else { upper })
                }
            };
            if result.slots[slot].is_some()
                || condition != Some(BoundaryCondition::Transmissive)
                || input.partial_densities.len() != n
            {
                return Err(Error::Invalid(
                    "reservoir needs a unique open volumetric side and original species axis",
                ));
            }
            let row = &mut result.states[index * width..(index + 1) * width];
            row[..n].copy_from_slice(input.partial_densities);
            row[n..n + 3].copy_from_slice(&input.velocity);
            row[n + 3] = input.pressure_pa;
            let original = gas
                .spatial_conserved_into(normal_primitive(row, n, 0), &mut conservative[..width])?;
            result.primitive[index] = Primitive {
                density: original.density,
                velocity: input.velocity,
                pressure: original.pressure,
                temperature: original.temperature,
                sound_speed: original.sound_speed,
            };
            result.reservoir_energy[slot] = conservative[width - 1];
            result.slots[slot] = Some(index);
        }
        stepping::poll::<Error>(cancelled)?;
        Ok(result)
    }

    pub(super) fn set_stagnation(
        &mut self,
        gas: &MixtureGas,
        sides: &[BoundarySide],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        stepping::poll::<Error>(cancelled)?;
        if sides.len() > 6 {
            return Err(Error::Invalid("too many compact reservoir sides"));
        }
        let n = gas.species_count();
        let mut candidate = std::array::from_fn(|slot| self.finite.has_side(slot));
        let mut selected = [false; 6];
        let mut flux = [0.0; WIDTH_LIMIT];
        for &side in sides {
            stepping::poll::<Error>(cancelled)?;
            let axis = side.axis();
            let sign = side.sign();
            let slot = side_index(axis, sign);
            let index = self.slots[slot].ok_or(Error::Invalid(
                "compact reservoir needs an already bound uniform exterior",
            ))?;
            if selected[slot] || self.primitive[index].velocity != [0.0; 3] {
                return Err(Error::Invalid(
                    "compact reservoir sides must be unique and exactly stationary",
                ));
            }
            let row = self.row(axis, sign, n).expect("uniform reservoir slot");
            let reservoir = fs_gas::multicomponent::Cell {
                partial_densities: &row[..n],
                momentum: 0.0,
                energy: self.reservoir_energy[slot],
            };
            gas.reservoir_spatial_flux_from_primitive(
                normal_primitive(row, n, axis),
                reservoir,
                sign,
                &mut flux[..n + 4],
                &mut *cancelled,
            )?;
            candidate[slot] = true;
            selected[slot] = true;
        }
        stepping::poll::<Error>(cancelled)?;
        self.stagnation = candidate;
        Ok(())
    }

    pub(super) fn is_stagnation(&self, axis: usize, sign: f64) -> bool {
        self.stagnation[side_index(axis, sign)]
    }

    pub(super) fn has_stagnation(&self) -> bool {
        self.stagnation.iter().any(|selected| *selected)
    }

    pub(super) fn stagnation_cell(
        &self,
        axis: usize,
        sign: f64,
        species: usize,
    ) -> Option<fs_gas::multicomponent::Cell<'_>> {
        self.stagnation_cell_at(axis, sign, species, false)
    }

    pub(super) fn stagnation_cell_at(
        &self,
        axis: usize,
        sign: f64,
        species: usize,
        trial_finite: bool,
    ) -> Option<fs_gas::multicomponent::Cell<'_>> {
        let slot = side_index(axis, sign);
        if let Some(state) = self.finite.state(slot, trial_finite) {
            return Some(state);
        }
        if !self.stagnation[slot] {
            return None;
        }
        let row = self
            .row(axis, sign, species)
            .expect("admitted compact reservoir");
        Some(fs_gas::multicomponent::Cell {
            partial_densities: &row[..species],
            momentum: 0.0,
            energy: self.reservoir_energy[slot],
        })
    }

    pub(super) fn stagnation_rate(
        &self,
        gas: &MixtureGas,
        grid: &CartesianGrid,
        states: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let n = gas.species_count();
        let width = n + 4;
        let mut flux = [0.0; WIDTH_LIMIT];
        let mut maximum: f64 = 0.0;
        for side in BoundarySide::ALL {
            let axis = side.axis();
            let sign = side.sign();
            let Some(reservoir) = self.stagnation_cell(axis, sign, n) else {
                continue;
            };
            let [a, b] = tangents(axis);
            let mut compact_rate = 0.0;
            let mut compact_correction = 0.0;
            for face in 0..grid.face_count(side) {
                let mut c = [0; 3];
                c[axis] = if sign < 0.0 { 0 } else { grid.shape[axis] - 1 };
                c[a] = face % grid.shape[a];
                c[b] = face / grid.shape[a];
                let k = grid.index(c);
                stepping::poll::<Error>(cancelled)?;
                let q = cell(&states[k * width..(k + 1) * width], n);
                let face = gas.reservoir_spatial_flux(
                    q.normal(axis),
                    reservoir,
                    sign,
                    &mut flux[..width],
                    &mut *cancelled,
                )?;
                let density = sum(q.partial_densities.iter().copied());
                let local = (0..3)
                    .map(|d| {
                        (face.signal_speed
                            + if d == axis {
                                0.0
                            } else {
                                (q.momentum[d] / density).abs()
                            })
                            / grid.spacing[d]
                    })
                    .sum::<f64>();
                maximum = maximum.max(local);
                if let Some(volume) = self.finite.volume(side_index(axis, sign)) {
                    let increment = 0.5 * grid.face_area(side) * face.signal_speed / volume;
                    let y = increment - compact_correction;
                    let next = compact_rate + y;
                    compact_correction = (next - compact_rate) - y;
                    compact_rate = next;
                }
            }
            if !compact_rate.is_finite() {
                return Err(Error::Invalid("finite volumetric reservoir CFL overflow"));
            }
            maximum = maximum.max(compact_rate);
        }
        Ok(maximum)
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        self.states.capacity() * size_of::<f64>()
            + self.primitive.capacity() * size_of::<Primitive>()
            + self.profiles.as_ref().map_or(0, |p| p.workspace_bytes())
            + self.finite.workspace_bytes()
    }

    pub(super) fn try_clone(
        &self,
        gas: &MixtureGas,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let mut states = allocated(self.states.len(), 0.0)?;
        copy_values(&self.states, &mut states, cancelled)?;
        let mut primitive = allocated(self.primitive.len(), Primitive::default())?;
        primitive.copy_from_slice(&self.primitive);
        let profiles = self
            .profiles
            .as_ref()
            .map(|profile| profile.try_clone(cancelled))
            .transpose()?;
        stepping::poll::<Error>(cancelled)?;
        Ok(Self {
            slots: self.slots,
            stagnation: self.stagnation,
            reservoir_energy: self.reservoir_energy,
            states,
            primitive,
            profiles,
            finite: self.finite.try_clone(gas, cancelled)?,
        })
    }

    pub(super) fn boundary(
        &self,
        side: BoundarySide,
        species: usize,
    ) -> Option<ReservoirBoundary<'_>> {
        let slot = side_index(side.axis(), side.sign());
        if let Some(state) = self.finite.state(slot, false) {
            return Some(ReservoirBoundary {
                side,
                partial_densities: state.partial_densities,
                velocity: [0.0; 3],
                pressure_pa: self.finite.primitive(slot, false)?.pressure,
            });
        }
        self.row(side.axis(), side.sign(), species)
            .map(|row| ReservoirBoundary {
                side,
                partial_densities: &row[..species],
                velocity: [row[species], row[species + 1], row[species + 2]],
                pressure_pa: row[species + 3],
            })
    }

    fn row(&self, axis: usize, sign: f64, species: usize) -> Option<&[f64]> {
        let width = species + 4;
        self.slots[side_index(axis, sign)]
            .map(|index| &self.states[index * width..(index + 1) * width])
    }

    pub(super) fn boundary_face(
        &self,
        side: BoundarySide,
        face: usize,
        species: usize,
    ) -> Option<ReservoirBoundary<'_>> {
        if let Some(profile) = &self.profiles
            && let Some(row) = profile.row(side, face, species)
        {
            return Some(ReservoirBoundary {
                side,
                partial_densities: &row[..species],
                velocity: [row[species], row[species + 1], row[species + 2]],
                pressure_pa: row[species + 3],
            });
        }
        self.boundary(side, species)
    }

    pub(super) fn primitive_face(&self, side: BoundarySide, face: usize) -> Option<Primitive> {
        self.primitive_face_at(side, face, false)
    }

    pub(super) fn primitive_face_at(
        &self,
        side: BoundarySide,
        face: usize,
        trial_profiles: bool,
    ) -> Option<Primitive> {
        self.profiles
            .as_ref()
            .and_then(|profile| profile.primitive_at(side, face, trial_profiles))
            .or_else(|| {
                let slot = side_index(side.axis(), side.sign());
                self.finite_primitive(slot)
                    .or_else(|| self.slots[slot].map(|index| self.primitive[index]))
            })
    }

    pub(super) fn normal_face(
        &self,
        axis: usize,
        sign: f64,
        face: usize,
        species: usize,
    ) -> Option<SpatialPrimitive<'_>> {
        self.normal_face_at(axis, sign, face, species, false)
    }

    pub(super) fn normal_face_at(
        &self,
        axis: usize,
        sign: f64,
        face: usize,
        species: usize,
        trial_profiles: bool,
    ) -> Option<SpatialPrimitive<'_>> {
        let side = BoundarySide::ALL[side_index(axis, sign)];
        if let Some(profile) = &self.profiles
            && let Some(row) = profile.row_at(side, face, species, trial_profiles)
        {
            return Some(normal_primitive(row, species, axis));
        }
        let slot = side_index(axis, sign);
        if let Some(state) = self.finite.state(slot, false) {
            return Some(SpatialPrimitive {
                partial_densities: state.partial_densities,
                normal_velocity: 0.0,
                transverse_velocity: [0.0; 2],
                pressure: self.finite.primitive(slot, false)?.pressure,
            });
        }
        self.row(axis, sign, species)
            .map(|row| normal_primitive(row, species, axis))
    }

    fn finite_primitive(&self, slot: usize) -> Option<Primitive> {
        self.finite.primitive(slot, false).map(|p| Primitive {
            density: p.density,
            velocity: [0.0; 3],
            pressure: p.pressure,
            temperature: p.temperature,
            sound_speed: p.sound_speed,
        })
    }

    /// Additional storage beyond the empty finite metadata already in the field.
    pub(super) fn extra_workspace_bytes(&self) -> usize {
        self.workspace_bytes() - size_of::<FiniteReservoirs<3, 6>>()
    }

    pub(super) fn maximum_rate(&self, grid: &CartesianGrid) -> f64 {
        self.maximum_rate_at(grid, false)
    }

    pub(super) fn maximum_rate_at(&self, grid: &CartesianGrid, trial_profiles: bool) -> f64 {
        self.slots
            .iter()
            .enumerate()
            .filter_map(|(slot, index)| {
                self.finite_primitive(slot)
                    .or_else(|| index.map(|index| self.primitive[index]))
            })
            .map(|p| {
                (0..3)
                    .map(|axis| (p.velocity[axis].abs() + p.sound_speed) / grid.spacing[axis])
                    .sum::<f64>()
            })
            .fold(0.0, f64::max)
            .max(
                self.profiles
                    .as_ref()
                    .map_or(0.0, |p| p.maximum_rate(trial_profiles)),
            )
    }
}

fn side_index(axis: usize, sign: f64) -> usize {
    2 * axis + usize::from(sign > 0.0)
}
