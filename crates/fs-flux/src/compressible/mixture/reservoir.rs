//! Immutable source-admitted boundary gas shared by advection and molecules.
use super::*;
use super::finite::FiniteReservoirs;
mod finite;
pub use finite::FiniteReservoir;
use crate::compressible::molecular::MolecularWall2d;
use crate::compressible::{AxisBoundaries, SIDES};
use fs_gas::multicomponent::dissipation::DissipationError;

/// Uniform externally maintained state on one physical boundary face plane.
/// Species densities use the field's exact canonical source order. The original
/// gas derives temperature, caloric precision and energy references; no second
/// EOS or user-supplied heat capacity/sound speed is admitted.
///
/// These are infinite controlled reservoirs, not evolving finite stores. Their
/// velocity is Cartesian; pressure is absolute. A selected side must be an open
/// (`Transmissive`) side of the Cartesian topology. The original Riemann solve
/// determines inward/outward advection from the actual pair of traces; this is
/// not a nonreflecting or characteristic-pressure-only boundary condition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirBoundary<'a> {
    pub side: BoundarySide,
    /// Original canonical partial mass densities [kg/m³].
    pub partial_densities: &'a [f64],
    /// Cartesian velocity [m/s].
    pub velocity: [f64; 2],
    /// Absolute pressure [Pa].
    pub pressure_pa: f64,
}

/// Flat source-admitted primitive rows; the immutable exterior numerical data
/// and primitive cache are included in the field's actual byte account.
pub(super) struct Reservoirs {
    slots: [Option<usize>; 4],
    stagnation: [bool; 4],
    reservoir_energy: [f64; 4],
    states: Vec<f64>,
    primitive: Vec<Primitive>,
    pub(super) finite: FiniteReservoirs<2, 4>,
}
impl Reservoirs {
    pub(super) fn required_bytes(gas: &MixtureGas, count: usize) -> Result<usize, Error> {
        if count > 4 {
            return Err(Error::Invalid("at most four Cartesian reservoir sides"));
        }
        count
            .checked_mul((gas.species_count() + 3) * size_of::<f64>() + size_of::<Primitive>())
            .ok_or(Error::Invalid(
                "Cartesian reservoir workspace size overflow",
            ))
    }

    pub(super) fn new(
        gas: &MixtureGas,
        topology: Boundaries2d,
        inputs: &[ReservoirBoundary<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::required_bytes(gas, inputs.len())?;
        let n = gas.species_count();
        let width = n + 3;
        let mut result = Self {
            slots: [None; 4],
            stagnation: [false; 4],
            reservoir_energy: [0.0; 4],
            states: allocated(inputs.len() * width, 0.0)?,
            primitive: allocated(inputs.len(), Primitive::default())?,
            finite: FiniteReservoirs::empty(),
        };
        let mut conservative = [0.0; WIDTH_LIMIT];
        for (index, input) in inputs.iter().enumerate() {
            cartesian::poll::<Error>(cancelled)?;
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
                    "reservoir needs a unique open side and original species axis",
                ));
            }
            let row = &mut result.states[index * width..(index + 1) * width];
            row[..n].copy_from_slice(input.partial_densities);
            row[n..n + 2].copy_from_slice(&input.velocity);
            row[n + 2] = input.pressure_pa;
            let original =
                gas.normal_conserved_into(normal_primitive(row, n, 0), &mut conservative[..width])?;
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
        cartesian::poll::<Error>(cancelled)?;
        Ok(result)
    }

    pub(super) fn set_stagnation(
        &mut self,
        gas: &MixtureGas,
        sides: &[BoundarySide],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        cartesian::poll::<Error>(cancelled)?;
        if sides.len() > 4 {
            return Err(Error::Invalid("too many compact reservoir sides"));
        }
        let n = gas.species_count();
        let mut candidate = std::array::from_fn(|slot| self.finite.has_side(slot));
        let mut seen = [false; 4];
        let mut flux = [0.0; WIDTH_LIMIT];
        for &side in sides {
            cartesian::poll::<Error>(cancelled)?;
            let axis = side.axis();
            let sign = side.sign();
            let slot = side_index(axis, sign);
            let index = self.slots[slot].ok_or(Error::Invalid(
                "compact reservoir needs an already bound uniform exterior",
            ))?;
            if seen[slot] || self.primitive[index].velocity != [0.0; 2] {
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
            gas.reservoir_normal_flux_from_primitive(
                normal_primitive(row, n, axis),
                reservoir,
                sign,
                &mut flux[..n + 3],
                &mut *cancelled,
            )?;
            candidate[slot] = true;
            seen[slot] = true;
        }
        cartesian::poll::<Error>(cancelled)?;
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
        trial: bool,
    ) -> Option<fs_gas::multicomponent::Cell<'_>> {
        let slot = side_index(axis, sign);
        if let Some(state) = self.finite.state(slot, trial) {
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
        let width = n + 3;
        let mut flux = [0.0; WIDTH_LIMIT];
        let mut maximum: f64 = 0.0;
        for side in SIDES {
            let axis = side.axis();
            let sign = side.sign();
            let Some(reservoir) = self.stagnation_cell(axis, sign, n) else {
                continue;
            };
            let position = if sign < 0.0 { 0 } else { grid.shape[axis] - 1 };
            let (mut wave_sum, mut wave_correction) = (0.0, 0.0);
            for face in 0..grid.face_count(side) {
                let k = if axis == 0 {
                    face * grid.shape[0] + position
                } else {
                    position * grid.shape[0] + face
                };
                cartesian::poll::<Error>(cancelled)?;
                let q = cell(&states[k * width..(k + 1) * width], n);
                let face = gas.reservoir_normal_flux(
                    q.normal(axis),
                    reservoir,
                    sign,
                    &mut flux[..width],
                    &mut *cancelled,
                )?;
                let wave = grid.face_area(side) * face.signal_speed;
                let y = wave - wave_correction;
                let next = wave_sum + y;
                wave_correction = (next - wave_sum) - y;
                wave_sum = next;
                let density = sum(q.partial_densities.iter().copied());
                let local = (0..2)
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
            }
            if let Some(volume) = self.finite.volume(side_index(axis, sign)) {
                let rate = 0.5 * wave_sum / volume;
                if !rate.is_finite() {
                    return Err(Error::Invalid("unrepresentable finite mixture aperture-wave rate"));
                }
                maximum = maximum.max(rate);
            }
        }
        Ok(maximum)
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        self.states.capacity() * size_of::<f64>()
            + self.primitive.capacity() * size_of::<Primitive>()
            + self.finite.workspace_bytes()
    }

    pub(super) fn try_clone(&self, gas: &MixtureGas, cancelled: &mut impl FnMut() -> bool) -> Result<Self, Error> {
        let mut states = allocated(self.states.len(), 0.0)?;
        copy_values(&self.states, &mut states, cancelled)?;
        let mut primitive = allocated(self.primitive.len(), Primitive::default())?;
        primitive.copy_from_slice(&self.primitive);
        cartesian::poll::<Error>(cancelled)?;
        Ok(Self {
            slots: self.slots,
            stagnation: self.stagnation,
            reservoir_energy: self.reservoir_energy,
            states,
            primitive,
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
                velocity: [0.0; 2],
                pressure_pa: self.finite.primitive(slot, false).expect("accepted finite primitive").pressure,
            });
        }
        self.row(side.axis(), side.sign(), species)
            .map(|row| ReservoirBoundary {
                side,
                partial_densities: &row[..species],
                velocity: [row[species], row[species + 1]],
                pressure_pa: row[species + 2],
            })
    }

    fn row(&self, axis: usize, sign: f64, species: usize) -> Option<&[f64]> {
        let width = species + 3;
        self.slots[side_index(axis, sign)]
            .map(|index| &self.states[index * width..(index + 1) * width])
    }

    pub(super) fn normal(
        &self,
        axis: usize,
        sign: f64,
        species: usize,
    ) -> Option<NormalPrimitive<'_>> {
        let slot = side_index(axis, sign);
        if let Some(state) = self.finite.state(slot, false) {
            return Some(NormalPrimitive {
                partial_densities: state.partial_densities,
                normal_velocity: 0.0,
                transverse_velocity: 0.0,
                pressure: self.finite.primitive(slot, false).expect("accepted finite primitive").pressure,
            });
        }
        self.row(axis, sign, species)
            .map(|row| normal_primitive(row, species, axis))
    }

    pub(super) fn maximum_rate(&self, grid: &CartesianGrid) -> f64 {
        self.slots.iter().enumerate().filter_map(|(slot, index)| {
            index.map(|index| {
                if let Some(p) = self.finite.primitive(slot, false) {
                    p.sound_speed / grid.spacing[0] + p.sound_speed / grid.spacing[1]
                } else {
                    let p = self.primitive[index];
                    (p.velocity[0].abs() + p.sound_speed) / grid.spacing[0]
                        + (p.velocity[1].abs() + p.sound_speed) / grid.spacing[1]
                }
            })
        }).fold(0.0, f64::max)
    }

    /// An Euler continuation or imposed pressure alone defines no molecular bath.
    pub(super) fn admit_molecular(&self, topology: Boundaries2d) -> Result<(), DissipationError> {
        if self.has_stagnation() {
            return Err(DissipationError::Invalid(
                "compact stagnation reservoirs do not specify a molecular face bath",
            ));
        }
        for side in SIDES {
            if let AxisBoundaries::Ends { lower, upper } = topology.axis(side.axis()) {
                let condition = if side.sign() < 0.0 { lower } else { upper };
                if condition != BoundaryCondition::Reflecting
                    && self.slots[side_index(side.axis(), side.sign())].is_none()
                {
                    return Err(DissipationError::Invalid(
                        "open molecular side needs an explicit original-source reservoir",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Cell-centered derivative ghost obtained by reflection about the SAME
    /// fixed physical-face velocity. No exterior thermodynamic ghost is made.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn velocity_neighbour(
        &self,
        grid: &CartesianGrid,
        topology: Boundaries2d,
        walls: &[MolecularWall2d; 4],
        primitive: &[Primitive],
        cell: usize,
        axis: usize,
        upper: bool,
    ) -> Result<[f64; 2], DissipationError> {
        let position = if axis == 0 {
            cell % grid.shape[0]
        } else {
            cell / grid.shape[0]
        };
        let edge = if upper {
            position + 1 == grid.shape[axis]
        } else {
            position == 0
        };
        let sign = if upper { 1.0 } else { -1.0 };
        let velocity = if edge && let Some(index) = self.slots[side_index(axis, sign)] {
            std::array::from_fn(|component| {
                2.0 * self.primitive[index].velocity[component]
                    - primitive[cell].velocity[component]
            })
        } else {
            let (other, reflect) = cartesian::neighbour(grid, topology, cell, axis, upper);
            if reflect {
                walls[side_index(axis, sign)]
                    .velocity_ghost(axis, 0.5 * grid.spacing[axis], primitive[cell].velocity)
                    .unwrap_or_else(|| {
                        std::array::from_fn(|component| {
                            primitive[other].velocity[component]
                                * if component == axis { -1.0 } else { 1.0 }
                        })
                    })
            } else {
                primitive[other].velocity
            }
        };
        if !velocity.iter().all(|value| value.is_finite()) {
            return Err(DissipationError::Invalid(
                "molecular wall/reservoir ghost velocity overflow",
            ));
        }
        Ok(velocity)
    }
}
fn side_index(axis: usize, sign: f64) -> usize {
    2 * axis + usize::from(sign > 0.0)
}
