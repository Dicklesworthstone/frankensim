//! Immutable source-admitted boundary gas shared by advection and molecules.
use super::*;
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
    states: Vec<f64>,
    primitive: Vec<Primitive>,
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
            states: allocated(inputs.len() * width, 0.0)?,
            primitive: allocated(inputs.len(), Primitive::default())?,
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
            result.slots[slot] = Some(index);
        }
        cartesian::poll::<Error>(cancelled)?;
        Ok(result)
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        self.states.capacity() * size_of::<f64>()
            + self.primitive.capacity() * size_of::<Primitive>()
    }

    pub(super) fn try_clone(&self, cancelled: &mut impl FnMut() -> bool) -> Result<Self, Error> {
        let mut states = allocated(self.states.len(), 0.0)?;
        copy_values(&self.states, &mut states, cancelled)?;
        let mut primitive = allocated(self.primitive.len(), Primitive::default())?;
        primitive.copy_from_slice(&self.primitive);
        cartesian::poll::<Error>(cancelled)?;
        Ok(Self {
            slots: self.slots,
            states,
            primitive,
        })
    }

    pub(super) fn boundary(
        &self,
        side: BoundarySide,
        species: usize,
    ) -> Option<ReservoirBoundary<'_>> {
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
        self.row(axis, sign, species)
            .map(|row| normal_primitive(row, species, axis))
    }

    pub(super) fn maximum_rate(&self, grid: &CartesianGrid) -> f64 {
        self.primitive
            .iter()
            .map(|p| {
                (p.velocity[0].abs() + p.sound_speed) / grid.spacing[0]
                    + (p.velocity[1].abs() + p.sound_speed) / grid.spacing[1]
            })
            .fold(0.0, f64::max)
    }

    /// An Euler zero-gradient continuation alone defines no molecular bath.
    pub(super) fn admit_molecular(&self, topology: Boundaries2d) -> Result<(), DissipationError> {
        for side in SIDES {
            if let AxisBoundaries::Ends { lower, upper } = topology.axis(side.axis()) {
                let condition = if side.sign() < 0.0 { lower } else { upper };
                if condition == BoundaryCondition::Transmissive
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
        primitive: &[Primitive],
        cell: usize,
        axis: usize,
        component: usize,
        upper: bool,
    ) -> f64 {
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
        if edge && let Some(index) = self.slots[side_index(axis, sign)] {
            return 2.0 * self.primitive[index].velocity[component]
                - primitive[cell].velocity[component];
        }
        let (other, reflect) = cartesian::neighbour(grid, topology, cell, axis, upper);
        primitive[other].velocity[component]
            * if reflect && component == axis {
                -1.0
            } else {
                1.0
            }
    }
}
fn side_index(axis: usize, sign: f64) -> usize {
    2 * axis + usize::from(sign > 0.0)
}
