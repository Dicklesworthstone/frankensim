//! Immutable source-admitted exterior composition on physical open sides.

use super::*;
use crate::compressible::AxisBoundaries;

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
    states: Vec<f64>,
    primitive: Vec<Primitive>,
}

impl Reservoirs {
    pub(super) fn required_bytes(gas: &MixtureGas, count: usize) -> Result<usize, Error> {
        if count > 6 {
            return Err(Error::Invalid("at most six volumetric mixture reservoir sides"));
        }
        count
            .checked_mul((gas.species_count() + 4) * size_of::<f64>() + size_of::<Primitive>())
            .ok_or(Error::Invalid("volumetric reservoir workspace size overflow"))
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
            states: allocated(inputs.len() * width, 0.0)?,
            primitive: allocated(inputs.len(), Primitive::default())?,
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
            let original = gas.spatial_conserved_into(
                normal_primitive(row, n, 0),
                &mut conservative[..width],
            )?;
            result.primitive[index] = Primitive {
                density: original.density,
                velocity: input.velocity,
                pressure: original.pressure,
                temperature: original.temperature,
                sound_speed: original.sound_speed,
            };
            result.slots[slot] = Some(index);
        }
        stepping::poll::<Error>(cancelled)?;
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
        stepping::poll::<Error>(cancelled)?;
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
                velocity: [row[species], row[species + 1], row[species + 2]],
                pressure_pa: row[species + 3],
            })
    }

    fn row(&self, axis: usize, sign: f64, species: usize) -> Option<&[f64]> {
        let width = species + 4;
        self.slots[side_index(axis, sign)]
            .map(|index| &self.states[index * width..(index + 1) * width])
    }

    pub(super) fn normal(
        &self,
        axis: usize,
        sign: f64,
        species: usize,
    ) -> Option<SpatialPrimitive<'_>> {
        self.row(axis, sign, species)
            .map(|row| normal_primitive(row, species, axis))
    }

    pub(super) fn maximum_rate(&self, grid: &CartesianGrid) -> f64 {
        self.primitive
            .iter()
            .map(|p| {
                (0..3)
                    .map(|axis| (p.velocity[axis].abs() + p.sound_speed) / grid.spacing[axis])
                    .sum::<f64>()
            })
            .fold(0.0, f64::max)
    }
}

fn side_index(axis: usize, sign: f64) -> usize {
    2 * axis + usize::from(sign > 0.0)
}
