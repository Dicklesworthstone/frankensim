//! Held exterior states on the original three-dimensional pressure-wave owner.

use super::*;

/// Uniform externally maintained gas on one open physical side.
///
/// The same source and frozen composition as the field supply temperature,
/// internal energy and sound speed. Density and absolute pressure must jointly
/// lie inside that source's caloric domain. Velocity is Cartesian, not outward
/// normal velocity; the face kernel rotates it into the original normal frame.
///
/// This is a complete-state Riemann boundary, with no finite reservoir storage,
/// characteristic-pressure-only condition or nonreflecting claim. The optional
/// molecular operator retains its existing zero molecular flux on open sides;
/// these advective traces do not introduce a conducting or viscous boundary bath.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirBoundary {
    /// One physical side whose topology is explicitly open.
    pub side: BoundarySide,
    /// Mass density [kg/m³].
    pub density: f64,
    /// Cartesian velocity [m/s].
    pub velocity: [f64; 3],
    /// Absolute pressure [Pa].
    pub pressure_pa: f64,
}

#[derive(Clone, Copy)]
pub(super) struct AdmittedReservoir {
    pub(super) input: ReservoirBoundary,
    pub(super) state: State,
    primitive: Primitive,
}
impl AdmittedReservoir {
    pub(super) fn normal(&self, axis: usize) -> SpatialPrimitive {
        let [a, b] = tangents(axis);
        SpatialPrimitive {
            density: self.primitive.density,
            normal_velocity: self.primitive.velocity[axis],
            transverse_velocity: [self.primitive.velocity[a], self.primitive.velocity[b]],
            pressure: self.primitive.pressure,
        }
    }
}

/// Fixed six-slot metadata; no face-sized storage or hot-loop allocation.
#[derive(Clone)]
pub(super) struct Reservoirs([Option<AdmittedReservoir>; 6]);

impl Reservoirs {
    pub(super) fn new(
        gas: &FrozenGas,
        topology: Boundaries3d,
        inputs: &[ReservoirBoundary],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if inputs.len() > 6 {
            return Err(Error::Invalid("at most six volumetric reservoir sides"));
        }
        let mut result = Self([None; 6]);
        for &input in inputs {
            poll(cancelled)?;
            let axis = input.side.axis();
            let sign = input.side.sign();
            let slot = &mut result.0[2 * axis + usize::from(sign > 0.0)];
            let condition = match topology.axis(axis) {
                AxisBoundaries::Periodic => None,
                AxisBoundaries::Ends { lower, upper } => {
                    Some(if sign < 0.0 { lower } else { upper })
                }
            };
            if slot.is_some() || condition != Some(BoundaryCondition::Transmissive) {
                return Err(Error::Invalid(
                    "reservoir needs a unique open volumetric side",
                ));
            }
            let state =
                State::from_primitive(gas, input.density, input.velocity, input.pressure_pa)?;
            // Cache the same decoded state that native cell reconstruction
            // uses, so an identical field/reservoir pair keeps exact consistency.
            let primitive = state.primitive(gas, &mut *cancelled)?;
            *slot = Some(AdmittedReservoir {
                input,
                state,
                primitive,
            });
        }
        poll(cancelled)?;
        Ok(result)
    }

    pub(super) fn get(&self, axis: usize, sign: f64) -> Option<&AdmittedReservoir> {
        self.0[2 * axis + usize::from(sign > 0.0)].as_ref()
    }

    pub(super) fn maximum_rate(&self, grid: &CartesianGrid) -> f64 {
        self.0
            .iter()
            .flatten()
            .map(|state| {
                (0..3)
                    .map(|axis| {
                        (state.primitive.velocity[axis].abs() + state.primitive.sound_speed)
                            / grid.spacing[axis]
                    })
                    .sum::<f64>()
            })
            .fold(0.0, f64::max)
    }
}
