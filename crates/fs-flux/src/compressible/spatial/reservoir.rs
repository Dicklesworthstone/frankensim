//! Prescribed open boundaries on the original three-dimensional pressure-wave owner.

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

/// Positive absolute backpressure on a side carrying outgoing gas.
///
/// The original gas source supplies the shock or expansion and samples it at
/// the stationary face. The interior supplies entropy and both shear velocities.
/// If the contact reverses, the step refuses because pressure alone cannot
/// specify incoming temperature or composition; use a reservoir for reversal.
/// A fixed backpressure is not a nonreflecting or molecular boundary condition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureOutletBoundary {
    /// One physical side whose topology is explicitly open.
    pub side: BoundarySide,
    /// Externally maintained absolute backpressure [Pa].
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
pub(super) struct Reservoirs {
    states: [Option<AdmittedReservoir>; 6],
    outlets: [Option<PressureOutletBoundary>; 6],
}

impl Reservoirs {
    pub(super) fn new(
        gas: &FrozenGas,
        topology: Boundaries3d,
        inputs: &[ReservoirBoundary],
        outlets: &[PressureOutletBoundary],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if inputs.len() > 6 || outlets.len() > 6 || inputs.len() + outlets.len() > 6 {
            return Err(Error::Invalid("at most six prescribed volumetric sides"));
        }
        let mut result = Self {
            states: [None; 6],
            outlets: [None; 6],
        };
        for &input in inputs {
            poll(cancelled)?;
            let axis = input.side.axis();
            let sign = input.side.sign();
            let slot = &mut result.states[2 * axis + usize::from(sign > 0.0)];
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
        for &input in outlets {
            poll(cancelled)?;
            let axis = input.side.axis();
            let sign = input.side.sign();
            let index = 2 * axis + usize::from(sign > 0.0);
            let condition = match topology.axis(axis) {
                AxisBoundaries::Periodic => None,
                AxisBoundaries::Ends { lower, upper } => {
                    Some(if sign < 0.0 { lower } else { upper })
                }
            };
            if result.states[index].is_some()
                || result.outlets[index].is_some()
                || condition != Some(BoundaryCondition::Transmissive)
                || !positive(input.pressure_pa)
            {
                return Err(Error::Invalid(
                    "pressure outlet needs a unique open side and positive absolute pressure",
                ));
            }
            result.outlets[index] = Some(input);
        }
        poll(cancelled)?;
        Ok(result)
    }

    pub(super) fn get(&self, axis: usize, sign: f64) -> Option<&AdmittedReservoir> {
        self.states[2 * axis + usize::from(sign > 0.0)].as_ref()
    }

    pub(super) fn outlet(&self, axis: usize, sign: f64) -> Option<PressureOutletBoundary> {
        self.outlets[2 * axis + usize::from(sign > 0.0)]
    }

    pub(super) fn maximum_rate(
        &self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        states: &[State],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let mut rate = self
            .states
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
            .fold(0.0, f64::max);
        for outlet in self.outlets.iter().flatten() {
            let axis = outlet.side.axis();
            let sign = outlet.side.sign();
            let [a, b] = tangents(axis);
            for face in 0..grid.face_count(outlet.side) {
                let mut c = [0; 3];
                c[axis] = if sign < 0.0 { 0 } else { grid.shape[axis] - 1 };
                c[a] = face % grid.shape[a];
                c[b] = face / grid.shape[a];
                let q = states[grid.index(c)];
                let flux = gas.pressure_outlet_spatial_flux(
                    q.normal(axis),
                    sign,
                    outlet.pressure_pa,
                    &mut *cancelled,
                )?;
                // The normal wave bound also bounds the face sound speed.
                // Include unchanged transverse advection so a pressure wave
                // bounds the entire unsplit interval before any split stage.
                let local = flux.signal_speed / grid.spacing[axis]
                    + ((q.momentum[a] / q.density).abs() + flux.signal_speed) / grid.spacing[a]
                    + ((q.momentum[b] / q.density).abs() + flux.signal_speed) / grid.spacing[b];
                rate = rate.max(local);
            }
        }
        Ok(rate)
    }
}
