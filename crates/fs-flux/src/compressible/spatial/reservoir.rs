//! Prescribed open boundaries on the original three-dimensional pressure-wave owner.

use super::*;

mod finite;
pub use finite::{FiniteReservoirBoundary, FiniteReservoirReport};
mod profiles;
pub use profiles::ReservoirFaceState;
pub(super) use profiles::TimedSource;

/// Uniform externally maintained gas on one open physical side.
///
/// The same source and frozen composition as the field supply temperature,
/// internal energy and sound speed. Density and absolute pressure must jointly
/// lie inside that source's caloric domain. Velocity is Cartesian, not outward
/// normal velocity; the face kernel rotates it into the original normal frame.
///
/// This is a complete-state Riemann boundary, with no finite reservoir storage,
/// characteristic-pressure-only condition or nonreflecting claim. The optional
/// molecular operator retains zero molecular flux on open sides by default;
/// its explicit `with_reservoir_exchange` option also uses this complete trace
/// as a conducting and viscous exterior bath.
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

/// Nonnegative absolute backpressure on a side carrying outgoing gas.
///
/// The original gas source supplies the shock or expansion and samples it at
/// the stationary face. The interior supplies entropy and both shear velocities.
/// If the contact reverses, the step refuses because pressure alone cannot
/// specify incoming temperature or composition; use a reservoir for reversal.
/// A fixed backpressure is not a nonreflecting or molecular boundary condition.
/// Zero denotes exterior vacuum with a source-supported positive-density sonic
/// face. Exterior expansion, internal vacuum and rarefied flow are not resolved.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureOutletBoundary {
    /// One physical side whose topology is explicitly open.
    pub side: BoundarySide,
    /// Externally maintained absolute backpressure [Pa].
    pub pressure_pa: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct AdmittedReservoir {
    pub(super) input: ReservoirFaceState,
    pub(super) state: State,
    primitive: Primitive,
}
impl AdmittedReservoir {
    pub(super) fn primitive(&self) -> Primitive {
        self.primitive
    }

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

/// Uniform sides keep fixed metadata; optional profiles have admitted reusable
/// face buffers. Neither representation allocates during a flow step.
pub(super) struct Reservoirs {
    states: [Option<AdmittedReservoir>; 6],
    outlets: [Option<PressureOutletBoundary>; 6],
    characteristic_outlets: [Option<(f64, f64)>; 6],
    stagnation: [Option<fs_gas::dynamics::Conserved>; 6],
    finite: finite::FiniteReservoirs,
    pub(super) profiles: profiles::Profiles,
}

/// Immutable view of one complete exterior stage. Uniform bindings never
/// change; only prepared face profiles can select uncommitted sample storage.
#[derive(Clone, Copy)]
pub(super) struct View<'a> {
    owner: &'a Reservoirs,
    staged_profiles: bool,
    staged_finite: bool,
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
            characteristic_outlets: [None; 6],
            stagnation: [None; 6],
            finite: finite::FiniteReservoirs::default(),
            profiles: profiles::Profiles::default(),
        };
        for side in BoundarySide::ALL {
            let pair = topology.axis(side.axis());
            pair.admit::<Error>()?;
            if let AxisBoundaries::Ends { lower, upper } = pair
                && let BoundaryCondition::PressureOutlet { pressure_pa } =
                    (if side.sign() < 0.0 { lower } else { upper })
            {
                result.outlets[2 * side.axis() + usize::from(side.sign() > 0.0)] =
                    Some(PressureOutletBoundary { side, pressure_pa });
            }
            if let AxisBoundaries::Ends { lower, upper } = pair
                && let BoundaryCondition::CharacteristicOutlet {
                    reference_pressure_pa,
                    reference_velocity_m_s,
                } = (if side.sign() < 0.0 { lower } else { upper })
            {
                result.characteristic_outlets[2 * side.axis() + usize::from(side.sign() > 0.0)] =
                    Some((reference_pressure_pa, reference_velocity_m_s));
            }
        }
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
                input: ReservoirFaceState {
                    density: input.density,
                    velocity: input.velocity,
                    pressure_pa: input.pressure_pa,
                },
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
                || !input.pressure_pa.is_finite()
                || input.pressure_pa < 0.0
            {
                return Err(Error::Invalid(
                    "pressure outlet needs a unique open side and nonnegative absolute pressure",
                ));
            }
            result.outlets[index] = Some(input);
        }
        poll(cancelled)?;
        Ok(result)
    }

    pub(super) fn set_stagnation(
        &mut self,
        gas: &FrozenGas,
        sides: &[BoundarySide],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if sides.len() > 6 {
            return Err(Error::Invalid("at most six compact reservoir sides"));
        }
        let mut candidate = [None; 6];
        for &side in sides {
            poll(cancelled)?;
            let axis = side.axis();
            let sign = side.sign();
            let slot = 2 * axis + usize::from(sign > 0.0);
            let exterior = self.get(axis, sign).ok_or(Error::Invalid(
                "compact reservoir needs an already bound uniform exterior",
            ))?;
            if candidate[slot].is_some() || exterior.input.velocity != [0.0; 3] {
                return Err(Error::Invalid(
                    "compact reservoir sides must be unique and exactly stationary",
                ));
            }
            let reservoir = fs_gas::dynamics::Conserved {
                density: exterior.state.density,
                momentum: 0.0,
                energy: exterior.state.energy,
            };
            gas.reservoir_spatial_flux(
                exterior.state.normal(axis),
                reservoir,
                sign,
                &mut *cancelled,
            )?;
            candidate[slot] = Some(reservoir);
        }
        if self
            .finite
            .accepted
            .iter()
            .zip(candidate)
            .any(|(tank, law)| tank.is_some() && law.is_none())
        {
            return Err(Error::Invalid(
                "finite reservoirs retain their compact gas law",
            ));
        }
        poll(cancelled)?;
        self.stagnation = candidate;
        Ok(())
    }

    pub(super) fn is_stagnation(&self, axis: usize, sign: f64) -> bool {
        self.stagnation[2 * axis + usize::from(sign > 0.0)].is_some()
    }

    pub(super) fn has_stagnation(&self) -> bool {
        self.stagnation.iter().any(Option::is_some)
    }

    pub(super) fn get(&self, axis: usize, sign: f64) -> Option<&AdmittedReservoir> {
        self.states[2 * axis + usize::from(sign > 0.0)].as_ref()
    }

    pub(super) fn get_at(
        &self,
        grid: &CartesianGrid,
        axis: usize,
        sign: f64,
        cell: usize,
    ) -> Option<&AdmittedReservoir> {
        self.accepted().get_at(grid, axis, sign, cell)
    }

    pub(super) fn accepted(&self) -> View<'_> {
        View {
            owner: self,
            staged_profiles: false,
            staged_finite: false,
        }
    }

    pub(super) fn staged(&self) -> View<'_> {
        View {
            owner: self,
            staged_profiles: true,
            staged_finite: false,
        }
    }

    pub(super) fn try_clone(&self, cancelled: &mut impl FnMut() -> bool) -> Result<Self, Error> {
        Ok(Self {
            states: self.states,
            outlets: self.outlets,
            characteristic_outlets: self.characteristic_outlets,
            stagnation: self.stagnation,
            finite: self.finite,
            profiles: self.profiles.try_clone(cancelled)?,
        })
    }

    pub(super) fn outlet(&self, axis: usize, sign: f64) -> Option<PressureOutletBoundary> {
        self.outlets[2 * axis + usize::from(sign > 0.0)]
    }

    pub(super) fn maximum_rate(
        &self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        states: &[State],
        hydrostatic: Option<&hydrostatic::HydrostaticBalance>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        self.accepted()
            .maximum_rate(gas, grid, states, hydrostatic, cancelled)
    }
}

impl<'a> View<'a> {
    pub(super) fn with_finite_trial(mut self) -> Self {
        self.staged_finite = true;
        self
    }

    pub(super) fn finite_volume(self, axis: usize, sign: f64) -> Option<f64> {
        self.owner.finite.accepted[2 * axis + usize::from(sign > 0.0)].map(|tank| tank.volume_m3())
    }

    pub(super) fn has_stagnation(self) -> bool {
        self.owner.has_stagnation()
    }

    pub(super) fn stagnation(self, axis: usize, sign: f64) -> Option<fs_gas::dynamics::Conserved> {
        let slot = 2 * axis + usize::from(sign > 0.0);
        if self.staged_finite
            && let Some(tank) = self.owner.finite.trial[slot]
        {
            return Some(tank.state());
        }
        self.owner.stagnation[slot]
    }

    pub(super) fn get_at(
        self,
        grid: &CartesianGrid,
        axis: usize,
        sign: f64,
        cell: usize,
    ) -> Option<&'a AdmittedReservoir> {
        let slot = 2 * axis + usize::from(sign > 0.0);
        if self.staged_finite
            && let Some(state) = self.owner.finite.traces[slot].as_ref()
        {
            return Some(state);
        }
        if let Some(state) = self.owner.states[slot].as_ref() {
            return Some(state);
        }
        let [a, b] = tangents(axis);
        let c = grid.coordinates(cell);
        let face = c[a] + grid.shape[a] * c[b];
        if self.staged_profiles {
            self.owner.profiles.get_staged(slot, face)
        } else {
            self.owner.profiles.get(slot, face)
        }
    }

    pub(super) fn outlet(self, axis: usize, sign: f64) -> Option<PressureOutletBoundary> {
        self.owner.outlet(axis, sign)
    }

    pub(super) fn maximum_rate(
        self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        states: &[State],
        hydrostatic: Option<&hydrostatic::HydrostaticBalance>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let mut rate = self
            .owner
            .states
            .iter()
            .enumerate()
            .filter_map(|(slot, state)| {
                if self.staged_finite {
                    self.owner.finite.traces[slot].as_ref().or(state.as_ref())
                } else {
                    state.as_ref()
                }
            })
            .map(|state| {
                (0..3)
                    .map(|axis| {
                        (state.primitive.velocity[axis].abs() + state.primitive.sound_speed)
                            / grid.spacing[axis]
                    })
                    .sum::<f64>()
            })
            .fold(0.0, f64::max);
        let profile_rate = if self.staged_profiles {
            self.owner.profiles.maximum_staged_rate(grid, cancelled)?
        } else {
            self.owner.profiles.maximum_rate(grid, cancelled)?
        };
        rate = rate.max(profile_rate);
        for side in BoundarySide::ALL {
            let axis = side.axis();
            let sign = side.sign();
            let slot = 2 * axis + usize::from(sign > 0.0);
            let compact = self.stagnation(axis, sign);
            let mut finite_rate = 0.0;
            let pressure_control = if let Some(outlet) = self.owner.outlets[slot] {
                Some((outlet.pressure_pa, None))
            } else {
                self.owner.characteristic_outlets[slot].map(|(p, v)| (p, Some(v)))
            };
            if compact.is_none() && pressure_control.is_none() {
                continue;
            }
            let [a, b] = tangents(axis);
            for face in 0..grid.face_count(side) {
                let mut c = [0; 3];
                c[axis] = if sign < 0.0 { 0 } else { grid.shape[axis] - 1 };
                c[a] = face % grid.shape[a];
                c[b] = face / grid.shape[a];
                let cell = grid.index(c);
                let q = states[cell];
                let (flux, shear) = if let Some(balance) = hydrostatic {
                    let p = balance.primitive(gas, grid, cell, q, cancelled)?;
                    let variables = balance.variables(grid, cell, p)?;
                    let reference = balance.trace(grid, cell, axis, sign);
                    let mut density = reference[0] * variables[0];
                    let mut pressure = reference[4] * variables[4];
                    let temperature = (pressure / density) / gas.gas_constant();
                    if !positive(density)
                        || !positive(pressure)
                        || temperature < gas.domain().minimum_temperature
                        || temperature > gas.domain().maximum_temperature
                        || !temperature.is_finite()
                    {
                        density = p.density;
                        pressure = p.pressure;
                    }
                    // The initial CFL/refusal must use the same physical face
                    // as the residual, not a stratified cell-average pressure.
                    let trace = SpatialPrimitive {
                        density,
                        normal_velocity: p.velocity[axis],
                        transverse_velocity: [p.velocity[a], p.velocity[b]],
                        pressure,
                    };
                    (
                        if let Some(reservoir) = compact {
                            gas.reservoir_spatial_flux_from_primitive(
                                trace,
                                reservoir,
                                sign,
                                &mut *cancelled,
                            )?
                        } else {
                            let (outlet_pressure, reference_velocity) =
                                pressure_control.expect("selected pressure control");
                            if let Some(velocity) = reference_velocity {
                                gas.characteristic_outlet_spatial_flux_from_primitive(
                                    trace,
                                    sign,
                                    outlet_pressure,
                                    velocity,
                                    &mut *cancelled,
                                )?
                            } else {
                                gas.pressure_outlet_spatial_flux_from_primitive(
                                    trace,
                                    sign,
                                    outlet_pressure,
                                    &mut *cancelled,
                                )?
                            }
                        },
                        trace.transverse_velocity,
                    )
                } else {
                    (
                        if let Some(reservoir) = compact {
                            gas.reservoir_spatial_flux(
                                q.normal(axis),
                                reservoir,
                                sign,
                                &mut *cancelled,
                            )?
                        } else {
                            let (outlet_pressure, reference_velocity) =
                                pressure_control.expect("selected pressure control");
                            if let Some(velocity) = reference_velocity {
                                gas.characteristic_outlet_spatial_flux(
                                    q.normal(axis),
                                    sign,
                                    outlet_pressure,
                                    velocity,
                                    &mut *cancelled,
                                )?
                            } else {
                                gas.pressure_outlet_spatial_flux(
                                    q.normal(axis),
                                    sign,
                                    outlet_pressure,
                                    &mut *cancelled,
                                )?
                            }
                        },
                        [q.momentum[a] / q.density, q.momentum[b] / q.density],
                    )
                };
                // The normal wave bound also bounds the face sound speed.
                // Include unchanged transverse advection so a pressure wave
                // bounds the entire unsplit interval before any split stage.
                let local = flux.signal_speed / grid.spacing[axis]
                    + (shear[0].abs() + flux.signal_speed) / grid.spacing[a]
                    + (shear[1].abs() + flux.signal_speed) / grid.spacing[b];
                rate = rate.max(local);
                if let Some(volume) = self.finite_volume(axis, sign) {
                    finite_rate += 0.5 * grid.face_area(side) * (flux.signal_speed / volume);
                    if !finite_rate.is_finite() {
                        return Err(Error::Invalid(
                            "unrepresentable finite reservoir signal rate",
                        ));
                    }
                }
            }
            rate = rate.max(finite_rate);
        }
        Ok(rate)
    }
}
