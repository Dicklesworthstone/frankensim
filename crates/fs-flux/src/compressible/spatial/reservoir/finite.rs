//! Finite rigid tanks on the original Cartesian gas transaction.
//!
//! A side transfers the exact native face mass and total-energy impulses to
//! one acoustically compact plenum. Its pressure is recovered from the same
//! frozen source at both SSP stages. The tank's stationary support receives
//! the Cartesian momentum impulse; it performs no work.

use super::*;
use fs_gas::dynamics::network::Plenum;

/// Replace one held stationary reservoir with an independent finite rigid tank.
///
/// The existing uniform binding supplies its initial density and pressure.
/// All faces on this side communicate with this same well-mixed volume.
/// Different entries represent independent tanks, not shared storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoirBoundary {
    /// Already bound stationary uniform exterior on a transmissive side.
    pub side: BoundarySide,
    /// Physical tank volume [m³], strictly positive.
    pub volume_m3: f64,
}

/// Accepted tank exchange over the same interval as the field's Euler report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoirReport {
    /// Source-compatible stationary tank before the accepted interval.
    pub before: Plenum,
    /// Source-compatible stationary tank after the accepted interval.
    pub after: Plenum,
    /// Outward field transfer integrated with the native SSP weights.
    ///
    /// Mass [kg] and referenced total energy [J] enter the tank. All incoming
    /// kinetic energy remains in its internal energy. Cartesian momentum
    /// [N s] goes to the stationary tank support, not a resolved tank velocity.
    pub transfer: BoundaryImpulse,
    /// `m_after - m_before - transfer.mass` [kg].
    pub mass_residual_kg: f64,
    /// `U_after - U_before - transfer.energy` in the source energy gauge [J].
    pub energy_residual_j: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct FiniteReservoirs {
    pub(super) accepted: [Option<Plenum>; 6],
    pub(super) trial: [Option<Plenum>; 6],
    pub(super) traces: [Option<AdmittedReservoir>; 6],
    reports: [Option<FiniteReservoirReport>; 6],
    trial_reports: [Option<FiniteReservoirReport>; 6],
}

impl Reservoirs {
    fn bind_finite(
        &mut self,
        gas: &FrozenGas,
        inputs: &[FiniteReservoirBoundary],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if inputs.len() > 6 || self.finite.accepted.iter().any(Option::is_some) {
            return Err(Error::Invalid(
                "finite reservoirs bind at most six independent tanks once",
            ));
        }
        let mut candidate = FiniteReservoirs::default();
        let mut compact = self.stagnation;
        for input in inputs {
            poll(cancelled)?;
            let axis = input.side.axis();
            let sign = input.side.sign();
            let slot = 2 * axis + usize::from(sign > 0.0);
            let exterior = self.get(axis, sign).ok_or(Error::Invalid(
                "finite tank needs an already bound uniform reservoir",
            ))?;
            if candidate.accepted[slot].is_some() || exterior.input.velocity != [0.0; 3] {
                return Err(Error::Invalid(
                    "finite reservoir sides must be unique and exactly stationary",
                ));
            }
            let state = fs_gas::dynamics::Conserved {
                density: exterior.state.density,
                momentum: 0.0,
                energy: exterior.state.energy,
            };
            let tank = Plenum::new(input.volume_m3, state)?.exchange_candidate(
                gas,
                0.0,
                0.0,
                &mut *cancelled,
            )?;
            // Admission includes the original compact boundary's source path.
            gas.reservoir_spatial_flux(exterior.state.normal(axis), state, sign, &mut *cancelled)?;
            candidate.accepted[slot] = Some(tank);
            candidate.trial[slot] = Some(tank);
            candidate.traces[slot] = Some(*exterior);
            compact[slot] = Some(state);
        }
        poll(cancelled)?;
        self.finite = candidate;
        self.stagnation = compact;
        Ok(())
    }

    /// The first native capture stores half of the forward-Euler transfer.
    /// Weight two creates the full forward tank stage; weight one creates the
    /// final SSP candidate from the accumulated half-weighted face impulses.
    pub(in crate::compressible::spatial) fn stage_finite(
        &mut self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        boundary: &[BoundaryImpulse],
        weight: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for side in BoundarySide::ALL {
            let slot = 2 * side.axis() + usize::from(side.sign() > 0.0);
            let Some(before) = self.finite.accepted[slot] else {
                continue;
            };
            poll(cancelled)?;
            let sum = sum_arrays(
                boundary[grid.boundary_range(side)]
                    .iter()
                    .map(|q| q.array()),
                cancelled,
            )?;
            let transfer = Totals::from_array(sum.map(|q| weight * q));
            if !transfer.array().into_iter().all(f64::is_finite) {
                return Err(Error::Invalid("unrepresentable finite reservoir transfer"));
            }
            // The field already includes enthalpy and kinetic-energy flux.
            // Adding p*dV or another enthalpy term here would double count it.
            let after =
                before.exchange_candidate(gas, transfer.mass, transfer.energy, &mut *cancelled)?;
            let q = after.state();
            let state = State {
                density: q.density,
                momentum: [0.0; 3],
                energy: q.energy,
            };
            let primitive = state.primitive(gas, &mut *cancelled)?;
            let report = FiniteReservoirReport {
                before,
                after,
                transfer,
                mass_residual_kg: q.density * after.volume_m3()
                    - before.state().density * before.volume_m3()
                    - transfer.mass,
                energy_residual_j: q.energy * after.volume_m3()
                    - before.state().energy * before.volume_m3()
                    - transfer.energy,
            };
            if !report.mass_residual_kg.is_finite() || !report.energy_residual_j.is_finite() {
                return Err(Error::Invalid("unrepresentable finite reservoir balance"));
            }
            self.finite.trial[slot] = Some(after);
            self.finite.traces[slot] = Some(AdmittedReservoir {
                input: ReservoirFaceState {
                    density: primitive.density,
                    velocity: [0.0; 3],
                    pressure_pa: primitive.pressure,
                },
                state,
                primitive,
            });
            self.finite.trial_reports[slot] = Some(report);
        }
        Ok(())
    }

    /// Called only after the field's last fallible operation and final poll.
    pub(in crate::compressible::spatial) fn publish_finite(&mut self) {
        for slot in 0..6 {
            if let Some(tank) = self.finite.trial[slot] {
                self.finite.accepted[slot] = Some(tank);
                self.states[slot] = self.finite.traces[slot];
                self.stagnation[slot] = Some(tank.state());
                self.finite.reports[slot] = self.finite.trial_reports[slot];
            }
        }
    }

    /// A complete split retry starts from accepted finite inventories as well
    /// as accepted field cells. This also restores the physical exterior cache.
    pub(in crate::compressible::spatial) fn copy_finite_from(&mut self, source: &Self) {
        self.finite = source.finite;
        for slot in 0..6 {
            if source.finite.accepted[slot].is_some() {
                self.states[slot] = source.states[slot];
                self.stagnation[slot] = source.stagnation[slot];
            }
        }
    }
}

impl Euler3d {
    /// Couple selected existing uniform reservoirs to finite rigid gas storage.
    ///
    /// Each binding becomes a compact, perfectly mixed tank using this field's
    /// original frozen composition and caloric law. Pressure and temperature
    /// evolve from the same signed mass and total-energy impulses as the field.
    /// Native SSP RK2 advances tank and field together; the intermediate tank
    /// pressure drives the second residual. The timestep also respects the
    /// tank's summed aperture-wave rate, so small volumes cannot evade CFL.
    ///
    /// Bind once before the first step. This selects the compact law on those
    /// sides; other boundary laws and profiles remain unchanged. Independent
    /// tanks may fill or drain, and incoming shear kinetic energy thermalizes.
    /// No tank heat transfer, volume motion, internal waves or mixing time is
    /// resolved. Explicit molecular bath exchange remains unavailable because
    /// a stagnation state does not supply a physical molecular face state.
    ///
    /// Accepted cells, tanks, clock, exterior caches and interval transfers
    /// publish together. All stage storage is fixed-size and allocated inline.
    pub fn with_finite_reservoirs(
        mut self,
        inputs: &[FiniteReservoirBoundary],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if self.time != 0.0 || self.last_step.is_some() {
            return Err(Error::Invalid(
                "finite reservoirs bind before the first flow step",
            ));
        }
        self.reservoirs
            .bind_finite(&self.gas, inputs, &mut cancelled)?;
        Ok(self)
    }

    /// Accepted finite gas inventory on this side, in the field's source gauge.
    pub fn finite_reservoir(&self, side: BoundarySide) -> Option<&Plenum> {
        self.reservoirs.finite.accepted[2 * side.axis() + usize::from(side.sign() > 0.0)].as_ref()
    }

    /// Last accepted tank and fixed-support exchange on this side.
    pub fn finite_reservoir_report(&self, side: BoundarySide) -> Option<&FiniteReservoirReport> {
        self.reservoirs.finite.reports[2 * side.axis() + usize::from(side.sign() > 0.0)].as_ref()
    }
}

#[cfg(test)]
mod tests;
