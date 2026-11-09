//! Finite rigid pressure reservoirs on the native planar Euler transaction.
//!
//! Each physical side communicates with one compact, perfectly mixed tank.
//! The original gas source supplies stagnation inflow, backpressure outflow,
//! choking and reversal. Its plenum receives the exact opposite mass and
//! total-energy transfer of the field; the fixed support receives momentum.

use super::*;
use fs_gas::dynamics::network::Plenum;

/// One independent source-compatible tank attached to a transmissive side.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoirBoundary {
    /// Physical side; all its faces share this same well-mixed inventory.
    pub side: BoundarySide,
    /// Finite stationary storage in this field's original frozen energy gauge.
    /// Construct its state through the same [`FrozenGas`] as the field.
    pub tank: Plenum,
}

/// Accepted tank and fixed-support exchange over the field's Euler interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoirReport {
    /// Original compact gas inventory.
    pub before: Plenum,
    /// Accepted compact gas inventory.
    pub after: Plenum,
    /// Outward field mass [kg], Cartesian support impulse [N s] and total
    /// energy [J], integrated with the native SSP weights. Energy includes
    /// incoming shear kinetic energy and flow work in the original gauge.
    pub transfer: BoundaryImpulse,
    /// `m_after - m_before - transfer.mass` [kg].
    pub mass_residual_kg: f64,
    /// `U_after - U_before - transfer.energy` [J] in the shared energy gauge.
    pub energy_residual_j: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct FiniteReservoirs {
    accepted: [Option<Plenum>; 4],
    trial: [Option<Plenum>; 4],
    reports: [Option<FiniteReservoirReport>; 4],
    trial_reports: [Option<FiniteReservoirReport>; 4],
}

impl FiniteReservoirs {
    pub(super) fn any(&self) -> bool {
        self.accepted.iter().any(Option::is_some)
    }

    pub(super) fn state(
        &self,
        axis: usize,
        sign: f64,
        trial: bool,
    ) -> Option<fs_gas::dynamics::Conserved> {
        let slot = 2 * axis + usize::from(sign > 0.0);
        (if trial {
            self.trial[slot]
        } else {
            self.accepted[slot]
        })
        .map(|tank| tank.state())
    }

    /// Include both the boundary-cell unsplit rate and the tank aperture rate
    /// before the first attempt, so tiny tanks do not spend the retry budget.
    pub(super) fn initial_rate(
        &self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        states: &[State],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let mut rate = 0.0_f64;
        let mut aperture = [0.0; 4];
        for (slot, side) in SIDES.into_iter().enumerate() {
            let Some(tank) = self.accepted[slot] else {
                continue;
            };
            let axis = side.axis();
            let end = if side.sign() < 0.0 {
                0
            } else {
                grid.shape[axis] - 1
            };
            let mut correction = 0.0;
            for line in 0..grid.face_count(side) {
                poll(cancelled)?;
                let cell = if axis == 0 {
                    line * grid.shape[0] + end
                } else {
                    end * grid.shape[0] + line
                };
                let q = states[cell];
                let p = q.primitive(gas, &mut *cancelled)?;
                let flux = gas.reservoir_normal_flux(
                    q.normal(axis),
                    tank.state(),
                    side.sign(),
                    &mut *cancelled,
                )?;
                rate = rate.max(
                    flux.signal_speed / grid.spacing[axis]
                        + (p.velocity[1 - axis].abs() + p.sound_speed) / grid.spacing[1 - axis],
                );
                let value = grid.face_area(side) * flux.signal_speed - correction;
                let next = aperture[slot] + value;
                correction = (next - aperture[slot]) - value;
                aperture[slot] = next;
            }
        }
        let rate = rate.max(self.aperture_rate(aperture, false)?);
        if !rate.is_finite() {
            return Err(Error::Invalid("unrepresentable finite reservoir face rate"));
        }
        Ok(rate)
    }

    pub(super) fn aperture_rate(&self, aperture: [f64; 4], trial: bool) -> Result<f64, Error> {
        let tanks = if trial { &self.trial } else { &self.accepted };
        let mut maximum = 0.0_f64;
        for (slot, tank) in tanks.iter().enumerate() {
            let Some(tank) = tank else { continue };
            // This is the compact source's 0.5 sum(A lambda)/V rate. The same
            // bound is recomputed with the evolved tank at the second stage.
            let rate = 0.5 * aperture[slot] / tank.volume_m3();
            if !rate.is_finite() || rate < 0.0 {
                return Err(Error::Invalid(
                    "unrepresentable finite reservoir aperture rate",
                ));
            }
            maximum = maximum.max(rate);
        }
        Ok(maximum)
    }

    /// First capture contains half a forward step; weight two creates the full
    /// tank stage. The final capture already contains both half-weighted fluxes.
    pub(super) fn stage(
        &mut self,
        gas: &FrozenGas,
        grid: &CartesianGrid,
        boundary: &[BoundaryImpulse],
        final_stage: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for (slot, side) in SIDES.into_iter().enumerate() {
            let Some(before) = self.accepted[slot] else {
                continue;
            };
            poll(cancelled)?;
            let sum = sum_arrays(
                boundary[grid.boundary_range(side)]
                    .iter()
                    .map(|q| q.array()),
                cancelled,
            )?;
            let weight = if final_stage { 1.0 } else { 2.0 };
            let transfer = Totals::from_array(sum.map(|q| weight * q));
            if !transfer.array().into_iter().all(f64::is_finite) {
                return Err(Error::Invalid("unrepresentable finite reservoir transfer"));
            }
            // The complete native energy flux already contains pressure work.
            // The stationary support does no work; it stores no bulk momentum.
            let after =
                before.exchange_candidate(gas, transfer.mass, transfer.energy, &mut *cancelled)?;
            let report = FiniteReservoirReport {
                before,
                after,
                transfer,
                mass_residual_kg: after.state().density * after.volume_m3()
                    - before.state().density * before.volume_m3()
                    - transfer.mass,
                energy_residual_j: after.state().energy * after.volume_m3()
                    - before.state().energy * before.volume_m3()
                    - transfer.energy,
            };
            if !report.mass_residual_kg.is_finite() || !report.energy_residual_j.is_finite() {
                return Err(Error::Invalid("unrepresentable finite reservoir balance"));
            }
            self.trial[slot] = Some(after);
            self.trial_reports[slot] = Some(report);
        }
        Ok(())
    }

    /// Infallible publication after the field's last cancellation poll.
    pub(super) fn publish(&mut self) {
        self.accepted = self.trial;
        self.reports = self.trial_reports;
    }
}

impl Euler2d {
    /// Bind independent finite tanks on selected physical transmissive sides.
    ///
    /// The original gas source recovers each tank's pressure and temperature
    /// from its evolving inventory at both native SSP RK2 stages. Stagnation
    /// inflow has zero shear; outflow carries the interior's complete kinetic
    /// energy. Choking and reversal use that same thermally perfect source,
    /// without a held-pressure approximation, fixed gamma or a second EOS.
    ///
    /// Bind once before evolution. A tank is rigid, adiabatic and perfectly
    /// mixed; different sides have independent storage. Internal tank waves,
    /// neck inertia, moving volume and thermal contact are not resolved.
    /// Retained hydrostatic references refuse until tank elevation and port
    /// potential transport are declared. Explicit body forcing remains valid.
    /// Molecular splitting retains zero molecular flux at these open faces.
    ///
    /// All stage storage is fixed-size. Accepted field, tank, clock and reports
    /// publish together; source refusal or cancellation preserves all of them.
    pub fn with_finite_reservoirs(
        mut self,
        inputs: &[FiniteReservoirBoundary],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if inputs.len() > 4 || self.reservoirs.any() || self.time != 0.0 || self.last_step.is_some()
        {
            return Err(Error::Invalid(
                "finite reservoirs bind at most four tanks once before evolution",
            ));
        }
        if !inputs.is_empty() && self.hydrostatic.is_some() {
            return Err(Error::Invalid(
                "finite tanks need declared hydrostatic port potentials",
            ));
        }
        let mut candidate = FiniteReservoirs::default();
        for input in inputs {
            poll(&mut cancelled)?;
            let axis = input.side.axis();
            let sign = input.side.sign();
            let slot = 2 * axis + usize::from(sign > 0.0);
            let AxisBoundaries::Ends { lower, upper } = self.boundaries.axis(axis) else {
                return Err(Error::Invalid(
                    "finite reservoir needs a physical transmissive side",
                ));
            };
            if [lower, upper][usize::from(sign > 0.0)] != BoundaryCondition::Transmissive
                || candidate.accepted[slot].is_some()
            {
                return Err(Error::Invalid(
                    "finite reservoir sides must be unique and transmissive",
                ));
            }
            let tank = input
                .tank
                .exchange_candidate(&self.gas, 0.0, 0.0, &mut cancelled)?;
            candidate.accepted[slot] = Some(tank);
            candidate.trial[slot] = Some(tank);
        }
        poll(&mut cancelled)?;
        self.reservoirs = candidate;
        Ok(self)
    }

    /// Accepted source-owned compact inventory; private retry stages are hidden.
    pub fn finite_reservoir(&self, side: BoundarySide) -> Option<&Plenum> {
        self.reservoirs.accepted[2 * side.axis() + usize::from(side.sign() > 0.0)].as_ref()
    }

    /// Last accepted transfer to the tank and its stationary support.
    pub fn finite_reservoir_report(&self, side: BoundarySide) -> Option<&FiniteReservoirReport> {
        self.reservoirs.reports[2 * side.axis() + usize::from(side.sign() > 0.0)].as_ref()
    }
}

#[cfg(test)]
mod tests;
