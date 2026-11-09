//! Planar bindings to the shared finite-mixture inventory transaction.
use super::*;
use fs_gas::multicomponent::network::Plenum;

/// Couple all physical faces on one bound side to one finite, well-mixed tank.
/// Different sides are independent stores, even when their initial states match.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoir {
    pub side: BoundarySide,
    /// Rigid tank volume [m³], strictly positive.
    pub volume_m3: f64,
}

impl Reservoirs {
    fn bind_finite(
        &mut self,
        gas: &MixtureGas,
        inputs: &[FiniteReservoir],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        cartesian::poll::<Error>(cancelled)?;
        if self.finite.has_any() || inputs.len() > 4 {
            return Err(Error::Invalid("finite mixture reservoirs bind at most four tanks once"));
        }
        let n = gas.species_count();
        let mut candidate = FiniteReservoirs::empty();
        let mut compact = self.stagnation;
        let mut flux = [0.0; WIDTH_LIMIT];
        for input in inputs {
            cartesian::poll::<Error>(cancelled)?;
            let axis = input.side.axis();
            let sign = input.side.sign();
            let slot = side_index(axis, sign);
            let exterior = self.boundary(input.side, n).ok_or(Error::Invalid(
                "finite mixture tank needs an already bound uniform reservoir",
            ))?;
            if exterior.velocity != [0.0; 2] {
                return Err(Error::Invalid("finite mixture tank must be exactly stationary"));
            }
            let state = fs_gas::multicomponent::Cell {
                partial_densities: exterior.partial_densities,
                momentum: 0.0,
                energy: self.reservoir_energy[slot],
            };
            candidate.bind(gas, slot, input.volume_m3, state, cancelled)?;
            // The same source compact-port path must admit the initial state.
            gas.reservoir_normal_flux_from_primitive(
                NormalPrimitive {
                    partial_densities: exterior.partial_densities,
                    normal_velocity: 0.0,
                    transverse_velocity: 0.0,
                    pressure: exterior.pressure_pa,
                },
                state,
                sign,
                &mut flux[..n + 3],
                &mut *cancelled,
            )?;
            compact[slot] = true;
        }
        cartesian::poll::<Error>(cancelled)?;
        self.finite = candidate;
        self.stagnation = compact;
        Ok(())
    }

    pub(in crate::compressible::mixture) fn copy_finite_from(
        &mut self,
        gas: &MixtureGas,
        source: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.finite.copy_accepted_from(gas, &source.finite, cancelled)
    }
}

impl MixtureEuler2d {
    /// Prepared total field capacity after binding this many finite tanks.
    ///
    /// Includes accepted and stage species inventories, fixed reports and all
    /// existing field/reservoir storage. Source-model storage is unchanged.
    pub fn required_workspace_bytes_with_finite_reservoirs(
        &self,
        count: usize,
    ) -> Result<usize, Error> {
        self.workspace_bytes
            .checked_sub(self.reservoirs.finite.workspace_bytes())
            .and_then(|bytes| bytes.checked_add(
                FiniteReservoirs::<2, 4>::required_bytes(&self.gas, count).ok()?,
            ))
            .ok_or(Error::Invalid("finite mixture reservoir workspace overflow"))
    }

    /// Evolve selected stationary reservoir compositions, pressure and energy
    /// as finite rigid tanks on this field's native SSP RK2 transaction.
    ///
    /// Each tank receives the exact native outward species and total-energy
    /// impulses. The full first-stage inventory drives the second residual;
    /// the final inventory receives the weighted accepted transfer. Incoming
    /// kinetic energy thermalizes, while the complete Cartesian impulse is
    /// reported to a stationary support. No extra enthalpy or p*dV is added.
    ///
    /// The source owns all calorics, species references, choking and reversal.
    /// Timestep admission includes each tank's summed aperture-wave rate, so
    /// a tiny tank cannot inherit an unsafe timestep from a large field cell.
    /// Warm stages, retries and publication allocate no storage.
    ///
    /// Bind once before the first step to existing uniform, zero-velocity
    /// reservoir sides. The supplied budget re-admits TOTAL numerical storage.
    /// Tanks do not resolve heat transfer, chemistry, internal waves or mixing
    /// times. A compact state does not provide a molecular boundary bath.
    /// Field cells, tanks, accepted exterior getters and reports commit together.
    pub fn with_finite_reservoirs(
        mut self,
        inputs: &[FiniteReservoir],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        cartesian::poll::<Error>(&mut cancelled)?;
        if self.time != 0.0 || self.last_step.is_some() {
            return Err(Error::Invalid("finite mixture reservoirs bind before the first step"));
        }
        let required = self.required_workspace_bytes_with_finite_reservoirs(inputs.len())?;
        if required > maximum_workspace_bytes || required > MAX_MIXTURE_WORKSPACE_BYTES {
            return Err(Error::Invalid("finite mixture reservoir total workspace budget"));
        }
        let before = self.reservoirs.finite.workspace_bytes();
        self.reservoirs.bind_finite(&self.gas, inputs, &mut cancelled)?;
        self.workspace_bytes = self.workspace_bytes
            .checked_sub(before)
            .and_then(|bytes| bytes.checked_add(self.reservoirs.finite.workspace_bytes()))
            .ok_or(Error::Invalid("finite mixture reservoir workspace overflow"))?;
        if self.workspace_bytes > maximum_workspace_bytes
            || self.workspace_bytes > MAX_MIXTURE_WORKSPACE_BYTES
        {
            return Err(Error::Invalid("finite mixture reservoir allocated workspace budget"));
        }
        cartesian::poll::<Error>(&mut cancelled)?;
        Ok(self)
    }

    /// Accepted species inventory and referenced internal energy for this tank.
    pub fn finite_reservoir(&self, side: BoundarySide) -> Option<&Plenum> {
        self.reservoirs.finite.inventory(side_index(side.axis(), side.sign()))
    }

    /// Accepted exchange, conservation residuals and complete support impulse.
    pub fn finite_reservoir_report(&self, side: BoundarySide) -> Option<&FiniteReservoirReport> {
        self.reservoirs.finite.report(side_index(side.axis(), side.sign()))
    }
}

#[cfg(test)]
mod tests;
