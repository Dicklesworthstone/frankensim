//! Finite rigid mixture tanks sharing the native volumetric SSP stages.

use super::*;
use fs_gas::multicomponent::network::Plenum;

/// Give one already-bound stationary uniform reservoir a finite rigid volume.
///
/// All faces on this side exchange with one perfectly mixed tank. Each side
/// names an independent tank, using the field's original canonical species,
/// calorics and energy references. The fixed support receives momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoir {
    pub side: BoundarySide,
    /// Positive physical tank volume [m³].
    pub volume_m3: f64,
}

/// Accepted species, pressure, energy and support exchange of one 3-D tank.
pub type FiniteReservoirReport = crate::compressible::mixture::finite::FiniteReservoirReport<3>;

impl MixtureEuler3d {
    /// Total numerical capacity after binding this many finite mixture tanks.
    ///
    /// Includes accepted and reusable stage inventories, without changing the
    /// field's already admitted capacities. Existing finite tanks cannot be
    /// rebound through this constructor.
    pub fn required_workspace_bytes_with_finite_reservoirs(
        &self,
        count: usize,
    ) -> Result<usize, Error> {
        if self.reservoirs.finite.has_any() {
            return Err(Error::Invalid("finite mixture tanks are already bound"));
        }
        self.workspace_bytes
            .checked_sub(self.reservoirs.finite.workspace_bytes())
            .and_then(|bytes| {
                FiniteReservoirs::<3, 6>::required_bytes(&self.gas, count)
                    .ok()
                    .and_then(|finite| bytes.checked_add(finite))
            })
            .ok_or(Error::Invalid("finite volumetric reservoir workspace size"))
    }

    /// Bind finite, adiabatic, well-mixed tanks to stationary uniform ports.
    ///
    /// The original compact pressure law receives the tank's current pressure
    /// and composition at both SSPRK2 stages. The same outward species and full
    /// total-energy face impulses enter the tank. Inflow kinetic energy becomes
    /// tank internal energy; the fixed support absorbs momentum without work.
    /// No additional enthalpy or pressure-work term is added.
    ///
    /// The explicit budget bounds the complete prepared field. No numerical
    /// storage is allocated while stepping. Tank aperture/volume wave rates
    /// join the Cartesian CFL bound, and source admission may shorten a trial
    /// before any tank, field state, clock, load or report becomes visible.
    ///
    /// These tanks have no separate reaction, heat-transfer or molecular-bath
    /// model. Held profiles on other sides and native field forcing can coexist.
    /// Each tank must start on an already-bound uniform side with exactly zero
    /// velocity; its original state is retained with the supplied rigid volume.
    pub fn with_finite_reservoirs(
        mut self,
        inputs: &[FiniteReservoir],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        let required = self.required_workspace_bytes_with_finite_reservoirs(inputs.len())?;
        if required > maximum_workspace_bytes || required > MAX_MIXTURE_WORKSPACE_BYTES {
            return Err(Error::Invalid(
                "finite volumetric reservoir workspace budget",
            ));
        }
        let mut candidate = FiniteReservoirs::<3, 6>::empty();
        let n = self.gas.species_count();
        for input in inputs {
            stepping::poll::<Error>(&mut cancelled)?;
            let axis = input.side.axis();
            let sign = input.side.sign();
            let slot = side_index(axis, sign);
            let index = self.reservoirs.slots[slot].ok_or(Error::Invalid(
                "finite tank requires an already-bound uniform volumetric reservoir",
            ))?;
            if self.reservoirs.primitive[index].velocity != [0.0; 3] {
                return Err(Error::Invalid("finite mixture tank must be stationary"));
            }
            let row = self
                .reservoirs
                .row(axis, sign, n)
                .expect("uniform reservoir slot");
            candidate.bind(
                &self.gas,
                slot,
                input.volume_m3,
                fs_gas::multicomponent::Cell {
                    partial_densities: &row[..n],
                    momentum: 0.0,
                    energy: self.reservoirs.reservoir_energy[slot],
                },
                &mut cancelled,
            )?;
        }
        let bytes = self
            .workspace_bytes
            .checked_sub(self.reservoirs.finite.workspace_bytes())
            .and_then(|bytes| bytes.checked_add(candidate.workspace_bytes()))
            .ok_or(Error::Invalid("finite volumetric reservoir workspace size"))?;
        if bytes > maximum_workspace_bytes || bytes > MAX_MIXTURE_WORKSPACE_BYTES {
            return Err(Error::Invalid(
                "finite volumetric reservoir workspace budget",
            ));
        }
        stepping::poll::<Error>(&mut cancelled)?;
        self.reservoirs.finite = candidate;
        for input in inputs {
            self.reservoirs.stagnation[side_index(input.side.axis(), input.side.sign())] = true;
        }
        self.workspace_bytes = bytes;
        Ok(self)
    }

    /// Accepted source-owned tank inventory, never a rejected private stage.
    pub fn finite_reservoir(&self, side: BoundarySide) -> Option<&Plenum> {
        self.reservoirs
            .finite
            .inventory(side_index(side.axis(), side.sign()))
    }

    /// Tank exchange over the latest accepted native field step.
    pub fn finite_reservoir_report(&self, side: BoundarySide) -> Option<&FiniteReservoirReport> {
        self.reservoirs
            .finite
            .report(side_index(side.axis(), side.sign()))
    }
}
