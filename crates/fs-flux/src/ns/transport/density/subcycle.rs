//! Conservative density subcycling for one frozen native volume advector.
//!
//! Equal dyadic subintervals reuse the original positive donor-cell Euler map.
//! The arithmetic mean of their PRE-step donor densities defines the integrated
//! mass flux because every original directional volume rate remains fixed.
//! Reconstructing that mean once retains the original constant cell divergence
//! and quadratic face moments needed by a native macrostep momentum balance.
//! The averaged density is an integration weight, not an accepted material state.
//!
//! This changes the density CFL cadence, not the order or physical accuracy of
//! an explicitly frozen advector. It does not move interfaces sharply, change
//! geometry, add a mass source, or supply a pressure/constitutive closure.

use super::super::mass_flux::PreparedDensityMassFlux;
use super::{
    DensityTransportCell, DensityTransportError, DensityTransportReport, DensityTransportWorkspace,
    EdgeDensityTransfer, PreparedDensityTransport, Result, Sum, admit_cell, admit_courant,
    admit_storage, allocate, finite, poll, positive, positive_value,
};

/// Physical interval and work admission for the original density Euler map.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensitySubcycleControls {
    /// Largest outgoing volume fraction in EACH density substep, in (0,1).
    pub maximum_courant: f64,
    /// Maximum admitted number of equal dyadic substeps. No global momentum
    /// solves occur in this owner, and no unbounded retry loop is hidden here.
    pub maximum_substeps: usize,
}

/// One accepted macro interval, retaining physical endpoint/boundary accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensitySubcycleReport<const N: usize> {
    /// The original inventory report over the COMPLETE macro interval.
    /// Its `maximum_courant` is the macro outgoing volume fraction and can
    /// exceed one. Positivity uses `maximum_substep_courant` below instead.
    pub transport: DensityTransportReport<N>,
    /// Actual power-of-two number of accepted density Euler steps.
    pub substeps: usize,
    /// Their equal physical duration [s]. Repeated binary splitting preserves
    /// the macro duration exactly; unrepresentable subdivisions are refused.
    pub substep_timestep_s: f64,
    /// Actual maximum outgoing volume fraction of an accepted substep.
    pub maximum_substep_courant: f64,
}

/// Reusable subcycle storage bound to ONE immutable mesh/volume-flow plan.
/// Warm advances allocate nothing; construction admits all retained elements,
/// including the original Euler workspace. A refused attempt leaves caller
/// output unchanged and invalidates the average/face-transfer views.
pub struct DensitySubcycleWorkspace<'p, 'm, const N: usize> {
    plan: &'p PreparedDensityTransport<'m>,
    euler: DensityTransportWorkspace<N>,
    current: Vec<DensityTransportCell<N>>,
    next: Vec<DensityTransportCell<N>>,
    average: Vec<Mean>,
    average_cells: Vec<DensityTransportCell<0>>,
    average_inflow: Vec<Option<DensityTransportCell<0>>>,
    transfers: Vec<TransferSum<N>>,
    accepted: bool,
}

#[derive(Clone, Copy, Default)]
pub(super) struct Mean {
    scale: f64,
    sum: Sum,
}
impl Mean {
    pub(super) fn add(&mut self, density: f64) -> Result<()> {
        if self.scale == 0.0 {
            self.scale = density;
            self.sum.add(1.0)?;
        } else if density > self.scale {
            let ratio = self.scale / density;
            // Keep both compensated parts while changing the positive scale.
            // Neither sum(rho) overflow nor sum(rho/n) underflow is required.
            self.sum.value *= ratio;
            self.sum.correction *= ratio;
            self.scale = density;
            self.sum.add(1.0)?;
        } else {
            self.sum.add(density / self.scale)?;
        }
        Ok(())
    }

    pub(super) fn value(self, count: usize) -> Result<f64> {
        positive_value(
            self.scale * (self.sum.value()? / count as f64),
            "time-averaged donor density",
        )
    }
}

#[derive(Clone, Copy)]
pub(super) struct TransferSum<const N: usize> {
    mass: [Sum; 2],
    carried: [[Sum; N]; 2],
}
impl<const N: usize> Default for TransferSum<N> {
    fn default() -> Self {
        Self {
            mass: [Sum::default(); 2],
            carried: [[Sum::default(); N]; 2],
        }
    }
}
impl<const N: usize> TransferSum<N> {
    pub(super) fn add(&mut self, transfer: &EdgeDensityTransfer<N>) -> Result<()> {
        for direction in 0..2 {
            self.mass[direction].add(transfer.mass_kg[direction])?;
            for component in 0..N {
                self.carried[direction][component].add(transfer.carried[direction][component])?;
            }
        }
        Ok(())
    }

    pub(super) fn value(self) -> Result<EdgeDensityTransfer<N>> {
        let mut value = EdgeDensityTransfer::default();
        for direction in 0..2 {
            value.mass_kg[direction] = self.mass[direction].value()?;
            for component in 0..N {
                value.carried[direction][component] = self.carried[direction][component].value()?;
            }
        }
        Ok(value)
    }
}

impl<'m> PreparedDensityTransport<'m> {
    /// Allocate reusable subcycle work for this frozen original plan. The byte
    /// cap covers BOTH stages, the original Euler workspace, compensated donor
    /// averages and cumulative face amounts. It excludes caller buffers, the
    /// plan itself and a subsequently prepared mass-flux operator.
    pub fn subcycle_workspace<const N: usize>(
        &self,
        max_storage_bytes: usize,
    ) -> Result<DensitySubcycleWorkspace<'_, 'm, N>> {
        let nc = self.volumes.len();
        let ne = self.volume_rates.len();
        admit_storage(
            &[
                // Keep generic-N sizes separate so every addition and
                // multiplication goes through the checked storage owner.
                (nc, size_of::<DensityTransportCell<N>>()),
                (nc, size_of::<DensityTransportCell<N>>()),
                (nc, size_of::<DensityTransportCell<N>>()),
                (nc, size_of::<f64>() + size_of::<Sum>()),
                (nc, size_of::<[Sum; N]>()),
                (nc, size_of::<Mean>() + size_of::<DensityTransportCell<0>>()),
                (ne, size_of::<EdgeDensityTransfer<N>>()),
                (ne, size_of::<TransferSum<N>>()),
                (ne, size_of::<Option<DensityTransportCell<0>>>()),
            ],
            max_storage_bytes,
        )?;
        Ok(DensitySubcycleWorkspace {
            plan: self,
            euler: self.workspace(max_storage_bytes)?,
            current: allocate(nc, DensityTransportCell::default())?,
            next: allocate(nc, DensityTransportCell::default())?,
            average: allocate(nc, Mean::default())?,
            average_cells: allocate(nc, DensityTransportCell::default())?,
            average_inflow: allocate(ne, None)?,
            transfers: allocate(ne, TransferSum::default())?,
            accepted: false,
        })
    }
}

impl<'m, const N: usize> DensitySubcycleWorkspace<'_, 'm, N> {
    /// Cumulative original directional mass/carried transfers over the entire
    /// latest accepted macro interval. A failed later attempt invalidates it.
    pub fn face_transfers(&self) -> Option<&[EdgeDensityTransfer<N>]> {
        self.accepted.then_some(self.euler.transfers.as_slice())
    }

    /// Equal-time average of the PRE-substep donor densities. The empty
    /// specific array makes explicit that these are mass-flux integration
    /// weights, not averaged enthalpy/species or accepted material states.
    pub fn average_donor_density(&self) -> Option<&[DensityTransportCell<0>]> {
        self.accepted.then_some(self.average_cells.as_slice())
    }

    /// Reconstruct the one time-averaged compatible mass flux of the accepted
    /// substeps using their original frozen advector. Reconstruction is linear
    /// in donor densities in exact arithmetic, so its divergence gives the
    /// macro density change and its quadratic face moments give the native
    /// momentum energy identity, up to conversion/averaging roundoff.
    /// Boundary donors are held constant over this macro interval.
    pub fn prepare_mass_flux(
        &self,
        max_storage_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDensityMassFlux<'m>> {
        if !self.accepted {
            return Err(DensityTransportError::Invalid(
                "an accepted density subcycle is required before reconstructing its mass flux",
            ));
        }
        self.plan.prepare_mass_flux(
            &self.average_cells,
            &self.average_inflow,
            max_storage_bytes,
            cancelled,
        )
    }

    /// Advance one macro interval with the original positive donor-cell Euler
    /// map. Density and ALL N extensive carried inventories use the actual
    /// microstep transfers. The advector, fixed volumes and exterior donor
    /// states stay frozen. No global momentum/pressure solve is done here.
    ///
    /// Substeps are equal dyadic subdivisions, at most twice the least integer
    /// CFL count in exact arithmetic. Their duration must split exactly in f64.
    /// The explicit work cap is admitted before touching any physical state.
    /// Preparation and every cell/edge pass poll cancellation. Failure leaves
    /// `output` unchanged; workspace remains reusable with invalidated views.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_into(
        &mut self,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensitySubcycleControls,
        output: &mut [DensityTransportCell<N>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensitySubcycleReport<N>> {
        self.accepted = false;
        self.euler.accepted = false;
        poll(&mut cancelled)?;
        if cells.len() != self.current.len()
            || output.len() != self.current.len()
            || inflow.len() != self.transfers.len()
            || !positive(timestep_s)
            || controls.maximum_substeps == 0
        {
            return Err(DensityTransportError::Invalid(
                "matching density subcycle dimensions, positive interval and nonzero work budget required",
            ));
        }
        admit_courant(controls.maximum_courant)?;
        let mut count = 1usize;
        let mut dt = timestep_s;
        loop {
            let courant = dt * self.plan.maximum_outgoing_rate_s;
            if courant.is_finite() && courant <= controls.maximum_courant {
                break;
            }
            count = count
                .checked_mul(2)
                .filter(|&n| n <= controls.maximum_substeps)
                .ok_or(DensityTransportError::SubstepLimit {
                    maximum: controls.maximum_substeps,
                })?;
            let half = positive_value(0.5 * dt, "density substep duration")?;
            if half + half != dt {
                return Err(DensityTransportError::Arithmetic(
                    "density macro interval does not split into representable equal substeps",
                ));
            }
            dt = half;
        }
        let macro_courant = finite(
            timestep_s * self.plan.maximum_outgoing_rate_s,
            "macro outgoing volume fraction",
        )?;
        for (index, cell) in cells.iter().enumerate() {
            poll(&mut cancelled)?;
            admit_cell(cell)?;
            self.current[index] = *cell;
            self.average[index] = Mean::default();
        }
        for (edge, state) in inflow.iter().enumerate() {
            poll(&mut cancelled)?;
            if self.plan.mesh.edges[edge].tris.1 != usize::MAX && state.is_some() {
                return Err(DensityTransportError::Invalid(
                    "inflow states belong only to boundary edges",
                ));
            }
            if let Some(state) = state {
                admit_cell(state)?;
            }
            self.average_inflow[edge] = state.map(|state| DensityTransportCell {
                density_kg_m3: state.density_kg_m3,
                specific: [],
            });
            self.transfers[edge] = TransferSum::default();
        }
        let mut maximum_substep_courant: f64 = 0.0;
        let mut absolute_boundary = [Sum::default(); N];
        for _ in 0..count {
            poll(&mut cancelled)?;
            let step = self.plan.advance_into(
                &self.current,
                inflow,
                dt,
                controls.maximum_courant,
                &mut self.next,
                &mut self.euler,
                &mut cancelled,
            )?;
            maximum_substep_courant = maximum_substep_courant.max(step.maximum_courant);
            for (mean, state) in self.average.iter_mut().zip(&self.current) {
                poll(&mut cancelled)?;
                mean.add(state.density_kg_m3)?;
            }
            for (edge, (sum, transfer)) in self
                .transfers
                .iter_mut()
                .zip(&self.euler.transfers)
                .enumerate()
            {
                poll(&mut cancelled)?;
                sum.add(transfer)?;
                if self.plan.mesh.edges[edge].tris.1 == usize::MAX {
                    for c in 0..N {
                        // Preserve the absolute scale even if a signed
                        // specific quantity reverses during this interval.
                        absolute_boundary[c].add(transfer.carried[0][c].abs())?;
                        absolute_boundary[c].add(transfer.carried[1][c].abs())?;
                    }
                }
            }
            core::mem::swap(&mut self.current, &mut self.next);
        }
        for (target, mean) in self.average_cells.iter_mut().zip(&self.average) {
            poll(&mut cancelled)?;
            target.density_kg_m3 = mean.value(count)?;
        }
        // Reuse the original Euler transfer storage for the accepted complete
        // interval. Its last microstep view is no longer needed.
        for (target, sum) in self.euler.transfers.iter_mut().zip(&self.transfers) {
            poll(&mut cancelled)?;
            *target = sum.value()?;
        }
        let transport = audit(
            self.plan,
            cells,
            &self.current,
            &self.euler.transfers,
            &absolute_boundary,
            timestep_s,
            macro_courant,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        output.copy_from_slice(&self.current);
        self.accepted = true;
        Ok(DensitySubcycleReport {
            transport,
            substeps: count,
            substep_timestep_s: dt,
            maximum_substep_courant,
        })
    }
}

// Re-audit the represented endpoints and actual cumulative face amounts once
// over the whole interval. Summing relative microstep residuals or subtracting
// two large total inventories would hide the physical macro accounting scale.
#[allow(clippy::too_many_arguments)]
pub(super) fn audit<const N: usize>(
    plan: &PreparedDensityTransport<'_>,
    old: &[DensityTransportCell<N>],
    new: &[DensityTransportCell<N>],
    transfers: &[EdgeDensityTransfer<N>],
    absolute_boundary: &[Sum; N],
    timestep_s: f64,
    maximum_courant: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DensityTransportReport<N>> {
    let mut mass_before = Sum::default();
    let mut mass_after = Sum::default();
    let mut mass_in = Sum::default();
    let mut mass_out = Sum::default();
    let mut mass_residual = Sum::default();
    let mut carried_before = [Sum::default(); N];
    let mut carried_after = [Sum::default(); N];
    let mut carried_in = [Sum::default(); N];
    let mut carried_out = [Sum::default(); N];
    let mut carried_residual = [Sum::default(); N];
    let mut absolute_before = [Sum::default(); N];
    let mut absolute_after = [Sum::default(); N];
    for cell in 0..old.len() {
        poll(cancelled)?;
        let before = positive_value(
            old[cell].density_kg_m3 * plan.volumes[cell],
            "old macro mass",
        )?;
        let after = positive_value(
            new[cell].density_kg_m3 * plan.volumes[cell],
            "new macro mass",
        )?;
        mass_before.add(before)?;
        mass_after.add(after)?;
        mass_residual.add(after - before)?;
        for c in 0..N {
            let a = finite(
                before * old[cell].specific[c],
                "old macro carried inventory",
            )?;
            let b = finite(after * new[cell].specific[c], "new macro carried inventory")?;
            carried_before[c].add(a)?;
            carried_after[c].add(b)?;
            absolute_before[c].add(a.abs())?;
            absolute_after[c].add(b.abs())?;
            carried_residual[c].add(b)?;
            carried_residual[c].add(-a)?;
            carried_residual[c].add(after.mul_add(new[cell].specific[c], -b))?;
            carried_residual[c].add(-before.mul_add(old[cell].specific[c], -a))?;
        }
    }
    for (edge, transfer) in plan.mesh.edges.iter().zip(transfers) {
        poll(cancelled)?;
        if edge.tris.1 != usize::MAX {
            continue;
        }
        mass_in.add(transfer.mass_kg[1])?;
        mass_out.add(transfer.mass_kg[0])?;
        mass_residual.add(-transfer.mass_kg[1])?;
        mass_residual.add(transfer.mass_kg[0])?;
        for c in 0..N {
            let incoming = transfer.carried[1][c];
            let outgoing = transfer.carried[0][c];
            carried_in[c].add(incoming)?;
            carried_out[c].add(outgoing)?;
            carried_residual[c].add(-incoming)?;
            carried_residual[c].add(outgoing)?;
        }
    }
    let before = mass_before.value()?;
    let after = mass_after.value()?;
    let incoming = mass_in.value()?;
    let outgoing = mass_out.value()?;
    let residual = mass_residual.value()?;
    let mut report = DensityTransportReport {
        timestep_s,
        maximum_courant,
        mass_before_kg: before,
        mass_after_kg: after,
        inlet_mass_kg: incoming,
        outlet_mass_kg: outgoing,
        mass_balance_residual_kg: residual,
        relative_mass_balance_residual: residual / before.max(after).max(incoming).max(outgoing),
        carried_before: [0.0; N],
        carried_after: [0.0; N],
        carried_in: [0.0; N],
        carried_out: [0.0; N],
        carried_balance_residual: [0.0; N],
        relative_carried_balance_residual: [0.0; N],
    };
    for c in 0..N {
        report.carried_before[c] = carried_before[c].value()?;
        report.carried_after[c] = carried_after[c].value()?;
        report.carried_in[c] = carried_in[c].value()?;
        report.carried_out[c] = carried_out[c].value()?;
        let residual = carried_residual[c].value()?;
        report.carried_balance_residual[c] = residual;
        let scale = absolute_before[c]
            .value()?
            .max(absolute_after[c].value()?)
            .max(absolute_boundary[c].value()?);
        report.relative_carried_balance_residual[c] =
            if scale == 0.0 { 0.0 } else { residual / scale };
    }
    Ok(report)
}
