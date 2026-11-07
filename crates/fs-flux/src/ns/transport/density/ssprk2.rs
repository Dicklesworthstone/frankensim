//! Conservative SSPRK2 density stages under one fixed original BDM advector.
//!
//! Each physical subinterval uses TWO positive conservative Euler maps,
//! then averages the initial and second-stage EXTENSIVE states. The accepted
//! face amounts are the half-sum of the two stage transfers; the compatible
//! mass flux uses the equal-time mean of their PRE-stage directional donors.
//! All stages use the same volume advector, geometry and exterior donors.
//!
//! This is second order in TIME for the fixed-advector spatially discrete
//! transport problem. Density uses the configured donor-cell or bounded linear
//! spatial reconstruction; specific inventories retain donor-cell mixing.
//! A native BDF1 momentum macrostep remains first
//! order; neither interface sharpening nor a nonlinear coupling is implied.

use super::super::mass_flux::PreparedDensityMassFlux;
use super::subcycle::{Mean, TransferSum, add_face_means, audit, finish_face_means};
use super::{
    DensityTransportCell, DensityTransportError, DensityTransportReport, DensityTransportWorkspace,
    EdgeDensityTransfer, PreparedDensityTransport, Result, Sum, admit_cell, admit_courant,
    admit_storage, allocate, finite, poll, positive, positive_value,
};

/// Explicit CFL and work admission for the SSPRK2 transport stages.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensitySspRk2Controls {
    /// Maximum outgoing-volume fraction of EACH original Euler stage, in (0,1).
    /// Bounded linear density uses one third of this sufficient volume allowance.
    pub maximum_courant: f64,
    /// Maximum equal dyadic physical subintervals. Each executes exactly TWO
    /// Euler evaluations; the accepted report records both work counts.
    pub maximum_substeps: usize,
}

/// Complete physical macro accounting and actual SSPRK2 work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensitySspRk2Report<const N: usize> {
    /// Original inventory report over the complete accepted physical interval.
    /// `maximum_courant` is the macro volume fraction and may exceed one.
    pub transport: DensityTransportReport<N>,
    /// Power-of-two number of accepted physical SSPRK2 subintervals.
    pub substeps: usize,
    /// Exactly twice `substeps`; stages are not extra physical elapsed time.
    pub stage_evaluations: usize,
    /// Equal physical duration of each SSPRK2 subinterval [s].
    pub substep_timestep_s: f64,
    /// Largest outgoing-volume fraction of either Euler stage.
    pub maximum_stage_courant: f64,
}

/// Reusable stage/accounting storage bound to one immutable original plan.
/// Warm advances allocate nothing. A failed attempt invalidates accepted
/// views and leaves caller output unchanged; the workspace can be retried.
pub struct DensitySspRk2Workspace<'p, 'm, const N: usize> {
    plan: &'p PreparedDensityTransport<'m>,
    euler: DensityTransportWorkspace<N>,
    current: Vec<DensityTransportCell<N>>,
    stage: Vec<DensityTransportCell<N>>,
    second: Vec<DensityTransportCell<N>>,
    first_transfers: Vec<EdgeDensityTransfer<N>>,
    average: Vec<Mean>,
    average_cells: Vec<DensityTransportCell<0>>,
    average_inflow: Vec<Option<DensityTransportCell<0>>>,
    face_average: Vec<[Mean; 2]>,
    transfers: Vec<TransferSum<N>>,
    accepted: bool,
}

impl<'m> PreparedDensityTransport<'m> {
    /// Allocate all retained SSPRK2 storage under one byte cap: original Euler
    /// scratch, three physical/stage buffers, first-stage face amounts,
    /// compensated donor means and cumulative accepted face amounts. Caller
    /// buffers, this plan and a later mass-flux reconstruction have separate
    /// ownership and allowances.
    pub fn ssprk2_workspace<const N: usize>(
        &self,
        max_storage_bytes: usize,
    ) -> Result<DensitySspRk2Workspace<'_, 'm, N>> {
        let nc = self.volumes.len();
        let ne = self.volume_rates.len();
        let linear_cells = if self.linear.is_some() { nc } else { 0 };
        let linear_edges = if self.linear.is_some() { ne } else { 0 };
        let inventory_cells = if self.reconstruct_inventories && N > 0 {
            nc
        } else {
            0
        };
        admit_storage(
            &[
                (nc, size_of::<DensityTransportCell<N>>()),
                (nc, size_of::<DensityTransportCell<N>>()),
                (nc, size_of::<DensityTransportCell<N>>()),
                (nc, size_of::<DensityTransportCell<N>>()),
                (nc, size_of::<f64>() + size_of::<Sum>()),
                (nc, size_of::<[Sum; N]>()),
                (nc, size_of::<Mean>() + size_of::<DensityTransportCell<0>>()),
                (ne, size_of::<EdgeDensityTransfer<N>>()),
                (ne, size_of::<EdgeDensityTransfer<N>>()),
                (ne, size_of::<TransferSum<N>>()),
                (ne, size_of::<Option<DensityTransportCell<0>>>()),
                (linear_cells, size_of::<[f64; 3]>()),
                (
                    if self.linear.is_some() {
                        self.mesh.verts.len()
                    } else {
                        0
                    },
                    size_of::<[f64; 2]>(),
                ),
                (linear_edges, size_of::<[f64; 2]>() + size_of::<[Mean; 2]>()),
                (inventory_cells, size_of::<[[f64; 3]; N]>()),
                (inventory_cells, size_of::<[[f64; 2]; N]>()),
                (inventory_cells, size_of::<f64>()),
            ],
            max_storage_bytes,
        )?;
        Ok(DensitySspRk2Workspace {
            plan: self,
            euler: self.workspace(max_storage_bytes)?,
            current: allocate(nc, DensityTransportCell::default())?,
            stage: allocate(nc, DensityTransportCell::default())?,
            second: allocate(nc, DensityTransportCell::default())?,
            first_transfers: allocate(ne, EdgeDensityTransfer::default())?,
            average: allocate(nc, Mean::default())?,
            average_cells: allocate(nc, DensityTransportCell::default())?,
            average_inflow: allocate(ne, None)?,
            face_average: allocate(linear_edges, [Mean::default(); 2])?,
            transfers: allocate(ne, TransferSum::default())?,
            accepted: false,
        })
    }
}

impl<'m, const N: usize> DensitySspRk2Workspace<'_, 'm, N> {
    /// RK-weighted directional mass and all carried amounts from the latest
    /// complete accepted physical interval. These are half-sums of actual
    /// stage transfers, not two full Euler intervals of boundary exchange.
    pub fn face_transfers(&self) -> Option<&[EdgeDensityTransfer<N>]> {
        self.accepted.then_some(self.euler.transfers.as_slice())
    }

    /// Equal-weight mean of the two PRE-Euler donor densities in every physical
    /// subinterval. These are flux integration weights, not material states.
    /// Bounded linear reconstruction separately averages the actual limited
    /// directional donors; reconstructing these averaged cells is not equivalent.
    pub fn average_donor_density(&self) -> Option<&[DensityTransportCell<0>]> {
        self.accepted.then_some(self.average_cells.as_slice())
    }

    /// Reconstruct one original compatible mass flux with the accepted RK
    /// weights. Under this frozen advector its integrated directional rates
    /// equal the accepted weighted stage amounts and its divergence gives the
    /// total density increment, up to reported representation/summation error.
    /// Original quadratic face moments and the native momentum owner are reused.
    pub fn prepare_mass_flux(
        &self,
        max_storage_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDensityMassFlux<'m>> {
        if !self.accepted {
            return Err(DensityTransportError::Invalid(
                "accepted SSPRK2 density stages are required before reconstructing their flux",
            ));
        }
        self.plan.prepare_mass_flux_with_donors(
            &self.average_cells,
            &self.average_inflow,
            self.plan
                .linear
                .as_ref()
                .map(|_| self.euler.face_densities.as_slice()),
            max_storage_bytes,
            cancelled,
        )
    }

    /// Advance one physical macro interval using bounded dyadic SSPRK2 steps.
    /// For the conservative state U=(mass, all extensive inventories),
    /// `U1=Euler_h(U0)`, `U2=Euler_h(U1)`, `U+=0.5*(U0+U2)`.
    /// Both stages obey the original outgoing-volume CFL condition, so their
    /// extensive convex combination preserves positive density and the original
    /// donor-specific bounds up to floating-point representation error.
    ///
    /// Each carried component is blended through physical mass TIMES specific
    /// value. Mass fractions are not formed first: a light donor can retain a
    /// finite inventory even when its fraction of total mass underflows. Equal
    /// specific values are copied exactly. No density or inventory is repaired.
    ///
    /// The advector, original volumes and boundary states remain frozen. All
    /// stage evaluations reuse the original Euler workspace. Work/storage are
    /// admitted explicitly; no output or accepted view is published on refusal
    /// or cancellation. No native pressure/momentum solve occurs here.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_into(
        &mut self,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensitySspRk2Controls,
        output: &mut [DensityTransportCell<N>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensitySspRk2Report<N>> {
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
                "matching SSPRK2 dimensions, positive interval and nonzero work budget required",
            ));
        }
        admit_courant(controls.maximum_courant)?;
        let mut count = 1usize;
        let mut dt = timestep_s;
        loop {
            let courant = dt * self.plan.maximum_outgoing_rate_s;
            if courant.is_finite()
                && courant <= self.plan.volume_courant_limit(controls.maximum_courant)
            {
                break;
            }
            count = count
                .checked_mul(2)
                .filter(|&n| n <= controls.maximum_substeps)
                .ok_or(DensityTransportError::SubstepLimit {
                    maximum: controls.maximum_substeps,
                })?;
            let half = positive_value(0.5 * dt, "SSPRK2 substep duration")?;
            if half + half != dt {
                return Err(DensityTransportError::Arithmetic(
                    "density macro interval does not split into representable equal SSPRK2 substeps",
                ));
            }
            dt = half;
        }
        let stage_evaluations = count
            .checked_mul(2)
            .ok_or(DensityTransportError::Arithmetic("SSPRK2 stage work count"))?;
        let macro_courant = finite(
            timestep_s * self.plan.maximum_outgoing_rate_s,
            "macro outgoing volume fraction",
        )?;
        for (index, state) in cells.iter().enumerate() {
            poll(&mut cancelled)?;
            admit_cell(state)?;
            self.current[index] = *state;
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
            if let Some(mean) = self.face_average.get_mut(edge) {
                *mean = [Mean::default(); 2];
            }
        }
        let mut maximum_stage_courant: f64 = 0.0;
        let mut absolute_boundary = [Sum::default(); N];
        for _ in 0..count {
            poll(&mut cancelled)?;
            let first = self.plan.advance_into(
                &self.current,
                inflow,
                dt,
                controls.maximum_courant,
                &mut self.stage,
                &mut self.euler,
                &mut cancelled,
            )?;
            maximum_stage_courant = maximum_stage_courant.max(first.maximum_courant);
            add_face_means(
                &mut self.face_average,
                &self.euler.face_densities,
                &mut cancelled,
            )?;
            for (saved, transfer) in self.first_transfers.iter_mut().zip(&self.euler.transfers) {
                poll(&mut cancelled)?;
                *saved = *transfer;
            }
            let second = self.plan.advance_into(
                &self.stage,
                inflow,
                dt,
                controls.maximum_courant,
                &mut self.second,
                &mut self.euler,
                &mut cancelled,
            )?;
            maximum_stage_courant = maximum_stage_courant.max(second.maximum_courant);
            add_face_means(
                &mut self.face_average,
                &self.euler.face_densities,
                &mut cancelled,
            )?;
            for cell in 0..self.current.len() {
                poll(&mut cancelled)?;
                self.average[cell].add(self.current[cell].density_kg_m3)?;
                self.average[cell].add(self.stage[cell].density_kg_m3)?;
                self.current[cell] = blend(
                    &self.current[cell],
                    &self.second[cell],
                    self.plan.volumes[cell],
                )?;
            }
            for edge in 0..self.transfers.len() {
                poll(&mut cancelled)?;
                let first = &self.first_transfers[edge];
                let second = &self.euler.transfers[edge];
                let mut accepted = EdgeDensityTransfer::default();
                for direction in 0..2 {
                    let a = first.mass_kg[direction];
                    let b = second.mass_kg[direction];
                    let mass = midpoint(a, b);
                    if mass == 0.0 && (a > 0.0 || b > 0.0) {
                        return Err(DensityTransportError::Arithmetic(
                            "SSPRK2 weighted face mass underflow",
                        ));
                    }
                    accepted.mass_kg[direction] = mass;
                    for component in 0..N {
                        let a = first.carried[direction][component];
                        let b = second.carried[direction][component];
                        accepted.carried[direction][component] = midpoint(a, b);
                        if self.plan.mesh.edges[edge].tris.1 == usize::MAX {
                            // Absolute RK-weighted exchange, even when signed
                            // carried flux changes sign between the two stages.
                            absolute_boundary[component].add(midpoint(a.abs(), b.abs()))?;
                        }
                    }
                }
                self.transfers[edge].add(&accepted)?;
            }
        }
        for (target, mean) in self.average_cells.iter_mut().zip(&self.average) {
            poll(&mut cancelled)?;
            target.density_kg_m3 = mean.value(stage_evaluations)?;
        }
        finish_face_means(
            &mut self.euler.face_densities,
            &self.face_average,
            stage_evaluations,
            &mut cancelled,
        )?;
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
        Ok(DensitySspRk2Report {
            transport,
            substeps: count,
            stage_evaluations,
            substep_timestep_s: dt,
            maximum_stage_courant,
        })
    }
}

// Preserve equal subnormal amounts and near-cancelling signed amounts before
// halving. Only halve first when the sum of equal-sign large values may overflow.
fn midpoint(a: f64, b: f64) -> f64 {
    if a.is_sign_positive() != b.is_sign_positive() || a.abs().max(b.abs()) <= 0.5 * f64::MAX {
        (a + b) * 0.5
    } else {
        a * 0.5 + b * 0.5
    }
}

fn blend<const N: usize>(
    initial: &DensityTransportCell<N>,
    second: &DensityTransportCell<N>,
    volume: f64,
) -> Result<DensityTransportCell<N>> {
    let a = positive_value(initial.density_kg_m3 * volume, "initial SSPRK2 blend mass")?;
    let b = positive_value(second.density_kg_m3 * volume, "second SSPRK2 blend mass")?;
    let mass = positive_value(midpoint(a, b), "SSPRK2 blended mass")?;
    let density = positive_value(mass / volume, "SSPRK2 blended density")?;
    let mut specific = initial.specific;
    for component in 0..N {
        let x = initial.specific[component];
        let y = second.specific[component];
        if x == y {
            continue;
        }
        let p = finite(a * x, "initial SSPRK2 carried inventory")?;
        let q = finite(b * y, "second SSPRK2 carried inventory")?;
        let dp = a.mul_add(x, -p);
        let dq = b.mul_add(y, -q);
        // Blend the actual extensive products BEFORE dividing by mass. A
        // small mass fraction can underflow while its carried product is finite.
        // When the unweighted sums fit, they also avoid prematurely halving
        // subnormal inventories. FMA remainders retain large-datum differences.
        let full_mass = a + b;
        let full_inventory = if full_mass.is_finite() {
            sum([p, q, dp, dq]).ok()
        } else {
            None
        };
        let (inventory, denominator) = if let Some(inventory) = full_inventory {
            (inventory, full_mass)
        } else {
            (sum([midpoint(p, q), midpoint(dp, dq)])?, mass)
        };
        let value = finite(inventory / denominator, "SSPRK2 blended specific quantity")?;
        if value == 0.0 && inventory != 0.0 {
            return Err(DensityTransportError::Arithmetic(
                "SSPRK2 blended specific quantity underflow",
            ));
        }
        specific[component] = value;
    }
    Ok(DensityTransportCell {
        density_kg_m3: density,
        specific,
    })
}

fn sum<const N: usize>(values: [f64; N]) -> Result<f64> {
    let mut sum = Sum::default();
    for value in values {
        sum.add(value)?;
    }
    sum.value()
}

#[cfg(test)]
mod tests {
    use super::{DensityTransportCell, blend};

    #[test]
    fn extensive_blend_retains_inventory_when_mass_fraction_underflows() {
        // Isolate the conservative stage blend: the original Euler increment
        // owner has a narrower admitted range for some incoming differences.
        let a = DensityTransportCell {
            density_kg_m3: 1e-300,
            specific: [1e300, -1e300],
        };
        let b = DensityTransportCell {
            density_kg_m3: 1e300,
            specific: [0.0; 2],
        };
        assert_eq!(a.density_kg_m3 / b.density_kg_m3, 0.0);
        let result = blend(&a, &b, 1.0).unwrap();
        for (value, expected) in result.specific.into_iter().zip([0.5, -0.5]) {
            assert!((result.density_kg_m3 * value - expected).abs() < 3e-16);
        }
    }

    #[test]
    fn extensive_blend_keeps_subnormal_amounts_and_exact_common_specifics() {
        let tiny = f64::from_bits(1);
        let a = DensityTransportCell {
            density_kg_m3: tiny,
            specific: [1.0, 7.0],
        };
        let b = DensityTransportCell {
            density_kg_m3: tiny,
            specific: [0.0, 7.0],
        };
        let result = blend(&a, &b, 1.0).unwrap();
        assert_eq!(result.density_kg_m3, tiny);
        // Average the extensive sum before division: no premature half(tiny).
        assert_eq!(result.specific, [0.5, 7.0]);
    }
}
