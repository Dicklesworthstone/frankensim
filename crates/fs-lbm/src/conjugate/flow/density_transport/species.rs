//! Constituent masses on the original accepted material-transfer routes.
//!
//! Species indices have the caller's fixed meaning. This transports a partition
//! of the carrier mass; it does not infer an EOS, partial volumes, reactions,
//! diffusion, mixing heat or a new velocity from the resulting composition.

use super::{
    AcceptedMassTransport, CancelGate, ChtError, InletLocation, MaterialAdvectionStep, Route, Sum,
    finite, finite_positive, poll, poll_every,
};
use crate::conjugate::krylov::bicgstab_ilu0_prepared;
use fs_sparse::precond::{IluControlledError, ilu0_controlled};

mod diffusion;
pub use diffusion::{
    AdaptiveSpeciesDiffusionControls, AdaptiveSpeciesDiffusionStep,
    BinarySpeciesDiffusionControls, BinarySpeciesDiffusionReport, BinarySpeciesDiffusionStep,
    SpeciesDiffusionControls, SpeciesDiffusionReport, SpeciesDiffusionStep,
    SpeciesEnthalpyDiffusionReport, SpeciesEnthalpyDiffusionStep,
};

/// Hard component bound; the caller also supplies solve and storage budgets.
pub const MAX_MATERIAL_SPECIES: usize = 128;

/// Authoritative constituent inventories and authored incoming compositions.
#[derive(Clone, Copy, Debug)]
pub struct MaterialSpeciesInput<'a, const N: usize> {
    /// One nonnegative constituent-mass array [kg] per original grid cell.
    /// Its sum must agree with the accepted transport's ORIGINAL cell mass.
    /// Solid cells contain exact zeros. Retain returned masses between steps.
    pub species_mass_kg: &'a [[f64; N]],
    /// Mass fractions of local injected material in original grid-cell order.
    /// Every actual source inlet requires a value; extraction uses its donor.
    pub source_mass_fractions: &'a [Option<[f64; N]>],
    /// Incoming/backflow fractions in x-, x+, y-, y+, z-, z+ order.
    /// An outlet becoming an inlet uses its own explicitly supplied composition.
    pub boundary_mass_fractions: [Option<[f64; N]>; 6],
}

/// Numerical work and independent physical-unit admission.
#[derive(Clone, Copy, Debug)]
pub struct MaterialSpeciesControls<const N: usize> {
    /// Native iteration ceiling per component; N is separately bounded.
    pub max_iterations: usize,
    pub relative_tolerance: f64,
    /// Per-species sum of local donor-equation, inventory-representation and
    /// global external-transfer defects [kg]. Trace species keep their OWN
    /// allowance rather than inheriting a bulk-mixture absolute tolerance.
    pub maximum_mass_defect_kg: [f64; N],
    /// Sum of absolute initial-cell, endpoint-cell and route partition
    /// defects [kg], all compared with the original accepted carrier masses.
    /// This permits declared roundoff; it never normalizes a composition.
    pub maximum_partition_defect_kg: f64,
    /// Conservative payload bound for output, one scalar transport workspace,
    /// original Krylov scratch and one shared ILU(0) factor. Existing accepted
    /// mass-transport storage and allocator overhead are additional.
    pub maximum_storage_bytes: usize,
}

/// Actual constituent inventory, external transfers and independent defects.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialSpeciesReport<const N: usize> {
    pub mass_before_kg: [f64; N],
    pub mass_after_kg: [f64; N],
    pub boundary_incoming_mass_kg: [f64; N],
    pub source_incoming_mass_kg: [f64; N],
    /// Boundary outflow and local extraction, using the accepted donor.
    pub outgoing_mass_kg: [f64; N],
    /// Direct sum of represented per-cell endpoint-minus-initial increments.
    pub mass_change_kg: [f64; N],
    pub mass_balance_defect_kg: [f64; N],
    pub absolute_mass_defect_kg: [f64; N],
    pub initial_partition_defect_kg: f64,
    pub endpoint_partition_defect_kg: f64,
    pub route_partition_defect_kg: f64,
    pub maximum_cell_partition_defect_kg: f64,
    pub iterations: [usize; N],
}

/// Complete constituent candidate. Neither the carrier step nor input changes.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialSpeciesStep<const N: usize> {
    /// Authoritative endpoint constituent masses [kg], in original cell order.
    /// Fractions are derived using the SAME MaterialAdvectionStep.mass_kg.
    pub species_mass_kg: Vec<[f64; N]>,
    pub report: MaterialSpeciesReport<N>,
}

fn refuse(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.material_species",
        reason: reason.into(),
    }
}

fn nonnegative(value: f64) -> Result<(), ChtError> {
    if !value.is_finite() || value < 0.0 {
        Err(refuse(
            "constituent mass or fraction is negative or nonfinite",
        ))
    } else {
        Ok(())
    }
}

fn fractions<const N: usize>(values: &[f64; N]) -> Result<(), ChtError> {
    let mut sum = Sum::default();
    sum.add(-1.0);
    for &value in values {
        nonnegative(value)?;
        if value > 1.0 {
            return Err(refuse("authored inlet mass fraction exceeds one"));
        }
        sum.add(value);
    }
    // Admission of authored unit fractions, not a physical mixing tolerance.
    if sum.value.abs() > 32.0 * N as f64 * f64::EPSILON {
        return Err(refuse("authored inlet fractions do not sum to one"));
    }
    Ok(())
}

fn inlet<'a, const N: usize>(
    route: &Route,
    input: &'a MaterialSpeciesInput<'_, N>,
) -> Result<&'a [f64; N], ChtError> {
    match route.inlet_location {
        Some(InletLocation::Boundary(face)) => input.boundary_mass_fractions[face].as_ref(),
        Some(InletLocation::Source(cell)) => input.source_mass_fractions[cell].as_ref(),
        None => None,
    }
    .ok_or_else(|| refuse("incoming mass route lacks its original inlet composition"))
}

fn ratio(mass: f64, carrier: f64) -> Result<f64, ChtError> {
    let value = mass / carrier;
    nonnegative(value)?;
    if mass > 0.0 && value == 0.0 {
        return Err(refuse("positive constituent fraction is not representable"));
    }
    // Do not impose an exact upper probability bound on a derived fraction:
    // the independently admitted carrier-mass partition can differ by ulps.
    Ok(value)
}

fn transfer(mass: f64, fraction: f64) -> Result<f64, ChtError> {
    let value = mass * fraction;
    nonnegative(value)?;
    if fraction > 0.0 && value == 0.0 {
        return Err(refuse("positive constituent transfer is not representable"));
    }
    Ok(value)
}

fn admit<const N: usize>(
    t: &AcceptedMassTransport,
    input: MaterialSpeciesInput<'_, N>,
    controls: MaterialSpeciesControls<N>,
    gate: &CancelGate,
) -> Result<(f64, f64), ChtError> {
    let cells = t.row_of.len();
    if N == 0
        || N > MAX_MATERIAL_SPECIES
        || input.species_mass_kg.len() != cells
        || input.source_mass_fractions.len() != cells
        || controls.max_iterations == 0
        || controls.max_iterations.checked_mul(N).is_none()
    {
        return Err(refuse(
            "invalid species dimensions, component count or work budget",
        ));
    }
    finite_positive("flow.species_tolerance", controls.relative_tolerance)?;
    if controls.relative_tolerance >= 1.0 {
        return Err(refuse("species relative tolerance must be less than one"));
    }
    for bound in controls
        .maximum_mass_defect_kg
        .into_iter()
        .chain([controls.maximum_partition_defect_kg])
    {
        nonnegative(bound)?;
    }
    // 160 bytes/active row bounds the scalar arrays, seven native Krylov
    // vectors, ILU row/diagonal/setup indices and their simultaneous payload.
    let bytes = cells
        .checked_mul(N)
        .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
        .and_then(|v| v.checked_add(t.cells.len().checked_mul(160)?))
        .and_then(|v| v.checked_add(t.routes.len().checked_mul(core::mem::size_of::<Sum>())?))
        .and_then(|v| {
            v.checked_add(
                t.operator
                    .nnz()
                    .checked_mul(core::mem::size_of::<usize>() + core::mem::size_of::<f64>())?,
            )
        })
        .and_then(|v| v.checked_add(core::mem::size_of::<MaterialSpeciesStep<N>>()))
        .and_then(|v| v.checked_add(2048))
        .ok_or_else(|| refuse("species storage extent is not representable"))?;
    if bytes > controls.maximum_storage_bytes {
        return Err(refuse("species storage budget exhausted"));
    }
    for values in input.boundary_mass_fractions.iter().flatten() {
        fractions(values)?;
    }
    let mut partition = Sum::default();
    let mut maximum: f64 = 0.0;
    for c in 0..cells {
        poll_every(gate, c)?;
        for &mass in &input.species_mass_kg[c] {
            nonnegative(mass)?;
        }
        if let Some(values) = &input.source_mass_fractions[c] {
            fractions(values)?;
        }
        let row = t.row_of[c];
        if row == usize::MAX {
            if input.species_mass_kg[c] != [0.0; N] || input.source_mass_fractions[c].is_some() {
                return Err(refuse(
                    "solid cells cannot carry constituent inventories or sources",
                ));
            }
            continue;
        }
        let mut sum = Sum::default();
        sum.add(-t.initial_mass_kg[row]);
        for &mass in &input.species_mass_kg[c] {
            ratio(mass, t.initial_mass_kg[row])?;
            sum.add(mass);
        }
        partition.add(sum.value.abs());
        maximum = maximum.max(sum.value.abs());
    }
    for (index, route) in t.routes.iter().enumerate() {
        poll_every(gate, index)?;
        if route.from.is_none() {
            inlet(route, &input)?;
        }
    }
    finite("flow.species_initial_partition", partition.value)?;
    if partition.value > controls.maximum_partition_defect_kg {
        return Err(refuse(
            "original constituent masses do not partition the original carrier mass",
        ));
    }
    Ok((partition.value, maximum))
}

impl MaterialAdvectionStep {
    /// Transport 1..=128 constituent masses on this step's EXACT accepted mass
    /// routes, including source-driven and phase-volume candidates. No second
    /// density solve, volume-flux reconstruction, or material injection occurs.
    ///
    /// The original mass*h donor matrix and native BiCGStab iteration are reused.
    /// One original ILU(0) factor serves all species. Each connected component's
    /// largest fraction scales one scalar solve, so an entirely trace species
    /// is not discarded by squaring tiny values in a Krylov norm. Absent species
    /// remain exact zeros; no seeded trace, clipping or renormalization is used.
    /// Paired transfers, not M*an independently rounded fraction, determine the
    /// endpoint constituent masses. Actual returned fractions then pass their
    /// unscaled extensive donor equations and global conservation audits.
    ///
    /// Initial, endpoint and route partitions are audited against the SAME
    /// retained carrier inventories/transfers under the declared kg allowance.
    /// A pure component can differ from carrier mass by admitted roundoff; its
    /// derived fraction is not forced into an exact probability interval.
    ///
    /// This retains first-order implicit donor advection and the caller's fixed
    /// species ordering. Material calorics, reactions, molecular diffusion and
    /// composition-dependent EOS/volume feedback remain separate owners.
    /// Invalid data, changed public carrier masses, nonrepresentable traces,
    /// negative candidates, source refusal, failed readback or cancellation
    /// returns no partial species result and leaves all inputs unchanged.
    pub fn advect_species<const N: usize>(
        &self,
        input: MaterialSpeciesInput<'_, N>,
        controls: MaterialSpeciesControls<N>,
        gate: &CancelGate,
    ) -> Result<MaterialSpeciesStep<N>, ChtError> {
        self.advect_species_checked(input, controls, gate, |_| poll(gate))
    }

    // The same gate is checked after each completed component, before any
    // complete candidate can escape. The callback permits deterministic late
    // cancellation tests without changing the executor gate or solver.
    #[allow(clippy::too_many_lines)]
    fn advect_species_checked<const N: usize>(
        &self,
        input: MaterialSpeciesInput<'_, N>,
        controls: MaterialSpeciesControls<N>,
        gate: &CancelGate,
        mut component_checkpoint: impl FnMut(usize) -> Result<(), ChtError>,
    ) -> Result<MaterialSpeciesStep<N>, ChtError> {
        poll(gate)?;
        let t = self.accepted_mass_transport(gate)?;
        let (initial_partition, mut maximum_partition) = admit(t, input, controls, gate)?;
        let factor = ilu0_controlled(&t.operator, || poll(gate)).map_err(|error| match error {
            IluControlledError::Breakdown(error) => ChtError::PreconditionerBreakdown {
                system: "material species transport",
                row: error.row,
            },
            IluControlledError::Capacity => refuse("species factor allocation refused"),
            IluControlledError::Checkpoint(error) => error,
        })?;
        let n = t.cells.len();
        let mut species_mass_kg = vec![[0.0; N]; t.row_of.len()];
        let mut rhs = vec![0.0; n];
        let mut donor = vec![0.0; n];
        let mut scales = vec![0.0_f64; n];
        let mut deltas = vec![Sum::default(); n];
        let mut route_partition = vec![Sum::default(); t.routes.len()];
        for (sum, route) in route_partition.iter_mut().zip(&t.routes) {
            sum.add(-route.mass);
        }
        let mut report = MaterialSpeciesReport {
            mass_before_kg: [0.0; N],
            mass_after_kg: [0.0; N],
            boundary_incoming_mass_kg: [0.0; N],
            source_incoming_mass_kg: [0.0; N],
            outgoing_mass_kg: [0.0; N],
            mass_change_kg: [0.0; N],
            mass_balance_defect_kg: [0.0; N],
            absolute_mass_defect_kg: [0.0; N],
            initial_partition_defect_kg: initial_partition,
            endpoint_partition_defect_kg: 0.0,
            route_partition_defect_kg: 0.0,
            maximum_cell_partition_defect_kg: 0.0,
            iterations: [0; N],
        };
        for k in 0..N {
            poll(gate)?;
            scales.fill(0.0);
            for (row, &c) in t.cells.iter().enumerate() {
                poll_every(gate, row)?;
                donor[row] = ratio(input.species_mass_kg[c][k], t.initial_mass_kg[row])?;
                scales[t.parents[row]] = scales[t.parents[row]].max(donor[row]);
            }
            for (index, route) in t.routes.iter().enumerate() {
                poll_every(gate, index)?;
                if route.from.is_none() {
                    let root = t.parents[route.to.expect("incoming receiver")];
                    scales[root] = scales[root].max(inlet(route, &input)?[k]);
                }
            }
            for (row, &c) in t.cells.iter().enumerate() {
                poll_every(gate, row)?;
                let scale = scales[t.parents[row]];
                if scale == 0.0 {
                    rhs[row] = 0.0;
                    donor[row] = 0.0;
                } else {
                    donor[row] /= scale;
                    rhs[row] = (input.species_mass_kg[c][k] / scale) / t.diagonal[row];
                    if input.species_mass_kg[c][k] > 0.0 && rhs[row] == 0.0 {
                        return Err(refuse("scaled constituent inventory is not representable"));
                    }
                    finite("flow.species_rhs", rhs[row])?;
                }
            }
            for (index, route) in t.routes.iter().enumerate() {
                poll_every(gate, index)?;
                if route.from.is_none() {
                    let to = route.to.expect("incoming receiver");
                    let scale = scales[t.parents[to]];
                    if scale > 0.0 {
                        let value = inlet(route, &input)?[k];
                        let increment = (route.mass / t.diagonal[to]) * (value / scale);
                        if value > 0.0 && increment == 0.0 {
                            return Err(refuse("scaled constituent inflow is not representable"));
                        }
                        rhs[to] += increment;
                        finite("flow.species_rhs", rhs[to])?;
                    }
                }
            }
            report.iterations[k] = bicgstab_ilu0_prepared(
                "material species transport",
                &t.operator,
                &rhs,
                &mut donor,
                controls.relative_tolerance,
                controls.max_iterations,
                gate,
                || Ok(&factor),
            )?
            .iterations;
            for row in 0..n {
                poll_every(gate, row)?;
                let scale = scales[t.parents[row]];
                let value = donor[row] * scale;
                nonnegative(value)?;
                if (donor[row] > 0.0 && value == 0.0) || (scale == 0.0 && donor[row] != 0.0) {
                    return Err(refuse("constituent donor fraction is not representable"));
                }
                donor[row] = value;
            }
            deltas.fill(Sum::default());
            let (mut boundary_in, mut source_in, mut outgoing) =
                (Sum::default(), Sum::default(), Sum::default());
            for (index, route) in t.routes.iter().enumerate() {
                poll_every(gate, index)?;
                let value = match route.from {
                    Some(from) => donor[from],
                    None => inlet(route, &input)?[k],
                };
                let amount = transfer(route.mass, value)?;
                route_partition[index].add(amount);
                if let Some(from) = route.from {
                    deltas[from].add(-amount);
                } else {
                    match route.inlet_location.expect("admitted original inlet") {
                        InletLocation::Boundary(_) => boundary_in.add(amount),
                        InletLocation::Source(_) => source_in.add(amount),
                    }
                }
                if let Some(to) = route.to {
                    deltas[to].add(amount);
                } else {
                    outgoing.add(amount);
                }
            }
            let (mut before, mut after, mut change, mut defect) = (
                Sum::default(),
                Sum::default(),
                Sum::default(),
                Sum::default(),
            );
            for (row, &c) in t.cells.iter().enumerate() {
                poll_every(gate, row)?;
                let old = input.species_mass_kg[c][k];
                let mass = old + deltas[row].value;
                nonnegative(mass)?;
                if mass == 0.0 && (old > 0.0 || donor[row] > 0.0) {
                    return Err(refuse(
                        "positive endpoint constituent inventory disappeared",
                    ));
                }
                species_mass_kg[c][k] = mass;
                donor[row] = ratio(mass, t.endpoint_mass_kg[row])?;
                before.add(old);
                after.add(mass);
                change.add(mass - old);
                defect.add((deltas[row].value - (mass - old)).abs());
                deltas[row] = Sum::default();
                deltas[row].add(mass - old);
            }
            // Independent extensive donor-equation readback uses the actual
            // returned masses, not the scaled solver iterate or its residual.
            for (index, route) in t.routes.iter().enumerate() {
                poll_every(gate, index)?;
                let value = match route.from {
                    Some(from) => donor[from],
                    None => inlet(route, &input)?[k],
                };
                let amount = transfer(route.mass, value)?;
                if let Some(from) = route.from {
                    deltas[from].add(amount);
                }
                if let Some(to) = route.to {
                    deltas[to].add(-amount);
                }
            }
            for (row, balance) in deltas.iter().enumerate() {
                poll_every(gate, row)?;
                defect.add(balance.value.abs());
            }
            let balance = change.value - boundary_in.value - source_in.value + outgoing.value;
            defect.add(balance.abs());
            for value in [
                before.value,
                after.value,
                boundary_in.value,
                source_in.value,
                outgoing.value,
                change.value,
                balance,
                defect.value,
            ] {
                finite("flow.species_account", value)?;
            }
            if defect.value > controls.maximum_mass_defect_kg[k] {
                return Err(refuse(
                    "constituent transport exceeds its independent mass budget",
                ));
            }
            report.mass_before_kg[k] = before.value;
            report.mass_after_kg[k] = after.value;
            report.boundary_incoming_mass_kg[k] = boundary_in.value;
            report.source_incoming_mass_kg[k] = source_in.value;
            report.outgoing_mass_kg[k] = outgoing.value;
            report.mass_change_kg[k] = change.value;
            report.mass_balance_defect_kg[k] = balance;
            report.absolute_mass_defect_kg[k] = defect.value;
            component_checkpoint(k + 1)?;
        }
        let mut endpoint_partition = Sum::default();
        for (row, &c) in t.cells.iter().enumerate() {
            poll_every(gate, row)?;
            let mut sum = Sum::default();
            sum.add(-t.endpoint_mass_kg[row]);
            for &mass in &species_mass_kg[c] {
                sum.add(mass);
            }
            endpoint_partition.add(sum.value.abs());
            maximum_partition = maximum_partition.max(sum.value.abs());
        }
        let mut transferred_partition = Sum::default();
        for (index, sum) in route_partition.iter().enumerate() {
            poll_every(gate, index)?;
            transferred_partition.add(sum.value.abs());
        }
        let partition = initial_partition + endpoint_partition.value + transferred_partition.value;
        finite("flow.species_partition", partition)?;
        if partition > controls.maximum_partition_defect_kg {
            return Err(refuse(
                "constituent inventories/transfers do not partition accepted carrier mass",
            ));
        }
        report.endpoint_partition_defect_kg = endpoint_partition.value;
        report.route_partition_defect_kg = transferred_partition.value;
        report.maximum_cell_partition_defect_kg = maximum_partition;
        poll(gate)?;
        Ok(MaterialSpeciesStep {
            species_mass_kg,
            report,
        })
    }
}

#[cfg(test)]
mod tests;
