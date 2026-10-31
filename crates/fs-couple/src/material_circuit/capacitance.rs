//! Nonlinear capacitor laws retained by the existing material-circuit owner.
//!
//! Resistor and contact updates change dissipation, not the stored charge law.
//! Replacing a charge law is a separate fixed-charge operation with explicit
//! parameter work integral_0^q (V_new - V_old) dq. The physical owner changing
//! that law must supply/receive this work when accepting the candidate.

use fs_phs::circuit::capacitance::{CapacitorChargeCurve, CapacitorCurve};
use fs_phs::circuit::{Branch, CircuitDae};

use super::storage::{MaterialStorageCandidate, StorageParameterWork};
use super::{MAX_STATES, MaterialCircuit, MaterialCircuitError, poll};

/// Original nonlinear dielectric material and fixed geometry in ordinary circuits.
pub mod dielectric;

fn invalid(message: &'static str) -> MaterialCircuitError {
    MaterialCircuitError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, MaterialCircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("capacitor storage/work is unrepresentable"))
    }
}
fn add(sum: &mut f64, correction: &mut f64, term: f64) -> Result<(), MaterialCircuitError> {
    let term = finite(term)?;
    let next = finite(*sum + term)?;
    *correction = finite(
        *correction
            + if sum.abs() >= term.abs() {
                (*sum - next) + term
            } else {
                (term - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}
fn law(curves: &[CapacitorCurve], branch: usize) -> Option<&CapacitorChargeCurve> {
    curves
        .binary_search_by_key(&branch, |law| law.branch)
        .ok()
        .map(|index| &curves[index].curve)
}
fn voltage(
    curve: Option<&CapacitorChargeCurve>,
    farads: f64,
    charge: f64,
) -> Result<f64, MaterialCircuitError> {
    match curve {
        Some(curve) => Ok(curve.voltage_v(charge)?),
        None => finite(charge / farads),
    }
}
fn energy(
    curve: Option<&CapacitorChargeCurve>,
    farads: f64,
    charge: f64,
) -> Result<f64, MaterialCircuitError> {
    let result = match curve {
        Some(curve) => curve.stored_energy_j(charge)?,
        None => finite(0.5 * charge * (charge / farads))?,
    };
    if result < 0.0 || (charge != 0.0 && result == 0.0) {
        return Err(invalid("nonzero capacitor energy disappeared"));
    }
    Ok(result)
}

// Integrate the voltage DIFFERENCE on the union of the original knots. Each
// cell is exactly linear, so trapezoidal integration is exact for the declared
// curves. Subtract voltages before multiplying by dq: biased endpoint energies
// can coincide after rounding while this parameter work is still resolvable.
pub(super) fn parameter_work(
    before: Option<&CapacitorChargeCurve>,
    after: Option<&CapacitorChargeCurve>,
    farads: f64,
    charge: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MaterialCircuitError> {
    if charge == 0.0 || before == after {
        return Ok(0.0);
    }
    let low = charge.min(0.0);
    let high = charge.max(0.0);
    let mut points = vec![0.0, charge];
    for curve in [before, after].into_iter().flatten() {
        // Each original curve has at most MAX_CAPACITOR_KNOTS knots. This
        // zero-time update is bounded independently of the time integrator.
        points.extend(
            curve
                .knots()
                .iter()
                .map(|&(q, _)| q)
                .filter(|&q| q > low && q < high),
        );
    }
    points.sort_by(f64::total_cmp);
    points.dedup();
    if charge < 0.0 {
        points.reverse();
    }
    let difference = |q| -> Result<f64, MaterialCircuitError> {
        finite(voltage(after, farads, q)? - voltage(before, farads, q)?)
    };
    let mut previous = difference(points[0])?;
    let (mut work, mut correction) = (0.0, 0.0);
    for cell in points.windows(2) {
        poll(cancelled)?;
        let next = difference(cell[1])?;
        let average = f64::midpoint(previous, next);
        let area = finite((cell[1] - cell[0]) * average)?;
        if average != 0.0 && area == 0.0 {
            return Err(invalid("capacitor parameter work underflowed"));
        }
        add(&mut work, &mut correction, area)?;
        previous = next;
    }
    finite(work + correction)
}

impl MaterialCircuit {
    /// Replace the complete set of explicitly authored capacitor curves while
    /// retaining charges, flux linkage, contacts, mutual inductance and topology.
    /// An omitted branch returns to its original graph capacitance; its energy
    /// change is included. Unchanged branches contribute exactly zero work.
    ///
    /// Each proposed curve must contain the accepted charge. The native owner
    /// checks ideal-source consistency without an unaccounted charge impulse.
    /// Linear material-derived capacitor records cannot simultaneously claim
    /// authority for a different curve and therefore refuse replacement here.
    /// This operation infers no temperature dependence or material provenance.
    /// Work is reversible parameter work, never extra resistor heat; publish
    /// the candidate with the physical owner that pays/receives that work.
    pub fn rebind_capacitor_curves_candidate(
        &self,
        curves: Vec<CapacitorCurve>,
        inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialStorageCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        if self.dielectric_capacitors.iter().any(|source| {
            curves
                .iter()
                .find(|curve| curve.branch == source.branch())
                .is_none_or(|curve| &curve.curve != source.curve())
        }) {
            return Err(invalid(
                "retained dielectric sources require their explicit source-work owner",
            ));
        }
        self.rebind_capacitor_curves_with_work(curves, inputs, cancelled)
    }

    // Both explicit source rebinds and authored-law rebinds use the SAME native
    // admission and union-of-knots parameter-work account. Public source-aware
    // entry points retain their respective material-authority checks.
    fn rebind_capacitor_curves_with_work(
        &self,
        curves: Vec<CapacitorCurve>,
        inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialStorageCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        if curves.len() > self.graph.branches.len() {
            return Err(invalid("capacitor laws exceed bounded branch capacity"));
        }
        if curves.iter().any(|curve| {
            self.storage
                .iter()
                .any(|source| source.branch == curve.branch)
        }) {
            return Err(invalid(
                "an authored capacitor curve cannot replace a retained linear material law",
            ));
        }
        let mut candidate = Self::new_with_capacitor_curves(
            self.graph.clone(),
            self.resistors.clone(),
            self.contacts.clone(),
            self.mutual.clone(),
            curves,
            self.fluxes_wb.clone(),
            self.charges_c.clone(),
            inputs,
        )?;
        candidate.storage = self.storage.clone();
        candidate.magnetic_couplings = self.magnetic_couplings.clone();
        candidate.dielectric_capacitors = self.dielectric_capacitors.clone();
        let mut work = Vec::new();
        let (mut total, mut correction) = (0.0, 0.0);
        let mut slot = 0;
        for (branch, &(_, _, element)) in self.graph.branches.iter().enumerate() {
            let Branch::Capacitor { farads } = element else {
                continue;
            };
            poll(&mut cancelled)?;
            let charge = self.charges_c[slot];
            slot += 1;
            let old = law(&self.capacitor_curves, branch);
            let new = law(&candidate.capacitor_curves, branch);
            if old.is_none() && new.is_none() {
                continue;
            }
            let stored_before_j = energy(old, farads, charge)?;
            let stored_after_j = energy(new, farads, charge)?;
            let supplied_work_j = parameter_work(old, new, farads, charge, &mut cancelled)?;
            add(&mut total, &mut correction, supplied_work_j)?;
            work.push(StorageParameterWork {
                branch,
                stored_before_j,
                stored_after_j,
                supplied_work_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(MaterialStorageCandidate {
            circuit: candidate,
            work,
            parameter_work_j: finite(total + correction)?,
        })
    }
}

/// Independent constitutive increment for the frozen ordinary circuit law.
/// Capacitors use their declared integral; coupled magnetic storage uses the
/// SAME native inverse-inductance gradient, including all shared-field terms.
/// Fixed buffers bound work and avoid an extra matrix clone or factorization.
pub(super) fn storage_change(
    graph: &fs_phs::circuit::CircuitGraph,
    curves: &[CapacitorCurve],
    dae: &CircuitDae,
    before: &[f64],
    after: &[f64],
) -> Result<f64, MaterialCircuitError> {
    let n = dae.system.state_dim();
    if n > MAX_STATES || before.len() != n || after.len() != n {
        return Err(invalid(
            "complete bounded circuit states required for storage change",
        ));
    }
    let (mut total, mut correction) = (0.0, 0.0);
    if !dae.flux_index.is_empty() {
        let mut midpoint = [0.0; MAX_STATES];
        let mut offset = [0.0; MAX_STATES];
        let mut effort = [0.0; MAX_STATES];
        let mut effort_offset = [0.0; MAX_STATES];
        for &i in &dae.flux_index {
            midpoint[i] = f64::midpoint(before[i], after[i]);
            offset[i] = 0.5 * (before[i] - midpoint[i]) + 0.5 * (after[i] - midpoint[i]);
        }
        // Capacitor charges are zero in these magnetic probes, inside every
        // admitted curve's support. Algebraic efforts have no stored energy.
        dae.system
            .storage()
            .gradient(&midpoint[..n], &mut effort[..n]);
        dae.system
            .storage()
            .gradient(&offset[..n], &mut effort_offset[..n]);
        for &i in &dae.flux_index {
            let change = finite(after[i] - before[i])?;
            add(&mut total, &mut correction, change * effort[i])?;
            add(&mut total, &mut correction, change * effort_offset[i])?;
        }
    }
    let mut slot = 0;
    for (branch, &(_, _, element)) in graph.branches.iter().enumerate() {
        let Branch::Capacitor { farads } = element else {
            continue;
        };
        let i = dae.charge_index[slot];
        slot += 1;
        if let Some(curve) = law(curves, branch) {
            add(
                &mut total,
                &mut correction,
                curve.energy_change_j(before[i], after[i])?,
            )?;
        } else {
            let change = finite(after[i] - before[i])?;
            let midpoint = f64::midpoint(before[i], after[i]);
            let offset = 0.5 * (before[i] - midpoint) + 0.5 * (after[i] - midpoint);
            add(&mut total, &mut correction, change * (midpoint / farads))?;
            add(&mut total, &mut correction, change * (offset / farads))?;
        }
    }
    finite(total + correction)
}
