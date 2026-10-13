//! Frictionless contact between independently deformable, history-aware roads.
//!
//! Each row measures `gap = reference_gap + w_upper - w_lower`. Concatenate
//! the signed curvature influence factors of ALL participating beams; their
//! Gram matrix is the relative flexibility. The existing nonnegative quadratic
//! owner solves one reaction per row. That same reaction enters the upper beam
//! upward and the lower beam downward. Neither beam is a rigid obstacle.
//!
//! This is fixed-pair, vertical contact at material-reference stations in
//! parallel x/z frames, not current-world contact discovery. Reference gaps
//! include the supplied initial surface offsets/thicknesses. Tangential slip,
//! evolving normals, area/pressure, adhesion, impact and drawing tension are
//! not represented. Built-in pin/clamp restraints remain bilateral. Contact
//! is pointwise; gaps between stations need separate spatial assessment.

use std::collections::BTreeMap;
use fs_material::thermovisco::Energy;
use fs_solver::orthant::{GramOrthantError, solve_gram_orthant};
use super::{Bridge, Config as BeamConfig, Support, supports};

/// Borrowed input for one physical road; all roads must share an accepted time.
#[derive(Clone, Copy)]
pub struct Body<'s, 'm> {
    pub id: u64,
    pub beam: &'s Bridge<'m>,
    pub temperatures_k: &'s BTreeMap<(u64, u64), f64>,
    pub gravity_m_s2: f64,
    pub additional_downward_n_m: &'s [f64],
    pub config: BeamConfig,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Station {
    pub body: u64,
    pub reference_x_m: f64,
}

/// One declared pair, with upward force on `upper` and its negative on `lower`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    pub id: u64,
    pub upper: Station,
    pub lower: Station,
    /// Signed undeformed surface separation, m, before adding w_upper-w_lower.
    /// Changing this at an endpoint is explicit imposed motion, not free work.
    pub reference_gap_m: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// Explicit retained-body ceiling, in 1..=64.
    pub max_bodies: usize,
    /// One shared allowance: columns count cells across ALL bodies; rows count
    /// ALL contacts. Gap and force-times-gap gates use the final actual shapes.
    pub contact: supports::Config,
    /// Whole-assembly work/storage gate, J, independent of each beam's gate.
    pub energy_tolerance_j: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reaction {
    pub contact: Contact,
    pub force_on_upper_n: f64,
    pub gap_m: f64,
}

#[derive(Debug, Clone)]
pub struct Report {
    pub beams: BTreeMap<u64, super::Report>,
    /// Stable contact-ID order, independent of body/contact input ordering.
    pub reactions: Vec<Reaction>,
    pub contact_sweeps: usize,
    pub maximum_gap_error_m: f64,
    pub maximum_complementarity_j: f64,
    /// Sum of gravity/additional-load endpoint work, excluding internal contact.
    pub external_work_j: f64,
    /// sum R_new*(delta_w_upper-delta_w_lower), J. Newly engaged contacts or
    /// changed reference gaps can give nonzero endpoint-quadrature work. This
    /// is not exact impact work, heat, or a whole-system passivity certificate.
    pub contact_work_j: f64,
    pub energy_j: Energy,
    pub energy_residual_j: f64,
}

/// Publish this whole result. No input body is changed on error/cancellation.
#[derive(Debug, Clone)]
pub struct Solution<'m> {
    pub bodies: BTreeMap<u64, Bridge<'m>>,
    pub report: Report,
}

#[derive(Debug, Clone)]
pub enum Error {
    Invalid(&'static str),
    Beam { body: u64, error: super::Error },
    Contact(GramOrthantError),
    Cancelled,
    Arithmetic,
    Residual { gap_m: f64, complementarity_j: f64 },
    Energy(f64),
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "neighbour contact refused: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<GramOrthantError> for Error {
    fn from(error: GramOrthantError) -> Self {
        if error == GramOrthantError::Cancelled { Self::Cancelled } else { Self::Contact(error) }
    }
}

/// Solve all relative contact forces before advancing any physical history.
/// Uses the original accepted Maxwell states in every compliance probe. The
/// complete assembly is staged: even failure while updating the last beam
/// returns no partial candidate. Pure immutable laws are required for replay.
///
/// At most 64 contacts and 4096 TOTAL curvature cells are admitted by the
/// existing numerical owner. No extra allowance is created per body or pair.
/// Dependent rows do not acquire unique reactions or an infeasibility proof.
/// Empty contacts give independent ordinary beam updates; removing a contact
/// does not erase prior stress. No thermal link or bond is inferred from force.
pub fn advance_controlled<'m>(bodies: &[Body<'_, 'm>], contacts: &[Contact],
    target_s: f64, config: Config, cancelled: &mut impl FnMut() -> bool,
) -> Result<Solution<'m>, Error> {
    poll(cancelled)?;
    let policy = config.contact.contact;
    solve_gram_orthant(&[], &[], policy, &mut *cancelled)?;
    if !(1..=64).contains(&config.max_bodies) || bodies.is_empty()
        || bodies.len() > config.max_bodies || !positive(config.energy_tolerance_j)
        || !positive(config.contact.gap_tolerance_m)
        || !positive(config.contact.complementarity_tolerance_j) {
        return Err(Error::Invalid("body count or energy/contact policy"));
    }
    let mut ordered = BTreeMap::new();
    let mut columns = 0usize;
    let time = bodies[0].beam.time_s();
    if !target_s.is_finite() || target_s <= time {
        return Err(Error::Invalid("target must advance the accepted physical time"));
    }
    for body in bodies {
        poll(cancelled)?;
        columns = columns.checked_add(body.beam.cells.len()).ok_or(Error::Arithmetic)?;
        if body.beam.time_s() != time || ordered.insert(body.id, body).is_some() {
            return Err(Error::Invalid("duplicate body ID or mismatched accepted times"));
        }
    }
    let rows = contacts.len();
    let products = rows.checked_mul(rows).and_then(|p| p.checked_mul(columns)).ok_or(Error::Arithmetic)?;
    if columns > policy.max_columns || rows > policy.max_variables || products > policy.max_setup_products {
        return Err(Error::Invalid("assembly exceeds shared contact setup/dimension budget"));
    }
    let mut pairs = BTreeMap::new();
    for &pair in contacts {
        poll(cancelled)?;
        if pair.upper.body == pair.lower.body || !pair.reference_gap_m.is_finite()
            || pairs.insert(pair.id, pair).is_some() {
            return Err(Error::Invalid("self-contact, duplicate contact ID or nonfinite gap"));
        }
        for station in [pair.upper, pair.lower] {
            let body = ordered.get(&station.body).ok_or(Error::Invalid("unknown contact body"))?;
            let x = station.reference_x_m;
            if !x.is_finite() || x <= 0.0 || x > body.beam.length_m
                || (body.beam.support == Support::SimplySupported && x == body.beam.length_m) {
                return Err(Error::Invalid("station outside beam or at a fixed endpoint"));
            }
        }
    }
    let mut predictions = BTreeMap::new();
    let mut offset = 0;
    for (&id, body) in &ordered {
        poll(cancelled)?;
        let prediction = supports::predict(body.beam, target_s, body.temperatures_k,
            body.gravity_m_s2, body.additional_downward_n_m, cancelled)
            .map_err(|error| Error::Beam { body: id, error })?;
        predictions.insert(id, (offset, prediction));
        offset += body.beam.cells.len();
    }
    let mut factor = Vec::with_capacity(rows);
    let mut gaps = Vec::with_capacity(rows);
    for pair in pairs.values() {
        poll(cancelled)?;
        let mut row = vec![0.0; columns];
        let mut gap = pair.reference_gap_m;
        for (station, sign) in [(pair.upper, 1.0), (pair.lower, -1.0)] {
            let beam = ordered[&station.body].beam;
            let (offset, prediction) = &predictions[&station.body];
            for (k, cell) in beam.cells.iter().enumerate() {
                if k % 256 == 0 { poll(cancelled)?; }
                let psi = sign*supports::curvature_basis(beam.support, beam.length_m,
                    cell.start_m, cell.length_m(), station.reference_x_m);
                row[*offset+k] = finite(psi*prediction.roots[k])?;
                gap = finite(gap+psi*prediction.free_curvature[k])?;
            }
        }
        factor.push(row); gaps.push(gap);
    }
    let solved = solve_gram_orthant(&factor, &gaps, policy, &mut *cancelled)?;
    let mut forces: BTreeMap<_, Vec<_>> = ordered.keys().map(|&id| (id, Vec::new())).collect();
    for (pair, &force) in pairs.values().zip(&solved.point) {
        poll(cancelled)?;
        for (station, sign) in [(pair.upper, 1.0), (pair.lower, -1.0)] {
            forces.get_mut(&station.body).expect("admitted body")
                .push((station.reference_x_m, sign*force));
        }
    }
    let mut next = BTreeMap::new();
    let mut reports = BTreeMap::new();
    let mut energy = Energy::default();
    let mut external_work = 0.0;
    let mut storage = 0.0;
    for (&id, body) in &ordered {
        poll(cancelled)?;
        let (beam, report) = body.beam.advance_with_point_forces(target_s, body.temperatures_k,
            body.gravity_m_s2, body.additional_downward_n_m, &forces[&id], body.config, cancelled)
            .map_err(|error| Error::Beam { body: id, error })?;
        super::add_energy(&mut energy, report.energy_j).map_err(|error| Error::Beam { body: id, error })?;
        storage = finite(storage+beam.stored_energy_j()-body.beam.stored_energy_j())?;
        // Independently integrate ONLY externally supplied loads, not a sum of
        // reports with internal forces silently classified as external work.
        for ((old, cell), &q_extra) in body.beam.cells.iter().zip(&beam.cells).zip(body.additional_downward_n_m) {
            poll(cancelled)?;
            let q = finite(body.gravity_m_s2*old.mass_per_length_kg_m+q_extra)?;
            external_work = finite(external_work-q*(cell.integral_height()-old.integral_height()))?;
        }
        next.insert(id, beam); reports.insert(id, report);
    }
    let mut reactions = Vec::with_capacity(rows);
    let mut gap_error = 0.0_f64;
    let mut complementarity = 0.0_f64;
    let mut contact_work = 0.0;
    for (pair, force) in pairs.into_values().zip(solved.point) {
        poll(cancelled)?;
        let mut gap = pair.reference_gap_m;
        let mut displacement = 0.0;
        for (station, sign) in [(pair.upper, 1.0), (pair.lower, -1.0)] {
            let old = height(ordered[&station.body].beam, station)?;
            let new = height(&next[&station.body], station)?;
            gap = finite(gap+sign*new)?;
            displacement = finite(displacement+sign*(new-old))?;
        }
        gap_error = gap_error.max(if force > 0.0 { gap.abs() } else { (-gap).max(0.0) });
        complementarity = complementarity.max(finite(force*gap)?.abs());
        contact_work = finite(contact_work+force*displacement)?;
        reactions.push(Reaction { contact: pair, force_on_upper_n: force, gap_m: gap });
    }
    if gap_error > config.contact.gap_tolerance_m || complementarity > config.contact.complementarity_tolerance_j {
        return Err(Error::Residual { gap_m: gap_error, complementarity_j: complementarity });
    }
    let work = finite(external_work+contact_work)?;
    let residual = finite(storage+energy.dissipated+energy.numerical_dissipation
        -work-energy.thermal_strain_work-energy.modulus_work)?;
    for defect in [residual, finite(work-energy.mechanical_work)?, finite(storage-energy.stored_change)?] {
        if defect.abs() > config.energy_tolerance_j { return Err(Error::Energy(defect)); }
    }
    energy.residual = residual;
    poll(cancelled)?;
    Ok(Solution { bodies: next, report: Report { beams: reports, reactions, contact_sweeps: solved.sweeps,
        maximum_gap_error_m: gap_error, maximum_complementarity_j: complementarity,
        external_work_j: external_work, contact_work_j: contact_work, energy_j: energy,
        energy_residual_j: residual } })
}

/// Executor-scoped form. Cancellation has the same whole-assembly semantics.
pub fn advance_to<'m>(cx: &fs_exec::Cx<'_>, bodies: &[Body<'_, 'm>], contacts: &[Contact],
    target_s: f64, config: Config) -> Result<Solution<'m>, Error> {
    advance_controlled(bodies, contacts, target_s, config, &mut || cx.checkpoint().is_err())
}
fn height(beam: &Bridge<'_>, station: Station) -> Result<f64, Error> {
    let (cell, s) = beam.locate(station.reference_x_m)
        .map_err(|error| Error::Beam { body: station.body, error })?;
    finite(cell.height(s))
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() { Ok(value) } else { Err(Error::Arithmetic) }
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests;
