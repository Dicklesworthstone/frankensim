//! Affine contact geometry on the original cohesive straight displacement path.
//!
//! This module samples actual native states. It never interpolates an accepted
//! motion, changes a gap, or certifies the first event of a continuous trajectory.
//! With fixed traces and affine prescribed motion, nonnegative endpoint gaps
//! admit the entire declared straight coordinate path. The caller owns complete
//! feature coverage, the fixed orthonormal frames, and the traces' geometry.

use super::super::impact::TraceTerm;
use super::super::impact::frictional::{FrictionalContact, NormalVelocityLaw};

/// A refused geometry trial leaves every original mechanical state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum GeometryError {
    Invalid(&'static str),
    Budget(&'static str),
    NonFinite(&'static str),
    Cancelled,
    Penetration { contact_id: u64, gap_m: f64 },
    NoProgress { step_s: f64 },
}
impl core::fmt::Display for GeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive affine contact geometry refused: {self:?}")
    }
}
impl std::error::Error for GeometryError {}

/// Original contact order. Negative proposed gaps remain available for retry.
pub(super) struct Sample {
    pub gaps_m: Vec<f64>,
    pub normal_velocities_m_s: Vec<f64>,
    /// Zero for an empty contact family; otherwise the smallest signed gap.
    pub minimum_gap_m: f64,
}

impl Sample {
    /// One active closing contact, selected by stable ID. The surrounding owner
    /// validates the finite nonnegative tolerance and rejects current penetration
    /// before requesting an event. All active contacts still enter one solve.
    pub(super) fn closing_contact(
        &self,
        contacts: &[FrictionalContact<'_>],
        velocity_tolerance_m_s: f64,
    ) -> Option<usize> {
        contacts
            .iter()
            .zip(&self.gaps_m)
            .zip(&self.normal_velocities_m_s)
            .enumerate()
            .filter(|(_, ((contact, gap), speed))| {
                **gap >= 0.0
                    && **gap <= contact.activation_tolerance_m
                    && **speed < -velocity_tolerance_m_s
            })
            .min_by_key(|(_, ((contact, _), _))| contact.id)
            .map(|(index, _)| index)
    }
}

/// Admit the fixed contact family before any native dynamics work.
///
/// Positive activation bands are required for bounded numerical approach. All
/// three conjugate trace rows use the original impact admission convention.
/// The workspace ceiling conservatively includes slice-table slots, two
/// retained samples, and scalar scratch; native body states are borrowed.
pub(super) fn admit(
    contacts: &[FrictionalContact<'_>],
    reference_positions: &[&[f64]],
    maximum_trace_terms: usize,
    maximum_workspace_entries: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GeometryError> {
    poll(cancelled)?;
    if contacts.len() > 21 {
        return Err(GeometryError::Budget("native Coulomb contact count"));
    }
    let workspace = reference_positions
        .len()
        .checked_mul(8)
        .and_then(|body| {
            contacts
                .len()
                .checked_mul(4)
                .and_then(|rows| body.checked_add(rows))
        })
        .and_then(|entries| entries.checked_add(32))
        .ok_or(GeometryError::Budget("geometry workspace overflow"))?;
    if workspace > maximum_workspace_entries {
        return Err(GeometryError::Budget("affine geometry workspace"));
    }
    validate_vectors(reference_positions, reference_positions, cancelled)?;
    let mut count = 0usize;
    for (index, contact) in contacts.iter().enumerate() {
        poll(cancelled)?;
        if !contact.gap_m.is_finite() {
            return Err(GeometryError::NonFinite("initial gap"));
        }
        if contact.gap_m < 0.0 {
            return Err(GeometryError::Penetration {
                contact_id: contact.id,
                gap_m: contact.gap_m,
            });
        }
        if contacts[..index].iter().any(|prior| prior.id == contact.id)
            || !contact.activation_tolerance_m.is_finite()
            || contact.activation_tolerance_m <= 0.0
            || !contact.friction_coefficient.is_finite()
            || contact.friction_coefficient < 0.0
            || contact
                .prescribed_velocity_m_s
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(GeometryError::Invalid(
                "unique contact IDs, positive activation bands and finite law data required",
            ));
        }
        match contact.normal_law {
            NormalVelocityLaw::NewtonRestitution { coefficient }
                if coefficient.is_finite() && (0.0..=1.0).contains(&coefficient) => {}
            NormalVelocityLaw::Target { speed } if speed.is_finite() && speed >= 0.0 => {}
            _ => return Err(GeometryError::Invalid("native normal velocity law")),
        }
        for terms in contact.trace {
            count = count
                .checked_add(terms.len())
                .ok_or(GeometryError::Budget("trace count overflow"))?;
            if count > maximum_trace_terms {
                return Err(GeometryError::Budget("contact trace terms"));
            }
            admit_trace(terms, reference_positions, cancelled)?;
        }
    }
    poll(cancelled)
}

/// Sample g = g_initial + B_normal (u - u_initial) + c_normal * elapsed,
/// and v_normal = B_normal v + c_normal. Differences precede trace summation,
/// avoiding subtraction of two large complete reference/current trace values.
///
/// Inputs must retain the admitted contact and body ordering. Proposed negative
/// gaps are returned without clamping so the owner can shorten and RECOMPUTE
/// native candidates. No material state from a rejected trial is reused.
pub(super) fn sample(
    contacts: &[FrictionalContact<'_>],
    reference_positions: &[&[f64]],
    positions: &[&[f64]],
    velocities: &[&[f64]],
    elapsed_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Sample, GeometryError> {
    poll(cancelled)?;
    if !elapsed_s.is_finite() || elapsed_s < 0.0 {
        return Err(GeometryError::Invalid(
            "finite nonnegative elapsed time required",
        ));
    }
    if contacts.len() > 21 {
        return Err(GeometryError::Budget("native Coulomb contact count"));
    }
    validate_vectors(positions, reference_positions, cancelled)?;
    validate_vectors(velocities, reference_positions, cancelled)?;
    let mut gaps_m = Vec::with_capacity(contacts.len());
    let mut normal_velocities_m_s = Vec::with_capacity(contacts.len());
    for contact in contacts {
        poll(cancelled)?;
        let mut gap = Sum::default();
        let mut speed = Sum::default();
        gap.add(contact.gap_m)?;
        gap.add(finite(
            elapsed_s * contact.prescribed_velocity_m_s[0],
            "prescribed gap increment",
        )?)?;
        speed.add(contact.prescribed_velocity_m_s[0])?;
        for (index, term) in contact.trace[0].iter().enumerate() {
            if index % 256 == 0 {
                poll(cancelled)?;
            }
            let before = coordinate(reference_positions, term)?;
            let after = coordinate(positions, term)?;
            let delta = finite(after - before, "native displacement increment")?;
            gap.add(finite(
                term.coefficient * delta,
                "traced displacement increment",
            )?)?;
            speed.add(finite(
                term.coefficient * coordinate(velocities, term)?,
                "traced normal velocity",
            )?)?;
        }
        gaps_m.push(gap.total()?);
        normal_velocities_m_s.push(speed.total()?);
    }
    let minimum_gap_m = gaps_m.iter().copied().reduce(f64::min).unwrap_or(0.0);
    poll(cancelled)?;
    Ok(Sample {
        gaps_m,
        normal_velocities_m_s,
        minimum_gap_m,
    })
}

/// Propose a strictly shorter positive duration after a penetrating trial.
///
/// Each crossing aims halfway into its declared activation band. The most
/// restrictive fraction wins, capped at one half. Already-near-contact starts
/// also halve: their recomputed dynamics may depart instead of penetrating.
/// No minimum-step floor or geometric projection is applied. The duration is
/// only a retry proposal; the nonlinear native endpoint must be solved again.
pub(super) fn shorter_step(
    contacts: &[FrictionalContact<'_>],
    start_gaps: &[f64],
    end_gaps: &[f64],
    dt_s: f64,
) -> Result<f64, GeometryError> {
    if contacts.len() > 21
        || contacts.len() != start_gaps.len()
        || contacts.len() != end_gaps.len()
        || !dt_s.is_finite()
        || dt_s <= 0.0
    {
        return Err(GeometryError::Invalid(
            "matching bounded gap arrays and positive finite step required",
        ));
    }
    let mut fraction = 0.5_f64;
    let mut crossing = false;
    for ((contact, &start), &end) in contacts.iter().zip(start_gaps).zip(end_gaps) {
        finite(start, "start gap")?;
        finite(end, "proposed gap")?;
        if start < 0.0 {
            return Err(GeometryError::Penetration {
                contact_id: contact.id,
                gap_m: start,
            });
        }
        if !contact.activation_tolerance_m.is_finite() || contact.activation_tolerance_m <= 0.0 {
            return Err(GeometryError::Invalid(
                "positive finite activation band required",
            ));
        }
        if end >= 0.0 {
            continue;
        }
        crossing = true;
        let target = 0.5 * contact.activation_tolerance_m;
        if target == 0.0 {
            return Err(GeometryError::NoProgress { step_s: dt_s });
        }
        if start > target {
            // Scaling avoids overflow in start-end for large opposite gaps.
            let scale = start.max(-end);
            let proposed = ((start - target) / scale) / (start / scale - end / scale);
            if !proposed.is_finite() || proposed <= 0.0 {
                return Err(GeometryError::NoProgress { step_s: dt_s });
            }
            fraction = fraction.min(proposed);
        }
    }
    if !crossing {
        return Err(GeometryError::Invalid(
            "shortening requires a penetrating proposed path",
        ));
    }
    let next = dt_s * fraction;
    if !next.is_finite() || next <= 0.0 || next >= dt_s {
        return Err(GeometryError::NoProgress { step_s: dt_s });
    }
    Ok(next)
}

pub(super) fn admit_trace(
    terms: &[TraceTerm],
    reference_positions: &[&[f64]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GeometryError> {
    let mut previous = None;
    for (index, term) in terms.iter().enumerate() {
        if index % 256 == 0 {
            poll(cancelled)?;
        }
        let key = (term.body, term.dof);
        if term.body >= reference_positions.len()
            || term.dof >= reference_positions[term.body].len()
            || !term.coefficient.is_finite()
            || term.coefficient == 0.0
            || previous.is_some_and(|old| old >= key)
        {
            return Err(GeometryError::Invalid(
                "strictly sorted unique finite nonzero trace entries required",
            ));
        }
        previous = Some(key);
    }
    Ok(())
}

pub(super) fn validate_vectors(
    values: &[&[f64]],
    reference_positions: &[&[f64]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GeometryError> {
    if values.is_empty() || values.len() != reference_positions.len() {
        return Err(GeometryError::Invalid(
            "original nonempty body ordering required",
        ));
    }
    for (body, reference) in values.iter().zip(reference_positions) {
        poll(cancelled)?;
        if body.is_empty() || body.len() != reference.len() {
            return Err(GeometryError::Invalid(
                "original body coordinate dimensions required",
            ));
        }
        for chunk in body.chunks(256) {
            poll(cancelled)?;
            for &value in chunk {
                finite(value, "body coordinate or velocity")?;
            }
        }
    }
    Ok(())
}

pub(super) fn coordinate(values: &[&[f64]], term: &TraceTerm) -> Result<f64, GeometryError> {
    values
        .get(term.body)
        .and_then(|body| body.get(term.dof))
        .copied()
        .ok_or(GeometryError::Invalid(
            "trace outside original body coordinates",
        ))
}

#[derive(Default)]
pub(super) struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    pub(super) fn add(&mut self, value: f64) -> Result<(), GeometryError> {
        finite(value, "affine trace contribution")?;
        let next = finite(self.value + value, "affine trace sum")?;
        let residual = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.correction = finite(self.correction + residual, "affine trace correction")?;
        self.value = next;
        Ok(())
    }
    pub(super) fn total(self) -> Result<f64, GeometryError> {
        finite(self.value + self.correction, "affine trace result")
    }
}

fn finite(value: f64, quantity: &'static str) -> Result<f64, GeometryError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GeometryError::NonFinite(quantity))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GeometryError> {
    if cancelled() {
        Err(GeometryError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "geometry/tests.rs"]
mod tests;
