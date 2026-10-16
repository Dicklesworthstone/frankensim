//! Resistance controls on the bound contact operator, not rebuilt trace geometry.
use super::{ConductionError, Cx, ThermalInterfaces};

/// One shared multiplier on every retained resistance of a named interface.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactResistanceGradient {
    /// Exact bound interface identity, in canonical name order.
    pub interface: String,
    /// dJ/ds at s=1 for R''_face -> s R''_face, also dJ/dln(s) there.
    /// Units are the objective's units; this is NOT dJ/dR'' for a mapped patch.
    pub derivative: f64,
    /// Original material authority, not a new or modified card claim.
    pub card_identity: fs_blake3::ContentHash,
    /// Whether the original interface carries one resistance per face pair.
    pub mapped: bool,
    /// Paired faces contracted for this named control.
    pub face_pairs: usize,
}

impl ThermalInterfaces {
    /// Contract complete nodal-load adjoints with all matching contact controls.
    ///
    /// For R''_face -> s R''_face at s=1, dK_contact/ds = -K_contact and
    /// dJ/ds = lambda^T K_contact T. Both traces remain independent. Each face
    /// uses its own retained resistance and the consistent P1 mass matrix,
    /// including nonuniform primal/dual jumps; products of patch means are not
    /// substituted. No geometry search, matrix assembly or solve is performed.
    ///
    /// The caller supplies the unchanged full temperature and TOTAL load dual
    /// from the same checked physical solve, including nonlinear conductivity,
    /// air and radiation feedback as applicable. Fixed-node dual entries must
    /// be zero; prescribed primal values must NOT be zeroed. This contraction
    /// cannot verify the caller's state/dual residual, topology or boundary
    /// binding. Add any explicit objective dependence on resistance separately.
    /// A transient consumer sums endpoint contractions without an extra dt.
    ///
    /// A mapped interface scales ALL its own values together, preserving their
    /// spatial pattern; separate interface names remain independent even when
    /// they share a card. This does not alter, sample or validate a card. These
    /// are Estimated local derivatives, not material uncertainty, moving-contact
    /// sensitivities, temperature-dependent contact laws or error certificates.
    ///
    /// This entry point deliberately requires a matching-only interface set.
    /// The existing nonmatching_log_resistance_pullback remains the separate
    /// common-refinement path; a delegated matching subpatch must never be
    /// silently published as the entire nonmatching interface's derivative.
    ///
    /// max_face_pairs bounds total trace work and output row count, not RSS or
    /// string storage. Canonical face/name order is retained. Cancellation is
    /// polled before work, between surfaces and every 512 paired faces.
    ///
    /// # Errors
    /// Nonmatching contacts, mismatched field lengths, missing/nonfinite trace
    /// entries, nonfinite arithmetic, exhausted trace work or allocation, or
    /// cancellation. No partial gradient list is returned.
    pub fn matching_resistance_scale_pullback(
        &self, cx: &Cx<'_>, temperature: &[f64], nodal_load_adjoint: &[f64],
        max_face_pairs: usize,
    ) -> Result<Vec<ContactResistanceGradient>, ConductionError> {
        poll(cx, 0)?;
        if !self.nonmatching.is_empty() {
            return Err(error("matching resistance controls cannot omit nonmatching contact pieces"));
        }
        self.matching.resistance_scale_pullback(cx, temperature, nodal_load_adjoint, max_face_pairs)
    }
}

pub(super) fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "contact-resistance-controls", at })
}

pub(super) fn error(what: impl Into<String>) -> ConductionError {
    ConductionError::Interface { interface: "<contact-resistance-controls>".into(), what: what.into(),
        fix: "supply the same retained contact operator, full physical field and complete nodal-load dual within the trace-work allowance".into() }
}

pub(super) fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) } else { Err(error("nonfinite contact-resistance contraction")) }
}
