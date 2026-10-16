//! Deformation hook for current-configuration queries (bead rjnd, E1
//! query upgrades, part 6).
//!
//! A deformable adapter owns a reference-configuration chart plus a
//! motion state. [`DeformationMap`] is the narrow contract it must
//! satisfy — a pull-back into the reference configuration with a
//! certified global Lipschitz bound — and [`DeformedChart`] is the
//! honest chart the rest of the query layer can then consume.
//!
//! The mathematics is deliberately modest. Signed DISTANCE does not
//! survive a general deformation, but sign and zero set do:
//! `f(x) = φ_ref(pull_back(x))` vanishes exactly on the deformed
//! boundary and keeps φ's sign. When the reference field is an exact
//! distance (1-Lipschitz) and the pull-back is `L`-Lipschitz, the
//! composed field is `L`-Lipschitz — exactly the
//! [`TraceStepClaim::LipschitzImplicit`] theorem, so steppers get a
//! certified no-tunneling radius `|f|/L` in the CURRENT configuration
//! without any distance overclaim. Gradients are declined (the chain
//! rule needs a certified Jacobian; a later surface), and a pull-back
//! that produces non-finite points yields a no-claim sample that every
//! typed consumer refuses.

use crate::QueryError;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim};

/// A certified pull-back from the current configuration into a
/// reference configuration.
///
/// Contract: `pull_back` is defined on the declared current support
/// box and satisfies `|pull_back(x) - pull_back(y)| ≤ L·|x - y|` with
/// `L = pull_back_lipschitz()` — a certified GLOBAL bound the adapter
/// derives from its motion state (e.g. an inverse-deformation-gradient
/// singular-value bound). The map must be a bijection onto the
/// reference domain for the sign/zero-set transfer to mean anything;
/// that bijectivity (no self-interpenetration of the deformation) is
/// the adapter's certificate, not this crate's.
pub trait DeformationMap: Send + Sync {
    /// Map a current-configuration point to the reference
    /// configuration.
    fn pull_back(&self, x: Point3) -> Point3;

    /// Certified global Lipschitz bound of [`Self::pull_back`].
    fn pull_back_lipschitz(&self) -> f64;

    /// Stable name for refusal messages.
    fn name(&self) -> &'static str;
}

/// A reference chart presented in the current configuration through a
/// [`DeformationMap`], claiming exactly what survives the composition:
/// sign, zero set, and an `L`-Lipschitz certified field.
pub struct DeformedChart<'a> {
    reference: &'a dyn Chart,
    map: &'a dyn DeformationMap,
    lipschitz: f64,
    support: Aabb,
}

impl<'a> DeformedChart<'a> {
    /// Wrap a reference chart and a pull-back map.
    ///
    /// `current_support` is the adapter's declared bound on the
    /// deformed region (the pull-back's domain); like every chart
    /// support it is a caller assertion.
    ///
    /// # Errors
    /// [`QueryError::DeformationRequiresExactDistance`] unless the
    /// reference chart claims `ExactDistance` (the composed Lipschitz
    /// theorem needs the reference field's certified 1-Lipschitz
    /// bound); [`QueryError::DeformationInvalidMap`] for a non-finite
    /// or non-positive pull-back Lipschitz bound.
    pub fn new(
        reference: &'a dyn Chart,
        map: &'a dyn DeformationMap,
        current_support: Aabb,
    ) -> Result<DeformedChart<'a>, QueryError> {
        let claim = reference.trace_step_claim();
        if claim != TraceStepClaim::ExactDistance {
            return Err(QueryError::DeformationRequiresExactDistance { claim });
        }
        let l = map.pull_back_lipschitz();
        if !l.is_finite() || l <= 0.0 {
            return Err(QueryError::DeformationInvalidMap {
                reason: "pull-back Lipschitz bound must be finite and positive",
            });
        }
        Ok(DeformedChart {
            reference,
            map,
            lipschitz: l.next_up(),
            support: current_support,
        })
    }

    fn no_claim_sample() -> ChartSample {
        ChartSample {
            signed_distance: f64::NAN,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        }
    }
}

impl Chart for DeformedChart<'_> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        if !(x.x.is_finite() && x.y.is_finite() && x.z.is_finite()) {
            return Self::no_claim_sample();
        }
        let p = self.map.pull_back(x);
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            // A broken pull-back cannot poison downstream consumers:
            // the sample refuses through the certificate, loudly.
            return Self::no_claim_sample();
        }
        let mut sample = self.reference.eval(p, cx);
        // The reference enclosure encloses φ_ref(p), which IS the
        // composed field value at x — it passes through unchanged.
        // The composed field's certified Lipschitz bound replaces the
        // reference's local one (which lived in reference coordinates).
        sample.gradient = None;
        sample.lipschitz = Some(self.lipschitz);
        sample
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        // The pass-through reference enclosure rigorously encloses the
        // composed field evaluation; anything weaker stays no-claim.
        let sound = matches!(
            sample.error.kind,
            NumericalKind::Exact | NumericalKind::Enclosure
        ) && sample.error.lo.is_finite()
            && sample.error.hi.is_finite()
            && sample.error.lo <= sample.error.hi;
        if sound {
            sample.error
        } else {
            NumericalCertificate::no_claim()
        }
    }

    fn name(&self) -> &'static str {
        "query/deformed"
    }
}
