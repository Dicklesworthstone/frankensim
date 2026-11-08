//! External gust incidence from total and explicitly declared reference air.
//!
//! Both fields use the native strip point, body pose, rigid angular velocity
//! and normal-section frame. Their geometric incidence difference supplies the
//! existing Kussner channel; the solved effective incidence minus that gust
//! supplies the existing motion channel. Twist, wing induction and prop wash
//! keep their original owners. Only the total field drives physical velocity,
//! dynamic pressure, propellers and drag. This is a declared decomposition,
//! not a second gust force or a chord-resolved aerodynamic model.

use super::super::{AerodynamicStage, Refusal, array, reject, vec, wing_error};
use fs_wing::nonlinear::SectionFrame;
use std::borrow::Cow;

#[derive(Clone, Copy)]
pub(in crate::spine::rotors::propulsion) enum GustInput<'a> {
    Angles(&'a [f64]),
    Reference(&'a dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>),
}

impl<'a> GustInput<'a> {
    pub(in crate::spine::rotors::propulsion) fn invalid_angles(self, count: usize) -> bool {
        match self {
            Self::Angles(gust) => gust.len() != count || !gust.iter().all(|x| x.is_finite()),
            Self::Reference(_) => false,
        }
    }

    pub(in crate::spine::rotors::propulsion) fn resolve(
        self,
        stage: &AerodynamicStage<'_>,
        count: usize,
        frame: impl Fn(usize) -> SectionFrame,
    ) -> Result<Cow<'a, [f64]>, Refusal> {
        let reference = match self {
            Self::Angles(gust) => return Ok(Cow::Borrowed(gust)),
            Self::Reference(reference) => reference,
        };
        if stage.strips.len() != count || stage.incident.strip_velocity_mps.len() != count {
            return Err(reject(
                "unsteady-gust-reference-invalid",
                "incomplete native strip geometry or incident air",
            ));
        }
        let mut gust = Vec::with_capacity(count);
        for (s, strip) in stage.strips.iter().enumerate() {
            let panel = &stage.panels[strip.panel_indices[0]];
            // Preserve the native material-air sample point and its arithmetic.
            let point = array((vec(panel.a) + vec(panel.b)).scale(0.5));
            let reference = stage
                .air
                .sample_reference(point, reference)?
                .relative_velocity_mps;
            let actual = stage.incident.strip_velocity_mps[s];
            let frame = frame(s);
            let actual = frame.flow(actual).map_err(wing_error)?;
            let reference = frame.flow(reference).map_err(wing_error)?;
            // The native flow gate rejects nonpositive chordwise components
            // and zero/underflowed section speed. Each atan2 is strictly in
            // (-pi/2, pi/2), so no +/-pi branch wrap or reversal floor is needed.
            gust.push(actual.incidence_rad - reference.incidence_rad);
        }
        Ok(Cow::Owned(gust))
    }
}
