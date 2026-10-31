//! Source-owned history after an independently admitted contact and drift.
//!
//! The declared first-order quadrature holds the actual post-contact tangent
//! velocity over the accepted interval: delta-slip = h v and path = |h v|.
//! It describes that discrete motion, not a continuously varying material-point
//! velocity or a curved-slip quadrature. Refine/partition the outer mechanical
//! path when that approximation is insufficient.
//!
//! Source speed and accumulated distance refresh from this accepted motion,
//! preserving the original claims and each source's accepted path origin.
//! Temperature and pressure come only from explicit physical conditions. The
//! native zero-storage history records resolved slip without another traction,
//! compliance or impulse solve. The mechanical owner's existing dissipation
//! receipt remains separate; no part of it is converted into heat here.

use super::motion::{MotionQueryError, query_for_distance};
use super::{
    FrictionQueryAxisKind, SourceFrictionConditions, SourceThermalCoulombError,
    SourceThermalEllipticCoulomb, SourceThermalEllipticState, principal_coefficients,
};
use fs_tribo::TriboError;

/// Independent absolute length gates on represented vector and source history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourcePostcontactMotionControls {
    /// Maximum roundoff in each accumulated native slip component [m].
    pub absolute_slip_tolerance_m: f64,
    /// Maximum roundoff in every owned accumulated source distance [m].
    pub absolute_path_tolerance_m: f64,
}

/// Candidate original history after the declared held-velocity finite motion.
///
/// There is no force, normal load, heat, or dissipation field. Those quantities
/// belong to the independently accepted mechanical/thermal owners. The native
/// source state has no clock; the caller advances its accepted scene clock once
/// alongside this candidate and the matching mechanical interval.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePostcontactMotionReport {
    pub state: SourceThermalEllipticState,
    pub duration_s: f64,
    /// Actual post-contact velocity times the accepted duration, in the supplied
    /// phase tangent frame. It excludes all instantaneous velocity jumps.
    pub slip_increment_m: [f64; 2],
    /// Unsigned length of this held-slip segment, including after reversal [m].
    pub path_increment_m: f64,
    /// Segment path / duration, used for every owned source speed coordinate.
    pub queried_speed_m_s: f64,
    /// Actual representable native slip increment minus `slip_increment_m`.
    pub slip_residual_m: [f64; 2],
    /// Largest absolute represented-minus-requested source path increment.
    /// Zero if no accumulated-distance source coordinate is declared.
    pub maximum_path_coordinate_residual_m: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SourcePostcontactMotionError {
    Source(SourceThermalCoulombError),
    Motion(MotionQueryError),
    Invalid(&'static str),
    PathResolution {
        principal_axis: usize,
        coordinate: String,
        requested_increment_m: f64,
        represented_increment_m: f64,
        tolerance_m: f64,
    },
    Cancelled,
}
impl From<SourceThermalCoulombError> for SourcePostcontactMotionError {
    fn from(error: SourceThermalCoulombError) -> Self {
        Self::Source(error)
    }
}
impl From<MotionQueryError> for SourcePostcontactMotionError {
    fn from(error: MotionQueryError) -> Self {
        Self::Motion(error)
    }
}
impl From<TriboError> for SourcePostcontactMotionError {
    fn from(error: TriboError) -> Self {
        Self::Source(error.into())
    }
}
impl core::fmt::Display for SourcePostcontactMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source directional accepted-motion history: {self:?}")
    }
}
impl std::error::Error for SourcePostcontactMotionError {}

impl SourceThermalEllipticCoulomb {
    /// Accept the finite post-contact phase into pinned source and slip history.
    ///
    /// Call this only after independently admitting the complete mechanical
    /// interval and its geometric path. Input and returned history share the
    /// supplied phase tangent frame. Existing objective transport must handle
    /// any resolved material-frame rotation before/after this operation; normal
    /// impulse, relative torsion or an arbitrary gauge angle cannot supply it.
    ///
    /// The original source owner replaces speed/path placeholders. All other
    /// coordinates retain their existing ownership and range checks. Endpoint
    /// coefficients refresh for subsequent motion; they do not retroactively
    /// change the admitted impulse. Fixed work is two original source refreshes
    /// and one native history update. Cancellation, a source refusal or an
    /// unrepresentable increment publishes no partial candidate.
    pub fn accept_postcontact_motion(
        &self,
        state: &SourceThermalEllipticState,
        velocity_m_s: [f64; 2],
        dt_s: f64,
        conditions: &SourceFrictionConditions,
        controls: SourcePostcontactMotionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourcePostcontactMotionReport, SourcePostcontactMotionError> {
        poll(cancelled)?;
        self.check(state)?;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || velocity_m_s.into_iter().any(|v| !v.is_finite())
            || [
                controls.absolute_slip_tolerance_m,
                controls.absolute_path_tolerance_m,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(SourcePostcontactMotionError::Invalid(
                "finite post-contact velocity, positive duration and nonnegative history gates required",
            ));
        }
        let slip_increment_m = velocity_m_s.map(|v| dt_s * v);
        if slip_increment_m
            .iter()
            .zip(velocity_m_s)
            .any(|(&slip, speed)| !slip.is_finite() || (speed != 0.0 && slip == 0.0))
        {
            return Err(SourcePostcontactMotionError::Invalid(
                "post-contact slip increment is unrepresentable",
            ));
        }
        let path_increment_m = fs_math::det::hypot(slip_increment_m[0], slip_increment_m[1]);
        let queried_speed_m_s = path_increment_m / dt_s;
        if !path_increment_m.is_finite()
            || !queried_speed_m_s.is_finite()
            || (path_increment_m > 0.0 && queried_speed_m_s == 0.0)
        {
            return Err(SourcePostcontactMotionError::Invalid(
                "post-contact path or source speed is unrepresentable",
            ));
        }
        let mut maximum_path_coordinate_residual_m = 0.0_f64;
        let mut refresh = |principal_axis: usize| {
            poll(cancelled)?;
            let previous = &state.sources[principal_axis];
            let source = &self.principal[principal_axis];
            let previous_point = previous.query();
            let point = query_for_distance(
                &conditions.query,
                &previous_point,
                source.query_axes(),
                path_increment_m,
                dt_s,
            )?;
            for axis in source.query_axes() {
                poll(cancelled)?;
                if axis.kind == FrictionQueryAxisKind::AccumulatedDistance {
                    let represented = point.axes()[&axis.name] - previous_point.axes()[&axis.name];
                    let residual = represented - path_increment_m;
                    if !represented.is_finite()
                        || !residual.is_finite()
                        || (path_increment_m > 0.0 && represented <= 0.0)
                        || residual.abs() > controls.absolute_path_tolerance_m
                    {
                        return Err(SourcePostcontactMotionError::PathResolution {
                            principal_axis,
                            coordinate: axis.name.clone(),
                            requested_increment_m: path_increment_m,
                            represented_increment_m: represented,
                            tolerance_m: controls.absolute_path_tolerance_m,
                        });
                    }
                    maximum_path_coordinate_residual_m =
                        maximum_path_coordinate_residual_m.max(residual.abs());
                }
            }
            let resolved = source.source_at_query_conditions(
                previous,
                conditions.contact_temperature_k,
                conditions.nominal_pressure_pa,
                &point,
            )?;
            poll(cancelled)?;
            Ok(resolved)
        };
        let sources = [refresh(0)?, refresh(1)?];
        let (static_mu, kinetic_mu) =
            principal_coefficients(&self.principal, [&sources[0], &sources[1]])?;
        let (law, continued) =
            state
                .law
                .continue_coefficients(&state.constitutive, static_mu, kinetic_mu)?;
        let native = law.accept_resolved_slip(
            &continued,
            slip_increment_m,
            conditions.contact_temperature_k,
            controls.absolute_slip_tolerance_m,
        )?;
        poll(cancelled)?;
        Ok(SourcePostcontactMotionReport {
            state: SourceThermalEllipticState {
                owner: self.identity,
                law,
                constitutive: native.state,
                sources,
            },
            duration_s: dt_s,
            slip_increment_m,
            path_increment_m,
            queried_speed_m_s,
            slip_residual_m: native.slip_residual_m,
            maximum_path_coordinate_residual_m,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SourcePostcontactMotionError> {
    if cancelled() {
        Err(SourcePostcontactMotionError::Cancelled)
    } else {
        Ok(())
    }
}
