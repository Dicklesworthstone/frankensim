//! Original speed-dependent source claims at every actual Hertz pressure node.
//!
//! A spinning patch has no single coefficient or source-refresh checkpoint.
//! Each local response therefore retains its own original use receipts. Native
//! history and accepted path coordinates remain untouched; this owner evaluates
//! an instantaneous prescribed velocity field, without inventing a distributed
//! slip history, normal work, thermal update or static capture decision.

use super::*;
use crate::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseSample;
use fs_tribo::sliding_patch::elliptic::resolved::{
    HertzLocalPoint, ResolvedEllipticHertzPatch, ResolvedHertzError, ResolvedHertzQuadrature,
    ResolvedHertzWrench,
};

/// One actual velocity/pressure node and its original selected principal claims.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzLocalSample {
    pub point: HertzLocalPoint,
    pub law: EllipticCoulomb,
    pub sources: [ResolvedInterfaceStatePoint; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceResolvedHertzReport {
    pub normal_model_identity: ContentHash,
    pub normal: PointNormalPatchReceipt,
    pub patch: Option<ResolvedEllipticHertzPatch>,
    pub response: ResolvedHertzWrench,
    /// Conservative full-footprint speed interval, admitted before integration.
    pub speed_envelope_m_s: [f64; 2],
    /// Source-domain probes at both envelope ends; not physical patch nodes.
    pub speed_envelope_samples: [SourceImpulseSample; 2],
    /// Canonically ordered actual nodes, bounded by maximum_evaluations.
    pub samples: Vec<SourceHertzLocalSample>,
}

#[derive(Debug)]
pub enum SourceResolvedHertzError {
    Admission(SourceEllipticHertzError),
    Native(ResolvedHertzError<SourceImpulseError>),
    Invalid(&'static str),
    Cancelled,
}
impl From<SourceEllipticHertzError> for SourceResolvedHertzError {
    fn from(error: SourceEllipticHertzError) -> Self {
        match error {
            SourceEllipticHertzError::Cancelled => Self::Cancelled,
            other => Self::Admission(other),
        }
    }
}
impl From<SourceImpulseError> for SourceResolvedHertzError {
    fn from(error: SourceImpulseError) -> Self {
        SourceEllipticHertzError::from(error).into()
    }
}
impl From<TriboError> for SourceResolvedHertzError {
    fn from(error: TriboError) -> Self {
        SourceEllipticHertzError::from(error).into()
    }
}
impl From<NormalPatchError> for SourceResolvedHertzError {
    fn from(error: NormalPatchError) -> Self {
        SourceEllipticHertzError::from(error).into()
    }
}
impl From<ResolvedHertzError<SourceImpulseError>> for SourceResolvedHertzError {
    fn from(error: ResolvedHertzError<SourceImpulseError>) -> Self {
        match error {
            ResolvedHertzError::Cancelled
            | ResolvedHertzError::Source(SourceImpulseError::Cancelled) => Self::Cancelled,
            other => Self::Native(other),
        }
    }
}
impl core::fmt::Display for SourceResolvedHertzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-resolved Hertz friction: {self:?}")
    }
}
impl std::error::Error for SourceResolvedHertzError {}

impl SourceThermalEllipticCoulomb {
    /// Integrate actual local-speed queries over the original normal footprint.
    ///
    /// Every node refreshes the SAME pinned static/kinetic claims from the
    /// accepted checkpoint, using the local norm of `u + omega J r`. The native
    /// associated ellipse supplies traction, and the same pressure measure owns
    /// force, centre moment and physical heat. Neither centre nor RMS speed is
    /// substituted for the field. The original uniform-source evaluator remains
    /// the cheaper analytic-radial path whenever its assumptions hold.
    ///
    /// The two original bulk temperatures and explicitly supplied nominal
    /// pressure keep their meanings. Local Hertz pressure is a mechanical
    /// weight, not a replacement source coordinate. The interval
    /// `[max(0, |u|-|omega|max(a,b)), |u|+|omega|max(a,b)]` is a conservative
    /// speed envelope; both endpoints must be source-admissible. It can refuse
    /// an otherwise admissible narrow ellipse rather than extrapolating a law
    /// between sampled nodes. Its two source probes are additional to the
    /// bounded native node count. Resolution must be refined by the caller;
    /// force/power residuals are not spatial convergence certificates.
    ///
    /// No node becomes the global history. The returned node receipts are
    /// instantaneous samples; future finite motion needs its own distributed
    /// history owner. Any source/normal/native/cancellation failure returns no
    /// patch response or source candidate, and the accepted checkpoint is intact.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_hertz_patch_resolved(
        &self,
        state: &SourceThermalEllipticState,
        normal: &SourcePointContact,
        input: SourceEllipticHertzInput<'_>,
        conditions: &SourceFrictionConditions,
        quadrature: ResolvedHertzQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceResolvedHertzReport, SourceResolvedHertzError> {
        quadrature.validate()?;
        let actual = self.hertz_conditions(state, normal, input, conditions, cancelled)?;
        let response = normal.evaluate(
            input.normal_state_id,
            input.indentation_m,
            input.indentation_rate_m_s,
            input.sample_duration_s,
        )?;
        let semi_axes = response
            .elliptic_patch_axes
            .map_or([response.patch_radius_m; 2], |axes| {
                [axes.semi_major_axis_m, axes.semi_minor_axis_m]
            });
        let speed = fs_math::det::hypot(input.slip_velocity_m_s[0], input.slip_velocity_m_s[1]);
        let edge = input.spin_rad_s.abs() * semi_axes[0].max(semi_axes[1]);
        let speed_envelope_m_s = [(speed - edge).max(0.0), speed + edge];
        if !speed.is_finite() || !edge.is_finite() || !speed_envelope_m_s[1].is_finite() {
            return Err(SourceResolvedHertzError::Invalid(
                "representable local-speed envelope required",
            ));
        }
        if input.spin_rad_s != 0.0 && semi_axes[0].max(semi_axes[1]) > 0.0 && edge == 0.0 {
            return Err(SourceResolvedHertzError::Invalid(
                "local spin speed is not representable",
            ));
        }
        if edge > 0.0 && speed_envelope_m_s[0] == speed_envelope_m_s[1] {
            return Err(SourceResolvedHertzError::Invalid(
                "local-speed variation is not representable",
            ));
        }
        let mut probe = |speed| {
            self.sample_preimpact_coefficients(
                state,
                [speed, 0.0],
                &actual,
                SourceImpulseCoefficientRole::Kinetic,
                cancelled,
            )
        };
        let speed_envelope_samples = [probe(speed_envelope_m_s[0])?, probe(speed_envelope_m_s[1])?];
        let mut samples = Vec::new();
        let (patch, integrated) = if input.indentation_m == 0.0 {
            (
                None,
                ResolvedHertzWrench {
                    wrench: HertzSlidingWrench {
                        force_n: [0.0; 2],
                        spin_moment_nm: 0.0,
                        dissipated_power_w: 0.0,
                        power_residual_w: 0.0,
                        at_rest: input.slip_velocity_m_s == [0.0; 2] && input.spin_rad_s == 0.0,
                    },
                    evaluations: 0,
                    integrated_normal_force_n: 0.0,
                    normal_force_residual_n: 0.0,
                },
            )
        } else {
            let patch = ResolvedEllipticHertzPatch::new(EllipticHertzPatchParameters {
                normal_force_n: response.normal_force_n,
                semi_axes_m: semi_axes,
                footprint_axis: input.footprint_axis,
            })?;
            samples
                .try_reserve_exact(quadrature.maximum_evaluations)
                .map_err(|_| {
                    SourceResolvedHertzError::Invalid("local source receipt allocation refused")
                })?;
            let integrated = patch.evaluate(
                input.slip_velocity_m_s,
                input.spin_rad_s,
                quadrature,
                |point, cancelled| {
                    let sampled = self.sample_preimpact_coefficients(
                        state,
                        point.slip_velocity_m_s,
                        &actual,
                        SourceImpulseCoefficientRole::Kinetic,
                        cancelled,
                    )?;
                    let law = EllipticCoulomb::new(
                        sampled.state.law.static_coefficients(),
                        sampled.state.law.kinetic_coefficients(),
                        sampled.state.constitutive.principal_axis(),
                    )?;
                    samples.push(SourceHertzLocalSample {
                        point,
                        law,
                        sources: sampled.state.sources,
                    });
                    Ok::<_, SourceImpulseError>(law)
                },
                cancelled,
            )?;
            (Some(patch), integrated)
        };
        poll(cancelled)?;
        Ok(SourceResolvedHertzReport {
            normal_model_identity: normal.model().identity(),
            normal: response,
            patch,
            response: integrated,
            speed_envelope_m_s,
            speed_envelope_samples,
            samples,
        })
    }
}
