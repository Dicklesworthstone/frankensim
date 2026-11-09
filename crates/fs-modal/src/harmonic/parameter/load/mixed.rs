//! Mixed material/attachment actions in original conjugate port coordinates.
//!
//! With L=B*Z*B^T, the extra mixed mechanical source is
//! -L_theta*v_omega - L_omega_theta*v. Both occurrences of B_theta
//! and the authored impedance mixed derivative remain factored through ports.

use super::super::super::HarmonicFrequencyResponse;
use super::{
    ComplexSum, HarmonicPortParameterDirection, HarmonicResponseError, PortAction, PowerSum,
    PreparedHarmonicResponse, Term, invalid, poll,
};
use fs_exec::Cx;
use fs_math::c64::C64;

/// Complete load direction for a mixed frequency/physical-parameter response.
///
/// B and B_theta are independent of angular frequency. The physical coordinate
/// and port dimensions stay fixed; all matrices use the original ordering.
/// Every impedance derivative is a partial derivative at this same frequency
/// and parameter, not a difference of neighboring impedance samples.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPortFrequencyParameterDirection<'a> {
    /// Complete B_theta and Z_theta in original force/flow and impedance units.
    pub parameter: HarmonicPortParameterDirection<'a>,
    /// Complete row-major Z_omega_theta, per rad/s per parameter unit.
    /// Signed, mutual and nonreciprocal terms are retained. An explicit zero
    /// matrix declares a zero mixed derivative; it is never inferred from
    /// either Z_theta or Z_omega at one sampled frequency.
    pub impedance_angular_parameter_derivative: &'a [C64],
}

impl PortAction<'_> {
    /// Use the original e_omega in the outer attachment contribution -B_theta*e_omega.
    /// The remaining complete inner derivative was retained in normalized form.
    pub(in super::super) fn visit_frequency_rhs(
        &self,
        owner: &PreparedHarmonicResponse,
        frequency: &HarmonicFrequencyResponse,
        cx: &Cx<'_>,
        visit: impl FnMut(usize, usize, Term, Term) -> Result<(), HarmonicResponseError>,
    ) -> Result<(), HarmonicResponseError> {
        let Some(port) = &frequency.port_load else {
            return Err(invalid(
                "mixed port action requires the original port frequency response",
            ));
        };
        self.visit_projected_rhs(owner, &port.effort_derivative, cx, visit)
    }
}

/// Prepare only direct load-law terms at fixed physical state parameter.
///
/// Let y=B^T*v, y_omega=B^T*v_omega, b=B_theta^T*v and
/// b_omega=B_theta^T*v_omega. The inner effort is
/// Z_theta*y_omega + Z_omega_theta*y + Z*b_omega + Z_omega*b.
/// The outer projection adds -B_theta*e_omega and -B times this effort.
/// The separate complete -A_omega*q_theta action stays with the existing
/// frequency owner; it must not be counted here a second time.
pub(in super::super) fn frequency_parameter_action<'a>(
    owner: &PreparedHarmonicResponse,
    frequency: &HarmonicFrequencyResponse,
    direction: Option<HarmonicPortFrequencyParameterDirection<'a>>,
    impedance_derivative: Option<&[C64]>,
    cx: &Cx<'_>,
) -> Result<Option<PortAction<'a>>, HarmonicResponseError> {
    let Some(direction) = direction else {
        return Ok(None);
    };
    let (Some(load), Some(primal), Some(slope), Some(impedance_slope)) = (
        &owner.port_load,
        &frequency.response.port_load,
        &frequency.port_load,
        impedance_derivative,
    ) else {
        return Err(invalid(
            "complete original load and frequency observations required",
        ));
    };
    let p = load.ports;
    let mut flow = Vec::with_capacity(p);
    let mut flow_slope = Vec::with_capacity(p);
    let mut flow_envelope = Vec::with_capacity(p);
    let mut flow_slope_envelope = Vec::with_capacity(p);
    for a in 0..p {
        poll(cx)?;
        let mut value = ComplexSum::default();
        let mut derivative = ComplexSum::default();
        let mut bound = PowerSum::default();
        let mut derivative_bound = PowerSum::default();
        for i in 0..owner.n {
            let db = direction.parameter.force_map[i * p + a];
            let v = frequency.response.velocity[i];
            let vw = frequency.velocity_derivative[i];
            for (axis, (v, vw)) in [(v.re, vw.re), (v.im, vw.im)].into_iter().enumerate() {
                let term = Term::product([db, v, 1.0, 1.0], 0);
                let slope_term = Term::product([db, vw, 1.0, 1.0], 0);
                value.0[axis].add(term)?;
                derivative.0[axis].add(slope_term)?;
                bound.add(term.abs())?;
                derivative_bound.add(slope_term.abs())?;
            }
        }
        flow.push(value.terms()?);
        flow_slope.push(derivative.terms()?);
        flow_envelope.push(bound.term()?);
        flow_slope_envelope.push(derivative_bound.term()?);
    }
    let mut effort = Vec::with_capacity(p);
    let mut envelope = Vec::with_capacity(p);
    for a in 0..p {
        poll(cx)?;
        let mut value = ComplexSum::default();
        let mut bound = PowerSum::default();
        for b in 0..p {
            let at = a * p + b;
            for (z, direct, direct_bound) in [
                (load.impedance[at], flow_slope[b], flow_slope_envelope[b]),
                (impedance_slope[at], flow[b], flow_envelope[b]),
            ] {
                value.scaled_product(z, direct)?;
                for part in [z.re, z.im] {
                    bound.add(direct_bound.times(part.abs(), 0))?;
                }
            }
            for (z, y) in [
                (direction.parameter.impedance[at], slope.flow_derivative[b]),
                (
                    direction.impedance_angular_parameter_derivative[at],
                    primal.flow[b],
                ),
            ] {
                value.complex_product(z, y)?;
                for z in [z.re, z.im] {
                    for y in [y.re, y.im] {
                        bound.add(Term::product([z, y, 1.0, 1.0], 0).abs())?;
                    }
                }
            }
        }
        effort.push(value.terms()?);
        envelope.push(bound.term()?);
    }
    Ok(Some(PortAction {
        map_direction: Some(direction.parameter.force_map),
        effort,
        envelope,
    }))
}
