//! Independent mixed source, viscous and boundary-work accounting.

use super::super::super::symmetric_entry;
use super::{
    HarmonicFrequencyCurvatureResponse, HarmonicParameterResponse,
    HarmonicPortFrequencyParameterResponse, HarmonicResponseError, PowerSum,
    PreparedHarmonicResponse, Term, invalid, poll,
};
use crate::quadratic_modes::QuadraticParameterDirection;
use fs_exec::Cx;
use fs_math::c64::C64;

pub(super) struct MixedPower {
    pub(super) input: f64,
    pub(super) viscous: f64,
    pub(super) balance: f64,
}

// Re(a^H*b)*c. Keep c as an original physical factor: first forming 2*c
// to account for a symmetric pair could overflow despite finite total work.
fn full_inner(
    sum: &mut PowerSum,
    a: C64,
    b: C64,
    coefficient: f64,
) -> Result<(), HarmonicResponseError> {
    sum.add(Term::product([coefficient, a.re, b.re, 1.0], 0))?;
    sum.add(Term::product([coefficient, a.im, b.im, 1.0], 0))
}

/// Contract admitted original phasors and physical directions. No allocation,
/// perturbed model, modal loss or inferred load heat/storage is introduced.
#[allow(clippy::too_many_arguments)]
pub(super) fn observe(
    owner: &PreparedHarmonicResponse,
    frequency: &HarmonicFrequencyCurvatureResponse,
    parameter: &HarmonicParameterResponse,
    velocity_mixed: &[C64],
    force: &[C64],
    force_frequency_derivative: Option<&[C64]>,
    force_parameter_derivative: Option<&[C64]>,
    force_mixed_derivative: Option<&[C64]>,
    direction: QuadraticParameterDirection<'_>,
    port_load: Option<&HarmonicPortFrequencyParameterResponse>,
    cx: &Cx<'_>,
) -> Result<MixedPower, HarmonicResponseError> {
    poll(cx)?;
    let first = &frequency.frequency_response;
    let mut input = PowerSum::default();
    let mut viscous = PowerSum::default();
    for i in 0..owner.n {
        poll(cx)?;
        let v = first.response.velocity[i];
        let vw = first.velocity_derivative[i];
        let vt = parameter.velocity_derivative[i];
        let vwt = velocity_mixed[i];
        let fw = force_frequency_derivative.map_or(C64::ZERO, |values| values[i]);
        let ft = force_parameter_derivative.map_or(C64::ZERO, |values| values[i]);
        let fwt = force_mixed_derivative.map_or(C64::ZERO, |values| values[i]);
        input.half_inner(fwt, v, 1.0)?;
        input.half_inner(fw, vt, 1.0)?;
        input.half_inner(ft, vw, 1.0)?;
        input.half_inner(force[i], vwt, 1.0)?;

        for j in 0..owner.n {
            let damping = if owner.layout.columns(i, owner.n).contains(&j) {
                owner.viscous_entry(i, j)
            } else {
                0.0
            };
            // Form each symmetric physical coefficient BEFORE contraction.
            // A changing gyroscope affects motion, but cannot add direct heat.
            let damping_parameter = symmetric_entry(direction.damping, owner.n, i, j);
            full_inner(&mut viscous, vwt, first.response.velocity[j], damping)?;
            full_inner(&mut viscous, vw, parameter.velocity_derivative[j], damping)?;
            full_inner(
                &mut viscous,
                vw,
                first.response.velocity[j],
                damping_parameter,
            )?;
        }
    }
    let input = input.value()?;
    let viscous = viscous.value()?;
    let load = port_load.map_or(0.0, |port| port.average_absorbed_power_mixed_derivative);
    let mut balance = PowerSum::default();
    for value in [input, -viscous, -load] {
        balance.add_scalar(value)?;
    }
    let balance = balance.value()?;
    if ![input, viscous, balance].into_iter().all(f64::is_finite) {
        return Err(invalid(
            "mixed derivative of physical work is not representable",
        ));
    }
    poll(cx)?;
    Ok(MixedPower {
        input,
        viscous,
        balance,
    })
}
