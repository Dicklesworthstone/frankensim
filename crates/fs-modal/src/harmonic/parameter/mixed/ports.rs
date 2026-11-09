//! Original conjugate port jets for mixed frequency/parameter observations.

use super::{
    HarmonicFrequencyCurvatureResponse, HarmonicParameterResponse,
    HarmonicPortFrequencyParameterDirection, HarmonicResponseError, PowerSum,
    PreparedHarmonicResponse, finite, invalid, poll,
};
use fs_exec::Cx;
use fs_math::c64::C64;

/// Mixed derivatives in the original physical flow/effort port coordinates.
///
/// These are partial derivatives per (rad/s) per caller-declared parameter
/// unit. B and B_theta are independent of frequency; their authored parameter
/// change remains in the flow. The impedance may depend on both variables.
#[derive(Clone, Debug)]
pub struct HarmonicPortFrequencyParameterResponse {
    /// B^T*v_omega_theta + B_theta^T*v_omega.
    pub flow_mixed_derivative: Vec<C64>,
    /// Z*flow_omega_theta + Z_omega*flow_theta
    /// + Z_theta*flow_omega + Z_omega_theta*flow.
    pub effort_mixed_derivative: Vec<C64>,
    /// Mixed derivative of signed cycle-mean Re(flow^H*effort)/2.
    /// Boundary work is not an inferred constitutive heat/storage split.
    pub average_absorbed_power_mixed_derivative: f64,
}

#[derive(Default)]
struct ComplexSum {
    real: PowerSum,
    imaginary: PowerSum,
}

impl ComplexSum {
    fn product(&mut self, a: C64, b: C64) -> Result<(), HarmonicResponseError> {
        // Every real factor stays binary-normalized until cancellation.
        // The exact factor 2 removes half_inner's 1/2 without first forming
        // an overflowing complex product or discarding a subnormal factor.
        self.real.half_inner(C64::new(a.re, -a.im), b, 2.0)?;
        self.imaginary.half_inner(C64::new(a.im, a.re), b, 2.0)
    }

    fn value(self) -> Result<C64, HarmonicResponseError> {
        let value = C64::new(self.real.value()?, self.imaginary.value()?);
        if !finite(value) {
            return Err(invalid(
                "mixed physical port observation is not representable",
            ));
        }
        Ok(value)
    }
}

/// Observe admitted, complete physical partial derivatives without another solve.
/// Work is O(n*p+p*p), with two retained length-p complex arrays and constant
/// scratch. The parent admits all shapes, finite inputs and aggregate resources.
#[allow(clippy::too_many_arguments)]
pub(super) fn observe(
    owner: &PreparedHarmonicResponse,
    frequency: &HarmonicFrequencyCurvatureResponse,
    parameter: &HarmonicParameterResponse,
    velocity_mixed: &[C64],
    impedance_derivative: Option<&[C64]>,
    direction: Option<HarmonicPortFrequencyParameterDirection<'_>>,
    cx: &Cx<'_>,
) -> Result<Option<HarmonicPortFrequencyParameterResponse>, HarmonicResponseError> {
    poll(cx)?;
    let Some(load) = &owner.port_load else {
        return Ok(None);
    };
    let first = &frequency.frequency_response;
    let (Some(primal), Some(angular), Some(parameter_port), Some(slope)) = (
        &first.response.port_load,
        &first.port_load,
        &parameter.port_load,
        impedance_derivative,
    ) else {
        return Err(invalid(
            "mixed port observation requires the complete original load jets",
        ));
    };
    let ports = load.ports;
    if velocity_mixed.len() != owner.n
        || first.velocity_derivative.len() != owner.n
        || [
            primal.flow.as_slice(),
            primal.effort.as_slice(),
            angular.flow_derivative.as_slice(),
            angular.effort_derivative.as_slice(),
            parameter_port.flow_derivative.as_slice(),
            parameter_port.effort_derivative.as_slice(),
        ]
        .into_iter()
        .any(|values| values.len() != ports)
        || slope.len() != load.impedance.len()
        || direction.is_some_and(|d| {
            d.parameter.force_map.len() != load.force_map.len()
                || d.parameter.impedance.len() != load.impedance.len()
                || d.impedance_angular_parameter_derivative.len() != load.impedance.len()
        })
    {
        return Err(invalid(
            "mixed port observations require matching original coordinate extents",
        ));
    }

    let mut flow = Vec::with_capacity(ports);
    for a in 0..ports {
        poll(cx)?;
        let mut sum = ComplexSum::default();
        for i in 0..owner.n {
            sum.product(
                C64::from_re(load.force_map[i * ports + a]),
                velocity_mixed[i],
            )?;
            if let Some(direction) = direction {
                sum.product(
                    C64::from_re(direction.parameter.force_map[i * ports + a]),
                    first.velocity_derivative[i],
                )?;
            }
        }
        flow.push(sum.value()?);
    }

    let mut effort = Vec::with_capacity(ports);
    let mut power = PowerSum::default();
    for a in 0..ports {
        poll(cx)?;
        let mut sum = ComplexSum::default();
        for b in 0..ports {
            let at = a * ports + b;
            sum.product(load.impedance[at], flow[b])?;
            sum.product(slope[at], parameter_port.flow_derivative[b])?;
            if let Some(direction) = direction {
                sum.product(
                    direction.parameter.impedance[at],
                    angular.flow_derivative[b],
                )?;
                sum.product(
                    direction.impedance_angular_parameter_derivative[at],
                    primal.flow[b],
                )?;
            }
        }
        let value = sum.value()?;
        // Differentiate BOTH phasors in Re(flow^H*effort)/2. Retain the two
        // distinct cross terms even for mutual or nonreciprocal impedances.
        power.half_inner(flow[a], primal.effort[a], 1.0)?;
        power.half_inner(
            angular.flow_derivative[a],
            parameter_port.effort_derivative[a],
            1.0,
        )?;
        power.half_inner(
            parameter_port.flow_derivative[a],
            angular.effort_derivative[a],
            1.0,
        )?;
        power.half_inner(primal.flow[a], value, 1.0)?;
        effort.push(value);
    }
    let power = power.value()?;
    if !power.is_finite() {
        return Err(invalid("mixed physical port work is not representable"));
    }
    poll(cx)?;
    Ok(Some(HarmonicPortFrequencyParameterResponse {
        flow_mixed_derivative: flow,
        effort_mixed_derivative: effort,
        average_absorbed_power_mixed_derivative: power,
    }))
}
