//! Original-port amplitudes and independently normalized wave-flux jets.

use super::{
    C64, HarmonicError, HarmonicFrequencyParameterResponse, HarmonicParameterDirection,
    ScatteringModel, SourceJet, Sum, invalid, physical, poll,
};

/// A complex value and local derivatives. Omega is angular frequency in rad/s;
/// theta is the caller's authored real physical parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringComplexJet {
    pub value: C64,
    pub frequency: C64,
    pub second_frequency: C64,
    pub parameter: C64,
    pub mixed: C64,
}

/// A power in watts and its signed local derivatives per rad/s and per
/// parameter unit. second_frequency is d²/domega²; mixed is d²/(domega dtheta).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringPowerJet {
    pub value: f64,
    pub frequency: f64,
    pub second_frequency: f64,
    pub parameter: f64,
    pub mixed: f64,
}

/// One selected original boundary, including the derivative of its watt
/// normalization. Outgoing amplitude derivatives remain defined at zero
/// reflection; no logarithm or division by the wave amplitude is used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringWaveFrequencyParameterJet {
    pub original_port: usize,
    /// Positive output/input reference, constant in omega.
    pub reference: f64,
    /// Authored local dZ/dtheta; no positivity is required of the slope.
    pub reference_derivative: f64,
    pub incident: ScatteringComplexJet,
    pub outgoing: ScatteringComplexJet,
    /// Independent derivatives of |a|²/(2Z).
    pub incident_power: ScatteringPowerJet,
    /// Independent derivatives of |b|²/(2Z).
    pub outgoing_power: ScatteringPowerJet,
    /// Original conjugate work Re(conj(u)*y)/2 and its derivatives.
    pub physical_power: ScatteringPowerJet,
    /// Incident minus outgoing minus physical work, independently evaluated.
    pub power_residual: ScatteringPowerJet,
}

impl ScatteringComplexJet {
    pub(super) fn from_source(source: SourceJet<'_>, port: usize) -> Self {
        Self {
            value: source.input[port],
            frequency: at(source.frequency, port),
            second_frequency: at(source.second_frequency, port),
            parameter: at(source.parameter, port),
            mixed: at(source.mixed, port),
        }
    }
}

impl ScatteringPowerJet {
    fn from_values(v: [f64; 5]) -> Self {
        Self {
            value: v[0],
            frequency: v[1],
            second_frequency: v[2],
            parameter: v[3],
            mixed: v[4],
        }
    }
    fn values(self) -> [f64; 5] {
        [
            self.value,
            self.frequency,
            self.second_frequency,
            self.parameter,
            self.mixed,
        ]
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn observe_waves(
    owner: &ScatteringModel,
    source: SourceJet<'_>,
    actual: SourceJet<'_>,
    direction: HarmonicParameterDirection<'_>,
    reference_derivative: Option<&[f64]>,
    physical: &HarmonicFrequencyParameterResponse,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<ScatteringWaveFrequencyParameterJet>, ScatteringPowerJet), HarmonicError> {
    let curvature = &physical.frequency_response;
    let first = &curvature.frequency_response;
    let original = &first.response;
    let parameter = &physical.parameter_response;
    let model = &owner.closed;
    let keep = |p: usize| !owner.selected[p];
    let mut supply = [
        crate::harmonic::physical_power_sum(
            model,
            actual.input,
            &original.effort,
            keep,
            cancelled,
        )?,
        crate::harmonic::physical_power_tangent_sum(
            model,
            actual.input,
            actual.frequency,
            &original.effort,
            &first.effort_derivative,
            None,
            keep,
            cancelled,
        )?,
        physical::physical_power_curvature_sum(
            model,
            actual,
            &original.effort,
            &first.effort_derivative,
            &curvature.effort_second_derivative,
            keep,
            cancelled,
        )?,
        crate::harmonic::physical_power_tangent_sum(
            model,
            actual.input,
            actual.parameter,
            &original.effort,
            &parameter.effort_derivative,
            direction.ports,
            keep,
            cancelled,
        )?,
        physical::physical_power_mixed_sum(
            model,
            actual,
            direction,
            &original.effort,
            &first.effort_derivative,
            &parameter.effort_derivative,
            &physical.effort_mixed_derivative,
            keep,
            cancelled,
        )?,
    ];
    let mut waves = Vec::with_capacity(owner.ports.len());
    for (index, &(p, z)) in owner.ports.iter().enumerate() {
        poll(cancelled)?;
        let dz = reference_derivative.map_or(0.0, |v| v[index]);
        let incident = ScatteringComplexJet::from_source(source, p);
        let outgoing = ScatteringComplexJet {
            value: original.output[p] - incident.value,
            frequency: first.output_derivative[p] - incident.frequency,
            second_frequency: curvature.output_second_derivative[p] - incident.second_frequency,
            parameter: parameter.output_derivative[p] - incident.parameter,
            mixed: physical.output_mixed_derivative[p] - incident.mixed,
        };
        let incident_power = wave_power(incident, z, dz)?;
        let outgoing_power = wave_power(outgoing, z, dz)?;
        let physical_power = ScatteringPowerJet::from_values([
            original.port_power_w[p],
            first.port_power_derivative[p],
            curvature.port_power_second_derivative[p],
            parameter.port_power_derivative[p],
            physical.port_power_mixed_derivative[p],
        ]);
        let mut residual = [0.0; 5];
        let (inc, out, work) = (
            incident_power.values(),
            outgoing_power.values(),
            physical_power.values(),
        );
        for order in 0..5 {
            let mut difference = Sum::default();
            difference
                .product(inc[order], 1.0)
                .and_then(|()| difference.product(out[order], -1.0))
                .and_then(|()| difference.product(work[order], -1.0))
                .and_then(|()| supply[order].product(inc[order], 1.0))
                .and_then(|()| supply[order].product(out[order], -1.0))
                .ok_or_else(|| invalid("independent loaded wave jet balance overflows"))?;
            residual[order] = power_value(difference)?;
        }
        waves.push(ScatteringWaveFrequencyParameterJet {
            original_port: p,
            reference: z,
            reference_derivative: dz,
            incident,
            outgoing,
            incident_power,
            outgoing_power,
            physical_power,
            power_residual: ScatteringPowerJet::from_values(residual),
        });
    }
    let loss = [
        original.dissipated_power_w,
        first.dissipated_power_derivative,
        curvature.dissipated_power_second_derivative,
        parameter.dissipated_power_derivative,
        physical.dissipated_power_mixed_derivative,
    ];
    let mut balance = [0.0; 5];
    for order in 0..5 {
        let mut residual = Sum::default();
        residual
            .product(loss[order], 1.0)
            .and_then(|()| residual.add_scaled(supply[order], -1.0))
            .ok_or_else(|| invalid("original loss minus wave source jet overflows"))?;
        balance[order] = power_value(residual)?;
    }
    Ok((waves, ScatteringPowerJet::from_values(balance)))
}

pub(in crate::harmonic::scattering) fn wave_power(
    a: ScatteringComplexJet,
    z: f64,
    dz: f64,
) -> Result<ScatteringPowerJet, HarmonicError> {
    let mut value = Sum::default();
    let frequency = wave_frequency_sum(a.value, a.frequency, z)?;
    let mut second = Sum::default();
    let mut parameter = Sum::default();
    let mut mixed = Sum::default();
    add_pair(&mut value, a.value, a.value, z)?;
    for (left, right) in [
        (a.second_frequency, a.value),
        (a.frequency, a.frequency),
        (a.frequency, a.frequency),
        (a.value, a.second_frequency),
    ] {
        add_pair(&mut second, left, right, z)?;
    }
    for (left, right) in [(a.parameter, a.value), (a.value, a.parameter)] {
        add_pair(&mut parameter, left, right, z)?;
    }
    for (left, right) in [
        (a.parameter, a.frequency),
        (a.value, a.mixed),
        (a.mixed, a.value),
        (a.frequency, a.parameter),
    ] {
        add_pair(&mut mixed, left, right, z)?;
    }
    parameter
        .add_scaled(value, -dz / z)
        .and_then(|()| mixed.add_scaled(frequency, -dz / z))
        .ok_or_else(|| invalid("reference-normalized wave power jet overflows"))?;
    Ok(ScatteringPowerJet::from_values([
        power_value(value)?,
        power_value(frequency)?,
        power_value(second)?,
        power_value(parameter)?,
        power_value(mixed)?,
    ]))
}

fn add_pair(total: &mut Sum, left: C64, right: C64, z: f64) -> Result<(), HarmonicError> {
    total
        .half_product(left.re / z, right.re)
        .and_then(|()| total.half_product(left.im / z, right.im))
        .ok_or_else(|| invalid("wave power jet is not representable"))
}

/// One frequency contraction shared by first-order and complete response jets.
/// Both conjugate terms retain the same normalized summation near a stationary
/// wave power; multiplying a separately rounded half-product by two would not
/// reproduce that small signed slope.
pub(in crate::harmonic::scattering) fn wave_power_frequency(
    amplitude: C64,
    derivative: C64,
    reference: f64,
) -> Result<f64, HarmonicError> {
    power_value(wave_frequency_sum(amplitude, derivative, reference)?)
}

fn wave_frequency_sum(amplitude: C64, derivative: C64, reference: f64) -> Result<Sum, HarmonicError> {
    let mut sum = Sum::default();
    add_pair(&mut sum, amplitude, derivative, reference)?;
    add_pair(&mut sum, derivative, amplitude, reference)?;
    Ok(sum)
}

fn power_value(sum: Sum) -> Result<f64, HarmonicError> {
    sum.value()
        .ok_or_else(|| invalid("wave power jet sum overflows"))
}
fn at(values: Option<&[C64]>, port: usize) -> C64 {
    values.map_or(C64::ZERO, |v| v[port])
}
