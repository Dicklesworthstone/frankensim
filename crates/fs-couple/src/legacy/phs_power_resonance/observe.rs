//! Select existing independent native accounts without recomputing physics.

use super::{
    HarmonicError, HarmonicPowerCurvatureSample, HarmonicPowerParameterJet, HarmonicPowerSample,
    PhsPowerObservation as Channel, invalid,
};
use fs_phs::harmonic::scattering::{
    ScatteringFrequencyCurvatureResponse, ScatteringFrequencyParameterResponse,
    ScatteringFrequencyResponse,
};
use fs_phs::harmonic::{
    HarmonicFrequencyCurvatureResponse, HarmonicFrequencyParameterResponse,
    HarmonicFrequencyResponse,
};

fn at(values: &[f64], p: usize) -> Result<f64, HarmonicError> {
    values
        .get(p)
        .copied()
        .ok_or_else(|| invalid("original power port is out of range"))
}
fn sample(value: f64, first: f64) -> HarmonicPowerSample {
    HarmonicPowerSample {
        average_power_w: value,
        angular_derivative: first,
    }
}

pub(super) fn original_sample(
    channel: Channel,
    r: &HarmonicFrequencyResponse,
) -> Result<HarmonicPowerSample, HarmonicError> {
    Ok(match channel {
        Channel::Dissipation => {
            sample(r.response.dissipated_power_w, r.dissipated_power_derivative)
        }
        Channel::Supplied => sample(r.response.supplied_power_w, r.supplied_power_derivative),
        Channel::Port(p) => sample(
            at(&r.response.port_power_w, p)?,
            at(&r.port_power_derivative, p)?,
        ),
        _ => {
            return Err(invalid(
                "selected wave power requires a loaded physical response",
            ));
        }
    })
}
pub(super) fn original_curvature(
    channel: Channel,
    r: &HarmonicFrequencyCurvatureResponse,
) -> Result<HarmonicPowerCurvatureSample, HarmonicError> {
    let sample = original_sample(channel, &r.frequency_response)?;
    let second = match channel {
        Channel::Dissipation => r.dissipated_power_second_derivative,
        Channel::Supplied => r.supplied_power_second_derivative,
        Channel::Port(p) => at(&r.port_power_second_derivative, p)?,
        _ => {
            return Err(invalid(
                "selected wave power requires a loaded physical response",
            ));
        }
    };
    Ok(HarmonicPowerCurvatureSample {
        sample,
        angular_second_derivative: second,
    })
}
pub(super) fn original_parameter(
    channel: Channel,
    r: &HarmonicFrequencyParameterResponse,
) -> Result<HarmonicPowerParameterJet, HarmonicError> {
    let curvature = original_curvature(channel, &r.frequency_response)?;
    let p = &r.parameter_response;
    let (parameter, mixed) = match channel {
        Channel::Dissipation => (
            p.dissipated_power_derivative,
            r.dissipated_power_mixed_derivative,
        ),
        Channel::Supplied => (
            p.supplied_power_derivative,
            r.supplied_power_mixed_derivative,
        ),
        Channel::Port(i) => (
            at(&p.port_power_derivative, i)?,
            at(&r.port_power_mixed_derivative, i)?,
        ),
        _ => {
            return Err(invalid(
                "selected wave power requires a loaded physical response",
            ));
        }
    };
    Ok(HarmonicPowerParameterJet {
        sample: curvature.sample,
        angular_second_derivative: curvature.angular_second_derivative,
        parameter_derivative: parameter,
        angular_parameter_derivative: mixed,
    })
}

pub(super) fn loaded_sample(
    channel: Channel,
    r: &ScatteringFrequencyResponse,
) -> Result<HarmonicPowerSample, HarmonicError> {
    Ok(match channel {
        Channel::Dissipation => sample(
            r.response.physical.dissipated_power_w,
            r.dissipated_power_derivative,
        ),
        Channel::Supplied => sample(
            r.response.physical.supplied_power_w,
            r.supplied_power_derivative,
        ),
        Channel::Port(p) => sample(
            at(&r.response.physical.port_power_w, p)?,
            at(&r.port_power_derivative, p)?,
        ),
        Channel::IncidentWave(p) | Channel::OutgoingWave(p) => {
            let wave = r
                .response
                .waves
                .get(p)
                .ok_or_else(|| invalid("selected wave port is out of range"))?;
            let first = &r.wave_derivatives[p];
            if matches!(channel, Channel::IncidentWave(_)) {
                sample(wave.incident_power_w, first.incident_power_w)
            } else {
                sample(wave.outgoing_power_w, first.outgoing_power_w)
            }
        }
    })
}
pub(super) fn loaded_curvature(
    channel: Channel,
    r: &ScatteringFrequencyCurvatureResponse,
) -> Result<HarmonicPowerCurvatureSample, HarmonicError> {
    match channel {
        Channel::IncidentWave(p) | Channel::OutgoingWave(p) => {
            let wave = r
                .waves
                .get(p)
                .ok_or_else(|| invalid("selected wave port is out of range"))?;
            let power = if matches!(channel, Channel::IncidentWave(_)) {
                wave.incident_power
            } else {
                wave.outgoing_power
            };
            Ok(HarmonicPowerCurvatureSample {
                sample: sample(power.value, power.frequency),
                angular_second_derivative: power.second_frequency,
            })
        }
        _ => original_curvature(channel, &r.physical),
    }
}
pub(super) fn loaded_parameter(
    channel: Channel,
    r: &ScatteringFrequencyParameterResponse,
) -> Result<HarmonicPowerParameterJet, HarmonicError> {
    match channel {
        Channel::IncidentWave(p) | Channel::OutgoingWave(p) => {
            let wave = r
                .waves
                .get(p)
                .ok_or_else(|| invalid("selected wave port is out of range"))?;
            let power = if matches!(channel, Channel::IncidentWave(_)) {
                wave.incident_power
            } else {
                wave.outgoing_power
            };
            Ok(HarmonicPowerParameterJet {
                sample: sample(power.value, power.frequency),
                angular_second_derivative: power.second_frequency,
                parameter_derivative: power.parameter,
                angular_parameter_derivative: power.mixed,
            })
        }
        _ => original_parameter(channel, &r.physical),
    }
}
