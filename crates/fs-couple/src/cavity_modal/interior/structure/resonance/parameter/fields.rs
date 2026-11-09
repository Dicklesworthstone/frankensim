//! Physical fields and work of the original complete structural/acoustic jet.
//! This recovery owns no factorization, discretization or resonance search.

use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicPowerParameterJet, HarmonicPowerSample, HarmonicReceiverParameterJet,
    HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet,
};
use fs_modal::harmonic::{
    HarmonicFrequencyParameterResponse, HarmonicResponse, HarmonicResponseError,
    harmonic_mean_power_parameter_jet_with_cancel,
};

type Error = HarmonicResponseError;

/// One derivative order of every original physical field. The containing
/// response states whether this is omega, omega-omega, theta or omega-theta.
#[derive(Clone, Debug)]
pub struct StructuralInteriorDerivativeFields {
    pub structural_displacement: Vec<C64>,
    pub structural_velocity: Vec<C64>,
    pub pressure_potential_pa_s: Vec<C64>,
    pub pressure_pa: Vec<C64>,
    pub fluid_force: Vec<C64>,
    pub displaced_flow_m3_s: Vec<C64>,
}

/// Select an original field without changing its units or coordinate ordering.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StructuralInteriorResponseField {
    StructuralDisplacement,
    StructuralVelocity,
    PressurePotential,
    Pressure,
    FluidForce,
    DisplacedFlow,
}

impl StructuralInteriorDerivativeFields {
    pub(super) fn field(&self, field: StructuralInteriorResponseField) -> &[C64] {
        match field {
            StructuralInteriorResponseField::StructuralDisplacement => {
                &self.structural_displacement
            }
            StructuralInteriorResponseField::StructuralVelocity => &self.structural_velocity,
            StructuralInteriorResponseField::PressurePotential => &self.pressure_potential_pa_s,
            StructuralInteriorResponseField::Pressure => &self.pressure_pa,
            StructuralInteriorResponseField::FluidForce => &self.fluid_force,
            StructuralInteriorResponseField::DisplacedFlow => &self.displaced_flow_m3_s,
        }
    }
}

/// Independent original work channels, including all source and coefficient
/// product rules. Values are cycle averages of peak phasors in watts; frequency
/// derivatives are per rad/s and parameter derivatives per authored theta unit.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorPowerParameterJet {
    pub structural_input: HarmonicPowerParameterJet,
    pub acoustic_input: HarmonicPowerParameterJet,
    /// Signed work from the interior fluid into the structure.
    pub fluid_to_structure: HarmonicPowerParameterJet,
    /// Symmetric structural damping only; gyroscopic work is exactly excluded.
    pub structural_dissipation: HarmonicPowerParameterJet,
    /// Relative-wall heat, including its directly authored parameter change.
    pub wall_dissipation: HarmonicPowerParameterJet,
    /// Whole-pencil input minus independently differentiated native heat.
    pub original_balance: HarmonicPowerParameterJet,
}

impl StructuralInteriorPowerParameterJet {
    /// Structural input plus reciprocal exchange minus structural heat.
    /// Cancellation is accumulated before recovering the physical range.
    pub fn structural_balance(&self) -> Result<HarmonicPowerParameterJet, Error> {
        signed_power_sum(
            [
                self.structural_input,
                self.fluid_to_structure,
                self.structural_dissipation,
            ],
            [1.0, 1.0, -1.0],
        )
    }

    /// Acoustic input minus the same exchange and independently observed wall heat.
    pub fn acoustic_balance(&self) -> Result<HarmonicPowerParameterJet, Error> {
        signed_power_sum(
            [
                self.acoustic_input,
                self.fluid_to_structure,
                self.wall_dissipation,
            ],
            [1.0, -1.0, -1.0],
        )
    }
}

fn signed_power_sum(
    powers: [HarmonicPowerParameterJet; 3],
    signs: [f64; 3],
) -> Result<HarmonicPowerParameterJet, Error> {
    let values = powers.map(|p| C64::from_re(p.sample.average_power_w));
    let angular = powers.map(|p| C64::from_re(p.sample.angular_derivative));
    let second = powers.map(|p| C64::from_re(p.angular_second_derivative));
    let parameter = powers.map(|p| C64::from_re(p.parameter_derivative));
    let mixed = powers.map(|p| C64::from_re(p.angular_parameter_derivative));
    let sum = HarmonicTransferParameterJet::project_real(
        &signs,
        HarmonicTransferFieldParameterJet {
            values: &values,
            angular_derivative: &angular,
            angular_second_derivative: &second,
            parameter_derivative: &parameter,
            angular_parameter_derivative: &mixed,
        },
    )?;
    Ok(HarmonicPowerParameterJet {
        sample: HarmonicPowerSample {
            average_power_w: sum.sample.amplitude.re,
            angular_derivative: sum.sample.angular_derivative.re,
        },
        angular_second_derivative: sum.angular_second_derivative.re,
        parameter_derivative: sum.parameter_derivative.re,
        angular_parameter_derivative: sum.angular_parameter_derivative.re,
    })
}

pub(super) struct Recovered {
    pub(super) response: HarmonicResponse,
    pub(super) angular_derivative: StructuralInteriorDerivativeFields,
    pub(super) angular_second_derivative: StructuralInteriorDerivativeFields,
    pub(super) parameter_derivative: StructuralInteriorDerivativeFields,
    pub(super) angular_parameter_derivative: StructuralInteriorDerivativeFields,
    pub(super) fluid_force: Vec<C64>,
    pub(super) displaced_flow: Vec<C64>,
    pub(super) power: StructuralInteriorPowerParameterJet,
    pub(super) relative_derivative_residual: f64,
    pub(super) relative_second_derivative_residual: f64,
    pub(super) relative_parameter_derivative_residual: f64,
    pub(super) relative_mixed_derivative_residual: f64,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

fn field(values: &[Vec<C64>; 5]) -> HarmonicTransferFieldParameterJet<'_> {
    HarmonicTransferFieldParameterJet {
        values: &values[0],
        angular_derivative: &values[1],
        angular_second_derivative: &values[2],
        parameter_derivative: &values[3],
        angular_parameter_derivative: &values[4],
    }
}

fn slice(
    values: HarmonicTransferFieldParameterJet<'_>,
    range: core::ops::Range<usize>,
) -> HarmonicTransferFieldParameterJet<'_> {
    HarmonicTransferFieldParameterJet {
        values: &values.values[range.clone()],
        angular_derivative: &values.angular_derivative[range.clone()],
        angular_second_derivative: &values.angular_second_derivative[range.clone()],
        parameter_derivative: &values.parameter_derivative[range.clone()],
        angular_parameter_derivative: &values.angular_parameter_derivative[range],
    }
}

fn complete(values: HarmonicTransferFieldParameterJet<'_>, n: usize) -> bool {
    [
        values.values,
        values.angular_derivative,
        values.angular_second_derivative,
        values.parameter_derivative,
        values.angular_parameter_derivative,
    ]
    .iter()
    .all(|v| v.len() == n)
}

/// A real physical map and its theta derivative act on the original field.
/// Both product-rule terms use the native normalized coherent projection.
/// Only two original rows and five output vectors are allocated; no full
/// complex derivative matrix is materialized.
fn map_field(
    values: HarmonicTransferFieldParameterJet<'_>,
    rows: usize,
    mut coefficient: impl FnMut(usize, usize) -> (f64, f64),
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Vec<C64>; 5], Error> {
    let n = values.values.len();
    let mut weights = vec![C64::ZERO; n];
    let mut direction = vec![C64::ZERO; n];
    let mut result: [Vec<C64>; 5] = core::array::from_fn(|_| Vec::with_capacity(rows));
    for i in 0..rows {
        poll(cancelled)?;
        for j in 0..n {
            let (a, da) = coefficient(i, j);
            weights[j] = C64::from_re(a);
            direction[j] = C64::from_re(da);
        }
        let value = HarmonicTransferParameterJet::project_complex_with_derivatives(
            HarmonicReceiverParameterJet {
                weights: &weights,
                angular_derivative: None,
                angular_second_derivative: None,
                parameter_derivative: Some(&direction),
                angular_parameter_derivative: None,
            },
            values,
        )?;
        for (out, value) in result.iter_mut().zip([
            value.sample.amplitude,
            value.sample.angular_derivative,
            value.angular_second_derivative,
            value.parameter_derivative,
            value.angular_parameter_derivative,
        ]) {
            out.push(value);
        }
    }
    poll(cancelled)?;
    Ok(result)
}

fn symmetric(matrix: &[f64], n: usize, i: usize, j: usize) -> f64 {
    let a = matrix[i * n + j];
    let b = matrix[j * n + i];
    if a == b { a } else { 0.5 * a + 0.5 * b }
}

fn heat(
    values: HarmonicTransferFieldParameterJet<'_>,
    damping: &[f64],
    direction: &[f64],
    n: usize,
    offset: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerParameterJet, Error> {
    let effort = map_field(
        values,
        values.values.len(),
        |i, j| {
            (
                symmetric(damping, n, offset + i, offset + j),
                symmetric(direction, n, offset + i, offset + j),
            )
        },
        cancelled,
    )?;
    harmonic_mean_power_parameter_jet_with_cancel(values, field(&effort), cancelled)
}

fn split(
    ns: usize,
    mut displacement: Vec<C64>,
    mut velocity: Vec<C64>,
    force: Vec<C64>,
    flow: Vec<C64>,
) -> StructuralInteriorDerivativeFields {
    let potential = displacement.split_off(ns);
    let pressure = velocity.split_off(ns);
    StructuralInteriorDerivativeFields {
        structural_displacement: displacement,
        structural_velocity: velocity,
        pressure_potential_pa_s: potential,
        pressure_pa: pressure,
        fluid_force: force,
        displaced_flow_m3_s: flow,
    }
}

/// Recover actual pressure, reciprocal force/flow and independently measured
/// physical work from one complete original solved jet. The acoustic coordinate
/// is pressure potential, so its native velocity is physical pressure and
/// already includes -i(psi + omega psi_omega) and every higher product rule.
#[allow(clippy::too_many_arguments)]
pub(super) fn recover(
    ns: usize,
    coupling: &[f64],
    damping: &[f64],
    damping_parameter: &[f64],
    coupling_parameter: &[f64],
    sources: HarmonicTransferFieldParameterJet<'_>,
    native: HarmonicFrequencyParameterResponse,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Recovered, Error> {
    poll(cancelled)?;
    let n = sources.values.len();
    if ns == 0 || ns >= n {
        return Err(Error::Invalid(
            "nonempty original structural and pressure coordinates required",
        ));
    }
    let na = n - ns;
    if n.checked_mul(n) != Some(damping.len())
        || damping_parameter.len() != damping.len()
        || ns.checked_mul(na) != Some(coupling.len())
        || coupling_parameter.len() != coupling.len()
        || !complete(sources, n)
        || !complete(native.displacement_field_parameter_jet(), n)
        || !complete(native.velocity_field_parameter_jet(), n)
    {
        return Err(Error::Invalid(
            "complete original physical field and coefficient jets required",
        ));
    }
    if native
        .frequency_response
        .frequency_response
        .response
        .port_load
        .is_some()
    {
        return Err(Error::Invalid(
            "structural-interior parameter recovery requires an unloaded whole pencil",
        ));
    }
    let velocity = native.velocity_field_parameter_jet();
    let structural_velocity = slice(velocity, 0..ns);
    let pressure = slice(velocity, ns..n);
    // Both uses of B and B_theta are retained. The reciprocal exchange remains
    // signed work and is never silently assigned to material or wall heat.
    let force = map_field(
        pressure,
        ns,
        |s, a| (coupling[a * ns + s], coupling_parameter[a * ns + s]),
        cancelled,
    )?;
    let flow = map_field(
        structural_velocity,
        na,
        |a, s| (coupling[a * ns + s], coupling_parameter[a * ns + s]),
        cancelled,
    )?;
    let structural_input = harmonic_mean_power_parameter_jet_with_cancel(
        structural_velocity,
        slice(sources, 0..ns),
        &mut *cancelled,
    )?;
    let acoustic_input = harmonic_mean_power_parameter_jet_with_cancel(
        slice(sources, ns..n),
        pressure,
        &mut *cancelled,
    )?;
    let fluid_to_structure = harmonic_mean_power_parameter_jet_with_cancel(
        structural_velocity,
        field(&force),
        &mut *cancelled,
    )?;
    let structural_dissipation = heat(
        structural_velocity,
        damping,
        damping_parameter,
        n,
        0,
        cancelled,
    )?;
    let wall_dissipation = heat(pressure, damping, damping_parameter, n, ns, cancelled)?;
    let first = &native.frequency_response.frequency_response;
    let power = StructuralInteriorPowerParameterJet {
        structural_input,
        acoustic_input,
        fluid_to_structure,
        structural_dissipation,
        wall_dissipation,
        original_balance: HarmonicPowerParameterJet {
            sample: HarmonicPowerSample {
                average_power_w: first.response.power_balance_residual_w,
                angular_derivative: first.power_balance_residual_derivative,
            },
            angular_second_derivative: native
                .frequency_response
                .power_balance_residual_second_derivative,
            parameter_derivative: native.parameter_response.power_balance_residual_derivative,
            angular_parameter_derivative: native.power_balance_residual_mixed_derivative,
        },
    };
    let relative_derivative_residual = first.relative_derivative_residual;
    let relative_second_derivative_residual = native
        .frequency_response
        .relative_second_derivative_residual;
    let relative_parameter_derivative_residual =
        native.parameter_response.relative_derivative_residual;
    let relative_mixed_derivative_residual = native.relative_mixed_derivative_residual;
    let [fluid_force, fw, fww, ft, fwt] = force;
    let [displaced_flow, uw, uww, ut, uwt] = flow;
    let first = native.frequency_response.frequency_response;
    let angular_derivative = split(
        ns,
        first.displacement_derivative,
        first.velocity_derivative,
        fw,
        uw,
    );
    let angular_second_derivative = split(
        ns,
        native.frequency_response.displacement_second_derivative,
        native.frequency_response.velocity_second_derivative,
        fww,
        uww,
    );
    let parameter_derivative = split(
        ns,
        native.parameter_response.displacement_derivative,
        native.parameter_response.velocity_derivative,
        ft,
        ut,
    );
    let angular_parameter_derivative = split(
        ns,
        native.displacement_mixed_derivative,
        native.velocity_mixed_derivative,
        fwt,
        uwt,
    );
    poll(cancelled)?;
    Ok(Recovered {
        response: first.response,
        angular_derivative,
        angular_second_derivative,
        parameter_derivative,
        angular_parameter_derivative,
        fluid_force,
        displaced_flow,
        power,
        relative_derivative_residual,
        relative_second_derivative_residual,
        relative_parameter_derivative_residual,
        relative_mixed_derivative_residual,
    })
}

#[cfg(test)]
#[path = "fields_tests.rs"]
mod tests;
