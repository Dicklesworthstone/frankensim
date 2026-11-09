//! First-order original work accounts, without curvature state solves.

use super::*;
use fs_modal::harmonic::resonance::{HarmonicPowerSample, HarmonicTransferSample};
use fs_modal::harmonic::{
    HarmonicResponseError, harmonic_mean_power_derivative_with_cancel,
    harmonic_quadratic_power_derivative_with_cancel,
};

/// Independently observed original work/heat and frequency slopes [W, W/(rad/s)].
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationPowerFrequency {
    pub structural_input: HarmonicPowerSample,
    pub acoustic_input: HarmonicPowerSample,
    pub fluid_to_structure: HarmonicPowerSample,
    pub structural_dissipation: HarmonicPowerSample,
    pub wall_dissipation: HarmonicPowerSample,
    pub relaxation_dissipation: HarmonicPowerSample,
    pub structural_to_relaxation: HarmonicPowerSample,
    pub acoustic_to_relaxation: HarmonicPowerSample,
}

/// One original branch's independent heat and signed work frequency slopes.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationMemoryFrequency {
    pub dissipation: HarmonicPowerSample,
    pub structural_work: HarmonicPowerSample,
    pub acoustic_work: HarmonicPowerSample,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn observe(
    p: &Prepared<'_>,
    base: &StructuralInteriorRelaxationResponse,
    first: &State,
    force_first: &[C64],
    input: &[C64],
    input_first: &[C64],
    effort: &[C64],
    effort_first: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<
    (
        StructuralInteriorRelaxationPowerFrequency,
        Vec<StructuralInteriorRelaxationMemoryFrequency>,
    ),
    Error,
> {
    let (n, ns) = (p.source.coordinates, p.source.structural);
    let velocity_first = &first.velocity[..ns];
    let pressure_first = &first.velocity[ns..];
    let structural_input = work(
        base.average_structural_input_power_w,
        [&base.structural_velocity, velocity_first],
        [&input[..ns], &input_first[..ns]],
        cancelled,
    )?;
    let acoustic_input = work(
        base.average_acoustic_input_power_w,
        [&base.pressure_pa, pressure_first],
        [&input[ns..], &input_first[ns..]],
        cancelled,
    )?;
    let fluid_to_structure = work(
        base.average_fluid_to_structure_power_w,
        [&base.structural_velocity, velocity_first],
        [&base.fluid_force, force_first],
        cancelled,
    )?;
    let block = |offset: usize, size: usize| {
        let mut matrix = Vec::with_capacity(size * size);
        for row in offset..offset + size {
            matrix.extend_from_slice(&p.source.damping[row * n + offset..row * n + offset + size]);
        }
        matrix
    };
    let structural_dissipation = quadratic(
        &base.structural_velocity,
        velocity_first,
        &block(0, ns),
        cancelled,
    )?;
    let wall_dissipation = quadratic(
        &base.pressure_pa,
        pressure_first,
        &block(ns, n - ns),
        cancelled,
    )?;
    let mut memory = Vec::with_capacity(p.branches.len());
    for (index, branch) in p.branches.iter().enumerate() {
        poll(cancelled)?;
        let canonical = 2 * n + index;
        let dissipation = quadratic(
            &[effort[canonical]],
            &[effort_first[canonical]],
            &[branch.relaxation_time_s.recip()],
            cancelled,
        )?;
        let structural = HarmonicTransferSample::project_real(
            &branch.strain[..ns],
            &base.structural_velocity,
            velocity_first,
        )
        .map_err(observation)?;
        let acoustic = HarmonicTransferSample::project_real(
            &branch.strain[ns..],
            &base.pressure_pa,
            pressure_first,
        )
        .map_err(observation)?;
        let branch_base = &base.memory[index];
        let stress = [[branch_base.stress], [first.stress[index]]];
        memory.push(StructuralInteriorRelaxationMemoryFrequency {
            dissipation,
            structural_work: work(
                branch_base.average_structural_work_w,
                [&[structural.amplitude], &[structural.angular_derivative]],
                [&stress[0], &stress[1]],
                cancelled,
            )?,
            acoustic_work: work(
                branch_base.average_acoustic_work_w,
                [&[acoustic.amplitude], &[acoustic.angular_derivative]],
                [&stress[0], &stress[1]],
                cancelled,
            )?,
        });
    }
    let mut total =
        |select: fn(&StructuralInteriorRelaxationMemoryFrequency) -> HarmonicPowerSample| {
            poll(cancelled)?;
            sum(&memory.iter().map(select).collect::<Vec<_>>())
        };
    let power = StructuralInteriorRelaxationPowerFrequency {
        structural_input,
        acoustic_input,
        fluid_to_structure,
        structural_dissipation,
        wall_dissipation,
        relaxation_dissipation: total(|m| m.dissipation)?,
        structural_to_relaxation: total(|m| m.structural_work)?,
        acoustic_to_relaxation: total(|m| m.acoustic_work)?,
    };
    poll(cancelled)?;
    Ok((power, memory))
}

fn quadratic(
    values: &[C64],
    first: &[C64],
    matrix: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerSample, Error> {
    harmonic_quadratic_power_derivative_with_cancel(values, first, matrix, cancelled)
        .map_err(observation)
}
fn work(
    value: f64,
    flow: [&[C64]; 2],
    effort: [&[C64]; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerSample, Error> {
    Ok(HarmonicPowerSample {
        average_power_w: value,
        angular_derivative: harmonic_mean_power_derivative_with_cancel(
            flow[0], flow[1], effort[0], effort[1], cancelled,
        )
        .map_err(observation)?,
    })
}
pub(in super::super) fn sum(values: &[HarmonicPowerSample]) -> Result<HarmonicPowerSample, Error> {
    if values.is_empty() {
        return Ok(HarmonicPowerSample {
            average_power_w: 0.0,
            angular_derivative: 0.0,
        });
    }
    let amplitude: Vec<_> = values
        .iter()
        .map(|p| C64::from_re(p.average_power_w))
        .collect();
    let first: Vec<_> = values
        .iter()
        .map(|p| C64::from_re(p.angular_derivative))
        .collect();
    let total = HarmonicTransferSample::project_real(&vec![1.0; values.len()], &amplitude, &first)
        .map_err(observation)?;
    Ok(HarmonicPowerSample {
        average_power_w: total.amplitude.re,
        angular_derivative: total.angular_derivative.re,
    })
}
fn observation(error: HarmonicResponseError) -> Error {
    match error {
        HarmonicResponseError::Cancelled => Error::Cancelled,
        _ => Error::Invalid("original physical power frequency observation is unrepresentable"),
    }
}
