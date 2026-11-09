//! Original material-memory, pressure and physical work curvature.

use super::*;
use fs_modal::harmonic::resonance::{
    HarmonicPowerCurvatureSample, HarmonicPowerSample, HarmonicTransferFieldParameterJet,
    HarmonicTransferParameterJet, HarmonicTransferSample,
};
use fs_modal::harmonic::{
    HarmonicResponseError, harmonic_mean_power_curvature_with_cancel,
    harmonic_mean_power_derivative_with_cancel, harmonic_quadratic_power_curvature_with_cancel,
};

/// Independent physical power channels and their two frequency derivatives.
/// Every entry uses original peak phasors, W, and derivatives per rad/s.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationPowerCurvature {
    pub structural_input: HarmonicPowerCurvatureSample,
    pub acoustic_input: HarmonicPowerCurvatureSample,
    pub fluid_to_structure: HarmonicPowerCurvatureSample,
    pub structural_dissipation: HarmonicPowerCurvatureSample,
    pub wall_dissipation: HarmonicPowerCurvatureSample,
    pub relaxation_dissipation: HarmonicPowerCurvatureSample,
    pub structural_to_relaxation: HarmonicPowerCurvatureSample,
    pub acoustic_to_relaxation: HarmonicPowerCurvatureSample,
}
impl StructuralInteriorRelaxationPowerCurvature {
    /// Structural source plus fluid exchange minus structural and memory work.
    pub fn structural_balance(&self) -> Result<HarmonicPowerCurvatureSample, Error> {
        sum(
            &[
                self.structural_input,
                self.fluid_to_structure,
                self.structural_dissipation,
                self.structural_to_relaxation,
            ],
            &[1.0, 1.0, -1.0, -1.0],
        )
    }
    /// Acoustic source minus the same interface exchange, wall and memory work.
    pub fn acoustic_balance(&self) -> Result<HarmonicPowerCurvatureSample, Error> {
        sum(
            &[
                self.acoustic_input,
                self.fluid_to_structure,
                self.wall_dissipation,
                self.acoustic_to_relaxation,
            ],
            &[1.0, -1.0, -1.0, -1.0],
        )
    }
    /// Independently observed branch work minus original branch heat.
    pub fn relaxation_balance(&self) -> Result<HarmonicPowerCurvatureSample, Error> {
        sum(
            &[
                self.structural_to_relaxation,
                self.acoustic_to_relaxation,
                self.relaxation_dissipation,
            ],
            &[1.0, 1.0, -1.0],
        )
    }
    /// Both physical source channels minus all independently observed heat.
    pub fn total_balance(&self) -> Result<HarmonicPowerCurvatureSample, Error> {
        sum(
            &[
                self.structural_input,
                self.acoustic_input,
                self.structural_dissipation,
                self.wall_dissipation,
                self.relaxation_dissipation,
            ],
            &[1.0, 1.0, -1.0, -1.0, -1.0],
        )
    }
}

/// One original Maxwell arm, retaining its separate storage and work accounts.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationMemoryCurvature {
    /// Original spring storage [J], including its first two frequency derivatives.
    pub mean_stored_energy_j: f64,
    pub mean_stored_energy_derivative: f64,
    pub mean_stored_energy_second_derivative: f64,
    pub dissipation: HarmonicPowerCurvatureSample,
    pub structural_work: HarmonicPowerCurvatureSample,
    pub acoustic_work: HarmonicPowerCurvatureSample,
}

/// Full original causal-memory response, tangent and frequency curvature.
#[derive(Clone, Debug)]
pub struct StructuralInteriorRelaxationCurvatureResponse {
    pub frequency_response: StructuralInteriorRelaxationFrequencyResponse,
    pub structural_displacement_second_derivative: Vec<C64>,
    pub structural_velocity_second_derivative: Vec<C64>,
    pub pressure_potential_second_derivative: Vec<C64>,
    pub pressure_second_derivative: Vec<C64>,
    pub fluid_force_second_derivative: Vec<C64>,
    pub displaced_flow_second_derivative: Vec<C64>,
    pub relaxation_stress_second_derivative: Vec<C64>,
    pub memory: Vec<StructuralInteriorRelaxationMemoryCurvature>,
    /// Derivatives of complete mean storage, including every material spring.
    pub mean_stored_energy_derivative: f64,
    pub mean_stored_energy_second_derivative: f64,
    pub power: StructuralInteriorRelaxationPowerCurvature,
    /// Maximum native and original second-equation backward error.
    pub relative_second_derivative_residual: f64,
    pub original_relative_second_derivative_residual: f64,
}

impl Factor<'_, '_> {
    pub(in super::super) fn solve_with_frequency_curvature(
        &self,
        structural: &[C64],
        acoustic: &[C64],
        source_derivative: Option<(&[C64], &[C64])>,
        source_second_derivative: Option<(&[C64], &[C64])>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StructuralInteriorRelaxationCurvatureResponse, Error> {
        poll(cancelled)?;
        let p = self.prepared;
        let (n, ns) = (p.source.coordinates, p.source.structural);
        // The existing complete preparation bound includes these O(n^2+b*n)
        // physical observations and retained arrays. Native solve limits remain
        // separately enforced by the original full-state harmonic owner.
        admit(n, p.branches.len(), p.limits)?;
        let input = p.load(structural, acoustic)?;
        let first_input = source_derivative.map(|(s, a)| p.load(s, a)).transpose()?;
        let second_input = source_second_derivative
            .map(|(s, a)| p.load(s, a))
            .transpose()?;
        let zero = vec![C64::ZERO; n];
        let f1 = first_input.as_deref().unwrap_or(&zero);
        let f2 = second_input.as_deref().unwrap_or(&zero);
        let native = self.harmonic.response_with_frequency_curvature(
            &input,
            first_input.as_deref(),
            second_input.as_deref(),
            &mut *cancelled,
        )?;
        let first = p.state(
            &native.frequency_response.state_derivative,
            &native.frequency_response.effort_derivative,
            &native.frequency_response.output_derivative,
            cancelled,
        )?;
        let mut second = p.state(
            &native.state_second_derivative,
            &native.effort_second_derivative,
            &native.output_second_derivative,
            cancelled,
        )?;
        let original_second_residual =
            p.residual(self.omega, &second, f2, Some((&first, 2.0)), cancelled)?;
        let mut memory = Vec::with_capacity(p.branches.len());
        for (j, branch) in p.branches.iter().enumerate() {
            poll(cancelled)?;
            // The original canonical branch effort is -sigma/sqrt(k).
            // Keep its independent positive coefficient inside the normalized
            // work contraction; never multiply fields by a huge resistance.
            let index = 2 * n + j;
            let e = [native.frequency_response.response.effort[index]];
            let e1 = [native.frequency_response.effort_derivative[index]];
            let e2 = [native.effort_second_derivative[index]];
            let stored = quadratic(&e, &e1, &e2, &[0.5], cancelled)?;
            let heat = quadratic(&e, &e1, &e2, &[branch.relaxation_time_s.recip()], cancelled)?;
            memory.push(StructuralInteriorRelaxationMemoryCurvature {
                mean_stored_energy_j: stored.sample.average_power_w,
                mean_stored_energy_derivative: stored.sample.angular_derivative,
                mean_stored_energy_second_derivative: stored.angular_second_derivative,
                dissipation: heat,
                structural_work: zero_power(),
                acoustic_work: zero_power(),
            });
        }
        let native_second_residual = native.relative_second_derivative_residual;
        let mean_stored_energy_derivative = native.frequency_response.mean_stored_energy_derivative;
        let mean_stored_energy_second_derivative = native.mean_stored_energy_second_derivative;
        let frequency_response = self.recover_frequency(
            native.frequency_response,
            &input,
            first_input.as_deref(),
            cancelled,
        )?;
        let potential_second = second.q.split_off(ns);
        let pressure_second = second.velocity.split_off(ns);
        let (force_second, flow_second) = response::interface(
            p.source.coupling,
            &second.velocity,
            &pressure_second,
            cancelled,
        )?;
        let r = &frequency_response;
        let base = &r.response;
        let structural_input = work(
            base.average_structural_input_power_w,
            [
                &base.structural_velocity,
                &r.structural_velocity_derivative_per_rad_s,
                &second.velocity,
            ],
            [&input[..ns], &f1[..ns], &f2[..ns]],
            cancelled,
        )?;
        let acoustic_input = work(
            base.average_acoustic_input_power_w,
            [
                &base.pressure_pa,
                &r.pressure_derivative_pa_per_rad_s,
                &pressure_second,
            ],
            [&input[ns..], &f1[ns..], &f2[ns..]],
            cancelled,
        )?;
        let fluid_to_structure = work(
            base.average_fluid_to_structure_power_w,
            [
                &base.structural_velocity,
                &r.structural_velocity_derivative_per_rad_s,
                &second.velocity,
            ],
            [
                &base.fluid_force,
                &r.fluid_force_derivative_per_rad_s,
                &force_second,
            ],
            cancelled,
        )?;
        let block = |offset: usize, size: usize| {
            let mut values = Vec::with_capacity(size * size);
            for i in 0..size {
                values.extend_from_slice(
                    &p.source.damping[(offset + i) * n + offset..(offset + i) * n + offset + size],
                );
            }
            values
        };
        let structural_dissipation = quadratic(
            &base.structural_velocity,
            &r.structural_velocity_derivative_per_rad_s,
            &second.velocity,
            &block(0, ns),
            cancelled,
        )?;
        let wall_dissipation = quadratic(
            &base.pressure_pa,
            &r.pressure_derivative_pa_per_rad_s,
            &pressure_second,
            &block(ns, n - ns),
            cancelled,
        )?;
        for (j, (branch, account)) in p.branches.iter().zip(&mut memory).enumerate() {
            poll(cancelled)?;
            let structural_flow = projection(
                &branch.strain[..ns],
                [
                    &base.structural_velocity,
                    &r.structural_velocity_derivative_per_rad_s,
                    &second.velocity,
                ],
            )?;
            let acoustic_flow = projection(
                &branch.strain[ns..],
                [
                    &base.pressure_pa,
                    &r.pressure_derivative_pa_per_rad_s,
                    &pressure_second,
                ],
            )?;
            let stress = [
                [base.memory[j].stress],
                [r.relaxation_stress_derivative_per_rad_s[j]],
                [second.stress[j]],
            ];
            account.structural_work = work(
                base.memory[j].average_structural_work_w,
                [
                    &structural_flow[0],
                    &structural_flow[1],
                    &structural_flow[2],
                ],
                [&stress[0], &stress[1], &stress[2]],
                cancelled,
            )?;
            account.acoustic_work = work(
                base.memory[j].average_acoustic_work_w,
                [&acoustic_flow[0], &acoustic_flow[1], &acoustic_flow[2]],
                [&stress[0], &stress[1], &stress[2]],
                cancelled,
            )?;
        }
        let mut total = |select: fn(
            &StructuralInteriorRelaxationMemoryCurvature,
        ) -> HarmonicPowerCurvatureSample| {
            poll(cancelled)?;
            let values: Vec<_> = memory.iter().map(select).collect();
            sum(&values, &vec![1.0; values.len()])
        };
        let power = StructuralInteriorRelaxationPowerCurvature {
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
        Ok(StructuralInteriorRelaxationCurvatureResponse {
            frequency_response,
            structural_displacement_second_derivative: second.q,
            structural_velocity_second_derivative: second.velocity,
            pressure_potential_second_derivative: potential_second,
            pressure_second_derivative: pressure_second,
            fluid_force_second_derivative: force_second,
            displaced_flow_second_derivative: flow_second,
            relaxation_stress_second_derivative: second.stress,
            memory,
            mean_stored_energy_derivative,
            mean_stored_energy_second_derivative,
            power,
            relative_second_derivative_residual: native_second_residual
                .max(original_second_residual),
            original_relative_second_derivative_residual: original_second_residual,
        })
    }
}

fn observation(error: HarmonicResponseError) -> Error {
    match error {
        HarmonicResponseError::Cancelled => Error::Cancelled,
        _ => Error::Invalid("original physical frequency observation is unrepresentable"),
    }
}
fn quadratic(
    values: &[C64],
    first: &[C64],
    second: &[C64],
    matrix: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerCurvatureSample, Error> {
    harmonic_quadratic_power_curvature_with_cancel(values, first, second, matrix, cancelled)
        .map_err(observation)
}
fn work(
    value: f64,
    flow: [&[C64]; 3],
    effort: [&[C64]; 3],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicPowerCurvatureSample, Error> {
    let first = harmonic_mean_power_derivative_with_cancel(
        flow[0],
        flow[1],
        effort[0],
        effort[1],
        &mut *cancelled,
    )
    .map_err(observation)?;
    let second = harmonic_mean_power_curvature_with_cancel(
        flow[0], flow[1], flow[2], effort[0], effort[1], effort[2], cancelled,
    )
    .map_err(observation)?;
    Ok(HarmonicPowerCurvatureSample {
        sample: HarmonicPowerSample {
            average_power_w: value,
            angular_derivative: first,
        },
        angular_second_derivative: second,
    })
}
fn projection(weights: &[f64], values: [&[C64]; 3]) -> Result<[[C64; 1]; 3], Error> {
    let first =
        HarmonicTransferSample::project_real(weights, values[0], values[1]).map_err(observation)?;
    let second = HarmonicTransferSample::project_real(weights, values[2], values[2])
        .map_err(observation)?
        .amplitude;
    Ok([[first.amplitude], [first.angular_derivative], [second]])
}
fn zero_power() -> HarmonicPowerCurvatureSample {
    HarmonicPowerCurvatureSample {
        sample: HarmonicPowerSample {
            average_power_w: 0.0,
            angular_derivative: 0.0,
        },
        angular_second_derivative: 0.0,
    }
}
fn sum(
    values: &[HarmonicPowerCurvatureSample],
    weights: &[f64],
) -> Result<HarmonicPowerCurvatureSample, Error> {
    if values.is_empty() {
        return Ok(zero_power());
    }
    let value: Vec<_> = values
        .iter()
        .map(|p| C64::from_re(p.sample.average_power_w))
        .collect();
    let first: Vec<_> = values
        .iter()
        .map(|p| C64::from_re(p.sample.angular_derivative))
        .collect();
    let second: Vec<_> = values
        .iter()
        .map(|p| C64::from_re(p.angular_second_derivative))
        .collect();
    let zero = vec![C64::ZERO; values.len()];
    let result = HarmonicTransferParameterJet::project_real(
        weights,
        HarmonicTransferFieldParameterJet {
            values: &value,
            angular_derivative: &first,
            angular_second_derivative: &second,
            parameter_derivative: &zero,
            angular_parameter_derivative: &zero,
        },
    )
    .map_err(observation)?;
    Ok(HarmonicPowerCurvatureSample {
        sample: HarmonicPowerSample {
            average_power_w: result.sample.amplitude.re,
            angular_derivative: result.sample.angular_derivative.re,
        },
        angular_second_derivative: result.angular_second_derivative.re,
    })
}
