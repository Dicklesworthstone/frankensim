//! Original quadrature history driven by a prescribed coherent strain history.

use super::{ViscoHarmonicMaterialBudget, ViscoHarmonicMaterialSpectrumPoint, cap};
use crate::linear3::visco::{
    ViscoElementError, ViscoElementState, ViscoElements, ViscoEnergy, ViscoHarmonicError,
    checkpoint,
};
use fs_exec::Cx;
use fs_fft::harmonic_window::phase_angular_interval_from_clock_dd;
use fs_material::visco::{
    MAX_HARMONIC_MATERIAL_LINES, TensorPronyHarmonicSpectrum, TensorPronyState,
};
use fs_math::{c64::C64, dd::Dd};

fn invalid(what: &'static str) -> ViscoElementError {
    ViscoElementError::InvalidHistory { what }
}
fn finite_clock(time: Dd) -> bool {
    time.hi.is_finite() && time.lo.is_finite()
}

/// An original accepted solid history with an explicit source-clock instant.
/// Binding assigns the instant of an existing history; it does not evolve,
/// relax, rephase or replace any of its original Maxwell internal variables.
#[derive(Clone)]
pub struct ViscoHarmonicMaterialHistory<'a> {
    history: ViscoElementState<'a>,
    time_s: Dd,
}
impl<'a> ViscoElementState<'a> {
    pub fn bind_harmonic_clock(
        self,
        time_s: Dd,
    ) -> Result<ViscoHarmonicMaterialHistory<'a>, ViscoElementError> {
        if !finite_clock(time_s) {
            return Err(invalid("finite original history clock required"));
        }
        let time_s = Dd::from_pair(time_s.hi, time_s.lo);
        if !finite_clock(time_s) {
            return Err(invalid("original history clock is unrepresentable"));
        }
        Ok(ViscoHarmonicMaterialHistory {
            history: self,
            time_s,
        })
    }
}
impl<'a> ViscoHarmonicMaterialHistory<'a> {
    pub const fn time_s_dd(&self) -> Dd {
        self.time_s
    }
    pub const fn history(&self) -> &ViscoElementState<'a> {
        &self.history
    }
    /// Continue with the original time-domain/dynamic owners without discarding
    /// the accepted material state or its exact represented clock.
    pub fn into_parts(self) -> (ViscoElementState<'a>, Dd) {
        (self.history, self.time_s)
    }
}

/// Independent physical interval account on every original material point.
/// The DD start and end retain short durations on a large source clock.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoHarmonicMaterialCausalField {
    pub start_time_s_dd: Dd,
    pub end_time_s_dd: Dd,
    pub duration_s: f64,
    pub angular_frequencies_rad_s: Vec<f64>,
    pub points: Vec<ViscoHarmonicMaterialSpectrumPoint>,
    pub start_stored_energy_j: f64,
    pub end_stored_energy_j: f64,
    pub stored_energy_change_j: f64,
    pub kelvin_heat_j: f64,
    pub maxwell_heat_j: f64,
    pub material_heat_j: f64,
    pub stress_work_j: f64,
    pub energy_balance_residual_j: f64,
}

/// Unpublished candidate. Its original per-point memory, endpoint stress,
/// full nodal displacement and cumulative ledgers can be accepted together.
pub struct ViscoHarmonicMaterialCausalStep<'a> {
    pub history: ViscoHarmonicMaterialHistory<'a>,
    pub field: ViscoHarmonicMaterialCausalField,
}

/// Frozen original geometry, laws and coherent frequency bank. Full node-major
/// input includes prescribed supports; no reduced coordinate or modal fit is
/// substituted for the actual material history.
pub struct ViscoHarmonicMaterialCausalObserver<'a> {
    source: &'a ViscoElements,
    frequencies: Vec<f64>,
    materials: Vec<TensorPronyHarmonicSpectrum<'a>>,
    points: usize,
}

impl ViscoElements {
    /// Admit causal observation and one complete candidate before evolution.
    /// Empty AC input admits free relaxation about the supplied static field.
    /// Payload includes retained spectra, the largest causal temporal cache,
    /// the returned field, endpoint Maxwell memory and full nodal displacement;
    /// caller-owned input histories, laws and phasors are excluded.
    /// Work includes each point's original material causal-window work bound,
    /// 4096*(F+1)*sum(branches+1) source preparation,
    /// 512*(F+1)*sum(points*local_dofs) original constitutive/compatible strain
    /// recovery, and F² + 256*(F+1)*full_dofs coherent nodal reconstruction,
    /// for authored F.
    pub fn prepare_harmonic_material_causal_observer(
        &self,
        frequencies: &[f64],
        budget: ViscoHarmonicMaterialBudget,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicMaterialCausalObserver<'_>, ViscoHarmonicError> {
        checkpoint(cx)?;
        if frequencies.len() > MAX_HARMONIC_MATERIAL_LINES
            || frequencies.iter().any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "causal material observation needs at most 64 positive finite frequencies",
            )
            .into());
        }
        let f = frequencies.len() as u128;
        let mut points = 0_u128;
        let mut bytes = core::mem::size_of::<ViscoHarmonicMaterialCausalObserver<'_>>() as u128
            + core::mem::size_of::<ViscoHarmonicMaterialCausalStep<'_>>() as u128
            + 128 * f
            + 8 * self.full_dofs as u128;
        let mut work = f * f + 256 * (f + 1) * self.full_dofs as u128;
        for element in &self.elements {
            checkpoint(cx)?;
            let b = element.law.branches().len() as u128;
            let p = element.points.len() as u128;
            points += p;
            bytes += core::mem::size_of::<TensorPronyHarmonicSpectrum<'_>>() as u128
                + 128 * f * (b + 1)
                + core::mem::size_of::<Vec<TensorPronyState>>() as u128;
            for point in &element.points {
                let fields = if point.compatible_gradients.is_some() {
                    2
                } else {
                    1
                };
                bytes += core::mem::size_of::<ViscoHarmonicMaterialSpectrumPoint>() as u128
                    + fields * f * core::mem::size_of::<[C64; 6]>() as u128
                    + core::mem::size_of::<TensorPronyState>() as u128
                    + 48 * b;
            }
            work += 4096 * (f + 1) * (b + 1) + 512 * (f + 1) * p * element.dofs.len() as u128;
        }
        cap(
            "causal harmonic material quadrature points",
            points,
            budget.maximum_points as u128,
        )?;
        cap(
            "causal harmonic material numerical payload bytes",
            bytes,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "causal harmonic material observation work",
            work,
            budget.maximum_work as u128,
        )?;
        let mut materials = Vec::with_capacity(self.elements.len());
        let mut temporal_bytes = 0;
        for element in &self.elements {
            checkpoint(cx)?;
            let material = element
                .law
                .prepare_harmonic_spectrum(frequencies)
                .map_err(ViscoElementError::Material)?;
            temporal_bytes = temporal_bytes.max(
                material
                    .causal_window_payload_bound_bytes()
                    .map_err(ViscoElementError::Material)?,
            );
            work += element.points.len() as u128
                * material
                    .causal_window_work_bound()
                    .map_err(ViscoElementError::Material)? as u128;
            materials.push(material);
        }
        cap(
            "causal harmonic material numerical payload bytes",
            bytes + temporal_bytes as u128,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "causal harmonic material observation work",
            work,
            budget.maximum_work as u128,
        )?;
        checkpoint(cx)?;
        Ok(ViscoHarmonicMaterialCausalObserver {
            source: self,
            frequencies: frequencies.to_vec(),
            materials,
            points: points as usize,
        })
    }
}

impl<'a> ViscoHarmonicMaterialCausalObserver<'a> {
    pub const fn source(&self) -> &'a ViscoElements {
        self.source
    }
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.frequencies
    }

    /// Advance the ORIGINAL accepted Maxwell memories through a prescribed
    /// coherent multitone interval. Every tone uses peak exp(-i omega t) on the
    /// bound source clock; `bias` and each tone contain ALL node-major xyz DOFs,
    /// including supports. The next interval begins at the returned DD clock.
    ///
    /// The initial displacement and original constitutive strain must match
    /// the waveform within scale-aware floating-point reconstruction error;
    /// this method never inserts an unaccounted strain jump. Full constitutive
    /// AND compatible B-bar waveforms pass original source admission, including
    /// DC coverage for nonzero strain bias. Initial memory need not be settled.
    /// Its natural decays are not a claim that startup is band limited.
    ///
    /// The source law and any accepted thermal eigenstrain stay frozen. The
    /// original thermal snapshot retains its source identities, raw samples
    /// and B-bar projection. Total and mechanical waveform envelopes must both
    /// satisfy their original material and thermal small-strain limits. A held
    /// free strain changes the mechanical bias but supplies no eigenstrain
    /// work during this interval. This is constitutive
    /// evolution, not a mechanical equilibrium, kinetic or boundary-work solve.
    /// Input history is borrowed; any refusal/cancellation publishes no partial
    /// next history, force, thermal deposit or cumulative-ledger change.
    pub fn advance_window(
        &self,
        accepted: &ViscoHarmonicMaterialHistory<'a>,
        displacements: &[&[C64]],
        bias: &[f64],
        duration_s: f64,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicMaterialCausalStep<'a>, ViscoHarmonicError> {
        checkpoint(cx)?;
        let old = &accepted.history;
        if !core::ptr::eq(self.source, old.system) {
            return Err(invalid(
                "causal observation requires the exact original solid history source",
            )
            .into());
        }
        let end_time = accepted.time_s + Dd::from_f64(duration_s);
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !finite_clock(end_time)
            || !accepted.time_s.lt(end_time)
            || displacements.len() != self.frequencies.len()
            || bias.len() != self.source.full_dofs
            || bias.iter().any(|x| !x.is_finite())
            || displacements.iter().any(|u| {
                u.len() != self.source.full_dofs
                    || u.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
        {
            return Err(invalid("causal material evolution needs finite full nodal fields and a positive advancing DD interval").into());
        }
        // A cancelled authored carrier must not enlarge the continuity
        // tolerance or erase the fixed bias before the cancellation resolves.
        // This also checks rigid/support motions, whose jumps cannot be found
        // by the constitutive strain checks at integration points.
        let mut coherent_frequencies = Vec::with_capacity(self.frequencies.len());
        let groups: Vec<_> = self
            .frequencies
            .iter()
            .map(|&omega| {
                coherent_frequencies
                    .iter()
                    .position(|&w| w == omega)
                    .unwrap_or_else(|| {
                        coherent_frequencies.push(omega);
                        coherent_frequencies.len() - 1
                    })
            })
            .collect();
        let phases = coherent_frequencies
            .iter()
            .map(|&omega| {
                phase_angular_interval_from_clock_dd(omega, accepted.time_s, duration_s)
                    .map_err(|_| invalid("causal nodal phase is unrepresentable"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut endpoint = Vec::with_capacity(self.source.full_dofs);
        for dof in 0..self.source.full_dofs {
            checkpoint(cx)?;
            let mut grouped = [[Dd::ZERO; 2]; MAX_HARMONIC_MATERIAL_LINES];
            for (tone, &group) in displacements.iter().zip(&groups) {
                grouped[group][0] = grouped[group][0] + Dd::from_f64(tone[dof].re);
                grouped[group][1] = grouped[group][1] + Dd::from_f64(tone[dof].im);
            }
            let mut start = Dd::ZERO;
            let mut end = Dd::ZERO;
            let mut scale = bias[dof].abs().max(old.displacement_m[dof].abs());
            for (z, phase) in grouped.iter().zip(&phases) {
                scale = scale
                    .max(z[0].hi.abs() + z[0].lo.abs())
                    .max(z[1].hi.abs() + z[1].lo.abs());
                start = start + z[0] * phase.start[0] - z[1] * phase.start[1];
                end = end + z[0] * phase.end[0] - z[1] * phase.end[1];
            }
            start = start + Dd::from_f64(bias[dof]);
            end = end + Dd::from_f64(bias[dof]);
            let difference = (start - Dd::from_f64(old.displacement_m[dof]))
                .to_f64()
                .abs();
            let tolerance = 128.0 * f64::EPSILON * (2 * phases.len() + 2) as f64;
            if !start.to_f64().is_finite()
                || !end.to_f64().is_finite()
                || !difference.is_finite()
                || !scale.is_finite()
                || (scale == 0.0 && difference != 0.0)
                || (scale > 0.0 && difference / scale > tolerance)
            {
                return Err(invalid(
                    "initial full nodal displacement is discontinuous with its harmonic history",
                )
                .into());
            }
            endpoint.push(end.to_f64());
        }
        let mut points = Vec::with_capacity(self.points);
        let mut histories = Vec::with_capacity(self.materials.len());
        let mut totals = [Dd::ZERO; 10];
        for (element, (source, material)) in
            self.source.elements.iter().zip(&self.materials).enumerate()
        {
            checkpoint(cx)?;
            let prepared = material
                .prepare_causal_window_from_clock_dd_with_checkpoint::<ViscoElementError>(
                    accepted.time_s,
                    duration_s,
                    || checkpoint(cx),
                )?;
            let mut states = Vec::with_capacity(source.points.len());
            for (quadrature_point, point) in source.points.iter().enumerate() {
                checkpoint(cx)?;
                let bias_at = |i: usize| bias[source.dofs[i]];
                let total_bias = point.strain_from_local(bias_at);
                let compatible_total_bias = point
                    .compatible_gradients
                    .as_ref()
                    .map(|_| point.compatible_strain_from_local(total_bias, bias_at));
                let mut strain = Vec::with_capacity(self.frequencies.len());
                let mut compatible = point
                    .compatible_gradients
                    .as_ref()
                    .map(|_| Vec::with_capacity(self.frequencies.len()));
                for tone in displacements {
                    let re = point.strain_from_local(|i| tone[source.dofs[i]].re);
                    let im = point.strain_from_local(|i| tone[source.dofs[i]].im);
                    strain.push(core::array::from_fn(|i| C64::new(re[i], im[i])));
                    if let Some(raw) = &mut compatible {
                        let re =
                            point.compatible_strain_from_local(re, |i| tone[source.dofs[i]].re);
                        let im =
                            point.compatible_strain_from_local(im, |i| tone[source.dofs[i]].im);
                        raw.push(core::array::from_fn(|i| C64::new(re[i], im[i])));
                    }
                }
                let mut bias_mandel = total_bias;
                let mut compatible_bias = compatible_total_bias;
                if let Some(thermal) = old.thermal_strains.as_deref() {
                    let sample = thermal.at(element, quadrature_point);
                    let ceiling = sample.maximum_absolute_mandel_strain();
                    let projected_free = thermal.free_strain(element, quadrature_point);
                    bias_mandel = core::array::from_fn(|i| total_bias[i] - projected_free[i]);
                    for offset in [&total_bias, &bias_mandel] {
                        material
                            .validate_kinematic_envelope_with_bias(&strain, offset, ceiling)
                            .map_err(ViscoElementError::Material)?;
                    }
                    if let (Some(raw), Some(total)) = (&compatible, compatible_total_bias) {
                        let raw_free = sample.free_strain_mandel();
                        let mechanical = core::array::from_fn(|i| total[i] - raw_free[i]);
                        for offset in [&total, &mechanical] {
                            material
                                .validate_kinematic_envelope_with_bias(raw, offset, ceiling)
                                .map_err(ViscoElementError::Material)?;
                        }
                        compatible_bias = Some(mechanical);
                    }
                }
                if let (Some(raw), Some(bias)) = (&compatible, &compatible_bias) {
                    material
                        .validate_strains_with_bias(raw, bias)
                        .map_err(ViscoElementError::Material)?;
                }
                let response = prepared.evaluate_with_bias_and_checkpoint::<ViscoElementError>(
                    &strain,
                    &bias_mandel,
                    &old.histories[element][quadrature_point],
                    || checkpoint(cx),
                )?;
                let energy = response.material;
                let volume = Dd::from_f64(point.volume_m3);
                for (sum, value) in totals.iter_mut().zip([
                    energy.start_stored_energy_j_m3,
                    energy.end_stored_energy_j_m3,
                    energy.stored_energy_change_j_m3,
                    energy.kelvin_heat_j_m3,
                    energy.maxwell_heat_j_m3,
                    energy.material_heat_j_m3,
                    energy.stress_work_j_m3,
                    response.state.stored,
                    response.state.work,
                    response.state.dissipated,
                ]) {
                    *sum = *sum + volume * Dd::from_f64(value);
                }
                states.push(response.state);
                points.push(ViscoHarmonicMaterialSpectrumPoint {
                    element,
                    quadrature_point,
                    volume_m3: point.volume_m3,
                    strain_mandel: strain,
                    compatible_strain_mandel: compatible,
                    bias_mandel,
                    compatible_bias_mandel: compatible_bias,
                    material: energy,
                });
            }
            histories.push(states);
        }
        let balance = (totals[6] - totals[2] - totals[5]).to_f64();
        let [
            start,
            end,
            change,
            kelvin,
            maxwell,
            heat,
            work,
            stored,
            cumulative_work,
            cumulative_heat,
        ] = totals.map(Dd::to_f64);
        if [
            start,
            end,
            change,
            kelvin,
            maxwell,
            heat,
            work,
            stored,
            cumulative_work,
            cumulative_heat,
            balance,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(ViscoElementError::NonFiniteResult.into());
        }
        checkpoint(cx)?;
        Ok(ViscoHarmonicMaterialCausalStep {
            history: ViscoHarmonicMaterialHistory {
                time_s: end_time,
                history: ViscoElementState {
                    system: self.source,
                    histories,
                    endpoint_stress_observed: true,
                    energy: ViscoEnergy {
                        stored_j: stored,
                        work_j: cumulative_work,
                        dissipated_j: cumulative_heat,
                    },
                    displacement_m: endpoint,
                    thermal_strains: old.thermal_strains.clone(),
                    eigenstrain_work_j: old.eigenstrain_work_j,
                    parameter_work_j: old.parameter_work_j,
                },
            },
            field: ViscoHarmonicMaterialCausalField {
                start_time_s_dd: accepted.time_s,
                end_time_s_dd: end_time,
                duration_s,
                angular_frequencies_rad_s: self.frequencies.clone(),
                points,
                start_stored_energy_j: start,
                end_stored_energy_j: end,
                stored_energy_change_j: change,
                kelvin_heat_j: kelvin,
                maxwell_heat_j: maxwell,
                material_heat_j: heat,
                stress_work_j: work,
                energy_balance_residual_j: balance,
            },
        })
    }
}
