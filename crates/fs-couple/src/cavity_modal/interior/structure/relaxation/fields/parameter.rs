//! Full original physical parameter response with every material history kept.

use super::*;
use crate::cavity_modal::interior::structure::resonance::{
    StructuralInteriorDerivativeFields, StructuralInteriorParameterDirection,
    StructuralInteriorResponseField, StructuralInteriorSourceParameterJet,
};
use fs_modal::harmonic::HarmonicResponseError;
use fs_modal::harmonic::resonance::{
    HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet, HarmonicTransferSample,
};
use fs_phs::mechanical::{
    MechanicalHarmonicDirection, MechanicalHarmonicDirectionControls,
    MechanicalRelaxationDirection, PreparedMechanicalHarmonicDirection,
};

mod observe;
mod residual;
pub use observe::{
    StructuralInteriorRelaxationMemoryParameterJet, StructuralInteriorRelaxationPowerParameterJet,
};

/// Authored local derivative of the original spatial and Maxwell laws.
/// Spatial K is relaxed. Branch order and coordinate topology remain fixed.
/// None declares zero branch directions; material/geometry derivatives are
/// never inferred from labels or a desired resonance shift.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationParameterDirection<'d> {
    pub spatial: StructuralInteriorParameterDirection<'d>,
    pub branches: Option<&'d [MechanicalRelaxationDirection<'d>]>,
}

struct BranchDirection {
    projection: Vec<f64>,
    stiffness: f64,
}

/// Cold mapping bound to this exact immutable original preparation.
/// Complete original direction blocks survive for independent physical audits.
pub(in super::super) struct PreparedDirection<'p, 's> {
    prepared: &'p Prepared<'s>,
    blocks: [Vec<f64>; 3],
    coupling: Vec<f64>,
    branches: Vec<BranchDirection>,
    native: PreparedMechanicalHarmonicDirection,
}

/// Full original memory response and its local physical parameter derivatives.
/// Angular derivatives are per rad/s; theta uses the caller's parameter unit.
/// These are partial derivatives at one frequency, not total derivatives along
/// a moving peak, finite material changes or transient parameter-change work.
#[derive(Clone, Debug)]
pub struct StructuralInteriorRelaxationFrequencyParameterResponse {
    pub frequency_response: StructuralInteriorRelaxationCurvatureResponse,
    pub parameter_derivative: StructuralInteriorDerivativeFields,
    pub angular_parameter_derivative: StructuralInteriorDerivativeFields,
    pub relaxation_stress_parameter_derivative: Vec<C64>,
    pub relaxation_stress_mixed_derivative: Vec<C64>,
    pub memory: Vec<StructuralInteriorRelaxationMemoryParameterJet>,
    pub power: StructuralInteriorRelaxationPowerParameterJet,
    pub mean_stored_energy_parameter_derivative: f64,
    pub mean_stored_energy_mixed_derivative: f64,
    /// Maximum original canonical and physical parameter-equation residual.
    pub relative_parameter_derivative_residual: f64,
    pub original_relative_parameter_derivative_residual: f64,
    /// Maximum original canonical and physical mixed-equation residual.
    pub relative_mixed_derivative_residual: f64,
    pub original_relative_mixed_derivative_residual: f64,
}

impl StructuralInteriorRelaxationFrequencyParameterResponse {
    /// Borrow the complete original field jet for coherent receivers and the
    /// existing peak/notch/linewidth parameter queries. Pressure remains the
    /// physical pressure, including every derivative of -i*omega*psi.
    pub fn field_parameter_jet(
        &self,
        which: StructuralInteriorResponseField,
    ) -> HarmonicTransferFieldParameterJet<'_> {
        let w = &self.frequency_response;
        let f = &w.frequency_response;
        let r = &f.response;
        let (values, first, second) = match which {
            StructuralInteriorResponseField::StructuralDisplacement => (
                &r.structural_displacement,
                &f.structural_displacement_derivative_per_rad_s,
                &w.structural_displacement_second_derivative,
            ),
            StructuralInteriorResponseField::StructuralVelocity => (
                &r.structural_velocity,
                &f.structural_velocity_derivative_per_rad_s,
                &w.structural_velocity_second_derivative,
            ),
            StructuralInteriorResponseField::PressurePotential => (
                &r.pressure_potential_pa_s,
                &f.pressure_potential_derivative_pa_s_per_rad_s,
                &w.pressure_potential_second_derivative,
            ),
            StructuralInteriorResponseField::Pressure => (
                &r.pressure_pa,
                &f.pressure_derivative_pa_per_rad_s,
                &w.pressure_second_derivative,
            ),
            StructuralInteriorResponseField::FluidForce => (
                &r.fluid_force,
                &f.fluid_force_derivative_per_rad_s,
                &w.fluid_force_second_derivative,
            ),
            StructuralInteriorResponseField::DisplacedFlow => (
                &r.displaced_flow_m3_s,
                &f.displaced_flow_derivative_m3_s_per_rad_s,
                &w.displaced_flow_second_derivative,
            ),
        };
        HarmonicTransferFieldParameterJet {
            values,
            angular_derivative: first,
            angular_second_derivative: second,
            parameter_derivative: field(&self.parameter_derivative, which),
            angular_parameter_derivative: field(&self.angular_parameter_derivative, which),
        }
    }

    /// One original material stress jet; a zero stress keeps ordinary complex
    /// derivatives, without dividing by an amplitude or inventing a phase.
    pub fn branch_stress_parameter_jet(
        &self,
        branch: usize,
    ) -> Result<HarmonicTransferParameterJet, Error> {
        let w = &self.frequency_response;
        let f = &w.frequency_response;
        let value = f
            .response
            .memory
            .get(branch)
            .ok_or(Error::Invalid(
                "original memory branch index is out of range",
            ))?
            .stress;
        let get = |v: &[C64]| {
            v.get(branch).copied().ok_or(Error::Invalid(
                "complete original branch derivative arrays required",
            ))
        };
        Ok(HarmonicTransferParameterJet {
            sample: HarmonicTransferSample {
                amplitude: value,
                angular_derivative: get(&f.relaxation_stress_derivative_per_rad_s)?,
            },
            angular_second_derivative: get(&w.relaxation_stress_second_derivative)?,
            parameter_derivative: get(&self.relaxation_stress_parameter_derivative)?,
            angular_parameter_derivative: get(&self.relaxation_stress_mixed_derivative)?,
        })
    }
}

fn field(
    fields: &StructuralInteriorDerivativeFields,
    which: StructuralInteriorResponseField,
) -> &[C64] {
    match which {
        StructuralInteriorResponseField::StructuralDisplacement => &fields.structural_displacement,
        StructuralInteriorResponseField::StructuralVelocity => &fields.structural_velocity,
        StructuralInteriorResponseField::PressurePotential => &fields.pressure_potential_pa_s,
        StructuralInteriorResponseField::Pressure => &fields.pressure_pa,
        StructuralInteriorResponseField::FluidForce => &fields.fluid_force,
        StructuralInteriorResponseField::DisplacedFlow => &fields.displaced_flow_m3_s,
    }
}

impl<'s> Prepared<'s> {
    pub(in super::super) fn prepare_parameter_direction<'p>(
        &'p self,
        direction: StructuralInteriorRelaxationParameterDirection<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedDirection<'p, 's>, Error> {
        poll(cancelled)?;
        admit_parameter(self)?;
        let n = self.source.coordinates;
        let blocks = direction.spatial.assemble(
            self.source.structural,
            n - self.source.structural,
            cancelled,
        )?;
        // Reuse the original admitted mechanical inverse, storage and arm map.
        // No model rebuild, factorization or frequency-dependent approximation.
        let native = self.mechanical.prepare_harmonic_parameter_direction(
            MechanicalHarmonicDirection {
                mass: Some(&blocks[0]),
                damping: Some(&blocks[1]),
                relaxed_stiffness: Some(&blocks[2]),
                force_map: None,
                branches: direction.branches,
            },
            MechanicalHarmonicDirectionControls {
                max_state_dimension: self.limits.numerical.max_dimension,
                max_dense_bytes: usize::try_from(self.limits.maximum_adapter_bytes)
                    .unwrap_or(usize::MAX),
                max_multiply_terms: u64::try_from(self.limits.maximum_adapter_work)
                    .unwrap_or(u64::MAX),
            },
            &mut *cancelled,
        )?;
        let mut branches = Vec::with_capacity(self.branches.len());
        for index in 0..self.branches.len() {
            poll(cancelled)?;
            let d = direction.branches.map(|v| v[index]).unwrap_or_default();
            branches.push(BranchDirection {
                projection: d
                    .projection
                    .map(<[f64]>::to_vec)
                    .unwrap_or_else(|| vec![0.0; n]),
                stiffness: d.stiffness,
            });
        }
        poll(cancelled)?;
        Ok(PreparedDirection {
            prepared: self,
            blocks,
            coupling: direction.spatial.coupling.to_vec(),
            branches,
            native,
        })
    }
}

impl Factor<'_, '_> {
    pub(in super::super) fn solve_with_frequency_parameter_derivative(
        &self,
        source: StructuralInteriorSourceParameterJet<'_>,
        direction: &PreparedDirection<'_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StructuralInteriorRelaxationFrequencyParameterResponse, Error> {
        poll(cancelled)?;
        let p = self.prepared;
        if !core::ptr::eq(p, direction.prepared) {
            return Err(Error::Invalid(
                "parameter direction belongs to another original preparation",
            ));
        }
        admit_parameter(p)?;
        let load =
            |s: crate::cavity_modal::interior::structure::resonance::StructuralInteriorSources<
                '_,
            >| p.load(s.structural_force, s.inward_flow_load);
        let input = load(source.values)?;
        let first = source.angular_derivative.map(load).transpose()?;
        let second = source.angular_second_derivative.map(load).transpose()?;
        let parameter = source.parameter_derivative.map(load).transpose()?;
        let mixed = source.angular_parameter_derivative.map(load).transpose()?;
        let native = self.harmonic.response_with_frequency_parameter_derivative(
            &input,
            direction.native.as_harmonic_direction(),
            first.as_deref(),
            second.as_deref(),
            parameter.as_deref(),
            mixed.as_deref(),
            &mut *cancelled,
        )?;
        let f = &native.frequency_response.frequency_response;
        let states = [
            p.state(
                &f.response.state,
                &f.response.effort,
                &f.response.output,
                cancelled,
            )?,
            p.state(
                &f.state_derivative,
                &f.effort_derivative,
                &f.output_derivative,
                cancelled,
            )?,
            p.state(
                &native.frequency_response.state_second_derivative,
                &native.frequency_response.effort_second_derivative,
                &native.frequency_response.output_second_derivative,
                cancelled,
            )?,
            parameter_state(
                p,
                direction,
                &native.parameter_response.state_derivative,
                &native.parameter_response.effort_derivative,
                &native.parameter_response.output_derivative,
                &f.response.effort,
                cancelled,
            )?,
            parameter_state(
                p,
                direction,
                &native.state_mixed_derivative,
                &native.effort_mixed_derivative,
                &native.output_mixed_derivative,
                &f.effort_derivative,
                cancelled,
            )?,
        ];
        let zero = vec![C64::ZERO; p.source.coordinates];
        let sources = HarmonicTransferFieldParameterJet {
            values: &input,
            angular_derivative: first.as_deref().unwrap_or(&zero),
            angular_second_derivative: second.as_deref().unwrap_or(&zero),
            parameter_derivative: parameter.as_deref().unwrap_or(&zero),
            angular_parameter_derivative: mixed.as_deref().unwrap_or(&zero),
        };
        let physical_parameter = residual::check(
            p,
            direction,
            self.omega,
            &states[3],
            &states[0],
            sources.parameter_derivative,
            None,
            cancelled,
        )?;
        let physical_mixed = residual::check(
            p,
            direction,
            self.omega,
            &states[4],
            &states[1],
            sources.angular_parameter_derivative,
            Some((&states[3], &states[0])),
            cancelled,
        )?;
        let observation = observe::observe(
            p,
            direction,
            &native,
            [&states[0], &states[1], &states[2], &states[3], &states[4]],
            sources,
            cancelled,
        )?;
        let mean_stored_energy_parameter_derivative =
            native.parameter_response.mean_stored_energy_derivative;
        let mean_stored_energy_mixed_derivative = native.mean_stored_energy_mixed_derivative;
        let relative_parameter_derivative_residual = native
            .parameter_response
            .relative_derivative_residual
            .max(physical_parameter);
        let relative_mixed_derivative_residual = native
            .relative_mixed_derivative_residual
            .max(physical_mixed);
        let [_, _, _, parameter, mixed] = states;
        let (parameter_derivative, relaxation_stress_parameter_derivative) = split(
            p.source.structural,
            parameter,
            observation.parameter_force,
            observation.parameter_flow,
        );
        let (angular_parameter_derivative, relaxation_stress_mixed_derivative) = split(
            p.source.structural,
            mixed,
            observation.mixed_force,
            observation.mixed_flow,
        );
        // Physical curvature recovery shares the already solved native image.
        let frequency_response = self.recover_curvature(
            native.frequency_response,
            &input,
            first.as_deref(),
            second.as_deref(),
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(StructuralInteriorRelaxationFrequencyParameterResponse {
            frequency_response,
            parameter_derivative,
            angular_parameter_derivative,
            relaxation_stress_parameter_derivative,
            relaxation_stress_mixed_derivative,
            memory: observation.memory,
            power: observation.power,
            mean_stored_energy_parameter_derivative,
            mean_stored_energy_mixed_derivative,
            relative_parameter_derivative_residual,
            original_relative_parameter_derivative_residual: physical_parameter,
            relative_mixed_derivative_residual,
            original_relative_mixed_derivative_residual: physical_mixed,
        })
    }
}

fn parameter_state(
    p: &Prepared<'_>,
    direction: &PreparedDirection<'_, '_>,
    state: &[C64],
    effort: &[C64],
    output: &[C64],
    base_effort: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<State, Error> {
    let n = p.source.coordinates;
    let mut stress = Vec::with_capacity(p.branches.len());
    for (j, (&root, &root_derivative)) in p
        .roots
        .iter()
        .zip(direction.native.branch_root_derivatives())
        .enumerate()
    {
        poll(cancelled)?;
        // sigma=-sqrt(k)*e_z. Both terms precede physical range recovery;
        // the native canonical effort derivative alone is not stress_theta.
        let value = HarmonicTransferSample::project_real(
            &[-root, -root_derivative],
            &[effort[2 * n + j], base_effort[2 * n + j]],
            &[C64::ZERO; 2],
        )
        .map_err(observation)?
        .amplitude;
        stress.push(value);
    }
    Ok(State {
        q: state[..n].to_vec(),
        velocity: output.to_vec(),
        stress,
    })
}
fn split(
    ns: usize,
    mut state: State,
    force: Vec<C64>,
    flow: Vec<C64>,
) -> (StructuralInteriorDerivativeFields, Vec<C64>) {
    let potential = state.q.split_off(ns);
    let pressure = state.velocity.split_off(ns);
    (
        StructuralInteriorDerivativeFields {
            structural_displacement: state.q,
            structural_velocity: state.velocity,
            pressure_potential_pa_s: potential,
            pressure_pa: pressure,
            fluid_force: force,
            displaced_flow_m3_s: flow,
        },
        state.stress,
    )
}
fn admit_parameter(p: &Prepared<'_>) -> Result<(), Error> {
    admit(p.source.coordinates, p.branches.len(), p.limits)?;
    let (n, b) = (p.source.coordinates as u128, p.branches.len() as u128);
    let d = 2 * n + b;
    for (quantity, requested, maximum) in [
        (
            "parameter adapter bytes",
            8 * (64 * d * d + 40 * n * n + 16 * b * n + 160 * d + 160 * n + 64 * b),
            p.limits.maximum_adapter_bytes,
        ),
        (
            "parameter adapter work",
            2048 * d * d * d + 512 * n * n + 1024 * b * n,
            p.limits.maximum_adapter_work,
        ),
    ] {
        if requested > maximum {
            return Err(Error::Limit {
                quantity,
                requested,
                maximum,
            });
        }
    }
    Ok(())
}
fn observation(error: HarmonicResponseError) -> Error {
    match error {
        HarmonicResponseError::Cancelled => Error::Cancelled,
        _ => Error::Invalid("original physical parameter observation is unrepresentable"),
    }
}
