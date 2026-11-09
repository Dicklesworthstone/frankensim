//! Complete physical forced-response jets of the same spatial structure/fluid.
//!
//! Implicit derivatives reuse the original factor. The source and coefficient
//! directions are authored in the original fixed coordinate basis; geometry
//! derivatives must already include trace, normal and area changes.

use super::{
    PreparedStructuralInterior, StructuralInteriorError, StructuralInteriorResponse,
    StructuralInteriorSources, invalid, poll, response,
};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::HarmonicTransferFieldParameterJet;
use fs_modal::harmonic::{HarmonicFrequencyParameterControls, HarmonicResponseError};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

/// The same complete physical direction used by the original coupled pole
/// owner. It assembles both opposite reciprocal interface blocks exactly once.
pub use super::super::poles::parameter::StructuralInteriorPoleDirection as StructuralInteriorParameterDirection;

mod fields;
pub use fields::{
    StructuralInteriorDerivativeFields, StructuralInteriorPowerParameterJet,
    StructuralInteriorResponseField,
};

/// Complete local source law at one physical frequency and parameter.
/// Each Some contains both force and inward-flow arrays, including explicit
/// zeros for an unchanged source. None declares that entire derivative zero.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorSourceParameterJet<'s> {
    pub values: StructuralInteriorSources<'s>,
    pub angular_derivative: Option<StructuralInteriorSources<'s>>,
    pub angular_second_derivative: Option<StructuralInteriorSources<'s>>,
    pub parameter_derivative: Option<StructuralInteriorSources<'s>>,
    pub angular_parameter_derivative: Option<StructuralInteriorSources<'s>>,
}

impl<'s> StructuralInteriorSourceParameterJet<'s> {
    /// Hold both original complex source phasors fixed in frequency and theta.
    #[must_use]
    pub fn fixed(values: StructuralInteriorSources<'s>) -> Self {
        Self {
            values,
            angular_derivative: None,
            angular_second_derivative: None,
            parameter_derivative: None,
            angular_parameter_derivative: None,
        }
    }
}

/// Complete solved fields for actual driven peaks, cancellation notches and
/// physical power linewidths. All derivatives are partials at the original
/// frequency/parameter, not total derivatives along a moving extremum.
#[derive(Clone, Debug)]
pub struct StructuralInteriorFrequencyParameterResponse {
    pub response: StructuralInteriorResponse,
    /// Physical derivatives per rad/s.
    pub angular_derivative: StructuralInteriorDerivativeFields,
    /// Physical second derivatives per (rad/s)^2.
    pub angular_second_derivative: StructuralInteriorDerivativeFields,
    /// Physical derivatives per caller-declared unit of theta.
    pub parameter_derivative: StructuralInteriorDerivativeFields,
    /// Physical mixed derivatives per rad/s per theta unit.
    pub angular_parameter_derivative: StructuralInteriorDerivativeFields,
    pub power: StructuralInteriorPowerParameterJet,
    /// Original-coordinate equation backward errors, not forward accuracy
    /// bounds near a poorly conditioned resonance.
    pub relative_derivative_residual: f64,
    pub relative_second_derivative_residual: f64,
    pub relative_parameter_derivative_residual: f64,
    pub relative_mixed_derivative_residual: f64,
}

impl StructuralInteriorFrequencyParameterResponse {
    /// Borrow any complete physical field for the shared coherent receiver and
    /// resonance-parameter queries. Pressure uses p=-i*omega*psi, never psi.
    /// A moving/interferometric receiver supplies its own complete receiver
    /// jet to HarmonicTransferParameterJet::project_complex_with_derivatives.
    #[must_use]
    pub fn field_parameter_jet(
        &self,
        field: StructuralInteriorResponseField,
    ) -> HarmonicTransferFieldParameterJet<'_> {
        let values = match field {
            StructuralInteriorResponseField::StructuralDisplacement => {
                &self.response.structural_displacement
            }
            StructuralInteriorResponseField::StructuralVelocity => {
                &self.response.structural_velocity
            }
            StructuralInteriorResponseField::PressurePotential => {
                &self.response.pressure_potential_pa_s
            }
            StructuralInteriorResponseField::Pressure => &self.response.pressure_pa,
            StructuralInteriorResponseField::FluidForce => &self.response.fluid_force,
            StructuralInteriorResponseField::DisplacedFlow => &self.response.displaced_flow_m3_s,
        };
        HarmonicTransferFieldParameterJet {
            values,
            angular_derivative: self.angular_derivative.field(field),
            angular_second_derivative: self.angular_second_derivative.field(field),
            parameter_derivative: self.parameter_derivative.field(field),
            angular_parameter_derivative: self.angular_parameter_derivative.field(field),
        }
    }
}

impl PreparedStructuralInterior<'_, '_, '_, '_> {
    /// Differentiate the complete coupled forced response with this SAME factor.
    ///
    /// The original structural and acoustic mass/stiffness, structural viscous
    /// and gyroscopic entries, relative-wall damping and reciprocal interface
    /// all participate. The native implicit equations retain A_theta*q_omega,
    /// A_omega*q_theta and A_omega_theta*q; no perturbed solve, eigenmode
    /// derivative, diagonal damping model or line-shape fit occurs.
    ///
    /// M/C/K and their theta directions are frequency independent. Source
    /// amplitude and phase laws may depend on both variables through all five
    /// supplied original phasor arrays. Coordinates and topology remain fixed.
    /// Constitutive/geometry derivatives are not inferred by this adapter.
    ///
    /// Physical pressure and velocity retain every derivative of -i*omega.
    /// Interface output includes B_theta^T*p and B_theta*v; mixed output also
    /// includes B_theta^T*p_omega and B_theta*v_omega. Independent work jets
    /// include every source product rule and direct symmetric damping change.
    /// Gyroscopic/interface skew coefficients affect motion, never direct heat.
    ///
    /// This entry covers the existing memoryless whole structural/interior
    /// pencil. An additional exterior port load refuses in the native owner;
    /// its complete dispersive/material law cannot be silently held constant.
    ///
    /// Beyond native full-jet resources, controls reserve 24*n*n+1536*n+2048
    /// bytes and work score 2048*n*n+4096*n for physical direction assembly,
    /// source mapping, original field/output recovery and independent work.
    /// These exclude the immutable model/factor, borrowed caller arrays and
    /// allocator metadata. Refusal/cancellation publishes no partial response.
    pub fn solve_with_frequency_parameter_derivative(
        &self,
        sources: StructuralInteriorSourceParameterJet<'_>,
        direction: StructuralInteriorParameterDirection<'_>,
        controls: HarmonicFrequencyParameterControls,
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorFrequencyParameterResponse, StructuralInteriorError> {
        poll(cx)?;
        let ns = self.model.structural_dofs();
        let na = self.model.pressure_dofs();
        let n = ns
            .checked_add(na)
            .ok_or_else(|| invalid("physical parameter extent overflow"))?;
        let controls = native_controls(n, controls)?;
        let blocks = direction.assemble(ns, na, &mut || cx.checkpoint().is_err())?;
        let load = |source: StructuralInteriorSources<'_>| {
            self.load(source.structural_force, source.inward_flow_load)
        };
        let force = load(sources.values)?;
        let angular = sources.angular_derivative.map(load).transpose()?;
        let second = sources.angular_second_derivative.map(load).transpose()?;
        let parameter = sources.parameter_derivative.map(load).transpose()?;
        let mixed = sources.angular_parameter_derivative.map(load).transpose()?;
        let native = self.harmonic.solve_with_frequency_parameter_derivative(
            &force,
            QuadraticParameterDirection {
                mass: &blocks[0],
                damping: &blocks[1],
                stiffness: &blocks[2],
            },
            angular.as_deref(),
            second.as_deref(),
            parameter.as_deref(),
            mixed.as_deref(),
            controls,
            cx,
        )?;
        let zero = vec![C64::ZERO; n];
        let mut cancelled = || cx.checkpoint().is_err();
        let recovered = fields::recover(
            ns,
            &self.model.coupling,
            &self.model.damping,
            &blocks[1],
            direction.coupling,
            HarmonicTransferFieldParameterJet {
                values: &force,
                angular_derivative: angular.as_deref().unwrap_or(&zero),
                angular_second_derivative: second.as_deref().unwrap_or(&zero),
                parameter_derivative: parameter.as_deref().unwrap_or(&zero),
                angular_parameter_derivative: mixed.as_deref().unwrap_or(&zero),
            },
            native,
            &mut cancelled,
        )?;
        let mut baseline = recovered.response;
        let potential = baseline.displacement.split_off(ns);
        let pressure = baseline.velocity.split_off(ns);
        let energy = response::energy(
            &self.model.mass,
            &self.model.stiffness,
            &baseline.displacement,
            &baseline.velocity,
            &potential,
            &pressure,
            &mut cancelled,
        )?;
        let power = recovered.power;
        let response = StructuralInteriorResponse {
            structural_displacement: baseline.displacement,
            structural_velocity: baseline.velocity,
            pressure_potential_pa_s: potential,
            pressure_pa: pressure,
            fluid_force: recovered.fluid_force,
            displaced_flow_m3_s: recovered.displaced_flow,
            energy,
            average_structural_input_power_w: power.structural_input.sample.average_power_w,
            average_acoustic_input_power_w: power.acoustic_input.sample.average_power_w,
            average_fluid_to_structure_power_w: power.fluid_to_structure.sample.average_power_w,
            average_structural_dissipation_w: power.structural_dissipation.sample.average_power_w,
            average_wall_dissipation_w: power.wall_dissipation.sample.average_power_w,
            structural_port_load: baseline.port_load,
            power_balance_residual_w: baseline.power_balance_residual_w,
            relative_residual: baseline.relative_residual,
        };
        poll(cx)?;
        Ok(StructuralInteriorFrequencyParameterResponse {
            response,
            angular_derivative: recovered.angular_derivative,
            angular_second_derivative: recovered.angular_second_derivative,
            parameter_derivative: recovered.parameter_derivative,
            angular_parameter_derivative: recovered.angular_parameter_derivative,
            power,
            relative_derivative_residual: recovered.relative_derivative_residual,
            relative_second_derivative_residual: recovered.relative_second_derivative_residual,
            relative_parameter_derivative_residual: recovered
                .relative_parameter_derivative_residual,
            relative_mixed_derivative_residual: recovered.relative_mixed_derivative_residual,
        })
    }
}

fn native_controls(
    n: usize,
    controls: HarmonicFrequencyParameterControls,
) -> Result<HarmonicFrequencyParameterControls, StructuralInteriorError> {
    let n = n as u128;
    let bytes = 24 * n * n + 1536 * n + 2048;
    let work = 2048 * n * n + 4096 * n;
    let bytes =
        usize::try_from(bytes).map_err(|_| invalid("physical parameter payload overflow"))?;
    let work = u64::try_from(work).map_err(|_| invalid("physical parameter work overflow"))?;
    let maximum_additional_payload_bytes = controls
        .maximum_additional_payload_bytes
        .checked_sub(bytes)
        .ok_or(HarmonicResponseError::Budget {
            quantity: "structural-interior parameter recovery payload bytes",
            requested: bytes,
            cap: controls.maximum_additional_payload_bytes,
        })?;
    let maximum_work =
        controls
            .maximum_work
            .checked_sub(work)
            .ok_or(HarmonicResponseError::WorkBudget {
                requested: work,
                cap: controls.maximum_work,
            })?;
    Ok(HarmonicFrequencyParameterControls {
        maximum_additional_payload_bytes,
        maximum_work,
    })
}
