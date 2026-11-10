//! A closed, spatially uniform ideal-gas cavity on the composition contract.
//!
//! The original material parcel owns calorics, pressure, source admission and
//! the discrete Gibbs heat/work solve. This adapter exposes its work pressure
//! and entropy-conjugate temperature as ports; it supplies no second equation
//! of state, constant-capacity approximation or independent pressure quadrature.
//!
//! Pressure port `p` receives a compression flow [m3/s], positive when the
//! cavity shrinks, and returns absolute pressure [Pa]. Heat port `heat`
//! receives heat rate [W] into the gas and returns its thermodynamic mean
//! temperature [K]. Both are `FlowIn`. A piston-area transformer can join `p`
//! to an inertia's effort-input port without adding piston state to the gas.
//!
//! This frozen-composition component has no open mass/species streams,
//! reactions, nonideal equation of state or internal wave field. It does not
//! know a connected piston's mass or area and therefore cannot enforce the
//! native piston owner's coupled acoustic-frequency timestep bound. Callers
//! must resolve their assembled dynamics; a closed ledger is not an accuracy
//! certificate or a claim to replace the complete native piston solver.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::parcel::discrete_gradient::GasDiscreteGradientCandidate;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

pub use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
pub use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};

fn invalid(message: &'static str) -> GasEnergyError {
    GasEnergyError::Invalid(message)
}

fn poll(cx: &Cx<'_>) -> Result<(), GasEnergyError> {
    cx.checkpoint().map_err(|_| GasEnergyError::Cancelled)
}

fn pair(values: &[f64]) -> Result<[f64; 2], GasEnergyError> {
    match values {
        [a, b] if a.is_finite() && b.is_finite() => Ok([*a, *b]),
        _ => Err(invalid("gas cavity requires two finite scalar port values")),
    }
}

/// Source-backed closed gas with conservative pressure-work and heat ports.
///
/// Its [`Component::State`] is the original [`ClosedGasParcel`]. The immutable
/// reference fixes only an energy datum and the admitted inventory/source
/// identity. Measured energy is the original caloric increment from that
/// reference, so large formation-energy constants cannot erase small work.
#[derive(Clone, Debug)]
pub struct GasCavity {
    reference: ClosedGasParcel,
    controls: GasEnergySolve,
    ports: [PortSpec; 2],
}

impl GasCavity {
    /// Bind the existing parcel and its original caloric solve controls.
    /// Pass a clone of `reference` as the initial state when adding the component.
    ///
    /// # Errors
    /// Original source, temperature-branch, derivative or control refusal at
    /// the unchanged state. No material parameters or source ranges are inferred.
    pub fn new(reference: ClosedGasParcel, controls: GasEnergySolve) -> Result<Self, GasEnergyError> {
        reference.discrete_gradient_heat_volume(reference.volume_m3(), 0.0, controls, || false)?;
        Ok(Self {
            reference,
            controls,
            ports: [
                PortSpec::power(
                    StableId::new("p").expect("static port name"),
                    PortKind::FluidPressureFlux,
                    PortRole::FlowIn,
                    NonZeroUsize::MIN,
                )
                .expect("pressure-volume power port"),
                PortSpec::heat(
                    StableId::new("heat").expect("static port name"),
                    PortRole::FlowIn,
                )
                .expect("heat input port"),
            ],
        })
    }

    /// Original parcel defining the source/inventory identity and energy datum.
    #[must_use]
    pub const fn reference(&self) -> &ClosedGasParcel {
        &self.reference
    }

    /// Original source-temperature bracket, absolute energy budget and work cap.
    #[must_use]
    pub const fn controls(&self) -> GasEnergySolve {
        self.controls
    }

    fn trial(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(GasDiscreteGradientCandidate, f64), GasEnergyError> {
        poll(cx)?;
        let [compression, heat_rate] = pair(inputs)?;
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(invalid("gas cavity requires a positive finite timestep"));
        }
        // This also checks unchanged mass, composition, complete source model
        // and auxiliary query coordinates before the native proposal.
        let stored_old = self.reference.stored_energy_change_to(old, || cx.is_cancel_requested())?;
        let displaced = dt_s * compression;
        let heat = dt_s * heat_rate;
        let volume = old.volume_m3() - displaced;
        if !(volume.is_finite() && volume > 0.0)
            || !heat.is_finite()
            || (compression != 0.0 && displaced == 0.0)
            || (heat_rate != 0.0 && heat == 0.0)
        {
            return Err(invalid("gas cavity volume or heat increment is unrepresentable"));
        }
        let trial = old.discrete_gradient_heat_volume(volume, heat, self.controls, || {
            cx.is_cancel_requested()
        })?;
        Ok((trial, stored_old))
    }

    fn jacobian(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<[[f64; 2]; 2], GasEnergyError> {
        let (trial, _) = self.trial(old, dt_s, inputs, cx)?;
        let step = trial.step();
        let theta = trial.thermodynamic_mean().state().temperature;
        let pressure = step.work_pressure_pa;
        let pressure_q = step.pressure_heat_derivative_pa_j;
        // The native ideal-gas residual is e(T1)+R ln(V1/V0) theta(T1)-Q/m.
        // With p=m R a theta, a=ln(V1/V0)/(V1-V0), implicit differentiation
        // gives theta_Q=(theta/p)p_Q and theta_V=-(m R theta/V1)theta_Q.
        // The owner's p_Q already includes its original caloric mean derivative,
        // including the coincident-state limit and source-region quadrature.
        // These identities are specific to ClosedGasParcel's ideal-gas law.
        let theta_q = (theta / pressure) * pressure_q;
        let theta_v = -(old.mass_kg() * old.state().state().specific_gas_constant * theta
            / step.final_volume_m3)
            * theta_q;
        let matrix = [
            [-dt_s * step.pressure_final_volume_derivative_pa_m3, dt_s * pressure_q],
            [-dt_s * theta_v, dt_s * theta_q],
        ];
        if matrix.iter().flatten().any(|entry| !entry.is_finite()) {
            return Err(invalid("gas cavity port tangent is unrepresentable"));
        }
        poll(cx)?;
        Ok(matrix)
    }
}

impl Component for GasCavity {
    type State = ClosedGasParcel;
    type Error = GasEnergyError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        // There are two coupled nonlinear inputs; the scalar monotone
        // pressure solve does not classify the complete heat/volume map.
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &ClosedGasParcel) -> Result<f64, GasEnergyError> {
        self.reference.stored_energy_change_to(state, || false)
    }

    fn propose(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        new: &mut ClosedGasParcel,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, GasEnergyError> {
        if outputs.len() != 2 {
            return Err(invalid("gas cavity requires two scalar port outputs"));
        }
        let (trial, stored_old) = self.trial(old, dt_s, inputs, cx)?;
        let step = trial.step();
        let pressure = step.work_pressure_pa;
        let theta = trial.thermodynamic_mean().state().temperature;
        let stored_new = self.reference.stored_energy_change_to(trial.endpoint(), || {
            cx.is_cancel_requested()
        })?;
        let budget = old.mass_kg() * self.controls.energy_tolerance_j_kg;
        // Include p*V: forming V1=V0-dt*f rounds the represented displacement,
        // while the kernel independently books the prescribed dt*f on the wire.
        // That representation defect stays visible in the ledger, not in heat.
        let rounding = 64.0 * f64::EPSILON;
        let allowance = budget
            + rounding * stored_old.abs()
            + rounding * stored_new.abs()
            + rounding * step.stored_energy_change_j.abs()
            + rounding * step.supplied_heat_j.abs()
            + rounding * step.work_by_gas_j.abs()
            + rounding * (pressure * old.volume_m3()).abs()
            + rounding * (pressure * step.final_volume_m3).abs();
        let residual = step.energy_residual_j.abs()
            .max(step.gibbs_residual_j.abs())
            .max(theta * step.entropy_residual_j_k.abs());
        if !(allowance.is_finite() && residual.is_finite()) {
            return Err(invalid("gas cavity energy account is unrepresentable"));
        }
        let receipt = Receipt {
            stored_change_j: step.stored_energy_change_j,
            dissipated_j: 0.0,
            numerical_loss_j: -step.energy_residual_j,
            // The gas-side heat entropy is Q/theta. The native Gibbs owner
            // separately bounds Delta S-Q/theta; a numerical entropy residual
            // is not physical entropy generation.
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: allowance,
            residual,
            // Converting the owner's entropy bound budget/theta back to
            // joules can round upward once, although that native gate passed.
            tolerance: budget.next_up().min(f64::MAX),
        };
        poll(cx)?;
        let (endpoint, _) = trial.into_parts();
        *new = endpoint;
        outputs.copy_from_slice(&[pressure, theta]);
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for GasCavity {
    fn jvp(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), GasEnergyError> {
        let direction = pair(d_inputs)?;
        if d_outputs.len() != 2 {
            return Err(invalid("gas cavity requires two scalar tangent outputs"));
        }
        let j = self.jacobian(old, dt_s, inputs, cx)?;
        let result = [
            j[0][0] * direction[0] + j[0][1] * direction[1],
            j[1][0] * direction[0] + j[1][1] * direction[1],
        ];
        pair(&result)?;
        poll(cx)?;
        d_outputs.copy_from_slice(&result);
        Ok(())
    }

    fn vjp(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), GasEnergyError> {
        let weights = pair(bar_outputs)?;
        if bar_inputs.len() != 2 {
            return Err(invalid("gas cavity requires two scalar tangent inputs"));
        }
        let j = self.jacobian(old, dt_s, inputs, cx)?;
        let result = [
            j[0][0] * weights[0] + j[1][0] * weights[1],
            j[0][1] * weights[0] + j[1][1] * weights[1],
        ];
        pair(&result)?;
        poll(cx)?;
        bar_inputs.copy_from_slice(&result);
        Ok(())
    }
}
