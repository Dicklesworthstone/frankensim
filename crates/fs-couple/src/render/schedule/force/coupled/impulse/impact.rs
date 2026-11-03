//! One frictionless impact through the original reciprocal modal attachments.
//!
//! The native mechanical coordinates have identity mass. The FULL relative
//! row is assembled before its mobility b^T b is evaluated: overlapping ports
//! share modes and are not independent masses. Springs, finite dashpots and
//! proper material/acoustic memories carry no zero-duration impulse.
//!
//! This is the existing fs-contact single-contact Newton law, not a simultaneous
//! contact selection, event locator, friction law, or compliant-contact limit.
//! Spatial momentum is conserved only for symmetries represented in the basis;
//! a constrained/truncated component does not become a free rigid body.

use super::super::{
    CoupledModalSystem, ModalAttachment, ModalCouplingError, finite, invalid, poll,
};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use fs_contact::impact::{
    NormalImpactError, NormalImpactInput, NormalImpactResult, solve_normal_impact,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Simultaneous normal contacts on this same continuing modal state.
pub mod multiple;

/// Simultaneous normal/tangential Coulomb contacts on the original modal state.
pub mod frictional;

/// Simultaneous reciprocal impacts between native rigid bodies and modal points.
pub mod rigid;

/// Physical normal ports in one common scalar direction.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalImpactContact {
    /// Receives +J; positive motion here increases separation.
    pub first: ModalAttachment,
    /// Receives -J; positive motion here decreases separation.
    pub second: ModalAttachment,
    /// Separation at zero modal displacement [m]. Actual gap is this value
    /// plus (first - second) displacement, evaluated on the accepted state.
    pub reference_gap_m: f64,
    /// Nonnegative proximity activation [m]; penetration is never repaired.
    pub activation_tolerance_m: f64,
    /// Newton restitution in [0,1], explicitly supplied by the contact model.
    pub restitution: f64,
}

/// Cold preparation of one original relative row; event application allocates
/// nothing and never rebuilds or restarts the coupled mechanical/acoustic state.
pub struct PreparedModalImpact {
    column: Vec<f64>,
    layout: Vec<usize>,
    components: [usize; 2],
    reference_gap_m: f64,
    activation_tolerance_m: f64,
    restitution: f64,
    inverse_effective_mass: f64,
}

/// Actual accepted reset and the original contact law's independently solved
/// response. Dissipated energy is an unresolved physical loss account, not an
/// automatic thermal credit or numerical-error repair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalImpactReceipt {
    pub sample: u64,
    pub gap_m: f64,
    pub relative_normal_speed_before_m_s: f64,
    pub relative_normal_speed_after_m_s: f64,
    pub inverse_effective_mass_kg_inv: f64,
    pub contact: NormalImpactResult,
    pub kinetic_energy_change_j: f64,
    /// Actual kinetic change plus the contact law's nonnegative loss [J].
    pub energy_residual_j: f64,
    /// Actual endpoint relative speed minus the Newton target [m/s].
    pub restitution_residual_m_s: f64,
    pub maximum_relative_momentum_residual: f64,
}

/// Refusal never publishes a partial modal/cavity reset.
#[derive(Debug)]
pub enum ModalImpactError {
    Coupling(ModalCouplingError),
    Contact(NormalImpactError),
    /// The represented endpoint velocities fail the original restitution target.
    Restitution {
        residual_m_s: f64,
        tolerance_m_s: f64,
    },
}
impl From<ModalCouplingError> for ModalImpactError {
    fn from(value: ModalCouplingError) -> Self {
        Self::Coupling(value)
    }
}
impl From<NormalImpactError> for ModalImpactError {
    fn from(value: NormalImpactError) -> Self {
        Self::Contact(value)
    }
}
impl core::fmt::Display for ModalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Coupling(error) => core::fmt::Display::fmt(error, f),
            Self::Contact(error) => core::fmt::Display::fmt(error, f),
            Self::Restitution {
                residual_m_s,
                tolerance_m_s,
            } => write!(
                f,
                "modal impact restitution residual {residual_m_s:e} m/s exceeds {tolerance_m_s:e} m/s"
            ),
        }
    }
}
impl std::error::Error for ModalImpactError {}

impl PreparedModalImpact {
    pub fn new(
        system: &CoupledModalSystem,
        contact: ModalImpactContact,
        gate: &CancelGate,
    ) -> Result<Self, ModalImpactError> {
        poll(Some(gate))?;
        finite(contact.reference_gap_m)?;
        let column = relative_column(
            system,
            (contact.first.component, &contact.first.shapes),
            (contact.second.component, &contact.second.shapes),
            gate,
        )?;
        let mut mobility = Dd::ZERO;
        for (index, &shape) in column.iter().enumerate() {
            if index % 64 == 0 {
                poll(Some(gate))?;
            }
            mobility = mobility + Dd::from_f64(shape) * Dd::from_f64(shape);
        }
        let inverse_effective_mass = checked(mobility)?;
        // Delegate original restitution, activation and positive-mobility
        // admission to the law owner, including invalid zero relative rows.
        let _ = solve_normal_impact(NormalImpactInput {
            gap: 0.0,
            activation_tolerance: contact.activation_tolerance_m,
            relative_normal_speed: 0.0,
            inverse_effective_mass,
            restitution: contact.restitution,
        })?;
        poll(Some(gate))?;
        Ok(Self {
            column,
            layout: system
                .models
                .iter()
                .map(|model| model.modes().len())
                .collect(),
            components: [contact.first.component, contact.second.component],
            reference_gap_m: contact.reference_gap_m,
            activation_tolerance_m: contact.activation_tolerance_m,
            restitution: contact.restitution,
            inverse_effective_mass,
        })
    }

    /// Original component incidence, also used by cavity surface admission.
    pub fn components(&self) -> [usize; 2] {
        self.components
    }

    /// Complete relative mobility, including shared-mode cross terms [1/kg].
    pub fn inverse_effective_mass_kg_inv(&self) -> f64 {
        self.inverse_effective_mass
    }

    fn sample(
        &self,
        models: &[ModalAcousticTimeModel],
        gate: &CancelGate,
    ) -> Result<(f64, Dd), ModalImpactError> {
        if models.len() != self.layout.len()
            || models
                .iter()
                .zip(&self.layout)
                .any(|(model, &count)| model.modes().len() != count)
        {
            return Err(
                invalid("prepared impact no longer matches the original component layout").into(),
            );
        }
        let mut gap = Dd::from_f64(self.reference_gap_m);
        let mut speed = Dd::ZERO;
        for (index, (&shape, state)) in self
            .column
            .iter()
            .zip(models.iter().flat_map(|model| model.states()))
            .enumerate()
        {
            if index % 64 == 0 {
                poll(Some(gate))?;
            }
            let shape = Dd::from_f64(shape);
            gap = gap + shape * Dd::from_f64(state.displacement_m_sqrt_kg);
            speed = speed + shape * Dd::from_f64(state.velocity_m_sqrt_kg_per_s);
        }
        checked(speed)?;
        Ok((checked(gap)?, speed))
    }
}

impl CoupledModalSystem {
    /// Resolve one current frictionless impact and continue this same network.
    /// The full original relative mobility supplies the native Newton law;
    /// its impulse is applied by the transpose of that same complete row.
    /// Frozen positions, memory, pressure storage, prepared propagators, held
    /// forces, clock and preceding frame/port records remain untouched.
    ///
    /// Actual represented velocities must satisfy restitution and the original
    /// momentum/energy/state budgets before publication. No simultaneous-contact
    /// law is inferred by repeatedly calling this single-contact operation.
    pub fn apply_normal_impact(
        &mut self,
        impact: &PreparedModalImpact,
        gate: &CancelGate,
    ) -> Result<ModalImpactReceipt, ModalImpactError> {
        self.apply_normal_impact_with_admission(impact, gate, |_| Ok(()))
    }

    pub(crate) fn apply_normal_impact_with_admission(
        &mut self,
        impact: &PreparedModalImpact,
        gate: &CancelGate,
        admit: impl FnOnce(&[ModalAcousticTimeModel]) -> Result<(), ModalCouplingError>,
    ) -> Result<ModalImpactReceipt, ModalImpactError> {
        poll(Some(gate))?;
        self.require_time_domain()?;
        let (gap, before_sum) = impact.sample(&self.models, gate)?;
        let before = checked(before_sum)?;
        let contact = solve_normal_impact(NormalImpactInput {
            gap,
            activation_tolerance: impact.activation_tolerance_m,
            relative_normal_speed: before,
            inverse_effective_mass: impact.inverse_effective_mass,
            restitution: impact.restitution,
        })?;
        let mut receipt = ModalImpactReceipt {
            sample: self.samples_rendered(),
            gap_m: gap,
            relative_normal_speed_before_m_s: before,
            relative_normal_speed_after_m_s: before,
            inverse_effective_mass_kg_inv: impact.inverse_effective_mass,
            contact,
            kinetic_energy_change_j: 0.0,
            energy_residual_j: 0.0,
            restitution_residual_m_s: 0.0,
            maximum_relative_momentum_residual: 0.0,
        };
        if !contact.impacted {
            admit(&self.models)?;
            poll(Some(gate))?;
            return Ok(receipt);
        }
        let jump =
            self.stage_scaled_modal_impulse(&impact.column, contact.normal_impulse, Some(gate))?;
        let (_, after_sum) = impact.sample(&self.candidates, gate)?;
        let after = checked(after_sum)?;
        let residual = checked(after_sum - Dd::from_f64(contact.relative_normal_speed_after))?;
        let velocity_tolerance = finite(
            self.config.solve_relative_tolerance
                * before
                    .abs()
                    .max(after.abs())
                    .max(contact.relative_normal_speed_after.abs()),
        )?;
        if residual.abs() > velocity_tolerance {
            return Err(ModalImpactError::Restitution {
                residual_m_s: residual,
                tolerance_m_s: velocity_tolerance,
            });
        }
        let energy_residual = checked(
            Dd::from_f64(jump.kinetic_energy_change_j) + Dd::from_f64(contact.dissipated_energy),
        )?;
        let energy_tolerance = finite(
            self.config.energy_absolute_tolerance_j
                + self.config.energy_relative_tolerance
                    * jump
                        .kinetic_energy_change_j
                        .abs()
                        .max(contact.dissipated_energy),
        )?;
        if energy_residual.abs() > energy_tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: energy_residual,
                tolerance_j: energy_tolerance,
            }
            .into());
        }
        receipt.relative_normal_speed_after_m_s = after;
        receipt.kinetic_energy_change_j = jump.kinetic_energy_change_j;
        receipt.energy_residual_j = energy_residual;
        receipt.restitution_residual_m_s = residual;
        receipt.maximum_relative_momentum_residual = jump.maximum_relative_momentum_residual;
        admit(&self.candidates)?;
        poll(Some(gate))?;
        self.publish_modal_impulse();
        Ok(receipt)
    }
}

fn checked(value: Dd) -> Result<f64, ModalCouplingError> {
    finite(value.hi)?;
    finite(value.lo)?;
    finite(value.to_f64())
}

// One actual displacement/velocity trace and its reciprocal impulse column.
// Tangential rows may be zero; the normal-law owner separately admits inertia.
fn relative_column(
    system: &CoupledModalSystem,
    first: (usize, &[f64]),
    second: (usize, &[f64]),
    gate: &CancelGate,
) -> Result<Vec<f64>, ModalCouplingError> {
    let mut column = vec![0.0; system.mode_count()];
    for ((component, shapes), sign) in [(first, 1.0), (second, -1.0)] {
        let Some(model) = system.models.get(component) else {
            return Err(invalid(
                "impact attachment component is outside the original network",
            ));
        };
        if shapes.len() != model.modes().len() || shapes.iter().any(|value| !value.is_finite()) {
            return Err(invalid(
                "impact attachment must cover its complete original modal basis",
            ));
        }
        let offset = system.offsets[component];
        for (index, &shape) in shapes.iter().enumerate() {
            if index % 64 == 0 {
                poll(Some(gate))?;
            }
            column[offset + index] = finite(column[offset + index] + sign * shape)?;
        }
    }
    Ok(column)
}

// Audit the FULL compensated reciprocal projection, including its rounding
// tail, against represented endpoint increments. Contact cancellation cannot
// excuse a swallowed small net momentum change.
fn audit_modal_momentum(
    system: &CoupledModalSystem,
    after_models: &[ModalAcousticTimeModel],
    requested: &[Dd],
    gate: &CancelGate,
) -> Result<f64, ModalCouplingError> {
    let mut maximum = 0.0_f64;
    for (mode, (before, after)) in system
        .models
        .iter()
        .flat_map(|model| model.states())
        .zip(after_models.iter().flat_map(|model| model.states()))
        .enumerate()
    {
        if mode % 64 == 0 {
            poll(Some(gate))?;
        }
        let actual = Dd::from_f64(after.velocity_m_sqrt_kg_per_s)
            - Dd::from_f64(before.velocity_m_sqrt_kg_per_s);
        let scale = checked(actual)?.abs().max(checked(requested[mode])?.abs());
        let residual = checked(actual - requested[mode])?.abs();
        let relative = if scale == 0.0 { 0.0 } else { residual / scale };
        if relative > system.config.solve_relative_tolerance {
            return Err(ModalCouplingError::ImpulseResidual {
                mode,
                relative,
                tolerance: system.config.solve_relative_tolerance,
            });
        }
        maximum = maximum.max(relative);
    }
    Ok(maximum)
}
