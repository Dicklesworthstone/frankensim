//! Original electrical interface laws in live material MHD.
//!
//! The shared interface owner retains ordered source responses and explicit
//! openings. The existing moving-current topology owner then supplies current,
//! contact heat and reciprocal Lorentz force to the unchanged MHD iteration.

use super::*;
use crate::conjugate::ElectricalContact;
use crate::conjugate::enthalpy::transport::interfaces::{
    ElectricalInterfaceState, SwitchingInterfaceEpoch, SwitchingInterfaceSample, admit_interfaces,
};
use std::borrow::Cow;
use std::sync::Arc;

/// Existing electrical/thermal/mechanical accounts with their original contacts.
#[derive(Debug)]
pub struct MhdInterfaceReport<J> {
    pub mhd: MhdEnthalpyReport,
    /// Contacts used by the last accepted interval, not rejected momentum trials.
    /// None for a zero-duration advance. Payloads do not need to implement Clone.
    pub last_applied: Option<Arc<SwitchingInterfaceEpoch<J>>>,
    /// Fresh endpoint support at the report's endpoint h/T and native velocity.
    /// This electrical evaluation deposits no additional heat or force impulse.
    pub endpoint: Arc<SwitchingInterfaceEpoch<J>>,
}

pub(super) struct ResolvedContacts<'a, E> {
    pub applied: Cow<'a, [ElectricalContact]>,
    pub open: Vec<[usize; 2]>,
    pub epoch: E,
}

// Share the complete MHD driver between unchanged authored contacts and original
// material interfaces. Only lowering the source response differs; native current,
// time integration, momentum iteration and all work/heat accounting stay shared.
pub(super) trait MhdInterfaces<F, T>: Clone {
    type Epoch: Clone;
    fn evaluate<'a>(
        &self,
        input: JouleEnthalpyInput<'a>,
        time_s: f64,
        samples: &[FlowMaterialSample<F, T>],
        gate: &CancelGate,
    ) -> Result<ResolvedContacts<'a, Self::Epoch>, JouleEnthalpyError>;
}

#[derive(Clone, Copy)]
pub(super) struct FixedInterfaces;
impl<F, T> MhdInterfaces<F, T> for FixedInterfaces {
    type Epoch = ();
    fn evaluate<'a>(
        &self,
        input: JouleEnthalpyInput<'a>,
        _: f64,
        _: &[FlowMaterialSample<F, T>],
        _: &CancelGate,
    ) -> Result<ResolvedContacts<'a, ()>, JouleEnthalpyError> {
        Ok(ResolvedContacts {
            applied: Cow::Borrowed(input.contacts),
            open: Vec::new(),
            epoch: (),
        })
    }
}

type InterfaceLaw<'a, F, T, J> = dyn Fn(
        usize,
        &ElectricalContact,
        [&FlowMaterialSample<F, T>; 2],
    ) -> Result<ElectricalInterfaceState<J>, JouleEnthalpyError>
    + 'a;

struct MaterialInterfaces<'a, F, T, J>(&'a InterfaceLaw<'a, F, T, J>);
impl<F, T, J> Clone for MaterialInterfaces<'_, F, T, J> {
    fn clone(&self) -> Self {
        Self(self.0)
    }
}
impl<F, T, J> MhdInterfaces<F, T> for MaterialInterfaces<'_, F, T, J> {
    type Epoch = Arc<SwitchingInterfaceEpoch<J>>;
    fn evaluate<'a>(
        &self,
        input: JouleEnthalpyInput<'a>,
        time_s: f64,
        samples: &[FlowMaterialSample<F, T>],
        gate: &CancelGate,
    ) -> Result<ResolvedContacts<'a, Self::Epoch>, JouleEnthalpyError> {
        // Use the original interface owner's bounded topology admission before
        // ordered indexing or calling any constitutive source.
        admit_interfaces(input, gate)?;
        let mut contacts = Vec::with_capacity(input.contacts.len());
        let mut applied = Vec::with_capacity(input.contacts.len());
        let mut open = Vec::new();
        for (index, &authored) in input.contacts.iter().enumerate() {
            checkpoint(gate)?;
            let sides = authored.cells.map(|c| &samples[c]);
            let response = (self.0)(index, &authored, sides)?;
            let sample = SwitchingInterfaceSample { authored, response };
            match sample.applied() {
                Some(contact) => applied.push(contact),
                None => open.push(authored.cells),
            }
            contacts.push(sample);
        }
        checkpoint(gate)?;
        Ok(ResolvedContacts {
            applied: Cow::Owned(applied),
            open,
            epoch: Arc::new(SwitchingInterfaceEpoch { time_s, contacts }),
        })
    }
}

/// Evolve original material contacts, conducting fluid, momentum and enthalpy.
///
/// `input.contacts` selects ordered face locations. The interface callback reads
/// the SAME retained adjacent material/phase samples as bulk current, viscosity,
/// conductivity and force coupling. It supplies a supported finite resistance
/// and explicit heat partition, or an explicit electrical opening. Its concrete
/// source payload survives in accepted and endpoint epochs without a Clone bound.
/// Use ElectricalInterfaceResponse::try_from_source to bind an original resolved
/// interface source to the actual ordered SurfaceSpecs; no temperature, pressure,
/// phase, surface condition, interface law or heat partition is guessed here.
///
/// Every momentum trial evaluates these contacts with its same material epoch
/// and velocity in ONE native moving-current solve. That solution supplies bulk
/// and contact heat once and its reciprocal Lorentz force. An opening removes
/// the face from current, motional voltage, contact heat and Lorentz transfer;
/// thermal conduction, material, mass and enthalpy remain owned independently.
/// A zero resistance is perfect contact, never an opening. Unsupported source
/// data must be an error. Fixed authored coefficients are reference data only
/// and are not added to the resolved resistance. Unlisted faces retain the
/// original perfect-contact semantics.
///
/// Callbacks must be pure with respect to physical history: rejected Picard
/// trials may repeat an evaluation. Source refusal, disconnected unreferenced
/// conductors, convergence, budgets or cancellation publish no candidate flow,
/// thermal or observer history. A zero-duration advance still checks endpoint
/// support, without depositing energy. Electrical openings do not split the
/// fluid domain or create a mechanical/thermal boundary.
///
/// All low-induction, fixed-field, creeping-flow and work gates of
/// advance_mhd_enthalpy_with_material_flow apply. Contact geometry is fixed and
/// carries no separate energy or momentum storage. Interfaces are sampled at
/// each actual thermal interval start; discontinuities require explicit advance
/// endpoints and timestep refinement. No event location or induction is added.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_mhd_enthalpy_with_material_interfaces<'t, 'f, F, R, O, T, J>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    config: StepConfig,
    input: JouleEnthalpyInput<'_>,
    motion: MhdEnthalpyInput<'_>,
    flow: &MaterialFlowFeedback<'f, F, R, O, T>,
    resistivity: impl Fn(&FlowMaterialSample<F, T>) -> Result<Option<f64>, ChtError>,
    interface: impl Fn(
        usize,
        &ElectricalContact,
        [&FlowMaterialSample<F, T>; 2],
    ) -> Result<ElectricalInterfaceState<J>, JouleEnthalpyError>,
) -> Result<
    (
        DepositionSimulation<'t>,
        MaterialFlowFeedback<'f, F, R, O, T>,
        MhdInterfaceReport<J>,
    ),
    MaterialFlowMhdError<O::Error>,
>
where
    F: Clone,
    T: FlowTransportSnapshot,
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
    O: DepositionObserver,
{
    let (next, accepted, mhd, last_applied, endpoint) = advance_mhd_enthalpy_impl(
        thermal,
        target_s,
        config,
        input,
        motion,
        flow,
        false,
        resistivity,
        MaterialInterfaces(&interface),
    )?;
    Ok((
        next,
        accepted,
        MhdInterfaceReport {
            mhd,
            last_applied,
            endpoint,
        },
    ))
}

/// The same source-resolved interfaces through externally anchored solid/mushy/
/// liquid phases. Frozen conductors retain their current and interface heat;
/// their Lorentz load belongs to the unchanged zero-work support account.
/// All anchored-phase admission and material-interface requirements still apply.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_mhd_enthalpy_with_anchored_interfaces<'t, 'f, F, R, O, J>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    config: StepConfig,
    input: JouleEnthalpyInput<'_>,
    motion: MhdEnthalpyInput<'_>,
    flow: &MaterialFlowFeedback<'f, F, R, O, PhaseFlowTransport>,
    resistivity: impl Fn(&FlowMaterialSample<F, PhaseFlowTransport>) -> Result<Option<f64>, ChtError>,
    interface: impl Fn(
        usize,
        &ElectricalContact,
        [&FlowMaterialSample<F, PhaseFlowTransport>; 2],
    ) -> Result<ElectricalInterfaceState<J>, JouleEnthalpyError>,
) -> Result<
    (
        DepositionSimulation<'t>,
        MaterialFlowFeedback<'f, F, R, O, PhaseFlowTransport>,
        MhdInterfaceReport<J>,
    ),
    MaterialFlowMhdError<O::Error>,
>
where
    F: Clone,
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
    O: DepositionObserver,
{
    let (next, accepted, mhd, last_applied, endpoint) = advance_mhd_enthalpy_impl(
        thermal,
        target_s,
        config,
        input,
        motion,
        flow,
        true,
        resistivity,
        MaterialInterfaces(&interface),
    )?;
    Ok((
        next,
        accepted,
        MhdInterfaceReport {
            mhd,
            last_applied,
            endpoint,
        },
    ))
}
