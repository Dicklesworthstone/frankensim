//! A common source-derived plate, not one independent sound filter per string.
//!
//! This compact image uses the uniform plate owner's existing force footprint
//! as a common bridge seat. Separate string spans can have different interior
//! bridge stations and connection laws. It is not a carved shell, distributed
//! bridge footprint, soundpost or broadband radiative back-loading model.
//! Optional body air reuses the existing lumped cavity/neck mechanics.

use super::*;
use crate::thin_plate::ResolvedPlateSpecimen;
use super::super::plate::{MaterialBodyCavity, MaterialBodyImage, reduce_material_body};

/// One string's bilateral support at the common plate force footprint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedPlateBridge {
    /// Interior station along this member's full fixed-end span.
    pub station_fraction: f64,
    /// Nonnegative coupling stiffness [N/m], independent of bulk material E.
    pub stiffness_n_m: f64,
    /// Nonnegative coupling resistance [N s/m], not inherited material damping.
    pub damping_n_s_m: f64,
    /// Signed stress-free string-minus-body displacement [m].
    pub rest_extension_m: f64,
}

/// One original body material/geometry with complete per-string bridge inputs.
#[derive(Clone, Debug)]
pub struct SharedMaterialPlateBody {
    /// Resolved plate material, thickness/mass constraint and support geometry.
    pub specimen: ResolvedPlateSpecimen,
    /// Exactly one connection per string, in original string declaration order.
    pub bridges: Vec<SharedPlateBridge>,
    /// Original modal time-domain admission for the body component.
    pub modal_budget: ModalAcousticTimeBudget,
    /// Compact baffled half-space observer distance [m]. No gain normalization.
    pub listener_m: f64,
}

impl ReciprocalStringBank {
    /// Compile all strings and ONE material plate in the same explicit gas.
    ///
    /// The original plate reducer determines frequencies, damping, masses,
    /// signed areas and mass-normalized bridge participation once. Every string
    /// connects to that same body, not to a cloned receiver. Pressure uses the
    /// body's actual endpoint acceleration, including all connection reactions.
    /// The same gas supplies every string's drag and this observer's density.
    /// All original string/layer/body material sources survive scheduling.
    ///
    /// This starts a new frozen-state realization at rest. It neither transfers
    /// thermal history nor morphs an already vibrating material into another.
    /// A Kelvin-Voigt string/plate model's validity is not full-instrument or
    /// measured timbre validation. The body image has no radiation back-loading.
    ///
    /// # Errors
    /// Incomplete bridge assignment, invalid supports/observer, unsupported body
    /// kinematics, material/loss/bandwidth refusal, setup limits or cancellation.
    pub fn from_material_plate(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateBody, law: FrictionLaw,
        options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::build_material_plate(cx, sources, gas, body, None, law, options, gate)
    }

    /// Add one shared body-air cavity and physical neck to ALL material strings.
    ///
    /// The original cavity law supplies compression, neck inertia and selected
    /// resistance. One compression spring loads the shared plate and neck;
    /// all string reactions see that same air load in the implicit bow solve.
    /// The compact pressure output includes actual body and outward neck
    /// accelerations. Cavity source/gas identity survives rendering/scheduling.
    /// This is a lumped, compact image, not a spatial f-hole or shell model.
    ///
    /// # Errors
    /// Stale cavity gas, invalid observer/geometry, budget and owner refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate_with_cavity(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateBody, cavity: MaterialBodyCavity,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::build_material_plate(cx, sources, gas, body, Some(cavity), law, options, gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn build_material_plate(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateBody, cavity: Option<MaterialBodyCavity>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        let string_modes = material_body_preflight(&sources, &body.bridges,
            body.specimen.plate().n_modes, cavity.as_ref(), &options, gate)?;
        let image = reduce_material_body(&body.specimen, &gas, options.sample_rate_hz,
            body.modal_budget, body.listener_m, cavity.as_ref(), gate)?;
        let mut bank = finish_material_body(cx, sources, &body.bridges,
            string_modes, gas, image, law, options, gate)?;
        bank.source.material_plate = Some(body);
        bank.source.material_cavity = cavity;
        poll(gate)?;
        Ok(bank)
    }
}

// Shared by uniform and regional geometry. Bound the complete retained network
// before either owner assembles and solves its finite-element eigenproblem.
pub(super) fn material_body_preflight(
    sources: &[MaterialStringSource], bridges: &[SharedPlateBridge], body_modes: usize,
    cavity: Option<&MaterialBodyCavity>, options: &StringBankOptions, gate: &CancelGate,
) -> Result<usize, ReciprocalStringError> {
    poll(gate)?;
    if sources.is_empty() || sources.len() != bridges.len() || body_modes == 0
        || sources.len() > options.network.max_connections || sources.len() > options.network.max_modes
        || options.bow.string >= sources.len() || !interior(options.bow.station_fraction) {
        return Err(invalid("material string bank needs a bounded, complete bridge assignment and an actual bow"));
    }
    let mut string_modes = 0usize;
    for (source, bridge) in sources.iter().zip(bridges) {
        poll(gate)?;
        if source.modes() == 0 || !interior(bridge.station_fraction)
            || !bridge.stiffness_n_m.is_finite() || bridge.stiffness_n_m < 0.0
            || !bridge.damping_n_s_m.is_finite() || bridge.damping_n_s_m < 0.0
            || (bridge.stiffness_n_m == 0.0 && bridge.damping_n_s_m == 0.0)
            || !bridge.rest_extension_m.is_finite() {
            return Err(invalid("every string needs modes and a finite nonzero interior bridge law"));
        }
        string_modes = string_modes.checked_add(source.modes())
            .ok_or_else(|| invalid("material string bank mode count overflow"))?;
    }
    let count = string_modes.checked_add(body_modes)
        .and_then(|n| n.checked_add(usize::from(cavity.is_some())))
        .ok_or_else(|| invalid("material string/body mode count overflow"))?;
    let links = sources.len().checked_add(cavity.map_or(0, |c|
        1 + usize::from(c.cavity.neck_resistance_n_s_m() != 0.0)))
        .ok_or_else(|| invalid("material string/body link count overflow"))?;
    if count > options.network.max_modes || links > options.network.max_connections {
        return Err(invalid("material string/body/cavity exceeds its complete-network setup budget"));
    }
    Ok(string_modes)
}

// One force/displacement binding and pressure-indexing owner for both shapes.
// The leading zero weights span ALL string modes, not the number of strings.
#[allow(clippy::too_many_arguments)]
pub(super) fn finish_material_body(
    cx: &Cx<'_>, sources: Vec<MaterialStringSource>, bridges: &[SharedPlateBridge],
    string_modes: usize, gas: GasState, image: MaterialBodyImage, law: FrictionLaw,
    options: StringBankOptions, gate: &CancelGate,
) -> Result<ReciprocalStringBank, ReciprocalStringError> {
    finish_material_body_with_ports(cx, sources, bridges, string_modes, gas, image,
        None, None, law, options, gate)
}

// Common and spatial bridge inputs differ only in the physical attachment map.
// Ports MUST belong to image's single eigenbasis; only the reducer calls here.
#[allow(clippy::too_many_arguments)]
pub(super) fn finish_material_body_with_ports(
    cx: &Cx<'_>, sources: Vec<MaterialStringSource>, bridges: &[SharedPlateBridge],
    string_modes: usize, gas: GasState, image: MaterialBodyImage,
    ports: Option<&[Vec<f64>]>, prepared: impl Into<ModalPreparation>, law: FrictionLaw,
    options: StringBankOptions, gate: &CancelGate,
) -> Result<ReciprocalStringBank, ReciprocalStringError> {
    poll(gate)?;
    if sources.len() != bridges.len() || ports.is_some_and(|p| p.len() != sources.len()) {
        return Err(invalid("every material string must retain exactly one body attachment"));
    }
    // All shape dimensions, nonzero support, finite values and constitutive
    // limits are checked by the original reciprocal network assembler below.
    let strings = sources.into_iter().zip(bridges).enumerate().map(|(i, (source, bridge))| MaterialStringMember {
        source, connections: vec![InteriorStringConnection {
            station_fraction: bridge.station_fraction,
            body: ModalAttachment { component: 0,
                shapes: ports.map_or_else(|| image.bridge.clone(), |p| p[i].clone()) },
            stiffness_n_m: bridge.stiffness_n_m, damping_n_s_m: bridge.damping_n_s_m,
            rest_extension_m: bridge.rest_extension_m,
        }],
    }).collect();
    let mut pressure = vec![0.0; string_modes];
    pressure.extend(image.pressure);
    let mut bank = ReciprocalStringBank::with_body_connections_prepared(cx, strings, gas, image.bodies,
        image.body_connections, prepared, law, options, gate)?;
    bank.system = bank.system.with_acceleration_pressure(pressure)
        .map_err(ReciprocalStringError::Friction)?;
    poll(gate)?;
    Ok(bank)
}
