//! A material-derived curved body in the SAME reciprocal string mechanics.
//!
//! Shear-deformable shell geometry, elastic spectrum, full material loss,
//! physical bridge footprints and signed displaced volume share one basis.
//! Exterior observation and reciprocal cavity pressure use independently owned
//! surfaces; internal structural members need not border either acoustic fluid.
//! This is an unloaded linear shell discretization: no assembled violin anatomy,
//! prestressed equilibrium, full exterior radiation or measured timbre claim.

use super::*;
use super::plate::{SharedPlateBridge, material_body_preflight, finish_material_body_with_ports};
use super::super::plate::{MaterialBodyCavity, reduce_radiating_body, validate_body_medium};
use crate::shear_shell::{ResolvedShearShell, ShellProjectionBudget, ShellReductionOptions};
use crate::shear_shell::assembly::{MaterialShellAssembly, ShellAssemblyForcePort,
    ShellAssemblyReductionOptions};
use fs_plate::shell::mitc::{MitcBudget, NodalDirectors, ShellConstraint, ShellMass};

mod point;
mod relaxation;
mod substrate;
pub use point::{MaterialShellPoint, ShellPointGeometry};

/// Complete source-bound curved body and physical connection geometry.
#[derive(Clone, Debug)]
pub struct SharedMaterialShellBody {
    /// Actual material/section data on the original oriented triangular mesh.
    pub specimen: ResolvedShearShell,
    /// Explicit unit shell directors, admitted against the original facets.
    pub directors: NodalDirectors,
    /// Homogeneous supports, with work-conjugate physical/local axes.
    pub constraints: Vec<ShellConstraint>,
    /// Caller-selected consistent or lumped translation/director inertia.
    pub mass: ShellMass,
    /// Positive drilling stabilization; not physical material dissipation.
    pub drilling_alpha: f64,
    /// Angular-frequency-squared spectral window [rad²/s²].
    pub eigenvalue_window: (f64,f64),
    /// Retained spatial mode count; all requested modes must fit.
    pub n_modes: usize,
    /// Exactly one physical nodal force/moment row per material string.
    /// Each has 6*nodes GLOBAL entries. Moment/N entries have units of meters.
    /// The resulting scalar displacement is conjugate to string force.
    pub force_ports: Vec<Vec<f64>>,
    /// Independent bilateral connection law per string, in source order.
    pub bridges: Vec<SharedPlateBridge>,
    /// Cold FE assembly extent, without an inferred mesh quality claim.
    pub assembly_budget: MitcBudget,
    /// Cold full loss/force/volume projection work and channel scale.
    pub projection_budget: ShellProjectionBudget,
    /// Existing physical modal time-owner admission.
    pub modal_budget: ModalAcousticTimeBudget,
    /// Compact baffled monopole listener distance [m].
    pub listener_m: f64,
}

/// Separate material component surfaces joined before a common body eigensolve.
/// Every string, cavity and observer receives the SAME assembled modal basis.
#[derive(Clone, Debug)]
pub struct SharedMaterialShellAssemblyBody {
    /// Original material parts, global supports, directors, seams and fluid roles.
    pub specimen: MaterialShellAssembly,
    /// Whole-body spectral slice and cumulative geometry/projection limits.
    pub reduction: ShellAssemblyReductionOptions,
    /// Exactly one physical force/moment footprint per string; a footprint may
    /// act on several parts and is never replaced by a normalized fitted gain.
    pub force_ports: Vec<ShellAssemblyForcePort>,
    /// One physical string/body attachment law in original string order.
    pub bridges: Vec<SharedPlateBridge>,
    /// Original body time-domain discretization/admission.
    pub modal_budget: ModalAcousticTimeBudget,
    /// Physical compact exterior listener distance [m].
    pub listener_m: f64,
}
impl ReciprocalStringBank {
    /// Connect every material string to one complete, tied material shell body.
    /// Component stiffness/mass are assembled before any modal truncation, and
    /// all bulk/interface loss, pressure and force maps use that common basis.
    /// One cavity compression law sees the sum of its authored wall-volume rows.
    /// When the specimen selects coating relaxation, its returned image uses
    /// equilibrium film storage. Complete the cold bank with
    /// [`Self::with_material_shell_coating_relaxation`] before playback, so the
    /// film's projected histories join that same mechanical owner. The material
    /// strings frontend performs this binding after its gas-inertia selection.
    /// Explicit substrate relaxation is attached atomically by this constructor,
    /// before coatings, and follows the same source-preserving cold rebuilds.
    ///
    /// This retains the original reciprocal bow, sympathetic strings, torsion,
    /// gas loading and scheduled performance. It is an unloaded linear assembly
    /// with explicit rigid seams, a compact observer and optional uniform air;
    /// it does not claim prestressed anatomy or measured instrument accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_shell_assembly(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialShellAssemblyBody, cavity: Option<MaterialBodyCavity>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::from_material_shell_assembly_prepared(cx, sources, gas, body, cavity, None, law, options, gate)
    }

    /// `from_material_shell_assembly` with explicit ordinary, prepared or
    /// harmonic-only construction. Existing `Option<PreparedModalBudget>`
    /// callers retain their playback choice; harmonic-only banks retain all
    /// physical bindings but have no time propagator and refuse playback.
    /// An assembled body whose modes, cavity, string planes and gas memory exceed the ordinary dense block (42 modal
    /// coordinates) is refused by its own constructor, so the record has to be present at construction. The budget is
    /// kept in the bank's source for every later rebuild and must not be bound again.
    ///
    /// # Errors
    /// Those of `from_material_shell_assembly`, plus a prepared budget the complete network does not fit.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_shell_assembly_prepared(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialShellAssemblyBody, cavity: Option<MaterialBodyCavity>,
        prepared: impl Into<ModalPreparation>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let string_modes = material_body_preflight(&sources, &body.bridges,
            body.reduction.n_modes, cavity.as_ref(), &options, gate)?;
        if body.force_ports.len() != sources.len() {
            return Err(invalid("material shell assembly needs one original force footprint per string"));
        }
        let has_loss = body.specimen.has_viscous_glue() || body.specimen.components.iter().any(|part|
            !part.specimen.coatings().is_empty()
                || part.specimen.regions().iter().any(|r| !r.viscosity().is_zero()));
        let interfaces = body.specimen.components.iter().map(|part|
            part.specimen.interfaces().iter().filter(|i| i.law().damping_n_s_m != 0.0).count())
            .try_fold(0usize, usize::checked_add)
            .ok_or_else(|| invalid("shell assembly interface connection extent overflow"))?;
        let cavity_connections = cavity.as_ref().map_or(0, |c|
            1 + usize::from(c.cavity.neck_resistance_n_s_m() != 0.0));
        let connections = sources.len().checked_add(if has_loss { body.reduction.n_modes } else { 0 })
            .and_then(|n| n.checked_add(interfaces))
            .and_then(|n| n.checked_add(cavity_connections))
            .ok_or_else(|| invalid("shell assembly connection extent overflow"))?;
        if connections > options.network.max_connections
            || connections > crate::render::schedule::force::coupled::MAX_MODAL_CONNECTIONS {
            return Err(invalid("shell assembly loss, bridges and cavity exceed connection allowance"));
        }
        validate_body_medium(&gas, options.sample_rate_hz, body.listener_m,
            cavity.as_ref(), gate)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        let projected = body.specimen.reduce(&body.reduction, &body.force_ports, gate)
            .map_err(ReciprocalStringError::Material)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        let mut image = reduce_radiating_body(projected.radiators(), &gas,
            options.sample_rate_hz, body.modal_budget, body.listener_m, cavity.as_ref(),
            Some(projected.cavity_volume_shapes_m2_per_sqrt_kg()), gate)?;
        image.body_connections.extend(projected.damping().connections(0));
        image.body_connections.extend(projected.interface_damping_connections(0));
        let mut bank = finish_material_body_with_ports(cx, sources, &body.bridges,
            string_modes, gas, image, Some(projected.force_ports()), prepared, law, options, gate)?;
        bank.source.material_shell_assembly = Some(body);
        bank.source.material_shell_assembly_image = Some(projected);
        bank.source.material_cavity = cavity;
        bank.finish_substrate_relaxation(gate)
    }

    /// Compile true curved shear-deformable shell mechanics into the existing
    /// all-string/body/friction/cavity network. Material swaps replace actual
    /// element density, membrane/bending/shear stiffness and full viscous C.
    /// All bridge directions and body volume are projected on ONE mode basis.
    /// The existing pressure owner uses the shared gas and accepted acceleration.
    ///
    /// Body tangential translations and rotations participate in mechanics even
    /// when their compact monopole projection vanishes. Neither material loss
    /// nor pressure is assembled from an instrument name or an output spectrum.
    /// Material/geometry sources survive optional torsion/fluid/prepared, stops
    /// and multi-bow scheduling through the original bank's immutable source.
    ///
    /// Starts a cold unloaded realization, including explicitly selected
    /// laminate B/I1 coupling. Finite strain, prestress and exterior radiation
    /// loading require their own physical owners. A specimen that selects
    /// coating relaxation uses equilibrium film storage: call
    /// [`Self::with_material_shell_coating_relaxation`] on the returned cold
    /// bank before playback to attach the complete projected film histories.
    /// The material strings frontend performs this binding automatically.
    /// Explicit substrate memory is attached here automatically, initially
    /// relaxed, and participates in every subsequent force/pressure solve.
    /// # Errors
    /// Incomplete/detached ports, material/geometry/cavity mismatch, bandwidth,
    /// budget, spectral/section refusal or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_shear_shell(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialShellBody, cavity: Option<MaterialBodyCavity>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::from_material_shear_shell_prepared(cx, sources, gas, body, cavity, None, law, options, gate)
    }

    /// `from_material_shear_shell` with explicit ordinary, prepared or
    /// harmonic-only construction. Existing `Option<PreparedModalBudget>`
    /// callers retain their playback choice; harmonic-only banks retain all
    /// physical bindings but have no time propagator and refuse playback.
    /// A shell body whose retained modes, cavity, string planes and gas memory exceed the ordinary dense block (42
    /// modal coordinates) is refused by its own constructor, so the record has to be present at construction. The
    /// budget is kept in the bank's source for every later rebuild and must not be bound again.
    ///
    /// # Errors
    /// Those of `from_material_shear_shell`, plus a prepared budget the complete network does not fit.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_shear_shell_prepared(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialShellBody, cavity: Option<MaterialBodyCavity>,
        prepared: impl Into<ModalPreparation>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let string_modes=material_body_preflight(&sources,&body.bridges,body.n_modes,
            cavity.as_ref(),&options,gate)?;
        if body.force_ports.len()!=sources.len() {
            return Err(invalid("material shell requires one original physical bridge footprint per string"));
        }
        let has_loss=!body.specimen.coatings().is_empty()
            || body.specimen.regions().iter().any(|r|!r.viscosity().is_zero());
        let links=sources.len().checked_add(if has_loss {body.n_modes} else {0})
            .and_then(|n|n.checked_add(body.specimen.interfaces().iter()
                .filter(|i|i.law().damping_n_s_m!=0.0).count()))
            .and_then(|n|n.checked_add(cavity.as_ref().map_or(0,|c|
                1+usize::from(c.cavity.neck_resistance_n_s_m()!=0.0))))
            .ok_or_else(||invalid("shell material damping/bridge/cavity count overflow"))?;
        if links>options.network.max_connections
            || links>crate::render::schedule::force::coupled::MAX_MODAL_CONNECTIONS {
            return Err(invalid("complete shell damping, bridge and cavity channels exceed connection allowance"));
        }
        validate_body_medium(&gas,options.sample_rate_hz,body.listener_m,cavity.as_ref(),gate)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        let reduction=ShellReductionOptions {
            constraints:body.constraints.clone(), mass:body.mass, drilling_alpha:body.drilling_alpha,
            eigenvalue_window:body.eigenvalue_window, n_modes:body.n_modes,
            assembly_budget:body.assembly_budget, projection_budget:body.projection_budget,
        };
        let projected=body.specimen.reduce(&body.directors,&reduction,&body.force_ports,gate)
            .map_err(ReciprocalStringError::Material)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        let mut image=reduce_radiating_body(projected.radiators(),&gas,options.sample_rate_hz,
            body.modal_budget,body.listener_m,cavity.as_ref(),
            Some(projected.cavity_volume_shapes_m2_per_sqrt_kg()),gate)?;
        image.body_connections.extend(projected.damping().connections(0));
        image.body_connections.extend(projected.interface_damping_connections(0));
        let mut bank=finish_material_body_with_ports(cx,sources,&body.bridges,string_modes,
            gas,image,Some(projected.force_ports()),prepared,law,options,gate)?;
        bank.source.material_shell=Some(body);
        bank.source.material_shell_image=Some(projected);
        bank.source.material_cavity=cavity;
        bank.finish_substrate_relaxation(gate)
    }
}
