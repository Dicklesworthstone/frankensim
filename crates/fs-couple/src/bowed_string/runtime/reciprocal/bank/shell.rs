//! A material-derived curved body in the SAME reciprocal string mechanics.
//!
//! Shear-deformable shell geometry, elastic spectrum, full material loss,
//! physical bridge footprints and signed displaced volume share one basis.
//! The compact baffled pressure/cavity owner observes accepted body acceleration.
//! This is an unloaded linear shell discretization: no assembled violin anatomy,
//! prestressed equilibrium, full exterior radiation or measured timbre claim.

use super::*;
use super::plate::{SharedPlateBridge, material_body_preflight, finish_material_body_with_ports};
use super::super::plate::{MaterialBodyCavity, reduce_radiating_body, validate_body_medium};
use crate::shear_shell::{ResolvedShearShell, ShellProjectionBudget, ShellReductionOptions};
use fs_plate::shell::mitc::{MitcBudget, NodalDirectors, ShellConstraint, ShellMass};

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
impl ReciprocalStringBank {
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
    /// Starts a cold unloaded realization. Finite strain, prestress, laminate
    /// B coupling, joints between separate parts and shell radiation loading
    /// require their own physical owners; they are not inferred by this adapter.
    /// # Errors
    /// Incomplete/detached ports, material/geometry/cavity mismatch, bandwidth,
    /// budget, spectral/section refusal or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_shear_shell(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialShellBody, cavity: Option<MaterialBodyCavity>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let string_modes=material_body_preflight(&sources,&body.bridges,body.n_modes,
            cavity.as_ref(),&options,gate)?;
        if body.force_ports.len()!=sources.len() {
            return Err(invalid("material shell requires one original physical bridge footprint per string"));
        }
        let has_loss=body.specimen.regions().iter().any(|r|!r.viscosity().is_zero());
        let links=sources.len().checked_add(if has_loss {body.n_modes} else {0})
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
            body.modal_budget,body.listener_m,cavity.as_ref(),gate)?;
        image.body_connections.extend(projected.damping().connections(0));
        let mut bank=finish_material_body_with_ports(cx,sources,&body.bridges,string_modes,
            gas,image,Some(projected.force_ports()),law,options,gate)?;
        bank.source.material_shell=Some(body);
        bank.source.material_shell_image=Some(projected);
        bank.source.material_cavity=cavity;
        poll(gate)?;
        Ok(bank)
    }
}
