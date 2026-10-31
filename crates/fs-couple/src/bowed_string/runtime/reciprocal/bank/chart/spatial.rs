//! Spatial string/body contacts, not independently synthesized voices.
//!
//! A nodal footprint w defines body displacement w^T u and applies force w F.
//! One shared plate reduction projects every footprint onto the SAME modes as
//! radiation, cavity volume and full tensor loss. Opposite forces at different
//! sites may therefore drive antisymmetric motion rather than cancel at one
//! fictitious common port. This does not add bridge mass, rotational inertia,
//! contact separation, carved geometry or an inferred lever-arm law.

use super::*;

/// Original physical unit-force footprints, in material-string order.
/// These are FE nodal weights, never precomputed or independently tuned modes.
#[derive(Clone, Debug, PartialEq)]
pub struct PlateBridgeFootprints {
    /// One complete nodal force distribution per string, in chart-node order.
    /// Each finite row must sum to one by the original plate projection rule.
    /// Signed entries are permitted; their physical justification belongs to
    /// the supplied attachment model. Repeated footprints remain distinct links.
    pub weights: Vec<Vec<f64>>,
    /// Original force-projection budget, separate from viscous projection work.
    /// Counts strings*(3*triangles + nodes*(1+requested_body_modes)).
    pub max_projection_terms: usize,
}

impl MaterialStringBankSource {
    /// Physical body footprints selected at construction. None retains the
    /// common reference footprint in material_chart(). Actual same-basis modal
    /// shapes are retained in strings()[i].connections[j].body.shapes.
    /// Sources survive both single- and multiple-bow scheduled performances.
    #[must_use]
    pub const fn material_chart_bridges(&self) -> Option<&PlateBridgeFootprints> {
        self.material_chart_bridges.as_ref()
    }
}

impl ReciprocalStringBank {
    /// Connect each string to its own physical body footprint without a second
    /// body solve. The common reference footprint still describes the reducer's
    /// reference drive only; it does not add another mechanical attachment.
    /// Existing isotropic source-loss selection and authored-loss behavior are
    /// unchanged. All original body, gas, support and modal limits still apply.
    ///
    /// A new frozen-state assembly starts at rest. Changing contact locations
    /// during playback requires a separate moving-contact work law; this method
    /// neither remaps a running state nor fabricates bridge inertia/rocking modes.
    /// The original compact pressure model and truncated flat DKT basis remain.
    ///
    /// # Errors
    /// Missing/extra footprints, invalid nodal weights, detached reduced ports,
    /// force-projection or original material/network budgets, and cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate_chart_with_force_ports(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, footprints: PlateBridgeFootprints,
        cavity: Option<MaterialBodyCavity>, law: FrictionLaw,
        options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::from_material_plate_chart_with_force_ports_prepared(cx, sources, gas,
            body, footprints, cavity, None, law, options, gate)
    }

    /// Spatial plate contacts with preparation selected before string history
    /// restoration. Each original force footprint shares the same full state.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plate_chart_with_force_ports_prepared(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, footprints: PlateBridgeFootprints,
        cavity: Option<MaterialBodyCavity>, prepared: impl Into<ModalPreparation>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::build_material_plate_chart(cx, sources, gas, body, cavity, prepared.into(),
            law, options, gate, None, Some(footprints))
    }

    /// The same spatial attachments WITH full, source-bound viscous cross terms.
    /// All force maps and the complete damping matrix share one original modal
    /// reduction. Damping and cavity reactions act back through each actual
    /// string footprint in the same simultaneous bow/bridge/cavity solve.
    /// Sources retain both nodal inputs and their resulting modal connections.
    ///
    /// # Errors
    /// All spatial-construction refusals plus stale/duplicate loss bindings,
    /// unsupported source bands and the original tensor/channel-work budgets.
    #[allow(clippy::too_many_arguments)]
    pub fn from_viscous_material_plate_chart_with_force_ports(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, footprints: PlateBridgeFootprints,
        viscosity: PlateViscosity, limits: PlateViscosityLimits,
        cavity: Option<MaterialBodyCavity>, law: FrictionLaw,
        options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::from_viscous_material_plate_chart_with_force_ports_prepared(cx, sources,
            gas, body, footprints, viscosity, limits, cavity, None, law, options, gate)
    }

    /// Spatial contacts and complete viscous cross terms with the propagator
    /// chosen before original material histories enlarge the connected state.
    #[allow(clippy::too_many_arguments)]
    pub fn from_viscous_material_plate_chart_with_force_ports_prepared(
        cx: &Cx<'_>, sources: Vec<MaterialStringSource>, gas: GasState,
        body: SharedMaterialPlateChartBody, footprints: PlateBridgeFootprints,
        viscosity: PlateViscosity, limits: PlateViscosityLimits,
        cavity: Option<MaterialBodyCavity>, prepared: impl Into<ModalPreparation>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::build_material_plate_chart(cx, sources, gas, body, cavity, prepared.into(),
            law, options, gate, Some((viscosity, limits)), Some(footprints))
    }
}

#[cfg(test)]
mod tests;
