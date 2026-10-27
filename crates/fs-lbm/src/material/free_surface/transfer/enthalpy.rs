//! Preserve spatial enthalpy while tetrahedral liquid enters native VOF motion.
//!
//! The parent geometry owner integrates P1 specific enthalpy over the SAME
//! tet/cell intersections as momentum. Native thermal transport then carries
//! mass*h over its original mass routes and resolves local source viscosity.
//! Mechanical projection error stays in the parent's report and is never heat.
//!
//! Input is already liquid and stress-free on one constant-density caloric
//! chart. This does not release a stressed solid, erase plastic history,
//! convert phases, introduce pressure work or supply capillary surface energy.

use std::collections::BTreeSet;

use fs_exec::{CancelGate, Cx};
use fs_matdb::ClaimSelection;
use fs_material::state_point::MaterialStatePointError;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_qty::Temperature;

use super::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceScaling, TetLiquidTransferError,
    TetLiquidTransferLimits, TetLiquidTransferReport, TetrahedralLiquidState,
    admit_limits, checkpoint, invalid,
};
use super::super::enthalpy::{
    FreeSurfaceEnthalpyControls, FreeSurfaceEnthalpyError,
    MaterialFreeSurfaceEnthalpy3,
};

/// Motion and thermal storage on the same original P1 tetrahedral field.
#[derive(Clone, Copy, Debug)]
pub struct TetrahedralLiquidEnthalpyState<'a> {
    /// The actual geometry and nodal velocity; no independent thermal mesh.
    pub motion: TetrahedralLiquidState<'a>,
    /// Nodal specific enthalpy [J/kg], affine within each original tet.
    pub specific_enthalpy_j_kg: &'a [f64],
}

/// Independent geometric/mechanical and extensive-thermal acceptance bounds.
#[derive(Clone, Copy, Debug)]
pub struct TetEnthalpyTransferLimits {
    /// Existing exact-overlap, momentum and mechanical projection gates.
    pub mechanical: TetLiquidTransferLimits,
    /// Bound on source-node, caloric-knot and source-property-line samples.
    /// The total caloric knot count must also fit this bound before scanning.
    pub max_source_queries: usize,
    /// Absolute allowed target minus source extensive enthalpy [J].
    pub enthalpy_tolerance_j: f64,
}

/// Separate thermal conservation and mechanical representation discrepancies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TetEnthalpyTransferReport {
    /// Original geometric transfer and its signed kinetic/potential changes.
    pub mechanical: TetLiquidTransferReport,
    /// Exact P1 source integral sum(rho*V*sum(h_i)/4) [J].
    pub source_enthalpy_j: f64,
    /// Actual native VOF sum(m*h), including its carry [J].
    pub target_enthalpy_j: f64,
    /// Target minus source enthalpy [J], without mechanical-error conversion.
    pub enthalpy_defect_j: f64,
    /// Actual source-node, caloric-knot and property-line sample queries.
    pub source_queries: usize,
}

/// Original geometric/momentum or thermal/material refusal.
#[derive(Debug)]
pub enum TetEnthalpyTransferError {
    Transfer(TetLiquidTransferError),
    Enthalpy(FreeSurfaceEnthalpyError),
}
impl core::fmt::Display for TetEnthalpyTransferError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Transfer(error) => core::fmt::Display::fmt(error, f),
            Self::Enthalpy(error) => core::fmt::Display::fmt(error, f),
        }
    }
}
impl std::error::Error for TetEnthalpyTransferError {}
impl From<TetLiquidTransferError> for TetEnthalpyTransferError {
    fn from(error: TetLiquidTransferError) -> Self { Self::Transfer(error) }
}
impl From<FreeSurfaceEnthalpyError> for TetEnthalpyTransferError {
    fn from(error: FreeSurfaceEnthalpyError) -> Self { Self::Enthalpy(error) }
}

impl<'s> MaterialFreeSurfaceScaling<'s> {
    /// Transfer actual P1 mass, motion and enthalpy, then attach the original
    /// evolving thermal/VOF owner. A returned candidate can immediately advance
    /// through its existing step, heating, conduction and viscosity interfaces.
    ///
    /// All referenced nodes and every caloric knot between their extrema must
    /// be fully liquid at the original density. Because P1 barycentric weights
    /// are nonnegative and the caloric chart interpolates these knots, source
    /// interiors cannot hide an unsupported phase/density interval. Source
    /// transport properties are resolved at these temperatures; the original
    /// thermal owner also admits each resulting cell at its actual averaged h.
    /// Original pinned source-line queries cover the complete temperature span,
    /// including each property's own breakpoints and sampled-surface topology;
    /// a source hole or interior density change refuses before geometric work.
    ///
    /// The same immutable card, pinned claims and complete query conventions
    /// must own the surface, temperature law and caloric chart, as enforced by
    /// the original thermal attachment. No material name selects coefficients.
    /// Native ingress retains exactly the supplied positive-mass wet cells.
    /// A missing thermal lane refuses; no fallback temperature is assigned.
    ///
    /// The caller supplies the execution context and its cancellation gate.
    /// Both are checked before publication. No input field or accepted flow is
    /// mutated on any source, projection, thermal, budget or cancellation refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn free_surface3_from_tetrahedra_with_enthalpy<'m>(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        grid_origin_m: [f64; 3],
        source: TetrahedralLiquidEnthalpyState<'_>,
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        limits: TetEnthalpyTransferLimits,
        controls: FreeSurfaceEnthalpyControls,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<(MaterialFreeSurfaceEnthalpy3<'s, 'm>, TetEnthalpyTransferReport), TetEnthalpyTransferError> {
        checkpoint(cx)?;
        if gate.is_requested() {
            return Err(TetLiquidTransferError::Cancelled.into());
        }
        admit_limits(limits.mechanical)?;
        if limits.max_source_queries == 0
            || caloric.knots().len() > limits.max_source_queries
            || !limits.enthalpy_tolerance_j.is_finite()
            || limits.enthalpy_tolerance_j < 0.0
            || source.motion.tetrahedra.is_empty()
            || source.motion.tetrahedra.len() > limits.mechanical.max_tetrahedra
            || source.specific_enthalpy_j_kg.len() != source.motion.positions_m.len()
        {
            return Err(invalid("invalid bounded tetrahedral enthalpy source").into());
        }
        if caloric.material_card_identity() != self.material().material().card_identity()
            || law.card().content_hash() != caloric.material_card_identity()
        {
            return Err(invalid("tetrahedral thermal fields require one original material card").into());
        }
        let density = self.material().density_kg_m3();
        let mut queried = BTreeSet::new();
        let mut source_queries = 0usize;
        let mut minimum_h = f64::INFINITY;
        let mut maximum_h = f64::NEG_INFINITY;
        for tet in source.motion.tetrahedra {
            for &index in tet {
                checkpoint(cx)?;
                let h = *source.specific_enthalpy_j_kg.get(index as usize)
                    .ok_or_else(|| invalid("tetrahedron index lies outside the enthalpy field"))?;
                if queried.insert(index) {
                    admit_source_sample(h, density, law, caloric, limits, &mut source_queries)?;
                    minimum_h = minimum_h.min(h);
                    maximum_h = maximum_h.max(h);
                }
            }
        }
        for knot in caloric.knots() {
            checkpoint(cx)?;
            if knot.specific_enthalpy_j_kg >= minimum_h
                && knot.specific_enthalpy_j_kg <= maximum_h
            {
                admit_source_sample(
                    knot.specific_enthalpy_j_kg, density, law, caloric,
                    limits, &mut source_queries,
                )?;
            }
        }
        admit_source_lines(
            minimum_h, maximum_h, density, law, caloric, limits,
            &mut source_queries, cx,
        )?;
        let (flow, mechanical, transferred) = self.transfer_tetrahedra(
            domain, grid_origin_m, source.motion, Some(source.specific_enthalpy_j_kg),
            limits.mechanical, cx,
        )?;
        let roundoff = 64.0 * f64::EPSILON * minimum_h.abs().max(maximum_h.abs()).max(1.0);
        for &specific in transferred.by_cell.values() {
            checkpoint(cx)?;
            if specific < minimum_h - roundoff || specific > maximum_h + roundoff {
                return Err(invalid("overlap projection left the source enthalpy convex hull").into());
            }
        }
        let thermal = flow.with_enthalpy(
            law, caloric, controls,
            |coordinates| transferred.by_cell.get(&coordinates).copied().unwrap_or(f64::NAN), gate,
        )?;
        let target_enthalpy_j = thermal.total_enthalpy_j()?;
        let enthalpy_defect_j = target_enthalpy_j - transferred.source_j;
        if !enthalpy_defect_j.is_finite()
            || enthalpy_defect_j.abs() > limits.enthalpy_tolerance_j
        {
            return Err(invalid("geometric transfer exceeds its extensive enthalpy bound").into());
        }
        checkpoint(cx)?;
        if gate.is_requested() {
            return Err(TetLiquidTransferError::Cancelled.into());
        }
        Ok((thermal, TetEnthalpyTransferReport {
            mechanical,
            source_enthalpy_j: transferred.source_j,
            target_enthalpy_j,
            enthalpy_defect_j,
            source_queries,
        }))
    }
}

fn admit_source_sample(
    h: f64,
    density: f64,
    law: &LiquidTemperatureLaw<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    limits: TetEnthalpyTransferLimits,
    queries: &mut usize,
) -> Result<(), TetEnthalpyTransferError> {
    *queries = queries.checked_add(1)
        .filter(|count| *count <= limits.max_source_queries)
        .ok_or_else(|| invalid("tetrahedral thermal source-query budget exhausted"))?;
    let phase = caloric.state_at_specific_enthalpy(h)
        .map_err(FreeSurfaceEnthalpyError::from)?;
    if phase.phase() != SolidLiquidPhase::Liquid
        || phase.bulk_density_kg_m3().to_bits() != density.to_bits()
    {
        return Err(invalid("source thermal field must remain liquid at the original density").into());
    }
    admit_transport_temperature(phase.temperature_k(), density, law)
}

fn admit_transport_temperature(
    temperature_k: f64,
    density: f64,
    law: &LiquidTemperatureLaw<'_>,
) -> Result<(), TetEnthalpyTransferError> {
    let material = law.at_temperature(Temperature::new(temperature_k))
        .map_err(FreeSurfaceEnthalpyError::from)?;
    if material.density_kg_m3().to_bits() != density.to_bits() {
        return Err(invalid("source thermal field requires a density/volume-change owner").into());
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn admit_source_lines(
    minimum_h: f64,
    maximum_h: f64,
    density: f64,
    law: &LiquidTemperatureLaw<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    limits: TetEnthalpyTransferLimits,
    queries: &mut usize,
    cx: &Cx<'_>,
) -> Result<(), TetEnthalpyTransferError> {
    let point_at = |h| -> Result<fs_matdb::QueryPoint, TetEnthalpyTransferError> {
        let temperature = caloric.state_at_specific_enthalpy(h)
            .map_err(FreeSurfaceEnthalpyError::from)?.temperature_k();
        law.reference_point().clone().with(law.temperature_axis(), temperature)
            .map_err(|_| invalid("invalid retained liquid temperature coordinate").into())
    };
    let lower = point_at(minimum_h)?;
    let upper = point_at(maximum_h)?;
    let card = law.card();
    for property in law.reference_state().material().properties() {
        checkpoint(cx)?;
        let selected = property.answer().receipt.selected;
        let claim = card.claims().claim(selected)
            .ok_or_else(|| invalid("original liquid property claim is absent"))?;
        if claim.value.axes().any(|(axis, dims)|
            axis == law.temperature_axis() && dims != Temperature::DIMS)
        {
            return Err(invalid("source temperature axis has incompatible dimensions").into());
        }
        let line = card.claims().query_line_typed(
            &claim.key, law.temperature_axis(), &lower, &upper,
            ClaimSelection::Pinned(selected),
        ).map_err(|source| FreeSurfaceEnthalpyError::Material(
            MaterialStatePointError::Support {
                property: claim.key.name().to_owned(), source: Box::new(source),
            },
        ))?;
        for (temperature, _) in line.knots {
            checkpoint(cx)?;
            *queries = queries.checked_add(1)
                .filter(|count| *count <= limits.max_source_queries)
                .ok_or_else(|| invalid("tetrahedral thermal source-query budget exhausted"))?;
            admit_transport_temperature(temperature, density, law)?;
        }
    }
    Ok(())
}
