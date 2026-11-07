//! Shared original-source cell Coulomb closure; no normal or wear integration.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_math::det;
use fs_tribo::{
    DryTangentialState, FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStep,
    TangentialStepInput,
};

use crate::interface_binding::BoundNormalInterfaceState;
use crate::interface_binding::thermal_friction::{SourceThermalCoulomb, SourceThermalCoulombError};

use super::{
    PlasticCellSliding, PlasticSlidingInput, SourcePlasticSlidingError, SourcePlasticSlidingWear,
    add_compensated, checkpoint,
};

pub(super) fn admit_query(
    normal: &BoundNormalInterfaceState,
    friction: &SourceThermalCoulomb,
    query: &QueryPoint,
) -> Result<(), SourcePlasticSlidingError> {
    if friction.interface() != normal.elastic().interface() || !friction.query_axes().is_empty() {
        return Err(SourcePlasticSlidingError::Invalid(
            "sliding requires the original ordered interface and only held temperature/local pressure",
        ));
    }
    let fixed = normal.interface_state().query();
    let pressure = friction.axes().nominal_pressure.as_deref();
    if !query
        .axes()
        .iter()
        .filter(|(name, _)| Some(name.as_str()) != pressure)
        .eq(fixed.axes().iter())
        || !query
            .axis_quantities()
            .iter()
            .filter(|(name, _)| Some(name.as_str()) != pressure)
            .eq(fixed.axis_quantities().iter())
    {
        return Err(SourcePlasticSlidingError::Invalid(
            "sliding must retain every original fixed query coordinate and quantity",
        ));
    }
    Ok(())
}

pub(super) fn query_at_temperature(
    friction: &SourceThermalCoulomb,
    temperature_k: f64,
) -> Result<QueryPoint, SourceThermalCoulombError> {
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(SourceThermalCoulombError::Invalid(
            "sliding requires a positive absolute interface temperature",
        ));
    }
    Ok(friction
        .initial_source()
        .query()
        .with(&friction.axes().temperature, temperature_k)?)
}

pub(super) fn coefficients_at_pressure(
    friction: &SourceThermalCoulomb,
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<(ResolvedInterfaceStatePoint, FrictionLaw), SourceThermalCoulombError> {
    let mut query = query_at_temperature(friction, temperature_k)?;
    if let Some(axis) = &friction.axes().nominal_pressure {
        query = query.with(axis, pressure_pa)?;
    }
    friction.coefficients_at_query(&query)
}

pub(super) struct CellSlidingAdvance {
    pub tangential: Vec<DryTangentialState>,
    pub cells: Vec<PlasticCellSliding>,
    pub friction_sources: Vec<ResolvedInterfaceStatePoint>,
    pub unsigned_slip_m: Vec<f64>,
    pub frictional_work_j: f64,
}

/// One original native solve per cell, with pressure-query sharing and no wear.
#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    friction: &SourceThermalCoulomb,
    initial_law: &FrictionLaw,
    temperature_k: f64,
    before: &[DryTangentialState],
    pressure: &[f64],
    area: f64,
    motion: &[PlasticSlidingInput],
    duration_s: f64,
    controls: TangentialSolveConfig,
    cx: &Cx<'_>,
) -> Result<CellSlidingAdvance, SourcePlasticSlidingError> {
    checkpoint(cx)?;
    if motion.len() != pressure.len() || before.len() != pressure.len() {
        return Err(SourcePlasticSlidingError::Invalid(
            "sliding requires one local motion for every original grid cell",
        ));
    }
    let mut sources = Vec::new();
    let mut laws = Vec::new();
    let mut sample_indices = BTreeMap::new();
    let mut tangential = Vec::with_capacity(pressure.len());
    let mut cells = Vec::with_capacity(pressure.len());
    let mut slip = Vec::with_capacity(pressure.len());
    let mut work = 0.0;
    let mut compensation = 0.0;
    for (index, (&p, &input)) in pressure.iter().zip(motion).enumerate() {
        checkpoint(cx)?;
        let error = |error| SourcePlasticSlidingError::Cell { index, error };
        let normal_n = SourcePlasticSlidingWear::cell_normal_force(p, area).map_err(error)?;
        // Detached cells need no source claim at an unsupported zero pressure.
        let source_index = if p == 0.0 {
            None
        } else {
            let key = if friction.axes().nominal_pressure.is_some() {
                p.to_bits()
            } else {
                0
            };
            let slot = if let Some(&slot) = sample_indices.get(&key) {
                slot
            } else {
                let (source, law) =
                    coefficients_at_pressure(friction, temperature_k, p).map_err(error)?;
                let slot = sources.len();
                sources.push(source);
                laws.push(law);
                sample_indices.insert(key, slot);
                slot
            };
            Some(slot)
        };
        let law = source_index.map_or(initial_law, |slot| &laws[slot]);
        let native = law
            .advance(
                &before[index],
                &TangentialStepInput {
                    normal_n,
                    slip_increment_m: input.released_increment_m,
                    dt_s: duration_s,
                    contact_temperature_k: Some(temperature_k),
                    nominal_pressure_pa: friction.axes().nominal_pressure.as_ref().map(|_| p),
                    compliance_m_per_n: input.compliance_m_per_n,
                    solve: controls,
                },
            )
            .map_err(|e| error(SourceThermalCoulombError::Tribo(e)))?;
        let distance = det::hypot(native.slip_increment_m[0], native.slip_increment_m[1]);
        if !distance.is_finite() {
            return Err(SourcePlasticSlidingError::Invalid(
                "finite actual local sliding distance",
            ));
        }
        slip.push(distance);
        add_compensated(&mut work, &mut compensation, native.dissipated_j);
        tangential.push(native.state);
        cells.push(PlasticCellSliding {
            input,
            pressure_pa: p,
            source_index,
            response: TangentialStep {
                state: (),
                traction_n: native.traction_n,
                slip_increment_m: native.slip_increment_m,
                regime: native.regime,
                stored_change_j: native.stored_change_j,
                dissipated_j: native.dissipated_j,
                closure_residual_j: native.closure_residual_j,
                compatibility_residual_m: native.compatibility_residual_m,
                iterations: native.iterations,
            },
        });
    }
    if !work.is_finite() || work < 0.0 {
        return Err(SourcePlasticSlidingError::Invalid(
            "finite summed native friction work",
        ));
    }
    checkpoint(cx)?;
    Ok(CellSlidingAdvance {
        tangential,
        cells,
        friction_sources: sources,
        unsigned_slip_m: slip,
        frictional_work_j: work,
    })
}
