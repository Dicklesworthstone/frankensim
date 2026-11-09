//! Explicit original equipment losses in the native regional heat inventory.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_qty::Pressure;
use fs_thermal::enthalpy::{CellId, DepositionSimulation};

use super::{ImplicitTransportError, TransientCarreauStep, finite, thermal_checkpoint};
use crate::graph::Sum;

/// Actual native cells receiving the original pump/drive routed fractions.
/// A positive authored fraction requires a recipient even at a zero-loss sample.
/// The fractions themselves remain owned by the original equipment source.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EquipmentHeatDestinations {
    pub pump: Option<CellId>,
    pub drive: Option<CellId>,
}

/// One branch's actual interval energies, all from its retained head evaluation.
/// The original power/routing value and source identity remain in the native
/// head report; these are integrated transfers, not another equipment model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquipmentEnergy {
    pub destinations: EquipmentHeatDestinations,
    pub hydraulic_work_j: f64,
    pub shaft_input_j: f64,
    pub electrical_input_j: f64,
    pub pump_loss_j: f64,
    pub drive_loss_j: f64,
    pub pump_heat_j: f64,
    pub drive_heat_j: f64,
    pub surroundings_heat_j: f64,
}

/// Electrical/source-boundary accounting for the same accepted caloric state.
/// Head-only sources stay explicit mechanical inputs; unknown input powers
/// remain unknown. Callback-authored heat is separately reported by the owner.
#[derive(Clone, Debug)]
pub struct EquipmentThermalReport {
    pub branches: Vec<Option<EquipmentEnergy>>,
    pub deposited_heat_j: f64,
    pub electrical_input_j: f64,
    pub surroundings_heat_j: f64,
    /// Held heads and original head sources lacking an input-power model.
    pub remaining_head_work_j: f64,
    /// Delta(U+K) - reservoir - injection - remaining head - electrical input
    /// + equipment surroundings - callback/native other heat - stream enthalpy.
    pub electrical_boundary_residual_j: f64,
    pub electrical_boundary_tolerance_j: f64,
}

pub(super) struct EquipmentDeposit {
    pub report: EquipmentThermalReport,
    pub sources_w: BTreeMap<CellId, f64>,
}

/// Read retained values only. No original source evaluation occurs here.
pub(super) fn prepare<E>(
    cx: &Cx<'_>,
    candidate: &DepositionSimulation<'_>,
    flow: &TransientCarreauStep,
    held_heads: Option<&[Pressure]>,
    destinations: Option<&[EquipmentHeatDestinations]>,
) -> Result<EquipmentDeposit, ImplicitTransportError<E>> {
    if flow.flow_heads.len() != flow.branches.len()
        || destinations.is_some_and(|d| d.len() != flow.branches.len())
        || held_heads.is_some_and(|h| h.len() != flow.branches.len())
    {
        return Err(ImplicitTransportError::Invalid(
            "one original equipment heat destination per branch required",
        ));
    }
    let mut branches = Vec::with_capacity(flow.branches.len());
    let mut sources = BTreeMap::<CellId, Sum>::new();
    let mut heat = Sum::default();
    let mut electrical = Sum::default();
    let mut surroundings = Sum::default();
    let mut remaining_head = Sum::default();
    for (branch, head) in flow.flow_heads.iter().enumerate() {
        thermal_checkpoint(cx)?;
        let dt = flow.branches[branch].duration_s;
        if let Some(heads) = held_heads {
            add(
                &mut remaining_head,
                finite(
                    (heads[branch].value() * flow.branches[branch].midpoint_volume_flow_m3_s) * dt,
                )?,
            )?;
        }
        let destination = destinations.map_or(EquipmentHeatDestinations::default(), |d| d[branch]);
        let Some((head, power)) = head
            .as_ref()
            .and_then(|head| head.response.power.map(|power| (head, power)))
        else {
            if destination.pump.is_some() || destination.drive.is_some() {
                return Err(ImplicitTransportError::Invalid(
                    "equipment heat requires original accepted input-power data",
                ));
            }
            if let Some(head) = head {
                add(&mut remaining_head, head.work_j)?;
            }
            branches.push(None);
            continue;
        };
        let r = power.routing;
        let state = power.state;
        let pump_w = finite(r.pump_to_reservoir * state.pump_loss_w)?;
        let drive_w = finite(r.drive_to_reservoir * state.drive_loss_w)?;
        for (fraction, cell, power_w) in [
            (r.pump_to_reservoir, destination.pump, pump_w),
            (r.drive_to_reservoir, destination.drive, drive_w),
        ] {
            if fraction > 0.0 && cell.is_none() {
                return Err(ImplicitTransportError::Invalid(
                    "positive original equipment heat routing requires an explicit recipient cell",
                ));
            }
            if let Some(cell) = cell {
                if !candidate.cells().contains_key(&cell) {
                    return Err(ImplicitTransportError::Invalid(
                        "equipment heat recipient is outside the native thermal inventory",
                    ));
                }
                add(sources.entry(cell).or_default(), power_w)?;
            }
        }
        let energy = EquipmentEnergy {
            destinations: destination,
            hydraulic_work_j: head.work_j,
            shaft_input_j: finite(state.shaft_input_w * dt)?,
            electrical_input_j: finite(state.electrical_input_w * dt)?,
            pump_loss_j: finite(state.pump_loss_w * dt)?,
            drive_loss_j: finite(state.drive_loss_w * dt)?,
            pump_heat_j: finite(pump_w * dt)?,
            drive_heat_j: finite(drive_w * dt)?,
            surroundings_heat_j: finite(state.surroundings_heat_w * dt)?,
        };
        add(&mut heat, energy.pump_heat_j)?;
        add(&mut heat, energy.drive_heat_j)?;
        add(&mut electrical, energy.electrical_input_j)?;
        add(&mut surroundings, energy.surroundings_heat_j)?;
        branches.push(Some(energy));
    }
    let mut sources_w = BTreeMap::new();
    for (cell, sum) in sources {
        thermal_checkpoint(cx)?;
        sources_w.insert(cell, value(sum)?);
    }
    Ok(EquipmentDeposit {
        report: EquipmentThermalReport {
            branches,
            deposited_heat_j: value(heat)?,
            electrical_input_j: value(electrical)?,
            surroundings_heat_j: value(surroundings)?,
            remaining_head_work_j: value(remaining_head)?,
            electrical_boundary_residual_j: 0.0,
            electrical_boundary_tolerance_j: 0.0,
        },
        sources_w,
    })
}

fn add<E>(sum: &mut Sum, value: f64) -> Result<(), ImplicitTransportError<E>> {
    sum.add(value).map_err(|_| {
        ImplicitTransportError::Invalid("unrepresentable equipment energy or deposited power")
    })
}
fn value<E>(sum: Sum) -> Result<f64, ImplicitTransportError<E>> {
    sum.value().map_err(|_| {
        ImplicitTransportError::Invalid("unrepresentable equipment energy or deposited power")
    })
}
