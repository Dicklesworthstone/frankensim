//! Conservative isothermal diffusion of insoluble surfactant between material faces.
//!
//! Each authored undirected edge supplies a nonnegative conductance `g` [m²/s]
//! for `dN_i/dt = g * (N_j/A_j - N_i/A_i)`. Areas, conductances and one common
//! constitutive law are held during this call. An exact two-face exponential
//! propagator preserves the area-weighted concentration mean and is a convex
//! averaging map. Thus every exchange preserves positivity and the concentration
//! hull and decreases the SAME convex ideal/Langmuir free energy. This would
//! not hold for arbitrary different temperatures or adsorption affinities.
//!
//! Palindromic half-step edge sweeps give second-order temporal splitting at
//! fixed geometry: Strang (1968), <https://doi.org/10.1137/0705041>. Positivity
//! is unconditional; accuracy still requires timestep refinement. The caller
//! owns the spatial discretization: no two-point consistency on arbitrary
//! nonorthogonal triangles is inferred from an edge list. Disconnected graph
//! components cannot exchange material. No adsorption, bulk exchange, heat,
//! remeshing or evolving-area interpolation occurs here.
//!
//! Work is O(faces + substeps * edges), with one candidate amount vector and
//! precomputed edge propagators. No allocation or exponential occurs inside a
//! sweep. Canonical caller edge order and deterministic elementary math fix the
//! arithmetic. Cancellation or any admission/balance failure publishes nothing.

use super::super::Sum;
use super::{SurfactantError, SurfactantFace, SurfactantLaw};
use fs_math::det;

/// One symmetric diffusion connection, supplied once in a fixed canonical order.
/// Repeated pairs represent additive parallel conductances, not duplicate removal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceDiffusionEdge {
    pub faces: [usize; 2],
    /// Surface diffusivity times the authored dual geometric factor [m²/s].
    pub conductance_m2_s: f64,
}

/// Explicit accuracy, storage, work and roundoff acceptance limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceDiffusionControls {
    /// Positive upper substep duration [s]; controls splitting error, not stability.
    pub maximum_step_s: f64,
    pub maximum_substeps: usize,
    pub maximum_faces: usize,
    pub maximum_edges: usize,
    /// Absolute accepted change of total surfactant amount [mol].
    pub amount_tolerance_mol: f64,
    /// Absolute positive free-energy change admitted as numerical roundoff [J].
    pub energy_increase_tolerance_j: f64,
}

/// The accepted fixed-area material exchange and its separate free-energy loss.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceDiffusionReport {
    pub duration_s: f64,
    pub substeps: usize,
    pub pair_exchanges: usize,
    pub amount_before_mol: f64,
    pub amount_after_mol: f64,
    /// Compensated sum of face amount changes, not subtraction of large totals.
    pub amount_defect_mol: f64,
    pub free_energy_before_j: f64,
    pub free_energy_after_j: f64,
    /// Actual constitutive changes at fixed area, with the clean baseline removed.
    pub free_energy_change_j: f64,
    /// Minus the free-energy change, never an inferred thermal deposit. It is
    /// nonnegative in exact arithmetic; admitted negative roundoff is bounded
    /// by `energy_increase_tolerance_j` and is reported without clipping.
    pub free_energy_dissipation_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SurfaceDiffusionError {
    InvalidInput(&'static str),
    InvalidEdge { edge: usize },
    Constitutive { face: usize, error: SurfactantError },
    WorkBudget,
    Allocation,
    Unrepresentable,
    Cancelled,
    AmountBalance { defect_mol: f64 },
    FreeEnergyIncrease { change_j: f64 },
}
impl core::fmt::Display for SurfaceDiffusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "insoluble surface diffusion: {self:?}")
    }
}
impl std::error::Error for SurfaceDiffusionError {}

struct Exchange {
    faces: [usize; 2],
    area_fractions: [f64; 2],
    relaxation: f64,
}

/// Diffuse the SAME material-face amounts used by the capillary constitutive law.
///
/// The amounts change only after complete success. A common law is required:
/// concentration-driven flux need not dissipate free energy for heterogeneous
/// chemical affinities, capacities or temperatures. Use actual material face
/// areas at this operator stage; moving-surface composition belongs to the
/// caller's mechanical clock. Empty edges and zero duration are valid no-ops.
#[allow(clippy::too_many_arguments)]
pub fn diffuse_surface_surfactant(
    law: SurfactantLaw,
    areas_m2: &[f64],
    amounts_mol: &mut [f64],
    edges: &[SurfaceDiffusionEdge],
    duration_s: f64,
    controls: SurfaceDiffusionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SurfaceDiffusionReport, SurfaceDiffusionError> {
    poll(&mut cancelled)?;
    if areas_m2.is_empty() || areas_m2.len() != amounts_mol.len() {
        return Err(SurfaceDiffusionError::InvalidInput(
            "nonempty matching areas and amounts required",
        ));
    }
    if !duration_s.is_finite()
        || duration_s < 0.0
        || !controls.maximum_step_s.is_finite()
        || controls.maximum_step_s <= 0.0
        || controls.maximum_substeps == 0
        || !controls.amount_tolerance_mol.is_finite()
        || controls.amount_tolerance_mol < 0.0
        || !controls.energy_increase_tolerance_j.is_finite()
        || controls.energy_increase_tolerance_j < 0.0
    {
        return Err(SurfaceDiffusionError::InvalidInput(
            "finite duration, positive time/work limits and nonnegative tolerances required",
        ));
    }
    if areas_m2.len() > controls.maximum_faces || edges.len() > controls.maximum_edges {
        return Err(SurfaceDiffusionError::WorkBudget);
    }
    let mut amount_before = Sum::default();
    let mut energy_before = Sum::default();
    for (face, (&area, &amount)) in areas_m2.iter().zip(amounts_mol.iter()).enumerate() {
        poll_at(&mut cancelled, face)?;
        let state = SurfactantFace {
            law,
            amount_mol: amount,
        }
        .at_area(area)
        .map_err(|error| SurfaceDiffusionError::Constitutive { face, error })?;
        amount_before.add(amount);
        energy_before.add(state.free_energy_j);
    }
    let amount_before = finite(amount_before.value())?;
    let energy_before = finite(energy_before.value())?;
    let substeps = if duration_s == 0.0 || edges.is_empty() {
        0
    } else {
        let required = (duration_s / controls.maximum_step_s).ceil().max(1.0);
        if required > controls.maximum_substeps as f64 {
            return Err(SurfaceDiffusionError::WorkBudget);
        }
        let count = required as usize;
        if count == 0 || count > controls.maximum_substeps {
            return Err(SurfaceDiffusionError::WorkBudget);
        }
        count
    };
    let pair_exchanges = if substeps == 0 {
        0
    } else {
        edges
            .len()
            .checked_mul(2)
            .and_then(|n| n.checked_sub(1))
            .and_then(|n| n.checked_mul(substeps))
            .ok_or(SurfaceDiffusionError::WorkBudget)?
    };
    let step = if substeps == 0 {
        0.0
    } else {
        duration_s / substeps as f64
    };
    if substeps != 0 && step <= 0.0 {
        return Err(SurfaceDiffusionError::Unrepresentable);
    }
    let mut propagators = Vec::new();
    propagators
        .try_reserve_exact(edges.len())
        .map_err(|_| SurfaceDiffusionError::Allocation)?;
    for (index, edge) in edges.iter().enumerate() {
        poll_at(&mut cancelled, index)?;
        let [a, b] = edge.faces;
        if a == b
            || a >= areas_m2.len()
            || b >= areas_m2.len()
            || !edge.conductance_m2_s.is_finite()
            || edge.conductance_m2_s < 0.0
        {
            return Err(SurfaceDiffusionError::InvalidEdge { edge: index });
        }
        let left = areas_m2[a];
        let right = areas_m2[b];
        // Scaled fractions avoid overflowing A_i+A_j. A ratio that vanishes
        // cannot safely carry material into its small face, so refuse it.
        let ratio = left.min(right) / left.max(right);
        if ratio == 0.0 {
            return Err(SurfaceDiffusionError::Unrepresentable);
        }
        let small = ratio / (1.0 + ratio);
        let large = 1.0 / (1.0 + ratio);
        let area_fractions = if left <= right {
            [small, large]
        } else {
            [large, small]
        };
        let interval = if index + 1 == edges.len() {
            step
        } else {
            0.5 * step
        };
        let relaxation = if edge.conductance_m2_s == 0.0 || interval == 0.0 {
            0.0
        } else {
            let rate = finite(edge.conductance_m2_s / left + edge.conductance_m2_s / right)?;
            let exponent = rate * interval;
            // A finite rate times a positive interval may overflow only into
            // complete equilibration; expm1 handles this stiff limit exactly.
            if exponent == 0.0 {
                return Err(SurfaceDiffusionError::Unrepresentable);
            }
            -det::expm1(-exponent)
        };
        if !relaxation.is_finite() || !(0.0..=1.0).contains(&relaxation) {
            return Err(SurfaceDiffusionError::Unrepresentable);
        }
        propagators.push(Exchange {
            faces: edge.faces,
            area_fractions,
            relaxation,
        });
    }
    let mut candidate = Vec::new();
    candidate
        .try_reserve_exact(amounts_mol.len())
        .map_err(|_| SurfaceDiffusionError::Allocation)?;
    candidate.extend_from_slice(amounts_mol);
    for _ in 0..substeps {
        poll(&mut cancelled)?;
        for (index, exchange) in propagators.iter().enumerate() {
            poll_at(&mut cancelled, index)?;
            exchange.apply(areas_m2, &mut candidate)?;
        }
        for (index, exchange) in propagators[..propagators.len() - 1]
            .iter()
            .rev()
            .enumerate()
        {
            poll_at(&mut cancelled, index)?;
            exchange.apply(areas_m2, &mut candidate)?;
        }
    }
    let mut amount_after = Sum::default();
    let mut amount_defect = Sum::default();
    let mut energy_after = Sum::default();
    let mut energy_change = Sum::default();
    for (face, ((&area, &before), &after)) in areas_m2
        .iter()
        .zip(amounts_mol.iter())
        .zip(&candidate)
        .enumerate()
    {
        poll_at(&mut cancelled, face)?;
        let initial = SurfactantFace {
            law,
            amount_mol: before,
        };
        let state = SurfactantFace {
            law,
            amount_mol: after,
        }
        .at_area(area)
        .map_err(|error| SurfaceDiffusionError::Constitutive { face, error })?;
        let change = initial
            .energy_change_at_fixed_area(area, after)
            .map_err(|error| SurfaceDiffusionError::Constitutive { face, error })?;
        amount_after.add(after);
        amount_defect.add(after - before);
        energy_after.add(state.free_energy_j);
        energy_change.add(change);
    }
    let amount_after = finite(amount_after.value())?;
    let amount_defect = finite(amount_defect.value())?;
    let energy_after = finite(energy_after.value())?;
    let energy_change = finite(energy_change.value())?;
    if amount_defect.abs() > controls.amount_tolerance_mol {
        return Err(SurfaceDiffusionError::AmountBalance {
            defect_mol: amount_defect,
        });
    }
    if energy_change > controls.energy_increase_tolerance_j {
        return Err(SurfaceDiffusionError::FreeEnergyIncrease {
            change_j: energy_change,
        });
    }
    poll(&mut cancelled)?;
    amounts_mol.copy_from_slice(&candidate);
    Ok(SurfaceDiffusionReport {
        duration_s,
        substeps,
        pair_exchanges,
        amount_before_mol: amount_before,
        amount_after_mol: amount_after,
        amount_defect_mol: amount_defect,
        free_energy_before_j: energy_before,
        free_energy_after_j: energy_after,
        free_energy_change_j: energy_change,
        free_energy_dissipation_j: -energy_change,
    })
}

impl Exchange {
    fn apply(&self, areas: &[f64], amounts: &mut [f64]) -> Result<(), SurfaceDiffusionError> {
        let [i, j] = self.faces;
        let [left, right] = [amounts[i], amounts[j]];
        // Preserve an exactly resolved uniform concentration without accumulating
        // a cancellation residue in its mathematically zero transfer.
        if self.relaxation == 0.0 || left / areas[i] == right / areas[j] {
            return Ok(());
        }
        let [wi, wj] = self.area_fractions;
        let transfer = self.relaxation * (wi * right - wj * left);
        let after_left = finite(left + transfer)?;
        let after_right = finite(right - transfer)?;
        if after_left < 0.0 || after_right < 0.0 {
            return Err(SurfaceDiffusionError::Unrepresentable);
        }
        amounts[i] = after_left;
        amounts[j] = after_right;
        Ok(())
    }
}

fn finite(value: f64) -> Result<f64, SurfaceDiffusionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SurfaceDiffusionError::Unrepresentable)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SurfaceDiffusionError> {
    if cancelled() {
        Err(SurfaceDiffusionError::Cancelled)
    } else {
        Ok(())
    }
}
fn poll_at(
    cancelled: &mut impl FnMut() -> bool,
    index: usize,
) -> Result<(), SurfaceDiffusionError> {
    if index % 64 == 0 {
        poll(cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
