//! Surface recession and fixed-grid material refresh with explicit normal work.
//!
//! Recession increases the relative undeformed separation and remains separate
//! from residual plastic indentation. Every continuation preserves the original
//! gap datum. Recession and material refresh normally hold the approach fixed;
//! load-controlled continuation instead accounts for constant-resultant normal
//! work. A finite-stiffness drive supplies its own normal-work integral.
//! Signed parameter work is the endpoint estimate
//! `delta(elastic energy) + incremental plastic work - normal work`.
//! Condition/geometry paths require refinement when intermediate equilibria
//! matter. No part of this bookkeeping is automatically converted to heat.

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::Cx;

use super::{
    FiniteGapContactError, FiniteGapPlasticAdvance, FiniteGapPlasticError, FiniteGapPlasticModel,
    FiniteGapPlasticReceipt, FiniteGapPlasticState, PlasticStepMetrics, box_residual, checkpoint,
    invalid,
};

/// Physical change advanced from an accepted contact history.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum FiniteGapPlasticContinuationKind {
    /// Nonnegative relative removal depth on the original contact grid.
    SurfaceRecession = 0,
    /// New indentation hardness and/or reduced modulus on the original grid.
    MaterialRefresh = 1,
    /// Relative surface recession with the original compressive resultant maintained.
    SurfaceRecessionAtLoad = 2,
    /// New material parameters with the original compressive resultant maintained.
    MaterialRefreshAtLoad = 3,
    /// Relative surface recession at a fixed finite-stiffness drive position.
    SurfaceRecessionWithDrive = 4,
    /// New material parameters at a fixed finite-stiffness drive position.
    MaterialRefreshWithDrive = 5,
}

/// Complete candidate from an explicitly authorized physical continuation.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticContinuation {
    /// Successor history; adopt only after accepting the complete coupled step.
    pub candidate: FiniteGapPlasticState,
    /// Source/model continuity, pressure response, and signed work bookkeeping.
    pub receipt: FiniteGapPlasticContinuationReceipt,
}

/// Continuation response with explicit predecessor/successor owners and normal work.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticContinuationReceipt {
    /// Content identity of the complete continuation and its physical response.
    pub identity: ContentHash,
    /// Parameter changed by this continuation.
    pub kind: FiniteGapPlasticContinuationKind,
    /// Compiled physical model that owns the accepted predecessor.
    pub predecessor_model_identity: ContentHash,
    /// Compiled physical model that owns the candidate.
    pub successor_model_identity: ContentHash,
    /// Exact accepted history from which this continuation starts.
    pub predecessor_state_identity: ContentHash,
    /// Exact candidate history after re-equilibration.
    pub successor_state_identity: ContentHash,
    /// Accepted predecessor approach relative to the original gap datum [m].
    pub previous_approach_m: f64,
    /// Candidate approach relative to the same original gap datum [m].
    pub approach_m: f64,
    /// Work of the controlled normal resultant during this continuation [J].
    ///
    /// Zero for held-approach operations. For maintained-load changes this is
    /// the prescribed constant force times the rigid approach change; the
    /// separate load receipt retains equilibrium force tolerances/residuals.
    /// At a held loading-device position, this is energy released by its
    /// spring. The stationary actuator performs no external work.
    pub normal_work_j: f64,
    /// Elastic energy before changing geometry or material parameters [J].
    pub previous_elastic_energy_j: f64,
    /// Elastic energy of the re-equilibrated candidate [J].
    pub elastic_energy_j: f64,
    /// Signed change in reversible elastic storage [J].
    pub elastic_energy_change_j: f64,
    /// Signed endpoint estimate `delta(U) + delta(W_plastic) - normal_work` [J].
    ///
    /// Negative values represent energy leaving through changed parameters;
    /// they need not imply a decrease in elastic storage.
    /// This is not wear heat, thermal heat capacity, or a complete external
    /// work integral along an unresolved condition/geometry path.
    pub parameter_work_j: f64,
    /// Representable increase in accumulated relative recession at each cell [m].
    ///
    /// Material refresh reports zeros. Requested positive increments that
    /// cannot change stored recession are refused rather than silently lost.
    pub recession_increment_m: Vec<f64>,
    /// Original native pressure, moments, elastic energy, and plastic response.
    /// Its predecessor identity is the accepted history, including on refresh.
    pub response: FiniteGapPlasticReceipt,
}

struct ContinuationAccounting {
    kind: FiniteGapPlasticContinuationKind,
    previous_elastic_energy_j: f64,
    previous_opening_m: Option<Vec<f64>>,
    normal_work_j: f64,
    recession_increment_m: Vec<f64>,
}

impl FiniteGapPlasticModel {
    /// Remove relative surface depth at the held approach and re-equilibrate.
    ///
    /// Every increment must be finite and nonnegative, on the original grid.
    /// Recession and plastic indentation remain independent histories. Loss of
    /// contact unloads elastic storage; pressure redistribution can also yield
    /// other cells. The accepted history is unchanged on every refusal.
    pub fn recede(
        &self,
        before: &FiniteGapPlasticState,
        relative_recession_increment_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticContinuation, FiniteGapPlasticError> {
        checkpoint(cx)?;
        self.check_state_owner(before)?;
        if relative_recession_increment_m.len() != self.gap_m.len()
            || relative_recession_increment_m
                .iter()
                .any(|increment| !increment.is_finite() || *increment < 0.0)
        {
            return Err(invalid("finite nonnegative surface recession field"));
        }
        let mut recession = before.recession_m.clone();
        let mut actual_increment = vec![0.0; recession.len()];
        for (i, increment) in relative_recession_increment_m.iter().enumerate() {
            recession[i] += increment;
            if !recession[i].is_finite() {
                return Err(invalid("finite accumulated surface recession"));
            }
            actual_increment[i] = recession[i] - before.recession_m[i];
            if *increment > 0.0 && actual_increment[i] == 0.0 {
                return Err(invalid("representable surface recession increment"));
            }
        }
        let (previous_elastic_energy_j, old_opening) =
            self.stored_energy_and_opening(before, cx)?;
        let (advance, previous_opening_m) = if actual_increment.iter().all(|value| *value == 0.0) {
            (self.unchanged_advance(before, old_opening, cx)?, None)
        } else {
            (
                self.advance_with_history(
                    before,
                    before.approach_m,
                    &recession,
                    &before.pressure_pa,
                    cx,
                )?,
                Some(old_opening),
            )
        };
        self.finish_continuation(
            self,
            before,
            advance,
            ContinuationAccounting {
                kind: FiniteGapPlasticContinuationKind::SurfaceRecession,
                previous_elastic_energy_j,
                previous_opening_m,
                normal_work_j: 0.0,
                recession_increment_m: actual_increment,
            },
            cx,
        )
    }

    /// Re-equilibrate an explicitly owned predecessor under this successor model.
    ///
    /// Only reduced modulus and the hardness field may change. The original
    /// grid, gap samples, datum, and numerical controls must match exactly.
    /// Plastic depth, recession, approach, and cumulative prior plastic work
    /// survive the transfer. A softer cap may cause additional plastic flow.
    ///
    /// This method checks the predecessor against `previous_model`; ordinary
    /// `advance` still refuses foreign histories. Identical model inputs return
    /// an identical physical state, without roundoff-induced plastic evolution.
    pub fn refresh_from(
        &self,
        previous_model: &Self,
        before: &FiniteGapPlasticState,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticContinuation, FiniteGapPlasticError> {
        checkpoint(cx)?;
        previous_model.check_state_owner(before)?;
        self.check_refresh_geometry(previous_model)?;
        let (previous_elastic_energy_j, old_opening) =
            previous_model.stored_energy_and_opening(before, cx)?;
        let (advance, previous_opening_m) = if self.identity == previous_model.identity {
            (self.unchanged_advance(before, old_opening, cx)?, None)
        } else {
            // Changing hardness can invalidate the old pressure bound. This
            // feasible warm guess is not a published or accepted intermediate.
            let pressure_guess: Vec<f64> = before
                .pressure_pa
                .iter()
                .zip(&self.config.hardness_pa)
                .map(|(pressure, hardness)| pressure.min(*hardness))
                .collect();
            (
                self.advance_with_history(
                    before,
                    before.approach_m,
                    &before.recession_m,
                    &pressure_guess,
                    cx,
                )?,
                Some(old_opening),
            )
        };
        self.finish_continuation(
            previous_model,
            before,
            advance,
            ContinuationAccounting {
                kind: FiniteGapPlasticContinuationKind::MaterialRefresh,
                previous_elastic_energy_j,
                previous_opening_m,
                normal_work_j: 0.0,
                recession_increment_m: vec![0.0; self.gap_m.len()],
            },
            cx,
        )
    }

    pub(super) fn check_refresh_geometry(
        &self,
        previous: &Self,
    ) -> Result<(), FiniteGapPlasticError> {
        let a = &self.config;
        let b = &previous.config;
        if self.linear_solver() != previous.linear_solver()
            || a.grid != b.grid
            || a.maximum_active_set_iterations != b.maximum_active_set_iterations
            || a.complementarity_tolerance_m.to_bits() != b.complementarity_tolerance_m.to_bits()
            || a.boundary_clearance_cells != b.boundary_clearance_cells
            || a.undeformed_gap_m.len() != b.undeformed_gap_m.len()
            || a.undeformed_gap_m
                .iter()
                .zip(&b.undeformed_gap_m)
                .any(|(left, right)| left.to_bits() != right.to_bits())
            || self.gap_datum_m.to_bits() != previous.gap_datum_m.to_bits()
        {
            return Err(invalid(
                "material refresh must retain original grid, gap, datum and solver controls",
            ));
        }
        Ok(())
    }

    pub(super) fn stored_energy_and_opening(
        &self,
        state: &FiniteGapPlasticState,
        cx: &Cx<'_>,
    ) -> Result<(f64, Vec<f64>), FiniteGapPlasticError> {
        let opening = self.apply_compliance(&state.pressure_pa, cx)?;
        let area = self.config.grid.cell_width_m * self.config.grid.cell_depth_m;
        // Keep the same reduction order as the original normal receipt.
        let energy = 0.5
            * area
            * state
                .pressure_pa
                .iter()
                .zip(&opening)
                .map(|(pressure, displacement)| pressure * displacement)
                .sum::<f64>();
        if !energy.is_finite() {
            return Err(invalid("finite predecessor elastic energy"));
        }
        Ok((energy, opening))
    }

    pub(super) fn unchanged_advance(
        &self,
        before: &FiniteGapPlasticState,
        opening: Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticAdvance, FiniteGapPlasticError> {
        self.check_state_owner(before)?;
        let target: Vec<f64> = self
            .gap_m
            .iter()
            .zip(&before.plastic_depth_m)
            .zip(&before.recession_m)
            .map(|((gap, plastic), recession)| before.approach_m - gap - plastic - recession)
            .collect();
        let mut physical_residual = 0.0_f64;
        for i in 0..target.len() {
            let separation = opening[i] - target[i];
            physical_residual = physical_residual.max(if before.pressure_pa[i] > 0.0 {
                separation.abs()
            } else {
                (-separation).max(0.0)
            });
        }
        if !physical_residual.is_finite()
            || physical_residual > self.config.complementarity_tolerance_m
        {
            return Err(FiniteGapContactError::ActiveSetDidNotConverge {
                iterations: 0,
                residual_m: physical_residual,
            }
            .into());
        }
        let kkt_residual_m = box_residual(
            &before.pressure_pa,
            &opening,
            &target,
            &self.config.hardness_pa,
        );
        let receipt = self.receipt_for_candidate(
            before.identity,
            before,
            opening,
            PlasticStepMetrics {
                iterations: 0,
                physical_residual_m: physical_residual,
                kkt_residual_m,
                plastic_increment_m: vec![0.0; target.len()],
                plastic_work_j: 0.0,
                yielded_cells: 0,
            },
            cx,
        )?;
        Ok(FiniteGapPlasticAdvance {
            candidate: before.clone(),
            receipt,
        })
    }

    pub(super) fn finish_load_recession(
        &self,
        before: &FiniteGapPlasticState,
        advance: FiniteGapPlasticAdvance,
        actual_recession_increment_m: Vec<f64>,
        maintained_force_n: f64,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticContinuation, FiniteGapPlasticError> {
        self.check_state_owner(before)?;
        let (previous_elastic_energy_j, opening) = self.stored_energy_and_opening(before, cx)?;
        let normal_work_j = maintained_force_n * (advance.candidate.approach_m - before.approach_m);
        if !normal_work_j.is_finite() {
            return Err(invalid("finite maintained-load normal work"));
        }
        self.finish_continuation(
            self,
            before,
            advance,
            ContinuationAccounting {
                kind: FiniteGapPlasticContinuationKind::SurfaceRecessionAtLoad,
                previous_elastic_energy_j,
                previous_opening_m: Some(opening),
                normal_work_j,
                recession_increment_m: actual_recession_increment_m,
            },
            cx,
        )
    }

    pub(super) fn finish_load_refresh(
        &self,
        previous_model: &Self,
        before: &FiniteGapPlasticState,
        advance: FiniteGapPlasticAdvance,
        maintained_force_n: f64,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticContinuation, FiniteGapPlasticError> {
        // The shared load oracle already admitted the original geometry and
        // both owners before constructing this successor response.
        previous_model.check_state_owner(before)?;
        let (previous_elastic_energy_j, opening) =
            previous_model.stored_energy_and_opening(before, cx)?;
        let normal_work_j = maintained_force_n * (advance.candidate.approach_m - before.approach_m);
        if !normal_work_j.is_finite() {
            return Err(invalid("finite maintained-load normal work"));
        }
        self.finish_continuation(
            previous_model,
            before,
            advance,
            ContinuationAccounting {
                kind: FiniteGapPlasticContinuationKind::MaterialRefreshAtLoad,
                previous_elastic_energy_j,
                previous_opening_m: Some(opening),
                normal_work_j,
                recession_increment_m: vec![0.0; self.gap_m.len()],
            },
            cx,
        )
    }

    pub(super) fn finish_drive_continuation(
        &self,
        previous_model: &Self,
        before: &FiniteGapPlasticState,
        advance: FiniteGapPlasticAdvance,
        actual_recession_increment_m: Option<Vec<f64>>,
        normal_work_j: f64,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticContinuation, FiniteGapPlasticError> {
        previous_model.check_state_owner(before)?;
        if !normal_work_j.is_finite() {
            return Err(invalid("finite drive normal work"));
        }
        let (previous_elastic_energy_j, opening) =
            previous_model.stored_energy_and_opening(before, cx)?;
        let (kind, recession_increment_m) = match actual_recession_increment_m {
            Some(increment) => (
                FiniteGapPlasticContinuationKind::SurfaceRecessionWithDrive,
                increment,
            ),
            None => (
                FiniteGapPlasticContinuationKind::MaterialRefreshWithDrive,
                vec![0.0; self.gap_m.len()],
            ),
        };
        self.finish_continuation(
            previous_model,
            before,
            advance,
            ContinuationAccounting {
                kind,
                previous_elastic_energy_j,
                previous_opening_m: Some(opening),
                normal_work_j,
                recession_increment_m,
            },
            cx,
        )
    }

    fn finish_continuation(
        &self,
        previous_model: &Self,
        before: &FiniteGapPlasticState,
        advance: FiniteGapPlasticAdvance,
        accounting: ContinuationAccounting,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticContinuation, FiniteGapPlasticError> {
        let elastic_energy_j = advance.receipt.normal.reversible_energy_j;
        let area = self.config.grid.cell_width_m * self.config.grid.cell_depth_m;
        let elastic_energy_change_j = match &accounting.previous_opening_m {
            Some(previous_opening) => compensated_energy_change(
                &before.pressure_pa,
                previous_opening,
                &advance.candidate.pressure_pa,
                &advance.receipt.normal.elastic_opening_m,
                area,
            ),
            None => 0.0,
        };
        let parameter_work_j =
            elastic_energy_change_j + advance.receipt.plastic_work_j - accounting.normal_work_j;
        if !elastic_energy_change_j.is_finite() || !parameter_work_j.is_finite() {
            return Err(invalid("finite continuation energy accounting"));
        }
        let mut bytes = Vec::with_capacity(240 + 8 * accounting.recession_increment_m.len());
        bytes.push(accounting.kind as u8);
        for identity in [
            previous_model.identity,
            self.identity,
            before.identity,
            advance.candidate.identity,
            advance.receipt.identity,
        ] {
            bytes.extend_from_slice(identity.as_bytes());
        }
        for value in [
            before.approach_m,
            advance.candidate.approach_m,
            accounting.normal_work_j,
            accounting.previous_elastic_energy_j,
            elastic_energy_j,
            elastic_energy_change_j,
            parameter_work_j,
        ] {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        for value in &accounting.recession_increment_m {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        let identity = hash_domain(
            "org.frankensim.fs-contact.finite-gap-plastic-continuation.v2",
            &bytes,
        );
        checkpoint(cx)?;
        Ok(FiniteGapPlasticContinuation {
            receipt: FiniteGapPlasticContinuationReceipt {
                identity,
                kind: accounting.kind,
                predecessor_model_identity: previous_model.identity,
                successor_model_identity: self.identity,
                predecessor_state_identity: before.identity,
                successor_state_identity: advance.candidate.identity,
                previous_approach_m: before.approach_m,
                approach_m: advance.candidate.approach_m,
                normal_work_j: accounting.normal_work_j,
                previous_elastic_energy_j: accounting.previous_elastic_energy_j,
                elastic_energy_j,
                elastic_energy_change_j,
                parameter_work_j,
                recession_increment_m: accounting.recession_increment_m,
                response: advance.receipt,
            },
            candidate: advance.candidate,
        })
    }
}

// Difference the two endpoint products before reducing their common storage.
// The symmetric identity is exact in real arithmetic; Neumaier compensation
// retains small cell contributions when pressure redistribution cancels terms.
fn compensated_energy_change(
    previous_pressure: &[f64],
    previous_opening: &[f64],
    pressure: &[f64],
    opening: &[f64],
    area: f64,
) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for i in 0..pressure.len() {
        let pressure_change = pressure[i] - previous_pressure[i];
        let opening_change = opening[i] - previous_opening[i];
        let mean_pressure = 0.5 * pressure[i] + 0.5 * previous_pressure[i];
        let mean_opening = 0.5 * opening[i] + 0.5 * previous_opening[i];
        for term in [
            pressure_change * mean_opening,
            opening_change * mean_pressure,
        ] {
            let next = sum + term;
            correction += if sum.abs() >= term.abs() {
                (sum - next) + term
            } else {
                (term - next) + sum
            };
            sum = next;
        }
    }
    0.5 * area * (sum + correction)
}

#[cfg(test)]
mod tests {
    use super::compensated_energy_change;

    #[test]
    fn small_cell_storage_change_survives_a_large_unchanged_energy() {
        // Two independently specified endpoint products differ by exactly one.
        // Subtracting separately rounded totals instead returns two here.
        let change =
            compensated_energy_change(&[1e16, 1.0], &[1.0, 1.0], &[1e16, 2.0], &[1.0, 1.0], 2.0);
        assert_eq!(change, 1.0);
    }
}
