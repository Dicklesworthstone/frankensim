//! Re-solve the original elastic load path before its first tensile crossing.
//!
//! The first lower event bracket is a proposed load fraction, not an accepted
//! FE state. This owner interpolates the retained physical loads and supports,
//! solves the original elastic problem again, and accepts only a native probe
//! with no tensile crossing. Bounded backoff may retreat from the bracket when
//! solve or onset resolution refuses. The reached fraction is reported exactly.
//!
//! Original target loads, solver controls and physical direction hints survive
//! in the remainder. Applied forces at constrained DOFs remain applied forces;
//! they are never reconstructed from support reactions. This is progress toward
//! first onset in a quasistatic load interval, not a physical timestep or an
//! event-aware damaged-solid solve after stress redistribution.

use fs_exec::Cx;
use fs_material::brittle::initiation::{TensileInitiationError, TensileInitiationSettings};
use fs_material::brittle::source::initiation::ResolvedTensileInitiationError;

use super::{TensileElasticTetSolid, TensileTetInitiationError, TensileTetOnset, TensileTetProbe};
use crate::brittle::solid::{BrittleTetError, checkpoint};
use crate::plastic3::{IncrementSettings, SmallStrainTetError};

/// Bounded native replay and minimum representable physical load progress.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensileFirstEventControls {
    /// Native solves, including rejected candidates; in 1..=128.
    pub maximum_replays: usize,
    /// Smallest admitted fraction of the ORIGINAL requested load interval.
    pub minimum_fraction: f64,
}

/// A native elastic endpoint strictly before the first resolved tensile event.
#[derive(Debug)]
pub struct TensileFirstEventStep<'a> {
    pub solid: TensileElasticTetSolid<'a>,
    /// Actual load fraction reached, which may precede the original bracket.
    pub reached_fraction: f64,
    pub remaining_fraction: f64,
    pub replays: usize,
    /// Exact reversible elastic work, equal to the original storage change.
    /// This is not an endpoint external-work quadrature or numerical heat.
    pub elastic_work_j: f64,
    /// Original target support entries, including their order and exact values.
    pub prescribed_end_m: Vec<(usize, f64)>,
    /// Original target force entries, including duplicates and support loads.
    pub nodal_forces_end_n: Vec<(usize, f64)>,
    /// Canonical remaining force interval per DOF: accepted start and original
    /// target. Duplicate force entries are added using native ordering.
    pub force_intervals_n: Vec<(usize, [f64; 2])>,
    solver_settings: IncrementSettings,
    initiation_settings: TensileInitiationSettings,
    direction_hints: Option<Vec<Option<[f64; 3]>>>,
}

impl<'a> TensileFirstEventStep<'a> {
    /// Probe the exact original target from this accepted elastic endpoint.
    /// This does not accept the predictor or choose post-fracture event order.
    pub fn probe_remainder(
        &self,
        cx: &Cx<'_>,
    ) -> Result<TensileTetProbe<'a>, TensileTetInitiationError> {
        self.solid.probe(
            &self.prescribed_end_m,
            &self.nodal_forces_end_n,
            self.solver_settings,
            self.initiation_settings,
            self.direction_hints.as_deref(),
            cx,
        )
    }
}

/// Every refusal preserves both the original accepted history and prediction.
#[derive(Clone, Debug, PartialEq)]
pub enum TensileFirstEventError {
    Invalid(&'static str),
    Allocation,
    Native(TensileTetInitiationError),
    /// No positive admissible load fraction or displacement update is resolved.
    NoProgress {
        fraction: f64,
    },
    MinimumFraction {
        attempted: f64,
        minimum: f64,
    },
    /// None denotes a still-crossing native probe; Some retains a refinable
    /// numerical refusal from the last attempt.
    ReplayLimit {
        replays: usize,
        last_fraction: f64,
        last_refusal: Option<TensileTetInitiationError>,
    },
}

impl core::fmt::Display for TensileFirstEventError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "first tensile-event approach refused: {self:?}")
    }
}
impl std::error::Error for TensileFirstEventError {}
impl From<TensileTetInitiationError> for TensileFirstEventError {
    fn from(error: TensileTetInitiationError) -> Self {
        Self::Native(error)
    }
}
impl From<BrittleTetError> for TensileFirstEventError {
    fn from(error: BrittleTetError) -> Self {
        Self::Native(error.into())
    }
}

impl<'a> TensileTetOnset<'a> {
    /// Advance to a checked elastic point before the FIRST predicted crossing.
    ///
    /// Original accepted dead forces and exact endpoint supports/forces were
    /// retained by probe. Affine interpolation uses those actual inputs, never
    /// reaction forces or a manufactured displacement field. The lower bracket
    /// is attempted first; each refinable refusal halves the original interval
    /// fraction. Every candidate passes the existing native solver and tensile
    /// detector. Cancellation, material/domain errors and invalid inputs are
    /// terminal. No zero-progress state or over-strength predictor is returned.
    pub fn advance_before_first_event(
        &self,
        controls: TensileFirstEventControls,
        cx: &Cx<'_>,
    ) -> Result<TensileFirstEventStep<'a>, TensileFirstEventError> {
        checkpoint(cx)?;
        if !(1..=128).contains(&controls.maximum_replays)
            || !controls.minimum_fraction.is_finite()
            || controls.minimum_fraction <= 0.0
            || controls.minimum_fraction >= 1.0
        {
            return Err(TensileFirstEventError::Invalid(
                "1..=128 replays and a finite minimum fraction in (0, 1) required",
            ));
        }
        let mut fraction = self.first_fraction[0];
        admit_fraction(fraction, controls.minimum_fraction)?;
        let old_u = self.initial.displacement_m.as_flattened();
        let mut force_intervals =
            if self.initial_forces_n.is_empty() && self.nodal_forces_end_n.is_empty() {
                Vec::new()
            } else {
                let mut start = buffer(old_u.len())?;
                let mut end = buffer(old_u.len())?;
                start.resize(old_u.len(), 0.0);
                end.resize(old_u.len(), 0.0);
                accumulate(&mut start, &self.initial_forces_n, cx)?;
                accumulate(&mut end, &self.nodal_forces_end_n, cx)?;
                let capacity = old_u.len().min(
                    self.initial_forces_n
                        .len()
                        .saturating_add(self.nodal_forces_end_n.len()),
                );
                let mut intervals = buffer(capacity)?;
                for (dof, (&a, &b)) in start.iter().zip(&end).enumerate() {
                    if dof % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    if a != 0.0 || b != 0.0 {
                        intervals.push((dof, [a, b]));
                    }
                }
                intervals
            };
        let mut forces = buffer(force_intervals.len())?;
        for &(dof, [value, _]) in &force_intervals {
            forces.push((dof, value));
        }
        let mut prescribed = buffer(self.prescribed_end_m.len())?;
        let mut supports = buffer(self.prescribed_end_m.len())?;
        for &(dof, _) in &self.prescribed_end_m {
            prescribed.push((dof, old_u[dof]));
            supports.push(dof);
        }
        supports.sort_unstable();
        let initial = TensileElasticTetSolid {
            geometry: self.geometry,
            materials: self.materials,
            budget: self.budget,
            solution: self.initial.clone(),
            stored_energy_j: self.elastic_storage_j,
            supports: Some(supports),
            nodal_forces_n: self.initial_forces_n.clone(),
        };
        let mut replays = 0;
        loop {
            checkpoint(cx)?;
            admit_fraction(fraction, controls.minimum_fraction)?;
            for (index, (entry, &(dof, target))) in prescribed
                .iter_mut()
                .zip(&self.prescribed_end_m)
                .enumerate()
            {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                entry.1 = interpolate(old_u[dof], target, fraction);
            }
            for (index, (entry, &(_, [a, b]))) in
                forces.iter_mut().zip(&force_intervals).enumerate()
            {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                entry.1 = interpolate(a, b, fraction);
            }
            replays += 1;
            let last_refusal = match initial.probe(
                &prescribed,
                &forces,
                self.solver_settings,
                self.settings,
                self.direction_hints.as_deref(),
                cx,
            ) {
                Ok(TensileTetProbe::Elastic(solid)) => {
                    if solid.solution.displacement_m == self.initial.displacement_m {
                        return Err(TensileFirstEventError::NoProgress { fraction });
                    }
                    let work = solid.stored_energy_j - self.elastic_storage_j;
                    if !work.is_finite() {
                        return Err(TensileFirstEventError::Invalid(
                            "elastic storage change is unrepresentable",
                        ));
                    }
                    for (entry, &(_, force)) in force_intervals.iter_mut().zip(&forces) {
                        entry.1[0] = force;
                    }
                    checkpoint(cx)?;
                    return Ok(TensileFirstEventStep {
                        solid,
                        reached_fraction: fraction,
                        remaining_fraction: 1.0 - fraction,
                        replays,
                        elastic_work_j: work,
                        prescribed_end_m: self.prescribed_end_m.clone(),
                        nodal_forces_end_n: self.nodal_forces_end_n.clone(),
                        force_intervals_n: force_intervals,
                        solver_settings: self.solver_settings,
                        initiation_settings: self.settings,
                        direction_hints: self.direction_hints.clone(),
                    });
                }
                Ok(TensileTetProbe::Onset(_)) => None,
                Err(error) if retryable(&error) => Some(error),
                Err(error) => return Err(error.into()),
            };
            if replays == controls.maximum_replays {
                return Err(TensileFirstEventError::ReplayLimit {
                    replays,
                    last_fraction: fraction,
                    last_refusal,
                });
            }
            let smaller = 0.5 * fraction;
            if smaller <= 0.0 || smaller == fraction {
                return Err(TensileFirstEventError::NoProgress { fraction: smaller });
            }
            fraction = smaller;
        }
    }
}

fn admit_fraction(fraction: f64, minimum: f64) -> Result<(), TensileFirstEventError> {
    if !fraction.is_finite() || fraction <= 0.0 || fraction >= 1.0 {
        return Err(TensileFirstEventError::NoProgress { fraction });
    }
    if fraction < minimum {
        return Err(TensileFirstEventError::MinimumFraction {
            attempted: fraction,
            minimum,
        });
    }
    Ok(())
}
fn buffer<T>(capacity: usize) -> Result<Vec<T>, TensileFirstEventError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| TensileFirstEventError::Allocation)?;
    Ok(values)
}
fn accumulate(
    total: &mut [f64],
    entries: &[(usize, f64)],
    cx: &Cx<'_>,
) -> Result<(), TensileFirstEventError> {
    for (index, &(dof, value)) in entries.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        let Some(sum) = total.get_mut(dof) else {
            return Err(TensileFirstEventError::Invalid(
                "retained force DOF outside mesh",
            ));
        };
        *sum += value;
        if !sum.is_finite() {
            return Err(TensileFirstEventError::Invalid(
                "retained applied-force sum must be finite",
            ));
        }
    }
    Ok(())
}
fn interpolate(start: f64, end: f64, fraction: f64) -> f64 {
    fraction.mul_add(end, (1.0 - fraction) * start)
}
fn retryable(error: &TensileTetInitiationError) -> bool {
    matches!(
        error,
        TensileTetInitiationError::Solid(BrittleTetError::Solve(
            SmallStrainTetError::NewtonStalled { .. }
                | SmallStrainTetError::LinearSolveFailed { .. }
        ))
    ) || matches!(
        error,
        TensileTetInitiationError::Element {
            source: ResolvedTensileInitiationError::Native(
                TensileInitiationError::UnresolvedSpectrum | TensileInitiationError::OnsetBudget
            ),
            ..
        }
    )
}
