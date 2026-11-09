//! Simultaneous original Coulomb slip and own-hardness material removal.
//!
//! The prepared full mechanical response retains cross-cell and cross-axis
//! recoil. Each law evaluation starts from the original shear history; solver
//! sweeps never accumulate fictitious slip, elapsed time, wear or heat.
//! Friction sources retain the existing held-temperature/local-pressure-only
//! admission. Speed/path-dependent claims are not frozen at released motion.
//!
//! Normal and tangential updates remain split at beginning pressure. The caller
//! owns the supplied mechanical response, fixed tangent axes and physical
//! interval. This does not add normal/shear cross-compliance or a heat partition.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_math::det;
use fs_tribo::tangential::multiple::{
    MultipleTangentialError, MultipleTangentialStep, PreparedTangentialContacts,
    TangentialContactInput,
};
use fs_tribo::{TangentialSolveConfig, TangentialStep, TangentialStepInput, WearAdvance};

use super::super::{SourcePlasticSlidingError, SourcePlasticSlidingWear, add_compensated, cells};
use super::{
    HardnessSlidingContinuation, HardnessSlidingControl, HardnessSlidingError, SlidingState,
    SourceHistorySlidingWear,
};
use crate::normal_patch::plastic_source::history_field::HistoryHardnessContactError;

/// Fixed ordered cell registration and the complete response on those sites.
/// Matrix entries retain their original [m/N] meaning, including every block.
#[derive(Clone, Debug, PartialEq)]
pub struct HardnessSlidingCoupling {
    cells: Arc<[usize]>,
    response: Arc<PreparedTangentialContacts>,
}
impl HardnessSlidingCoupling {
    /// Register each original grid cell at most once, preserving caller order.
    /// Grid-range checks are made against the actual contact owner at each step.
    pub fn new(
        cell_indices: Vec<usize>,
        response: PreparedTangentialContacts,
        cx: &Cx<'_>,
    ) -> Result<Self, CoupledHardnessSlidingError> {
        check(cx)?;
        if cell_indices.len() != response.contact_count() {
            return Err(CoupledHardnessSlidingError::Invalid(
                "one original cell per prepared site",
            ));
        }
        let mut unique = BTreeSet::new();
        for &cell in &cell_indices {
            check(cx)?;
            if !unique.insert(cell) {
                return Err(CoupledHardnessSlidingError::Invalid(
                    "duplicate original contact cell",
                ));
            }
        }
        Ok(Self {
            cells: cell_indices.into(),
            response: Arc::new(response),
        })
    }
    /// Original grid cells in the response's unchanged site order.
    pub fn cell_indices(&self) -> &[usize] {
        &self.cells
    }
    /// Complete immutable physical matrix and native numerical policy.
    pub fn mechanical_response(&self) -> &PreparedTangentialContacts {
        &self.response
    }
}

/// Physical/source refusal never publishes partial shear, wear or normal history.
#[derive(Clone, Debug, PartialEq)]
pub enum CoupledHardnessSlidingError {
    Invalid(&'static str),
    /// Original hardness/contact/wear or selected Coulomb source refusal.
    Source(HardnessSlidingError),
    /// Full-matrix closure refusal; contact indices are registered site indices.
    Shear(MultipleTangentialError),
}
impl From<HardnessSlidingError> for CoupledHardnessSlidingError {
    fn from(error: HardnessSlidingError) -> Self {
        Self::Source(error)
    }
}
impl From<SourcePlasticSlidingError> for CoupledHardnessSlidingError {
    fn from(error: SourcePlasticSlidingError) -> Self {
        Self::Source(error.into())
    }
}
impl From<HistoryHardnessContactError> for CoupledHardnessSlidingError {
    fn from(error: HistoryHardnessContactError) -> Self {
        Self::Source(error.into())
    }
}
impl From<MultipleTangentialError> for CoupledHardnessSlidingError {
    fn from(error: MultipleTangentialError) -> Self {
        Self::Shear(error)
    }
}
impl core::fmt::Display for CoupledHardnessSlidingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled original-hardness sliding/wear: {self:?}")
    }
}
impl std::error::Error for CoupledHardnessSlidingError {}

/// Original source samples, full mechanical closure and actual slip, in fixed axes.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledHardnessSlidingReceipt {
    /// Shared exact registration, complete physical response and solver policy.
    pub coupling: HardnessSlidingCoupling,
    /// Original zero-traction motion at each registered site [m].
    pub released_increment_m: Vec<[f64; 2]>,
    /// Beginning normal pressures in registered site order [Pa].
    pub pressure_pa: Vec<f64>,
    pub duration_s: f64,
    pub friction_controls: TangentialSolveConfig,
    /// Native simultaneous response; candidate histories live in the successor.
    pub shear: MultipleTangentialStep<()>,
    /// Original selected T/p samples shared across identical loaded queries.
    pub friction_sources: Vec<ResolvedInterfaceStatePoint>,
    /// Per-site sample index; detached sites need no zero-pressure source claim.
    pub cell_source_indices: Vec<Option<usize>>,
    /// Actual unsigned wear distances on the complete original normal grid [m].
    /// Unregistered cells receive zero wear input; no zero-slip reaction is
    /// solved there.
    pub unsigned_slip_m: Vec<f64>,
    /// Sum p*A*actual_distance; distinct from shear dissipation [J].
    pub pressure_weighted_sliding_j: f64,
}

/// Complete coupled successor and distinct shear, wear and normal accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledHardnessSlidingAdvance<C> {
    pub source_identity: ContentHash,
    pub candidate: C,
    pub sliding: CoupledHardnessSlidingReceipt,
    pub wear: [WearAdvance; 2],
    pub recession_increment_m: [Vec<f64>; 2],
    pub normal: HardnessSlidingContinuation,
}

impl SourceHistorySlidingWear {
    /// Resolve full-matrix Coulomb stick/slip before removing either body.
    ///
    /// Only registered sites participate in this tangential solve. Omitted
    /// loaded cells supply no shear reaction here and retain their histories.
    /// Their influence is present only if already eliminated consistently into
    /// the supplied response and released motion. Zero wear input at an omitted
    /// cell does not establish that it physically sticks. This method does not
    /// claim to resolve every shear reaction on the normal grid.
    ///
    /// Pressure and original T/p coefficient samples are held during the native
    /// solve; actual solved irreversible slip drives each body's mapped H law.
    ///
    /// Wear then re-equilibrates the original normal owner at the supplied held
    /// approach, maintained load or finite spring. Refine physical intervals
    /// when this split changes pressure appreciably. The original source owner
    /// rejects extra speed/path axes; this method does not introduce a second
    /// nonlinear coefficient closure or infer thermal evolution from friction.
    pub fn advance_coupled(
        &self,
        coupling: &HardnessSlidingCoupling,
        released_increment_m: &[[f64; 2]],
        duration_s: f64,
        friction_controls: TangentialSolveConfig,
        control: HardnessSlidingControl,
        cx: &Cx<'_>,
    ) -> Result<CoupledHardnessSlidingAdvance<Self>, CoupledHardnessSlidingError> {
        let grid = self.wear.model().config().grid;
        let (sliding, receipt) = self.sliding.advance_coupled(
            self.wear.contact_state().pressure_pa(),
            grid.cell_width_m * grid.cell_depth_m,
            coupling,
            released_increment_m,
            duration_s,
            friction_controls,
            cx,
        )?;
        let wear = self
            .wear
            .advance_wear(&receipt.unsigned_slip_m, control, cx)?;
        check(cx)?;
        Ok(CoupledHardnessSlidingAdvance {
            source_identity: self.source_identity(),
            candidate: Self {
                wear: wear.candidate,
                sliding,
            },
            sliding: receipt,
            wear: wear.wear,
            recession_increment_m: wear.recession_increment_m,
            normal: HardnessSlidingContinuation {
                normal: wear.receipt,
                boundary: wear.control,
            },
        })
    }
}

impl SlidingState {
    #[allow(clippy::too_many_arguments)]
    fn advance_coupled(
        &self,
        pressure: &[f64],
        area: f64,
        coupling: &HardnessSlidingCoupling,
        released: &[[f64; 2]],
        duration_s: f64,
        controls: TangentialSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, CoupledHardnessSlidingReceipt), CoupledHardnessSlidingError> {
        check(cx)?;
        if released.len() != coupling.cells.len()
            || self.tangential.len() != pressure.len()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || released.iter().flatten().any(|x| !x.is_finite())
            || !self.friction.query_axes().is_empty()
        {
            return Err(CoupledHardnessSlidingError::Invalid(
                "finite registered motion, positive duration and held T/p-only friction required",
            ));
        }
        let mut sources = Vec::new();
        let mut laws = Vec::new();
        let mut sample_indices = BTreeMap::new();
        let mut source_indices = Vec::with_capacity(coupling.cells.len());
        let mut pressures = Vec::with_capacity(coupling.cells.len());
        let mut normal = Vec::with_capacity(coupling.cells.len());
        for &cell in coupling.cells.iter() {
            check(cx)?;
            let &p = pressure
                .get(cell)
                .ok_or(CoupledHardnessSlidingError::Invalid(
                    "registered cell lies outside the original contact grid",
                ))?;
            let error = |error| SourcePlasticSlidingError::Cell { index: cell, error };
            normal.push(SourcePlasticSlidingWear::cell_normal_force(p, area).map_err(error)?);
            pressures.push(p);
            let slot = if p == 0.0 {
                None
            } else {
                let key = if self.friction.axes().nominal_pressure.is_some() {
                    p.to_bits()
                } else {
                    0
                };
                let slot = if let Some(&slot) = sample_indices.get(&key) {
                    slot
                } else {
                    let (source, law) =
                        cells::coefficients_at_pressure(&self.friction, self.temperature_k, p)
                            .map_err(error)?;
                    let slot = sources.len();
                    sources.push(source);
                    laws.push(law);
                    sample_indices.insert(key, slot);
                    slot
                };
                Some(slot)
            };
            source_indices.push(slot);
        }
        let inputs: Vec<_> = coupling
            .cells
            .iter()
            .enumerate()
            .map(|(site, &cell)| TangentialContactInput {
                law: source_indices[site].map_or(&self.initial_law, |slot| &laws[slot]),
                state: &self.tangential[cell],
                input: TangentialStepInput {
                    normal_n: normal[site],
                    slip_increment_m: released[site],
                    dt_s: duration_s,
                    contact_temperature_k: Some(self.temperature_k),
                    nominal_pressure_pa: self
                        .friction
                        .axes()
                        .nominal_pressure
                        .as_ref()
                        .map(|_| pressures[site]),
                    compliance_m_per_n: 0.0,
                    solve: controls,
                },
            })
            .collect();
        let native = coupling
            .response
            .advance(&inputs, || cx.checkpoint().is_err())?;
        let mut tangential = self.tangential.clone();
        let mut unsigned_slip_m = vec![0.0; pressure.len()];
        let mut contacts = Vec::with_capacity(coupling.cells.len());
        let mut measure = 0.0;
        let mut correction = 0.0;
        for (site, step) in native.contacts.into_iter().enumerate() {
            check(cx)?;
            let cell = coupling.cells[site];
            let distance = det::hypot(step.slip_increment_m[0], step.slip_increment_m[1]);
            let local_measure = normal[site] * distance;
            if !distance.is_finite() || !local_measure.is_finite() {
                return Err(CoupledHardnessSlidingError::Invalid(
                    "finite actual local sliding",
                ));
            }
            unsigned_slip_m[cell] = distance;
            add_compensated(&mut measure, &mut correction, local_measure);
            tangential[cell] = step.state;
            contacts.push(TangentialStep {
                state: (),
                traction_n: step.traction_n,
                slip_increment_m: step.slip_increment_m,
                regime: step.regime,
                stored_change_j: step.stored_change_j,
                dissipated_j: step.dissipated_j,
                closure_residual_j: step.closure_residual_j,
                compatibility_residual_m: step.compatibility_residual_m,
                iterations: step.iterations,
            });
        }
        if !measure.is_finite() {
            return Err(CoupledHardnessSlidingError::Invalid(
                "finite total pressure-sliding measure",
            ));
        }
        check(cx)?;
        Ok((
            Self {
                friction: Arc::clone(&self.friction),
                initial_law: self.initial_law.clone(),
                temperature_k: self.temperature_k,
                tangential,
            },
            CoupledHardnessSlidingReceipt {
                coupling: coupling.clone(),
                released_increment_m: released.to_vec(),
                pressure_pa: pressures,
                duration_s,
                friction_controls: controls,
                shear: MultipleTangentialStep {
                    contacts,
                    contact_work_residuals_j: native.contact_work_residuals_j,
                    sweeps: native.sweeps,
                    law_evaluations: native.law_evaluations,
                    maximum_compatibility_residual_m: native.maximum_compatibility_residual_m,
                    stored_change_j: native.stored_change_j,
                    dissipated_j: native.dissipated_j,
                    body_work_residual_j: native.body_work_residual_j,
                },
                friction_sources: sources,
                cell_source_indices: source_indices,
                unsigned_slip_m,
                pressure_weighted_sliding_j: measure,
            },
        ))
    }
}

fn check(cx: &Cx<'_>) -> Result<(), CoupledHardnessSlidingError> {
    cx.checkpoint()
        .map_err(|_| MultipleTangentialError::Cancelled.into())
}
