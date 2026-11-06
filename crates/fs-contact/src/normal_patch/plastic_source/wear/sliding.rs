//! Source Coulomb sliding on the actual hardness-limited pressure field.
//!
//! Each original grid cell receives its beginning pressure times its area as
//! normal force. The existing tangential owner solves `du = du_free - C F`;
//! only the length of that actual irreversible slip drives original Archard
//! removal. Sticking therefore removes no material. Recession then changes
//! the normal equilibrium under the chosen approach, force, or loading spring.
//!
//! This is a first-order beginning-pressure split. Refine the physical motion
//! when wear changes pressure appreciably; neither released motion nor the
//! norm of accumulated signed displacement is a wear-distance substitute.
//! Tangent axes and material-to-grid registration stay fixed. Supplied local
//! body compliance and its storage/work belong to the external mechanical
//! owner. They are not Mindlin contact springs or cross-cell Cerruti coupling.
//! Frictional dissipation, normal plastic work and recession parameter work
//! remain separate. No thermal partition, heat evolution, or debris law is
//! inferred, and no combined mechanical-energy theorem is claimed.

pub mod temperature;

use std::collections::BTreeMap;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_math::det;
use fs_tribo::{
    DryTangentialState, FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStep,
    TangentialStepInput, WearAdvance,
};

use crate::interface_binding::thermal_friction::{SourceThermalCoulomb, SourceThermalCoulombError};
use crate::normal_patch::{
    FiniteGapPlasticContinuationReceipt, FiniteGapPlasticDrive, FiniteGapPlasticDriveReceipt,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadReceipt, FiniteGapPlasticReceipt,
};

use super::{
    SourcePlasticArchardWear, SourcePlasticWearAdvance, SourcePlasticWearError,
    SourcePlasticWearState, add_compensated, checkpoint,
};

/// Released local motion and explicit diagonal body response in unchanged axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlasticSlidingInput {
    /// Zero-traction A-minus-B displacement increment [m].
    pub released_increment_m: [f64; 2],
    /// Nonnegative isotropic held-force body response [m/N].
    /// Zero prescribes actual motion; positive response can admit sticking.
    pub compliance_m_per_n: f64,
}

/// Original interface, numerical, or physical refusal, without partial history.
#[derive(Clone, Debug, PartialEq)]
pub enum SourcePlasticSlidingError {
    /// Incompatible original sources, unsupported coordinates, or invalid shape.
    Invalid(&'static str),
    /// Original friction binding failed before any cell was advanced.
    Source(SourceThermalCoulombError),
    /// Original source or native Coulomb failure at this exact grid cell.
    Cell {
        index: usize,
        error: SourceThermalCoulombError,
    },
    /// Original wear, normal equilibrium, physical envelope, or cancellation.
    Wear(SourcePlasticWearError),
}
impl From<SourcePlasticWearError> for SourcePlasticSlidingError {
    fn from(error: SourcePlasticWearError) -> Self {
        Self::Wear(error)
    }
}
impl core::fmt::Display for SourcePlasticSlidingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source hardness-limited sliding/wear: {self:?}")
    }
}
impl std::error::Error for SourcePlasticSlidingError {}

/// One immutable original pressure/wear owner and its sourced Coulomb law.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticSlidingWear {
    wear: SourcePlasticArchardWear,
    friction: SourceThermalCoulomb,
    initial_law: FrictionLaw,
    temperature_k: f64,
    identity: ContentHash,
}

/// Normal history, removed material, and actual signed shear travel accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticSlidingWearState {
    owner: ContentHash,
    wear: SourcePlasticWearState,
    tangential: Vec<DryTangentialState>,
}
impl SourcePlasticSlidingWearState {
    /// Complete original normal history, both profiles, and removed volumes.
    pub const fn wear_state(&self) -> &SourcePlasticWearState {
        &self.wear
    }
    /// Native actual slip histories in the original fixed cell/tangent order.
    pub fn tangential_states(&self) -> &[DryTangentialState] {
        &self.tangential
    }
}

/// Native per-cell result without duplicating the candidate's retained state.
#[derive(Clone, Debug, PartialEq)]
pub struct PlasticCellSliding {
    /// Actual supplied released motion and external body response.
    pub input: PlasticSlidingInput,
    /// Beginning normal pressure used for friction and wear [Pa].
    pub pressure_pa: f64,
    /// Index into the returned original source samples. Unloaded cells have
    /// no coefficient query or source claim at zero pressure.
    pub source_index: Option<usize>,
    /// Original native traction, actual slip, branch, dissipation and residuals.
    pub response: TangentialStep<()>,
}

/// Atomic sliding/removal candidate and separate original physical accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticSlidingWearAdvance<B = ()> {
    /// Original complete normal, hardness, wear, friction and registration owner.
    pub source_identity: ContentHash,
    /// All accepted normal, shear, and removed-material histories.
    pub candidate: SourcePlasticSlidingWearState,
    /// One native response for every original cell, including unloaded cells.
    pub cells: Vec<PlasticCellSliding>,
    /// Exact selected claims at the actual loaded pressures, shared when equal.
    pub friction_sources: Vec<ResolvedInterfaceStatePoint>,
    /// Supplied positive held-force interval [s].
    pub duration_s: f64,
    /// Explicit original per-cell numerical and effort budgets.
    pub friction_controls: TangentialSolveConfig,
    /// Compensated sum of the native nonnegative frictional work [J].
    /// No heat allocation is implied.
    pub frictional_work_j: f64,
    /// Original Archard volume increments in immutable A/B order.
    pub wear: [WearAdvance; 2],
    /// Actual represented recession increments for each original body [m].
    pub recession_increment_m: [Vec<f64>; 2],
    /// Beginning-pressure-weighted actual unsigned sliding [N m].
    pub pressure_weighted_sliding_j: f64,
    /// Contact storage, irreversible indentation, normal and parameter work.
    pub normal: FiniteGapPlasticContinuationReceipt,
    /// Unit at held approach, or the native load/drive receipt for that boundary.
    pub boundary: B,
}

/// Ordinary normal loading with the complete shear and removal history retained.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticSlidingNormalAdvance<B = ()> {
    /// Original complete normal, wear and friction owner.
    pub source_identity: ContentHash,
    /// Updated normal equilibrium and unchanged actual shear travel/inventories.
    pub candidate: SourcePlasticSlidingWearState,
    /// Original pressure, elastic storage and irreversible indentation response.
    pub receipt: FiniteGapPlasticReceipt,
    /// Unit at prescribed approach, or the native load/drive receipt.
    pub boundary: B,
}

impl SourcePlasticSlidingWear {
    /// Bind exact original ordered sources at one held interface temperature.
    ///
    /// Friction may additionally own a pressure coordinate absent from the
    /// fixed normal/wear query. All other coordinates and quantities must match
    /// exactly. Speed/path-dependent coefficients need an implicit physical
    /// source closure and are not silently frozen by this adapter.
    pub fn new(
        wear: SourcePlasticArchardWear,
        friction: SourceThermalCoulomb,
    ) -> Result<Self, SourcePlasticSlidingError> {
        let query = friction.initial_source().query();
        Self::admit_query(&wear, &friction, &query)?;
        let temperature_k = query.axes()[&friction.axes().temperature];
        let (_, initial_law) = friction
            .coefficients_at_query(&query)
            .map_err(SourcePlasticSlidingError::Source)?;
        let identity = Self::source_identity(&wear, &friction);
        Ok(Self {
            wear,
            friction,
            initial_law,
            temperature_k,
            identity,
        })
    }

    fn admit_query(
        wear: &SourcePlasticArchardWear,
        friction: &SourceThermalCoulomb,
        query: &QueryPoint,
    ) -> Result<(), SourcePlasticSlidingError> {
        let normal = wear.contact().normal_source();
        if friction.interface() != normal.elastic().interface() || !friction.query_axes().is_empty()
        {
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

    fn source_identity(
        wear: &SourcePlasticArchardWear,
        friction: &SourceThermalCoulomb,
    ) -> ContentHash {
        let mut hash =
            DomainHasher::new("org.frankensim.fs-contact.source-plastic-sliding-wear.v1");
        hash.update(wear.identity().as_bytes());
        hash.update(friction.identity().as_bytes());
        hash.finalize()
    }

    fn query_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<QueryPoint, SourceThermalCoulombError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(SourceThermalCoulombError::Invalid(
                "sliding requires a positive absolute interface temperature",
            ));
        }
        Ok(self
            .friction
            .initial_source()
            .query()
            .with(&self.friction.axes().temperature, temperature_k)?)
    }

    fn coefficients_at_pressure(
        &self,
        pressure_pa: f64,
    ) -> Result<(ResolvedInterfaceStatePoint, FrictionLaw), SourceThermalCoulombError> {
        let mut query = self.query_at_temperature(self.temperature_k)?;
        if let Some(axis) = &self.friction.axes().nominal_pressure {
            query = query.with(axis, pressure_pa)?;
        }
        self.friction.coefficients_at_query(&query)
    }

    fn cell_normal_force(pressure_pa: f64, area_m2: f64) -> Result<f64, SourceThermalCoulombError> {
        let normal_n = pressure_pa * area_m2;
        if !normal_n.is_finite() || (pressure_pa > 0.0 && normal_n == 0.0) {
            return Err(SourceThermalCoulombError::Invalid(
                "a positive cell pressure must produce a representable positive normal force",
            ));
        }
        Ok(normal_n)
    }

    /// Original pressure, indentation-hardness and per-body removal owner.
    pub const fn wear_source(&self) -> &SourcePlasticArchardWear {
        &self.wear
    }
    /// Exact original coefficient selection and pressure/temperature semantics.
    pub const fn friction_source(&self) -> &SourceThermalCoulomb {
        &self.friction
    }
    /// Current held interface temperature [K], including a temperature-only continuation.
    /// Native tangential states separately retain the temperature of their last shear step.
    pub const fn interface_temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Start both normal/wear and signed tangential travel from their original zero states.
    pub fn initial_state(&self) -> SourcePlasticSlidingWearState {
        self.bind_initial_shear(self.wear.initial_state())
    }

    /// Start shear bookkeeping on an already loaded original wear checkpoint.
    ///
    /// Dents, recession and inventories remain intact. The supplied checkpoint
    /// establishes a new zero of signed shear travel; it does not transfer or
    /// relabel a pre-existing tangential history from another owner.
    pub fn state_from_wear(
        &self,
        wear: SourcePlasticWearState,
    ) -> Result<SourcePlasticSlidingWearState, SourcePlasticSlidingError> {
        self.wear.admit(&wear)?;
        Ok(self.bind_initial_shear(wear))
    }

    fn bind_initial_shear(&self, wear: SourcePlasticWearState) -> SourcePlasticSlidingWearState {
        let count = wear.contact_state().contact_state().pressure_pa().len();
        SourcePlasticSlidingWearState {
            owner: self.identity,
            wear,
            tangential: vec![self.initial_law.initial_state(); count],
        }
    }

    /// Change normal approach without resetting signed shear travel or removed material.
    pub fn advance_normal(
        &self,
        before: &SourcePlasticSlidingWearState,
        approach_m: f64,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticSlidingNormalAdvance, SourcePlasticSlidingError> {
        self.advance_normal_with(before, cx, |state| {
            self.wear
                .advance_normal(state, approach_m, cx)
                .map(|next| (next.candidate, next.receipt, ()))
        })
    }

    /// Change applied normal force while retaining all actual shear and wear history.
    pub fn advance_normal_at_load(
        &self,
        before: &SourcePlasticSlidingWearState,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        SourcePlasticSlidingNormalAdvance<FiniteGapPlasticLoadReceipt>,
        SourcePlasticSlidingError,
    > {
        self.advance_normal_with(before, cx, |state| {
            self.wear
                .advance_normal_at_load(state, target_force_n, controls, cx)
                .map(|next| (next.candidate, next.receipt, next.load))
        })
    }

    /// Change the loading-device endpoint with original shear and removal history.
    /// No actuator-work path between device endpoints is inferred.
    pub fn advance_normal_with_drive(
        &self,
        before: &SourcePlasticSlidingWearState,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        SourcePlasticSlidingNormalAdvance<FiniteGapPlasticDriveReceipt>,
        SourcePlasticSlidingError,
    > {
        self.advance_normal_with(before, cx, |state| {
            self.wear
                .advance_normal_with_drive(state, drive, controls, cx)
                .map(|next| (next.candidate, next.receipt, next.drive))
        })
    }

    fn advance_normal_with<B>(
        &self,
        before: &SourcePlasticSlidingWearState,
        cx: &Cx<'_>,
        load: impl FnOnce(
            &SourcePlasticWearState,
        ) -> Result<
            (SourcePlasticWearState, FiniteGapPlasticReceipt, B),
            SourcePlasticWearError,
        >,
    ) -> Result<SourcePlasticSlidingNormalAdvance<B>, SourcePlasticSlidingError> {
        checkpoint(cx)?;
        self.admit(before)?;
        let (wear, receipt, boundary) = load(&before.wear)?;
        let candidate = SourcePlasticSlidingWearState {
            owner: self.identity,
            wear,
            tangential: before.tangential.clone(),
        };
        checkpoint(cx)?;
        Ok(SourcePlasticSlidingNormalAdvance {
            source_identity: self.identity,
            candidate,
            receipt,
            boundary,
        })
    }

    /// Resolve actual sliding, remove material, and re-equilibrate at held approach.
    pub fn advance(
        &self,
        before: &SourcePlasticSlidingWearState,
        motion: &[PlasticSlidingInput],
        duration_s: f64,
        controls: TangentialSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticSlidingWearAdvance, SourcePlasticSlidingError> {
        self.advance_with(before, motion, duration_s, controls, cx, |slip| {
            self.wear
                .advance_wear(&before.wear, slip, cx)
                .map(|wear| (wear, ()))
        })
    }

    /// Resolve sliding and recession while the already equilibrated force stays fixed.
    pub fn advance_at_load(
        &self,
        before: &SourcePlasticSlidingWearState,
        motion: &[PlasticSlidingInput],
        duration_s: f64,
        controls: TangentialSolveConfig,
        target_force_n: f64,
        normal_controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        SourcePlasticSlidingWearAdvance<FiniteGapPlasticLoadReceipt>,
        SourcePlasticSlidingError,
    > {
        self.advance_with(before, motion, duration_s, controls, cx, |slip| {
            self.wear
                .advance_wear_at_load(&before.wear, slip, target_force_n, normal_controls, cx)
                .map(|next| (next.wear, next.load))
        })
    }

    /// Resolve sliding and recession against the same held, equilibrated loading spring.
    pub fn advance_with_drive(
        &self,
        before: &SourcePlasticSlidingWearState,
        motion: &[PlasticSlidingInput],
        duration_s: f64,
        controls: TangentialSolveConfig,
        drive: FiniteGapPlasticDrive,
        normal_controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        SourcePlasticSlidingWearAdvance<FiniteGapPlasticDriveReceipt>,
        SourcePlasticSlidingError,
    > {
        self.advance_with(before, motion, duration_s, controls, cx, |slip| {
            self.wear
                .advance_wear_with_drive(&before.wear, slip, drive, normal_controls, cx)
                .map(|next| (next.wear, next.drive))
        })
    }

    fn advance_with<B>(
        &self,
        before: &SourcePlasticSlidingWearState,
        motion: &[PlasticSlidingInput],
        duration_s: f64,
        controls: TangentialSolveConfig,
        cx: &Cx<'_>,
        remove: impl FnOnce(&[f64]) -> Result<(SourcePlasticWearAdvance, B), SourcePlasticWearError>,
    ) -> Result<SourcePlasticSlidingWearAdvance<B>, SourcePlasticSlidingError> {
        checkpoint(cx)?;
        self.admit(before)?;
        let pressure = before.wear.contact_state().contact_state().pressure_pa();
        if motion.len() != pressure.len() {
            return Err(SourcePlasticSlidingError::Invalid(
                "sliding requires one local motion for every original grid cell",
            ));
        }
        let grid = self.wear.contact().model().config().grid;
        let area = grid.cell_width_m * grid.cell_depth_m;
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
            let normal_n = Self::cell_normal_force(p, area).map_err(error)?;
            // A detached cell needs no coefficient at an unsupported zero
            // pressure. The native unloaded branch has zero effort and work;
            // its receipt deliberately carries no source sample.
            let source_index = if p == 0.0 {
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
                    let (source, law) = self.coefficients_at_pressure(p).map_err(error)?;
                    let slot = sources.len();
                    sources.push(source);
                    laws.push(law);
                    sample_indices.insert(key, slot);
                    slot
                };
                Some(slot)
            };
            let law = source_index.map_or(&self.initial_law, |slot| &laws[slot]);
            let native = law
                .advance(
                    &before.tangential[index],
                    &TangentialStepInput {
                        normal_n,
                        slip_increment_m: input.released_increment_m,
                        dt_s: duration_s,
                        contact_temperature_k: Some(self.temperature_k),
                        nominal_pressure_pa: self
                            .friction
                            .axes()
                            .nominal_pressure
                            .as_ref()
                            .map(|_| p),
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
        let (removed, boundary) = remove(&slip)?;
        checkpoint(cx)?;
        Ok(SourcePlasticSlidingWearAdvance {
            source_identity: self.identity,
            candidate: SourcePlasticSlidingWearState {
                owner: self.identity,
                wear: removed.candidate,
                tangential,
            },
            cells,
            friction_sources: sources,
            duration_s,
            friction_controls: controls,
            frictional_work_j: work,
            wear: removed.wear,
            recession_increment_m: removed.recession_increment_m,
            pressure_weighted_sliding_j: removed.pressure_weighted_sliding_j,
            normal: removed.receipt,
            boundary,
        })
    }

    fn admit(
        &self,
        before: &SourcePlasticSlidingWearState,
    ) -> Result<(), SourcePlasticSlidingError> {
        self.wear.admit(&before.wear)?;
        if before.owner != self.identity
            || before.tangential.len()
                != before
                    .wear
                    .contact_state()
                    .contact_state()
                    .pressure_pa()
                    .len()
        {
            return Err(SourcePlasticSlidingError::Invalid(
                "sliding requires its complete original fixed-grid history",
            ));
        }
        Ok(())
    }
}
