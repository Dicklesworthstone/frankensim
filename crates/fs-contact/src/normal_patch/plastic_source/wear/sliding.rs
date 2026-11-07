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

mod cells;
pub mod hardness;
pub mod temperature;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_tribo::{
    DryTangentialState, FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStep,
    WearAdvance,
};

use crate::interface_binding::thermal_friction::{SourceThermalCoulomb, SourceThermalCoulombError};
use crate::normal_patch::{
    FiniteGapPlasticContinuationReceipt, FiniteGapPlasticDrive, FiniteGapPlasticDriveReceipt,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadReceipt, FiniteGapPlasticReceipt,
};

use super::coefficients::LocalWearCoefficientSamples;
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
    /// Original wear claims at actual local pressure/irreversible mean speed.
    /// Present only when the wear owner evolves local coefficient coordinates.
    pub coefficient_samples: Option<LocalWearCoefficientSamples>,
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
    /// fixed normal query. All other friction coordinates and quantities must
    /// match exactly. Local wear coefficients may follow the actual native
    /// irreversible sliding speed after friction resolves its physical closure.
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
        cells::admit_query(wear.contact().normal_source(), friction, query)
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
        cells::query_at_temperature(&self.friction, temperature_k)
    }

    fn coefficients_at_pressure(
        &self,
        pressure_pa: f64,
    ) -> Result<(ResolvedInterfaceStatePoint, FrictionLaw), SourceThermalCoulombError> {
        cells::coefficients_at_pressure(&self.friction, self.temperature_k, pressure_pa)
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
                .advance_wear_during(&before.wear, slip, duration_s, cx)
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
                .advance_wear_at_load_during(
                    &before.wear,
                    slip,
                    duration_s,
                    target_force_n,
                    normal_controls,
                    cx,
                )
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
                .advance_wear_with_drive_during(
                    &before.wear,
                    slip,
                    duration_s,
                    drive,
                    normal_controls,
                    cx,
                )
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
        let grid = self.wear.contact().model().config().grid;
        let area = grid.cell_width_m * grid.cell_depth_m;
        let sliding = cells::advance(
            &self.friction,
            &self.initial_law,
            self.temperature_k,
            &before.tangential,
            pressure,
            area,
            motion,
            duration_s,
            controls,
            cx,
        )?;
        let (removed, boundary) = remove(&sliding.unsigned_slip_m)?;
        checkpoint(cx)?;
        Ok(SourcePlasticSlidingWearAdvance {
            source_identity: self.identity,
            candidate: SourcePlasticSlidingWearState {
                owner: self.identity,
                wear: removed.candidate,
                tangential: sliding.tangential,
            },
            cells: sliding.cells,
            friction_sources: sliding.friction_sources,
            duration_s,
            friction_controls: controls,
            frictional_work_j: sliding.frictional_work_j,
            wear: removed.wear,
            recession_increment_m: removed.recession_increment_m,
            pressure_weighted_sliding_j: removed.pressure_weighted_sliding_j,
            coefficient_samples: removed.coefficient_samples,
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
