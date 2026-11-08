//! Retained thin-layer shear on the original plastically flattened contact grid.
//!
//! Normal loading and held-normal shear are separate physical substeps. The
//! original half-space model owns all normal indentation; independent layer
//! potentials own shear only. This is not a pressure-dependent bulk yield law,
//! normal/shear cross-compliance, a moving registration, or a thermal/wear model.
//! Every layer occupies one complete original cell, with extensive coefficients
//! supplied for that cell's area. Unmapped cells have no layer shear in this model.
//!
//! All mapped cells must be compressive during a shear interval. Normal opening
//! with recoverable shear storage refuses: no spring energy, Maxwell memory, or
//! permanent slip is silently discarded or relabelled as heat.

use std::sync::Arc;

use fs_exec::Cx;
use fs_tribo::layer_contact::{
    ThirdBodyLayerContact, ThirdBodyLayerDissipation, ThirdBodyLayerState,
};
use fs_tribo::tangential::multiple::{
    MultipleTangentialError, MultipleTangentialStep, PreparedTangentialContacts,
    TangentialContactInput,
};
use fs_tribo::{
    TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput, TriboError,
};

use super::{
    FiniteGapPlasticAdvance, FiniteGapPlasticContinuation, FiniteGapPlasticContinuationReceipt,
    FiniteGapPlasticDrive, FiniteGapPlasticDriveError, FiniteGapPlasticDriveReceipt,
    FiniteGapPlasticError, FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError,
    FiniteGapPlasticLoadReceipt, FiniteGapPlasticModel, FiniteGapPlasticReceipt,
    FiniteGapPlasticState, checkpoint,
};

/// One fixed cell registration and its original extensive shear law.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLayerSite {
    /// Row-major cell in the original finite-gap grid.
    pub cell_index: usize,
    /// Integrated cell stiffness, yield force and resistance, not per-area moduli.
    /// Scalar and full tensor laws retain their original Maxwell/Kelvin branches.
    pub law: ThirdBodyLayerContact,
}

#[derive(Clone, Debug, PartialEq)]
struct LayerOwner {
    normal: FiniteGapPlasticModel,
    sites: Vec<FiniteGapPlasticLayerSite>,
    response: PreparedTangentialContacts,
}

/// Fixed normal model, cell registration, layer laws and full mechanical response.
///
/// The prepared response contains all cross-cell and two-axis shear recoil.
/// Normal/body shear compliance is not added to the layer's own stored potential.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLayerContact {
    owner: Arc<LayerOwner>,
}

/// Joint original normal and irreversible/recoverable layer history.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLayerState {
    owner: Arc<LayerOwner>,
    normal: Arc<FiniteGapPlasticState>,
    shear: Vec<ThirdBodyLayerState>,
}

impl FiniteGapPlasticLayerState {
    /// Original pressure, permanent indentation, recession and plastic work.
    pub fn normal_state(&self) -> &FiniteGapPlasticState {
        &self.normal
    }

    /// Retained native layer histories in the model's fixed site order.
    pub fn layer_states(&self) -> &[ThirdBodyLayerState] {
        &self.shear
    }
}

/// Normal endpoint update with every original shear history retained.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLayerNormalAdvance<B = ()> {
    /// Adopt only after the enclosing mechanical operation accepts.
    pub candidate: FiniteGapPlasticLayerState,
    /// Original normal storage and plastic work; no normal path work is inferred.
    pub normal: FiniteGapPlasticReceipt,
    /// Unit at prescribed approach, or the original load/spring receipt.
    pub boundary: B,
}

/// Explicit material-coefficient continuation at fixed retained layer strains.
///
/// Only the normal owner re-equilibrates during this substep. Layer spring
/// forces can change while elastic/plastic/Maxwell strains remain fixed. The
/// next positive-time shear interval resolves actual mechanical recoil.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLayerRefreshAdvance<B = ()> {
    /// All original histories bound to the successor coefficients together.
    pub candidate: FiniteGapPlasticLayerState,
    /// Native normal re-equilibration, plastic work, boundary work and parameter work.
    pub normal: FiniteGapPlasticContinuationReceipt,
    /// Native signed change in each layer's mechanical free energy at fixed
    /// retained strains [J], in the original site order. This is not heat,
    /// shear dissipation, or a body's mechanical work at fixed external motion.
    pub layer_coefficient_free_energy_change_j: Vec<f64>,
    /// Compensated sum of the preceding layer potential changes [J].
    /// Caloric/entropy accounting, if requested, belongs to the thermal owner.
    pub total_layer_coefficient_free_energy_change_j: f64,
    /// Unit at held approach, or the original maintained-load/spring receipt.
    pub boundary: B,
}

/// Held-normal simultaneous shear candidate and independent physical accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLayerShearAdvance {
    /// Normal history is unchanged; all layer histories advance together.
    pub candidate: FiniteGapPlasticLayerState,
    /// Native full-response result, without duplicating candidate histories.
    /// Its work uses actual resolved shear motion, including cross-cell recoil.
    pub shear: MultipleTangentialStep<()>,
    /// Original local pressures used by the native layers, in fixed site order [Pa].
    pub pressure_pa: Vec<f64>,
    /// Independently evaluated original plastic, Maxwell and Kelvin losses [J].
    pub dissipation: Vec<ThirdBodyLayerDissipation>,
    /// Sum of cell shear efforts entering the interface in its fixed axes [N].
    pub shear_force_n: [f64; 2],
    /// Sum of x*F_y-y*F_x about the original grid centre [N m].
    /// Body reaction moments have the opposite sign.
    pub spin_moment_n_m: f64,
}

/// Refusal without changing any accepted normal or layer history.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteGapPlasticLayerError {
    /// Invalid registration, shape, arithmetic, or a foreign composite owner.
    Invalid(&'static str),
    /// Original normal constitutive failure.
    Normal(FiniteGapPlasticError),
    /// Original force-controlled normal failure.
    Load(FiniteGapPlasticLoadError),
    /// Original finite-spring normal failure.
    Drive(FiniteGapPlasticDriveError),
    /// Original full-matrix simultaneous shear failure.
    Shear(MultipleTangentialError),
    /// Original layer dissipation evaluation at this original grid cell failed.
    Cell {
        cell_index: usize,
        error: TriboError,
    },
    /// A held-normal shear step requires every registered layer to be compressed.
    DetachedLayer { cell_index: usize },
    /// Opening requires a separate physical release/transfer of this stored energy.
    /// A coefficient refresh reports the larger predecessor/successor storage.
    OpeningWithStoredShear {
        cell_index: usize,
        stored_energy_j: f64,
    },
}

impl From<FiniteGapPlasticError> for FiniteGapPlasticLayerError {
    fn from(value: FiniteGapPlasticError) -> Self {
        Self::Normal(value)
    }
}
impl From<FiniteGapPlasticLoadError> for FiniteGapPlasticLayerError {
    fn from(value: FiniteGapPlasticLoadError) -> Self {
        Self::Load(value)
    }
}
impl From<FiniteGapPlasticDriveError> for FiniteGapPlasticLayerError {
    fn from(value: FiniteGapPlasticDriveError) -> Self {
        Self::Drive(value)
    }
}
impl From<MultipleTangentialError> for FiniteGapPlasticLayerError {
    fn from(value: MultipleTangentialError) -> Self {
        Self::Shear(value)
    }
}
impl core::fmt::Display for FiniteGapPlasticLayerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "retained shear on plastic finite-gap contact: {self:?}")
    }
}
impl std::error::Error for FiniteGapPlasticLayerError {}

impl FiniteGapPlasticLayerContact {
    /// Bind layer sites in strictly increasing original cell order.
    ///
    /// The prepared matrix must use exactly that same site/axis order and the
    /// original mechanical body's compliance, excluding these layer springs.
    /// Its existing contact, storage, evaluation and cancellation limits apply.
    pub fn new(
        normal: FiniteGapPlasticModel,
        sites: Vec<FiniteGapPlasticLayerSite>,
        response: PreparedTangentialContacts,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteGapPlasticLayerError> {
        checkpoint(cx)?;
        if sites.len() != response.contact_count() {
            return Err(FiniteGapPlasticLayerError::Invalid(
                "layer/matrix site count",
            ));
        }
        let count = normal.config().undeformed_gap_m.len();
        let mut previous = None;
        for site in &sites {
            checkpoint(cx)?;
            if site.cell_index >= count || previous.is_some_and(|index| index >= site.cell_index) {
                return Err(FiniteGapPlasticLayerError::Invalid(
                    "unique increasing original layer cell indices",
                ));
            }
            previous = Some(site.cell_index);
        }
        Ok(Self {
            owner: Arc::new(LayerOwner {
                normal,
                sites,
                response,
            }),
        })
    }

    /// Immutable original normal owner and prepared half-space kernel.
    pub fn normal_model(&self) -> &FiniteGapPlasticModel {
        &self.owner.normal
    }

    /// Fixed registration and native law at each simultaneous shear site.
    pub fn sites(&self) -> &[FiniteGapPlasticLayerSite] {
        &self.owner.sites
    }

    /// Complete physical response used by every held-normal shear solve.
    pub fn mechanical_response(&self) -> &PreparedTangentialContacts {
        &self.owner.response
    }

    /// Virgin normal and shear history. Load the registered cells before shearing.
    pub fn initial_state(&self) -> FiniteGapPlasticLayerState {
        self.bind_initial_shear(self.owner.normal.initial_state())
    }

    /// Begin an independently unstrained layer on an existing normal history.
    ///
    /// Accepted dents, recession and cumulative normal plastic work survive.
    /// This does not transfer or reset an existing composite shear checkpoint.
    pub fn state_from_normal(
        &self,
        normal: FiniteGapPlasticState,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLayerState, FiniteGapPlasticLayerError> {
        checkpoint(cx)?;
        self.owner.normal.check_state_owner(&normal)?;
        Ok(self.bind_initial_shear(normal))
    }

    fn bind_initial_shear(&self, normal: FiniteGapPlasticState) -> FiniteGapPlasticLayerState {
        FiniteGapPlasticLayerState {
            owner: Arc::clone(&self.owner),
            normal: Arc::new(normal),
            shear: self
                .owner
                .sites
                .iter()
                .map(|site| site.law.initial_state())
                .collect(),
        }
    }

    /// Change approach without advancing shear time, loss, or plastic slip.
    pub fn advance_normal(
        &self,
        before: &FiniteGapPlasticLayerState,
        approach_m: f64,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLayerNormalAdvance, FiniteGapPlasticLayerError> {
        self.advance_normal_with(before, cx, |state| {
            Ok((self.owner.normal.advance(state, approach_m, cx)?, ()))
        })
    }

    /// Change normal force through the original bounded equilibrium owner.
    pub fn advance_normal_at_load(
        &self,
        before: &FiniteGapPlasticLayerState,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteGapPlasticLayerNormalAdvance<FiniteGapPlasticLoadReceipt>,
        FiniteGapPlasticLayerError,
    > {
        self.advance_normal_with(before, cx, |state| {
            let next = self
                .owner
                .normal
                .advance_at_load(state, target_force_n, controls, cx)?;
            Ok((next.contact, next.load))
        })
    }

    /// Change the original finite loading-spring endpoint with shear retained.
    pub fn advance_normal_with_drive(
        &self,
        before: &FiniteGapPlasticLayerState,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteGapPlasticLayerNormalAdvance<FiniteGapPlasticDriveReceipt>,
        FiniteGapPlasticLayerError,
    > {
        self.advance_normal_with(before, cx, |state| {
            let next = self
                .owner
                .normal
                .advance_with_drive(state, drive, controls, cx)?;
            Ok((next.contact, next.drive))
        })
    }

    fn advance_normal_with<B>(
        &self,
        before: &FiniteGapPlasticLayerState,
        cx: &Cx<'_>,
        load: impl FnOnce(
            &FiniteGapPlasticState,
        ) -> Result<(FiniteGapPlasticAdvance, B), FiniteGapPlasticLayerError>,
    ) -> Result<FiniteGapPlasticLayerNormalAdvance<B>, FiniteGapPlasticLayerError> {
        checkpoint(cx)?;
        self.admit(before)?;
        let (next, boundary) = load(&before.normal)?;
        for (site, state) in self.owner.sites.iter().zip(&before.shear) {
            checkpoint(cx)?;
            if next.candidate.pressure_pa()[site.cell_index] == 0.0 {
                let energy = site.law.stored_energy_j(state);
                if !energy.is_finite() || energy < 0.0 {
                    return Err(FiniteGapPlasticLayerError::Invalid("native layer storage"));
                }
                if energy > 0.0 {
                    return Err(FiniteGapPlasticLayerError::OpeningWithStoredShear {
                        cell_index: site.cell_index,
                        stored_energy_j: energy,
                    });
                }
            }
        }
        checkpoint(cx)?;
        Ok(FiniteGapPlasticLayerNormalAdvance {
            candidate: FiniteGapPlasticLayerState {
                owner: Arc::clone(&self.owner),
                normal: Arc::new(next.candidate),
                shear: before.shear.clone(),
            },
            normal: next.receipt,
            boundary,
        })
    }

    /// Advance every original layer against the complete body response.
    ///
    /// Pressure and permanent normal indentation are held for this positive
    /// interval. Native layer temperatures and geometry are fixed material
    /// inputs. Each trial borrows the same original shear state; full-matrix
    /// sweeps never advance time or dissipate work more than once.
    pub fn advance_shear(
        &self,
        before: &FiniteGapPlasticLayerState,
        released_increment_m: &[[f64; 2]],
        duration_s: f64,
        controls: TangentialSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLayerShearAdvance, FiniteGapPlasticLayerError> {
        checkpoint(cx)?;
        self.admit(before)?;
        if released_increment_m.len() != self.owner.sites.len() {
            return Err(FiniteGapPlasticLayerError::Invalid(
                "one shear increment per layer site",
            ));
        }
        let grid = self.owner.normal.config().grid;
        let area = grid.cell_width_m * grid.cell_depth_m;
        let mut pressure_pa = Vec::with_capacity(self.owner.sites.len());
        let mut inputs = Vec::with_capacity(self.owner.sites.len());
        for (index, site) in self.owner.sites.iter().enumerate() {
            checkpoint(cx)?;
            let pressure = before.normal.pressure_pa()[site.cell_index];
            if pressure == 0.0 {
                return Err(FiniteGapPlasticLayerError::DetachedLayer {
                    cell_index: site.cell_index,
                });
            }
            let normal_n = pressure * area;
            if !normal_n.is_finite() || normal_n <= 0.0 {
                return Err(FiniteGapPlasticLayerError::Invalid(
                    "finite positive cell normal force",
                ));
            }
            pressure_pa.push(pressure);
            inputs.push(TangentialContactInput {
                law: &site.law,
                state: &before.shear[index],
                input: TangentialStepInput {
                    normal_n,
                    slip_increment_m: released_increment_m[index],
                    dt_s: duration_s,
                    contact_temperature_k: Some(site.law.temperature_kelvin()),
                    nominal_pressure_pa: Some(pressure),
                    compliance_m_per_n: 0.0,
                    solve: controls,
                },
            });
        }
        let result = self
            .owner
            .response
            .advance(&inputs, || cx.checkpoint().is_err())?;
        let mut shear = Vec::with_capacity(self.owner.sites.len());
        let mut contacts = Vec::with_capacity(self.owner.sites.len());
        let mut dissipation = Vec::with_capacity(self.owner.sites.len());
        let mut wrench = [0.0; 3];
        let mut correction = [0.0; 3];
        for (index, step) in result.contacts.into_iter().enumerate() {
            checkpoint(cx)?;
            let site = &self.owner.sites[index];
            dissipation.push(
                site.law
                    .step_dissipation(&before.shear[index], &step.state, duration_s)
                    .map_err(|error| FiniteGapPlasticLayerError::Cell {
                        cell_index: site.cell_index,
                        error,
                    })?,
            );
            let (x, y) = grid.center(site.cell_index);
            let terms = [
                step.traction_n[0],
                step.traction_n[1],
                x * step.traction_n[1] - y * step.traction_n[0],
            ];
            for axis in 0..3 {
                let value = terms[axis] - correction[axis];
                let next = wrench[axis] + value;
                correction[axis] = (next - wrench[axis]) - value;
                wrench[axis] = next;
            }
            shear.push(step.state);
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
        if wrench.iter().any(|value| !value.is_finite()) {
            return Err(FiniteGapPlasticLayerError::Invalid(
                "finite layer resultant wrench",
            ));
        }
        checkpoint(cx)?;
        Ok(FiniteGapPlasticLayerShearAdvance {
            candidate: FiniteGapPlasticLayerState {
                owner: Arc::clone(&self.owner),
                normal: Arc::clone(&before.normal),
                shear,
            },
            shear: MultipleTangentialStep {
                contacts,
                contact_work_residuals_j: result.contact_work_residuals_j,
                sweeps: result.sweeps,
                law_evaluations: result.law_evaluations,
                maximum_compatibility_residual_m: result.maximum_compatibility_residual_m,
                stored_change_j: result.stored_change_j,
                dissipated_j: result.dissipated_j,
                body_work_residual_j: result.body_work_residual_j,
            },
            pressure_pa,
            dissipation,
            shear_force_n: [wrench[0], wrench[1]],
            spin_moment_n_m: wrench[2],
        })
    }

    /// Continue explicit material coefficients at the accepted normal approach.
    ///
    /// Native normal continuation retains dents and recession while changed
    /// modulus/hardness may create additional normal plastic flow. Native
    /// layer continuation keeps every elastic, plastic and Maxwell strain
    /// fixed and reports the signed coefficient potential change. Changing
    /// viscosity or yield force alone contributes no stored-energy change.
    ///
    /// This is not shear equilibrium at held external displacement: the
    /// prepared response maps the next released-motion/mean-force interval.
    /// No shear time, slip increment, loss, heat or thermal field is inferred.
    /// Original grid, cell registration, layer domains/branch order and full
    /// mechanical response stay fixed. Prescribed temperature is part of the
    /// explicitly supplied successor laws, not an autonomous temperature solve.
    pub fn refresh_from(
        &self,
        previous: &Self,
        before: &FiniteGapPlasticLayerState,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLayerRefreshAdvance, FiniteGapPlasticLayerError> {
        self.refresh_with(previous, before, cx, |state| {
            Ok((
                self.owner
                    .normal
                    .refresh_from(&previous.owner.normal, state, cx)?,
                (),
            ))
        })
    }

    /// Continue coefficients while maintaining the already accepted normal force.
    ///
    /// The native force owner acts directly on accepted normal history, without
    /// a held-approach plastic intermediate. Layer strains remain fixed under
    /// the same coefficient-substep convention as the held-approach refresh.
    pub fn refresh_from_at_load(
        &self,
        previous: &Self,
        before: &FiniteGapPlasticLayerState,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteGapPlasticLayerRefreshAdvance<FiniteGapPlasticLoadReceipt>,
        FiniteGapPlasticLayerError,
    > {
        self.refresh_with(previous, before, cx, |state| {
            let next = self.owner.normal.refresh_from_at_load(
                &previous.owner.normal,
                state,
                target_force_n,
                controls,
                cx,
            )?;
            Ok((next.continuation, next.load))
        })
    }

    /// Continue coefficients with the original loading spring held at its position.
    ///
    /// Native normal work is energy released by that spring; the stationary
    /// actuator does no work. It stays separate from the fixed-strain layer
    /// coefficient potential change and from any subsequent shear loss.
    pub fn refresh_from_with_drive(
        &self,
        previous: &Self,
        before: &FiniteGapPlasticLayerState,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteGapPlasticLayerRefreshAdvance<FiniteGapPlasticDriveReceipt>,
        FiniteGapPlasticLayerError,
    > {
        self.refresh_with(previous, before, cx, |state| {
            let next = self.owner.normal.refresh_from_with_drive(
                &previous.owner.normal,
                state,
                drive,
                controls,
                cx,
            )?;
            Ok((next.continuation, next.drive))
        })
    }

    fn refresh_with<B>(
        &self,
        previous: &Self,
        before: &FiniteGapPlasticLayerState,
        cx: &Cx<'_>,
        refresh_normal: impl FnOnce(
            &FiniteGapPlasticState,
        ) -> Result<
            (FiniteGapPlasticContinuation, B),
            FiniteGapPlasticLayerError,
        >,
    ) -> Result<FiniteGapPlasticLayerRefreshAdvance<B>, FiniteGapPlasticLayerError> {
        checkpoint(cx)?;
        previous.admit(before)?;
        if self.owner.sites.len() != previous.owner.sites.len()
            || self
                .owner
                .sites
                .iter()
                .zip(&previous.owner.sites)
                .any(|(new, old)| new.cell_index != old.cell_index)
            || (!Arc::ptr_eq(&self.owner, &previous.owner)
                && self.owner.response != previous.owner.response)
        {
            return Err(FiniteGapPlasticLayerError::Invalid(
                "coefficient refresh requires original layer registration and mechanical response",
            ));
        }

        // Stage native layer transfers before the expensive normal equilibrium.
        // Each retains its original branch correspondence and finite strain domain.
        let mut shear = Vec::with_capacity(self.owner.sites.len());
        let mut layer_work = Vec::with_capacity(self.owner.sites.len());
        let mut total_work = 0.0;
        let mut correction = 0.0;
        for (site, old_state) in self.owner.sites.iter().zip(&before.shear) {
            checkpoint(cx)?;
            let next = site
                .law
                .continue_state_with_work(old_state)
                .map_err(|error| FiniteGapPlasticLayerError::Cell {
                    cell_index: site.cell_index,
                    error,
                })?;
            let value = next.coefficient_free_energy_change_j - correction;
            let sum = total_work + value;
            correction = (sum - total_work) - value;
            total_work = sum;
            layer_work.push(next.coefficient_free_energy_change_j);
            shear.push(next.state);
        }
        if !total_work.is_finite() {
            return Err(FiniteGapPlasticLayerError::Invalid(
                "finite layer coefficient work",
            ));
        }
        let (normal, boundary) = refresh_normal(&before.normal)?;
        for (index, site) in self.owner.sites.iter().enumerate() {
            checkpoint(cx)?;
            if normal.candidate.pressure_pa()[site.cell_index] == 0.0 {
                let old_energy = previous.owner.sites[index]
                    .law
                    .stored_energy_j(&before.shear[index]);
                let new_energy = site.law.stored_energy_j(&shear[index]);
                if !old_energy.is_finite()
                    || !new_energy.is_finite()
                    || old_energy < 0.0
                    || new_energy < 0.0
                {
                    return Err(FiniteGapPlasticLayerError::Invalid("native layer storage"));
                }
                let energy = old_energy.max(new_energy);
                if energy > 0.0 {
                    return Err(FiniteGapPlasticLayerError::OpeningWithStoredShear {
                        cell_index: site.cell_index,
                        stored_energy_j: energy,
                    });
                }
            }
        }
        checkpoint(cx)?;
        let normal_state = if normal.candidate.identity() == before.normal.identity() {
            Arc::clone(&before.normal)
        } else {
            Arc::new(normal.candidate)
        };
        Ok(FiniteGapPlasticLayerRefreshAdvance {
            candidate: FiniteGapPlasticLayerState {
                owner: Arc::clone(&self.owner),
                normal: normal_state,
                shear,
            },
            normal: normal.receipt,
            layer_coefficient_free_energy_change_j: layer_work,
            total_layer_coefficient_free_energy_change_j: total_work,
            boundary,
        })
    }

    fn admit(&self, before: &FiniteGapPlasticLayerState) -> Result<(), FiniteGapPlasticLayerError> {
        if (!Arc::ptr_eq(&self.owner, &before.owner) && self.owner != before.owner)
            || before.shear.len() != self.owner.sites.len()
        {
            return Err(FiniteGapPlasticLayerError::Invalid(
                "original layer/normal/matrix owner",
            ));
        }
        self.owner.normal.check_state_owner(&before.normal)?;
        Ok(())
    }
}
