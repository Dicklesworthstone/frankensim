//! Local Archard recession of a retained hardness-limited contact patch.
//!
//! Ordinary updates integrate the *initial* equilibrated cell pressure against the
//! supplied unsigned local sliding distance, holding each sampled local k/H
//! fixed over the step. Original k claims may follow actual cell pressure and
//! irreversible mean sliding speed; separate body hardness stays independent.
//! It removes volume from each original body and re-equilibrates the increased
//! relative gap at held approach. Refine sliding increments when wear changes
//! pressure appreciably. This bounded surface law does not transport debris or
//! infer local sliding from macroscopic displacement, and it assigns no heat.
//! The implicit adapter instead closes the same quadrature at endpoint pressure,
//! retaining the original normal, inventory and parameter-work owners.

/// Original local pressure/speed-dependent wear coefficients and query receipts.
pub mod coefficients;
/// Wear under an explicitly finite loading-device stiffness.
pub mod drive;
/// Implicit endpoint-pressure integration of the original local wear law.
pub mod implicit;
/// Adaptive integration of pressure-dependent local surface recession.
pub mod integration;
/// Prescribed-load wear without intermediate indentation histories.
pub mod load;
/// Actual source Coulomb sliding drives local recession and pressure feedback.
pub mod sliding;
/// Prescribed-temperature continuation with both original wear inventories.
pub mod temperature;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::{MaterialCard, PropertyKey};
use fs_material::hardness::IndentationHardnessStatePoint;
use fs_material::state_point::{MaterialStatePointError, ResolvedInterfaceStatePoint};
use fs_qty::{Dims, QuantitySpec};
use fs_tribo::{ArchardLaw, LocalArchardSample, TriboError, WearAdvance, WearState};

use coefficients::{LocalWearAxes, LocalWearCoefficientSamples};

use super::drive::SourcePlasticDriveError;
use super::load::SourcePlasticLoadError;
use super::temperature::SourcePlasticTemperatureError;
use super::{SourcePlasticContact, SourcePlasticContactError, SourcePlasticContactState};
use crate::normal_patch::{
    FiniteGapContactError, FiniteGapPlasticContinuation, FiniteGapPlasticContinuationReceipt,
    FiniteGapPlasticError, FiniteGapPlasticReceipt,
};

/// Refusal leaves contact history, removed volumes, and recession unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum SourcePlasticWearError {
    /// Incompatible source ownership, unsupported coordinates, or invalid slip.
    Invalid(&'static str),
    /// Original normal-contact applicability or numerical refusal.
    Contact(SourcePlasticContactError),
    /// The prescribed normal load or original source load envelope refused.
    Load(SourcePlasticLoadError),
    /// Finite-stiffness device equilibrium or original physical bounds refused.
    Drive(SourcePlasticDriveError),
    /// Native Archard arithmetic or provenance refusal.
    Tribology(TriboError),
    /// Original normal and pressure-cap sources refused temperature continuation.
    Temperature(SourcePlasticTemperatureError),
    /// The original selected wear coefficient or per-body hardness refused refresh.
    Material(MaterialStatePointError),
}
impl From<SourcePlasticContactError> for SourcePlasticWearError {
    fn from(error: SourcePlasticContactError) -> Self {
        Self::Contact(error)
    }
}
impl From<SourcePlasticLoadError> for SourcePlasticWearError {
    fn from(error: SourcePlasticLoadError) -> Self {
        Self::Load(error)
    }
}
impl From<SourcePlasticDriveError> for SourcePlasticWearError {
    fn from(error: SourcePlasticDriveError) -> Self {
        Self::Drive(error)
    }
}
impl From<TriboError> for SourcePlasticWearError {
    fn from(error: TriboError) -> Self {
        Self::Tribology(error)
    }
}
impl From<SourcePlasticTemperatureError> for SourcePlasticWearError {
    fn from(error: SourcePlasticTemperatureError) -> Self {
        Self::Temperature(error)
    }
}
impl From<MaterialStatePointError> for SourcePlasticWearError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for SourcePlasticWearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source finite-gap Archard wear: {self:?}")
    }
}
impl std::error::Error for SourcePlasticWearError {}

/// Exact original interface coefficients and separate body hardness.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticArchardWear {
    contact: SourcePlasticContact,
    coefficients: ResolvedInterfaceStatePoint,
    coefficient_keys: [PropertyKey; 2],
    local_axes: Option<LocalWearAxes>,
    hardness: [IndentationHardnessStatePoint; 2],
    laws: [ArchardLaw; 2],
    maximum_recession_m: [f64; 2],
    identity: ContentHash,
}

/// Contact and removed material are accepted as one immutable candidate.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearState {
    owner: ContentHash,
    contact: SourcePlasticContactState,
    volumes: [WearState; 2],
    recession_m: [Vec<f64>; 2],
}
impl SourcePlasticWearState {
    /// Retained source-bound normal history.
    pub const fn contact_state(&self) -> &SourcePlasticContactState {
        &self.contact
    }
    /// Cumulative removed volume in original body A/B order [m^3].
    pub fn volumes_m3(&self) -> [f64; 2] {
        self.volumes.map(|state| state.volume_m3())
    }
    /// Cumulative recession at each original grid cell in body A/B order [m].
    pub fn recession_m(&self) -> [&[f64]; 2] {
        [&self.recession_m[0], &self.recession_m[1]]
    }
}

/// A normal loading step that retains both wear inventories and surface profiles.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearNormalAdvance {
    /// Complete candidate normal and wear state.
    pub candidate: SourcePlasticWearState,
    /// Original normal pressure, elastic storage, and plastic work.
    pub receipt: FiniteGapPlasticReceipt,
}

/// Local removal and its held-approach normal re-equilibration.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearAdvance {
    /// Complete owner of coefficient, hardness, and geometry inputs.
    pub source_identity: ContentHash,
    /// Normal history, both wear inventories, and both recession profiles.
    pub candidate: SourcePlasticWearState,
    /// Native Archard volumes and retained interface provenance in A/B order.
    pub wear: [WearAdvance; 2],
    /// Recession added to each original body at every original grid cell [m].
    pub recession_increment_m: [Vec<f64>; 2],
    /// Pressure quadrature `sum(p * unsigned_local_slip * area)` [N m]. Ordinary
    /// advances use beginning pressure; the implicit adapter declares its
    /// endpoint-pressure approximation and residuals in its enclosing receipt.
    /// This is a sliding measure and has no automatic heat interpretation.
    pub pressure_weighted_sliding_j: f64,
    /// Actual local wear-claim queries, present for a pressure/speed-dependent
    /// source. Static-coefficient steps reuse their original source directly.
    pub coefficient_samples: Option<LocalWearCoefficientSamples>,
    /// Native normal response and old/new storage plus signed parameter work.
    /// Material removal at held approach may release stored elastic energy;
    /// the caller must retain this account in its coupled energy balance.
    pub receipt: FiniteGapPlasticContinuationReceipt,
}

/// Prescribed temperatures with both original wear inventories retained.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearThermalRefresh {
    source: SourcePlasticArchardWear,
    candidate: SourcePlasticWearState,
    receipt: FiniteGapPlasticContinuationReceipt,
}
impl SourcePlasticWearThermalRefresh {
    /// Refreshed original normal, coefficient, and separate body-hardness sources.
    pub const fn source(&self) -> &SourcePlasticArchardWear {
        &self.source
    }
    /// Retained per-body removal and irreversible contact history.
    pub const fn candidate(&self) -> &SourcePlasticWearState {
        &self.candidate
    }
    /// Native normal endpoint and signed parameter-work account.
    pub const fn receipt(&self) -> &FiniteGapPlasticContinuationReceipt {
        &self.receipt
    }
    /// Adopt the complete continuation without copying its cell arrays.
    pub fn into_parts(
        self,
    ) -> (
        SourcePlasticArchardWear,
        SourcePlasticWearState,
        FiniteGapPlasticContinuationReceipt,
    ) {
        (self.source, self.candidate, self.receipt)
    }
}

impl SourcePlasticArchardWear {
    /// Bind original k_A/k_B and typed indentation hardness at fixed conditions.
    ///
    /// The pressure cap's minimum hardness is not substituted for the separate
    /// per-body wear hardness. Exact selected test keys remain part of identity.
    /// Each recession limit is a caller-declared physical applicability bound;
    /// zero admits only a surface with no wear during the attempted update.
    pub fn new(
        contact: SourcePlasticContact,
        coefficients: ResolvedInterfaceStatePoint,
        coefficient_keys: [PropertyKey; 2],
        hardness: [IndentationHardnessStatePoint; 2],
        maximum_recession_m: [f64; 2],
    ) -> Result<Self, SourcePlasticWearError> {
        Self::build(
            contact,
            coefficients,
            coefficient_keys,
            hardness,
            maximum_recession_m,
            None,
        )
    }

    /// Bind original wear coefficients to actual local pressure and/or speed.
    ///
    /// Additional typed coordinates belong to this contact/sliding owner. The
    /// same exact claims and original fixed conditions are retained on every
    /// query; unsupported actual conditions refuse the whole candidate. Local
    /// speed requires an actual duration and unsigned irreversible sliding, not
    /// released body motion. No coefficient or hardness calibration is inferred.
    pub fn new_with_local_coefficients(
        contact: SourcePlasticContact,
        coefficients: ResolvedInterfaceStatePoint,
        coefficient_keys: [PropertyKey; 2],
        hardness: [IndentationHardnessStatePoint; 2],
        maximum_recession_m: [f64; 2],
        axes: LocalWearAxes,
    ) -> Result<Self, SourcePlasticWearError> {
        Self::build(
            contact,
            coefficients,
            coefficient_keys,
            hardness,
            maximum_recession_m,
            Some(axes),
        )
    }

    fn build(
        contact: SourcePlasticContact,
        coefficients: ResolvedInterfaceStatePoint,
        coefficient_keys: [PropertyKey; 2],
        hardness: [IndentationHardnessStatePoint; 2],
        maximum_recession_m: [f64; 2],
        local_axes: Option<LocalWearAxes>,
    ) -> Result<Self, SourcePlasticWearError> {
        let original = contact.normal_source().interface_state();
        if coefficients.card_identity() != original.card_identity()
            || coefficients.surface_materials() != original.surface_materials()
        {
            return Err(invalid(
                "wear requires the original ordered interface and fixed query",
            ));
        }
        let owned = if let Some(axes) = &local_axes {
            axes.admit(&coefficients, original)?;
            axes.owned()
        } else {
            if coefficients.query_point() != original.query_point()
                || coefficients.axis_quantities() != original.axis_quantities()
            {
                return Err(invalid("wear requires the original fixed interface query"));
            }
            Vec::new()
        };
        let elastic = contact.normal_source().elastic();
        let bulk = [
            elastic.surface_a_state().resolved(),
            elastic.surface_b_state().resolved(),
        ];
        let mut laws = [ArchardLaw {
            wear_coefficient: 0.0,
        }; 2];
        let mut selected = [0_usize; 2];
        for side in 0..2 {
            let k = coefficients
                .property_by_key(&coefficient_keys[side])
                .ok_or_else(|| invalid("missing exact original wear coefficient"))?;
            let h = hardness[side].resolved();
            if coefficient_keys[side].quantity() != QuantitySpec::dimensional(Dims::NONE)
                || !k.value_si().is_finite()
                || k.value_si() < 0.0
                || !k.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
                || !hardness[side]
                    .property()
                    .is_static_at_fixed_external_conditions()
                || h.card_identity() != bulk[side].card_identity()
                || h.material() != bulk[side].material()
                || h.query_point() != bulk[side].query_point()
                || h.axis_quantities() != bulk[side].axis_quantities()
            {
                return Err(invalid(
                    "wear requires original k with only owned coordinates and ordered typed indentation hardness",
                ));
            }
            if !maximum_recession_m[side].is_finite() || maximum_recession_m[side] < 0.0 {
                return Err(invalid("finite nonnegative per-body recession bound"));
            }
            selected[side] = coefficients
                .properties()
                .iter()
                .position(|property| property == k)
                .ok_or_else(|| {
                    invalid("wear coefficient is not retained by its original bundle")
                })?;
            laws[side].wear_coefficient = k.value_si();
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.source-plastic-archard.v1");
        hash.update(contact.identity().as_bytes());
        hash.update(coefficients.identity().as_bytes());
        for side in 0..2 {
            hash.update(&(selected[side] as u64).to_le_bytes());
            hash.update(hardness[side].identity().as_bytes());
            hash.update(&maximum_recession_m[side].to_bits().to_le_bytes());
        }
        if let Some(axes) = &local_axes {
            hash.update(b"local-pressure-speed-coefficients-v1");
            for axis in [&axes.pressure, &axes.sliding_speed] {
                hash.update(&[u8::from(axis.is_some())]);
                let name = axis.as_deref().unwrap_or("");
                hash.update(&(name.len() as u64).to_le_bytes());
                hash.update(name.as_bytes());
            }
        }
        Ok(Self {
            contact,
            coefficients,
            coefficient_keys,
            local_axes,
            hardness,
            laws,
            maximum_recession_m,
            identity: hash.finalize(),
        })
    }

    /// Original hardness-limited normal owner.
    pub const fn contact(&self) -> &SourcePlasticContact {
        &self.contact
    }
    /// Exact local coordinates this owner evolves when querying wear claims.
    pub const fn local_coefficient_axes(&self) -> Option<&LocalWearAxes> {
        self.local_axes.as_ref()
    }
    /// Complete original coefficient receipt.
    pub const fn coefficient_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.coefficients
    }
    /// Exact selected original A/B coefficient keys.
    pub const fn coefficient_keys(&self) -> &[PropertyKey; 2] {
        &self.coefficient_keys
    }
    /// Complete strict per-body indentation-hardness receipts.
    pub const fn hardness_sources(&self) -> &[IndentationHardnessStatePoint; 2] {
        &self.hardness
    }
    /// Identity of all original sources, selected roles, and physical bounds.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Start without material removal or normal-contact history.
    pub fn initial_state(&self) -> SourcePlasticWearState {
        self.zero_wear_state(self.contact.initial_state())
    }

    /// Start wear from an original loaded contact that has no prior recession.
    /// Existing residual plastic indentation is retained. A profile already worn
    /// through another owner cannot be relabelled as unworn material.
    pub fn state_from_contact(
        &self,
        contact: SourcePlasticContactState,
    ) -> Result<SourcePlasticWearState, SourcePlasticWearError> {
        if contact.owner != self.contact.identity()
            || contact.contact_state().owner_identity() != self.contact.model().identity()
            || contact
                .contact_state()
                .recession_m()
                .iter()
                .any(|depth| *depth != 0.0)
        {
            return Err(invalid(
                "wear checkpoint requires original contact without existing recession",
            ));
        }
        Ok(self.zero_wear_state(contact))
    }

    /// Change approach while retaining cumulative recession and removed volumes.
    pub fn advance_normal(
        &self,
        before: &SourcePlasticWearState,
        approach_m: f64,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearNormalAdvance, SourcePlasticWearError> {
        self.admit(before)?;
        let next = self.contact.advance(&before.contact, approach_m, cx)?;
        self.admit_receded_profile(&next.candidate, &before.recession_m, cx)?;
        Ok(SourcePlasticWearNormalAdvance {
            candidate: SourcePlasticWearState {
                owner: self.identity,
                contact: next.candidate,
                volumes: before.volumes,
                recession_m: before.recession_m.clone(),
            },
            receipt: next.receipt,
        })
    }

    /// Remove material using the currently equilibrated pressure and local slip.
    ///
    /// Slip is unsigned *irreversible local relative sliding*, one distance per
    /// original cell. Elastic tangential motion is not admitted as sliding by
    /// implication. This first-order pressure/distance quadrature should be
    /// refined whenever removal appreciably redistributes contact pressure.
    /// All sources and approach stay fixed, plastic depth is retained, and no
    /// removed volume is relabelled as heat or plastic dissipation.
    pub fn advance_wear(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearAdvance, SourcePlasticWearError> {
        self.advance_wear_over_interval(before, local_unsigned_slip_m, None, cx)
    }

    /// Remove material over a supplied positive physical interval.
    ///
    /// Speed-dependent claims use `unsigned actual slip / duration` separately
    /// at every loaded sliding cell. This held mean-speed quadrature requires
    /// refinement when local pressure, speed, or temperature changes within a
    /// step. Detached and sticking cells remove nothing without a source query.
    pub fn advance_wear_during(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        duration_s: f64,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearAdvance, SourcePlasticWearError> {
        self.advance_wear_over_interval(before, local_unsigned_slip_m, Some(duration_s), cx)
    }

    fn advance_wear_over_interval(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        duration_s: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearAdvance, SourcePlasticWearError> {
        self.advance_wear_during_with_recession(
            before,
            local_unsigned_slip_m,
            duration_s,
            cx,
            |increment| {
                let geometry = self
                    .contact
                    .model()
                    .recede(&before.contact.contact, increment, cx)
                    .map_err(SourcePlasticContactError::from)?;
                Ok((geometry, ()))
            },
        )
        .map(|(advance, ())| advance)
    }

    fn advance_wear_at_speeds(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        duration_s: Option<f64>,
        prescribed_speeds_m_per_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearAdvance, SourcePlasticWearError> {
        self.advance_wear_at_pressure_and_speeds_with_recession(
            before,
            before.contact.contact_state().pressure_pa(),
            local_unsigned_slip_m,
            duration_s,
            Some(prescribed_speeds_m_per_s),
            cx,
            |increment| {
                let geometry = self
                    .contact
                    .model()
                    .recede(&before.contact.contact, increment, cx)
                    .map_err(SourcePlasticContactError::from)?;
                Ok((geometry, ()))
            },
        )
        .map(|(advance, ())| advance)
    }

    // The native boundary condition changes the equilibrium solve, not the
    // Archard law, source ownership, volume/profile accounting, or admissibility.
    // A trial geometry is accepted only with the complete per-body inventory.
    fn advance_wear_during_with_recession<T>(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        duration_s: Option<f64>,
        cx: &Cx<'_>,
        equilibrate: impl FnOnce(
            &[f64],
        )
            -> Result<(FiniteGapPlasticContinuation, T), SourcePlasticWearError>,
    ) -> Result<(SourcePlasticWearAdvance, T), SourcePlasticWearError> {
        self.advance_wear_at_pressure_with_recession(
            before,
            before.contact.contact_state().pressure_pa(),
            local_unsigned_slip_m,
            duration_s,
            cx,
            equilibrate,
        )
    }

    // The pressure is a QUADRATURE input, never a manufactured accepted contact
    // state. The implicit adapter gates it against the resulting native endpoint.
    // Every geometry solve still receives the original accepted predecessor.
    fn advance_wear_at_pressure_with_recession<T>(
        &self,
        before: &SourcePlasticWearState,
        pressure: &[f64],
        local_unsigned_slip_m: &[f64],
        duration_s: Option<f64>,
        cx: &Cx<'_>,
        equilibrate: impl FnOnce(
            &[f64],
        )
            -> Result<(FiniteGapPlasticContinuation, T), SourcePlasticWearError>,
    ) -> Result<(SourcePlasticWearAdvance, T), SourcePlasticWearError> {
        self.advance_wear_at_pressure_and_speeds_with_recession(
            before,
            pressure,
            local_unsigned_slip_m,
            duration_s,
            None,
            cx,
            equilibrate,
        )
    }

    fn advance_wear_at_pressure_and_speeds_with_recession<T>(
        &self,
        before: &SourcePlasticWearState,
        pressure: &[f64],
        local_unsigned_slip_m: &[f64],
        duration_s: Option<f64>,
        prescribed_speeds_m_per_s: Option<&[f64]>,
        cx: &Cx<'_>,
        equilibrate: impl FnOnce(
            &[f64],
        )
            -> Result<(FiniteGapPlasticContinuation, T), SourcePlasticWearError>,
    ) -> Result<(SourcePlasticWearAdvance, T), SourcePlasticWearError> {
        self.admit(before)?;
        checkpoint(cx)?;
        if pressure.len() != before.contact.contact_state().pressure_pa().len()
            || pressure
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
            || local_unsigned_slip_m.len() != pressure.len()
            || local_unsigned_slip_m
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "finite nonnegative local slip at every original cell",
            ));
        }
        let sampled = match prescribed_speeds_m_per_s {
            Some(speeds) => self.sample_local_coefficients_at_speeds(
                pressure,
                local_unsigned_slip_m,
                duration_s,
                Some(speeds),
                cx,
            )?,
            None => {
                self.sample_local_coefficients(pressure, local_unsigned_slip_m, duration_s, cx)?
            }
        };
        let grid = self.contact.model().config().grid;
        let area = grid.cell_width_m * grid.cell_depth_m;
        let mut increment = [vec![0.0; pressure.len()], vec![0.0; pressure.len()]];
        let mut recession = before.recession_m.clone();
        let mut relative_increment = vec![0.0; pressure.len()];
        let mut measure = 0.0;
        let mut compensation = 0.0;
        let mut geometric_volume = [0.0; 2];
        let mut geometric_compensation = [0.0; 2];
        let mut profile_roundoff_volume = [0.0; 2];
        for (index, (&p, &slip)) in pressure.iter().zip(local_unsigned_slip_m).enumerate() {
            if index % 64 == 0 {
                checkpoint(cx)?;
            }
            let cell_measure = (p * area) * slip;
            if !cell_measure.is_finite() || cell_measure < 0.0 {
                return Err(invalid("finite pressure-weighted local sliding"));
            }
            // Compensated summation protects heterogeneous pressure/sliding
            // fields. Equal local source queries share their original receipts.
            add_compensated(&mut measure, &mut compensation, cell_measure);
            for side in 0..2 {
                let coefficient = sampled
                    .as_ref()
                    .map_or(self.laws[side].wear_coefficient, |local| {
                        local.coefficient(side, index)
                    });
                let depth = (coefficient * (p / self.hardness[side].hardness_pa())) * slip;
                let old = recession[side][index];
                let total = old + depth;
                let retained_increment = total - old;
                if !depth.is_finite() || !total.is_finite() || depth < 0.0 {
                    return Err(invalid("finite nonnegative Archard recession"));
                }
                if depth > 0.0 && retained_increment == 0.0 {
                    return Err(invalid(
                        "wear increment is below the retained profile resolution",
                    ));
                }
                if total > self.maximum_recession_m[side] {
                    return Err(SourcePlasticContactError::OutsideApplicability {
                        quantity: if side == 0 {
                            "body A cumulative recession [m]"
                        } else {
                            "body B cumulative recession [m]"
                        },
                        value: total,
                        limit: self.maximum_recession_m[side],
                    }
                    .into());
                }
                let cell_volume = retained_increment * area;
                if !cell_volume.is_finite() || (retained_increment > 0.0 && cell_volume == 0.0) {
                    return Err(invalid(
                        "wear increment is outside volume arithmetic resolution",
                    ));
                }
                add_compensated(
                    &mut geometric_volume[side],
                    &mut geometric_compensation[side],
                    cell_volume,
                );
                // Four roundoff units cover addition/subtraction of depth and
                // multiplication by cell area; this is an arithmetic bound.
                profile_roundoff_volume[side] +=
                    4.0 * f64::EPSILON * ((old * area).abs() + (total * area).abs());
                increment[side][index] = retained_increment;
                recession[side][index] = total;
            }
            relative_increment[index] = increment[0][index] + increment[1][index];
            if !relative_increment[index].is_finite() {
                return Err(invalid("finite relative recession"));
            }
        }
        let interface = self
            .contact
            .normal_source()
            .elastic()
            .interface()
            .interface();
        let mut volumes = before.volumes;
        let mut advance_volume = |side: usize| -> Result<WearAdvance, SourcePlasticWearError> {
            if let Some(local) = &sampled {
                let mut refusal = None;
                let samples = pressure
                    .iter()
                    .zip(local_unsigned_slip_m)
                    .enumerate()
                    .scan((), |_, (cell, (&p, &slip))| {
                        if cell % 64 == 0 {
                            if let Err(error) = checkpoint(cx) {
                                refusal = Some(error);
                                return None;
                            }
                        }
                        Some(local.receipt.cell_source_indices[cell].map(|slot| {
                            LocalArchardSample {
                                wear_coefficient: local.laws[slot][side].wear_coefficient,
                                pressure_weighted_sliding_j: (p * area) * slip,
                            }
                        }))
                    })
                    .flatten();
                let result = ArchardLaw::advance_local_pressure_weighted_sliding(
                    interface,
                    &mut volumes[side],
                    samples,
                    self.hardness[side].hardness_pa(),
                );
                // The inventory is still private. A cancelled iterator may
                // not publish its partial native quadrature as an accepted step.
                if let Some(error) = refusal {
                    Err(error)
                } else {
                    result.map_err(Into::into)
                }
            } else {
                self.laws[side]
                    .advance_pressure_weighted_sliding(
                        interface,
                        &mut volumes[side],
                        measure,
                        self.hardness[side].hardness_pa(),
                    )
                    .map_err(Into::into)
            }
        };
        let wear = [advance_volume(0)?, advance_volume(1)?];
        for side in 0..2 {
            let scalar = wear[side].increment_m3;
            let retained_volume = volumes[side].volume_m3() - before.volumes[side].volume_m3();
            if scalar > 0.0 && retained_volume == 0.0 {
                return Err(invalid(
                    "wear increment is below the retained inventory resolution",
                ));
            }
            // Both formulas describe the same k*p*ds/H integral. Allow only
            // floating arithmetic error (product/division, compensated sums,
            // and represented profile/inventory additions), never a tunable
            // physical tolerance that could hide removed-volume inconsistency.
            let arithmetic_allowance = profile_roundoff_volume[side]
                + 16.0 * f64::EPSILON * scalar.abs().max(geometric_volume[side].abs());
            let inventory_allowance = arithmetic_allowance
                + 4.0
                    * f64::EPSILON
                    * (before.volumes[side].volume_m3().abs() + volumes[side].volume_m3().abs());
            if !arithmetic_allowance.is_finite()
                || !inventory_allowance.is_finite()
                || !geometric_volume[side].is_finite()
                || (scalar - geometric_volume[side]).abs() > arithmetic_allowance
                || (retained_volume - geometric_volume[side]).abs() > inventory_allowance
            {
                return Err(invalid(
                    "Archard inventory and retained recession disagree beyond roundoff",
                ));
            }
        }
        let (geometry, boundary_response) = equilibrate(&relative_increment)?;
        let contact = SourcePlasticContactState {
            owner: self.contact.identity(),
            contact: geometry.candidate,
        };
        self.admit_receded_profile(&contact, &recession, cx)?;
        let advance = SourcePlasticWearAdvance {
            source_identity: self.identity,
            candidate: SourcePlasticWearState {
                owner: self.identity,
                contact,
                volumes,
                recession_m: recession,
            },
            wear,
            recession_increment_m: increment,
            pressure_weighted_sliding_j: measure,
            coefficient_samples: sampled.map(|local| local.receipt),
            receipt: geometry.receipt,
        };
        Ok((advance, boundary_response))
    }

    /// Continue prescribed temperatures without resetting either body's wear.
    ///
    /// Every k and hardness key remains pinned to its exact original claim and
    /// test context. Interface k follows the original interface T; the two wear
    /// hardnesses follow their original bulk temperatures independently. All
    /// other external conditions stay fixed. This endpoint supplies no heat law.
    pub fn refresh_at_temperatures(
        &self,
        before: &SourcePlasticWearState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearThermalRefresh, SourcePlasticWearError> {
        self.refresh_temperature_with(before, material_cards, bulk_temperatures_k, cx, || {
            let refreshed = self.contact.refresh_at_temperatures(
                &before.contact,
                material_cards,
                bulk_temperatures_k,
                interface_temperature_k,
                cx,
            )?;
            let (contact, state, receipt) = refreshed.into_parts();
            Ok((contact, state, receipt, ()))
        })
        .map(|(refresh, ())| refresh)
    }

    // A fixed original source may only describe material still within its
    // declared homogeneous depth after removal. Inspect the farthest loaded-cell
    // corner exactly as the parent half-space envelope does.
    fn admit_receded_profile(
        &self,
        candidate: &SourcePlasticContactState,
        recession: &[Vec<f64>; 2],
        cx: &Cx<'_>,
    ) -> Result<(), SourcePlasticWearError> {
        self.contact
            .admit_candidate(candidate.contact_state(), cx)?;
        let grid = self.contact.model().config().grid;
        let mut radius = 0.0_f64;
        for (index, pressure) in candidate.contact_state().pressure_pa().iter().enumerate() {
            if *pressure > 0.0 {
                let x = (index % grid.cells_x) as f64 + 0.5 - 0.5 * grid.cells_x as f64;
                let y = (index / grid.cells_x) as f64 + 0.5 - 0.5 * grid.cells_y as f64;
                radius = radius.max(
                    ((x.abs() + 0.5) * grid.cell_width_m)
                        .hypot((y.abs() + 0.5) * grid.cell_depth_m),
                );
            }
        }
        for side in 0..2 {
            let maximum = recession[side].iter().copied().fold(0.0_f64, f64::max);
            let remaining = self.contact.envelope.half_space_depth_m[side] - maximum;
            if remaining <= 0.0 {
                return Err(invalid(
                    "recession exhausted the original homogeneous material depth",
                ));
            }
            super::within(
                "patch radius / remaining homogeneous depth",
                radius / remaining,
                self.contact.envelope.maximum_patch_to_depth[side],
            )?;
        }
        checkpoint(cx)
    }

    fn zero_wear_state(&self, contact: SourcePlasticContactState) -> SourcePlasticWearState {
        let count = self.contact.model().config().undeformed_gap_m.len();
        SourcePlasticWearState {
            owner: self.identity,
            contact,
            volumes: [WearState::default(); 2],
            recession_m: [vec![0.0; count], vec![0.0; count]],
        }
    }

    fn admit(&self, before: &SourcePlasticWearState) -> Result<(), SourcePlasticWearError> {
        if before.owner != self.identity || before.contact.owner != self.contact.identity() {
            return Err(invalid(
                "wear history requires its original ordered sources",
            ));
        }
        Ok(())
    }
}

fn invalid(message: &'static str) -> SourcePlasticWearError {
    SourcePlasticWearError::Invalid(message)
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SourcePlasticWearError> {
    cx.checkpoint().map_err(|_| {
        SourcePlasticContactError::from(FiniteGapPlasticError::Contact(
            FiniteGapContactError::Cancelled,
        ))
        .into()
    })
}

fn add_compensated(sum: &mut f64, compensation: &mut f64, value: f64) {
    let corrected = value - *compensation;
    let next = *sum + corrected;
    *compensation = (next - *sum) - corrected;
    *sum = next;
}
