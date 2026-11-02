//! Irreversible local sliding and original-source Archard fretting wear.
//!
//! The fixed-load scalar contact supplies I = integral p(r) ds(r) dA from
//! actual local sliding. At constant positive mu this is dissipated work / mu;
//! the explicitly frictionless limit instead uses N times actual travel.
//! Source wear is k_i I / H_i for each original ordered surface. Elastic
//! displacement and recoverable work never enter I. This is a small-removal,
//! fixed-geometry diagnostic, not profile recession, debris or damage energy.
//!
//! The local Archard/constant-Coulomb equivalence is described by Willert
//! (2021), section 2: <https://doi.org/10.3390/sym13091582>.

use super::*;
use fs_matdb::PropertyKey;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::Pressure;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_tribo::mindlin::surface::sliding::{MindlinSurfaceAdvance, MindlinSurfaceSliding};
use fs_tribo::{ArchardLaw, WearAdvance, WearState};

/// Simultaneous changing normal load, shear history and original surface wear.
#[path = "fretting/loading.rs"]
pub mod loading;

/// A single admitted source step, retaining its exact native sliding path.
/// The borrowed original state prevents an observation from silently spanning
/// several reversals; consume this receipt to obtain the source-bound candidate.
pub struct SourceMindlinSurfaceAdvance<'a> {
    source: &'a SourceMindlinContact,
    native: MindlinSurfaceAdvance<'a>,
}

impl SourceMindlinContact {
    /// Advance once under original source admission, with local sliding fields.
    pub fn advance_with_surface_sliding<'a>(
        &'a self,
        before: &'a SourceMindlinState,
        input: &TangentialStepInput,
    ) -> Result<SourceMindlinSurfaceAdvance<'a>, SourceMindlinError> {
        self.admit_step(before, input)?;
        Ok(SourceMindlinSurfaceAdvance {
            source: self,
            native: self
                .law
                .advance_with_surface_sliding(&before.constitutive, input)?,
        })
    }
}

impl SourceMindlinSurfaceAdvance<'_> {
    /// Read-only native effort, energy and solved motion for this exact step.
    pub fn tangential_step(&self) -> &TangentialStep<MindlinState> {
        self.native.tangential_step()
    }

    /// Pressure-weighted unsigned local sliding integral [N m].
    pub fn pressure_weighted_sliding_j(&self) -> f64 {
        self.native.pressure_weighted_sliding_j()
    }

    /// Local irreversible sliding and dissipated work per contact area.
    pub fn surface_sliding(
        &self,
        radius_m: f64,
    ) -> Result<MindlinSurfaceSliding, SourceMindlinError> {
        Ok(self.native.surface_sliding(radius_m)?)
    }

    /// Transfer the accepted native history back into its original source owner.
    pub fn into_tangential_step(self) -> TangentialStep<SourceMindlinState> {
        self.source.bind_step(self.native.into_tangential_step())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum SourceFrettingError {
    Source(SourceMindlinError),
    Path(super::loading::path::SourceMindlinPathError),
    Cancelled,
}
impl From<SourceMindlinError> for SourceFrettingError {
    fn from(error: SourceMindlinError) -> Self {
        Self::Source(error)
    }
}
impl From<TriboError> for SourceFrettingError {
    fn from(error: TriboError) -> Self {
        Self::Source(SourceMindlinError::Tribo(error))
    }
}
impl From<super::loading::path::SourceMindlinPathError> for SourceFrettingError {
    fn from(error: super::loading::path::SourceMindlinPathError) -> Self {
        use super::loading::path::SourceMindlinPathError;
        use fs_tribo::mindlin::loading::path::MindlinPathError;
        match error {
            SourceMindlinPathError::Source(error) => Self::Source(error),
            SourceMindlinPathError::Path(MindlinPathError::Cancelled) => Self::Cancelled,
            error => Self::Path(error),
        }
    }
}
impl fmt::Display for SourceFrettingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "source Mindlin fretting: {self:?}")
    }
}
impl std::error::Error for SourceFrettingError {}

// Bundle identities pin every original claim and quantity context. The small
// canonical property indices additionally pin the chosen A/B roles within each
// bundle; no property-key strings or source cards are cloned per step.
#[derive(Clone, Copy, Debug, PartialEq)]
struct WearOwner {
    coefficients: ContentHash,
    hardness: [ContentHash; 2],
    selected: [usize; 4],
}

/// Retained contact history and cumulative wear accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceMindlinWearState {
    owner: WearOwner,
    contact: SourceMindlinState,
    volumes: [WearState; 2],
}
impl SourceMindlinWearState {
    pub const fn contact_state(&self) -> &SourceMindlinState {
        &self.contact
    }
    /// Cumulative removed volumes in original A/B order [m^3].
    pub fn volumes_m3(&self) -> [f64; 2] {
        self.volumes.map(|state| state.volume_m3())
    }
}

/// Source-selected, spatially constant k/H around a fixed-load Mindlin contact.
/// The complete interface and both original material cards and queries remain
/// fixed, including humidity, roughness/history and any hardness test context.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceMindlinWear {
    contact: SourceMindlinContact,
    coefficients: ResolvedInterfaceStatePoint,
    hardness: [ResolvedMaterialStatePoint; 2],
    laws: [ArchardLaw; 2],
    hardness_pa: [f64; 2],
    owner: WearOwner,
}

/// One candidate contains complete contact memory and both wear inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceMindlinWearAdvance {
    pub tangential: TangentialStep<SourceMindlinWearState>,
    pub wear: [WearAdvance; 2],
    /// Integral of actual unsigned local sliding weighted by pressure [N m].
    pub pressure_weighted_sliding_j: f64,
}

impl SourceMindlinWear {
    /// Bind dimensionless Archard k and positive indentation hardness in Pa.
    /// Exact property keys retain source test contexts and ordered surface roles.
    /// No coefficient, hardness conversion or material calibration is inferred.
    pub fn new(
        contact: SourceMindlinContact,
        coefficients: ResolvedInterfaceStatePoint,
        coefficient_keys: [PropertyKey; 2],
        hardness: [ResolvedMaterialStatePoint; 2],
        hardness_keys: [PropertyKey; 2],
    ) -> Result<Self, SourceFrettingError> {
        let friction = contact.friction_source();
        if coefficients.card_identity() != friction.card_identity()
            || coefficients.surface_materials() != friction.surface_materials()
            || coefficients.query_point() != friction.query_point()
            || coefficients.axis_quantities() != friction.axis_quantities()
        {
            return Err(invalid(
                "wear requires the original ordered interface and frozen query",
            ));
        }
        let elastic = contact.normal.model().normal_state().elastic();
        let bulk = [
            elastic.surface_a_state().resolved(),
            elastic.surface_b_state().resolved(),
        ];
        let mut laws = [ArchardLaw {
            wear_coefficient: 0.0,
        }; 2];
        let mut hardness_pa = [0.0; 2];
        let mut selected = [0; 4];
        for side in 0..2 {
            let k = coefficients
                .property_by_key(&coefficient_keys[side])
                .ok_or_else(|| invalid("missing exact original wear coefficient"))?;
            let h = hardness[side]
                .property_by_key(&hardness_keys[side])
                .ok_or_else(|| invalid("missing exact original hardness and test context"))?;
            let quantity = hardness_keys[side].quantity();
            let physical_hardness = quantity.dims() == Pressure::DIMS
                && quantity.semantic_type().is_none_or(|kind| {
                    matches!(kind.kind(), QuantityKind::Hardness(_))
                        && kind.form() == ValueForm::Static
                });
            if coefficient_keys[side].quantity() != QuantitySpec::dimensional(Dims::NONE)
                || !k.value_si().is_finite()
                || k.value_si() < 0.0
                || !k.is_static_at_fixed_external_conditions()
                || !physical_hardness
                || !h.value_si().is_finite()
                || h.value_si() <= 0.0
                || !h.is_static_at_fixed_external_conditions()
                || hardness[side].card_identity() != bulk[side].card_identity()
                || hardness[side].material() != bulk[side].material()
                || hardness[side].query_point() != bulk[side].query_point()
                || hardness[side].axis_quantities() != bulk[side].axis_quantities()
            {
                return Err(invalid(
                    "wear requires static nonnegative k and original positive Pa hardness",
                ));
            }
            selected[side] = coefficients
                .properties()
                .iter()
                .position(|property| property == k)
                .ok_or_else(|| {
                    invalid("wear coefficient is not retained by its original bundle")
                })?;
            selected[2 + side] = hardness[side]
                .properties()
                .iter()
                .position(|property| property == h)
                .ok_or_else(|| invalid("hardness is not retained by its original bundle"))?;
            laws[side].wear_coefficient = k.value_si();
            hardness_pa[side] = h.value_si();
        }
        let owner = WearOwner {
            coefficients: coefficients.identity(),
            hardness: hardness.each_ref().map(|source| source.identity()),
            selected,
        };
        Ok(Self {
            contact,
            coefficients,
            hardness,
            laws,
            hardness_pa,
            owner,
        })
    }

    pub const fn contact(&self) -> &SourceMindlinContact {
        &self.contact
    }
    pub const fn coefficient_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.coefficients
    }
    pub const fn hardness_sources(&self) -> &[ResolvedMaterialStatePoint; 2] {
        &self.hardness
    }

    /// Start from zero shear memory and zero removed volume.
    pub fn initial_state(&self) -> SourceMindlinWearState {
        SourceMindlinWearState {
            owner: self.owner,
            contact: self.contact.initial_state(),
            volumes: [WearState::default(); 2],
        }
    }

    /// Attach explicitly supplied existing inventories to an original checkpoint.
    pub fn state_from_contact(
        &self,
        contact: SourceMindlinState,
        volumes: [WearState; 2],
    ) -> Result<SourceMindlinWearState, SourceFrettingError> {
        if !self.contact.stored_energy_j(&contact).is_finite() {
            return Err(invalid(
                "fretting checkpoint requires original contact history",
            ));
        }
        Ok(SourceMindlinWearState {
            owner: self.owner,
            contact,
            volumes,
        })
    }

    /// Solve contact once and accumulate wear from its actual local sliding.
    /// Frozen N, mu, k/H and monotone scalar motion are the admitted phase.
    /// Arbitrary reversals are successive calls. Cancellation or either surface
    /// refusing leaves the complete accepted history and both volumes untouched.
    /// Wear is a volume diagnostic and does not subtract energy from the native
    /// dissipation/heat port or alter the retained Hertz geometry.
    pub fn advance_candidate(
        &self,
        state: &SourceMindlinWearState,
        input: &TangentialStepInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceMindlinWearAdvance, SourceFrettingError> {
        poll(cancelled)?;
        if state.owner != self.owner {
            return Err(invalid(
                "fretting checkpoint requires its original ordered k/H claims",
            ));
        }
        let phase = self
            .contact
            .advance_with_surface_sliding(&state.contact, input)?;
        let measure = phase.pressure_weighted_sliding_j();
        poll(cancelled)?;
        let (volumes, wear) = self.advance_volumes(state.volumes, measure)?;
        let candidate = phase.into_tangential_step();
        poll(cancelled)?;
        Ok(SourceMindlinWearAdvance {
            tangential: TangentialStep {
                traction_n: candidate.traction_n,
                slip_increment_m: candidate.slip_increment_m,
                state: SourceMindlinWearState {
                    owner: self.owner,
                    contact: candidate.state,
                    volumes,
                },
                regime: candidate.regime,
                stored_change_j: candidate.stored_change_j,
                dissipated_j: candidate.dissipated_j,
                closure_residual_j: candidate.closure_residual_j,
                compatibility_residual_m: candidate.compatibility_residual_m,
                iterations: candidate.iterations,
            },
            wear,
            pressure_weighted_sliding_j: measure,
        })
    }

    // Both fixed and changing-load owners pass an admitted physical sliding
    // integral. This common leaf keeps source roles and cumulative arithmetic
    // identical, without performing a second constitutive solve or material query.
    fn advance_volumes(
        &self,
        mut volumes: [WearState; 2],
        measure: f64,
    ) -> Result<([WearState; 2], [WearAdvance; 2]), SourceFrettingError> {
        let interface = self
            .contact
            .normal
            .model()
            .normal_state()
            .elastic()
            .interface()
            .interface();
        let first = self.laws[0].advance_pressure_weighted_sliding(
            interface,
            &mut volumes[0],
            measure,
            self.hardness_pa[0],
        )?;
        let second = self.laws[1].advance_pressure_weighted_sliding(
            interface,
            &mut volumes[1],
            measure,
            self.hardness_pa[1],
        )?;
        Ok((volumes, [first, second]))
    }
}

fn invalid(message: &'static str) -> SourceFrettingError {
    SourceMindlinError::Invalid(message).into()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SourceFrettingError> {
    if cancelled() {
        Err(SourceFrettingError::Cancelled)
    } else {
        Ok(())
    }
}
