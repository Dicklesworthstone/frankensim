//! Source-backed, constant-load Cattaneo--Mindlin tangential contact.
//!
//! The original normal owner supplies the circular Hertz footprint and load.
//! Both original bulk states supply G = E / (2 (1 + nu)); the local Coulomb
//! coefficient is explicitly selected from the ORDERED interface card. No bulk
//! material name or normal-contact coefficient supplies a friction default.
//!
//! This explicitly selects the classical decoupled, isotropic Hertz/Mindlin
//! approximation. Normal/tangential elastic coupling, changing load, elliptic
//! contacts, adhesion, thermal evolution and nonproportional slip are excluded.
//! The whole source query is frozen, including its external conditions. A
//! different temperature or load needs a physical history-transfer owner;
//! constructing a new unloaded contact does not perform such a transfer.

use core::fmt;

use fs_blake3::ContentHash;
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_qty::{Dims, QuantitySpec};
use fs_tribo::mindlin::{MindlinDeresiewiczPoint, MindlinPointParameters, MindlinState};
use fs_tribo::{TangentialLaw, TangentialStep, TangentialStepInput, TriboError};

use super::BoundNormalContactLaw;
use crate::normal_patch::source::{PointContactGeometry, SourcePointContact};
use crate::normal_patch::{NormalPatchError, PointNormalPatchReceipt};

/// Sourced normal loading, unloading and recontact with retained shear history.
#[path = "mindlin/loading.rs"]
pub mod loading;

/// Exact source-backed gross sliding and normal-axis spin on the circular patch.
#[path = "mindlin/sliding.rs"]
pub mod sliding;

/// Original bulk, interface and patch sources around one stateful friction law.
///
/// The scalar shear coordinate is tangent axis 0. Its fixed physical direction
/// belongs to the mechanical attachment; axis 1 must have zero released motion.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceMindlinContact {
    normal: SourcePointContact,
    patch: PointNormalPatchReceipt,
    friction: ResolvedInterfaceStatePoint,
    law: MindlinDeresiewiczPoint,
    temperature_kelvin: f64,
}

/// Retained tangential memory bound to the original normal and friction sources.
/// Equal numerical coefficients do not authorize transfer to a different card,
/// ordered interface history, query, or selected property-use receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceMindlinState {
    normal_model_identity: ContentHash,
    friction_state_identity: ContentHash,
    constitutive: MindlinState,
}

impl SourceMindlinState {
    /// Read-only constitutive checkpoint, including the complete analytic memory.
    #[must_use]
    pub const fn constitutive_state(&self) -> &MindlinState {
        &self.constitutive
    }

    /// Accepted endpoint force [N], distinct from the held path-average effort.
    #[must_use]
    pub const fn endpoint_force_n(&self) -> f64 {
        self.constitutive.endpoint_force_n()
    }

    /// Accepted displacement along the fixed first tangent [m].
    #[must_use]
    pub const fn slip_m(&self) -> f64 {
        self.constitutive.slip_m()
    }

    /// Number of retained analytic history intervals.
    #[must_use]
    pub fn history_intervals(&self) -> usize {
        self.constitutive.history_intervals()
    }
}

/// A source, normal geometry or constitutive admission refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceMindlinError {
    /// The selected interface property could not be resolved.
    Material(MaterialStatePointError),
    /// The original normal contact refused the footprint or operating point.
    Normal(NormalPatchError),
    /// The native tangential law refused its parameters.
    Tribo(TriboError),
    /// A requested source or model lies outside this fixed-load rung.
    Invalid(&'static str),
}

impl fmt::Display for SourceMindlinError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "source-backed Mindlin contact: {self:?}")
    }
}
impl std::error::Error for SourceMindlinError {}
impl From<MaterialStatePointError> for SourceMindlinError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<NormalPatchError> for SourceMindlinError {
    fn from(error: NormalPatchError) -> Self {
        Self::Normal(error)
    }
}
impl From<TriboError> for SourceMindlinError {
    fn from(error: TriboError) -> Self {
        Self::Tribo(error)
    }
}

impl SourceMindlinContact {
    /// Resolve an explicit single Coulomb coefficient and retain its full source.
    ///
    /// `friction_requirement` selects the local coefficient used both at the
    /// microslip boundary and during gross sliding. This does not combine
    /// unequal static/kinetic coefficients or silently select one of them.
    /// The property must be dimensionless and valid at fixed external source
    /// conditions; rate-, displacement- and history-bounded claims refuse.
    /// `indentation_m` is the accepted, positive normal indentation to hold.
    pub fn resolve(
        normal: SourcePointContact,
        indentation_m: f64,
        friction_requirement: ScalarPropertyRequirement,
        selection: MaterialPropertySelection,
        maximum_history_intervals: usize,
    ) -> Result<Self, SourceMindlinError> {
        if !matches!(normal.geometry(), PointContactGeometry::SpherePlane { .. })
            || normal.model().law() != BoundNormalContactLaw::ElasticHertz
        {
            return Err(SourceMindlinError::Invalid(
                "Mindlin source requires an elastic circular Hertz point contact",
            ));
        }
        if friction_requirement.quantity() != QuantitySpec::dimensional(Dims::NONE) {
            return Err(SourceMindlinError::Invalid(
                "local Coulomb coefficient must have an explicit dimensionless requirement",
            ));
        }
        let state = normal.model().normal_state();
        let elastic = state.elastic();
        let friction = resolve_interface_state_point(
            elastic.interface().card(),
            &state.interface_state().query(),
            &[friction_requirement],
            selection,
        )?;
        if friction
            .properties()
            .iter()
            .any(|p| !p.is_static_at_fixed_external_conditions())
            || elastic
                .surface_a_state()
                .resolved()
                .properties()
                .iter()
                .chain(elastic.surface_b_state().resolved().properties())
                .any(|p| !p.is_static_at_fixed_external_conditions())
        {
            return Err(SourceMindlinError::Invalid(
                "Mindlin material and friction sources must admit fixed external conditions",
            ));
        }
        let temperature_kelvin =
            elastic
                .state_coordinate("T")
                .ok_or(SourceMindlinError::Invalid(
                    "original normal source has no contact temperature",
                ))?;
        let patch = normal.evaluate(
            "source-mindlin/fixed-normal-patch",
            indentation_m,
            0.0,
            normal.envelope().maximum_step_s,
        )?;
        let solids = [elastic.surface_a_state(), elastic.surface_b_state()];
        let poisson_ratios = solids.map(|s| s.poisson_ratio());
        let shear_moduli_pa =
            solids.map(|s| (0.5 * s.young_modulus_pa()) / (1.0 + s.poisson_ratio()));
        let law = MindlinDeresiewiczPoint::new(MindlinPointParameters {
            normal_force_n: patch.normal_force_n,
            contact_radius_m: patch.patch_radius_m,
            shear_moduli_pa,
            poisson_ratios,
            friction_coefficient: friction.properties()[0].value_si(),
            maximum_history_intervals,
        })?;
        Ok(Self {
            normal,
            patch,
            friction,
            law,
            temperature_kelvin,
        })
    }

    /// Complete original normal model, both bulk states and ordered interface.
    #[must_use]
    pub const fn normal_source(&self) -> &SourcePointContact {
        &self.normal
    }

    /// Original Hertz response at the frozen normal indentation.
    #[must_use]
    pub const fn normal_response(&self) -> &PointNormalPatchReceipt {
        &self.patch
    }

    /// Selected coefficient, query, original uncertainty and property-use receipt.
    #[must_use]
    pub const fn friction_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.friction
    }

    /// Native stateful law with material-derived tangential stiffness.
    #[must_use]
    pub const fn law(&self) -> &MindlinDeresiewiczPoint {
        &self.law
    }

    /// Frozen original contact temperature required at each step [K].
    #[must_use]
    pub const fn temperature_kelvin(&self) -> f64 {
        self.temperature_kelvin
    }

    fn owns_state(&self, state: &SourceMindlinState) -> bool {
        state.normal_model_identity == self.normal.model().identity()
            && state.friction_state_identity == self.friction.identity()
    }
}

impl TangentialLaw for SourceMindlinContact {
    type State = SourceMindlinState;

    fn model_id(&self) -> &str {
        "source-cattaneo-mindlin-fixed-load"
    }
    fn initial_state(&self) -> Self::State {
        SourceMindlinState {
            normal_model_identity: self.normal.model().identity(),
            friction_state_identity: self.friction.identity(),
            constitutive: self.law.initial_state(),
        }
    }
    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        if !self.owns_state(state) {
            f64::NAN
        } else {
            self.law.stored_energy_j(&state.constitutive)
        }
    }
    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        if !self.owns_state(state) {
            return Err(TriboError::InvalidInput {
                field: "Mindlin checkpoint requires its original normal and friction sources",
            });
        }
        if input.contact_temperature_k != Some(self.temperature_kelvin)
            || input.nominal_pressure_pa.is_some()
        {
            return Err(TriboError::InvalidInput {
                field: "Mindlin source requires its frozen temperature and no nominal-pressure refresh",
            });
        }
        self.normal
            .validate_kinematics(self.patch.approach_m, 0.0, input.dt_s)
            .map_err(|_| TriboError::InvalidInput {
                field: "Mindlin step exceeds the original normal source envelope",
            })?;
        let candidate = self.law.advance(&state.constitutive, input)?;
        Ok(TangentialStep {
            traction_n: candidate.traction_n,
            slip_increment_m: candidate.slip_increment_m,
            state: SourceMindlinState {
                normal_model_identity: state.normal_model_identity,
                friction_state_identity: state.friction_state_identity,
                constitutive: candidate.state,
            },
            regime: candidate.regime,
            stored_change_j: candidate.stored_change_j,
            dissipated_j: candidate.dissipated_j,
            closure_residual_j: candidate.closure_residual_j,
            compatibility_residual_m: candidate.compatibility_residual_m,
            iterations: candidate.iterations,
        })
    }
}
