//! Original Hertz/material sources around explicit changing-load shear memory.
//!
//! The fixed source admits geometry, both bulk cards and the ordered interface
//! once. This owner never refreshes coefficients while retaining history.
//! Pressure, force, footprint, stress, rate and unclassified query coordinates
//! are excluded because their validity cannot silently remain frozen when the
//! contact changes. The canonical material policy still admits fixed absolute
//! temperature, composition, moisture and concentration conditions.

use fs_blake3::ContentHash;
use fs_material::state_point::is_fixed_external_coordinate;
use fs_qty::semantic::QuantityKind;
use fs_tribo::mindlin::loading::{HertzMindlinContact, HertzMindlinParameters, HertzMindlinState};
use fs_tribo::{TangentialLaw, TangentialStep, TangentialStepInput, TriboError};

use super::{SourceMindlinContact, SourceMindlinError};
use crate::normal_patch::PointNormalPatchReceipt;
use crate::normal_patch::source::PointContactGeometry;

/// Isothermal changing-load contact retaining the complete original sources.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzMindlinContact {
    source: SourceMindlinContact,
    law: HertzMindlinContact,
    initial: SourceHertzMindlinState,
}

/// Source-bound shear memory and original Hertz normal initial/accepted energy.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzMindlinState {
    normal_model_identity: ContentHash,
    friction_state_identity: ContentHash,
    constitutive: HertzMindlinState,
    indentation_m: f64,
    normal_energy_j: f64,
}

impl SourceHertzMindlinState {
    /// Original analytic history, exposed without a mutation bypass.
    #[must_use]
    pub const fn constitutive_state(&self) -> &HertzMindlinState {
        &self.constitutive
    }
    /// Normal load from the original source response [N].
    #[must_use]
    pub fn normal_force_n(&self) -> f64 {
        self.constitutive.normal_force_n()
    }
    /// Original accepted geometric indentation [m].
    #[must_use]
    pub const fn indentation_m(&self) -> f64 {
        self.indentation_m
    }
    /// Normal energy from the original Hertz evaluator [J], separate from shear.
    #[must_use]
    pub const fn normal_energy_j(&self) -> f64 {
        self.normal_energy_j
    }
    /// Current contact radius [m], independently checked against the source.
    #[must_use]
    pub fn contact_radius_m(&self) -> f64 {
        self.constitutive.contact_radius_m()
    }
    /// Accepted tangential endpoint force [N].
    #[must_use]
    pub fn endpoint_force_n(&self) -> f64 {
        self.constitutive.endpoint_force_n()
    }
    /// Accumulated scalar shear displacement [m], including motion while open.
    #[must_use]
    pub const fn slip_m(&self) -> f64 {
        self.constitutive.slip_m()
    }
    /// Retained analytic intervals; zero after separation has paid for release.
    #[must_use]
    pub fn history_intervals(&self) -> usize {
        self.constitutive.history_intervals()
    }
}

/// One original sourced normal path and its physical shear-energy transfer.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceMindlinNormalStep {
    /// Candidate state; admission failures never publish it.
    pub state: SourceHertzMindlinState,
    /// Original Hertz footprint, pressure, authority and uncertainty receipt.
    pub normal_response: PointNormalPatchReceipt,
    /// Signed normal work from the original Hertz potential [J].
    pub normal_work_j: f64,
    /// Change in that same original normal potential [J].
    pub normal_stored_change_j: f64,
    /// Shear-storage change while macroscopic shear is held fixed [J].
    pub tangential_stored_change_j: f64,
    /// Physical microslip/release heat, excluding reversible normal work [J].
    pub dissipated_j: f64,
    /// Both storage changes plus heat minus original normal work [J].
    pub closure_residual_j: f64,
}

impl SourceHertzMindlinContact {
    /// Reuse an original admitted fixed source as the immutable material owner.
    ///
    /// The original positive indentation supplies the initial normal condition.
    /// Its Hertz energy is retained as initial storage, not newly supplied work;
    /// initial shear is zero. No previously loaded shear checkpoint is consumed
    /// or reset by this constructor. Runtime load changes use `normal_step`.
    pub fn from_fixed_source(source: SourceMindlinContact) -> Result<Self, SourceMindlinError> {
        let normal = source.normal_source().model().normal_state();
        let elastic = normal.elastic();
        let properties = elastic
            .surface_a_state()
            .resolved()
            .properties()
            .iter()
            .chain(elastic.surface_b_state().resolved().properties())
            .chain(normal.interface_state().properties())
            .chain(source.friction_source().properties());
        for property in properties {
            if !property.is_static_at_fixed_external_conditions() {
                return Err(invalid(
                    "changing-load Mindlin requires fixed external source conditions",
                ));
            }
            let receipt = &property.answer().receipt;
            for (axis, value) in &receipt.query_point {
                let quantity = receipt.axis_quantities.get(axis).copied();
                let pressure = quantity
                    .and_then(|q| q.semantic_type())
                    .is_some_and(|ty| ty.kind() == QuantityKind::Pressure);
                if pressure || !is_fixed_external_coordinate(axis, quantity, Some(*value)) {
                    return Err(invalid(
                        "changing-load Mindlin cannot freeze pressure or unowned load/footprint query conditions",
                    ));
                }
            }
        }
        let PointContactGeometry::SpherePlane { effective_radius_m } =
            source.normal_source().geometry()
        else {
            return Err(invalid(
                "changing-load Mindlin requires the original circular Hertz geometry",
            ));
        };
        let fixed = source.law().parameters();
        let law = HertzMindlinContact::new(HertzMindlinParameters {
            effective_radius_m,
            shear_moduli_pa: fixed.shear_moduli_pa,
            poisson_ratios: fixed.poisson_ratios,
            friction_coefficient: fixed.friction_coefficient,
            maximum_history_intervals: fixed.maximum_history_intervals,
        })?;
        let initialized = law.normal_step(
            &law.initial_state(),
            source.normal_response().normal_force_n,
        )?;
        match_original(&initialized.state, source.normal_response())?;
        let initial = SourceHertzMindlinState {
            normal_model_identity: source.normal_source().model().identity(),
            friction_state_identity: source.friction_source().identity(),
            constitutive: initialized.state,
            indentation_m: source.normal_response().approach_m,
            normal_energy_j: source.normal_response().reversible_energy_j,
        };
        Ok(Self {
            source,
            law,
            initial,
        })
    }

    /// All original immutable normal, bulk and ordered-interface evidence.
    #[must_use]
    pub const fn fixed_source(&self) -> &SourceMindlinContact {
        &self.source
    }
    /// Original analytic changing-load owner; no private material equation.
    #[must_use]
    pub const fn law(&self) -> &HertzMindlinContact {
        &self.law
    }
    /// Exact original isothermal condition required by each shear interval [K].
    #[must_use]
    pub const fn temperature_kelvin(&self) -> f64 {
        self.source.temperature_kelvin()
    }

    /// Move monotonically to an explicit normal indentation while shear is held.
    ///
    /// The declared duration admits the constant segment rate against the same
    /// original operating envelope. Its endpoint load, pressure and footprint
    /// come from the original normal evaluator. Normal and shear path ordering
    /// belongs to the caller; this is not a simultaneous body-motion solver.
    pub fn normal_step(
        &self,
        state: &SourceHertzMindlinState,
        indentation_m: f64,
        dt_s: f64,
    ) -> Result<SourceMindlinNormalStep, SourceMindlinError> {
        if !self.owns_state(state) {
            return Err(invalid(
                "changing-load Mindlin checkpoint requires original source identities",
            ));
        }
        let normal = self.source.normal_source();
        normal.validate_kinematics(state.indentation_m, 0.0, dt_s)?;
        let rate = (indentation_m - state.indentation_m) / dt_s;
        let response = normal.evaluate(
            "source-mindlin/changing-normal-patch",
            indentation_m,
            rate,
            dt_s,
        )?;
        let step = self
            .law
            .normal_step(&state.constitutive, response.normal_force_n)?;
        match_original(&step.state, &response)?;
        let change = response.reversible_energy_j - state.normal_energy_j;
        let state = SourceHertzMindlinState {
            normal_model_identity: state.normal_model_identity,
            friction_state_identity: state.friction_state_identity,
            constitutive: step.state,
            indentation_m,
            normal_energy_j: response.reversible_energy_j,
        };
        Ok(SourceMindlinNormalStep {
            state,
            normal_response: response,
            normal_work_j: change,
            normal_stored_change_j: change,
            tangential_stored_change_j: step.tangential_stored_change_j,
            dissipated_j: step.dissipated_j,
            closure_residual_j: step.closure_residual_j,
        })
    }

    fn owns_state(&self, state: &SourceHertzMindlinState) -> bool {
        state.normal_model_identity == self.source.normal_source().model().identity()
            && state.friction_state_identity == self.source.friction_source().identity()
            && self.law.stored_energy_j(&state.constitutive).is_finite()
    }
}

impl TangentialLaw for SourceHertzMindlinContact {
    type State = SourceHertzMindlinState;

    fn model_id(&self) -> &str {
        "source-hertz-mindlin-explicit-normal-path-v1"
    }
    fn initial_state(&self) -> Self::State {
        self.initial.clone()
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
                field: "changing-load Mindlin requires original source history",
            });
        }
        if input.contact_temperature_k != Some(self.source.temperature_kelvin())
            || input.nominal_pressure_pa.is_some()
        {
            return Err(TriboError::InvalidInput {
                field: "changing-load Mindlin requires frozen temperature and no pressure refresh",
            });
        }
        self.source
            .normal_source()
            .validate_kinematics(state.indentation_m, 0.0, input.dt_s)
            .map_err(|_| TriboError::InvalidInput {
                field: "changing-load Mindlin shear exceeds original normal envelope",
            })?;
        let candidate = self.law.advance(&state.constitutive, input)?;
        Ok(TangentialStep {
            traction_n: candidate.traction_n,
            slip_increment_m: candidate.slip_increment_m,
            state: SourceHertzMindlinState {
                normal_model_identity: state.normal_model_identity,
                friction_state_identity: state.friction_state_identity,
                constitutive: candidate.state,
                indentation_m: state.indentation_m,
                normal_energy_j: state.normal_energy_j,
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

fn invalid(message: &'static str) -> SourceMindlinError {
    SourceMindlinError::Invalid(message)
}

// These are two evaluations of the SAME admitted elastic law, not independent
// uncertain material estimates. Permit only arithmetic/reduction roundoff;
// original domain/pressure/rate/geometry guards are never relaxed.
fn match_original(
    state: &HertzMindlinState,
    original: &PointNormalPatchReceipt,
) -> Result<(), SourceMindlinError> {
    for (actual, expected) in [
        (state.normal_force_n(), original.normal_force_n),
        (state.contact_radius_m(), original.patch_radius_m),
        (state.indentation_m(), original.approach_m),
        (state.normal_energy_j(), original.reversible_energy_j),
    ] {
        let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
        if !actual.is_finite()
            || !expected.is_finite()
            || (actual - expected).abs() > 512.0 * f64::EPSILON * scale
        {
            return Err(invalid(
                "changing-load Mindlin must retain the original Hertz response",
            ));
        }
    }
    Ok(())
}
