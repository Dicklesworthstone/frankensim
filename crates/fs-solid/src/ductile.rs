//! Ductile-initiation gates on the original inelastic continuum solvers.
//!
//! Each call proposes exactly one physical increment with the existing Newton
//! and return-map owners. Mechanics and ductility are returned together only
//! when every material point remains below initiation. A refused candidate
//! leaves the borrowed displacement, plastic history and ductility unchanged;
//! callers can refine the load/time interval from that same accepted state.
//!
//! The criterion uses accumulated equivalent PLASTIC strain, including
//! reversed plastic flow, and physical endpoint stress. Creep accumulation,
//! total strain and effective stress relative to a backstress are not
//! substitutes. Endpoint-held quadrature is first order for varying conditions.
//! A crossing brackets initiation; continuing needs a regularized fracture law.
//! No stress degradation, fracture energy, element deletion or topology is
//! inferred. Timed small-strain laws must have been prepared for the same
//! physical duration supplied to this adapter.

use fs_exec::Cx;
use fs_material::creep::plastic::NortonJ2State;
use fs_material::creep::thermal::expansion::ThermalCreepState;
use fs_material::ductile::{DuctileError, DuctileInitiationLaw, DuctileState, DuctileStep};
use fs_material::plastic::chaboche::ChabocheState;
use fs_material::plastic::finite::FiniteJ2State;
use fs_material::plastic::reduction::ReducedState;
use fs_material::plastic::J2State;
use fs_material::{SmallStrainLaw, Voigt};

use crate::plastic2::{SmallStrainPlaneError, SmallStrainPlaneProblem, SmallStrainPlaneSolution};
use crate::plastic3::dynamics::BackwardEulerStep;
use crate::plastic3::hex::{SmallStrainHexProblem, SmallStrainHexSolution};
use crate::plastic3::{
    SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem, SmallStrainTetSolution,
};

/// Finite-deformation, thermoplastic and inertial consumers of the same gate.
pub mod finite;

/// Progressive J2 softening with adaptive native tetrahedral continuation.
pub mod evolution;

/// Read original equivalent plastic history without re-running a law.
/// Implementations must exclude creep and retain accumulated reverse flow.
pub trait PlasticStrainHistory {
    /// Dimensionless accumulated equivalent plastic strain.
    fn plastic_strain_accumulation(&self) -> f64;
}

/// Explicit material admission for a ductile law used with source-resolved
/// mechanics. Authored numerical laws retain their declared assumptions;
/// sourced laws must match the original card, named state and complete fixed
/// context. There is no blanket unchecked implementation for sourced laws.
pub trait DuctileMaterialBinding: DuctileInitiationLaw {
    /// Admit the complete source identity and fixed quantity-typed context.
    /// Authored numerical laws may explicitly carry no material-card binding.
    fn matches_material(
        &self, material: &fs_material::state_point::ResolvedMaterialStatePoint,
    ) -> bool;
}

impl DuctileMaterialBinding for fs_material::ductile::JohnsonCookInitiation {
    fn matches_material(&self, _: &fs_material::state_point::ResolvedMaterialStatePoint) -> bool {
        true
    }
}
impl DuctileMaterialBinding for fs_material::ductile::TabulatedDuctileLocus {
    fn matches_material(&self, _: &fs_material::state_point::ResolvedMaterialStatePoint) -> bool {
        true
    }
}
impl DuctileMaterialBinding for fs_material::ductile::source::ResolvedDuctileLocus {
    fn matches_material(&self, material: &fs_material::state_point::ResolvedMaterialStatePoint) -> bool {
        self.source().matches_material(material)
    }
}
impl<L: DuctileMaterialBinding + ?Sized> DuctileMaterialBinding for &L {
    fn matches_material(&self, material: &fs_material::state_point::ResolvedMaterialStatePoint) -> bool {
        (**self).matches_material(material)
    }
}
impl PlasticStrainHistory for J2State {
    fn plastic_strain_accumulation(&self) -> f64 { self.alpha }
}
impl<const N: usize> PlasticStrainHistory for ChabocheState<N> {
    fn plastic_strain_accumulation(&self) -> f64 { self.plastic.alpha }
}
impl PlasticStrainHistory for NortonJ2State {
    fn plastic_strain_accumulation(&self) -> f64 { self.plastic.alpha }
}
impl PlasticStrainHistory for FiniteJ2State {
    fn plastic_strain_accumulation(&self) -> f64 { self.alpha }
}
impl<S: PlasticStrainHistory> PlasticStrainHistory for ReducedState<S> {
    fn plastic_strain_accumulation(&self) -> f64 { self.history.plastic_strain_accumulation() }
}
impl<S: PlasticStrainHistory> PlasticStrainHistory for ThermalCreepState<S> {
    fn plastic_strain_accumulation(&self) -> f64 { self.mechanical.plastic_strain_accumulation() }
}

/// Immutable calibration assignment in native integration-point order.
/// Source-backed models retain their original cards without queries or card
/// clones in the point loop. Heterogeneous fields cover every point.
#[derive(Debug)]
pub enum DuctileLawField<'a, L: DuctileInitiationLaw> {
    /// One original calibration for all integration points.
    Uniform(&'a L),
    /// A calibration per point in the native solution's ordering.
    PerPoint(&'a [L]),
}
impl<L: DuctileInitiationLaw> Copy for DuctileLawField<'_, L> {}
impl<L: DuctileInitiationLaw> Clone for DuctileLawField<'_, L> {
    fn clone(&self) -> Self { *self }
}
impl<'a, L: DuctileInitiationLaw> DuctileLawField<'a, L> {
    fn at(self, point: usize) -> &'a L {
        match self { Self::Uniform(law) => law, Self::PerPoint(laws) => &laws[point] }
    }
}

/// Actual endpoint temperatures [K], in the same point order as mechanics.
#[derive(Clone, Copy, Debug)]
pub enum DuctileTemperatures<'a> {
    /// Spatially uniform physical temperature.
    Uniform(f64),
    /// A temperature per native integration point.
    PerPoint(&'a [f64]),
    /// One temperature per native eight-point hexahedron, in cell order.
    /// This broadcasts without expanding an element field into a point array.
    PerHexElement(&'a [f64]),
}
impl DuctileTemperatures<'_> {
    fn at(self, point: usize) -> f64 {
        match self {
            Self::Uniform(t) => t,
            Self::PerPoint(t) => t[point],
            Self::PerHexElement(t) => t[point / 8],
        }
    }
}

/// Physical interval used to form delta_alpha / duration_s. A quasistatic load
/// factor is not a time; even a static solve needs this explicit clock.
#[derive(Clone, Copy, Debug)]
pub struct DuctileStepConditions<'a> {
    /// Positive physical duration of exactly one mechanical increment, seconds.
    pub duration_s: f64,
    /// Actual endpoint temperatures, not ambient reservoir temperatures.
    pub temperatures_k: DuctileTemperatures<'a>,
}

/// Retained initiation history tied to an immutable calibration field.
/// Clone for restart or alternative proposals. Material, mesh, quadrature or
/// point-order changes require explicit history transfer, never a reset.
#[derive(Debug)]
pub struct DuctileHistory<'a, L: DuctileInitiationLaw> {
    laws: DuctileLawField<'a, L>,
    points: Vec<DuctileState>,
}
impl<L: DuctileInitiationLaw> Clone for DuctileHistory<'_, L> {
    fn clone(&self) -> Self { Self { laws: self.laws, points: self.points.clone() } }
}

/// One accepted native interval and its corresponding initiation history.
#[derive(Debug)]
pub struct DuctileAccepted<'a, M, L: DuctileInitiationLaw> {
    /// Original native result, including reactions and existing energy accounts.
    pub mechanics: M,
    /// Initiation history at exactly that native endpoint.
    pub ductility: DuctileHistory<'a, L>,
}

/// A ductility refusal never publishes any part of a proposed increment.
#[derive(Clone, Debug, PartialEq)]
pub enum DuctileSolidError {
    /// Shape, clock, calibration or old-history mismatch.
    Input(&'static str),
    /// Original calibrated material criterion refusal at one point.
    Point { point: usize, source: DuctileError },
    /// Earliest held-endpoint initiation estimate, with point-order tie breaking.
    /// The fraction is a fraction of PLASTIC accumulation, not resolved time.
    /// Refine the mechanical path to localize the actual initiation event.
    Initiated { point: usize, candidate: DuctileStep },
    /// Cancellation before the complete mechanics/history pair was returned.
    Cancelled,
}
impl core::fmt::Display for DuctileSolidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ductile-initiation gate refused: {self:?}")
    }
}
impl std::error::Error for DuctileSolidError {}

/// Preserve native mechanical errors separately from the initiation gate.
#[derive(Debug)]
pub enum DuctileAdvanceError<E> {
    /// Original geometry, material, equilibrium, energy or cancellation error.
    Mechanics(E),
    /// Original state and calibration remain unchanged on this refusal.
    Ductility(DuctileSolidError),
}
impl<E> From<DuctileSolidError> for DuctileAdvanceError<E> {
    fn from(value: DuctileSolidError) -> Self { Self::Ductility(value) }
}
impl<E: core::fmt::Display> core::fmt::Display for DuctileAdvanceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self { Self::Mechanics(e) => e.fmt(f), Self::Ductility(e) => e.fmt(f) }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DuctileAdvanceError<E> {}

fn checkpoint(cx: &Cx<'_>) -> Result<(), DuctileSolidError> {
    cx.checkpoint().map_err(|_| DuctileSolidError::Cancelled)
}

impl<'a, L: DuctileInitiationLaw> DuctileHistory<'a, L> {
    /// Allocate bounded virgin history. Old mechanical alpha must be zero;
    /// resume strained specimens with their retained returned history.
    pub fn virgin(
        laws: DuctileLawField<'a, L>, point_count: usize, maximum_points: usize,
    ) -> Result<Self, DuctileSolidError> {
        if point_count == 0 || point_count > maximum_points {
            return Err(DuctileSolidError::Input("nonempty point count within the caller budget required"));
        }
        if let DuctileLawField::PerPoint(entries) = laws {
            if entries.len() != point_count {
                return Err(DuctileSolidError::Input("calibration field must cover every integration point"));
            }
        }
        Ok(Self { laws, points: vec![DuctileState::virgin(); point_count] })
    }

    /// Original accepted initiation states in native point order.
    pub fn points(&self) -> &[DuctileState] { &self.points }
    /// Immutable original calibration assignment.
    pub fn laws(&self) -> DuctileLawField<'a, L> { self.laws }

    fn prepare<S: PlasticStrainHistory>(
        &self, states: &[S], conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<(), DuctileSolidError> {
        checkpoint(cx)?;
        if states.len() != self.points.len() {
            return Err(DuctileSolidError::Input("native and ductile point counts differ"));
        }
        if !conditions.duration_s.is_finite() || conditions.duration_s <= 0.0 {
            return Err(DuctileSolidError::Input("finite positive physical duration required"));
        }
        let temperatures_cover_points = match conditions.temperatures_k {
            DuctileTemperatures::Uniform(_) => true,
            DuctileTemperatures::PerPoint(temperatures) => temperatures.len() == states.len(),
            DuctileTemperatures::PerHexElement(temperatures) =>
                temperatures.len().checked_mul(8) == Some(states.len()),
        };
        if !temperatures_cover_points {
            return Err(DuctileSolidError::Input("temperature field must cover every integration point"));
        }
        for (point, (native, history)) in states.iter().zip(&self.points).enumerate() {
            if point % 256 == 0 { checkpoint(cx)?; }
            if native.plastic_strain_accumulation() != history.accumulated_plastic_strain() {
                return Err(DuctileSolidError::Input("native plastic history differs from retained ductility; no history reset or skipped increments allowed"));
            }
            let temperature = conditions.temperatures_k.at(point);
            let domain = self.laws.at(point).domain().temperature_k;
            if !temperature.is_finite() || temperature <= 0.0
                || temperature < domain[0] || temperature > domain[1]
            {
                return Err(DuctileSolidError::Input("actual point temperature lies outside the ductility calibration"));
            }
        }
        Ok(())
    }

    fn trial<S: PlasticStrainHistory>(
        &self, states: &[S], stress: &[Voigt],
        conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<Self, DuctileSolidError> {
        if states.len() != self.points.len() || stress.len() != states.len() {
            return Err(DuctileSolidError::Input("proposed native point arrays do not match retained history"));
        }
        self.trial_stress(states, |point| Ok(stress[point]), conditions, cx)
    }

    fn trial_stress<S: PlasticStrainHistory>(
        &self, states: &[S], stress: impl Fn(usize) -> Result<Voigt, DuctileSolidError>,
        conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<Self, DuctileSolidError> {
        if states.len() != self.points.len() {
            return Err(DuctileSolidError::Input("proposed native point count differs from retained history"));
        }
        let mut points = Vec::with_capacity(self.points.len());
        let mut initiation: Option<(usize, DuctileStep)> = None;
        for (point, (native, history)) in states.iter().zip(&self.points).enumerate() {
            if point % 256 == 0 { checkpoint(cx)?; }
            let candidate = history.trial(
                self.laws.at(point), native.plastic_strain_accumulation(),
                conditions.duration_s, conditions.temperatures_k.at(point), &stress(point)?,
            ).map_err(|source| DuctileSolidError::Point { point, source })?;
            points.push(candidate.state);
            if candidate.state.initiated() && initiation.as_ref().is_none_or(|(_, prior)| {
                candidate.initiation_fraction.unwrap_or(0.0) < prior.initiation_fraction.unwrap_or(0.0)
            }) {
                initiation = Some((point, candidate));
            }
        }
        checkpoint(cx)?;
        if let Some((point, candidate)) = initiation {
            return Err(DuctileSolidError::Initiated { point, candidate });
        }
        Ok(Self { laws: self.laws, points })
    }

    /// One original P1 tetrahedral quasistatic solve, gating every element.
    /// The problem must contain exactly one absolute physical load target.
    pub fn advance_tet<M: SmallStrainLaw>(
        &self, problem: &SmallStrainTetProblem<'_, M>,
        displacement_m: &[f64], states: &[M::State],
        conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, SmallStrainTetSolution<M::State>, L>, DuctileAdvanceError<SmallStrainTetError>>
    where M::State: PlasticStrainHistory {
        if problem.load_path.len() != 1 {
            return Err(DuctileSolidError::Input("one physical load target per ductility update required").into());
        }
        self.prepare(states, conditions, cx)?;
        let mechanics = problem.solve_from(displacement_m, states, cx).map_err(DuctileAdvanceError::Mechanics)?;
        let ductility = self.trial(&mechanics.element_states, &mechanics.element_stress, conditions, cx)?;
        Ok(DuctileAccepted { mechanics, ductility })
    }

    /// Original tetrahedral backward-Euler motion. One duration drives both
    /// inertia and the fracture calibration's plastic strain rate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_tet_backward_euler<M: SmallStrainLaw>(
        &self, problem: &SmallStrainTetProblem<'_, M>, displacement_m: &[f64],
        velocity_mps: &[[f64; 3]], states: &[M::State], reference_density_kg_m3: &[f64],
        conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, BackwardEulerStep<M::State>, L>, DuctileAdvanceError<SmallStrainTetError>>
    where M::State: PlasticStrainHistory {
        self.prepare(states, conditions, cx)?;
        let mechanics = problem.advance_backward_euler(displacement_m, velocity_mps, states,
            reference_density_kg_m3, conditions.duration_s, cx).map_err(DuctileAdvanceError::Mechanics)?;
        let ductility = self.trial(&mechanics.solution.element_states,
            &mechanics.solution.element_stress, conditions, cx)?;
        Ok(DuctileAccepted { mechanics, ductility })
    }

    /// Original plane-stress/plane-strain P1/Q1 increment using full physical
    /// three-dimensional stresses and separate native quadrature histories.
    pub fn advance_plane<M: SmallStrainLaw>(
        &self, problem: &SmallStrainPlaneProblem<'_, M>, displacement_m: &[f64],
        states: &[ReducedState<M::State>], conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, SmallStrainPlaneSolution<M::State>, L>, DuctileAdvanceError<SmallStrainPlaneError>>
    where M::State: PlasticStrainHistory {
        if problem.load_path.len() != 1 {
            return Err(DuctileSolidError::Input("one physical load target per ductility update required").into());
        }
        self.prepare(states, conditions, cx)?;
        let mechanics = problem.solve_from(displacement_m, states, cx).map_err(DuctileAdvanceError::Mechanics)?;
        let ductility = self.trial(&mechanics.quadrature_states,
            &mechanics.quadrature_stress_pa, conditions, cx)?;
        Ok(DuctileAccepted { mechanics, ductility })
    }

    /// Original fully integrated/B-bar hex increment. All eight Gauss-point
    /// histories are retained; ductility is never element-averaged.
    pub fn advance_hex<M: SmallStrainLaw>(
        &self, problem: &SmallStrainHexProblem<'_, M>, displacement_m: &[f64],
        states: &[M::State], conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, SmallStrainHexSolution<M::State>, L>, DuctileAdvanceError<SmallStrainTetError>>
    where M::State: PlasticStrainHistory {
        if problem.load_path.len() != 1 {
            return Err(DuctileSolidError::Input("one physical load target per ductility update required").into());
        }
        self.prepare(states, conditions, cx)?;
        let mechanics = problem.solve_from(displacement_m, states, cx).map_err(DuctileAdvanceError::Mechanics)?;
        let ductility = self.trial(&mechanics.quadrature_states,
            &mechanics.quadrature_stress_pa, conditions, cx)?;
        Ok(DuctileAccepted { mechanics, ductility })
    }
}

impl<'a, L: DuctileMaterialBinding> DuctileHistory<'a, L> {
    fn check_materials(
        &self, materials: &[fs_material::plastic::source::ResolvedJ2Plasticity], cx: &Cx<'_>,
    ) -> Result<(), DuctileSolidError> {
        checkpoint(cx)?;
        if materials.len() != self.points.len() {
            return Err(DuctileSolidError::Input("material sources must cover every ductile point"));
        }
        for (point, material) in materials.iter().enumerate() {
            if point % 256 == 0 { checkpoint(cx)?; }
            if !self.laws.at(point).matches_material(material.material()) {
                return Err(DuctileSolidError::Input(
                    "fracture calibration does not match native material card, state and context",
                ));
            }
        }
        checkpoint(cx)
    }

    /// Admit accepted material/history/temperature before an external coupled
    /// owner performs an expensive proposal. This reads the native history;
    /// it never re-runs a constitutive update or changes initiation history.
    pub fn validate_material_history<S: PlasticStrainHistory>(
        &self, states: &[S], materials: &[fs_material::plastic::source::ResolvedJ2Plasticity],
        conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<(), DuctileSolidError> {
        self.check_materials(materials, cx)?;
        self.prepare(states, conditions, cx)
    }

    /// Gate one already converged native tetrahedral interval, preserving the
    /// complete source binding at both endpoints. Thermal/coupled owners use
    /// their actual endpoint temperatures and physical duration; neither a
    /// second Newton solve nor a second constitutive increment is performed.
    /// The caller retains the complete coupled candidate until this succeeds.
    #[allow(clippy::too_many_arguments)]
    pub fn check_material_tet_candidate<S: PlasticStrainHistory>(
        &self, old_states: &[S],
        old_materials: &[fs_material::plastic::source::ResolvedJ2Plasticity],
        candidate: &SmallStrainTetSolution<S>,
        candidate_materials: &[fs_material::plastic::source::ResolvedJ2Plasticity],
        conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<Self, DuctileSolidError> {
        self.validate_material_history(old_states, old_materials, conditions, cx)?;
        self.check_materials(candidate_materials, cx)?;
        if candidate.increments.len() != 1 || candidate.element_strain.len() != self.points.len() {
            return Err(DuctileSolidError::Input("one complete native tetrahedral interval required"));
        }
        self.trial(&candidate.element_states, &candidate.element_stress, conditions, cx)
    }
}

impl<'a> DuctileHistory<'a, fs_material::ductile::source::ResolvedDuctileLocus> {
    fn check_j2_sources(
        &self,
        problem: &SmallStrainTetProblem<'_, fs_material::plastic::source::ResolvedJ2Plasticity>,
        cx: &Cx<'_>,
    ) -> Result<(), DuctileSolidError> {
        checkpoint(cx)?;
        if problem.tetrahedra.len() != self.points.len() {
            return Err(DuctileSolidError::Input("source-bound native and ductile point counts differ"));
        }
        if let SmallStrainLawField::PerElement(native) = problem.laws {
            if native.len() != self.points.len() {
                return Err(DuctileSolidError::Input("native J2 field does not cover every tetrahedron"));
            }
        }
        for point in 0..self.points.len() {
            if point % 256 == 0 { checkpoint(cx)?; }
            let native = match problem.laws {
                SmallStrainLawField::Uniform(law) => law,
                SmallStrainLawField::PerElement(entries) => &entries[point],
            };
            if !self.laws.at(point).source().matches_material(native.material()) {
                return Err(DuctileSolidError::Input("fracture calibration does not match native material card, state and context"));
            }
        }
        checkpoint(cx)
    }

    /// Bind each fracture calibration to the original card, named state and
    /// fixed external context of its native J2 tetrahedron.
    pub fn from_source_bound_j2_tet(
        problem: &SmallStrainTetProblem<'_, fs_material::plastic::source::ResolvedJ2Plasticity>,
        laws: DuctileLawField<'a, fs_material::ductile::source::ResolvedDuctileLocus>,
        cx: &Cx<'_>,
    ) -> Result<Self, DuctileSolidError> {
        checkpoint(cx)?;
        let history = Self::virgin(laws, problem.tetrahedra.len(), problem.budget.maximum_tetrahedra)?;
        history.check_j2_sources(problem, cx)?;
        Ok(history)
    }

    /// Source-bound native J2 continuation, rechecking original material/card
    /// association on EVERY increment so a later grade or context swap refuses.
    pub fn advance_source_bound_j2_tet(
        &self,
        problem: &SmallStrainTetProblem<'_, fs_material::plastic::source::ResolvedJ2Plasticity>,
        displacement_m: &[f64], states: &[J2State],
        conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, SmallStrainTetSolution<J2State>,
        fs_material::ductile::source::ResolvedDuctileLocus>, DuctileAdvanceError<SmallStrainTetError>> {
        self.check_j2_sources(problem, cx)?;
        self.advance_tet(problem, displacement_m, states, conditions, cx)
    }
}
