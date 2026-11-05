//! Elastic prehistory, tensile-plane prediction and conservative brittle handoff.
//!
//! The original per-element elastic backbones own every accepted elastic solve.
//! A proposed affine load increment is accepted as elastic only when the native
//! tensile detector finds no crossing. Otherwise the prediction is disposable:
//! source-bound events select planes, reference geometry supplies V/A widths,
//! and handoff initializes native crack bands at the OLD accepted strain.
//! Existing elastic storage is transferred, never charged as new external work.
//! The caller replays the requested endpoint with the retained brittle solver.
//!
//! This is an elastic-predictor, fixed-plane model. Predicted crossings after
//! the first event do not account for fracture-induced stress redistribution;
//! their brackets are not accepted fracture times. Elements without a predicted
//! crossing retain their explicitly supplied source plane. Resolve physical
//! increments and any plane sensitivity. No rotating cracks, all-plane mixed-
//! mode detection, topology change, contact or dynamic fracture is supplied.

/// Native replay toward the first tensile event and an explicit load remainder.
pub mod first_event;

use fs_exec::Cx;
use fs_material::brittle::initiation::TensileInitiationSettings;
use fs_material::brittle::source::ResolvedBrittleCrackBand;
use fs_material::brittle::source::initiation::{
    ResolvedTensileInitiationError, ResolvedTensileOnset,
};
use fs_material::tensor::contract;

use super::{BrittleTetError, BrittleTetSolid, checkpoint, energy_gate};
use crate::brittle::TetCrackBandGeometry;
use crate::linear3::TetAssemblyBudget;
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    SmallStrainTetSolution,
};

/// A retained history that has stayed below the original tensile criterion.
///
/// It is constructed virgin and can advance only through checked elastic
/// endpoints. The original source bundles remain immutable. This history is
/// tensile-only: it is not permission to erase prior shear/plastic damage.
#[derive(Debug)]
pub struct TensileElasticTetSolid<'a> {
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedBrittleCrackBand],
    budget: TetAssemblyBudget,
    solution: SmallStrainTetSolution<()>,
    stored_energy_j: f64,
    supports: Option<Vec<usize>>,
    nodal_forces_n: Vec<(usize, f64)>,
}

/// A checked elastic endpoint or a disposable onset prediction.
#[derive(Debug)]
pub enum TensileTetProbe<'a> {
    Elastic(TensileElasticTetSolid<'a>),
    Onset(TensileTetOnset<'a>),
}

/// The original elastic state and a proposed elastic endpoint with crossings.
///
/// Original source owners are borrowed. Receiving bands are resolved separately
/// so the eventual brittle histories can borrow them without self references.
#[derive(Debug)]
pub struct TensileTetOnset<'a> {
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedBrittleCrackBand],
    budget: TetAssemblyBudget,
    initial: SmallStrainTetSolution<()>,
    predictor: SmallStrainTetSolution<()>,
    elastic_storage_j: f64,
    events: Vec<Option<ResolvedTensileOnset<'a>>>,
    first_fraction: [f64; 2],
    first_candidate: usize,
    settings: TensileInitiationSettings,
    initial_forces_n: Vec<(usize, f64)>,
    prescribed_end_m: Vec<(usize, f64)>,
    nodal_forces_end_n: Vec<(usize, f64)>,
    solver_settings: IncrementSettings,
    direction_hints: Option<Vec<Option<[f64; 3]>>>,
}

/// A change of constitutive representation at the original accepted endpoint.
#[derive(Debug)]
pub struct BrittleTetHandoff<'a> {
    pub solid: BrittleTetSolid<'a>,
    pub elastic_storage_j: f64,
    pub storage_defect_j: f64,
    pub maximum_stress_defect_pa: f64,
}

/// Native solve, source/onset admission, geometry or continuity refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum TensileTetInitiationError {
    Invalid(&'static str),
    Element {
        element: usize,
        source: ResolvedTensileInitiationError,
    },
    Solid(BrittleTetError),
}

impl core::fmt::Display for TensileTetInitiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "elastic-to-brittle solid handoff refused: {self:?}")
    }
}
impl std::error::Error for TensileTetInitiationError {}
impl From<BrittleTetError> for TensileTetInitiationError {
    fn from(error: BrittleTetError) -> Self {
        Self::Solid(error)
    }
}
impl From<SmallStrainTetError> for TensileTetInitiationError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Solid(error.into())
    }
}

impl<'a> TensileElasticTetSolid<'a> {
    /// Retain original geometry and one resolved source per element, virgin.
    ///
    /// Source planes are explicit dormant modeling inputs, not direction hints.
    /// Their eventual widths are derived from this actual geometry at handoff.
    pub fn try_new(
        geometry: &'a TetLoadGeometry3<'a>,
        materials: &'a [ResolvedBrittleCrackBand],
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, TensileTetInitiationError> {
        checkpoint(cx)?;
        if materials.len() != geometry.tetrahedra().len() {
            return Err(BrittleTetError::MaterialCount {
                expected: geometry.tetrahedra().len(),
                got: materials.len(),
            }
            .into());
        }
        if materials.is_empty()
            || geometry.nodes_m().is_empty()
            || geometry.nodes_m().len() > budget.maximum_nodes
            || materials.len() > budget.maximum_tetrahedra
            || geometry.nodes_m().len().checked_mul(3).is_none()
            || !budget.minimum_scaled_jacobian.is_finite()
            || budget.minimum_scaled_jacobian <= 0.0
        {
            return Err(TensileTetInitiationError::Invalid(
                "nonempty original geometry must fit the native assembly budget",
            ));
        }
        let solution = SmallStrainTetSolution {
            displacement_m: vec![[0.0; 3]; geometry.nodes_m().len()],
            element_strain: vec![[0.0; 6]; materials.len()],
            element_stress: vec![[0.0; 6]; materials.len()],
            element_states: vec![(); materials.len()],
            increments: Vec::new(),
        };
        checkpoint(cx)?;
        Ok(Self {
            geometry,
            materials,
            budget,
            solution,
            stored_energy_j: 0.0,
            supports: None,
            nodal_forces_n: Vec::new(),
        })
    }

    #[must_use]
    pub const fn geometry(&self) -> &'a TetLoadGeometry3<'a> {
        self.geometry
    }
    #[must_use]
    pub const fn materials(&self) -> &'a [ResolvedBrittleCrackBand] {
        self.materials
    }
    #[must_use]
    pub const fn solution(&self) -> &SmallStrainTetSolution<()> {
        &self.solution
    }
    #[must_use]
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }

    /// Applied dead-force entries at the accepted elastic endpoint, including
    /// forces on supports. These are original inputs, not inferred reactions.
    #[must_use]
    pub fn nodal_forces_n(&self) -> &[(usize, f64)] {
        &self.nodal_forces_n
    }

    /// Solve a proposed native elastic endpoint and inspect every actual strain.
    ///
    /// Prescribed and force values are absolute targets. The prescribed DOF set
    /// stays fixed after the first endpoint, so affine loads/support motion give
    /// an affine elastic strain path. All reversals require separate endpoints.
    /// Direction hints are in element order and resolve only repeated principal
    /// tension; a source plane is never silently substituted for a missing hint.
    /// No crossing makes the returned elastic endpoint admissible. Any crossing
    /// returns a disposable predictor and leaves this accepted state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn probe(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        solver_settings: IncrementSettings,
        initiation_settings: TensileInitiationSettings,
        direction_hints: Option<&[Option<[f64; 3]>]>,
        cx: &Cx<'_>,
    ) -> Result<TensileTetProbe<'a>, TensileTetInitiationError> {
        checkpoint(cx)?;
        let dofs = self.solution.displacement_m.len() * 3;
        if prescribed_m.len() > dofs
            || nodal_forces_n.len() > dofs
            || direction_hints.is_some_and(|h| h.len() != self.materials.len())
        {
            return Err(TensileTetInitiationError::Invalid(
                "endpoint arrays must fit the original DOF and element counts",
            ));
        }
        let mut supports: Vec<_> = prescribed_m.iter().map(|&(dof, _)| dof).collect();
        supports.sort_unstable();
        if supports.iter().any(|&dof| dof >= dofs)
            || supports.windows(2).any(|pair| pair[0] == pair[1])
            || self.supports.as_ref().is_some_and(|old| old != &supports)
        {
            return Err(TensileTetInitiationError::Invalid(
                "elastic history requires the same unique prescribed DOF set",
            ));
        }
        let mut laws = Vec::with_capacity(self.materials.len());
        for (element, source) in self.materials.iter().enumerate() {
            checkpoint(cx)?;
            source
                .tensile_initiator(initiation_settings)
                .map_err(|source| TensileTetInitiationError::Element { element, source })?;
            laws.push(source.law().elastic());
        }
        let problem = SmallStrainTetProblem {
            nodes_m: self.geometry.nodes_m(),
            tetrahedra: self.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            budget: self.budget,
            settings: solver_settings,
        };
        let solution = problem.solve_from(
            self.solution.displacement_m.as_flattened(),
            &self.solution.element_states,
            cx,
        )?;
        let mut events = Vec::with_capacity(self.materials.len());
        let mut first_fraction = [1.0_f64; 2];
        let mut first_candidate = None;
        let mut smallest_upper = f64::INFINITY;
        let mut stored_energy_j = 0.0;
        for element in 0..self.materials.len() {
            checkpoint(cx)?;
            let source: &'a ResolvedBrittleCrackBand = &self.materials[element];
            let detector = source
                .tensile_initiator(initiation_settings)
                .map_err(|source| TensileTetInitiationError::Element { element, source })?;
            let event = detector
                .between_controlled(
                    &self.solution.element_strain[element],
                    &solution.element_strain[element],
                    direction_hints.and_then(|h| h[element]),
                    &mut || cx.checkpoint().is_err(),
                )
                .map_err(|source| TensileTetInitiationError::Element { element, source })?;
            if let Some(onset) = &event {
                let bracket = onset.native().fraction_bracket();
                first_fraction[0] = first_fraction[0].min(bracket[0]);
                first_fraction[1] = first_fraction[1].min(bracket[1]);
                if bracket[1] < smallest_upper {
                    smallest_upper = bracket[1];
                    first_candidate = Some(element);
                }
            }
            let density = 0.5
                * contract(
                    &solution.element_stress[element],
                    &solution.element_strain[element],
                );
            if !density.is_finite() || density < 0.0 {
                return Err(TensileTetInitiationError::Invalid(
                    "native elastic storage must be finite and nonnegative",
                ));
            }
            stored_energy_j =
                self.geometry.element_volumes_m3()[element].mul_add(density, stored_energy_j);
            events.push(event);
        }
        if !stored_energy_j.is_finite() {
            return Err(TensileTetInitiationError::Invalid(
                "whole-solid elastic storage is unrepresentable",
            ));
        }
        checkpoint(cx)?;
        match first_candidate {
            None => Ok(TensileTetProbe::Elastic(Self {
                geometry: self.geometry,
                materials: self.materials,
                budget: self.budget,
                solution,
                stored_energy_j,
                supports: Some(supports),
                nodal_forces_n: nodal_forces_n.to_vec(),
            })),
            Some(first_candidate) => Ok(TensileTetProbe::Onset(TensileTetOnset {
                geometry: self.geometry,
                materials: self.materials,
                budget: self.budget,
                initial: self.solution.clone(),
                predictor: solution,
                elastic_storage_j: self.stored_energy_j,
                events,
                first_fraction,
                first_candidate,
                settings: initiation_settings,
                initial_forces_n: self.nodal_forces_n.clone(),
                prescribed_end_m: prescribed_m.to_vec(),
                nodal_forces_end_n: nodal_forces_n.to_vec(),
                solver_settings,
                direction_hints: direction_hints.map(|hints| hints.to_vec()),
            })),
        }
    }
}

impl<'a> TensileTetOnset<'a> {
    #[must_use]
    pub fn events(&self) -> &[Option<ResolvedTensileOnset<'a>>] {
        &self.events
    }
    #[must_use]
    pub const fn initial_solution(&self) -> &SmallStrainTetSolution<()> {
        &self.initial
    }
    /// Proposed original elastic response; not an accepted damaged endpoint.
    #[must_use]
    pub const fn elastic_predictor(&self) -> &SmallStrainTetSolution<()> {
        &self.predictor
    }
    /// Numerical bracket for the first crossing of the elastic predictor.
    #[must_use]
    pub const fn first_fraction_bracket(&self) -> [f64; 2] {
        self.first_fraction
    }
    /// First element attaining the smallest upper bracket. Overlapping brackets
    /// need not identify a unique exact event order.
    #[must_use]
    pub const fn first_candidate_element(&self) -> usize {
        self.first_candidate
    }

    /// Resolve actual projected widths on the original geometry.
    ///
    /// Predicted events use their discovered principal planes. Every other
    /// element retains its explicitly supplied plane, without inventing a
    /// direction from a nonunique or noncritical stress. These fixed planes
    /// are an inspectable modeling choice, not an adaptive crack-path solution.
    pub fn resolve_bands(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<ResolvedBrittleCrackBand>, TensileTetInitiationError> {
        checkpoint(cx)?;
        let mut bands = Vec::with_capacity(self.materials.len());
        for (element, source) in self.materials.iter().enumerate() {
            checkpoint(cx)?;
            let event = self.events[element].as_ref();
            let normal = event.map_or_else(|| source.law().normal(), |e| e.native().normal());
            let geometry = TetCrackBandGeometry::from_element(self.geometry, element, normal)
                .map_err(BrittleTetError::Geometry)?;
            let resolved = if let Some(event) = event {
                event.resolve_band(geometry.band_width_m())
            } else {
                let mut settings = source.law().settings();
                settings.band_width_m = geometry.band_width_m();
                ResolvedBrittleCrackBand::try_from_material(source.material().clone(), settings)
                    .map_err(ResolvedTensileInitiationError::Source)
            }
            .map_err(|source| TensileTetInitiationError::Element { element, source })?;
            bands.push(resolved);
        }
        checkpoint(cx)?;
        Ok(bands)
    }

    /// Preserve accepted elastic strain, stress and storage in receiving bands.
    ///
    /// This returns the OLD accepted endpoint under the brittle representation.
    /// Replay the requested load with BrittleTetSolid::advance or advance_path.
    /// No predictor work/damage is accepted. Every receiving source, geometric
    /// width, plane, undamaged seed and continuity gate is checked before return.
    pub fn handoff<'next>(
        &self,
        materials: &'next [ResolvedBrittleCrackBand],
        energy_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetHandoff<'next>, TensileTetInitiationError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        if !energy_tolerance_j.is_finite() || energy_tolerance_j < 0.0 {
            return Err(TensileTetInitiationError::Invalid(
                "handoff requires a finite nonnegative whole-solid energy tolerance",
            ));
        }
        let mut solid = BrittleTetSolid::try_new(self.geometry, materials, self.budget, cx)?;
        let mut states = Vec::with_capacity(materials.len());
        let mut stresses = Vec::with_capacity(materials.len());
        let mut storage = 0.0;
        let mut maximum_stress_defect_pa = 0.0_f64;
        for (element, receiving) in materials.iter().enumerate() {
            checkpoint(cx)?;
            let original = &self.materials[element];
            if receiving.material() != original.material()
                || receiving.law().elastic() != original.law().elastic()
                || receiving.law().parameters() != original.law().parameters()
            {
                return Err(TensileTetInitiationError::Element {
                    element,
                    source: ResolvedTensileInitiationError::SourceMismatch,
                });
            }
            if receiving.law().settings().maximum_strain != original.law().settings().maximum_strain
                || receiving.law().settings().energy_tolerance_j_m3
                    != original.law().settings().energy_tolerance_j_m3
            {
                return Err(TensileTetInitiationError::Element {
                    element,
                    source: ResolvedTensileInitiationError::ControlsMismatch,
                });
            }
            let strain = &self.initial.element_strain[element];
            let state = if let Some(event) = &self.events[element] {
                event
                    .seed(receiving)
                    .map_err(|source| TensileTetInitiationError::Element { element, source })?
                    .state
            } else {
                let a = original.law().normal();
                let b = receiving.law().normal();
                let sign = if a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>() < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                if (0..3).any(|i| (a[i] - sign * b[i]).abs() > 128.0 * f64::EPSILON) {
                    return Err(TensileTetInitiationError::Invalid(
                        "an untriggered element must retain its declared source plane",
                    ));
                }
                receiving
                    .law()
                    .trial(strain, &receiving.law().initial_state())
                    .map_err(|source| BrittleTetError::Material { element, source })?
                    .state
            };
            if state.damage() != 0.0
                || state.dissipated_j_m3() != 0.0
                || state.cohesive_point().mixed_mode_history().is_some()
            {
                return Err(TensileTetInitiationError::Invalid(
                    "handoff cannot erase prior mixed-mode damage",
                ));
            }
            let held = receiving
                .law()
                .trial(strain, &state)
                .map_err(|source| BrittleTetError::Material { element, source })?;
            if held.state != state || held.energy.work != 0.0 || held.energy.dissipated != 0.0 {
                return Err(TensileTetInitiationError::Invalid(
                    "seeded history must be stationary at its accepted strain",
                ));
            }
            let stress_defect = held
                .stress_pa
                .iter()
                .zip(self.initial.element_stress[element])
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            let stress_tolerance = self.settings.relative_stress_tolerance
                * original.law().parameters().tensile_strength_pa;
            if !stress_defect.is_finite() || stress_defect > stress_tolerance {
                return Err(TensileTetInitiationError::Invalid(
                    "receiving band changes the original elastic stress",
                ));
            }
            maximum_stress_defect_pa = maximum_stress_defect_pa.max(stress_defect);
            let old_density = 0.5 * contract(&self.initial.element_stress[element], strain);
            let density = state.stored_energy_j_m3();
            if !density.is_finite()
                || density < 0.0
                || (density - old_density).abs() > receiving.law().settings().energy_tolerance_j_m3
            {
                return Err(TensileTetInitiationError::Invalid(
                    "receiving band changes the original elastic storage",
                ));
            }
            storage = self.geometry.element_volumes_m3()[element].mul_add(density, storage);
            stresses.push(held.stress_pa);
            states.push(state);
        }
        let storage_defect_j = storage - self.elastic_storage_j;
        energy_gate(
            "elastic-to-brittle storage transfer",
            storage_defect_j,
            energy_tolerance_j,
        )?;
        solid.solution = SmallStrainTetSolution {
            displacement_m: self.initial.displacement_m.clone(),
            element_strain: self.initial.element_strain.clone(),
            element_stress: stresses,
            element_states: states,
            increments: Vec::new(),
        };
        solid.stored_energy_j = storage;
        checkpoint(cx)?;
        Ok(BrittleTetHandoff {
            solid,
            elastic_storage_j: self.elastic_storage_j,
            storage_defect_j,
            maximum_stress_defect_pa,
        })
    }
}
