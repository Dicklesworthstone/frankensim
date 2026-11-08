//! Ductile initiation on the original CST/DKT shell material points.
//!
//! Every facet, DKT surface point and thickness point retains its own criterion,
//! in that order. Membrane and bending resultants cannot replace local stress:
//! tensile and compressive layers can yield while their membrane forces cancel.
//! The original source-bound J2, Voce or Chaboche law recovers full physical
//! stress from accepted elastic strain through its read-only observation seam.
//! No second return map advances history, and backstress is never subtracted.
//!
//! Mechanical sources and their explicitly typed absolute temperatures remain
//! fixed. Each supplied point temperature must match its original source query;
//! changing-temperature, creep and viscous laws need their own stress owner.
//! The same physical duration drives inertia and the initiation rate. All
//! points gate one native increment atomically, including during reversal.
//! A refusal leaves the borrowed shell, velocity and ductility unchanged.
//! This supplies initiation only, without degraded stiffness, fracture energy,
//! crack topology or a finite-strain forming model.

use fs_material::fatigue::cyclic::SourceBoundFatigueLaw;
use fs_material::plastic::section::LayeredPlaneStressSection;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::*;
use crate::plastic_shell::dynamics::{PlasticShellDynamicStep, PlasticShellReferenceInertia};
use crate::plastic_shell::{
    PlasticShellSections, SmallStrainShellProblem, SmallStrainShellSolution,
};

fn section_at<'m, M: SmallStrainLaw>(
    sections: &PlasticShellSections<'m, M>,
    element: usize,
) -> &'m LayeredPlaneStressSection<M> {
    match sections {
        PlasticShellSections::Uniform(section) => section,
        PlasticShellSections::PerElement(sections) => &sections[element],
        PlasticShellSections::PerElementBorrowed(sections) => sections[element],
    }
}

// Admit the complete nested layout before indexing or allocating a point bank.
fn point_count<M: SmallStrainLaw>(
    problem: &SmallStrainShellProblem<'_, M>,
    cx: &Cx<'_>,
) -> Result<usize, DuctileSolidError> {
    checkpoint(cx)?;
    let elements = problem.mesh.tris.len();
    let sections = match &problem.sections {
        PlasticShellSections::Uniform(_) => elements,
        PlasticShellSections::PerElement(sections) => sections.len(),
        PlasticShellSections::PerElementBorrowed(sections) => sections.len(),
    };
    if elements == 0 || elements > problem.budget.maximum_elements || sections != elements {
        return Err(DuctileSolidError::Input(
            "invalid shell section assignment or element budget",
        ));
    }
    let mut count = 0usize;
    for element in 0..elements {
        if element % 256 == 0 {
            checkpoint(cx)?;
        }
        count = section_at(&problem.sections, element)
            .point_count()
            .checked_mul(3)
            .and_then(|n| count.checked_add(n))
            .ok_or(DuctileSolidError::Input(
                "shell material point count overflow",
            ))?;
        if count > problem.budget.maximum_material_points {
            return Err(DuctileSolidError::Input(
                "shell ductility exceeds the material point budget",
            ));
        }
    }
    Ok(count)
}

fn source_temperature_matches(material: &ResolvedMaterialStatePoint, temperature_k: f64) -> bool {
    let mut found = false;
    for (axis, quantity) in material.axis_quantities() {
        let absolute = *quantity == QuantitySpec::dimensional(Temperature::DIMS)
            || quantity.semantic_type().is_some_and(|semantic| {
                semantic.kind() == QuantityKind::AbsoluteTemperature
                    && semantic.form() == ValueForm::Static
            });
        if absolute {
            found = true;
            if material.state_coordinate(axis) != Some(temperature_k) {
                return false;
            }
        }
    }
    found
}

fn admit_binding<L: DuctileMaterialBinding>(
    law: &L,
    material: &ResolvedMaterialStatePoint,
    temperature: f64,
) -> Result<(), DuctileSolidError> {
    let domain = law.domain().temperature_k;
    if !temperature.is_finite()
        || temperature <= 0.0
        || temperature < domain[0]
        || temperature > domain[1]
        || !source_temperature_matches(material, temperature)
    {
        return Err(DuctileSolidError::Input(
            "actual shell temperature differs from the fixed mechanical source or ductility calibration",
        ));
    }
    if !law.matches_material(material) {
        return Err(DuctileSolidError::Input(
            "shell ductility calibration differs from the original material card, state or fixed context",
        ));
    }
    Ok(())
}

impl<'a, L: DuctileMaterialBinding> DuctileHistory<'a, L> {
    /// Allocate virgin initiation history in the original facet/surface/thickness
    /// order, using the native shell's material-point budget. All three surface
    /// points and every layer quadrature point are included automatically.
    pub fn virgin_shell<M: SourceBoundFatigueLaw>(
        laws: DuctileLawField<'a, L>,
        problem: &SmallStrainShellProblem<'_, M>,
        cx: &Cx<'_>,
    ) -> Result<Self, DuctileSolidError> {
        Self::virgin(
            laws,
            point_count(problem, cx)?,
            problem.budget.maximum_material_points,
        )
    }

    fn prepare_shell<M: SourceBoundFatigueLaw>(
        &self,
        problem: &SmallStrainShellProblem<'_, M>,
        initial: Option<&SmallStrainShellSolution<M::State>>,
        conditions: DuctileStepConditions<'_>,
        cx: &Cx<'_>,
    ) -> Result<(), DuctileSolidError>
    where
        M::State: PlasticStrainHistory,
    {
        if point_count(problem, cx)? != self.points.len()
            || problem.load_path.len() != 1
            || !conditions.duration_s.is_finite()
            || conditions.duration_s <= 0.0
            || initial.is_some_and(|old| {
                old.coordinates.len() != problem.mesh.nodes.len()
                    || old.section_states.len() != 3 * problem.mesh.tris.len()
            })
        {
            return Err(DuctileSolidError::Input(
                "shell ductility needs matching point history and one positive-duration physical increment",
            ));
        }
        match conditions.temperatures_k {
            DuctileTemperatures::Uniform(_) => {}
            DuctileTemperatures::PerPoint(values) if values.len() == self.points.len() => {}
            _ => {
                return Err(DuctileSolidError::Input(
                    "shell temperature must be uniform or cover every facet/surface/thickness point",
                ));
            }
        }
        let uniform = match (self.laws, conditions.temperatures_k) {
            (DuctileLawField::Uniform(law), DuctileTemperatures::Uniform(t)) => Some((law, t)),
            _ => None,
        };
        let mut index = 0;
        for element in 0..problem.mesh.tris.len() {
            let section = section_at(&problem.sections, element);
            // The common uniform field binds each physical layer once instead
            // of reconstructing its source query at every thickness/surface
            // point. Heterogeneous fields still bind each original point.
            if let Some((law, temperature)) = uniform {
                for (layer, native) in section.layers().iter().enumerate() {
                    if layer % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    admit_binding(law, native.law.fatigue_source(), temperature)?;
                }
            }
            for surface in 0..3 {
                let old = initial.map(|old| &old.section_states[3 * element + surface]);
                if old.is_some_and(|state| state.points.len() != section.point_count()) {
                    return Err(DuctileSolidError::Input(
                        "shell thickness history does not match the original quadrature",
                    ));
                }
                for (thickness, point) in section.points().iter().enumerate() {
                    if index % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    let alpha = old.map_or(0.0, |state| {
                        state.points[thickness].plastic_strain_accumulation()
                    });
                    if alpha != self.points[index].accumulated_plastic_strain() {
                        return Err(DuctileSolidError::Input(
                            "shell plastic history differs from retained ductility; no reset or skipped increments allowed",
                        ));
                    }
                    if uniform.is_none() {
                        admit_binding(
                            self.laws.at(index),
                            section.layers()[point.layer].law.fatigue_source(),
                            conditions.temperatures_k.at(index),
                        )?;
                    }
                    index += 1;
                }
            }
        }
        checkpoint(cx)
    }

    fn trial_shell<M: SourceBoundFatigueLaw>(
        &self,
        problem: &SmallStrainShellProblem<'_, M>,
        solution: &SmallStrainShellSolution<M::State>,
        conditions: DuctileStepConditions<'_>,
        cx: &Cx<'_>,
    ) -> Result<Self, DuctileAdvanceError<SmallStrainTetError>>
    where
        M::State: PlasticStrainHistory,
    {
        let mut points = Vec::with_capacity(self.points.len());
        let mut initiation: Option<(usize, DuctileStep)> = None;
        for (surface, state) in solution.section_states.iter().enumerate() {
            let section = section_at(&problem.sections, surface / 3);
            for (thickness, native) in state.points.iter().enumerate() {
                let index = points.len();
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                let layer = section.points()[thickness].layer;
                let stress = section.layers()[layer].law
                    .accepted_fatigue_sample(&native.strain, &native.history)
                    .map_err(|error| DuctileAdvanceError::Mechanics(SmallStrainTetError::InvalidInput {
                        what: format!("shell ductility element {}, surface {}, thickness {thickness}: {error}", surface / 3, surface % 3),
                    }))?.stress_pa;
                let candidate = self.points[index]
                    .trial(
                        self.laws.at(index),
                        native.plastic_strain_accumulation(),
                        conditions.duration_s,
                        conditions.temperatures_k.at(index),
                        &stress,
                    )
                    .map_err(|source| DuctileSolidError::Point {
                        point: index,
                        source,
                    })?;
                points.push(candidate.state);
                if candidate.state.initiated()
                    && initiation.as_ref().is_none_or(|(_, prior)| {
                        candidate.initiation_fraction.unwrap_or(0.0)
                            < prior.initiation_fraction.unwrap_or(0.0)
                    })
                {
                    initiation = Some((index, candidate));
                }
            }
        }
        checkpoint(cx)?;
        if let Some((point, candidate)) = initiation {
            return Err(DuctileSolidError::Initiated { point, candidate }.into());
        }
        Ok(Self {
            laws: self.laws,
            points,
        })
    }

    /// Propose one original quasistatic shell increment. `None` starts from
    /// virgin zero coordinates; subsequent increments require the complete
    /// previously accepted shell solution and its matching initiation history.
    /// Source binding covers every layer, including unyielded layers. Physical
    /// point temperatures use the same flattened order as [`Self::points`].
    pub fn advance_shell<M: SourceBoundFatigueLaw>(
        &self,
        problem: &SmallStrainShellProblem<'_, M>,
        initial: Option<&SmallStrainShellSolution<M::State>>,
        conditions: DuctileStepConditions<'_>,
        cx: &Cx<'_>,
    ) -> Result<
        DuctileAccepted<'a, SmallStrainShellSolution<M::State>, L>,
        DuctileAdvanceError<SmallStrainTetError>,
    >
    where
        M::State: PlasticStrainHistory,
    {
        self.prepare_shell(problem, initial, conditions, cx)?;
        let mechanics = match initial {
            Some(old) => problem.solve_from(&old.coordinates, &old.section_states, cx),
            None => problem.solve(cx),
        }
        .map_err(DuctileAdvanceError::Mechanics)?;
        let ductility = self.trial_shell(problem, &mechanics, conditions, cx)?;
        Ok(DuctileAccepted {
            mechanics,
            ductility,
        })
    }

    /// Propose one original backward-Euler shell step with retained reference
    /// inertia. Native momentum and the ductility rate use this same duration;
    /// physical and numerical kinetic-energy accounts remain with mechanics.
    pub fn advance_shell_backward_euler<M: SourceBoundFatigueLaw>(
        &self,
        problem: &SmallStrainShellProblem<'_, M>,
        initial: Option<&SmallStrainShellSolution<M::State>>,
        velocity: &[[f64; 6]],
        reference: &PlasticShellReferenceInertia<'_>,
        conditions: DuctileStepConditions<'_>,
        cx: &Cx<'_>,
    ) -> Result<
        DuctileAccepted<'a, PlasticShellDynamicStep<M::State>, L>,
        DuctileAdvanceError<SmallStrainTetError>,
    >
    where
        M::State: PlasticStrainHistory,
    {
        self.prepare_shell(problem, initial, conditions, cx)?;
        let mechanics = problem
            .advance_backward_euler(initial, velocity, reference, conditions.duration_s, cx)
            .map_err(DuctileAdvanceError::Mechanics)?;
        let ductility = self.trial_shell(problem, &mechanics.solution, conditions, cx)?;
        Ok(DuctileAccepted {
            mechanics,
            ductility,
        })
    }
}
