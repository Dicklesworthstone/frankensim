//! Directional closed-loop circulation on the original thermal/pressure owners.
//!
//! A reverse flow has a different thermal trajectory: visit the original regions,
//! walls, resistances and gravity segments in reverse order, then solve again.
//! Retained thermal/hydraulic receipts use that actual traversal order. Signed
//! accessors express the result in the author's original direction without
//! negating heat or dissipation. The common inlet belongs to the closing reservoir.

use fs_exec::Cx;
use fs_qty::{MassFlowRate, Pressure, Temperature};

use super::super::{
    LiquidFeedbackBudget, LiquidHeatMarch, LiquidHeatPath, LiquidHeatSection, checkpoint,
};
use super::{
    LiquidResistanceSection, PressureDrivenBudget, PressureDrivenError, PressureDrivenSolution,
    SeriesGravity,
};

/// Direction relative to the authored section and vertex order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LiquidFlowDirection {
    Forward,
    Reverse,
}
impl LiquidFlowDirection {
    const fn sign(self) -> f64 {
        match self {
            Self::Forward => 1.0,
            Self::Reverse => -1.0,
        }
    }
}

/// How the caller selects a quasi-steady circulation branch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LiquidFlowSelection {
    /// Select one orientation and solve within the positive magnitude bracket.
    Direction(LiquidFlowDirection),
    /// Evaluate both orientation brackets and accept only if exactly one closes.
    /// This says nothing about roots outside those brackets, multiplicity inside
    /// a bracket, or dynamical stability. All source failures remain failures.
    OnlyBracketed,
}

/// Original positive-magnitude solve in ACTUAL flow order, with authored signs.
#[derive(Debug)]
pub struct DirectedLiquidSolution<'a> {
    direction: LiquidFlowDirection,
    operating_point: PressureDrivenSolution<LiquidHeatMarch<'a>>,
    total_flow_evaluations: usize,
}
impl<'a> DirectedLiquidSolution<'a> {
    pub const fn direction(&self) -> LiquidFlowDirection {
        self.direction
    }
    /// All source, region, heat and loss receipts in actual traversal order.
    pub const fn operating_point(&self) -> &PressureDrivenSolution<LiquidHeatMarch<'a>> {
        &self.operating_point
    }
    pub fn signed_mass_flow_kg_per_s(&self) -> f64 {
        self.direction.sign() * self.operating_point.mass_flow_kg_per_s
    }
    pub fn signed_driving_pressure_pa(&self) -> f64 {
        self.direction.sign() * self.operating_point.driving_pressure_pa
    }
    pub fn signed_pressure_drop_pa(&self) -> f64 {
        self.direction.sign() * self.operating_point.hydraulics.total_pressure_drop_pa
    }
    pub fn signed_hydrostatic_gain_pa(&self) -> f64 {
        self.direction.sign()
            * self
                .operating_point
                .hydrostatic
                .as_ref()
                .map_or(0.0, |h| h.pressure_gain_pa())
    }
    pub fn signed_pressure_residual_pa(&self) -> f64 {
        self.direction.sign() * self.operating_point.pressure_residual_pa
    }
    /// Includes the opposite orientation's actual endpoint/root evaluations.
    pub const fn total_flow_evaluations(&self) -> usize {
        self.total_flow_evaluations
    }
}

impl<'a, 'card> LiquidHeatPath<'a, 'card> {
    /// Solve either direction without changing the authored path or its source.
    ///
    /// Geometry must close at the same vertex, so either traversal takes supply
    /// from the SAME inlet reservoir. Bounds remain positive mass magnitudes;
    /// reverse traversal negates the imposed differential pressure, not gravity.
    /// Regional bulk material is recomputed on its new thermal trajectory.
    ///
    /// `budget.max_flow_evaluations` is shared across both orientations when
    /// selecting `OnlyBracketed`. Only an evaluated `NotBracketed` counts as an
    /// absent candidate; unsupported material/correlation states are never
    /// mistaken for an absent root. No bracket expansion, zero-flow heat model,
    /// inertial continuation or globally unique/stable branch is supplied.
    pub fn solve_gravity_circulation(
        &self,
        cx: &Cx<'_>,
        walls: &[f64],
        resistances: &[LiquidResistanceSection<'_>],
        drive: Pressure,
        gravity: &SeriesGravity,
        selection: LiquidFlowSelection,
        local: LiquidFeedbackBudget,
        budget: PressureDrivenBudget,
    ) -> Result<DirectedLiquidSolution<'a>, PressureDrivenError> {
        checkpoint(cx)?;
        budget.admit_with_gravity(drive)?;
        if gravity.section_count() != self.sections.len()
            || gravity.vertices_m().first() != gravity.vertices_m().last()
            || walls.len() != self.sections.len()
        {
            return Err(PressureDrivenError::Invalid(
                "circulation needs matching regions and coincident path endpoints",
            ));
        }
        self.admit_resistances(resistances)?;
        if let LiquidFlowSelection::Direction(direction) = selection {
            let operating_point = self.solve_circulation_direction(
                cx,
                walls,
                resistances,
                drive,
                gravity,
                direction,
                local,
                budget,
            )?;
            let total_flow_evaluations = operating_point.flow_evaluations;
            return Ok(DirectedLiquidSolution {
                direction,
                operating_point,
                total_flow_evaluations,
            });
        }
        let mut performed = 0;
        let mut found: Option<DirectedLiquidSolution<'a>> = None;
        let mut failures = [[0.0; 2]; 2];
        for (index, direction) in [LiquidFlowDirection::Forward, LiquidFlowDirection::Reverse]
            .into_iter()
            .enumerate()
        {
            checkpoint(cx)?;
            let remaining = budget.max_flow_evaluations.saturating_sub(performed);
            if remaining < 2 {
                return Err(PressureDrivenError::CirculationWorkBudget {
                    performed,
                    maximum: budget.max_flow_evaluations,
                });
            }
            match self.solve_circulation_direction(
                cx,
                walls,
                resistances,
                drive,
                gravity,
                direction,
                local,
                PressureDrivenBudget {
                    max_flow_evaluations: remaining,
                    ..budget
                },
            ) {
                Ok(operating_point) => {
                    performed += operating_point.flow_evaluations;
                    if let Some(previous) = &found {
                        return Err(PressureDrivenError::AmbiguousCirculationBrackets {
                            forward_mass_flow_kg_per_s: previous.signed_mass_flow_kg_per_s(),
                            reverse_mass_flow_kg_per_s: -operating_point.mass_flow_kg_per_s,
                        });
                    }
                    found = Some(DirectedLiquidSolution {
                        direction,
                        operating_point,
                        total_flow_evaluations: performed,
                    });
                }
                Err(PressureDrivenError::NotBracketed {
                    lower_residual_pa,
                    upper_residual_pa,
                }) => {
                    // The existing root evaluates both endpoints before this refusal.
                    performed += 2;
                    failures[index] = [lower_residual_pa, upper_residual_pa];
                }
                Err(error) => return Err(error),
            }
        }
        checkpoint(cx)?;
        if let Some(mut solution) = found {
            solution.total_flow_evaluations = performed;
            Ok(solution)
        } else {
            Err(PressureDrivenError::NoCirculationBracket {
                forward_residuals_pa: failures[0],
                reverse_residuals_pa: failures[1],
            })
        }
    }

    fn solve_circulation_direction(
        &self,
        cx: &Cx<'_>,
        walls: &[f64],
        resistances: &[LiquidResistanceSection<'_>],
        drive: Pressure,
        gravity: &SeriesGravity,
        direction: LiquidFlowDirection,
        local: LiquidFeedbackBudget,
        budget: PressureDrivenBudget,
    ) -> Result<PressureDrivenSolution<LiquidHeatMarch<'a>>, PressureDrivenError> {
        if direction == LiquidFlowDirection::Forward {
            return self.solve_pressure_gravity_driven(
                cx,
                walls,
                resistances,
                drive,
                gravity,
                local,
                budget,
            );
        }
        let sections: Vec<_> = self
            .sections
            .iter()
            .rev()
            .map(|s| LiquidHeatSection {
                region: &s.name,
                geometry: s.geometry,
                correlation: s.correlation,
            })
            .collect();
        let path = Self::try_new(
            cx,
            self.law,
            self.source,
            MassFlowRate::new(self.mass_flow),
            Temperature::new(self.inlet.temperature_k()),
            &sections,
            self.property_temperature,
        )?;
        let walls: Vec<_> = walls.iter().rev().copied().collect();
        let resistances: Vec<_> = resistances.iter().rev().copied().collect();
        let gravity = SeriesGravity::try_new(
            gravity.gravity_m_s2(),
            gravity.vertices_m().iter().rev().copied().collect(),
        )
        .map_err(PressureDrivenError::Gravity)?;
        path.solve_pressure_gravity_driven(
            cx,
            &walls,
            &resistances,
            Pressure::new(-drive.value()),
            &gravity,
            local,
            budget,
        )
    }
}
