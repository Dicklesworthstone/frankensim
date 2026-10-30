//! Source density drives loads through the existing layout and sizing owners.
//! Uniform member mass is lumped half to each endpoint. This equivalent nodal
//! load preserves total force and first moment; it resolves neither distributed
//! member bending nor the within-member axial stress gradient under gravity.

use super::{
    LayoutCase, MaterialCatalogAudit, MaterialLayout, MaterialLayoutSolution,
    MaterialLayoutSolveError, PdhgSettings, TrussConstructionError, checkpoint, invalid,
};
use fs_exec::Cx;

/// Original per-member masses and their equivalent applied nodal body forces.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialBodyLoads {
    /// In original ground-member order; zero area produces zero mass.
    pub member_mass_kg: Vec<f64>,
    pub mass_kg: f64,
    /// In original node order, including every supported node.
    pub nodal_forces_n: Vec<[f64; 2]>,
    pub total_force_n: [f64; 2],
}

/// A catalog design whose own material mass generated its actual solved loads.
#[derive(Debug)]
pub struct MaterialBodyLoadDesign {
    pub load_case: LayoutCase,
    pub body_loads: MaterialBodyLoads,
    /// Original member order; zero denotes a pruned/absent member.
    pub areas_m2: Vec<f64>,
    pub design: MaterialCatalogAudit,
    pub solution: MaterialLayoutSolution,
    /// Body-load solves after the initial external-load design.
    pub iterations: usize,
}

/// Refusal returns no inconsistent or unconverged design.
#[derive(Debug)]
pub enum BodyLoadDesignError {
    Construction(TrussConstructionError),
    Solve(MaterialLayoutSolveError),
    /// The common PDHG owner did not meet its caller-selected diagnostics.
    SolverNotConverged {
        iteration: usize,
        gap: f64,
        equilibrium_residual: f64,
    },
    /// Catalog, equilibrium, yield or Euler checks failed at this iteration.
    InvalidSizing {
        iteration: usize,
    },
    /// Exact catalog areas/topology did not settle within the admitted bound.
    IterationLimit {
        maximum_iterations: usize,
    },
}
impl From<TrussConstructionError> for BodyLoadDesignError {
    fn from(error: TrussConstructionError) -> Self {
        Self::Construction(error)
    }
}
impl From<MaterialLayoutSolveError> for BodyLoadDesignError {
    fn from(error: MaterialLayoutSolveError) -> Self {
        Self::Solve(error)
    }
}
impl core::fmt::Display for BodyLoadDesignError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material body-load design: {self:?}")
    }
}
impl std::error::Error for BodyLoadDesignError {}

fn sum(a: f64, b: f64) -> Result<f64, TrussConstructionError> {
    let next = a + b;
    if !next.is_finite() || (b != 0.0 && next == a) || (a != 0.0 && next == b) {
        return Err(invalid(
            "body-load accumulation",
            "must retain every finite nonzero contribution",
        ));
    }
    Ok(next)
}

impl MaterialLayout<'_> {
    /// Integrate original `density * area * length` and apply half each member's
    /// mass times `acceleration_m_s2` to each endpoint. Areas are actual current
    /// sections in original member order; zero means the member is absent.
    /// The explicit acceleration may represent gravity or an equivalent inertial
    /// load; no sign convention, gravity constant or second density is inferred.
    ///
    /// # Errors
    /// Malformed sections, nonfinite acceleration, unrepresentable mass/force,
    /// lost nonzero contributions or cancellation refuse before publication.
    pub fn body_loads(
        &self,
        areas_m2: &[f64],
        acceleration_m_s2: [f64; 2],
        cx: &Cx<'_>,
    ) -> Result<MaterialBodyLoads, TrussConstructionError> {
        checkpoint(cx)?;
        if areas_m2.len() != self.assignments.len() {
            return Err(TrussConstructionError::VectorLength {
                field: "body-load member areas",
                expected: self.assignments.len(),
                actual: areas_m2.len(),
            });
        }
        if acceleration_m_s2.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "body acceleration",
                "must contain finite components in m/s2",
            ));
        }
        let mut member_mass_kg = Vec::with_capacity(areas_m2.len());
        let mut nodal_forces_n = vec![[0.0; 2]; self.ground.nodes().len()];
        let mut mass_kg = 0.0;
        for (member, &area) in areas_m2.iter().enumerate() {
            if member % 256 == 0 {
                checkpoint(cx)?;
            }
            if !area.is_finite() || area < 0.0 {
                return Err(invalid(
                    "body-load member area",
                    "must be finite and nonnegative",
                ));
            }
            if area == 0.0 {
                member_mass_kg.push(0.0);
                continue;
            }
            let volume = area * self.ground.lengths()[member];
            let mass = volume * self.materials[self.assignments[member]].density;
            if !volume.is_finite() || volume <= 0.0 || !mass.is_finite() || mass <= 0.0 {
                return Err(invalid(
                    "member body mass",
                    "positive sections require representable positive volume and mass",
                ));
            }
            mass_kg = sum(mass_kg, mass)?;
            member_mass_kg.push(mass);
            let (a, b) = self.ground.members()[member];
            for (axis, acceleration) in acceleration_m_s2.into_iter().enumerate() {
                let force = mass * acceleration;
                let half = 0.5 * force;
                if !force.is_finite() || (acceleration != 0.0 && half == 0.0) {
                    return Err(invalid(
                        "member body force",
                        "nonzero acceleration requires finite representable endpoint forces",
                    ));
                }
                nodal_forces_n[a][axis] = sum(nodal_forces_n[a][axis], half)?;
                nodal_forces_n[b][axis] = sum(nodal_forces_n[b][axis], half)?;
            }
        }
        let mut total_force_n = [0.0; 2];
        for (node, force) in nodal_forces_n.iter().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for axis in 0..2 {
                total_force_n[axis] = sum(total_force_n[axis], force[axis])?;
            }
        }
        checkpoint(cx)?;
        Ok(MaterialBodyLoads {
            member_mass_kg,
            mass_kg,
            nodal_forces_n,
            total_force_n,
        })
    }

    /// Close catalog sizing against its own material-derived body loads.
    ///
    /// First solve and size the unchanged external-load case. Then derive body
    /// forces from those actual snapped sections, add them ONCE to the original
    /// loads, and use the same LP/refit/yield/Euler/catalog owners again. Only
    /// exact equality of the complete section vector (including absent members)
    /// accepts. The returned loads, mass and design therefore describe one state.
    /// Original support loads, material assignments and assembly limits persist.
    ///
    /// This bounded fixed-point iteration can cycle or lack a feasible catalog;
    /// it does not establish a global self-weight optimum. It retains the
    /// existing synchronous PDHG/refit cancellation limitation, with checkpoints
    /// before/after those owners and throughout body-load reduction.
    ///
    /// # Errors
    /// Zero iteration allowance, source/numerical refusal, unconverged PDHG,
    /// failed sizing, unsettled sections or cancellation publish no design.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_body_loads(
        &self,
        acceleration_m_s2: [f64; 2],
        catalog: &[f64],
        prune_frac: f64,
        settings: PdhgSettings,
        maximum_iterations: usize,
        cx: &Cx<'_>,
    ) -> Result<MaterialBodyLoadDesign, BodyLoadDesignError> {
        checkpoint(cx)?;
        if maximum_iterations == 0 || acceleration_m_s2.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "body-load solve controls",
                "need a positive iteration bound and finite acceleration",
            )
            .into());
        }
        let (_, _, mut areas) = self.solve_sized(catalog, prune_frac, settings, 0, cx)?;
        for iteration in 1..=maximum_iterations {
            checkpoint(cx)?;
            let body_loads = self.body_loads(&areas, acceleration_m_s2, cx)?;
            let mut loads = self.case.loads().to_vec();
            for (node, (load, body)) in loads.iter_mut().zip(&body_loads.nodal_forces_n).enumerate()
            {
                if node % 256 == 0 {
                    checkpoint(cx)?;
                }
                for axis in 0..2 {
                    load[axis] = sum(load[axis], body[axis])?;
                }
            }
            let load_case = LayoutCase::try_new(
                self.case.supported().to_vec(),
                loads,
                self.case.node_count(),
            )?;
            let loaded = MaterialLayout::try_new(
                self.ground,
                &load_case,
                self.materials,
                self.assignments,
                self.limits,
                cx,
            )?;
            let (solution, design, next_areas) =
                loaded.solve_sized(catalog, prune_frac, settings, iteration, cx)?;
            if areas == next_areas {
                checkpoint(cx)?;
                return Ok(MaterialBodyLoadDesign {
                    load_case,
                    body_loads,
                    areas_m2: next_areas,
                    design,
                    solution,
                    iterations: iteration,
                });
            }
            areas = next_areas;
        }
        Err(BodyLoadDesignError::IterationLimit { maximum_iterations })
    }

    fn solve_sized(
        &self,
        catalog: &[f64],
        prune_frac: f64,
        settings: PdhgSettings,
        iteration: usize,
        cx: &Cx<'_>,
    ) -> Result<(MaterialLayoutSolution, MaterialCatalogAudit, Vec<f64>), BodyLoadDesignError> {
        checkpoint(cx)?;
        let solution = self.solve(settings)?;
        checkpoint(cx)?;
        let numerical = &solution.numerical;
        if !(numerical.gap < settings.gap_tol && numerical.eq_residual < settings.gap_tol) {
            return Err(BodyLoadDesignError::SolverNotConverged {
                iteration,
                gap: numerical.gap,
                equilibrium_residual: numerical.eq_residual,
            });
        }
        let design = self.size_and_snap(&solution.forces, catalog, prune_frac, cx)?;
        if !design.audit.all_pass || design.mass_kg.is_none() || design.volume_m3.is_none() {
            return Err(BodyLoadDesignError::InvalidSizing { iteration });
        }
        let mut areas = vec![0.0; self.assignments.len()];
        for (index, member) in design.audit.members.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            areas[member.index] = member.area_catalog;
        }
        checkpoint(cx)?;
        Ok((solution, design, areas))
    }
}
