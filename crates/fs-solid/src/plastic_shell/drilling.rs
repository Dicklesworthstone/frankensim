//! Retained numerical drilling stiffness for shells with evolving material laws.
//!
//! Temperature-dependent elastic moduli and thermal strains may change a
//! section's virgin tangent. They must not change the numerical spin gauge in
//! the middle of a physical trajectory. This owner retains the original gauge
//! weights; force, tangent and work still use the same native shell rows.
//! Gauge storage and endpoint integration loss are numerical, never heat.

use super::dynamics::{PlasticShellReferenceInertia, ShellInertialEnergy};
use super::{
    Element, Operator, SmallStrainShellProblem, SmallStrainShellSolution, checkpoint, dot, invalid,
};
use crate::plastic3::SmallStrainTetError;
use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::plastic::reduction::PlaneStressWorkspace;
use fs_material::plastic::section::LayeredPlaneStressState;
use fs_plate::shell::ShellMesh;
use std::cell::RefCell;

#[derive(Clone, Debug)]
struct FacetGauge {
    bounds: Vec<[f64; 2]>,
    weight: f64,
}

/// Original numerical gauge, bound to its reference mesh, physical layer bounds,
/// coordinate length and drilling coefficient. Material coefficients may change;
/// these reference data may not. No material history or stiffness matrix is copied.
#[derive(Clone, Debug)]
pub struct PlasticShellReferenceDrilling<'m> {
    mesh: &'m ShellMesh,
    coordinate_length_m: f64,
    drilling_alpha: f64,
    facets: Vec<FacetGauge>,
}

impl<'m> PlasticShellReferenceDrilling<'m> {
    /// Bind the original gauge to a mesh whose lifetime is independent of the
    /// temporary section-law borrows in `problem`. The mesh must match exactly.
    pub fn prepare<L: SmallStrainLaw>(
        mesh: &'m ShellMesh,
        problem: &SmallStrainShellProblem<'_, L>,
        cx: &Cx<'_>,
    ) -> Result<Self, SmallStrainTetError> {
        checkpoint(cx)?;
        if !std::ptr::eq(mesh, problem.mesh) && mesh != problem.mesh {
            return Err(invalid(
                "reference drilling mesh differs from its original shell",
            ));
        }
        let elements = problem.prepare(cx)?;
        let mut facets = Vec::with_capacity(elements.len());
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            facets.push(FacetGauge {
                bounds: problem
                    .sections
                    .at(e)
                    .layers()
                    .iter()
                    .map(|layer| [layer.bottom_z_m, layer.top_z_m])
                    .collect(),
                weight: element.drilling_weight,
            });
        }
        Ok(Self {
            mesh,
            coordinate_length_m: problem.coordinate_length_m,
            drilling_alpha: problem.drilling_alpha,
            facets,
        })
    }

    pub(super) fn admit<L: SmallStrainLaw>(
        &self,
        problem: &SmallStrainShellProblem<'_, L>,
        cx: &Cx<'_>,
    ) -> Result<(), SmallStrainTetError> {
        checkpoint(cx)?;
        if (!std::ptr::eq(self.mesh, problem.mesh) && self.mesh != problem.mesh)
            || self.coordinate_length_m != problem.coordinate_length_m
            || self.drilling_alpha != problem.drilling_alpha
            || self.facets.len() != problem.mesh.tris.len()
        {
            return Err(invalid(
                "retained drilling belongs to different reference geometry, coordinates or gauge",
            ));
        }
        // The common preparation admits section counts before entering here.
        for (e, facet) in self.facets.iter().enumerate() {
            checkpoint(cx)?;
            let layers = problem.sections.at(e).layers();
            if layers.len() != facet.bounds.len()
                || layers
                    .iter()
                    .zip(&facet.bounds)
                    .any(|(layer, bounds)| [layer.bottom_z_m, layer.top_z_m] != *bounds)
            {
                return Err(invalid(
                    "retained drilling belongs to different physical layer bounds",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn weight(&self, element: usize) -> f64 {
        self.facets[element].weight
    }
}

/// Numerical drilling storage and work [J]. None of these fields is material heat.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellDrillingEnergy {
    pub initial_j: f64,
    pub final_j: f64,
    pub work_j: f64,
    pub numerical_loss_j: f64,
}

impl ShellDrillingEnergy {
    /// Shared kernel for retained-reference and fixed-material shell audits.
    pub(super) fn add(&mut self, element: &Element, old: &[f64; 18], new: &[f64; 18]) {
        for row in &element.drilling {
            let before = dot(row, old);
            let after = dot(row, new);
            let delta = after - before;
            let weight = element.drilling_weight;
            self.initial_j += 0.5 * weight * before * before;
            self.final_j += 0.5 * weight * after * after;
            self.work_j += weight * after * delta;
            self.numerical_loss_j += 0.5 * weight * delta * delta;
        }
    }

    fn admit(self) -> Result<Self, SmallStrainTetError> {
        let scale = self
            .initial_j
            .max(self.final_j)
            .max(self.work_j.abs())
            .max(self.numerical_loss_j);
        if [self.initial_j, self.final_j, self.numerical_loss_j]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
            || !self.work_j.is_finite()
            || (self.work_j - (self.final_j - self.initial_j) - self.numerical_loss_j).abs()
                > 4096.0 * f64::EPSILON * scale
        {
            return Err(invalid(
                "retained drilling misses its numerical storage/work identity",
            ));
        }
        Ok(self)
    }
}

/// A single accepted endpoint with its independently retained numerical gauge.
#[derive(Debug)]
pub struct PlasticShellReferenceDrillingStep<S> {
    pub solution: SmallStrainShellSolution<S>,
    pub drilling: ShellDrillingEnergy,
}

/// Motion on the original reference mass and drilling stiffness. Material heat
/// remains the responsibility of the native material/thermal energy owner.
#[derive(Debug)]
pub struct PlasticShellReferenceDrillingDynamicStep<S> {
    pub solution: SmallStrainShellSolution<S>,
    pub velocity: Vec<[f64; 6]>,
    pub inertia: ShellInertialEnergy,
    pub drilling: ShellDrillingEnergy,
}

impl<'m, L: SmallStrainLaw> SmallStrainShellProblem<'m, L> {
    /// Capture the virgin numerical gauge once, before material coefficients evolve.
    pub fn prepare_reference_drilling(
        &self,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellReferenceDrilling<'m>, SmallStrainTetError> {
        PlasticShellReferenceDrilling::prepare(self.mesh, self, cx)
    }

    /// Continue one factor-one endpoint with the original numerical stiffness.
    /// Temperature-dependent laws retain their own frozen old histories. No
    /// thermodynamic work is inferred from this purely numerical gauge audit.
    pub fn solve_from_with_reference_drilling(
        &self,
        coordinates: &[[f64; 6]],
        states: &[LayeredPlaneStressState<L::State>],
        reference: &PlasticShellReferenceDrilling<'_>,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellReferenceDrillingStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        if self.load_path != [1.0] {
            return Err(invalid(
                "retained drilling requires one factor-one shell endpoint",
            ));
        }
        let elements = self.prepare_with_reference_drilling(Some(reference), cx)?;
        let solution = self.solve_prepared(&elements, Some((coordinates, states)), None, cx)?;
        let drilling =
            self.audit_reference_drilling(&elements, Some(coordinates), &solution.coordinates, cx)?;
        Ok(PlasticShellReferenceDrillingStep { solution, drilling })
    }

    /// Advance with unchanged reference mass and numerical drilling stiffness.
    /// The native material laws must represent this same physical `dt_s`.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_backward_euler_with_reference_drilling(
        &self,
        initial: Option<&SmallStrainShellSolution<L::State>>,
        velocity: &[[f64; 6]],
        inertia: &PlasticShellReferenceInertia<'_>,
        drilling: &PlasticShellReferenceDrilling<'_>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellReferenceDrillingDynamicStep<L::State>, SmallStrainTetError> {
        let elements = self.prepare_with_reference_drilling(Some(drilling), cx)?;
        let step =
            self.advance_backward_euler_prepared(&elements, initial, velocity, inertia, dt_s, cx)?;
        let drilling = self.audit_reference_drilling(
            &elements,
            initial.map(|old| old.coordinates.as_slice()),
            &step.solution.coordinates,
            cx,
        )?;
        Ok(PlasticShellReferenceDrillingDynamicStep {
            solution: step.solution,
            velocity: step.velocity,
            inertia: step.inertia,
            drilling,
        })
    }

    pub(super) fn audit_reference_drilling(
        &self,
        elements: &[Element],
        old: Option<&[[f64; 6]]>,
        new: &[[f64; 6]],
        cx: &Cx<'_>,
    ) -> Result<ShellDrillingEnergy, SmallStrainTetError> {
        checkpoint(cx)?;
        let n = self.mesh.nodes.len();
        if new.len() != n || old.is_some_and(|old| old.len() != n) {
            return Err(invalid(
                "complete drilling-work endpoint coordinates required",
            ));
        }
        let length = self.coordinate_length_m;
        let scaled = |value, component| value * if component < 3 { 1.0 } else { length };
        let before: Vec<_> = (0..6 * n)
            .map(|d| scaled(old.map_or(0.0, |q| q[d / 6][d % 6]), d % 6))
            .collect();
        let after: Vec<_> = new
            .iter()
            .flat_map(|q| q.iter().enumerate().map(|(c, &value)| scaled(value, c)))
            .collect();
        if before.iter().chain(&after).any(|value| !value.is_finite()) {
            return Err(invalid(
                "scaled drilling-work coordinates are unrepresentable",
            ));
        }
        let operator = Operator {
            problem: self,
            elements,
            workspace: RefCell::new(
                PlaneStressWorkspace::new().map_err(|e| invalid(e.to_string()))?,
            ),
            inertia: None,
        };
        let mut account = ShellDrillingEnergy::default();
        for element in elements {
            checkpoint(cx)?;
            let old_local = operator.local_coordinates(element, &before)?;
            let new_local = operator.local_coordinates(element, &after)?;
            account.add(element, &old_local, &new_local);
        }
        checkpoint(cx)?;
        account.admit()
    }
}
