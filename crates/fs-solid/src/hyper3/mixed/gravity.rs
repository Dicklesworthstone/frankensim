//! Conservative spatial gravity on the original mixed Tet10 configuration.
//!
//! The same cached quadratic shapes, positive reference weights and element
//! densities supply gravity and the native consistent modal mass. At each
//! current point, U=sum rho0 w phi(x), f_ai=sum rho0 w N_a g_i(x), and
//! df_ai/du_bj=sum rho0 w N_a N_b dg_i/dx_j. No current-volume mass update,
//! corner lumping or cell-centroid substitution is used. Curved geometry and
//! nonlinear fields retain the caller's explicit quadrature approximation.
//!
//! The stationary field must satisfy the original GravityField3 contract.
//! Its domain and errors remain intact. Gravity acts only on displacement;
//! the original pressure constraint, gauge and constitutive split remain in
//! force. The optional uniform pressure acts on the complete closed material
//! boundary. Both loads ramp with the native static load factor and are held
//! during modes. This adapter supplies equilibrium and vibration, not motion.

use super::pressure::{self, UniformBoundaryPressure};
use super::vibration::incompressible::IncompressibleMixedVibration;
use super::vibration::{MixedUpVibration, MixedVibrationError, MixedVibrationSettings};
use super::{
    MixedPoint, MixedUpError, MixedUpEvaluation, MixedUpProblem, MixedUpSolution,
    PreparedMixedElement, PreparedMixedUpProblem, checkpoint, invalid, point_deformation,
};
use crate::hyper3::HyperTetError;
use crate::hyper3::gravity::{GravityField3, GravitySample3};
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::MixedFiniteStrainPotential;

/// Original stationary field and per-element reference mass density. Modal
/// inertia uses this same density; the field is never replaced by dead forces.
#[derive(Clone, Copy)]
pub struct HeldMixedGravity<'a> {
    reference_density_kg_m3: &'a [f64],
    field: &'a GravityField3<'a>,
}

impl core::fmt::Debug for HeldMixedGravity<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("HeldMixedGravity")
            .field("reference_density_kg_m3", &self.reference_density_kg_m3)
            .finish_non_exhaustive()
    }
}

impl<'a> HeldMixedGravity<'a> {
    pub const fn reference_density_kg_m3(self) -> &'a [f64] {
        self.reference_density_kg_m3
    }

    pub const fn field(self) -> &'a GravityField3<'a> {
        self.field
    }

    pub(super) fn sample(
        self,
        element_index: usize,
        point_index: usize,
        element: &PreparedMixedElement,
        point: &MixedPoint,
        nodes_m: &[[f64; 3]],
        displacement_m: &[[f64; 3]],
    ) -> Result<(f64, GravitySample3), MixedUpError> {
        let origin = element.connectivity[0];
        // Use the original Tet10 map's translation-safe relative coordinates;
        // a rounded partition-of-unity residual must not multiply a huge origin.
        let finite = |value: f64| {
            if value.is_finite() {
                Ok(value)
            } else {
                Err(invalid("mixed gravity position is unrepresentable"))
            }
        };
        let mut position = [0.0; 3];
        for i in 0..3 {
            position[i] = finite(nodes_m[origin][i] + displacement_m[origin][i])?;
            for a in 1..10 {
                let node = element.connectivity[a];
                let relative = finite(
                    (nodes_m[node][i] - nodes_m[origin][i])
                        + (displacement_m[node][i] - displacement_m[origin][i]),
                )?;
                position[i] = finite(point.shapes[a].mul_add(relative, position[i]))?;
            }
        }
        let mass = self.reference_density_kg_m3[element_index] * point.weight_m3;
        let sample = (self.field)(position)
            .and_then(GravitySample3::validate)
            .map_err(|error| match error {
                HyperTetError::Cancelled => MixedUpError::Cancelled,
                error => MixedUpError::Gravity {
                    element: element_index,
                    point: point_index,
                    error,
                },
            })?;
        Ok((mass, sample))
    }
}

/// External loading potentials, separate from the native material saddle
/// potential. The original field fixes the gravitational potential datum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedGravityPotentials {
    pub gravitational_j: f64,
    /// q V for the optional uniform boundary pressure; zero when absent.
    pub boundary_pressure_j: f64,
}

/// One immutable prepared mixed solid with its original conservative loads.
pub struct PreparedMixedGravityProblem<'a, L: MixedFiniteStrainPotential = Hyperelastic> {
    prepared: PreparedMixedUpProblem<'a, L>,
}

#[derive(Debug)]
pub struct MixedGravitySolution<'a, L: MixedFiniteStrainPotential = Hyperelastic> {
    problem: MixedUpProblem<'a, L>,
    gravity: HeldMixedGravity<'a>,
    boundary_pressure: Option<UniformBoundaryPressure>,
    equilibrium: MixedUpSolution,
    loading_potentials: MixedGravityPotentials,
}

impl<'a, L: MixedFiniteStrainPotential> MixedUpProblem<'a, L> {
    /// Sample the original field at each CURRENT quadratic material point.
    /// Densities are positive mass per original reference volume, one per
    /// element, and also define the consistent mass used by both mode routes.
    pub fn prepare_with_gravity(
        &self,
        reference_density_kg_m3: &'a [f64],
        field: &'a GravityField3<'a>,
        boundary_pressure: Option<UniformBoundaryPressure>,
        cx: &Cx<'_>,
    ) -> Result<PreparedMixedGravityProblem<'a, L>, MixedUpError> {
        checkpoint(cx)?;
        if reference_density_kg_m3.len() != self.elements.len() {
            return Err(invalid(
                "mixed gravity requires one reference density per element",
            ));
        }
        for (element, &density) in reference_density_kg_m3.iter().enumerate() {
            if element % 128 == 0 {
                checkpoint(cx)?;
            }
            if !density.is_finite() || density <= 0.0 {
                return Err(invalid(
                    "mixed gravity reference density must be finite and positive",
                ));
            }
        }
        let mut prepared = self.prepare(cx)?;
        // Admit the same positive point masses and accumulation order as the
        // native reference inertia, without allocating a second mass matrix.
        let mut reference_mass_kg = 0.0;
        for (element, &density) in prepared.mesh.elements.iter().zip(reference_density_kg_m3) {
            for point in &element.points {
                checkpoint(cx)?;
                let mass = density * point.weight_m3;
                reference_mass_kg += mass;
                if !mass.is_finite() || mass <= 0.0 || !reference_mass_kg.is_finite() {
                    return Err(invalid(
                        "mixed gravity reference point mass or total is unrepresentable",
                    ));
                }
            }
        }
        if boundary_pressure.is_some() {
            pressure::admit_closed_domain(self, cx)?;
        }
        prepared.boundary_pressure = boundary_pressure;
        prepared.gravity = Some(HeldMixedGravity {
            reference_density_kg_m3,
            field,
        });
        checkpoint(cx)?;
        Ok(PreparedMixedGravityProblem { prepared })
    }
}

impl<'a, L: MixedFiniteStrainPotential> PreparedMixedGravityProblem<'a, L> {
    pub const fn problem(&self) -> &MixedUpProblem<'a, L> {
        &self.prepared.problem
    }

    pub fn gravity(&self) -> HeldMixedGravity<'a> {
        self.prepared
            .gravity
            .expect("prepared with original gravity")
    }

    pub const fn boundary_pressure(&self) -> Option<UniformBoundaryPressure> {
        self.prepared.boundary_pressure
    }

    pub fn pressure_nodes(&self) -> &[usize] {
        &self.prepared.mesh.pressure_nodes
    }

    /// Full-load material force MINUS gravity and optional pressure force,
    /// with the complete consistent mixed Hessian when requested. The pressure
    /// equations and the material-only saddle potential retain their meanings.
    pub fn evaluate(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        with_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<MixedUpEvaluation, MixedUpError> {
        self.prepared
            .evaluate(displacement_m, pressure_pa, with_tangent, cx)
    }

    /// Integrate full-load potentials with the SAME original reference rule.
    /// Positive sampled Jacobians are required; this does not certify global
    /// injectivity or replace a constitutive-domain check by the native solve.
    pub fn loading_potentials(
        &self,
        displacement_m: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<MixedGravityPotentials, MixedUpError> {
        checkpoint(cx)?;
        if displacement_m.len() != self.prepared.mesh.nodes_count
            || displacement_m.as_flattened().iter().any(|u| !u.is_finite())
        {
            return Err(invalid("invalid mixed gravity potential displacement"));
        }
        let mut gravitational_j = 0.0;
        let mut volume_m3 = 0.0;
        for (e, element) in self.prepared.mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let f = point_deformation(element, point, displacement_m);
                let cof = pressure::cofactor(&f);
                let j = (0..3).map(|i| f[i] * cof[i]).sum::<f64>();
                if !j.is_finite() || j <= 0.0 {
                    return Err(invalid(
                        "mixed loading potential requires positive sampled det(F)",
                    ));
                }
                let (mass, sample) = self.gravity().sample(
                    e,
                    q,
                    element,
                    point,
                    self.prepared.problem.nodes_m,
                    displacement_m,
                )?;
                gravitational_j += mass * sample.potential_j_kg;
                volume_m3 += point.weight_m3 * j;
            }
        }
        let boundary_pressure_j = self
            .boundary_pressure()
            .map_or(0.0, |p| p.pressure_pa() * volume_m3);
        if !gravitational_j.is_finite()
            || !volume_m3.is_finite()
            || volume_m3 <= 0.0
            || !boundary_pressure_j.is_finite()
        {
            return Err(invalid(
                "mixed gravity loading potential is unrepresentable",
            ));
        }
        checkpoint(cx)?;
        Ok(MixedGravityPotentials {
            gravitational_j,
            boundary_pressure_j,
        })
    }

    pub fn solve(&self, cx: &Cx<'_>) -> Result<MixedGravitySolution<'a, L>, MixedUpError> {
        self.retain(self.prepared.solve(cx)?, cx)
    }

    /// Re-equilibrate without changing either supplied iterate. Native absolute
    /// supports, pressure gauge and all stable-KKT admission gates still apply.
    pub fn solve_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<MixedGravitySolution<'a, L>, MixedUpError> {
        self.retain(
            self.prepared.solve_from(displacement_m, pressure_pa, cx)?,
            cx,
        )
    }

    fn retain(
        &self,
        equilibrium: MixedUpSolution,
        cx: &Cx<'_>,
    ) -> Result<MixedGravitySolution<'a, L>, MixedUpError> {
        let loading_potentials = self.loading_potentials(&equilibrium.displacement_m, cx)?;
        Ok(MixedGravitySolution {
            problem: self.prepared.problem,
            gravity: self.gravity(),
            boundary_pressure: self.boundary_pressure(),
            equilibrium,
            loading_potentials,
        })
    }
}

impl<L: MixedFiniteStrainPotential> MixedGravitySolution<'_, L> {
    pub const fn problem(&self) -> &MixedUpProblem<'_, L> {
        &self.problem
    }
    pub const fn gravity(&self) -> HeldMixedGravity<'_> {
        self.gravity
    }
    pub const fn boundary_pressure(&self) -> Option<UniformBoundaryPressure> {
        self.boundary_pressure
    }
    pub const fn equilibrium(&self) -> &MixedUpSolution {
        &self.equilibrium
    }
    pub const fn loading_potentials(&self) -> MixedGravityPotentials {
        self.loading_potentials
    }
}

impl<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential> PreparedMixedGravityProblem<'a, L> {
    /// Follow the same mixed equilibrium with its original held loads.
    /// The path parameter scales current pressure/gravity, dead forces and
    /// the original support/gauge targets together.
    pub fn into_continuation(
        self,
        cx: &Cx<'_>,
    ) -> Result<
        super::continuation::MixedContinuation<'a, L>,
        super::continuation::MixedContinuationError,
    > {
        self.prepared.into_continuation(cx)
    }

    /// Finite-compliance modes of the loaded equilibrium. Density for inertia
    /// is the SAME retained reference density that supplies gravitational mass.
    pub fn solve_vibration(
        &self,
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.prepared.solve_vibration_impl(
            None,
            self.gravity().reference_density_kg_m3(),
            settings,
            cx,
        )
    }

    pub fn solve_vibration_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.prepared.solve_vibration_impl(
            Some((displacement_m, pressure_pa)),
            self.gravity().reference_density_kg_m3(),
            settings,
            cx,
        )
    }

    /// Exact-constraint modes with native pressure recovery, gauge handling and
    /// full original weak-volume residual checks. Pressure still has no mass.
    pub fn solve_incompressible_vibration(
        &self,
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<IncompressibleMixedVibration<'a, L>, MixedVibrationError> {
        self.prepared.solve_incompressible_vibration_impl(
            None,
            self.gravity().reference_density_kg_m3(),
            settings,
            cx,
        )
    }

    pub fn solve_incompressible_vibration_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<IncompressibleMixedVibration<'a, L>, MixedVibrationError> {
        self.prepared.solve_incompressible_vibration_impl(
            Some((displacement_m, pressure_pa)),
            self.gravity().reference_density_kg_m3(),
            settings,
            cx,
        )
    }
}

/// Add the negative external force and its displacement derivative to the
/// existing element block. Pressure rows receive no gravity or fictitious mass.
pub(super) fn add_point(
    shapes: &[f64; 10],
    mass: f64,
    sample: GravitySample3,
    residual: &mut [f64; 34],
    tangent: Option<&mut [[f64; 34]; 34]>,
) {
    for (a, &na) in shapes.iter().enumerate() {
        for i in 0..3 {
            residual[3 * a + i] -= mass * na * sample.acceleration_m_s2[i];
        }
    }
    if let Some(tangent) = tangent {
        for (a, &na) in shapes.iter().enumerate() {
            for (b, &nb) in shapes.iter().enumerate() {
                let weight = mass * na * nb;
                for i in 0..3 {
                    for j in 0..3 {
                        tangent[3 * a + i][3 * b + j] -=
                            weight * sample.acceleration_gradient_s2[i][j];
                    }
                }
            }
        }
    }
}
