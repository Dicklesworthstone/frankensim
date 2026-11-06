//! Local derivatives of one hardness-capped contact return map.
//!
//! With accepted plastic history fixed and a strictly separated active set,
//! lower pressures have derivative zero, capped pressures follow hardness,
//! and free pressures solve the original compliance principal system. This is
//! implicit differentiation of equilibrium, not differentiation of iterations.
//! Opening, cap activation and plastic-flow onset have no classical tangent here.

use fs_blake3::ContentHash;
use fs_exec::Cx;

use super::{
    FiniteGapPlasticAdvance, FiniteGapPlasticError, FiniteGapPlasticModel, FiniteGapPlasticState,
    box_residual, checkpoint,
};

mod control;
pub use control::*;

/// Explicit active-set separation and directional-equation accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGapPlasticSensitivityControls {
    /// Required distance of every free pressure from both bounds [Pa].
    pub minimum_pressure_margin_pa: f64,
    /// Required open gap or positive plastic-flow closure [m].
    /// Must exceed the original primal complementarity tolerance.
    pub minimum_displacement_margin_m: f64,
    /// Maximum recomputed free-equation residual [m per dimensionless parameter].
    /// Scale this with the chosen direction; it is not a forward-error bound.
    pub maximum_linear_residual_m: f64,
}

/// Signed input derivatives with respect to one dimensionless parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGapPlasticDirection<'a> {
    /// Rigid approach derivative [m].
    pub approach_m: f64,
    /// Homogeneous reduced-modulus derivative [Pa].
    pub reduced_modulus_pa: f64,
    /// Pressure-hardness derivative at each original cell [Pa].
    pub hardness_pa: &'a [f64],
    /// Relative recession derivative at each original cell [m].
    pub recession_m: &'a [f64],
}

/// Measured active-set slack, not a certified admissible step radius.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FiniteGapPlasticSensitivityMargins {
    /// Number of separated zero-pressure cells.
    pub open_cells: usize,
    /// Number of strictly interior pressures.
    pub elastic_cells: usize,
    /// Number of hardness-capped cells with positive incremental flow.
    pub flowing_cells: usize,
    /// Smallest distance of an elastic pressure from either bound [Pa].
    /// None when no elastic cells exist.
    pub minimum_elastic_pressure_margin_pa: Option<f64>,
    /// Smallest open gap or plastic-flow closure [m].
    /// None when every cell has an interior pressure.
    pub minimum_bound_displacement_margin_m: Option<f64>,
}

/// A refused tangent never changes the accepted or candidate contact history.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteGapPlasticSensitivityError {
    /// Invalid direction, controls, or mismatched primal response.
    Invalid(&'static str),
    /// Original owner, cancellation, factorization, or iterative-solve refusal.
    Contact(FiniteGapPlasticError),
    /// A classical fixed-active-set derivative is unresolved at this cell.
    SwitchingCell {
        /// Original grid-cell index.
        cell: usize,
        /// Margin and its units.
        quantity: &'static str,
        /// Measured active-set separation.
        margin: f64,
        /// Required strict lower bound.
        minimum: f64,
    },
    /// Direct application of the original kernel failed the derivative budget.
    LinearResidual {
        /// Largest free-equation defect [m per dimensionless parameter].
        residual_m: f64,
        /// Requested maximum defect.
        maximum_m: f64,
    },
    /// The prescribed force does not resolve a unique local approach.
    UnresolvedLoadStiffness {
        /// Computed contact force per unit approach [N/m].
        stiffness_n_per_m: f64,
        /// Required strict positive lower bound [N/m].
        minimum_n_per_m: f64,
    },
    /// The differentiated physical loading boundary failed its force budget.
    BoundaryResidual {
        /// Absolute force-balance defect [N per dimensionless parameter].
        residual_n: f64,
        /// Requested maximum defect.
        maximum_n: f64,
    },
}

impl From<FiniteGapPlasticError> for FiniteGapPlasticSensitivityError {
    fn from(error: FiniteGapPlasticError) -> Self {
        Self::Contact(error)
    }
}

impl core::fmt::Display for FiniteGapPlasticSensitivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gap plastic sensitivity: {self:?}")
    }
}

impl std::error::Error for FiniteGapPlasticSensitivityError {}

/// Borrowed linearization of one original equilibrium and its fixed predecessor.
///
/// Repeated directions reuse the prepared active set and half-space operator.
/// Each direction needs at most one free linear solve, with the original solver
/// and work limits. No full Jacobian or perturbed contact histories are stored.
#[derive(Debug)]
pub struct FiniteGapPlasticLinearization<'a> {
    model: &'a FiniteGapPlasticModel,
    before: &'a FiniteGapPlasticState,
    candidate: &'a FiniteGapPlasticState,
    opening_m: Vec<f64>,
    free: Vec<usize>,
    upper: Vec<usize>,
    controls: FiniteGapPlasticSensitivityControls,
    margins: FiniteGapPlasticSensitivityMargins,
}

/// Jacobian-vector product of the admitted local constitutive map.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticSensitivity {
    /// Exact original physical and numerical owner.
    pub model_identity: ContentHash,
    /// Fixed accepted plastic history.
    pub predecessor_identity: ContentHash,
    /// Equilibrium whose active set was differentiated.
    pub candidate_identity: ContentHash,
    /// Signed pressure derivative at every original cell [Pa].
    pub pressure_pa: Vec<f64>,
    /// Signed elastic opening derivative [m].
    pub elastic_opening_m: Vec<f64>,
    /// Signed new plastic-depth derivative, with accepted history fixed [m].
    pub plastic_depth_m: Vec<f64>,
    /// Signed resultant-force derivative [N].
    pub normal_force_n: f64,
    /// Signed elastic-storage derivative [J].
    pub elastic_energy_j: f64,
    /// Signed incremental plastic-work derivative [J], not heat.
    pub plastic_work_j: f64,
    /// Measured free-equation defect [m per dimensionless parameter].
    pub linear_residual_m: f64,
}

impl FiniteGapPlasticModel {
    /// Prepare derivatives of an ordinary advance without resolving contact.
    ///
    /// The original grid, datum and accepted plastic history are fixed.
    /// Recession and material directions describe the corresponding local
    /// parameterized equilibrium equations, not a finite history-transfer API.
    /// The caller still owns physical/source applicability of a proposed step.
    ///
    /// A capped cell must be actively acquiring indentation in this advance.
    /// An already yielded cell at its unloading corner is refused: loading and
    /// unloading do not have the same stiffness. Source selection (including
    /// the minimum of two equal body hardnesses) is not differentiated here.
    pub fn linearize_advance<'a>(
        &'a self,
        before: &'a FiniteGapPlasticState,
        advance: &'a FiniteGapPlasticAdvance,
        controls: FiniteGapPlasticSensitivityControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLinearization<'a>, FiniteGapPlasticSensitivityError> {
        checkpoint(cx)?;
        self.check_state_owner(before)?;
        self.check_state_owner(&advance.candidate)?;
        if !controls.minimum_pressure_margin_pa.is_finite()
            || controls.minimum_pressure_margin_pa <= 0.0
            || !controls.minimum_displacement_margin_m.is_finite()
            || controls.minimum_displacement_margin_m <= self.config.complementarity_tolerance_m
            || !controls.maximum_linear_residual_m.is_finite()
            || controls.maximum_linear_residual_m <= 0.0
        {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "positive derivative tolerance and resolved active-set margins",
            ));
        }
        let candidate = &advance.candidate;
        if advance.receipt.predecessor_identity != before.identity
            || candidate.recession_m != before.recession_m
        {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "linearization requires the original ordinary advance and predecessor",
            ));
        }
        // Read the opaque candidate, not caller-editable receipt arrays.
        let opening_m = self.apply_compliance(&candidate.pressure_pa, cx)?;
        let target = self.trial_closure(before, candidate.approach_m, &before.recession_m)?;
        let residual = box_residual(
            &candidate.pressure_pa,
            &opening_m,
            &target,
            &self.config.hardness_pa,
        );
        if residual > self.config.complementarity_tolerance_m {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "candidate does not satisfy the original predecessor equations",
            ));
        }
        let mut free = Vec::new();
        let mut upper = Vec::new();
        let mut margins = FiniteGapPlasticSensitivityMargins::default();
        for i in 0..target.len() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            let p = candidate.pressure_pa[i];
            let h = self.config.hardness_pa[i];
            let closure = target[i] - opening_m[i];
            let increment = candidate.plastic_depth_m[i] - before.plastic_depth_m[i];
            let expected = if p == h { closure.max(0.0) } else { 0.0 };
            if increment < 0.0
                || (increment - expected).abs() > self.config.complementarity_tolerance_m
            {
                return Err(FiniteGapPlasticSensitivityError::Invalid(
                    "candidate plastic history is not this predecessor's return map",
                ));
            }
            let (margin, minimum, quantity) = if p == 0.0 || p == h {
                let margin = if p == 0.0 {
                    margins.open_cells += 1;
                    -closure
                } else {
                    margins.flowing_cells += 1;
                    upper.push(i);
                    closure
                };
                margins.minimum_bound_displacement_margin_m = Some(
                    margins
                        .minimum_bound_displacement_margin_m
                        .map_or(margin, |v| v.min(margin)),
                );
                (
                    margin,
                    controls.minimum_displacement_margin_m,
                    "open gap or plastic-flow closure [m]",
                )
            } else {
                free.push(i);
                margins.elastic_cells += 1;
                let margin = p.min(h - p);
                margins.minimum_elastic_pressure_margin_pa = Some(
                    margins
                        .minimum_elastic_pressure_margin_pa
                        .map_or(margin, |v| v.min(margin)),
                );
                (
                    margin,
                    controls.minimum_pressure_margin_pa,
                    "elastic pressure distance from either bound [Pa]",
                )
            };
            if margin <= minimum {
                return Err(FiniteGapPlasticSensitivityError::SwitchingCell {
                    cell: i,
                    quantity,
                    margin,
                    minimum,
                });
            }
        }
        checkpoint(cx)?;
        Ok(FiniteGapPlasticLinearization {
            model: self,
            before,
            candidate,
            opening_m,
            free,
            upper,
            controls,
            margins,
        })
    }
}

impl FiniteGapPlasticLinearization<'_> {
    /// Actual opening/yield separation of the differentiated equilibrium.
    pub const fn margins(&self) -> FiniteGapPlasticSensitivityMargins {
        self.margins
    }

    /// Evaluate one signed direction through the fixed active set.
    ///
    /// For free cells, C_FF dp_F = da - dr + (dE/E) u - C_FU dH_U.
    /// The capped pressures follow dH; open pressures stay zero. Since
    /// dC = -(dE/E) C, du = C dp - (dE/E) u. Flowing cells then have
    /// dz = da - dr - du. Work follows d(H dz_increment), retaining the
    /// original increment as well as the change in new plastic depth.
    ///
    /// The numerical residual is measured against the original signed kernel.
    /// These are local derivatives of the pressure-saturation approximation,
    /// not finite-step admissibility, an adjoint, or a derivative across switches.
    pub fn apply(
        &self,
        direction: FiniteGapPlasticDirection<'_>,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticSensitivity, FiniteGapPlasticSensitivityError> {
        let pressure_pa = self.solve_direction(direction, 0.0, cx)?;
        self.finish_direction(direction, pressure_pa, cx)
    }

    // The finite loading spring uses the same principal operator with a
    // positive rank-one compliance; its resolved approach is recovered before
    // evaluating the physical opening, storage and plastic-work derivatives.
    fn solve_direction(
        &self,
        direction: FiniteGapPlasticDirection<'_>,
        rank_one: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteGapPlasticSensitivityError> {
        checkpoint(cx)?;
        let count = self.opening_m.len();
        if !direction.approach_m.is_finite()
            || !direction.reduced_modulus_pa.is_finite()
            || direction.hardness_pa.len() != count
            || direction.recession_m.len() != count
            || direction
                .hardness_pa
                .iter()
                .chain(direction.recession_m)
                .any(|v| !v.is_finite())
        {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "finite signed direction on the complete original cell grid",
            ));
        }
        let modulus_fraction = direction.reduced_modulus_pa / self.model.config.reduced_modulus_pa;
        let mut target = Vec::with_capacity(count);
        for (i, &opening) in self.opening_m.iter().enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            let value =
                modulus_fraction.mul_add(opening, direction.approach_m - direction.recession_m[i]);
            if !value.is_finite() {
                return Err(FiniteGapPlasticSensitivityError::Invalid(
                    "representable directional equilibrium",
                ));
            }
            target.push(value);
        }
        let mut pressure_pa = vec![0.0; count];
        for (slot, &i) in self.upper.iter().enumerate() {
            if slot % 64 == 0 {
                checkpoint(cx)?;
            }
            pressure_pa[i] = direction.hardness_pa[i];
        }
        if !self.free.is_empty() {
            let solution = self.model.solve_free_with_tolerance(
                &self.free,
                &self.upper,
                &pressure_pa,
                &target,
                rank_one,
                self.controls.maximum_linear_residual_m,
                cx,
            )?;
            for (slot, (&i, value)) in self.free.iter().zip(solution).enumerate() {
                if slot % 64 == 0 {
                    checkpoint(cx)?;
                }
                pressure_pa[i] = value;
            }
        }
        Ok(pressure_pa)
    }

    fn finish_direction(
        &self,
        direction: FiniteGapPlasticDirection<'_>,
        pressure_pa: Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticSensitivity, FiniteGapPlasticSensitivityError> {
        let count = self.opening_m.len();
        let modulus_fraction = direction.reduced_modulus_pa / self.model.config.reduced_modulus_pa;
        let mut elastic_opening_m = self.model.apply_compliance(&pressure_pa, cx)?;
        let mut plastic_depth_m = vec![0.0; count];
        for i in 0..count {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            elastic_opening_m[i] =
                (-modulus_fraction).mul_add(self.opening_m[i], elastic_opening_m[i]);
        }
        let mut linear_residual_m = 0.0_f64;
        for (slot, &i) in self.free.iter().enumerate() {
            if slot % 64 == 0 {
                checkpoint(cx)?;
            }
            let defect = elastic_opening_m[i] - (direction.approach_m - direction.recession_m[i]);
            if !defect.is_finite() {
                return Err(FiniteGapPlasticSensitivityError::Invalid(
                    "finite directional residual",
                ));
            }
            linear_residual_m = linear_residual_m.max(defect.abs());
        }
        if linear_residual_m > self.controls.maximum_linear_residual_m {
            return Err(FiniteGapPlasticSensitivityError::LinearResidual {
                residual_m: linear_residual_m,
                maximum_m: self.controls.maximum_linear_residual_m,
            });
        }
        for (slot, &i) in self.upper.iter().enumerate() {
            if slot % 64 == 0 {
                checkpoint(cx)?;
            }
            plastic_depth_m[i] =
                direction.approach_m - direction.recession_m[i] - elastic_opening_m[i];
        }
        let mut sums = [0.0; 3];
        let mut corrections = [0.0; 3];
        for i in 0..count {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            let old_increment = self.candidate.plastic_depth_m[i] - self.before.plastic_depth_m[i];
            let values = [
                pressure_pa[i],
                0.5 * pressure_pa[i].mul_add(
                    self.opening_m[i],
                    self.candidate.pressure_pa[i] * elastic_opening_m[i],
                ),
                direction.hardness_pa[i].mul_add(
                    old_increment,
                    self.model.config.hardness_pa[i] * plastic_depth_m[i],
                ),
            ];
            for k in 0..3 {
                let corrected = values[k] - corrections[k];
                let next = sums[k] + corrected;
                corrections[k] = (next - sums[k]) - corrected;
                sums[k] = next;
            }
        }
        let area = self.model.config.grid.cell_width_m * self.model.config.grid.cell_depth_m;
        for sum in &mut sums {
            *sum *= area;
        }
        if sums
            .iter()
            .chain(&elastic_opening_m)
            .chain(&plastic_depth_m)
            .any(|v| !v.is_finite())
        {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "finite directional pressure, storage and plastic work",
            ));
        }
        checkpoint(cx)?;
        Ok(FiniteGapPlasticSensitivity {
            model_identity: self.model.identity,
            predecessor_identity: self.before.identity,
            candidate_identity: self.candidate.identity,
            pressure_pa,
            elastic_opening_m,
            plastic_depth_m,
            normal_force_n: sums[0],
            elastic_energy_j: sums[1],
            plastic_work_j: sums[2],
            linear_residual_m,
        })
    }
}
