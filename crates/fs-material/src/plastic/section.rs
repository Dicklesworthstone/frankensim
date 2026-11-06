//! Retained plane-stress material histories through a thin section.
//!
//! The generalized engineering strain is `[ex, ey, gamma_xy, kx, ky, kxy]`.
//! At signed depth `z`, the native three-dimensional law sees the in-plane
//! tensor strain `[ex + z kx, ey + z ky, (gamma_xy + z kxy)/2]`. The original
//! [`PlaneStressReduction`] supplies transverse equilibrium, the complete
//! material history, and the algorithmic Schur-complement tangent. Integrating
//! stress with weights `1` and `z` gives work-conjugate membrane forces and
//! bending moments, including membrane/bending coupling after uneven yield.
//!
//! Every thickness point retains its own accepted history. All evaluations
//! borrow that history; only the returned complete candidate may be committed.
//! A timed material advances its declared interval once at acceptance, never
//! once per local Newton iteration. No constitutive update is duplicated here.
//!
//! Layer bounds share one caller-chosen reference surface and one material
//! frame. Layers must be contiguous; voids and independent layer rotations are
//! not inferred. This is a small-strain plane-stress section, with no transverse
//! shear, finite-rotation kinematics, delamination or thermal transport.
//! Gauss integration is exact for elastic `A/B/D`, but a moving yield front
//! requires thickness refinement; finite sampling does not certify its error.

use super::reduction::{
    PlaneStressReduction, PlaneStressWorkspace, ReducedState, ReductionConfig, ReductionError,
};
use crate::{SmallStrainLaw, tensor};
use core::fmt;

/// Native mechanical storage and physical/numerical loss of accepted sections.
pub mod energy;

/// One ordered physical layer and its original native material law.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneStressLayer<L> {
    /// Signed bottom coordinate from the common reference surface [m].
    pub bottom_z_m: f64,
    /// Signed top coordinate, strictly above the bottom [m].
    pub top_z_m: f64,
    /// Uniform positive number of thickness integration intervals in this layer.
    pub subdivisions: usize,
    /// Original immutable law, in the section axes; references are supported.
    pub law: L,
}

/// Explicit integration rule on every declared layer subdivision.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThicknessQuadrature {
    /// Two-point Gauss integration: exact through cubic thickness polynomials.
    Gauss2,
    /// Three-point Gauss integration: exact through quintic polynomials.
    Gauss3,
}

impl ThicknessQuadrature {
    /// Number of independent histories per subdivision.
    #[must_use]
    pub const fn points_per_subdivision(self) -> usize {
        match self {
            Self::Gauss2 => 2,
            Self::Gauss3 => 3,
        }
    }

    fn rule(self) -> ([(f64, f64); 3], usize) {
        match self {
            Self::Gauss2 => {
                let x = fs_math::det::sqrt(1.0 / 3.0);
                ([(-x, 1.0), (x, 1.0), (0.0, 0.0)], 2)
            }
            Self::Gauss3 => {
                let x = fs_math::det::sqrt(3.0 / 5.0);
                ([(-x, 5.0 / 9.0), (0.0, 8.0 / 9.0), (x, 5.0 / 9.0)], 3)
            }
        }
    }
}

/// Cold allocation bounds, independent of material applicability and accuracy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SectionBudget {
    /// Maximum number of original material layers.
    pub maximum_layers: usize,
    /// Maximum total retained thickness points per section state.
    pub maximum_quadrature_points: usize,
}

/// Prepared physical thickness quadrature, in layer and increasing-depth order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThicknessPoint {
    /// Index of the original layer; material identities are never averaged.
    pub layer: usize,
    /// Signed depth from the original reference surface [m].
    pub z_m: f64,
    /// Strictly positive reference thickness integration weight [m].
    pub weight_m: f64,
}

/// Complete accepted generalized strain and all native material histories.
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredPlaneStressState<S> {
    /// `[ex, ey, gamma_xy, kx, ky, kxy]`; curvatures are in inverse metres.
    pub generalized_strain: [f64; 6],
    /// Full tensor strains and native histories in [`ThicknessPoint`] order.
    pub points: Vec<ReducedState<S>>,
}

/// One frozen-history section response. Acceptance belongs to the caller.
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredPlaneStressResponse<S> {
    /// Complete candidate, never written into the borrowed accepted state.
    pub state: LayeredPlaneStressState<S>,
    /// `[Nx, Ny, Nxy, Mx, My, Mxy]`, forces in N/m and moments in N.
    /// Their ordinary dot product with generalized strain is work per area.
    pub resultants: [f64; 6],
    /// Exact derivative of `resultants` at the same incoming native history.
    /// Engineering shear is already included. Nonsymmetric native tangents
    /// remain nonsymmetric; no averaging or symmetrization is applied.
    pub tangent: [[f64; 6]; 6],
}

/// Geometry, state, numerical or native transverse-equilibrium refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SectionError {
    /// Invalid geometry, generalized strain, state binding or arithmetic.
    Invalid(&'static str),
    /// A caller-declared allocation ceiling was exceeded before allocation.
    Budget {
        /// Resource being bounded.
        quantity: &'static str,
        /// Requested count.
        requested: usize,
        /// Caller-declared ceiling.
        limit: usize,
    },
    /// The shared plane-stress reduction refused this integration point.
    Reduction {
        /// Global thickness point index, or zero during initial configuration.
        point: usize,
        /// Original reduction error, without replacing its meaning.
        source: ReductionError,
    },
}

impl fmt::Display for SectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "plastic plane-stress section: {what}"),
            Self::Budget {
                quantity,
                requested,
                limit,
            } => write!(
                f,
                "plastic plane-stress section {quantity}: requested {requested}, limit {limit}"
            ),
            Self::Reduction { point, source } => {
                write!(f, "plastic plane-stress section point {point}: {source}")
            }
        }
    }
}

impl std::error::Error for SectionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Reduction { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Immutable layer geometry, quadrature and original material bindings.
///
/// Geometry is prepared once. Each evaluation allocates one candidate state
/// vector; [`Self::evaluate_with_workspace`] also reuses the original local
/// factorization storage. Native history types retain their own allocation
/// behavior. The section imposes no extra heap allocation inside local Newton.
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredPlaneStressSection<L> {
    layers: Vec<PlaneStressLayer<L>>,
    quadrature: ThicknessQuadrature,
    points: Vec<ThicknessPoint>,
    reduction: ReductionConfig,
}

impl<L: SmallStrainLaw> LayeredPlaneStressSection<L> {
    /// Admit ordered contiguous layers, positive quadrature weights, explicit
    /// point budgets and the original bounded plane-stress solve controls.
    ///
    /// Adjacent layers must share exactly the same interface coordinate. The
    /// common reference surface may lie inside or outside the physical stack.
    pub fn new(
        layers: Vec<PlaneStressLayer<L>>,
        quadrature: ThicknessQuadrature,
        reduction: ReductionConfig,
        budget: SectionBudget,
    ) -> Result<Self, SectionError> {
        if layers.is_empty() {
            return Err(SectionError::Invalid(
                "at least one material layer is required",
            ));
        }
        if layers.len() > budget.maximum_layers {
            return Err(SectionError::Budget {
                quantity: "layers",
                requested: layers.len(),
                limit: budget.maximum_layers,
            });
        }
        // Admission is owned by the same reducer that will evaluate the law.
        PlaneStressReduction::new(&layers[0].law, reduction)
            .map_err(|source| SectionError::Reduction { point: 0, source })?;
        let order = quadrature.points_per_subdivision();
        let mut point_count = 0usize;
        for (index, layer) in layers.iter().enumerate() {
            let thickness = layer.top_z_m - layer.bottom_z_m;
            if !layer.bottom_z_m.is_finite()
                || !layer.top_z_m.is_finite()
                || !thickness.is_finite()
                || thickness <= 0.0
                || layer.subdivisions == 0
            {
                return Err(SectionError::Invalid(
                    "layer bounds and subdivision count must define positive finite thickness",
                ));
            }
            if index > 0 && layer.bottom_z_m != layers[index - 1].top_z_m {
                return Err(SectionError::Invalid(
                    "layers must be ordered and contiguous at identical interface coordinates",
                ));
            }
            point_count = layer
                .subdivisions
                .checked_mul(order)
                .and_then(|count| point_count.checked_add(count))
                .ok_or(SectionError::Invalid("thickness point count overflow"))?;
            if point_count > budget.maximum_quadrature_points {
                return Err(SectionError::Budget {
                    quantity: "quadrature points",
                    requested: point_count,
                    limit: budget.maximum_quadrature_points,
                });
            }
        }
        let span = layers[layers.len() - 1].top_z_m - layers[0].bottom_z_m;
        if !span.is_finite() || span <= 0.0 {
            return Err(SectionError::Invalid(
                "total section thickness is unrepresentable",
            ));
        }
        let (rule, _) = quadrature.rule();
        let mut points = Vec::with_capacity(point_count);
        for (index, layer) in layers.iter().enumerate() {
            let thickness = layer.top_z_m - layer.bottom_z_m;
            for part in 0..layer.subdivisions {
                let bottom =
                    layer.bottom_z_m + thickness * (part as f64 / layer.subdivisions as f64);
                let top = if part + 1 == layer.subdivisions {
                    layer.top_z_m
                } else {
                    layer.bottom_z_m + thickness * ((part + 1) as f64 / layer.subdivisions as f64)
                };
                let half = 0.5 * (top - bottom);
                let center = bottom + half;
                for &(xi, weight) in &rule[..order] {
                    let z_m = center + half * xi;
                    let weight_m = half * weight;
                    if !z_m.is_finite()
                        || z_m <= bottom
                        || z_m >= top
                        || !weight_m.is_finite()
                        || weight_m <= 0.0
                    {
                        return Err(SectionError::Invalid(
                            "thickness quadrature is not representable inside its interval",
                        ));
                    }
                    points.push(ThicknessPoint {
                        layer: index,
                        z_m,
                        weight_m,
                    });
                }
            }
        }
        Ok(Self {
            layers,
            quadrature,
            points,
            reduction,
        })
    }

    /// Original laws and layer bounds, without averaging or property extraction.
    #[must_use]
    pub fn layers(&self) -> &[PlaneStressLayer<L>] {
        &self.layers
    }

    /// Immutable physical quadrature and native history ordering.
    #[must_use]
    pub fn points(&self) -> &[ThicknessPoint] {
        &self.points
    }

    /// Total number of independent three-dimensional material histories.
    #[must_use]
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// Original local strain domain and transverse-equilibrium controls.
    #[must_use]
    pub const fn reduction_config(&self) -> ReductionConfig {
        self.reduction
    }

    /// Whether every original law declares a symmetric work-conjugate tangent.
    #[must_use]
    pub fn tangent_symmetric(&self) -> bool {
        self.layers
            .iter()
            .all(|layer| layer.law.admissibility().tangent_symmetric)
    }

    /// Virgin native histories at zero full strain. A law with eigenstrain can
    /// require an initial equilibrium solve; virgin does not imply stress-free.
    #[must_use]
    pub fn initial_state(&self) -> LayeredPlaneStressState<L::State> {
        LayeredPlaneStressState {
            generalized_strain: [0.0; 6],
            points: self
                .points
                .iter()
                .map(|point| ReducedState {
                    strain: [0.0; 6],
                    history: self.layers[point.layer].law.initial_state(),
                })
                .collect(),
        }
    }

    /// Validate the old geometric strain binding without advancing material
    /// time or projecting an accepted viscous state onto an equilibrium surface.
    /// Native internal-variable admissibility remains the law owner's contract.
    pub fn validate_state(
        &self,
        old: &LayeredPlaneStressState<L::State>,
    ) -> Result<(), SectionError> {
        self.validate_generalized(&old.generalized_strain)?;
        if old.points.len() != self.points.len() {
            return Err(SectionError::Invalid(
                "section history count does not match quadrature",
            ));
        }
        for (point, state) in self.points.iter().zip(&old.points) {
            let expected = strain_at(&old.generalized_strain, point.z_m);
            if [state.strain[0], state.strain[1], state.strain[3]] != expected
                || state.strain.iter().any(|value| !value.is_finite())
                || tensor::contract(&state.strain, &state.strain)
                    > self.reduction.strain_limit * self.reduction.strain_limit
            {
                return Err(SectionError::Invalid(
                    "retained full strain is out of domain or inconsistent with section kinematics",
                ));
            }
        }
        Ok(())
    }

    fn validate_generalized(&self, strain: &[f64; 6]) -> Result<(), SectionError> {
        if strain.iter().any(|value| !value.is_finite()) {
            return Err(SectionError::Invalid("generalized strain must be finite"));
        }
        // In-plane strain is affine through thickness. Bound its norm at every
        // material face too, so interior Gauss nodes cannot hide a large face
        // strain. Full transverse strain is admitted by the native point solve.
        let limit_squared = self.reduction.strain_limit * self.reduction.strain_limit;
        for layer in &self.layers {
            for z in [layer.bottom_z_m, layer.top_z_m] {
                let s = strain_at(strain, z);
                let norm_squared = s[0] * s[0] + s[1] * s[1] + 2.0 * s[2] * s[2];
                if !norm_squared.is_finite() || norm_squared > limit_squared {
                    return Err(SectionError::Invalid(
                        "in-plane layer-face strain exceeds the declared small-strain domain",
                    ));
                }
            }
        }
        Ok(())
    }

    /// One immutable candidate using a fresh fixed-size transverse workspace.
    pub fn evaluate(
        &self,
        strain: [f64; 6],
        old: &LayeredPlaneStressState<L::State>,
    ) -> Result<LayeredPlaneStressResponse<L::State>, SectionError> {
        let mut workspace = PlaneStressWorkspace::new()
            .map_err(|source| SectionError::Reduction { point: 0, source })?;
        self.evaluate_with_workspace(strain, old, &mut workspace)
    }

    /// Integrate native plane stress and its consistent derivative using the
    /// same incoming history at every point. A late refusal discards the whole
    /// candidate; the borrowed accepted state is unchanged.
    pub fn evaluate_with_workspace(
        &self,
        strain: [f64; 6],
        old: &LayeredPlaneStressState<L::State>,
        workspace: &mut PlaneStressWorkspace,
    ) -> Result<LayeredPlaneStressResponse<L::State>, SectionError> {
        self.validate_state(old)?;
        self.validate_generalized(&strain)?;
        let mut resultants = [0.0; 6];
        let mut tangent = [[0.0; 6]; 6];
        let mut states = Vec::with_capacity(self.points.len());
        let mut offset = 0;
        for layer in &self.layers {
            let reduction =
                PlaneStressReduction::new(&layer.law, self.reduction).map_err(|source| {
                    SectionError::Reduction {
                        point: offset,
                        source,
                    }
                })?;
            let end = offset + layer.subdivisions * self.quadrature.points_per_subdivision();
            for point_index in offset..end {
                let point = self.points[point_index];
                let response = reduction
                    .evaluate_with_workspace(
                        strain_at(&strain, point.z_m),
                        &old.points[point_index],
                        workspace,
                    )
                    .map_err(|source| SectionError::Reduction {
                        point: point_index,
                        source,
                    })?;
                let w = point.weight_m;
                let wz = w * point.z_m;
                let wzz = wz * point.z_m;
                for i in 0..3 {
                    resultants[i] += w * response.stress_pa[i];
                    resultants[i + 3] += wz * response.stress_pa[i];
                    for j in 0..3 {
                        // Only the strain column changes from tensor to
                        // engineering shear. S^T diag(1,1,2) is the identity,
                        // so physical stress remains the conjugate resultant.
                        let q = response.tangent_pa[i][j] * if j == 2 { 0.5 } else { 1.0 };
                        tangent[i][j] += w * q;
                        tangent[i][j + 3] += wz * q;
                        tangent[i + 3][j] += wz * q;
                        tangent[i + 3][j + 3] += wzz * q;
                    }
                }
                states.push(response.state);
            }
            offset = end;
        }
        if resultants.iter().any(|value| !value.is_finite())
            || tangent.iter().flatten().any(|value| !value.is_finite())
        {
            return Err(SectionError::Invalid(
                "section force, moment or tangent overflow",
            ));
        }
        Ok(LayeredPlaneStressResponse {
            state: LayeredPlaneStressState {
                generalized_strain: strain,
                points: states,
            },
            resultants,
            tangent,
        })
    }
}

fn strain_at(generalized: &[f64; 6], z_m: f64) -> [f64; 3] {
    [
        generalized[0] + z_m * generalized[3],
        generalized[1] + z_m * generalized[4],
        0.5 * (generalized[2] + z_m * generalized[5]),
    ]
}
