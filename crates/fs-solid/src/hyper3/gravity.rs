//! Conservative gravity sampled on the CURRENT P1 material configuration.
//!
//! Cellwise density uses four equal-weight points exact for quadratics. A
//! positive nodal P1 reference density uses eight positive points exact for
//! cubics, retaining density moments instead of replacing them by a cell mean.
//! With fixed reference mass weights `w`, the SAME discrete potential gives
//! `U = sum w phi(x_q)`, `f_ai = sum w N_a g_i(x_q)` and
//! `df_ai/du_bj = sum w N_a N_b dg_i/dx_j`. Thus uniform and affine gravity
//! reproduce the native consistent mass projection, including supported DOFs.
//! General smooth fields require spatial refinement; these fixed rules do
//! not resolve a nearby singular source or certify a quadrature error bound.
//!
//! Gravity is an external conservative body field, not a current-area follower
//! pressure and not a force density multiplied by the changing material volume.
//! The callback must be deterministic and stationary, and must supply mutually
//! consistent potential, acceleration and gradient in the mesh's world axes.

use fs_exec::Cx;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_math::dd::Dd;
use fs_sparse::{Coo, Csr};

use super::{Element, HyperTetError, HyperTetProblem, checkpoint};
use crate::linear3::p1_tet_mass_scale;

/// One conservative field evaluation in fixed world coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GravitySample3 {
    /// Specific potential [J/kg], with one fixed datum for the entire field.
    pub potential_j_kg: f64,
    /// Acceleration [m/s²], equal to minus the gradient of `potential_j_kg`.
    pub acceleration_m_s2: [f64; 3],
    /// Spatial acceleration derivative [1/s²]: entry `(i,j)` is `dg_i/dx_j`.
    /// It must be symmetric to roundoff for this conservative field interface.
    pub acceleration_gradient_s2: [[f64; 3]; 3],
}

/// A borrowed stationary field at a current world position [m].
///
/// Refuse outside the original field's domain; do not clamp a singularity or
/// extrapolate source data. Derivative consistency is the field owner's
/// obligation, not a property established by the Rust callback type.
pub type GravityField3<'a> = dyn Fn([f64; 3]) -> Result<GravitySample3, HyperTetError> + 'a;

/// Specific potential at the second point minus that at the first [J/kg].
///
/// A field-owned analytic difference can cancel its arbitrary datum before
/// rounding. It describes the original field at the two endpoint source
/// snapshots: the same field for stationary gravity, or explicitly retained
/// source poses for reciprocal motion. The field owner admits those snapshots
/// and retains its original parameters; the quadrature retains reference mass.
/// Substituting midpoint force work changes the law for nonlinear fields.
pub type GravityPotentialIncrement3<'a> =
    dyn Fn([f64; 3], [f64; 3]) -> Result<f64, HyperTetError> + 'a;

/// One stationary conservative field and its optional direct potential change.
///
/// Both callbacks must describe the SAME source, parameters, domain and fixed
/// datum. The increment is `phi(after) - phi(before)` from that original law,
/// not midpoint acceleration work. It may cancel an arbitrary datum before
/// rounding, while the field continues to admit and supply endpoint loads.
/// Changed source poses or parameters require a separate source-work owner;
/// they are not conservative static loads represented by this binding.
#[derive(Clone, Copy)]
pub struct StationaryGravity3<'a> {
    field: &'a GravityField3<'a>,
    potential_increment: Option<&'a GravityPotentialIncrement3<'a>>,
}

impl<'a> StationaryGravity3<'a> {
    /// Retain an existing endpoint field without inventing a change formula.
    /// Native equilibrium keeps its endpoint-roundoff safeguard for this case.
    pub const fn from_samples(field: &'a GravityField3<'a>) -> Self {
        Self {
            field,
            potential_increment: None,
        }
    }

    /// Retain the field owner's direct difference alongside its endpoint field.
    ///
    /// The increment callback must admit both original endpoint positions and
    /// propagate domain refusals. A failed increment is never replaced by
    /// endpoint subtraction, a clamped field or a different force law.
    pub const fn with_potential_increment(
        field: &'a GravityField3<'a>,
        potential_increment: &'a GravityPotentialIncrement3<'a>,
    ) -> Self {
        Self {
            field,
            potential_increment: Some(potential_increment),
        }
    }

    /// The original endpoint field, including its domain and derivative data.
    pub const fn field(self) -> &'a GravityField3<'a> {
        self.field
    }
}

/// Gravity integrated with the original tetrahedral reference mass inventory.
#[derive(Debug)]
pub struct TetGravityLoad3 {
    /// Full node-major world forces [N], including prescribed coordinates.
    pub force_n: Vec<f64>,
    /// Derivative of the EXTERNAL force [N/m]; subtract from internal stiffness.
    pub force_tangent_n_m: Csr,
    /// Gravitational potential [J], separate from elastic strain energy.
    pub potential_energy_j: f64,
}

pub(super) struct GravityEvaluation {
    pub force_n: Vec<f64>,
    pub tangent: Option<Csr>,
    pub potential_energy_j: f64,
    /// Direct source-owned change, retained before final reduction rounding.
    /// `None` means no before configuration or no original increment callback.
    pub potential_energy_change_j: Option<Dd>,
}

/// Positive mass density per original reference volume, never current volume.
/// Nodal input retains the complete P1 field and its reference mass moments.
#[derive(Debug, Clone, Copy)]
pub enum ReferenceGravityDensity<'a> {
    /// One fixed density [kg/m³] per original tetrahedron, in connectivity order.
    PerElement(&'a [f64]),
    /// One fixed P1 density [kg/m³] per original mesh node.
    Nodal(&'a [f64]),
}

pub(super) enum ReferenceGravityWeights {
    PerElement(Vec<f64>),
    Nodal(Vec<[f64; 8]>),
}

fn invalid(what: &str) -> HyperTetError {
    HyperTetError::InvalidInput { what: what.into() }
}

fn finite(value: f64) -> Result<f64, HyperTetError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("gravity sample or integral is not representable"))
    }
}

impl GravitySample3 {
    pub(super) fn validate(self) -> Result<Self, HyperTetError> {
        finite(self.potential_j_kg)?;
        for value in self.acceleration_m_s2 {
            finite(value)?;
        }
        let mut scale = f64::MIN_POSITIVE;
        for row in self.acceleration_gradient_s2 {
            for value in row {
                finite(value)?;
                scale = scale.max(value.abs());
            }
        }
        for i in 0..3 {
            for j in 0..i {
                if (self.acceleration_gradient_s2[i][j] - self.acceleration_gradient_s2[j][i]).abs()
                    > 64.0 * f64::EPSILON * scale
                {
                    return Err(invalid(
                        "conservative gravity needs a symmetric acceleration gradient",
                    ));
                }
            }
        }
        Ok(self)
    }
}

/// One quarter of each element's fixed reference mass, shared by all four
/// integration points. The scalar is the native consistent P1 mass kernel.
pub(super) fn reference_weights(
    elements: &[Element],
    reference_density_kg_m3: &[f64],
    cx: &Cx<'_>,
) -> Result<Vec<f64>, HyperTetError> {
    checkpoint(cx)?;
    if reference_density_kg_m3.len() != elements.len() {
        return Err(invalid(
            "gravity needs one reference density per tetrahedron",
        ));
    }
    let mut weights = Vec::with_capacity(elements.len());
    for (element, &density) in elements.iter().zip(reference_density_kg_m3) {
        checkpoint(cx)?;
        if !density.is_finite() || density <= 0.0 {
            return Err(invalid("gravity needs finite positive reference densities"));
        }
        let weight = finite(5.0 * p1_tet_mass_scale(element.volume, density))?;
        if weight <= 0.0 {
            return Err(invalid("gravity reference mass is unrepresentable"));
        }
        weights.push(weight);
    }
    Ok(weights)
}

/// Bind the original material inventory once, before Newton or line search.
/// For nodal rho, the positive cubic rule integrates rho*N_a*N_b exactly.
pub(super) fn bind_reference_weights(
    elements: &[Element],
    node_count: usize,
    density: ReferenceGravityDensity<'_>,
    cx: &Cx<'_>,
) -> Result<ReferenceGravityWeights, HyperTetError> {
    let densities = match density {
        ReferenceGravityDensity::PerElement(densities) => {
            return reference_weights(elements, densities, cx)
                .map(ReferenceGravityWeights::PerElement);
        }
        ReferenceGravityDensity::Nodal(densities) => densities,
    };
    checkpoint(cx)?;
    if node_count == 0
        || densities.len() != node_count
        || densities.iter().any(|rho| !rho.is_finite() || *rho <= 0.0)
    {
        return Err(invalid(
            "gravity needs one finite positive reference density per node",
        ));
    }
    // Keep the original scalar-density rule, sampling and accumulation bits
    // when the new interface receives a globally uniform material inventory.
    if densities
        .iter()
        .all(|rho| rho.to_bits() == densities[0].to_bits())
    {
        let cells = vec![densities[0]; elements.len()];
        return reference_weights(elements, &cells, cx).map(ReferenceGravityWeights::PerElement);
    }
    let mut weights = Vec::with_capacity(elements.len());
    for element in elements {
        checkpoint(cx)?;
        let rho = element.nodes.map(|node| densities[node]);
        let scale = rho.into_iter().fold(0.0_f64, f64::max);
        let normalized = rho.map(|value| value / scale);
        let sum: f64 = normalized.iter().sum();
        let mass = finite(element.volume * (scale * (0.25 * sum)))?;
        let fractions = normalized.map(|value| value / sum);
        let mut local = [0.0; 8];
        for a in 0..4 {
            // Vertices carry V*rho_a/40. Opposite face centroids carry
            // (3*V/40)*sum_{b!=a} rho_b. Summing the OTHER fractions avoids
            // cancellation when one vertex carries almost all the density.
            local[a] = finite(mass * (0.1 * fractions[a]))?;
            let others: f64 = (0..4).filter(|&b| b != a).map(|b| fractions[b]).sum();
            local[4 + a] = finite(mass * (0.3 * others))?;
        }
        if local.iter().any(|weight| *weight <= 0.0) {
            return Err(invalid(
                "graded gravity reference quadrature mass is unrepresentable",
            ));
        }
        weights.push(local);
    }
    Ok(ReferenceGravityWeights::Nodal(weights))
}

// N_i values of the symmetric degree-two tetrahedral rule. Each point's
// weight is V/4: sum w N_i = V/4 and sum w N_i N_j = V(1+delta_ij)/20.
const LARGE: f64 = 0.585_410_196_624_968_5;
const SMALL: f64 = 0.138_196_601_125_010_5;

// A positive cubic tetrahedral rule: four vertices of weight V/40 and four
// opposite face centroids of weight 9V/40. Its moments are
// integral N_i^3 = V/20, N_i^2*N_j = V/60, N_i*N_j*N_k = V/120.
const CUBIC_SHAPES: [[f64; 4]; 8] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
    [0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0],
];

impl ReferenceGravityWeights {
    fn quadrature_size(&self) -> (usize, usize) {
        match self {
            Self::PerElement(weights) => (weights.len(), 4),
            Self::Nodal(weights) => (weights.len(), 8),
        }
    }

    fn point(&self, element: usize, point: usize) -> ([f64; 4], f64) {
        match self {
            Self::PerElement(weights) => (
                std::array::from_fn(|a| if a == point { LARGE } else { SMALL }),
                weights[element],
            ),
            Self::Nodal(weights) => (CUBIC_SHAPES[point], weights[element][point]),
        }
    }
}

#[allow(clippy::needless_range_loop)]
fn position_at(
    nodes_m: &[[f64; 3]],
    element: &Element,
    displacement_m: &[f64],
    shape: [f64; 4],
) -> Result<[f64; 3], HyperTetError> {
    // Relative interpolation preserves common translation in the geometry
    // instead of adding four large absolute coordinates. Energy increments
    // and force evaluations must sample exactly the same material points.
    let origin = element.nodes[0];
    let mut position = [0.0; 3];
    for i in 0..3 {
        position[i] = finite(nodes_m[origin][i] + displacement_m[3 * origin + i])?;
        for a in 1..4 {
            let node = element.nodes[a];
            let relative = finite(
                (nodes_m[node][i] - nodes_m[origin][i])
                    + (displacement_m[3 * node + i] - displacement_m[3 * origin + i]),
            )?;
            position[i] = finite(shape[a].mul_add(relative, position[i]))?;
        }
    }
    Ok(position)
}

/// Integrate field-owned differences before rounding or summing absolute
/// potentials. Fixed material weights retain exactly the original inventory.
pub(super) fn potential_increment(
    nodes_m: &[[f64; 3]],
    elements: &[Element],
    mass_weights: &ReferenceGravityWeights,
    before_m: &[f64],
    after_m: &[f64],
    field_increment: &GravityPotentialIncrement3<'_>,
    cx: &Cx<'_>,
) -> Result<f64, HyperTetError> {
    checkpoint(cx)?;
    let dofs = 3 * nodes_m.len();
    for displacement in [before_m, after_m] {
        if displacement.len() != dofs || displacement.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "gravity needs one finite xyz displacement per reference node",
            ));
        }
    }
    let (element_count, point_count) = mass_weights.quadrature_size();
    if element_count != elements.len() {
        return Err(invalid(
            "gravity reference weights do not cover the original mesh",
        ));
    }
    let mut change = Dd::from_f64(0.0);
    for (element_index, element) in elements.iter().enumerate() {
        checkpoint(cx)?;
        for point in 0..point_count {
            let (shape, weight) = mass_weights.point(element_index, point);
            let before = position_at(nodes_m, element, before_m, shape)?;
            let after = position_at(nodes_m, element, after_m, shape)?;
            let specific_change = finite(field_increment(before, after)?)?;
            change = change + Dd::from_f64(weight) * Dd::from_f64(specific_change);
        }
    }
    checkpoint(cx)?;
    finite(change.to_f64())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate(
    nodes_m: &[[f64; 3]],
    elements: &[Element],
    mass_weights: &ReferenceGravityWeights,
    displacement_m: &[f64],
    field: &GravityField3<'_>,
    with_tangent: bool,
    cx: &Cx<'_>,
) -> Result<GravityEvaluation, HyperTetError> {
    evaluate_with_increment(
        nodes_m,
        elements,
        mass_weights,
        displacement_m,
        StationaryGravity3::from_samples(field),
        None,
        with_tangent,
        cx,
    )
}

/// Evaluate the trial load and, when available, its direct source work in the
/// SAME reference-mass quadrature. The caller has already admitted the before
/// configuration through the original field. Trial endpoint sample admission
/// always precedes its increment; no extra before sample or force vector is
/// needed. A source without a direct increment keeps the original sample count
/// and endpoint arithmetic exactly.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub(super) fn evaluate_with_increment(
    nodes_m: &[[f64; 3]],
    elements: &[Element],
    mass_weights: &ReferenceGravityWeights,
    displacement_m: &[f64],
    source: StationaryGravity3<'_>,
    before_m: Option<&[f64]>,
    with_tangent: bool,
    cx: &Cx<'_>,
) -> Result<GravityEvaluation, HyperTetError> {
    checkpoint(cx)?;
    let dofs = 3 * nodes_m.len();
    for displacement in std::iter::once(displacement_m).chain(before_m) {
        if displacement.len() != dofs || displacement.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "gravity needs one finite xyz displacement per reference node",
            ));
        }
    }
    let (element_count, point_count) = mass_weights.quadrature_size();
    if element_count != elements.len() {
        return Err(invalid(
            "gravity reference weights do not cover the original mesh",
        ));
    }
    let mut force_n = vec![0.0; dofs];
    let mut potential_energy_j = 0.0;
    let interval = before_m.zip(source.potential_increment);
    let mut potential_energy_change_j = interval.map(|_| Dd::from_f64(0.0));
    let mut entries = with_tangent.then(|| Coo::new(dofs, dofs));
    for (element_index, element) in elements.iter().enumerate() {
        checkpoint(cx)?;
        let mut local = [[0.0; 12]; 12];
        for point in 0..point_count {
            let (shape, weight) = mass_weights.point(element_index, point);
            let position = position_at(nodes_m, element, displacement_m, shape)?;
            let sample = (source.field)(position)?.validate()?;
            potential_energy_j = finite(weight.mul_add(sample.potential_j_kg, potential_energy_j))?;
            if let Some((before, increment)) = interval {
                let before_position = position_at(nodes_m, element, before, shape)?;
                let specific_change = finite(increment(before_position, position)?)?;
                if let Some(change) = &mut potential_energy_change_j {
                    *change = *change + Dd::from_f64(weight) * Dd::from_f64(specific_change);
                }
            }
            for a in 0..4 {
                let share = weight * shape[a];
                for i in 0..3 {
                    let dof = 3 * element.nodes[a] + i;
                    force_n[dof] =
                        finite(share.mul_add(sample.acceleration_m_s2[i], force_n[dof]))?;
                }
                if with_tangent {
                    for b in 0..4 {
                        let product = share * shape[b];
                        for i in 0..3 {
                            for j in 0..3 {
                                local[3 * a + i][3 * b + j] = finite(product.mul_add(
                                    sample.acceleration_gradient_s2[i][j],
                                    local[3 * a + i][3 * b + j],
                                ))?;
                            }
                        }
                    }
                }
            }
        }
        if let Some(entries) = &mut entries {
            for a in 0..4 {
                for b in 0..4 {
                    for i in 0..3 {
                        for j in 0..3 {
                            let value = local[3 * a + i][3 * b + j];
                            if value != 0.0 {
                                entries.push(
                                    3 * element.nodes[a] + i,
                                    3 * element.nodes[b] + j,
                                    value,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    let tangent = entries.map(|entries| entries.assemble());
    if let Some(matrix) = &tangent {
        for row in 0..dofs {
            if row % 4096 == 0 {
                checkpoint(cx)?;
            }
            for &value in matrix.row(row).1 {
                finite(value)?;
            }
        }
    }
    if let Some(change) = potential_energy_change_j {
        finite(change.to_f64())?;
    }
    checkpoint(cx)?;
    Ok(GravityEvaluation {
        force_n,
        tangent,
        potential_energy_j,
        potential_energy_change_j,
    })
}

impl<L: HyperelasticMaterial> HyperTetProblem<'_, L> {
    /// Sample current-configuration gravity with the native reference mass.
    ///
    /// Densities are exactly the element-ordered reference densities used by
    /// [`Self::reference_mass_matrix`]; they must remain fixed during motion.
    /// The field has no implicit Earth preset. This operation assembles a load;
    /// native equilibrium separately admits the current constitutive state.
    ///
    /// # Errors
    /// Invalid reference problem, density/displacement/field, overflow or cancellation.
    pub fn gravity_load(
        &self,
        reference_density_kg_m3: &[f64],
        displacement_m: &[f64],
        field: &GravityField3<'_>,
        cx: &Cx<'_>,
    ) -> Result<TetGravityLoad3, HyperTetError> {
        self.prepare(cx)?
            .gravity_load(reference_density_kg_m3, displacement_m, field, cx)
    }

    /// Current-configuration conservative gravity for a positive P1 reference
    /// density, supplied once per original mesh node in kg/m³.
    ///
    /// The cubic rule retains the graded inventory's center of mass and full
    /// consistent inertia moments. During admitted positive-J deformation its fixed reference
    /// weights represent rho_current = rho_reference / det F. No current-volume
    /// multiplier or cell-mean density substitution enters force or potential.
    /// Uniform nodal density retains the original four-point rule exactly.
    /// Native equilibrium separately admits the current constitutive state.
    ///
    /// # Errors
    /// Invalid reference problem, density/displacement/field, overflow or cancellation.
    pub fn gravity_load_with_nodal_density(
        &self,
        reference_nodal_density_kg_m3: &[f64],
        displacement_m: &[f64],
        field: &GravityField3<'_>,
        cx: &Cx<'_>,
    ) -> Result<TetGravityLoad3, HyperTetError> {
        self.prepare(cx)?.gravity_load_with_nodal_density(
            reference_nodal_density_kg_m3,
            displacement_m,
            field,
            cx,
        )
    }
}

#[cfg(test)]
mod nodal_density_tests {
    use super::*;
    use crate::inertia::ReferenceMass;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    const NODES: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    const ELEMENT: Element = Element {
        nodes: [0, 1, 2, 3],
        gradients: [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        volume: 1.0 / 6.0,
    };
    const DENSITY: [f64; 4] = [3.0, 15.0, 3.0, 3.0];
    const SHARES: [f64; 4] = [0.225, 0.325, 0.225, 0.225];
    const GRADIENT: [[f64; 3]; 3] = [[-0.4, 0.1, 0.0], [0.1, -0.7, 0.05], [0.0, 0.05, -1.2]];
    const CONSTANT: [f64; 3] = [0.3, -0.5, -2.0];

    fn with_cx(run: impl FnOnce(&Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            run(&Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ));
        });
    }

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 4e-13 * expected.abs().max(1.0),
            "{actual:.17e} != {expected:.17e}"
        );
    }

    fn affine(x: [f64; 3]) -> Result<GravitySample3, HyperTetError> {
        let acceleration = std::array::from_fn(|i| {
            CONSTANT[i] + (0..3).map(|j| GRADIENT[i][j] * x[j]).sum::<f64>()
        });
        Ok(GravitySample3 {
            potential_j_kg: -(0..3)
                .map(|i| 0.5 * x[i] * (acceleration[i] + CONSTANT[i]))
                .sum::<f64>(),
            acceleration_m_s2: acceleration,
            acceleration_gradient_s2: GRADIENT,
        })
    }

    fn displacement() -> [f64; 12] {
        std::array::from_fn(|d| {
            let x = NODES[d / 3];
            [
                0.7 + 0.7 * x[0] + 0.2 * x[1] - 0.1 * x[2],
                -0.2 + 0.1 * x[0] - 0.2 * x[1] + 0.3 * x[2],
                0.3 + 0.15 * x[0] - 0.05 * x[1] + 0.2 * x[2],
            ][d % 3]
        })
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    #[test]
    fn cubic_gravity_preserves_native_graded_mass_force_tangent_energy_and_moment() {
        // Independent simplex monomial integral, normalized by volume:
        // 6*product(p_i!)/(sum(p_i)+3)!. Check ALL 35 monomials through degree 3.
        let factorial = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0];
        for a in 0..=3 {
            for b in 0..=3 - a {
                for c in 0..=3 - a - b {
                    for d in 0..=3 - a - b - c {
                        let powers = [a, b, c, d];
                        let exact = 6.0 * powers.map(|p| factorial[p]).into_iter().product::<f64>()
                            / factorial[3 + a + b + c + d];
                        let rule: f64 = CUBIC_SHAPES
                            .iter()
                            .enumerate()
                            .map(|(q, shape)| {
                                let value: f64 =
                                    (0..4).map(|i| shape[i].powi(powers[i] as i32)).product();
                                value * if q < 4 { 1.0 / 40.0 } else { 9.0 / 40.0 }
                            })
                            .sum();
                        close(rule, exact);
                    }
                }
            }
        }
        with_cx(|cx| {
            let weights =
                bind_reference_weights(&[ELEMENT], 4, ReferenceGravityDensity::Nodal(&DENSITY), cx)
                    .unwrap();
            let mass = ReferenceMass::from_nodal_density(
                [(ELEMENT.nodes, ELEMENT.volume)].into_iter(),
                &DENSITY,
                12,
                cx,
            )
            .unwrap()
            .0;
            let u = displacement();
            let positions: [[f64; 3]; 4] =
                std::array::from_fn(|a| std::array::from_fn(|i| NODES[a][i] + u[3 * a + i]));
            let acceleration = positions.map(|x| affine(x).unwrap().acceleration_m_s2);
            let expected = mass.action(acceleration.as_flattened(), cx).unwrap();
            let load = evaluate(&NODES, &[ELEMENT], &weights, &u, &affine, true, cx).unwrap();
            for (actual, expected) in load.force_n.iter().zip(&expected) {
                close(*actual, *expected);
            }
            let tangent = load.tangent.as_ref().unwrap();
            for column in 0..12 {
                let mut basis = [0.0; 12];
                basis[column] = 1.0;
                let column_mass = mass.action(&basis, cx).unwrap();
                for row in 0..12 {
                    close(
                        tangent.get(row, column),
                        column_mass[3 * (row / 3) + column % 3] * GRADIENT[row % 3][column % 3],
                    );
                }
            }
            let potential: f64 = (0..12)
                .map(|d| {
                    -0.5 * positions[d / 3][d % 3] * (expected[d] + SHARES[d / 3] * CONSTANT[d % 3])
                })
                .sum();
            close(load.potential_energy_j, potential);
            close(mass.reference_mass_kg(), 1.0);

            // A changed current Jacobian changes rho_current, not inventory.
            // Uniform acceleration retains total weight and the MOVED graded COM.
            let uniform = |x: [f64; 3]| {
                Ok(GravitySample3 {
                    potential_j_kg: -(0..3).map(|i| CONSTANT[i] * x[i]).sum::<f64>(),
                    acceleration_m_s2: CONSTANT,
                    acceleration_gradient_s2: [[0.0; 3]; 3],
                })
            };
            let undeformed = evaluate(
                &NODES,
                &[ELEMENT],
                &weights,
                &[0.0; 12],
                &uniform,
                false,
                cx,
            )
            .unwrap();
            let deformed = evaluate(&NODES, &[ELEMENT], &weights, &u, &uniform, false, cx).unwrap();
            assert_eq!(undeformed.force_n, deformed.force_n);
            let com = std::array::from_fn(|i| (0..4).map(|a| SHARES[a] * positions[a][i]).sum());
            let expected_moment = cross(com, CONSTANT);
            let mut moment = [0.0; 3];
            for a in 0..4 {
                let force = std::array::from_fn(|i| deformed.force_n[3 * a + i]);
                for i in 0..3 {
                    close(force[i], SHARES[a] * CONSTANT[i]);
                }
                let arm = cross(positions[a], force);
                for i in 0..3 {
                    moment[i] += arm[i];
                }
            }
            for i in 0..3 {
                close(moment[i], expected_moment[i]);
            }
            let centroid = std::array::from_fn(|i| positions.iter().map(|x| 0.25 * x[i]).sum());
            let averaged_moment = cross(centroid, CONSTANT);
            assert!((moment[1] - averaged_moment[1]).abs() > 0.01);
        });
    }

    #[test]
    fn graded_nonlinear_gravity_shares_one_discrete_potential_force_and_tangent() {
        with_cx(|cx| {
            let weights =
                bind_reference_weights(&[ELEMENT], 4, ReferenceGravityDensity::Nodal(&DENSITY), cx)
                    .unwrap();
            let field = |x: [f64; 3]| {
                Ok(GravitySample3 {
                    potential_j_kg: x.into_iter().map(|v| 0.25 * v.powi(4)).sum(),
                    acceleration_m_s2: x.map(|v| -v.powi(3)),
                    acceleration_gradient_s2: std::array::from_fn(|i| {
                        std::array::from_fn(|j| if i == j { -3.0 * x[i] * x[i] } else { 0.0 })
                    }),
                })
            };
            let u = displacement();
            let load = evaluate(&NODES, &[ELEMENT], &weights, &u, &field, true, cx).unwrap();
            let h = 2e-6;
            for column in 0..12 {
                let mut plus = u;
                let mut minus = u;
                plus[column] += h;
                minus[column] -= h;
                let up = evaluate(&NODES, &[ELEMENT], &weights, &plus, &field, false, cx).unwrap();
                let down =
                    evaluate(&NODES, &[ELEMENT], &weights, &minus, &field, false, cx).unwrap();
                let force = -(up.potential_energy_j - down.potential_energy_j) / (2.0 * h);
                assert!((force - load.force_n[column]).abs() < 2e-9);
                for row in 0..12 {
                    let tangent = (up.force_n[row] - down.force_n[row]) / (2.0 * h);
                    assert!(
                        (tangent - load.tangent.as_ref().unwrap().get(row, column)).abs() < 2e-9
                    );
                }
            }
        });
    }

    #[test]
    fn uniform_nodal_gravity_keeps_scalar_bits_and_invalid_inventory_refuses() {
        with_cx(|cx| {
            let scalar = bind_reference_weights(
                &[ELEMENT],
                4,
                ReferenceGravityDensity::PerElement(&[7.3]),
                cx,
            )
            .unwrap();
            let nodal = bind_reference_weights(
                &[ELEMENT],
                4,
                ReferenceGravityDensity::Nodal(&[7.3; 4]),
                cx,
            )
            .unwrap();
            assert!(matches!(&nodal, ReferenceGravityWeights::PerElement(_)));
            let u = displacement();
            let original = evaluate(&NODES, &[ELEMENT], &scalar, &u, &affine, true, cx).unwrap();
            let uniform = evaluate(&NODES, &[ELEMENT], &nodal, &u, &affine, true, cx).unwrap();
            assert_eq!(original.force_n, uniform.force_n);
            assert_eq!(original.tangent, uniform.tangent);
            assert_eq!(
                original.potential_energy_j.to_bits(),
                uniform.potential_energy_j.to_bits()
            );
            for rho in [
                &[1.0, 2.0, 3.0][..],
                &[1.0, 0.0, 2.0, 3.0],
                &[1.0, f64::NAN, 2.0, 3.0],
            ] {
                assert!(
                    bind_reference_weights(&[ELEMENT], 4, ReferenceGravityDensity::Nodal(rho), cx)
                        .is_err()
                );
            }
            let outside = |_| {
                Err(HyperTetError::InvalidInput {
                    what: "outside original field".into(),
                })
            };
            assert!(evaluate(&NODES, &[ELEMENT], &nodal, &u, &outside, true, cx).is_err());
            let cancelled = |_| Err(HyperTetError::Cancelled);
            assert!(matches!(
                evaluate(&NODES, &[ELEMENT], &nodal, &u, &cancelled, true, cx),
                Err(HyperTetError::Cancelled)
            ));
        });
    }

    #[test]
    fn potential_increments_preserve_small_work_and_graded_mass_moments() {
        with_cx(|cx| {
            let before = [0.0; 12];
            let after = displacement();
            let field = |x: [f64; 3]| {
                Ok(GravitySample3 {
                    potential_j_kg: 1e20 + 8.0 * x[1],
                    acceleration_m_s2: [0.0, -8.0, 0.0],
                    acceleration_gradient_s2: [[0.0; 3]; 3],
                })
            };
            for (density, shares) in [
                (ReferenceGravityDensity::PerElement(&[6.0]), [0.25; 4]),
                (ReferenceGravityDensity::Nodal(&DENSITY), SHARES),
            ] {
                let weights = bind_reference_weights(&[ELEMENT], 4, density, cx).unwrap();
                let a = evaluate(&NODES, &[ELEMENT], &weights, &before, &field, false, cx).unwrap();
                let b = evaluate(&NODES, &[ELEMENT], &weights, &after, &field, false, cx).unwrap();
                // The datum has erased the motion from the rounded diagnostics.
                assert_eq!(a.potential_energy_j, b.potential_energy_j);
                let increment = |a: [f64; 3], b: [f64; 3]| Ok(8.0 * (b[1] - a[1]));
                let change = potential_increment(
                    &NODES,
                    &[ELEMENT],
                    &weights,
                    &before,
                    &after,
                    &increment,
                    cx,
                )
                .unwrap();
                // Independent exact P1 rho*N moments, not a quadrature oracle.
                let expected = 8.0 * (0..4).map(|a| shares[a] * after[3 * a + 1]).sum::<f64>();
                close(change, expected);
                let reverse = potential_increment(
                    &NODES,
                    &[ELEMENT],
                    &weights,
                    &after,
                    &before,
                    &increment,
                    cx,
                )
                .unwrap();
                close(reverse, -expected);
                assert!(change.abs() > 0.1);

                let fused = evaluate_with_increment(
                    &NODES,
                    &[ELEMENT],
                    &weights,
                    &after,
                    StationaryGravity3::with_potential_increment(&field, &increment),
                    Some(&before),
                    false,
                    cx,
                )
                .unwrap();
                close(fused.potential_energy_change_j.unwrap().to_f64(), expected);
                assert_eq!(fused.force_n, b.force_n);
                assert_eq!(fused.potential_energy_j, b.potential_energy_j);
            }
        });
    }

    #[test]
    fn fused_gravity_work_keeps_one_sample_per_original_mass_point() {
        with_cx(|cx| {
            let sample_count = std::cell::Cell::new(0);
            let increment_count = std::cell::Cell::new(0);
            let field = |x| {
                sample_count.set(sample_count.get() + 1);
                affine(x)
            };
            // This exact polarization belongs to the original QUADRATIC
            // potential; the nonlinear regression below uses its own quartic.
            let increment = |a: [f64; 3], b: [f64; 3]| {
                increment_count.set(increment_count.get() + 1);
                let ga = affine(a)?.acceleration_m_s2;
                let gb = affine(b)?.acceleration_m_s2;
                Ok(-(0..3)
                    .map(|i| (b[i] - a[i]) * (0.5 * ga[i] + 0.5 * gb[i]))
                    .sum::<f64>())
            };
            let before = [0.0; 12];
            let after = displacement();
            for density in [
                ReferenceGravityDensity::PerElement(&[6.0]),
                ReferenceGravityDensity::Nodal(&[6.0; 4]),
                ReferenceGravityDensity::Nodal(&DENSITY),
            ] {
                let weights = bind_reference_weights(&[ELEMENT], 4, density, cx).unwrap();
                let points = weights.quadrature_size().1;
                let original =
                    evaluate(&NODES, &[ELEMENT], &weights, &after, &affine, true, cx).unwrap();
                for (source, previous, expected_increments) in [
                    (
                        StationaryGravity3::from_samples(&field),
                        Some(before.as_slice()),
                        0,
                    ),
                    (
                        StationaryGravity3::with_potential_increment(&field, &increment),
                        None,
                        0,
                    ),
                    (
                        StationaryGravity3::with_potential_increment(&field, &increment),
                        Some(before.as_slice()),
                        points,
                    ),
                ] {
                    sample_count.set(0);
                    increment_count.set(0);
                    let evaluated = evaluate_with_increment(
                        &NODES,
                        &[ELEMENT],
                        &weights,
                        &after,
                        source,
                        previous,
                        true,
                        cx,
                    )
                    .unwrap();
                    assert_eq!(sample_count.get(), points);
                    assert_eq!(increment_count.get(), expected_increments);
                    assert_eq!(
                        evaluated.potential_energy_change_j.is_some(),
                        expected_increments > 0
                    );
                    assert_eq!(evaluated.force_n, original.force_n);
                    assert_eq!(evaluated.tangent, original.tangent);
                    assert_eq!(
                        evaluated.potential_energy_j.to_bits(),
                        original.potential_energy_j.to_bits()
                    );
                }
            }
        });
    }

    #[test]
    fn potential_increment_retains_nonlinear_storage_and_source_refusals() {
        with_cx(|cx| {
            let weights =
                bind_reference_weights(&[ELEMENT], 4, ReferenceGravityDensity::Nodal(&DENSITY), cx)
                    .unwrap();
            let before = [0.0; 12];
            let after = displacement();
            let field = |x: [f64; 3]| {
                Ok(GravitySample3 {
                    potential_j_kg: x.into_iter().map(|v| 0.25 * v.powi(4)).sum(),
                    acceleration_m_s2: x.map(|v| -v.powi(3)),
                    acceleration_gradient_s2: std::array::from_fn(|i| {
                        std::array::from_fn(|j| if i == j { -3.0 * x[i] * x[i] } else { 0.0 })
                    }),
                })
            };
            let increment = |a: [f64; 3], b: [f64; 3]| {
                Ok((0..3)
                    .map(|i| 0.25 * (b[i] - a[i]) * (b[i] + a[i]) * (b[i] * b[i] + a[i] * a[i]))
                    .sum())
            };
            let change = potential_increment(
                &NODES,
                &[ELEMENT],
                &weights,
                &before,
                &after,
                &increment,
                cx,
            )
            .unwrap();
            let first = evaluate(&NODES, &[ELEMENT], &weights, &before, &field, false, cx).unwrap();
            let last = evaluate(&NODES, &[ELEMENT], &weights, &after, &field, false, cx).unwrap();
            close(change, last.potential_energy_j - first.potential_energy_j);
            let shifted_field = |x| {
                let mut sample = field(x)?;
                sample.potential_j_kg += 1e20;
                Ok(sample)
            };
            let fused = evaluate_with_increment(
                &NODES,
                &[ELEMENT],
                &weights,
                &after,
                StationaryGravity3::with_potential_increment(&shifted_field, &increment),
                Some(&before),
                true,
                cx,
            )
            .unwrap();
            close(fused.potential_energy_change_j.unwrap().to_f64(), change);
            assert_eq!(fused.force_n, last.force_n);
            // Domain/sample validation is still mandatory, even if an
            // increment callback would produce a finite number.
            let outside = |_| Err(invalid("outside original field domain"));
            assert!(
                evaluate_with_increment(
                    &NODES,
                    &[ELEMENT],
                    &weights,
                    &after,
                    StationaryGravity3::with_potential_increment(&outside, &increment),
                    Some(&before),
                    false,
                    cx,
                )
                .is_err()
            );
            for refused in [f64::NAN, f64::INFINITY] {
                assert!(
                    evaluate_with_increment(
                        &NODES,
                        &[ELEMENT],
                        &weights,
                        &after,
                        StationaryGravity3::with_potential_increment(&field, &|_, _| Ok(refused)),
                        Some(&before),
                        false,
                        cx,
                    )
                    .is_err()
                );
            }
            assert!(matches!(
                evaluate_with_increment(
                    &NODES,
                    &[ELEMENT],
                    &weights,
                    &after,
                    StationaryGravity3::with_potential_increment(&field, &|_, _| Err(
                        HyperTetError::Cancelled
                    )),
                    Some(&before),
                    false,
                    cx,
                ),
                Err(HyperTetError::Cancelled)
            ));
            assert!(
                potential_increment(
                    &NODES,
                    &[ELEMENT],
                    &weights,
                    &before[..11],
                    &after,
                    &increment,
                    cx,
                )
                .is_err()
            );
            assert!(
                potential_increment(
                    &NODES,
                    &[ELEMENT],
                    &weights,
                    &before,
                    &after,
                    &|_, _| Ok(f64::NAN),
                    cx,
                )
                .is_err()
            );
            assert!(matches!(
                potential_increment(
                    &NODES,
                    &[ELEMENT],
                    &weights,
                    &before,
                    &after,
                    &|_, _| Err(HyperTetError::Cancelled),
                    cx,
                ),
                Err(HyperTetError::Cancelled)
            ));
        });
    }
}
