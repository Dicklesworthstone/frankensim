//! Insulated-slab thermoelastic factors on the original MITC trial space.
//!
//! A homogeneous symmetric slab has two different thermal mechanisms. Its
//! uniform through-thickness temperature mode contributes reversible membrane
//! stiffness; odd diffusion modes contribute bending storage and dissipation.
//! This module projects their ORIGINAL constitutive square factors separately,
//! without forming a dense material or nodal Gram matrix first. The caller owns
//! the thermal source, finite positive diffusion spectrum and time integration.
//!
//! Directors, physical rotations, membrane strain, curvature and elastic
//! director-bubble condensation all come from the existing MITC element owner.
//! The same positive seven-point quadrature integrates the squared quadratic
//! curvature exactly in the represented trial space. No transverse-shear or
//! numerical-drilling loss is invented. Static condensation remains the
//! original fixed elastic trial-space approximation, not a dynamic elimination
//! of additional thermal/director coordinates.

use super::{MitcBudget, NodalDirectors, PlateError, ShearShellSection, ShellMesh, element};
use fs_math::det;

/// Complete retained mode count; the caller's work and storage ceilings still
/// apply. No accepted mode or small constitutive channel is truncated.
pub const MAX_THERMOELASTIC_MODES: usize = super::viscosity::MAX_VISCOSITY_MODES;

/// Original rank-one thermal resultants in the facet's engineering-strain frame.
///
/// For a homogeneous symmetric slab, obtain bending from the shared
/// SlabThermoelasticBending::new_coupled(...).relaxation_factor() and set
/// membrane = sqrt(12) / h * bending using the SAME physical thickness.
/// The source owner must supply isothermal stiffness and consistent fixed-state
/// thermal data; these numeric factors do not establish material authority.
/// A coupled laminate or independently resolved thermal field needs its own
/// original thermal model and must not be silently replaced by this slab law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellThermoelasticFactors {
    /// Square factor of the adiabatic membrane increment [sqrt(N/m)].
    /// Its zero thermal-diffusion mode stores energy but produces no heat.
    pub membrane: [f64; 3],
    /// Square factor of the full adiabatic bending increment [sqrt(N m)].
    /// Positive dimensionless thermal-spectrum weights multiply its Gram.
    pub bending: [f64; 3],
}

/// Original coordinate maps, before any optional orthogonal row compression.
/// Both arrays are row-major row_count * mode_count; rows follow selected
/// facets in mesh order and the owner's seven-point rule. Exact zero rows pad
/// the end when there are fewer quadrature points than retained modes.
///
/// For original mass-normalized physical modes these are modal stiffness
/// square factors [1/s]. Their COMPLETE Grams, including every off-diagonal
/// entry, are the corresponding thermal stiffness increments [1/s²].
/// Arbitrary physical test vectors are also admitted without rescaling them.
#[derive(Clone, Debug)]
pub struct ProjectedThermoelasticRows {
    /// Row-major conservative membrane factors; no relaxation weight is applied.
    pub membrane: Vec<f64>,
    /// Row-major bending factors before the source's thermal spectrum weights.
    pub bending: Vec<f64>,
    /// Seven rows per selected facet, padded to at least mode_count.
    pub row_count: usize,
    /// Complete supplied physical basis size, with no renormalization.
    pub mode_count: usize,
    /// Original element for each quadrature row; None identifies zero padding.
    /// Different regional thermal rates can share one physical FE projection.
    pub row_elements: Vec<Option<usize>>,
}

/// Conservative cold-work allowance for source projection, excluding downstream
/// row compression and spectrum expansion. Includes whole-input admission and
/// the original elastic condensation on every selected facet. These are scalar
/// admission units, not measured timings or a continuum error bound.
pub fn projection_terms(
    nodes: usize,
    triangles: usize,
    selected_triangles: usize,
    modes: usize,
) -> Option<usize> {
    if selected_triangles > triangles { return None; }
    let selected_work = selected_triangles.checked_mul(
        1_048_576usize.checked_add(8192usize.checked_mul(modes)?)?,
    )?;
    // The original director admission compares incident facet pairs. Include
    // its worst-case valence cost, even when only one region is selected.
    triangles.checked_mul(triangles)?.checked_mul(9)?
        .checked_add(triangles.checked_mul(1024)?)?
        .checked_add(nodes.checked_mul(12)?)?
        .checked_add(6usize.checked_mul(nodes)?.checked_mul(modes)?)?
        .checked_add(selected_work)
}

fn bad(what: &'static str) -> PlateError { PlateError::BadSection { what } }

/// Project both thermal mechanisms on the COMPLETE supplied physical basis.
///
/// Modes must be expanded by the original MitcShellModel::expand, retaining
/// the support frames, mode signs and mass normalization. Sections must be
/// the exact elastic sections used by that model: their original bubble map is
/// reused. Each Some selects one facet; None permits separate regional
/// projection when different source materials have different thermal rates.
/// An explicit zero factor remains exactly zero, without a damping floor.
///
/// max_factor_entries bounds the COMBINED number of f64 entries in both
/// output arrays, including QR padding. Geometry, mode, work and storage
/// admission precedes output allocation. Cancellation is checked before work,
/// at each facet/quadrature/mode and before publication. Inputs are immutable.
///
/// Arithmetic is deterministic ordinary floating point with compensated short
/// source contractions. This is neither an eigensolve, a source-applicability
/// check nor a spatial/roundoff certificate. The caller retains the spectrum's
/// separate finite-band bound and must not infer a forcing bandlimit from it.
#[allow(clippy::too_many_arguments)]
pub fn project_thermoelastic_factors(
    mesh: &ShellMesh,
    directors: &NodalDirectors,
    sections: &[ShearShellSection],
    factors: &[Option<ShellThermoelasticFactors>],
    modes: &[Vec<f64>],
    budget: MitcBudget,
    max_projection_terms: usize,
    max_factor_entries: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ProjectedThermoelasticRows, PlateError> {
    if cancelled() { return Err(PlateError::Cancelled); }
    let n = modes.len();
    let nodes = mesh.nodes.len();
    let triangles = mesh.tris.len();
    let dofs = nodes.checked_mul(6)
        .ok_or_else(|| bad("thermal shell physical coordinate extent overflow"))?;
    if nodes > budget.max_nodes || triangles > budget.max_triangles
        || n == 0 || n > MAX_THERMOELASTIC_MODES
        || sections.len() != triangles || factors.len() != triangles {
        return Err(bad("thermal shell requires complete bounded sections, factors and physical modes"));
    }
    let selected = factors.iter().filter(|f| f.is_some()).count();
    let row_count = selected.checked_mul(7)
        .ok_or_else(|| bad("thermal shell quadrature extent overflow"))?.max(n);
    let entries = row_count.checked_mul(n)
        .ok_or_else(|| bad("thermal shell factor extent overflow"))?;
    if entries.checked_mul(2).is_none_or(|count| count > max_factor_entries)
        || projection_terms(nodes, triangles, selected, n)
            .is_none_or(|work| work > max_projection_terms) {
        return Err(bad("complete thermal shell projection exceeds its explicit work or factor allowance"));
    }
    if modes.iter().any(|p| p.len() != dofs || p.iter().any(|x| !x.is_finite()))
        || factors.iter().flatten().any(|f|
            f.membrane.iter().chain(&f.bending).any(|v| !v.is_finite())) {
        return Err(bad("thermal shell requires finite complete physical vectors and source factors"));
    }
    directors.validate(mesh)?;
    for section in sections {
        if cancelled() { return Err(PlateError::Cancelled); }
        section.validate()?;
    }
    let mut out = ProjectedThermoelasticRows {
        membrane: vec![0.0; entries], bending: vec![0.0; entries],
        row_count, mode_count: n, row_elements: vec![None; row_count],
    };
    let rule = element::stiffness_rule();
    let mut row = 0;
    for (e, tri) in mesh.tris.iter().enumerate() {
        if cancelled() { return Err(PlateError::Cancelled); }
        let Some(factor) = factors[e] else { continue; };
        out.row_elements[row..row + rule.len()].fill(Some(e));
        if factor.membrane.iter().chain(&factor.bending).all(|v| *v == 0.0) {
            row += rule.len();
            continue;
        }
        let geometry = mesh.facet(e)?;
        let field = element::PreparedStrainField::new(
            mesh, tri, geometry, tri.map(|i| directors.vectors[i]),
            &sections[e], &mut cancelled,
        )?;
        // Read the ORIGINAL elastic bubble coefficients through its unit
        // expansion. Compensated contractions below preserve weak represented
        // motion when larger physical coordinates cancel; no new condensation.
        let bubble: [[f64; 2]; 18] = core::array::from_fn(|j| {
            let mut unit = [0.0; 18];
            unit[j] = 1.0;
            let expanded = field.expand(&unit);
            [expanded[18], expanded[19]]
        });
        let mut local = Vec::with_capacity(n);
        for mode in modes {
            if cancelled() { return Err(PlateError::Cancelled); }
            let nodal: [f64; 18] =
                core::array::from_fn(|j| mode[6 * tri[j / 6] + j % 6]);
            let mut expanded = [0.0; 20];
            expanded[..18].copy_from_slice(&nodal);
            for k in 0..2 {
                expanded[18 + k] = compensated_sum(
                    (0..18).map(|j| bubble[j][k] * nodal[j]),
                );
            }
            if expanded.iter().any(|v| !v.is_finite()) {
                return Err(bad("thermal shell elastic expansion is unrepresentable"));
            }
            local.push(expanded);
        }
        for &(lambda, weight) in &rule {
            if cancelled() { return Err(PlateError::Cancelled); }
            let scale = det::sqrt(geometry.area_m2) * det::sqrt(weight);
            if !scale.is_finite() || scale <= 0.0 {
                return Err(bad("thermal shell requires positive finite original quadrature"));
            }
            // Evaluate each original kinematic basis once, then form each
            // actual strain BEFORE contracting its constitutive factor.
            let basis: [([f64; 3], [f64; 3], [f64; 2]); 20] =
                core::array::from_fn(|j| {
                    let mut unit = [0.0; 20];
                    unit[j] = 1.0;
                    field.eval(&unit, lambda[1], lambda[2])
                });
            for (column, displacement) in local.iter().enumerate() {
                if cancelled() { return Err(PlateError::Cancelled); }
                let membrane: [f64; 3] = core::array::from_fn(|k|
                    compensated_sum((0..20).map(|j| basis[j].0[k] * displacement[j])));
                let bending: [f64; 3] = core::array::from_fn(|k|
                    compensated_sum((0..20).map(|j| basis[j].1[k] * displacement[j])));
                let a = scale * compensated_sum(
                    (0..3).map(|k| factor.membrane[k] * membrane[k]));
                let d = scale * compensated_sum(
                    (0..3).map(|k| factor.bending[k] * bending[k]));
                if !a.is_finite() || !d.is_finite() {
                    return Err(bad("thermal shell strain projection is unrepresentable"));
                }
                out.membrane[row * n + column] = a;
                out.bending[row * n + column] = d;
            }
            row += 1;
        }
    }
    if cancelled() { return Err(PlateError::Cancelled); }
    Ok(out)
}

// Fixed, short 18/20/3-term contractions; no magnitude threshold or discarded
// residual channel. In particular [1,2^-55,-1] retains the middle contribution.
fn compensated_sum(values: impl Iterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

#[cfg(test)]
mod tests;
