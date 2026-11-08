//! Passive hereditary A/B/D and transverse shear on the original MITC space.
//!
//! Every facet can carry independently oriented excess tensors and times.
//! The original equilibrium section owns its director-bubble elimination;
//! material memory neither changes that map nor introduces another time step.
//! This is the same fixed condensed trial-space approximation as the original
//! modal shell, not frequency-dependent elimination of dynamic bubbles.

use super::viscosity::{ShearDissipation, ShearShellViscosity, projected_strains};
use super::{MitcBudget, MitcShellModel, PlateError, element};
use crate::shell::relaxation::ProjectedShellRelaxationBranch;
use fs_math::{dd::Dd, det};
use std::collections::{BTreeMap, BTreeSet};

fn bad(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}

/// One passive excess resultant tensor above the relaxed elastic section.
/// A and transverse shear have units N/m, B has units N, D has units N m.
/// Engineering shear and the facet tangent axes match the original elastic
/// section. Any physical shear correction is included explicitly by the
/// source compiler, just as for instantaneous MITC viscosity.
#[derive(Clone, Debug)]
pub struct ShearShellRelaxationBranch {
    excess: ShearShellViscosity,
    factor: ShearDissipation,
    relaxation_time_s: f64,
}

impl ShearShellRelaxationBranch {
    /// Uncoupled membrane/bending excess with full symmetric transverse shear
    /// `[S13_13, S13_23, S23_23]` and its own positive material time.
    pub fn new(
        membrane: [f64; 9],
        bending: [f64; 9],
        shear: [f64; 3],
        relaxation_time_s: f64,
    ) -> Result<Self, PlateError> {
        Self::new_coupled(membrane, [0.0; 9], bending, shear, relaxation_time_s)
    }

    /// Full passive excess A/B/D block and transverse shear. Signed B is
    /// admitted through the original complete resultant factor; B need not be
    /// positive by itself. Directional and exact-zero increments are legal.
    pub fn new_coupled(
        membrane: [f64; 9],
        coupling: [f64; 9],
        bending: [f64; 9],
        shear: [f64; 3],
        relaxation_time_s: f64,
    ) -> Result<Self, PlateError> {
        if !relaxation_time_s.is_finite()
            || relaxation_time_s <= 0.0
            || !(1.0 / relaxation_time_s).is_finite()
        {
            return Err(bad(
                "MITC relaxation requires finite positive material time and rate",
            ));
        }
        let excess = ShearShellViscosity::new_coupled(membrane, coupling, bending, shear)?;
        let factor = excess.dissipation()?;
        Ok(Self {
            excess,
            factor,
            relaxation_time_s,
        })
    }

    /// Original A/B/D/shear stiffness increments, without spectral averaging.
    pub const fn coupled_matrices(&self) -> ([f64; 9], [f64; 9], [f64; 9], [[f64; 2]; 2]) {
        self.excess.coupled_matrices()
    }

    pub const fn relaxation_time_s(&self) -> f64 {
        self.relaxation_time_s
    }

    pub(super) const fn excess(&self) -> &ShearShellViscosity {
        &self.excess
    }

    pub(super) const fn factor(&self) -> &ShearDissipation {
        &self.factor
    }

    pub fn is_zero(&self) -> bool {
        self.excess.is_zero()
    }
}

/// Complete cold spatial/projection allocation and work ceilings.
#[derive(Clone, Copy, Debug)]
pub struct ShearShellRelaxationBudget {
    pub mesh: MitcBudget,
    pub max_modes: usize,
    /// Total supplied facet branches, including explicitly zero increments.
    pub max_element_branches: usize,
    /// Distinct positive times among nonzero increments; bounds the number
    /// of complete output modal squares before allocation.
    pub max_relaxation_times: usize,
    pub max_projection_terms: usize,
}

/// Conservative scalar cold-work count. Original physical expansion and
/// equilibrium condensation occur once per mode/facet, while passive channel
/// contractions and complete Gram accumulation are charged per source branch.
pub fn projection_terms(
    nodes: usize,
    triangles: usize,
    modes: usize,
    element_branches: usize,
    relaxation_times: usize,
) -> Option<usize> {
    let square = modes.checked_mul(modes)?;
    let original = super::viscosity::projection_terms(nodes, triangles, modes)?;
    let per_branch = 4096usize
        .checked_add(4096usize.checked_mul(modes)?)?
        .checked_add(2048usize.checked_mul(square)?)?;
    original
        .checked_add(nodes.checked_mul(modes)?.checked_mul(72)?)?
        .checked_add(element_branches.checked_mul(per_branch)?)?
        .checked_add(relaxation_times.checked_mul(square)?.checked_mul(4)?)
}

/// Project complete per-facet material spectra from the original supported
/// MITC model. Input vectors are in its reduced support coordinates. The
/// immutable model supplies geometry, directors, elastic sections and the
/// original work-conjugate expansion; its bare `dof_map` is never substituted
/// for oblique support frames.
///
/// The same PreparedStrainField, original elastic bubble map and positive
/// seven-point quadrature as viscosity contract membrane, signed curvature
/// and transverse shear before passive material factors. Numerical drilling,
/// prestress and mass contribute no hereditary channel. Geometry/strains are
/// evaluated once per facet across the entire spectrum.
///
/// Exactly equal positive f64 relaxation times aggregate in double-double
/// Gram sums. Distinct times and all modal cross terms remain present.
/// Arbitrary finite test vectors are allowed; only the consuming modal owner
/// may establish mass normalization and stiffness units of 1/s².
/// Source bands and the interpretation as excess above relaxed equilibrium
/// remain caller data. Tied assemblies and independent dynamic bubbles have
/// different coordinate owners and cannot be supplied as this model.
pub fn project_relaxation(
    model: &MitcShellModel,
    spectra: &[Vec<ShearShellRelaxationBranch>],
    modes: &[Vec<f64>],
    budget: ShearShellRelaxationBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<ProjectedShellRelaxationBranch>, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    let mesh = &model.mesh;
    let n = modes.len();
    if mesh.nodes.len() > budget.mesh.max_nodes
        || mesh.tris.len() > budget.mesh.max_triangles
        || spectra.len() != mesh.tris.len()
        || n == 0
        || n > budget.max_modes
        || n > super::viscosity::MAX_VISCOSITY_MODES
        || modes.iter().any(|mode| mode.len() != model.pencil.free)
    {
        return Err(bad(
            "MITC relaxation requires its complete bounded model, facet spectra and reduced vectors",
        ));
    }
    let mut branches = 0usize;
    let mut times = BTreeSet::new();
    for section in spectra {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        branches = branches
            .checked_add(section.len())
            .ok_or_else(|| bad("MITC relaxation source extent overflows"))?;
        if branches > budget.max_element_branches {
            return Err(PlateError::Budget {
                what: "MITC relaxation facet branches",
                requested: branches,
                limit: budget.max_element_branches,
            });
        }
        for branch in section {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            if !branch.is_zero() {
                times.insert(branch.relaxation_time_s.to_bits());
                if times.len() > budget.max_relaxation_times {
                    return Err(PlateError::Budget {
                        what: "MITC relaxation distinct times",
                        requested: times.len(),
                        limit: budget.max_relaxation_times,
                    });
                }
            }
        }
    }
    let work = projection_terms(mesh.nodes.len(), mesh.tris.len(), n, branches, times.len())
        .ok_or_else(|| bad("MITC relaxation projection work overflows"))?;
    if work > budget.max_projection_terms {
        return Err(PlateError::Budget {
            what: "MITC relaxation projection work",
            requested: work,
            limit: budget.max_projection_terms,
        });
    }
    let mut physical = Vec::with_capacity(n);
    for mode in modes {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        physical.push(model.expand(mode)?);
    }
    let mut forms: BTreeMap<_, _> = times
        .into_iter()
        .map(|time| (time, vec![Dd::ZERO; n * n]))
        .collect();
    let mut strains = vec![([0.0; 3], [0.0; 3], [0.0; 2]); n];
    let mut channels = vec![[0.0; 12]; n];
    let rule = element::stiffness_rule();
    for (element, (tri, section)) in mesh.tris.iter().zip(spectra).enumerate() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if section.iter().all(ShearShellRelaxationBranch::is_zero) {
            continue;
        }
        let geometry = mesh.facet(element)?;
        let field = element::PreparedStrainField::new(
            mesh,
            tri,
            geometry,
            tri.map(|i| model.directors.vectors[i]),
            &model.sections[element],
            &mut cancelled,
        )?;
        let mut local = Vec::with_capacity(n);
        for mode in &physical {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            let nodal = core::array::from_fn(|i| mode[6 * tri[i / 6] + i % 6]);
            let expanded = field.expand(&nodal);
            if expanded.iter().any(|value| !value.is_finite()) {
                return Err(bad("MITC relaxation original elastic expansion overflows"));
            }
            local.push(expanded);
        }
        for &(lambda, weight) in &rule {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            projected_strains(&field, &local, lambda, &mut strains, &mut cancelled)?;
            let scale = det::sqrt(geometry.area_m2) * det::sqrt(weight);
            for branch in section {
                if cancelled() {
                    return Err(PlateError::Cancelled);
                }
                if branch.is_zero() {
                    continue;
                }
                for (&(membrane, bending, shear), channel) in strains.iter().zip(&mut channels) {
                    *channel = branch
                        .factor
                        .channels(membrane, bending, shear)
                        .map(|value| scale * value);
                    if channel.iter().any(|value| !value.is_finite()) {
                        return Err(bad("MITC relaxation physical channels overflow"));
                    }
                }
                let form = forms
                    .get_mut(&branch.relaxation_time_s.to_bits())
                    .ok_or_else(|| bad("MITC relaxation source time layout changed"))?;
                for i in 0..n {
                    if cancelled() {
                        return Err(PlateError::Cancelled);
                    }
                    for j in i..n {
                        let mut value = form[n * i + j];
                        for k in 0..12 {
                            value =
                                value + Dd::from_f64(channels[i][k]) * Dd::from_f64(channels[j][k]);
                        }
                        if !value.hi.is_finite() || !value.lo.is_finite() {
                            return Err(bad("MITC relaxation projected stiffness overflows"));
                        }
                        form[n * i + j] = value;
                    }
                }
            }
        }
    }
    let mut projected = Vec::with_capacity(forms.len());
    for (time, form) in forms {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let mut matrix = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                matrix[n * i + j] = form[n * i + j].to_f64();
                matrix[n * j + i] = matrix[n * i + j];
            }
        }
        projected.push(ProjectedShellRelaxationBranch::from_projected(
            f64::from_bits(time),
            matrix,
        ));
    }
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    Ok(projected)
}

#[cfg(test)]
mod tests;
