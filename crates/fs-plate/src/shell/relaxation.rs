//! Passive hereditary resultants on the original CST/DKT physical strains.
//!
//! Each excess A/B/D tensor has its own relaxation time. The equilibrium
//! elastic pencil remains unchanged. Equal times are combined exactly at
//! projection; distinct times are never averaged. This is a cold spatial
//! projection, with no local time step, fictitious mass or modal damping fit.

use super::viscosity::{ShellViscosity, ShellViscosityBudget, admit, projected_strains};
use super::{ShellMesh, ShellModel};
use crate::PlateError;
use fs_math::{dd::Dd, det};
use std::collections::{BTreeMap, BTreeSet};

fn bad(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}

/// A facet's passive Maxwell excess stiffness in its original tangent axes.
/// A [N/m], B [N] and D [N m] use engineering shear in their third entries.
/// These are increments above the equilibrium pencil, not viscosity or a
/// frequency-sampled total storage tensor. Source applicability is explicit
/// caller data; this type does not fit or infer a material spectrum.
#[derive(Clone, Debug)]
pub struct ShellRelaxationBranch {
    // Reuse the one admitted physical resultant Gram factor. The numeric
    // contraction is identical; this wrapper assigns stiffness rather than
    // viscosity units and owns the distinct first-order constitutive time.
    excess: ShellViscosity,
    relaxation_time_s: f64,
}

impl ShellRelaxationBranch {
    /// Independent excess membrane and bending stiffness; B is exactly zero.
    pub fn new(
        membrane: [f64; 9],
        bending: [f64; 9],
        relaxation_time_s: f64,
    ) -> Result<Self, PlateError> {
        Self::new_coupled(membrane, [0.0; 9], bending, relaxation_time_s)
    }

    /// Complete passive excess A/B/D block, including signed membrane/bending
    /// coupling. Reversing physical layer order reverses B. Directional and
    /// exact-zero branches are legal; no positive stiffness floor is added.
    pub fn new_coupled(
        membrane: [f64; 9],
        coupling: [f64; 9],
        bending: [f64; 9],
        relaxation_time_s: f64,
    ) -> Result<Self, PlateError> {
        if !relaxation_time_s.is_finite()
            || relaxation_time_s <= 0.0
            || !(1.0 / relaxation_time_s).is_finite()
        {
            return Err(bad(
                "shell relaxation requires a finite positive time and rate",
            ));
        }
        Ok(Self {
            excess: ShellViscosity::new_coupled(membrane, coupling, bending)?,
            relaxation_time_s,
        })
    }

    /// Original excess stiffness resultants, without rotation or averaging.
    pub const fn matrices(&self) -> ([f64; 9], [f64; 9], [f64; 9]) {
        self.excess.matrices()
    }

    pub const fn relaxation_time_s(&self) -> f64 {
        self.relaxation_time_s
    }

    pub fn is_zero(&self) -> bool {
        self.excess.is_zero()
    }
}

/// Explicit spatial, source-spectrum and cold projection allocation bounds.
#[derive(Clone, Copy, Debug)]
pub struct ShellRelaxationBudget {
    pub projection: ShellViscosityBudget,
    /// Total supplied facet branches, including declared zero branches.
    pub max_element_branches: usize,
    /// Distinct times among nonzero source branches. Output storage is bounded
    /// by this count times the complete modal square, before allocation.
    pub max_relaxation_times: usize,
}

/// One complete excess modal stiffness at an unchanged constitutive time.
#[derive(Clone, Debug)]
pub struct ProjectedShellRelaxationBranch {
    relaxation_time_s: f64,
    matrix: Vec<f64>,
}

impl ProjectedShellRelaxationBranch {
    // Both original shell element owners publish already admitted complete
    // modal Grams through this one immutable result type.
    pub(super) fn from_projected(relaxation_time_s: f64, matrix: Vec<f64>) -> Self {
        Self { relaxation_time_s, matrix }
    }

    pub const fn relaxation_time_s(&self) -> f64 {
        self.relaxation_time_s
    }

    /// Row-major Phi^T K_excess Phi. Units are 1/s² only when the supplied
    /// modes are mass-normalized. Every cross-mode term remains present.
    pub fn matrix(&self) -> &[f64] {
        &self.matrix
    }
}

/// Conservative scalar-work estimate: kinematics are evaluated once per
/// facet, regardless of branch count. Material contractions and complete Gram
/// accumulation are charged per branch, and output per distinct time.
pub fn projection_terms(
    full_dofs: usize,
    elements: usize,
    modes: usize,
    element_branches: usize,
    relaxation_times: usize,
) -> Option<usize> {
    let square = modes.checked_mul(modes)?;
    let base = super::viscosity::projection_terms(full_dofs, elements, modes)?;
    let per_branch = 128usize
        .checked_add(256usize.checked_mul(modes)?)?
        .checked_add(512usize.checked_mul(square)?)?;
    base.checked_add(element_branches.checked_mul(per_branch)?)?
        .checked_add(relaxation_times.checked_mul(square)?.checked_mul(4)?)
}

/// Project every facet's complete excess spectrum on the SAME supported
/// CST/DKT coordinates. Strains are formed before material channels using the
/// instantaneous-viscosity owner's original CST and signed DKT kinematics.
/// Thus rigid motion and drilling stabilization acquire no material memory.
///
/// Outputs are ordered by increasing positive relaxation time; only exactly
/// equal f64 times aggregate. An all-zero source produces no output. A
/// nonzero branch whose supplied modes all have zero strain retains its zero
/// projected matrix and needs no memory in the consuming factorization.
///
/// Mesh, pencil and basis must refer to the same bare CST/DKT shell. Its
/// original global-axis support map is checked; geometry identity cannot be
/// recovered from public `ShellModel`. MITC, tied and oblique constraints have
/// their own spatial projection owners. No mass normalization is inferred.
pub fn project_relaxation(
    mesh: &ShellMesh,
    sections: &[Vec<ShellRelaxationBranch>],
    model: &ShellModel,
    modes: &[Vec<f64>],
    budget: ShellRelaxationBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<ProjectedShellRelaxationBranch>, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    let n = modes.len();
    if n == 0 || n > budget.projection.max_modes || n > super::mitc::viscosity::MAX_VISCOSITY_MODES
    {
        return Err(bad(
            "shell relaxation requires a nonempty bounded modal basis",
        ));
    }
    admit(mesh, sections.len(), model, budget.projection, n)?;
    if modes
        .iter()
        .any(|mode| mode.len() != model.free || mode.iter().any(|value| !value.is_finite()))
    {
        return Err(bad(
            "shell relaxation requires finite original-coordinate vectors",
        ));
    }
    let mut branches = 0usize;
    let mut times = BTreeSet::new();
    for section in sections {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        branches = branches
            .checked_add(section.len())
            .ok_or_else(|| bad("shell relaxation source extent overflows"))?;
        if branches > budget.max_element_branches {
            return Err(PlateError::Budget {
                what: "shell relaxation facet branches",
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
                        what: "shell relaxation distinct times",
                        requested: times.len(),
                        limit: budget.max_relaxation_times,
                    });
                }
            }
        }
    }
    let work = projection_terms(
        mesh.nodes.len() * 6,
        mesh.tris.len(),
        n,
        branches,
        times.len(),
    )
    .ok_or_else(|| bad("shell relaxation projection work overflows"))?;
    if work > budget.projection.max_work {
        return Err(PlateError::Budget {
            what: "shell relaxation projection work",
            requested: work,
            limit: budget.projection.max_work,
        });
    }
    let mut forms: BTreeMap<_, _> = times
        .into_iter()
        .map(|time| (time, vec![Dd::ZERO; n * n]))
        .collect();
    let mut strains = vec![[[0.0; 6]; 3]; n];
    let mut channels = vec![[0.0; 6]; n];
    for (element, (tri, section)) in mesh.tris.iter().zip(sections).enumerate() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let geometry = mesh.facet(element)?;
        if section.iter().all(ShellRelaxationBranch::is_zero) {
            continue;
        }
        projected_strains(
            &geometry,
            element,
            tri,
            model,
            modes,
            &mut strains,
            &mut cancelled,
        )?;
        let scale = det::sqrt(geometry.area_m2) * det::sqrt(1.0 / 3.0);
        for branch in section {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            if branch.is_zero() {
                continue;
            }
            let form = forms
                .get_mut(&branch.relaxation_time_s.to_bits())
                .ok_or_else(|| bad("shell relaxation source time layout changed"))?;
            for point in 0..3 {
                for (strain, channel) in strains.iter().zip(&mut channels) {
                    *channel = branch.excess.channels(&strain[point]).map(|x| scale * x);
                    if channel.iter().any(|value| !value.is_finite()) {
                        return Err(bad("shell relaxation projected channels overflow"));
                    }
                }
                for i in 0..n {
                    if cancelled() {
                        return Err(PlateError::Cancelled);
                    }
                    for j in i..n {
                        let mut value = form[n * i + j];
                        for k in 0..6 {
                            value =
                                value + Dd::from_f64(channels[i][k]) * Dd::from_f64(channels[j][k]);
                        }
                        if !value.hi.is_finite() || !value.lo.is_finite() {
                            return Err(bad("shell relaxation projected stiffness overflow"));
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
            f64::from_bits(time), matrix,
        ));
    }
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    Ok(projected)
}
