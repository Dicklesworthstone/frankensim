//! Fixed-volume nonideal solution diffusion on the existing adaptive timeline.
//!
//! Every connection consumes the original source-owned solvent-frame MS flux.
//! Conservative explicit pair exchanges are swept in canonical endpoint order;
//! the enclosing full/two-half comparison controls their first-order temporal
//! error. This mode admits no water/membrane connections, so a closed mixing
//! free-energy decrease does not omit mechanical pressure/volume work.

use super::{
    MAX_COMPARTMENTS, MAX_MEMBRANES, MAX_OSMOTIC_SPECIES, OsmoticCompartment, OsmoticNetwork,
    OsmoticNetworkError, compensated_add, poll, positive,
};
use fs_material::osmosis::{
    OsmosisError,
    diffusion::{
        SolutionDiffusionError, SolutionDiffusionSolve, SolutionDiffusionWorkspace,
        SolutionMaxwellStefan,
    },
};

/// Physical conservation/passivity and withdrawal admission, independent of
/// the enclosing adaptive timeline's temporal error tolerances.
#[derive(Clone, Copy, Debug)]
pub struct SolutionDiffusionControls {
    pub constitutive: SolutionDiffusionSolve,
    /// Largest withdrawal / original individual donor amount at each split pair
    /// stage, in (0,1). This is not a whole-window withdrawal fraction.
    pub maximum_species_fraction: f64,
    /// Absolute closed-window individual-species balance allowance [mol].
    pub amount_balance_tolerance_mol: f64,
    /// Absolute original-to-final mixing-free-energy increase allowance [J].
    /// This same limit also applies to every accepted physical pair step.
    pub free_energy_tolerance_j: f64,
}

/// Recoverable finite-step failure; source/domain refusals are different errors.
#[derive(Clone, Debug, PartialEq)]
pub enum SolutionDiffusionStepRefusal {
    Withdrawal {
        species: usize,
        fraction: f64,
        maximum: f64,
    },
    FreeEnergy {
        change_j: f64,
        maximum_j: f64,
    },
}

/// One physical diffusion face with complete original material/friction data.
#[derive(Clone, Debug)]
pub struct SolutionDiffusionLink {
    left: usize,
    right: usize,
    area_m2: f64,
    distance_m: f64,
    material: SolutionMaxwellStefan,
    workspace: SolutionDiffusionWorkspace,
}
impl SolutionDiffusionLink {
    /// Endpoints are canonicalized because diffusion has no imposed direction.
    /// Explicit zero area closes a face; distance must remain positive.
    pub fn new(
        left: usize,
        right: usize,
        area_m2: f64,
        distance_m: f64,
        material: SolutionMaxwellStefan,
    ) -> Result<Self, OsmoticNetworkError> {
        if left == right || !area_m2.is_finite() || area_m2 < 0.0 || !positive(distance_m) {
            return Err(OsmoticNetworkError::Invalid(
                "distinct endpoints, nonnegative area and positive diffusion distance required",
            ));
        }
        let workspace = material
            .workspace()
            .map_err(OsmoticNetworkError::SolutionDiffusion)?;
        Ok(Self {
            left: left.min(right),
            right: left.max(right),
            area_m2,
            distance_m,
            material,
            workspace,
        })
    }
    pub const fn endpoints(&self) -> [usize; 2] {
        [self.left, self.right]
    }
    pub fn material(&self) -> &SolutionMaxwellStefan {
        &self.material
    }
    fn heap_bytes(&self) -> usize {
        self.material.heap_bytes() + self.workspace.workspace_bytes()
    }
}

impl OsmoticNetwork {
    /// Prepare a diffusion-only network of the SAME original convex solution
    /// model, species order and fixed solvent volumes. All solutes must be
    /// strictly present; this initial law does not invent zero-log limits.
    ///
    /// Reuses the existing bounded allocation-free adaptive timeline. Explicit
    /// conservative diffusion stages replace membrane water-root stages in this
    /// mode; the method remains first order and may need small steps. Mechanical
    /// pressure laws are retained but perform no work because volumes are fixed.
    /// Distinct physical endpoint pairs are required; duplicate faces refuse.
    pub fn with_solution_diffusion(
        species: Vec<String>,
        cells: Vec<OsmoticCompartment>,
        mut links: Vec<SolutionDiffusionLink>,
        controls: SolutionDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        let friction = controls.constitutive.friction;
        if cells.is_empty()
            || cells.len() > MAX_COMPARTMENTS
            || links.len() > MAX_MEMBRANES
            || !positive(controls.maximum_species_fraction)
            || controls.maximum_species_fraction >= 1.0
            || !positive(controls.amount_balance_tolerance_mol)
            || !positive(controls.free_energy_tolerance_j)
            || !positive(controls.constitutive.entropy_tolerance_w_m3_k)
            || !positive(friction.relative_tolerance)
            || friction.relative_tolerance > 1e-3
            || !positive(friction.absolute_gradient_tolerance_per_m)
            || !positive(friction.absolute_molar_flux_tolerance_mol_m2_s)
        {
            return Err(OsmoticNetworkError::Invalid(
                "bounded diffusion graph and explicit positive physical budgets required",
            ));
        }
        let source = cells[0].solution.nonideal_model().ok_or(OsmoticNetworkError::Invalid(
            "solution diffusion requires the explicit original virial model (zero interactions admit its ideal limit)"))?;
        if cells.iter().any(|cell| {
            cell.solution.nonideal_model() != Some(source)
                || cell.solution.solute_moles().iter().any(|&n| !positive(n))
        }) {
            return Err(OsmoticNetworkError::Invalid(
                "one original thermodynamic model and positive inventories required",
            ));
        }
        links.sort_by_key(|link| (link.left, link.right));
        if links
            .windows(2)
            .any(|pair| pair[0].endpoints() == pair[1].endpoints())
        {
            return Err(OsmoticNetworkError::Invalid(
                "duplicate physical solution-diffusion face",
            ));
        }
        for link in &links {
            if link.right >= cells.len()
                || link.material.species().len() != species.len()
                || !link
                    .material
                    .species()
                    .iter()
                    .zip(&species)
                    .all(|(id, name)| id.as_str() == name)
            {
                return Err(OsmoticNetworkError::Invalid(
                    "diffusion faces must retain the exact original species order and existing endpoints",
                ));
            }
            link.material
                .admit_solution(&cells[link.left].solution)
                .map_err(OsmoticNetworkError::SolutionDiffusion)?;
            link.material
                .admit_solution(&cells[link.right].solution)
                .map_err(OsmoticNetworkError::SolutionDiffusion)?;
        }
        let required = links
            .iter()
            .fold(
                links
                    .capacity()
                    .checked_mul(core::mem::size_of::<SolutionDiffusionLink>()),
                |total, link| total.and_then(|bytes| bytes.checked_add(link.heap_bytes())),
            )
            .ok_or(OsmoticNetworkError::Invalid(
                "diffusion workspace byte count overflow",
            ))?;
        let remaining = maximum_workspace_bytes.checked_sub(required).ok_or(
            OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_workspace_bytes,
            },
        )?;
        let mut network = Self::new(species, cells, Vec::new(), remaining)?;
        network.workspace_bytes += required;
        network.diffusion = links;
        network.diffusion_controls = Some(controls);
        Ok(network)
    }
}

fn source(error: OsmosisError) -> OsmoticNetworkError {
    OsmoticNetworkError::SolutionDiffusion(SolutionDiffusionError::Thermodynamics(error))
}
fn finite(value: f64) -> Result<f64, OsmoticNetworkError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OsmoticNetworkError::Invalid(
            "unrepresentable finite solution diffusion account",
        ))
    }
}

pub(super) fn sweep(
    cells: &mut [OsmoticCompartment],
    links: &mut [SolutionDiffusionLink],
    dt: f64,
    controls: SolutionDiffusionControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    for link in links {
        poll(cancelled)?;
        if link.area_m2 == 0.0 {
            continue;
        }
        let left = &cells[link.left].solution;
        let right = &cells[link.right].solution;
        let n = left.solute_moles().len();
        let mut flux = [0.0; MAX_OSMOTIC_SPECIES];
        *evaluations = evaluations
            .checked_add(1)
            .ok_or(OsmoticNetworkError::WorkBudget)?;
        link.material
            .flux_into(
                left,
                right,
                link.distance_m,
                controls.constitutive,
                &mut link.workspace,
                &mut flux[..n],
                &mut *cancelled,
            )
            .map_err(|error| {
                if error == SolutionDiffusionError::Cancelled {
                    OsmoticNetworkError::Cancelled
                } else {
                    OsmoticNetworkError::SolutionDiffusion(error)
                }
            })?;
        let area_time = dt * link.area_m2;
        if !positive(area_time) {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        let mut next_left = [0.0; MAX_OSMOTIC_SPECIES];
        let mut next_right = [0.0; MAX_OSMOTIC_SPECIES];
        for i in 0..n {
            let transfer = finite(area_time * flux[i])?;
            if flux[i] != 0.0 && transfer == 0.0 {
                return Err(OsmoticNetworkError::StepUnderflow);
            }
            let (a, b) = (left.solute_moles()[i], right.solute_moles()[i]);
            let donor = if transfer >= 0.0 { a } else { b };
            let fraction = transfer.abs() / donor;
            if !fraction.is_finite() || fraction > controls.maximum_species_fraction {
                return Err(OsmoticNetworkError::DiffusionStep(
                    SolutionDiffusionStepRefusal::Withdrawal {
                        species: i,
                        fraction,
                        maximum: controls.maximum_species_fraction,
                    },
                ));
            }
            next_left[i] = finite(a - transfer)?;
            next_right[i] = finite(b + transfer)?;
            if transfer != 0.0 && next_left[i] == a && next_right[i] == b {
                return Err(OsmoticNetworkError::StepUnderflow);
            }
            if !positive(next_left[i]) || !positive(next_right[i]) {
                return Err(OsmoticNetworkError::Invalid(
                    "positive solution endpoint is unrepresentable",
                ));
            }
            let residual = (next_left[i] - a) + (next_right[i] - b);
            if residual.abs() > controls.amount_balance_tolerance_mol {
                return Err(OsmoticNetworkError::DiffusionBalance {
                    species: i,
                    residual_mol: residual,
                });
            }
        }
        // Admit BOTH original thermodynamic endpoints before modifying either.
        left.osmotic_pressure_at(left.solvent_volume_m3(), &next_left[..n])
            .map_err(source)?;
        right
            .osmotic_pressure_at(right.solvent_volume_m3(), &next_right[..n])
            .map_err(source)?;
        let mut change = 0.0;
        let mut correction = 0.0;
        for (old, next) in [(left, &next_left[..n]), (right, &next_right[..n])] {
            let mut before = [0.0; MAX_OSMOTIC_SPECIES];
            let mut after = [0.0; MAX_OSMOTIC_SPECIES];
            old.concentrations_into(&mut before[..n]).map_err(source)?;
            for i in 0..n {
                after[i] = next[i] / old.solvent_volume_m3();
            }
            let density = link
                .material
                .thermodynamics()
                .reference_cancelled_free_energy_change_j_m3(&before[..n], &after[..n])
                .map_err(OsmosisError::from)
                .map_err(source)?;
            compensated_add(
                &mut change,
                &mut correction,
                finite(density * old.solvent_volume_m3())?,
            );
        }
        finite(change)?;
        if change > controls.free_energy_tolerance_j {
            return Err(OsmoticNetworkError::DiffusionStep(
                SolutionDiffusionStepRefusal::FreeEnergy {
                    change_j: change,
                    maximum_j: controls.free_energy_tolerance_j,
                },
            ));
        }
        poll(cancelled)?;
        let volumes = [left.solvent_volume_m3(), right.solvent_volume_m3()];
        cells[link.left]
            .solution
            .replace_inventories(volumes[0], &next_left[..n])
            .map_err(source)?;
        cells[link.right]
            .solution
            .replace_inventories(volumes[1], &next_right[..n])
            .map_err(source)?;
    }
    Ok(())
}

pub(super) fn free_energy_change(
    original: &[OsmoticCompartment],
    candidate: &[OsmoticCompartment],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, OsmoticNetworkError> {
    let mut sum = 0.0;
    let mut correction = 0.0;
    for (before, after) in original.iter().zip(candidate) {
        poll(cancelled)?;
        let (before, after) = (&before.solution, &after.solution);
        if before.solvent_volume_m3().to_bits() != after.solvent_volume_m3().to_bits()
            || before.nonideal_model() != after.nonideal_model()
        {
            return Err(OsmoticNetworkError::Invalid(
                "closed diffusion changed solvent volume or original free energy",
            ));
        }
        let n = before.solute_moles().len();
        let mut a = [0.0; MAX_OSMOTIC_SPECIES];
        let mut b = [0.0; MAX_OSMOTIC_SPECIES];
        before.concentrations_into(&mut a[..n]).map_err(source)?;
        after.concentrations_into(&mut b[..n]).map_err(source)?;
        let law = before.nonideal_model().ok_or(OsmoticNetworkError::Invalid(
            "missing original solution free energy",
        ))?;
        let density = law
            .reference_cancelled_free_energy_change_j_m3(&a[..n], &b[..n])
            .map_err(OsmosisError::from)
            .map_err(source)?;
        compensated_add(
            &mut sum,
            &mut correction,
            finite(density * before.solvent_volume_m3())?,
        );
    }
    finite(sum)
}

#[cfg(test)]
mod tests;
