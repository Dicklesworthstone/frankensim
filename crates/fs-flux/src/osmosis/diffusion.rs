//! Fixed-volume nonideal solution diffusion on the existing adaptive timeline.
//!
//! Every connection consumes the original source-owned solvent-frame MS flux.
//! Conservative pair exchanges are swept in canonical endpoint order;
//! the enclosing full/two-half comparison controls their first-order temporal
//! error. This mode admits no water/membrane connections, so a closed mixing
//! free-energy decrease does not omit mechanical pressure/volume work.
//! A zero-current option evolves neutral ionic inventories on a forest, with
//! optional charge-constrained implicit binary-salt exchange. General closed
//! graphs can select a native global potential solve and simultaneous molecular
//! updates. These ionic paths share the original source-qualified face and the
//! same adaptive timeline, conservation, charge and free-energy gates.

mod electrical;
mod electroneutral;
mod implicit;
mod implicit_salt;
mod multicomponent;
pub(crate) use electroneutral::ElectroneutralDiffusion;
pub use electroneutral::{
    ElectroneutralDiffusionControls, ElectroneutralDiffusionReport, IonicDiffusionConnection,
    IonicDiffusionSnapshot,
};
pub use fs_material::osmosis::diffusion::electrical::ZeroCurrentDiffusion;
pub(crate) use implicit::ImplicitBinaryDiffusion;
pub use implicit::ImplicitBinaryDiffusionControls;
pub(crate) use implicit_salt::ImplicitSaltDiffusion;
pub use implicit_salt::ImplicitSaltDiffusionControls;
pub(crate) use multicomponent::ImplicitMulticomponentDiffusion;
pub use multicomponent::ImplicitMulticomponentDiffusionControls;

use super::{
    MAX_COMPARTMENTS, MAX_MEMBRANES, MAX_OSMOTIC_SPECIES, OsmoticCompartment, OsmoticNetwork,
    OsmoticNetworkError, compensated_add, nonnegative, poll, positive,
};
use fs_material::osmosis::{
    AqueousSolution, OsmosisError,
    diffusion::{
        SolutionDiffusionError, SolutionDiffusionSolve, SolutionDiffusionWorkspace,
        SolutionMaxwellStefan, source::MaterialSolutionMaxwellStefan,
    },
};
use std::sync::Arc;

/// Physical conservation/passivity and withdrawal admission, independent of
/// the enclosing adaptive timeline's temporal error tolerances.
#[derive(Clone, Copy, Debug)]
pub struct SolutionDiffusionControls {
    pub constitutive: SolutionDiffusionSolve,
    /// Largest withdrawal / original individual donor amount at each split pair
    /// stage, in (0,1). Global electrical stages instead sum ALL outgoing
    /// withdrawals per donor species before checking this same bound.
    /// This is not a whole-window withdrawal fraction.
    pub maximum_species_fraction: f64,
    /// Absolute closed-window individual-species balance allowance [mol].
    pub amount_balance_tolerance_mol: f64,
    /// Absolute original-to-final mixing-free-energy increase allowance [J].
    /// This same limit also applies to every accepted physical pair step, or
    /// the simultaneous whole-network stage under a global electrical closure.
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
    material: DiffusionMaterial,
    workspace: SolutionDiffusionWorkspace,
    electrical: Option<Arc<ZeroCurrentDiffusion>>,
}
#[derive(Clone, Debug)]
enum DiffusionMaterial {
    Constant(SolutionMaxwellStefan),
    Source(Box<MaterialSolutionMaxwellStefan>),
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
        Self::prepare(
            left,
            right,
            area_m2,
            distance_m,
            DiffusionMaterial::Constant(material),
        )
    }
    /// State-dependent original material pairs. The same finite-volume
    /// controller and shared friction solver are used; source queries allocate.
    /// Every pair keeps its original selected claim across the complete window.
    pub fn from_material_sources(
        left: usize,
        right: usize,
        area_m2: f64,
        distance_m: f64,
        material: MaterialSolutionMaxwellStefan,
    ) -> Result<Self, OsmoticNetworkError> {
        Self::prepare(
            left,
            right,
            area_m2,
            distance_m,
            DiffusionMaterial::Source(Box::new(material)),
        )
    }
    fn prepare(
        left: usize,
        right: usize,
        area_m2: f64,
        distance_m: f64,
        material: DiffusionMaterial,
    ) -> Result<Self, OsmoticNetworkError> {
        if left == right || !area_m2.is_finite() || area_m2 < 0.0 || !positive(distance_m) {
            return Err(OsmoticNetworkError::Invalid(
                "distinct endpoints, nonnegative area and positive diffusion distance required",
            ));
        }
        let workspace = match &material {
            DiffusionMaterial::Constant(law) => law.workspace(),
            DiffusionMaterial::Source(law) => law.workspace(),
        }
        .map_err(OsmoticNetworkError::SolutionDiffusion)?;
        Ok(Self {
            left: left.min(right),
            right: left.max(right),
            area_m2,
            distance_m,
            material,
            workspace,
            electrical: None,
        })
    }
    pub const fn endpoints(&self) -> [usize; 2] {
        [self.left, self.right]
    }
    /// Original physical owner. For source-backed links its stored binary
    /// matrix describes only the reference state; advances always refresh it.
    pub fn material(&self) -> &SolutionMaxwellStefan {
        match &self.material {
            DiffusionMaterial::Constant(law) => law,
            DiffusionMaterial::Source(law) => law.reference_law(),
        }
    }
    /// Original immutable cards, complete pair assignments and selected IDs.
    /// Runtime source receipts are transient, not a last-successful-state cache.
    pub fn material_sources(&self) -> Option<&MaterialSolutionMaxwellStefan> {
        match &self.material {
            DiffusionMaterial::Constant(_) => None,
            DiffusionMaterial::Source(law) => Some(law),
        }
    }
    /// Common ionic assignment selected by the enclosing network, if any.
    pub fn zero_current_diffusion(&self) -> Option<&ZeroCurrentDiffusion> {
        self.electrical.as_deref()
    }
    fn heap_bytes(&self) -> usize {
        let material = match &self.material {
            DiffusionMaterial::Constant(law) => law.heap_bytes(),
            DiffusionMaterial::Source(law) => {
                core::mem::size_of::<MaterialSolutionMaxwellStefan>()
                    + law.estimated_heap_payload_bytes()
            }
        };
        material.saturating_add(self.workspace.workspace_bytes())
    }
    fn qualify_sources(
        &self,
        volumes: [f64; 2],
        amounts: [&[f64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        if let Some(electrical) = &self.electrical {
            for inventory in amounts {
                electrical
                    .admit_inventories(inventory)
                    .map_err(diffusion_source)?;
            }
        }
        if let Some(source) = self.material_sources() {
            if let Some(electrical) = &self.electrical {
                source.resolve_zero_current_inventories(volumes, amounts, electrical, cancelled)
            } else {
                source.resolve_inventories(volumes, amounts, cancelled)
            }
            .map_err(diffusion_source)?;
        }
        Ok(())
    }

    /// The original constitutive path for explicit and implicit trial states.
    /// Every call, including root search and rejected time trials, is charged.
    fn flux_into(
        &mut self,
        left: &AqueousSolution,
        right: &AqueousSolution,
        controls: SolutionDiffusionSolve,
        output: &mut [f64],
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        *evaluations = evaluations
            .checked_add(1)
            .ok_or(OsmoticNetworkError::WorkBudget)?;
        if let Some(electrical) = &self.electrical {
            return match &self.material {
                DiffusionMaterial::Constant(material) => material
                    .flux_zero_current_into(
                        left,
                        right,
                        self.distance_m,
                        controls,
                        electrical,
                        &mut self.workspace,
                        output,
                        cancelled,
                    )
                    .map(|_| ()),
                DiffusionMaterial::Source(material) => material
                    .flux_zero_current_into(
                        left,
                        right,
                        self.distance_m,
                        controls,
                        electrical,
                        &mut self.workspace,
                        output,
                        cancelled,
                    )
                    .map(|_| ()),
            }
            .map_err(diffusion_source);
        }
        match &self.material {
            DiffusionMaterial::Constant(material) => material
                .flux_into(
                    left,
                    right,
                    self.distance_m,
                    controls,
                    &mut self.workspace,
                    output,
                    cancelled,
                )
                .map(|_| ()),
            DiffusionMaterial::Source(material) => material
                .flux_into(
                    left,
                    right,
                    self.distance_m,
                    controls,
                    &mut self.workspace,
                    output,
                    cancelled,
                )
                .map(|_| ()),
        }
        .map_err(diffusion_source)
    }
}

impl OsmoticNetwork {
    /// Prepared diffusion links, including original source cards and pinned
    /// claim identities. Current endpoint receipts can be queried explicitly
    /// through `material_sources`; the timeline retains no query cache.
    pub fn solution_diffusion_links(&self) -> &[SolutionDiffusionLink] {
        &self.diffusion
    }

    /// Prepare a diffusion-only network of the SAME original convex solution
    /// model, species order and fixed solvent volumes. Exact absent solutes are
    /// retained: the original logarithmic-mean face supplies a finite limiting
    /// flux into pure solvent without adding a concentration floor or seed.
    ///
    /// Reuses the existing bounded adaptive timeline. Constant-data links retain
    /// allocation-free advancement. Material-source links allocate original
    /// queries and receipts at every stage; these transient allocations are not
    /// covered by `maximum_workspace_bytes`. That limit charges retained arrays,
    /// models and the estimated FULL card payload per source pair, including
    /// unselected claims/observations/models and excluding allocator overhead.
    /// No runtime source receipts are retained in the prepared network. Explicit
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
        let required = prepare_solution_diffusion(
            &species,
            &cells,
            &mut links,
            controls,
            maximum_workspace_bytes,
        )?;
        let remaining = maximum_workspace_bytes - required;
        let mut network = Self::new(species, cells, Vec::new(), remaining)?;
        network.workspace_bytes += required;
        network.diffusion = links;
        network.diffusion_controls = Some(controls);
        Ok(network)
    }
}

/// Shared preparation for fixed-volume and coupled-water networks.
pub(super) fn prepare_solution_diffusion(
    species: &[String],
    cells: &[OsmoticCompartment],
    links: &mut Vec<SolutionDiffusionLink>,
    controls: SolutionDiffusionControls,
    maximum_workspace_bytes: usize,
) -> Result<usize, OsmoticNetworkError> {
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
            || cell
                .solution
                .solute_moles()
                .iter()
                .any(|&n| !nonnegative(n))
    }) {
        return Err(OsmoticNetworkError::Invalid(
            "one original thermodynamic model and nonnegative inventories required",
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
    for link in links.iter() {
        if link.right >= cells.len()
            || link.material().species().len() != species.len()
            || !link
                .material()
                .species()
                .iter()
                .zip(species)
                .all(|(id, name)| id.as_str() == name)
        {
            return Err(OsmoticNetworkError::Invalid(
                "diffusion faces must retain the exact original species order and existing endpoints",
            ));
        }
        link.material()
            .admit_solution(&cells[link.left].solution)
            .map_err(OsmoticNetworkError::SolutionDiffusion)?;
        link.material()
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
    maximum_workspace_bytes
        .checked_sub(required)
        .ok_or(OsmoticNetworkError::Workspace {
            required_bytes: required,
            budget_bytes: maximum_workspace_bytes,
        })?;
    for link in links.iter() {
        link.qualify_sources(
            [
                cells[link.left].solution.solvent_volume_m3(),
                cells[link.right].solution.solvent_volume_m3(),
            ],
            [
                cells[link.left].solution.solute_moles(),
                cells[link.right].solution.solute_moles(),
            ],
            &mut || false,
        )?;
    }
    Ok(required)
}

fn source(error: OsmosisError) -> OsmoticNetworkError {
    OsmoticNetworkError::SolutionDiffusion(SolutionDiffusionError::Thermodynamics(error))
}
fn diffusion_source(error: SolutionDiffusionError) -> OsmoticNetworkError {
    if error == SolutionDiffusionError::Cancelled {
        OsmoticNetworkError::Cancelled
    } else {
        OsmoticNetworkError::SolutionDiffusion(error)
    }
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

#[allow(clippy::too_many_arguments)]
pub(super) fn sweep(
    cells: &mut [OsmoticCompartment],
    links: &mut [SolutionDiffusionLink],
    mut implicit: Option<&mut ImplicitBinaryDiffusion>,
    mut multicomponent: Option<&mut ImplicitMulticomponentDiffusion>,
    mut salt: Option<&mut ImplicitSaltDiffusion>,
    electrical: Option<&mut ElectroneutralDiffusion>,
    dt: f64,
    controls: SolutionDiffusionControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    if let Some(electrical) = electrical {
        return electrical.sweep(cells, links, dt, controls, evaluations, cancelled);
    }
    for link in links.iter_mut() {
        poll(cancelled)?;
        if link.area_m2 == 0.0 {
            continue;
        }
        let left = &cells[link.left].solution;
        let right = &cells[link.right].solution;
        let n = left.solute_moles().len();
        let mut flux = [0.0; MAX_OSMOTIC_SPECIES];
        let mut area_time = 0.0;
        let mut transfers = [0.0; MAX_OSMOTIC_SPECIES];
        let implicit_transfer = if let Some(implicit) = implicit.as_deref_mut() {
            transfers[0] = implicit.transfer(
                left,
                right,
                link,
                dt,
                controls.constitutive,
                evaluations,
                cancelled,
            )?;
            true
        } else if let Some(implicit) = multicomponent.as_deref_mut() {
            implicit.transfers_into(
                left,
                right,
                link,
                dt,
                controls.constitutive,
                evaluations,
                &mut transfers[..n],
                cancelled,
            )?;
            true
        } else if let Some(implicit) = salt.as_deref_mut() {
            implicit.transfers_into(
                left,
                right,
                link,
                dt,
                controls.constitutive,
                evaluations,
                cancelled,
                &mut transfers[..n],
            )?;
            true
        } else {
            link.flux_into(
                left,
                right,
                controls.constitutive,
                &mut flux[..n],
                evaluations,
                cancelled,
            )?;
            area_time = dt * link.area_m2;
            if !positive(area_time) {
                return Err(OsmoticNetworkError::StepUnderflow);
            }
            false
        };
        let mut next_left = [0.0; MAX_OSMOTIC_SPECIES];
        let mut next_right = [0.0; MAX_OSMOTIC_SPECIES];
        for i in 0..n {
            // Preserve the explicit owner's original per-species arithmetic
            // and refusal order, including recoverable withdrawal checks.
            let transfer = if implicit_transfer {
                transfers[i]
            } else {
                let transfer = finite(area_time * flux[i])?;
                if flux[i] != 0.0 && transfer == 0.0 {
                    return Err(OsmoticNetworkError::StepUnderflow);
                }
                transfer
            };
            let (a, b) = (left.solute_moles()[i], right.solute_moles()[i]);
            let donor = if transfer >= 0.0 { a } else { b };
            // An absent species with no driving force must stay exactly absent.
            // Its zero transfer withdraws nothing; 0/0 is not a failed step.
            // A nonzero withdrawal from an empty donor still refuses below.
            let fraction = if transfer == 0.0 {
                0.0
            } else {
                transfer.abs() / donor
            };
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
            if !nonnegative(next_left[i]) || !nonnegative(next_right[i]) {
                return Err(OsmoticNetworkError::Invalid(
                    "nonnegative solution endpoint is unrepresentable",
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
                .material()
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
        let volumes = [left.solvent_volume_m3(), right.solvent_volume_m3()];
        // Old flux samples do not qualify a newly proposed physical state.
        link.qualify_sources(volumes, [&next_left[..n], &next_right[..n]], cancelled)?;
        poll(cancelled)?;
        cells[link.left]
            .solution
            .replace_inventories(volumes[0], &next_left[..n])
            .map_err(source)?;
        cells[link.right]
            .solution
            .replace_inventories(volumes[1], &next_right[..n])
            .map_err(source)?;
    }
    qualify_sources(cells, links, cancelled)
}

/// A later water or species connection can move an earlier source endpoint.
pub(super) fn qualify_sources(
    cells: &[OsmoticCompartment],
    links: &[SolutionDiffusionLink],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    for link in links.iter() {
        if link.material_sources().is_none() && link.electrical.is_none() {
            continue;
        }
        poll(cancelled)?;
        let (left, right) = (&cells[link.left].solution, &cells[link.right].solution);
        link.qualify_sources(
            [left.solvent_volume_m3(), right.solvent_volume_m3()],
            [left.solute_moles(), right.solute_moles()],
            cancelled,
        )?;
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
