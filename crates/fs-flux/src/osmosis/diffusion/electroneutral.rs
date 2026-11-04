//! Closed ionic networks with globally compatible diffusion potentials.
//!
//! At a frozen physical state the original Maxwell--Stefan owner supplies
//! `q = q_chemical - G*(psi_right-psi_left)`, where psi=F*phi/(RT).
//! Charge continuity gives a weighted graph Laplacian. One gauge per conducting
//! component removes its nullspace; the original checked-CG owner solves the
//! admitted SPD problem. Actual combined-force fluxes, not the affine estimate,
//! must then satisfy every node's absolute charge-rate budget.
//!
//! All species transfers are simultaneous. A pair may gain chemical free
//! energy while transferring electrical work to another pair, so the finite
//! step checks the GLOBAL original mixing free energy. The existing temporal
//! controller, source qualification and atomic whole-window publication remain
//! in charge. This is a bulk electroneutral limit, not a Poisson/double-layer,
//! electrode, imposed-current or caloric energy model. The conjugate current
//! and entropy structure follows <https://arxiv.org/abs/2211.13505>, section 6.

mod stage;

use super::{
    Arc, OsmoticCompartment, OsmoticNetwork, OsmoticNetworkError, SolutionDiffusionLink,
    ZeroCurrentDiffusion, diffusion_source, positive,
};
use crate::osmosis::{OsmoticMembrane, coupled::CoupledSolutionControls};

/// Accuracy and complete-window work limits for the global electrical closure.
/// These supplement the original amount, energy and temporal error budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectroneutralDiffusionControls {
    /// Recomputed Euclidean relative graph residual, in `(0, 1e-3]`.
    /// Actual node charge rates are checked independently after native replay.
    pub potential_relative_tolerance: f64,
    /// Absolute `sum_outgoing area*sum_i(z_i*N_i)` at EVERY node [mol/s].
    /// Charge is in moles of elementary-charge equivalents, not coulombs.
    pub maximum_node_charge_rate_mol_s: f64,
    /// Native CG iterations per physical potential solve, including corrections.
    pub maximum_iterations_per_solve: usize,
    /// Total CG iterations over the COMPLETE requested advance, including
    /// discarded coarse/fine probes and temporal rejections.
    pub maximum_potential_iterations: usize,
    /// Original face preparations and combined-force replays over the complete
    /// advance. A preparation uses at most two native friction solves; a replay
    /// uses one. Original source qualification queries are additional work.
    pub maximum_face_evaluations: usize,
}

impl Default for ElectroneutralDiffusionControls {
    fn default() -> Self {
        Self {
            potential_relative_tolerance: 1e-10,
            maximum_node_charge_rate_mol_s: 1e-16,
            maximum_iterations_per_solve: 1024,
            maximum_potential_iterations: 1_000_000,
            maximum_face_evaluations: 1_000_000,
        }
    }
}

impl ElectroneutralDiffusionControls {
    fn admit(self) -> Result<(), OsmoticNetworkError> {
        if !positive(self.potential_relative_tolerance)
            || self.potential_relative_tolerance > 1e-3
            || !positive(self.maximum_node_charge_rate_mol_s)
            || self.maximum_iterations_per_solve == 0
            || self.maximum_iterations_per_solve > 65_536
            || self.maximum_potential_iterations == 0
            || self.maximum_face_evaluations == 0
        {
            return Err(OsmoticNetworkError::Invalid(
                "electroneutral diffusion needs explicit charge accuracy and bounded native potential/face work",
            ));
        }
        Ok(())
    }
}

/// Work and physical diagnostics over all probes in a successful advance.
/// Maxima include valid physical stages discarded by temporal error control;
/// they do not purport to be a final-state potential or edge-current snapshot.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ElectroneutralDiffusionReport {
    pub potential_solves: usize,
    pub potential_iterations: usize,
    pub face_evaluations: usize,
    /// Largest actual signed-node-balance magnitude after original replay.
    pub maximum_node_charge_rate_mol_s: f64,
    /// Largest actual edge-current magnitude, including circulating currents.
    pub maximum_edge_charge_rate_mol_s: f64,
    /// Largest |F*(phi_right-phi_left)/(RT)| on a conducting molecular edge.
    pub maximum_dimensionless_potential_difference: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct ElectroneutralDiffusion {
    controls: ElectroneutralDiffusionControls,
    report: ElectroneutralDiffusionReport,
    parents: Vec<usize>,
    rhs: Vec<f64>,
    diagonal: Vec<f64>,
    pub(super) potentials: Vec<f64>,
    charge_rate: Vec<f64>,
    charge_correction: Vec<f64>,
    // During assembly these are paired net transfers. After admission the
    // first array holds complete candidate inventories for simultaneous commit.
    amounts: Vec<f64>,
    amount_correction: Vec<f64>,
    withdrawn: Vec<f64>,
    withdrawal_correction: Vec<f64>,
}

impl OsmoticNetwork {
    /// Prepare general closed ionic graphs, including cycles, with correct
    /// neutral material-source faces from the FIRST source query. Optional
    /// water paths use the original bare, fully retaining membrane owner.
    /// Empty water paths give a fixed-volume ionic diffusion network.
    #[allow(clippy::too_many_arguments)]
    pub fn with_electroneutral_solution_transport(
        species: Vec<String>,
        cells: Vec<OsmoticCompartment>,
        membranes: Vec<OsmoticMembrane>,
        mut links: Vec<SolutionDiffusionLink>,
        physical: CoupledSolutionControls,
        ionic: ZeroCurrentDiffusion,
        electrical: ElectroneutralDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        electrical.admit()?;
        let ionic = Arc::new(ionic);
        for link in &mut links {
            link.electrical = Some(Arc::clone(&ionic));
        }
        Self::with_coupled_solution_transport(
            species,
            cells,
            membranes,
            links,
            physical,
            maximum_workspace_bytes,
        )?
        .prepare_electroneutral(ionic, electrical, maximum_workspace_bytes)
    }

    /// Select a global floating potential for an unadvanced molecular network.
    /// The policy supplies original valences and an absolute inventory-charge
    /// budget; its per-face zero-current tolerance is NOT imposed on cycles.
    /// Every component is electrically isolated and may carry loop currents.
    /// No external voltage/current port or independent ionic leak is inferred.
    ///
    /// This mode uses simultaneous explicit molecular stages, controlled by the
    /// existing full/two-half timeline. Independent-transfer and pairwise salt
    /// implicit modes are incompatible with this global charge constraint.
    ///
    /// Storage admission includes retained arrays, transient bound face
    /// responses, resolved binary matrices and a conservative allowance for
    /// the original native CG numerical vectors/history. Native solves and
    /// original source queries allocate transient storage. As with sourced
    /// diffusion, this is a numerical payload estimate, not an allocator/RSS
    /// bound; original query/receipt strings and allocator overhead are extra.
    pub fn with_electroneutral_diffusion(
        self,
        ionic: ZeroCurrentDiffusion,
        controls: ElectroneutralDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        if self.diffusion.iter().any(|link| link.electrical.is_some()) {
            return Err(OsmoticNetworkError::Invalid(
                "ionic assignment is already selected",
            ));
        }
        self.prepare_electroneutral(Arc::new(ionic), controls, maximum_workspace_bytes)
    }

    fn prepare_electroneutral(
        mut self,
        ionic: Arc<ZeroCurrentDiffusion>,
        controls: ElectroneutralDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        controls.admit()?;
        if self.time_s != 0.0
            || self.diffusion.is_empty()
            || self.diffusion_controls.is_none()
            || self.implicit_diffusion.is_some()
            || self.implicit_multicomponent.is_some()
            || self.implicit_salt.is_some()
            || self.electroneutral_diffusion.is_some()
            || ionic.valences().len() != self.species.len()
        {
            return Err(OsmoticNetworkError::Invalid(
                "global electroneutral transport needs an unadvanced explicit molecular network and one original valence per solute",
            ));
        }
        let nodes = self.cells.len();
        let entries = nodes
            .checked_mul(self.species.len())
            .ok_or(OsmoticNetworkError::Invalid(
                "electroneutral inventory storage overflow",
            ))?;
        let retained = nodes
            .checked_mul(core::mem::size_of::<usize>() + 5 * core::mem::size_of::<f64>())
            .and_then(|bytes| {
                bytes.checked_add(entries.checked_mul(4 * core::mem::size_of::<f64>())?)
            })
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<ElectroneutralDiffusion>()))
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<ZeroCurrentDiffusion>()))
            .and_then(|bytes| bytes.checked_add(ionic.heap_bytes()))
            .and_then(|bytes| bytes.checked_add(2 * core::mem::size_of::<usize>()));
        let transient = stage::workspace_bytes(nodes, &self.diffusion, controls)?;
        let required = retained
            .and_then(|bytes| bytes.checked_add(transient))
            .and_then(|bytes| bytes.checked_add(self.workspace_bytes))
            .ok_or(OsmoticNetworkError::Invalid(
                "electroneutral workspace byte count overflow",
            ))?;
        if required > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        for cell in &self.cells {
            ionic
                .admit_inventories(cell.solution.solute_moles())
                .map_err(diffusion_source)?;
        }
        for link in &mut self.diffusion {
            link.electrical = Some(Arc::clone(&ionic));
        }
        super::qualify_sources(&self.cells, &self.diffusion, &mut || false)?;
        self.electroneutral_diffusion = Some(ElectroneutralDiffusion {
            controls,
            report: ElectroneutralDiffusionReport::default(),
            parents: filled(nodes, 0usize)?,
            rhs: filled(nodes, 0.0)?,
            diagonal: filled(nodes, 0.0)?,
            potentials: filled(nodes, 0.0)?,
            charge_rate: filled(nodes, 0.0)?,
            charge_correction: filled(nodes, 0.0)?,
            amounts: filled(entries, 0.0)?,
            amount_correction: filled(entries, 0.0)?,
            withdrawn: filled(entries, 0.0)?,
            withdrawal_correction: filled(entries, 0.0)?,
        });
        self.workspace_bytes = required;
        Ok(self)
    }

    pub fn electroneutral_diffusion_controls(&self) -> Option<ElectroneutralDiffusionControls> {
        self.electroneutral_diffusion
            .as_ref()
            .map(|state| state.controls)
    }
}

impl ElectroneutralDiffusion {
    pub(crate) fn begin_window(&mut self) {
        self.report = ElectroneutralDiffusionReport::default();
    }
    pub(crate) fn report(&self) -> ElectroneutralDiffusionReport {
        self.report
    }
}

fn filled<T: Clone>(length: usize, value: T) -> Result<Vec<T>, OsmoticNetworkError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(length)
        .map_err(|_| OsmoticNetworkError::Invalid("electroneutral workspace allocation failed"))?;
    result.resize(length, value);
    Ok(result)
}
