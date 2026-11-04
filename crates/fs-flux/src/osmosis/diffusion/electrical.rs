//! Charge-conserving molecular paths on an electrically isolated forest.

use super::{
    Arc, OsmoticCompartment, OsmoticNetwork, OsmoticNetworkError, SolutionDiffusionLink,
    ZeroCurrentDiffusion, diffusion_source,
};
use crate::osmosis::{OsmoticMembrane, coupled::CoupledSolutionControls};

impl OsmoticNetwork {
    /// Prepare ionic molecular paths and optional fully retaining water paths
    /// with the correct neutral face from the first network source query.
    /// An empty water path set gives fixed-volume ionic diffusion. This entry
    /// point supports sources whose joint domain excludes off-neutral states;
    /// preparing a neutral molecular network first could query the wrong face.
    #[allow(clippy::too_many_arguments)]
    pub fn with_zero_current_solution_transport(
        species: Vec<String>,
        cells: Vec<OsmoticCompartment>,
        membranes: Vec<OsmoticMembrane>,
        mut links: Vec<SolutionDiffusionLink>,
        controls: CoupledSolutionControls,
        electrical: ZeroCurrentDiffusion,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        let electrical = Arc::new(electrical);
        for link in &mut links {
            link.electrical = Some(Arc::clone(&electrical));
        }
        Self::with_coupled_solution_transport(
            species,
            cells,
            membranes,
            links,
            controls,
            maximum_workspace_bytes,
        )?
        .prepare_zero_current(electrical, maximum_workspace_bytes)
    }

    /// Bind explicit ion valences to ALL molecular paths before first advance.
    /// Compartments must already be neutral within the absolute molar budget.
    /// Existing coupled water paths retain every ion; no independent ionic leak
    /// or solvent-drag law is inferred. Constant-data advances allocate nothing.
    ///
    /// Conducting molecular links must form a forest. Independent zero-current
    /// face potentials can then be integrated to a potential on each component;
    /// a cycle needs a global current/potential solve that this option excludes.
    /// Zero-area closed links impose no electrical constraint.
    ///
    /// Uses the original explicit conservative sweep and adaptive full/two-half
    /// timeline. The current independent-transfer implicit solvers are refused:
    /// their nonlinear trial coordinates do not preserve electroneutrality.
    /// The byte budget includes the temporary forest check as well as the one
    /// shared retained ionic assignment, never one copy per molecular link.
    pub fn with_zero_current_diffusion(
        self,
        electrical: ZeroCurrentDiffusion,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        if self.zero_current_diffusion().is_some() {
            return Err(OsmoticNetworkError::Invalid(
                "ionic assignment is already selected",
            ));
        }
        self.prepare_zero_current(Arc::new(electrical), maximum_workspace_bytes)
    }

    fn prepare_zero_current(
        mut self,
        electrical: Arc<ZeroCurrentDiffusion>,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        if self.time_s != 0.0
            || self.diffusion.is_empty()
            || self.diffusion_controls.is_none()
            || electrical.valences().len() != self.species.len()
            || self.implicit_diffusion.is_some()
            || self.implicit_multicomponent.is_some()
        {
            return Err(OsmoticNetworkError::Invalid(
                "zero-current transport needs an unadvanced explicit molecular network and one original valence per solute",
            ));
        }
        let required = self
            .workspace_bytes
            .checked_add(core::mem::size_of::<ZeroCurrentDiffusion>())
            .and_then(|bytes| bytes.checked_add(electrical.heap_bytes()))
            .and_then(|bytes| bytes.checked_add(2 * core::mem::size_of::<usize>()))
            .ok_or(OsmoticNetworkError::Invalid("ionic storage size overflow"))?;
        let preparation = self
            .cells
            .len()
            .checked_mul(core::mem::size_of::<usize>())
            .and_then(|bytes| bytes.checked_add(required))
            .ok_or(OsmoticNetworkError::Invalid(
                "ionic forest workspace size overflow",
            ))?;
        if preparation > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: preparation,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        for cell in &self.cells {
            electrical
                .admit_inventories(cell.solution.solute_moles())
                .map_err(diffusion_source)?;
        }
        let mut components = Vec::new();
        components
            .try_reserve_exact(self.cells.len())
            .map_err(|_| {
                OsmoticNetworkError::Invalid("ionic forest workspace allocation failed")
            })?;
        components.extend(0..self.cells.len());
        for link in &self.diffusion {
            if link.area_m2 == 0.0 {
                continue;
            }
            let left = root(&mut components, link.left);
            let right = root(&mut components, link.right);
            if left == right {
                return Err(OsmoticNetworkError::Invalid(
                    "cyclic ionic paths require a global electrochemical current/potential solve",
                ));
            }
            components[left] = right;
        }
        for link in &mut self.diffusion {
            link.electrical = Some(Arc::clone(&electrical));
        }
        super::qualify_sources(&self.cells, &self.diffusion, &mut || false)?;
        self.workspace_bytes = required;
        Ok(self)
    }

    /// The one common immutable ionic assignment, if selected.
    pub fn zero_current_diffusion(&self) -> Option<&ZeroCurrentDiffusion> {
        self.diffusion
            .first()
            .and_then(|link| link.zero_current_diffusion())
    }
}

fn root(parents: &mut [usize], mut node: usize) -> usize {
    // Path halving avoids a recursive traversal and its stack-depth limit.
    while parents[node] != node {
        parents[node] = parents[parents[node]];
        node = parents[node];
    }
    node
}
