//! Continuing static load paths on one immutable, originally unloaded assembly.
//!
//! A new endpoint changes dead forces only. Geometry, material sections, ties,
//! supports and axial connection laws retain their original reference. The
//! original Newton owner starts at the last accepted deformation; tethers and
//! struts stay fully active while dead loads interpolate between endpoints.
//! This is quasi-static equilibrium, not elapsed-time dynamics or load work.

use super::{
    AssemblyStrut, AssemblyTether, EdgeTie, FaceTie, PreloadOptions, PreloadProblem,
    PreloadedConnectedShellAssembly, PreparedAssembly, ShellPreloadComponent, TieBudget,
    TiedShellModel, prepare_assembly, validate_options,
};

/// A fixed physical assembly with a continuing, privately retained equilibrium.
///
/// Reference matrices, sparse component/connection maps, element laws and the
/// reference preconditioner are prepared once. Only the accepted tangent is
/// rebuilt at each nonlinear load increment. No displaced geometry is admitted
/// as stress free, and no accepted tangent is reused as an elastic reference.
pub struct ShellAssemblyLoadPath {
    prepared: PreparedAssembly,
    accepted: PreloadedConnectedShellAssembly,
    component_loads: Vec<Vec<f64>>,
    retained_load: Vec<f64>,
}

impl ShellAssemblyLoadPath {
    /// Prepare the original sources and establish their first loaded equilibrium.
    ///
    /// The first solve has exactly the one-shot assembly's load-ramp semantics:
    /// dead loads and complete axial potentials grow from zero together. Later
    /// `advance` calls change only dead loads and retain full axial potentials.
    /// Inputs are copied into the prepared physical kernels, so this path never
    /// borrows mutable caller geometry or load storage.
    ///
    /// # Errors
    /// The original assembly's geometry, source, connection, work, deformation,
    /// force-balance, tangent-inertia and cancellation refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        components: &[ShellPreloadComponent<'_>],
        edge_ties: &[EdgeTie],
        face_ties: &[FaceTie],
        tie_budget: TieBudget,
        tethers: &[AssemblyTether],
        struts: &[AssemblyStrut],
        options: PreloadOptions,
        solve: &mut impl FnMut(
            &PreloadProblem<'_>,
            Vec<f64>,
            usize,
        ) -> Result<(Vec<f64>, usize), String>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, String> {
        validate_options(options)?;
        let prepared = prepare_assembly(
            components,
            edge_ties,
            face_ties,
            tie_budget,
            tethers,
            struts,
            &mut cancelled,
        )?;
        let component_loads: Vec<_> = components.iter().map(|c| c.loads.to_vec()).collect();
        let retained_load = prepared
            .assembly
            .project_loads(
                &component_loads
                    .iter()
                    .map(Vec::as_slice)
                    .collect::<Vec<_>>(),
            )
            .map_err(|e| e.to_string())?;
        let accepted =
            prepared.solve(&retained_load, None, None, options, solve, &mut cancelled)?;
        if cancelled() {
            return Err("assembled shell load path cancelled".into());
        }
        Ok(Self {
            prepared,
            accepted,
            component_loads,
            retained_load,
        })
    }

    /// The complete last accepted loaded body, suitable for the existing
    /// monolithic vibration/response consumers. Its source cannot be edited.
    pub const fn accepted(&self) -> &PreloadedConnectedShellAssembly {
        &self.accepted
    }

    /// Original unloaded elastic pencil, mass and constraint maps. Axial
    /// connection potentials are additional prepared sources, not in this K.
    pub const fn reference_assembly(&self) -> &TiedShellModel {
        &self.prepared.assembly
    }

    /// Accepted full component forces [N] and moments [N m], before the same
    /// work-conjugate tie projection. Loads on supports are retained here.
    pub fn component_loads(&self) -> &[Vec<f64>] {
        &self.component_loads
    }

    /// Continue to complete new dead loads on the same original components.
    ///
    /// Each load increment uses `(1-lambda) f_old + lambda f_new` and the
    /// complete original connection potentials. All intermediate candidates
    /// pass the original physical-force, deformation and tangent-inertia gates.
    /// The original residual/Jacobian budget covers the whole requested change.
    /// A solver failure, invalid load, exhausted budget or final cancellation
    /// leaves the accepted loads, coordinates, tangent and port states unchanged.
    ///
    /// Repeating the exact full load vectors returns the existing body without
    /// another solve or changing its bits. Its iteration/evaluation fields still
    /// describe the last actual equilibrium solve; this is not recertification
    /// under newly supplied tolerances. Parameters are nevertheless validated.
    ///
    /// No elapsed time, inertia, damping, quasi-static work integral, material
    /// history, follower force or continuation through unstable branches is
    /// inferred. Connection opening/closing is permitted by the original laws;
    /// a subsequent vibration model must retain its selected active branches.
    ///
    /// # Errors
    /// Original preload refusals, mismatched/nonfinite load vectors, or a new
    /// deformation/force tolerance that refuses the accepted starting state.
    pub fn advance(
        &mut self,
        component_loads: &[&[f64]],
        options: PreloadOptions,
        solve: &mut impl FnMut(
            &PreloadProblem<'_>,
            Vec<f64>,
            usize,
        ) -> Result<(Vec<f64>, usize), String>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&PreloadedConnectedShellAssembly, String> {
        validate_options(options)?;
        if cancelled() {
            return Err("assembled shell load path cancelled".into());
        }
        let load = self
            .prepared
            .assembly
            .project_loads(component_loads)
            .map_err(|e| e.to_string())?;
        let same = component_loads.len() == self.component_loads.len()
            && component_loads
                .iter()
                .zip(&self.component_loads)
                .all(|(new, old)| new.iter().zip(old).all(|(a, b)| a.to_bits() == b.to_bits()));
        if same {
            if cancelled() {
                return Err("assembled shell load path cancelled".into());
            }
            return Ok(&self.accepted);
        }
        let candidate = self.prepared.solve(
            &load,
            Some(&self.retained_load),
            Some(&self.accepted.equilibrium.coordinates),
            options,
            solve,
            &mut cancelled,
        )?;
        let loads = component_loads.iter().map(|p| p.to_vec()).collect();
        if cancelled() {
            return Err("assembled shell load path cancelled".into());
        }
        self.accepted = candidate;
        self.component_loads = loads;
        self.retained_load = load;
        Ok(&self.accepted)
    }
}

#[cfg(test)]
mod tests;
