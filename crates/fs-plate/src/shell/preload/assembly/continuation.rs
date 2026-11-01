//! Continuing static load paths on one immutable, originally unloaded assembly.
//!
//! Endpoints may change dead forces, tether reference tension and strut free
//! length. Geometry, sections, ties, supports, tether stiffness and strut EA
//! retain their original reference. The original Newton owner starts at the
//! last accepted deformation and evaluates the original unilateral axial laws
//! at each interpolated setting. This is quasi-static equilibrium, not elapsed
//! time, material history, or an actuation/dead-load work integral.

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
    tether_reference_tensions_n: Vec<f64>,
    strut_free_lengths_m: Vec<f64>,
}

impl ShellAssemblyLoadPath {
    /// Prepare the original sources and establish their first loaded equilibrium.
    ///
    /// The first solve has exactly the one-shot assembly's load-ramp semantics:
    /// dead loads and complete axial potentials grow from zero together. Later
    /// `advance` calls change only dead loads and retain full axial potentials;
    /// `advance_with_axial_settings` also changes their explicit settings.
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
            tether_reference_tensions_n: tethers.iter().map(|t| t.reference_tension_n).collect(),
            strut_free_lengths_m: struts.iter().map(|s| s.free_length_m).collect(),
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

    /// Accepted tensions at the ORIGINAL unloaded geometry [N], in tether order.
    /// Actual loaded tensions remain in `accepted().tether_responses`.
    pub fn tether_reference_tensions_n(&self) -> &[f64] {
        &self.tether_reference_tensions_n
    }

    /// Accepted uncompressed lengths [m], in original strut order. Original EA
    /// remains fixed; the original law computes each axial stiffness as EA/L.
    pub fn strut_free_lengths_m(&self) -> &[f64] {
        &self.strut_free_lengths_m
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

    /// Continue dead loads and the original axial settings together.
    ///
    /// Tether reference tension and strut free length interpolate linearly in
    /// the existing load-continuation fraction. Each intermediate strut uses
    /// `EA / interpolated_free_length`; interpolating stiffness or blending
    /// endpoint potentials would describe a different member. The original
    /// tension-only/compression-only branches, direction-updating forces and
    /// complete material/geometric tangents remain authoritative.
    ///
    /// Endpoints, ties, supports, shell reference geometry, tether stiffness and
    /// strut EA cannot change. Prepared connection maps share their immutable
    /// columns; no geometry or reference matrix is rebuilt. Every intermediate
    /// equilibrium retains the existing deformation, force, tangent-inertia and
    /// whole-call work gates. All loads, settings and state publish atomically.
    ///
    /// Repeating the exact settings and full loads is the existing exact no-op.
    /// Connection energy reports use the NEW law's original-geometry datum;
    /// their difference across changing settings is not an actuator-work audit.
    /// No physical duration, damping, actuator work, kinetics or changing
    /// material history is inferred from this numerical continuation.
    ///
    /// # Errors
    /// Original preload refusals, setting-count mismatch, nonfinite/negative
    /// tension, free length below the original 1e-9 m floor, or unrepresentable
    /// EA/L. Any refusal or cancellation preserves the complete accepted state.
    pub fn advance_with_axial_settings(
        &mut self,
        component_loads: &[&[f64]],
        tether_reference_tensions_n: &[f64],
        strut_free_lengths_m: &[f64],
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
            return Err("assembled axial setting path cancelled".into());
        }
        if tether_reference_tensions_n.len() != self.tether_reference_tensions_n.len()
            || strut_free_lengths_m.len() != self.strut_free_lengths_m.len()
        {
            return Err("one setting per original tether and strut is required".into());
        }
        let same = tether_reference_tensions_n
            .iter()
            .zip(&self.tether_reference_tensions_n)
            .chain(strut_free_lengths_m.iter().zip(&self.strut_free_lengths_m))
            .all(|(a, b)| a.to_bits() == b.to_bits());
        if same {
            return self.advance(component_loads, options, solve, cancelled);
        }
        let load = self
            .prepared
            .assembly
            .project_loads(component_loads)
            .map_err(|e| e.to_string())?;
        let mut connections = Vec::with_capacity(self.prepared.tethers.len());
        for (i, original) in self.prepared.tethers.iter().enumerate() {
            if cancelled() {
                return Err("assembled axial setting path cancelled".into());
            }
            connections.push(if i < self.prepared.tether_count {
                original.with_reference_tension(tether_reference_tensions_n[i])?
            } else {
                original.with_free_length(strut_free_lengths_m[i - self.prepared.tether_count])?
            });
        }
        let candidate = self.prepared.solve_with_connections(
            &connections,
            &load,
            Some(&self.retained_load),
            Some(&self.accepted.equilibrium.coordinates),
            options,
            solve,
            &mut cancelled,
        )?;
        for connection in &mut connections {
            connection.accept_setting();
        }
        let loads = component_loads.iter().map(|p| p.to_vec()).collect();
        let tensions = tether_reference_tensions_n.to_vec();
        let lengths = strut_free_lengths_m.to_vec();
        if cancelled() {
            return Err("assembled axial setting path cancelled".into());
        }
        self.accepted = candidate;
        self.component_loads = loads;
        self.retained_load = load;
        self.prepared.tethers = connections;
        self.tether_reference_tensions_n = tensions;
        self.strut_free_lengths_m = lengths;
        Ok(&self.accepted)
    }
}

#[cfg(test)]
mod axial_tests;
#[cfg(test)]
mod tests;
