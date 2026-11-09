//! Prescribed translations and physical rotations on the original material shell.

use super::{
    C64, Cx, HarmonicPortControls, HarmonicResponseControls, HarmonicResponseError, PreparedFacet,
    PreparedSection, PreparedTensorPronyShellOperator, TensorPronyShell, TensorShellHarmonicError,
    admit_material_displacement, checkpoint, dense, invalid,
};
use fs_modal::harmonic::{PreparedSupportedHarmonicResponse, SupportedHarmonicResponse};

/// Reusable material-aware base excitation on complete original shell coordinates.
///
/// Each prescribed translation or physical axial rotation remains in the original
/// CST/DKT strain field and consistent inertia. Original homogeneous coordinates
/// removed during shell assembly stay fixed; an absent moving support refuses.
/// This is a particular periodic response, not a causal startup or actuator model.
pub struct PreparedSupportedTensorPronyShell<'a> {
    source: &'a TensorPronyShell,
    omega_rad_s: f64,
    materials: Vec<PreparedSection>,
    facets: Vec<PreparedFacet>,
    represented_dofs: Vec<usize>,
    free_dofs: Vec<usize>,
    support_dofs: Vec<usize>,
    response: PreparedSupportedHarmonicResponse,
}

impl TensorPronyShell {
    /// Prepare prescribed support motion using the complete dispersive shell.
    ///
    /// Support labels are original global DOFs `6*node + component`, where the
    /// components are `(u,v,w,theta_x,theta_y,theta_z)` in global axes. Retain
    /// every moving coordinate when assembling the source; this method applies
    /// their nonzero phasors after complete material and inertia assembly.
    /// Supply strictly increasing labels; empty and all-represented maps are
    /// admitted. All original layers must admit the positive physical frequency.
    ///
    /// The shared native solver factors only A_ff and solves
    /// A_ff*q_f = F_f - A_fc*q_c. Full cross inertia, membrane/bending coupling,
    /// every dispersive layer and the original lossless drilling term remain.
    /// No new mechanical solver, penalty support or artificial damping is used.
    ///
    /// Conservative dense payload is `344*n*n + 544*n + section_bytes`, where
    /// section_bytes is the original material/facet/resultant allowance used by
    /// `prepare_harmonic_operator`. Native sparse assembly retains the source's
    /// spatial budget. All represented coordinates count toward admission.
    pub fn prepare_supported_harmonic_response(
        &self,
        support_dofs: &[usize],
        omega_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedSupportedTensorPronyShell<'_>, TensorShellHarmonicError> {
        self.prepare_supported(support_dofs, omega_rad_s, controls, None, cx)
    }

    /// Add a complete physical effort/flow impedance before support elimination.
    ///
    /// Real B rows follow `represented_dofs()`, including moving support rows;
    /// opposing effort is Z*(B^T*v), with all self, mutual, reactive and
    /// nonreciprocal entries retained. The full pencil is
    /// K'(omega) - omega^2*M - i*K''(omega) - i*omega*B*Z*B^T.
    /// Only the free principal block needs to be invertible. The shared physical
    /// port owner separately admits its payload, projection and Hermitian-PSD
    /// impedance. The caller supplies source frequency validity and units.
    ///
    /// Boundary absorption, signed support work and original material heat stay
    /// separate; the impedance's owner determines radiation/heat/storage meaning.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_supported_harmonic_response_with_port_impedance(
        &self,
        support_dofs: &[usize],
        omega_rad_s: f64,
        force_map: &[f64],
        ports: usize,
        impedance: &[C64],
        controls: HarmonicResponseControls,
        port_controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedSupportedTensorPronyShell<'_>, TensorShellHarmonicError> {
        self.prepare_supported(
            support_dofs,
            omega_rad_s,
            controls,
            Some((force_map, ports, impedance, port_controls)),
            cx,
        )
    }

    fn prepare_supported(
        &self,
        support_dofs: &[usize],
        omega_rad_s: f64,
        controls: HarmonicResponseControls,
        boundary: Option<(&[f64], usize, &[C64], HarmonicPortControls)>,
        cx: &Cx<'_>,
    ) -> Result<PreparedSupportedTensorPronyShell<'_>, TensorShellHarmonicError> {
        let section_bytes = self.harmonic_section_payload(omega_rad_s, controls, cx)?;
        let n = self.model.free;
        if n == 0
            || support_dofs.len() > n
            || support_dofs.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(invalid(
                "ordered original support labels and a nonempty shell required",
            ));
        }
        let source_bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(24))
            .and_then(|v| n.checked_mul(32)?.checked_add(v))
            .and_then(|v| v.checked_add(section_bytes))
            .ok_or_else(|| invalid("supported shell source payload overflows"))?;
        let total = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(344))
            .and_then(|v| n.checked_mul(544)?.checked_add(v))
            .and_then(|v| v.checked_add(section_bytes))
            .ok_or_else(|| invalid("supported shell total payload overflows"))?;
        if total > controls.maximum_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "supported shell numerical payload bytes",
                requested: total,
                cap: controls.maximum_payload_bytes,
            }
            .into());
        }
        let mut support_rows = Vec::with_capacity(support_dofs.len());
        for &label in support_dofs {
            checkpoint(cx)?;
            let row = self
                .model
                .dof_map
                .get(label)
                .copied()
                .flatten()
                .ok_or_else(|| invalid("moving shell support is absent or already eliminated"))?;
            support_rows.push(row);
        }
        let mut represented_dofs = Vec::with_capacity(n);
        for (label, &row) in self.model.dof_map.iter().enumerate() {
            if label % 256 == 0 {
                checkpoint(cx)?;
            }
            if let Some(row) = row {
                if row != represented_dofs.len() {
                    return Err(invalid("original shell coordinate map is inconsistent"));
                }
                represented_dofs.push(label);
            }
        }
        if represented_dofs.len() != n {
            return Err(invalid("original shell coordinate map is incomplete"));
        }
        let PreparedTensorPronyShellOperator {
            materials,
            facets,
            stiffness,
            damping,
            ..
        } = self.prepare_harmonic_operator(omega_rad_s, controls, cx)?;
        let mass = dense(&self.model.m, n, cx)?;
        let storage = dense(&stiffness, n, cx)?;
        let viscosity = dense(&damping, n, cx)?;
        drop(stiffness);
        drop(damping);
        let solve_controls = HarmonicResponseControls {
            maximum_payload_bytes: controls.maximum_payload_bytes - source_bytes,
            ..controls
        };
        let response = if let Some((b, ports, z, load_controls)) = boundary {
            PreparedSupportedHarmonicResponse::try_with_port_impedance(
                &mass,
                &viscosity,
                &storage,
                n,
                &support_rows,
                b,
                ports,
                z,
                omega_rad_s,
                solve_controls,
                load_controls,
                cx,
            )?
        } else {
            PreparedSupportedHarmonicResponse::try_new(
                &mass,
                &viscosity,
                &storage,
                n,
                &support_rows,
                omega_rad_s,
                solve_controls,
                cx,
            )?
        };
        let free_dofs = response
            .free_dofs()
            .iter()
            .map(|&i| represented_dofs[i])
            .collect();
        checkpoint(cx)?;
        Ok(PreparedSupportedTensorPronyShell {
            source: self,
            omega_rad_s,
            materials,
            facets,
            represented_dofs,
            free_dofs,
            support_dofs: support_dofs.to_vec(),
            response,
        })
    }
}

impl<'a> PreparedSupportedTensorPronyShell<'a> {
    pub const fn source(&self) -> &'a TensorPronyShell {
        self.source
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Original labels of all applied-force and full displacement/velocity rows.
    pub fn represented_dofs(&self) -> &[usize] {
        &self.represented_dofs
    }
    /// Original labels of coordinates solved by their physical force balance.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    /// Original labels of supplied support phasors and returned reactions.
    pub fn support_dofs(&self) -> &[usize] {
        &self.support_dofs
    }

    /// Solve peak force/torque and support translation/rotation under exp(-i omega t).
    ///
    /// Inputs follow `represented_dofs()` and `support_dofs()` respectively.
    /// The result retains every represented coordinate, including the supports;
    /// its reaction is the force/torque exerted BY the support ON the structure.
    /// `average_viscous_power_w` is only the original tensor-Prony/Kelvin heat.
    /// All layers admit their complete strain cycle at every facet vertex and
    /// layer face before publication. No factor, source or history is mutated
    /// on failure/cancellation. The returned complete displacement can enter
    /// this same source's coherent material spectrum observer directly.
    pub fn solve(
        &self,
        applied_force: &[C64],
        support_displacement: &[C64],
        cx: &Cx<'_>,
    ) -> Result<SupportedHarmonicResponse, TensorShellHarmonicError> {
        let response = self
            .response
            .solve(applied_force, support_displacement, cx)?;
        admit_material_displacement(
            self.source,
            &self.materials,
            &self.facets,
            &response.displacement,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(response)
    }
}
