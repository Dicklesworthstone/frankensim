//! Stress and support reactions from the last accepted inertial endpoint.
use super::*;
use crate::linear3::TetThermalStrainState;

/// Actual elastic stress and full-node equation terms at a single accepted time.
/// `inertia_force_n` is +M*a on the equation's left side, not -M*a.
/// Reactions are forces exerted by the admitted supports ON the solid:
/// R = M*a + C*v + integral(B' sigma) - f_external at constrained components.
/// Free residuals are measured from those full rows, never reported as zero by
/// construction. These are discrete fields, not continuum error estimates.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientResponse3 {
    pub time_s: f64,
    pub step: usize,
    pub stress: TetThermalStressReport,
    pub inertia_force_n: Vec<[f64; 3]>,
    pub damping_force_n: Vec<[f64; 3]>,
    pub external_force_n: Vec<[f64; 3]>,
    pub support_reactions_n: Vec<[f64; 3]>,
    pub free_force_residual_n: f64,
    /// 1/2 v'Mv + 1/2 integral (strain-strain_free)' C (strain-strain_free).
    /// This differs from the stepper's unshifted quadratic energy when heated.
    /// No thermal storage or damping-heat field is inferred.
    pub kinetic_plus_elastic_energy_j: f64,
}

#[derive(Clone)]
pub(super) enum FrozenThermal3 {
    Uniform(TetThermalStrainState),
    PerElement(Arc<[TetThermalStrainState]>),
}
impl FrozenThermal3 {
    pub(super) fn capture(field: TetThermalStrainField<'_>, cx: &Cx<'_>)
        -> Result<Self, TransientError3> {
        check(cx)?;
        Ok(match field {
            TetThermalStrainField::Uniform(v) => Self::Uniform(*v),
            TetThermalStrainField::PerElement(values) => {
                let mut saved = Vec::with_capacity(values.len());
                for value in values { check(cx)?; saved.push(*value); }
                Self::PerElement(Arc::from(saved))
            }
        })
    }
    pub(super) fn field(&self) -> TetThermalStrainField<'_> {
        match self {
            Self::Uniform(v) => TetThermalStrainField::Uniform(v),
            Self::PerElement(v) => TetThermalStrainField::PerElement(v),
        }
    }
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Frozen initial-stress pencil from this actual accepted inertial state,
    /// using an explicitly selected incremental material tangent.
    ///
    /// Native response recovery retains the run's exact motion, original
    /// material/thermal field and supports; no quasistatic replacement solve
    /// or external stress array enters the path. The replacement tangent must
    /// preserve each reference density. Signed eigenvalues remain signed.
    ///
    /// This is a local conservative snapshot, not a normal-mode interpretation
    /// of a generally time-varying trajectory. Damping, moving-support work,
    /// follower loads and further coupled-state evolution are not added here.
    pub fn assemble_prestressed_constitutive_tangent(
        &self, constitutive_materials: crate::linear3::TetMaterialField<'_>,
        mesh_frame: ContentHash, maximum_strain: f64, cx: &Cx<'_>,
    ) -> Result<crate::prestress::PrestressedTetAssembly, TransientError3> {
        let response = self.recover_response(mesh_frame, maximum_strain, cx)?;
        self.problem.assemble_constitutive_tangent_from_stress(
            response.stress, constitutive_materials, cx).map_err(TransientError3::from)
    }

    /// Recover stress, inertia, damping and physical support forces WITHOUT a
    /// static solve, new constraints, force resampling or global reassembly.
    /// The accepted full external load (including constrained-node loads) and
    /// the exact frozen thermal field are retained by the run itself. Callers
    /// cannot accidentally combine yesterday's load with today's displacement.
    /// The existing stress report supports its normal frame-transform methods.
    /// Finite total AND elastic Mandel strains must meet `maximum_strain`.
    /// Recovery/refusal/cancellation never changes the accepted trajectory.
    pub fn recover_response(&self, mesh_frame: ContentHash, maximum_strain: f64,
        cx: &Cx<'_>) -> Result<TransientResponse3, TransientError3> {
        check(cx)?;
        let motion = self.nodal_motion(cx)?;
        let fields = self.problem.recover_kinematic_response(&motion,
            self.thermal_strains.as_ref().map(FrozenThermal3::field), self.state.time_s,
            self.state.step, self.config.mass_damping_per_s, self.config.stiffness_damping_s,
            mesh_frame, maximum_strain, cx)?;
        let mut external = vec![[0.0; 3]; self.geometry.nodes.len()];
        for &(d, f) in &self.full_force { check(cx)?; external[d/3][d%3] = f; }
        let mut imbalance = external.clone();
        for (n, row) in imbalance.iter_mut().enumerate() {
            check(cx)?;
            for c in 0..3 {
                row[c] = number(fields.inertia_n[n][c] + fields.damping_n[n][c]
                    + fields.stress.nodal_internal_force_n()[n][c] - external[n][c])?;
            }
        }
        let mut reactions = vec![[0.0; 3]; self.geometry.nodes.len()];
        for &d in self.problem.fixed_dofs {
            check(cx)?; reactions[d/3][d%3] = imbalance[d/3][d%3];
        }
        let mut residual = 0.0_f64;
        for &d in self.free_dofs() {
            check(cx)?; residual = number(residual.hypot(imbalance[d/3][d%3]))?;
        }
        let energy = number(self.state.kinetic_energy_j + fields.stress.elastic_energy_j())?;
        check(cx)?;
        Ok(TransientResponse3 { time_s: self.state.time_s, step: self.state.step,
            stress: fields.stress, inertia_force_n: fields.inertia_n,
            damping_force_n: fields.damping_n, external_force_n: external,
            support_reactions_n: reactions, free_force_residual_n: residual,
            kinetic_plus_elastic_energy_j: energy })
    }
}
