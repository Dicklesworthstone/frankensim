//! Conservative native force/pressure impulses returned to transported cells.

use super::{
    BoundaryPressure, CancelGate, ChtError, DensityForceControls, DensityForceInput,
    DensityMomentumReconstruction, DensityMomentumReconstructionControls,
    DensityMomentumReconstructionReport, DensityProjection, DensityProjectionControls,
    DensityResistanceReport, FlowFace, FlowField, ForcedDensityProjection, PressureReadback,
    ReservoirPressure, ResistiveDensityProjection, StaggeredBoundaryPressure, Sum, VoxelDomain,
    face_dual_mass, face_resistance_kg_s, finite, finite_positive, nonnegative, poll, poll_every,
    refuse,
};
use fs_math::dd::Dd;

/// Original extensive material state and authored sources for one fresh kick.
#[derive(Debug, Clone, Copy)]
pub struct MaterialDensityMomentumInput<'a> {
    /// Authoritative fluid mass per original cell, kg; zero in solids.
    pub mass_kg: &'a [f64],
    /// Transported extensive momentum in world axes, kg m/s; zero in solids.
    pub momentum_kg_m_s: &'a [[f64; 3]],
    /// Uniform physical acceleration, m/s^2, including gravity if present.
    pub acceleration_m_s2: [f64; 3],
    /// Additional body force, N/m^3; excludes pressure, and is zero in solids.
    pub force_density_n_m3: &'a [[f64; 3]],
    /// Authored net outward volume flux, m^3/s, for this fresh correction.
    /// Its associated material mass/enthalpy/EOS evolution remains caller-owned.
    pub target_outflow_m3_s: &'a [f64],
}

/// Original native stage controls and independent cell-feedback readback gates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialDensityMomentumControls {
    pub reconstruction: DensityMomentumReconstructionControls,
    pub forcing: DensityForceControls,
    pub projection: DensityProjectionControls,
    /// Sum of absolute local and global impulse defects, kg m/s. Includes
    /// native Q/A readback and actual stored cell-momentum increments.
    pub maximum_feedback_momentum_defect_kg_m_s: f64,
    /// Corresponding angular-impulse defect about the domain origin, kg m^2/s.
    pub maximum_feedback_angular_impulse_defect_kg_m2_s: f64,
    /// Absolute mass-representation work plus absolute local/global energy
    /// defects, J. Naming numerical representation work never waives this gate.
    pub maximum_feedback_energy_defect_j: f64,
}

/// Actual returned cell momentum against the independently read native stages.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialDensityMomentumReport {
    /// Held reconstruction reaction plus fresh body/pressure/resistance impulse,
    /// kg m/s. Resistance is present only in the resistance-aware entry points.
    /// The component accounts remain in `reconstruction` and `evolution`.
    pub face_impulse_kg_m_s: [f64; 3],
    /// Actual `sum(P_returned-P_input)`, kg m/s.
    pub cell_momentum_change_kg_m_s: [f64; 3],
    pub absolute_momentum_defect_kg_m_s: f64,
    /// `sum(r_face cross I_face)`, about the original domain origin, kg m^2/s.
    pub face_angular_impulse_kg_m2_s: [f64; 3],
    /// Actual `sum(r_cell cross (P_returned-P_input))`, kg m^2/s.
    pub cell_angular_momentum_change_kg_m2_s: [f64; 3],
    pub absolute_angular_impulse_defect_kg_m2_s: f64,
    pub kinetic_energy_before_j: f64,
    pub kinetic_energy_after_j: f64,
    /// Cell kinetic change evaluated from actual momentum increments, J.
    pub kinetic_energy_change_j: f64,
    /// Native unconstrained-face work `sum((p/Mn+J/(2Mn))*J)`, J, where
    /// p=sum(P_axis/2), Mn is original pressure inertia, and J includes the
    /// held-face reconstruction reaction. This is not another physical force.
    pub native_impulse_work_j: f64,
    /// Nonnegative `sum(M_cell*(delta_v_plus-delta_v_minus)^2/8)`, J, with
    /// delta_v=J/sum(M/2). Numerical scatter loss never heats material.
    pub numerical_scatter_loss_j: f64,
    /// Signed `sum((Mn-Ma)*(p/Ma+J/(2Ma))*(J/Mn))`, J, Ma=sum(M/2).
    /// Native M/V representation is numerical; its absolute work is budgeted.
    pub mass_representation_work_j: f64,
    /// `cell_delta_K + scatter_loss - native_impulse_work
    /// - mass_representation_work`, evaluated in centered local accounts, J.
    pub feedback_energy_defect_j: f64,
    /// `cell_delta_K + scatter_loss + reconstruction_constraint_loss
    /// + pressure_projection_loss + force_pressure_split_work
    /// - reconstruction_reservoir_work - body_endpoint_work
    /// - force_constraint_endpoint_work - pressure_source_work
    /// + fixed_boundary_pressure_work + open_boundary_pressure_work
    /// - mass_representation_work`, J.
    /// Initial reconstruction Jensen loss cancels; it is NOT debited again.
    /// For the implicit resistance entries, replace pressure projection loss
    /// and force/pressure split work by native backward-Euler loss plus physical
    /// resistance heat. Their body and constraint work use final face velocity.
    pub composed_energy_defect_j: f64,
    pub absolute_energy_defect_j: f64,
}

/// One accepted native force/pressure kick and conservative cell feedback.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialDensityMomentumStep {
    /// Updated extensive momentum, suitable for the next actual material route.
    pub momentum_kg_m_s: Vec<[f64; 3]>,
    /// The exact admitted M/V values used by all native stages, kg/m^3.
    pub density_kg_m3: Vec<f64>,
    pub reconstruction: DensityMomentumReconstructionReport,
    /// Actual native pressure-corrected face field and original physical reports.
    pub evolution: ForcedDensityProjection,
    pub report: MaterialDensityMomentumReport,
}

/// Actual-density implicit stationary resistance and conservative cell feedback.
#[derive(Debug, Clone, PartialEq)]
pub struct ResistiveMaterialDensityMomentumStep {
    pub momentum_kg_m_s: Vec<[f64; 3]>,
    /// The original current M/V, never the resistance-modified pressure metric.
    pub density_kg_m3: Vec<f64>,
    pub reconstruction: DensityMomentumReconstructionReport,
    pub evolution: ResistiveDensityProjection,
    pub report: MaterialDensityMomentumReport,
}

impl FlowField {
    /// Reconstruct transported momentum, apply one native force/pressure kick,
    /// and return its impulses to the ORIGINAL cell momenta by half-cell mass.
    /// `self` supplies held boundary fluxes; other starting fluxes are rebuilt.
    /// All stages use the same M/V and original geometry. Every fluid half-cell
    /// participates, including at walls and either orientation of a solid face.
    ///
    /// A face impulse J is shared with each contributor in proportion to its
    /// authoritative M/2. The updated momentum is original P plus these shares;
    /// it is not an inverse reconstruction or an average of final velocities.
    /// Zero impulses retain original sub-face momentum exactly. Held-face
    /// reconstruction reactions are included, so supported momentum cannot
    /// silently reappear on the next transport step. Normal impulse transfer
    /// preserves both linear and orbital angular impulse on this rectilinear
    /// grid, subject to the explicit readback gates.
    ///
    /// Fresh stages are computed here; mutable public reports are never accepted
    /// as proof that arbitrary fields belong together. Actual starting/final
    /// face fluxes are checked against half momenta, authored force, returned
    /// pressure and volume constraints. Representational Q/A roundoff is a
    /// numerical defect, never an extra force scattered into the cells.
    /// Independently public mass/momentum arrays are validated as supplied data;
    /// no private association with a previous material result is invented.
    ///
    /// The volume target must be for a FRESH pressure correction. Do not consume
    /// the pressure impulse of an already evolved phase-volume result again.
    /// The enclosing owner still supplies compatible transport, EOS/enthalpy,
    /// phase mobility and source work. This bounded split component does not
    /// evolve staggered advective inertia or define a complete Navier-Stokes
    /// timestep. Reconstructing the returned P need not reproduce the final
    /// divergence-constrained face field. Numerical scatter/constraint/projection
    /// losses and mass-representation work never alter material enthalpy.
    ///
    /// # Errors
    /// Any native stage refusal, invalid feedback budget, mismatched field,
    /// excessive independent momentum/angular/work defect, or cancellation.
    /// No input or accepted material state changes on refusal.
    pub fn force_project_material_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: MaterialDensityMomentumInput<'_>,
        controls: MaterialDensityMomentumControls,
        gate: &CancelGate,
    ) -> Result<MaterialDensityMomentumStep, ChtError> {
        self.force_project_material_momentum_with_optional_pressure(
            domain, faces, None, input, controls, gate,
        )
    }

    /// Return a fresh body-force and reservoir-pressure impulse to the original
    /// material momentum. The borrowed profile has the same native face layout,
    /// open-face admission and physical half-cell inertia as
    /// [`Self::force_and_project_density_with_pressure_profile`].
    ///
    /// Every open-face impulse uses the actual cell-to-reservoir pressure
    /// difference. Its half-cell share reaches material momentum once, including
    /// hydrostatic support and either flow direction. The existing independent
    /// angular/energy balances include the original reservoir work; numerical
    /// scatter and projection losses remain separate. This method supplies no
    /// reservoir mass, enthalpy, constitutive law or additional transport solve.
    ///
    /// # Errors
    /// Invalid pressure data or any original reconstruction, forcing, projection,
    /// feedback or cancellation refusal leaves all supplied state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn force_project_material_momentum_with_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        input: MaterialDensityMomentumInput<'_>,
        controls: MaterialDensityMomentumControls,
        gate: &CancelGate,
    ) -> Result<MaterialDensityMomentumStep, ChtError> {
        self.force_project_material_momentum_with_reservoir_pressure(
            domain, faces, profile, input, controls, gate,
        )
    }

    /// Return one native force/pressure impulse to original material momentum
    /// using either six optional plane pressures in Pa or a borrowed spatial
    /// profile. Uniform reservoirs need no staggered pressure buffers. Both
    /// inputs reuse the existing component-relative pressure and work readback.
    /// Reservoir work remains mechanical and does not alter material enthalpy.
    ///
    /// # Errors
    /// Nonfinite pressure, pressure on a held plane, invalid profile, or any
    /// refusal of [`Self::force_project_material_momentum_with_pressure_profile`].
    #[allow(clippy::too_many_arguments)]
    pub fn force_project_material_momentum_with_reservoir_pressure<'p>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: impl Into<ReservoirPressure<'p>>,
        input: MaterialDensityMomentumInput<'_>,
        controls: MaterialDensityMomentumControls,
        gate: &CancelGate,
    ) -> Result<MaterialDensityMomentumStep, ChtError> {
        self.force_project_material_momentum_with_optional_pressure(
            domain,
            faces,
            Some(reservoir_pressure.into().admit(faces, gate)?),
            input,
            controls,
            gate,
        )
    }

    /// Apply stationary Darcy resistance implicitly in the native pressure
    /// solve, then scatter its accepted impulse to the original material P.
    /// `resistance_pa_s_m2` is one nonnegative frozen beta per fluid cell and
    /// zero in solids. Face resistance integrates the original half-cell
    /// volumes. Pressure mobility uses M + dt B; inertia, kinetic energy and
    /// impulse shares continue to use the actual M/V and authoritative M/2.
    ///
    /// Drag impulse is -dt B u_final, with equal support for held velocities.
    /// Physical resistance heat remains in `evolution.report`; this momentum
    /// entry does not insert it into material enthalpy. Numerical reconstruction,
    /// backward-Euler and scatter losses are separately audited and never heat.
    ///
    /// # Errors
    /// Invalid resistance or any fresh reconstruction, native solve, impulse,
    /// work, scatter or cancellation refusal leaves all supplied inputs unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn force_project_material_momentum_with_resistance(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: MaterialDensityMomentumInput<'_>,
        resistance_pa_s_m2: &[f64],
        controls: MaterialDensityMomentumControls,
        gate: &CancelGate,
    ) -> Result<ResistiveMaterialDensityMomentumStep, ChtError> {
        self.force_project_material_momentum_with_optional_resistance_pressure(
            domain,
            faces,
            None,
            input,
            resistance_pa_s_m2,
            controls,
            gate,
        )
    }

    /// The same implicit-resistance transfer with authored plane or spatial
    /// reservoir pressures. Component-relative pressure differences reach each
    /// boundary half-cell once. Opening work remains mechanical; it is not
    /// added to resistance heat or represented as a material source.
    ///
    /// # Errors
    /// Invalid pressure or any refusal of
    /// [`Self::force_project_material_momentum_with_resistance`].
    #[allow(clippy::too_many_arguments)]
    pub fn force_project_material_momentum_with_resistance_and_reservoir_pressure<'p>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: impl Into<ReservoirPressure<'p>>,
        input: MaterialDensityMomentumInput<'_>,
        resistance_pa_s_m2: &[f64],
        controls: MaterialDensityMomentumControls,
        gate: &CancelGate,
    ) -> Result<ResistiveMaterialDensityMomentumStep, ChtError> {
        self.force_project_material_momentum_with_optional_resistance_pressure(
            domain,
            faces,
            Some(reservoir_pressure.into().admit(faces, gate)?),
            input,
            resistance_pa_s_m2,
            controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn force_project_material_momentum_with_optional_resistance_pressure(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure: Option<BoundaryPressure<'_>>,
        input: MaterialDensityMomentumInput<'_>,
        resistance_pa_s_m2: &[f64],
        controls: MaterialDensityMomentumControls,
        gate: &CancelGate,
    ) -> Result<ResistiveMaterialDensityMomentumStep, ChtError> {
        poll(gate)?;
        admit_feedback_controls(controls)?;
        let reconstructed = self.reconstruct_material_momentum(
            domain,
            faces,
            input.mass_kg,
            input.momentum_kg_m_s,
            controls.reconstruction,
            gate,
        )?;
        let (evolution, pressure_readback) = reconstructed
            .field
            .force_and_project_density_with_resistance_readback(
                domain,
                faces,
                boundary_pressure,
                DensityForceInput {
                    density_kg_m3: &reconstructed.density_kg_m3,
                    acceleration_m_s2: input.acceleration_m_s2,
                    force_density_n_m3: input.force_density_n_m3,
                },
                resistance_pa_s_m2,
                input.target_outflow_m3_s,
                controls.projection,
                controls.forcing,
                gate,
            )?;
        let (momentum, report) = scatter_resistive_projected_material_momentum(
            self,
            domain,
            faces,
            boundary_pressure,
            &pressure_readback,
            input,
            controls,
            &reconstructed,
            &evolution.projection,
            &evolution.report,
            resistance_pa_s_m2,
            gate,
        )?;
        Ok(ResistiveMaterialDensityMomentumStep {
            momentum_kg_m_s: momentum,
            density_kg_m3: reconstructed.density_kg_m3,
            reconstruction: reconstructed.report,
            evolution,
            report,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn force_project_material_momentum_with_optional_pressure(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure: Option<BoundaryPressure<'_>>,
        input: MaterialDensityMomentumInput<'_>,
        controls: MaterialDensityMomentumControls,
        gate: &CancelGate,
    ) -> Result<MaterialDensityMomentumStep, ChtError> {
        poll(gate)?;
        admit_feedback_controls(controls)?;
        let reconstructed = self.reconstruct_material_momentum(
            domain,
            faces,
            input.mass_kg,
            input.momentum_kg_m_s,
            controls.reconstruction,
            gate,
        )?;
        let (evolution, pressure_readback) = reconstructed
            .field
            .force_and_project_density_with_optional_pressure(
                domain,
                faces,
                boundary_pressure,
                DensityForceInput {
                    density_kg_m3: &reconstructed.density_kg_m3,
                    acceleration_m_s2: input.acceleration_m_s2,
                    force_density_n_m3: input.force_density_n_m3,
                },
                input.target_outflow_m3_s,
                controls.projection,
                controls.forcing,
                gate,
            )?;
        let (momentum, report) = scatter_projected_material_momentum(
            self,
            domain,
            faces,
            boundary_pressure,
            pressure_readback.as_ref(),
            input,
            controls,
            &reconstructed,
            &evolution.projection,
            gate,
        )?;
        poll(gate)?;
        Ok(MaterialDensityMomentumStep {
            momentum_kg_m_s: momentum,
            density_kg_m3: reconstructed.density_kg_m3,
            reconstruction: reconstructed.report,
            evolution,
            report,
        })
    }
}

pub(in crate::conjugate::flow) fn admit_feedback_controls(
    controls: MaterialDensityMomentumControls,
) -> Result<(), ChtError> {
    for (name, maximum) in [
        (
            "flow.feedback_momentum_budget",
            controls.maximum_feedback_momentum_defect_kg_m_s,
        ),
        (
            "flow.feedback_angular_budget",
            controls.maximum_feedback_angular_impulse_defect_kg_m2_s,
        ),
        (
            "flow.feedback_energy_budget",
            controls.maximum_feedback_energy_defect_j,
        ),
    ] {
        nonnegative(name, maximum)?;
    }
    Ok(())
}

/// Consume only an internally computed projection. Borrowing lets the thermal
/// volume owner keep its accepted mass routes without copying their pressure
/// field or applying the pressure impulse a second time.
#[allow(clippy::too_many_arguments)]
pub(in crate::conjugate::flow) fn scatter_projected_material_momentum(
    held: &FlowField,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    boundary_pressure: Option<BoundaryPressure<'_>>,
    pressure_readback: Option<&PressureReadback>,
    input: MaterialDensityMomentumInput<'_>,
    controls: MaterialDensityMomentumControls,
    reconstructed: &DensityMomentumReconstruction,
    projection: &DensityProjection,
    gate: &CancelGate,
) -> Result<(Vec<[f64; 3]>, MaterialDensityMomentumReport), ChtError> {
    let exchange = Exchange {
        held,
        domain,
        faces,
        boundary_pressure,
        pressure_readback,
        input,
        controls,
        reconstructed,
        projection,
        resistance_pa_s_m2: None,
    };
    exchange.admit(gate)?;
    let momentum = exchange.scatter(gate)?;
    let report = exchange.audit(&momentum, gate)?;
    poll(gate)?;
    Ok((momentum, report))
}

/// Return the accepted internally computed resistance impulse to original P.
/// The projection, relative pressure readback and physical resistance account
/// must come from the same fresh native solve. Borrowing lets thermal-volume
/// trials retain only their accepted result without another reconstruction,
/// pressure solve or copy of the accepted material transport field.
#[allow(clippy::too_many_arguments)]
pub(in crate::conjugate::flow) fn scatter_resistive_projected_material_momentum(
    held: &FlowField,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    boundary_pressure: Option<BoundaryPressure<'_>>,
    pressure_readback: &PressureReadback,
    input: MaterialDensityMomentumInput<'_>,
    controls: MaterialDensityMomentumControls,
    reconstructed: &DensityMomentumReconstruction,
    projection: &DensityProjection,
    native: &DensityResistanceReport,
    resistance_pa_s_m2: &[f64],
    gate: &CancelGate,
) -> Result<(Vec<[f64; 3]>, MaterialDensityMomentumReport), ChtError> {
    let exchange = Exchange {
        held,
        domain,
        faces,
        boundary_pressure,
        pressure_readback: Some(pressure_readback),
        input,
        controls,
        reconstructed,
        projection,
        resistance_pa_s_m2: Some(resistance_pa_s_m2),
    };
    exchange.admit(gate)?;
    let momentum = exchange.scatter(gate)?;
    let mut report = exchange.audit(&momentum, gate)?;
    // Independently combine actual cell energy with the fresh physical
    // resistance account. The exchange's local centered identity above
    // treats body plus accepted drag as one algebraic impulse; that is not
    // a claim that physical drag is a midpoint-work body-force predictor.
    let pressure = &projection.report;
    let terms = [
        report.kinetic_energy_change_j,
        report.numerical_scatter_loss_j,
        reconstructed.report.numerical_constraint_loss_j,
        native.backward_euler_loss_j,
        native.resistance_heat_j,
        -reconstructed.report.constraint_reservoir_work_j,
        -native.body_endpoint_work_j,
        -native.constraint_endpoint_work_j,
        -pressure.pressure_source_work_j,
        pressure.fixed_boundary_pressure_work_j,
        pressure.open_boundary_pressure_work_j,
        -report.mass_representation_work_j,
    ];
    let composed = terms
        .into_iter()
        .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value))
        .to_f64();
    finite("flow.resistance_feedback_composed_energy", composed)?;
    let absolute =
        (Dd::from_f64(report.absolute_energy_defect_j) + Dd::from_f64(composed.abs())).to_f64();
    check_budget(
        "flow.resistance_feedback_energy_budget",
        absolute,
        controls.maximum_feedback_energy_defect_j,
    )?;
    report.composed_energy_defect_j = composed;
    report.absolute_energy_defect_j = absolute;
    poll(gate)?;
    Ok((momentum, report))
}

struct Exchange<'a> {
    held: &'a FlowField,
    domain: &'a VoxelDomain,
    faces: [FlowFace; 6],
    boundary_pressure: Option<BoundaryPressure<'a>>,
    pressure_readback: Option<&'a PressureReadback>,
    input: MaterialDensityMomentumInput<'a>,
    controls: MaterialDensityMomentumControls,
    reconstructed: &'a DensityMomentumReconstruction,
    projection: &'a DensityProjection,
    resistance_pa_s_m2: Option<&'a [f64]>,
}

struct FaceImpulse {
    neighbor: Option<usize>,
    mass: f64,
    native_mass: f64,
    momentum: f64,
    unconstrained: f64,
    constraint: f64,
    force: f64,
    pressure: f64,
    impulse: f64,
    kick: f64,
    area: f64,
    q0: f64,
    q1: f64,
}

impl Exchange<'_> {
    fn admit(&self, gate: &CancelGate) -> Result<(), ChtError> {
        let n = self.domain.cell_count();
        if self.input.mass_kg.len() != n
            || self.input.momentum_kg_m_s.len() != n
            || self.input.force_density_n_m3.len() != n
            || self.input.target_outflow_m3_s.len() != n
            || self.reconstructed.density_kg_m3.len() != n
            || self.projection.pressure_correction_pa.len() != n
            || self.resistance_pa_s_m2.is_some_and(|beta| beta.len() != n)
        {
            return Err(refuse(
                "flow.feedback_shape",
                "all stages need the same original cells",
            ));
        }
        for c in 0..n {
            poll_every(gate, c)?;
            let expected = if self.domain.is_fluid(c) {
                self.input.mass_kg[c] / self.domain.volume(c)
            } else {
                0.0
            };
            if self.reconstructed.density_kg_m3[c] != expected {
                return Err(refuse(
                    "flow.feedback_density",
                    "native stages must retain the supplied actual M/V",
                ));
            }
            finite(
                "flow.feedback_pressure",
                self.projection.pressure_correction_pa[c],
            )?;
        }
        for field in [self.held, &self.reconstructed.field, &self.projection.field] {
            field.admit_density_field(
                self.domain,
                self.faces,
                &self.reconstructed.density_kg_m3,
                Some(self.input.target_outflow_m3_s),
                gate,
            )?;
        }
        Ok(())
    }

    fn face(&self, c: usize, f: usize) -> Result<FaceImpulse, ChtError> {
        let axis = f / 2;
        let neighbor = self
            .domain
            .neighbor(c, f)
            .filter(|&n| self.domain.is_fluid(n));
        let held = self
            .domain
            .neighbor(c, f)
            .map_or(self.faces[f] != FlowFace::Free, |n| {
                !self.domain.is_fluid(n)
            });
        let mass =
            0.5 * self.input.mass_kg[c] + neighbor.map_or(0.0, |n| 0.5 * self.input.mass_kg[n]);
        let native_mass = face_dual_mass(
            self.domain,
            &self.reconstructed.density_kg_m3,
            c,
            f,
            neighbor,
        );
        finite_positive("flow.feedback_half_mass", mass)?;
        finite_positive("flow.feedback_native_mass", native_mass)?;
        let momentum = 0.5 * self.input.momentum_kg_m_s[c][axis]
            + neighbor.map_or(0.0, |n| 0.5 * self.input.momentum_kg_m_s[n][axis]);
        let unconstrained = momentum / native_mass;
        let area = self.domain.face_area(c, axis);
        let (slot, sign) = self.held.slot(self.domain, c, f);
        let q0 = self.reconstructed.field.face(axis, slot);
        let q1 = self.projection.field.face(axis, slot);
        let expected_q0 = if held {
            self.held.face(axis, slot)
        } else {
            area * unconstrained
        };
        if q0 != expected_q0 || (held && q1 != q0) {
            return Err(refuse(
                "flow.feedback_starting_field",
                "native fields disagree with original half momenta or held fluxes",
            ));
        }
        let constraint = if held {
            native_mass * (q0 / area - unconstrained)
        } else {
            0.0
        };
        let dt = self.controls.projection.duration_s;
        let (force, pressure) = if held {
            // Native body force and drag are balanced by the held reaction.
            (0.0, 0.0)
        } else {
            let force_per_area =
                0.5 * self.domain.widths(c)[axis] * self.input.force_density_n_m3[c][axis]
                    + neighbor.map_or(0.0, |n| {
                        0.5 * self.domain.widths(n)[axis] * self.input.force_density_n_m3[n][axis]
                    });
            let body =
                dt * (native_mass * self.input.acceleration_m_s2[axis] + area * force_per_area);
            // This is the accepted native drag, not a second drag predictor.
            // Its exact impulse joins the same authoritative half-mass shares.
            let force = self.resistance_pa_s_m2.map_or(body, |beta| {
                body - dt * face_resistance_kg_s(self.domain, beta, c, f, neighbor) * (q1 / area)
            });
            let difference = self.pressure_readback.map_or_else(
                || {
                    self.projection.pressure_correction_pa[c]
                        - neighbor.map_or(0.0, |n| self.projection.pressure_correction_pa[n])
                },
                |readback| {
                    readback.difference(
                        c,
                        neighbor,
                        self.boundary_pressure.map_or(0.0, |p| p.value(f, slot)),
                    )
                },
            );
            let pressure = sign * dt * area * difference;
            (force, pressure)
        };
        let impulse = constraint + (force + pressure);
        let kick = impulse / mass;
        for value in [
            momentum,
            unconstrained,
            constraint,
            force,
            pressure,
            impulse,
            kick,
        ] {
            finite("flow.feedback_face_impulse", value)?;
        }
        Ok(FaceImpulse {
            neighbor,
            mass,
            native_mass,
            momentum,
            unconstrained,
            constraint,
            force,
            pressure,
            impulse,
            kick,
            area,
            q0,
            q1,
        })
    }

    fn scatter(&self, gate: &CancelGate) -> Result<Vec<[f64; 3]>, ChtError> {
        let mut momentum = self.input.momentum_kg_m_s.to_vec();
        for (c, candidate) in momentum.iter_mut().enumerate() {
            poll_every(gate, c)?;
            if !self.domain.is_fluid(c) {
                continue;
            }
            let half_mass = 0.5 * self.input.mass_kg[c];
            for (axis, value) in candidate.iter_mut().enumerate() {
                let left = self.face(c, 2 * axis)?;
                let right = self.face(c, 2 * axis + 1)?;
                let impulse = half_mass * left.kick + half_mass * right.kick;
                // Preserve original sub-face momentum, including signed zero.
                if impulse != 0.0 {
                    *value += impulse;
                }
                finite("flow.feedback_cell_momentum", *value)?;
            }
        }
        Ok(momentum)
    }

    #[allow(clippy::too_many_lines)]
    fn audit(
        &self,
        momentum: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<MaterialDensityMomentumReport, ChtError> {
        if momentum.len() != self.domain.cell_count() {
            return Err(refuse(
                "flow.feedback_result_shape",
                "one returned momentum per original cell is required",
            ));
        }
        let mut face_p: [Sum; 3] = std::array::from_fn(|_| Sum::default());
        let mut cell_p: [Sum; 3] = std::array::from_fn(|_| Sum::default());
        let mut face_l: [Sum; 3] = std::array::from_fn(|_| Sum::default());
        let mut cell_l: [Sum; 3] = std::array::from_fn(|_| Sum::default());
        let mut impulse_defect = Sum::default();
        let mut angular_defect = Sum::default();
        let mut before = Sum::default();
        let mut after = Sum::default();
        let mut change = Sum::default();
        let mut native_work = Sum::default();
        let mut scatter_loss = Sum::default();
        let mut mass_work = Sum::default();
        let mut feedback_defect = Sum::default();
        let mut composition_defect = Sum::default();
        let mut energy_defect = Sum::default();
        for c in 0..self.domain.cell_count() {
            poll_every(gate, c)?;
            if !self.domain.is_fluid(c) {
                if momentum[c] != [0.0; 3] {
                    return Err(refuse(
                        "flow.feedback_solid_momentum",
                        "solid cells cannot acquire fluid momentum",
                    ));
                }
                continue;
            }
            let mass = self.input.mass_kg[c];
            let half = 0.5 * mass;
            let [x, y, z] = self.domain.coords(c);
            let center = self.domain.center(x, y, z);
            let mut delta = [0.0; 3];
            let mut local_p_error = [0.0; 3];
            for axis in 0..3 {
                let left = self.face(c, 2 * axis)?;
                let right = self.face(c, 2 * axis + 1)?;
                let expected = half * left.kick + half * right.kick;
                let p0 = self.input.momentum_kg_m_s[c][axis];
                let p1 = momentum[c][axis];
                finite("flow.feedback_returned_momentum", p1)?;
                let actual = p1 - p0;
                delta[axis] = actual;
                local_p_error[axis] = actual - expected;
                impulse_defect.add(local_p_error[axis].abs());
                cell_p[axis].add(actual);
                let u0 = p0 / mass;
                let du = actual / mass;
                before.add(0.5 * p0 * u0);
                after.add(0.5 * p1 * (p1 / mass));
                change.add((u0 + 0.5 * du) * actual);
                let relative = right.kick - left.kick;
                let loss = 0.125 * mass * relative * relative;
                scatter_loss.add(loss);
                // Only kick energies are subtracted. A common original velocity
                // appears solely multiplied by the actual momentum defect.
                let local = u0.mul_add(
                    local_p_error[axis],
                    0.5 * actual * du
                        - 0.5 * half * (left.kick * left.kick + right.kick * right.kick)
                        + loss,
                );
                feedback_defect.add(local);
                energy_defect.add(local.abs());
            }
            add_vector(&mut cell_l, cross(center, delta));
            for component in cross(center, local_p_error) {
                angular_defect.add(component.abs());
            }
            let volume_defect = self.projection.field.net_outflow(self.domain, c)
                - self.input.target_outflow_m3_s[c];
            check_budget(
                "flow.feedback_volume_budget",
                volume_defect.abs(),
                self.controls.projection.maximum_volume_defect_m3_s,
            )?;
            let pressure_volume_defect = self.controls.projection.duration_s
                * self
                    .pressure_readback
                    .map_or(self.projection.pressure_correction_pa[c], |readback| {
                        readback.cell_pressure(c)
                    })
                * volume_defect;
            composition_defect.add(pressure_volume_defect);
            energy_defect.add(pressure_volume_defect.abs());
            for f in 0..6 {
                let face = self.face(c, f)?;
                if face.neighbor.is_some() && f % 2 == 0 {
                    continue;
                }
                let axis = f / 2;
                face_p[axis].add(face.impulse);
                let mut position = center;
                position[axis] +=
                    (if f % 2 == 0 { -0.5 } else { 0.5 }) * self.domain.widths(c)[axis];
                let mut normal_impulse = [0.0; 3];
                normal_impulse[axis] = face.impulse;
                add_vector(&mut face_l, cross(position, normal_impulse));
                let stage_actual = face.native_mass * ((face.q1 - face.q0) / face.area);
                let stage_error = stage_actual - (face.force + face.pressure);
                let final_du = face.q1 / face.area - face.unconstrained;
                let final_error = face.native_mass * final_du - face.impulse;
                let other_half = face.neighbor.map_or(0.0, |n| 0.5 * self.input.mass_kg[n]);
                let distributed = half * face.kick + other_half * face.kick;
                let share_error = distributed - face.impulse;
                for defect in [stage_error, final_error, share_error] {
                    impulse_defect.add(defect.abs());
                    let mut error = [0.0; 3];
                    error[axis] = defect;
                    for component in cross(position, error) {
                        angular_defect.add(component.abs());
                    }
                }
                let native_du = face.impulse / face.native_mass;
                native_work.add((face.unconstrained + 0.5 * native_du) * face.impulse);
                let representation = (face.native_mass - face.mass)
                    * (face.momentum / face.mass + 0.5 * face.kick)
                    * native_du;
                mass_work.add(representation);
                energy_defect.add(representation.abs());
                // Compare authoritative P/2 with actual P/M at both donors;
                // neither this representation nor final Q/A rounding is force.
                let represented_p = half * (self.input.momentum_kg_m_s[c][axis] / mass)
                    + face.neighbor.map_or(0.0, |n| {
                        other_half * (self.input.momentum_kg_m_s[n][axis] / self.input.mass_kg[n])
                    });
                let local =
                    (represented_p - face.momentum) * face.kick + 0.5 * face.kick * share_error;
                feedback_defect.add(local);
                energy_defect.add(local.abs());
                energy_defect
                    .add(((face.unconstrained + 0.5 * (final_du + native_du)) * final_error).abs());
                let q_predicted = face.q0 + face.area * (face.force / face.native_mass);
                let pressure_du = (face.q1 - q_predicted) / face.area;
                let constraint_du = if face.constraint == 0.0 {
                    0.0
                } else {
                    face.q0 / face.area - face.unconstrained
                };
                let endpoint_relative = face.q0 / face.area - face.unconstrained
                    + 0.5 * ((face.q1 - face.q0) / face.area);
                // Center all actual work at the original unconstrained face
                // velocity. Applied and held-body forces cancel on held faces.
                let local_composed = 0.5 * face.impulse * native_du
                    + 0.5 * face.native_mass * constraint_du * constraint_du
                    + 0.5 * face.native_mass * pressure_du * pressure_du
                    + 0.5 * face.force * pressure_du
                    - face.constraint * constraint_du
                    - face.force * endpoint_relative
                    - face.pressure * final_du;
                composition_defect.add(local_composed);
                energy_defect.add(local_composed.abs());
            }
        }
        let face_p = face_p.map(|s| s.value);
        let cell_p = cell_p.map(|s| s.value);
        let face_l = face_l.map(|s| s.value);
        let cell_l = cell_l.map(|s| s.value);
        for axis in 0..3 {
            for value in [face_p[axis], cell_p[axis], face_l[axis], cell_l[axis]] {
                finite("flow.feedback_inventory", value)?;
            }
            impulse_defect.add((cell_p[axis] - face_p[axis]).abs());
            angular_defect.add((cell_l[axis] - face_l[axis]).abs());
        }
        let composed = feedback_defect.value + composition_defect.value;
        energy_defect.add(feedback_defect.value.abs());
        energy_defect.add(composed.abs());
        for value in [
            before.value,
            after.value,
            change.value,
            native_work.value,
            scatter_loss.value,
            mass_work.value,
            feedback_defect.value,
            composed,
            energy_defect.value,
        ] {
            finite("flow.feedback_energy", value)?;
        }
        check_budget(
            "flow.feedback_momentum_budget",
            impulse_defect.value,
            self.controls.maximum_feedback_momentum_defect_kg_m_s,
        )?;
        check_budget(
            "flow.feedback_angular_budget",
            angular_defect.value,
            self.controls
                .maximum_feedback_angular_impulse_defect_kg_m2_s,
        )?;
        check_budget(
            "flow.feedback_energy_budget",
            energy_defect.value,
            self.controls.maximum_feedback_energy_defect_j,
        )?;
        Ok(MaterialDensityMomentumReport {
            face_impulse_kg_m_s: face_p,
            cell_momentum_change_kg_m_s: cell_p,
            absolute_momentum_defect_kg_m_s: impulse_defect.value,
            face_angular_impulse_kg_m2_s: face_l,
            cell_angular_momentum_change_kg_m2_s: cell_l,
            absolute_angular_impulse_defect_kg_m2_s: angular_defect.value,
            kinetic_energy_before_j: before.value,
            kinetic_energy_after_j: after.value,
            kinetic_energy_change_j: change.value,
            native_impulse_work_j: native_work.value,
            numerical_scatter_loss_j: scatter_loss.value,
            mass_representation_work_j: mass_work.value,
            feedback_energy_defect_j: feedback_defect.value,
            composed_energy_defect_j: composed,
            absolute_energy_defect_j: energy_defect.value,
        })
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn add_vector(sum: &mut [Sum; 3], values: [f64; 3]) {
    for (sum, value) in sum.iter_mut().zip(values) {
        sum.add(value);
    }
}

fn check_budget(field: &'static str, defect: f64, maximum: f64) -> Result<(), ChtError> {
    finite(field, defect)?;
    if defect > maximum {
        return Err(refuse(
            field,
            &format!("actual absolute defect {defect:e} exceeds {maximum:e}"),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
