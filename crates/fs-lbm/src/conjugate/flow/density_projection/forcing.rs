//! Physical body-force impulse on the same frozen face masses as pressure.

use super::{
    BoundaryPressure, CancelGate, ChtError, DensityProjection, DensityProjectionControls, FlowFace,
    FlowField, PressureReadback, ReservoirPressure, StaggeredBoundaryPressure, Sum, VoxelDomain,
    face_dual_mass, finite, finite_positive, nonnegative, poll, poll_every, refuse,
};

/// Frozen material density and prescribed physical body force for one step.
#[derive(Debug, Clone, Copy)]
pub struct DensityForceInput<'a> {
    /// Positive fluid density, kg/m^3; entries in solid cells are ignored.
    pub density_kg_m3: &'a [f64],
    /// Uniform acceleration in world axes, m/s^2, including gravity if present.
    pub acceleration_m_s2: [f64; 3],
    /// Additional cell force per physical volume in world axes, N/m^3.
    /// Supply one entry per cell, with exactly zero entries in solids. This
    /// excludes pressure, which the composed projection applies once itself.
    pub force_density_n_m3: &'a [[f64; 3]],
}

/// Allocation and independent physical readback budgets for body forcing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DensityForceControls {
    /// Maximum domain cell count, checked before cloning a face field.
    pub max_cells: usize,
    /// Absolute sum of local impulse defects plus the global vector L1 defect,
    /// kg m/s. Applies independently to the predictor and the composed result.
    pub maximum_impulse_defect_kg_m_s: f64,
    /// Absolute sum of local force-work defects plus the global work defect,
    /// J. Also bounds the composed force/pressure energy identity.
    pub maximum_energy_defect_j: f64,
}

/// Independent readback of a body-force predictor and its held-face reactions.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityForceReport {
    /// `dt sum(V (rho acceleration + force_density))` in world axes, kg m/s.
    /// Includes the fluid half-dual volumes adjacent to walls and solids.
    pub applied_impulse_kg_m_s: [f64; 3],
    /// Reaction on the fluid that cancels forcing on held `Fixed`, `Wall`,
    /// and fluid-solid face degrees of freedom, in world axes, kg m/s.
    pub constraint_impulse_kg_m_s: [f64; 3],
    /// Actual change of staggered `sum(M_face u_face)` in world axes, kg m/s.
    pub momentum_change_kg_m_s: [f64; 3],
    /// Sum of local absolute defects and the absolute global vector defect.
    pub absolute_impulse_defect_kg_m_s: f64,
    /// Face-dual kinetic energy before forcing, J.
    pub kinetic_energy_before_j: f64,
    /// Face-dual kinetic energy after forcing, J.
    pub kinetic_energy_after_j: f64,
    /// Stable `sum(M (u_before + delta_u/2) delta_u)`, J.
    pub kinetic_energy_change_j: f64,
    /// Applied impulse dotted with the face midpoint velocity, J. This is
    /// discrete force-stage work, not a gravitational-potential-energy update.
    pub applied_work_j: f64,
    /// Held-face reaction work, J. It cancels the applied work on a nonzero
    /// `Fixed` flow; stationary wall and fluid-solid faces do zero work.
    pub constraint_work_j: f64,
    /// `delta_K - applied_work - constraint_work`, J.
    pub energy_defect_j: f64,
    /// Sum of local absolute work defects and the absolute global defect, J.
    pub absolute_energy_defect_j: f64,
}

/// Unprojected physical body-force predictor, suitable for a material owner
/// that will derive and apply its own compatible volume-source projection.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityForcePrediction {
    /// Predicted native volumetric face fluxes, m^3/s.
    pub field: FlowField,
    /// Accepted force and constraint work evidence.
    pub report: DensityForceReport,
}

/// Independent readback spanning the original field and both accepted stages.
#[derive(Debug, Clone, PartialEq)]
pub struct ForcedDensityProjectionReport {
    /// Pressure impulse on the pressure-corrected face degrees of freedom,
    /// reconstructed from returned cell pressure and physical areas, kg m/s.
    pub pressure_impulse_kg_m_s: [f64; 3],
    /// Pressure impulse plus the predictor's held-face constraint impulse,
    /// kg m/s. In a sealed hydrostatic column this balances its total weight.
    /// This is a net discrete support account, not resolved wall tractions.
    pub pressure_and_constraint_impulse_kg_m_s: [f64; 3],
    /// Actual original-to-final staggered momentum change, kg m/s.
    pub momentum_change_kg_m_s: [f64; 3],
    /// Predictor defects, local pressure-impulse defects, and the absolute
    /// original-to-final global vector defect, kg m/s.
    pub absolute_impulse_defect_kg_m_s: f64,
    /// Stable original-to-final face-dual kinetic-energy increment, J.
    pub kinetic_energy_change_j: f64,
    /// Applied body impulse dotted with original-to-final midpoint velocity,
    /// J. This differs from the intermediate predictor-stage force work.
    pub applied_endpoint_work_j: f64,
    /// Held-face constraint impulse dotted with original-to-final midpoint
    /// velocity, J; cancels applied endpoint work on held nonzero fluxes.
    pub constraint_endpoint_work_j: f64,
    /// Signed splitting cross term `sum((I_force + I_constraint)
    /// (u_final - u_predicted) / 2)`, J. Endpoint force/constraint work equals
    /// predictor work plus this term. At exact hydrostatic rest this term
    /// cancels projection loss, while endpoint gravity work is zero.
    pub force_pressure_split_work_j: f64,
    /// `delta_K + projection_loss + split_work - applied_endpoint_work
    /// - constraint_endpoint_work - pressure_source_work
    /// + fixed_boundary_pressure_work + open_boundary_pressure_work`, J.
    pub energy_defect_j: f64,
}

/// One body-force kick followed by one native density/volume projection.
#[derive(Debug, Clone, PartialEq)]
pub struct ForcedDensityProjection {
    /// Corrected fluxes, pressure in Pa, and native continuity/work evidence.
    pub projection: DensityProjection,
    /// Independent force-stage readback, including held-face reactions.
    pub forcing: DensityForceReport,
    /// Independent original-to-final pressure and force readback.
    pub report: ForcedDensityProjectionReport,
}

impl FlowField {
    /// Apply prescribed physical body force with frozen native face inertia.
    /// Each internal face uses the two actual half-cell masses and forces;
    /// open faces use the adjacent half-cell. In particular gravity is
    /// `rho g`, with no reference-density or Boussinesq approximation.
    ///
    /// Every fluid half-cell is counted once per axis, including those beside
    /// walls or solids. `Wall`, fluid-solid and `Fixed` faces stay unchanged;
    /// their cancelled impulse and reaction work are reported explicitly.
    /// The independent audit reads back the returned fluxes and compares them
    /// with the cell-integrated applied force, rather than trusting the kick.
    ///
    /// This is a force predictor, not a continuity solve, advective momentum
    /// update, moving-density integrator, or cell-to-face momentum coupling.
    /// Call a compatible density projection once afterward, or pass its field
    /// to the material/phase owner that constructs that projection. Density
    /// and geometry must match that owner's frozen input. Pressure force must
    /// not also be supplied through `force_density_n_m3`.
    ///
    /// # Errors
    /// Invalid shapes/data, nonpositive fluid density or duration, blocked
    /// fluxes, force in a solid cell, exceeded cell/impulse/work budget, or
    /// cancellation. The original field is unchanged on every refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn predict_density_body_force(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: DensityForceInput<'_>,
        duration_s: f64,
        controls: DensityForceControls,
        gate: &CancelGate,
    ) -> Result<DensityForcePrediction, ChtError> {
        admit(self, domain, faces, input, duration_s, controls, gate)?;
        let mut field = self.clone();
        for c in 0..domain.cell_count() {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                continue;
            }
            for f in 0..6 {
                let nb = fluid_neighbor(domain, c, f);
                if !owns_face(f, nb) || !corrected_face(domain, faces, c, f) {
                    continue;
                }
                let axis = f / 2;
                let area = domain.face_area(c, axis);
                let mass = face_dual_mass(domain, input.density_kg_m3, c, f, nb);
                finite_positive("flow.force_face_mass_kg", mass)?;
                let impulse = face_impulse(domain, input, duration_s, c, f, nb, mass);
                finite("flow.body_force_impulse", impulse)?;
                let (slot, _) = self.slot(domain, c, f);
                let q = self.face(axis, slot) + area * (impulse / mass);
                finite("flow.forced_face_flux", q)?;
                *field.face_mut(axis, slot) = q;
            }
        }
        let report = audit_force(
            self, &field, domain, faces, input, duration_s, controls, gate,
        )?;
        poll(gate)?;
        Ok(DensityForcePrediction { field, report })
    }

    /// Apply one density-aware force predictor followed by one native pressure
    /// correction onto `target_outflow_m3_s`. Both stages use the identical
    /// frozen face dual mass and duration. This supports exact discrete
    /// hydrostatic balance on graded, sharply stratified columns and genuine
    /// density-driven circulation where pressure and density gradients differ.
    ///
    /// This bounded split stage does not advect staggered inertia, reconstruct
    /// faces from conservative cell momentum, evolve density, or update
    /// gravitational potential energy. The returned projection loss remains
    /// numerical loss and never heats material. See the predictor's contract
    /// for force ownership and the native projection's source compatibility.
    ///
    /// # Errors
    /// Any predictor or projection refusal, an excessive independently read
    /// back combined impulse/work defect, or cancellation; no input mutation.
    #[allow(clippy::too_many_arguments)]
    pub fn force_and_project_density(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: DensityForceInput<'_>,
        target_outflow_m3_s: &[f64],
        projection_controls: DensityProjectionControls,
        force_controls: DensityForceControls,
        gate: &CancelGate,
    ) -> Result<ForcedDensityProjection, ChtError> {
        self.force_and_project_density_with_optional_pressure(
            domain,
            faces,
            None,
            input,
            target_outflow_m3_s,
            projection_controls,
            force_controls,
            gate,
        )
        .map(|(evolution, _)| evolution)
    }

    /// Apply body forcing and a spatial reservoir-pressure profile through the
    /// same density-weighted projection. The profile supplies pressure on actual
    /// open staggered faces; their impulses use the adjacent half-cell inertia.
    /// A stratified hydrostatic head can balance gravity at an extended opening,
    /// and unequal reservoir pressures can drive inflow or outflow.
    ///
    /// Pressure impulses and reservoir work are independently read from the
    /// returned native field. Reservoir work retains the projection's component
    /// gauge and is included in the combined mechanical energy balance. The
    /// borrowed profile is neither copied nor retained. Reservoir material,
    /// thermodynamic state and any subsequent mass transport remain caller-owned.
    ///
    /// # Errors
    /// Any original force/projection refusal, invalid native pressure profile,
    /// exceeded impulse/work budget or cancellation leaves the input unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn force_and_project_density_with_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        input: DensityForceInput<'_>,
        target_outflow_m3_s: &[f64],
        projection_controls: DensityProjectionControls,
        force_controls: DensityForceControls,
        gate: &CancelGate,
    ) -> Result<ForcedDensityProjection, ChtError> {
        self.force_and_project_density_with_reservoir_pressure(
            domain,
            faces,
            profile,
            input,
            target_outflow_m3_s,
            projection_controls,
            force_controls,
            gate,
        )
    }

    /// Apply the original force/pressure stage with either six optional plane
    /// pressures in Pa or a borrowed [`StaggeredBoundaryPressure`] profile.
    /// Uniform reservoirs need no staggered pressure buffers. Both inputs use
    /// the same native relative-pressure readback and per-face opening work;
    /// large common pressure references do not degrade the impulse calculation.
    /// Pressure must not also be included in the authored body force.
    ///
    /// # Errors
    /// Nonfinite pressure, pressure on a held plane, invalid profile, or any
    /// refusal of [`Self::force_and_project_density_with_pressure_profile`].
    #[allow(clippy::too_many_arguments)]
    pub fn force_and_project_density_with_reservoir_pressure<'p>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: impl Into<ReservoirPressure<'p>>,
        input: DensityForceInput<'_>,
        target_outflow_m3_s: &[f64],
        projection_controls: DensityProjectionControls,
        force_controls: DensityForceControls,
        gate: &CancelGate,
    ) -> Result<ForcedDensityProjection, ChtError> {
        self.force_and_project_density_with_optional_pressure(
            domain,
            faces,
            Some(reservoir_pressure.into().admit(faces, gate)?),
            input,
            target_outflow_m3_s,
            projection_controls,
            force_controls,
            gate,
        )
        .map(|(evolution, _)| evolution)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn force_and_project_density_with_optional_pressure(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure: Option<BoundaryPressure<'_>>,
        input: DensityForceInput<'_>,
        target_outflow_m3_s: &[f64],
        projection_controls: DensityProjectionControls,
        force_controls: DensityForceControls,
        gate: &CancelGate,
    ) -> Result<(ForcedDensityProjection, Option<PressureReadback>), ChtError> {
        let predicted = self.predict_density_body_force(
            domain,
            faces,
            input,
            projection_controls.duration_s,
            force_controls,
            gate,
        )?;
        let (projection, pressure_readback) = match boundary_pressure {
            None => (
                predicted.field.project_density_volume_source(
                    domain,
                    faces,
                    input.density_kg_m3,
                    target_outflow_m3_s,
                    projection_controls,
                    gate,
                )?,
                None,
            ),
            Some(boundary_pressure) => {
                let (projection, readback) = predicted
                    .field
                    .project_density_volume_source_with_pressure_readback(
                        domain,
                        faces,
                        boundary_pressure,
                        input.density_kg_m3,
                        target_outflow_m3_s,
                        projection_controls,
                        gate,
                    )?;
                (projection, Some(readback))
            }
        };
        let report = audit_combined(
            self,
            &predicted,
            &projection,
            domain,
            faces,
            boundary_pressure,
            pressure_readback.as_ref(),
            input,
            projection_controls.duration_s,
            force_controls,
            gate,
        )?;
        poll(gate)?;
        Ok((
            ForcedDensityProjection {
                projection,
                forcing: predicted.report,
                report,
            },
            pressure_readback,
        ))
    }
}

#[allow(clippy::too_many_arguments)]
fn admit(
    field: &FlowField,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    input: DensityForceInput<'_>,
    duration_s: f64,
    controls: DensityForceControls,
    gate: &CancelGate,
) -> Result<(), ChtError> {
    poll(gate)?;
    finite_positive("flow.force_duration_s", duration_s)?;
    nonnegative(
        "flow.maximum_impulse_defect_kg_m_s",
        controls.maximum_impulse_defect_kg_m_s,
    )?;
    nonnegative(
        "flow.maximum_force_energy_defect_j",
        controls.maximum_energy_defect_j,
    )?;
    if domain.cell_count() > controls.max_cells {
        return Err(refuse(
            "flow.force_cell_budget",
            "domain exceeds the explicit force cell budget",
        ));
    }
    field.admit_density_field(domain, faces, input.density_kg_m3, None, gate)?;
    if input.force_density_n_m3.len() != domain.cell_count() {
        return Err(refuse(
            "flow.body_force_shape",
            "force density must have one vector per domain cell",
        ));
    }
    for acceleration in input.acceleration_m_s2 {
        finite("flow.acceleration_m_s2", acceleration)?;
    }
    for (c, force) in input.force_density_n_m3.iter().enumerate() {
        poll_every(gate, c)?;
        for &component in force {
            finite("flow.force_density_n_m3", component)?;
            if !domain.is_fluid(c) && component != 0.0 {
                return Err(refuse(
                    "flow.solid_body_force",
                    "a solid cell cannot carry fluid body force",
                ));
            }
        }
    }
    Ok(())
}

fn fluid_neighbor(domain: &VoxelDomain, c: usize, f: usize) -> Option<usize> {
    domain.neighbor(c, f).filter(|&n| domain.is_fluid(n))
}

fn owns_face(f: usize, fluid_neighbor: Option<usize>) -> bool {
    fluid_neighbor.is_none() || f % 2 == 1
}

fn corrected_face(domain: &VoxelDomain, faces: [FlowFace; 6], c: usize, f: usize) -> bool {
    domain
        .neighbor(c, f)
        .map_or(faces[f] == FlowFace::Free, |n| domain.is_fluid(n))
}

fn face_impulse(
    domain: &VoxelDomain,
    input: DensityForceInput<'_>,
    duration_s: f64,
    c: usize,
    f: usize,
    nb: Option<usize>,
    mass: f64,
) -> f64 {
    let axis = f / 2;
    let force_per_area = 0.5 * domain.widths(c)[axis] * input.force_density_n_m3[c][axis]
        + nb.map_or(0.0, |n| {
            0.5 * domain.widths(n)[axis] * input.force_density_n_m3[n][axis]
        });
    duration_s * (mass * input.acceleration_m_s2[axis] + domain.face_area(c, axis) * force_per_area)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn audit_force(
    before: &FlowField,
    after: &FlowField,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    input: DensityForceInput<'_>,
    duration_s: f64,
    controls: DensityForceControls,
    gate: &CancelGate,
) -> Result<DensityForceReport, ChtError> {
    let mut applied: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut constraint: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut momentum: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut impulse_defect = Sum::default();
    let mut initial = Sum::default();
    let mut final_energy = Sum::default();
    let mut change = Sum::default();
    let mut work = Sum::default();
    let mut constraint_work = Sum::default();
    let mut work_defect = Sum::default();
    for c in 0..domain.cell_count() {
        poll_every(gate, c)?;
        if !domain.is_fluid(c) {
            continue;
        }
        // This cell-integrated input total is independent of face ownership
        // and catches missing/duplicated half-duals, including solid contacts.
        for (axis, total) in applied.iter_mut().enumerate() {
            total.add(
                duration_s
                    * domain.volume(c)
                    * (input.density_kg_m3[c] * input.acceleration_m_s2[axis]
                        + input.force_density_n_m3[c][axis]),
            );
        }
        for f in 0..6 {
            let nb = fluid_neighbor(domain, c, f);
            if !owns_face(f, nb) {
                continue;
            }
            let axis = f / 2;
            let mass = face_dual_mass(domain, input.density_kg_m3, c, f, nb);
            finite_positive("flow.force_face_mass_kg", mass)?;
            let impulse = face_impulse(domain, input, duration_s, c, f, nb, mass);
            let area = domain.face_area(c, axis);
            let (slot, _) = before.slot(domain, c, f);
            let q0 = before.face(axis, slot);
            let q1 = after.face(axis, slot);
            let u0 = q0 / area;
            let u1 = q1 / area;
            let du = (q1 - q0) / area;
            let midpoint = u0 + 0.5 * du;
            let reaction = if corrected_face(domain, faces, c, f) {
                0.0
            } else {
                -impulse
            };
            let delta_p = mass * du;
            let delta_k = mass * midpoint * du;
            let body_work = impulse * midpoint;
            let held_work = reaction * midpoint;
            constraint[axis].add(reaction);
            momentum[axis].add(delta_p);
            impulse_defect.add((delta_p - impulse - reaction).abs());
            initial.add(0.5 * mass * u0 * u0);
            final_energy.add(0.5 * mass * u1 * u1);
            change.add(delta_k);
            work.add(body_work);
            constraint_work.add(held_work);
            work_defect.add((delta_k - body_work - held_work).abs());
        }
    }
    let applied = applied.map(|s| s.value);
    let constraint = constraint.map(|s| s.value);
    let momentum = momentum.map(|s| s.value);
    for axis in 0..3 {
        finite("flow.applied_impulse", applied[axis])?;
        finite("flow.constraint_impulse", constraint[axis])?;
        finite("flow.momentum_change", momentum[axis])?;
        impulse_defect.add((momentum[axis] - applied[axis] - constraint[axis]).abs());
    }
    let energy_defect = change.value - work.value - constraint_work.value;
    work_defect.add(energy_defect.abs());
    for value in [
        initial.value,
        final_energy.value,
        change.value,
        work.value,
        constraint_work.value,
        energy_defect,
        work_defect.value,
    ] {
        finite("flow.force_energy", value)?;
    }
    check_budget(
        "flow.force_impulse_budget",
        impulse_defect.value,
        controls.maximum_impulse_defect_kg_m_s,
    )?;
    check_budget(
        "flow.force_energy_budget",
        work_defect.value,
        controls.maximum_energy_defect_j,
    )?;
    Ok(DensityForceReport {
        applied_impulse_kg_m_s: applied,
        constraint_impulse_kg_m_s: constraint,
        momentum_change_kg_m_s: momentum,
        absolute_impulse_defect_kg_m_s: impulse_defect.value,
        kinetic_energy_before_j: initial.value,
        kinetic_energy_after_j: final_energy.value,
        kinetic_energy_change_j: change.value,
        applied_work_j: work.value,
        constraint_work_j: constraint_work.value,
        energy_defect_j: energy_defect,
        absolute_energy_defect_j: work_defect.value,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn audit_combined(
    before: &FlowField,
    predicted: &DensityForcePrediction,
    projected: &DensityProjection,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    boundary_pressure: Option<BoundaryPressure<'_>>,
    pressure_readback: Option<&PressureReadback>,
    input: DensityForceInput<'_>,
    duration_s: f64,
    controls: DensityForceControls,
    gate: &CancelGate,
) -> Result<ForcedDensityProjectionReport, ChtError> {
    let mut pressure: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut momentum: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut impulse_defect = Sum::default();
    impulse_defect.add(predicted.report.absolute_impulse_defect_kg_m_s);
    let mut change = Sum::default();
    let mut endpoint_work = Sum::default();
    let mut constraint_work = Sum::default();
    let mut split_work = Sum::default();
    for c in 0..domain.cell_count() {
        poll_every(gate, c)?;
        if !domain.is_fluid(c) {
            continue;
        }
        for f in 0..6 {
            let nb = fluid_neighbor(domain, c, f);
            if !owns_face(f, nb) {
                continue;
            }
            let axis = f / 2;
            let area = domain.face_area(c, axis);
            let mass = face_dual_mass(domain, input.density_kg_m3, c, f, nb);
            let (slot, sign) = before.slot(domain, c, f);
            let q0 = before.face(axis, slot);
            let q_predicted = predicted.field.face(axis, slot);
            let q1 = projected.field.face(axis, slot);
            let expected = if corrected_face(domain, faces, c, f) {
                let difference = pressure_readback.map_or_else(
                    || {
                        projected.pressure_correction_pa[c]
                            - nb.map_or(0.0, |n| projected.pressure_correction_pa[n])
                    },
                    |readback| {
                        readback.difference(
                            c,
                            nb,
                            boundary_pressure.map_or(0.0, |p| p.value(f, slot)),
                        )
                    },
                );
                sign * duration_s * area * difference
            } else {
                0.0
            };
            pressure[axis].add(expected);
            impulse_defect.add((mass * ((q1 - q_predicted) / area) - expected).abs());
            let du = (q1 - q0) / area;
            momentum[axis].add(mass * du);
            let midpoint = q0 / area + 0.5 * du;
            change.add(mass * midpoint * du);
            let applied = face_impulse(domain, input, duration_s, c, f, nb, mass);
            let reaction = if corrected_face(domain, faces, c, f) {
                0.0
            } else {
                -applied
            };
            endpoint_work.add(applied * midpoint);
            constraint_work.add(reaction * midpoint);
            split_work.add(0.5 * (applied + reaction) * ((q1 - q_predicted) / area));
        }
    }
    let pressure = pressure.map(|s| s.value);
    let momentum = momentum.map(|s| s.value);
    let support = std::array::from_fn(|axis| {
        pressure[axis] + predicted.report.constraint_impulse_kg_m_s[axis]
    });
    for axis in 0..3 {
        finite("flow.combined_pressure_impulse", pressure[axis])?;
        finite("flow.combined_support_impulse", support[axis])?;
        finite("flow.combined_momentum_change", momentum[axis])?;
        impulse_defect.add(
            (momentum[axis] - predicted.report.applied_impulse_kg_m_s[axis] - support[axis]).abs(),
        );
    }
    finite("flow.combined_energy_change", change.value)?;
    finite("flow.endpoint_body_work", endpoint_work.value)?;
    finite("flow.endpoint_constraint_work", constraint_work.value)?;
    finite("flow.force_pressure_split_work", split_work.value)?;
    let mut energy_defect = change.value + projected.report.projection_loss_j + split_work.value
        - endpoint_work.value
        - constraint_work.value
        - projected.report.pressure_source_work_j
        + projected.report.fixed_boundary_pressure_work_j;
    if projected.report.open_boundary_pressure_work_j != 0.0 {
        energy_defect += projected.report.open_boundary_pressure_work_j;
    }
    check_budget(
        "flow.combined_impulse_budget",
        impulse_defect.value,
        controls.maximum_impulse_defect_kg_m_s,
    )?;
    check_budget(
        "flow.combined_energy_budget",
        energy_defect.abs(),
        controls.maximum_energy_defect_j,
    )?;
    Ok(ForcedDensityProjectionReport {
        pressure_impulse_kg_m_s: pressure,
        pressure_and_constraint_impulse_kg_m_s: support,
        momentum_change_kg_m_s: momentum,
        absolute_impulse_defect_kg_m_s: impulse_defect.value,
        kinetic_energy_change_j: change.value,
        applied_endpoint_work_j: endpoint_work.value,
        constraint_endpoint_work_j: constraint_work.value,
        force_pressure_split_work_j: split_work.value,
        energy_defect_j: energy_defect,
    })
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
