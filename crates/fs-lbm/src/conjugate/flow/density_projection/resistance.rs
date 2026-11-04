//! Stationary Darcy resistance in the original density/pressure momentum equation.

use super::forcing::{admit, corrected_face, face_impulse, fluid_neighbor, owns_face};
use super::*;
use fs_math::dd::Dd;

/// Physical impulse, heat and work of one implicit stationary-resistance step.
#[derive(Clone, Debug, PartialEq)]
pub struct DensityResistanceReport {
    pub applied_impulse_kg_m_s: [f64; 3],
    pub pressure_impulse_kg_m_s: [f64; 3],
    /// Impulse exerted on the fluid by the prescribed stationary matrix.
    pub resistance_impulse_kg_m_s: [f64; 3],
    /// Reaction needed to preserve held normal velocities, including drag.
    pub constraint_impulse_kg_m_s: [f64; 3],
    pub momentum_change_kg_m_s: [f64; 3],
    pub absolute_impulse_defect_kg_m_s: f64,
    /// Actual native face masses, including each fluid half-cell once per axis.
    pub kinetic_energy_before_j: f64,
    pub kinetic_energy_after_j: f64,
    pub kinetic_energy_change_j: f64,
    /// Applied body impulse dotted with the accepted endpoint face velocity.
    pub body_endpoint_work_j: f64,
    /// Held-velocity controller impulse dotted with that same endpoint velocity.
    pub constraint_endpoint_work_j: f64,
    /// Nonnegative sum M_face*(u_final-u_initial)^2/2. This is numerical loss.
    pub backward_euler_loss_j: f64,
    /// Physical dt*beta*|u_final|^2 integrated on each original cell's half-duals.
    /// This step reports heat; only an enclosing caloric owner may deposit it.
    pub cell_resistance_heat_j: Vec<f64>,
    /// Sum of the actually returned cell heat, with its partition defect gated.
    pub resistance_heat_j: f64,
    pub absolute_heat_partition_defect_j: f64,
    /// delta_K + BE_loss + resistance_heat - body_endpoint_work
    /// - constraint_endpoint_work - pressure_source_work
    /// + fixed_boundary_pressure_work + open_boundary_pressure_work.
    pub energy_defect_j: f64,
    /// Non-cancelling local work defects, heat-partition defect and global defect.
    pub absolute_energy_defect_j: f64,
}

/// One native pressure/volume endpoint with physical stationary Darcy response.
#[derive(Clone, Debug, PartialEq)]
pub struct ResistiveDensityProjection {
    pub projection: DensityProjection,
    pub report: DensityResistanceReport,
}

impl FlowField {
    /// Apply body forcing, implicit Darcy resistance and pressure consistently.
    ///
    /// For every active native face, the same equation is used in prediction,
    /// pressure conductance and independent readback:
    /// `M*(u1-u0) = I_body + dt*A*delta_p - dt*B*u1`.
    /// M integrates ACTUAL density over the two adjacent half-cells; B integrates
    /// the caller's nonnegative `beta` [Pa s/m2] over those same volumes. Thus
    /// pressure conductance is `dt*A^2/(M+dt*B)`, while kinetic energy and all
    /// material inertia still use M. Zero resistance preserves the original
    /// force/projection face arithmetic.
    ///
    /// This supports sharp density/resistance jumps, graded cells, hydrostatic
    /// balance and the stiff Darcy limit without an explicit drag timestep
    /// restriction. Finite representable coefficients and the native residual,
    /// volume, impulse and work budgets remain required. The model is a frozen
    /// isotropic linear resistance against a stationary matrix; it supplies no
    /// constitutive permeability, moving skeleton, solid locking or viscosity.
    ///
    /// Physical resistance heat uses the ACCEPTED velocity and is split by each
    /// cell's own resistance/half-volume contribution. Body and controller work
    /// use that same endpoint. The separate backward-Euler loss is numerical
    /// and never becomes heat. Volume sources and mechanical pressure retain
    /// their original work conventions; density and enthalpy do not evolve here.
    ///
    /// # Errors
    /// Invalid source/geometry/resistance, incompatible volume, unrepresentable
    /// arithmetic, any native solve or balance refusal, or cancellation returns
    /// no candidate. All caller-owned fields remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn force_and_project_density_with_resistance(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: DensityForceInput<'_>,
        resistance_pa_s_m2: &[f64],
        target_outflow_m3_s: &[f64],
        projection_controls: DensityProjectionControls,
        force_controls: DensityForceControls,
        gate: &CancelGate,
    ) -> Result<ResistiveDensityProjection, ChtError> {
        self.force_and_project_density_with_resistance_readback(
            domain,
            faces,
            None,
            input,
            resistance_pa_s_m2,
            target_outflow_m3_s,
            projection_controls,
            force_controls,
            gate,
        )
        .map(|(step, _)| step)
    }

    /// The same implicit response with uniform planes or a borrowed spatial
    /// reservoir-pressure profile. Pressure acts through the original boundary
    /// half-cell mass AND resistance; either flow direction is admitted.
    /// Relative component gauges preserve small heads under a large reference.
    /// Prescribed nonzero Fixed flux requires controller work against drag,
    /// which is retained explicitly instead of deleting its physical heat.
    ///
    /// # Errors
    /// Conflicting pressure/velocity boundary roles or any original resistance
    /// step refusal returns no endpoint and preserves all supplied state.
    #[allow(clippy::too_many_arguments)]
    pub fn force_and_project_density_with_resistance_and_reservoir_pressure<'p>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: impl Into<ReservoirPressure<'p>>,
        input: DensityForceInput<'_>,
        resistance_pa_s_m2: &[f64],
        target_outflow_m3_s: &[f64],
        projection_controls: DensityProjectionControls,
        force_controls: DensityForceControls,
        gate: &CancelGate,
    ) -> Result<ResistiveDensityProjection, ChtError> {
        let pressure = reservoir_pressure.into().admit(faces, gate)?;
        self.force_and_project_density_with_resistance_readback(
            domain,
            faces,
            Some(pressure),
            input,
            resistance_pa_s_m2,
            target_outflow_m3_s,
            projection_controls,
            force_controls,
            gate,
        )
        .map(|(step, _)| step)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::conjugate::flow) fn force_and_project_density_with_resistance_readback(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure: Option<BoundaryPressure<'_>>,
        input: DensityForceInput<'_>,
        resistance_pa_s_m2: &[f64],
        target_outflow_m3_s: &[f64],
        projection_controls: DensityProjectionControls,
        force_controls: DensityForceControls,
        gate: &CancelGate,
    ) -> Result<(ResistiveDensityProjection, PressureReadback), ChtError> {
        let dt = projection_controls.duration_s;
        admit(self, domain, faces, input, dt, force_controls, gate)?;
        admit_resistance_field(domain, resistance_pa_s_m2, gate)?;
        if target_outflow_m3_s.len() != domain.cell_count() {
            return Err(refuse(
                "flow.resistance_shape",
                "one volume target per cell is required",
            ));
        }
        let mut predictor = self.clone();
        for cell in 0..domain.cell_count() {
            poll_every(gate, cell)?;
            if !domain.is_fluid(cell) {
                continue;
            }
            for face in 0..6 {
                let neighbor = fluid_neighbor(domain, cell, face);
                if !owns_face(face, neighbor) || !corrected_face(domain, faces, cell, face) {
                    continue;
                }
                let axis = face / 2;
                let area = domain.face_area(cell, axis);
                let mass = face_dual_mass(domain, input.density_kg_m3, cell, face, neighbor);
                finite_positive("flow.resistance_face_mass", mass)?;
                let body = face_impulse(domain, input, dt, cell, face, neighbor, mass);
                finite("flow.resistance_body_impulse", body)?;
                let resistance =
                    face_resistance_kg_s(domain, resistance_pa_s_m2, cell, face, neighbor);
                nonnegative("flow.resistance_face_coefficient", resistance)?;
                let (slot, _) = self.slot(domain, cell, face);
                let flux = if resistance == 0.0 {
                    // Preserve the original undamped predictor's arithmetic.
                    self.face(axis, slot) + area * (body / mass)
                } else {
                    let metric = mass + dt * resistance;
                    finite_positive("flow.resistance_momentum_metric", metric)?;
                    // Avoid constructing a potentially enormous unresisted
                    // velocity before applying stiff, representable damping.
                    (mass / metric) * self.face(axis, slot) + area * (body / metric)
                };
                finite("flow.resistance_predictor_flux", flux)?;
                *predictor.face_mut(axis, slot) = flux;
            }
        }
        let pressure = boundary_pressure.unwrap_or(BoundaryPressure::Planes([0.0; 6]));
        let (projection, readback) = predictor
            .project_density_volume_source_with_resistance_readback(
                domain,
                faces,
                pressure,
                input.density_kg_m3,
                target_outflow_m3_s,
                resistance_pa_s_m2,
                projection_controls,
                gate,
            )?;
        let report = audit_resistance(
            self,
            &projection,
            &readback,
            domain,
            faces,
            pressure,
            input,
            resistance_pa_s_m2,
            dt,
            force_controls,
            gate,
        )?;
        poll(gate)?;
        Ok((ResistiveDensityProjection { projection, report }, readback))
    }
}

fn add_heat(cell_heat: &mut f64, heat: f64, defect: &mut Dd) -> Result<(), ChtError> {
    nonnegative("flow.resistance_cell_heat", heat)?;
    let next = *cell_heat + heat;
    finite("flow.resistance_cell_heat_sum", next)?;
    *defect = *defect
        + Dd::from_f64(
            (Dd::from_f64(*cell_heat) + Dd::from_f64(heat) - Dd::from_f64(next))
                .to_f64()
                .abs(),
        );
    *cell_heat = next;
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn audit_resistance(
    before: &FlowField,
    projected: &DensityProjection,
    readback: &PressureReadback,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    pressure_boundary: BoundaryPressure<'_>,
    input: DensityForceInput<'_>,
    beta: &[f64],
    dt: f64,
    controls: DensityForceControls,
    gate: &CancelGate,
) -> Result<DensityResistanceReport, ChtError> {
    let mut applied = [Dd::ZERO; 3];
    let mut pressure = [Dd::ZERO; 3];
    let mut resistance = [Dd::ZERO; 3];
    let mut constraint = [Dd::ZERO; 3];
    let mut momentum = [Dd::ZERO; 3];
    let mut impulse_defect = Dd::ZERO;
    let mut initial = Dd::ZERO;
    let mut final_energy = Dd::ZERO;
    let mut change = Dd::ZERO;
    let mut body_work = Dd::ZERO;
    let mut constraint_work = Dd::ZERO;
    let mut numerical_loss = Dd::ZERO;
    let mut physical_heat = Dd::ZERO;
    let mut cell_heat = vec![0.0; domain.cell_count()];
    let mut heat_defect = Dd::ZERO;
    let mut energy_defect = Dd::ZERO;
    for cell in 0..domain.cell_count() {
        poll_every(gate, cell)?;
        if !domain.is_fluid(cell) {
            continue;
        }
        for axis in 0..3 {
            let impulse = dt
                * domain.volume(cell)
                * (input.density_kg_m3[cell] * input.acceleration_m_s2[axis]
                    + input.force_density_n_m3[cell][axis]);
            finite("flow.resistance_applied_impulse", impulse)?;
            applied[axis] = applied[axis] + Dd::from_f64(impulse);
        }
        for face in 0..6 {
            let neighbor = fluid_neighbor(domain, cell, face);
            if !owns_face(face, neighbor) {
                continue;
            }
            let axis = face / 2;
            let area = domain.face_area(cell, axis);
            let mass = face_dual_mass(domain, input.density_kg_m3, cell, face, neighbor);
            finite_positive("flow.resistance_audit_mass", mass)?;
            let coefficient = face_resistance_kg_s(domain, beta, cell, face, neighbor);
            nonnegative("flow.resistance_audit_coefficient", coefficient)?;
            let (slot, sign) = before.slot(domain, cell, face);
            let q0 = before.face(axis, slot);
            let q1 = projected.field.face(axis, slot);
            let u0 = q0 / area;
            let u1 = q1 / area;
            let du = (q1 - q0) / area;
            let body = face_impulse(domain, input, dt, cell, face, neighbor, mass);
            let drag = -(dt * coefficient) * u1;
            let active = corrected_face(domain, faces, cell, face);
            let pressure_impulse = if active {
                sign * dt
                    * area
                    * readback.difference(cell, neighbor, pressure_boundary.value(face, slot))
            } else {
                0.0
            };
            let held = if active { 0.0 } else { -(body + drag) };
            let delta_p = mass * du;
            let delta_k = mass * (u0 + 0.5 * du) * du;
            let loss = 0.5 * mass * du * du;
            let heat = -drag * u1;
            let applied_work = body * u1;
            let held_work = held * u1;
            let pressure_work = pressure_impulse * u1;
            for value in [
                u0,
                u1,
                du,
                body,
                drag,
                pressure_impulse,
                held,
                delta_p,
                delta_k,
                loss,
                heat,
                applied_work,
                held_work,
                pressure_work,
            ] {
                finite("flow.resistance_face_readback", value)?;
            }
            pressure[axis] = pressure[axis] + Dd::from_f64(pressure_impulse);
            resistance[axis] = resistance[axis] + Dd::from_f64(drag);
            constraint[axis] = constraint[axis] + Dd::from_f64(held);
            momentum[axis] = momentum[axis] + Dd::from_f64(delta_p);
            impulse_defect = impulse_defect
                + Dd::from_f64(
                    (Dd::from_f64(delta_p)
                        - Dd::from_f64(body)
                        - Dd::from_f64(drag)
                        - Dd::from_f64(pressure_impulse)
                        - Dd::from_f64(held))
                    .to_f64()
                    .abs(),
                );
            initial = initial + Dd::from_f64(0.5 * mass * u0 * u0);
            final_energy = final_energy + Dd::from_f64(0.5 * mass * u1 * u1);
            change = change + Dd::from_f64(delta_k);
            numerical_loss = numerical_loss + Dd::from_f64(loss);
            physical_heat = physical_heat + Dd::from_f64(heat);
            body_work = body_work + Dd::from_f64(applied_work);
            constraint_work = constraint_work + Dd::from_f64(held_work);
            let local = Dd::from_f64(delta_k) + Dd::from_f64(loss) + Dd::from_f64(heat)
                - Dd::from_f64(applied_work)
                - Dd::from_f64(held_work)
                - Dd::from_f64(pressure_work);
            energy_defect = energy_defect + Dd::from_f64(local.to_f64().abs());
            let mut partition = Dd::ZERO;
            for contributor in [Some(cell), neighbor].into_iter().flatten() {
                let b = area * (0.5 * domain.widths(contributor)[axis] * beta[contributor]);
                let heat = ((dt * b) * u1) * u1;
                add_heat(&mut cell_heat[contributor], heat, &mut heat_defect)?;
                partition = partition + Dd::from_f64(heat);
            }
            heat_defect =
                heat_defect + Dd::from_f64((partition - Dd::from_f64(heat)).to_f64().abs());
        }
    }
    for axis in 0..3 {
        impulse_defect = impulse_defect
            + Dd::from_f64(
                (momentum[axis]
                    - applied[axis]
                    - pressure[axis]
                    - resistance[axis]
                    - constraint[axis])
                    .to_f64()
                    .abs(),
            );
    }
    let mut represented_heat = Dd::ZERO;
    for (cell, &heat) in cell_heat.iter().enumerate() {
        poll_every(gate, cell)?;
        represented_heat = represented_heat + Dd::from_f64(heat);
    }
    heat_defect = heat_defect + Dd::from_f64((represented_heat - physical_heat).to_f64().abs());
    let report = &projected.report;
    let balance = change + numerical_loss + represented_heat
        - body_work
        - constraint_work
        - Dd::from_f64(report.pressure_source_work_j)
        + Dd::from_f64(report.fixed_boundary_pressure_work_j)
        + Dd::from_f64(report.open_boundary_pressure_work_j);
    energy_defect = energy_defect + heat_defect + Dd::from_f64(balance.to_f64().abs());
    for value in [
        initial,
        final_energy,
        change,
        numerical_loss,
        represented_heat,
        body_work,
        constraint_work,
        balance,
        impulse_defect,
        heat_defect,
        energy_defect,
    ] {
        finite("flow.resistance_account", value.to_f64())?;
    }
    if impulse_defect.to_f64() > controls.maximum_impulse_defect_kg_m_s {
        return Err(refuse(
            "flow.resistance_impulse_budget",
            "implicit resistance exceeds the impulse defect budget",
        ));
    }
    if energy_defect.to_f64() > controls.maximum_energy_defect_j {
        return Err(refuse(
            "flow.resistance_energy_budget",
            "implicit resistance exceeds the physical work/heat budget",
        ));
    }
    Ok(DensityResistanceReport {
        applied_impulse_kg_m_s: applied.map(Dd::to_f64),
        pressure_impulse_kg_m_s: pressure.map(Dd::to_f64),
        resistance_impulse_kg_m_s: resistance.map(Dd::to_f64),
        constraint_impulse_kg_m_s: constraint.map(Dd::to_f64),
        momentum_change_kg_m_s: momentum.map(Dd::to_f64),
        absolute_impulse_defect_kg_m_s: impulse_defect.to_f64(),
        kinetic_energy_before_j: initial.to_f64(),
        kinetic_energy_after_j: final_energy.to_f64(),
        kinetic_energy_change_j: change.to_f64(),
        body_endpoint_work_j: body_work.to_f64(),
        constraint_endpoint_work_j: constraint_work.to_f64(),
        backward_euler_loss_j: numerical_loss.to_f64(),
        cell_resistance_heat_j: cell_heat,
        resistance_heat_j: represented_heat.to_f64(),
        absolute_heat_partition_defect_j: heat_defect.to_f64(),
        energy_defect_j: balance.to_f64(),
        absolute_energy_defect_j: energy_defect.to_f64(),
    })
}

#[cfg(test)]
mod tests;
