//! Extensive cell momentum transferred onto the original pressure face masses.

use super::{
    CancelGate, ChtError, FlowFace, FlowField, Sum, VoxelDomain, face_dual_mass, finite,
    finite_positive, nonnegative, poll, poll_every, refuse,
};

/// Allocation and independent physical readback limits for momentum transfer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DensityMomentumReconstructionControls {
    /// Maximum domain cells, including solids, before allocation.
    pub max_cells: usize,
    /// Sum of absolute local native-versus-half-cell mass differences and
    /// global per-axis native-versus-cell mass differences, kg.
    pub maximum_mass_defect_kg: f64,
    /// Sum of local impulse defects and the global vector L1 defect, kg m/s.
    pub maximum_momentum_defect_kg_m_s: f64,
    /// Sum of absolute native-mass representation energy, local corrected
    /// work-identity defects, and the global corrected defect, J. Numerical
    /// mass representation cannot be accepted merely by naming its account.
    pub maximum_energy_defect_j: f64,
}

/// Readback of the returned native fluxes against authoritative cell M and P.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityMomentumReconstructionReport {
    /// Original fluid mass, counted once over cells, kg.
    pub cell_mass_kg: f64,
    /// Sum of native face dual masses separately for x, y and z, kg. Each
    /// component includes every fluid half-cell once, including held faces.
    pub native_face_mass_kg: [f64; 3],
    pub absolute_mass_defect_kg: f64,
    /// Original extensive cell momentum in world axes, kg m/s.
    pub cell_momentum_kg_m_s: [f64; 3],
    /// Actual `sum(M_native Q/A)` in world axes, kg m/s.
    pub face_momentum_kg_m_s: [f64; 3],
    /// Impulse on the fluid needed to retain `Fixed`, `Wall` and solid-face
    /// velocities, kg m/s. This is a discrete support/controller account,
    /// not a spatially resolved wall traction or a force over a known time.
    pub constraint_impulse_kg_m_s: [f64; 3],
    pub absolute_momentum_defect_kg_m_s: f64,
    /// Original `sum(P dot (P/M))/2`, J.
    pub kinetic_energy_before_j: f64,
    /// Native face-dual kinetic energy read from the returned field, J.
    pub kinetic_energy_after_j: f64,
    /// Face minus cell energy, evaluated about local velocity datums, J.
    pub kinetic_energy_change_j: f64,
    /// Nonnegative Jensen loss from averaging the two half-cell velocities,
    /// J. Uses authoritative M/2 weights and relative donor velocities.
    pub numerical_mixing_loss_j: f64,
    /// `sum(I_constraint u_held)`, J, supplied by the prescribed-velocity
    /// reservoirs. Stationary walls and solid interfaces perform zero work.
    pub constraint_reservoir_work_j: f64,
    /// Nonnegative `sum(M_native (u_held-u_unconstrained)^2/2)`, J.
    pub numerical_constraint_loss_j: f64,
    /// Signed native-inertia representation change at fixed face momentum:
    /// `sum((M_half-M_native) (p/M_native) (p/M_half)/2)`, J. This accounts
    /// for deriving rho=M/V and then using the original pressure dual mass;
    /// it is roundoff, never mixing, material dissipation or heat.
    pub mass_representation_energy_j: f64,
    /// `delta_K + mixing_loss + constraint_loss - reservoir_work
    /// - mass_representation_energy`, J. Evaluated with centered local
    /// velocities so a large common velocity does not hide the residual.
    pub energy_defect_j: f64,
    /// Includes absolute mass-representation energy as well as remaining
    /// local/global kinetic defects, J; see the corresponding control.
    pub absolute_energy_defect_j: f64,
}

/// Conservative cell-to-face transfer, before the physical owner's projection.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityMomentumReconstruction {
    /// Reconstructed native volumetric face fluxes, m^3/s. This field has NOT
    /// been pressure-projected and need not satisfy a volume constraint yet.
    pub field: FlowField,
    /// Exact admitted M/V values used to compute native face inertia, kg/m^3;
    /// zero in solids. Reuse these values for the subsequent density stage.
    pub density_kg_m3: Vec<f64>,
    pub report: DensityMomentumReconstructionReport,
}

impl FlowField {
    /// Transfer authoritative cell mass and momentum onto native face inertia.
    /// Each cell contributes M/2 and P_axis/2 to each of its two faces in that
    /// axis. An internal face merges those extensive contributions; it does
    /// not geometrically interpolate velocities across a density jump.
    /// Donor velocity is piecewise constant within each original cell.
    /// The returned velocity divides that momentum by the SAME dual mass used
    /// by native pressure, including graded grids and boundary half-cells.
    ///
    /// `self` supplies held boundary fluxes. Its internal and `Free` fluxes
    /// are replaced; `Fixed`, wall and fluid-solid velocities stay unchanged.
    /// Their reactions, prescribed-reservoir work and constraint-projection
    /// losses are explicit. Every fluid half-cell participates, including at
    /// walls and on either orientation of a solid interface.
    ///
    /// Supply actual returned material masses and extensive momenta, such as
    /// `MaterialAdvectionStep::mass_kg` and `MaterialMomentumStep::momentum_kg_m_s`.
    /// This method validates the supplied inventories; it does not infer a
    /// private association between those independently public values. The
    /// supplied domain must be their original geometry and occupancy.
    ///
    /// The output is unprojected so the existing force, material or phase owner
    /// can apply its one compatible density pressure solve. No pressure impulse
    /// is applied to the cell momenta, and no new density evolution is invented.
    /// This is a one-way conservative transfer, not full staggered advection,
    /// a bidirectional momentum coupling, or a Navier-Stokes time integrator.
    /// Repeated reconstruction introduces the reported averaging diffusion;
    /// all numerical losses and mass-representation defects remain out of heat.
    ///
    /// # Errors
    /// Invalid shapes, nonpositive fluid mass/density, solid inventory, invalid
    /// held fluxes, non-finite or unrepresentable arithmetic, exceeded explicit
    /// mass/momentum/energy budgets, or cancellation. Inputs never change.
    #[allow(clippy::too_many_arguments)]
    pub fn reconstruct_material_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        mass_kg: &[f64],
        momentum_kg_m_s: &[[f64; 3]],
        controls: DensityMomentumReconstructionControls,
        gate: &CancelGate,
    ) -> Result<DensityMomentumReconstruction, ChtError> {
        let density = admit(
            self,
            domain,
            faces,
            mass_kg,
            momentum_kg_m_s,
            controls,
            gate,
        )?;
        let mut field = self.clone();
        for c in 0..domain.cell_count() {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                continue;
            }
            for f in 0..6 {
                let nb = fluid_neighbor(domain, c, f);
                if !owns_face(f, nb) || held_face(domain, faces, c, f) {
                    continue;
                }
                let axis = f / 2;
                let native_mass = face_dual_mass(domain, &density, c, f, nb);
                finite_positive("flow.reconstruction_face_mass_kg", native_mass)?;
                let p = half_momentum(momentum_kg_m_s, c, axis, nb);
                finite("flow.reconstruction_face_momentum", p)?;
                let q = domain.face_area(c, axis) * (p / native_mass);
                finite("flow.reconstructed_face_flux", q)?;
                let (slot, _) = self.slot(domain, c, f);
                *field.face_mut(axis, slot) = q;
            }
        }
        let report = audit(
            self,
            &field,
            domain,
            faces,
            mass_kg,
            momentum_kg_m_s,
            &density,
            controls,
            gate,
        )?;
        poll(gate)?;
        Ok(DensityMomentumReconstruction {
            field,
            density_kg_m3: density,
            report,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn admit(
    field: &FlowField,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    mass_kg: &[f64],
    momentum: &[[f64; 3]],
    controls: DensityMomentumReconstructionControls,
    gate: &CancelGate,
) -> Result<Vec<f64>, ChtError> {
    poll(gate)?;
    if domain.cell_count() > controls.max_cells {
        return Err(refuse(
            "flow.reconstruction_cell_budget",
            "domain exceeds the explicit reconstruction cell budget",
        ));
    }
    for (name, value) in [
        (
            "flow.reconstruction_mass_budget",
            controls.maximum_mass_defect_kg,
        ),
        (
            "flow.reconstruction_momentum_budget",
            controls.maximum_momentum_defect_kg_m_s,
        ),
        (
            "flow.reconstruction_energy_budget",
            controls.maximum_energy_defect_j,
        ),
    ] {
        nonnegative(name, value)?;
    }
    if mass_kg.len() != domain.cell_count() || momentum.len() != domain.cell_count() {
        return Err(refuse(
            "flow.reconstruction_shape",
            "one mass and momentum vector per original grid cell are required",
        ));
    }
    let mut density = vec![0.0; domain.cell_count()];
    for c in 0..domain.cell_count() {
        poll_every(gate, c)?;
        finite("flow.reconstruction_cell_mass_kg", mass_kg[c])?;
        for &p in &momentum[c] {
            finite("flow.reconstruction_cell_momentum", p)?;
        }
        if domain.is_fluid(c) {
            finite_positive("flow.reconstruction_half_mass_kg", 0.5 * mass_kg[c])?;
            density[c] = mass_kg[c] / domain.volume(c);
            finite_positive("flow.reconstruction_density_kg_m3", density[c])?;
            for &p in &momentum[c] {
                finite("flow.reconstruction_cell_velocity", p / mass_kg[c])?;
            }
        } else if mass_kg[c] != 0.0 || momentum[c] != [0.0; 3] {
            return Err(refuse(
                "flow.reconstruction_solid_inventory",
                "solid cells cannot carry fluid mass or momentum",
            ));
        }
    }
    field.admit_density_field(domain, faces, &density, None, gate)?;
    Ok(density)
}

fn fluid_neighbor(domain: &VoxelDomain, c: usize, f: usize) -> Option<usize> {
    domain.neighbor(c, f).filter(|&n| domain.is_fluid(n))
}

fn owns_face(f: usize, nb: Option<usize>) -> bool {
    nb.is_none() || f % 2 == 1
}

fn held_face(domain: &VoxelDomain, faces: [FlowFace; 6], c: usize, f: usize) -> bool {
    domain
        .neighbor(c, f)
        .map_or(faces[f] != FlowFace::Free, |n| !domain.is_fluid(n))
}

fn half_momentum(p: &[[f64; 3]], c: usize, axis: usize, nb: Option<usize>) -> f64 {
    0.5 * p[c][axis] + nb.map_or(0.0, |n| 0.5 * p[n][axis])
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn audit(
    held: &FlowField,
    field: &FlowField,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    mass: &[f64],
    momentum: &[[f64; 3]],
    density: &[f64],
    controls: DensityMomentumReconstructionControls,
    gate: &CancelGate,
) -> Result<DensityMomentumReconstructionReport, ChtError> {
    let mut cell_mass = Sum::default();
    let mut native_mass: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut mass_defect = Sum::default();
    let mut cell_p: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut face_p: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut reaction: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut momentum_defect = Sum::default();
    let mut before = Sum::default();
    let mut after = Sum::default();
    let mut change = Sum::default();
    let mut mixing = Sum::default();
    let mut work = Sum::default();
    let mut constraint_loss = Sum::default();
    let mut representation = Sum::default();
    let mut energy_defect = Sum::default();
    let mut absolute_energy_defect = Sum::default();
    for c in 0..domain.cell_count() {
        poll_every(gate, c)?;
        if !domain.is_fluid(c) {
            continue;
        }
        // Independent cell inventory totals detect missing/duplicate halves.
        cell_mass.add(mass[c]);
        for axis in 0..3 {
            cell_p[axis].add(momentum[c][axis]);
            before.add(0.5 * momentum[c][axis] * (momentum[c][axis] / mass[c]));
        }
        for f in 0..6 {
            let nb = fluid_neighbor(domain, c, f);
            if !owns_face(f, nb) {
                continue;
            }
            let axis = f / 2;
            let ma = 0.5 * mass[c];
            let mb = nb.map_or(0.0, |n| 0.5 * mass[n]);
            let original_mass = ma + mb;
            let represented_mass = face_dual_mass(domain, density, c, f, nb);
            finite_positive("flow.reconstruction_half_sum_kg", original_mass)?;
            finite_positive("flow.reconstruction_face_mass_kg", represented_mass)?;
            let p = half_momentum(momentum, c, axis, nb);
            let original_velocity = p / original_mass;
            let unconstrained = p / represented_mass;
            let reference = momentum[c][axis] / mass[c];
            let other = nb.map_or(reference, |n| momentum[n][axis] / mass[n]);
            let (slot, _) = field.slot(domain, c, f);
            let constrained = held_face(domain, faces, c, f);
            if constrained && field.face(axis, slot) != held.face(axis, slot) {
                return Err(refuse(
                    "flow.reconstruction_held_flux",
                    "reconstruction changed a held boundary or solid flux",
                ));
            }
            let actual = field.face(axis, slot) / domain.face_area(c, axis);
            let actual_p = represented_mass * actual;
            let delta_p = actual_p - p;
            let delta_mass = represented_mass - original_mass;
            // Smaller half mass times a ratio >= 1/2 avoids an unnecessary
            // product of masses; differences preserve a common velocity datum.
            let reduced_mass = ma.min(mb) * (ma.max(mb) / original_mass);
            let difference = other - reference;
            let mix = 0.5 * reduced_mass * difference * difference;
            let constraint_du = if constrained {
                actual - unconstrained
            } else {
                0.0
            };
            let impulse = represented_mass * constraint_du;
            let loss = 0.5 * represented_mass * constraint_du * constraint_du;
            let mass_energy = -0.5 * delta_mass * unconstrained * original_velocity;
            let local_momentum_defect = delta_p - impulse;

            // K_after-K_before = K_relative + u_ref*delta_P
            //                         - u_ref^2*delta_M/2.
            // Apply the same centering to reservoir work before taking the
            // residual. Subtracting absolute kinetic totals would erase small
            // mixing losses at a large common velocity.
            let relative = actual - reference;
            let relative_change =
                0.5 * represented_mass * relative * relative - 0.5 * mb * difference * difference;
            let mass_datum = 0.5 * delta_mass * reference * reference;
            let physical_change = reference.mul_add(delta_p, relative_change - mass_datum);
            let local_defect = reference.mul_add(
                local_momentum_defect,
                relative_change - impulse * relative + mix + loss - mass_datum - mass_energy,
            );
            native_mass[axis].add(represented_mass);
            mass_defect.add(delta_mass.abs());
            face_p[axis].add(actual_p);
            reaction[axis].add(impulse);
            momentum_defect.add(local_momentum_defect.abs());
            after.add(0.5 * represented_mass * actual * actual);
            change.add(physical_change);
            mixing.add(mix);
            work.add(impulse * actual);
            constraint_loss.add(loss);
            representation.add(mass_energy);
            energy_defect.add(local_defect);
            absolute_energy_defect.add(mass_energy.abs());
            absolute_energy_defect.add(local_defect.abs());
        }
    }
    let native_mass = native_mass.map(|s| s.value);
    let cell_p = cell_p.map(|s| s.value);
    let face_p = face_p.map(|s| s.value);
    let reaction = reaction.map(|s| s.value);
    for axis in 0..3 {
        for value in [
            native_mass[axis],
            cell_p[axis],
            face_p[axis],
            reaction[axis],
        ] {
            finite("flow.reconstruction_inventory", value)?;
        }
        mass_defect.add((native_mass[axis] - cell_mass.value).abs());
        momentum_defect.add((face_p[axis] - cell_p[axis] - reaction[axis]).abs());
    }
    absolute_energy_defect.add(energy_defect.value.abs());
    for value in [
        cell_mass.value,
        before.value,
        after.value,
        change.value,
        mixing.value,
        work.value,
        constraint_loss.value,
        representation.value,
        energy_defect.value,
        absolute_energy_defect.value,
    ] {
        finite("flow.reconstruction_energy", value)?;
    }
    check_budget(
        "flow.reconstruction_mass_budget",
        mass_defect.value,
        controls.maximum_mass_defect_kg,
    )?;
    check_budget(
        "flow.reconstruction_momentum_budget",
        momentum_defect.value,
        controls.maximum_momentum_defect_kg_m_s,
    )?;
    check_budget(
        "flow.reconstruction_energy_budget",
        absolute_energy_defect.value,
        controls.maximum_energy_defect_j,
    )?;
    Ok(DensityMomentumReconstructionReport {
        cell_mass_kg: cell_mass.value,
        native_face_mass_kg: native_mass,
        absolute_mass_defect_kg: mass_defect.value,
        cell_momentum_kg_m_s: cell_p,
        face_momentum_kg_m_s: face_p,
        constraint_impulse_kg_m_s: reaction,
        absolute_momentum_defect_kg_m_s: momentum_defect.value,
        kinetic_energy_before_j: before.value,
        kinetic_energy_after_j: after.value,
        kinetic_energy_change_j: change.value,
        numerical_mixing_loss_j: mixing.value,
        constraint_reservoir_work_j: work.value,
        numerical_constraint_loss_j: constraint_loss.value,
        mass_representation_energy_j: representation.value,
        energy_defect_j: energy_defect.value,
        absolute_energy_defect_j: absolute_energy_defect.value,
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
