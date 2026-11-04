//! Molecular momentum diffusion on the original MAC stress and Krylov owners.
//!
//! Pressure, convection and calorics remain separate physical substeps. This
//! predictor uses actual M/V face inertia and the existing full symmetric
//! molecular stress, including harmonic material jumps and half-cell walls.

use super::*;
use fs_math::dd::Dd;

/// Work and independent physical-unit gates for one native viscosity substep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialViscosityControls {
    pub max_cells: usize,
    /// Each sweep uses the original three momentum solves, each bounded by
    /// its native 20,000 Krylov iterations, and one frozen residual readback.
    pub max_sweeps: usize,
    pub relative_tolerance: f64,
    pub maximum_impulse_defect_kg_m_s: f64,
    pub maximum_angular_impulse_defect_kg_m2_s: f64,
    pub maximum_energy_defect_j: f64,
}

/// Physical molecular heat and the distinct numerical/mechanical accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialViscosityReport {
    pub sweeps: usize,
    pub krylov_iterations: usize,
    pub relative_residual: f64,
    /// Native extensive impulse and its actual feedback into original cells.
    pub face_impulse_kg_m_s: [f64; 3],
    pub cell_momentum_change_kg_m_s: [f64; 3],
    pub absolute_impulse_defect_kg_m_s: f64,
    pub face_angular_impulse_kg_m2_s: [f64; 3],
    pub cell_angular_momentum_change_kg_m2_s: [f64; 3],
    pub absolute_angular_impulse_defect_kg_m2_s: f64,
    /// Exactly duration times the original molecular-dissipation quadrature.
    /// The composition deposits only this physical channel into enthalpy.
    pub cell_molecular_heat_j: Vec<f64>,
    pub molecular_heat_j: f64,
    pub wall_work_by_face_j: [f64; 6],
    pub wall_work_j: f64,
    pub native_kinetic_energy_change_j: f64,
    pub cell_kinetic_energy_change_j: f64,
    /// Backward-Euler numerical loss on actual native face masses.
    pub backward_euler_loss_j: f64,
    /// Half-cell scatter loss; neither numerical loss is material heat.
    pub numerical_scatter_loss_j: f64,
    /// Numerical cell/native-inertia and representable-impulse work. Its
    /// absolute face contributions remain inside the energy defect budget.
    pub representation_work_j: f64,
    /// Cell delta_K + molecular_heat + BE_loss + scatter_loss
    /// - wall_work - representation_work.
    pub energy_defect_j: f64,
    pub absolute_energy_defect_j: f64,
}

pub(in crate::conjugate) struct NativeMaterialViscosityInput<'a> {
    pub field: &'a FlowField,
    pub mass_kg: &'a [f64],
    pub momentum_kg_m_s: &'a [[f64; 3]],
    pub viscosity_pa_s: &'a [f64],
    /// Produced by the ORIGINAL density projection's face_dual_mass.
    pub face_mass_kg: [Vec<f64>; 3],
    pub boundaries: [FvBoundary; 6],
    /// One opposed normal-traction pair, with zero viscous-stage pressure and
    /// zero tangential trace. The following projection owns reservoir loading.
    pub pressure_mouth_axis: Option<usize>,
    pub duration_s: f64,
}

pub(in crate::conjugate) struct NativeMaterialViscosityStep {
    pub momentum_kg_m_s: Vec<[f64; 3]>,
    pub report: MaterialViscosityReport,
}

fn refuse(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.material_viscosity",
        reason: reason.into(),
    }
}

fn budget(value: Dd, maximum: f64, reason: &str) -> Result<f64, ChtError> {
    let value = value.to_f64();
    finite("flow.material_viscosity_readback", value)?;
    if value.abs() > maximum {
        return Err(refuse(reason));
    }
    Ok(value)
}

fn magnitude(value: Dd) -> Dd {
    Dd::from_f64(value.to_f64().abs())
}

fn torque(position: [f64; 3], impulse: [Dd; 3]) -> [Dd; 3] {
    let r = position.map(Dd::from_f64);
    [
        r[1] * impulse[2] - r[2] * impulse[1],
        r[2] * impulse[0] - r[0] * impulse[2],
        r[0] * impulse[1] - r[1] * impulse[0],
    ]
}

#[allow(clippy::too_many_lines)] // one bounded native solve and acceptance
pub(in crate::conjugate) fn kick(
    domain: &VoxelDomain,
    input: NativeMaterialViscosityInput<'_>,
    controls: MaterialViscosityControls,
    gate: &CancelGate,
) -> Result<NativeMaterialViscosityStep, ChtError> {
    poll(gate)?;
    let n = domain.cell_count();
    if n > controls.max_cells
        || controls.max_sweeps == 0
        || input.mass_kg.len() != n
        || input.momentum_kg_m_s.len() != n
        || input.viscosity_pa_s.len() != n
    {
        return Err(refuse("viscosity requires bounded original cells"));
    }
    if input.pressure_mouth_axis.is_some_and(|a| a >= 3)
        || input.boundaries.iter().enumerate().any(|(f, b)| {
            if input.pressure_mouth_axis == Some(f / 2) {
                !matches!(b, FvBoundary::Outlet)
            } else {
                !matches!(b, FvBoundary::Wall { .. } | FvBoundary::Symmetry)
            }
        })
    {
        return Err(refuse(
            "viscosity needs closed walls/symmetry or one opposed zero-pressure traction pair",
        ));
    }
    finite_positive("flow.viscosity_duration", input.duration_s)?;
    finite_positive("flow.viscosity_tolerance", controls.relative_tolerance)?;
    for bound in [
        controls.maximum_impulse_defect_kg_m_s,
        controls.maximum_angular_impulse_defect_kg_m2_s,
        controls.maximum_energy_defect_j,
    ] {
        if !bound.is_finite() || bound < 0.0 {
            return Err(refuse("invalid viscosity balance budget"));
        }
    }
    let first = (0..n)
        .find(|&c| domain.is_fluid(c))
        .ok_or_else(|| refuse("viscosity requires mobile material"))?;
    let reference_density = input.mass_kg[first] / domain.volume(first);
    finite_positive("flow.viscosity_reference_density", reference_density)?;
    let mut config = SimpleConfig::new(input.boundaries);
    config.max_iterations = controls.max_sweeps;
    config.tolerance = controls.relative_tolerance;
    admit_molecular_boundaries(domain, &config)?;
    let mut solver = Solver::with_molecular_coefficients(
        domain,
        reference_density,
        input.viscosity_pa_s[first],
        &config,
    );
    solver.set_material_transport(input.viscosity_pa_s, None, gate)?;
    solver.pressure_work_axis = input.pressure_mouth_axis;
    // This is only molecular diffusion. The enclosing owner subsequently
    // applies body force, pressure and material convection, each once.
    solver.creeping = true;
    let mut face_mass = input.face_mass_kg;
    for (a, fluxes) in input.field.face_arrays().into_iter().enumerate() {
        if fluxes.len() != solver.vel[a].len() || face_mass[a].len() != fluxes.len() {
            return Err(refuse(
                "viscosity face arrays differ from the original grid",
            ));
        }
        for (i, &flux) in fluxes.iter().enumerate() {
            if i % 256 == 0 {
                poll(gate)?;
            }
            let velocity = flux / solver.comp_face_area(a, i);
            finite("flow.viscosity_initial_velocity", velocity)?;
            match solver.comps[a].kind[i] {
                Kind::Fixed(required) if velocity != required => {
                    return Err(refuse(
                        "viscosity cannot change an initial held normal velocity",
                    ));
                }
                Kind::Unknown(_) => {
                    finite_positive("flow.viscosity_face_mass", face_mass[a][i])?;
                }
                Kind::Outlet(_) => {
                    finite_positive("flow.viscosity_face_mass", face_mass[a][i])?;
                    // Existing traction rows mirror the interior cell: every
                    // stress coefficient is twice its physical half-cell value.
                    // Mirror inertia by the SAME factor; audit uses weight 1/2.
                    face_mass[a][i] *= 2.0;
                    finite_positive("flow.viscosity_mirrored_face_mass", face_mass[a][i])?;
                }
                Kind::Fixed(_) => {}
            }
            solver.vel[a][i] = velocity;
        }
    }
    solver.time = Some(TimeLevels {
        dt: input.duration_s,
        scheme: TimeScheme::BackwardEuler,
        previous: solver.vel.clone(),
        older: None,
        face_mass_kg: Some(face_mass),
    });
    let mut relative = f64::INFINITY;
    let mut residual_impulse = f64::INFINITY;
    for sweep in 1..=controls.max_sweeps {
        poll(gate)?;
        solver.sweeps = sweep;
        for a in 0..3 {
            solver.momentum(a, true, gate)?;
        }
        solver.guard_finite()?;
        let mut rms = 0.0_f64;
        let mut absolute = Dd::ZERO;
        // All three residuals see the SAME final coupled velocity, including
        // reciprocal cross derivatives. No pressure correction is performed.
        for a in 0..3 {
            let residual = solver.momentum(a, false, gate)?;
            rms = rms.max(residual.rms);
            absolute = absolute + Dd::from_f64(residual.absolute_force_n);
        }
        let scale = solver.momentum_rms.iter().copied().fold(0.0_f64, f64::max);
        relative = if scale > 0.0 { rms / scale } else { rms };
        residual_impulse = (Dd::from_f64(input.duration_s) * absolute).to_f64();
        finite("flow.viscosity_frozen_residual", relative)?;
        finite("flow.viscosity_frozen_impulse", residual_impulse)?;
        if relative <= controls.relative_tolerance
            && residual_impulse <= controls.maximum_impulse_defect_kg_m_s
        {
            break;
        }
    }
    if relative > controls.relative_tolerance
        || residual_impulse > controls.maximum_impulse_defect_kg_m_s
    {
        return Err(ChtError::SolverNotConverged {
            system: "material molecular momentum",
            iterations: solver.sweeps,
            relative_residual: relative,
            tolerance: controls.relative_tolerance,
        });
    }
    audit(
        &solver,
        input.mass_kg,
        input.momentum_kg_m_s,
        controls,
        relative,
        residual_impulse,
        gate,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn audit(
    solver: &Solver<'_>,
    mass: &[f64],
    original: &[[f64; 3]],
    controls: MaterialViscosityControls,
    relative: f64,
    residual_impulse: f64,
    gate: &CancelGate,
) -> Result<NativeMaterialViscosityStep, ChtError> {
    let domain = solver.domain;
    let time = solver.time.as_ref().expect("native viscosity time level");
    let inertia = time.face_mass_kg.as_ref().expect("actual native inertia");
    let duration = Dd::from_f64(time.dt);
    let half = Dd::from_f64(0.5);
    let power = solver.sample_molecular_power(gate)?;
    let mut cell_heat = Vec::with_capacity(domain.cell_count());
    let mut heat = Dd::ZERO;
    for (c, watts) in power.heat.iter().copied().enumerate() {
        if c % 256 == 0 {
            poll(gate)?;
        }
        let joules = (duration * Dd::from_f64(watts)).to_f64();
        finite("flow.viscosity_heat", joules)?;
        if joules < 0.0 || (!domain.is_fluid(c) && joules != 0.0) {
            return Err(refuse(
                "molecular heat must be nonnegative and mobile-cell owned",
            ));
        }
        cell_heat.push(joules);
        heat = heat + Dd::from_f64(joules);
    }
    let wall_by_face = power
        .wall_power_w
        .map(|w| (duration * Dd::from_f64(w)).to_f64());
    let wall = wall_by_face
        .into_iter()
        .fold(Dd::ZERO, |s, w| s + Dd::from_f64(w));
    let partition = magnitude(heat - duration * Dd::from_f64(power.quadrature_w))
        + magnitude(heat - duration * Dd::from_f64(power.molecular_heat_w));
    let mut impulses: [Vec<f64>; 3] = std::array::from_fn(|a| vec![0.0; solver.vel[a].len()]);
    let mut face_impulse = [Dd::ZERO; 3];
    let mut face_angular = [Dd::ZERO; 3];
    let mut native_change = Dd::ZERO;
    let mut time_loss = Dd::ZERO;
    let mut representation = Dd::ZERO;
    let mut absolute_representation = Dd::ZERO;
    let mut impulse_defect = Dd::from_f64(residual_impulse);
    for a in 0..3 {
        for &i in &solver.comps[a].unknowns {
            poll(gate)?;
            let at = solver.comps[a].coords(i);
            let (minus, plus, weight) = solver.power_face_cells(a, i);
            let m = cell_of(solver.n, minus);
            let p = cell_of(solver.n, plus);
            let old = Dd::from_f64(time.previous[a][i]);
            let delta = Dd::from_f64(solver.vel[a][i]) - old;
            let native_mass = Dd::from_f64(weight * inertia[a][i]);
            let physical_impulse = native_mass * delta;
            let represented = physical_impulse.to_f64();
            finite("flow.viscosity_impulse", represented)?;
            impulses[a][i] = represented;
            let impulse = Dd::from_f64(represented);
            impulse_defect = impulse_defect + magnitude(impulse - physical_impulse);
            face_impulse[a] = face_impulse[a] + impulse;
            let mut position = domain.center(minus[0], minus[1], minus[2]);
            position[a] = domain.face_coord(a, at[a]);
            let mut vector = [Dd::ZERO; 3];
            vector[a] = impulse;
            for (sum, value) in face_angular.iter_mut().zip(torque(position, vector)) {
                *sum = *sum + value;
            }
            let change = native_mass * delta * (old + half * delta);
            native_change = native_change + change;
            time_loss = time_loss + half * native_mass * delta * delta;
            let (half_mass, half_momentum) = if weight == 0.5 {
                (
                    Dd::from_f64(0.5 * mass[m]),
                    half * Dd::from_f64(original[m][a]),
                )
            } else {
                (
                    Dd::from_f64(0.5 * mass[m] + 0.5 * mass[p]),
                    half * (Dd::from_f64(original[m][a]) + Dd::from_f64(original[p][a])),
                )
            };
            // Virtual common face work minus actual native kinetic change.
            // This keeps M/V and Q/A/impulse representation work explicit,
            // including low components before a large velocity datum cancels.
            let difference = (half_momentum + half * impulse) * (impulse / half_mass) - change;
            representation = representation + difference;
            absolute_representation = absolute_representation + magnitude(difference);
        }
    }
    let mut momentum = original.to_vec();
    let mut cell_change = [Dd::ZERO; 3];
    let mut cell_angular = [Dd::ZERO; 3];
    let mut cell_energy = Dd::ZERO;
    let mut scatter = Dd::ZERO;
    for c in 0..domain.cell_count() {
        if c % 256 == 0 {
            poll(gate)?;
        }
        if !domain.is_fluid(c) {
            continue;
        }
        let at = domain.coords(c);
        let mut increment = [Dd::ZERO; 3];
        for a in 0..3 {
            let kick = |plus: bool| {
                let f = 2 * a + usize::from(plus);
                let index = solver.cell_face(a, at, plus);
                if impulses[a][index] == 0.0 {
                    0.0
                } else {
                    let half_mass = domain.neighbor(c, f).map_or(0.5 * mass[c], |neighbor| {
                        0.5 * mass[c] + 0.5 * mass[neighbor]
                    });
                    impulses[a][index] / half_mass
                }
            };
            let left = kick(false);
            let right = kick(true);
            let intended = (0.5 * mass[c]) * left + (0.5 * mass[c]) * right;
            if intended != 0.0 {
                momentum[c][a] += intended;
            }
            finite("flow.viscous_cell_momentum", momentum[c][a])?;
            let old = Dd::from_f64(original[c][a]);
            let new = Dd::from_f64(momentum[c][a]);
            let actual = new - old;
            increment[a] = actual;
            impulse_defect = impulse_defect + magnitude(actual - Dd::from_f64(intended));
            cell_change[a] = cell_change[a] + actual;
            cell_energy = cell_energy + half * (new + old) * (actual / Dd::from_f64(mass[c]));
            let difference = Dd::from_f64(left) - Dd::from_f64(right);
            scatter = scatter + Dd::from_f64(mass[c] / 8.0) * difference * difference;
        }
        for (sum, value) in cell_angular
            .iter_mut()
            .zip(torque(domain.center(at[0], at[1], at[2]), increment))
        {
            *sum = *sum + value;
        }
    }
    let mut angular_defect = Dd::ZERO;
    for a in 0..3 {
        impulse_defect = impulse_defect + magnitude(cell_change[a] - face_impulse[a]);
        angular_defect = angular_defect + magnitude(cell_angular[a] - face_angular[a]);
        // Bound the original row residual's possible moment about the domain
        // origin as well as checking actual transfer of angular impulse.
        angular_defect = angular_defect
            + Dd::from_f64(domain.face_coord(a, domain.dims()[a])) * Dd::from_f64(residual_impulse);
    }
    let native_defect = native_change + time_loss + heat - wall;
    let feedback_defect = cell_energy + scatter - native_change - representation;
    let defect = cell_energy + time_loss + scatter + heat - wall - representation;
    let energy_defect = magnitude(native_defect)
        + magnitude(feedback_defect)
        + magnitude(defect)
        + absolute_representation
        + partition;
    let impulse_defect = budget(
        impulse_defect,
        controls.maximum_impulse_defect_kg_m_s,
        "viscous impulse feedback exceeds its independent budget",
    )?;
    let angular_defect = budget(
        angular_defect,
        controls.maximum_angular_impulse_defect_kg_m2_s,
        "viscous angular impulse exceeds its independent budget",
    )?;
    let energy_defect = budget(
        energy_defect,
        controls.maximum_energy_defect_j,
        "molecular work, numerical losses or representation exceed their independent budget",
    )?;
    let report = MaterialViscosityReport {
        sweeps: solver.sweeps,
        krylov_iterations: solver.momentum_krylov,
        relative_residual: relative,
        face_impulse_kg_m_s: face_impulse.map(Dd::to_f64),
        cell_momentum_change_kg_m_s: cell_change.map(Dd::to_f64),
        absolute_impulse_defect_kg_m_s: impulse_defect,
        face_angular_impulse_kg_m2_s: face_angular.map(Dd::to_f64),
        cell_angular_momentum_change_kg_m2_s: cell_angular.map(Dd::to_f64),
        absolute_angular_impulse_defect_kg_m2_s: angular_defect,
        cell_molecular_heat_j: cell_heat,
        molecular_heat_j: heat.to_f64(),
        wall_work_by_face_j: wall_by_face,
        wall_work_j: wall.to_f64(),
        native_kinetic_energy_change_j: native_change.to_f64(),
        cell_kinetic_energy_change_j: cell_energy.to_f64(),
        backward_euler_loss_j: time_loss.to_f64(),
        numerical_scatter_loss_j: scatter.to_f64(),
        representation_work_j: representation.to_f64(),
        energy_defect_j: defect.to_f64(),
        absolute_energy_defect_j: energy_defect,
    };
    for value in report
        .face_impulse_kg_m_s
        .into_iter()
        .chain(report.cell_momentum_change_kg_m_s)
        .chain(report.face_angular_impulse_kg_m2_s)
        .chain(report.cell_angular_momentum_change_kg_m2_s)
        .chain(report.wall_work_by_face_j)
        .chain([
            report.molecular_heat_j,
            report.wall_work_j,
            report.native_kinetic_energy_change_j,
            report.cell_kinetic_energy_change_j,
            report.backward_euler_loss_j,
            report.numerical_scatter_loss_j,
            report.representation_work_j,
            report.energy_defect_j,
        ])
    {
        finite("flow.viscosity_report", value)?;
    }
    poll(gate)?;
    Ok(NativeMaterialViscosityStep {
        momentum_kg_m_s: momentum,
        report,
    })
}
