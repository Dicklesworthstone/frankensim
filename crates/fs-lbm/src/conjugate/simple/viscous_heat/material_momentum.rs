//! Molecular momentum diffusion on the original MAC stress and Krylov owners.
//!
//! Pressure, convection and calorics remain separate physical substeps. This
//! predictor uses actual M/V face inertia and either the original 2 mu D stress
//! or the explicit independent-bulk response 2 mu dev(D) + zeta div(u) I.
//! Both retain harmonic material jumps, half-cell walls and matching heat.
//! The source opt-in resolves viscosity at the accepted native endpoint shear.
//! The SGS opt-in adds a deviatoric modeled stress using actual cell density
//! and reports its energy transfer separately from molecular heating. Reservoir
//! mouths retain normal traction and a zero tangential velocity trace; their
//! boundary closure and half-edge power use the same native stress coefficient.

mod rheology;
pub(in crate::conjugate) use rheology::{
    NativeMaterialRheologyInput, NativeMaterialRheologyStep, kick_with_source,
};

use super::*;
use fs_material::fluid_transport::viscous::NewtonianViscosity;
use fs_math::dd::Dd;

/// Work and independent physical-unit gates for one native viscosity substep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialViscosityControls {
    pub max_cells: usize,
    /// Each sweep uses the original three momentum solves, each bounded by
    /// its native 20,000 Krylov iterations, and one same-velocity residual
    /// readback. The endpoint-rheology entry refreshes its exact source
    /// coefficients before that readback; this bounds the entire native solve.
    pub max_sweeps: usize,
    pub relative_tolerance: f64,
    pub maximum_impulse_defect_kg_m_s: f64,
    pub maximum_angular_impulse_defect_kg_m2_s: f64,
    pub maximum_energy_defect_j: f64,
}

/// Physical molecular heat, modeled SGS transfer and distinct numerical work.
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
    /// Modeled transfer out of resolved kinetic energy, owned by each original
    /// cell. This is NOT molecular heat or an evolved subgrid-energy store.
    /// Empty on the original molecular path, including exact zero SGS.
    pub cell_sgs_transfer_j: Vec<f64>,
    pub sgs_transfer_j: f64,
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
    /// Cell delta_K + molecular_heat + sgs_transfer + BE_loss + scatter_loss
    /// - wall_work - representation_work.
    pub energy_defect_j: f64,
    pub absolute_energy_defect_j: f64,
}

pub(in crate::conjugate) struct NativeMaterialViscosityInput<'a> {
    pub field: &'a FlowField,
    pub mass_kg: &'a [f64],
    pub momentum_kg_m_s: &'a [[f64; 3]],
    pub viscosity_pa_s: &'a [f64],
    /// None retains 2 mu D; Some selects independent physical bulk response.
    pub bulk_viscosity_pa_s: Option<&'a [f64]>,
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

fn admit_split_boundaries(
    boundaries: [FvBoundary; 6],
    pressure_mouth_axis: Option<usize>,
) -> Result<(), ChtError> {
    if pressure_mouth_axis.is_some_and(|a| a >= 3)
        || boundaries.iter().enumerate().any(|(f, b)| {
            if pressure_mouth_axis == Some(f / 2) {
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
    Ok(())
}

/// Observe the physical shear invariant on an original reconstructed field.
/// The unit molecular operator retains the existing symmetric-stress cell,
/// shear-edge, half-cell wall and reservoir-mouth quadratures. No momentum
/// solve, inertia, material coefficient or thermal source enters this sample.
pub(in crate::conjugate) fn physical_shear_rates(
    domain: &VoxelDomain,
    field: &FlowField,
    boundaries: [FvBoundary; 6],
    pressure_mouth_axis: Option<usize>,
    max_cells: usize,
    gate: &CancelGate,
) -> Result<Vec<f64>, ChtError> {
    poll(gate)?;
    if domain.cell_count() > max_cells || field.dims() != domain.dims() {
        return Err(refuse("physical shear needs the bounded original grid"));
    }
    admit_split_boundaries(boundaries, pressure_mouth_axis)?;
    let config = SimpleConfig::new(boundaries);
    admit_molecular_boundaries(domain, &config)?;
    poll(gate)?;
    let mut sampler = Solver::with_molecular_coefficients(domain, 1.0, 1.0, &config);
    sampler.molecular_viscosity = vec![1.0; domain.cell_count()];
    sampler.pressure_work_axis = pressure_mouth_axis;
    for (axis, fluxes) in field.face_arrays().into_iter().enumerate() {
        if fluxes.len() != sampler.vel[axis].len() {
            return Err(refuse(
                "physical shear face arrays differ from the original grid",
            ));
        }
        for (slot, &flux) in fluxes.iter().enumerate() {
            if slot % 256 == 0 {
                poll(gate)?;
            }
            let velocity = flux / sampler.comp_face_area(axis, slot);
            finite("flow.rheology_physical_velocity", velocity)?;
            if let Kind::Fixed(required) = sampler.comps[axis].kind[slot] {
                if velocity != required {
                    return Err(refuse(
                        "physical shear cannot change a held normal velocity",
                    ));
                }
            }
            sampler.vel[axis][slot] = velocity;
        }
    }
    sampler.molecular_shear_rates(gate)
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

pub(in crate::conjugate) fn kick(
    domain: &VoxelDomain,
    input: NativeMaterialViscosityInput<'_>,
    controls: MaterialViscosityControls,
    gate: &CancelGate,
) -> Result<NativeMaterialViscosityStep, ChtError> {
    kick_with_refresh(domain, input, controls, gate, |_| Ok(()))
}

/// Add resolved-gradient SGS on the original material inertia and stress owner.
/// Closed walls/symmetry and one opposed reservoir-traction pair retain their
/// original boundary kinematics. Source molecular heat and modeled SGS transfer
/// stay separate, including the two-cell half-edges at reservoir mouths.
pub(in crate::conjugate) fn kick_with_sgs(
    domain: &VoxelDomain,
    input: NativeMaterialViscosityInput<'_>,
    sgs: crate::conjugate::FvSgsConfig,
    controls: MaterialViscosityControls,
    gate: &CancelGate,
) -> Result<NativeMaterialViscosityStep, ChtError> {
    poll(gate)?;
    sgs.validate()?;
    if sgs.model.coefficient() == 0.0 {
        return kick(domain, input, controls, gate);
    }
    admit_split_boundaries(input.boundaries, input.pressure_mouth_axis)?;
    let mass = input.mass_kg;
    let mut initialized = false;
    kick_with_refresh(domain, input, controls, gate, |solver| {
        if !initialized {
            solver.set_sgs_density(mass, gate)?;
            if solver.molecular_bulk_viscosity.is_empty() {
                // Preserve the ORIGINAL molecular 2 mu D law while adding
                // only 2 rho nu_t dev(D). Its existing isotropic part is
                // exactly zeta = 2 mu / 3, not an inferred new material law.
                let n = solver.domain.cell_count();
                solver.molecular_bulk_viscosity.resize(n, 0.0);
                for cell in 0..n {
                    if cell % 256 == 0 {
                        poll(gate)?;
                    }
                    if solver.domain.is_fluid(cell) {
                        let bulk = (2.0 / 3.0) * solver.molecular_viscosity[cell];
                        finite_positive("flow.sgs_original_isotropic_viscosity", bulk)?;
                        solver.molecular_bulk_viscosity[cell] = bulk;
                    }
                }
            }
            // Explicit original bulk coefficients have already been admitted
            // unchanged by the shared native setup. SGS adds no bulk loss.
            solver.enable_sgs(sgs, gate)?;
            initialized = true;
        }
        // Called before the first solve and after every coupled sweep, BEFORE
        // all three residuals. No stale closure enters an accepted power audit.
        solver.update_sgs(gate)
    })
}

/// Retain one original native time level through every coupled sweep.
/// The frozen owner supplies a no-op; a constitutive owner refreshes exact
/// endpoint coefficients before all three residual rows and the same audit.
#[allow(clippy::too_many_lines)]
fn kick_with_refresh(
    domain: &VoxelDomain,
    input: NativeMaterialViscosityInput<'_>,
    controls: MaterialViscosityControls,
    gate: &CancelGate,
    mut refresh: impl FnMut(&mut Solver<'_>) -> Result<(), ChtError>,
) -> Result<NativeMaterialViscosityStep, ChtError> {
    poll(gate)?;
    let n = domain.cell_count();
    if n > controls.max_cells
        || controls.max_sweeps == 0
        || input.mass_kg.len() != n
        || input.momentum_kg_m_s.len() != n
        || input.viscosity_pa_s.len() != n
        || input
            .bulk_viscosity_pa_s
            .is_some_and(|bulk| bulk.len() != n)
    {
        return Err(refuse("viscosity requires bounded original cells"));
    }
    admit_split_boundaries(input.boundaries, input.pressure_mouth_axis)?;
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
    if let Some(bulk) = input.bulk_viscosity_pa_s {
        set_bulk_viscosity(&mut solver, input.viscosity_pa_s, bulk, gate)?;
    }
    solver.pressure_work_axis = input.pressure_mouth_axis;
    // This is only stress diffusion. The enclosing owner subsequently
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
    refresh(&mut solver)?;
    let mut relative = f64::INFINITY;
    let mut residual_impulse = f64::INFINITY;
    for sweep in 1..=controls.max_sweeps {
        poll(gate)?;
        solver.sweeps = sweep;
        for a in 0..3 {
            solver.momentum(a, true, gate)?;
        }
        solver.guard_finite()?;
        refresh(&mut solver)?;
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

/// Refresh independent bulk response while retaining the original time level
/// and the complete existing shear/volumetric operator. Solid slots are inert.
fn set_bulk_viscosity(
    solver: &mut Solver<'_>,
    shear: &[f64],
    bulk: &[f64],
    gate: &CancelGate,
) -> Result<(), ChtError> {
    let n = solver.domain.cell_count();
    if shear.len() != n || bulk.len() != n {
        return Err(refuse(
            "shear and bulk viscosity must match the original grid",
        ));
    }
    solver.molecular_bulk_viscosity.resize(n, 0.0);
    for (cell, &value) in bulk.iter().enumerate() {
        if cell % 256 == 0 {
            poll(gate)?;
        }
        if !solver.domain.is_fluid(cell) {
            solver.molecular_bulk_viscosity[cell] = 0.0;
            continue;
        }
        let law = NewtonianViscosity::new(shear[cell], value).map_err(|error| {
            ChtError::InvalidInput {
                field: "flow.bulk_viscosity_source",
                reason: format!("cell {cell}: {error}"),
            }
        })?;
        finite_positive(
            "flow.longitudinal_viscosity",
            (4.0 / 3.0) * law.shear_viscosity_pa_s() + law.bulk_viscosity_pa_s(),
        )?;
        solver.molecular_bulk_viscosity[cell] = law.bulk_viscosity_pa_s();
    }
    Ok(())
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
    let power = if solver.has_sgs() {
        solver.sample_material_sgs_power(gate)?
    } else {
        solver.sample_molecular_power::<false>(gate)?
    };
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
    let mut cell_sgs_transfer = Vec::with_capacity(power.sgs_transfer.len());
    let mut sgs_transfer = Dd::ZERO;
    for (cell, watts) in power.sgs_transfer.iter().copied().enumerate() {
        if cell % 256 == 0 {
            poll(gate)?;
        }
        let joules = (duration * Dd::from_f64(watts)).to_f64();
        finite("flow.sgs_energy_transfer", joules)?;
        if joules < 0.0 || (!domain.is_fluid(cell) && joules != 0.0) {
            return Err(refuse(
                "SGS transfer must be nonnegative and mobile-cell owned",
            ));
        }
        cell_sgs_transfer.push(joules);
        sgs_transfer = sgs_transfer + Dd::from_f64(joules);
    }
    // Preserve the original molecular arithmetic exactly on the empty branch.
    let dissipation = if cell_sgs_transfer.is_empty() {
        heat
    } else {
        heat + sgs_transfer
    };
    let wall_by_face = power
        .wall_power_w
        .map(|w| (duration * Dd::from_f64(w)).to_f64());
    let wall = wall_by_face
        .into_iter()
        .fold(Dd::ZERO, |s, w| s + Dd::from_f64(w));
    let mut partition = magnitude(dissipation - duration * Dd::from_f64(power.quadrature_w))
        + magnitude(heat - duration * Dd::from_f64(power.molecular_heat_w));
    if !cell_sgs_transfer.is_empty() {
        partition =
            partition + magnitude(sgs_transfer - duration * Dd::from_f64(power.sgs_transfer_w));
    }
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
    let native_defect = native_change + time_loss + dissipation - wall;
    let feedback_defect = cell_energy + scatter - native_change - representation;
    let defect = cell_energy + time_loss + scatter + dissipation - wall - representation;
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
        if solver.has_sgs() {
            "molecular/SGS work, numerical losses or representation exceed their independent budget"
        } else {
            "molecular work, numerical losses or representation exceed their independent budget"
        },
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
        cell_sgs_transfer_j: cell_sgs_transfer,
        sgs_transfer_j: sgs_transfer.to_f64(),
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
            report.sgs_transfer_j,
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
