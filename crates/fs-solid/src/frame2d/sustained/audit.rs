//! Independent measured-state, force, mode and work gates for native trials.
use super::*;
use crate::frame2d::{DampingPreparation, mass::MassMatrix};

const ROUND: f64 = 512.0 * f64::EPSILON;

pub(super) struct Measured {
    pub displacement_change: Vec<f64>,
    pub velocity_change: Vec<f64>,
    pub generalized_momentum_change: Vec<f64>,
    pub support_impulse: Vec<f64>,
    pub kinetic_energy_change_j: f64,
    pub strain_energy_change_j: f64,
    pub material_damping_energy_j: f64,
    pub rayleigh_damping_energy_j: f64,
    pub external_work_j: f64,
    pub obstacle_work_j: f64,
    pub normal_contact_work_j: f64,
    pub frictional_work_j: f64,
    pub energy_balance_residual_j: f64,
    pub native_energy_balance_residual_j: f64,
}
fn gate<E>(
    quantity: &'static str,
    residual: Dd,
    allowance: Dd,
) -> Result<(), FrameSustainedError<E>> {
    let residual = finite(residual)?;
    let allowed = finite(allowance)?;
    if allowed < 0.0 || residual.abs() > allowed {
        Err(FrameSustainedError::Balance {
            quantity,
            residual,
            allowed,
        })
    } else {
        Ok(())
    }
}
fn rounding(terms: &[Dd]) -> Dd {
    terms
        .iter()
        .fold(Dd::ZERO, |sum, term| sum + term.abs() * Dd::from_f64(ROUND))
}
fn action(mass: &MassMatrix, row: usize, x: &[f64]) -> Dd {
    x.iter().enumerate().fold(Dd::ZERO, |sum, (j, &x)| {
        sum + Dd::from_f64(mass.entry(row, j)) * Dd::from_f64(x)
    })
}
fn dot_midpoint(a: &[f64], b: &[Dd]) -> Dd {
    a.iter()
        .zip(b)
        .fold(Dd::ZERO, |sum, (&a, &b)| sum + Dd::from_f64(a) * b)
}
fn negative(value: Dd) -> bool {
    value.hi < 0.0 || (value.hi == 0.0 && value.lo < 0.0)
}
fn mode<E>(id: u64, reason: &'static str) -> FrameSustainedError<E> {
    FrameSustainedError::ContactMode {
        contact_id: id,
        reason,
    }
}

pub(super) fn check<E>(
    before: &Frame2d,
    trial: &Trial,
    contacts: &[FrameSustainedContact<'_>],
    external: &[f64],
    controls: FrameSustainedControls,
    mass: &MassMatrix,
    damping: &DampingPreparation,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<Measured, FrameSustainedError<E>> {
    let after = &trial.frame;
    let n = before.ndof();
    let h = Dd::from_f64(controls.duration_s);
    let energy_allowance = Dd::from_f64(controls.energy_tolerance_j);
    if after
        .u
        .iter()
        .chain(&after.vel)
        .chain(&after.acc)
        .any(|v| !v.is_finite())
        || trial.history.u.len() != 2
        || trial.history.t.len() != 2
        || [
            &trial.history.kinetic_energy,
            &trial.history.internal_work,
            &trial.history.input_energy,
            &trial.history.damping_energy,
            &trial.history.material_damping_energy,
            &trial.history.rayleigh_damping_energy,
        ]
        .iter()
        .any(|values| values.len() != 2 || values.iter().any(|v| !v.is_finite()))
    {
        return Err(FrameSustainedError::Invalid(
            "nonfinite or incomplete native interval",
        ));
    }
    let du: Vec<_> = after.u.iter().zip(&before.u).map(|(a, b)| a - b).collect();
    let dv: Vec<_> = after
        .vel
        .iter()
        .zip(&before.vel)
        .map(|(a, b)| a - b)
        .collect();
    let vm: Vec<_> = before
        .vel
        .iter()
        .zip(&after.vel)
        .map(|(&a, &b)| Dd::from_f64(0.5 * a) + Dd::from_f64(0.5 * b))
        .collect();
    for i in 0..n {
        // Physical metre/radian allowances have no absolute world-position
        // floor: resolved motion lost on large u cannot become contact work.
        let left = h * Dd::from_f64(0.5 * before.vel[i]);
        let right = h * Dd::from_f64(0.5 * after.vel[i]);
        gate(
            "actual Newmark displacement",
            Dd::from_f64(du[i]) - left - right,
            Dd::from_f64(if i % 3 == 2 {
                controls.rotation_tolerance_rad
            } else {
                controls.displacement_tolerance_m
            }) + rounding(&[Dd::from_f64(du[i]), left, right]),
        )?;
        if before.fixed[i] && (du[i] != 0.0 || after.vel[i] != 0.0 || after.acc[i] != 0.0) {
            return Err(FrameSustainedError::Invalid(
                "stationary support moved during native interval",
            ));
        }
    }
    poll(cancelled)?;
    let (f0, _, _) = before.assemble(&before.u)?;
    let (f1, _, _) = after.assemble(&after.u)?;
    poll(cancelled)?;
    let kinetic_before = mass.kinetic_energy(&before.vel);
    let kinetic_after = mass.kinetic_energy(&after.vel);
    if !kinetic_before.is_finite()
        || !kinetic_after.is_finite()
        || kinetic_before < 0.0
        || kinetic_after < 0.0
    {
        return Err(FrameSustainedError::Invalid(
            "native kinetic energy is unrepresentable",
        ));
    }
    let mut momentum = Vec::with_capacity(n);
    let mut support = vec![0.0; n];
    let mut kinetic = Dd::ZERO;
    let mut elastic = Dd::ZERO;
    let mut rayleigh = Dd::ZERO;
    let mut material = Dd::ZERO;
    let mut external_work = Dd::ZERO;
    let mut applied_work = Dd::ZERO;
    let mut work_roundoff = Dd::ZERO;
    for i in 0..n {
        poll(cancelled)?;
        let change = action(mass, i, &dv);
        momentum.push(finite(change)?);
        let internal = Dd::from_f64(0.5 * f0[i]) + Dd::from_f64(0.5 * f1[i]);
        let viscous = dot_midpoint(&damping.cmat[i * n..(i + 1) * n], &vm);
        let load = Dd::from_f64(trial.load[i]);
        let expected = h * (load - internal - viscous);
        let tol = if i % 3 == 2 {
            controls.moment_tolerance_nm
        } else {
            controls.force_tolerance_n
        };
        if before.fixed[i] {
            support[i] = finite(change - expected)?;
        } else {
            gate(
                "original generalized momentum",
                change - expected,
                h * Dd::from_f64(tol) + rounding(&[change, h * load, h * internal, h * viscous]),
            )?;
            let inertia = action(mass, i, &after.acc);
            let endpoint_viscous = dot(&damping.cmat[i * n..(i + 1) * n], &after.vel);
            gate(
                "native finite endpoint force",
                inertia - load + Dd::from_f64(f1[i]) + endpoint_viscous,
                Dd::from_f64(tol)
                    + rounding(&[inertia, load, Dd::from_f64(f1[i]), endpoint_viscous]),
            )?;
        }
        let displacement = Dd::from_f64(du[i]);
        let e_kin = vm[i] * change;
        let e_elastic = internal * displacement;
        let e_rayleigh =
            dot_midpoint(&damping.rayleigh_matrix[i * n..(i + 1) * n], &vm) * displacement;
        let e_material = damping.material_damping.as_ref().map_or(Dd::ZERO, |op| {
            dot_midpoint(&op.entries()[i * n..(i + 1) * n], &vm) * displacement
        });
        let e_external = Dd::from_f64(external[i]) * displacement;
        let e_applied = load * displacement;
        kinetic = kinetic + e_kin;
        elastic = elastic + e_elastic;
        rayleigh = rayleigh + e_rayleigh;
        material = material + e_material;
        external_work = external_work + e_external;
        applied_work = applied_work + e_applied;
        work_roundoff = work_roundoff
            + rounding(&[
                e_kin, e_elastic, e_rayleigh, e_material, e_external, e_applied,
            ]);
    }
    let rayleigh_value = finite(rayleigh)?;
    let material_value = finite(material)?;
    if rayleigh_value < -controls.energy_tolerance_j
        || material_value < -controls.energy_tolerance_j
    {
        return Err(FrameSustainedError::Invalid(
            "native damping produced negative work beyond the allowance",
        ));
    }
    let mut obstacle = Dd::ZERO;
    let mut normal_work = Dd::ZERO;
    let mut friction_work = Dd::ZERO;
    let mut contact_world = Dd::ZERO;
    let mut friction_roundoff = Dd::ZERO;
    for (i, c) in contacts.iter().enumerate() {
        poll(cancelled)?;
        let force = trial.forces[i];
        let pn = Dd::from_f64(force[0]);
        let pt = Dd::from_f64(force[1]);
        if force[0] < 0.0 || force.iter().any(|v| !v.is_finite()) {
            return Err(mode(c.id, "nonnegative finite normal force required"));
        }
        let relative_n = dot(c.normal_covector, &du) - h * Dd::from_f64(c.obstacle_velocity[0]);
        let slip = dot(c.tangent_covector, &du) - h * Dd::from_f64(c.obstacle_velocity[1]);
        let gap = Dd::from_f64(c.gap_m) + relative_n;
        let gap_value = finite(gap)?;
        let slip_value = finite(slip)?;
        if negative(gap) {
            return Err(mode(
                c.id,
                "negative actual endpoint gap; position repair is forbidden",
            ));
        }
        let capacity = Dd::from_f64(c.friction_coefficient) * pn;
        let capacity_value = finite(capacity)?;
        let cone = finite(pt.abs() - capacity)?;
        if cone > controls.force_tolerance_n + finite(rounding(&[pt, capacity]))? {
            return Err(mode(c.id, "force exceeds the original Coulomb cone"));
        }
        gate(
            "normal endpoint complementarity",
            pn * gap,
            Dd::from_f64(controls.complementarity_tolerance_j) + rounding(&[pn * gap]),
        )?;
        gate(
            "tangent displacement duality",
            capacity * slip.abs() + pt * slip,
            Dd::from_f64(controls.complementarity_tolerance_j)
                + rounding(&[capacity * slip.abs(), pt * slip]),
        )?;
        let vn0 =
            finite(dot(c.normal_covector, &before.vel) - Dd::from_f64(c.obstacle_velocity[0]))?;
        let vn1 =
            finite(dot(c.normal_covector, &after.vel) - Dd::from_f64(c.obstacle_velocity[0]))?;
        let vt0 =
            finite(dot(c.tangent_covector, &before.vel) - Dd::from_f64(c.obstacle_velocity[1]))?;
        let vt1 =
            finite(dot(c.tangent_covector, &after.vel) - Dd::from_f64(c.obstacle_velocity[1]))?;
        if force[0] > controls.force_tolerance_n {
            if c.gap_m > controls.displacement_tolerance_m
                || gap_value > controls.displacement_tolerance_m
                || vn0.abs() > controls.velocity_tolerance_m_s
                || vn1.abs() > controls.velocity_tolerance_m_s
            {
                return Err(mode(
                    c.id,
                    "loaded point does not dwell; arrival/detachment needs separate event resolution",
                ));
            }
            gate(
                "normal initial complementarity",
                pn * Dd::from_f64(c.gap_m),
                Dd::from_f64(controls.complementarity_tolerance_j)
                    + rounding(&[pn * Dd::from_f64(c.gap_m)]),
            )?;
            if capacity_value > controls.force_tolerance_n {
                let inside = force[1].abs() < capacity_value - controls.force_tolerance_n;
                if inside
                    && (slip_value.abs() > controls.displacement_tolerance_m
                        || vt0.abs() > controls.velocity_tolerance_m_s
                        || vt1.abs() > controls.velocity_tolerance_m_s)
                {
                    return Err(mode(
                        c.id,
                        "interior-cone force requires sticking throughout the admitted interval",
                    ));
                }
                // Opposing endpoint signs prove a reversal even when net slip
                // retains its old sign. Same-sign endpoints alone never prove
                // absence of an interior reversal: that is the path owner's job.
                if (vt0 > controls.velocity_tolerance_m_s && vt1 < -controls.velocity_tolerance_m_s)
                    || (vt0 < -controls.velocity_tolerance_m_s
                        && vt1 > controls.velocity_tolerance_m_s)
                    || (force[1] > controls.force_tolerance_n
                        && (vt0 > controls.velocity_tolerance_m_s
                            || vt1 > controls.velocity_tolerance_m_s))
                    || (force[1] < -controls.force_tolerance_n
                        && (vt0 < -controls.velocity_tolerance_m_s
                            || vt1 < -controls.velocity_tolerance_m_s))
                {
                    return Err(mode(
                        c.id,
                        "held friction opposes the wrong endpoint velocity or crosses a stop",
                    ));
                }
            }
        }
        let e_obstacle = h
            * (pn * Dd::from_f64(c.obstacle_velocity[0])
                + pt * Dd::from_f64(c.obstacle_velocity[1]));
        let e_normal = pn * relative_n;
        let e_friction = pt * slip;
        let e_world = pn * dot(c.normal_covector, &du) + pt * dot(c.tangent_covector, &du);
        obstacle = obstacle + e_obstacle;
        normal_work = normal_work + e_normal;
        friction_work = friction_work + e_friction;
        friction_roundoff = friction_roundoff + rounding(&[e_friction]);
        contact_world = contact_world + e_world;
        work_roundoff = work_roundoff + rounding(&[e_obstacle, e_normal, e_friction, e_world]);
    }
    let friction_value = finite(friction_work)?;
    if friction_value > finite(energy_allowance + friction_roundoff)? {
        return Err(FrameSustainedError::Balance {
            quantity: "positive relative friction work",
            residual: friction_value,
            allowed: finite(energy_allowance + friction_roundoff)?,
        });
    }
    let work_allowance = energy_allowance + work_roundoff;
    gate(
        "native applied load versus conjugate contact work",
        applied_work - external_work - contact_world,
        work_allowance,
    )?;
    gate(
        "contact source work partition",
        contact_world - obstacle - normal_work - friction_work,
        work_allowance,
    )?;
    let balance = kinetic + elastic + rayleigh + material
        - external_work
        - obstacle
        - normal_work
        - friction_work;
    gate("whole native/contact energy", balance, work_allowance)?;
    let native_input = Dd::from_f64(trial.history.input_energy[1]);
    let native_internal = Dd::from_f64(trial.history.internal_work[1]);
    let native_damping = Dd::from_f64(trial.history.damping_energy[1]);
    gate(
        "actual versus native input work",
        applied_work - native_input,
        work_allowance + rounding(&[native_input]),
    )?;
    gate(
        "actual versus native elastic work",
        elastic - native_internal,
        work_allowance + rounding(&[native_internal]),
    )?;
    gate(
        "actual versus native material work",
        material - Dd::from_f64(trial.history.material_damping_energy[1]),
        work_allowance + rounding(&[material]),
    )?;
    gate(
        "actual versus native Rayleigh work",
        rayleigh - Dd::from_f64(trial.history.rayleigh_damping_energy[1]),
        work_allowance + rounding(&[rayleigh]),
    )?;
    let native_balance = kinetic + native_internal + native_damping - native_input;
    gate(
        "original native energy ledger",
        native_balance,
        work_allowance + rounding(&[native_input, native_internal, native_damping]),
    )?;
    Ok(Measured {
        displacement_change: du,
        velocity_change: dv,
        generalized_momentum_change: momentum,
        support_impulse: support,
        kinetic_energy_change_j: finite(kinetic)?,
        strain_energy_change_j: finite(elastic)?,
        material_damping_energy_j: material_value,
        rayleigh_damping_energy_j: rayleigh_value,
        external_work_j: finite(external_work)?,
        obstacle_work_j: finite(obstacle)?,
        normal_contact_work_j: finite(normal_work)?,
        frictional_work_j: friction_value,
        energy_balance_residual_j: finite(balance)?,
        native_energy_balance_residual_j: finite(native_balance)?,
    })
}
