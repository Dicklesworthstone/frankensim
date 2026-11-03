//! Independent acceptance of actual retained-body sustained contact paths.
//!
//! Inputs are sampled from the original native successors and native work
//! receipts. This audit supplies no motion, contact force, state projection or
//! heat allocation. The supported path begins in closed, normal-stationary
//! contact; release is permitted, while an unresolved sliding reversal refuses.

/// Actual contact coordinates in one fixed orthonormal normal-first frame.
#[derive(Clone, Copy, Debug)]
pub struct ContactAuditInput {
    pub id: u64,
    pub gap_before_m: f64,
    pub gap_after_m: f64,
    /// Full relative displacement increment, including prescribed motion.
    pub slip_increment_m: [f64; 3],
    pub relative_velocity_before_m_s: [f64; 3],
    pub relative_velocity_after_m_s: [f64; 3],
    /// Held force on the retained bodies, conjugate to the positive trace.
    pub force_n: [f64; 3],
    pub friction_coefficient: f64,
    /// c in relative velocity = B v + c, constant during this interval.
    pub prescribed_velocity_m_s: [f64; 3],
    pub activation_tolerance_m: f64,
}

/// Original body increments. Native total external work includes the supplied
/// contact forces; authored work is evaluated independently from caller loads.
#[derive(Clone, Copy, Debug)]
pub struct BodyWork {
    pub kinetic_change_j: f64,
    pub bulk_change_j: f64,
    pub cohesive_change_j: f64,
    /// Original fracture, viscosity and internal-friction losses, once each.
    pub material_loss_j: f64,
    pub authored_work_j: f64,
    pub prescribed_jump_work_j: f64,
    pub native_total_external_work_j: f64,
    pub native_energy_defect_j: f64,
}

/// Physical acceptance tolerances; none changes a supplied force or position.
#[derive(Clone, Copy, Debug)]
pub struct AuditControls {
    pub force_tolerance_n: f64,
    pub gap_tolerance_m: f64,
    pub slip_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct AuditReport {
    pub friction_loss_j: f64,
    /// Work delivered through the prescribed surface: -h sum(c dot lambda).
    pub surface_work_j: f64,
    /// Signed normal reaction work. This is a constraint defect, not a loss.
    pub normal_constraint_work_j: f64,
    pub absolute_normal_constraint_work_j: f64,
    pub relative_contact_work_j: f64,
    /// Actual native external work less independently authored body work.
    pub contact_body_work_j: f64,
    /// Actual body contact work minus relative work and prescribed work.
    pub contact_work_residual_j: f64,
    /// Authored + prescribed work minus storage and all physical losses.
    /// Normal constraint work remains visible; it is never debited as heat.
    pub energy_residual_j: f64,
    /// Sum of the larger absolute original/reconstructed defect of each body.
    pub absolute_native_energy_defect_j: f64,
    pub absolute_friction_work_defect_j: f64,
    pub maximum_kinematic_residual_m: f64,
    pub maximum_complementarity_work_j: f64,
    pub maximum_force_residual_n: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AuditError {
    Invalid(&'static str),
    NonFinite(&'static str),
    Cancelled,
    Contact {
        id: u64,
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
    Energy {
        quantity: &'static str,
        residual_j: f64,
        tolerance_j: f64,
    },
    UnresolvedReversal {
        contact_id: u64,
    },
}
impl core::fmt::Display for AuditError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive sustained-contact audit refused: {self:?}")
    }
}
impl std::error::Error for AuditError {}

/// Audit a complete private candidate before any native successor publishes.
/// The original owner bounds body work; this allocation-free audit polls each
/// input record and accepts at most the native 21-contact family.
pub(super) fn audit(
    contacts: &[ContactAuditInput],
    bodies: &[BodyWork],
    dt_s: f64,
    controls: AuditControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AuditReport, AuditError> {
    poll(cancelled)?;
    if bodies.is_empty()
        || contacts.len() > 21
        || !positive(dt_s)
        || !positive(controls.force_tolerance_n)
        || !nonnegative(controls.gap_tolerance_m)
        || !positive(controls.slip_tolerance_m)
        || !positive(controls.velocity_tolerance_m_s)
        || !positive(controls.energy_tolerance_j)
    {
        return Err(AuditError::Invalid(
            "finite bounded contact audit inputs required",
        ));
    }
    let mut physical = Sum::default();
    let mut contact_body = Sum::default();
    let mut contact_closure = Sum::default();
    let mut native_defect = Sum::default();
    for body in bodies {
        poll(cancelled)?;
        if !nonnegative(body.material_loss_j) {
            return Err(AuditError::Invalid(
                "nonnegative original material losses required",
            ));
        }
        let mut reconstructed = Sum::default();
        for value in [
            body.native_total_external_work_j,
            body.prescribed_jump_work_j,
            -body.kinetic_change_j,
            -body.bulk_change_j,
            -body.cohesive_change_j,
            -body.material_loss_j,
        ] {
            reconstructed.add(value)?;
        }
        let original = finite(body.native_energy_defect_j, "original body energy defect")?;
        native_defect.add(reconstructed.total()?.abs().max(original.abs()))?;
        for value in [
            body.authored_work_j,
            body.prescribed_jump_work_j,
            -body.kinetic_change_j,
            -body.bulk_change_j,
            -body.cohesive_change_j,
            -body.material_loss_j,
        ] {
            physical.add(value)?;
        }
        for value in [body.native_total_external_work_j, -body.authored_work_j] {
            contact_body.add(value)?;
            contact_closure.add(value)?;
        }
    }

    let mut friction_loss = Sum::default();
    let mut surface_work = Sum::default();
    let mut normal_work = Sum::default();
    let mut absolute_normal_work = Sum::default();
    let mut relative_work = Sum::default();
    let mut friction_defect = Sum::default();
    let mut maximum_kinematic_residual_m = 0.0_f64;
    let mut maximum_complementarity_work_j = 0.0_f64;
    let mut maximum_force_residual_n = 0.0_f64;

    for (index, contact) in contacts.iter().enumerate() {
        poll(cancelled)?;
        if contacts[..index].iter().any(|old| old.id == contact.id) {
            return Err(AuditError::Invalid("unique contact audit IDs required"));
        }
        validate_contact(contact)?;
        let id = contact.id;
        gate(id, "initial penetration [m]", -contact.gap_before_m, 0.0)?;
        gate(id, "endpoint penetration [m]", -contact.gap_after_m, 0.0)?;
        gate(id, "tensile normal reaction [N]", -contact.force_n[0], 0.0)?;
        gate(
            id,
            "initial contact activation [m]",
            contact.gap_before_m - contact.activation_tolerance_m,
            0.0,
        )?;
        gate(
            id,
            "initial normal velocity [m/s]",
            contact.relative_velocity_before_m_s[0].abs(),
            controls.velocity_tolerance_m_s,
        )?;

        let gap_defect = finite(
            finite(
                contact.gap_after_m - contact.gap_before_m,
                "normal gap increment",
            )? - contact.slip_increment_m[0],
            "normal gap/trace consistency",
        )?
        .abs();
        maximum_kinematic_residual_m = maximum_kinematic_residual_m.max(gap_defect);
        gate(
            id,
            "normal gap/trace consistency [m]",
            gap_defect,
            controls.gap_tolerance_m,
        )?;
        for axis in 0..3 {
            let expected = finite(
                dt_s * f64::midpoint(
                    contact.relative_velocity_before_m_s[axis],
                    contact.relative_velocity_after_m_s[axis],
                ),
                "native midpoint trace increment",
            )?;
            let defect = finite(
                contact.slip_increment_m[axis] - expected,
                "native midpoint trace consistency",
            )?
            .abs();
            maximum_kinematic_residual_m = maximum_kinematic_residual_m.max(defect);
            let tolerance = if axis == 0 {
                controls.gap_tolerance_m
            } else {
                controls.slip_tolerance_m
            };
            gate(id, "actual midpoint kinematics [m]", defect, tolerance)?;
        }
        gate(
            id,
            "inward endpoint normal velocity [m/s]",
            -contact.relative_velocity_after_m_s[0],
            controls.velocity_tolerance_m_s,
        )?;
        if contact.force_n[0] > controls.force_tolerance_n {
            gate(
                id,
                "loaded endpoint normal velocity [m/s]",
                contact.relative_velocity_after_m_s[0].abs(),
                controls.velocity_tolerance_m_s,
            )?;
        }
        let complementarity = finite(
            contact.force_n[0] * contact.gap_after_m,
            "normal complementarity work",
        )?;
        maximum_complementarity_work_j = maximum_complementarity_work_j.max(complementarity);
        gate(
            id,
            "normal complementarity [J]",
            complementarity,
            controls.energy_tolerance_j,
        )?;

        let normal = finite(
            contact.force_n[0] * contact.slip_increment_m[0],
            "normal constraint work",
        )?;
        normal_work.add(normal)?;
        absolute_normal_work.add(normal.abs())?;
        let mut tangent_work = Sum::default();
        for axis in 0..3 {
            let work = finite(
                contact.force_n[axis] * contact.slip_increment_m[axis],
                "relative contact work",
            )?;
            relative_work.add(work)?;
            contact_closure.add(-work)?;
            if axis != 0 {
                tangent_work.add(work)?;
            }
            let displacement = finite(
                dt_s * contact.prescribed_velocity_m_s[axis],
                "prescribed surface displacement",
            )?;
            let supplied = finite(
                -displacement * contact.force_n[axis],
                "prescribed surface work",
            )?;
            surface_work.add(supplied)?;
            contact_closure.add(-supplied)?;
            physical.add(supplied)?;
        }

        let force = finite(
            contact.force_n[1].hypot(contact.force_n[2]),
            "tangential force magnitude",
        )?;
        let cap = finite(
            contact.friction_coefficient * contact.force_n[0],
            "Coulomb force capacity",
        )?;
        let excess = finite(force - cap, "Coulomb disk excess")?.max(0.0);
        maximum_force_residual_n = maximum_force_residual_n.max(excess);
        gate(
            id,
            "Coulomb force disk [N]",
            excess,
            controls.force_tolerance_n,
        )?;
        let distance = finite(
            contact.slip_increment_m[1].hypot(contact.slip_increment_m[2]),
            "actual tangential slip distance",
        )?;
        let v0 = &contact.relative_velocity_before_m_s;
        let v1 = &contact.relative_velocity_after_m_s;
        let speed0 = finite(v0[1].hypot(v0[2]), "initial tangential speed")?;
        let speed1 = finite(v1[1].hypot(v1[2]), "final tangential speed")?;
        let velocity_tol = controls.velocity_tolerance_m_s;
        if cap > 0.0 && speed0 > velocity_tol && speed1 > velocity_tol {
            // Normalize before products so large valid speeds do not overflow.
            let direction =
                (v0[1] / speed0) * (v1[1] / speed1) + (v0[2] / speed0) * (v1[2] / speed1);
            let allowance = (velocity_tol / speed0) * (velocity_tol / speed1);
            if direction < -allowance {
                return Err(AuditError::UnresolvedReversal { contact_id: id });
            }
        }
        let static_path = distance <= controls.slip_tolerance_m
            && speed0 <= velocity_tol
            && speed1 <= velocity_tol;
        // A released or frictionless trace has no dry-friction direction to
        // continue. Its native turnaround remains valid; the zero-radius disk
        // and independent work gate still constrain any residual force.
        let loss = if static_path || cap == 0.0 {
            0.0
        } else {
            if distance == 0.0 {
                return Err(AuditError::UnresolvedReversal { contact_id: id });
            }
            for axis in 1..3 {
                let expected = finite(
                    -cap * (contact.slip_increment_m[axis] / distance),
                    "sliding Coulomb force",
                )?;
                let residual =
                    finite(contact.force_n[axis] - expected, "sliding force residual")?.abs();
                maximum_force_residual_n = maximum_force_residual_n.max(residual);
                gate(
                    id,
                    "actual sliding force direction [N]",
                    residual,
                    controls.force_tolerance_n,
                )?;
            }
            finite(cap * distance, "Coulomb friction loss")?
        };
        let mut law_work = tangent_work;
        law_work.add(loss)?;
        let defect = law_work.total()?.abs();
        gate(
            id,
            "actual Coulomb work [J]",
            defect,
            controls.energy_tolerance_j,
        )?;
        friction_defect.add(defect)?;
        friction_loss.add(loss)?;
        physical.add(-loss)?;
    }

    let report = AuditReport {
        friction_loss_j: friction_loss.total()?,
        surface_work_j: surface_work.total()?,
        normal_constraint_work_j: normal_work.total()?,
        absolute_normal_constraint_work_j: absolute_normal_work.total()?,
        relative_contact_work_j: relative_work.total()?,
        contact_body_work_j: contact_body.total()?,
        contact_work_residual_j: contact_closure.total()?,
        energy_residual_j: physical.total()?,
        absolute_native_energy_defect_j: native_defect.total()?,
        absolute_friction_work_defect_j: friction_defect.total()?,
        maximum_kinematic_residual_m,
        maximum_complementarity_work_j,
        maximum_force_residual_n,
    };
    for (quantity, residual_j) in [
        (
            "normal constraint work",
            report.absolute_normal_constraint_work_j,
        ),
        (
            "actual contact/body work",
            report.contact_work_residual_j.abs(),
        ),
        ("full sustained first law", report.energy_residual_j.abs()),
        (
            "absolute native body defects",
            report.absolute_native_energy_defect_j,
        ),
        (
            "absolute Coulomb work defects",
            report.absolute_friction_work_defect_j,
        ),
    ] {
        if residual_j > controls.energy_tolerance_j {
            return Err(AuditError::Energy {
                quantity,
                residual_j,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
    }
    poll(cancelled)?;
    Ok(report)
}

fn validate_contact(contact: &ContactAuditInput) -> Result<(), AuditError> {
    for value in [contact.gap_before_m, contact.gap_after_m] {
        finite(value, "actual contact gap")?;
    }
    if !nonnegative(contact.friction_coefficient) || !nonnegative(contact.activation_tolerance_m) {
        return Err(AuditError::Invalid(
            "finite nonnegative contact law data required",
        ));
    }
    for vector in [
        contact.slip_increment_m,
        contact.relative_velocity_before_m_s,
        contact.relative_velocity_after_m_s,
        contact.force_n,
        contact.prescribed_velocity_m_s,
    ] {
        for value in vector {
            finite(value, "actual contact vector")?;
        }
    }
    Ok(())
}

fn gate(id: u64, quantity: &'static str, residual: f64, tolerance: f64) -> Result<(), AuditError> {
    finite(residual, quantity)?;
    if residual > tolerance {
        Err(AuditError::Contact {
            id,
            quantity,
            residual,
            tolerance,
        })
    } else {
        Ok(())
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), AuditError> {
        finite(value, "physical work contribution")?;
        let next = finite(self.value + value, "physical work sum")?;
        let residual = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.correction = finite(self.correction + residual, "physical work correction")?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, AuditError> {
        finite(self.value + self.correction, "physical work result")
    }
}
fn finite(value: f64, quantity: &'static str) -> Result<f64, AuditError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AuditError::NonFinite(quantity))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AuditError> {
    if cancelled() {
        Err(AuditError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "audit/tests.rs"]
mod tests;
