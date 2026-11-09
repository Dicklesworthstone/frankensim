//! Cold-to-hot material endpoints in the original undamped Newmark owner.
//!
//! The initial internal force and consistent acceleration belong to the OLD
//! temperature field. The original Newmark acceleration solve evaluates the
//! NEW field at its endpoint. This distinction is essential: refreshing first
//! and then reinitializing at the hot force invents a different thermal history.
//! Average acceleration supplies no numerical damping in this linear-reference
//! lane; an energy defect is checked, never assigned to heat or dissipation.

use super::{Frame2d, checkpoint, invalid};
use crate::SolidError;
use crate::frame2d::{DynamicHistory, DynamicInitialState, Geometry, Rayleigh, mass::MassMatrix};
use fs_exec::Cx;
use fs_math::dd::Dd;

/// Native size bounds and independently dimensioned accepted-state allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameLaminateNewmarkControls {
    pub maximum_dofs: usize,
    pub maximum_members: usize,
    pub maximum_force_residual_n: f64,
    pub maximum_moment_residual_n_m: f64,
    /// Allowed loss when adding the resolved increment to stored positions.
    pub maximum_displacement_defect_m: f64,
    pub maximum_rotation_defect_rad: f64,
    pub maximum_energy_defect_j: f64,
}

/// One unpublished original frame endpoint and its independent native account.
pub struct FrameLaminateNewmarkStep {
    pub frame: Frame2d,
    pub history: DynamicHistory,
    /// Resolved native Newmark increment before addition to absolute positions.
    pub displacement_increment: Vec<f64>,
    pub kinetic_energy_change_j: f64,
    /// Original endpoint-average nodal internal force dotted with the increment.
    pub internal_work_j: f64,
    pub external_work_j: f64,
    /// Delta K + internal work - external work; never a physical heat source.
    pub energy_balance_residual_j: f64,
    pub maximum_force_residual_n: f64,
    pub maximum_moment_residual_n_m: f64,
    pub maximum_displacement_defect_m: f64,
    pub maximum_rotation_defect_rad: f64,
    retained_kinetic_change: Dd,
}

impl FrameLaminateNewmarkStep {
    /// Retained polarization of the same original kinetic energy.
    pub const fn retained_kinetic_energy_change_j(&self) -> Dd {
        self.retained_kinetic_change
    }
}

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn finite(x: Dd) -> Result<f64, SolidError> {
    let value = x.to_f64();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "native laminate dynamic account is unrepresentable",
        ))
    }
}

// The diagonal route stays O(n); the consistent route contracts its original
// entries without forming another mass matrix or subtracting stored energies.
fn mass_action(mass: &MassMatrix, row: usize, x: impl Fn(usize) -> Dd) -> Dd {
    match mass {
        MassMatrix::Lumped(values) => dd(values[row]) * x(row),
        MassMatrix::Consistent { n, entries } => {
            (0..*n).fold(Dd::ZERO, |sum, j| sum + dd(entries[row * n + j]) * x(j))
        }
    }
}

impl Frame2d {
    /// Advance one undamped, source-bound laminate interval at held nodal load.
    ///
    /// The complete original frame is cloned once by the native temperature
    /// refresh. Its initial displacement and velocity are retained, while the
    /// initial acceleration is resolved from the OLD source force and original
    /// mass. The same existing average-acceleration Newton/line-search kernel
    /// solves the endpoint with the new temperatures. No second integrator,
    /// fitted damping, interpolated material law or hot-density mass is used.
    ///
    /// This supplies prescribed thermal endpoints, not a thermal time solver.
    /// A reciprocal caloric owner must account for material internal-energy
    /// change and stress work. Linear reference geometry, original beam strain
    /// domains, stationary supports and quasistatic cross-section relaxation
    /// apply. Temperature-path/time refinement is still required. Any source,
    /// force, kinematic, work or cancellation refusal returns no changed frame.
    pub fn laminate_newmark_candidate(
        &self,
        temperatures_k: &[Vec<Vec<f64>>],
        held_load: &[f64],
        duration_s: f64,
        controls: FrameLaminateNewmarkControls,
        cx: &Cx<'_>,
    ) -> Result<FrameLaminateNewmarkStep, SolidError> {
        checkpoint(cx)?;
        let n = self.ndof();
        if n == 0
            || n > controls.maximum_dofs
            || n.checked_mul(n).is_none()
            || self.elements.is_empty()
            || self.elements.len() > controls.maximum_members
            || self.geometry != Geometry::Linear
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || held_load.len() != n
            || held_load
                .iter()
                .chain(&self.u)
                .chain(&self.vel)
                .any(|x| !x.is_finite())
            || [
                controls.maximum_force_residual_n,
                controls.maximum_moment_residual_n_m,
                controls.maximum_displacement_defect_m,
                controls.maximum_rotation_defect_rad,
                controls.maximum_energy_defect_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
            || self
                .fixed
                .iter()
                .zip(&self.vel)
                .any(|(&fixed, &v)| fixed && v != 0.0)
        {
            return Err(invalid(
                "bounded linear laminate dynamics, finite inputs and stationary supports required",
            ));
        }
        let dt2 = duration_s * duration_s;
        if !dt2.is_finite() || dt2 == 0.0 {
            return Err(invalid("native Newmark position weight is unrepresentable"));
        }
        for element in 0..self.elements.len() {
            checkpoint(cx)?;
            self.reference_strain_traces(element)?;
        }
        // This also refuses separate viscous/shear operators and validates all
        // original source, reference mass, strength and potential bindings.
        let mut frame = self.with_laminate_temperature_field(temperatures_k, cx)?;
        let mass = self.reference_mass()?;
        let (before_force, _, _) = self.assemble_controlled(&self.u, &mut || checkpoint(cx))?;
        let initial_force: Vec<_> = held_load
            .iter()
            .zip(&before_force)
            .map(|(&p, &f)| p - f)
            .collect();
        let initial_acceleration = mass.acceleration(&initial_force, &self.fixed)?;
        // Massless coordinates must already meet their algebraic equilibrium;
        // a source update cannot silently repair an invalid accepted prefix.
        let mut maximum_force: f64 = 0.0;
        let mut maximum_moment: f64 = 0.0;
        let mut check_force = |force: &[f64], acceleration: &[f64]| -> Result<(), SolidError> {
            for i in 0..n {
                checkpoint(cx)?;
                if self.fixed[i] {
                    continue;
                }
                let residual = finite(
                    mass_action(&mass, i, |j| dd(acceleration[j])) + dd(force[i])
                        - dd(held_load[i]),
                )?
                .abs();
                let (maximum, allowance) = if i % 3 == 2 {
                    (&mut maximum_moment, controls.maximum_moment_residual_n_m)
                } else {
                    (&mut maximum_force, controls.maximum_force_residual_n)
                };
                *maximum = maximum.max(residual);
                if residual > allowance {
                    return Err(invalid(
                        "original laminate dynamic force or moment exceeds its physical allowance",
                    ));
                }
            }
            Ok(())
        };
        check_force(&before_force, &initial_acceleration)?;
        let history = frame.newmark_with_load_initialized(
            &[0.0, 0.0],
            duration_s,
            (0.0, 0.0),
            Rayleigh { a0: 0.0, a1: 0.0 },
            held_load,
            Some(DynamicInitialState {
                internal_force: &before_force,
                acceleration: &initial_acceleration,
            }),
            &mut || checkpoint(cx),
            &mut |_| Ok(()),
        )?;
        let (after_force, _, _) = frame.assemble_controlled(&frame.u, &mut || checkpoint(cx))?;
        check_force(&after_force, &frame.acc)?;
        let mut displacement_increment = Vec::with_capacity(n);
        let (mut kinetic, mut internal, mut external) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        let (mut displacement_defect, mut rotation_defect) = (0.0_f64, 0.0_f64);
        for i in 0..n {
            checkpoint(cx)?;
            // Exactly the same kinematic evaluation as the native kernel.
            let du = duration_s.mul_add(
                self.vel[i],
                dt2 * 0.25_f64.mul_add(frame.acc[i], 0.25 * initial_acceleration[i]),
            );
            let defect = finite(dd(frame.u[i]) - dd(self.u[i]) - dd(du))?.abs();
            let (maximum, allowance) = if i % 3 == 2 {
                (&mut rotation_defect, controls.maximum_rotation_defect_rad)
            } else {
                (
                    &mut displacement_defect,
                    controls.maximum_displacement_defect_m,
                )
            };
            *maximum = maximum.max(defect);
            if defect > allowance {
                return Err(invalid(
                    "resolved laminate motion is not represented within its physical displacement or rotation allowance",
                ));
            }
            displacement_increment.push(du);
            kinetic = kinetic
                + dd(0.5)
                    * (dd(frame.vel[i]) + dd(self.vel[i]))
                    * mass_action(&mass, i, |j| dd(frame.vel[j]) - dd(self.vel[j]));
            internal = internal + dd(0.5) * (dd(before_force[i]) + dd(after_force[i])) * dd(du);
            external = external + dd(held_load[i]) * dd(du);
        }
        let residual = finite(kinetic + internal - external)?;
        let internal_work_j = finite(internal)?;
        let external_work_j = finite(external)?;
        if residual.abs() > controls.maximum_energy_defect_j
            || finite(internal - dd(history.internal_work[1]))?.abs()
                > controls.maximum_energy_defect_j
            || finite(external - dd(history.input_energy[1]))?.abs()
                > controls.maximum_energy_defect_j
        {
            return Err(invalid(
                "native laminate Newmark work and kinetic energy do not close within the original joule allowance",
            ));
        }
        checkpoint(cx)?;
        Ok(FrameLaminateNewmarkStep {
            frame,
            history,
            displacement_increment,
            kinetic_energy_change_j: finite(kinetic)?,
            internal_work_j,
            external_work_j,
            energy_balance_residual_j: residual,
            maximum_force_residual_n: maximum_force,
            maximum_moment_residual_n_m: maximum_moment,
            maximum_displacement_defect_m: displacement_defect,
            maximum_rotation_defect_rad: rotation_defect,
            retained_kinetic_change: kinetic,
        })
    }
}
