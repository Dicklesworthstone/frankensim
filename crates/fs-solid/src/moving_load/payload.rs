//! Finite rigid payload inertia at the existing reciprocal shell ports.
//!
//! A payload follows the interpolated small rotation and the velocity at its
//! offset center of mass. Its kinetic form adds B^T diag(m I, I_com) B to the
//! original full mass, before support restriction. No new timestepper, guessed
//! nodal mass lumping, or payload damping is introduced.

use super::{MovingLoadError, MovingShellModel, ShellLoadPoint, cross, finite, invalid, poll};
use crate::shell::{ShellMassFormulation, ShellPlate};

/// One permanently attached rigid payload in the shell's reference SI frame.
///
/// The point's offset locates the payload center of mass relative to the P1
/// midsurface location. Its COM translation and small Cartesian rotation use
/// the same material motion as force/couple ports. Membership and orientation
/// remain fixed; this does not describe sliding, impacts, or finite rotations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellRigidPayload {
    /// Material attachment and reference arm to the center of mass, m.
    pub point: ShellLoadPoint,
    /// Positive finite payload mass, kg.
    pub mass_kg: f64,
    /// Exactly symmetric COM rotational inertia in the model frame, kg m².
    /// The physical second-moment tensor trace(I)/2 * identity - I must be
    /// positive semidefinite. An all-zero tensor represents a point mass.
    pub inertia_com_kg_m2: [[f64; 3]; 3],
}

impl ShellRigidPayload {
    fn check(&self) -> Result<(), MovingLoadError> {
        if !self.mass_kg.is_finite() || self.mass_kg <= 0.0 {
            return Err(invalid("payload mass must be finite and positive"));
        }
        let inertia = &self.inertia_com_kg_m2;
        for row in inertia {
            finite(row)?;
        }
        for i in 0..3 {
            for j in 0..i {
                if inertia[i][j] != inertia[j][i] {
                    return Err(invalid("payload COM inertia must be exactly symmetric"));
                }
            }
        }
        let scale = inertia.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
        if scale == 0.0 {
            return Ok(());
        }
        // Normalize before taking the trace or products, so physical units and
        // magnitude do not impose an absolute pivot cutoff. A PSD second moment
        // enforces both nonnegative principal inertias and their triangle laws.
        let normalized = inertia.map(|row| row.map(|value| value / scale));
        let half_trace = 0.5 * (normalized[0][0] + normalized[1][1] + normalized[2][2]);
        let second: [[f64; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                if i == j {
                    half_trace - normalized[i][j]
                } else {
                    -normalized[i][j]
                }
            })
        });
        let mut factor = [[0.0; 3]; 3];
        for j in 0..3 {
            let pivot = second[j][j] - (0..j).map(|k| factor[j][k] * factor[j][k]).sum::<f64>();
            if !pivot.is_finite() || pivot < 0.0 {
                return Err(invalid("payload COM inertia has an inadmissible second moment"));
            }
            factor[j][j] = pivot.sqrt();
            for i in j + 1..3 {
                let residual =
                    second[i][j] - (0..j).map(|k| factor[i][k] * factor[j][k]).sum::<f64>();
                if pivot == 0.0 {
                    if residual != 0.0 {
                        return Err(invalid("payload COM inertia has an unresolved zero pivot"));
                    }
                } else {
                    factor[i][j] = residual / factor[j][j];
                    finite(&[factor[i][j]])?;
                }
            }
        }
        Ok(())
    }

    fn spatial_inertia(&self) -> Result<[[f64; 6]; 6], MovingLoadError> {
        // Unit generalized motions at the midsurface produce precisely the
        // port velocity u_dot + theta_dot cross offset, not its negative.
        let mut translation = [[0.0; 3]; 6];
        for axis in 0..3 {
            translation[axis][axis] = 1.0;
            let mut rotation = [0.0; 3];
            rotation[axis] = 1.0;
            translation[axis + 3] = cross(rotation, self.point.offset_m);
        }
        let mut spatial = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in i..6 {
                let value = self.mass_kg
                    * (0..3).map(|a| translation[i][a] * translation[j][a]).sum::<f64>()
                    + if i >= 3 && j >= 3 {
                        self.inertia_com_kg_m2[i - 3][j - 3]
                    } else {
                        0.0
                    };
                finite(&[value])?;
                spatial[i][j] = value;
                spatial[j][i] = value;
            }
        }
        Ok(spatial)
    }
}

impl MovingShellModel {
    /// Assemble native shell dynamics with permanently attached finite inertia.
    ///
    /// The same explicit point bound caps each endpoint load set and the
    /// independently stored payload set. Assembly adds each payload's full
    /// 18-by-18 P1 inertial coupling once, including offset/rotation cross terms.
    /// The original shell stiffness and damping stay unchanged: a payload does
    /// not acquire a fictitious Rayleigh drag. All native solves, prescribed
    /// pin terms, reactions, momentum and energy then use the augmented mass.
    ///
    /// Payload locations cannot change after assembly. Supply uniform gravity
    /// explicitly with the body-acceleration entry points, which include shell
    /// and payload weight once, or author external forces at the COM port.
    /// Numerically unresolved semidefinite inertia is refused, never repaired
    /// by clipping a negative pivot. No new finite-rotation claim is made.
    pub fn new_with_payloads(
        plate: ShellPlate,
        max_load_points: usize,
        payloads: &[ShellRigidPayload],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MovingLoadError> {
        Self::new_with_payloads_and_mass_formulation(
            plate,
            ShellMassFormulation::Lumped,
            max_load_points,
            payloads,
            cancelled,
        )
    }

    /// Attach the same physical payloads to a selected reference shell mass.
    /// Payload inertia is added once to the full mass before support reduction;
    /// Rayleigh damping still belongs only to the selected bare shell inertia.
    /// The original [`Self::new_with_payloads`] entry selects lumped shell mass.
    pub fn new_with_payloads_and_mass_formulation(
        plate: ShellPlate,
        mass_formulation: ShellMassFormulation,
        max_load_points: usize,
        payloads: &[ShellRigidPayload],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MovingLoadError> {
        check_payloads(payloads, max_load_points, &mut cancelled)?;
        let model = Self::new_with_mass_formulation(
            plate,
            mass_formulation,
            max_load_points,
            &mut cancelled,
        )?;
        model.attach_payloads(payloads, cancelled)
    }

    /// Combine a graded reference-density shell with permanently attached
    /// payloads. Each triangle's three values are bulk densities [kg/m³] at
    /// its original corners; the native P1 density integration owns shell
    /// mass, while each payload supplies its independent mass and COM inertia.
    ///
    /// The original density field is retained without averaging it into a
    /// uniform material. Shell and payload inertia are assembled before any
    /// support reduction, so contact, impulses, prescribed pins and Newmark
    /// all use the same total mass. Rayleigh damping uses only the graded
    /// bare shell mass; adding payloads does not manufacture material loss.
    /// No density evolution, mass transfer or finite rotation is inferred.
    pub fn new_with_payloads_and_reference_density(
        plate: ShellPlate,
        reference_density_kg_m3: Vec<[f64; 3]>,
        mass_formulation: ShellMassFormulation,
        max_load_points: usize,
        payloads: &[ShellRigidPayload],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MovingLoadError> {
        check_payloads(payloads, max_load_points, &mut cancelled)?;
        let model = Self::new_with_reference_density(
            plate,
            reference_density_kg_m3,
            mass_formulation,
            max_load_points,
            &mut cancelled,
        )?;
        model.attach_payloads(payloads, cancelled)
    }

    fn attach_payloads(
        mut self,
        payloads: &[ShellRigidPayload],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MovingLoadError> {
        if payloads.is_empty() {
            return Ok(self);
        }
        let n = self.assembly.full_mass.dimension();
        let mut mass = self.assembly.full_mass.values().to_vec();
        for payload in payloads {
            poll(&mut cancelled)?;
            let (nodes, _) = self.point(payload.point)?;
            let spatial = payload.spatial_inertia()?;
            for i in 0..18 {
                let row = 6 * nodes[i / 6] + i % 6;
                for j in i..18 {
                    let column = 6 * nodes[j / 6] + j % 6;
                    let value = payload.point.barycentric[i / 6]
                        * payload.point.barycentric[j / 6]
                        * spatial[i % 6][j % 6];
                    mass[row * n + column] += value;
                    if row != column {
                        mass[column * n + row] += value;
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        self
            .assembly
            .replace_mass_preserving_damping(mass, self.plate.budget.max_conditioning_dofs)
            .map_err(MovingLoadError::Spatial)?;
        self.body_acceleration_map =
            Self::prepare_body_acceleration_map(&self.assembly, &mut cancelled)?;
        self.payloads = payloads.to_vec();
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Fixed physical payload descriptors in their assembly order.
    pub fn payloads(&self) -> &[ShellRigidPayload] {
        &self.payloads
    }
}

fn check_payloads(
    payloads: &[ShellRigidPayload],
    max_load_points: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MovingLoadError> {
    poll(cancelled)?;
    if payloads.len() > max_load_points || payloads.len() > 4096 {
        return Err(invalid("attached shell payload point budget exceeded"));
    }
    for payload in payloads {
        poll(cancelled)?;
        payload.check()?;
    }
    Ok(())
}
