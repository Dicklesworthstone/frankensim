//! Finite rigid payload inertia at the existing reciprocal shell ports.
//!
//! A payload follows the interpolated small rotation and the velocity at its
//! offset center of mass. Its kinetic form adds B^T diag(m I, I_com) B to the
//! original full mass, before support restriction. No new timestepper, guessed
//! nodal mass lumping, or payload damping is introduced.

use super::{MovingLoadError, MovingShellModel, ShellLoadPoint, cross, finite, invalid, poll};
use crate::shell::ShellPlate;

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
    /// Payload locations cannot change after assembly. Apply any gravity or
    /// other external payload force explicitly through its existing COM port.
    /// Numerically unresolved semidefinite inertia is refused, never repaired
    /// by clipping a negative pivot. No new finite-rotation claim is made.
    pub fn new_with_payloads(
        plate: ShellPlate,
        max_load_points: usize,
        payloads: &[ShellRigidPayload],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MovingLoadError> {
        poll(&mut cancelled)?;
        if payloads.len() > max_load_points || payloads.len() > 4096 {
            return Err(invalid("attached shell payload point budget exceeded"));
        }
        for payload in payloads {
            poll(&mut cancelled)?;
            payload.check()?;
        }
        let mut model = Self::new(plate, max_load_points, &mut cancelled)?;
        if payloads.is_empty() {
            return Ok(model);
        }
        let n = model.assembly.full_mass.dimension();
        let mut mass = model.assembly.full_mass.values().to_vec();
        for payload in payloads {
            poll(&mut cancelled)?;
            let (nodes, _) = model.point(payload.point)?;
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
        model
            .assembly
            .replace_mass_preserving_damping(mass, model.plate.budget.max_conditioning_dofs)
            .map_err(MovingLoadError::Spatial)?;
        model.payloads = payloads.to_vec();
        poll(&mut cancelled)?;
        Ok(model)
    }

    /// Fixed physical payload descriptors in their assembly order.
    pub fn payloads(&self) -> &[ShellRigidPayload] {
        &self.payloads
    }
}
