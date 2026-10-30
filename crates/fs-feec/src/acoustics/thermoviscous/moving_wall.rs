//! Full-vector no-slip structural motion in the ORIGINAL resolved fluid solve.
//!
//! A boundary lifting v=Phi*u is imposed exactly in the hierarchical P2 trace.
//! The congruence P^T A P gives the reciprocal generalized fluid reaction,
//! including pressure, shear and dilatational stress. Original mass-normalized
//! structural inertia, stiffness and viscous damping enter those same rows;
//! no iteration, acoustic spring, pressure-only force or second loss is added.
//!
//! Shapes are vectors in the mesh's global Cartesian frame, with their actual
//! signs. No face-normal magnitude or automatic sign flip replaces them.
//! This is linear harmonic motion on a FIXED reference domain. It is not ALE,
//! finite deformation, a time advance, or a modal/mesh convergence certificate.
//! The no-slip coupling is full displacement/velocity continuity, unlike a
//! normal-only slip boundary; see
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.aco/aco_ug_multiphysics_couplings.14.06.html>.

use super::*;
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeModel};

/// Hard cap on retained structural coordinates in this local dense owner.
pub const MAX_MOVING_WALL_MODES: usize = 64;

/// One retained ORIGINAL structural model and its conforming velocity traces.
///
/// Mode order must be the source model's order. Every shape has one vector
/// coefficient [1/sqrt(kg)] per GLOBAL P2 scalar DOF: vertices first, then the
/// hierarchical edge bubbles, exactly as in ThermoviscousField::velocity_m_s.
/// Edge coefficients are hierarchical coefficients, NOT midpoint samples.
/// Only DOFs on `faces` may be nonzero; nodes shared with stationary walls
/// must be zero. All selected boundary DOFs follow the sum of these shapes,
/// including tangential components and any identically zero components.
///
/// The source's modal damping must exclude this resolved fluid's loading and
/// losses. Its observer transfers and current time states are retained as
/// source metadata but neither observed nor advanced by this harmonic solve.
/// Source-mode accuracy and the trace transfer remain the caller's obligations;
/// this owner does not infer a nonmatching mortar or an unprovided solid mesh.
#[derive(Clone, Debug)]
pub struct MovingWall {
    /// Moving subset of the original no-slip faces, without duplicates.
    pub faces: Vec<[u32; 3]>,
    /// Original global Cartesian P2 shapes, one per source structural mode.
    pub velocity_shapes: Vec<Vec<[f64; 3]>>,
    /// Original material/geometry-derived mass-normalized structural modes.
    pub structure: ModalAcousticTimeModel,
}

/// Harmonic structural response and its work-conjugate fluid reaction.
#[derive(Clone, Debug, Default)]
pub struct MovingWallField {
    /// Generalized structural velocities [m sqrt(kg)/s], exp(+i omega t).
    pub modal_velocity_m_sqrt_kg_per_s: Vec<C64>,
    /// Total fluid force ON the solid [N/sqrt(kg)], including viscous traction.
    /// Computed from the projected original fluid equations, not from an
    /// independently sampled pressure or a force-balance definition.
    pub fluid_force_n_per_sqrt_kg: Vec<C64>,
    /// External structural drive work [W], already included in total input_w.
    pub external_input_w: f64,
    /// Structural damping [W], already included in structural_damping_w.
    pub damping_w: f64,
    /// Solid-to-fluid traction work [W]. This is an INTERNAL exchange diagnostic,
    /// not an additional loss. Positive means the moving solid supplies fluid work.
    pub power_into_fluid_w: f64,
}

pub(super) fn modal_impedance(mode: &ModalAcousticMode, omega: f64) -> C64 {
    let natural = mode.angular_frequency_rad_s;
    C64::new(
        2.0 * mode.damping_ratio * natural,
        omega - natural * natural / omega,
    )
}

impl<'c> ThermoviscousFe<'c> {
    /// Assemble a general-fluid/structural problem with exact conforming
    /// no-slip modal motion. For a gas, use the original ideal_gas derivatives.
    /// Empty pressure ports are admitted: a closed cavity may be driven entirely
    /// by external structural force. Isothermal and finite conductive patches
    /// retain the original temperature owner on the same moving reference faces.
    ///
    /// The preparation budget covers the released P2 space AND the constrained
    /// mixed system. No boundary velocities are condensed before their original
    /// full viscous, inertial and pressure-work operators have been assembled.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_moving_wall(
        complex: &'c TetComplex,
        positions: &[[f64; 3]],
        medium: ThermoviscousMedium,
        thermodynamics: ThermoviscousThermodynamics,
        boundary: &ThermoviscousBoundary,
        wall: MovingWall,
        limits: ThermoviscousLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousError> {
        poll(cancelled)?;
        let modes = wall.structure.modes().len();
        if wall.faces.is_empty()
            || modes == 0
            || modes > MAX_MOVING_WALL_MODES
            || wall.velocity_shapes.len() != modes
        {
            return Err(ThermoviscousError::InvalidInput(
                "moving wall faces/mode count",
            ));
        }
        let mut faces = BTreeSet::new();
        for &face in &wall.faces {
            if !faces.insert(canonical(face)) {
                return Err(ThermoviscousError::InvalidInput(
                    "duplicate moving wall face",
                ));
            }
        }
        let mut model = Self::assemble_medium(
            complex,
            positions,
            medium,
            thermodynamics,
            false,
            boundary,
            &faces,
            limits,
            cancelled,
        )?;
        model.constrain_moving_wall(wall, &faces, limits, cancelled)?;
        Ok(model)
    }

    /// Original retained structural source and exact trace coefficients.
    pub fn moving_wall(&self) -> Option<&MovingWall> {
        self.moving_wall.as_ref()
    }

    fn constrain_moving_wall(
        &mut self,
        wall: MovingWall,
        faces: &BTreeSet<[u32; 3]>,
        limits: ThermoviscousLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ThermoviscousError> {
        let mut moving = vec![false; self.space.ndof];
        for face in faces {
            for &v in face {
                moving[v as usize] = true;
            }
            for a in 0..3 {
                for b in a + 1..3 {
                    let edge = self
                        .space
                        .complex
                        .edges
                        .binary_search(&[face[a], face[b]])
                        .expect("validated canonical face edge");
                    moving[self.space.edge_off + edge] = true;
                }
            }
        }
        for shape in &wall.velocity_shapes {
            poll(cancelled)?;
            if shape.len() != self.space.ndof {
                return Err(ThermoviscousError::InvalidInput(
                    "global P2 wall shape length",
                ));
            }
            let mut nonzero = false;
            for (i, phi) in shape.iter().enumerate() {
                if phi.iter().any(|v| !v.is_finite())
                    || ((!moving[i] || self.velocity_map[i].is_none())
                        && phi.iter().any(|v| *v != 0.0))
                {
                    return Err(ThermoviscousError::InvalidInput(
                        "wall shape must vanish outside its moving trace and on stationary intersections",
                    ));
                }
                nonzero |= phi.iter().any(|v| *v != 0.0);
            }
            if !nonzero {
                return Err(ThermoviscousError::InvalidInput(
                    "unresolved zero wall shape",
                ));
            }
        }
        let mut count = 0;
        let velocity_map = self
            .velocity_map
            .iter()
            .enumerate()
            .map(|(i, old)| {
                if old.is_some() && !moving[i] {
                    let index = count;
                    count += 1;
                    Some(index)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        let modes = wall.structure.modes().len();
        let pressure_start = 3 * count + modes;
        let temperature_start = pressure_start + self.space.complex.vertex_count;
        let n = temperature_start + self.n - self.temperature_start;
        let cap = limits.max_dofs.min(MAX_LNS_DOFS);
        if n > cap {
            return Err(ThermoviscousError::Budget {
                what: "moving-wall mixed unknowns",
                requested: n,
                limit: cap,
            });
        }
        let peak = n.max(self.n);
        let bytes = peak
            .checked_mul(peak)
            .and_then(|v| v.checked_mul(48))
            .ok_or(ThermoviscousError::NonFinite)?;
        if bytes > limits.max_dense_bytes {
            return Err(ThermoviscousError::Budget {
                what: "moving-wall dense bytes",
                requested: bytes,
                limit: limits.max_dense_bytes,
            });
        }
        // x_full=P*x_constrained. Every row retains its actual Cartesian sign.
        let mut projection = vec![Vec::<(usize, f64)>::new(); self.n];
        for (i, &old) in self.velocity_map.iter().enumerate() {
            if let Some(old) = old {
                for c in 0..3 {
                    let row = &mut projection[c * self.velocity_count + old];
                    if let Some(new) = velocity_map[i] {
                        row.push((c * count + new, 1.0));
                    } else {
                        for (mode, shape) in wall.velocity_shapes.iter().enumerate() {
                            if shape[i][c] != 0.0 {
                                row.push((3 * count + mode, shape[i][c]));
                            }
                        }
                    }
                }
            }
        }
        for i in self.pressure_start..self.n {
            projection[i].push((pressure_start + i - self.pressure_start, 1.0));
        }
        let constant = congruence(&self.constant, self.n, n, &projection, cancelled)?;
        let mass = congruence(&self.mass, self.n, n, &projection, cancelled)?;
        let mut loads = Vec::with_capacity(self.loads.len());
        for old in &self.loads {
            poll(cancelled)?;
            let mut load = vec![0.0; n];
            for (row, terms) in projection.iter().enumerate() {
                for &(col, coefficient) in terms {
                    load[col] += coefficient * old[row];
                }
            }
            loads.push(load);
        }
        if constant
            .iter()
            .chain(&mass)
            .chain(loads.iter().flatten())
            .any(|x| !x.is_finite())
        {
            return Err(ThermoviscousError::NonFinite);
        }
        self.constant = constant;
        self.mass = mass;
        self.loads = loads;
        self.velocity_map = velocity_map;
        self.velocity_count = count;
        self.pressure_start = pressure_start;
        self.temperature_start = temperature_start;
        self.n = n;
        self.moving_wall = Some(wall);
        poll(cancelled)
    }
}

/// Two sparse/dense multiplies avoid a quadratic expansion in the number of
/// shapes at each boundary node. Peak dense storage stays below 48*max(n)^2.
fn congruence(
    original: &[f64],
    old_n: usize,
    n: usize,
    projection: &[Vec<(usize, f64)>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, ThermoviscousError> {
    let mut right = vec![0.0; old_n * n];
    for i in 0..old_n {
        poll(cancelled)?;
        for j in 0..old_n {
            let a = original[i * old_n + j];
            if a != 0.0 {
                for &(col, b) in &projection[j] {
                    right[i * n + col] += a * b;
                }
            }
        }
    }
    let mut result = vec![0.0; n * n];
    for i in 0..old_n {
        poll(cancelled)?;
        for &(row, a) in &projection[i] {
            for j in 0..n {
                result[row * n + j] += a * right[i * n + j];
            }
        }
    }
    Ok(result)
}

impl PreparedThermoviscous<'_, '_> {
    /// Solve external structural forces and pressure tractions simultaneously.
    /// Forces are generalized [N/sqrt(kg)], in the retained source's mode order
    /// and exp(+i omega t). The same prepared LU serves every combined drive.
    pub fn response_with_wall_forces(
        &self,
        drive: &[C64],
        wall_forces: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousField, ThermoviscousError> {
        poll(cancelled)?;
        let model = self.model;
        let modes = model.pressure_start - 3 * model.velocity_count;
        if drive.len() != model.loads.len()
            || wall_forces.len() != modes
            || drive.iter().chain(wall_forces).any(|z| !finite(*z))
        {
            return Err(ThermoviscousError::InvalidInput(
                "traction/structural force drive",
            ));
        }
        let mut rhs = vec![C64::ZERO; model.n];
        for (load, &p) in model.loads.iter().zip(drive) {
            for (entry, &b) in rhs.iter_mut().zip(load) {
                *entry = *entry + p.scale(b);
            }
        }
        for (i, &force) in wall_forces.iter().enumerate() {
            rhs[3 * model.velocity_count + i] = rhs[3 * model.velocity_count + i] + force;
        }
        let x = model.solve_system(
            &self.factor,
            &self.rows,
            &self.columns,
            std::f64::consts::TAU * self.frequency_hz,
            self.cp,
            &self.wall_admittances,
            None,
            &rhs,
            cancelled,
        )?;
        self.field(drive, wall_forces, &x, cancelled)
    }

    pub(super) fn wall_field(
        &self,
        drive: &[C64],
        forces: &[C64],
        x: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingWallField, ThermoviscousError> {
        let mut result = MovingWallField::default();
        let model = self.model;
        let Some(wall) = &model.moving_wall else {
            return Ok(result);
        };
        let omega = std::f64::consts::TAU * self.frequency_hz;
        for (i, mode) in wall.structure.modes().iter().enumerate() {
            poll(cancelled)?;
            let row = 3 * model.velocity_count + i;
            let u = x[row];
            // -Phi^T (A_fluid*x-B_pressure*P) is the discrete total traction
            // ON the solid. Interior equations enforce the same fluid motion.
            // Near static compression the conservative terms are much larger
            // than traction work. Retain the same extended precision as the
            // original physical-equation residual, rather than subtracting
            // rounded pressure and inertia samples.
            let mut re = Dd::ZERO;
            let mut im = Dd::ZERO;
            for (col, &value) in x.iter().enumerate() {
                let a = model.fluid_operator_entry(row, col, omega, self.cp);
                re = re
                    - (Dd::from_f64(a.re) * Dd::from_f64(value.re)
                        - Dd::from_f64(a.im) * Dd::from_f64(value.im));
                im = im
                    - (Dd::from_f64(a.re) * Dd::from_f64(value.im)
                        + Dd::from_f64(a.im) * Dd::from_f64(value.re));
            }
            for (load, &p) in model.loads.iter().zip(drive) {
                re = re + Dd::from_f64(p.re) * Dd::from_f64(load[row]);
                im = im + Dd::from_f64(p.im) * Dd::from_f64(load[row]);
            }
            let force = C64::new(re.to_f64(), im.to_f64());
            result.external_input_w +=
                0.5 * (forces.get(i).copied().unwrap_or(C64::ZERO).conj() * u).re;
            result.damping_w += mode.damping_ratio * mode.angular_frequency_rad_s * u.norm_sq();
            result.power_into_fluid_w -= 0.5 * (force.conj() * u).re;
            result.modal_velocity_m_sqrt_kg_per_s.push(u);
            result.fluid_force_n_per_sqrt_kg.push(force);
        }
        if result.fluid_force_n_per_sqrt_kg.iter().any(|z| !finite(*z))
            || [
                result.external_input_w,
                result.damping_w,
                result.power_into_fluid_w,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(result)
    }
}
