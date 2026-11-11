//! Full-vector no-slip structural motion in the ORIGINAL resolved fluid solve.
//!
//! A boundary lifting v=Phi*u is imposed exactly in the hierarchical P2 trace.
//! The congruence P^T A P gives the reciprocal generalized fluid reaction,
//! including pressure, shear and dilatational stress. Original mass-normalized
//! structural inertia, stiffness and viscous damping enter those same rows;
//! no iteration, acoustic spring, pressure-only force or second loss is added.
//! Optional original passive damping channels retain the complete cross-mode
//! material/joint operator in harmonic and causal solves. Their signed channel
//! work also retains undamped common motion and other exact loss null spaces.
//!
//! Shapes are vectors in the mesh's global Cartesian frame, with their actual
//! signs. No face-normal magnitude or automatic sign flip replaces them.
//! This is linear motion on a FIXED reference domain. A causal advance uses
//! the separately prepared time cache. The attachment is not ALE, finite
//! deformation, or a modal/mesh convergence certificate.
//! The no-slip coupling is full displacement/velocity continuity, unlike a
//! normal-only slip boundary; see
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.aco/aco_ug_multiphysics_couplings.14.06.html>.

use super::*;
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeModel};

/// Original passive material histories in harmonic and causal wall motion.
pub mod relaxation;

/// Hard cap on retained structural coordinates in this local dense owner.
pub const MAX_MOVING_WALL_MODES: usize = 64;

/// Original passive channels and their complete cached Gram matrix [1/s].
/// Both harmonic and causal operators read this one matrix; harmonic heat
/// retains the channels to avoid cancellation in a quadratic matrix contraction.
pub(super) struct MovingWallDamping {
    channels: Vec<Vec<f64>>,
    matrix: Vec<f64>,
}

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
    /// Original material energy-coordinate peak phasors [sqrt(J)].
    pub relaxation_coordinates_sqrt_j: Vec<C64>,
    /// Cycle-mean stored structural-memory energy [J].
    pub relaxation_mean_energy_j: f64,
    /// Structural-memory heat [W], included in structural_relaxation_w.
    pub relaxation_w: f64,
    /// Solid-to-fluid traction work [W]. This is an INTERNAL exchange diagnostic,
    /// not an additional loss. Positive means the moving solid supplies fluid work.
    pub power_into_fluid_w: f64,
}

pub(super) fn modal_impedance(
    mode: &ModalAcousticMode,
    damping_rate_per_s: f64,
    omega: f64,
) -> C64 {
    let natural = mode.angular_frequency_rad_s;
    C64::new(damping_rate_per_s, omega - natural * natural / omega)
}

fn moving_faces(wall: &MovingWall) -> Result<BTreeSet<[u32; 3]>, ThermoviscousError> {
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
    Ok(faces)
}

impl<'c> ThermoviscousFe<'c> {
    /// Attach full material/joint damping in the ORIGINAL mass-normalized
    /// structural coordinates, before preparing harmonic or causal motion.
    /// Channel k supplies `l_k` [1/sqrt(s)] and defines
    /// `C = sum_k l_k l_k^T` [1/s]. Its resisting generalized force is `C v`;
    /// instantaneous heat is `sum_k (l_k dot v)^2`. Every signed cross term
    /// survives, including rank-deficient damping and undamped common motion.
    ///
    /// The same original moving wall, fluid, shapes, inertia and stiffness
    /// remain. Existing source damping rates must ALL be zero, and a second
    /// attachment refuses. Channels must exclude this resolved fluid's losses.
    /// Empty channels explicitly select zero structural dissipation. A missing
    /// attachment preserves the source's original diagonal damping unchanged.
    ///
    /// At most 64 channels are admitted, each covering all original modes.
    /// Matrix-only source owners can supply their existing passive factor;
    /// this method does not fit, diagonalize, clip or invent a material law.
    /// Matching the source's physical basis, not merely its count, is the
    /// caller's responsibility. Cancellation publishes no partial attachment.
    pub fn with_moving_wall_damping_channels(
        mut self,
        channels_per_sqrt_s: Vec<Vec<f64>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousError> {
        poll(cancelled)?;
        let wall = self
            .moving_wall
            .as_ref()
            .ok_or(ThermoviscousError::InvalidInput(
                "structural damping requires an original moving wall",
            ))?;
        let modes = wall.structure.modes().len();
        if self.moving_wall_damping.is_some()
            || wall
                .structure
                .damping_rates_per_s()
                .iter()
                .any(|&c| c != 0.0)
        {
            return Err(ThermoviscousError::InvalidInput(
                "full structural damping and existing damping owners are exclusive",
            ));
        }
        if channels_per_sqrt_s.len() > MAX_MOVING_WALL_MODES {
            return Err(ThermoviscousError::Budget {
                what: "moving-wall damping channels",
                requested: channels_per_sqrt_s.len(),
                limit: MAX_MOVING_WALL_MODES,
            });
        }
        for channel in &channels_per_sqrt_s {
            poll(cancelled)?;
            if channel.len() != modes || channel.iter().any(|x| !x.is_finite()) {
                return Err(ThermoviscousError::InvalidInput(
                    "finite damping channels in the complete original wall basis required",
                ));
            }
        }
        let mut matrix = vec![0.0; modes * modes];
        for i in 0..modes {
            poll(cancelled)?;
            for j in 0..=i {
                let mut value = Dd::ZERO;
                for channel in &channels_per_sqrt_s {
                    value = value + Dd::from_f64(channel[i]) * Dd::from_f64(channel[j]);
                }
                let value = value.to_f64();
                if !value.is_finite()
                    || (i == j
                        && value == 0.0
                        && channels_per_sqrt_s.iter().any(|channel| channel[i] != 0.0))
                {
                    return Err(ThermoviscousError::NonFinite);
                }
                matrix[i * modes + j] = value;
                matrix[j * modes + i] = value;
            }
        }
        poll(cancelled)?;
        self.moving_wall_damping = Some(MovingWallDamping {
            channels: channels_per_sqrt_s,
            matrix,
        });
        Ok(self)
    }

    /// Original attached passive channels [1/sqrt(s)], in structural order.
    /// None means the source's original per-coordinate damping remains active.
    pub fn moving_wall_damping_channels(&self) -> Option<&[Vec<f64>]> {
        self.moving_wall_damping
            .as_ref()
            .map(|owner| owner.channels.as_slice())
    }

    /// Complete attached structural C [1/s], row-major, without diagonalization.
    pub fn moving_wall_damping_matrix_per_s(&self) -> Option<&[f64]> {
        self.moving_wall_damping
            .as_ref()
            .map(|owner| owner.matrix.as_slice())
    }

    pub(super) fn wall_damping_entry(&self, row: usize, column: usize) -> f64 {
        let structure = &self
            .moving_wall
            .as_ref()
            .expect("moving wall coordinate")
            .structure;
        if let Some(owner) = &self.moving_wall_damping {
            owner.matrix[row * structure.modes().len() + column]
        } else if row == column {
            structure.damping_rate_per_s(row)
        } else {
            0.0
        }
    }

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
        let faces = moving_faces(&wall)?;
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

    /// Consume the original admitted geometry with full-vector conforming
    /// structural motion. The same mapped inertial, pressure and viscous forms
    /// are assembled before the original exact congruence imposes the wall
    /// velocity. Thus curved walls retain tangential no-slip work and the full
    /// reciprocal fluid reaction, without replacing vector shapes by normals.
    ///
    /// The P2 shapes remain hierarchical coefficients in global Cartesian
    /// axes on a fixed reference domain. This does not infer solid geometry or
    /// validate its modal transfer; it is not a moving-mesh/ALE solve.
    #[allow(clippy::too_many_arguments)]
    pub fn from_geometry_with_moving_wall(
        geometry: HelmholtzFe<'c>,
        medium: ThermoviscousMedium,
        thermodynamics: ThermoviscousThermodynamics,
        boundary: &ThermoviscousBoundary,
        wall: MovingWall,
        limits: ThermoviscousLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousError> {
        poll(cancelled)?;
        let faces = moving_faces(&wall)?;
        let mut model = Self::assemble_geometry(
            geometry,
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
        let mut mass = congruence(&self.mass, self.n, n, &projection, cancelled)?;
        // P^T M P is symmetric in exact arithmetic when M is, and the two sparse passes of `congruence` sum the pair (a, b) and
        // the pair (b, a) in different orders, so the computed triangles are two roundings of one value: a pair whose exact
        // value is 0 (a wall shape that is antisymmetric about the centre against a symmetric interior function) came out as
        // a residue of 2e-22 in one triangle (2e-17 after the normalization by the diagonal) and exactly 0 in the other.
        // Keeping the lower triangle, the one the Cholesky factorization of this matrix reads, is exact bookkeeping of one
        // symmetric value, not a tolerance: the stored and the factored matrix are the same. `constant` is not symmetric and
        // keeps the plain product.
        for i in 0..n {
            for j in 0..i {
                mass[j * n + i] = mass[i * n + j];
            }
        }
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
            &self.wall_relaxation,
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
        for (i, &damping) in wall.structure.damping_rates_per_s().iter().enumerate() {
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
            if model.moving_wall_damping.is_none() {
                result.damping_w += 0.5 * damping * u.norm_sq();
            }
            result.power_into_fluid_w -= 0.5 * (force.conj() * u).re;
            result.modal_velocity_m_sqrt_kg_per_s.push(u);
            result.fluid_force_n_per_sqrt_kg.push(force);
        }
        if let Some(owner) = &model.moving_wall_damping {
            for channel in &owner.channels {
                poll(cancelled)?;
                let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
                for (&coefficient, &velocity) in
                    channel.iter().zip(&result.modal_velocity_m_sqrt_kg_per_s)
                {
                    re = re + Dd::from_f64(coefficient) * Dd::from_f64(velocity.re);
                    im = im + Dd::from_f64(coefficient) * Dd::from_f64(velocity.im);
                }
                result.damping_w += 0.5 * C64::new(re.to_f64(), im.to_f64()).norm_sq();
            }
        }
        self.wall_relaxation_field(&mut result, cancelled)?;
        if result.fluid_force_n_per_sqrt_kg.iter().any(|z| !finite(*z))
            || [
                result.external_input_w,
                result.damping_w,
                result.relaxation_w,
                result.relaxation_mean_energy_j,
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
