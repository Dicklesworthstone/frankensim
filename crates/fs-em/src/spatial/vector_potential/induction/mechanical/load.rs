//! Consistent P1 load integration and canonical mechanical-frame transfer.

use fs_ga::{GaError, Se3, Vec3, Wrench};

/// Explicit allocation and selected-cell work admission for a Lorentz load.
#[derive(Clone, Copy, Debug)]
pub struct SpatialLorentzLoadLimits {
    /// Maximum original mesh vertices (one nodal load per original vertex).
    pub max_vertices: usize,
    /// Maximum selected conducting tetrahedra.
    pub max_cells: usize,
}

/// Transactional refusal during mechanical load extraction or transformation.
#[derive(Debug)]
pub enum SpatialLorentzLoadError {
    /// Invalid selection, budget, geometry or nonfinite arithmetic.
    Invalid(&'static str),
    /// Cancellation requested before extraction completes.
    Cancelled,
    /// A receiving rigid-frame transformation refused.
    Geometry(GaError),
}

impl core::fmt::Display for SpatialLorentzLoadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "invalid Lorentz load: {reason}"),
            Self::Cancelled => f.write_str("Lorentz load extraction cancelled"),
            Self::Geometry(error) => write!(f, "Lorentz load frame: {error}"),
        }
    }
}
impl std::error::Error for SpatialLorentzLoadError {}
impl From<GaError> for SpatialLorentzLoadError {
    fn from(error: GaError) -> Self {
        Self::Geometry(error)
    }
}

type Result<T> = core::result::Result<T, SpatialLorentzLoadError>;

/// Source-bound mechanical forces from one selected induction stage.
///
/// Nodal forces preserve all P1 virtual work, not only rigid-body work. The
/// reduced wrench retains both net force and moment; extracting it is O(1)
/// after the selected elements have been integrated. Loads are immutable.
#[derive(Debug)]
pub struct SpatialLorentzLoads {
    nodal_forces_n: Vec<[f64; 3]>,
    reference_m: [f64; 3],
    torque_nm: [f64; 3],
    force_n: [f64; 3],
}

impl SpatialLorentzLoads {
    /// Forces [N] in the original magnetic mesh's coordinates and vertex order.
    pub fn nodal_forces_n(&self) -> &[[f64; 3]] {
        &self.nodal_forces_n
    }

    /// Canonical `[torque, force]` wrench in the receiving mechanical frame.
    ///
    /// `placement` maps the magnetic mesh into that frame; `moment_origin_m`
    /// is expressed in the SAME receiving frame. Choose the receiving body's
    /// origin and axes for native `fs-mbd` body-wrench or impulse inputs.
    /// Forces are N and torques N m: multiply both by a held-load duration
    /// before supplying a momentum impulse. This method does not advance time.
    pub fn wrench(&self, placement: &Se3, moment_origin_m: [f64; 3]) -> Result<Wrench> {
        for value in moment_origin_m {
            finite(value)?;
        }
        let force = placement.transform_vector(vector(self.force_n))?;
        let moment = placement.transform_vector(vector(self.torque_nm))?;
        let reference = placement.transform_point(vector(self.reference_m))?;
        let arm = reference - vector(moment_origin_m);
        let wrench = Wrench::new(moment + arm.cross(force), force);
        for value in wrench.to_array() {
            finite(value)?;
        }
        Ok(wrench)
    }

    /// Exact discrete mechanical power [W] for an original-mesh P1 velocity.
    ///
    /// Velocities use the same original coordinates and vertex order as the
    /// nodal loads. This also admits virtual velocities independently of the
    /// prescribed velocity that generated the induction stage.
    pub fn power_w(
        &self,
        nodal_velocity_m_s: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        if nodal_velocity_m_s.len() != self.nodal_forces_n.len() {
            return Err(SpatialLorentzLoadError::Invalid(
                "velocity length differs from original mesh",
            ));
        }
        let (mut sum, mut correction) = (0.0, 0.0);
        for (velocity, force) in nodal_velocity_m_s.iter().zip(&self.nodal_forces_n) {
            poll(cancelled)?;
            for axis in 0..3 {
                finite(velocity[axis])?;
                add(&mut sum, &mut correction, velocity[axis] * force[axis])?;
            }
        }
        finite(sum + correction)
    }
}

pub(crate) struct CellLoad<'a> {
    pub vertices: [u32; 4],
    pub volume_m3: f64,
    pub force_density_n_m3: &'a [[f64; 3]; 4],
}

pub(crate) fn integrate<'a>(
    positions_m: &[[f64; 3]],
    cells: impl ExactSizeIterator<Item = CellLoad<'a>>,
    limits: SpatialLorentzLoadLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialLorentzLoads> {
    use SpatialLorentzLoadError::Invalid;
    poll(cancelled)?;
    if positions_m.is_empty()
        || positions_m.len() > limits.max_vertices
        || cells.len() == 0
        || cells.len() > limits.max_cells
    {
        return Err(Invalid("empty or over-budget Lorentz load integration"));
    }
    let mut nodal_forces_n = vec![[0.0; 3]; positions_m.len()];
    let mut corrections = vec![[0.0; 3]; positions_m.len()];
    let mut reference_m = None;
    for cell in cells {
        poll(cancelled)?;
        if !cell.volume_m3.is_finite() || cell.volume_m3 <= 0.0 {
            return Err(Invalid("Lorentz cell volume must be finite and positive"));
        }
        let weight = cell.volume_m3 / 20.0;
        if weight == 0.0 {
            return Err(Invalid("Lorentz cell mass weights are not representable"));
        }
        for vertex in cell.vertices {
            let position = positions_m
                .get(vertex as usize)
                .ok_or(Invalid("Lorentz cell vertex is outside its original mesh"))?;
            for &value in position {
                finite(value)?;
            }
            reference_m.get_or_insert(*position);
        }
        for a in 0..4 {
            let node = cell.vertices[a] as usize;
            for b in 0..4 {
                // Integral N_a N_b dV = V*(1 + delta_ab)/20.
                // Multiply the smaller mass weight first, avoiding the sum of
                // possibly large force samples before their volume scaling.
                let mass = if a == b { 2.0 * weight } else { weight };
                for axis in 0..3 {
                    let density = finite(cell.force_density_n_m3[b][axis])?;
                    add(
                        &mut nodal_forces_n[node][axis],
                        &mut corrections[node][axis],
                        mass * density,
                    )?;
                }
            }
        }
    }
    let reference_m = reference_m.ok_or(Invalid("missing Lorentz load reference"))?;
    let (mut total, mut correction) = ([0.0; 6], [0.0; 6]);
    for (node, force) in nodal_forces_n.iter_mut().enumerate() {
        poll(cancelled)?;
        for axis in 0..3 {
            force[axis] = finite(force[axis] + corrections[node][axis])?;
            add(&mut total[axis + 3], &mut correction[axis + 3], force[axis])?;
        }
        // Unselected mesh vertices have no load; their coordinates need not
        // enter the selected body's moment reduction.
        if *force != [0.0; 3] {
            let moment = (vector(positions_m[node]) - vector(reference_m)).cross(vector(*force));
            for (axis, value) in [moment.x, moment.y, moment.z].into_iter().enumerate() {
                add(&mut total[axis], &mut correction[axis], value)?;
            }
        }
    }
    for axis in 0..6 {
        total[axis] = finite(total[axis] + correction[axis])?;
    }
    Ok(SpatialLorentzLoads {
        nodal_forces_n,
        reference_m,
        torque_nm: [total[0], total[1], total[2]],
        force_n: [total[3], total[4], total[5]],
    })
}

fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SpatialLorentzLoadError::Invalid(
            "nonfinite Lorentz load arithmetic",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialLorentzLoadError::Cancelled)
    } else {
        Ok(())
    }
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<()> {
    finite(value)?;
    let next = finite(*sum + value)?;
    *correction = finite(
        *correction
            + if sum.abs() >= value.abs() {
                (*sum - next) + value
            } else {
                (value - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}

#[cfg(test)]
#[path = "load/tests.rs"]
mod tests;
