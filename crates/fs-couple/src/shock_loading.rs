//! Reflected shock loading on stationary, explicitly matched structural surfaces.
//!
//! The original Euler step owns wall pressure and its time integral. This
//! adapter consumes that accepted impulse instead of sampling an endpoint
//! pressure or replacing the shock with a linear acoustic wave. The original
//! shell pressure port owns the spatial P1 load and its conjugate volume row.
//! Cold preparation compresses those rows; interval projection visits only
//! retained nodes. No gas or structural state is advanced by this adapter.
//!
//! Forces are the held interval means. Their impulse and spatial virtual work
//! are preserved; unresolved time variation still requires a timestep ladder.
//! A stationary wall does no gas work. Driving a deformable body with these
//! loads is one-way prescribed loading, not a reciprocal moving-wall FSI solve.

use fs_gas::dynamics::{Boundaries, Boundary, Euler1d, StepReport};
use fs_plate::shell::{
    ShellMesh,
    pressure::{OrientedShellFacet, ShellFace, ShellPressurePort},
    profile::ProfileBudget,
};
use std::collections::BTreeSet;

#[cfg(feature = "experimental-compressible")]
mod cartesian;

/// Input, geometry and cancellation refusals publish no partial load.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// Inconsistent geometry, boundary, clock or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The original shell pressure owner refused the surface.
    Surface(fs_plate::PlateError),
    /// Cancellation before returning a complete result.
    Cancelled,
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Euler wall/structural loading: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<fs_plate::PlateError> for Error {
    fn from(error: fs_plate::PlateError) -> Self {
        if error == fs_plate::PlateError::Cancelled {
            Self::Cancelled
        } else {
            Self::Surface(error)
        }
    }
}

/// End of a one-dimensional branch, in increasing axial-coordinate order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    /// Smallest axial coordinate.
    Left,
    /// Largest axial coordinate.
    Right,
}

/// Actual accepted pressure-time integral on one rigid fluid boundary.
/// Fields are private so a transported momentum flux cannot masquerade as
/// reflected-wall pressure. Exterior pressure is absolute, held over the SAME
/// accepted interval; it is never subtracted from the gas owner's flux.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallPressureImpulse {
    duration_s: f64,
    end_time_s: f64,
    area_m2: f64,
    normal: [f64; 3],
    gas_pa_s: f64,
    exterior_pa_s: f64,
}
impl WallPressureImpulse {
    /// Consume the latest accepted frozen-gas report at a reflecting end.
    /// `axis` is the unit vector of increasing branch coordinate in the same
    /// global frame as the structural surface. Both SSP RK stages contribute.
    /// The caller supplies its latest report; current extensive state is
    /// checked, but these public reports are not authenticated history handles.
    pub fn from_euler1d(
        flow: &Euler1d,
        report: &StepReport,
        end: End,
        axis: [f64; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let boundary = match (flow.boundaries(), end) {
            (Boundaries::Ends { left, .. }, End::Left) => left,
            (Boundaries::Ends { right, .. }, End::Right) => right,
            (Boundaries::Periodic, _) => {
                return Err(Error::Invalid("periodic faces have no wall load"));
            }
        };
        if boundary != Boundary::Reflecting || report.after != flow.totals() {
            return Err(Error::Invalid(
                "latest accepted reflecting-wall report required",
            ));
        }
        let (change, face, sign) = match end {
            End::Left => (report.left_boundary_change, 0, -1.0),
            End::Right => (report.right_boundary_change, flow.mesh().cell_count(), 1.0),
        };
        if change.mass != 0.0 || change.energy != 0.0 {
            return Err(Error::Invalid(
                "stationary adiabatic wall must transfer only momentum",
            ));
        }
        Self::new(
            report.elapsed.value(),
            flow.time().value(),
            flow.mesh().face_areas()[face],
            axis.map(|value| sign * value),
            -sign * change.momentum,
            exterior_pressure_pa,
        )
    }

    /// Same accepted wall account for every species in the original mixture.
    /// No constant-gamma substitution or separate species/energy transport is
    /// introduced. The composition-aware wall solver already owns reflection.
    pub fn from_mixture_euler1d(
        flow: &fs_gas::multicomponent::MixtureEuler1d,
        report: &fs_gas::multicomponent::StepReport,
        end: End,
        axis: [f64; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let boundary = match (flow.boundaries(), end) {
            (Boundaries::Ends { left, .. }, End::Left) => left,
            (Boundaries::Ends { right, .. }, End::Right) => right,
            (Boundaries::Periodic, _) => {
                return Err(Error::Invalid("periodic faces have no wall load"));
            }
        };
        if boundary != Boundary::Reflecting || report.after != flow.totals() {
            return Err(Error::Invalid(
                "latest accepted reflecting-wall mixture report required",
            ));
        }
        let (change, face, sign) = match end {
            End::Left => (report.left_boundary_change, 0, -1.0),
            End::Right => (report.right_boundary_change, flow.mesh().cell_count(), 1.0),
        };
        if change.species_masses.iter().any(|mass| *mass != 0.0) || change.energy != 0.0 {
            return Err(Error::Invalid(
                "stationary adiabatic wall must transfer only momentum",
            ));
        }
        Self::new(
            report.elapsed.value(),
            flow.time().value(),
            flow.mesh().face_areas()[face],
            axis.map(|value| sign * value),
            -sign * change.momentum,
            exterior_pressure_pa,
        )
    }

    fn new(
        duration_s: f64,
        end_time_s: f64,
        area_m2: f64,
        normal: [f64; 3],
        gas_impulse_n_s: f64,
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        if !positive(duration_s)
            || !end_time_s.is_finite()
            || end_time_s < duration_s
            || !positive(area_m2)
            || !positive(gas_impulse_n_s)
            || !unit(normal)
            || !exterior_pressure_pa.is_finite()
            || exterior_pressure_pa < 0.0
        {
            return Err(Error::Invalid(
                "finite positive wall interval, area, pressure and unit frame required",
            ));
        }
        let gas_pa_s = gas_impulse_n_s / area_m2;
        let exterior_pa_s = exterior_pressure_pa * duration_s;
        if !positive(gas_pa_s)
            || !exterior_pa_s.is_finite()
            || (exterior_pressure_pa > 0.0 && exterior_pa_s == 0.0)
            || !(gas_pa_s - exterior_pa_s).is_finite()
        {
            return Err(Error::Invalid(
                "wall pressure-time integral is unrepresentable",
            ));
        }
        Ok(Self {
            duration_s,
            end_time_s,
            area_m2,
            normal,
            gas_pa_s,
            exterior_pa_s,
        })
    }

    /// Actual accepted duration [s].
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Accepted gas clock at the end of the interval [s].
    pub const fn end_time_s(&self) -> f64 {
        self.end_time_s
    }
    /// Physical gas face area [m²].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Positive unit normal from the gas into the support.
    pub const fn outward_normal(&self) -> [f64; 3] {
        self.normal
    }
    /// Original wall pressure integrated over both accepted RK stages [Pa s].
    pub const fn gas_pressure_impulse_pa_s(&self) -> f64 {
        self.gas_pa_s
    }
    /// Opposing held exterior absolute pressure times accepted duration [Pa s].
    pub const fn exterior_pressure_impulse_pa_s(&self) -> f64 {
        self.exterior_pa_s
    }
    /// Signed gas-minus-exterior pressure-time integral [Pa s].
    pub fn net_pressure_impulse_pa_s(&self) -> f64 {
        self.gas_pa_s - self.exterior_pa_s
    }
    /// Actual interval-mean reflected gas pressure [Pa].
    pub fn mean_gas_pressure_pa(&self) -> f64 {
        self.gas_pa_s / self.duration_s
    }
}

/// Cold preparation limits include the original global indexing and all dense
/// pressure rows visited before their sparse compression. No mesh is truncated.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceBudget {
    /// Original shell geometry limits, passed unchanged to its pressure owner.
    pub geometry: ProfileBudget,
    /// Maximum separately loaded fluid faces.
    pub maximum_patches: usize,
    /// Maximum dense-row and selected-geometry scalar visits during preparation.
    pub maximum_preparation_terms: usize,
}

#[derive(Clone, Debug)]
struct Patch {
    area_m2: f64,
    normal: [f64; 3],
    row: Vec<(usize, [f64; 3])>,
}

/// One planar oriented structural patch per fluid face, preserving original
/// node indices. Selected triangles may not be repeated across fluid faces.
#[derive(Clone, Debug)]
pub struct ShockLoadSurface {
    positions: Vec<[f64; 3]>,
    patches: Vec<Patch>,
}
impl ShockLoadSurface {
    /// Prepare the existing P1 pressure rows, then retain only nonzero nodes.
    /// Selected normals point OUT of the fluid; facet side is explicit. A patch
    /// must be planar with one consistent normal. Its area must later match
    /// the gas face to roundoff, so mesh mismatch is never hidden by rescaling.
    pub fn new(
        mesh: &ShellMesh,
        selections: Vec<Vec<OrientedShellFacet>>,
        budget: SurfaceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        let terms = mesh
            .nodes
            .len()
            .checked_mul(6)
            .and_then(|width| width.checked_mul(selections.len()))
            .and_then(|rows| {
                mesh.tris
                    .len()
                    .checked_mul(12)
                    .and_then(|faces| rows.checked_add(faces))
            });
        if selections.is_empty()
            || selections.len() > budget.maximum_patches
            || terms.is_none_or(|n| n > budget.maximum_preparation_terms)
        {
            return Err(Error::Invalid(
                "nonempty pressure patches exceed the preparation allowance",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut patches = Vec::with_capacity(selections.len());
        for selected in selections {
            poll(&mut cancelled)?;
            if selected.is_empty() {
                return Err(Error::Invalid("empty fluid-wall patch"));
            }
            let first = mesh.facet(selected[0].triangle)?;
            let sign = orientation(selected[0].side);
            let normal = first.frame[2].map(|x| sign * x);
            let base = mesh.nodes[mesh.tris[selected[0].triangle][0]];
            for selection in &selected {
                poll(&mut cancelled)?;
                if !seen.insert(selection.triangle) {
                    return Err(Error::Invalid(
                        "an original structural facet belongs to two fluid faces",
                    ));
                }
                let facet = mesh.facet(selection.triangle)?;
                let n = facet.frame[2].map(|x| orientation(selection.side) * x);
                if !same_normal(normal, n) {
                    return Err(Error::Invalid(
                        "fluid face requires one oriented planar patch",
                    ));
                }
                for &node in &mesh.tris[selection.triangle] {
                    let delta = core::array::from_fn(|i| mesh.nodes[node][i] - base[i]);
                    let scale = delta.iter().fold(0.0_f64, |m, value| m.max(value.abs()));
                    if dot(delta, normal).abs() > ROUNDING * scale {
                        return Err(Error::Invalid("fluid-wall patch is not planar"));
                    }
                }
            }
            let port = ShellPressurePort::new(mesh, selected, budget.geometry, &mut cancelled)?;
            let mut row = Vec::new();
            for (node, weights) in port.physical_row_m2().chunks_exact(6).enumerate() {
                if node % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                let area = [weights[0], weights[1], weights[2]];
                if area != [0.0; 3] {
                    row.push((node, area));
                }
            }
            patches.push(Patch {
                area_m2: port.area_m2(),
                normal,
                row,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            positions: mesh.nodes.clone(),
            patches,
        })
    }

    /// Original full structural node count, including unloaded nodes.
    pub fn node_count(&self) -> usize {
        self.positions.len()
    }
    /// Number of separately matched fluid-wall pressure patches.
    pub fn patch_count(&self) -> usize {
        self.patches.len()
    }

    /// Project one simultaneous accepted impulse per original fluid face.
    /// The exterior pressure acts on the opposite side. The resultant and
    /// torque are independently summed from the returned original-node loads.
    /// Supplied `origin_m` changes only the reported moment reference.
    pub fn project(
        &self,
        impulses: &[WallPressureImpulse],
        origin_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShockLoad, Error> {
        poll(&mut cancelled)?;
        if impulses.len() != self.patches.len() || origin_m.iter().any(|x| !x.is_finite()) {
            return Err(Error::Invalid(
                "one wall impulse per patch and finite moment origin required",
            ));
        }
        let first = impulses[0];
        let mut nodal_impulses = vec![[0.0; 3]; self.positions.len()];
        let mut gas_impulse = [0.0; 3];
        let mut exterior_impulse = [0.0; 3];
        for (patch, impulse) in self.patches.iter().zip(impulses) {
            poll(&mut cancelled)?;
            if impulse.duration_s != first.duration_s
                || impulse.end_time_s != first.end_time_s
                || !same_normal(patch.normal, impulse.normal)
                || (patch.area_m2 - impulse.area_m2).abs()
                    > ROUNDING * patch.area_m2.max(impulse.area_m2)
            {
                return Err(Error::Invalid(
                    "wall/surface area, outward normal or accepted interval mismatch",
                ));
            }
            for axis in 0..3 {
                gas_impulse[axis] += impulse.normal[axis] * impulse.gas_pa_s * impulse.area_m2;
                exterior_impulse[axis] -=
                    impulse.normal[axis] * impulse.exterior_pa_s * impulse.area_m2;
            }
            for (i, &(node, row)) in patch.row.iter().enumerate() {
                if i % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                for axis in 0..3 {
                    nodal_impulses[node][axis] += row[axis] * impulse.net_pressure_impulse_pa_s();
                }
            }
        }
        let mut forces = Vec::with_capacity(nodal_impulses.len());
        let mut resultant = [0.0; 3];
        let mut moment = [0.0; 3];
        for (i, (position, impulse)) in self.positions.iter().zip(&nodal_impulses).enumerate() {
            if i % 256 == 0 {
                poll(&mut cancelled)?;
            }
            let lever = core::array::from_fn(|axis| position[axis] - origin_m[axis]);
            let torque = cross(lever, *impulse);
            for axis in 0..3 {
                resultant[axis] += impulse[axis];
                moment[axis] += torque[axis];
            }
            forces.push(impulse.map(|x| x / first.duration_s));
        }
        if forces
            .iter()
            .flatten()
            .chain(nodal_impulses.iter().flatten())
            .chain(&resultant)
            .chain(&moment)
            .chain(&gas_impulse)
            .chain(&exterior_impulse)
            .any(|x| !x.is_finite())
        {
            return Err(Error::Invalid(
                "structural impulse, load or moment overflow",
            ));
        }
        poll(&mut cancelled)?;
        Ok(ShockLoad {
            duration_s: first.duration_s,
            end_time_s: first.end_time_s,
            nodal_impulses_n_s: nodal_impulses,
            mean_nodal_forces_n: forces,
            gas_impulse_n_s: gas_impulse,
            exterior_impulse_n_s: exterior_impulse,
            resultant_impulse_n_s: resultant,
            angular_impulse_n_m_s: moment,
        })
    }
}

/// Original structural node order, including supported nodes. A constraint
/// owner may retain their reactions; this adapter never redistributes them.
#[derive(Clone, Debug, PartialEq)]
pub struct ShockLoad {
    /// Common accepted physical duration [s].
    pub duration_s: f64,
    /// Common accepted gas endpoint clock [s].
    pub end_time_s: f64,
    /// Integrated force on each original structural node [N s].
    pub nodal_impulses_n_s: Vec<[f64; 3]>,
    /// Corresponding interval-mean nodal forces [N], with no constraint projection.
    pub mean_nodal_forces_n: Vec<[f64; 3]>,
    /// Gas-side support impulse before exterior pressure subtraction [N s].
    pub gas_impulse_n_s: [f64; 3],
    /// Opposite impulse from the explicitly held exterior pressure [N s].
    pub exterior_impulse_n_s: [f64; 3],
    /// Sum of the returned structural nodal impulses [N s].
    pub resultant_impulse_n_s: [f64; 3],
    /// Moment of returned nodal impulses about the requested origin [N m s].
    pub angular_impulse_n_m_s: [f64; 3],
}
impl ShockLoad {
    /// Spatial virtual work of the interval-mean traction on original nodal
    /// displacements. This is not a claim about unresolved pressure-time work.
    pub fn virtual_work_j(&self, displacement_m: &[[f64; 3]]) -> Result<f64, Error> {
        if displacement_m.len() != self.mean_nodal_forces_n.len()
            || displacement_m.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(Error::Invalid(
                "finite displacement at each original structural node required",
            ));
        }
        let work: f64 = self
            .mean_nodal_forces_n
            .iter()
            .zip(displacement_m)
            .map(|(force, u)| dot(*force, *u))
            .sum();
        if !work.is_finite() {
            return Err(Error::Invalid("virtual work overflow"));
        }
        Ok(work)
    }
}

const ROUNDING: f64 = 1024.0 * f64::EPSILON;
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn unit(x: [f64; 3]) -> bool {
    x.iter().all(|x| x.is_finite()) && (dot(x, x) - 1.0).abs() <= ROUNDING
}
fn same_normal(a: [f64; 3], b: [f64; 3]) -> bool {
    a.into_iter().zip(b).all(|(a, b)| (a - b).abs() <= ROUNDING)
}
fn orientation(side: ShellFace) -> f64 {
    match side {
        ShellFace::Positive => 1.0,
        ShellFace::Negative => -1.0,
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
