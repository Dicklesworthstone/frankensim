//! Original JKR/Maugis pressure on actual exterior P1 tetrahedral faces.
//! Mixed-mode JKR profiles retain their actual shear-reduced footprint and
//! normal opening. Their normal pressure follows this same projection path;
//! the original mixed-mode owner separately retains tangential efforts and
//! storage. Prescribing that field is not a self-consistent finite-body
//! normal/shear contact solve.
//!
//! Each polar cell belongs wholly to one coplanar exterior face. Its three
//! barycentric functions come from the ORIGINAL mesh, not caller-authored
//! weights. Integrating N_i p uses the original force and first moments; the
//! JKR edge singularity is never sampled or clipped. Positive source pressure
//! acts inward, so its tensile contribution pulls the receiver outward.
//!
//! The fixed grid must fit the selected facets, including every curved arc.
//! This supports an interior patch or a fan of faces around a contact centre;
//! it does not silently clip arbitrary triangle/polar intersections. Loads are
//! prescribed fields on reference geometry, not an implicit adhesion solve or
//! a follower-pressure law. The existing native solver owns deformation and
//! any explicitly prescribed support motion. Pressure work uses that same
//! full nodal motion; support-reaction work remains in its own native account.

mod transient;
pub use transient::{AdhesiveP1Step3, AdhesiveP1SupportEndpoint3, AdhesiveP1Transient3};

use super::{LoadError3, TetLoadCase3, TetLoadGeometry3, cross, dot, norm, sub};
use core::f64::consts::{PI, TAU};
use fs_contact::normal_patch::traction::{
    AdhesiveTractionProfile,
    sectors::{
        MAX_POLAR_TRACTION_SECTORS, MAX_RADIAL_MOMENT_EVALUATIONS, PolarRingScratch,
        PolarTractionTransfer, PolarTransferBudget, SectorAdhesiveTraction, SectorTractionError,
    },
};
use fs_exec::{CancelGate, Cx};
use fs_math::dd::Dd;

/// Fixed SI contact frame. The two unit axes define a right-handed inward
/// normal and the SAME local x,y coordinates used by the original profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveP1Frame3 {
    pub centre_m: [f64; 3],
    pub x_axis: [f64; 3],
    pub y_axis: [f64; 3],
}

/// Complete radial-major/angular-minor polar partition. Each face index names
/// the canonical exterior-face order of the original TetLoadGeometry3.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveP1Grid3 {
    pub radii_m: Vec<f64>,
    pub angles_rad: Vec<f64>,
    pub faces: Vec<usize>,
}

/// Explicit geometry, storage and small-work gates. Moment quadrature carries
/// an ERROR ESTIMATE, not a certified enclosure or material uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveP1Controls3 {
    pub maximum_nodes: usize,
    pub geometry_tolerance_m: f64,
    pub maximum_frame_error: f64,
    pub transfer: PolarTransferBudget,
    pub projection_force_tolerance_n: f64,
    pub projection_moment_tolerance_n_m: f64,
    pub work_absolute_tolerance_j: f64,
    pub work_relative_tolerance: f64,
    pub maximum_quadrature_work_estimate_j: f64,
}

#[derive(Debug)]
pub enum AdhesiveP1Error3 {
    Invalid(&'static str),
    Loads(LoadError3),
    Spatial(SectorTractionError),
    Native(super::transient::TransientError3),
    /// One entire sector leaves its selected original exterior triangle.
    Coverage {
        cell: usize,
        face: usize,
    },
    Projection {
        force_residual_n: [f64; 3],
        moment_residual_n_m: [f64; 3],
    },
    Work {
        residual_j: f64,
        tolerance_j: f64,
    },
    QuadratureWork {
        estimate_j: f64,
        maximum_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for AdhesiveP1Error3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native P1 adhesive pressure refused: {self:?}")
    }
}
impl std::error::Error for AdhesiveP1Error3 {}
impl From<LoadError3> for AdhesiveP1Error3 {
    fn from(value: LoadError3) -> Self {
        Self::Loads(value)
    }
}
impl From<SectorTractionError> for AdhesiveP1Error3 {
    fn from(value: SectorTractionError) -> Self {
        Self::Spatial(value)
    }
}
impl From<super::transient::TransientError3> for AdhesiveP1Error3 {
    fn from(value: super::transient::TransientError3) -> Self {
        Self::Native(value)
    }
}

/// Original source and independently recovered nodal resultant, after adding
/// other external forces. Supported-node efforts are included in this account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveP1Transfer3 {
    pub source: PolarTractionTransfer,
    pub force_world_n: [f64; 3],
    pub moment_about_centre_n_m: [f64; 3],
    pub force_residual_n: [f64; 3],
    pub moment_residual_n_m: [f64; 3],
}

#[derive(Clone, Copy)]
struct Trace {
    nodes: [usize; 3],
    // N_i(x,y) = coefficient[i][0] + coefficient[i][1] x + coefficient[i][2] y.
    coefficient: [[f64; 3]; 3],
}

/// Immutable original geometry/traces with bounded reusable projection scratch.
/// A native load case owns only its canonical sparse nodal-force output.
pub struct AdhesiveP1Patch3<'g, 'm> {
    geometry: &'g TetLoadGeometry3<'m>,
    frame: AdhesiveP1Frame3,
    normal: [f64; 3],
    grid: AdhesiveP1Grid3,
    traces: Vec<Trace>,
    controls: AdhesiveP1Controls3,
    cells: Vec<SectorAdhesiveTraction>,
    rings: Vec<PolarRingScratch>,
    sums: Vec<Dd>,
    represented: Vec<Dd>,
    background: Vec<f64>,
}

impl<'g, 'm> AdhesiveP1Patch3<'g, 'm> {
    /// Admit the complete pressure support against the true exterior mesh.
    /// Bounds are analytic extrema of each affine barycentric function over
    /// the entire radial/angular cell. No corner-only or centroid-only test.
    pub fn new(
        geometry: &'g TetLoadGeometry3<'m>,
        frame: AdhesiveP1Frame3,
        grid: AdhesiveP1Grid3,
        controls: AdhesiveP1Controls3,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<Self, AdhesiveP1Error3> {
        poll(gate, cx)?;
        validate(controls)?;
        let nodes = geometry.nodes_m().len();
        let dofs = nodes
            .checked_mul(3)
            .ok_or(AdhesiveP1Error3::Invalid("nodal extent overflow"))?;
        if nodes > controls.maximum_nodes {
            return Err(AdhesiveP1Error3::Invalid(
                "original mesh exceeds node budget",
            ));
        }
        for value in frame
            .centre_m
            .into_iter()
            .chain(frame.x_axis)
            .chain(frame.y_axis)
        {
            finite(value)?;
        }
        let normal = cross(frame.x_axis, frame.y_axis);
        if (norm(frame.x_axis) - 1.0).abs() > controls.maximum_frame_error
            || (norm(frame.y_axis) - 1.0).abs() > controls.maximum_frame_error
            || dot(frame.x_axis, frame.y_axis).abs() > controls.maximum_frame_error
            || (norm(normal) - 1.0).abs() > controls.maximum_frame_error
        {
            return Err(AdhesiveP1Error3::Invalid(
                "unit orthogonal contact axes required",
            ));
        }
        let rings = grid
            .radii_m
            .len()
            .checked_sub(1)
            .ok_or(AdhesiveP1Error3::Invalid("empty polar grid"))?;
        let angles = grid
            .angles_rad
            .len()
            .checked_sub(1)
            .ok_or(AdhesiveP1Error3::Invalid("empty polar grid"))?;
        let count = rings
            .checked_mul(angles)
            .ok_or(AdhesiveP1Error3::Invalid("polar extent overflow"))?;
        if rings == 0
            || angles == 0
            || count > MAX_POLAR_TRACTION_SECTORS
            || grid.faces.len() != count
            || grid.radii_m[0] != 0.0
            || grid.angles_rad[0] != 0.0
            || *grid.angles_rad.last().unwrap() != TAU
        {
            return Err(AdhesiveP1Error3::Invalid(
                "complete bounded polar partition required",
            ));
        }
        for edges in grid.radii_m.windows(2).chain(grid.angles_rad.windows(2)) {
            poll(gate, cx)?;
            if !edges[0].is_finite() || !edges[1].is_finite() || edges[1] <= edges[0] {
                return Err(AdhesiveP1Error3::Invalid(
                    "finite increasing polar edges required",
                ));
            }
        }
        let mut traces = Vec::with_capacity(count);
        for (cell, &id) in grid.faces.iter().enumerate() {
            poll(gate, cx)?;
            let face = geometry
                .faces()
                .get(id)
                .ok_or(AdhesiveP1Error3::Invalid("unknown exterior face"))?;
            if norm(std::array::from_fn(|i| normal[i] + face.outward_normal[i]))
                > controls.maximum_frame_error
            {
                return Err(AdhesiveP1Error3::Invalid(
                    "pressure normal must point into each selected face",
                ));
            }
            let mut p = [[0.0; 2]; 3];
            for (i, &node) in face.nodes.iter().enumerate() {
                let d = sub(geometry.nodes_m()[node], frame.centre_m);
                if finite(dot(d, normal))?.abs() > controls.geometry_tolerance_m {
                    return Err(AdhesiveP1Error3::Invalid(
                        "patch and selected faces must be coplanar",
                    ));
                }
                p[i] = [finite(dot(d, frame.x_axis))?, finite(dot(d, frame.y_axis))?];
            }
            let coefficient = barycentric(p)?;
            let radial = cell / angles;
            let angular = cell % angles;
            for [c, x, y] in coefficient {
                let minimum = sector_minimum(
                    [c, x, y],
                    [grid.radii_m[radial], grid.radii_m[radial + 1]],
                    [grid.angles_rad[angular], grid.angles_rad[angular + 1]],
                )?;
                let scale = finite(c.abs() + grid.radii_m[radial + 1] * x.hypot(y))?;
                let allowance = finite(
                    controls.geometry_tolerance_m * x.hypot(y) + 64.0 * f64::EPSILON * scale,
                )?;
                if minimum < -allowance {
                    return Err(AdhesiveP1Error3::Coverage { cell, face: id });
                }
            }
            traces.push(Trace {
                nodes: face.nodes,
                coefficient,
            });
        }
        poll(gate, cx)?;
        Ok(Self {
            geometry,
            frame,
            normal,
            grid,
            traces,
            controls,
            cells: vec![SectorAdhesiveTraction::default(); count],
            rings: vec![PolarRingScratch::default(); rings],
            sums: vec![Dd::ZERO; dofs],
            represented: vec![Dd::ZERO; dofs],
            background: vec![0.0; dofs],
        })
    }

    /// Integrate the original immutable source and superpose an existing load
    /// case on the SAME reference mesh. The returned native load preserves
    /// constrained efforts, so static and transient reaction owners see them.
    /// Scratch may change on refusal; no caller-owned accepted state does.
    pub fn assemble(
        &mut self,
        profile: &AdhesiveTractionProfile,
        other: &TetLoadCase3<'_, '_>,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<(TetLoadCase3<'g, 'm>, AdhesiveP1Transfer3), AdhesiveP1Error3> {
        poll(gate, cx)?;
        if other.geometry.nodes != self.geometry.nodes || other.geometry.tets != self.geometry.tets
        {
            return Err(AdhesiveP1Error3::Invalid(
                "background loads belong to another reference mesh",
            ));
        }
        let source = profile.transfer_polar_grid(
            &self.grid.radii_m,
            &self.grid.angles_rad,
            &mut self.cells,
            &mut self.rings,
            self.controls.transfer,
            gate,
        )?;
        self.sums.fill(Dd::ZERO);
        self.background.fill(0.0);
        for (&Trace { nodes, coefficient }, cell) in self.traces.iter().zip(&self.cells) {
            poll(gate, cx)?;
            for (node, b) in nodes.into_iter().zip(coefficient) {
                let value = product(b[0], cell.normal_force_n)
                    + product(b[1], cell.first_moment_n_m[0])
                    + product(b[2], cell.first_moment_n_m[1]);
                for (axis, &normal) in self.normal.iter().enumerate() {
                    let d = 3 * node + axis;
                    self.sums[d] = self.sums[d] + value * Dd::from_f64(normal);
                }
            }
        }
        for &(d, force) in other.nodal_forces_n() {
            poll(gate, cx)?;
            self.background[d] = force;
        }
        let mut nodal_forces_n = Vec::new();
        let mut force = [Dd::ZERO; 3];
        let mut moment = [Dd::ZERO; 3];
        for node in 0..self.geometry.nodes.len() {
            poll(gate, cx)?;
            let r = sub(self.geometry.nodes[node], self.frame.centre_m);
            let mut actual = [Dd::ZERO; 3];
            for axis in 0..3 {
                let d = 3 * node + axis;
                let applied = finite(self.background[d] + finite(self.sums[d].to_f64())?)?;
                if applied != 0.0 {
                    nodal_forces_n.push((d, applied));
                }
                // Charge only the correction represented after the final f64
                // addition. Large unrelated loads cannot hide erased adhesion.
                actual[axis] = Dd::from_f64(applied) - Dd::from_f64(self.background[d]);
                self.represented[d] = actual[axis];
                force[axis] = force[axis] + actual[axis];
            }
            for axis in 0..3 {
                let j = (axis + 1) % 3;
                let k = (axis + 2) % 3;
                moment[axis] =
                    moment[axis] + Dd::from_f64(r[j]) * actual[k] - Dd::from_f64(r[k]) * actual[j];
            }
        }
        let residual = std::array::from_fn(|i| {
            (force[i] - product(profile.normal_force_n(), self.normal[i])).to_f64()
        });
        let moment = moment.map(Dd::to_f64);
        for value in residual.into_iter().chain(moment) {
            finite(value)?;
        }
        if norm(residual) > self.controls.projection_force_tolerance_n
            || norm(moment) > self.controls.projection_moment_tolerance_n_m
        {
            return Err(AdhesiveP1Error3::Projection {
                force_residual_n: residual,
                moment_residual_n_m: moment,
            });
        }
        let transfer = AdhesiveP1Transfer3 {
            source,
            force_world_n: force.map(Dd::to_f64),
            moment_about_centre_n_m: moment,
            force_residual_n: residual,
            moment_residual_n_m: moment,
        };
        poll(gate, cx)?;
        Ok((
            TetLoadCase3 {
                geometry: self.geometry,
                nodal_forces_n,
            },
            transfer,
        ))
    }
}

fn barycentric(p: [[f64; 2]; 3]) -> Result<[[f64; 3]; 3], AdhesiveP1Error3> {
    let a = [p[1][0] - p[0][0], p[1][1] - p[0][1]];
    let b = [p[2][0] - p[0][0], p[2][1] - p[0][1]];
    let h = finite(a[0].hypot(a[1]).max(b[0].hypot(b[1])))?;
    if h == 0.0 {
        return Err(AdhesiveP1Error3::Invalid("degenerate face trace"));
    }
    let a = a.map(|v| v / h);
    let b = b.map(|v| v / h);
    let det = finite(a[0] * b[1] - a[1] * b[0])?;
    if det == 0.0 {
        return Err(AdhesiveP1Error3::Invalid("degenerate projected face"));
    }
    let mut c = [[0.0; 3]; 3];
    c[1][1] = (b[1] / det) / h;
    c[1][2] = (-b[0] / det) / h;
    c[2][1] = (-a[1] / det) / h;
    c[2][2] = (a[0] / det) / h;
    for i in [1, 2] {
        c[i][0] = -c[i][1] * p[0][0] - c[i][2] * p[0][1];
    }
    c[0] = [
        1.0 - c[1][0] - c[2][0],
        -c[1][1] - c[2][1],
        -c[1][2] - c[2][2],
    ];
    for value in c.into_iter().flatten() {
        finite(value)?;
    }
    Ok(c)
}

// A linear barycentric function attains its radial minimum at an endpoint.
// Its angular minimum is an endpoint or the opposite gradient direction.
fn sector_minimum(b: [f64; 3], radii: [f64; 2], angles: [f64; 2]) -> Result<f64, AdhesiveP1Error3> {
    let evaluate = |angle: f64| {
        let (s, c) = angle.sin_cos();
        b[1] * c + b[2] * s
    };
    let mut minimum = evaluate(angles[0]).min(evaluate(angles[1]));
    let direction = (b[2].atan2(b[1]) + PI).rem_euclid(TAU);
    if angles[0] <= direction && direction <= angles[1] {
        minimum = minimum.min(-b[1].hypot(b[2]));
    }
    finite(b[0] + (radii[0] * minimum).min(radii[1] * minimum))
}

fn product(a: f64, b: f64) -> Dd {
    Dd::from_f64(a) * Dd::from_f64(b)
}
fn finite(value: f64) -> Result<f64, AdhesiveP1Error3> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AdhesiveP1Error3::Invalid(
            "nonfinite adhesive projection arithmetic",
        ))
    }
}
fn poll(gate: &CancelGate, cx: &Cx<'_>) -> Result<(), AdhesiveP1Error3> {
    if gate.is_requested() {
        return Err(AdhesiveP1Error3::Cancelled);
    }
    cx.checkpoint().map_err(|_| AdhesiveP1Error3::Cancelled)
}
fn validate(c: AdhesiveP1Controls3) -> Result<(), AdhesiveP1Error3> {
    if c.maximum_nodes == 0
        || c.transfer.radial_moment.maximum_evaluations > MAX_RADIAL_MOMENT_EVALUATIONS
        || [
            c.geometry_tolerance_m,
            c.projection_force_tolerance_n,
            c.projection_moment_tolerance_n_m,
            c.work_absolute_tolerance_j,
            c.maximum_quadrature_work_estimate_j,
            c.transfer.force_absolute_tolerance_n,
            c.transfer.moment_absolute_tolerance_n_m,
            c.transfer.radial_moment.absolute_tolerance_n_m,
            c.transfer.radial_moment.relative_tolerance,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        || !c.maximum_frame_error.is_finite()
        || !(0.0..1e-4).contains(&c.maximum_frame_error)
        || !c.work_relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&c.work_relative_tolerance)
    {
        return Err(AdhesiveP1Error3::Invalid(
            "bounded mesh and finite geometry/accuracy allowances required",
        ));
    }
    Ok(())
}
