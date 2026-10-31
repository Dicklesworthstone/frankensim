//! Plane-clipped displacement geometry. Material mass geometry is deliberately
//! separate: a sealed hollow body displaces its envelope, not just its walls.
use fs_ga::{Se3, Vec3};
use super::{HydrostaticError, finite, finite_vec, poll};
use crate::VolumeMoments;

/// Maximum tetrahedra in a displacement mesh and one interface-clipping pass.
/// Layered loading can traverse the mesh once per nonzero density jump.
pub const MAX_DISPLACEMENT_CELLS: usize = 65_536;

/// Caller-supplied, non-overlapping tetrahedral partition of the volume which
/// excludes fluid. Sealed voids belong in this volume; flooded voids do not.
/// Numeric admission does not certify mesh coverage, overlap or watertightness.
#[derive(Clone, Debug)]
pub struct DisplacementMesh {
    cells: Vec<[Vec3; 4]>,
    volume_m3: f64,
    centroid_body_m: Vec3,
}
impl DisplacementMesh {
    /// Admit an immutable displacement envelope in the original body frame.
    /// The original volume-moment owner validates each nondegenerate cell.
    /// Vertex orientation does not select a sign for displaced volume.
    pub fn new(cells: Vec<[Vec3; 4]>, max_cells: usize,
        mut cancelled: impl FnMut() -> bool) -> Result<Self, HydrostaticError> {
        poll(&mut cancelled)?;
        if max_cells == 0 || max_cells > MAX_DISPLACEMENT_CELLS
            || cells.is_empty() || cells.len() > max_cells {
            return Err(HydrostaticError::Invalid("displacement cell budget exceeded or empty mesh"));
        }
        let mut volume_m3 = 0.0;
        let mut centroid_body_m = Vec3::new(0.0, 0.0, 0.0);
        for cell in &cells {
            poll(&mut cancelled)?;
            let moments = VolumeMoments::tetrahedron(*cell).map_err(HydrostaticError::Geometry)?;
            volume_m3 = finite(volume_m3 + moments.volume_m3())?;
            centroid_body_m = finite_vec(centroid_body_m
                + (moments.centroid_m()-centroid_body_m).scale(moments.volume_m3()/volume_m3))?;
        }
        poll(&mut cancelled)?;
        Ok(Self { cells, volume_m3, centroid_body_m })
    }

    /// A centered cuboid envelope, partitioned into six straight tetrahedra.
    /// Side lengths are full lengths [m]. No density or material is selected.
    pub fn cuboid(sides: Vec3, cancelled: impl FnMut() -> bool) -> Result<Self, HydrostaticError> {
        VolumeMoments::cuboid(sides).map_err(HydrostaticError::Geometry)?;
        let start = sides.scale(-0.5);
        let axes = [Vec3::new(sides.x,0.0,0.0), Vec3::new(0.0,sides.y,0.0), Vec3::new(0.0,0.0,sides.z)];
        let cells = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]
            .map(|p| [start, start+axes[p[0]], start+axes[p[0]]+axes[p[1]], sides.scale(0.5)]);
        Self::new(cells.to_vec(), 6, cancelled)
    }

    /// Complete undeformed fluid-excluding volume [m3], not material volume.
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    /// Original body-frame cells in the deterministic integration order.
    pub fn cells(&self) -> &[[Vec3; 4]] { &self.cells }

    pub(super) fn submerged(&self, world_from_body: Se3, down: Vec3,
        surface: Vec3, cancelled: &mut impl FnMut() -> bool) -> Result<DisplacedVolume, HydrostaticError> {
        self.submerged_below(world_from_body, down, surface, 0.0, cancelled)
    }

    pub(super) fn whole(&self, world_from_body: Se3, down: Vec3,
        datum: Vec3) -> Result<DisplacedVolume, HydrostaticError> {
        let direction = world_from_body.rotation().inverse().rotate(down)?;
        let mean_depth = finite(down.dot(world_from_body.translation()-datum)
            + direction.dot(self.centroid_body_m))?;
        Ok(DisplacedVolume { volume_m3: self.volume_m3,
            centroid_body_m: Some(self.centroid_body_m),
            depth_integral_m4: finite(self.volume_m3*mean_depth)? })
    }

    pub(super) fn submerged_below(&self, world_from_body: Se3, down: Vec3,
        datum: Vec3, interface_depth_m: f64,
        cancelled: &mut impl FnMut() -> bool) -> Result<DisplacedVolume, HydrostaticError> {
        poll(cancelled)?;
        // Transform the plane once, not every vertex into a large world origin.
        let direction = world_from_body.rotation().inverse().rotate(down)?;
        let origin_depth = finite(down.dot(world_from_body.translation()-datum)-interface_depth_m)?;
        let mut sum = Accumulator::default();
        for cell in &self.cells {
            poll(cancelled)?;
            let mut vertices = [Vertex { point: Vec3::new(0.0,0.0,0.0), depth: 0.0 };4];
            for (out, &point) in vertices.iter_mut().zip(cell) {
                *out = Vertex { point, depth: finite(origin_depth+direction.dot(point))? };
            }
            clip(vertices, &mut sum)?;
        }
        let depth_integral_m4 = finite(sum.volume * sum.mean_depth)?;
        if sum.volume > 0.0 && depth_integral_m4 <= 0.0 {
            return Err(HydrostaticError::Invalid("submerged depth integral is unrepresentable"));
        }
        poll(cancelled)?;
        Ok(DisplacedVolume { volume_m3: sum.volume,
            centroid_body_m: (sum.volume > 0.0).then_some(sum.center), depth_integral_m4 })
    }
}

/// Analytically integrated clipped volume and first moments. Floating-point
/// arithmetic is measured numerics, not a certified geometric enclosure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplacedVolume {
    volume_m3: f64,
    centroid_body_m: Option<Vec3>,
    pub(super) depth_integral_m4: f64,
}
impl DisplacedVolume {
    pub(super) const fn dry() -> Self {
        Self { volume_m3: 0.0, centroid_body_m: None, depth_integral_m4: 0.0 }
    }
    /// Immersed volume [m3]; exactly zero for a dry envelope.
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    /// Geometric centroid in the ORIGINAL geometry frame; None when dry.
    /// Equals the center of buoyancy only for a homogeneous fluid density.
    pub const fn centroid_body_m(&self) -> Option<Vec3> { self.centroid_body_m }
}

#[derive(Clone, Copy)]
struct Vertex { point: Vec3, depth: f64 }
impl Vertex {
    fn intersection(self, dry: Self) -> Result<Self, HydrostaticError> {
        // Opposite signs, with a possible exactly-on-plane dry endpoint.
        let scale = self.depth.abs().max(dry.depth.abs());
        let fraction = (self.depth/scale) / (self.depth/scale-dry.depth/scale);
        let point = finite_vec(self.point+(dry.point-self.point).scale(fraction))?;
        Ok(Self { point, depth: 0.0 })
    }
}

struct Accumulator { volume: f64, center: Vec3, mean_depth: f64 }
impl Default for Accumulator {
    fn default() -> Self { Self { volume:0.0, center:Vec3::new(0.0,0.0,0.0), mean_depth:0.0 } }
}
impl Accumulator {
    fn add(&mut self, vertices: [Vertex;4]) -> Result<(), HydrostaticError> {
        let p = vertices.map(|v| v.point);
        let offset = p.map(|v| v-p[0]);
        let volume = finite(offset[1].dot(offset[2].cross(offset[3])).abs()/6.0)?;
        // On-plane vertices can make some pieces exactly degenerate. They
        // contribute zero; no epsilon discards a small positive wet volume.
        if volume == 0.0 { return Ok(()); }
        let center = finite_vec(p[0]+(offset[1]+offset[2]+offset[3]).scale(0.25))?;
        let mean_depth = finite(vertices.iter().map(|v| 0.25*v.depth).sum())?;
        let next = finite(self.volume+volume)?;
        let weight = volume/next;
        if weight <= 0.0 || mean_depth < 0.0 {
            return Err(HydrostaticError::Invalid("clipped volume moments are unrepresentable"));
        }
        self.center = finite_vec(self.center+(center-self.center).scale(weight))?;
        self.mean_depth = finite(self.mean_depth+(mean_depth-self.mean_depth)*weight)?;
        self.volume = next;
        Ok(())
    }
}

fn clip(v: [Vertex;4], sum: &mut Accumulator) -> Result<(), HydrostaticError> {
    let mut wet = [0usize;4]; let mut dry = [0usize;4]; let mut n = 0; let mut m = 0;
    for (i, vertex) in v.iter().enumerate() {
        if vertex.depth > 0.0 { wet[n]=i; n+=1; } else { dry[m]=i; m+=1; }
    }
    match n {
        0 => Ok(()),
        4 => sum.add(v),
        1 => {
            let a=v[wet[0]];
            sum.add([a,a.intersection(v[dry[0]])?,a.intersection(v[dry[1]])?,a.intersection(v[dry[2]])?])
        }
        2 => {
            let (a,b,c,d)=(v[wet[0]],v[wet[1]],v[dry[0]],v[dry[1]]);
            let (ac,ad,bc,bd)=(a.intersection(c)?,a.intersection(d)?,b.intersection(c)?,b.intersection(d)?);
            sum.add([a,b,ac,ad])?; sum.add([b,ac,ad,bc])?; sum.add([b,ad,bc,bd])
        }
        _ => {
            let (a,b,c,d)=(v[wet[0]],v[wet[1]],v[wet[2]],v[dry[0]]);
            let (ad,bd,cd)=(a.intersection(d)?,b.intersection(d)?,c.intersection(d)?);
            sum.add([a,b,c,ad])?; sum.add([b,c,ad,bd])?; sum.add([c,ad,bd,cd])
        }
    }
}
