//! A declared elastic barrel, clamped at rigid hoops, sharing the head/air clock.
//! The full linear shell pencil owns its frequencies. Inner pressure and outer
//! radiation use the same finite-thickness translation/rotation kinematics.
use std::io::Read;
use fs_exec::CancelGate;
use fs_plate::shell::{ShellMesh, ShellSupport, assemble_shell, assemble_shell_sections, modes_shell};
use fs_plate::shell::reduction::{ReductionBudget, ShellReduction};
use fs_plate::shell::reduction::radiation::ShellFace;
use fs_couple::render::plate::impact::cavity::cylinder::CylindricalCavity;
use fs_couple::vibroacoustic::{StructuralModes, assemble_coupling};
use super::{BodyPotential, Error, ImpactBody, drum_spec};

pub const HEADER: &str = "frankensim-drum-barrel-v1";
pub const ORTHOTROPIC_HEADER: &str = "frankensim-drum-barrel-v2";

/// Homogeneous through-thickness constitutive data. A laminate's independent
/// membrane, coupling and bending tensors cannot be inferred from these values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Material {
    Isotropic { young_pa: f64, poisson: f64 },
    Orthotropic {
        young_1_pa: f64, young_2_pa: f64, poisson_12: f64, shear_12_pa: f64,
        /// Material axis 1, from the upward barrel axis toward increasing
        /// azimuth, on the actual facet tangent plane. Radians, in [-pi, pi].
        grain_angle_rad: f64,
    },
}
impl Material {
    fn section(self, thickness: f64, density: f64, local_angle: f64)
        -> Result<fs_plate::PlateSection, Error>
    {
        Ok(match self {
            Self::Isotropic { young_pa, poisson } =>
                fs_plate::PlateSection::isotropic(young_pa, poisson, thickness, density)?,
            Self::Orthotropic { young_1_pa, young_2_pa, poisson_12, shear_12_pa, .. } =>
                fs_plate::PlateSection::orthotropic_plane_stress_at_angle(
                    young_1_pa, young_2_pa, poisson_12, shear_12_pa, thickness, density, local_angle)?,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spec {
    pub material: Material,
    pub density_kg_m3: f64,
    pub damping_ratio: f64,
    pub axial_intervals: usize,
    pub band_hz: [f64; 2],
}

pub fn option(args: &mut Vec<String>) -> Result<Option<Spec>, Error> {
    let mut positions = args.iter().enumerate().filter(|(_, a)| a.as_str() == "--elastic-barrel");
    let Some((at, _)) = positions.next() else { return Ok(None); };
    if positions.next().is_some() { return Err("--elastic-barrel may be supplied only once".into()); }
    let path = args.get(at + 1).ok_or("--elastic-barrel needs a supplied barrel specification")?;
    let mut text = String::new();
    std::fs::File::open(path)?.take(16_385).read_to_string(&mut text)?;
    let spec = Spec::read(&text)?;
    args.drain(at..at + 2);
    Ok(Some(spec))
}

pub fn admit_command(spec: Option<&Spec>, command: &str) -> Result<(), Error> {
    if spec.is_some() && !matches!(command, "drum" | "drum-wav" | "drum-mic" |
        "drum-modal" | "drum-modal-wav" | "drum-modal-mic" |
        "drum-stretch" | "drum-stretch-wav" | "drum-stretch-mic" |
        "snare" | "snare-wav" | "snare-mic" | "snare-off" | "snare-off-wav" | "snare-off-mic") {
        return Err("--elastic-barrel requires a drum or snare command".into());
    }
    Ok(())
}

impl Spec {
    pub fn read(text: &str) -> Result<Self, Error> {
        if text.len() > 16_384 { return Err("barrel specification exceeds 16 KiB".into()); }
        let (mut header, mut material, mut mesh, mut band) = (None, None, None, None);
        for (line, raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() { continue; }
            if header.is_none() {
                header = Some(match row {
                    HEADER => false, ORTHOTROPIC_HEADER => true,
                    _ => return Err("expected frankensim-drum-barrel-v1 or frankensim-drum-barrel-v2".into()),
                });
                continue;
            }
            let fields: Vec<_> = row.split(',').map(str::trim).collect();
            match fields.as_slice() {
                ["material", e, nu, rho, zeta] if header == Some(false) && material.is_none() => {
                    material = Some((Material::Isotropic { young_pa: e.parse()?, poisson: nu.parse()? },
                        rho.parse()?, zeta.parse()?));
                }
                ["orthotropic", e1, e2, nu12, g12, rho, zeta, angle]
                    if header == Some(true) && material.is_none() => {
                    material = Some((Material::Orthotropic { young_1_pa: e1.parse()?, young_2_pa: e2.parse()?,
                        poisson_12: nu12.parse()?, shear_12_pa: g12.parse()?, grain_angle_rad: angle.parse()? },
                        rho.parse()?, zeta.parse()?));
                }
                ["axial_intervals", n] if mesh.is_none() => { mesh = Some(n.parse()?); }
                ["band_hz", low, high] if band.is_none() => { band = Some([low.parse()?, high.parse()?]); }
                _ => return Err(format!("barrel specification line {}: unknown, repeated or malformed record", line + 1).into()),
            }
        }
        let (material, density_kg_m3, damping_ratio) = material.ok_or("barrel material record is required")?;
        let spec = Self { material, density_kg_m3, damping_ratio,
            axial_intervals: mesh.ok_or("barrel axial_intervals record is required")?,
            band_hz: band.ok_or("barrel band_hz record is required")? };
        spec.validate()?;
        Ok(spec)
    }

    fn validate(self) -> Result<(), Error> {
        self.material.section(0.001, self.density_kg_m3, 0.0)?;
        if let Material::Orthotropic { grain_angle_rad, .. } = self.material {
            if !grain_angle_rad.is_finite() || grain_angle_rad.abs() > core::f64::consts::PI {
                return Err("barrel grain angle must be finite radians in [-pi, pi]".into());
            }
        }
        if !self.damping_ratio.is_finite() || self.damping_ratio < 0.0
            || !(2..=32).contains(&self.axial_intervals)
            || self.band_hz.iter().any(|v| !v.is_finite()) || self.band_hz[0] < 0.0
            || self.band_hz[1] <= self.band_hz[0]
            || !(core::f64::consts::TAU * self.band_hz[1]).powi(2).is_finite() {
            return Err("barrel needs passive finite damping, 2..32 axial intervals and an increasing finite frequency window".into());
        }
        Ok(())
    }

    /// Rotate the material into EACH actual element frame. The two triangles
    /// on a rectangular barrel strip generally have different first edges;
    /// reusing one local D would rotate the grain at every diagonal.
    fn sections(&self, mesh: &ShellMesh, thickness: f64) -> Result<Vec<fs_plate::PlateSection>, Error> {
        let Material::Orthotropic { grain_angle_rad, .. } = self.material else {
            return Ok(vec![self.material.section(thickness, self.density_kg_m3, 0.0)?; mesh.tris.len()]);
        };
        let (sin, cos) = grain_angle_rad.sin_cos();
        let mut sections = Vec::with_capacity(mesh.tris.len());
        for facet in 0..mesh.tris.len() {
            let frame = mesh.facet(facet)?.frame;
            let normal = frame[2];
            let mut axial: [f64; 3] = core::array::from_fn(|i|
                if i == 2 { 1.0 - normal[2]*normal[i] } else { -normal[2]*normal[i] });
            let length = dot(axial, axial).sqrt();
            if !length.is_finite() || length <= 1e-12 {
                return Err("barrel grain needs a resolved axial direction in every facet".into());
            }
            for x in &mut axial { *x /= length; }
            // z x outward-normal follows increasing cylinder azimuth. The
            // angle is defined in physical space, independent of edge order.
            let hoop = cross(axial, normal);
            let grain = core::array::from_fn(|i| cos*axial[i] + sin*hoop[i]);
            let angle = dot(grain, frame[1]).atan2(dot(grain, frame[0]));
            sections.push(self.material.section(thickness, self.density_kg_m3, angle)?);
        }
        Ok(sections)
    }

    pub fn prepare(&self, drum: &drum_spec::Spec, dt_s: f64, audio: bool) -> Result<Prepared, Error> {
        self.validate()?; drum.validate()?;
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s * self.band_hz[1] >= 0.45
            || audio && (self.band_hz[0] < 40.0 || self.band_hz[1] > 1640.0) {
            return Err("barrel window exceeds the mechanical Nyquist guard or the 40..1640 Hz audio band".into());
        }
        let azimuths = drum.azimuths;
        let count = (self.axial_intervals + 1).checked_mul(azimuths).ok_or("barrel mesh size overflow")?;
        let facets = 2 * self.axial_intervals * azimuths;
        if count > 512 || facets > 1024 { return Err("barrel exceeds its 512-node/1024-facet mesh budget".into()); }
        let thickness = drum.outer_radius_m - drum.radius_m;
        let radius = drum.radius_m + 0.5 * thickness;
        let mut nodes = Vec::with_capacity(count);
        for level in 0..=self.axial_intervals {
            let z = drum.depth_m * (level as f64 / self.axial_intervals as f64 - 0.5);
            for i in 0..azimuths {
                let theta = core::f64::consts::TAU * i as f64 / azimuths as f64;
                nodes.push([radius * theta.cos(), radius * theta.sin(), z]);
            }
        }
        let mut triangles = Vec::with_capacity(facets);
        for level in 0..self.axial_intervals { for i in 0..azimuths {
            let j = (i + 1) % azimuths;
            let (a, b, c, d) = (level * azimuths + i, level * azimuths + j,
                (level + 1) * azimuths + i, (level + 1) * azimuths + j);
            triangles.extend([[a, b, d], [a, d, c]]);
        }}
        let mesh = ShellMesh::new(nodes, triangles)?;
        let sections = self.sections(&mesh, thickness)?;
        let clamps: Vec<_> = (0..azimuths).chain(count - azimuths..count).collect();
        let model = match self.material {
            Material::Isotropic { .. } => assemble_shell(&mesh, &sections[0], &clamps, ShellSupport::Clamped)?,
            Material::Orthotropic { .. } => assemble_shell_sections(&mesh, &sections, &clamps, ShellSupport::Clamped)?,
        };
        let omega = self.band_hz.map(|f| core::f64::consts::TAU * f);
        let modes = modes_shell(&model, (omega[0].powi(2), omega[1].powi(2)),
            &fs_plate::SliceOptions::default())?.modes;
        if modes.is_empty() || modes.len() > 63 {
            return Err("barrel window must contain 1..63 modes; no missing modes are invented or returned modes dropped".into());
        }
        let reduction = ShellReduction::new(&mesh, &sections, &model, &modes,
            ReductionBudget { max_modes: 63, max_facet_modes: 65_536, relative_tolerance: 1e-5 })?;
        let thicknesses = vec![thickness; count];
        let outer_positions = reduction.surface_positions(&thicknesses, ShellFace::Positive)?;
        let inner_positions = reduction.surface_positions(&thicknesses, ShellFace::Negative)?;
        let mut outer_weights = vec![vec![0.0; facets]; modes.len()];
        let mut inner_points = Vec::with_capacity(3 * facets);
        let mut inner_areas = Vec::with_capacity(3 * facets);
        let mut inner_weights = vec![Vec::with_capacity(3 * facets); modes.len()];
        for (f, triangle) in mesh.tris.iter().enumerate() {
            let (out_normal, _) = normal_area(triangle.map(|i| outer_positions[i]))?;
            let port = reduction.surface_point_port(&thicknesses, f, [1.0 / 3.0; 3], ShellFace::Positive, out_normal)?;
            for (row, weight) in outer_weights.iter_mut().zip(port.weights) { row[f] = weight; }
            // Original triangle winding points OUTWARD FROM THE GAS on its
            // negative skin. The solid's negative-face winding is the reverse.
            let (gas_normal, area) = normal_area(triangle.map(|i| inner_positions[i]))?;
            for bary in [[2.0/3.0,1.0/6.0,1.0/6.0], [1.0/6.0,2.0/3.0,1.0/6.0], [1.0/6.0,1.0/6.0,2.0/3.0]] {
                let port = reduction.surface_point_port(&thicknesses, f, bary, ShellFace::Negative, gas_normal)?;
                let [x, y, z] = port.position_m;
                let point = [x, y, 0.5 * drum.depth_m - z];
                if x.hypot(y) > drum.radius_m * (1.0 + 32.0 * f64::EPSILON)
                    || point[2] < 0.0 || point[2] > drum.depth_m {
                    return Err("barrel inner-skin quadrature lies outside the declared cylindrical cavity; refine the shell geometry".into());
                }
                inner_points.push(point); inner_areas.push(area / 3.0);
                for (row, weight) in inner_weights.iter_mut().zip(port.weights) { row.push(weight); }
            }
        }
        let compression = inner_weights.iter().map(|row| -row.iter().zip(&inner_areas).map(|(w, a)| w * a).sum::<f64>()).collect();
        Ok(Prepared { omegas: reduction.omegas().to_vec(), damping_ratio: self.damping_ratio,
            azimuths, axial_intervals: self.axial_intervals, outer_positions, triangles: mesh.tris,
            outer_weights, inner_points, inner_areas, inner_weights, compression,
            reduction, nodal_thicknesses: thicknesses,
            radius_m: drum.radius_m, outer_radius_m: drum.outer_radius_m, depth_m: drum.depth_m })
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}

/// Cached physical projections; no independent dynamic state or force clock.
pub struct Prepared {
    omegas: Vec<f64>, damping_ratio: f64,
    // Keep the admitted operator for physical point contacts on the same skin.
    // Its translations AND rotations remain in the pressure/radiation basis.
    pub(super) reduction: ShellReduction,
    pub(super) nodal_thicknesses: Vec<f64>,
    pub(super) azimuths: usize,
    pub(super) axial_intervals: usize,
    pub(super) outer_positions: Vec<[f64; 3]>,
    pub(super) triangles: Vec<[usize; 3]>,
    pub(super) outer_weights: Vec<Vec<f64>>,
    inner_points: Vec<[f64; 3]>, inner_areas: Vec<f64>, inner_weights: Vec<Vec<f64>>,
    compression: Vec<f64>,
    pub(super) radius_m: f64,
    pub(super) outer_radius_m: f64,
    pub(super) depth_m: f64,
}
impl Prepared {
    pub fn mode_count(&self) -> usize { self.omegas.len() }
    pub fn body(&self) -> ImpactBody {
        let mut body = super::zero_body(BodyPotential::Linear(self.omegas.clone()), &self.omegas);
        body.damping_per_s = self.omegas.iter().map(|w| 2.0 * self.damping_ratio * w).collect();
        body
    }
    pub fn compression_areas(&self) -> &[f64] { &self.compression }
    pub(super) fn coupling(&self, air: &CylindricalCavity, gate: &CancelGate) -> Result<Vec<f64>, Error> {
        if gate.is_requested() { return Err("barrel cavity assembly cancelled".into()); }
        if air.spec().radius_m != self.radius_m || air.spec().depth_m != self.depth_m {
            return Err("barrel and cavity dimensions disagree".into());
        }
        let sampled = air.sample(&self.inner_points, 24_576)?;
        let structure = StructuralModes { omegas: self.omegas.clone(), shapes: self.inner_weights.clone(), loss_factor: 0.0 };
        let coupling = assemble_coupling(&structure, &sampled, &self.inner_areas)?;
        if gate.is_requested() { return Err("barrel cavity assembly cancelled".into()); }
        Ok(coupling)
    }
}

pub(super) fn normal_area(p: [[f64; 3]; 3]) -> Result<([f64; 3], f64), Error> {
    let a: [f64; 3] = core::array::from_fn(|i| p[1][i] - p[0][i]);
    let b: [f64; 3] = core::array::from_fn(|i| p[2][i] - p[0][i]);
    let n = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
    let length = n.iter().map(|v| v*v).sum::<f64>().sqrt();
    if !length.is_finite() || length <= 0.0 { return Err("barrel skin has a collapsed or nonfinite facet".into()); }
    Ok((n.map(|v| v / length), 0.5 * length))
}

#[cfg(test)]
#[path = "barrel_tests.rs"]
mod tests;
