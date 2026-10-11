//! A normal point impact on the actual outer skin of an elastic barrel.
//! The existing tip contact and time owner supply force, storage and reaction.
//! Azimuth locates a skin point; the local facet normal defines strike motion.
use fs_plate::shell::reduction::radiation::{ShellFace, ShellSurfacePort};
use super::{Error, barrel};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spec {
    pub azimuth_rad: f64,
    pub height_below_batter_m: f64,
}

pub fn option(args: &mut Vec<String>) -> Result<Option<Spec>, Error> {
    let mut positions = args.iter().enumerate().filter(|(_, a)| a.as_str() == "--barrel-strike");
    let Some((at, _)) = positions.next() else { return Ok(None); };
    if positions.next().is_some() { return Err("--barrel-strike may be supplied only once".into()); }
    let azimuth = args.get(at + 1).ok_or("--barrel-strike needs AZIMUTH_RAD HEIGHT_BELOW_BATTER_M")?;
    let height = args.get(at + 2).ok_or("--barrel-strike needs AZIMUTH_RAD HEIGHT_BELOW_BATTER_M")?;
    let spec = Spec { azimuth_rad: azimuth.parse()?, height_below_batter_m: height.parse()? };
    spec.validate()?;
    args.drain(at..at + 3);
    Ok(Some(spec))
}

pub fn admit_command(spec: Option<&Spec>, command: &str) -> Result<(), Error> {
    if spec.is_some() && !matches!(command, "drum" | "drum-wav" | "drum-mic" |
        "drum-modal" | "drum-modal-wav" | "drum-modal-mic" |
        "drum-stretch" | "drum-stretch-wav" | "drum-stretch-mic" |
        "snare" | "snare-wav" | "snare-mic" | "snare-off" | "snare-off-wav" | "snare-off-mic") {
        return Err("--barrel-strike requires a drum or snare command with --elastic-barrel".into());
    }
    Ok(())
}

impl Spec {
    fn validate(self) -> Result<(), Error> {
        if !self.azimuth_rad.is_finite() || self.azimuth_rad.abs() > core::f64::consts::TAU
            || !self.height_below_batter_m.is_finite() || self.height_below_batter_m <= 0.0 {
            return Err("barrel strike needs azimuth in [-2*pi,2*pi] rad and positive finite height below the batter [m]".into());
        }
        Ok(())
    }

    /// Positive OUTWARD skin displacement, including the moment of a force at
    /// its thickness offset. For an inward-positive stick row s, contact uses
    /// s + this row: an outward-moving wall closes the exterior tip gap.
    /// No amplitude normalization, state change, or nearest-node strike occurs.
    pub fn port(&self, barrel: &barrel::Prepared) -> Result<ShellSurfacePort, Error> {
        let hit = self.locate(barrel)?;
        let port = barrel.reduction.surface_point_port(&barrel.nodal_thicknesses,
            hit.facet, hit.barycentric, ShellFace::Positive, hit.normal)?;
        if port.weights.iter().all(|w| *w == 0.0) {
            return Err("barrel strike has no retained surface motion; choose a resolved interior station or enrich the declared window".into());
        }
        Ok(port)
    }

    fn locate(&self, barrel: &barrel::Prepared) -> Result<Hit, Error> {
        self.validate()?;
        let half_depth = 0.5 * barrel.depth_m;
        let z = half_depth - self.height_below_batter_m;
        if self.height_below_batter_m >= barrel.depth_m || z <= -half_depth || z >= half_depth {
            return Err("barrel strike must lie strictly between the rigid hoops; no clamp or axial snapping".into());
        }
        let origin = [0.0, 0.0, z];
        let direction = [self.azimuth_rad.cos(), self.azimuth_rad.sin(), 0.0];
        // This is floating-point edge recognition, not a physical contact
        // radius. Inputs outside the axial domain are refused above, not moved.
        let edge_roundoff = 128.0 * f64::EPSILON;
        let distance_roundoff = edge_roundoff * barrel.outer_radius_m.max(barrel.depth_m);
        let mut nearest: Option<Hit> = None;
        for (facet, triangle) in barrel.triangles.iter().enumerate() {
            let points = triangle.map(|i| barrel.outer_positions[i]);
            let (normal, area) = barrel::normal_area(points)?;
            let denominator = dot(normal, direction);
            if denominator <= edge_roundoff { continue; }
            let distance = dot(normal, sub(points[0], origin)) / denominator;
            if !distance.is_finite() || distance <= 0.0 { continue; }
            let point = core::array::from_fn(|i| origin[i] + distance * direction[i]);
            let a = sub(points[1], points[0]); let b = sub(points[2], points[0]);
            let p = sub(point, points[0]);
            let u = dot(cross(p, b), normal) / (2.0 * area);
            let v = dot(cross(a, p), normal) / (2.0 * area);
            let mut barycentric = [1.0 - u - v, u, v];
            if barycentric.iter().any(|b| !b.is_finite() || *b < -edge_roundoff || *b > 1.0 + edge_roundoff) {
                continue;
            }
            for b in &mut barycentric { *b = b.clamp(0.0, 1.0); }
            let sum: f64 = barycentric.iter().sum();
            for b in &mut barycentric { *b /= sum; }
            let hit = Hit { facet, barycentric, normal, distance };
            if let Some(old) = &nearest {
                if (distance - old.distance).abs() <= distance_roundoff {
                    // A ray on a shared edge/vertex has one physical point.
                    // Retain the lowest facet index and that facet's normal;
                    // never average normals or duplicate the contact force.
                    if !triangle.iter().any(|i| barrel.triangles[old.facet].contains(i))
                        || !barycentric.iter().any(|b| *b <= edge_roundoff)
                        || !old.barycentric.iter().any(|b| *b <= edge_roundoff) {
                        return Err("barrel strike has ambiguous coincident outer-skin intersections".into());
                    }
                    continue;
                }
                if distance >= old.distance { continue; }
            }
            nearest = Some(hit);
        }
        nearest.ok_or_else(|| "barrel strike ray misses the admitted outer skin".into())
    }
}

struct Hit { facet: usize, barycentric: [f64; 3], normal: [f64; 3], distance: f64 }
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { core::array::from_fn(|i| a[i] - b[i]) }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}

#[cfg(test)]
#[path = "barrel_strike_tests.rs"]
mod tests;
