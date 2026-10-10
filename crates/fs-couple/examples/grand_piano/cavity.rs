//! Explicit sealed rectangular gas beneath a planar soundboard interface.
//! Positive board z is away from the gas. The supplied upper-face origin and
//! dimensions declare rigid walls outside the actual structural patch; no
//! closed piano enclosure, geometry, gas state or damping is inferred.
//! The existing pressure-mode and physical-motion owners provide the basis.
use std::io::Read;
use fs_couple::{render::plate::impact::cavity::CavityCoupling,
    vibroacoustic::{AcousticMedium, CavityModes, rectangular_cavity_modes}};
use fs_material::gas::{GasSpec, GasState};
use super::{board_geometry::{PreparedBoard, motion::{MotionSurface, SourceBridgePort,
    EDGE_CUBIC_QUADRATURE}}, linear::{Bank, BoardMode, MAX_BOARD_MODES}};

pub const HEADER: &str = "frankensim-piano-cavity-si-v1";
const MAX_BYTES: u64 = 64 * 1024;
const MAX_CAVITY_MODES: usize = 8;
const MAX_INTERFACE_SAMPLES: usize = 262_144;
// Count overlap accumulation products. Motion evaluation adds bounded work
// per structural coordinate; this is not a wall-time or full-operation claim.
const MAX_INTEGRATION_TERMS: usize = 250_000_000;
const MAX_REFINEMENT: usize = 6;
const MAX_RELATIVE_CHANGE: f64 = 1e-6;

#[derive(Clone, Debug)]
pub struct Specification {
    source: String,
    origin_m: [f64; 3],
    dimensions_m: [f64; 3],
    count: usize,
    damping_ratio: f64,
    medium: AcousticMedium,
    temperature_k: f64,
    pressure_pa: f64,
}

#[derive(Debug)]
pub struct Projected {
    // The same mode normalization sampled at the interface origin; complete
    // quadrature tables are cold scratch, not retained acoustic state.
    cavity: CavityModes,
    structural: usize,
    overlaps: Vec<f64>, // row-major structural-by-acoustic, bare board basis
    damping_per_s: Vec<f64>,
    description: String,
}

fn finite(text: &str) -> Result<f64, String> {
    text.parse::<f64>().ok().filter(|v| v.is_finite())
        .ok_or_else(|| "cavity row requires finite SI values".into())
}
fn once<T>(slot: &mut Option<T>, value: T, row: &str) -> Result<(), String> {
    if slot.is_some() { return Err(format!("duplicate cavity {row} row")); }
    *slot = Some(value); Ok(())
}

impl Specification {
    pub fn load(path: &str) -> Result<Self, String> {
        let mut text = String::new();
        std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?
            .take(MAX_BYTES + 1).read_to_string(&mut text).map_err(|e| format!("{path}: {e}"))?;
        Self::read(&text)
    }

    /// All rows are required. `interface-origin-m` is the upper-face corner;
    /// gas occupies x..x+Lx, y..y+Ly, z-Lz..z. The damping ratio specifies
    /// causal acoustic momentum drag, not a hysteretic material loss factor.
    pub fn read(text: &str) -> Result<Self, String> {
        if text.len() as u64 > MAX_BYTES { return Err("cavity specification exceeds 64 KiB".into()); }
        let mut rows = text.lines().map(|line| line.split('#').next().unwrap_or("").trim())
            .filter(|row| !row.is_empty());
        if rows.next() != Some(HEADER) { return Err(format!("expected {HEADER}")); }
        let mut source = None; let mut origin = None; let mut dimensions = None;
        let mut count = None; let mut damping = None; let mut gas = None;
        for row in rows {
            let fields: Vec<_> = row.split(',').map(str::trim).collect();
            match fields.as_slice() {
                ["source", authority, attribution @ ..] if !attribution.is_empty()
                    && ["estimated", "mixed", "published", "measured"].contains(authority)
                    && attribution.iter().any(|s| !s.is_empty()) => {
                    once(&mut source, format!("{authority}: {}", attribution.join(",")), "source")?;
                }
                ["interface-origin-m", x, y, z] => {
                    once(&mut origin, [finite(x)?, finite(y)?, finite(z)?], "interface-origin-m")?;
                }
                ["dimensions-m", x, y, z] => {
                    once(&mut dimensions, [finite(x)?, finite(y)?, finite(z)?], "dimensions-m")?;
                }
                ["modes", n] => {
                    once(&mut count, n.parse::<usize>().map_err(|_| "invalid cavity mode count")?, "modes")?;
                }
                ["damping-ratio", zeta] => { once(&mut damping, finite(zeta)?, "damping-ratio")?; }
                ["gas", "dry-air-ussa1976", temperature, pressure] => {
                    let temperature = finite(temperature)?; let pressure = finite(pressure)?;
                    let state = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure)
                        .map_err(|e| format!("cavity gas: {e}"))?;
                    once(&mut gas, (AcousticMedium { rho0: state.density, c0: state.sound_speed },
                        temperature, pressure), "gas")?;
                }
                _ => return Err(format!("unknown cavity row or wrong field count: {row}")),
            }
        }
        let origin_m = origin.ok_or("missing cavity interface-origin-m row")?;
        let dimensions_m = dimensions.ok_or("missing cavity dimensions-m row")?;
        let count = count.ok_or("missing cavity modes row")?;
        let damping_ratio = damping.ok_or("missing cavity damping-ratio row")?;
        let (medium, temperature_k, pressure_pa) = gas.ok_or("missing cavity gas row")?;
        if origin_m.iter().any(|v| v.abs() > 100.)
            || dimensions_m.iter().any(|v| !(1e-4..=20.).contains(v))
            || !(1..=MAX_CAVITY_MODES).contains(&count)
            || !(0.0..1.0).contains(&damping_ratio) {
            return Err("cavity requires origin within 100 m, dimensions in [0.0001,20] m, 1..=8 pressure modes and damping ratio in [0,1)".into());
        }
        Ok(Self { source: source.ok_or("missing cavity source attribution")?, origin_m,
            dimensions_m, count, damping_ratio, medium, temperature_k, pressure_pa })
    }

    pub fn project(&self, board: &PreparedBoard) -> Result<Projected, String> {
        self.project_motion(&board.modes, board.motion.as_ref()
            .ok_or("cavity coupling requires the complete geometric board motion")?)
    }

    pub fn project_motion(&self, modes: &[BoardMode], motion: &MotionSurface)
        -> Result<Projected, String> {
        if !(1..=MAX_BOARD_MODES).contains(&modes.len()) || motion.shapes.len() != modes.len()
            || motion.mesh.nodes.is_empty() || motion.mesh.tris.is_empty()
            || motion.mesh.nodes.len() > 20_000 || motion.mesh.tris.len() > 40_000
            || motion.shapes.iter().any(|shape| shape.len() != motion.mesh.nodes.len()) {
            return Err("cavity interface needs the complete finite geometric board basis".into());
        }
        // Check the complete interface before sampling any overlap. A crown or
        // a patch beyond the declared wall must not silently become a box lid.
        let tolerance = 1e-9 * self.dimensions_m.iter().copied().fold(1.0_f64, f64::max);
        for point in &motion.mesh.nodes {
            if point.iter().any(|x| !x.is_finite())
                || (point[2] - self.origin_m[2]).abs() > tolerance
                || (0..2).any(|i| point[i] < self.origin_m[i] - tolerance
                    || point[i] > self.origin_m[i] + self.dimensions_m[i] + tolerance) {
                return Err("cavity interface must lie on the declared planar upper face and inside its rectangle; no crown flattening or wall extrapolation".into());
            }
        }
        for shape in &motion.shapes {
            if shape.iter().flatten().any(|v| !v.is_finite()) {
                return Err("cavity interface has nonfinite structural motion".into());
            }
        }
        let cavity = self.basis(&[[0., 0.]])?;
        let damping_per_s: Vec<_> = cavity.omegas.iter().map(|w| 2. * self.damping_ratio * w).collect();
        let mut cells = vec![[[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]];
        let mut previous: Option<(Vec<f64>, Vec<f64>)> = None;
        let mut work = 0usize;
        for level in 0..=MAX_REFINEMENT {
            let samples = motion.mesh.tris.len().checked_mul(4 * cells.len())
                .ok_or("cavity interface sample count overflow")?;
            let terms = samples.checked_mul(modes.len() * self.count)
                .ok_or("cavity overlap work count overflow")?;
            work = work.checked_add(terms).ok_or("cavity overlap work count overflow")?;
            if samples > MAX_INTERFACE_SAMPLES || work > MAX_INTEGRATION_TERMS {
                return Err("cavity interface quadrature exceeds its explicit sample/work budget; refine the structural interface or reduce the retained cavity/board basis".into());
            }
            let (overlaps, scales) = self.integrate(motion, &cells, samples)?;
            if let Some((old, old_scales)) = &previous {
                let relative = overlaps.iter().zip(old).zip(scales.iter().zip(old_scales))
                    .map(|((a,b),(scale,old_scale))| (a-b).abs()
                        / scale.max(*old_scale).max(f64::MIN_POSITIVE)).fold(0.0_f64, f64::max);
                if relative <= MAX_RELATIVE_CHANGE {
                    let description = format!("explicit sealed rectangular cavity beneath planar board: source {}; upper-face origin {:?} m, dimensions {:?} m, {} pressure modes including uniform compression; dry-air-ussa1976 at {} K and {} Pa, rho={} kg/m3, c={} m/s; viscous acoustic damping ratio {}; {} interface samples at refinement {}, relative overlap-change estimate {:.3e} (not an integration or modal-truncation certificate); rigid remaining walls, no inferred grand-piano enclosure",
                        self.source, self.origin_m, self.dimensions_m, self.count, self.temperature_k,
                        self.pressure_pa, self.medium.rho0, self.medium.c0, self.damping_ratio,
                        samples, level, relative);
                    return Ok(Projected { cavity, structural: modes.len(), overlaps, damping_per_s, description });
                }
            }
            previous = Some((overlaps, scales));
            if level < MAX_REFINEMENT {
                let mut children = Vec::with_capacity(4 * cells.len());
                for [a,b,c] in cells {
                    let midpoint = |x: [f64;3], y: [f64;3]| std::array::from_fn(|i| (x[i]+y[i])*0.5);
                    let ab = midpoint(a,b); let bc = midpoint(b,c); let ca = midpoint(c,a);
                    children.extend([[a,ab,ca], [ab,b,bc], [ca,bc,c], [ab,bc,ca]]);
                }
                cells = children;
            }
        }
        Err("cavity interface quadrature did not resolve the requested cosine/structural overlap; no unresolved coupling admitted".into())
    }

    fn basis(&self, points: &[[f64;2]]) -> Result<CavityModes, String> {
        rectangular_cavity_modes(self.dimensions_m[0], self.dimensions_m[1], self.dimensions_m[2],
            self.medium, 0.0, self.count, points).map_err(|e| e.to_string())
    }

    fn integrate(&self, motion: &MotionSurface, cells: &[[[f64;3];3]], samples: usize)
        -> Result<(Vec<f64>, Vec<f64>), String> {
        let mut points = Vec::with_capacity(samples); let mut sites = Vec::with_capacity(samples);
        for (element, tri) in motion.mesh.tris.iter().enumerate() {
            let area = motion.mesh.facet(element).map_err(|e| e.to_string())?.area_m2;
            for cell in cells { for (barycentric, weight) in EDGE_CUBIC_QUADRATURE {
                let weights: [f64;3] = std::array::from_fn(|i|
                    (0..3).map(|j| barycentric[j] * cell[j][i]).sum());
                let point = std::array::from_fn(|coordinate|
                    (0..3).map(|j| weights[j] * (motion.mesh.nodes[tri[j]][coordinate]
                        - self.origin_m[coordinate])).sum());
                points.push(point); sites.push((element, weights, area * weight / cells.len() as f64));
            } }
        }
        let sampled = self.basis(&points)?;
        let mut overlaps: Vec<f64> = vec![0.; motion.shapes.len() * self.count];
        let mut scales = vec![0.; overlaps.len()];
        for (i, &(triangle, weights, area)) in sites.iter().enumerate() {
            let projection = SourceBridgePort { triangle, weights, arm_m: [0.;3], direction: [0.,0.,1.] }
                .prepare(&motion.mesh, motion.is_edge_cubic())?;
            for (r, shape) in motion.shapes.iter().enumerate() {
                let (phi, _) = projection.project_nodal(projection.nodes().map(|node| shape[node]))?;
                for j in 0..self.count {
                    let term = area * phi * sampled.interface[j][i];
                    overlaps[r*self.count+j] += term; scales[r*self.count+j] += term.abs();
                }
            }
        }
        if overlaps.iter().chain(&scales).any(|v| !v.is_finite()) {
            return Err("cavity interface overlap overflow".into());
        }
        Ok((overlaps, scales))
    }
}

impl Projected {
    /// Bind this same physical interface to the actual string-loaded board.
    /// Cavity normalization and full wood damping stay with their owners.
    pub fn loaded(&self, bank: &Bank) -> Result<CavityCoupling, String> {
        if bank.board_count != self.structural {
            return Err("cavity projection and loaded board bases differ".into());
        }
        let count = self.cavity.omegas.len();
        let mut loaded = vec![0.; self.overlaps.len()];
        for j in 0..count {
            let bare: Vec<_> = (0..self.structural).map(|r| self.overlaps[r*count+j]).collect();
            let row = bank.project_board_shape(&bare)?;
            for (r, value) in row.into_iter().enumerate() { loaded[r*count+j] = value; }
        }
        CavityCoupling::new_with_mode_budget(&self.cavity, self.structural, &loaded,
            &self.damping_per_s, MAX_BOARD_MODES + MAX_CAVITY_MODES).map_err(|e| e.to_string())
    }
    pub fn report(&self) -> String { self.description.clone() }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn card(dimensions: &str, count: usize) -> String {
        format!("{HEADER}\nsource,estimated,authored sealed-box coupon\ninterface-origin-m,0,0,0\ndimensions-m,{dimensions}\nmodes,{count}\ndamping-ratio,0.02\ngas,dry-air-ussa1976,293.15,101325\n")
    }
    fn modes() -> Vec<BoardMode> {
        vec![BoardMode { frequency_hz: 150., damping_ratio: 0.01, bridge: [1.;88], volume: 0.4 },
            BoardMode { frequency_hz: 250., damping_ratio: 0.01, bridge: [-0.4;88], volume: 0.2 }]
    }
    fn piston() -> MotionSurface {
        let mesh = fs_plate::ShellMesh::new(vec![[0.,0.,0.], [1.,0.,0.], [1.,0.4,0.], [0.,0.4,0.]],
            vec![[0,1,2], [0,2,3]]).unwrap();
        let uniform = vec![[0.,0.,1.,0.,0.,0.];4];
        let linear = mesh.nodes.iter().map(|p| [0.,0.,p[0],0.,-1.,0.]).collect();
        MotionSurface::new(mesh, vec![uniform, linear]).unwrap()
    }

    #[test]
    fn geometric_overlap_and_cavity_norm_match_independent_rectangle_integrals() {
        let spec = Specification::read(&card("1,0.4,0.2", 3)).unwrap();
        let p = spec.project_motion(&modes(), &piston()).unwrap();
        let pi = std::f64::consts::PI;
        assert_eq!(p.cavity.omegas[0], 0.);
        assert!((p.cavity.omegas[1] - spec.medium.c0*pi).abs() < 1e-10);
        assert!((p.cavity.omegas[2] - spec.medium.c0*2.*pi).abs() < 1e-10);
        assert_eq!(p.damping_per_s[0], 0.);
        for j in 1..3 { assert_eq!(p.damping_per_s[j], 0.04*p.cavity.omegas[j]); }
        for (actual, expected) in p.cavity.lambdas.iter().zip([0.08,0.04,0.04]) {
            assert!((actual-expected).abs() < 1e-15);
        }
        // Integrals over [0,1] x [0,.4]: constant; x; and x*cos(pi*x).
        let expected = [0.4,0.,0.,0.2,-0.8/(pi*pi),0.];
        for (actual, expected) in p.overlaps.iter().zip(expected) {
            assert!((actual-expected).abs() < 2e-7, "{actual} versus {expected}");
        }
        assert!(p.report().contains("sealed rectangular cavity"));
        assert!(p.report().contains("not an integration or modal-truncation certificate"));
    }

    #[test]
    fn cubic_interface_retains_signed_volume_and_actual_loaded_pressure_normalization() {
        let mesh = fs_plate::ShellMesh::new(vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]],
            vec![[0,1,2]]).unwrap();
        // Existing reconstruction: w=x^3-(1-x-y)*x*y, whose integral is
        // 1/20 - 1/120 = 1/24. P1 interpolation of these nodes gives 1/6.
        let q = vec![[0.;6], [0.,0.,1.,0.,-3.,0.], [0.;6]];
        let negative = q.iter().map(|node| node.map(|v| -0.5*v)).collect();
        let motion = MotionSurface::new_edge_cubic(mesh, vec![q,negative]).unwrap();
        let spec = Specification::read(&card("1,1,0.2", 1)).unwrap();
        let p = spec.project_motion(&modes(), &motion).unwrap();
        assert!((p.overlaps[0]-1./24.).abs() < 1e-14);
        assert!((p.overlaps[1]+1./48.).abs() < 1e-14);
        let course = super::super::geometry::demonstration_scale().unwrap()[48];
        let bank = Bank::new(&[course], &modes(), 192_000, 1000., 3, false).unwrap();
        let loaded = p.loaded(&bank).unwrap();
        assert_eq!(loaded.cavity_modes(), 1);
        assert_eq!(loaded.total_modes(), 2, "uniform compression has no invented air momentum");
        let state = [0.003,0.,-0.002,0.];
        let (basis, _) = bank.board_pressure_basis().unwrap();
        let bare_q: [f64;2] = std::array::from_fn(|i| basis[2*i]*state[0] + basis[2*i+1]*state[2]);
        let volume: f64 = p.overlaps.iter().zip(bare_q).map(|(area,q)| area*q).sum();
        let expected = -spec.medium.rho0*spec.medium.c0.powi(2)/0.2*volume;
        let mut pressure = [0.]; loaded.pressures_into(&state, &mut pressure).unwrap();
        assert!((pressure[0]-expected).abs() < 1e-11*expected.abs());
        assert!(p.loaded(&Bank::new(&[course], &modes()[..1], 192_000, 1000., 3, false).unwrap()).is_err());
    }

    #[test]
    fn incomplete_cards_nonplanar_interfaces_and_outside_patches_refuse() {
        let good = card("1,0.4,0.2", 3);
        for invalid in [good.replace("modes,3", "modes,0"), good.replace("modes,3", "modes,9"),
            good.replace("dimensions-m,1,0.4,0.2", "dimensions-m,1,0.4,0"),
            good.replace("damping-ratio,0.02", "damping-ratio,1"),
            good.replace("dry-air-ussa1976", "unknown-gas"), good.replace("293.15", "NaN"),
            good.replace("source,estimated,authored sealed-box coupon\n", ""),
            format!("{good}modes,3\n")] {
            assert!(Specification::read(&invalid).is_err());
        }
        let spec = Specification::read(&good).unwrap();
        let mut crown = piston(); crown.mesh.nodes[2][2] = 0.001;
        assert!(spec.project_motion(&modes(), &crown).unwrap_err().contains("planar upper face"));
        let mut outside = piston(); outside.mesh.nodes[1][0] = 1.01;
        assert!(spec.project_motion(&modes(), &outside).unwrap_err().contains("inside its rectangle"));
        assert!(spec.project_motion(&modes()[..1], &piston()).is_err());
    }
}
