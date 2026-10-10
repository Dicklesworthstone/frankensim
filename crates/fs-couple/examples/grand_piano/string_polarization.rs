//! Two transverse string directions from the played board's actual motion.
//! The supplied bridge site, arm and orthonormal string/hammer frame are
//! physical input. Neither a lateral bridge coefficient nor a damper ratio is
//! inferred from the primary displacement, a material name or a piano preset.
use std::{collections::BTreeMap, io::Read};
use super::{board_geometry::motion::{MotionSurface, SourceBridgeFrame, SourceBridgePort},
    geometry::Course, linear::BoardMode};

pub const HEADER: &str = "frankensim-piano-string-polarization-v1";
const MAX_BYTES: u64 = 64 * 1024;

#[derive(Clone, Debug)]
struct Site {
    triangle: usize,
    weights: [f64; 3],
    arm_m: [f64; 3],
    string_axis: [f64; 3],
    hammer_axis: [f64; 3],
    lateral_damper_ratio: f64,
}

#[derive(Clone, Debug)]
pub struct Specification {
    sites: BTreeMap<u8, Site>,
    pub source: String,
}

pub struct Prepared {
    lateral: Vec<Vec<f64>>,
    damper_ratios: Vec<f64>,
    pub source: String,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a.iter().zip(b).map(|(a, b)| a * b).sum() }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

impl Specification {
    pub fn load(path: &str, courses: &[Course]) -> Result<Self, String> {
        let mut text = String::new();
        std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?
            .take(MAX_BYTES + 1).read_to_string(&mut text).map_err(|e| format!("{path}: {e}"))?;
        Self::read(&text, courses)
    }

    pub fn read(text: &str, courses: &[Course]) -> Result<Self, String> {
        if text.len() as u64 > MAX_BYTES { return Err("string polarization file exceeds 64 KiB".into()); }
        let mut header = false;
        let mut source = None;
        let mut sites = BTreeMap::new();
        for (line, raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() { continue; }
            let bad = || format!("string polarization line {}: invalid source or complete bridge-frame row", line + 1);
            if !header {
                if row != HEADER { return Err(bad()); }
                header = true;
                continue;
            }
            let fields: Vec<_> = row.split(',').map(str::trim).collect();
            if fields.first() == Some(&"source") {
                if fields.len() != 3 || source.is_some() || fields[2].is_empty()
                    || !["estimated", "mixed", "published", "measured"].contains(&fields[1]) {
                    return Err(bad());
                }
                source = Some(format!("{}: {}", fields[1], fields[2]));
                continue;
            }
            if fields.len() != 16 || fields[0] != "course" { return Err(bad()); }
            let key = fields[1].parse::<u8>().map_err(|_| bad())?;
            let triangle = fields[2].parse::<usize>().map_err(|_| bad())?;
            let values: Vec<f64> = fields[3..].iter().map(|v| v.parse().map_err(|_| bad()))
                .collect::<Result<_, _>>()?;
            if values.iter().any(|v| !v.is_finite()) { return Err(bad()); }
            let site = Site {
                triangle, weights: [values[0], values[1], values[2]],
                arm_m: [values[3], values[4], values[5]],
                string_axis: [values[6], values[7], values[8]],
                hammer_axis: [values[9], values[10], values[11]],
                lateral_damper_ratio: values[12],
            };
            if site.weights.iter().any(|w| !(0.0..=1.0).contains(w))
                || (site.weights.iter().sum::<f64>() - 1.0).abs() > 1e-10
                || (dot(site.string_axis, site.string_axis) - 1.0).abs() > 1e-10
                || (dot(site.hammer_axis, site.hammer_axis) - 1.0).abs() > 1e-10
                || dot(site.string_axis, site.hammer_axis).abs() > 1e-10
                || !(0.0..=10.0).contains(&site.lateral_damper_ratio) {
                return Err(format!("string polarization key {key}: require in-triangle weights, orthonormal string/hammer axes and lateral damper ratio in [0,10]"));
            }
            if sites.insert(key, site).is_some() { return Err(format!("duplicate string polarization key {key}")); }
        }
        let source = source.ok_or("string polarization requires a source attribution row")?;
        let spec = Self { sites, source };
        spec.validate_scale(courses)?;
        Ok(spec)
    }

    fn validate_scale(&self, courses: &[Course]) -> Result<(), String> {
        if courses.is_empty() || courses.len() > 88 || self.sites.len() != courses.len()
            || courses.iter().enumerate().any(|(i, c)| courses[..i].iter().any(|p| p.midi == c.midi)) {
            return Err("string polarization must cover the complete unique scale, including silent keys".into());
        }
        for c in courses {
            c.validate()?;
            if !self.sites.contains_key(&c.midi) { return Err(format!("missing string polarization key {}", c.midi)); }
        }
        Ok(())
    }

    /// Expose the complete admitted source frames before modal reduction.
    /// These are physical probes, not coefficients copied from a retained
    /// basis: a lateral-only source mode must participate in Ritz selection.
    /// Rows follow the supplied scale order, including its silent courses.
    pub fn source_ports(&self, courses: &[Course]) -> Result<Vec<SourceBridgeFrame>, String> {
        self.validate_scale(courses)?;
        Ok(courses.iter().map(|course| {
            let site = &self.sites[&course.midi];
            let primary = SourceBridgePort { triangle: site.triangle, weights: site.weights,
                arm_m: site.arm_m, direction: site.hammer_axis };
            let secondary = SourceBridgePort { direction: cross(site.string_axis, site.hammer_axis),
                ..primary };
            SourceBridgeFrame { midi: course.midi, primary, secondary }
        }).collect())
    }

    /// Project u + theta x arm in the SAME bare-board basis as primary bridge
    /// mechanics. The primary check refuses incompatible sites or frames that
    /// would silently change the struck plane; it never repairs supplied input.
    pub fn project(&self, courses: &[Course], board: &[BoardMode], motion: Option<&MotionSurface>)
        -> Result<Prepared, String> {
        let frames = self.source_ports(courses)?;
        let motion = motion.ok_or("string polarization requires full-vector geometric board motion; modal CSV has no lateral-motion source")?;
        if board.is_empty() || motion.shapes.len() != board.len()
            || motion.shapes.iter().any(|m| m.len() != motion.mesh.nodes.len()) {
            return Err("string polarization requires the same complete board motion basis as primary mechanics".into());
        }
        let mut lateral = Vec::with_capacity(courses.len());
        let mut damper_ratios = Vec::with_capacity(courses.len());
        for (course, frame) in courses.iter().zip(frames) {
            let site = &self.sites[&course.midi];
            let (primary, scales) = motion.project_at(frame.primary.triangle, frame.primary.weights,
                frame.primary.arm_m, frame.primary.direction)
                .map_err(|e| format!("string polarization key {}: {e}", course.midi))?;
            let (secondary, _) = motion.project_at(frame.secondary.triangle, frame.secondary.weights,
                frame.secondary.arm_m, frame.secondary.direction)
                .map_err(|e| format!("string polarization key {}: {e}", course.midi))?;
            let mut row = Vec::with_capacity(board.len());
            for (mode_index, mode) in board.iter().enumerate() {
                let primary = primary[mode_index];
                let secondary = secondary[mode_index];
                let scale = scales[mode_index];
                let expected = mode.bridge[usize::from(course.midi - 21)];
                if !primary.is_finite() || !secondary.is_finite() || !expected.is_finite() || !scale.is_finite() {
                    return Err("string polarization motion projection is not finite".into());
                }
                if (primary - expected).abs() > 1e-10 * scale.max(expected.abs()).max(f64::MIN_POSITIVE) {
                    return Err(format!("string polarization key {}, mode {mode_index}: supplied hammer-plane projection does not match the existing bridge site/basis", course.midi));
                }
                row.push(secondary);
            }
            lateral.push(row);
            damper_ratios.push(site.lateral_damper_ratio);
        }
        Ok(Prepared { lateral, damper_ratios, source: self.source.clone() })
    }
}

impl Prepared {
    pub fn secondary(&self) -> (&[Vec<f64>], &[f64]) { (&self.lateral, &self.damper_ratios) }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn course() -> Course { super::super::geometry::demonstration_scale().unwrap()[48] }
    fn card(row: &str) -> String { format!("{HEADER}\nsource,estimated,analytic supplied geometry\n{row}\n") }
    const ROW: &str = "course,69,0,0.2,0.3,0.5,0,0,0.02,0,1,0,0,0,1,0.4";
    fn motion() -> MotionSurface {
        let mesh = fs_plate::ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]], vec![[0,1,2]]).unwrap();
        MotionSurface::new(mesh, vec![vec![[1.,2.,3.,4.,5.,6.]; 3]]).unwrap()
    }
    fn board() -> Vec<BoardMode> {
        vec![BoardMode { frequency_hz: 150., damping_ratio: 0., bridge: [3.; 88], volume: 1. }]
    }
    #[test]
    fn full_vector_motion_and_physical_arm_determine_both_transverse_ports() {
        let spec = Specification::read(&card(ROW), &[course()]).unwrap();
        let projected = spec.project(&[course()], &board(), Some(&motion())).unwrap();
        // y-string, z-hammer => lateral +x. theta_y * arm_z adds 0.1 m
        // to x; primary z stays 3. No fit of a lateral/vertical coefficient.
        assert!((projected.secondary().0[0][0] - 1.1).abs() < 1e-14);
        assert_eq!(projected.secondary().1, &[0.4]);
        let reversed = Specification::read(&card(&ROW.replace("0,1,0,0,0,1", "0,-1,0,0,0,1")), &[course()]).unwrap();
        assert!((reversed.project(&[course()], &board(), Some(&motion())).unwrap().secondary().0[0][0] + 1.1).abs() < 1e-14);
    }
    #[test]
    fn source_frames_retain_a_primary_invisible_mode_and_mutual_bridge_compliance() {
        use super::super::board_geometry::ritz::{RitzOptions, bridge_basis};
        let spec = Specification::read(&card(ROW), &[course()]).unwrap();
        let frames = spec.source_ports(&[course()]).unwrap();
        assert_eq!(frames.len(), 1); assert_eq!(frames[0].midi, 69);
        let mesh = fs_plate::ShellMesh::new(vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]],
            vec![[0,1,2]]).unwrap();
        // Source mode 1 moves only along the supplied lateral direction. A
        // primary-only Ritz basis cannot see it, even with spare capacity.
        let source = [[[0.,0.,1.,0.,2.,0.]; 3], [[1.,0.,0.,0.,0.,0.]; 3]];
        let ports: Vec<Vec<f64>> = [frames[0].primary, frames[0].secondary].iter().map(|port| {
            let projection = port.prepare(&mesh, false).unwrap();
            source.iter().map(|mode| projection.project_nodal(*mode).unwrap().0).collect()
        }).collect();
        assert_eq!(ports[0], [1., 0.]);
        assert!((ports[1][0]-0.04).abs() < 1e-15);
        assert_eq!(ports[1][1], 1.);
        let lambda = [100_f64, 200.].map(|hz| (std::f64::consts::TAU*hz).powi(2));
        let damping = [0.3, 0.4];
        let options = RitzOptions::parse("2,0,150").unwrap();
        let primary_only = bridge_basis(&lambda, &damping, &ports[..1], &options).unwrap();
        assert_eq!(primary_only.columns, vec![vec![1.,0.]]);
        let both = bridge_basis(&lambda, &damping, &ports, &options).unwrap();
        assert_eq!(both.columns.len(), 2);
        assert!(both.max_relative_snapshot_error < 1e-14);
        let reduced: Vec<Vec<f64>> = ports.iter().map(|port| both.columns.iter().map(|q|
            q.iter().zip(port).map(|(q,b)| q*b).sum()).collect()).collect();
        let k = &both.stiffness;
        let determinant = k[0]*k[3]-k[1]*k[2];
        assert!(determinant > 0.);
        let solve = |force: &[f64]| [(k[3]*force[0]-k[1]*force[1])/determinant,
            (k[0]*force[1]-k[2]*force[0])/determinant];
        // Independent 2x2 solve checks both directions of reciprocal transfer
        // and lateral self-compliance, including the previously unseen mode.
        for drive in 0..2 { for receive in 0..2 {
            let displacement = solve(&reduced[drive]);
            let response: f64 = reduced[receive].iter().zip(displacement).map(|(b,q)| b*q).sum();
            let exact: f64 = (0..2).map(|i| ports[drive][i]*ports[receive][i]/lambda[i]).sum();
            assert!((response-exact).abs() < 1e-12*exact.abs());
        } }
    }
    #[test]
    fn incomplete_or_nonphysical_frames_and_unavailable_motion_refuse_without_fallback() {
        for row in [ROW.replace("0.2,0.3,0.5", "0.2,0.3,0.6"),
            ROW.replace("0,1,0,0,0,1", "0,2,0,0,0,1"),
            ROW.replace("0,1,0,0,0,1", "0,1,0,0,1,0"),
            ROW.replace("0.4", "NaN"), ROW.replace("0.4", "10.1"), format!("{ROW}\n{ROW}")] {
            assert!(Specification::read(&card(&row), &[course()]).is_err(), "accepted {row}");
        }
        let mut other = course(); other.midi = 70;
        assert!(Specification::read(&card(ROW), &[course(), other]).is_err());
        let spec = Specification::read(&card(ROW), &[course()]).unwrap();
        assert!(spec.project(&[course()], &board(), None).err().unwrap().contains("full-vector"));
        let mut wrong = board(); wrong[0].bridge[48] = 4.;
        assert!(spec.project(&[course()], &wrong, Some(&motion())).err().unwrap().contains("does not match"));
        // A near-node coefficient cannot evade the primary geometry check by
        // falling below a fixed absolute displacement threshold.
        let mut tiny = motion();
        for value in tiny.shapes.iter_mut().flatten().flatten() { *value *= 1e-12; }
        wrong[0].bridge.fill(3e-12);
        assert!(spec.project(&[course()], &wrong, Some(&tiny)).is_ok());
        wrong[0].bridge[48] = 4e-12;
        assert!(spec.project(&[course()], &wrong, Some(&tiny)).is_err());
        let wrong_site = Specification::read(&card(&ROW.replace("69,0,", "69,1,")), &[course()]).unwrap();
        assert!(wrong_site.project(&[course()], &board(), Some(&motion())).is_err());
    }
}
