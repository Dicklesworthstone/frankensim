use fs_math::det;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::{ConstraintFrame, NodalDirectors, ShearShellSection, ShellConstraint};

pub const DISTORTION_KEY: u64 = 0x4952_5053_4845_4c4c;

#[derive(Clone, Copy, Debug)]
pub enum Case {
    Roof,
    Cylinder,
    Hole,
    Full,
}
impl Case {
    pub fn name(self) -> &'static str {
        match self {
            Self::Roof => "scordelis-lo",
            Self::Cylinder => "pinched-cylinder",
            Self::Hole => "hemisphere-18deg-hole",
            Self::Full => "full-hemisphere",
        }
    }
    pub fn reference(self) -> f64 {
        match self {
            Self::Roof => 0.3024,
            Self::Cylinder => 1.8248e-5,
            Self::Hole => 0.094,
            Self::Full => 0.0924,
        }
    }
    pub fn rungs(self) -> &'static [usize] {
        match self {
            Self::Cylinder => &[4, 8, 16, 32, 64],
            _ => &[4, 8, 16, 32],
        }
    }
}
pub struct Fixture {
    pub mesh: ShellMesh,
    pub directors: NodalDirectors,
    pub sections: Vec<ShearShellSection>,
    pub constraints: Vec<ShellConstraint>,
    pub load: Vec<f64>,
    pub observation: usize,
    pub component: usize,
    pub sign: f64,
    pub e: f64,
    pub nu: f64,
    pub h: f64,
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|k| a[k] - b[k])
}

fn parameter(i: usize, j: usize, n: usize, distorted: bool) -> (f64, f64) {
    let mut a = i as f64 / n as f64;
    let mut b = j as f64 / n as f64;
    if distorted && i > 0 && i < n && j > 0 && j < n {
        let draw = fs_rand::Stream::at(
            fs_rand::StreamKey {
                seed: DISTORTION_KEY,
                kernel: 47,
                tile: 0,
            },
            (j * (n + 1) + i) as u64,
        );
        a += 0.15 / n as f64 * (draw[0] as f64 / u32::MAX as f64 - 0.5);
        b += 0.15 / n as f64 * (draw[1] as f64 / u32::MAX as f64 - 0.5);
    }
    (a, b)
}
fn triangle(
    tris: &mut Vec<[usize; 3]>,
    nodes: &[[f64; 3]],
    normals: &[[f64; 3]],
    mut t: [usize; 3],
) {
    if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
        return;
    }
    let n = cross(sub(nodes[t[1]], nodes[t[0]]), sub(nodes[t[2]], nodes[t[0]]));
    let direction: [f64; 3] =
        core::array::from_fn(|k| normals[t[0]][k] + normals[t[1]][k] + normals[t[2]][k]);
    if dot(n, direction) < 0.0 {
        t.swap(1, 2);
    }
    tris.push(t);
}

pub fn build(case: Case, n: usize, distorted: bool) -> Fixture {
    match case {
        Case::Roof | Case::Cylinder => cylinder(case, n, distorted),
        _ => hemisphere(case, n, distorted),
    }
}
fn finish(
    nodes: Vec<[f64; 3]>,
    normals: Vec<[f64; 3]>,
    tris: Vec<[usize; 3]>,
    masks: Vec<[bool; 6]>,
    load: Vec<f64>,
    observation: usize,
    component: usize,
    sign: f64,
    e: f64,
    nu: f64,
    h: f64,
) -> Fixture {
    let mesh = ShellMesh::new(nodes, tris).unwrap();
    let directors = NodalDirectors::new(&mesh, normals, 0.5).unwrap();
    let shear = e / (2.0 * (1.0 + nu));
    // rho=1 is an explicit inertia-only analytic fixture choice. Published
    // static loads below are supplied independently; this is not a density
    // measurement or a modal validation claim.
    let section = ShearShellSection::from_orthotropic_solid(
        e,
        e,
        nu,
        shear,
        shear,
        shear,
        h,
        1.0,
        0.0,
        5.0 / 6.0,
    )
    .unwrap();
    let sections = vec![section; mesh.tris.len()];
    let constraints = masks
        .into_iter()
        .enumerate()
        .filter(|(_, mask)| mask.iter().any(|v| *v))
        .map(|(node, mask)| ShellConstraint {
            node,
            mask,
            frame: ConstraintFrame::Global,
        })
        .collect();
    Fixture {
        mesh,
        directors,
        sections,
        constraints,
        load,
        observation,
        component,
        sign,
        e,
        nu,
        h,
    }
}
fn cylinder(case: Case, n: usize, distorted: bool) -> Fixture {
    let roof = matches!(case, Case::Roof);
    let (radius, length, half_angle, e, nu, h) = if roof {
        (
            25.0,
            50.0,
            40.0 * core::f64::consts::PI / 180.0,
            4.32e8,
            0.0,
            0.25,
        )
    } else {
        (300.0, 600.0, core::f64::consts::FRAC_PI_2, 3e6, 0.3, 3.0)
    };
    let mut nodes = Vec::new();
    let mut normals = Vec::new();
    let mut masks = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let (a, b) = parameter(i, j, n, distorted);
            let angle = a * half_angle;
            let normal = [det::sin(angle), 0.0, det::cos(angle)];
            nodes.push([radius * normal[0], 0.5 * length * b, radius * normal[2]]);
            normals.push(normal);
            let mut mask = [false; 6];
            // x=0 symmetry: ux=theta_y=theta_z=0.
            if i == 0 {
                mask[0] = true;
                mask[4] = true;
                mask[5] = true;
            }
            // y=0 mid-length symmetry: uy=theta_x=theta_z=0.
            if j == 0 {
                mask[1] = true;
                mask[3] = true;
                mask[5] = true;
            }
            // End diaphragm restrains the two radial-plane displacements and
            // rotation about the cylinder axis, while axial displacement is free.
            if j == n {
                mask[0] = true;
                mask[2] = true;
                mask[4] = true;
            }
            if !roof && i == n {
                mask[2] = true;
                mask[3] = true;
                mask[4] = true;
            }
            masks.push(mask);
        }
    }
    let mut tris = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let a = j * (n + 1) + i;
            let b = a + 1;
            let c = a + n + 1;
            let d = c + 1;
            triangle(&mut tris, &nodes, &normals, [a, b, d]);
            triangle(&mut tris, &nodes, &normals, [a, d, c]);
        }
    }
    let mut load = vec![0.0; 6 * nodes.len()];
    let observation = if roof { n } else { 0 };
    if roof {
        for t in &tris {
            let normal = cross(sub(nodes[t[1]], nodes[t[0]]), sub(nodes[t[2]], nodes[t[0]]));
            let area = 0.5 * det::sqrt(dot(normal, normal));
            for &i in t {
                load[6 * i + 2] -= 90.0 * area / 3.0;
            }
        }
    } else {
        // The load point lies on two symmetry planes: one quarter of the
        // full physical P=1 top pinching force enters this one-eighth model.
        load[6 * observation + 2] = -0.25;
    }
    finish(
        nodes,
        normals,
        tris,
        masks,
        load,
        observation,
        2,
        -1.0,
        e,
        nu,
        h,
    )
}

fn hemisphere(case: Case, n: usize, distorted: bool) -> Fixture {
    let hole = matches!(case, Case::Hole);
    let radius = 10.0;
    let polar_start = if hole {
        18.0 * core::f64::consts::PI / 180.0
    } else {
        0.0
    };
    let mut nodes = Vec::new();
    let mut normals = Vec::new();
    let mut masks = Vec::<[bool; 6]>::new();
    let mut ids = vec![vec![0; n + 1]; n + 1];
    for j in 0..=n {
        for i in 0..=n {
            let (a, b) = parameter(i, j, n, distorted);
            let phi = core::f64::consts::FRAC_PI_2 * a;
            let theta = polar_start + (core::f64::consts::FRAC_PI_2 - polar_start) * b;
            let normal = [
                det::sin(theta) * det::cos(phi),
                det::sin(theta) * det::sin(phi),
                det::cos(theta),
            ];
            let id = if !hole && j == 0 && i > 0 {
                0
            } else {
                let id = nodes.len();
                nodes.push(normal.map(|v| radius * v));
                normals.push(normal);
                masks.push([false; 6]);
                id
            };
            ids[j][i] = id;
            if i == 0 {
                masks[id][1] = true;
                masks[id][3] = true;
                masks[id][5] = true;
            }
            if i == n {
                masks[id][0] = true;
                masks[id][4] = true;
                masks[id][5] = true;
            }
        }
    }
    let mut tris = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (ids[j][i], ids[j][i + 1], ids[j + 1][i], ids[j + 1][i + 1]);
            triangle(&mut tris, &nodes, &normals, [a, b, d]);
            triangle(&mut tris, &nodes, &normals, [a, d, c]);
        }
    }
    let observation = ids[n][0];
    let other = ids[n][n];
    masks[observation][2] = true;
    let mut load = vec![0.0; 6 * nodes.len()];
    // Each equatorial load is bisected by one symmetry plane; F=2 on the
    // complete body gives unit forces on the two quarter-model boundary nodes.
    load[6 * observation] = 1.0;
    load[6 * other + 1] = -1.0;
    finish(
        nodes,
        normals,
        tris,
        masks,
        load,
        observation,
        0,
        1.0,
        6.825e7,
        0.3,
        0.04,
    )
}
