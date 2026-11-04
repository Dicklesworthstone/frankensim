//! Axial stretch of the original eccentric, cubic Euler--Bernoulli beam field.
//!
//! The original linear pencil already contains EA*u'^2/2. Add only the
//! polynomial difference to EA*(u' + (u'^2+v'^2+w'^2)/2)^2/2, so the beam
//! participates in the same equilibrium as its shell and its tension or
//! compression enters the accepted tangent exactly once. Reference-frame
//! bending, torsion and rigid arms retain their original moderate-rotation law.

use super::{CoordinateRows, ShellBeam, ShellMesh, ShellModel};
use crate::shell::stiffened::element_kinematics;
use fs_sparse::Coo;

pub(super) const MAX_BONDED_BEAMS: usize = 40_000;

// Degree-eight axial storage and its derivatives are integrated exactly.
const GAUSS: [(f64, f64); 5] = [
    (0.046_910_077_030_668, 0.118_463_442_528_094_54),
    (0.230_765_344_947_158_45, 0.239_314_335_249_683_24),
    (0.5, 0.284_444_444_444_444_44),
    (0.769_234_655_052_841_5, 0.239_314_335_249_683_24),
    (0.953_089_922_969_332, 0.118_463_442_528_094_54),
];

struct Sample {
    // Derivatives of the local u,v,w fields with respect to shell-node DOFs.
    rows: [[f64; 12]; 3],
    ea_weight: f64,
}

pub(super) struct BeamStretch {
    dofs: [Option<usize>; 12],
    pub(super) transform: Option<CoordinateRows>,
    samples: [Sample; 5],
    // Quadratic Bernstein controls bound the gradient on the WHOLE segment,
    // including endpoints and extrema between the integration points.
    controls: [[[f64; 12]; 3]; 3],
}

fn apply(rows: &[[f64; 12]; 3], u: &[f64; 12]) -> [f64; 3] {
    std::array::from_fn(|c| rows[c].iter().zip(u).map(|(a, b)| a * b).sum())
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

impl BeamStretch {
    pub(super) fn new(
        mesh: &ShellMesh,
        model: &ShellModel,
        beam: &ShellBeam,
    ) -> Result<Self, String> {
        let (length, frame) = element_kinematics(mesh, beam).map_err(|e| e.to_string())?;
        let ea_length = beam.section.young_pa * beam.section.area_m2 * length;
        if !ea_length.is_finite() || ea_length <= 0. {
            return Err("unrepresentable bonded-beam axial storage".into());
        }
        let gradient = |s: f64| {
            let mut local = [[0.; 12]; 3];
            local[0][0] = -1. / length;
            local[0][6] = 1. / length;
            let h = [
                (6. * s * s - 6. * s) / length,
                3. * s * s - 4. * s + 1.,
                (-6. * s * s + 6. * s) / length,
                3. * s * s - 2. * s,
            ];
            for i in 0..4 {
                local[1][[1, 5, 7, 11][i]] = h[i];
                local[2][[2, 4, 8, 10][i]] = h[i] * [1., -1., 1., -1.][i];
            }
            std::array::from_fn(|c| {
                std::array::from_fn(|j| (0..12).map(|i| local[c][i] * frame[12 * i + j]).sum())
            })
        };
        let ends: [[[f64; 12]; 3]; 3] = [gradient(0.), gradient(0.5), gradient(1.)];
        let controls = [
            ends[0],
            std::array::from_fn(|c| {
                std::array::from_fn(|i| 2. * ends[1][c][i] - 0.5 * (ends[0][c][i] + ends[2][c][i]))
            }),
            ends[2],
        ];
        let samples = GAUSS.map(|(s, weight)| Sample {
            rows: gradient(s),
            ea_weight: ea_length * weight,
        });
        if controls.iter().flatten().flatten().any(|v| !v.is_finite())
            || samples
                .iter()
                .any(|s| s.rows.iter().flatten().any(|v| !v.is_finite()))
        {
            return Err("bonded-beam strain map overflow".into());
        }
        Ok(Self {
            dofs: std::array::from_fn(|i| model.dof_map[6 * beam.nodes[i / 6] + i % 6]),
            transform: None,
            samples,
            controls,
        })
    }

    fn for_dof(&self, i: usize, mut f: impl FnMut(usize, f64)) {
        if let Some(local) = self.dofs[i] {
            if let Some(rows) = &self.transform {
                for &(global, weight) in &rows[local] {
                    f(global, weight);
                }
            } else {
                f(local, 1.);
            }
        }
    }

    fn values(&self, u: &[f64]) -> [f64; 12] {
        std::array::from_fn(|i| {
            let mut value = 0.;
            self.for_dof(i, |d, w| value += w * u[d]);
            value
        })
    }

    fn admitted_values(&self, u: &[f64], limit: f64) -> Option<[f64; 12]> {
        let values = self.values(u);
        self.controls
            .iter()
            .all(|rows| {
                let g = apply(rows, &values);
                let bound = g[0].hypot(g[1]).hypot(g[2]);
                bound.is_finite() && bound <= limit
            })
            .then_some(values)
    }

    pub(super) fn add_force_energy(&self, u: &[f64], force: &mut [f64], limit: f64) -> Option<f64> {
        let values = self.admitted_values(u, limit)?;
        let mut energy = 0.;
        let mut local_force = [0.; 12];
        for sample in &self.samples {
            let g = apply(&sample.rows, &values);
            let quadratic = 0.5 * dot(g, g);
            let strain = g[0] + quadratic;
            energy += sample.ea_weight * (g[0] * quadratic + 0.5 * quadratic * quadratic);
            for (i, force) in local_force.iter_mut().enumerate() {
                let nonlinear = (0..3).map(|c| g[c] * sample.rows[c][i]).sum::<f64>();
                *force += sample.ea_weight * (sample.rows[0][i] * quadratic + nonlinear * strain);
            }
        }
        for (i, f) in local_force.into_iter().enumerate() {
            self.for_dof(i, |d, w| force[d] += w * f);
        }
        energy.is_finite().then_some(energy)
    }

    // O(5*12) local action: Newton--Krylov never builds a 12x12 beam tangent or
    // expands a dense mortar transform just to multiply it by a direction.
    pub(super) fn add_action(
        &self,
        u: &[f64],
        v: &[f64],
        out: &mut [f64],
        limit: f64,
    ) -> Option<()> {
        let values = self.admitted_values(u, limit)?;
        let direction = self.values(v);
        let mut local_action = [0.; 12];
        for sample in &self.samples {
            let g = apply(&sample.rows, &values);
            let dg = apply(&sample.rows, &direction);
            let strain = g[0] + 0.5 * dot(g, g);
            let dq = dot(g, dg);
            let de = dg[0] + dq;
            for (i, action) in local_action.iter_mut().enumerate() {
                let nonlinear = (0..3).map(|c| g[c] * sample.rows[c][i]).sum::<f64>();
                let derivative = (0..3).map(|c| dg[c] * sample.rows[c][i]).sum::<f64>();
                *action += sample.ea_weight
                    * (sample.rows[0][i] * dq + derivative * strain + nonlinear * de);
            }
        }
        if local_action.iter().any(|v| !v.is_finite()) {
            return None;
        }
        for (i, f) in local_action.into_iter().enumerate() {
            self.for_dof(i, |d, w| out[d] += w * f);
        }
        Some(())
    }

    pub(super) fn add_tangent(&self, u: &[f64], k: &mut Coo, limit: f64) -> Option<()> {
        let values = self.admitted_values(u, limit)?;
        let mut tangent = [0.; 144];
        for sample in &self.samples {
            let g = apply(&sample.rows, &values);
            let strain = g[0] + 0.5 * dot(g, g);
            let nonlinear: [f64; 12] =
                std::array::from_fn(|i| (0..3).map(|c| g[c] * sample.rows[c][i]).sum());
            for i in 0..12 {
                for j in i..12 {
                    let metric = (0..3)
                        .map(|c| sample.rows[c][i] * sample.rows[c][j])
                        .sum::<f64>();
                    tangent[12 * i + j] += sample.ea_weight
                        * (sample.rows[0][i] * nonlinear[j]
                            + nonlinear[i] * sample.rows[0][j]
                            + nonlinear[i] * nonlinear[j]
                            + strain * metric);
                }
            }
        }
        if tangent.iter().any(|v| !v.is_finite()) {
            return None;
        }
        for i in 0..12 {
            for j in i..12 {
                let value = tangent[12 * i + j];
                if value == 0. {
                    continue;
                }
                // A shared retained coordinate receives both off-diagonal
                // local terms. Mirror only identical floating-point addends.
                self.for_dof(i, |a, wa| {
                    self.for_dof(j, |b, wb| {
                        if i == j && a > b {
                            return;
                        }
                        let entry = wa * value * wb;
                        if a == b {
                            k.push(a, a, if i == j { entry } else { 2. * entry });
                        } else {
                            k.push(a, b, entry);
                            k.push(b, a, entry);
                        }
                    })
                });
            }
        }
        Some(())
    }
}

#[cfg(test)]
mod tests;
