//! SGS stress at the original reciprocal MAC shear edges.
//!
//! Cell-centred differentiation can annihilate an alternating shear mode.
//! At each four-fluid-cell edge, the two reciprocal shear derivatives instead
//! use the actual adjacent staggered velocities. Other gradient entries and
//! the authored filter-width field are interpolated to the same point using
//! physical, nonuniform centre-to-edge distances. One nonnegative viscosity
//! is cached for both momentum components; wall and incomplete stencils keep
//! the original boundary treatment.

use super::{
    SgsModel, Solver, State, admit_local_viscosity, eddy_viscosity, refuse, vreman_eddy_viscosity,
};
use crate::conjugate::{ChtError, poll};
use fs_exec::CancelGate;

fn axes(third: usize) -> (usize, usize) {
    match third {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    }
}

fn dimensions(n: [usize; 3], third: usize) -> [usize; 3] {
    core::array::from_fn(|axis| n[axis] - usize::from(axis != third))
}

fn index(n: [usize; 3], third: usize, upper: [usize; 3]) -> usize {
    let dims = dimensions(n, third);
    let at: [usize; 3] = core::array::from_fn(|axis| upper[axis] - usize::from(axis != third));
    (at[2] * dims[1] + at[1]) * dims[0] + at[0]
}

fn coordinates(n: [usize; 3], third: usize, slot: usize) -> [usize; 3] {
    let dims = dimensions(n, third);
    let yz = slot / dims[0];
    let at = [slot % dims[0], yz % dims[1], yz / dims[1]];
    core::array::from_fn(|axis| at[axis] + usize::from(axis != third))
}

fn cells(upper: [usize; 3], a: usize, d: usize) -> [[usize; 3]; 4] {
    let mut low = upper;
    low[a] -= 1;
    low[d] -= 1;
    let mut along_a = low;
    along_a[a] += 1;
    let mut along_d = low;
    along_d[d] += 1;
    [low, along_a, along_d, upper]
}

pub(super) fn storage(n: [usize; 3]) -> [Vec<f64>; 3] {
    core::array::from_fn(|third| {
        let dims = dimensions(n, third);
        // Negative is an absent stencil, never a physical viscosity.
        vec![-1.0; dims[0] * dims[1] * dims[2]]
    })
}

impl Solver<'_> {
    pub(super) fn update_sgs_edges(
        &self,
        state: &mut State,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        for third in 0..3 {
            let (a, d) = axes(third);
            for slot in 0..state.shear_eddy[third].len() {
                if slot % 64 == 0 {
                    poll(gate)?;
                }
                let upper = coordinates(self.n, third, slot);
                let at = cells(upper, a, d);
                let ids = at.map(|c| super::super::cell_of(self.n, c));
                if ids.iter().any(|&cell| !self.domain.is_fluid(cell)) {
                    state.shear_eddy[third][slot] = -1.0;
                    continue;
                }
                let wa = [
                    self.domain.width(a, upper[a] - 1),
                    self.domain.width(a, upper[a]),
                ];
                let wd = [
                    self.domain.width(d, upper[d] - 1),
                    self.domain.width(d, upper[d]),
                ];
                // Linear interpolation to the edge weights each cell by the
                // OPPOSITE centre distance, including on a graded mesh.
                let fa = wa[1] / (wa[0] + wa[1]);
                let fd = wd[1] / (wd[0] + wd[1]);
                let weights = [
                    fa * fd,
                    (1.0 - fa) * fd,
                    fa * (1.0 - fd),
                    (1.0 - fa) * (1.0 - fd),
                ];
                let mut gradient = [[0.0; 3]; 3];
                for (cell, weight) in ids.into_iter().zip(weights) {
                    for (row, source) in gradient.iter_mut().zip(state.gradients[cell]) {
                        for (entry, value) in row.iter_mut().zip(source) {
                            *entry += weight * value;
                        }
                    }
                }
                let mut below = upper;
                below[d] -= 1;
                gradient[a][d] = (self.vel[a][self.comps[a].index(upper)]
                    - self.vel[a][self.comps[a].index(below)])
                    / (0.5 * wd[0] + 0.5 * wd[1]);
                let mut left = upper;
                left[a] -= 1;
                gradient[d][a] = (self.vel[d][self.comps[d].index(upper)]
                    - self.vel[d][self.comps[d].index(left)])
                    / (0.5 * wa[0] + 0.5 * wa[1]);
                let eddy = match state.config.model {
                    SgsModel::Vreman { cv } => {
                        let widths: [f64; 3] = core::array::from_fn(|axis| {
                            ids.into_iter()
                                .zip(weights)
                                .map(|(cell, weight)| {
                                    weight
                                        * (self.domain.widths(cell)[axis]
                                            * state.config.filter_width_factor)
                                })
                                .sum()
                        });
                        vreman_eddy_viscosity(cv, widths, gradient)
                    }
                    model => {
                        let width = ids
                            .into_iter()
                            .zip(weights)
                            .map(|(cell, weight)| weight * state.filter_width[cell])
                            .sum();
                        eddy_viscosity(model, width, gradient)
                    }
                }
                .map_err(|error| refuse(format!("shear edge {upper:?}, axes {a}/{d}: {error}")))?;
                // Each adjacent material supplies a local physical cap. A
                // harmonic average cannot hide a too-small molecular value.
                for cell in ids {
                    let molecular = self
                        .molecular_viscosity
                        .get(cell)
                        .copied()
                        .unwrap_or(self.mu);
                    admit_local_viscosity(
                        cell,
                        eddy,
                        molecular,
                        self.rho,
                        state.config.maximum_eddy_to_molecular,
                    )?;
                }
                state.shear_eddy[third][slot] = eddy;
            }
        }
        Ok(())
    }

    /// Borrow the same coefficient from either reciprocal component row.
    pub(in crate::conjugate::simple) fn sgs_edge_eddy(
        &self,
        minus: [usize; 3],
        plus: [usize; 3],
        d: usize,
        up: bool,
    ) -> Option<f64> {
        let state = self.sgs.as_ref()?;
        let a = (0..3).find(|&axis| minus[axis] != plus[axis])?;
        if a == d {
            return None;
        }
        let third = 3 - a - d;
        let mut upper = plus;
        upper[d] += usize::from(up);
        if upper[a] == 0 || upper[a] >= self.n[a] || upper[d] == 0 || upper[d] >= self.n[d] {
            return None;
        }
        let eddy = state.shear_eddy[third][index(self.n, third, upper)];
        (eddy >= 0.0).then_some(eddy)
    }

    pub(super) fn admit_sgs_edge_material(
        &self,
        state: &State,
        molecular: &[f64],
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        for third in 0..3 {
            let (a, d) = axes(third);
            for (slot, &eddy) in state.shear_eddy[third].iter().enumerate() {
                if slot % 64 == 0 {
                    poll(gate)?;
                }
                if eddy < 0.0 {
                    continue;
                }
                let upper = coordinates(self.n, third, slot);
                for at in cells(upper, a, d) {
                    let cell = super::super::cell_of(self.n, at);
                    admit_local_viscosity(
                        cell,
                        eddy,
                        molecular[cell],
                        self.rho,
                        state.config.maximum_eddy_to_molecular,
                    )?;
                }
            }
        }
        Ok(())
    }
}
