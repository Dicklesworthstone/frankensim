//! SGS at reservoir shear half-edges, using the native traction-mouth stencil.
//!
//! Each side has its own cache even for a single layer of cells. The normal
//! velocity difference along the mouth and the zero tangential trace at half
//! a cell supply the two reciprocal derivatives used by momentum and power.
//! Other tensor entries and the authored cell filter widths use nonuniform
//! centre-to-edge interpolation over the physical half-edge control volume,
//! as for the original interior MAC edges. This is a volume reconstruction,
//! not an exact gradient of the tangential trace on the reservoir surface.

use super::{
    SgsModel, Solver, State, admit_local_viscosity, eddy_viscosity, refuse, vreman_eddy_viscosity,
};
use crate::conjugate::{ChtError, poll};
use fs_exec::CancelGate;

pub(super) fn storage(n: [usize; 3], normal: Option<usize>) -> Result<[Vec<f64>; 3], ChtError> {
    let mut values = core::array::from_fn(|_| Vec::new());
    if let Some(normal) = normal {
        for tangent in 0..3 {
            if tangent == normal {
                continue;
            }
            let third = 3 - normal - tangent;
            let count = (n[tangent] - 1)
                .checked_mul(n[third])
                .and_then(|value| value.checked_mul(2))
                .ok_or_else(|| refuse("reservoir SGS edge count overflow"))?;
            values[tangent]
                .try_reserve_exact(count)
                .map_err(|_| refuse("reservoir SGS edge allocation failed"))?;
            // Negative marks an absent fluid stencil, never a viscosity.
            values[tangent].resize(count, -1.0);
        }
    }
    Ok(values)
}

fn cells(n: [usize; 3], normal: usize, tangent: usize, high: bool, slot: usize) -> [[usize; 3]; 2] {
    let third = 3 - normal - tangent;
    let mut plus = [0; 3];
    plus[normal] = if high { n[normal] - 1 } else { 0 };
    plus[tangent] = slot % (n[tangent] - 1) + 1;
    plus[third] = slot / (n[tangent] - 1);
    let mut minus = plus;
    minus[tangent] -= 1;
    [minus, plus]
}

impl Solver<'_> {
    pub(super) fn update_sgs_pressure_mouths(
        &self,
        state: &mut State,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        let Some(normal) = self.pressure_work_axis else {
            return Ok(());
        };
        for tangent in 0..3 {
            if tangent == normal {
                continue;
            }
            let per_side = state.mouth_eddy[tangent].len() / 2;
            for side in 0..2 {
                let high = side == 1;
                for slot in 0..per_side {
                    if slot % 64 == 0 {
                        poll(gate)?;
                    }
                    let [minus, plus] = cells(self.n, normal, tangent, high, slot);
                    let ids = [minus, plus].map(|at| super::super::cell_of(self.n, at));
                    let index = side * per_side + slot;
                    if ids.iter().any(|&cell| !self.domain.is_fluid(cell)) {
                        state.mouth_eddy[tangent][index] = -1.0;
                        continue;
                    }
                    let wm = self.domain.width(tangent, minus[tangent]);
                    let wp = self.domain.width(tangent, plus[tangent]);
                    let fraction = wp / (wm + wp);
                    let weights = [fraction, 1.0 - fraction];
                    let mut gradient = [[0.0; 3]; 3];
                    for (cell, weight) in ids.into_iter().zip(weights) {
                        for (row, source) in gradient.iter_mut().zip(state.gradients[cell]) {
                            for (entry, value) in row.iter_mut().zip(source) {
                                *entry += weight * value;
                            }
                        }
                    }
                    gradient[normal][tangent] = (self.vel[normal]
                        [self.cell_face(normal, plus, high)]
                        - self.vel[normal][self.cell_face(normal, minus, high)])
                        / (0.5 * (wm + wp));
                    gradient[tangent][normal] = (if high { -1.0 } else { 1.0 })
                        * self.vel[tangent][self.cell_face(tangent, minus, true)]
                        / (0.5 * self.domain.width(normal, plus[normal]));
                    let eddy = match state.config.model {
                        SgsModel::Vreman { cv } => {
                            let widths = core::array::from_fn(|axis| {
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
                    .map_err(|error| {
                        refuse(format!(
                            "reservoir edge {plus:?}, axes {normal}/{tangent}, high={high}: {error}"
                        ))
                    })?;
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
                            self.sgs_density(cell),
                            state.config.maximum_eddy_to_molecular,
                        )?;
                    }
                    state.mouth_eddy[tangent][index] = eddy;
                }
            }
        }
        Ok(())
    }

    /// One cached coefficient for both reciprocal momentum rows and their
    /// power partition. The side is explicit because one cell may touch both
    /// reservoir mouths with different velocity gradients.
    pub(in crate::conjugate::simple) fn pressure_mouth_sgs_eddy(
        &self,
        tangent: usize,
        minus: [usize; 3],
        plus: [usize; 3],
        high: bool,
    ) -> Option<f64> {
        let normal = self.pressure_work_axis?;
        let state = self.sgs.as_ref()?;
        if tangent == normal || plus[tangent] == 0 || plus[tangent] >= self.n[tangent] {
            return None;
        }
        let third = 3 - normal - tangent;
        let plane = if high { self.n[normal] - 1 } else { 0 };
        if minus[normal] != plane
            || plus[normal] != plane
            || minus[tangent] + 1 != plus[tangent]
            || minus[third] != plus[third]
        {
            return None;
        }
        let per_side = state.mouth_eddy[tangent].len() / 2;
        let index =
            usize::from(high) * per_side + plus[third] * (self.n[tangent] - 1) + plus[tangent] - 1;
        let eddy = *state.mouth_eddy[tangent].get(index)?;
        (eddy >= 0.0).then_some(eddy)
    }

    pub(super) fn admit_sgs_pressure_material(
        &self,
        state: &State,
        molecular: &[f64],
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        let Some(normal) = self.pressure_work_axis else {
            return Ok(());
        };
        for tangent in 0..3 {
            if tangent == normal {
                continue;
            }
            let per_side = state.mouth_eddy[tangent].len() / 2;
            for side in 0..2 {
                for slot in 0..per_side {
                    if slot % 64 == 0 {
                        poll(gate)?;
                    }
                    let eddy = state.mouth_eddy[tangent][side * per_side + slot];
                    if eddy < 0.0 {
                        continue;
                    }
                    for at in cells(self.n, normal, tangent, side == 1, slot) {
                        let cell = super::super::cell_of(self.n, at);
                        admit_local_viscosity(
                            cell,
                            eddy,
                            molecular[cell],
                            self.sgs_density(cell),
                            state.config.maximum_eddy_to_molecular,
                        )?;
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
