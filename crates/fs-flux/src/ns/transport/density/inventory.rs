//! Bounded specific traces on the existing conservative directional mass flux.

use super::{
    DensityTransportCell, DensityTransportError, DensityTransportWorkspace,
    PreparedDensityTransport, Result, Sum, finite, linear, poll, positive_value,
};

impl PreparedDensityTransport<'_> {
    pub(super) fn reconstruct_carried_inventories<const N: usize>(
        &self,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        work: &mut DensityTransportWorkspace<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        if work.specific_deviations.len() != cells.len()
            || work.specific_margins.len() != cells.len()
            || work.specific_limiter.len() != cells.len()
        {
            return Err(DensityTransportError::Invalid(
                "matching bounded inventory reconstruction workspace required",
            ));
        }
        let linear = self.linear.as_ref().expect("bounded density geometry");
        for limiter in &mut work.specific_limiter {
            poll(cancelled)?;
            *limiter = 1.0;
        }
        // Reuse the scalar density-bound buffer component by component. Slopes,
        // retained-inventory margins and one COMMON limiter use warm storage.
        for component in 0..N {
            for bound in &mut work.vertex_bounds {
                poll(cancelled)?;
                *bound = [f64::INFINITY, f64::NEG_INFINITY];
            }
            for (cell, state) in cells.iter().enumerate() {
                poll(cancelled)?;
                for &vertex in &self.mesh.tris[cell] {
                    let bound = &mut work.vertex_bounds[vertex];
                    bound[0] = bound[0].min(state.specific[component]);
                    bound[1] = bound[1].max(state.specific[component]);
                }
            }
            for (edge, exterior) in inflow.iter().enumerate() {
                poll(cancelled)?;
                if self.volume_rates[edge].against_normal_m2_s() > 0.0 {
                    if let Some(state) = exterior {
                        for vertex in [self.mesh.edges[edge].verts.0, self.mesh.edges[edge].verts.1]
                        {
                            let bound = &mut work.vertex_bounds[vertex];
                            bound[0] = bound[0].min(state.specific[component]);
                            bound[1] = bound[1].max(state.specific[component]);
                        }
                    }
                }
            }
            for (cell, state) in cells.iter().enumerate() {
                poll(cancelled)?;
                let value = state.specific[component];
                let bounds = self.mesh.tris[cell].map(|vertex| work.vertex_bounds[vertex]);
                let lower = bounds.iter().fold(value, |a, b| a.min(b[0]));
                let upper = bounds.iter().fold(value, |a, b| a.max(b[1]));
                work.specific_deviations[cell][component] = [0.0; 3];
                work.specific_margins[cell][component] = [0.0; 2];
                if lower == upper || work.specific_limiter[cell] == 0.0 {
                    continue;
                }
                let scale = lower.abs().max(upper.abs());
                let mut differences = [0.0; 3];
                for (difference, &(edge, _)) in
                    differences.iter_mut().zip(&self.mesh.tri_edges[cell])
                {
                    let other = linear::neighbor(self.mesh, cell, edge);
                    let sample = if other != usize::MAX {
                        cells[other].specific[component]
                    } else if self.volume_rates[edge].against_normal_m2_s() > 0.0 {
                        inflow[edge].as_ref().expect("admitted inlet").specific[component]
                    } else {
                        value
                    };
                    let delta = sample - value;
                    *difference = if delta.is_finite() {
                        delta / scale
                    } else {
                        sample / scale - value / scale
                    };
                }
                let variation = linear
                    .specific_deviations(cell, differences, work.vertex_densities[cell])
                    .map(|delta| delta * scale);
                if variation.iter().any(|delta| !delta.is_finite()) {
                    // An unrepresentable slope keeps the already admitted
                    // constant donor for every coupled component in this cell.
                    work.specific_limiter[cell] = 0.0;
                    continue;
                }
                work.specific_deviations[cell][component] = variation;
                let mut theta = work.specific_limiter[cell];
                for (&delta, bound) in variation.iter().zip(bounds) {
                    if delta > 0.0 {
                        theta = theta.min((bound[1] - value) / delta);
                    } else if delta < 0.0 {
                        theta = theta.min((value - bound[0]) / -delta);
                    }
                }
                let mut outgoing = Sum::default();
                let mut correction = Sum::default();
                for &(edge, sign) in &self.mesh.tri_edges[cell] {
                    let direction = usize::from(sign < 0.0);
                    let mass = work.transfers[edge].mass_kg[direction];
                    if mass == 0.0 {
                        continue;
                    }
                    outgoing.add(mass)?;
                    let delta = self.specific_face_deviation(
                        cell,
                        edge,
                        direction,
                        variation,
                        work.vertex_densities[cell],
                    )?;
                    correction.add(finite(mass * delta, "linear carried inventory correction")?)?;
                }
                let retained = positive_value(
                    work.masses[cell] - outgoing.value()?,
                    "positive retained donor mass for reconstruction",
                )?;
                let correction = correction.value()?;
                // The retained specific state is s - theta*C/m_retained.
                // Bound it before transport. Incoming traces are already
                // bounded, so the final state is a positive convex mixture.
                let margins = [
                    margin(retained, value, lower),
                    margin(retained, upper, value),
                ];
                work.specific_margins[cell][component] = margins;
                let allowance = margins[usize::from(correction < 0.0)];
                if correction != 0.0 {
                    theta = theta.min(allowance / correction.abs());
                }
                work.specific_limiter[cell] = theta;
            }
        }
        for cell in 0..cells.len() {
            poll(cancelled)?;
            let mut theta = work.specific_limiter[cell];
            if theta == 0.0 {
                continue;
            }
            if theta < 1.0 {
                theta *= 1.0 - 64.0 * f64::EPSILON;
            }
            let mut carried_faces = [[0.0; N]; 3];
            let mut accepted = false;
            // Rounding a large enthalpy datum can exceed a tiny retained
            // inventory margin even when the ideal slope fits. Check the
            // REPRESENTED corrections before emitting any of this donor's
            // faces. A bounded shared contraction preserves composition;
            // unresolved corrections fall back to the original donor flux.
            for _ in 0..8 {
                let mut fits = true;
                for component in 0..N {
                    let mut total = Sum::default();
                    for (local, &(edge, sign)) in self.mesh.tri_edges[cell].iter().enumerate() {
                        let direction = usize::from(sign < 0.0);
                        let mass = work.transfers[edge].mass_kg[direction];
                        let original = work.transfers[edge].carried[direction][component];
                        let delta = if mass == 0.0 {
                            0.0
                        } else {
                            self.specific_face_deviation(
                                cell,
                                edge,
                                direction,
                                work.specific_deviations[cell][component],
                                work.vertex_densities[cell],
                            )?
                        };
                        let carried = (theta * delta).mul_add(mass, original);
                        if !carried.is_finite() || total.add(carried - original).is_err() {
                            fits = false;
                            break;
                        }
                        carried_faces[local][component] = carried;
                    }
                    if !fits {
                        break;
                    }
                    let correction = total.value()?;
                    let [lower, upper] = work.specific_margins[cell][component];
                    if correction > lower || correction < -upper {
                        fits = false;
                        break;
                    }
                }
                if fits {
                    accepted = true;
                    break;
                }
                theta *= 0.5;
            }
            if !accepted {
                continue;
            }
            for (local, &(edge, sign)) in self.mesh.tri_edges[cell].iter().enumerate() {
                let direction = usize::from(sign < 0.0);
                let mass = work.transfers[edge].mass_kg[direction];
                if mass == 0.0 {
                    continue;
                }
                let receiver = linear::neighbor(self.mesh, cell, edge);
                for component in 0..N {
                    let original = work.transfers[edge].carried[direction][component];
                    let carried = carried_faces[local][component];
                    let correction = finite(carried - original, "accepted inventory correction")?;
                    // Use the ACTUAL represented face correction on both sides.
                    // Common limiting preserves linear composition constraints;
                    // neither final-state clipping nor species normalization occurs.
                    work.transfers[edge].carried[direction][component] = carried;
                    work.specific_change[cell][component].add(-correction)?;
                    if receiver != usize::MAX {
                        work.specific_change[receiver][component].add(correction)?;
                    }
                }
            }
        }
        Ok(())
    }

    fn specific_face_deviation(
        &self,
        cell: usize,
        edge: usize,
        direction: usize,
        variation: [f64; 3],
        density: [f64; 3],
    ) -> Result<f64> {
        let face = &self.mesh.edges[edge];
        let triangle = self.mesh.tris[cell];
        let a = triangle
            .iter()
            .position(|&v| v == face.verts.0)
            .expect("admitted face vertex");
        let b = triangle
            .iter()
            .position(|&v| v == face.verts.1)
            .expect("admitted face vertex");
        let s = self
            .linear
            .as_ref()
            .expect("bounded density geometry")
            .mass_centroid(edge, direction, [density[a], density[b]]);
        if !s.is_finite() || !(0.0..=1.0).contains(&s) {
            return Err(DensityTransportError::Arithmetic(
                "directional mass centroid",
            ));
        }
        let delta = variation[b] - variation[a];
        finite(
            if delta.is_finite() {
                delta.mul_add(s, variation[a])
            } else {
                (1.0 - s) * variation[a] + s * variation[b]
            },
            "mass-weighted specific trace",
        )
    }
}

// Difference first retains close large datums; product first avoids overflowing
// the subtraction of opposite finite values when the physical mass is small.
fn margin(mass: f64, upper: f64, lower: f64) -> f64 {
    let difference = upper - lower;
    if difference.is_finite() {
        mass * difference
    } else {
        mass * upper - mass * lower
    }
}
