//! Prescribed piecewise-constant fluid density on horizontal depth layers.
//!
//! With q = g dot (x - surface), hydrostatic gauge pressure is
//! p(x) = sum_j (rho_j-rho_{j-1}) max(q-|g| depth_j, 0), rho_{-1}=0.
//! The same continuous pressure potential supplies surface traction, its
//! follower tangent, exact endpoint work, and the native discrete gradient.
//! Density may decrease with depth; this prescribes a reservoir, not its
//! stability or fluid motion. No local material density is substituted for
//! the displaced fluid's density.

use fs_exec::Cx;

use super::algorithmic::{self, PressureMean};
use super::free_surface::{self, FaceMoments};
use super::{
    ClosedTetPressureSurface, HyperTetError, TetPressure, checkpoint, finite, invalid,
    signed_volume,
};

/// One fluid layer, extending from this depth to the next layer's top.
/// Depth is measured along gravity from the prescribed free surface.
#[derive(Debug, Clone, Copy)]
pub struct HydrostaticDensityLayer {
    /// Nonnegative top depth [m]; the first layer starts at zero.
    pub top_depth_m: f64,
    /// Positive physical fluid density [kg/m^3].
    pub density_kg_m3: f64,
}

#[derive(Debug, Clone, Copy)]
struct Ramp {
    offset_m2_s2: f64,
    density_jump_kg_m3: f64,
}

/// Immutable hydrostatic loading in a prescribed stratified reservoir.
///
/// The last layer extends indefinitely. Each layer has constant positive
/// density; pressure is continuous across every interface and zero above the
/// free surface. Both increasing and decreasing density profiles are admitted.
/// The profile does not evolve a free surface, include fluid inertia, or
/// certify the stability of a density inversion.
///
/// Preparation compresses equal neighboring densities. Evaluation reuses the
/// original clipped-simplex moments, with no volume mesh or material query
/// per layer. The same reference mesh and current deformation remain owned by
/// `ClosedTetPressureSurface`.
#[derive(Debug, Clone)]
pub struct StratifiedHydrostaticPressure {
    surface_point_m: [f64; 3],
    gravity_m_s2: [f64; 3],
    ramps: Vec<Ramp>,
}

impl StratifiedHydrostaticPressure {
    /// Admit SI data and a caller-supplied layer-count budget before allocation.
    /// Zero gravity is admitted and produces zero gauge pressure.
    pub fn new(
        surface_point_m: [f64; 3],
        gravity_m_s2: [f64; 3],
        layers: &[HydrostaticDensityLayer],
        maximum_layers: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, HyperTetError> {
        checkpoint(cx)?;
        if layers.is_empty() || layers.len() > maximum_layers {
            return Err(invalid(
                "hydrostatic density layers are empty or exceed their budget",
            ));
        }
        for value in surface_point_m.into_iter().chain(gravity_m_s2) {
            finite(value)?;
        }
        let gravity = finite(
            gravity_m_s2[0]
                .hypot(gravity_m_s2[1])
                .hypot(gravity_m_s2[2]),
        )?;
        let mut ramps = Vec::with_capacity(layers.len());
        let mut previous_density = 0.0;
        let mut previous_depth = 0.0;
        for (index, layer) in layers.iter().enumerate() {
            checkpoint(cx)?;
            if !layer.top_depth_m.is_finite()
                || !layer.density_kg_m3.is_finite()
                || layer.density_kg_m3 <= 0.0
                || (index == 0 && layer.top_depth_m != 0.0)
                || (index > 0 && layer.top_depth_m <= previous_depth)
            {
                return Err(invalid(
                    "hydrostatic layers need positive densities and increasing depths starting at zero",
                ));
            }
            let offset_m2_s2 = finite(gravity * layer.top_depth_m)?;
            let density_jump_kg_m3 = finite(layer.density_kg_m3 - previous_density)?;
            if density_jump_kg_m3 != 0.0 {
                ramps.push(Ramp {
                    offset_m2_s2,
                    density_jump_kg_m3,
                });
            }
            previous_density = layer.density_kg_m3;
            previous_depth = layer.top_depth_m;
        }
        checkpoint(cx)?;
        Ok(Self {
            surface_point_m,
            gravity_m_s2,
            ramps,
        })
    }

    /// Apply this physical reservoir through the native pressure-load APIs.
    pub fn loading(&self) -> TetPressure<'_> {
        TetPressure::Stratified {
            profile: self,
            load_factor: 1.0,
        }
    }

    /// Gauge pressure at a world point [Pa], with no interpolation of density
    /// across a layer interface.
    pub fn pressure_pa(&self, point_m: [f64; 3], cx: &Cx<'_>) -> Result<f64, HyperTetError> {
        checkpoint(cx)?;
        let mut q = 0.0;
        for axis in 0..3 {
            let offset = finite(finite(point_m[axis])? - self.surface_point_m[axis])?;
            q = finite(self.gravity_m_s2[axis].mul_add(offset, q))?;
        }
        let mut pressure = Sum::default();
        for ramp in &self.ramps {
            checkpoint(cx)?;
            pressure.add(ramp.density_jump_kg_m3 * finite(q - ramp.offset_m2_s2)?.max(0.0))?;
        }
        checkpoint(cx)?;
        finite(pressure.value)
    }

    pub(in crate::hyper3) fn surface_point_m(&self) -> [f64; 3] {
        self.surface_point_m
    }

    pub(in crate::hyper3) fn gravity_m_s2(&self) -> [f64; 3] {
        self.gravity_m_s2
    }

    pub(super) fn is_zero(&self) -> bool {
        self.gravity_m_s2 == [0.0; 3]
    }

    pub(super) fn face_moments<const TANGENT: bool>(
        &self,
        surface: &ClosedTetPressureSurface,
        displacement_m: &[f64],
        face_index: usize,
        load_factor: f64,
        cx: &Cx<'_>,
    ) -> Result<FaceMoments, HyperTetError> {
        let mut q = [0.0; 3];
        for (local, &node) in surface.faces[face_index].iter().enumerate() {
            q[local] = surface.spatial_pressure_at(
                displacement_m,
                node,
                self.surface_point_m,
                0.0,
                self.gravity_m_s2,
            )?;
        }
        let mut weights = [Sum::default(); 3];
        let mut basis = [[Sum::default(); 3]; 3];
        for ramp in &self.ramps {
            checkpoint(cx)?;
            let values = [
                finite(q[0] - ramp.offset_m2_s2)?,
                finite(q[1] - ramp.offset_m2_s2)?,
                finite(q[2] - ramp.offset_m2_s2)?,
            ];
            if values.iter().all(|&value| value <= 0.0) {
                continue;
            }
            let moments = free_surface::face_moments::<TANGENT>(values)?;
            let coefficient = finite(load_factor * ramp.density_jump_kg_m3)?;
            for row in 0..3 {
                weights[row].add(coefficient * moments.weights_pa[row])?;
                if TANGENT {
                    for column in 0..3 {
                        basis[row][column].add(coefficient * moments.basis[row][column])?;
                    }
                }
            }
        }
        Ok(FaceMoments {
            weights_pa: weights.map(|sum| sum.value),
            // Here q has units m^2/s^2, so this basis includes dp/dq.
            // Multiplication by the shared gravity vector gives dp/dx.
            basis: basis.map(|row| row.map(|sum| sum.value)),
        })
    }

    /// Use one geometry admission and one pair of nodal depth fields for all
    /// layers. Sum element increments rather than large global potentials.
    pub(super) fn work_j(
        &self,
        surface: &ClosedTetPressureSurface,
        before_m: &[f64],
        after_m: &[f64],
        load_factor: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        let (before, _) = surface.positions_and_volume(before_m, cx)?;
        let (after, _) = surface.positions_and_volume(after_m, cx)?;
        let mut before_q = Vec::with_capacity(surface.nodes_m.len());
        let mut after_q = Vec::with_capacity(surface.nodes_m.len());
        for node in 0..surface.nodes_m.len() {
            checkpoint(cx)?;
            before_q.push(surface.spatial_pressure_at(
                before_m,
                node,
                self.surface_point_m,
                0.0,
                self.gravity_m_s2,
            )?);
            after_q.push(surface.spatial_pressure_at(
                after_m,
                node,
                self.surface_point_m,
                0.0,
                self.gravity_m_s2,
            )?);
        }
        let mut work = Sum::default();
        for (tet, &orientation) in surface.tetrahedra.iter().zip(&surface.orientations) {
            checkpoint(cx)?;
            let volume_before = orientation * signed_volume(tet.map(|node| before[node]))?;
            let volume_after = orientation * signed_volume(tet.map(|node| after[node]))?;
            let volume_mean = finite(0.5 * volume_before + 0.5 * volume_after)?;
            let volume_change = finite(volume_after - volume_before)?;
            let mut full_mean_change = 0.0;
            for &node in tet {
                for axis in 0..3 {
                    full_mean_change = finite((self.gravity_m_s2[axis] / 4.0).mul_add(
                        finite(after_m[3 * node + axis] - before_m[3 * node + axis])?,
                        full_mean_change,
                    ))?;
                }
            }
            for ramp in &self.ramps {
                checkpoint(cx)?;
                let first = tet.map(|node| before_q[node] - ramp.offset_m2_s2);
                let last = tet.map(|node| after_q[node] - ramp.offset_m2_s2);
                for value in first.into_iter().chain(last) {
                    finite(value)?;
                }
                let mean_before = free_surface::tet_mean(first)?;
                let mean_after = free_surface::tet_mean(last)?;
                if mean_before == 0.0 && mean_after == 0.0 {
                    continue;
                }
                let mean_change = if first.iter().chain(&last).all(|&value| value >= 0.0) {
                    // Cancel the layer datum before subtraction in deep water.
                    full_mean_change
                } else {
                    finite(mean_after - mean_before)?
                };
                let mean_pressure = finite(0.5 * mean_before + 0.5 * mean_after)?;
                let increment = finite(
                    -mean_pressure.mul_add(volume_change, finite(volume_mean * mean_change)?),
                )?;
                let coefficient = finite(load_factor * ramp.density_jump_kg_m3)?;
                work.add(coefficient * increment)?;
            }
        }
        checkpoint(cx)?;
        finite(work.value)
    }

    /// Exact linear combination of the native clipped-volume discrete
    /// gradients. Nodal arguments are q=g dot (x-surface), not pressure.
    pub(in crate::hyper3) fn algorithmic_mean(
        &self,
        before_q: [f64; 4],
        after_q: [f64; 4],
        load_factor: f64,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<PressureMean, HyperTetError> {
        let mut before = Sum::default();
        let mut after = Sum::default();
        let mut derivative = [Sum::default(); 4];
        let mut weights = [Sum::default(); 4];
        let mut matrix = [[Sum::default(); 4]; 4];
        for ramp in &self.ramps {
            checkpoint(cx)?;
            let first = before_q.map(|q| q - ramp.offset_m2_s2);
            let last = after_q.map(|q| q - ramp.offset_m2_s2);
            for value in first.into_iter().chain(last) {
                finite(value)?;
            }
            let mean = algorithmic::evaluate(first, last, true, tangent)?;
            let coefficient = finite(load_factor * ramp.density_jump_kg_m3)?;
            before.add(coefficient * mean.before_pa)?;
            after.add(coefficient * mean.after_pa)?;
            for i in 0..4 {
                derivative[i].add(coefficient * mean.after_derivative[i])?;
                weights[i].add(coefficient * mean.weights[i])?;
                if tangent {
                    for j in 0..4 {
                        matrix[i][j].add(coefficient * mean.tangent_after[i][j])?;
                    }
                }
            }
        }
        Ok(PressureMean {
            before_pa: before.value,
            after_pa: after.value,
            after_derivative: derivative.map(|sum| sum.value),
            weights: weights.map(|sum| sum.value),
            tangent_after: matrix.map(|row| row.map(|sum| sum.value)),
        })
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) -> Result<(), HyperTetError> {
        let corrected = finite(finite(value)? - self.correction)?;
        let next = finite(self.value + corrected)?;
        self.correction = finite((next - self.value) - corrected)?;
        self.value = next;
        Ok(())
    }
}
