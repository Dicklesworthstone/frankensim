//! Prescribed stratified density with an optional continuous depth gradient.
//!
//! Depth d is measured along the unit gravity direction. Hydrostatic pressure
//! is |g| sum_j (rho_j-rho_{j-1}) max(d-depth_j,0) + |g| beta max(d,0)^2/2.
//! Here beta is a nonnegative density gradient in kg/m^4, and rho_{-1}=0.
//! The same continuous pressure potential supplies surface traction, its
//! follower tangent, exact endpoint work, and the native discrete gradient.
//! Density may decrease with depth; this prescribes a reservoir, not its
//! stability or fluid motion. No local material density is substituted for
//! the displaced fluid's density.

use fs_exec::Cx;
use fs_math::dd::Dd;

use super::algorithmic::{self, PressureMean};
use super::free_surface::{self, FaceMoments};
use super::{
    ClosedTetPressureSurface, HyperTetError, TetPressure, checkpoint, finite, invalid,
    signed_volume, signed_volume_increment,
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
    depth_m: f64,
    pressure_gradient_pa_m: f64,
    density_kg_m3: f64,
}

/// Immutable hydrostatic loading in a prescribed stratified reservoir.
///
/// The last layer extends indefinitely. The base layers have positive constant
/// densities; an optional nonnegative gradient adds beta*depth to each layer.
/// Pressure is continuous across every interface and zero above the free
/// surface. Both increasing and decreasing base density profiles are admitted.
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
    depth_gradient: [f64; 3],
    density_gradient_kg_m4: f64,
    quadratic_pressure_pa_m2: f64,
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
            let density_jump = finite(layer.density_kg_m3 - previous_density)?;
            if density_jump != 0.0 {
                ramps.push(Ramp {
                    depth_m: layer.top_depth_m,
                    pressure_gradient_pa_m: finite(gravity * density_jump)?,
                    density_kg_m3: layer.density_kg_m3,
                });
            }
            previous_density = layer.density_kg_m3;
            previous_depth = layer.top_depth_m;
        }
        checkpoint(cx)?;
        Ok(Self {
            surface_point_m,
            gravity_m_s2,
            depth_gradient: if gravity == 0.0 {
                [0.0; 3]
            } else {
                gravity_m_s2.map(|value| value / gravity)
            },
            density_gradient_kg_m4: 0.0,
            quadratic_pressure_pa_m2: 0.0,
            ramps,
        })
    }

    /// Add or replace a continuous density increase beta*depth [kg/m^3]
    /// on top of the prepared base layers. With one layer this gives the
    /// unbounded law rho(depth)=rho_surface+beta*depth.
    ///
    /// A negative slope is refused because it eventually produces negative
    /// density in this unbounded reservoir. Zero slope retains the prepared
    /// piecewise-constant path exactly. Geometry, material state and the
    /// native solver are unchanged; the added quadratic pressure is integrated
    /// analytically, including its force derivative and endpoint work.
    pub fn with_density_gradient(mut self, gradient_kg_m4: f64) -> Result<Self, HyperTetError> {
        if !gradient_kg_m4.is_finite() || gradient_kg_m4 < 0.0 {
            return Err(invalid(
                "unbounded hydrostatic density gradient must be finite and nonnegative",
            ));
        }
        let gravity = self.gravity_m_s2[0]
            .hypot(self.gravity_m_s2[1])
            .hypot(self.gravity_m_s2[2]);
        self.density_gradient_kg_m4 = gradient_kg_m4;
        self.quadratic_pressure_pa_m2 = finite(gravity * gradient_kg_m4)? * 0.5;
        Ok(self)
    }

    /// Physical density at nonnegative downward depth, or zero above the
    /// free surface. A layer interface takes the density immediately below it.
    pub fn density_at_depth_kg_m3(&self, depth_m: f64) -> Result<f64, HyperTetError> {
        finite(depth_m)?;
        if depth_m < 0.0 {
            return Ok(0.0);
        }
        let index = self.ramps.partition_point(|ramp| ramp.depth_m <= depth_m) - 1;
        finite(
            self.density_gradient_kg_m4
                .mul_add(depth_m, self.ramps[index].density_kg_m3),
        )
    }

    // Admit the whole physical density column down to the deepest sampled
    // point. A density decrease must not conceal overflow in an earlier layer.
    fn admit_depth(&self, deepest_m: f64, cx: &Cx<'_>) -> Result<(), HyperTetError> {
        if self.density_gradient_kg_m4 == 0.0 || deepest_m <= 0.0 {
            return Ok(());
        }
        for (index, ramp) in self.ramps.iter().enumerate() {
            checkpoint(cx)?;
            if ramp.depth_m > deepest_m {
                break;
            }
            let end = self
                .ramps
                .get(index + 1)
                .map_or(deepest_m, |next| next.depth_m.min(deepest_m));
            finite(self.density_gradient_kg_m4.mul_add(end, ramp.density_kg_m3))?;
        }
        Ok(())
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
            q = finite(self.depth_gradient[axis].mul_add(offset, q))?;
        }
        self.admit_depth(q, cx)?;
        let mut pressure = Sum::default();
        if self.quadratic_pressure_pa_m2 != 0.0 {
            let depth = q.max(0.0);
            pressure.add((self.quadratic_pressure_pa_m2 * depth) * depth)?;
        }
        for ramp in &self.ramps {
            checkpoint(cx)?;
            pressure.add(ramp.pressure_gradient_pa_m * finite(q - ramp.depth_m)?.max(0.0))?;
        }
        checkpoint(cx)?;
        finite(pressure.value)
    }

    pub(in crate::hyper3) fn surface_point_m(&self) -> [f64; 3] {
        self.surface_point_m
    }

    pub(in crate::hyper3) fn depth_gradient(&self) -> [f64; 3] {
        self.depth_gradient
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
                self.depth_gradient,
            )?;
        }
        self.admit_depth(q.into_iter().fold(0.0, f64::max), cx)?;
        let mut weights = [Sum::default(); 3];
        let mut basis = [[Sum::default(); 3]; 3];
        if self.quadratic_pressure_pa_m2 != 0.0 {
            let moments = free_surface::face_square_moments::<TANGENT>(q)?;
            let coefficient = finite(load_factor * self.quadratic_pressure_pa_m2)?;
            for row in 0..3 {
                weights[row].add(coefficient * moments.weights_pa[row])?;
                if TANGENT {
                    for column in 0..3 {
                        basis[row][column].add(coefficient * moments.basis[row][column])?;
                    }
                }
            }
        }
        for ramp in &self.ramps {
            checkpoint(cx)?;
            let values = [
                finite(q[0] - ramp.depth_m)?,
                finite(q[1] - ramp.depth_m)?,
                finite(q[2] - ramp.depth_m)?,
            ];
            if values.iter().all(|&value| value <= 0.0) {
                continue;
            }
            let moments = free_surface::face_moments::<TANGENT>(values)?;
            let coefficient = finite(load_factor * ramp.pressure_gradient_pa_m)?;
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
            // The basis includes dp/d(depth). Multiplication by the shared
            // unit depth gradient supplies the complete spatial derivative.
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
                self.depth_gradient,
            )?);
            after_q.push(surface.spatial_pressure_at(
                after_m,
                node,
                self.surface_point_m,
                0.0,
                self.depth_gradient,
            )?);
        }
        self.admit_depth(
            before_q.iter().chain(&after_q).copied().fold(0.0, f64::max),
            cx,
        )?;
        let mut work = Sum::default();
        for (tet, &orientation) in surface.tetrahedra.iter().zip(&surface.orientations) {
            checkpoint(cx)?;
            let volume_before = orientation * signed_volume(tet.map(|node| before[node]))?;
            let volume_after = orientation * signed_volume(tet.map(|node| after[node]))?;
            let volume_mean = finite(0.5 * volume_before + 0.5 * volume_after)?;
            let volume_change = Dd::from_f64(orientation)
                * signed_volume_increment(
                    tet.map(|node| before[node]),
                    tet.map(|node| after[node]),
                )?;
            let mut depth_changes = [0.0; 4];
            for (local, &node) in tet.iter().enumerate() {
                for axis in 0..3 {
                    depth_changes[local] = finite(self.depth_gradient[axis].mul_add(
                        finite(after_m[3 * node + axis] - before_m[3 * node + axis])?,
                        depth_changes[local],
                    ))?;
                }
            }
            let full_mean_change = finite(depth_changes.into_iter().map(|v| v / 4.0).sum())?;
            if self.quadratic_pressure_pa_m2 != 0.0 {
                let first = tet.map(|node| before_q[node]);
                let last = tet.map(|node| after_q[node]);
                let mean_before = free_surface::tet_square_mean(first)?;
                let mean_after = free_surface::tet_square_mean(last)?;
                let mean_change = if first.iter().chain(&last).all(|&v| v >= 0.0) {
                    // Difference of quadratic products from displacement
                    // increments, retaining tiny drift at a large depth.
                    let midpoint: [f64; 4] =
                        core::array::from_fn(|i| 0.5 * first[i] + 0.5 * last[i]);
                    let mut change = Sum::default();
                    for i in 0..4 {
                        for j in i..4 {
                            change.add(
                                (midpoint[i] / 10.0) * depth_changes[j]
                                    + (midpoint[j] / 10.0) * depth_changes[i],
                            )?;
                        }
                    }
                    change.value
                } else {
                    finite(mean_after - mean_before)?
                };
                let mean = finite(0.5 * mean_before + 0.5 * mean_after)?;
                let increment = finite(
                    (-Dd::from_f64(mean) * volume_change
                        - Dd::from_f64(volume_mean) * Dd::from_f64(mean_change))
                    .to_f64(),
                )?;
                work.add(finite(load_factor * self.quadratic_pressure_pa_m2)? * increment)?;
            }
            for ramp in &self.ramps {
                checkpoint(cx)?;
                let first = tet.map(|node| before_q[node] - ramp.depth_m);
                let last = tet.map(|node| after_q[node] - ramp.depth_m);
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
                    (-Dd::from_f64(mean_pressure) * volume_change
                        - Dd::from_f64(volume_mean) * Dd::from_f64(mean_change))
                    .to_f64(),
                )?;
                let coefficient = finite(load_factor * ramp.pressure_gradient_pa_m)?;
                work.add(coefficient * increment)?;
            }
        }
        checkpoint(cx)?;
        finite(work.value)
    }

    /// Exact linear combination of the native clipped-volume discrete
    /// gradients. Nodal arguments are downward depth in metres, not pressure.
    /// The optional quadratic law uses the same native Secant/Dual algebra.
    pub(in crate::hyper3) fn algorithmic_mean(
        &self,
        before_q: [f64; 4],
        after_q: [f64; 4],
        load_factor: f64,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<PressureMean, HyperTetError> {
        self.admit_depth(before_q.into_iter().chain(after_q).fold(0.0, f64::max), cx)?;
        let mut before = Sum::default();
        let mut after = Sum::default();
        let mut derivative = [Sum::default(); 4];
        let mut weights = [Sum::default(); 4];
        let mut matrix = [[Sum::default(); 4]; 4];
        if self.quadratic_pressure_pa_m2 != 0.0 {
            checkpoint(cx)?;
            let mean = algorithmic::evaluate_square(before_q, after_q, tangent)?;
            let coefficient = finite(load_factor * self.quadratic_pressure_pa_m2)?;
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
        for ramp in &self.ramps {
            checkpoint(cx)?;
            let first = before_q.map(|q| q - ramp.depth_m);
            let last = after_q.map(|q| q - ramp.depth_m);
            for value in first.into_iter().chain(last) {
                finite(value)?;
            }
            let mean = algorithmic::evaluate(first, last, true, tangent)?;
            let coefficient = finite(load_factor * ramp.pressure_gradient_pa_m)?;
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
