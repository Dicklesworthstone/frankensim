//! Full-wave diffraction by finite pressure-release (sound-soft) screens.
//!
//! With G=exp(ikr)/(4*pi*r), write p_scattered=S*sigma and solve
//! integral(test_i S sigma)=-integral(test_i p_incident). Pressure is
//! continuous and zero on both faces; sigma=-[d_n p] is a single-layer
//! density [Pa/m], NOT a pressure jump or a prescribed surface velocity.
//! Piecewise constants are retained on every triangle, including free edges:
//! the physical inverse-square-root edge singularity must not be clamped to
//! zero. Caller-supplied edge grading is an independent convergence axis.
//!
//! The original analytic Laplace triangle moments, outgoing Green remainder,
//! positive angular Gram identity and cancellable LU are reused. This is a
//! Galerkin discretization of the single-layer screen equation, not a local
//! shadow, Kirchhoff-aperture or infinite-admittance approximation. See
//! Chandler-Wilde et al., Numer. Math. 147 (2021), doi:10.1007/s00211-021-01182-y.
use super::*;
use crate::{ComplexSum, helmholtz::bc::IncidentField};
use fs_la::eigen_complex::LuComplex;

mod receivers;
pub use receivers::PreparedSoftScreenReceivers;

/// Original open geometry and frequency-independent P0 single-layer moments.
/// Multiple disconnected open components retain every mutual interaction.
pub struct PreparedSoftScreen<'s> {
    surface: &'s OpenPanels,
    facets: Vec<Facet>,
    samples: Vec<Vec<Sample>>,
    static_single: Vec<f64>,
    options: OpenOptions,
    highest_k: f64,
    bytes: usize,
}

/// One reusable full-wave matrix and LU at the declared frequency and medium.
pub struct SoftScreenFrequency<'p, 's> {
    prepared: &'p PreparedSoftScreen<'s>,
    k: f64,
    medium: Medium,
    matrix: Vec<C64>,
    factor: LuComplex,
}

/// Complete coherent incident-bank limits, checked before any source solve.
#[derive(Clone, Copy, Debug)]
pub struct SoftScreenScatteringLimits {
    pub maximum_cases: usize,
    /// Every retained case/triangle density; no source or edge is truncated.
    pub maximum_density_entries: usize,
    /// Analytic source terms times (quadrature points + triangle-clearance
    /// visits). This ceiling also bounds source sampling at receivers.
    pub maximum_source_point_products: u64,
    pub minimum_source_clearance_m: f64,
}

/// The original incident field and solved scattered single-layer density.
/// The prepared frequency's factorization need not outlive this result.
pub struct SoftScreenScattering<'p, 's, 'i> {
    prepared: &'p PreparedSoftScreen<'s>,
    incident: IncidentField<'i>,
    density: Vec<C64>,
    k: f64,
    medium: Medium,
    maximum_source_point_products: u64,
    pub backward_residual: f64,
    /// Outgoing scattered power from the ORIGINAL angular Gram operator [W].
    /// This is not independent evidence of mesh or angular convergence. A
    /// stationary ideal pressure-release sheet absorbs no acoustic power.
    pub scattered_power_w: f64,
}

// Pair tiles bound the interval between cancellation polls independently of
// the dense row length. One worker runs directly; more use the original
// deterministic row executor on pair indices, then scatter in serial order.
fn pair_tiles<T: Send>(
    n: usize,
    workers: usize,
    cancelled: &mut impl FnMut() -> bool,
    value: impl Fn(usize, usize) -> T + Sync,
    mut record: impl FnMut(usize, usize, T),
) -> Result<(), HelmholtzError> {
    let mut pairs = (0..n).flat_map(|i| (i..n).map(move |j| (i, j)));
    loop {
        poll(cancelled)?;
        let tile: Vec<_> = pairs.by_ref().take(32 * workers).collect();
        if tile.is_empty() {
            break;
        }
        if workers == 1 {
            for (i, j) in tile {
                record(i, j, value(i, j));
            }
        } else {
            let values = rows_on_workers(0..tile.len(), workers, |index| {
                let (i, j) = tile[index];
                value(i, j)
            });
            for ((i, j), v) in tile.into_iter().zip(values) {
                record(i, j, v);
            }
        }
    }
    Ok(())
}

impl<'s> PreparedSoftScreen<'s> {
    /// Retain the supplied mesh and its declared triangle/angular/edge rungs.
    /// `maximum_unknowns` counts P0 TRIANGLES here, not the hard-screen P1
    /// vertices. All facets, including those with only free-edge vertices,
    /// have a density unknown. Logical cache bounds exclude allocator overhead
    /// and temporary pair tiles. At most 64 declared workers execute at once.
    pub fn new(
        surface: &'s OpenPanels,
        highest_k: f64,
        options: OpenOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HelmholtzError> {
        poll(&mut cancelled)?;
        let n = surface.triangles.len();
        if n == 0 || n > MAX_DENSE_PANELS || n > options.maximum_unknowns {
            return Err(HelmholtzError::WorkCap { panels: n });
        }
        if !highest_k.is_finite()
            || highest_k <= 0.0
            || ![4, 8].contains(&options.triangle_order)
            || options.angular_refinement > 2
            || options.workers == 0
            || options.edge_level > options.maximum_edge_level
            || options.maximum_work == 0
        {
            return Err(bad(
                "invalid soft-screen frequency, quadrature or work declaration",
            ));
        }
        let size = surface.resolution_length();
        if highest_k > TAU / (MIN_PANELS_PER_WAVELENGTH * size) {
            return Err(HelmholtzError::TooCoarse {
                panels_per_wavelength: TAU / highest_k / size,
            });
        }
        let samples_per_facet = options.triangle_order.pow(2);
        let bytes = n
            .checked_mul(n)
            .and_then(|m| m.checked_mul(8))
            .and_then(|m| {
                m.checked_add(
                    n * (core::mem::size_of::<Facet>()
                        + core::mem::size_of::<Vec<Sample>>()
                        + samples_per_facet * core::mem::size_of::<Sample>()),
                )
            })
            .and_then(|m| m.checked_add(core::mem::size_of::<Self>()))
            .ok_or_else(|| bad("soft-screen cache extent overflow"))?;
        if bytes > options.maximum_cache_bytes {
            return Err(HelmholtzError::CacheBudget {
                required_bytes: bytes,
                maximum_bytes: options.maximum_cache_bytes,
            });
        }
        if work(n, n, options) > u128::from(options.maximum_work) {
            return Err(bad("soft-screen frequency work cap exceeded"));
        }
        // A bounding-box diagonal bounds all source-pair separations without
        // an O(vertices²) scan before admission. Its conservatism is explicit:
        // this is a sufficient angular-support condition, not an error bound.
        let mut lower = surface.vertices[0];
        let mut upper = lower;
        for (index, point) in surface.vertices.iter().enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            for axis in 0..3 {
                lower[axis] = lower[axis].min(point[axis]);
                upper[axis] = upper[axis].max(point[axis]);
            }
        }
        let diameter_bound = norm(core::array::from_fn(|axis| {
            (upper[axis] - lower[axis]).next_up()
        }))
        .next_up();
        if highest_k * diameter_bound
            > ANGULAR_REACH * f64::from(1_u32 << options.angular_refinement)
        {
            return Err(bad("soft-screen band exceeds its angular-resolution chart"));
        }
        let mut facets = Vec::with_capacity(n);
        for &ids in &surface.triangles {
            poll(&mut cancelled)?;
            let triangle = ids.map(|i| surface.vertices[i]);
            let cross = cross(sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]));
            let size = norm(cross);
            facets.push(Facet::new(
                triangle,
                ids,
                0.5 * size,
                cross.map(|v| v / size),
            )?);
        }
        let samples = facets
            .iter()
            .map(|f| samples(f, options.triangle_order))
            .collect();
        let mut out = Self {
            surface,
            facets,
            samples,
            static_single: zeroed(n * n, 0.0)?,
            options,
            highest_k,
            bytes,
        };
        let cached_samples = &out.samples;
        let cached_facets = &out.facets;
        let pair = |i: usize, j: usize| {
            let integrate = |a: usize, b: usize| {
                let mut sum = fs_math::dd::Dd::ZERO;
                for x in &cached_samples[a] {
                    sum = sum
                        + fs_math::dd::Dd::from_f64(x.weight)
                            * fs_math::dd::Dd::from_f64(
                                static_moments(x.point, &cached_facets[b], a == b).0,
                            );
                }
                sum.to_f64()
            };
            let forward = integrate(i, j);
            if i == j {
                forward
            } else {
                0.5 * forward + 0.5 * integrate(j, i)
            }
        };
        let workers = options.workers.min(n).min(64);
        pair_tiles(n, workers, &mut cancelled, pair, |i, j, value| {
            out.static_single[i * n + j] = value;
            out.static_single[j * n + i] = value;
        })?;
        if out.static_single.iter().any(|v| !v.is_finite()) {
            return Err(bad("nonfinite soft-screen static moment"));
        }
        poll(&mut cancelled)?;
        Ok(out)
    }

    pub fn unknowns(&self) -> usize {
        self.facets.len()
    }
    pub fn cache_bytes(&self) -> usize {
        self.bytes
    }

    /// Retain and factor the original P0 single-layer matrix once. The byte
    /// ceiling covers matrix, LU, worst pair-tile scratch, angular row and LU
    /// vectors; it excludes the separately admitted geometry cache and runtime
    /// thread/allocator metadata. Source solves reuse this frequency owner.
    pub fn prepare_frequency(
        &self,
        k: f64,
        medium: Medium,
        maximum_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SoftScreenFrequency<'_, 's>, HelmholtzError> {
        poll(&mut cancelled)?;
        let rho_omega = k * medium.sound_speed * medium.density;
        if !k.is_finite()
            || k <= 0.0
            || k > self.highest_k
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !rho_omega.is_finite()
            || rho_omega <= 0.0
            || !rho_omega.recip().is_finite()
            || !(medium.density * medium.sound_speed).is_finite()
        {
            return Err(bad(
                "soft-screen frequency requires its admitted band and finite positive medium",
            ));
        }
        let n = self.unknowns();
        let workers = self.options.workers.min(n).min(64);
        let lu_scratch = if workers == 1 {
            0
        } else {
            fs_la::eigen_complex::lu_complex_blocked_workspace_bytes(n, workers)
                .ok_or_else(|| bad("soft-screen LU workspace extent overflow"))?
        };
        let required = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(32))
            .and_then(|v| v.checked_add(n * (64 + core::mem::size_of::<usize>())))
            .and_then(|v| v.checked_add(32 * workers * (32 + 2 * core::mem::size_of::<usize>())))
            .and_then(|v| v.checked_add(lu_scratch))
            .and_then(|v| v.checked_add(core::mem::size_of::<SoftScreenFrequency<'_, '_>>()))
            .ok_or_else(|| bad("soft-screen frequency storage overflow"))?;
        if required > maximum_bytes {
            return Err(HelmholtzError::CacheBudget {
                required_bytes: required,
                maximum_bytes,
            });
        }
        let mut matrix = zeroed(n * n, C64::ZERO)?;
        let pair = |i: usize, j: usize| {
            let mut sum = fs_math::dd::Dd::ZERO;
            for x in &self.samples[i] {
                for y in &self.samples[j] {
                    sum = sum
                        + fs_math::dd::Dd::from_f64(
                            green_remainder(k, norm(sub(x.point, y.point))).re,
                        ) * fs_math::dd::Dd::from_f64(x.weight)
                            * fs_math::dd::Dd::from_f64(y.weight);
                }
            }
            self.static_single[i * n + j] + sum.to_f64()
        };
        pair_tiles(n, workers, &mut cancelled, pair, |i, j, value| {
            matrix[i * n + j].re = value;
            matrix[j * n + i].re = value;
        })?;
        // Im G = k/(16*pi²) integral_S² exp(ik d.(x-y)) dOmega.
        // The original rule stores one partner of each antipodal pair.
        // This Gram form retains nonnegative outgoing power before solving.
        for (direction, weight) in directions(self.options.angular_refinement) {
            poll(&mut cancelled)?;
            let c = self.coefficients(k, direction, &mut cancelled)?;
            let scale = k * weight / (8.0 * PI * PI);
            for i in 0..n {
                if i % 32 == 0 {
                    poll(&mut cancelled)?;
                }
                for j in i..n {
                    matrix[i * n + j].im += scale * (c[i].conj() * c[j]).re;
                    matrix[j * n + i].im = matrix[i * n + j].im;
                }
            }
        }
        if matrix.iter().any(|&v| !finite(v)) {
            return Err(bad("nonfinite soft-screen operator"));
        }
        let factor = fs_la::eigen_complex::with_lu_workers(workers, || {
            lu_complex_with_cancel(&matrix, n, &mut cancelled)
        })
        .map_err(lu_error)?;
        poll(&mut cancelled)?;
        Ok(SoftScreenFrequency {
            prepared: self,
            k,
            medium,
            matrix,
            factor,
        })
    }

    fn coefficients(
        &self,
        k: f64,
        direction: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        // A common phase cancels from the radiation Gram form. Removing the
        // original mesh origin avoids large translation-dependent phases.
        let origin = self.surface.vertices[0];
        let mut out = Vec::with_capacity(self.unknowns());
        for (index, samples) in self.samples.iter().enumerate() {
            if index % 8 == 0 {
                poll(cancelled)?;
            }
            let mut sum = ComplexSum::default();
            for sample in samples {
                let phase = -k * dot(direction, sub(sample.point, origin));
                sum.add_product(
                    C64::new(det::cos(phase), det::sin(phase)),
                    C64::from_re(sample.weight),
                );
            }
            out.push(sum.finish());
        }
        Ok(out)
    }

    fn admit_points(
        &self,
        points: &[Point],
        separation: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HelmholtzError> {
        let mut scale = self
            .surface
            .vertices
            .iter()
            .flatten()
            .fold(1.0_f64, |s, x| s.max(x.abs()));
        for point in points {
            if point.iter().any(|x| !x.is_finite() || x.abs() > 1e9) {
                return Err(bad(
                    "soft-screen source/receiver coordinates are unresolved",
                ));
            }
            for x in point {
                scale = scale.max(x.abs());
            }
        }
        if !separation.is_finite() || separation <= 64.0 * f64::EPSILON * scale {
            return Err(bad("soft-screen clearance is below coordinate resolution"));
        }
        for &point in points {
            for (i, facet) in self.facets.iter().enumerate() {
                if i % 32 == 0 {
                    poll(cancelled)?;
                }
                let distance = crate::near_field::distance(point, facet.triangle);
                if !distance.is_finite() || distance <= separation {
                    return Err(bad("soft-screen source or receiver touches a face or edge"));
                }
            }
        }
        Ok(())
    }
}

impl<'p, 's> SoftScreenFrequency<'p, 's> {
    pub fn scatter<'i>(
        &self,
        incident: IncidentField<'i>,
        limits: SoftScreenScatteringLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SoftScreenScattering<'p, 's, 'i>, HelmholtzError> {
        self.scatter_many(std::slice::from_ref(&incident), limits, cancelled)?
            .pop()
            .ok_or_else(|| bad("empty soft-screen scattering result"))
    }

    /// Complete coherent cases share the retained factor. The incident field
    /// is sampled at every original Duffy point, never only at a centroid.
    /// Refusal or cancellation publishes no partial source bank.
    pub fn scatter_many<'i>(
        &self,
        incidents: &[IncidentField<'i>],
        limits: SoftScreenScatteringLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<SoftScreenScattering<'p, 's, 'i>>, HelmholtzError> {
        poll(&mut cancelled)?;
        let n = self.prepared.unknowns();
        if incidents.is_empty()
            || incidents.len() > limits.maximum_cases
            || limits.maximum_cases == 0
            || limits.maximum_cases > MAX_RADIATION_FIELDS_PER_BATCH
            || limits.maximum_source_point_products == 0
        {
            return Err(bad(
                "soft-screen incident bank requires finite nonempty budgets",
            ));
        }
        let mut terms = 0_u128;
        for &incident in incidents {
            poll(&mut cancelled)?;
            incident.validate()?;
            terms += incident.source_count() as u128;
        }
        let points = n as u128 * (self.prepared.options.triangle_order.pow(2) + 1) as u128;
        let source_work = terms * points;
        if incidents.len() as u128 * n as u128 > limits.maximum_density_entries as u128
            || source_work > limits.maximum_source_point_products as u128
            || work(n, n, self.prepared.options)
                + 256 * source_work
                + 128 * incidents.len() as u128 * n as u128 * n as u128
                > self.prepared.options.maximum_work as u128
        {
            return Err(bad(
                "complete soft-screen source bank exceeds density or work budget",
            ));
        }
        for &incident in incidents {
            self.prepared.admit_points(
                &incident.source_positions(),
                limits.minimum_source_clearance_m,
                &mut cancelled,
            )?;
        }
        let mut out = Vec::with_capacity(incidents.len());
        let rho_omega = self.k * self.medium.sound_speed * self.medium.density;
        for &incident in incidents {
            let mut rhs = Vec::with_capacity(n);
            for samples in &self.prepared.samples {
                poll(&mut cancelled)?;
                let mut sum = ComplexSum::default();
                for sample in samples {
                    let field = incident.sample_wave(
                        self.k,
                        rho_omega,
                        Some(self.medium.density * self.medium.sound_speed),
                        sample.point,
                        &mut cancelled,
                    )?;
                    sum.add_product(field.pressure_pa, C64::from_re(-sample.weight));
                }
                rhs.push(sum.finish());
            }
            if rhs.iter().any(|&v| !finite(v)) {
                return Err(bad("soft-screen incident weak pressure overflows"));
            }
            let mut density = rhs.clone();
            self.factor
                .solve_with_cancel(&mut density, &mut cancelled)
                .map_err(lu_error)?;
            let mut residual = 0.0_f64;
            let mut scale = 0.0_f64;
            let mut power = ComplexSum::default();
            let mut absolute_power = 0.0_f64;
            for i in 0..n {
                poll(&mut cancelled)?;
                let mut value = ComplexSum::default();
                let mut radiating = ComplexSum::default();
                let mut row = rhs[i].abs();
                for j in 0..n {
                    let a = self.matrix[i * n + j];
                    value.add_product(a, density[j]);
                    radiating.add_product(density[j], C64::from_re(a.im));
                    row += a.abs() * density[j].abs();
                }
                let flux = density[i].conj() * radiating.finish();
                power.add(flux);
                absolute_power += flux.abs();
                residual = residual.max((value.finish() - rhs[i]).abs());
                scale = scale.max(row);
            }
            let backward_residual = residual / scale.max(f64::MIN_POSITIVE);
            let power = power.finish().re;
            let scattered_power_w = power / (2.0 * rho_omega);
            if !backward_residual.is_finite()
                || backward_residual > 1e-10
                || density.iter().any(|&v| !finite(v))
                || !scattered_power_w.is_finite()
                || !absolute_power.is_finite()
                || power < -1e-10 * absolute_power.max(f64::MIN_POSITIVE)
            {
                return Err(bad("soft-screen solve residual or outgoing-power refusal"));
            }
            out.push(SoftScreenScattering {
                prepared: self.prepared,
                incident,
                density,
                k: self.k,
                medium: self.medium,
                maximum_source_point_products: limits.maximum_source_point_products,
                backward_residual,
                scattered_power_w,
            });
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
}

impl SoftScreenScattering<'_, '_, '_> {
    /// Per-triangle single-layer density [Pa/m], sigma=-[d_n p]. Pressure
    /// itself has no jump. No free-edge density is clamped or discarded.
    pub fn single_layer_density_pa_m(&self) -> &[C64] {
        &self.density
    }

    /// Scattered far-field amplitude r*exp(-ikr)*p_scattered(r*d) [Pa m].
    /// This remains an outgoing field even when the total incident plane wave
    /// has no finite outgoing far field. Directions must already be unit.
    pub fn scattered_far_field(
        &self,
        direction: Point,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<C64, HelmholtzError> {
        poll(&mut cancelled)?;
        if direction.iter().any(|x| !x.is_finite()) || (norm(direction) - 1.0).abs() > 1e-10 {
            return Err(bad(
                "soft-screen far field requires a finite unit direction",
            ));
        }
        let coefficients = self
            .prepared
            .coefficients(self.k, direction, &mut cancelled)?;
        let mut sum = ComplexSum::default();
        for (i, (&c, &density)) in coefficients.iter().zip(&self.density).enumerate() {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            sum.add_product(c, density);
        }
        let dot = direction
            .into_iter()
            .zip(self.prepared.surface.vertices[0])
            .fold(fs_math::dd::Dd::ZERO, |sum, (a, b)| {
                sum + fs_math::dd::Dd::from_f64(a) * fs_math::dd::Dd::from_f64(b)
            });
        let (sine, cosine) = det::sin_cos_dd(-fs_math::dd::Dd::from_f64(self.k) * dot);
        let result =
            sum.finish() * C64::new(cosine.to_f64(), sine.to_f64()).scale(1.0 / (4.0 * PI));
        if !finite(result) {
            return Err(bad("soft-screen far field overflows"));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
