//! Physical point-force diffraction by the original spherical half-plane Green function.
//!
//! A force ON the fluid satisfies (Laplacian+k²)p=div(F delta), hence
//! p=F.grad_source G and u=grad_receiver(p)/(i omega rho), away from the source.
//! Both derivatives are analytic derivatives of the SAME Biot--Tolstoy integral
//! documented by the parent (Ouis 2002, doi:10.1006/jsvi.2000.3553). The mixed
//! receiver/source derivative supplies the complete near-field velocity tensor.
//! No displaced source pair or pressure/impedance approximation is used.
//!
//! Differentiation precedes the dimensionless contour substitution. At an exact
//! geometrical boundary q=0, the second derivative has the half-weight of the
//! visible Green term: its jump cancels the one-sided limit of the contour's
//! atan(q/s) derivative. This retains the smooth physical mixed derivative.

use super::*;

pub use crate::curle3d::refraction::full_wave::FullWaveForceLimits as HalfPlaneForceLimits;

const SCOPE: &str = "Prescribed stationary physical point force ON fluid [N] \
diffracted by one infinite straight-edged rigid or pressure-release sheet in a \
stationary homogeneous lossless fluid. Analytic source/receiver Green derivatives \
retain full spherical pressure and Euler velocity. Estimated quadrature errors; \
no finite barrier, flow, absorption, source feedback or measured-SPL qualification.";

/// Embedded quadrature estimates for the original Cartesian force transfer.
/// For an actual force, sum each component error times the corresponding |F_j|.
/// Coordinate, geometry and floating-point errors are outside these estimates.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneForceDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub minimum_edge_path_m: f64,
    pub pressure_quadrature_error_pa_per_n: [f64; 3],
    /// Rows are world velocity components; columns are world force components.
    pub velocity_quadrature_error_m_s_per_n: [[f64; 3]; 3],
}

/// One original physical force location, sheet, receiver, material and clock.
/// All three force directions share every quadrature node. Subsequent force
/// phasors require twelve complex products and DD sums, without allocation,
/// quadrature, transcendentals or a new solve.
#[derive(Clone, Copy, Debug)]
pub struct PreparedHalfPlaneForce {
    acoustic: FluidAcousticState,
    frequency_hz: f64,
    frame: HalfPlaneFrame,
    boundary: IdealPlaneBoundary,
    source_position_m: [f64; 3],
    receiver: HalfPlaneReceiver,
    pressure: [C64; 3],
    velocity: [[C64; 3]; 3],
    diagnostics: HalfPlaneForceDiagnostics,
}

impl PreparedHalfPlaneForce {
    /// Physical force ON the fluid in world axes, with exp(-i omega t) peak
    /// phasors on the original source clock. The source lies in open fluid;
    /// explicit receiver sheet faces retain their original limiting traces.
    ///
    /// Limits apply to the pressure-transfer Euclidean norm and the complete
    /// velocity-transfer Frobenius norm. They are not relative-error claims
    /// at a directional or coherent force null. All twelve physical channels
    /// share the existing bounded adaptive quadrature and its cancellation.
    ///
    /// # Errors
    /// Invalid geometry/frequency, source or edge singularities, unresolved
    /// physical transfers, exhausted work or cancellation publish no observer.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        frame: HalfPlaneFrame,
        boundary: IdealPlaneBoundary,
        source_position_m: [f64; 3],
        receiver: HalfPlaneReceiver,
        limits: HalfPlaneForceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let controls = QuadratureLimits::from(limits);
        controls.admit()?;
        let omega = TAU * frequency_hz;
        let k = omega / acoustic.sound_speed_m_s();
        let rho_omega = omega * acoustic.density_kg_m3();
        if !frequency_hz.is_finite()
            || frequency_hz <= 0.0
            || !k.is_finite()
            || k <= 0.0
            || !rho_omega.is_finite()
            || rho_omega <= 0.0
        {
            return Err(invalid(
                "half-plane physical force needs a resolved positive frequency",
            ));
        }
        let source = frame.local_receiver(HalfPlaneReceiver::WorldPosition(source_position_m))?;
        let point = frame.local_receiver(receiver)?;
        let kernel = ForceKernel::new(
            source,
            point,
            k,
            rho_omega,
            boundary,
            frame.axes(),
            controls,
        )?;
        let result = integrate_spectrum::<12, 3>(&kernel, [0.0; 12], controls, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            frequency_hz,
            frame,
            boundary,
            source_position_m,
            receiver,
            pressure: core::array::from_fn(|j| result.value[j]),
            velocity: core::array::from_fn(|i| {
                core::array::from_fn(|j| result.value[3 + 3 * i + j])
            }),
            diagnostics: HalfPlaneForceDiagnostics {
                evaluations: result.evaluations,
                maximum_depth_reached: result.maximum_depth_reached,
                retained_panels: result.retained_panels,
                minimum_edge_path_m: kernel.point.distance,
                pressure_quadrature_error_pa_per_n: core::array::from_fn(|j| result.error[j]),
                velocity_quadrature_error_m_s_per_n: core::array::from_fn(|i| {
                    core::array::from_fn(|j| result.error[3 + 3 * i + j])
                }),
            },
        })
    }

    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn frame(&self) -> HalfPlaneFrame {
        self.frame
    }
    pub const fn boundary(&self) -> IdealPlaneBoundary {
        self.boundary
    }
    pub const fn source_position_m(&self) -> [f64; 3] {
        self.source_position_m
    }
    pub const fn receiver(&self) -> HalfPlaneReceiver {
        self.receiver
    }
    pub const fn diagnostics(&self) -> &HalfPlaneForceDiagnostics {
        &self.diagnostics
    }
    pub const fn pressure_transfer_pa_per_n(&self) -> &[C64; 3] {
        &self.pressure
    }
    pub const fn velocity_transfer_m_s_per_n(&self) -> &[[C64; 3]; 3] {
        &self.velocity
    }

    pub fn evaluate_harmonic_line(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        if force_on_fluid_n.into_iter().any(|v| !finite(v)) {
            return Err(invalid(
                "half-plane physical force must have finite Cartesian components",
            ));
        }
        let pressure_pa = contract(self.pressure, force_on_fluid_n);
        let particle_velocity_m_s = self.velocity.map(|row| contract(row, force_on_fluid_n));
        if !finite(pressure_pa) || particle_velocity_m_s.into_iter().any(|v| !finite(v)) {
            return Err(invalid(
                "half-plane physical force field is unrepresentable",
            ));
        }
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa,
            particle_velocity_m_s,
        })
    }

    pub fn evaluate(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        let line = self.evaluate_harmonic_line(force_on_fluid_n)?;
        let flux = line
            .particle_velocity_m_s
            .map(|v| line.pressure_pa.scale(0.5) * v.conj());
        if flux.into_iter().any(|v| !finite(v)) {
            return Err(invalid("half-plane force intensity is unrepresentable"));
        }
        Ok(CoherentCompactField {
            pressure_pa: line.pressure_pa,
            particle_velocity_m_s: line.particle_velocity_m_s,
            active_intensity_w_m2: flux.map(|v| v.re),
            reactive_intensity_w_m2: flux.map(|v| v.im),
            scope: SCOPE,
        })
    }

    pub fn evaluate_energy_transport(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<crate::curle3d::energy::AcousticEnergyTransport, AeroacError> {
        super::super::energy_transport(
            self.acoustic,
            [0.0; 3],
            self.evaluate_harmonic_line(force_on_fluid_n)?,
        )
    }
}

fn contract(row: [C64; 3], force: [C64; 3]) -> C64 {
    let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
    for j in 0..3 {
        let (a, b, c, d) = (
            Dd::from_f64(row[j].re),
            Dd::from_f64(row[j].im),
            Dd::from_f64(force[j].re),
            Dd::from_f64(force[j].im),
        );
        re = re + a * c - b * d;
        im = im + a * d + b * c;
    }
    C64::new(re.to_f64(), im.to_f64())
}

fn rotate_tensor(axes: [[f64; 3]; 3], tensor: [[C64; 3]; 3]) -> [[C64; 3]; 3] {
    let columns: [[C64; 3]; 3] =
        core::array::from_fn(|j| rotate(axes, [tensor[0][j], tensor[1][j], tensor[2][j]]));
    core::array::from_fn(|i| rotate(axes, [columns[0][i], columns[1][i], columns[2][i]]))
}

struct ForceKernel {
    point: PointKernel,
    source_edge: [f64; 3],
    source_q: [[f64; 3]; 2],
    mixed_edge: [[f64; 3]; 3],
    mixed_q: [[[f64; 3]; 3]; 2],
    pressure_scale: f64,
    velocity_scale: f64,
    base: [C64; 12],
    breaks: [f64; 4],
    count: usize,
}

impl ForceKernel {
    #[allow(clippy::too_many_arguments)]
    fn new(
        source: [f64; 3],
        receiver: [f64; 3],
        k: f64,
        rho_omega: f64,
        boundary: IdealPlaneBoundary,
        axes: [[f64; 3]; 3],
        controls: QuadratureLimits,
    ) -> Result<Self, AeroacError> {
        // Reuse the original geometry, contour coordinates, phase, and exact
        // first-derivative geometrical terms. Reciprocity gives the source
        // derivative by exchanging coordinates; no displaced samples occur.
        let point = PointKernel::new(source, receiver, k, 1.0, boundary, axes, controls)?;
        let reciprocal = PointKernel::new(receiver, source, k, 1.0, boundary, axes, controls)?;
        let l = point.distance;
        let pressure_scale = point.velocity_scale;
        let velocity_scale = pressure_scale / l / rho_omega;
        if !velocity_scale.is_finite() || velocity_scale <= 0.0 {
            return Err(invalid(
                "half-plane force velocity scale is unrepresentable",
            ));
        }
        let r = det::hypot(receiver[0], receiver[1]);
        let rs = det::hypot(source[0], source[1]);
        let radial = [receiver[0] / r, receiver[1] / r, 0.0];
        let tangent = [-radial[1], radial[0], 0.0];
        let source_radial = [source[0] / rs, source[1] / rs, 0.0];
        let source_tangent = [-source_radial[1], source_radial[0], 0.0];
        let [ct, st] = half_angle(radial[0], radial[1]);
        let [cs, ss] = half_angle(source_radial[0], source_radial[1]);
        let cosine = [ct * cs + st * ss, ct * cs - st * ss];
        let sine = [st * cs - ct * ss, st * cs + ct * ss];
        let mixed_edge = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                radial[i] * source_radial[j]
                    - (if i == 2 && j == 2 { 1.0 } else { 0.0 })
                    - point.edge_gradient[i] * reciprocal.edge_gradient[j]
            })
        });
        let coefficient = (0.5 * l / det::sqrt(r)) / det::sqrt(rs);
        let mixed_q = core::array::from_fn(|term| {
            let angular_sign = if term == 0 { -1.0 } else { 1.0 };
            core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    coefficient
                        * ((cosine[term] * radial[i] - sine[term] * tangent[i]) * source_radial[j]
                            - angular_sign
                                * (sine[term] * radial[i] + cosine[term] * tangent[i])
                                * source_tangent[j])
                })
            })
        });
        let edge = point.phase.scale(1.0 / l);
        let edge_first = edge * C64::new(-1.0 / l, k);
        let edge_second = edge * C64::new(2.0 / l / l - k * k, -2.0 * k / l);
        let edge_mixed: [[C64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                edge_second.scale(point.edge_gradient[i] * reciprocal.edge_gradient[j])
                    + edge_first.scale(mixed_edge[i][j] / l)
            })
        });
        let mut base_mixed = [[C64::ZERO; 3]; 3];
        for term in 0..2 {
            let sign = if term == 0 { 1.0 } else { point.sign };
            let reflection = [1.0, if term == 0 { 1.0 } else { -1.0 }, 1.0];
            let delta: [f64; 3] = core::array::from_fn(|i| {
                (Dd::from_f64(receiver[i]) - Dd::from_f64(source[i] * reflection[i])).to_f64()
            });
            let d = det::hypot(det::hypot(delta[0], delta[1]), delta[2]);
            let visible = if point.eta[term] > 0.0 {
                1.0
            } else if point.eta[term] < 0.0 {
                0.0
            } else {
                0.5
            };
            let mut addition = edge_mixed.map(|row| row.map(|v| v.scale(1.0 / (8.0 * PI))));
            if visible > 0.0 {
                if d <= 0.0 || !d.is_finite() {
                    return Err(invalid("half-plane force coincides with a visible source"));
                }
                let direct = cis(Dd::from_f64(k) * Dd::from_f64(d)).scale(1.0 / d);
                let first = direct * C64::new(-1.0 / d, k);
                let second = direct * C64::new(2.0 / d / d - k * k, -2.0 * k / d);
                for i in 0..3 {
                    for j in 0..3 {
                        let normal_product = (delta[i] / d) * (delta[j] / d);
                        let mixed = (second.scale(normal_product)
                            + first.scale(((if i == j { 1.0 } else { 0.0 }) - normal_product) / d))
                        .scale(-reflection[j]);
                        addition[i][j] = super::super::add(
                            addition[i][j],
                            (mixed - edge_mixed[i][j]).scale(visible / (4.0 * PI)),
                        );
                    }
                }
            }
            for i in 0..3 {
                for j in 0..3 {
                    base_mixed[i][j] =
                        super::super::add(base_mixed[i][j], addition[i][j].scale(sign));
                }
            }
        }
        let world = rotate_tensor(axes, base_mixed);
        let mut base = [C64::ZERO; 12];
        for j in 0..3 {
            base[j] = -reciprocal.base[j + 1];
        }
        for i in 0..3 {
            for j in 0..3 {
                let value = world[i][j];
                base[3 + 3 * i + j] = C64::new(value.im, -value.re).scale(1.0 / rho_omega);
            }
        }
        if base.into_iter().any(|v| !finite(v))
            || mixed_edge.iter().flatten().any(|v| !v.is_finite())
            || mixed_q.iter().flatten().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("half-plane force derivatives are unrepresentable"));
        }
        // The second angular derivative has finite integrated weight even
        // as its transition narrows. Logarithmic finite regions resolve BOTH
        // sides of each width; a lone linear [q,1] panel could miss its peak.
        // Include the decay width so the final rational tail never starts
        // below every retained geometric/oscillatory scale.
        let mut candidates = [0.0, point.eta[0].abs(), point.eta[1].abs(), point.scale];
        candidates.sort_by(f64::total_cmp);
        let mut breaks = [0.0; 4];
        let mut count = 1;
        for value in candidates {
            if value > breaks[count - 1] {
                breaks[count] = value;
                count += 1;
            }
        }
        Ok(Self {
            source_edge: reciprocal.edge_gradient,
            source_q: reciprocal.q_gradient,
            mixed_edge,
            mixed_q,
            pressure_scale,
            velocity_scale,
            base,
            breaks,
            count,
            point,
        })
    }

    fn coordinate(&self, region: usize, x: f64) -> (f64, f64) {
        let lower = self.breaks[region];
        if region + 1 < self.count {
            let upper = self.breaks[region + 1];
            if lower == 0.0 || upper <= 2.0 * lower {
                let width = upper - lower;
                (lower + width * x, width)
            } else {
                let start = det::ln(lower);
                let width = det::ln(upper) - start;
                let value = det::exp(start + width * x);
                (value, width * value)
            }
        } else {
            let rest = 1.0 - x;
            (
                lower + self.point.scale * x / rest,
                self.point.scale / rest / rest,
            )
        }
    }

    // Algebraic contour envelopes decide phase refinement only; no tail is
    // removed. The Hessian needs its own bound, beyond the monopole gradient.
    fn tail_is_small(&self, w: f64) -> bool {
        if !w.is_finite() || w < 1.0 {
            return false;
        }
        let inv = 1.0 / w;
        let k = self.point.kappa.to_f64();
        let g = k * inv * inv + inv.powi(3);
        let dg = 3.0 * inv.powi(5) + 3.0 * k * inv.powi(4) + k * k * inv.powi(3);
        let ddg = 15.0 * inv.powi(7)
            + 15.0 * k * inv.powi(6)
            + 6.0 * k * k * inv.powi(5)
            + k * k * k * inv.powi(4);
        let mut p = 0.0;
        let mut u = 0.0;
        for term in 0..2 {
            let a = self.point.eta[term].abs();
            let norm = |v: [f64; 3]| v.into_iter().fold(0.0, det::hypot);
            let qx = norm(self.point.q_gradient[term]);
            let qs = norm(self.source_q[term]);
            let lm = self
                .mixed_edge
                .iter()
                .flatten()
                .fold(0.0_f64, |a, b| a.max(b.abs()));
            let qm = self.mixed_q[term]
                .iter()
                .flatten()
                .fold(0.0_f64, |a, b| a.max(b.abs()));
            p += 2.0 * a * dg + g * qs;
            u += 2.0 * a * dg * lm
                + g * qm
                + 2.0 * a * (dg + ddg)
                + dg * (qx + qs)
                + 2.0 * a * g * inv * inv * qx * qs;
        }
        let attenuation = self.point.phase_extent(w).1;
        // Every algebraic term decays at least w^-2. Four covers the three
        // pressure components and nine tensor components in their norms.
        let bound = |scale: f64, value: f64| {
            log_weighted(C64::ONE, det::ln(scale) - attenuation, 4.0 * w * value).re
        };
        bound(self.pressure_scale, p) <= 0.01 * self.point.absolute[0]
            && bound(self.velocity_scale, u) <= 0.01 * self.point.absolute[1]
    }
}

impl Spectrum<12> for ForceKernel {
    fn regions(&self) -> usize {
        self.count
    }
    fn extent(&self) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let (a, _) = self.coordinate(region, lower);
        let (b, _) = self.coordinate(region, upper);
        self.tail_is_small(a)
            || self.point.phase_extent(b).0 - self.point.phase_extent(a).0 <= FRAC_PI_2
    }

    fn sample(&self, region: usize, x: f64) -> Result<[C64; 12], AeroacError> {
        let (w, jacobian) = self.coordinate(region, x);
        let r = C64::new(1.0, w * w).sqrt();
        let delta = C64::new(0.0, w * w) / (r + C64::ONE);
        let attenuation = (self.point.kappa * Dd::from_f64(r.im)).to_f64();
        if !w.is_finite() || !jacobian.is_finite() || !attenuation.is_finite() {
            return Err(invalid(
                "half-plane force contour coordinate is unrepresentable",
            ));
        }
        let logarithmic = attenuation > 700.0;
        let decay = if logarithmic {
            1.0
        } else {
            det::exp(-attenuation)
        };
        let e = cis(self.point.kappa * Dd::from_f64(delta.re)).scale(decay);
        let kr = r.scale(self.point.kappa.to_f64());
        let z = C64::new(-kr.im, kr.re);
        let inverse = r.recip();
        let inverse3 = inverse * inverse * inverse;
        let g = e * (z - C64::ONE) * inverse3;
        let dg = e * (C64::from_re(3.0) - z.scale(3.0) + z * z) * inverse3 * inverse * inverse;
        let ddg = e
            * (z * z * z - (z * z).scale(6.0) + z.scale(15.0) - C64::from_re(15.0))
            * inverse3
            * inverse
            * inverse
            * inverse
            * inverse;
        let e45 = C64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2);
        let iw = C64::new(0.0, w);
        let mut gradient = [C64::ZERO; 3];
        let mut mixed = [[C64::ZERO; 3]; 3];
        for term in 0..2 {
            let sign = if term == 0 { 1.0 } else { self.point.sign };
            let a = self.point.eta[term];
            let angle = contour_arctangent(a, w);
            let angular_scale = a.abs().max(w);
            let aw = a / angular_scale;
            let ww = w / angular_scale;
            let denominator = C64::new(aw * aw, ww * ww);
            let ratio = C64::new(0.0, ww * ww) / denominator;
            let cl = (iw * dg * angle).scale(jacobian);
            let cq = (e45 * g * ratio).scale(jacobian);
            let cll = (iw * (dg + ddg) * angle).scale(jacobian);
            let clq = (e45 * dg * ratio).scale(jacobian);
            // Combine jacobian/width before forming the second angular
            // derivative: arbitrarily narrow shadow panels need no 1/q
            // intermediate that could overflow before its finite integral.
            let cqq = (e45 * g * ratio / denominator).scale(-2.0 * aw * (jacobian / angular_scale));
            for j in 0..3 {
                gradient[j] = super::super::add(
                    gradient[j],
                    (cl.scale(self.source_edge[j]) + cq.scale(self.source_q[term][j])).scale(sign),
                );
                for i in 0..3 {
                    let value = cl.scale(self.mixed_edge[i][j])
                        + cq.scale(self.mixed_q[term][i][j])
                        + cll.scale(self.point.edge_gradient[i] * self.source_edge[j])
                        + clq.scale(
                            self.point.edge_gradient[i] * self.source_q[term][j]
                                + self.point.q_gradient[term][i] * self.source_edge[j],
                        )
                        + cqq.scale(self.point.q_gradient[term][i] * self.source_q[term][j]);
                    mixed[i][j] = super::super::add(mixed[i][j], value.scale(sign));
                }
            }
        }
        let scale = |value: C64, physical: f64| {
            if logarithmic {
                log_weighted(value, det::ln(physical) - attenuation, 1.0)
            } else {
                value.scale(physical)
            }
        };
        let gradient = rotate(self.point.axes, gradient);
        let tensor = rotate_tensor(self.point.axes, mixed);
        let mut result = self.base.map(|v| v.scale(1.0 / self.count as f64));
        for j in 0..3 {
            result[j] = super::super::add(
                result[j],
                scale(-self.point.phase * gradient[j], self.pressure_scale),
            );
            for i in 0..3 {
                let v = self.point.phase * tensor[i][j];
                result[3 + 3 * i + j] = super::super::add(
                    result[3 + 3 * i + j],
                    scale(C64::new(-v.im, v.re), self.velocity_scale),
                );
            }
        }
        if result.into_iter().any(|v| !finite(v)) {
            return Err(invalid(
                "half-plane force outgoing contour is unrepresentable",
            ));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
