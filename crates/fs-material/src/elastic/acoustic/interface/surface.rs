//! Coherent reflection by an infinite isotropic solid surface and passive support.
//!
//! The original solid occupies z < 0; +z points out of it. The wave, traction,
//! phase and modal velocity conventions are those of the bonded interface.
//! A free surface sets all three tractions to zero; a stationary clamp sets
//! all three velocities to zero. P and SV generally convert into each other,
//! including an evanescent P field for supercritical SV incidence. SH remains
//! independent for axisymmetric boundaries. A passive local support imposes
//! `traction = -Z velocity`; a reciprocal tensor can couple all three modes.
//! Its resistance absorbs work; its reactance preserves coherent reflection phase.
//!
//! The same original tensor supplies every wave trace. A 2x2 boundary solve
//! prepares axisymmetric reflection; a 6x6 solve handles a full support tensor.
//! Both are cached per tangential slowness, then coherent drives allocate nothing.
//! This is a bulk plane-wave boundary, not a finite plate, flexural edge,
//! Rayleigh-wave source, nonlocal support or surface-roughness model. Supplied
//! impedances belong to the same frequency as the drive. No causal time-domain
//! realization or experimental qualification is inferred from one frequency.
//! Free-surface analytical coefficients: Sánchez-Sesma et al. (2023), Appendix 2,
//! <https://doi.org/10.1186/s40623-023-01780-0>.

use super::*;

pub(super) mod tensor;
pub use tensor::ElasticSurfaceImpedanceTensor;

/// Passive locally reacting support in the original `exp(-i omega t)` convention.
///
/// Each impedance is traction opposing motion divided by velocity [Pa s/m].
/// Tangential response is isotropic in the surface plane. Nonnegative real
/// parts dissipate mechanical work; either sign of reactance is admitted.
/// The caller evaluates any dispersive material/support law at the drive's
/// original frequency. Zero values recover a free surface exactly.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElasticSurfaceImpedance {
    normal_pa_s_m: C64,
    tangential_pa_s_m: C64,
}

impl ElasticSurfaceImpedance {
    /// Admit finite normal/tangential impedances with nonnegative resistance.
    pub fn new(normal_pa_s_m: C64, tangential_pa_s_m: C64) -> Result<Self, MaterialError> {
        if [normal_pa_s_m, tangential_pa_s_m]
            .into_iter()
            .any(|z| !finite(z) || z.re < 0.0)
        {
            return Err(invalid(
                "elastic surface impedance must be finite and passive",
            ));
        }
        Ok(Self {
            normal_pa_s_m,
            tangential_pa_s_m,
        })
    }

    /// Specific impedance opposing normal motion [Pa s/m].
    pub const fn normal_pa_s_m(self) -> C64 {
        self.normal_pa_s_m
    }

    /// Specific impedance opposing either tangential component [Pa s/m].
    pub const fn tangential_pa_s_m(self) -> C64 {
        self.tangential_pa_s_m
    }

    fn component(self, axis: usize) -> C64 {
        if axis == 2 {
            self.normal_pa_s_m
        } else {
            self.tangential_pa_s_m
        }
    }

    // Bounded coefficients of a*t/Zref + b*v=0. They retain the exact zero
    // and stiff limits without constructing a huge complex norm or inverse.
    fn coordinates(self, axis: usize, reference: f64) -> (f64, C64) {
        let z = self.component(axis);
        let scale = reference.max(z.re).max(z.im.abs());
        (reference / scale, divided(z, scale))
    }
}

/// Explicit stationary mechanical boundary condition.
#[derive(Clone, Copy, Debug, PartialEq)]
#[allow(clippy::large_enum_variant)] // Fixed storage keeps prepared drives allocation-free.
pub enum ElasticSurfaceCondition {
    /// Vacuum/free boundary: `sigma.e_z = 0`, with unconstrained velocity.
    TractionFree,
    /// Ideal immovable support: `v = 0`, with unconstrained reaction traction.
    Clamped,
    /// Passive support with `traction = -Z velocity` on the outward normal.
    Impedance(ElasticSurfaceImpedance),
    /// Reciprocal passive tensor support, including normal/tangential coupling.
    TensorImpedance(ElasticSurfaceImpedanceTensor),
}

/// Original isotropic material and its cached P/shear speeds and impedances.
#[derive(Debug)]
pub struct IsotropicElasticSurface<'a> {
    state: &'a ElasticAcousticState,
    condition: ElasticSurfaceCondition,
    speeds: [f64; 2],
    impedances: [f64; 2],
}

/// Fixed-size reflection operator at a conserved tangential slowness.
#[derive(Debug)]
pub struct PreparedElasticSurface<'a> {
    state: &'a ElasticAcousticState,
    condition: ElasticSurfaceCondition,
    tangential_slowness_s_m: f64,
    incoming: [ElasticWaveTrace; 3],
    reflected: [ElasticWaveTrace; 3],
    scattering: [[C64; 3]; 3],
    impedance_reference_pa_s_m: f64,
    // A coupled support solves wall motion directly, retaining weak motion
    // beside stiff axes without subtracting nearly opposite modal velocities.
    tensor_boundary_velocity: Option<[[C64; 3]; 3]>,
}

/// Coherent physical fields and independently reconstructed surface work.
/// Modal arrays are ordered P, SV, SH; powers are per physical surface area.
#[derive(Clone, Debug)]
pub struct ElasticSurfaceResponse {
    pub material_state_identity: ContentHash,
    pub incoming_velocity_m_s: [C64; 3],
    pub reflected_velocity_m_s: [C64; 3],
    pub incoming_normal_power_w_m2: [f64; 3],
    pub reflected_normal_power_w_m2: [f64; 3],
    pub total_incoming_power_w_m2: f64,
    pub total_reflected_power_w_m2: f64,
    /// Positive local support work `Re(v* . Z v)/2` [W/m²], evaluated through
    /// the admitted resistance factor for a tensor rather than a power deficit.
    /// Reactive supports have exactly zero loss.
    pub absorbed_normal_power_w_m2: f64,
    /// Coherent boundary velocity, retained through both traction/velocity
    /// coordinates or a coupled support solve to preserve weak wall motion.
    pub boundary_velocity_m_s: [C64; 3],
    /// Traction `sigma.e_z` on the outward normal, in Pa.
    pub boundary_traction_pa: [C64; 3],
    /// Independently evaluated coherent `-Re(v*.traction)/2` before boundary
    /// reconstruction; zero for ideal or purely reactive conditions.
    pub boundary_normal_flux_w_m2: f64,
    /// Conservative Mandel-strain norm bound throughout the solid half-space.
    pub strain_norm_bound: f64,
}

impl<'a> IsotropicElasticSurface<'a> {
    /// Retain the original stiffness, density, identity and strain domain.
    /// An anisotropic source is never projected to an isotropic approximation.
    pub fn new(
        state: &'a ElasticAcousticState,
        condition: ElasticSurfaceCondition,
    ) -> Result<Self, MaterialError> {
        if !state.is_explicitly_isotropic() {
            return Err(invalid(
                "elastic surface requires an explicitly isotropic state",
            ));
        }
        let modes = state.directional_waves([0.0, 0.0, 1.0])?.modes;
        Ok(Self {
            state,
            condition,
            speeds: [modes[2].phase_speed_m_s, modes[0].phase_speed_m_s],
            impedances: [modes[2].impedance_pa_s_m, modes[0].impedance_pa_s_m],
        })
    }

    pub const fn state(&self) -> &'a ElasticAcousticState {
        self.state
    }

    /// P and shear phase speeds [m/s], in that order.
    pub const fn phase_speeds_m_s(&self) -> [f64; 2] {
        self.speeds
    }

    /// Fix tangential slowness using a propagating incident mode. Every other
    /// propagating channel at that slowness can subsequently be driven too.
    pub fn prepare_incidence(
        &self,
        mode: ElasticMode,
        cosine: f64,
    ) -> Result<PreparedElasticSurface<'a>, MaterialError> {
        self.prepare_incidence_cancellable(mode, cosine, &mut || false)
    }

    /// Bounded preparation with cancellation before traces and each solve.
    /// Exact grazing incidence and singular boundary systems refuse.
    pub fn prepare_incidence_cancellable(
        &self,
        mode: ElasticMode,
        cosine: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedElasticSurface<'a>, MaterialError> {
        poll(cancelled)?;
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(invalid("elastic surface incident cosine must be in (0,1]"));
        }
        let incident_speed = self.speeds[mode.index().min(1)];
        let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine));
        let tangential = sine / incident_speed;
        if !tangential.is_finite() || (sine > 0.0 && tangential == 0.0) {
            return Err(invalid(
                "elastic surface tangential slowness is unrepresentable",
            ));
        }
        let mut result = PreparedElasticSurface {
            state: self.state,
            condition: self.condition,
            tangential_slowness_s_m: tangential,
            incoming: [empty_trace(); 3],
            reflected: [empty_trace(); 3],
            scattering: [[C64::ZERO; 3]; 3],
            impedance_reference_pa_s_m: self.impedances[0].max(self.impedances[1]),
            tensor_boundary_velocity: None,
        };
        for kind in 0..3 {
            poll(cancelled)?;
            let speed = self.speeds[kind.min(1)];
            let normal = normal_ratio(incident_speed, speed, cosine)?;
            let ratio = speed / incident_speed;
            let wave_cosine = normal.scale(ratio);
            let wave_sine = sine * ratio;
            let normal_slowness = divided(normal, incident_speed);
            if !finite(wave_cosine)
                || !wave_sine.is_finite()
                || (normal != C64::ZERO && normal_slowness == C64::ZERO)
            {
                return Err(invalid("elastic surface wave direction is unrepresentable"));
            }
            for (direction, traces) in [(1.0, &mut result.incoming), (-1.0, &mut result.reflected)]
            {
                traces[kind] = trace(
                    self.state,
                    self.impedances,
                    kind,
                    direction,
                    wave_sine,
                    wave_cosine,
                    tangential,
                    normal_slowness,
                )?;
            }
        }
        if let ElasticSurfaceCondition::TensorImpedance(support) = self.condition {
            result.prepare_tensor(support, cancelled)?;
            poll(cancelled)?;
            return Ok(result);
        }
        result.prepare_psv(cancelled)?;
        if result.incoming[2].normal_power_factor_pa_s_m > 0.0 {
            result.scattering[2][2] = match self.condition {
                ElasticSurfaceCondition::TractionFree => C64::ONE,
                ElasticSurfaceCondition::Clamped => -C64::ONE,
                ElasticSurfaceCondition::Impedance(_) => {
                    let incoming = result.boundary_entry(result.incoming[2], 1);
                    let reflected = result.boundary_entry(result.reflected[2], 1);
                    let coefficient = -incoming / reflected;
                    if !finite(coefficient) {
                        return Err(invalid("elastic SH support reflection is unrepresentable"));
                    }
                    coefficient
                }
                ElasticSurfaceCondition::TensorImpedance(_) => {
                    unreachable!("tensor reflection uses the coupled preparation")
                }
            };
        }
        poll(cancelled)?;
        Ok(result)
    }
}

impl PreparedElasticSurface<'_> {
    pub const fn condition(&self) -> ElasticSurfaceCondition {
        self.condition
    }

    pub const fn tangential_slowness_s_m(&self) -> f64 {
        self.tangential_slowness_s_m
    }

    /// Evanescent reflected modes remain physical fields with zero normal power.
    pub fn reflected_mode(&self, mode: ElasticMode) -> ElasticWaveTrace {
        self.reflected[mode.index()]
    }

    /// Only propagating incoming modes are available for external excitation.
    pub fn incoming_mode(&self, mode: ElasticMode) -> Option<ElasticWaveTrace> {
        let wave = self.incoming[mode.index()];
        (wave.normal_power_factor_pa_s_m > 0.0).then_some(wave)
    }

    /// Reflected P/SV/SH modal velocities per unit incident modal velocity.
    pub fn velocity_coefficients(&self, mode: ElasticMode) -> Option<[C64; 3]> {
        self.incoming_mode(mode)?;
        Some(core::array::from_fn(|i| self.scattering[i][mode.index()]))
    }

    /// Apply coherent peak modal velocities [m/s] at a common frequency/phase
    /// reference. No allocation, factorization, clipping or phase averaging.
    pub fn drive(&self, incoming: [C64; 3]) -> Result<ElasticSurfaceResponse, MaterialError> {
        for (i, &amplitude) in incoming.iter().enumerate() {
            if !finite(amplitude)
                || (amplitude != C64::ZERO && self.incoming[i].normal_power_factor_pa_s_m == 0.0)
            {
                return Err(invalid(
                    "elastic surface drive needs finite propagating amplitudes",
                ));
            }
        }
        let reflected: [C64; 3] = core::array::from_fn(|i| {
            complex_sum((0..3).map(|j| self.scattering[i][j] * incoming[j]))
        });
        if reflected.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "elastic surface reflected amplitude is unrepresentable",
            ));
        }
        let mut input_power = [0.0; 3];
        let mut reflected_power = [0.0; 3];
        let mut strain = Dd::ZERO;
        for mode in 0..3 {
            input_power[mode] = power(self.incoming[mode], incoming[mode])?;
            reflected_power[mode] = power(self.reflected[mode], reflected[mode])?;
            for (wave, amplitude) in [
                (self.incoming[mode], incoming[mode]),
                (self.reflected[mode], reflected[mode]),
            ] {
                let norm = wave
                    .strain_mandel_s_m
                    .iter()
                    .fold(0.0, |n, e| det::hypot(n, e.abs()));
                strain = strain + Dd::from_f64(norm * amplitude.abs());
            }
        }
        let strain_bound = strain.to_f64();
        if !strain_bound.is_finite() || strain_bound > self.state.maximum_absolute_mandel_strain() {
            return Err(invalid(
                "elastic surface drive exceeds original strain domain",
            ));
        }
        let mut velocity: [C64; 3] = core::array::from_fn(|component| {
            complex_sum((0..3).flat_map(|m| {
                [
                    self.incoming[m].velocity[component] * incoming[m],
                    self.reflected[m].velocity[component] * reflected[m],
                ]
            }))
        });
        let mut traction: [C64; 3] = core::array::from_fn(|component| {
            complex_sum((0..3).flat_map(|m| {
                [
                    self.incoming[m].traction_pa_s_m[component] * incoming[m],
                    self.reflected[m].traction_pa_s_m[component] * reflected[m],
                ]
            }))
        });
        let total_incoming = real_sum(input_power);
        let total_reflected = real_sum(reflected_power);
        let flux = -0.5 * real_sum((0..3).map(|i| (velocity[i].conj() * traction[i]).re));
        let mut absorbed = 0.0;
        if let ElasticSurfaceCondition::Impedance(support) = self.condition {
            let reference = self.impedance_reference_pa_s_m;
            let mut losses = [0.0; 3];
            for axis in 0..3 {
                let (a, b) = support.coordinates(axis, reference);
                // Project both independently summed traces onto their original
                // constitutive coordinates. At a stiff support the traction
                // retains the small wall velocity that v_in+v_ref would lose;
                // a weak support instead retains the well-resolved velocity.
                let amplitude = (velocity[axis].scale(a)
                    - b.conj() * divided(traction[axis], reference))
                    / C64::from_re(a * a + b.norm_sq());
                velocity[axis] = amplitude.scale(a);
                traction[axis] = -(b * amplitude).scale(reference);
                let resistance = support.component(axis).re;
                let root =
                    det::sqrt(resistance) * core::f64::consts::FRAC_1_SQRT_2 * velocity[axis].abs();
                losses[axis] = root * root;
                if resistance > 0.0 && velocity[axis] != C64::ZERO && losses[axis] == 0.0 {
                    return Err(invalid("elastic support absorption is unrepresentable"));
                }
            }
            absorbed = real_sum(losses);
        }
        if let ElasticSurfaceCondition::TensorImpedance(support) = self.condition {
            let map = self
                .tensor_boundary_velocity
                .as_ref()
                .expect("tensor boundary was prepared together with reflection");
            velocity = core::array::from_fn(|axis| {
                tensor::complex_dot((0..3).map(|mode| (map[axis][mode], incoming[mode])))
            });
            traction = support.traction(velocity);
            absorbed = support.absorbed_power(velocity)?;
        }
        if velocity.iter().chain(&traction).any(|&x| !finite(x))
            || [total_incoming, total_reflected, flux, absorbed]
                .into_iter()
                .any(|x| !x.is_finite())
            || !balanced(total_incoming, total_reflected + absorbed, total_incoming)
            || !balanced(flux, absorbed, total_incoming)
        {
            return Err(invalid(
                "elastic surface power is unrepresentable or does not close",
            ));
        }
        Ok(ElasticSurfaceResponse {
            material_state_identity: self.state.material_state_identity(),
            incoming_velocity_m_s: incoming,
            reflected_velocity_m_s: reflected,
            incoming_normal_power_w_m2: input_power,
            reflected_normal_power_w_m2: reflected_power,
            total_incoming_power_w_m2: total_incoming,
            total_reflected_power_w_m2: total_reflected,
            absorbed_normal_power_w_m2: absorbed,
            boundary_velocity_m_s: velocity,
            boundary_traction_pa: traction,
            boundary_normal_flux_w_m2: flux,
            strain_norm_bound: strain_bound,
        })
    }

    fn boundary_entry(&self, wave: ElasticWaveTrace, component: usize) -> C64 {
        match self.condition {
            ElasticSurfaceCondition::TractionFree => wave.traction_pa_s_m[component],
            ElasticSurfaceCondition::Clamped => wave.velocity[component],
            ElasticSurfaceCondition::Impedance(support) => {
                let (a, b) = support.coordinates(component, self.impedance_reference_pa_s_m);
                divided(
                    wave.traction_pa_s_m[component],
                    self.impedance_reference_pa_s_m,
                )
                .scale(a)
                    + b * wave.velocity[component]
            }
            ElasticSurfaceCondition::TensorImpedance(_) => {
                unreachable!("tensor support retains coupled velocity coordinates")
            }
        }
    }

    fn prepare_psv(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialError> {
        let mut matrix = [C64::ZERO; 4];
        let mut scales = [0.0_f64; 2];
        for row in 0..2 {
            for column in 0..2 {
                let entry = self.boundary_entry(self.reflected[column], 2 * row);
                matrix[2 * row + column] = entry;
                scales[row] = scales[row].max(entry.abs());
            }
            if scales[row] == 0.0 || !scales[row].is_finite() {
                return Err(invalid(
                    "elastic surface boundary row is singular or unrepresentable",
                ));
            }
            for column in 0..2 {
                matrix[2 * row + column] = divided(matrix[2 * row + column], scales[row]);
            }
        }
        let factor = lu_complex_with_cancel(&matrix, 2, &mut *cancelled).map_err(solve_error)?;
        for incident in 0..2 {
            if self.incoming[incident].normal_power_factor_pa_s_m == 0.0 {
                continue;
            }
            let rhs: [C64; 2] = core::array::from_fn(|row| {
                -divided(
                    self.boundary_entry(self.incoming[incident], 2 * row),
                    scales[row],
                )
            });
            let mut amplitudes = rhs;
            factor
                .solve_with_cancel(&mut amplitudes, &mut *cancelled)
                .map_err(solve_error)?;
            if amplitudes.iter().any(|&a| !finite(a)) {
                return Err(invalid(
                    "elastic surface reflection solve is unrepresentable",
                ));
            }
            for row in 0..2 {
                let terms: [C64; 2] = core::array::from_fn(|j| matrix[2 * row + j] * amplitudes[j]);
                let residual = (complex_sum(terms) - rhs[row]).abs();
                let magnitude = real_sum(terms.map(C64::abs)) + rhs[row].abs();
                if !residual.is_finite() || (magnitude > 0.0 && residual / magnitude > 1e-11) {
                    return Err(invalid("elastic surface boundary solve did not converge"));
                }
                self.scattering[row][incident] = amplitudes[row];
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
