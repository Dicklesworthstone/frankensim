//! Coherent plane waves at a welded, isotropic elastic interface.
//!
//! The local frame is `(x, y, z)`: x is the conserved tangential propagation
//! direction, y the SH direction, and z points from the first solid into the
//! second. Phasors are peak amplitudes with `exp(i omega (s.x - t))`.
//! Continuity of ALL velocity and traction components supplies the exact
//! Zoeppritz system. P/SV share one 4x4 factorization; SH is independent.
//! No small-angle or weak-contrast approximation is made.
//!
//! Outgoing normal slowness has positive imaginary part into its own solid.
//! Evanescent waves retain complex polarization and traction but carry zero
//! individual normal power. They are not absorption. Incoming excitation must
//! propagate; finite-source near fields and head waves need a spatial solver.
//! See Zoeppritz, *On the Reflection and Transmission of Seismic Waves at
//! Surfaces of Discontinuity*, <https://doi.org/10.1190/1.9781560801931.ch3n>,
//! and the CREWES energy-flux discussion:
//! <https://www.crewes.org/ResearchLinks/ExplorerPrograms/EnergyEx/>.

use super::{ElasticAcousticState, invalid};
use crate::MaterialError;
use fs_blake3::ContentHash;
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};

/// Coherent scattering and interior fields of a finite bonded elastic slab.
pub mod layer;

/// Physical side of the interface. The common normal points first to second.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElasticSide {
    First,
    Second,
}
impl ElasticSide {
    const fn index(self) -> usize {
        match self {
            Self::First => 0,
            Self::Second => 1,
        }
    }
}

/// Local velocity polarizations, in the order used by three-element arrays.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElasticMode {
    /// Longitudinal, parallel to the (possibly complex) wave normal.
    P,
    /// Shear in the xz incidence plane.
    Sv,
    /// Shear in the y direction, uncoupled from P/SV.
    Sh,
}
impl ElasticMode {
    const fn index(self) -> usize {
        match self {
            Self::P => 0,
            Self::Sv => 1,
            Self::Sh => 2,
        }
    }
}

/// One wave at the common interface, per unit modal peak velocity [m/s].
/// For propagation sign d, P polarization is `(sin(theta),0,d*cos(theta))`,
/// SV is `(cos(theta),0,-d*sin(theta))`, and SH is `(0,1,0)`. The same formulas
/// continue analytically for evanescence; their complex norms need not be one.
#[derive(Clone, Copy, Debug)]
pub struct ElasticWaveTrace {
    /// Signed z slowness [s/m], including the outgoing evanescent branch.
    pub normal_slowness_s_m: C64,
    /// Local Cartesian velocity divided by modal velocity amplitude.
    pub velocity: [C64; 3],
    /// `sigma . e_z` divided by modal velocity amplitude [Pa s/m].
    pub traction_pa_s_m: [C64; 3],
    /// Mandel strain divided by modal velocity amplitude [s/m].
    pub strain_mandel_s_m: [C64; 6],
    /// Magnitude of normal mean flux per squared modal velocity,
    /// `rho*c*Re(cos(theta))/2` [Pa s/m]. Zero for grazing/evanescent modes.
    pub normal_power_factor_pa_s_m: f64,
}

/// Original lossless isotropic states, with their Christoffel speeds cached.
/// No material source is reselected, projected to isotropy, or assigned damping.
#[derive(Debug)]
pub struct IsotropicElasticInterface<'a> {
    states: [&'a ElasticAcousticState; 2],
    speeds: [[f64; 2]; 2],
    impedances: [[f64; 2]; 2],
}

/// Six-channel response at one conserved tangential slowness. Channel ordering
/// is `[first P, first SV, first SH, second P, second SV, second SH]`.
/// Preparation uses bounded small LU storage; subsequent `drive` calls allocate
/// nothing and use the retained amplitudes, fields and original strain limits.
#[derive(Debug)]
pub struct PreparedElasticInterface<'a> {
    states: [&'a ElasticAcousticState; 2],
    tangential_slowness_s_m: f64,
    incoming: [ElasticWaveTrace; 6],
    outgoing: [ElasticWaveTrace; 6],
    // Columns for grazing/evanescent incoming channels are not populated and
    // never exposed as admitted scattering columns.
    scattering: [[C64; 6]; 6],
}

/// Physical coherent fields for the two specified incoming P/SV/SH triples.
/// All powers are per unit interface area, and both boundary fluxes use +z.
#[derive(Clone, Debug)]
pub struct ElasticInterfaceResponse {
    pub material_state_identities: [ContentHash; 2],
    pub incoming_velocity_m_s: [[C64; 3]; 2],
    pub outgoing_velocity_m_s: [[C64; 3]; 2],
    pub incoming_normal_power_w_m2: [[f64; 3]; 2],
    pub outgoing_normal_power_w_m2: [[f64; 3]; 2],
    pub total_incoming_power_w_m2: f64,
    pub total_outgoing_power_w_m2: f64,
    pub boundary_velocity_m_s: [[C64; 3]; 2],
    /// Traction on the SAME +z normal for both sides, not outward tractions.
    pub boundary_traction_pa: [[C64; 3]; 2],
    /// Independently reconstructed `-Re(v* . traction)/2` on each side.
    pub boundary_normal_flux_w_m2: [f64; 2],
    /// Triangle bound on the Mandel strain norm over each entire half-space.
    /// This conservatively enforces the source's componentwise strain ceiling,
    /// including phases away from the interface and arbitrary frame rotations.
    pub strain_norm_bound: [f64; 2],
}

impl<'a> IsotropicElasticInterface<'a> {
    pub fn new(
        first: &'a ElasticAcousticState,
        second: &'a ElasticAcousticState,
    ) -> Result<Self, MaterialError> {
        if !first.is_explicitly_isotropic() || !second.is_explicitly_isotropic() {
            return Err(invalid(
                "elastic interface requires two explicitly isotropic states",
            ));
        }
        let states = [first, second];
        let mut speeds = [[0.0; 2]; 2];
        let mut impedances = [[0.0; 2]; 2];
        for side in 0..2 {
            let modes = states[side].directional_waves([0.0, 0.0, 1.0])?.modes;
            speeds[side] = [modes[2].phase_speed_m_s, modes[0].phase_speed_m_s];
            impedances[side] = [modes[2].impedance_pa_s_m, modes[0].impedance_pa_s_m];
        }
        Ok(Self {
            states,
            speeds,
            impedances,
        })
    }

    pub fn state(&self, side: ElasticSide) -> &'a ElasticAcousticState {
        self.states[side.index()]
    }

    /// Original P and shear speeds [m/s], in that order.
    pub const fn phase_speeds_m_s(&self, side: ElasticSide) -> [f64; 2] {
        self.speeds[side.index()]
    }

    /// Specify tangential slowness through one physically propagating incident
    /// mode and its cosine in `(0,1]`. The resulting operator admits every
    /// propagating incoming channel on BOTH sides at that same slowness.
    /// Reflected/transmitted critical modes remain finite traction states;
    /// exact grazing INCIDENT power is outside this normalization.
    pub fn prepare_incidence(
        &self,
        side: ElasticSide,
        mode: ElasticMode,
        cosine: f64,
    ) -> Result<PreparedElasticInterface<'a>, MaterialError> {
        self.prepare_incidence_cancellable(side, mode, cosine, &mut || false)
    }

    pub fn prepare_incidence_cancellable(
        &self,
        side: ElasticSide,
        mode: ElasticMode,
        cosine: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedElasticInterface<'a>, MaterialError> {
        let mut result = self.prepare_wave_traces(side, mode, cosine, cancelled)?;
        // A fictitious split of the identical mechanical medium scatters no
        // wave. This exact identity also removes the undriven homogeneous
        // grazing channel that makes the unreduced boundary matrix singular.
        if self.states[0].density_kg_m3() == self.states[1].density_kg_m3()
            && self.states[0].stiffness_mandel_pa() == self.states[1].stiffness_mandel_pa()
        {
            for channel in 0..6 {
                if result.incoming[channel].normal_power_factor_pa_s_m > 0.0 {
                    result.scattering[(channel + 3) % 6][channel] = C64::ONE;
                }
            }
            poll(cancelled)?;
            return Ok(result);
        }
        result.prepare_psv(cancelled)?;
        result.prepare_sh()?;
        poll(cancelled)?;
        Ok(result)
    }

    // Finite layers use the same exterior traces, including critical branches,
    // without first solving a fictitious direct contact between their faces.
    fn prepare_wave_traces(
        &self,
        side: ElasticSide,
        mode: ElasticMode,
        cosine: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedElasticInterface<'a>, MaterialError> {
        poll(cancelled)?;
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(invalid("elastic incident cosine must be in (0,1]"));
        }
        let incident_speed = self.speeds[side.index()][mode.index().min(1)];
        let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine));
        let tangential = sine / incident_speed;
        if !tangential.is_finite() || (sine > 0.0 && tangential == 0.0) {
            return Err(invalid("elastic tangential slowness is unrepresentable"));
        }
        let mut incoming = [empty_trace(); 6];
        let mut outgoing = [empty_trace(); 6];
        for material in 0..2 {
            poll(cancelled)?;
            for kind in 0..3 {
                let speed = self.speeds[material][kind.min(1)];
                let normal = normal_ratio(incident_speed, speed, cosine)?;
                let ratio = speed / incident_speed;
                let wave_cosine = normal.scale(ratio);
                let wave_sine = sine * ratio;
                let normal_slowness = divided(normal, incident_speed);
                if !finite(wave_cosine)
                    || !wave_sine.is_finite()
                    || (normal != C64::ZERO && normal_slowness == C64::ZERO)
                {
                    return Err(invalid("elastic wave direction is unrepresentable"));
                }
                let index = material * 3 + kind;
                let direction = if material == 0 { -1.0 } else { 1.0 };
                outgoing[index] = trace(
                    self.states[material],
                    self.impedances[material],
                    kind,
                    direction,
                    wave_sine,
                    wave_cosine,
                    tangential,
                    normal_slowness,
                )?;
                incoming[index] = trace(
                    self.states[material],
                    self.impedances[material],
                    kind,
                    -direction,
                    wave_sine,
                    wave_cosine,
                    tangential,
                    normal_slowness,
                )?;
            }
        }
        Ok(PreparedElasticInterface {
            states: self.states,
            tangential_slowness_s_m: tangential,
            incoming,
            outgoing,
            scattering: [[C64::ZERO; 6]; 6],
        })
    }
}

impl PreparedElasticInterface<'_> {
    pub const fn tangential_slowness_s_m(&self) -> f64 {
        self.tangential_slowness_s_m
    }

    pub fn material_state_identities(&self) -> [ContentHash; 2] {
        self.states
            .map(ElasticAcousticState::material_state_identity)
    }

    /// Every outgoing mode exists as a field, including zero-flux evanescence.
    pub fn outgoing_mode(&self, side: ElasticSide, mode: ElasticMode) -> ElasticWaveTrace {
        self.outgoing[3 * side.index() + mode.index()]
    }

    /// An incoming channel is admitted only if it transports positive normal
    /// power toward the interface. No evanescent illumination is inferred.
    pub fn incoming_mode(&self, side: ElasticSide, mode: ElasticMode) -> Option<ElasticWaveTrace> {
        let trace = self.incoming[3 * side.index() + mode.index()];
        (trace.normal_power_factor_pa_s_m > 0.0).then_some(trace)
    }

    /// Outgoing modal velocity coefficients for unit incident modal velocity.
    /// Arrays are `[first P,SV,SH]` and `[second P,SV,SH]`. These infinitesimal
    /// coefficients are independent of a finite excitation's strain admission.
    pub fn velocity_coefficients(
        &self,
        side: ElasticSide,
        mode: ElasticMode,
    ) -> Option<[[C64; 3]; 2]> {
        self.incoming_mode(side, mode)?;
        let column = 3 * side.index() + mode.index();
        Some(core::array::from_fn(|s| {
            core::array::from_fn(|m| self.scattering[3 * s + m][column])
        }))
    }

    /// Apply six coherent peak velocity amplitudes [m/s]. Zero excitation in
    /// a grazing/evanescent incoming channel is allowed; nonzero excitation
    /// refuses. Amplitudes share one frequency and the same interface phase
    /// reference. There is no hidden phase averaging or power renormalization.
    pub fn drive(
        &self,
        incoming: [[C64; 3]; 2],
    ) -> Result<ElasticInterfaceResponse, MaterialError> {
        let amplitudes: [C64; 6] = core::array::from_fn(|i| incoming[i / 3][i % 3]);
        for (i, &a) in amplitudes.iter().enumerate() {
            if !finite(a) || (a != C64::ZERO && self.incoming[i].normal_power_factor_pa_s_m == 0.0)
            {
                return Err(invalid(
                    "elastic drive needs finite propagating incident amplitudes",
                ));
            }
        }
        let outgoing: [C64; 6] = core::array::from_fn(|i| {
            complex_sum((0..6).map(|j| self.scattering[i][j] * amplitudes[j]))
        });
        if outgoing.iter().any(|&a| !finite(a)) {
            return Err(invalid("elastic outgoing amplitudes are unrepresentable"));
        }
        let mut input_power = [[0.0; 3]; 2];
        let mut output_power = [[0.0; 3]; 2];
        let mut velocity = [[C64::ZERO; 3]; 2];
        let mut traction = [[C64::ZERO; 3]; 2];
        let mut strain_bound = [0.0; 2];
        for side in 0..2 {
            let mut strain_sum = Dd::from_f64(0.0);
            for mode in 0..3 {
                let i = 3 * side + mode;
                input_power[side][mode] = power(self.incoming[i], amplitudes[i])?;
                output_power[side][mode] = power(self.outgoing[i], outgoing[i])?;
                for (trace, a) in [
                    (self.incoming[i], amplitudes[i]),
                    (self.outgoing[i], outgoing[i]),
                ] {
                    let norm = trace
                        .strain_mandel_s_m
                        .iter()
                        .fold(0.0, |n, e| det::hypot(n, e.abs()));
                    strain_sum = strain_sum + Dd::from_f64(norm * a.abs());
                }
            }
            strain_bound[side] = strain_sum.to_f64();
            if !strain_bound[side].is_finite()
                || strain_bound[side] > self.states[side].maximum_absolute_mandel_strain()
            {
                return Err(invalid(
                    "elastic interface drive exceeds original strain domain",
                ));
            }
            for component in 0..3 {
                velocity[side][component] = complex_sum((0..3).flat_map(|m| {
                    let i = 3 * side + m;
                    [
                        self.incoming[i].velocity[component] * amplitudes[i],
                        self.outgoing[i].velocity[component] * outgoing[i],
                    ]
                }));
                traction[side][component] = complex_sum((0..3).flat_map(|m| {
                    let i = 3 * side + m;
                    [
                        self.incoming[i].traction_pa_s_m[component] * amplitudes[i],
                        self.outgoing[i].traction_pa_s_m[component] * outgoing[i],
                    ]
                }));
            }
        }
        let total_incoming = real_sum(input_power.into_iter().flatten());
        let total_outgoing = real_sum(output_power.into_iter().flatten());
        let flux = core::array::from_fn(|s| {
            -0.5 * real_sum((0..3).map(|i| (velocity[s][i].conj() * traction[s][i]).re))
        });
        if velocity
            .iter()
            .chain(&traction)
            .flatten()
            .any(|&v| !finite(v))
            || [total_incoming, total_outgoing]
                .into_iter()
                .chain(flux)
                .any(|x| !x.is_finite())
            || !balanced(total_incoming, total_outgoing, total_incoming)
            || !balanced(
                flux[0],
                input_power[0].iter().sum::<f64>() - output_power[0].iter().sum::<f64>(),
                total_incoming,
            )
            || !balanced(
                flux[1],
                output_power[1].iter().sum::<f64>() - input_power[1].iter().sum::<f64>(),
                total_incoming,
            )
        {
            return Err(invalid(
                "elastic interface power is unrepresentable or does not close",
            ));
        }
        Ok(ElasticInterfaceResponse {
            material_state_identities: self.material_state_identities(),
            incoming_velocity_m_s: incoming,
            outgoing_velocity_m_s: core::array::from_fn(|s| {
                core::array::from_fn(|m| outgoing[3 * s + m])
            }),
            incoming_normal_power_w_m2: input_power,
            outgoing_normal_power_w_m2: output_power,
            total_incoming_power_w_m2: total_incoming,
            total_outgoing_power_w_m2: total_outgoing,
            boundary_velocity_m_s: velocity,
            boundary_traction_pa: traction,
            boundary_normal_flux_w_m2: flux,
            strain_norm_bound: strain_bound,
        })
    }

    fn prepare_psv(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialError> {
        const CHANNELS: [usize; 4] = [0, 1, 3, 4];
        let mut matrix = [C64::ZERO; 16];
        let mut scales = [0.0_f64; 4];
        for row in 0..4 {
            for (column, &channel) in CHANNELS.iter().enumerate() {
                let entry = boundary_entry(self.outgoing[channel], row);
                matrix[4 * row + column] = if channel < 3 { entry } else { -entry };
                scales[row] = scales[row].max(entry.abs());
            }
            if scales[row] == 0.0 || !scales[row].is_finite() {
                return Err(invalid("elastic interface boundary row is unrepresentable"));
            }
            for column in 0..4 {
                matrix[4 * row + column] = divided(matrix[4 * row + column], scales[row]);
            }
        }
        let factor = lu_complex_with_cancel(&matrix, 4, &mut *cancelled).map_err(solve_error)?;
        for &incident in &CHANNELS {
            if self.incoming[incident].normal_power_factor_pa_s_m == 0.0 {
                continue;
            }
            let rhs: [C64; 4] = core::array::from_fn(|row| {
                divided(boundary_entry(self.incoming[incident], row), scales[row])
                    .scale(if incident < 3 { -1.0 } else { 1.0 })
            });
            let mut amplitudes = rhs;
            factor
                .solve_with_cancel(&mut amplitudes, &mut *cancelled)
                .map_err(solve_error)?;
            if amplitudes.iter().any(|&a| !finite(a)) {
                return Err(invalid("elastic interface solve is unrepresentable"));
            }
            for row in 0..4 {
                let terms: [C64; 4] = core::array::from_fn(|j| matrix[4 * row + j] * amplitudes[j]);
                let residual = (complex_sum(terms) - rhs[row]).abs();
                let magnitude = real_sum(terms.map(C64::abs)) + rhs[row].abs();
                if !residual.is_finite() || (magnitude > 0.0 && residual / magnitude > 1e-11) {
                    return Err(invalid(
                        "elastic interface continuity solve did not converge",
                    ));
                }
            }
            for (i, &out) in CHANNELS.iter().enumerate() {
                self.scattering[out][incident] = amplitudes[i];
            }
        }
        Ok(())
    }

    fn prepare_sh(&mut self) -> Result<(), MaterialError> {
        let first = self.outgoing[2].traction_pa_s_m[1];
        let second = -self.outgoing[5].traction_pa_s_m[1];
        let scale = first.abs().max(second.abs());
        let a = divided(first, scale);
        let b = divided(second, scale);
        let denominator = a + b;
        if !scale.is_finite() || scale == 0.0 || !finite(denominator) || denominator == C64::ZERO {
            return Err(invalid("elastic SH boundary is unrepresentable"));
        }
        if self.incoming[2].normal_power_factor_pa_s_m > 0.0 {
            self.scattering[2][2] = (a - b) / denominator;
            self.scattering[5][2] = a.scale(2.0) / denominator;
        }
        if self.incoming[5].normal_power_factor_pa_s_m > 0.0 {
            self.scattering[5][5] = (b - a) / denominator;
            self.scattering[2][5] = b.scale(2.0) / denominator;
        }
        Ok(())
    }
}

// All fields follow sigma = -C:sym(s tensor v), with the SAME original
// stiffness/density as the Christoffel calculation. Analytic isotropy keeps the
// critical normal-traction states without division by a cosine or impedance.
#[allow(clippy::too_many_arguments)]
fn trace(
    state: &ElasticAcousticState,
    impedances: [f64; 2],
    kind: usize,
    direction: f64,
    sine: f64,
    cosine: C64,
    tangential: f64,
    normal: C64,
) -> Result<ElasticWaveTrace, MaterialError> {
    let z = cosine.scale(direction);
    let zero = C64::ZERO;
    let velocity = match kind {
        0 => [C64::from_re(sine), zero, z],
        1 => [cosine, zero, C64::from_re(-direction * sine)],
        _ => [zero, C64::ONE, zero],
    };
    let normal = normal.scale(direction);
    let strain = [
        velocity[0].scale(-tangential),
        zero,
        -(normal * velocity[2]),
        velocity[1].scale(-tangential * core::f64::consts::FRAC_1_SQRT_2),
        (normal * velocity[1]).scale(-core::f64::consts::FRAC_1_SQRT_2),
        (normal * velocity[0] + velocity[2].scale(tangential))
            .scale(-core::f64::consts::FRAC_1_SQRT_2),
    ];
    // Use the retained stiffness directly. For example, lambda=0 at grazing
    // P incidence has O(cosine^2) normal traction; subtracting rounded speed
    // ratios in rho*c_P*(1-2*sin(theta_S)^2) invents a much larger stress.
    let stress_component = |row: usize| {
        complex_sum((0..6).map(|j| strain[j].scale(state.stiffness_mandel_pa()[row][j])))
    };
    let traction = [
        stress_component(5).scale(core::f64::consts::FRAC_1_SQRT_2),
        stress_component(4).scale(core::f64::consts::FRAC_1_SQRT_2),
        stress_component(2),
    ];
    let power_factor = (0.5 * impedances[kind.min(1)]) * cosine.re;
    if velocity
        .iter()
        .chain(&traction)
        .chain(&strain)
        .any(|&x| !finite(x))
        || !finite(normal)
        || !power_factor.is_finite()
        || power_factor < 0.0
        || (cosine.re > 0.0 && power_factor == 0.0)
    {
        return Err(invalid("elastic wave trace is unrepresentable"));
    }
    Ok(ElasticWaveTrace {
        normal_slowness_s_m: normal,
        velocity,
        traction_pa_s_m: traction,
        strain_mandel_s_m: strain,
        normal_power_factor_pa_s_m: power_factor,
    })
}

// k_z/k_incident, retaining the ORIGINAL speed pair at near-critical and
// nearly matched grazing incidence. Only a cancelling critical subtraction
// needs double-double refinement; equal speeds preserve the authored cosine.
fn normal_ratio(incident: f64, transmitted: f64, cosine: f64) -> Result<C64, MaterialError> {
    let ratio = incident / transmitted;
    if !ratio.is_finite() || ratio == 0.0 {
        return Err(invalid("elastic speed contrast is unrepresentable"));
    }
    if cosine == 1.0 {
        return Ok(C64::from_re(ratio));
    }
    if incident == transmitted {
        return Ok(C64::from_re(cosine));
    }
    if ratio > 2.0 {
        let s = det::sqrt((1.0 - cosine) * (1.0 + cosine)) / ratio;
        return Ok(C64::from_re(ratio * det::sqrt(1.0 - s * s)));
    }
    let (mut squared, cancellation_scale) = if cosine > 0.5 {
        let tangential_squared = (1.0 - cosine) * (1.0 + cosine);
        (
            ratio * ratio - tangential_squared,
            ratio * ratio + tangential_squared,
        )
    } else {
        let contrast = ((incident - transmitted) / transmitted) * (ratio + 1.0);
        (
            cosine.mul_add(cosine, contrast),
            cosine * cosine + contrast.abs(),
        )
    };
    if incident < transmitted && squared.abs() <= 64.0 * f64::EPSILON * cancellation_scale {
        let scale =
            f64::from_bits(transmitted.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE);
        let a = Dd::from_f64(incident / scale);
        let b = Dd::from_f64(transmitted / scale);
        let normal = b * Dd::from_f64(cosine);
        squared = (((a - b) * (a + b) + normal * normal) / (b * b)).to_f64();
    }
    Ok(if squared >= 0.0 {
        C64::from_re(det::sqrt(squared))
    } else {
        C64::new(0.0, det::sqrt(-squared))
    })
}

fn empty_trace() -> ElasticWaveTrace {
    ElasticWaveTrace {
        normal_slowness_s_m: C64::ZERO,
        velocity: [C64::ZERO; 3],
        traction_pa_s_m: [C64::ZERO; 3],
        strain_mandel_s_m: [C64::ZERO; 6],
        normal_power_factor_pa_s_m: 0.0,
    }
}
fn boundary_entry(trace: ElasticWaveTrace, row: usize) -> C64 {
    match row {
        0 => trace.velocity[0],
        1 => trace.velocity[2],
        2 => trace.traction_pa_s_m[0],
        _ => trace.traction_pa_s_m[2],
    }
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn divided(value: C64, divisor: f64) -> C64 {
    C64::new(value.re / divisor, value.im / divisor)
}
fn real_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::from_f64(0.0), |s, v| s + Dd::from_f64(v))
        .to_f64()
}
fn complex_sum(values: impl IntoIterator<Item = C64>) -> C64 {
    let (mut re, mut im) = (Dd::from_f64(0.0), Dd::from_f64(0.0));
    for value in values {
        re = re + Dd::from_f64(value.re);
        im = im + Dd::from_f64(value.im);
    }
    C64::new(re.to_f64(), im.to_f64())
}
fn power(trace: ElasticWaveTrace, amplitude: C64) -> Result<f64, MaterialError> {
    if amplitude == C64::ZERO || trace.normal_power_factor_pa_s_m == 0.0 {
        return Ok(0.0);
    }
    let root = det::sqrt(trace.normal_power_factor_pa_s_m) * amplitude.abs();
    let value = root * root;
    if !value.is_finite() || value == 0.0 {
        return Err(invalid("elastic wave power is unrepresentable"));
    }
    Ok(value)
}
fn balanced(first: f64, second: f64, scale: f64) -> bool {
    if scale == 0.0 {
        first == second
    } else {
        (first - second).abs() / scale <= 1e-10
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialError> {
    if cancelled() {
        Err(invalid("elastic interface preparation cancelled"))
    } else {
        Ok(())
    }
}
fn solve_error(error: LuComplexError) -> MaterialError {
    match error {
        LuComplexError::Cancelled => invalid("elastic interface preparation cancelled"),
        LuComplexError::Singular { .. } => invalid("elastic interface boundary system is singular"),
    }
}

#[cfg(test)]
mod tests;
