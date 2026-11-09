//! Mixed frequency/physical-parameter response jets of the original pencil.
//!
//! The source derivatives are authored physical derivatives, not perturbed
//! samples. All solves reuse the admitted factor at the same frequency.

use super::super::{
    HarmonicFrequencyCurvatureControls, HarmonicFrequencyCurvatureResponse, HarmonicResponseError,
    PreparedHarmonicResponse, finite, invalid, poll, scale_complex,
};
use super::{HarmonicParameterControls, HarmonicParameterResponse, PowerSum, Term};
use crate::quadratic_modes::QuadraticParameterDirection;
use fs_exec::Cx;
use fs_math::{c64::C64, dd::Dd};

/// Total resources for the original response, both frequency derivatives,
/// one physical-parameter derivative and their mixed derivative.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicFrequencyParameterControls {
    /// Conservative added numerical payload: 3072*n bytes when unloaded;
    /// 8*n*p+3072*n+1536*p+128*k for p ports and compact-factor dimension k.
    /// Excludes the immutable factor, borrowed source/direction arrays and
    /// allocator metadata. Dense derivative matrices are never allocated.
    pub maximum_additional_payload_bytes: usize,
    /// Conservative work score: 1024*n*n+2048*n when unloaded; loading adds
    /// 512*(n*p+p*p+k*k)+2048*(p+k), plus 512*n*n*p for a compact factor.
    /// Includes all same-factor solves, complete physical directions and
    /// original-equation audits. Owners with a cached LU/compact factor reserve
    /// three times this score for bounded same-factor residual corrections.
    /// Direct scalar/diagonal owners keep the baseline score. A dense physical
    /// direction can still couple an exactly diagonal base.
    pub maximum_work: u64,
}

/// Physical response jet needed to differentiate a driven peak or notch.
///
/// Frequency derivatives are per rad/s. Parameter derivatives are per the
/// caller's declared unit of theta, at fixed physical frequency. The original
/// coordinate basis, mass, damping and stiffness directions remain fixed.
#[derive(Clone, Debug)]
pub struct HarmonicFrequencyParameterResponse {
    /// Original response, frequency slope and frequency curvature.
    pub frequency_response: HarmonicFrequencyCurvatureResponse,
    /// Original response and the fixed-frequency physical parameter tangent.
    pub parameter_response: HarmonicParameterResponse,
    /// d^2 q/(domega dtheta) in the original retained coordinates.
    pub displacement_mixed_derivative: Vec<C64>,
    /// -i*(q_theta + omega*q_omega_theta), including the DC kinematic term.
    pub velocity_mixed_derivative: Vec<C64>,
    /// Maximum same-factor solve and published-coordinate mixed-equation
    /// backward residual. This is not a forward-error bound near a pole.
    pub relative_mixed_derivative_residual: f64,
}

impl PreparedHarmonicResponse {
    /// Evaluate physical frequency/parameter motion jets with the cached factor.
    ///
    /// For A=K-omega^2*M-i*omega*C under exp(-i*omega*t), the mixed equation is
    /// A*q_omega_theta = F_omega_theta - A_omega_theta*q
    ///                  - A_theta*q_omega - A_omega*q_theta.
    /// This retains the full original M/C/K directions, including nonproportional
    /// damping and signed gyroscopic coupling. The mass and stiffness directions
    /// must be symmetric; all directions are finite complete row-major matrices.
    /// Neither a derivative matrix nor another dynamic factor is constructed.
    ///
    /// M/C/K and their parameter directions are independent of frequency, and
    /// theta does not change the drive frequency or physical coordinate basis.
    /// The five source arrays are the physical partial derivatives F, F_omega,
    /// F_omegaomega, F_theta and F_omega_theta at this same frequency/parameter.
    /// None explicitly declares the corresponding source derivative zero.
    /// Source amplitude and phase changes therefore remain in every derivative.
    ///
    /// This producer currently admits unloaded mechanical pencils only. Any
    /// attached port load refuses explicitly: frequency/parameter derivatives
    /// of dispersive impedances and changing port maps cannot be omitted.
    ///
    /// The native curvature and parameter owners provide their existing motion,
    /// power and residual diagnostics. One further normalized solve produces the
    /// mixed tangent, with independent amplitude scaling and an uncancelled
    /// original-equation envelope. Resource, cancellation or numerical refusal
    /// publishes no partial jet and leaves this immutable owner reusable.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_frequency_parameter_derivative(
        &self,
        force: &[C64],
        direction: QuadraticParameterDirection<'_>,
        force_frequency_derivative: Option<&[C64]>,
        force_frequency_second_derivative: Option<&[C64]>,
        force_parameter_derivative: Option<&[C64]>,
        force_frequency_parameter_derivative: Option<&[C64]>,
        controls: HarmonicFrequencyParameterControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicFrequencyParameterResponse, HarmonicResponseError> {
        poll(cx)?;
        if self.port_load.is_some() {
            return Err(invalid(
                "mixed frequency/parameter derivatives require an unloaded mechanical pencil",
            ));
        }
        self.solve_frequency_parameter_jet(
            force,
            direction,
            force_frequency_derivative,
            force_frequency_second_derivative,
            force_parameter_derivative,
            force_frequency_parameter_derivative,
            None,
            None,
            controls,
            cx,
        )
    }

    /// Evaluate physical frequency/parameter jets with a dispersive fixed load.
    ///
    /// The source and mechanical direction arguments have the same meaning as
    /// solve_with_frequency_parameter_derivative. The original attachment map B
    /// and the entire impedance law Z(omega) are held fixed with respect to theta:
    /// B_theta=0 and Z_theta=Z_omega_theta=0. Their frequency dependence remains.
    /// Complete finite row-major Z_omega and Z_omegaomega are BOTH required;
    /// explicit zero matrices declare the corresponding frequency derivatives zero.
    /// An unloaded owner refuses this entry; use the unloaded entry instead.
    ///
    /// The mixed equation retains the full action
    /// -A_omega*q_theta=(2*omega*M+i*C+i*B*Z*B^T
    ///                   +i*omega*B*Z_omega*B^T)*q_theta.
    /// It reuses the original frequency-action owner, including its normalized
    /// factored port map, full mutual/nonreciprocal entries and uncancelled
    /// envelopes. No dense load-derivative matrix or additional factor is formed.
    /// The existing curvature owner also retains Z_omegaomega in q_omegaomega.
    ///
    /// These are partial derivatives at fixed physical frequency/parameter,
    /// never derivatives of a frozen impedance sample along frequency. This
    /// entry does not represent changes of the attachment, impedance parameters,
    /// port units, moving coordinates or frequency-dependent M/C/K directions.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_loaded_frequency_parameter_derivative(
        &self,
        force: &[C64],
        direction: QuadraticParameterDirection<'_>,
        force_frequency_derivative: Option<&[C64]>,
        force_frequency_second_derivative: Option<&[C64]>,
        force_parameter_derivative: Option<&[C64]>,
        force_frequency_parameter_derivative: Option<&[C64]>,
        impedance_derivative: &[C64],
        impedance_second_derivative: &[C64],
        controls: HarmonicFrequencyParameterControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicFrequencyParameterResponse, HarmonicResponseError> {
        poll(cx)?;
        if self.port_load.is_none() {
            return Err(invalid(
                "loaded mixed frequency/parameter derivatives require an attached port load",
            ));
        }
        self.solve_frequency_parameter_jet(
            force,
            direction,
            force_frequency_derivative,
            force_frequency_second_derivative,
            force_parameter_derivative,
            force_frequency_parameter_derivative,
            Some(impedance_derivative),
            Some(impedance_second_derivative),
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_frequency_parameter_jet(
        &self,
        force: &[C64],
        direction: QuadraticParameterDirection<'_>,
        force_frequency_derivative: Option<&[C64]>,
        force_frequency_second_derivative: Option<&[C64]>,
        force_parameter_derivative: Option<&[C64]>,
        force_frequency_parameter_derivative: Option<&[C64]>,
        impedance_derivative: Option<&[C64]>,
        impedance_second_derivative: Option<&[C64]>,
        controls: HarmonicFrequencyParameterControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicFrequencyParameterResponse, HarmonicResponseError> {
        self.admit_frequency_parameter(controls)?;
        match (
            &self.port_load,
            impedance_derivative,
            impedance_second_derivative,
        ) {
            (None, None, None) => {}
            (Some(load), Some(first), Some(second)) => {
                let extent = load
                    .ports
                    .checked_mul(load.ports)
                    .ok_or_else(|| invalid("mixed impedance derivative extent overflow"))?;
                for values in [first, second] {
                    poll(cx)?;
                    if values.len() != extent || values.iter().any(|&value| !finite(value)) {
                        return Err(invalid(
                            "complete finite first and second impedance frequency derivatives required",
                        ));
                    }
                }
            }
            _ => {
                return Err(invalid(
                    "loaded mixed jets require both impedance frequency derivatives; unloaded jets require none",
                ));
            }
        }
        for values in [
            Some(force),
            force_frequency_derivative,
            force_frequency_second_derivative,
            force_parameter_derivative,
            force_frequency_parameter_derivative,
        ]
        .into_iter()
        .flatten()
        {
            poll(cx)?;
            if values.len() != self.n || values.iter().any(|&value| !finite(value)) {
                return Err(invalid(
                    "one finite source-jet value per physical coordinate required",
                ));
            }
        }
        let parameter_controls = HarmonicParameterControls {
            maximum_additional_payload_bytes: controls.maximum_additional_payload_bytes,
            maximum_work: controls.maximum_work,
        };
        // Validate every physical direction before the first solve/allocation.
        self.admit_parameter(
            direction,
            force_parameter_derivative,
            None,
            0.0,
            parameter_controls,
            cx,
        )?;
        let frequency_response = self.solve_with_frequency_curvature(
            force,
            force_frequency_derivative,
            force_frequency_second_derivative,
            impedance_derivative,
            impedance_second_derivative,
            HarmonicFrequencyCurvatureControls {
                maximum_additional_payload_bytes: controls.maximum_additional_payload_bytes,
                maximum_work: controls.maximum_work,
            },
            cx,
        )?;
        let parameter_response = self.solve_with_parameter_derivative(
            force,
            direction,
            force_parameter_derivative,
            parameter_controls,
            cx,
        )?;
        let frequency_action = if self.port_load.is_some() {
            Some(self.physical_frequency_action(
                &parameter_response.displacement_derivative,
                impedance_derivative,
                cx,
            )?)
        } else {
            None
        };
        let mut rhs_exponent = i32::MIN;
        self.frequency_parameter_terms(
            &frequency_response,
            &parameter_response,
            direction,
            force_frequency_parameter_derivative,
            cx,
            |_, _, term| {
                if term.fraction.hi != 0.0 {
                    rhs_exponent = rhs_exponent.max(term.exponent);
                }
                Ok(())
            },
        )?;
        if let Some((values, envelopes, exponent)) = &frequency_action {
            for i in 0..self.n {
                poll(cx)?;
                for component in [values[i].re, values[i].im, envelopes[i]] {
                    if component != 0.0 {
                        rhs_exponent =
                            rhs_exponent.max(crate::quadratic_fraction(component).1 + exponent);
                    }
                }
            }
        }
        rhs_exponent = if rhs_exponent == i32::MIN {
            0
        } else {
            // Room for every binary-fraction product and its uncancelled
            // row envelope, including both original complex quadratures.
            rhs_exponent + self.n.next_power_of_two().ilog2() as i32 + 8
        };
        let mut sums = vec![[Dd::ZERO; 2]; self.n];
        let mut envelopes = vec![Dd::ZERO; self.n];
        self.frequency_parameter_terms(
            &frequency_response,
            &parameter_response,
            direction,
            force_frequency_parameter_derivative,
            cx,
            |row, axis, term| {
                let value = term.scaled(rhs_exponent)?;
                sums[row][axis] = sums[row][axis] + value;
                envelopes[row] = envelopes[row] + value.abs();
                Ok(())
            },
        )?;
        if let Some((values, bounds, exponent)) = &frequency_action {
            for i in 0..self.n {
                poll(cx)?;
                for (axis, component) in [values[i].re, values[i].im].into_iter().enumerate() {
                    sums[i][axis] = sums[i][axis]
                        + Term::product([component, 1.0, 1.0, 1.0], *exponent)
                            .scaled(rhs_exponent)?;
                }
                envelopes[i] = envelopes[i]
                    + Term::product([bounds[i], 1.0, 1.0, 1.0], *exponent).scaled(rhs_exponent)?;
            }
        }
        let rhs: Vec<_> = sums
            .iter()
            .map(|value| C64::new(value[0].to_f64(), value[1].to_f64()))
            .collect();
        let (mixed, solve_residual) = self.solve_normalized(&rhs, cx)?;
        let mut displacement_mixed_derivative = Vec::with_capacity(self.n);
        let mut velocity_mixed_derivative = Vec::with_capacity(self.n);
        for (i, &value) in mixed.iter().enumerate() {
            poll(cx)?;
            let dq = scale_complex(
                value,
                rhs_exponent - self.coordinate_exponents[i] - 2 * self.time_exponent,
            )?;
            displacement_mixed_derivative.push(dq);
            // Sum before rounding: a changing source phase may cancel a large
            // kinematic contribution. At DC the q_theta term still remains.
            let mut real = PowerSum::default();
            let mut imaginary = PowerSum::default();
            for (rate, motion) in [
                (1.0, parameter_response.displacement_derivative[i]),
                (self.angular_frequency_rad_s, dq),
            ] {
                real.add(Term::product([rate, motion.im, 1.0, 1.0], 0))?;
                imaginary.add(Term::product([-rate, motion.re, 1.0, 1.0], 0))?;
            }
            velocity_mixed_derivative.push(C64::new(real.value()?, imaginary.value()?));
        }
        let residual = self.parameter_residual(
            &displacement_mixed_derivative,
            &rhs,
            &envelopes,
            rhs_exponent,
            cx,
        )?;
        poll(cx)?;
        Ok(HarmonicFrequencyParameterResponse {
            frequency_response,
            parameter_response,
            displacement_mixed_derivative,
            velocity_mixed_derivative,
            relative_mixed_derivative_residual: solve_residual.max(residual),
        })
    }

    fn admit_frequency_parameter(
        &self,
        controls: HarmonicFrequencyParameterControls,
    ) -> Result<(), HarmonicResponseError> {
        self.n
            .checked_mul(self.n)
            .ok_or_else(|| invalid("mixed parameter matrix extent overflow"))?;
        let n = self.n as u128;
        let p = self.port_load.as_ref().map_or(0, |load| load.ports) as u128;
        let k = self
            .compact_port_factor
            .as_ref()
            .map_or(0, |factor| factor.dimension()) as u128;
        let bytes = usize::try_from(8 * n * p + 3072 * n + 1536 * p + 128 * k)
            .map_err(|_| invalid("mixed frequency/parameter payload extent overflow"))?;
        if bytes > controls.maximum_additional_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "mixed frequency/parameter additional payload bytes",
                requested: bytes,
                cap: controls.maximum_additional_payload_bytes,
            });
        }
        let work = 1024 * n * n
            + 2048 * n
            + 512 * (n * p + p * p + k * k)
            + 2048 * (p + k)
            + if k == 0 { 0 } else { 512 * n * n * p };
        // Coupled LU/compact solves reserve the original solve plus at most
        // two same-factor residual corrections; scalar owners stay direct.
        let work = if self.factor.is_some() || self.compact_port_factor.is_some() {
            3 * work
        } else {
            work
        };
        let requested = u64::try_from(work)
            .map_err(|_| invalid("mixed frequency/parameter work extent overflow"))?;
        if requested > controls.maximum_work {
            return Err(HarmonicResponseError::WorkBudget {
                requested,
                cap: controls.maximum_work,
            });
        }
        Ok(())
    }

    // Original physical RHS products, divided by the fixed coordinate scale
    // of their row. Do not differentiate that numerical coordinate/time map.
    #[allow(clippy::too_many_arguments)]
    fn frequency_parameter_terms(
        &self,
        frequency: &HarmonicFrequencyCurvatureResponse,
        parameter: &HarmonicParameterResponse,
        direction: QuadraticParameterDirection<'_>,
        force_mixed: Option<&[C64]>,
        cx: &Cx<'_>,
        mut visit: impl FnMut(usize, usize, Term) -> Result<(), HarmonicResponseError>,
    ) -> Result<(), HarmonicResponseError> {
        let w = self.angular_frequency_rad_s;
        for i in 0..self.n {
            poll(cx)?;
            let shift = -self.coordinate_exponents[i];
            let force = force_mixed.map_or(C64::ZERO, |values| values[i]);
            for (axis, value) in [force.re, force.im].into_iter().enumerate() {
                visit(i, axis, Term::product([value, 1.0, 1.0, 1.0], shift))?;
            }
            // -K_theta*q_omega + M_theta*(w*w*q_omega+2*w*q)
            // -C_theta*v_omega. The last term retains BOTH damping terms
            // i*C_theta*q and i*w*C_theta*q_omega, including at DC.
            for j in 0..self.n {
                let at = i * self.n + j;
                let q = frequency.frequency_response.response.displacement[j];
                let qw = frequency.frequency_response.displacement_derivative[j];
                let vw = frequency.frequency_response.velocity_derivative[j];
                for (axis, (q, qw, vw)) in [(q.re, qw.re, vw.re), (q.im, qw.im, vw.im)]
                    .into_iter()
                    .enumerate()
                {
                    visit(
                        i,
                        axis,
                        Term::product([-direction.stiffness[at], qw, 1.0, 1.0], shift),
                    )?;
                    visit(
                        i,
                        axis,
                        Term::product([direction.mass[at], qw, w, w], shift),
                    )?;
                    visit(
                        i,
                        axis,
                        Term::product([direction.mass[at], q, w, 1.0], shift + 1),
                    )?;
                    visit(
                        i,
                        axis,
                        Term::product([-direction.damping[at], vw, 1.0, 1.0], shift),
                    )?;
                }
            }
            // An unloaded owner retains its original normalized mechanical
            // action. Loaded owners append the COMPLETE shared frequency
            // action once, including Z and Z_omega, outside this traversal.
            if self.port_load.is_none() {
                for j in self.layout.columns(i, self.n) {
                    let at = self.layout.index(i, j, self.n);
                    let q = parameter.displacement_derivative[j];
                    for (axis, value) in [q.re, q.im].into_iter().enumerate() {
                        visit(
                            i,
                            axis,
                            Term::product(
                                [self.scaled_mass[at], value, w, 1.0],
                                self.coordinate_exponents[j] + 1,
                            ),
                        )?;
                    }
                    visit(
                        i,
                        0,
                        Term::product([-self.damping[at], q.im, 1.0, 1.0], shift),
                    )?;
                    visit(
                        i,
                        1,
                        Term::product([self.damping[at], q.re, 1.0, 1.0], shift),
                    )?;
                }
            }
        }
        Ok(())
    }
}
