//! Resolved fields of the original rectangular thermoviscous modes.
//!
//! For short side a, long side b and lambda=m*pi/a (odd m), the existing
//! +i omega diffusion solution is
//! `h = sum 4/(m*pi) sin(lambda*y) s/(s+D*lambda²) f(z)`,
//! `f(z)=1-cosh(q*(z-b/2))/cosh(q*b/2)`, `q²=lambda²+s/D`.
//! This is the separated solution of Stinson (1991), equations 55--69,
//! <https://doi.org/10.1121/1.400379>. Its same retained modes supply the
//! mean wave, both transverse gradients and the positive gradient integrals.
//! Only decaying exponentials are formed, even for very thin wall layers.
//!
//! Peak exp(-i omega t) fields are `u=(U/A) conj(h_v/<h_v>)` and a thermal
//! pressure response `p conj(h_t)`. Shear and thermal-conduction dissipation
//! are `mu |grad u|²/2` and `(Cp/Cv-1) D_t |p grad h_t|²/(2 K_s)`.
//! The latter uses the original fluid's compressibility identity without
//! guessing its expansion coefficient or temperature. The original uniform-
//! pressure, fully developed, rigid/no-slip/isothermal-wall assumptions apply.
//! Entrance effects, transverse velocity, axial-gradient/bulk/molecular losses
//! and wall/fluid heat partition are not supplied by these observations.

use super::{
    FluidDuctError, PI, RectangularAdmission, RectangularSection, RectangularWave, bad,
    finite, mode_response, poll, positive,
};
use crate::DuctError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::gas::GasState;
use fs_math::{c64::C64, dd::Dd, det};

/// Actual transverse fields. Coordinates and gradient components are ordered
/// [short side, long side], independently of constructor argument order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularField {
    /// Distance from the first wall along the shorter side [m].
    pub short_coordinate_m: f64,
    /// Distance from the adjacent wall along the longer side [m].
    pub long_coordinate_m: f64,
    /// Cross-section-uniform pressure, peak exp(-i omega t) phasor [Pa].
    pub pressure_pa: C64,
    /// Actual axial particle velocity, rather than the section mean U/A [m/s].
    pub axial_velocity_m_s: C64,
    /// Transverse derivatives of axial velocity [1/s].
    pub transverse_velocity_gradient_per_s: [C64; 2],
    /// p*h_t [Pa]; the original adiabatic dT/dp converts it to T' [K].
    /// Absent exactly when Cp=Cv, with no guessed liquid expansion law.
    pub thermal_pressure_response_pa: Option<C64>,
    /// Both transverse derivatives of the pressure-normalized thermal field.
    pub thermal_pressure_gradient_pa_m: Option<[C64; 2]>,
    /// Irreversible cycle-mean shear work per fluid volume [W/m³].
    pub viscous_dissipation_w_m3: f64,
    /// Irreversible cycle-mean heat-conduction work per fluid volume [W/m³].
    pub thermal_dissipation_w_m3: f64,
}

/// Independently positive local cross-sectional work, per metre of duct.
/// This is neither a finite-length integral nor a boundary-flux difference.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RectangularDissipationPerLength {
    /// Original transverse shear work [W/m].
    pub viscous_w_m: f64,
    /// Original transverse heat-conduction work [W/m].
    pub thermal_w_m: f64,
    /// Sum of the two independently positive channels [W/m].
    pub total_w_m: f64,
}

#[derive(Debug)]
struct FieldMode {
    phase: f64,
    transverse: f64,
    q: C64,
    amplitude: C64,
    denominator_inverse: C64,
}

#[derive(Debug)]
struct DiffusionField {
    modes: Vec<FieldMode>,
    /// sqrt(integral |grad h|² dA), evaluated from original modal gradients.
    gradient_norm: f64,
}

/// Source-bound fields of the same finite Fourier section as the mean wave.
/// Preparation is O(N); a point query is O(N) without allocation or material
/// evaluation, and integrated power is O(1). At most two N-mode arrays are
/// retained, bounded by the original 4096-term cap and fallible allocation.
/// The existing tail admission bounds section means, not pointwise gradients:
/// refine the authored term count for local fields, especially near corners.
#[derive(Debug)]
pub struct PreparedRectangularFields<'geometry, 'material, M: ThermoviscousFluid + ?Sized> {
    section: &'geometry RectangularSection,
    source: &'material M,
    state: ThermoviscousFluidState,
    omega: f64,
    wave: RectangularWave,
    viscous: DiffusionField,
    thermal: Option<DiffusionField>,
    velocity_scale: C64,
    loss_scales: [f64; 2],
    total_loss_roots: [f64; 2],
}

impl RectangularSection {
    /// Prepare actual no-slip velocity and isothermal thermal fields from ONE
    /// coherent original material snapshot. The native wave is unchanged.
    pub fn prepare_fields<'g, 'm, M: ThermoviscousFluid + ?Sized>(
        &'g self,
        material: &'m M,
        omega: f64,
        admission: RectangularAdmission,
    ) -> Result<PreparedRectangularFields<'g, 'm, M>, FluidDuctError> {
        self.prepare_fields_controlled(material, omega, admission, || false)
    }

    /// As `prepare_fields`, with bounded polls at the original diffusion
    /// boundaries and every 32 modes. Refusal publishes no partial fields.
    pub fn prepare_fields_controlled<'g, 'm, M: ThermoviscousFluid + ?Sized>(
        &'g self,
        material: &'m M,
        omega: f64,
        admission: RectangularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRectangularFields<'g, 'm, M>, FluidDuctError> {
        poll(&mut cancelled, 0)?;
        if !positive(omega)
            || !positive(admission.maximum_relative_tail)
            || admission.maximum_relative_tail >= 1.0
            || !positive(admission.maximum_transverse_phase)
            || admission.maximum_transverse_phase >= PI
        {
            return Err(
                bad("rectangular frequency, tail tolerance or phase envelope is invalid").into(),
            );
        }
        let state = material.thermoviscous_state()?;
        // Evaluate the immutable local snapshot, never re-query the source.
        let wave = self.evaluate_controlled(&state, omega, admission, &mut cancelled)?;
        let acoustic = state.acoustic();
        let transport = state.transport();
        let viscous = self.prepare_diffusion_fields(
            omega / transport.kinematic_viscosity_m2_s(),
            &mut cancelled,
        )?;
        let thermal = if wave.thermal.is_some() {
            Some(self.prepare_diffusion_fields(
                omega / transport.thermal_diffusivity_m2_s(),
                &mut cancelled,
            )?)
        } else {
            None
        };
        let velocity_scale = wave.viscous.mean.conj().recip() / C64::from_re(self.area_m2);
        let loss_scales = [
            det::sqrt(transport.dynamic_viscosity_pa_s()) * core::f64::consts::FRAC_1_SQRT_2,
            (det::sqrt(state.heat_capacity_ratio() - 1.0)
                / det::sqrt(acoustic.adiabatic_bulk_modulus_pa()))
                * det::sqrt(transport.thermal_diffusivity_m2_s())
                * core::f64::consts::FRAC_1_SQRT_2,
        ];
        let total_loss_roots = [
            (velocity_scale.scale(loss_scales[0]) * C64::from_re(viscous.gradient_norm)).abs(),
            thermal
                .as_ref()
                .map_or(0.0, |t| loss_scales[1] * t.gradient_norm),
        ];
        if !finite(velocity_scale)
            || velocity_scale == C64::ZERO
            || !positive(loss_scales[0])
            || !positive(total_loss_roots[0])
            || !loss_scales[1].is_finite()
            || !total_loss_roots[1].is_finite()
            || (thermal.is_some() && (loss_scales[1] <= 0.0 || total_loss_roots[1] <= 0.0))
        {
            return Err(bad("rectangular field normalization or work is unrepresentable").into());
        }
        poll(&mut cancelled, self.terms())?;
        Ok(PreparedRectangularFields {
            section: self,
            source: material,
            state,
            omega,
            wave,
            viscous,
            thermal,
            velocity_scale,
            loss_scales,
            total_loss_roots,
        })
    }

    fn prepare_diffusion_fields(
        &self,
        frequency: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DiffusionField, DuctError> {
        let mut modes = Vec::new();
        modes
            .try_reserve_exact(self.terms())
            .map_err(|_| bad("rectangular field allocation exceeds available memory"))?;
        let mut energy = Dd::ZERO;
        for (index, mode) in self.modes.iter().enumerate() {
            if index % 32 == 0 {
                poll(cancelled, index)?;
            }
            let response = mode_response(mode, self.long_m, C64::new(0.0, frequency))?;
            let phase = PI * (2 * index + 1) as f64;
            let amplitude = response.factor.scale(4.0 / phase);
            let q = response.q;
            // Integral |cosh(q(z-b/2))/cosh(qb/2)|² and its derivative.
            // Scale hyperbolic functions by exp(-Re(q)b); Re(q)b >= pi.
            let a = q.re * self.long_m;
            let b = q.im * self.long_m;
            let decay = det::exp(-a);
            let (sine, cosine) = det::sin_cos(b);
            let denominator = 1.0 + decay * decay + 2.0 * decay * cosine;
            let hyperbolic = (1.0 - decay * decay) / (denominator * q.re);
            let sinc = if b == 0.0 { 1.0 } else { sine / b };
            let circular = 2.0 * decay * self.long_m * sinc / denominator;
            let field_integral =
                self.long_m - 4.0 * (response.tanh / q).re + hyperbolic + circular;
            let derivative_integral = q.norm_sq() * (hyperbolic - circular);
            // Orthogonality removes every cross-mode term in both transverse
            // gradient directions. Each retained modal contribution is positive.
            let norm = det::sqrt(0.5 * self.short_m)
                * amplitude.abs()
                * det::sqrt(mode.eigenvalue * field_integral + derivative_integral);
            if !positive(field_integral) || !positive(derivative_integral) || !positive(norm) {
                return Err(bad("rectangular modal gradient integral is unrepresentable"));
            }
            energy = energy + Dd::from_f64(norm * norm);
            modes.push(FieldMode {
                phase,
                transverse: phase / self.short_m,
                q,
                amplitude,
                denominator_inverse: (C64::ONE + response.decay).recip(),
            });
        }
        let gradient_norm = det::sqrt(energy.to_f64());
        if !positive(gradient_norm) {
            return Err(bad("rectangular total gradient integral is unrepresentable"));
        }
        Ok(DiffusionField {
            modes,
            gradient_norm,
        })
    }
}

impl<'g, 'm, M: ThermoviscousFluid + ?Sized> PreparedRectangularFields<'g, 'm, M> {
    /// Retained original material owner, including its own source receipts.
    pub const fn source(&self) -> &'m M {
        self.source
    }
    /// Original geometry and declared transverse resolution.
    pub const fn section(&self) -> &'g RectangularSection {
        self.section
    }
    /// The single coherent snapshot used by the wave and every observation.
    pub const fn state(&self) -> &ThermoviscousFluidState {
        &self.state
    }
    /// Angular frequency binding the local pressure/flow amplitudes [rad/s].
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// The same section-mean response as `RectangularSection::evaluate`.
    pub const fn wave(&self) -> &RectangularWave {
        &self.wave
    }

    /// Physical coordinates measured from the two adjacent walls, along the
    /// short and long sides. All four no-slip/isothermal walls are retained.
    /// Supply coherent section pressure and volume flow at this frequency.
    pub fn field_at(
        &self,
        short_coordinate_m: f64,
        long_coordinate_m: f64,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<RectangularField, DuctError> {
        self.field_at_controlled(
            short_coordinate_m,
            long_coordinate_m,
            pressure_pa,
            volume_flow_m3_s,
            || false,
        )
    }

    /// As `field_at`, with cancellation polled every 32 retained modes.
    pub fn field_at_controlled(
        &self,
        short_coordinate_m: f64,
        long_coordinate_m: f64,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RectangularField, DuctError> {
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let a = self.section.short_m;
        let b = self.section.long_m;
        if !short_coordinate_m.is_finite()
            || !(0.0..=a).contains(&short_coordinate_m)
            || !long_coordinate_m.is_finite()
            || !(0.0..=b).contains(&long_coordinate_m)
        {
            return Err(bad("rectangular field point is outside the original section"));
        }
        let v = self
            .viscous
            .at(short_coordinate_m / a, long_coordinate_m, b, &mut cancelled)?
            .map(|value| value.conj() * self.velocity_scale * volume_flow_m3_s);
        let thermal = self
            .thermal
            .as_ref()
            .map(|t| {
                t.at(short_coordinate_m / a, long_coordinate_m, b, &mut cancelled)
                    .map(|values| values.map(|value| value.conj() * pressure_pa))
            })
            .transpose()?;
        let viscous_work: f64 = v[1..]
            .iter()
            .map(|g| g.scale(self.loss_scales[0]).norm_sq())
            .sum();
        let thermal_work = thermal.map_or(0.0, |t| {
            t[1..]
                .iter()
                .map(|g| g.scale(self.loss_scales[1]).norm_sq())
                .sum()
        });
        if !v.into_iter().all(finite)
            || thermal.is_some_and(|t| !t.into_iter().all(finite))
            || ![viscous_work, thermal_work].into_iter().all(f64::is_finite)
        {
            return Err(bad("rectangular physical field or local work is unrepresentable"));
        }
        Ok(RectangularField {
            short_coordinate_m,
            long_coordinate_m,
            pressure_pa,
            axial_velocity_m_s: v[0],
            transverse_velocity_gradient_per_s: [v[1], v[2]],
            thermal_pressure_response_pa: thermal.map(|t| t[0]),
            thermal_pressure_gradient_pa_m: thermal.map(|t| [t[1], t[2]]),
            viscous_dissipation_w_m3: viscous_work,
            thermal_dissipation_w_m3: thermal_work,
        })
    }

    /// Whole-section positive physical work, from the analytic integrals of
    /// the same original gradients. No wave-power subtraction supplies it.
    pub fn dissipation_per_axial_length(
        &self,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<RectangularDissipationPerLength, DuctError> {
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let viscous = volume_flow_m3_s.scale(self.total_loss_roots[0]).norm_sq();
        let thermal = pressure_pa.scale(self.total_loss_roots[1]).norm_sq();
        let total = viscous + thermal;
        if ![viscous, thermal, total].into_iter().all(f64::is_finite) {
            return Err(bad("rectangular integrated physical work is unrepresentable"));
        }
        Ok(RectangularDissipationPerLength {
            viscous_w_m: viscous,
            thermal_w_m: thermal,
            total_w_m: total,
        })
    }
}

impl PreparedRectangularFields<'_, '_, GasState> {
    /// Excess temperature [K] using this original ideal gas's 1/(rho Cp).
    /// A general fluid's thermal expansion or reference temperature is not guessed.
    pub fn gas_temperature_at(
        &self,
        short_coordinate_m: f64,
        long_coordinate_m: f64,
        pressure_pa: C64,
    ) -> Result<C64, DuctError> {
        let field = self.field_at(short_coordinate_m, long_coordinate_m, pressure_pa, C64::ZERO)?;
        let temperature = field.thermal_pressure_response_pa.unwrap_or(C64::ZERO)
            / C64::from_re(self.state.acoustic().density_kg_m3())
            / C64::from_re(self.state.transport().specific_heat_j_kg_k());
        if !finite(temperature) {
            return Err(bad("rectangular gas temperature is unrepresentable"));
        }
        Ok(temperature)
    }
}

impl DiffusionField {
    fn at(
        &self,
        fraction: f64,
        z: f64,
        long_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[C64; 3], DuctError> {
        let (mut re, mut im) = ([Dd::ZERO; 3], [Dd::ZERO; 3]);
        for (index, mode) in self.modes.iter().enumerate() {
            if index % 32 == 0 {
                poll(cancelled, index)?;
            }
            let (sine, cosine) = if fraction == 0.0 {
                (0.0, 1.0)
            } else if fraction == 1.0 {
                (0.0, -1.0)
            } else if fraction == 0.5 {
                (det::sin_cos(0.5 * mode.phase).0, 0.0)
            } else if fraction < 0.5 {
                det::sin_cos(mode.phase * fraction)
            } else {
                // Odd sine modes are symmetric about the mid-plane; their
                // short-direction derivatives are antisymmetric.
                let (s, c) = det::sin_cos(mode.phase * (1.0 - fraction));
                (s, -c)
            };
            let (left, one_minus_left) = decay_pair(mode.q.scale(z));
            let (right, one_minus_right) = decay_pair(mode.q.scale(long_m - z));
            let shape = one_minus_left * one_minus_right * mode.denominator_inverse;
            let derivative = mode.q * (left - right) * mode.denominator_inverse;
            let values = [
                shape.scale(sine),
                shape.scale(mode.transverse * cosine),
                derivative.scale(sine),
            ];
            for (i, value) in values.map(|v| v * mode.amplitude).into_iter().enumerate() {
                re[i] = re[i] + Dd::from_f64(value.re);
                im[i] = im[i] + Dd::from_f64(value.im);
            }
        }
        Ok(core::array::from_fn(|i| C64::new(re[i].to_f64(), im[i].to_f64())))
    }
}

// exp(-z) and 1-exp(-z), retaining tiny wall distances and exact wall zeros.
fn decay_pair(z: C64) -> (C64, C64) {
    let decay = det::exp(-z.re);
    if decay == 0.0 {
        return (C64::ZERO, C64::ONE);
    }
    let (sine, cosine) = det::sin_cos(z.im);
    let half_sine = det::sin_cos(0.5 * z.im).0;
    (
        C64::new(cosine, -sine).scale(decay),
        C64::new(
            -det::expm1(-z.re) + 2.0 * decay * half_sine * half_sine,
            decay * sine,
        ),
    )
}

fn amplitudes(pressure: C64, flow: C64) -> Result<(), DuctError> {
    if !finite(pressure) || !finite(flow) {
        Err(bad("rectangular pressure and volume flow must be finite"))
    } else {
        Ok(())
    }
}
