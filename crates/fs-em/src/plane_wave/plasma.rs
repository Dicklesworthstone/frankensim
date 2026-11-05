//! Collisional plasma in the original plane-wave half-space owner.
//!
//! With exp(+i omega t), Maxwell-Ampere uses
//! epsilon_eff = epsilon_background - i sigma(omega)/omega. The free-carrier
//! conductivity is evaluated by the original prepared material owner; there
//! is no refitted Drude pole or DC conductivity substitution. Negative real
//! epsilon, passive attenuation and the outgoing evanescent branch remain
//! those of the existing interface kernel.
//!
//! Reference: R. Fitzpatrick, equations 1204, 1205 and 1208 (opposite phasor):
//! <https://farside.ph.utexas.edu/teaching/em/lectures/node103.html>.
//! A single outgoing wave decays away from the interface. Its surface field
//! is therefore the largest peak drift anywhere in the homogeneous half-space.
//! Validate that actual three-dimensional field once; depth queries then reuse
//! its original carrier currents, friction partitions and kinetic storage.
//! Frozen populations/background, unmagnetized weak drift and local Markovian
//! collisions are inherited. No sheath, changing ionization, spatial kinetic
//! response, transient turn-on or measured broadband accuracy is supplied.

use super::{
    DielectricResponse, InterfaceReflection, PlanarDielectricInterface, Polarization,
    ReflectionError, finite, positive,
};
use fs_material::{
    field_storage::{linear_electric_permittivity, linear_magnetic_permeability},
    ionized_gas::{IonizedGasError, IonizedGasHarmonic, IonizedGasHarmonicResponse},
    state_point::ResolvedMaterialStatePoint,
};
use fs_math::{c64::C64, dd::Dd, det};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};

/// One frequency of the original collisional transport plus an explicit
/// lossless bound-charge background. Both complete sources remain borrowed.
/// Preparation evaluates the material once at zero field; normalized interface
/// queries subsequently allocate nothing and cannot silently change frequency.
#[derive(Clone, Copy, Debug)]
pub struct IonizedGasWaveResponse<'a> {
    transport: &'a IonizedGasHarmonic<'a>,
    background: &'a ResolvedMaterialStatePoint,
    omega_rad_s: f64,
    conductivity_s_m: C64,
    permittivity_f_m: C64,
    background_permittivity_f_m: f64,
    permeability_h_m: f64,
}

impl<'a> IonizedGasWaveResponse<'a> {
    /// The caller explicitly declares that background epsilon excludes these
    /// free carriers and applies to this frozen mixture. No vacuum constant,
    /// species polarizability or chemistry-name mapping is inferred.
    ///
    /// Background T/T_h, T_e and typed pressure p coordinates, when supplied,
    /// must match the original heavy/electron temperatures and pressure.
    /// Other condition coordinates refuse because no mapping to this gas
    /// state is supplied. An unconditioned background remains unconditioned.
    pub fn try_new(
        transport: &'a IonizedGasHarmonic<'a>,
        background: &'a ResolvedMaterialStatePoint,
        omega_rad_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReflectionError> {
        poll(&mut cancelled)?;
        let omega = positive(omega_rad_s)?;
        let epsilon_background = linear_electric_permittivity(background)?;
        let permeability = linear_magnetic_permeability(background)?;
        admit_background(transport, background)?;
        let response = transport.response(omega, [C64::ZERO; 3], &mut cancelled)?;
        let sigma = response.conductivity_s_m;
        // Retain the small real part around the plasma cutoff after division.
        // This improves arithmetic; it is not a material-error enclosure.
        let epsilon = C64::new(
            (Dd::from_f64(epsilon_background) + Dd::from_f64(sigma.im) / Dd::from_f64(omega))
                .to_f64(),
            -positive(sigma.re / omega)?,
        );
        if !finite(epsilon) {
            return Err(bad("plasma effective permittivity is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            transport,
            background,
            omega_rad_s: omega,
            conductivity_s_m: sigma,
            permittivity_f_m: epsilon,
            background_permittivity_f_m: epsilon_background,
            permeability_h_m: permeability,
        })
    }

    pub const fn transport(self) -> &'a IonizedGasHarmonic<'a> {
        self.transport
    }
    /// This is the bound-charge background only. The separate transport()
    /// retains the actual gas populations and all original collision sources.
    pub const fn material(self) -> &'a ResolvedMaterialStatePoint {
        self.background
    }
    pub const fn omega_rad_s(self) -> f64 {
        self.omega_rad_s
    }
    pub const fn conductivity_s_m(self) -> C64 {
        self.conductivity_s_m
    }
    pub const fn effective_permittivity_f_m(self) -> C64 {
        self.permittivity_f_m
    }
    pub(super) fn permittivity(self, omega: f64) -> Result<C64, ReflectionError> {
        if omega != self.omega_rad_s {
            return Err(bad(
                "plasma wave source was prepared at a different angular frequency",
            ));
        }
        Ok(self.permittivity_f_m)
    }
}

/// One admitted physical outgoing wave, retaining both electromagnetic sources
/// and the complete original material response at its maximum surface field.
/// Preparation costs one bounded material response. Any depth query is O(1),
/// allocation-free, and never repeats the carrier eigenspace calculation.
pub struct PreparedIonizedGasWave<'a> {
    incident: DielectricResponse<'a>,
    plasma: IonizedGasWaveResponse<'a>,
    reflection: InterfaceReflection,
    incident_tangential_e_v_m: C64,
    incidence_cosine: f64,
    polarization: Polarization,
    surface: IonizedGasHarmonicResponse<'a>,
    magnetic_a_m: [C64; 3],
    poynting_w_m2: [C64; 3],
    background_energy_j_m3: f64,
    power_balance_relative_residual: f64,
}

/// Actual peak vector fields at x=y=0, z=depth_m in the original half-space.
/// Carrier rows and signed thermal partitions are read from the retained
/// surface response and scaled by the same passive propagation factor.
/// Wavelength and attenuation scales of the admitted normal propagation.
/// These come from the same complex normal wave number used by `field`; no
/// second material solve, vacuum-index substitution or fitted skin depth is
/// introduced. `None` means that the corresponding finite positive scale is
/// absent (for example a purely evanescent phase or a lossless amplitude).
/// Group velocity is deliberately not inferred from one monochromatic sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonizedGasWaveScales {
    /// 2*pi/Re(k_z) [m], when the transmitted normal phase propagates.
    pub normal_phase_wavelength_m: Option<f64>,
    /// omega/Re(k_z) [m/s], a phase velocity, never a group velocity.
    pub normal_phase_velocity_m_s: Option<f64>,
    /// 1/(-Im(k_z)) [m], where field amplitude falls by e^-1.
    pub amplitude_skin_depth_m: Option<f64>,
    /// 1/(-2 Im(k_z)) [m], where intensity/mean power falls by e^-1.
    pub intensity_skin_depth_m: Option<f64>,
    /// Re(k_z)/(2*pi*(-Im(k_z))): phase cycles per amplitude e-fold.
    pub phase_cycles_per_amplitude_efold: Option<f64>,
}

pub struct IonizedGasWavePoint<'w, 'a> {
    wave: &'w PreparedIonizedGasWave<'a>,
    depth_m: f64,
    propagation: C64,
    attenuation: f64,
    pub electric_v_m: [C64; 3],
    pub magnetic_a_m: [C64; 3],
    pub current_density_a_m2: [C64; 3],
    /// E cross conj(H)/2 [W/m²], including reactive flux.
    pub complex_poynting_w_m2: [C64; 3],
    /// Physical cycle-mean heat integrated from the interface to this depth
    /// [W/m²]. This excludes flux still travelling into deeper material.
    pub absorbed_power_to_depth_w_m2: f64,
}

impl<'a> PlanarDielectricInterface<'a> {
    /// Prepare a physical wave transmitted into an original ionized gas.
    /// The actual peak tangential incident E is in V/m, including source phase.
    /// Its complete transmitted surface vector (including longitudinal TM E)
    /// must pass the original weak-drift and Coulomb-relative-drift limits.
    ///
    /// A single outgoing homogeneous wave has no interior enhancement; every
    /// deeper point is a common complex phase times this admitted surface.
    /// Finite films/standing waves require a separate whole-field admission.
    pub fn prepare_ionized_gas_wave(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        incident_tangential_e_v_m: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedIonizedGasWave<'a>, ReflectionError> {
        poll(&mut cancelled)?;
        let DielectricResponse::IonizedGas(plasma) = self.transmitted else {
            return Err(bad(
                "transmitted wave does not retain an ionized-gas material",
            ));
        };
        if !finite(incident_tangential_e_v_m) {
            return Err(bad("plasma incident peak field must be finite"));
        }
        let omega = plasma.omega_rad_s;
        let reflection = self.response(omega, incidence_cosine, polarization)?;
        let epsilon_i = self.incident_response.lossless_permittivity(omega)?;
        let root_epsilon = det::sqrt(epsilon_i);
        let root_mu = det::sqrt(self.incident_mu_h_m);
        let base_admittance = positive(root_epsilon / root_mu)?;
        let admittance = positive(match polarization {
            Polarization::Te => base_admittance * incidence_cosine,
            Polarization::Tm => base_admittance / incidence_cosine,
        })?;
        let sine = det::sqrt((1.0 - incidence_cosine) * (1.0 + incidence_cosine));
        let parallel_over_omega = root_epsilon * (root_mu * sine);
        let electric = reflection.transmitted_tangential_e * incident_tangential_e_v_m;
        // Retain the independent native H numerator. Forming 1-r would erase
        // a small physical H after a reflection rounded to one.
        let magnetic =
            (reflection.transmitted_conjugate_h * incident_tangential_e_v_m).scale(admittance);
        let (e, h) = match polarization {
            Polarization::Te => (
                [C64::ZERO, electric, C64::ZERO],
                [
                    -magnetic,
                    C64::ZERO,
                    electric.scale(parallel_over_omega / plasma.permeability_h_m),
                ],
            ),
            Polarization::Tm => (
                [
                    electric,
                    C64::ZERO,
                    (C64::from_re(-parallel_over_omega) / plasma.permittivity_f_m) * magnetic,
                ],
                [C64::ZERO, magnetic, C64::ZERO],
            ),
        };
        if !parallel_over_omega.is_finite() || e.iter().chain(&h).any(|&z| !finite(z)) {
            return Err(bad("transmitted plasma vector field is unrepresentable"));
        }
        let surface = plasma.transport.response(omega, e, &mut cancelled)?;
        let poynting = complex_poynting(e, h);
        let incident_power = incident_tangential_e_v_m
            .scale(det::sqrt(admittance) * core::f64::consts::FRAC_1_SQRT_2)
            .norm_sq();
        let transmitted_power = incident_power * reflection.transmitted_power_fraction;
        let flux_loss = (-reflection.transmitted_normal_wavenumber_rad_m.im * poynting[2].re) * 2.0;
        let residual = balance(poynting[2].re, transmitted_power)?.max(balance(
            flux_loss,
            surface.mean_heating_w_m3.total_power_density_w_m3(),
        )?);
        let electric_scale = 0.5 * det::sqrt(plasma.background_permittivity_f_m);
        let magnetic_scale = 0.5 * det::sqrt(plasma.permeability_h_m);
        let background_energy = e
            .iter()
            .map(|z| z.scale(electric_scale).norm_sq())
            .chain(h.iter().map(|z| z.scale(magnetic_scale).norm_sq()))
            .fold(Dd::ZERO, |s, x| s + Dd::from_f64(x))
            .to_f64();
        if poynting.iter().any(|&z| !finite(z))
            || !background_energy.is_finite()
            || background_energy < 0.0
        {
            return Err(bad("plasma wave flux or field storage is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(PreparedIonizedGasWave {
            incident: self.incident_response,
            plasma,
            reflection,
            incident_tangential_e_v_m,
            incidence_cosine,
            polarization,
            surface,
            magnetic_a_m: h,
            poynting_w_m2: poynting,
            background_energy_j_m3: background_energy,
            power_balance_relative_residual: residual,
        })
    }
}

impl<'a> PreparedIonizedGasWave<'a> {
    pub const fn incident_response(&self) -> DielectricResponse<'a> {
        self.incident
    }
    pub const fn plasma_response(&self) -> IonizedGasWaveResponse<'a> {
        self.plasma
    }
    pub const fn interface_response(&self) -> &InterfaceReflection {
        &self.reflection
    }
    pub const fn surface_response(&self) -> &IonizedGasHarmonicResponse<'a> {
        &self.surface
    }
    pub const fn incident_tangential_e_v_m(&self) -> C64 {
        self.incident_tangential_e_v_m
    }
    pub const fn incidence_cosine(&self) -> f64 {
        self.incidence_cosine
    }
    pub const fn polarization(&self) -> Polarization {
        self.polarization
    }
    /// Surface normal-flux admission and local Maxwell loss versus original
    /// carrier friction. This is a residual, not a continuum error bound.
    pub const fn power_balance_relative_residual(&self) -> f64 {
        self.power_balance_relative_residual
    }

    /// Physical normal wavelength, phase speed and attenuation lengths of this
    /// already-admitted wave. Every value is derived from the retained complex
    /// `k_z`; the original carrier conductivity is not reevaluated.
    #[must_use]
    pub fn propagation_scales(&self) -> IonizedGasWaveScales {
        let kz = self.reflection.transmitted_normal_wavenumber_rad_m;
        let phase_wavelength = positive_scale(core::f64::consts::TAU / kz.re);
        let phase_velocity = positive_scale(self.plasma.omega_rad_s / kz.re);
        let amplitude_depth = positive_scale(-1.0 / kz.im);
        let intensity_depth = positive_scale(-0.5 / kz.im);
        let cycles = positive_scale(kz.re / (-core::f64::consts::TAU * kz.im));
        IonizedGasWaveScales {
            normal_phase_wavelength_m: phase_wavelength,
            normal_phase_velocity_m_s: phase_velocity,
            amplitude_skin_depth_m: amplitude_depth,
            intensity_skin_depth_m: intensity_depth,
            phase_cycles_per_amplitude_efold: cycles,
        }
    }

    /// Evaluate an admitted outgoing wave at finite depth >=0 [m].
    /// Fully attenuated tails need no otherwise unrepresentable phase.
    /// Thermal terms are scaled as (surface_power*a)*a rather than multiplying
    /// by a prematurely underflowed attenuation squared.
    pub fn field(&self, depth_m: f64) -> Result<IonizedGasWavePoint<'_, 'a>, ReflectionError> {
        if !depth_m.is_finite() || depth_m < 0.0 {
            return Err(bad("plasma depth must be finite and nonnegative"));
        }
        let kz = self.reflection.transmitted_normal_wavenumber_rad_m;
        let decay = kz.im * depth_m;
        let attenuation = det::exp(decay);
        let propagation = if attenuation == 0.0 {
            C64::ZERO
        } else {
            let phase = Dd::from_f64(kz.re) * Dd::from_f64(depth_m);
            if !phase.to_f64().is_finite() {
                return Err(bad("plasma propagation phase is unrepresentable"));
            }
            let (sine, cosine) = det::sin_cos_dd(phase);
            C64::new(attenuation * cosine.to_f64(), -attenuation * sine.to_f64())
        };
        if !finite(propagation) || !attenuation.is_finite() || attenuation > 1.0 {
            return Err(bad("plasma outgoing attenuation is unrepresentable"));
        }
        let electric = self.surface.electric_field_v_m.map(|e| e * propagation);
        let magnetic = self.magnetic_a_m.map(|h| h * propagation);
        let current = self.surface.current_density_a_m2.map(|j| j * propagation);
        let poynting = self.poynting_w_m2.map(|s| {
            C64::new(
                scaled_power(s.re, attenuation),
                scaled_power(s.im, attenuation),
            )
        });
        let absorbed = self.poynting_w_m2[2].re * -det::expm1(2.0 * decay);
        if electric
            .iter()
            .chain(&magnetic)
            .chain(&current)
            .chain(&poynting)
            .any(|&z| !finite(z))
            || !absorbed.is_finite()
            || absorbed < 0.0
        {
            return Err(bad("plasma field or accumulated heat is unrepresentable"));
        }
        Ok(IonizedGasWavePoint {
            wave: self,
            depth_m,
            propagation,
            attenuation,
            electric_v_m: electric,
            magnetic_a_m: magnetic,
            current_density_a_m2: current,
            complex_poynting_w_m2: poynting,
            absorbed_power_to_depth_w_m2: absorbed,
        })
    }
}

impl<'w, 'a> IonizedGasWavePoint<'w, 'a> {
    pub const fn wave(&self) -> &'w PreparedIonizedGasWave<'a> {
        self.wave
    }
    pub const fn depth_m(&self) -> f64 {
        self.depth_m
    }
    pub fn carrier_velocity_m_s(&self, carrier: usize) -> Option<[C64; 3]> {
        self.wave
            .surface
            .carrier_velocity_m_s
            .get(carrier)
            .map(|row| row.map(|v| v * self.propagation))
    }
    /// Total physical friction heat [W/m³], deposited only once by a consumer.
    pub fn mean_heating_w_m3(&self) -> f64 {
        scaled_power(
            self.wave
                .surface
                .mean_heating_w_m3
                .total_power_density_w_m3(),
            self.attenuation,
        )
    }
    pub fn mean_electron_heating_w_m3(&self) -> f64 {
        scaled_power(
            self.wave
                .surface
                .mean_heating_w_m3
                .electron_power_density_w_m3(),
            self.attenuation,
        )
    }
    pub fn mean_heavy_heating_w_m3(&self) -> f64 {
        scaled_power(
            self.wave
                .surface
                .mean_heating_w_m3
                .heavy_power_density_w_m3(),
            self.attenuation,
        )
    }
    /// Signed drift correction, already included in the two pool powers.
    /// Baseline zero-drift thermal exchange remains a separate process.
    pub fn mean_exchange_correction_to_electrons_w_m3(&self) -> f64 {
        scaled_power(
            self.wave
                .surface
                .mean_heating_w_m3
                .thermal_exchange_correction_to_electrons_w_m3(),
            self.attenuation,
        )
    }
    pub fn mean_drift_kinetic_energy_j_m3(&self) -> f64 {
        scaled_power(
            self.wave.surface.mean_drift_kinetic_energy_j_m3,
            self.attenuation,
        )
    }
    /// Lossless background electric plus magnetic energy [J/m³]. Carrier
    /// kinetic storage is separate; complex epsilon is not an energy density.
    pub fn mean_background_field_energy_j_m3(&self) -> f64 {
        scaled_power(self.wave.background_energy_j_m3, self.attenuation)
    }
    pub fn reactive_carrier_power_var_m3(&self) -> f64 {
        scaled_power(self.wave.surface.reactive_power_var_m3, self.attenuation)
    }
}

fn positive_scale(value: f64) -> Option<f64> {
    (value.is_finite() && value > 0.0).then_some(value)
}

fn admit_background(
    transport: &IonizedGasHarmonic<'_>,
    background: &ResolvedMaterialStatePoint,
) -> Result<(), ReflectionError> {
    let source = transport.source();
    let state = source.state();
    let pressure = source
        .equilibrium()
        .map_or_else(|| state.pressure().value(), |s| s.pressure().value());
    for (axis, value) in background.query_point() {
        let value = *value;
        let (expected, kind) = match axis.as_str() {
            "T" | "T_h" => (
                state.heavy_temperature().value(),
                QuantityKind::AbsoluteTemperature,
            ),
            "T_e" => (
                state.electron_temperature().value(),
                QuantityKind::AbsoluteTemperature,
            ),
            "p" => (pressure, QuantityKind::Pressure),
            _ => {
                return Err(bad(
                    "plasma background has an unmapped condition coordinate",
                ));
            }
        };
        let quantity = background.axis_quantities().get(axis).copied();
        let typed = quantity
            .and_then(QuantitySpec::semantic_type)
            .is_some_and(|ty| ty.kind() == kind && ty.form() == ValueForm::Static);
        // Preserve the material owner's existing legacy absolute-T convention.
        let legacy_t = axis == "T"
            && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if !(typed || legacy_t)
            || !value.is_finite()
            || value <= 0.0
            || (value != expected
                && (axis != "p" || (value / expected - 1.0).abs() > 512.0 * f64::EPSILON))
        {
            return Err(bad(
                "plasma background and carrier source conditions disagree",
            ));
        }
    }
    Ok(())
}

fn complex_poynting(e: [C64; 3], h: [C64; 3]) -> [C64; 3] {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        (e[j] * h[k].conj() - e[k] * h[j].conj()).scale(0.5)
    })
}
fn scaled_power(value: f64, attenuation: f64) -> f64 {
    (value * attenuation) * attenuation
}
fn balance(a: f64, b: f64) -> Result<f64, ReflectionError> {
    let scale = a.max(b);
    let residual = if scale == 0.0 {
        0.0
    } else {
        (a - b).abs() / scale
    };
    if !a.is_finite()
        || !b.is_finite()
        || a < 0.0
        || b < 0.0
        || !residual.is_finite()
        || residual > 1e-10
    {
        return Err(bad(
            "plasma wave Poynting flux and original friction disagree",
        ));
    }
    Ok(residual)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReflectionError> {
    if cancelled() {
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    } else {
        Ok(())
    }
}
fn bad(message: &'static str) -> ReflectionError {
    ReflectionError::Invalid(message)
}
