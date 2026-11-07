//! Source-bound circular waves at a magnetized plasma half-space.
//!
//! Normal incidence is along or opposite the original static B. For any
//! right-handed transverse frame, e_plus/minus = (e1 +/- i e2)/sqrt(2).
//! With exp(+i omega t) and forward exp(-i k depth), each circular eigenmode
//! uses its ORIGINAL sigma, epsilon = epsilon_background - i sigma/omega,
//! and the existing passive Fresnel kernel. Different complex k values give
//! Faraday rotation and circular dichroism; no scalar Hall surrogate is made.
//!
//! The material is local, classical, weak-drift and linear about its frozen
//! state. The wave magnetic force on the perturbed velocity is second order
//! and is excluded. No oblique propagation, longitudinal plasma mode, finite
//! film, sheath, kinetic spatial dispersion or time-domain pulse is supplied.
//!
//! Circular eigenmodes and the rotation mechanism: Fitzpatrick, equations
//! 5.88--5.91 and 1175--1182 (the references use the opposite time convention):
//! <https://farside.ph.utexas.edu/teaching/plasma/lectures/node74.html>
//! <https://farside.ph.utexas.edu/teaching/em/lectures/node101.html>

use super::super::{Polarization, WaveMedium, interface_response};
use super::{
    DielectricResponse, InterfaceReflection, ReflectionError, admit_background, bad, balance,
    complex_poynting, finite, linear_electric_permittivity, linear_magnetic_permeability, poll,
    positive, scaled_power,
};
use fs_material::{
    ionized_gas::{IonizedGasMagnetizedHarmonic, IonizedGasMagnetizedHarmonicResponse},
    state_point::ResolvedMaterialStatePoint,
};
use fs_math::{c64::C64, dd::Dd, det};

type Vector = [C64; 3];
type Result<T> = core::result::Result<T, ReflectionError>;

mod absorption;
pub use absorption::MagnetizedIonizedGasSlabAbsorption;

/// Original fixed-frequency transport, explicit bound-charge background and
/// lossless incident medium. Circular labels are relative to the outgoing
/// interface normal, so reversing B at fixed normal swaps the two responses.
#[derive(Clone, Copy, Debug)]
pub struct MagnetizedIonizedGasInterface<'a> {
    incident: DielectricResponse<'a>,
    transport: &'a IonizedGasMagnetizedHarmonic<'a>,
    background: &'a ResolvedMaterialStatePoint,
    normal: [f64; 3],
    incident_admittance_s: f64,
    background_permittivity_f_m: f64,
    permeability_h_m: f64,
    conductivity_s_m: [C64; 2],
    modes: [InterfaceReflection; 2],
}

impl<'a> MagnetizedIonizedGasInterface<'a> {
    /// Bind one homogeneous half-space at the transport's prepared frequency.
    /// Background epsilon excludes the retained free carriers; the caller
    /// declares its applicability to this mixture. Its supplied T/T_h, T_e
    /// and typed pressure conditions must match the original gas owner.
    ///
    /// The nonzero finite outgoing normal is normalized. Nonzero B must be
    /// parallel or antiparallel to it, to floating-point angular roundoff.
    /// With zero B any normal is admitted and both modes recover the original
    /// unmagnetized response. The incident source must be positive/lossless.
    pub fn try_new(
        incident: DielectricResponse<'a>,
        transport: &'a IonizedGasMagnetizedHarmonic<'a>,
        background: &'a ResolvedMaterialStatePoint,
        outgoing_normal: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let normal = unit(outgoing_normal)?;
        let epsilon = linear_electric_permittivity(background)?;
        let mu = linear_magnetic_permeability(background)?;
        admit_background(transport.source(), background)?;
        let omega = transport.omega_rad_s();
        let incident_epsilon = incident.lossless_permittivity(omega)?;
        let incident_mu = linear_magnetic_permeability(incident.material())?;
        let admittance = positive(det::sqrt(incident_epsilon) / det::sqrt(incident_mu))?;
        let wave_number = positive(omega * det::sqrt(incident_epsilon) * det::sqrt(incident_mu))?;
        let b = transport.magnetic_field_t();
        let magnitude = real_norm(b);
        let reversed = if magnitude == 0.0 {
            false
        } else {
            let axis = b.map(|v| v / magnitude);
            let cross = real_cross(normal, axis);
            if !magnitude.is_finite() || real_norm(cross) > 512.0 * f64::EPSILON {
                return Err(bad(
                    "magnetized half-space needs propagation parallel to the original B",
                ));
            }
            real_dot(normal, axis) < 0.0
        };
        let sigma = transport.circular_conductivities_s_m();
        let conductivity = if reversed {
            [sigma[2], sigma[1]]
        } else {
            [sigma[1], sigma[2]]
        };
        let circular = |sigma: C64| -> Result<InterfaceReflection> {
            // Keep cancellation near the cutoff in the existing double-double
            // arithmetic. The original carrier drag supplies all electric loss.
            let epsilon = C64::new(
                (Dd::from_f64(epsilon) + Dd::from_f64(sigma.im) / Dd::from_f64(omega)).to_f64(),
                -positive(sigma.re / omega)?,
            );
            let medium = WaveMedium::new(epsilon, mu, incident_epsilon, incident_mu)?;
            interface_response(medium, wave_number, 1.0, Polarization::Te)
        };
        let modes = [circular(conductivity[0])?, circular(conductivity[1])?];
        poll(&mut cancelled)?;
        Ok(Self {
            incident,
            transport,
            background,
            normal,
            incident_admittance_s: admittance,
            background_permittivity_f_m: epsilon,
            permeability_h_m: mu,
            conductivity_s_m: conductivity,
            modes,
        })
    }

    pub const fn transport(&self) -> &'a IonizedGasMagnetizedHarmonic<'a> {
        self.transport
    }
    /// The explicit bound-charge background, separate from the carrier source.
    pub const fn background(&self) -> &'a ResolvedMaterialStatePoint {
        self.background
    }
    pub const fn incident_response(&self) -> DielectricResponse<'a> {
        self.incident
    }
    pub const fn outgoing_normal(&self) -> [f64; 3] {
        self.normal
    }
    /// Plus/minus conductivity in the outgoing-normal circular basis [S/m].
    pub const fn circular_conductivities_s_m(&self) -> [C64; 2] {
        self.conductivity_s_m
    }
    /// Plus/minus native Fresnel amplitudes, passive k and physical port powers.
    /// These are normalized mode data, not an admission of a physical amplitude.
    pub const fn circular_modes(&self) -> &[InterfaceReflection; 2] {
        &self.modes
    }

    /// Admit the actual Cartesian peak incident E [V/m], with any transverse
    /// Jones polarization. A longitudinal component beyond roundoff refuses.
    ///
    /// For each carrier, u = a_plus e_plus + a_minus e_minus has instantaneous
    /// peak speed (|a_plus|+|a_minus|)/sqrt(2), independent of relative phase.
    /// Both modal magnitudes decay with depth, so the original combined surface
    /// response bounds ALL deeper individual and pair-relative drift speeds.
    /// A zero-field conductivity probe alone would not establish this bound.
    pub fn prepare_wave(
        &self,
        incident_electric_v_m: Vector,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMagnetizedIonizedGasWave<'a>> {
        poll(&mut cancelled)?;
        let electric = transverse(self.normal, incident_electric_v_m)?;
        let circular = project(self.normal, electric);
        let transmitted: [Vector; 2] = std::array::from_fn(|m| {
            circular[m].map(|e| e * self.modes[m].transmitted_tangential_e)
        });
        let combined = add(transmitted[0], transmitted[1]);
        // Validate the physical superposition before accepting either mode.
        let surface = self.transport.response(combined, &mut cancelled)?;
        let modes = [
            self.prepare_mode(0, circular[0], transmitted[0], &mut cancelled)?,
            self.prepare_mode(1, circular[1], transmitted[1], &mut cancelled)?,
        ];
        let reflected = add(
            circular[0].map(|e| e * self.modes[0].reflected_tangential_e),
            circular[1].map(|e| e * self.modes[1].reflected_tangential_e),
        );
        let incident_power = sum(modes.iter().map(|m| m.incident_power_w_m2));
        let reflected_power = sum(modes.iter().map(|m| m.reflected_power_w_m2));
        let transmitted_power = sum(modes.iter().map(|m| m.normal_flux_w_m2));
        let heat = sum(modes
            .iter()
            .map(|m| m.material.mean_heating_w_m3.total_power_density_w_m3()));
        let kinetic = sum(modes
            .iter()
            .map(|m| m.material.mean_drift_kinetic_energy_j_m3));
        let residual = balance(incident_power, reflected_power + transmitted_power)?
            .max(balance(
                surface.mean_heating_w_m3.total_power_density_w_m3(),
                heat,
            )?)
            .max(balance(surface.mean_drift_kinetic_energy_j_m3, kinetic)?)
            .max(modes[0].relative_residual)
            .max(modes[1].relative_residual);
        if reflected.iter().any(|&e| !finite(e)) {
            return Err(bad("magnetized reflected field is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(PreparedMagnetizedIonizedGasWave {
            interface: *self,
            incident_electric_v_m: electric,
            reflected_electric_v_m: reflected,
            incident_power_w_m2: incident_power,
            reflected_power_w_m2: reflected_power,
            transmitted_power_w_m2: transmitted_power,
            surface,
            modes,
            power_balance_relative_residual: residual,
        })
    }

    fn prepare_mode(
        &self,
        mode: usize,
        incident: Vector,
        electric: Vector,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ModeWave<'a>> {
        poll(cancelled)?;
        let reflection = self.modes[mode];
        // The independent Fresnel H numerator preserves small magnetic fields
        // when reflection rounds to one; never recover it by subtracting 1-r.
        let magnetic = cross(self.normal, incident)
            .map(|e| (e * reflection.transmitted_conjugate_h).scale(self.incident_admittance_s));
        let material = self.transport.response(electric, &mut *cancelled)?;
        let poynting = complex_poynting(electric, magnetic);
        let normal_flux = dot(self.normal, poynting).re;
        let incident_scale =
            det::sqrt(self.incident_admittance_s) * core::f64::consts::FRAC_1_SQRT_2;
        let incident_power = sum(incident.iter().map(|e| e.scale(incident_scale).norm_sq()));
        let reflected_power = incident_power * reflection.reflected_power_fraction;
        let transmitted_power = incident_power * reflection.transmitted_power_fraction;
        let k = reflection.transmitted_normal_wavenumber_rad_m;
        let heat_from_flux = (-k.im * normal_flux) * 2.0;
        let omega = self.transport.omega_rad_s();
        let k_over_omega = C64::new(k.re / omega, k.im / omega);
        let faraday = cross(self.normal, electric).map(|e| e * k_over_omega);
        let magnetic_flux = magnetic.map(|h| h.scale(self.permeability_h_m));
        let curl_h = cross(self.normal, magnetic).map(|h| h * k_over_omega);
        let displacement = electric.map(|e| e.scale(self.background_permittivity_f_m));
        let free_current = material
            .current_density_a_m2
            .map(|j| C64::new(-j.im / omega, j.re / omega));
        let relative_residual = balance(normal_flux, transmitted_power)?
            .max(balance(
                heat_from_flux,
                material.mean_heating_w_m3.total_power_density_w_m3(),
            )?)
            .max(vector_residual(
                subtract(faraday, magnetic_flux),
                norm(faraday) + norm(magnetic_flux),
            )?)
            .max(vector_residual(
                subtract(add(curl_h, displacement), free_current),
                norm(curl_h) + norm(displacement) + norm(free_current),
            )?);
        let electric_scale = 0.5 * det::sqrt(self.background_permittivity_f_m);
        let magnetic_scale = 0.5 * det::sqrt(self.permeability_h_m);
        let energy = sum(electric
            .iter()
            .map(|e| e.scale(electric_scale).norm_sq())
            .chain(magnetic.iter().map(|h| h.scale(magnetic_scale).norm_sq())));
        if magnetic.iter().chain(&poynting).any(|&z| !finite(z))
            || !energy.is_finite()
            || energy < 0.0
        {
            return Err(bad(
                "magnetized wave flux or background storage is unrepresentable",
            ));
        }
        Ok(ModeWave {
            material,
            magnetic_a_m: magnetic,
            poynting_w_m2: poynting,
            normal_flux_w_m2: normal_flux,
            background_energy_j_m3: energy,
            incident_power_w_m2: incident_power,
            reflected_power_w_m2: reflected_power,
            relative_residual,
        })
    }
}

struct ModeWave<'a> {
    material: IonizedGasMagnetizedHarmonicResponse<'a>,
    magnetic_a_m: Vector,
    poynting_w_m2: Vector,
    normal_flux_w_m2: f64,
    background_energy_j_m3: f64,
    incident_power_w_m2: f64,
    reflected_power_w_m2: f64,
    relative_residual: f64,
}

/// One physically admitted superposition of the two outgoing circular waves.
/// Preparation evaluates the original transport three times (combined surface
/// plus its modes). Depth queries are allocation-free O(1), retain signed
/// thermal partitions, and never repeat carrier factorization or transport.
pub struct PreparedMagnetizedIonizedGasWave<'a> {
    interface: MagnetizedIonizedGasInterface<'a>,
    incident_electric_v_m: Vector,
    reflected_electric_v_m: Vector,
    incident_power_w_m2: f64,
    reflected_power_w_m2: f64,
    transmitted_power_w_m2: f64,
    surface: IonizedGasMagnetizedHarmonicResponse<'a>,
    modes: [ModeWave<'a>; 2],
    power_balance_relative_residual: f64,
}

impl<'a> PreparedMagnetizedIonizedGasWave<'a> {
    pub const fn interface(&self) -> &MagnetizedIonizedGasInterface<'a> {
        &self.interface
    }
    pub const fn surface_response(&self) -> &IonizedGasMagnetizedHarmonicResponse<'a> {
        &self.surface
    }
    pub const fn incident_electric_v_m(&self) -> Vector {
        self.incident_electric_v_m
    }
    pub const fn reflected_electric_v_m(&self) -> Vector {
        self.reflected_electric_v_m
    }
    pub fn incident_magnetic_a_m(&self) -> Vector {
        cross(self.interface.normal, self.incident_electric_v_m)
            .map(|e| e.scale(self.interface.incident_admittance_s))
    }
    pub fn reflected_magnetic_a_m(&self) -> Vector {
        cross(self.interface.normal, self.reflected_electric_v_m)
            .map(|e| e.scale(-self.interface.incident_admittance_s))
    }
    pub const fn incident_power_w_m2(&self) -> f64 {
        self.incident_power_w_m2
    }
    pub const fn reflected_power_w_m2(&self) -> f64 {
        self.reflected_power_w_m2
    }
    /// Flux entering the plasma, not power surviving at infinity.
    pub const fn transmitted_power_w_m2(&self) -> f64 {
        self.transmitted_power_w_m2
    }
    /// Native port, local Maxwell and original material-power residual.
    /// A numerical residual is not a continuum or material error bound.
    pub const fn power_balance_relative_residual(&self) -> f64 {
        self.power_balance_relative_residual
    }

    /// Evaluate depth >= 0 along the admitted outgoing normal [m].
    /// Each mode retains its own complex phase and attenuation. Orthogonal
    /// circular modes have no cross term in cycle-mean carrier/pair speed
    /// squared, storage or normal Poynting flux, even at different complex k.
    pub fn field(&self, depth_m: f64) -> Result<MagnetizedIonizedGasWavePoint<'_, 'a>> {
        if !depth_m.is_finite() || depth_m < 0.0 {
            return Err(bad("magnetized wave depth must be finite and nonnegative"));
        }
        let mut attenuation = [0.0; 2];
        let mut propagation = [C64::ZERO; 2];
        let mut absorbed = [0.0; 2];
        for mode in 0..2 {
            let k = self.interface.modes[mode].transmitted_normal_wavenumber_rad_m;
            let decay = k.im * depth_m;
            let a = det::exp(decay);
            if !a.is_finite() || !(0.0..=1.0).contains(&a) {
                return Err(bad("magnetized outgoing attenuation is unrepresentable"));
            }
            attenuation[mode] = a;
            let active = self.modes[mode]
                .material
                .electric_field_v_m
                .iter()
                .any(|&e| e != C64::ZERO);
            if a > 0.0 && active {
                let phase = Dd::from_f64(k.re) * Dd::from_f64(depth_m);
                if !phase.to_f64().is_finite() {
                    return Err(bad("magnetized propagation phase is unrepresentable"));
                }
                let (sine, cosine) = det::sin_cos_dd(phase);
                propagation[mode] = C64::new(a * cosine.to_f64(), -a * sine.to_f64());
            }
            // Preserve weak absorption without subtracting nearby port powers.
            absorbed[mode] = self.modes[mode].normal_flux_w_m2 * -det::expm1(2.0 * decay);
        }
        let combine = |a: Vector, b: Vector| -> Vector {
            std::array::from_fn(|d| a[d] * propagation[0] + b[d] * propagation[1])
        };
        let electric = combine(
            self.modes[0].material.electric_field_v_m,
            self.modes[1].material.electric_field_v_m,
        );
        let magnetic = combine(self.modes[0].magnetic_a_m, self.modes[1].magnetic_a_m);
        let current = combine(
            self.modes[0].material.current_density_a_m2,
            self.modes[1].material.current_density_a_m2,
        );
        let poynting = std::array::from_fn(|d| {
            C64::new(
                sum((0..2).map(|m| scaled_power(self.modes[m].poynting_w_m2[d].re, attenuation[m]))),
                sum((0..2).map(|m| scaled_power(self.modes[m].poynting_w_m2[d].im, attenuation[m]))),
            )
        });
        let absorbed_power = sum(absorbed.into_iter());
        if electric
            .iter()
            .chain(&magnetic)
            .chain(&current)
            .chain(&poynting)
            .any(|&z| !finite(z))
            || !absorbed_power.is_finite()
            || absorbed_power < 0.0
        {
            return Err(bad(
                "magnetized field or accumulated heat is unrepresentable",
            ));
        }
        Ok(MagnetizedIonizedGasWavePoint {
            wave: self,
            depth_m,
            propagation,
            attenuation,
            electric_v_m: electric,
            magnetic_a_m: magnetic,
            current_density_a_m2: current,
            complex_poynting_w_m2: poynting,
            absorbed_power_to_depth_w_m2: absorbed_power,
        })
    }
}

/// Physical Cartesian fields and cycle means at one point on the normal ray.
pub struct MagnetizedIonizedGasWavePoint<'w, 'a> {
    wave: &'w PreparedMagnetizedIonizedGasWave<'a>,
    depth_m: f64,
    propagation: [C64; 2],
    attenuation: [f64; 2],
    pub electric_v_m: Vector,
    pub magnetic_a_m: Vector,
    /// Free-carrier current, excluding the bound-charge displacement current.
    pub current_density_a_m2: Vector,
    /// E cross conj(H)/2, including reactive flux [W/m²].
    pub complex_poynting_w_m2: Vector,
    /// Physical total friction heat integrated from the interface [W/m²].
    pub absorbed_power_to_depth_w_m2: f64,
}

impl<'w, 'a> MagnetizedIonizedGasWavePoint<'w, 'a> {
    pub const fn wave(&self) -> &'w PreparedMagnetizedIonizedGasWave<'a> {
        self.wave
    }
    pub const fn depth_m(&self) -> f64 {
        self.depth_m
    }
    pub fn carrier_velocity_m_s(&self, carrier: usize) -> Option<Vector> {
        let first = self.wave.modes[0]
            .material
            .carrier_velocity_m_s
            .get(carrier)?;
        let second = self.wave.modes[1]
            .material
            .carrier_velocity_m_s
            .get(carrier)?;
        Some(std::array::from_fn(|d| {
            first[d] * self.propagation[0] + second[d] * self.propagation[1]
        }))
    }
    fn scaled_response(
        &self,
        value: impl Fn(&IonizedGasMagnetizedHarmonicResponse<'a>) -> f64,
    ) -> f64 {
        sum((0..2).map(|m| scaled_power(value(&self.wave.modes[m].material), self.attenuation[m])))
    }
    /// Total physical friction heat; thermal consumers deposit it only once.
    pub fn mean_heating_w_m3(&self) -> f64 {
        self.scaled_response(|m| m.mean_heating_w_m3.total_power_density_w_m3())
    }
    /// Signed electron-pool power, including the original drift correction.
    pub fn mean_electron_heating_w_m3(&self) -> f64 {
        self.scaled_response(|m| m.mean_heating_w_m3.electron_power_density_w_m3())
    }
    pub fn mean_heavy_heating_w_m3(&self) -> f64 {
        self.scaled_response(|m| m.mean_heating_w_m3.heavy_power_density_w_m3())
    }
    /// Already included in the pool powers; baseline thermal exchange is separate.
    pub fn mean_exchange_correction_to_electrons_w_m3(&self) -> f64 {
        self.scaled_response(|m| {
            m.mean_heating_w_m3
                .thermal_exchange_correction_to_electrons_w_m3()
        })
    }
    pub fn mean_drift_kinetic_energy_j_m3(&self) -> f64 {
        self.scaled_response(|m| m.mean_drift_kinetic_energy_j_m3)
    }
    /// Lossless bound-charge electric plus wave magnetic storage. Complex
    /// permittivity is not used as an energy density and excludes static B.
    pub fn mean_background_field_energy_j_m3(&self) -> f64 {
        sum((0..2).map(|m| {
            scaled_power(
                self.wave.modes[m].background_energy_j_m3,
                self.attenuation[m],
            )
        }))
    }
    pub fn reactive_carrier_power_var_m3(&self) -> f64 {
        self.scaled_response(|m| m.reactive_power_var_m3)
    }
    /// Reactive electrical power equals 2*omega*kinetic storage plus this signed
    /// gyroscopic term. Neither contribution is deposited as material heat.
    pub fn gyroscopic_reactive_power_var_m3(&self) -> f64 {
        self.scaled_response(|m| m.gyroscopic_reactive_power_var_m3)
    }
    /// The prescribed static B does no real work, to the original residual.
    pub fn mean_magnetic_power_w_m3(&self) -> f64 {
        self.scaled_response(|m| m.mean_magnetic_power_w_m3)
    }
}

fn sum(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(Dd::ZERO, |s, x| s + Dd::from_f64(x)).to_f64()
}
fn dot(n: [f64; 3], v: Vector) -> C64 {
    C64::new(
        sum((0..3).map(|d| n[d] * v[d].re)),
        sum((0..3).map(|d| n[d] * v[d].im)),
    )
}
fn real_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    sum((0..3).map(|d| a[d] * b[d]))
}
fn norm(v: Vector) -> f64 {
    det::hypot(det::hypot(v[0].abs(), v[1].abs()), v[2].abs())
}
fn real_norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn unit(v: [f64; 3]) -> Result<[f64; 3]> {
    if v.iter().any(|x| !x.is_finite()) {
        return Err(bad("magnetized outgoing normal must be finite and nonzero"));
    }
    let length = positive(real_norm(v))?;
    Ok(v.map(|x| x / length))
}
fn real_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        a[j] * b[k] - a[k] * b[j]
    })
}
fn cross(a: [f64; 3], b: Vector) -> Vector {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        b[k].scale(a[j]) - b[j].scale(a[k])
    })
}
fn add(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|d| a[d] + b[d])
}
fn subtract(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|d| a[d] - b[d])
}
fn transverse(normal: [f64; 3], field: Vector) -> Result<Vector> {
    if field.iter().any(|&e| !finite(e)) {
        return Err(bad("magnetized incident field must be finite"));
    }
    let longitudinal = dot(normal, field);
    if longitudinal.abs() > 512.0 * f64::EPSILON * norm(field) {
        return Err(bad(
            "normal-incidence magnetized wave requires transverse electric field",
        ));
    }
    Ok(subtract(field, normal.map(|n| longitudinal.scale(n))))
}
fn project(normal: [f64; 3], field: Vector) -> [Vector; 2] {
    let rotated = cross(normal, field).map(|z| C64::new(-z.im, z.re));
    [
        add(field.map(|z| z.scale(0.5)), rotated.map(|z| z.scale(0.5))),
        subtract(field.map(|z| z.scale(0.5)), rotated.map(|z| z.scale(0.5))),
    ]
}
fn vector_residual(residual: Vector, scale: f64) -> Result<f64> {
    let magnitude = norm(residual);
    let relative = if scale == 0.0 && magnitude == 0.0 {
        0.0
    } else {
        magnitude / scale
    };
    if !scale.is_finite() || !relative.is_finite() || relative > 1e-10 {
        return Err(bad("magnetized mode violates original Maxwell equations"));
    }
    Ok(relative)
}
