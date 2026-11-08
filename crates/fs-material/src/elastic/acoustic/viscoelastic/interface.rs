//! Oblique sound reflection from an original causal isotropic solid half-space.
//!
//! The lossless fluid occupies z<0 and the solid z>0. Peak phasors use
//! exp(+i omega t), exp(-i omega (p*x+q*z)). The real tangential slowness p is
//! conserved; q_P²=rho/(K*+4G*/3)-p² and q_S²=rho/G*-p² use outgoing roots
//! with Im(q)<=0. Normal velocity/traction are continuous and shear traction
//! vanishes. Tangential velocity may slip. No real refracted angle is assigned
//! to an inhomogeneous attenuating wave.
//!
//! Strain is -sym(s tensor v), and stress is K* tr(strain) I+2G* dev(strain).
//! Normal power is -Re(v.conj().dot(sigma.e_z))/2. Total P/S fields are combined
//! BEFORE power and constitutive heat are evaluated: lossy converted waves
//! have interference fluxes. See Carcione (1997), GJI 129, 669-680,
//! <https://doi.org/10.1111/j.1365-246X.1997.tb04502.x>, and Cervený & Pšenčík
//! (2006), GJI 166, 1299-1317, §§2-3,
//! <https://doi.org/10.1111/j.1365-246X.2006.03057.x>. The latter uses the
//! conjugate time convention. Here positive Im(K*,G*) dissipates energy.
//!
//! This fixed homogeneous planar contact retains the original bulk/shear
//! Maxwell sources, frequency band, density and small-strain support. It does
//! not infer material authority, anisotropy, roughness, finite thickness,
//! thermal evolution, an incident solid wave or a moving interface.

use super::{ViscoelasticSlabError as Error, ViscoelasticSlabSource, poll, source_moduli};
use crate::elastic::acoustic::interface::fluid::FluidDriveLimits;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};

#[derive(Clone, Copy)]
struct Mode {
    q: C64,
    velocity: [C64; 3],
    strain: [C64; 6],
    trace: C64,
    stress: [C64; 6],
}

/// Prepared at one frequency and fluid incidence angle. Only preparation
/// allocates (bounded three-by-three LU); drives and depth queries do not.
pub struct FluidViscoelasticInterface<'a> {
    source: ViscoelasticSlabSource<'a>,
    fluid_source: &'a dyn FluidAcoustics,
    fluid: FluidAcousticState,
    moduli: [C64; 2],
    omega: f64,
    cosine: f64,
    sine: f64,
    tangential: f64,
    modes: [Mode; 2],
    // Reflected fluid, transmitted P, transmitted SV modal velocities,
    // divided by the incoming fluid's longitudinal particle velocity.
    coefficients: [C64; 3],
}

/// Peak fields at x=0 and a caller-selected solid depth. The same common
/// tangential phase multiplies all fields at other x; mean power is unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticHalfSpaceField {
    pub velocity_m_s: [C64; 3],
    /// [xx, yy, zz, sqrt(2)xy, sqrt(2)yz, sqrt(2)zx].
    pub strain_mandel: [C64; 6],
    /// Retained analytically, avoiding cancellation between large opposite
    /// normal strains in evanescent P/S fields.
    pub volumetric_strain: C64,
    pub stress_mandel_pa: [C64; 6],
    /// sigma.e_z, with tensile stress positive.
    pub traction_pa: [C64; 3],
    pub normal_power_w_m2: f64,
    /// omega Im(K*) |tr(strain)|² / 2.
    pub bulk_heat_density_w_m3: f64,
    /// omega Im(G*) dev(strain):dev(strain)*.
    pub shear_heat_density_w_m3: f64,
}

/// Independent fluid traveling powers and total coherent boundary fluxes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticInterfacePower {
    pub incident_power_w_m2: f64,
    pub reflected_power_w_m2: f64,
    pub fluid_boundary_power_w_m2: f64,
    /// Flux entering the solid; not heat deposited at the mathematical contact.
    /// A lossless transmitted wave can carry this power to infinity.
    pub solid_boundary_power_w_m2: f64,
    pub balance_residual_w_m2: f64,
    pub solid_strain_norm_bound: f64,
    pub fluid_relative_density_bound: f64,
    pub fluid_particle_mach_bound: f64,
}

/// An admitted actual pressure excitation of the same prepared interface.
pub struct DrivenFluidViscoelasticInterface<'s, 'a> {
    interface: &'s FluidViscoelasticInterface<'a>,
    amplitudes: [C64; 2],
    reflected_pressure: C64,
    fluid_pressure: C64,
    fluid_velocity: [C64; 3],
    power: ViscoelasticInterfacePower,
    boundary: ViscoelasticHalfSpaceField,
}

impl<'a> FluidViscoelasticInterface<'a> {
    /// Prepare from one original fluid snapshot and the original solid laws.
    /// The cosine must be in (0,1]; exact grazing has no incoming normal power.
    /// The branch allowance counts bulk and shear together, as for the slab.
    pub fn try_new(
        source: ViscoelasticSlabSource<'a>,
        fluid_source: &'a dyn FluidAcoustics,
        angular_frequency_rad_s: f64,
        incidence_cosine: f64,
        maximum_branches: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        let moduli = source_moduli(source, angular_frequency_rad_s, maximum_branches, cancelled)?;
        if !source.density_kg_m3.is_finite()
            || source.density_kg_m3 <= 0.0
            || moduli
                .iter()
                .any(|m| !finite(*m) || m.re <= 0.0 || m.im < 0.0)
            || !incidence_cosine.is_finite()
            || incidence_cosine <= 0.0
            || incidence_cosine > 1.0
        {
            return Err(Error::Invalid(
                "invalid viscoelastic interface material or angle",
            ));
        }
        let fluid = fluid_source.acoustic_state()?;
        let sine = det::sqrt((1.0 - incidence_cosine) * (1.0 + incidence_cosine));
        let tangential = sine / fluid.sound_speed_m_s();
        if !tangential.is_finite() || (sine > 0.0 && tangential == 0.0) {
            return Err(Error::Invalid("tangential slowness is unrepresentable"));
        }
        let modes = [
            mode(source.density_kg_m3, moduli, tangential, false)?,
            mode(source.density_kg_m3, moduli, tangential, true)?,
        ];
        let z = fluid.characteristic_impedance_pa_s_m();
        let mut matrix = [
            C64::from_re(incidence_cosine),
            modes[0].velocity[2],
            modes[1].velocity[2],
            C64::from_re(z),
            modes[0].stress[2],
            modes[1].stress[2],
            C64::ZERO,
            modes[0].stress[5],
            modes[1].stress[5],
        ];
        let mut rhs = [C64::from_re(incidence_cosine), C64::from_re(-z), C64::ZERO];
        for row in 0..3 {
            let scale = matrix[3 * row..3 * row + 3]
                .iter()
                .fold(0.0_f64, |s, x| s.max(x.abs()));
            if !scale.is_finite() || scale == 0.0 {
                return Err(Error::Invalid(
                    "viscoelastic contact row is unrepresentable",
                ));
            }
            for value in &mut matrix[3 * row..3 * row + 3] {
                *value = divide(*value, scale);
            }
            rhs[row] = divide(rhs[row], scale);
        }
        let factor = lu_complex_with_cancel(&matrix, 3, &mut *cancelled).map_err(solve_error)?;
        let mut coefficients = rhs;
        factor
            .solve_with_cancel(&mut coefficients, &mut *cancelled)
            .map_err(solve_error)?;
        if coefficients.iter().any(|x| !finite(*x)) {
            return Err(Error::Invalid("viscoelastic reflection is unrepresentable"));
        }
        for row in 0..3 {
            let terms: [C64; 3] = core::array::from_fn(|c| matrix[3 * row + c] * coefficients[c]);
            let residual = (terms.into_iter().fold(C64::ZERO, |a, b| a + b) - rhs[row]).abs();
            let scale = terms.iter().map(|x| x.abs()).sum::<f64>() + rhs[row].abs();
            if !residual.is_finite() || (scale > 0.0 && residual / scale > 1e-11) {
                return Err(Error::Invalid(
                    "viscoelastic contact continuity did not converge",
                ));
            }
        }
        poll(cancelled)?;
        Ok(Self {
            source,
            fluid_source,
            fluid,
            moduli,
            omega: angular_frequency_rad_s,
            cosine: incidence_cosine,
            sine,
            tangential,
            modes,
            coefficients,
        })
    }

    pub const fn source(&self) -> ViscoelasticSlabSource<'a> {
        self.source
    }
    pub fn fluid_source(&self) -> &'a dyn FluidAcoustics {
        self.fluid_source
    }
    pub const fn fluid_state(&self) -> FluidAcousticState {
        self.fluid
    }
    pub const fn moduli_pa(&self) -> [C64; 2] {
        self.moduli
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn tangential_slowness_s_m(&self) -> f64 {
        self.tangential
    }
    /// Pressure reflection at z=0; complex phase is retained.
    pub const fn pressure_reflection(&self) -> C64 {
        self.coefficients[0]
    }
    /// Outgoing P and SV normal slownesses [s/m], in that order.
    pub fn normal_slownesses_s_m(&self) -> [C64; 2] {
        self.modes.map(|m| m.q)
    }

    /// Admit the actual complex peak pressure and caller-declared fluid
    /// small-signal limits. The triangle strain bound holds at every z>=0,
    /// including converted evanescent waves and P/S interference.
    pub fn drive<'s>(
        &'s self,
        incident_pressure_pa: C64,
        limits: FluidDriveLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DrivenFluidViscoelasticInterface<'s, 'a>, Error> {
        poll(cancelled)?;
        if !finite(incident_pressure_pa) {
            return Err(Error::Invalid("incident pressure must be finite"));
        }
        let z = self.fluid.characteristic_impedance_pa_s_m();
        let incoming = divide(incident_pressure_pa, z);
        let reflected = incident_pressure_pa * self.coefficients[0];
        let amplitudes = [
            self.coefficients[1] * incoming,
            self.coefficients[2] * incoming,
        ];
        let fluid_bound =
            (incoming.abs() + divide(reflected, z).abs()) / self.fluid.sound_speed_m_s();
        let strain_bound = self
            .modes
            .iter()
            .zip(amplitudes)
            .map(|(m, a)| m.strain.iter().fold(0.0, |n, e| det::hypot(n, e.abs())) * a.abs())
            .sum::<f64>();
        if !fluid_bound.is_finite()
            || !strain_bound.is_finite()
            || fluid_bound > limits.maximum_relative_density()
            || fluid_bound > limits.maximum_particle_mach()
            || strain_bound > self.source.maximum_absolute_strain
        {
            return Err(Error::OutsideSourceSupport);
        }
        let boundary = self.field(amplitudes, 0.0)?;
        let fluid_pressure = incident_pressure_pa + reflected;
        let fluid_velocity = [
            divide(fluid_pressure, z).scale(self.sine),
            C64::ZERO,
            divide(incident_pressure_pa - reflected, z).scale(self.cosine),
        ];
        let incoming_power = quadratic(0.5 * self.cosine * z, incoming);
        let reflected_power = quadratic(0.5 * self.cosine * z, divide(reflected, z));
        let fluid_power = (product(fluid_pressure, fluid_velocity[2]) * Dd::from_f64(0.5)).to_f64();
        let residual =
            (incoming_power - reflected_power - Dd::from_f64(boundary.normal_power_w_m2)).to_f64();
        let incident = incoming_power.to_f64();
        let reflected_w = reflected_power.to_f64();
        if !finite(fluid_pressure)
            || fluid_velocity.iter().any(|v| !finite(*v))
            || [incident, reflected_w, fluid_power, residual]
                .iter()
                .any(|x| !x.is_finite())
            || (incident_pressure_pa != C64::ZERO && incident == 0.0)
            || (incident > 0.0
                && (residual.abs() / incident > 1e-10
                    || (fluid_power - boundary.normal_power_w_m2).abs() / incident > 1e-10
                    || boundary.normal_power_w_m2 < -1e-10 * incident))
        {
            return Err(Error::Invalid(
                "viscoelastic interface power is unrepresentable or unbalanced",
            ));
        }
        poll(cancelled)?;
        Ok(DrivenFluidViscoelasticInterface {
            interface: self,
            amplitudes,
            reflected_pressure: reflected,
            fluid_pressure,
            fluid_velocity,
            boundary,
            power: ViscoelasticInterfacePower {
                incident_power_w_m2: incident,
                reflected_power_w_m2: reflected_w,
                fluid_boundary_power_w_m2: fluid_power,
                solid_boundary_power_w_m2: boundary.normal_power_w_m2,
                balance_residual_w_m2: residual,
                solid_strain_norm_bound: strain_bound,
                fluid_relative_density_bound: fluid_bound,
                fluid_particle_mach_bound: fluid_bound,
            },
        })
    }

    fn field(&self, amplitudes: [C64; 2], depth: f64) -> Result<ViscoelasticHalfSpaceField, Error> {
        let a = [
            amplitudes[0] * propagation(self.modes[0].q, self.omega, depth)?,
            amplitudes[1] * propagation(self.modes[1].q, self.omega, depth)?,
        ];
        let velocity = core::array::from_fn(|i| {
            self.modes[0].velocity[i] * a[0] + self.modes[1].velocity[i] * a[1]
        });
        let strain: [C64; 6] = core::array::from_fn(|i| {
            self.modes[0].strain[i] * a[0] + self.modes[1].strain[i] * a[1]
        });
        let stress: [C64; 6] = core::array::from_fn(|i| {
            self.modes[0].stress[i] * a[0] + self.modes[1].stress[i] * a[1]
        });
        let traction = [
            stress[5].scale(core::f64::consts::FRAC_1_SQRT_2),
            C64::ZERO,
            stress[2],
        ];
        let trace = self.modes[0].trace * a[0] + self.modes[1].trace * a[1];
        let mean = trace.scale(1.0 / 3.0);
        let bulk = (quadratic(0.5 * self.moduli[0].im, trace) * Dd::from_f64(self.omega)).to_f64();
        let shear = (strain.iter().enumerate().fold(Dd::ZERO, |s, (i, &e)| {
            s + quadratic(self.moduli[1].im, if i < 3 { e - mean } else { e })
        }) * Dd::from_f64(self.omega))
        .to_f64();
        let flux = ((0..3).fold(Dd::ZERO, |s, i| s - product(velocity[i], traction[i]))
            * Dd::from_f64(0.5))
        .to_f64();
        if velocity
            .iter()
            .chain(&strain)
            .chain(&stress)
            .any(|v| !finite(*v))
            || [bulk, shear, flux].iter().any(|x| !x.is_finite())
            || bulk < 0.0
            || shear < 0.0
        {
            return Err(Error::Invalid(
                "viscoelastic half-space field is unrepresentable",
            ));
        }
        Ok(ViscoelasticHalfSpaceField {
            velocity_m_s: velocity,
            strain_mandel: strain,
            volumetric_strain: trace,
            stress_mandel_pa: stress,
            traction_pa: traction,
            normal_power_w_m2: flux,
            bulk_heat_density_w_m3: bulk,
            shear_heat_density_w_m3: shear,
        })
    }
}

impl DrivenFluidViscoelasticInterface<'_, '_> {
    pub const fn power(&self) -> ViscoelasticInterfacePower {
        self.power
    }
    pub const fn reflected_pressure_pa(&self) -> C64 {
        self.reflected_pressure
    }
    pub const fn fluid_pressure_pa(&self) -> C64 {
        self.fluid_pressure
    }
    pub const fn fluid_velocity_m_s(&self) -> [C64; 3] {
        self.fluid_velocity
    }
    pub const fn solid_boundary_field(&self) -> ViscoelasticHalfSpaceField {
        self.boundary
    }
    /// Physical coherent P/S fields at a nonnegative finite depth in the solid.
    pub fn field(
        &self,
        depth_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticHalfSpaceField, Error> {
        poll(cancelled)?;
        if !depth_m.is_finite() || depth_m < 0.0 {
            return Err(Error::Invalid("solid depth must be finite and nonnegative"));
        }
        let field = self.interface.field(self.amplitudes, depth_m)?;
        poll(cancelled)?;
        Ok(field)
    }
}

fn mode(density: f64, moduli: [C64; 2], p: f64, shear: bool) -> Result<Mode, Error> {
    let [bulk, g] = moduli;
    let modulus = if shear { g } else { bulk + g.scale(4.0 / 3.0) };
    let root_modulus = modulus.sqrt();
    let speed = divide(root_modulus, det::sqrt(density));
    // sqrt(rho-M*p²)/sqrt(M) avoids first squaring the tiny slowness or
    // losing a representable attenuation term in rho/M. Evaluate M*p² with
    // the dimensional coefficient first, and compensate critical subtraction.
    let p_dd = Dd::from_f64(p);
    let real = (Dd::from_f64(density) - (Dd::from_f64(modulus.re) * p_dd) * p_dd).to_f64();
    let imaginary = -((Dd::from_f64(modulus.im) * p_dd) * p_dd).to_f64();
    let mut q = C64::new(real, imaginary).sqrt() / root_modulus;
    // Principal sqrt chooses +i for a negative lossless argument; the outgoing
    // evanescent branch in this convention is -i. Lossy roots already have -Im.
    if q.im > 0.0 {
        q = -q;
    }
    let velocity = if shear {
        [q * speed, C64::ZERO, -speed.scale(p)]
    } else {
        [speed.scale(p), C64::ZERO, q * speed]
    };
    let strain = [
        -velocity[0].scale(p),
        C64::ZERO,
        -(q * velocity[2]),
        C64::ZERO,
        C64::ZERO,
        -(q * velocity[0] + velocity[2].scale(p)).scale(core::f64::consts::FRAC_1_SQRT_2),
    ];
    // q²+p²=rho/M exactly for a P wave; an SV wave has identically zero
    // divergence. Keep those identities even when evanescent strains cancel.
    let trace = if shear { C64::ZERO } else { -C64::ONE / speed };
    let mean = trace.scale(1.0 / 3.0);
    let stress = core::array::from_fn(|i| {
        if i < 3 {
            bulk * trace + g.scale(2.0) * (strain[i] - mean)
        } else {
            g.scale(2.0) * strain[i]
        }
    });
    if !finite(speed)
        || speed == C64::ZERO
        || !finite(q)
        || q.re < 0.0
        || q.im > 0.0
        || (modulus.im > 0.0 && q.im == 0.0)
        || velocity
            .iter()
            .chain(&strain)
            .chain(&stress)
            .any(|v| !finite(*v))
    {
        return Err(Error::Invalid(
            "viscoelastic outgoing mode is unrepresentable",
        ));
    }
    Ok(Mode {
        q,
        velocity,
        strain,
        trace,
        stress,
    })
}

fn propagation(q: C64, omega: f64, depth: f64) -> Result<C64, Error> {
    let attenuation = (omega * q.im) * depth;
    if attenuation < -745.0 {
        return Ok(C64::ZERO);
    }
    let phase = (omega * q.re) * depth;
    if !attenuation.is_finite() || !phase.is_finite() {
        return Err(Error::Invalid(
            "viscoelastic propagation phase is unrepresentable",
        ));
    }
    let (s, c) = det::sin_cos(phase);
    Ok(C64::new(c, -s).scale(det::exp(attenuation)))
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn divide(z: C64, s: f64) -> C64 {
    C64::new(z.re / s, z.im / s)
}
fn quadratic(weight: f64, z: C64) -> Dd {
    let w = Dd::from_f64(weight);
    let re = Dd::from_f64(z.re);
    let im = Dd::from_f64(z.im);
    (w * re) * re + (w * im) * im
}
fn product(a: C64, b: C64) -> Dd {
    Dd::from_f64(a.re) * Dd::from_f64(b.re) + Dd::from_f64(a.im) * Dd::from_f64(b.im)
}
fn solve_error(e: LuComplexError) -> Error {
    match e {
        LuComplexError::Cancelled => Error::Cancelled,
        _ => Error::Invalid("viscoelastic contact system is singular or unrepresentable"),
    }
}
