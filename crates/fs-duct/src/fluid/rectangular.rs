//! Thermoviscous propagation through the actual uniform rectangular section.
//!
//! In +i omega telegraph bookkeeping, transverse diffusion solves
//! `(s - D Laplacian) h = s`, with h=0 on all four walls. Its mean is
//! `H_D = sum_(m,n odd) 64/(pi^4 m^2 n^2) s/(s+D lambda_mn)`.
//! One transverse sum is evaluated analytically with tanh; the remaining
//! Fourier prefix has explicit componentwise tail bounds. No hydraulic or
//! equal-area radius replaces this geometry. The original material supplies
//! both momentum and thermal diffusivity, density, sound speed and Cp/Cv.
//!
//! This is the uniform-pressure, fully developed, linear continuum model with
//! rigid no-slip isothermal walls: Stinson, JASA 89 (1991), 550-558,
//! <https://doi.org/10.1121/1.400379>. It does not resolve entrance/edge fields,
//! mean flow, moving walls, rarefaction, bulk viscosity or axial heat diffusion.
//! The explicit transverse-phase ceiling is a model envelope, not a 3-D error
//! certificate. Inlet/outlet and shape-transition equivalence remain the
//! caller's responsibility; the source material retains its own validity.

use super::{FluidDuctError, bad, positive};
use crate::{DuctError, SegmentWave, compose::TwoPort};
use core::f64::consts::PI;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_math::{c64::C64, dd::Dd, det};

/// Molecular reservoirs and wall heat diffusion in the same energy equation.
pub mod molecular;
/// Resolved two-dimensional thermal fields with independently authored walls.
pub mod thermal_wall;
mod field;
pub use field::{
    PreparedRectangularFields, RectangularDissipationPerLength, RectangularField,
};

/// Hard bound on the one-dimensional Fourier work and cached geometry.
pub const MAX_RECTANGULAR_TERMS: usize = 4096;

/// Explicit numerical and plane-pressure admission; neither has a hidden default.
#[derive(Clone, Copy, Debug)]
pub struct RectangularAdmission {
    /// Relative omitted-series error in both diffusion components and the
    /// real wall fraction. Must lie in (0,1); excludes floating-point rounding.
    pub maximum_relative_tail: f64,
    /// Upper bound on max(omega/c, |k|) times the longer side, in (0,pi).
    /// The pi ceiling also excludes the first lossless transverse cutoff.
    pub maximum_transverse_phase: f64,
}

/// Mean transverse diffusion and its uncomputed infinite Fourier tail.
#[derive(Clone, Copy, Debug)]
pub struct RectangularDiffusion {
    /// H_D(s), under +i omega bookkeeping, with positive real/imaginary parts.
    pub mean: C64,
    /// 1-H_D, formed before rounding the accumulated real part to f64.
    pub wall: C64,
    /// Upper bounds on the omitted positive real and imaginary components.
    /// These are analytic truncation bounds, not interval roundoff enclosures.
    pub tail_real: f64,
    pub tail_imaginary: f64,
}

/// Native exp(-i omega t) wave, retaining the actual area and resolution evidence.
#[derive(Clone, Copy, Debug)]
pub struct RectangularWave {
    /// Pressure/mean-velocity and pressure/volume-flow impedances. The shear
    /// diagnostic uses HALF the shorter side, not an equivalent circular radius.
    pub wave: SegmentWave,
    pub viscous: RectangularDiffusion,
    /// None exactly when Cp=Cv, so uncoupled thermal data do no numerical work.
    pub thermal: Option<RectangularDiffusion>,
    pub transverse_phase: f64,
    area_m2: f64,
}

impl RectangularWave {
    /// Stable finite-length pressure/volume-flow block for the original area.
    pub fn two_port(&self, length_m: f64, reference: f64) -> Result<TwoPort, DuctError> {
        TwoPort::from_uniform_wave(&self.wave, self.area_m2, length_m, reference)
    }
}

#[derive(Clone, Copy, Debug)]
struct Mode {
    eigenvalue: f64,
    weight: f64,
}

/// One analytically summed long-direction solution, shared by the mean wave
/// and its resolved transverse field. All quantities use +i omega time.
struct ModeResponse {
    q: C64,
    decay: C64,
    factor: C64,
    tanh: C64,
    mean: C64,
}

fn mode_response(mode: &Mode, long_m: f64, argument: C64) -> Result<ModeResponse, DuctError> {
    let denominator = C64::new(mode.eigenvalue + argument.re, argument.im);
    let q = denominator.sqrt();
    let z = q.scale(0.5 * long_m);
    if !finite(z) {
        return Err(bad(
            "rectangular transverse diffusion argument is unrepresentable",
        ));
    }
    // Re(z)>=pi/2 because the shorter side owns the Fourier series.
    // Only a decaying exponential is formed, including for opaque layers.
    let decay = det::exp(-2.0 * z.re);
    let decay = if decay == 0.0 {
        C64::ZERO
    } else {
        let (sine, cosine) = det::sin_cos(-2.0 * z.im);
        C64::new(cosine, sine).scale(decay)
    };
    let tanh = if decay == C64::ZERO {
        C64::ONE
    } else {
        (C64::ONE - decay) / (C64::ONE + decay)
    };
    let factor = argument / denominator;
    let mean = (C64::ONE - tanh / z) * factor;
    Ok(ModeResponse {
        q,
        decay,
        factor,
        tanh,
        mean,
    })
}

/// Prepared geometry of a uniform duct; width/height exchange is identical.
/// Preparation is O(N); a frequency evaluation is O(N), allocation-free and
/// polls cancellation every 32 Fourier terms. The longer direction is summed
/// analytically, so a very thin slit needs no aspect-ratio-sized 2-D grid.
#[derive(Clone, Debug)]
pub struct RectangularSection {
    short_m: f64,
    long_m: f64,
    area_m2: f64,
    modes: Vec<Mode>,
}

impl RectangularSection {
    pub fn new(width_m: f64, height_m: f64, terms: usize) -> Result<Self, DuctError> {
        let area_m2 = width_m * height_m;
        if !positive(width_m)
            || !positive(height_m)
            || !positive(area_m2)
            || !(1..=MAX_RECTANGULAR_TERMS).contains(&terms)
        {
            return Err(bad(
                "rectangular section needs positive dimensions and 1..=4096 terms",
            ));
        }
        let short_m = width_m.min(height_m);
        let long_m = width_m.max(height_m);
        let mut modes = Vec::new();
        modes
            .try_reserve_exact(terms)
            .map_err(|_| bad("rectangular Fourier allocation exceeds available memory"))?;
        for index in 0..terms {
            let phase = PI * (2 * index + 1) as f64;
            let transverse = phase / short_m;
            let eigenvalue = transverse * transverse;
            if !positive(eigenvalue) {
                return Err(bad("rectangular transverse eigenvalue is unrepresentable"));
            }
            modes.push(Mode {
                eigenvalue,
                weight: 8.0 / (phase * phase),
            });
        }
        Ok(Self {
            short_m,
            long_m,
            area_m2,
            modes,
        })
    }

    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub fn terms(&self) -> usize {
        self.modes.len()
    }

    pub fn evaluate<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        admission: RectangularAdmission,
    ) -> Result<RectangularWave, FluidDuctError> {
        self.evaluate_controlled(material, omega, admission, || false)
    }

    /// Evaluate ONE original material snapshot. Unresolved Fourier tails or
    /// departure from the authored transverse-phase envelope refuse outright.
    pub fn evaluate_controlled<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        admission: RectangularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RectangularWave, FluidDuctError> {
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
        let acoustic = state.acoustic();
        let transport = state.transport();
        let k0 = omega / acoustic.sound_speed_m_s();
        if !positive(k0) || k0 * self.long_m > admission.maximum_transverse_phase {
            return Err(bad(
                "rectangular section is outside its plane-pressure frequency envelope",
            )
            .into());
        }
        let viscous =
            self.diffusion(omega, transport.kinematic_viscosity_m2_s(), &mut cancelled)?;
        admit_tail(viscous, admission.maximum_relative_tail)?;
        let gamma = state.heat_capacity_ratio();
        let thermal = if gamma == 1.0 {
            None
        } else {
            let response =
                self.diffusion(omega, transport.thermal_diffusivity_m2_s(), &mut cancelled)?;
            admit_tail(response, admission.maximum_relative_tail)?;
            Some(response)
        };
        // Exactly the existing +i omega telegraph convention, then its shared
        // branch conversion to exp(-i omega t). Thermal loss is counted once.
        let s = C64::new(0.0, omega);
        let series = s.scale(acoustic.density_kg_m3() / self.area_m2) / viscous.mean;
        let thermal_factor = C64::ONE + thermal.map_or(C64::ZERO, |t| t.wall.scale(gamma - 1.0));
        let shunt = s.scale(self.area_m2 / acoustic.adiabatic_bulk_modulus_pa()) * thermal_factor;
        let (wavenumber, characteristic_impedance) =
            crate::helmholtz_wave_from_telegraph(series, shunt);
        let specific_impedance = characteristic_impedance.scale(self.area_m2);
        let transverse_phase = k0.max(wavenumber.abs()) * self.long_m;
        if [wavenumber, characteristic_impedance, specific_impedance]
            .iter()
            .any(|v| !finite(*v))
            || wavenumber.re <= 0.0
            || wavenumber.im < 0.0
            || characteristic_impedance.re <= 0.0
            || !transverse_phase.is_finite()
            || transverse_phase > admission.maximum_transverse_phase
        {
            return Err(bad(
                "rectangular wave is unrepresentable or outside its transverse-phase envelope",
            )
            .into());
        }
        poll(&mut cancelled, self.modes.len())?;
        Ok(RectangularWave {
            wave: SegmentWave {
                wavenumber,
                characteristic_impedance,
                specific_impedance,
                shear_number: 0.5
                    * self.short_m
                    * det::sqrt(omega / transport.kinematic_viscosity_m2_s()),
            },
            viscous,
            thermal,
            transverse_phase,
            area_m2: self.area_m2,
        })
    }

    fn diffusion(
        &self,
        omega: f64,
        diffusivity: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RectangularDiffusion, DuctError> {
        let frequency = omega / diffusivity;
        if !positive(frequency) {
            return Err(bad("rectangular diffusion frequency is unrepresentable"));
        }
        self.diffusion_argument(C64::new(0.0, frequency), cancelled)
    }

    // H(z) = sum w_mn z/(z+lambda_mn). A relaxing heat capacity makes
    // z = s*rho*cp(s)/conductivity genuinely complex with Re(z)>=0.
    // The original pure-imaginary path retains its summation arithmetic.
    fn diffusion_argument(
        &self,
        argument: C64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RectangularDiffusion, DuctError> {
        if !finite(argument) || argument.re < 0.0 || argument.im <= 0.0 {
            return Err(bad(
                "rectangular diffusion needs a finite first-quadrant argument",
            ));
        }
        let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
        for (index, mode) in self.modes.iter().enumerate() {
            if index % 32 == 0 {
                poll(cancelled, index)?;
            }
            let value = mode_response(mode, self.long_m, argument)?.mean;
            re = re + Dd::from_f64(value.re) * Dd::from_f64(mode.weight);
            im = im + Dd::from_f64(value.im) * Dd::from_f64(mode.weight);
        }
        // Omitted modes have positive weights. For z=x+iy, x>=0, y>0:
        // Re[z/(z+lambda)]<=min(1,x/lambda+(y/lambda)^2),
        // Im[z/(z+lambda)]<=min(1/2,y/lambda). The midpoint integral
        // bounds the remaining odd sums: sum_(m>=2N+1) m^-p <=
        // (2N)^(1-p)/(2(p-1)). The other transverse weights sum to one.
        let n = self.modes.len() as f64;
        let ratio = argument.im / self.modes[0].eigenvalue;
        let real_ratio = argument.re / self.modes[0].eigenvalue;
        let mass = 2.0 / (PI * PI * n);
        let reduced = ratio / (n * n);
        let response = RectangularDiffusion {
            mean: C64::new(re.to_f64(), im.to_f64()),
            wall: C64::new((Dd::ONE - re).to_f64(), -im.to_f64()),
            tail_real: mass.min(
                reduced * reduced / (40.0 * PI * PI * n)
                    + (real_ratio / (n * n)) / (6.0 * PI * PI * n),
            ),
            tail_imaginary: (0.5 * mass).min(reduced / (6.0 * PI * PI * n)),
        };
        if !finite(response.mean)
            || !finite(response.wall)
            || !ratio.is_finite()
            || !real_ratio.is_finite()
            || response.mean.re <= 0.0
            || response.mean.im <= 0.0
            || response.wall.re <= 0.0
        {
            return Err(bad("rectangular diffusion response is unrepresentable"));
        }
        Ok(response)
    }
}

fn admit_tail(response: RectangularDiffusion, tolerance: f64) -> Result<(), DuctError> {
    if response.tail_real > tolerance * response.mean.re
        || response.tail_imaginary > tolerance * response.mean.im
        || response.tail_real > tolerance * (response.wall.re - response.tail_real)
    {
        Err(bad(
            "rectangular Fourier tail is unresolved; increase the declared terms",
        ))
    } else {
        Ok(())
    }
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn poll(cancelled: &mut impl FnMut() -> bool, completed_steps: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps })
    } else {
        Ok(())
    }
}
