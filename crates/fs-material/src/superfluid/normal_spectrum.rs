//! Passive transverse Stokes memory of the normal component at fixed walls.
//!
//! For uniform axial acceleration `g_n`, each normalized velocity obeys
//! `dU_m/dt + lambda_m*U_m = g_n`. The section mean is `sum(w_m*U_m)`,
//! normal kinetic energy per fluid volume is `rho_n/2 * sum(w_m*U_m^2)`,
//! and viscous heating is `rho_n * sum(w_m*lambda_m*U_m^2)`. Positive
//! weights and rates retain passivity. Superfluid velocity has no wall drag.
//!
//! Circular modes use the positive zeros `j_m` of **J0**, with
//! `w_m=4/j_m^2`, `lambda_m=(eta/rho_n)*j_m^2/R^2`. Parallel plates use
//! the FULL gap `H`: `w_m=8/((2m+1)*pi)^2` and
//! `lambda_m=(eta/rho_n)*((2m+1)*pi/H)^2`, starting at `m=0`.
//! A finite prefix is followed by one positive aggregate tail. It matches
//! both omitted weight and static mobility, so total normal inertia and the
//! exact Poiseuille resistance survive truncation. The aggregate tail is an
//! approximation; [`NormalViscousSpectrum::harmonic_mobility_error_bound_s`]
//! bounds its frequency-response error against the infinite Stokes spectrum.
//!
//! This linear, fully developed, continuum/no-slip model resolves unsteady
//! transverse momentum diffusion without requiring a small Womersley number.
//! Material properties and geometry remain fixed. Entrance effects, transverse
//! thermodynamics, moving walls, nonlinear normal turbulence, and pointwise
//! vortex drag across the profile are outside this closure. Initializing all
//! `U_m` to the same velocity represents uniform initial normal flow, not an
//! already developed parabolic profile. A spatial owner must keep the full
//! modal kinetic energy, external wall impulse, and heat exactly once.
//!
//! Primary references: P. Szymanski, J. Math. Pures Appl. 11 (1932), 67-107,
//! <https://www.numdam.org/item/JMPA_1932_9_11__67_0.pdf>, derives the
//! cylindrical Fourier-Bessel expansion. J. R. Womersley, J. Physiol. 127
//! (1955), 553-563, DOI: 10.1113/jphysiol.1955.sp005276, treats oscillatory
//! tube flow and viscous drag. The planar eigenfunction solution is given in
//! *Unified Theory of Unsteady Planar Laminar Flow in the Presence of Arbitrary
//! Pressure Gradients and Boundary Movement*, DOI: 10.3390/sym14040757.
//! Root conventions and Rayleigh sums are documented in
//! <https://dlmf.nist.gov/10.21>. The positive moment-matched tail and its bound
//! below are a numerical reduction of these spectra, not a measured helium fit.

use super::{NormalViscousWall, TwoFluidState, invalid, nonnegative, positive};
use crate::MaterialError;
use fs_math::{c64::C64, det};

/// Maximum number of exact geometry modes retained before the aggregate tail.
pub const MAX_RESOLVED_NORMAL_VISCOUS_MODES: usize = 32;
/// Fixed storage capacity, including the positive aggregate tail.
pub const MAX_NORMAL_VISCOUS_MODES: usize = MAX_RESOLVED_NORMAL_VISCOUS_MODES + 1;

/// A normalized transverse velocity carries `weight` of the normal inertia.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalViscousMode {
    pub weight: f64,
    /// Decay rate [1/s]; zero only for the admitted zero-viscosity limit.
    pub decay_rate_s_inv: f64,
}

impl NormalViscousMode {
    const EMPTY: Self = Self {
        weight: 0.0,
        decay_rate_s_inv: 0.0,
    };
}

/// Explicit information about the unresolved infinite tail.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalViscousTruncation {
    pub omitted_weight: f64,
    /// Omitted `sum(w/lambda)` divided by the full DC mobility `k*rho_n/eta`.
    /// This dimensionless moment remains defined in the zero-viscosity limit.
    pub omitted_dc_mobility_fraction: f64,
    pub first_omitted_decay_rate_s_inv: f64,
    pub modeled_tail_decay_rate_s_inv: f64,
}

/// Frozen geometry, normal inertia and viscosity; fixed storage, no allocation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalViscousSpectrum {
    material: TwoFluidState,
    wall: NormalViscousWall,
    modes: [NormalViscousMode; MAX_NORMAL_VISCOUS_MODES],
    resolved_modes: usize,
    darcy_decay_rate_s_inv: f64,
    truncation: NormalViscousTruncation,
}

impl NormalViscousSpectrum {
    /// Circular capillary with stationary no-slip wall and radius in metres.
    /// Retains `1..=32` exact J0 modes plus one moment-matched aggregate tail.
    pub fn circular_capillary(
        material: TwoFluidState,
        radius_m: f64,
        resolved_modes: usize,
    ) -> Result<Self, MaterialError> {
        Self::build(
            material,
            NormalViscousWall::circular_capillary(radius_m)?,
            resolved_modes,
            true,
        )
    }

    /// Wide parallel-plate slit with stationary walls and FULL gap in metres.
    /// Retains `1..=32` exact odd-sine modes plus one aggregate tail.
    pub fn parallel_plates(
        material: TwoFluidState,
        full_gap_m: f64,
        resolved_modes: usize,
    ) -> Result<Self, MaterialError> {
        Self::build(
            material,
            NormalViscousWall::parallel_plates(full_gap_m)?,
            resolved_modes,
            false,
        )
    }

    fn build(
        material: TwoFluidState,
        wall: NormalViscousWall,
        resolved_modes: usize,
        circle: bool,
    ) -> Result<Self, MaterialError> {
        if !(1..=MAX_RESOLVED_NORMAL_VISCOUS_MODES).contains(&resolved_modes) {
            return Err(invalid(
                "normal Stokes spectrum requires 1 through 32 resolved modes",
            ));
        }
        let darcy_rate = material.normal_wall_decay_rate_s_inv(wall)?;
        let mut modes = [NormalViscousMode::EMPTY; MAX_NORMAL_VISCOUS_MODES];
        let mut weight_sum = 0.0;
        let mut mobility_sum = 0.0;
        let mut mobility_correction = 0.0;
        for (index, mode) in modes[..resolved_modes].iter_mut().enumerate() {
            let (weight, relative_rate) = geometry_mode(circle, index);
            let rate = checked_rate(darcy_rate, relative_rate)?;
            *mode = NormalViscousMode {
                weight,
                decay_rate_s_inv: rate,
            };
            // The tail uses precisely this ordinary sum: reducing the entire
            // published weight slice in order then recovers one exactly.
            weight_sum += weight;
            add(
                &mut mobility_sum,
                &mut mobility_correction,
                weight / relative_rate,
            );
        }
        let tail_weight = 1.0 - weight_sum;
        let tail_mobility = 1.0 - mobility_sum;
        if !positive(tail_weight) || !positive(tail_mobility) {
            return Err(invalid("normal Stokes tail moments are unresolved"));
        }
        let tail_rate = checked_rate(darcy_rate, tail_weight / tail_mobility)?;
        let first_omitted_rate = checked_rate(darcy_rate, geometry_mode(circle, resolved_modes).1)?;
        if tail_rate < first_omitted_rate {
            return Err(invalid(
                "normal Stokes tail rate violates the omitted spectral support",
            ));
        }
        modes[resolved_modes] = NormalViscousMode {
            weight: tail_weight,
            decay_rate_s_inv: tail_rate,
        };
        Ok(Self {
            material,
            wall,
            modes,
            resolved_modes,
            darcy_decay_rate_s_inv: darcy_rate,
            truncation: NormalViscousTruncation {
                omitted_weight: tail_weight,
                omitted_dc_mobility_fraction: tail_mobility,
                first_omitted_decay_rate_s_inv: first_omitted_rate,
                modeled_tail_decay_rate_s_inv: tail_rate,
            },
        })
    }

    pub const fn material(&self) -> &TwoFluidState {
        &self.material
    }
    pub const fn wall(&self) -> NormalViscousWall {
        self.wall
    }
    pub const fn resolved_mode_count(&self) -> usize {
        self.resolved_modes
    }
    pub fn modes(&self) -> &[NormalViscousMode] {
        &self.modes[..=self.resolved_modes]
    }
    pub const fn truncation(&self) -> NormalViscousTruncation {
        self.truncation
    }

    /// Static mean velocity divided by normal acceleration [s].
    /// Zero viscosity has no finite steady wall-balanced response.
    pub fn dc_mobility_s(&self) -> Result<f64, MaterialError> {
        let mobility = self.darcy_decay_rate_s_inv.recip();
        if !positive(mobility) {
            return Err(invalid(
                "normal Stokes DC mobility is not finite and positive",
            ));
        }
        Ok(mobility)
    }

    /// Mean velocity / uniform acceleration for the `exp(+i*omega*t)`
    /// convention, in seconds. The imaginary part is nonpositive (phase lag).
    /// `omega` is a finite nonnegative angular frequency in radians/second.
    pub fn harmonic_mobility_s(&self, omega_rad_s: f64) -> Result<C64, MaterialError> {
        if !nonnegative(omega_rad_s) {
            return Err(invalid(
                "normal Stokes angular frequency must be nonnegative finite",
            ));
        }
        if omega_rad_s == 0.0 {
            return self.dc_mobility_s().map(C64::from_re);
        }
        let mut real = 0.0;
        let mut imag = 0.0;
        let mut real_correction = 0.0;
        let mut imag_correction = 0.0;
        for mode in self.modes() {
            // Scale by the larger frequency instead of squaring either one.
            let (re, im) = if mode.decay_rate_s_inv >= omega_rad_s {
                let ratio = omega_rad_s / mode.decay_rate_s_inv;
                let re = (mode.weight / mode.decay_rate_s_inv) / (1.0 + ratio * ratio);
                (re, -ratio * re)
            } else {
                let ratio = mode.decay_rate_s_inv / omega_rad_s;
                let im = -(mode.weight / omega_rad_s) / (1.0 + ratio * ratio);
                (-ratio * im, im)
            };
            add(&mut real, &mut real_correction, re);
            add(&mut imag, &mut imag_correction, im);
        }
        if !real.is_finite() || !imag.is_finite() {
            return Err(invalid(
                "normal Stokes harmonic mobility is unrepresentable",
            ));
        }
        Ok(C64::new(real, imag))
    }

    /// Absolute complex-mobility error bound [s] from replacing ONLY the
    /// omitted infinite Stokes tail with its single matched pole.
    ///
    /// Let `W=sum_tail(w)`, `mu=sum_tail(w/lambda)`, `tau_bar=mu/W`
    /// and `tau_max=1/lambda_first_omitted`. Taylor expansion of
    /// `tau/(1+s*tau)` about `tau_bar` cancels the matched first moment.
    /// At `s=i*omega`, the remaining magnitude is at most
    /// `omega*mu*(tau_max-tau_bar)/(1+(omega*tau_bar)^2)`.
    /// The variance inequality uses only `0<=tau<=tau_max`; no subtraction
    /// of tiny high-order Rayleigh tails is needed. Zero frequency is exact.
    ///
    /// This analytical truncation bound does not enclose floating-point
    /// rounding, spatial/time discretization, or departures from the stated
    /// linear fully developed continuum physics.
    pub fn harmonic_mobility_error_bound_s(&self, omega_rad_s: f64) -> Result<f64, MaterialError> {
        if !nonnegative(omega_rad_s) {
            return Err(invalid(
                "normal Stokes angular frequency must be nonnegative finite",
            ));
        }
        if omega_rad_s == 0.0 || self.darcy_decay_rate_s_inv == 0.0 {
            return Ok(0.0);
        }
        let tail = self.truncation;
        let spread = tail.modeled_tail_decay_rate_s_inv / tail.first_omitted_decay_rate_s_inv - 1.0;
        let bound = if omega_rad_s <= tail.modeled_tail_decay_rate_s_inv {
            let ratio = omega_rad_s / tail.modeled_tail_decay_rate_s_inv;
            // Divide only after attenuating the dimensionless numerator so
            // a very large DC mobility need not itself be representable.
            (tail.omitted_weight * spread * ratio / (1.0 + ratio * ratio))
                / tail.modeled_tail_decay_rate_s_inv
        } else {
            let ratio = tail.modeled_tail_decay_rate_s_inv / omega_rad_s;
            (tail.omitted_weight * spread / (1.0 + ratio * ratio)) / omega_rad_s
        };
        if !nonnegative(bound) {
            return Err(invalid(
                "normal Stokes harmonic tail bound is unrepresentable",
            ));
        }
        Ok(bound)
    }

    /// Prepare all modal exponentials once for an allocation-free spatial batch.
    pub fn prepare(&self, duration_s: f64) -> Result<PreparedNormalViscousSpectrum, MaterialError> {
        let mut modes = [PreparedNormalViscousMode {
            weight: 0.0,
            factors: NormalViscousFactors::FREE,
        }; MAX_NORMAL_VISCOUS_MODES];
        // Admit even an empty-duration preparation before filling the batch.
        NormalViscousFactors::try_new(0.0, duration_s)?;
        for (prepared, mode) in modes.iter_mut().zip(self.modes()) {
            *prepared = PreparedNormalViscousMode {
                weight: mode.weight,
                factors: NormalViscousFactors::try_new(mode.decay_rate_s_inv, duration_s)?,
            };
        }
        Ok(PreparedNormalViscousSpectrum {
            modes,
            mode_count: self.resolved_modes + 1,
            duration_s,
        })
    }
}

/// Exponential factors for `z' = -rate*z + g` with constant `g` over `dt`.
///
/// With `w=dt*g`, the final value is `decay*z0+phi1*w`, and the interval
/// mean is `phi1*z0+phi2*w`. Wall impulse per unit mass is
/// `-(loss*z0+drag_forced*w)`. Nonnegative heat per unit mass is the sum
/// of squares `(heat_initial_root*(z0+heat_shift*w))^2 +
/// (heat_forced_root*w)^2`. Thus `delta(z^2/2)=w*mean-heat` without
/// cancellation-prone kinetic-energy subtraction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalViscousFactors {
    pub decay: f64,
    pub loss: f64,
    pub phi1: f64,
    pub phi2: f64,
    pub drag_forced: f64,
    pub heat_initial_root: f64,
    pub heat_shift: f64,
    pub heat_forced_root: f64,
}

impl NormalViscousFactors {
    const FREE: Self = Self {
        decay: 1.0,
        loss: 0.0,
        phi1: 1.0,
        phi2: 0.5,
        drag_forced: 0.0,
        heat_initial_root: 0.0,
        heat_shift: 0.5,
        heat_forced_root: 0.0,
    };

    /// Nonnegative finite rate [1/s] and duration [s]. Zero rate or duration
    /// gives the exact free-motion factors. An overflowed product `rate*dt`
    /// is admitted using reciprocal scaling in the fully stiff limit.
    pub fn try_new(rate_s_inv: f64, duration_s: f64) -> Result<Self, MaterialError> {
        if !nonnegative(rate_s_inv) || !nonnegative(duration_s) {
            return Err(invalid(
                "normal Stokes rate and duration must be nonnegative finite",
            ));
        }
        if rate_s_inv == 0.0 || duration_s == 0.0 {
            return Ok(Self::FREE);
        }
        let s = rate_s_inv * duration_s;
        let decay = det::exp(-s);
        let loss = -det::expm1(-s);
        let (phi1, phi2, drag_forced, heat_shift, heat_forced_root) = if s < 0.125 {
            let mut term = 0.5;
            let mut phi2 = term;
            for denominator in 3..=13 {
                term *= -s / f64::from(denominator);
                phi2 += term;
            }
            let drag_forced = s * phi2;
            let square = s * s;
            let heat_shift = 0.5
                + square
                    * (-1.0 / 24.0
                        + square
                            * (1.0 / 240.0
                                + square
                                    * (-17.0 / 40320.0
                                        + square
                                            * (31.0 / 725760.0 - square * 691.0 / 159667200.0))));
            let heat_forced_per_s = 1.0 / 12.0
                + square
                    * (-1.0 / 120.0
                        + square
                            * (17.0 / 20160.0
                                + square
                                    * (-31.0 / 362880.0
                                        + square
                                            * (691.0 / 79833600.0
                                                - square * 5461.0 / 6227020800.0))));
            // Taking the roots separately retains a representable coefficient
            // even when s/12 would underflow as a subnormal product.
            let root_s = if s > 0.0 {
                s.sqrt()
            } else {
                rate_s_inv.sqrt() * duration_s.sqrt()
            };
            (
                1.0 - drag_forced,
                phi2,
                drag_forced,
                heat_shift,
                root_s * heat_forced_per_s.sqrt(),
            )
        } else {
            let inverse = if s.is_finite() {
                s.recip()
            } else {
                rate_s_inv.recip() / duration_s
            };
            let phi1 = loss * inverse;
            let drag_forced = 1.0 - phi1;
            let heat_shift = det::tanh(0.5 * s) * inverse;
            let inverse_root = if inverse > 0.0 {
                inverse.sqrt()
            } else {
                rate_s_inv.sqrt().recip() / duration_s.sqrt()
            };
            let heat_forced_root = (1.0 - 2.0 * heat_shift).sqrt() * inverse_root;
            (
                phi1,
                drag_forced * inverse,
                drag_forced,
                heat_shift,
                heat_forced_root,
            )
        };
        Ok(Self {
            decay,
            loss,
            phi1,
            phi2,
            drag_forced,
            heat_initial_root: if s > 0.0 {
                (-0.5 * det::expm1(-2.0 * s)).sqrt()
            } else {
                rate_s_inv.sqrt() * duration_s.sqrt()
            },
            heat_shift,
            heat_forced_root,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedNormalViscousMode {
    pub weight: f64,
    pub factors: NormalViscousFactors,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedNormalViscousSpectrum {
    modes: [PreparedNormalViscousMode; MAX_NORMAL_VISCOUS_MODES],
    mode_count: usize,
    duration_s: f64,
}

impl PreparedNormalViscousSpectrum {
    pub fn modes(&self) -> &[PreparedNormalViscousMode] {
        &self.modes[..self.mode_count]
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
}

fn geometry_mode(circle: bool, index: usize) -> (f64, f64) {
    if circle {
        let square = J0_ZEROS[index] * J0_ZEROS[index];
        (4.0 / square, square / 8.0)
    } else {
        let root = (2 * index + 1) as f64 * core::f64::consts::PI;
        let square = root * root;
        (8.0 / square, square / 12.0)
    }
}

fn checked_rate(darcy_rate: f64, multiplier: f64) -> Result<f64, MaterialError> {
    let rate = darcy_rate * multiplier;
    if !nonnegative(rate) || (darcy_rate > 0.0 && rate == 0.0) {
        return Err(invalid("normal Stokes modal decay rate is unrepresentable"));
    }
    Ok(rate)
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let adjusted = value - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

// Positive J0 Dirichlet roots, rounded to binary64 from 110-digit evaluation
// of the defining J0/J1 power series. The 33rd root bounds the omitted tail
// when 32 modes are retained. These are mathematical constants, not an
// empirical geometry fit; tests independently check the J0 integral identity.
const J0_ZEROS: [f64; MAX_NORMAL_VISCOUS_MODES] = [
    2.404825557695773,
    5.520078110286311,
    8.653727912911013,
    11.791534439014281,
    14.930917708487787,
    18.071063967910924,
    21.21163662987926,
    24.352471530749302,
    27.493479132040253,
    30.634606468431976,
    33.77582021357357,
    36.917098353664045,
    40.05842576462824,
    43.19979171317673,
    46.341188371661815,
    49.482609897397815,
    52.624051841115,
    55.76551075501998,
    58.90698392608094,
    62.048469190227166,
    65.18996480020687,
    68.3314693298568,
    71.47298160359374,
    74.61450064370183,
    77.75602563038805,
    80.89755587113763,
    84.0390907769382,
    87.18062984364116,
    90.32217263721049,
    93.46371878194478,
    96.60526795099626,
    99.7468198586806,
    102.8883742541948,
];
