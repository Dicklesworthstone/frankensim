//! Smooth plane-wave Webster propagation. Geometry is data; local losses
//! come exclusively from `segment_wave`. Fourth-order two-node Magnus
//! integration acts on [p, Z_ref U], so no scalar impedance pole is crossed.
//! This converges to the plane-wave model, not the three-dimensional field.

use crate::{DuctError, LossModel, mul2, segment_wave};
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_phs::WallPin;
use crate::WallModel;

/// Continuous primitive gas fields on the same plane-wave propagation owner.
pub mod gas_field;

/// Named admission caps, never silent resolution clamps.
pub const MAX_RICCATI_STEPS: usize = 65_536;
/// Largest admitted spherical-slice axis; larger requests refuse.
pub const MAX_CONE_SLICES: usize = 65_536;
/// Largest admitted measured-station data set.
pub const MAX_PROFILE_STATIONS: usize = 4_096;

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Fixed, refinable integration resolution; admission occurs at the call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RiccatiSteps(pub usize);
impl RiccatiSteps {
    /// Admit a nonzero resolution without silently reducing it.
    pub fn validate(self) -> Result<(), DuctError> {
        if self.0 == 0 || self.0 > MAX_RICCATI_STEPS {
            return Err(bad("Riccati steps outside 1..=MAX_RICCATI_STEPS"));
        }
        Ok(())
    }
}

/// Lossy-cone resolution. Auto retains the old radius-change rule with
/// its minimum of two slices, and refuses rather than clamping at the cap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConeSlices {
    /// Existing radius-change heuristic, with a named upper-cap refusal.
    Auto,
    /// Exactly the requested number of spherical slices.
    Fixed(usize),
}
impl ConeSlices {
    /// Resolve the axis for positive finite endpoint radii [m].
    pub fn count(self, inlet: f64, outlet: f64) -> Result<usize, DuctError> {
        if !positive(inlet) || !positive(outlet) {
            return Err(bad("cone radii must be positive and finite"));
        }
        let n = match self {
            Self::Auto => {
                let n = ((outlet - inlet).abs() / (0.2 * inlet.min(outlet).max(1.0e-9)))
                    .ceil()
                    .max(2.0);
                if !n.is_finite() || n > MAX_CONE_SLICES as f64 {
                    return Err(bad("automatic cone resolution exceeds MAX_CONE_SLICES"));
                }
                n as usize
            }
            Self::Fixed(n) => n,
        };
        if n == 0 || n > MAX_CONE_SLICES {
            return Err(bad("cone slices outside 1..=MAX_CONE_SLICES"));
        }
        Ok(n)
    }
}

#[derive(Debug, Clone, PartialEq)]
enum Shape {
    Exponential {
        r0: f64,
        area_flare: f64,
    },
    Power {
        r0: f64,
        x0: f64,
        exponent: f64,
    },
    Spline {
        stations: Vec<(f64, f64)>,
        second: Vec<f64>,
    },
}

/// Validated positive radius on a finite local axial interval [0,length].
/// Spline stations retain their original SI abscissae and values; an
/// upstream measured-data owner retains source/uncertainty authority.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiusProfile {
    shape: Shape,
    length: f64,
}

/// Geometric measures of one admitted axial interval of a radius profile.
/// The lateral area follows the actual wall slope and excludes end caps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileGeometry {
    /// Enclosed gas volume [m3].
    pub volume_m3: f64,
    /// Gas-facing lateral wall area [m2].
    pub lateral_area_m2: f64,
    /// Actual two-node Gauss panels, including splits at measured spline knots.
    pub integration_steps: usize,
}

impl RadiusProfile {
    /// r(x)=r0 exp(area_flare*x/2); S(x)=S0 exp(area_flare*x).
    pub fn exponential(r0: f64, area_flare: f64, length: f64) -> Result<Self, DuctError> {
        if !positive(r0)
            || !positive(length)
            || !area_flare.is_finite()
            || !positive(r0 * det::exp(0.5 * area_flare * length))
        {
            return Err(bad("invalid exponential radius profile"));
        }
        Ok(Self {
            shape: Shape::Exponential { r0, area_flare },
            length,
        })
    }

    /// r(x)=r0 ((x0+x)/x0)^exponent, x0>0. Exponents 0,1,2 are
    /// the cylinder, linear cone and elementary order-3/2 Bessel horn.
    pub fn power(r0: f64, x0: f64, exponent: f64, length: f64) -> Result<Self, DuctError> {
        if !positive(r0)
            || !positive(x0)
            || !positive(length)
            || !exponent.is_finite()
            || !positive(r0 * det::exp(exponent * det::ln(1.0 + length / x0)))
        {
            return Err(bad("invalid power-law radius profile"));
        }
        Ok(Self {
            shape: Shape::Power { r0, x0, exponent },
            length,
        })
    }

    /// C2 cubic through SI (axial position,radius) stations. None selects
    /// natural ends; Some([left,right]) supplies dr/dx at both ends.
    /// Positivity is admitted by the cubic Bernstein convex hull on every
    /// span; data failing that sufficient condition refuse, never clip.
    pub fn spline(stations: Vec<(f64, f64)>, slopes: Option<[f64; 2]>) -> Result<Self, DuctError> {
        let n = stations.len();
        if !(2..=MAX_PROFILE_STATIONS).contains(&n)
            || stations
                .iter()
                .any(|&(x, r)| !x.is_finite() || !positive(r))
            || stations.windows(2).any(|w| !positive(w[1].0 - w[0].0))
            || slopes.is_some_and(|s| s.iter().any(|v| !v.is_finite()))
        {
            return Err(bad(
                "invalid spline stations, endpoint slopes or station cap",
            ));
        }
        let length = stations[n - 1].0 - stations[0].0;
        if !positive(length) {
            return Err(bad("invalid spline length"));
        }
        let mut lower = vec![0.0; n];
        let mut diag = vec![1.0; n];
        let mut upper = vec![0.0; n];
        let mut second = vec![0.0; n];
        for i in 1..n - 1 {
            let hl = stations[i].0 - stations[i - 1].0;
            let hr = stations[i + 1].0 - stations[i].0;
            lower[i] = hl;
            diag[i] = 2.0 * (hl + hr);
            upper[i] = hr;
            second[i] = 6.0
                * ((stations[i + 1].1 - stations[i].1) / hr
                    - (stations[i].1 - stations[i - 1].1) / hl);
        }
        if let Some([left, right]) = slopes {
            let hl = stations[1].0 - stations[0].0;
            let hr = stations[n - 1].0 - stations[n - 2].0;
            diag[0] = 2.0 * hl;
            upper[0] = hl;
            second[0] = 6.0 * ((stations[1].1 - stations[0].1) / hl - left);
            lower[n - 1] = hr;
            diag[n - 1] = 2.0 * hr;
            second[n - 1] = 6.0 * (right - (stations[n - 1].1 - stations[n - 2].1) / hr);
        }
        for i in 1..n {
            let factor = lower[i] / diag[i - 1];
            diag[i] -= factor * upper[i - 1];
            second[i] -= factor * second[i - 1];
        }
        for i in (0..n).rev() {
            if i + 1 < n {
                second[i] -= upper[i] * second[i + 1];
            }
            second[i] /= diag[i];
        }
        if second.iter().any(|v| !v.is_finite()) {
            return Err(bad("spline solve is not finite"));
        }
        for i in 0..n - 1 {
            let h = stations[i + 1].0 - stations[i].0;
            let secant = (stations[i + 1].1 - stations[i].1) / h;
            let dl = secant - h * (2.0 * second[i] + second[i + 1]) / 6.0;
            let dr = secant + h * (second[i] + 2.0 * second[i + 1]) / 6.0;
            if !positive(stations[i].1 + h * dl / 3.0)
                || !positive(stations[i + 1].1 - h * dr / 3.0)
            {
                return Err(bad(
                    "spline radius positivity not enclosed by its cubic control hull",
                ));
            }
        }
        Ok(Self {
            shape: Shape::Spline { stations, second },
            length,
        })
    }

    pub fn length(&self) -> f64 {
        self.length
    }
    pub fn inlet_radius(&self) -> f64 {
        self.eval(0.0)
    }
    pub fn outlet_radius(&self) -> f64 {
        self.eval(self.length)
    }
    pub fn stations(&self) -> Option<&[(f64, f64)]> {
        match &self.shape {
            Shape::Spline { stations, .. } => Some(stations),
            _ => None,
        }
    }
    pub fn radius(&self, x: f64) -> Result<f64, DuctError> {
        if !x.is_finite() || !(0.0..=self.length).contains(&x) {
            return Err(bad("profile coordinate outside its axial support"));
        }
        let r = self.eval(x);
        if !positive(r) {
            return Err(bad("profile radius evaluation is not positive and finite"));
        }
        Ok(r)
    }
    /// dr/dx on the axial support. Every shape is at least C1 (the spline is
    /// C2), so the continuous multimodal coupling `R'/R` is defined everywhere.
    pub fn slope(&self, x: f64) -> Result<f64, DuctError> {
        if !x.is_finite() || !(0.0..=self.length).contains(&x) {
            return Err(bad("profile coordinate outside its axial support"));
        }
        let s = match &self.shape {
            Shape::Exponential { area_flare, .. } => 0.5 * area_flare * self.eval(x),
            Shape::Power { x0, exponent, .. } => exponent * self.eval(x) / (x0 + x),
            Shape::Spline { stations, second } => {
                let x = x + stations[0].0;
                let i = stations
                    .partition_point(|&(xi, _)| xi <= x)
                    .saturating_sub(1)
                    .min(stations.len() - 2);
                let h = stations[i + 1].0 - stations[i].0;
                let b = (x - stations[i].0) / h;
                let a = 1.0 - b;
                (stations[i + 1].1 - stations[i].1) / h
                    + ((1.0 - 3.0 * a * a) * second[i] + (3.0 * b * b - 1.0) * second[i + 1]) * h
                        / 6.0
            }
        };
        if !s.is_finite() {
            return Err(bad("profile slope evaluation is not finite"));
        }
        Ok(s)
    }

    /// Integrate volume and lateral wall area over `start_m..end_m`.
    ///
    /// Composite two-node Gauss quadrature is fourth order on each smooth span.
    /// The explicit step count applies to EACH intersected spline span, with a
    /// total cap of `MAX_RICCATI_STEPS`; no knot is skipped or decimated. The
    /// caller refines this axis to control geometric integration error, which is
    /// distinct from thermal-cell and acoustic-field discretization error.
    pub fn geometry(
        &self, start_m: f64, end_m: f64, steps: RiccatiSteps,
    ) -> Result<ProfileGeometry, DuctError> {
        self.geometry_controlled(start_m, end_m, steps, &mut || false)
    }

    /// The same bounded geometry integral with cancellation at each panel.
    /// Refusal or cancellation publishes no partial measures.
    pub fn geometry_controlled(
        &self, start_m: f64, end_m: f64, steps: RiccatiSteps,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ProfileGeometry, DuctError> {
        steps.validate()?;
        if !start_m.is_finite() || !end_m.is_finite()
            || start_m < 0.0 || end_m > self.length || end_m <= start_m
        {
            return Err(bad("profile geometry interval outside its axial support"));
        }
        let mut cuts = vec![start_m];
        if let Some(stations) = self.stations() {
            for &(position, _) in stations {
                let x = position - stations[0].0;
                if x > start_m && x < end_m {
                    cuts.push(x);
                }
            }
        }
        cuts.push(end_m);
        let count = (cuts.len() - 1).checked_mul(steps.0)
            .filter(|count| *count <= MAX_RICCATI_STEPS)
            .ok_or_else(|| bad("profile geometry exceeds MAX_RICCATI_STEPS"))?;
        let root3 = det::sqrt(3.0);
        let nodes = [0.5 - root3 / 6.0, 0.5 + root3 / 6.0];
        let mut result = ProfileGeometry {
            volume_m3: 0.0, lateral_area_m2: 0.0, integration_steps: 0,
        };
        for span in cuts.windows(2) {
            let dx = (span[1] - span[0]) / steps.0 as f64;
            for index in 0..steps.0 {
                if cancelled() {
                    return Err(DuctError::Cancelled { completed_steps: result.integration_steps });
                }
                let left = span[0] + index as f64 * dx;
                let right = if index + 1 == steps.0 { span[1] }
                    else { span[0] + (index + 1) as f64 * dx };
                let width = right - left;
                if !positive(width) {
                    return Err(bad("profile geometry panel is numerically unresolved"));
                }
                for node in nodes {
                    let x = left + node * width;
                    let radius = self.radius(x)?;
                    let slope = self.slope(x)?;
                    let half_width = 0.5 * width;
                    result.volume_m3 += core::f64::consts::PI * radius * radius * half_width;
                    result.lateral_area_m2 += core::f64::consts::TAU * radius
                        * slope.hypot(1.0) * half_width;
                }
                result.integration_steps += 1;
            }
        }
        debug_assert_eq!(result.integration_steps, count);
        if !positive(result.volume_m3) || !positive(result.lateral_area_m2) {
            return Err(bad("profile geometric measures are not positive and finite"));
        }
        Ok(result)
    }

    pub(super) fn eval(&self, x: f64) -> f64 {
        match &self.shape {
            Shape::Exponential { r0, area_flare } => r0 * det::exp(0.5 * area_flare * x),
            Shape::Power { r0, x0, exponent } => r0 * det::exp(exponent * det::ln(1.0 + x / x0)),
            Shape::Spline { stations, second } => {
                let x = x + stations[0].0;
                let i = stations
                    .partition_point(|&(xi, _)| xi <= x)
                    .saturating_sub(1)
                    .min(stations.len() - 2);
                let h = stations[i + 1].0 - stations[i].0;
                let b = (x - stations[i].0) / h;
                let a = 1.0 - b;
                a * stations[i].1
                    + b * stations[i + 1].1
                    + ((a * a * a - a) * second[i] + (b * b * b - b) * second[i + 1]) * h * h / 6.0
            }
        }
    }
}

/// Physical [p_in,U_in]=matrix*[p_out,U_out] and the actual axis used.
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileTransfer {
    /// Physical row-major ABCD two-port in pressure/volume-flow units.
    pub matrix: [C64; 4],
    /// Minimum shear number at the material quadrature points.
    pub min_shear_number: f64,
    /// Actual integration steps, or spherical slices for cone_matrix.
    pub steps_used: usize,
}
/// Mouth-load response in the exp(-i omega t) convention.
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileResponse {
    /// Inlet pressure / inlet volume flow [Pa s/m3].
    pub z_in: C64,
    /// Mouth pressure / inlet volume flow [Pa s/m3].
    pub p_out_over_u_in: C64,
    /// Mouth volume flow / inlet volume flow [dimensionless].
    pub u_out_over_u_in: C64,
    /// Minimum shear number at the material quadrature points.
    pub min_shear_number: f64,
    /// Actual integration steps.
    pub steps_used: usize,
}

/// exp of the traceless 2x2 Magnus generator. Entire series avoid
/// branch choices and the removable sinh(delta)/delta singularity.
fn exp_traceless(d: C64, a: C64, b: C64) -> Result<[C64; 4], DuctError> {
    let q = d * d + a * b;
    let (mut c, mut s, mut tc, mut ts) = (C64::ONE, C64::ONE, C64::ONE, C64::ONE);
    for j in 1..=64_u32 {
        let n = f64::from(2 * j);
        tc = tc * q.scale(1.0 / ((n - 1.0) * n));
        ts = ts * q.scale(1.0 / (n * (n + 1.0)));
        c = c + tc;
        s = s + ts;
        if finite(c)
            && finite(s)
            && tc.abs() <= f64::EPSILON * c.abs()
            && ts.abs() <= f64::EPSILON * s.abs()
        {
            return Ok([c + s * d, s * a, s * b, c - s * d]);
        }
    }
    Err(bad("Magnus exponential exhausted; refine Riccati steps"))
}

// Shared local constitutive and integration owners for both uniform gas and
// source-backed fields. The latter changes where the gas is sampled, not the
// acoustic equation, loss law or Magnus exponential.
fn local_generator(
    state: &GasState,
    radius: f64,
    omega: f64,
    loss: LossModel,
    zref: f64,
    wall: WallModel<'_>,
    slant: f64,
) -> Result<((C64, C64), f64), DuctError> {
    let wave = wall.apply(segment_wave(state, radius, omega, loss)?, radius, omega, slant)?;
    let ik = C64::new(0.0, 1.0) * wave.wavenumber;
    Ok((
        (
            ik * wave.characteristic_impedance.scale(1.0 / zref),
            ik * wave.characteristic_impedance.recip().scale(zref),
        ),
        wave.shear_number,
    ))
}

pub(crate) fn magnus_step(
    (a1, b1): (C64, C64),
    (a2, b2): (C64, C64),
    h: f64,
    root3: f64,
) -> Result<[C64; 4], DuctError> {
    let d = (a1 * b2 - a2 * b1).scale(-root3 * h * h / 12.0);
    exp_traceless(d, (a1 + a2).scale(h / 2.0), (b1 + b2).scale(h / 2.0))
}

/// Transfer through validated geometry, with material and resolution explicit.
pub fn profile_matrix(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    steps: RiccatiSteps,
) -> Result<ProfileTransfer, DuctError> {
    profile_matrix_controlled(profile, state, omega, loss, steps, &mut || false)
}

/// Cancellation is polled once per step, before its two local material
/// evaluations. A cancellation refuses without publishing a partial port.
pub fn profile_matrix_controlled<F: FnMut() -> bool>(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    steps: RiccatiSteps,
    cancel: &mut F,
) -> Result<ProfileTransfer, DuctError> {
    profile_matrix_wall_controlled(profile, state, omega, loss, steps, None, cancel)
}

/// Smooth-profile propagation with the existing locally reacting wall law.
/// The supplied per-area mass, stiffness and resistance remain constant; the
/// actual radius and surface slope set its acoustic shunt at every quadrature
/// point. This is the same wall model as cylinder/cone chains, not a shell solve.
pub fn profile_matrix_wall(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    steps: RiccatiSteps,
    wall: Option<&WallPin>,
) -> Result<ProfileTransfer, DuctError> {
    profile_matrix_wall_controlled(profile, state, omega, loss, steps, wall, &mut || false)
}

/// Wall-aware propagation with the same cancellation and publication contract.
pub fn profile_matrix_wall_controlled<F: FnMut() -> bool>(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    steps: RiccatiSteps,
    wall: Option<&WallPin>,
    cancel: &mut F,
) -> Result<ProfileTransfer, DuctError> {
    profile_matrix_with_wall_controlled(profile, state, omega, loss, steps,
        WallModel::from_pin(wall), cancel)
}

/// Smooth geometry and a source-owned wall law at the same quadrature radii.
pub fn profile_matrix_with_wall(
    profile: &RadiusProfile, state: &GasState, omega: f64, loss: LossModel,
    steps: RiccatiSteps, wall: WallModel<'_>,
) -> Result<ProfileTransfer, DuctError> {
    profile_matrix_with_wall_controlled(profile, state, omega, loss, steps, wall, &mut || false)
}

/// Complete wall response with the existing cancellation/publication contract.
pub fn profile_matrix_with_wall_controlled<F: FnMut() -> bool>(
    profile: &RadiusProfile, state: &GasState, omega: f64, loss: LossModel,
    steps: RiccatiSteps, wall: WallModel<'_>, cancel: &mut F,
) -> Result<ProfileTransfer, DuctError> {
    let mut m = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    let (min_shear, zref) =
        for_each_profile_step(profile, state, omega, loss, steps, wall, cancel, |e, _| {
            m = mul2(e, m);
            if m.iter().any(|&v| !finite(v)) {
                Err(DuctError::Singular)
            } else {
                Ok(())
            }
        })?;
    m[1] = m[1].scale(zref);
    m[2] = m[2].scale(1.0 / zref);
    if m.iter().any(|&v| !finite(v)) {
        return Err(DuctError::Singular);
    }
    Ok(ProfileTransfer {
        matrix: m,
        min_shear_number: min_shear,
        steps_used: steps.0,
    })
}

/// Same local Magnus owner, exposing resolved physical two-ports in inlet to
/// outlet order. Stable composition never has to form the full profile ABCD.
pub(crate) fn profile_blocks_controlled<F: FnMut() -> bool>(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    steps: RiccatiSteps,
    wall: WallModel<'_>,
    cancel: &mut F,
) -> Result<(Vec<[C64; 4]>, f64), DuctError> {
    steps.validate()?;
    let mut blocks = Vec::with_capacity(steps.0);
    let (rv, _) =
        for_each_profile_step(profile, state, omega, loss, steps, wall, cancel, |mut e, zref| {
            e[1] = e[1].scale(zref);
            e[2] = e[2].scale(1.0 / zref);
            if e.iter().any(|&v| !finite(v)) {
                return Err(DuctError::Singular);
            }
            blocks.push(e);
            Ok(())
        })?;
    blocks.reverse();
    Ok((blocks, rv))
}

#[allow(clippy::too_many_arguments)]
fn for_each_profile_step<F: FnMut() -> bool>(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    steps: RiccatiSteps,
    wall: WallModel<'_>,
    cancel: &mut F,
    mut visit: impl FnMut([C64; 4], f64) -> Result<(), DuctError>,
) -> Result<(f64, f64), DuctError> {
    steps.validate()?;
    if !positive(omega) || !positive(state.density) || !positive(state.sound_speed) {
        return Err(bad("profile gas/frequency must be positive and finite"));
    }
    let r0 = profile.inlet_radius();
    let zref = state.density * state.sound_speed / (core::f64::consts::PI * r0 * r0);
    if !positive(zref) {
        return Err(bad("profile reference impedance is not finite"));
    }
    let h = -profile.length / steps.0 as f64;
    let root3 = det::sqrt(3.0);
    let (c1, c2) = (0.5 - root3 / 6.0, 0.5 + root3 / 6.0);
    let mut min_shear = f64::INFINITY;
    let mut local = |x: f64| -> Result<(C64, C64), DuctError> {
        let slant = if !wall.is_rigid() {
            let slope = profile.slope(x)?;
            det::sqrt(1.0 + slope * slope)
        } else { 1.0 };
        let (generator, shear) = local_generator(
            state, profile.radius(x)?, omega, loss, zref, wall, slant,
        )?;
        min_shear = min_shear.min(shear);
        Ok(generator)
    };
    for j in 0..steps.0 {
        if cancel() {
            return Err(DuctError::Cancelled { completed_steps: j });
        }
        let x = profile.length + j as f64 * h;
        let (a1, b1) = local(x + c1 * h)?;
        let (a2, b2) = local(x + c2 * h)?;
        let e = magnus_step((a1, b1), (a2, b2), h, root3)?;
        visit(e, zref)?;
    }
    Ok((min_shear, zref))
}

/// Finite p/U mouth load, or (+infinity,0) for a rigid cap. Zero is
/// ideal pressure release. Ratios carry pressure/flow units explicitly.
pub fn profile_impedance(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    z_load: C64,
    steps: RiccatiSteps,
) -> Result<ProfileResponse, DuctError> {
    let closed = z_load.re == f64::INFINITY && z_load.im == 0.0;
    if !closed && !finite(z_load) {
        return Err(bad("invalid profile mouth impedance"));
    }
    let t = profile_matrix(profile, state, omega, loss, steps)?;
    let [a, b, c, d] = t.matrix;
    let (pout, uout) = if closed {
        (C64::ONE, C64::ZERO)
    } else {
        (z_load, C64::ONE)
    };
    let pin = a * pout + b * uout;
    let uin = c * pout + d * uout;
    if uin.abs() == 0.0 {
        return Err(DuctError::Singular);
    }
    let out = ProfileResponse {
        z_in: pin / uin,
        p_out_over_u_in: pout / uin,
        u_out_over_u_in: uout / uin,
        min_shear_number: t.min_shear_number,
        steps_used: t.steps_used,
    };
    if !finite(out.z_in) || !finite(out.p_out_over_u_in) || !finite(out.u_out_over_u_in) {
        return Err(DuctError::Singular);
    }
    Ok(out)
}

/// Refinable existing spherical-cone slices, sharing the same local ZK
/// owner as smooth propagation. No second constitutive loss formula.
pub fn cone_matrix(
    inlet: f64,
    outlet: f64,
    length: f64,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    slices: ConeSlices,
) -> Result<ProfileTransfer, DuctError> {
    if !positive(length) {
        return Err(bad("cone length must be positive and finite"));
    }
    let n = slices.count(inlet, outlet)?;
    let (matrix, min_shear_number) =
        crate::cone_lossy_matrix_axis(inlet, outlet, length, state, omega, loss, None, slices)?;
    Ok(ProfileTransfer {
        matrix,
        min_shear_number,
        steps_used: n,
    })
}

#[cfg(test)]
mod tests;
