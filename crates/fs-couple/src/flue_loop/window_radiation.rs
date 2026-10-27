//! Radiation load of a flue window on the head's outer cylinder.
//!
//! The opening is the window's outer face on the head's rigid outer cylinder
//! `r = A` (bore radius plus wall), as [`fs_bem::wall_window`] solves it: the
//! face's cells and the flux through each. As a bracket it is one uniform
//! piston over the face's angle and the window's cut-up. Each cell moves as a
//! uniform piston at its own flux. The cylinder is infinite, rigid and in
//! free space, and its gas is the declared exterior.
//!
//! Resistance (Junger and Feit, Sound, Structures, and Their Interaction,
//! ch. 6: the cylinder's modal radiation impedance; `e^{-i omega t}`). With
//! `k_z = k sin(theta)`, `k_r = k cos(theta)`, cell `p` at angle `phi_p` and
//! axial centre `z_p`, spanning `dphi_p` and `l_p`, carrying flux `q_p` of
//! the total `U`:
//! - `F_m(k_z) = sum_p (q_p/U)/(dphi_p l_p) e^{-i m phi_p} s_m(dphi_p)
//!   e^{-i k_z z_p} b(k_z, l_p)`, with `s_m(d) = 2 sin(m d/2)/m` (`s_0 = d`)
//!   and `b(k_z, l) = 2 sin(k_z l/2)/k_z` (`b(0, l) = l`);
//! - `R = omega rho/(pi^3 A^2) int_0^{pi/2} (1/k_r) sum_m |F_m|^2 /
//!   |H_m'(k_r A)|^2 dtheta`, the sum over all integer `m`, folded onto
//!   `m >= 0` by `|F_m|^2 + |F_-m|^2 = 2 (|C_m|^2 + |S_m|^2)`, where `C_m`
//!   and `S_m` carry `cos(m phi_p)` and `sin(m phi_p)`.
//!
//! For a uniform piston of area `S` this is the series
//! `omega rho/(pi^3 S^2) sum eps_m a_m^2 int b^2/(k_r |H_m'|^2) dtheta`. As
//! `kA -> 0` it tends to the free monopole `rho omega^2/(4 pi c)`; at high
//! frequency, to `R_inf = rho c sum_p (q_p/U)^2/S_p`. The series stops at
//! `m = k_r A + 40`; `theta` takes `4 k max(span, A) + 64` Gauss-Legendre
//! nodes, rounded up to a multiple of 64.
//!
//! Static mass (Kramers-Kronig): `M0 = (2/pi) int_0^inf R/omega^2 domega`,
//! by 100 Gauss-Legendre nodes in `ln omega` from `0.01 c/reach` to
//! `omega_top = 20 c/min(a_eq, A)`. Below the bottom the monopole limit
//! holds; above the top `R_inf` is taken, and its share of `M0` is reported.
//!
//! Realization: a series chain of parallel cells ([`ParallelCellLoad`]),
//! passive by construction. The dictionary is R-L poles at `k a_eq = 0.5 *
//! 2^i` up to 256, and R-L-C resonances of `Q = 1.8` at `k a_eq = 0.05 *
//! 1.45^i` up to 10. (A chain of R-L cells alone has `Re Z/omega^2`
//! non-increasing; the cylinder's rises from the free toward the baffled
//! value near `kA ~ 1`, which needs the resonances.) Nonnegative least
//! squares fits the relative resistance at 150 log-spaced frequencies from
//! 20 Hz to `omega_top` or the mechanics Nyquist, whichever is higher:
//! weight 1 to 24 kHz and 0.3 above. Two heavy rows hold the cells' static masses to `M0` and
//! the poles' high-frequency resistance to `R_inf`. Lawson reweighting moves
//! the fit toward minimax. At most 16 cells are kept: while more are active,
//! the resonance whose removal costs least is dropped. Nothing is tuned to a
//! note. A realization missing `R` by more than 5% anywhere from 20 Hz to
//! 24 kHz is refused. The error from 24 kHz to the mechanics Nyquist is
//! reported, as are the realized (bilinear) band error and the `kA` range.
//!
//! The voice's series inertance becomes the window's computed mass less the
//! cells' static mass, so the static mass from junction to far field stays
//! the computed one.
//!
//! A face beside the window ([`fs_bem::wall_window::CoverSide`]) enters
//! through the computed mass only: it closes the static exterior the solve
//! sees, while the resistance and the cells' static mass stay the bare
//! cylinder's, so the series inertance carries the face's static mass.
//!
//! Not represented: the face's effect on the resistance, and a lip plate or
//! baffle near the window; the window's mutual radiation with the foot or
//! the open holes; flow through the opening (the radiation is linear); a
//! finite head length.

use super::{FlueError, FlueVoice};
use core::f64::consts::{FRAC_PI_4, PI, TAU};
use fs_bem::wall_window::OuterFaceCell;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::det;
use fs_vfit::relaxation::{
    MAX_PARALLEL_CELLS, ParallelCell, ParallelCellFrame, ParallelCellLoad, ParallelCellSpec,
};
use std::collections::BTreeMap;

/// The band the realized resistance is checked on [Hz].
pub const WINDOW_RADIATION_BAND_HZ: [f64; 2] = [20.0, 24_000.0];
/// A realization missing the resistance by more than this in the band is refused.
pub const MAX_WINDOW_RADIATION_ERROR: f64 = 0.05;
/// Bounded cold construction: an opening with more cells is refused.
pub const MAX_WINDOW_RADIATION_CELLS: usize = 65_536;

/// Fit samples, Kramers-Kronig nodes, Lawson passes and weights.
const FIT_SAMPLES: usize = 150;
const STATIC_MASS_NODES: usize = 100;
const LAWSON_PASSES: usize = 60;
const PRUNED_PASSES: usize = 40;
const TAIL_WEIGHT: f64 = 0.3;
const EQUALITY_WEIGHT: f64 = 30.0;
const RESONANCE_Q: f64 = 1.8;
/// `k a_eq` at the fit's and the mass integral's top.
const TOP_KA: f64 = 20.0;

fn invalid(what: &'static str) -> FlueError {
    FlueError::Invalid { what }
}
fn physics(error: impl core::fmt::Display) -> FlueError {
    FlueError::Physics(crate::acoustic_realize::AcousticRealizeError::Nonlinear(
        error.to_string(),
    ))
}

/// How the flux spreads over the opening.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WindowFlux {
    /// The wall window's own solved outer-face cells.
    WallWindow,
    /// One uniform piston over the outer face's angle and the window's
    /// cut-up, centred on the opening: a bracket.
    CutUpPiston,
}

impl WindowFlux {
    /// Receipt label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::WallWindow => "wall-window-outer-face-flux",
            Self::CutUpPiston => "uniform-cut-up-piston",
        }
    }
}

/// The uniform piston over the angle the outer `cells` span and the
/// window's `cut_up_m`, centred on that span (on the window when the cells
/// are symmetric), carrying unit flux.
#[must_use]
pub fn cut_up_piston(cells: &[OuterFaceCell], cut_up_m: f64) -> OuterFaceCell {
    let reach = cells
        .iter()
        .map(|cell| cell.theta.abs() + 0.5 * cell.arc)
        .fold(0.0, f64::max);
    // A one-sided opening's span does not straddle the mid-plane.
    let lo = cells
        .iter()
        .map(|cell| cell.theta - 0.5 * cell.arc)
        .fold(f64::INFINITY, f64::min);
    let hi = cells
        .iter()
        .map(|cell| cell.theta + 0.5 * cell.arc)
        .fold(f64::NEG_INFINITY, f64::max);
    let (theta, arc) = if lo.is_finite() && (lo + hi).abs() > 1e-12 * reach {
        (0.5 * (lo + hi), hi - lo)
    } else {
        (0.0, 2.0 * reach)
    };
    OuterFaceCell {
        theta,
        arc,
        axial: 0.0,
        length: cut_up_m,
        flux: 1.0,
    }
}

/// The opening, its passive realization and how well it holds.
#[derive(Clone, Copy, Debug)]
pub struct WindowRadiationLoad {
    flux: WindowFlux,
    outer_radius_m: f64,
    cells: usize,
    open_area_m2: f64,
    exterior: GasState,
    time_step_s: f64,
    static_mass: f64,
    static_mass_tail_share: f64,
    high_frequency_resistance: f64,
    ka_band: [f64; 2],
    band_error: f64,
    discrete_band_error: f64,
    tail_error: f64,
    tail_top_hz: f64,
    load: ParallelCellSpec,
}

impl WindowRadiationLoad {
    /// Fit the opening `cells` on the outer cylinder of `outer_radius_m` in
    /// `exterior`, for a voice stepping at `time_step_s`.
    ///
    /// # Errors
    /// Invalid gas, radius, cells or clock; a clock whose Nyquist does not
    /// exceed the band; cancellation; a realization missing the resistance
    /// by more than 5% in the band; any refusal of the cell owner.
    pub fn new(
        outer_radius_m: f64,
        cells: &[OuterFaceCell],
        flux: WindowFlux,
        exterior: GasState,
        time_step_s: f64,
        gate: &CancelGate,
    ) -> Result<Self, FlueError> {
        exterior.validate_coherence().map_err(physics)?;
        let (rho, c) = (exterior.density, exterior.sound_speed);
        if !(time_step_s.is_finite() && time_step_s > 0.0) {
            return Err(invalid("window radiation requires a positive finite time step"));
        }
        let nyquist_hz = 0.5 / time_step_s;
        if !(nyquist_hz > WINDOW_RADIATION_BAND_HZ[1]) {
            return Err(invalid(
                "window radiation is checked to 24 kHz; the mechanics clock must resolve it",
            ));
        }
        let opening = Opening::new(outer_radius_m, cells)?;
        let a_eq = det::sqrt(opening.area / PI);
        let omega_top = TOP_KA * c / a_eq.min(outer_radius_m);
        let high_frequency_resistance = rho * c * opening.concentration;
        let mut quadrature = Quadrature::default();
        let mut resistance = |omega: f64| {
            if gate.is_requested() {
                return Err(invalid("window radiation construction cancelled"));
            }
            opening.resistance(omega, rho, c, &mut quadrature)
        };

        // Kramers-Kronig static mass.
        let reach = outer_radius_m
            .max(outer_radius_m * opening.angular_span)
            .max(opening.axial_span);
        let omega_lo = 0.01 * c / reach;
        let span = det::ln(omega_top / omega_lo);
        let (nodes, weights) = fs_feec::highorder::quad1d::gauss_legendre(STATIC_MASS_NODES);
        let mut integral = 0.0;
        for (&t, &w) in nodes.iter().zip(&weights) {
            let omega = det::exp(det::ln(omega_lo) + 0.5 * (t + 1.0) * span);
            integral += 0.5 * w * span * resistance(omega)? / omega;
        }
        let low = rho / (4.0 * PI * c) * omega_lo;
        let tail = high_frequency_resistance / omega_top;
        let static_mass = 2.0 / PI * (integral + low + tail);
        let static_mass_tail_share = tail / (integral + low + tail);

        // Samples to the integral's top or the mechanics Nyquist, whichever
        // is higher, and the fit.
        let omega_bottom = TAU * WINDOW_RADIATION_BAND_HZ[0];
        let nyquist = PI / time_step_s;
        let ratio = det::ln(omega_top.max(nyquist) / omega_bottom) / (FIT_SAMPLES - 1) as f64;
        let mut samples = Vec::with_capacity(FIT_SAMPLES);
        for i in 0..FIT_SAMPLES {
            let omega = omega_bottom * det::exp(ratio * i as f64);
            samples.push((omega, resistance(omega)?));
        }
        let band_top = TAU * WINDOW_RADIATION_BAND_HZ[1];
        let fitted = fit(&samples, band_top, static_mass, high_frequency_resistance, a_eq, c)?;
        let load = ParallelCellSpec::new(&fitted).map_err(physics)?;
        ParallelCellLoad::new(load, time_step_s).map_err(physics)?;

        let realized = |omega: f64| realized_resistance(&load, omega);
        let warped = |omega: f64| 2.0 * det::tan(0.5 * omega * time_step_s) / time_step_s;
        let (mut band_error, mut discrete_band_error, mut tail_error) = (0.0f64, 0.0f64, 0.0f64);
        let mut tail_top = 0.0f64;
        for &(omega, r) in &samples {
            if omega <= band_top {
                band_error = band_error.max((realized(omega) / r - 1.0).abs());
                discrete_band_error = discrete_band_error.max((realized(warped(omega)) / r - 1.0).abs());
            } else if omega <= nyquist {
                tail_error = tail_error.max((realized(omega) / r - 1.0).abs());
                tail_top = tail_top.max(omega);
            }
        }
        if !(band_error.is_finite() && band_error <= MAX_WINDOW_RADIATION_ERROR) {
            return Err(invalid(
                "window radiation realization misses the resistance by more than 5% from 20 Hz to 24 kHz",
            ));
        }
        Ok(Self {
            flux,
            outer_radius_m,
            cells: cells.len(),
            open_area_m2: opening.area,
            exterior,
            time_step_s,
            static_mass,
            static_mass_tail_share,
            high_frequency_resistance,
            ka_band: WINDOW_RADIATION_BAND_HZ.map(|f| TAU * f / c * outer_radius_m),
            band_error,
            discrete_band_error,
            tail_error,
            tail_top_hz: tail_top / TAU,
            load,
        })
    }

    /// The radiation resistance of the opening at `frequency_hz`, from the
    /// series itself rather than its realization [Pa s/m^3].
    ///
    /// # Errors
    /// A frequency that is not finite and positive, or a cell refusal.
    pub fn exact_resistance(&self, cells: &[OuterFaceCell], frequency_hz: f64) -> Result<f64, FlueError> {
        if !(frequency_hz.is_finite() && frequency_hz > 0.0) {
            return Err(invalid("window radiation query requires a positive finite frequency"));
        }
        Opening::new(self.outer_radius_m, cells)?.resistance(
            TAU * frequency_hz,
            self.exterior.density,
            self.exterior.sound_speed,
            &mut Quadrature::default(),
        )
    }

    /// The realized resistance at `frequency_hz`, continuous or as the
    /// bilinear midpoint step realizes it [Pa s/m^3].
    #[must_use]
    pub fn realized_resistance(&self, frequency_hz: f64, discrete: bool) -> f64 {
        let mut omega = TAU * frequency_hz;
        if discrete {
            omega = 2.0 * det::tan(0.5 * omega * self.time_step_s) / self.time_step_s;
        }
        realized_resistance(&self.load, omega)
    }

    /// The passive cells.
    #[must_use]
    pub const fn load(&self) -> &ParallelCellSpec { &self.load }
    /// Which flux distribution was radiated.
    #[must_use]
    pub const fn flux(&self) -> WindowFlux { self.flux }
    /// The outer cylinder's radius `A` [m].
    #[must_use]
    pub const fn outer_radius_m(&self) -> f64 { self.outer_radius_m }
    /// Cells of the opening.
    #[must_use]
    pub const fn cells(&self) -> usize { self.cells }
    /// Open area on the outer cylinder [m^2].
    #[must_use]
    pub const fn open_area_m2(&self) -> f64 { self.open_area_m2 }
    /// The exterior gas.
    #[must_use]
    pub const fn exterior(&self) -> &GasState { &self.exterior }
    /// Kramers-Kronig static mass `M0` [Pa s^2/m^3].
    #[must_use]
    pub const fn static_mass_pa_s2_m3(&self) -> f64 { self.static_mass }
    /// The share of `M0` above the integral's top, taken at `R_inf`.
    #[must_use]
    pub const fn static_mass_tail_share(&self) -> f64 { self.static_mass_tail_share }
    /// `R_inf` [Pa s/m^3].
    #[must_use]
    pub const fn high_frequency_resistance_pa_s_m3(&self) -> f64 { self.high_frequency_resistance }
    /// `kA` at 20 Hz and at 24 kHz.
    #[must_use]
    pub const fn ka_band(&self) -> [f64; 2] { self.ka_band }
    /// Largest relative resistance error of the continuous realization, 20 Hz to 24 kHz.
    #[must_use]
    pub const fn band_error(&self) -> f64 { self.band_error }
    /// The same for the bilinear midpoint step.
    #[must_use]
    pub const fn discrete_band_error(&self) -> f64 { self.discrete_band_error }
    /// Largest relative resistance error from 24 kHz to [`Self::tail_top_hz`].
    #[must_use]
    pub const fn tail_error(&self) -> f64 { self.tail_error }
    /// The last fit sample at or below the mechanics Nyquist [Hz].
    #[must_use]
    pub const fn tail_top_hz(&self) -> f64 { self.tail_top_hz }
}

/// `Re Z` of the cell chain at `omega`.
fn realized_resistance(load: &ParallelCellSpec, omega: f64) -> f64 {
    load.cells()
        .iter()
        .map(|cell| {
            // 1/Z = 1/R + i (omega C - 1/(omega L)).
            let susceptance = omega * cell.compliance_m3_pa - 1.0 / (omega * cell.inertance_pa_s2_m3);
            let r = cell.resistance_pa_s_m3;
            r / (1.0 + r * r * susceptance * susceptance)
        })
        .sum()
}

/// The opening grouped for the series: cells of one angle and span share a
/// column, cells of one axial centre and length a strip.
struct Opening {
    radius: f64,
    /// (angle, span) [rad].
    columns: Vec<(f64, f64)>,
    /// (centre, length) [m].
    strips: Vec<(f64, f64)>,
    /// (column, strip, `(q_p/U)/(dphi_p l_p)`).
    cells: Vec<(usize, usize, f64)>,
    area: f64,
    /// `sum_p (q_p/U)^2/S_p` [1/m^2].
    concentration: f64,
    /// Twice the largest axial reach from the window centre [m].
    axial_span: f64,
    /// Angle the cells span [rad].
    angular_span: f64,
}

impl Opening {
    fn new(radius: f64, cells: &[OuterFaceCell]) -> Result<Self, FlueError> {
        if !(radius.is_finite() && radius > 0.0) || cells.is_empty()
            || cells.len() > MAX_WINDOW_RADIATION_CELLS
        {
            return Err(invalid(
                "window radiation requires a positive outer radius and 1 to 65536 cells",
            ));
        }
        let total: f64 = cells.iter().map(|cell| cell.flux).sum();
        if cells.iter().any(|cell| {
            !(cell.arc > 0.0 && cell.arc <= TAU && cell.length > 0.0 && cell.length.is_finite()
                && cell.theta.is_finite() && cell.axial.is_finite() && cell.flux.is_finite())
        }) || !(total.is_finite() && total > 0.0)
        {
            return Err(invalid(
                "window radiation cells need finite positions, positive spans and a positive total flux",
            ));
        }
        let mut columns = Vec::new();
        let mut strips = Vec::new();
        let mut column_index = BTreeMap::new();
        let mut strip_index = BTreeMap::new();
        let mut grouped = Vec::with_capacity(cells.len());
        let (mut area, mut concentration, mut axial_span) = (0.0, 0.0, 0.0f64);
        let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
        for cell in cells {
            let column = *column_index
                .entry((cell.theta.to_bits(), cell.arc.to_bits()))
                .or_insert_with(|| {
                    columns.push((cell.theta, cell.arc));
                    columns.len() - 1
                });
            let strip = *strip_index
                .entry((cell.axial.to_bits(), cell.length.to_bits()))
                .or_insert_with(|| {
                    strips.push((cell.axial, cell.length));
                    strips.len() - 1
                });
            let share = cell.flux / total;
            let cell_area = radius * cell.arc * cell.length;
            grouped.push((column, strip, share / (cell.arc * cell.length)));
            area += cell_area;
            concentration += share * share / cell_area;
            axial_span = axial_span.max(2.0 * (cell.axial.abs() + 0.5 * cell.length));
            low = low.min(cell.theta - 0.5 * cell.arc);
            high = high.max(cell.theta + 0.5 * cell.arc);
        }
        Ok(Self {
            radius,
            columns,
            strips,
            cells: grouped,
            area,
            concentration,
            axial_span,
            angular_span: high - low,
        })
    }

    /// The series resistance at `omega` [Pa s/m^3].
    fn resistance(&self, omega: f64, rho: f64, c: f64, quadrature: &mut Quadrature) -> Result<f64, FlueError> {
        let (a, k) = (self.radius, omega / c);
        if !(k.is_finite() && k > 0.0) {
            return Err(invalid("window radiation requires a positive finite wavenumber"));
        }
        let nodes = (4.0 * k * self.axial_span.max(a)) as usize + 64;
        let (t, w) = quadrature.rule(nodes.div_ceil(64) * 64);
        let top = (k * a).ceil() as usize + 40;
        let trig = self.trig(top);
        let mut phasors = vec![(0.0, 0.0); self.strips.len()];
        let mut sums = vec![(0.0, 0.0); self.columns.len()];
        let mut bessel = Bessel::default();
        let mut total = 0.0;
        for (&ti, &wi) in t.iter().zip(w) {
            let (sin, cos) = det::sin_cos(FRAC_PI_4 * (ti + 1.0));
            let (kz, kr) = (k * sin, k * cos);
            for (phasor, &(z, l)) in phasors.iter_mut().zip(&self.strips) {
                let half = 0.5 * kz * l;
                let b = if half.abs() < 1e-8 { l } else { l * det::sin(half) / half };
                let (s, co) = det::sin_cos(kz * z);
                *phasor = (b * co, -b * s);
            }
            sums.fill((0.0, 0.0));
            for &(column, strip, weight) in &self.cells {
                let (re, im) = phasors[strip];
                sums[column].0 += weight * re;
                sums[column].1 += weight * im;
            }
            let derivatives = bessel.hankel_derivatives(kr * a, ((kr * a).ceil() as usize + 40).min(top));
            let mut sum = 0.0;
            for (m, &h2) in derivatives.iter().enumerate() {
                let row = &trig[m * self.columns.len()..(m + 1) * self.columns.len()];
                let (mut cr, mut ci, mut sr, mut si) = (0.0, 0.0, 0.0, 0.0);
                for (&(cm, sm), &(re, im)) in row.iter().zip(&sums) {
                    cr += cm * re;
                    ci += cm * im;
                    sr += sm * re;
                    si += sm * im;
                }
                let folded = if m == 0 { 1.0 } else { 2.0 };
                sum += folded * (cr * cr + ci * ci + sr * sr + si * si) / h2;
            }
            total += wi * sum / kr;
        }
        let r = omega * rho / (PI * PI * PI * a * a) * FRAC_PI_4 * total;
        if !(r.is_finite() && r > 0.0) {
            return Err(invalid("window radiation resistance is not finite and positive"));
        }
        Ok(r)
    }

    /// `(cos(m phi) s_m(dphi), sin(m phi) s_m(dphi))` per order `m <= top`
    /// and column, row-major by order.
    fn trig(&self, top: usize) -> Vec<(f64, f64)> {
        let mut out = Vec::with_capacity((top + 1) * self.columns.len());
        for m in 0..=top {
            for &(phi, span) in &self.columns {
                let s = if m == 0 { span } else { 2.0 * det::sin(0.5 * m as f64 * span) / m as f64 };
                let (sin, cos) = det::sin_cos(m as f64 * phi);
                out.push((cos * s, sin * s));
            }
        }
        out
    }
}

/// Gauss-Legendre rules on `[-1, 1]` by node count, built once each.
#[derive(Default)]
struct Quadrature {
    rules: BTreeMap<usize, (Vec<f64>, Vec<f64>)>,
}

impl Quadrature {
    fn rule(&mut self, nodes: usize) -> (&[f64], &[f64]) {
        let (t, w) = self
            .rules
            .entry(nodes)
            .or_insert_with(|| fs_feec::highorder::quad1d::gauss_legendre(nodes));
        (t, w)
    }
}

/// Workspace for `|H_m'(x)|^2`.
#[derive(Default)]
struct Bessel {
    j: Vec<f64>,
    y: Vec<f64>,
    out: Vec<f64>,
}

impl Bessel {
    /// `|H_m'(x)|^2 = J_m'(x)^2 + Y_m'(x)^2` for `m = 0..=top`, `x > 0`,
    /// stopping early where `Y_m` passes 1e150: those terms are below 1e-300
    /// of the rest. `Y` recurs upward from `Y_0`, `Y_1`; `J` by Miller's
    /// downward recurrence, normalized by `J_0 + 2 sum J_2k = 1`.
    fn hankel_derivatives(&mut self, x: f64, top: usize) -> &[f64] {
        let orders = top + 2;
        self.y.clear();
        self.y.push(fs_aeroac::bessel::y0(x));
        self.y.push(fs_aeroac::bessel::y1(x));
        while self.y.len() < orders {
            let m = self.y.len() - 1;
            let next = 2.0 * m as f64 / x * self.y[m] - self.y[m - 1];
            if !(next.abs() <= 1e150) {
                break;
            }
            self.y.push(next);
        }
        let orders = self.y.len();
        let start = orders + det::sqrt(200.0 * orders as f64) as usize + 20;
        let start = start + start % 2;
        self.j.clear();
        self.j.resize(orders, 0.0);
        let (mut above, mut here, mut norm) = (0.0f64, 1e-300f64, 0.0f64);
        for n in (0..=start).rev() {
            if n < orders {
                self.j[n] = here;
            }
            if n % 2 == 0 {
                norm += if n == 0 { here } else { 2.0 * here };
            }
            if n > 0 {
                let below = 2.0 * n as f64 / x * here - above;
                above = here;
                here = below;
                if here.abs() > 1e250 {
                    for v in &mut self.j {
                        *v *= 1e-250;
                    }
                    norm *= 1e-250;
                    above *= 1e-250;
                    here *= 1e-250;
                }
            }
        }
        for v in &mut self.j {
            *v /= norm;
        }
        self.out.clear();
        for m in 0..orders - 1 {
            let (dj, dy) = if m == 0 {
                (-self.j[1], -self.y[1])
            } else {
                (
                    0.5 * (self.j[m - 1] - self.j[m + 1]),
                    0.5 * (self.y[m - 1] - self.y[m + 1]),
                )
            };
            self.out.push(dj * dj + dy * dy);
        }
        &self.out
    }
}

/// One element of the dictionary: an R-L pole at its rate, or an R-L-C
/// resonance at its angular frequency.
#[derive(Clone, Copy)]
enum Element {
    Pole(f64),
    Resonance(f64),
}

impl Element {
    /// `Re Z` per unit cell resistance.
    fn shape(self, omega: f64) -> f64 {
        match self {
            Self::Pole(p) => omega * omega / (omega * omega + p * p),
            Self::Resonance(w0) => {
                let d = RESONANCE_Q * (omega / w0 - w0 / omega);
                1.0 / (1.0 + d * d)
            }
        }
    }
    /// Static mass per unit cell resistance.
    fn mass(self) -> f64 {
        match self {
            Self::Pole(p) => 1.0 / p,
            Self::Resonance(w0) => 1.0 / (RESONANCE_Q * w0),
        }
    }
    /// High-frequency resistance per unit cell resistance.
    fn high(self) -> f64 {
        match self {
            Self::Pole(_) => 1.0,
            Self::Resonance(_) => 0.0,
        }
    }
    fn cell(self, resistance: f64) -> ParallelCell {
        match self {
            Self::Pole(p) => ParallelCell {
                resistance_pa_s_m3: resistance,
                inertance_pa_s2_m3: resistance / p,
                compliance_m3_pa: 0.0,
            },
            Self::Resonance(w0) => ParallelCell {
                resistance_pa_s_m3: resistance,
                inertance_pa_s2_m3: resistance / (RESONANCE_Q * w0),
                compliance_m3_pa: RESONANCE_Q / (resistance * w0),
            },
        }
    }
}

/// The weighted fit on the `active` elements over `passes` Lawson passes:
/// the best score, the resistances and the relative errors per sample.
struct Fitted {
    score: f64,
    resistances: Vec<f64>,
}

#[allow(clippy::too_many_arguments)]
fn weighted_fit(
    dictionary: &[Element],
    active: &[usize],
    samples: &[(f64, f64)],
    base: &[f64],
    band_top: f64,
    static_mass: f64,
    high_resistance: f64,
    passes: usize,
) -> Fitted {
    let (n, cols) = (samples.len(), active.len());
    let rows = n + 2;
    let mut weights = base.to_vec();
    let mut best: Option<Fitted> = None;
    let mut matrix = vec![0.0; rows * cols];
    for _ in 0..passes {
        for (i, &(omega, r)) in samples.iter().enumerate() {
            for (j, &e) in active.iter().enumerate() {
                matrix[i * cols + j] = dictionary[e].shape(omega) * weights[i] / r;
            }
        }
        for (j, &e) in active.iter().enumerate() {
            matrix[n * cols + j] = EQUALITY_WEIGHT * dictionary[e].mass() / static_mass;
            matrix[(n + 1) * cols + j] = EQUALITY_WEIGHT * dictionary[e].high() / high_resistance;
        }
        // Unit columns, so the solver's fixed gradient threshold is relative.
        let norms: Vec<f64> = (0..cols)
            .map(|j| det::sqrt((0..rows).map(|i| matrix[i * cols + j] * matrix[i * cols + j]).sum()))
            .collect();
        let scaled: Vec<f64> = matrix
            .iter()
            .enumerate()
            .map(|(k, &v)| v / norms[k % cols])
            .collect();
        let mut rhs = weights.clone();
        rhs.extend([EQUALITY_WEIGHT, EQUALITY_WEIGHT]);
        let x = fs_material::visco::nnls(&scaled, &rhs, rows, cols);
        let resistances: Vec<f64> = x.iter().zip(&norms).map(|(x, norm)| x / norm).collect();
        let errors: Vec<f64> = samples
            .iter()
            .map(|&(omega, r)| {
                let fit: f64 = active
                    .iter()
                    .zip(&resistances)
                    .map(|(&e, g)| dictionary[e].shape(omega) * g)
                    .sum();
                (fit / r - 1.0).abs()
            })
            .collect();
        let (mut in_band, mut above) = (0.0f64, 0.0f64);
        for (&(omega, _), &e) in samples.iter().zip(&errors) {
            if omega <= band_top {
                in_band = in_band.max(e);
            } else {
                above = above.max(e);
            }
        }
        let score = in_band.max(TAIL_WEIGHT * above);
        if best.as_ref().is_none_or(|b| score < b.score) {
            best = Some(Fitted { score, resistances });
        }
        let largest = errors.iter().copied().fold(0.0, f64::max);
        if !(largest > 0.0) {
            break;
        }
        for ((w, &e), &b) in weights.iter_mut().zip(&errors).zip(base) {
            *w *= (e / largest + 0.02) * b;
        }
        let mean = weights.iter().sum::<f64>() / n as f64;
        for w in &mut weights {
            *w /= mean;
        }
    }
    best.expect("at least one Lawson pass")
}

/// The dictionary fit and its pruning to at most [`MAX_PARALLEL_CELLS`].
fn fit(
    samples: &[(f64, f64)],
    band_top: f64,
    static_mass: f64,
    high_resistance: f64,
    a_eq: f64,
    c: f64,
) -> Result<Vec<ParallelCell>, FlueError> {
    let mut dictionary = Vec::new();
    let mut ka = 0.5;
    while ka <= 256.0 * 1.0001 {
        dictionary.push(Element::Pole(ka * c / a_eq));
        ka *= 2.0;
    }
    let mut ka = 0.05;
    while ka <= 10.0 * 1.0001 {
        dictionary.push(Element::Resonance(ka * c / a_eq));
        ka *= 1.45;
    }
    let base: Vec<f64> = samples
        .iter()
        .map(|&(omega, _)| if omega <= band_top { 1.0 } else { TAIL_WEIGHT })
        .collect();
    let solve = |active: &[usize], passes: usize| {
        weighted_fit(&dictionary, active, samples, &base, band_top, static_mass, high_resistance, passes)
    };
    let all: Vec<usize> = (0..dictionary.len()).collect();
    let first = solve(&all, LAWSON_PASSES);
    let mut active: Vec<usize> = all.into_iter().filter(|&e| first.resistances[e] > 0.0).collect();
    while active.len() > MAX_PARALLEL_CELLS {
        let mut trials: Vec<(f64, Vec<usize>)> = active
            .iter()
            .filter(|&&e| matches!(dictionary[e], Element::Resonance(_)))
            .map(|&drop| {
                let trial: Vec<usize> = active.iter().copied().filter(|&e| e != drop).collect();
                (solve(&trial, PRUNED_PASSES).score, trial)
            })
            .collect();
        if trials.is_empty() {
            return Err(invalid("window radiation needs more than 16 poles"));
        }
        trials.sort_by(|p, q| p.0.total_cmp(&q.0));
        let trial = trials.swap_remove(0).1;
        let refit = solve(&trial, PRUNED_PASSES);
        active = trial
            .iter()
            .zip(&refit.resistances)
            .filter(|&(_, &g)| g > 0.0)
            .map(|(&e, _)| e)
            .collect();
    }
    let last = solve(&active, PRUNED_PASSES);
    Ok(active
        .iter()
        .zip(&last.resistances)
        .filter(|&(_, &g)| g > 0.0)
        .map(|(&e, &g)| dictionary[e].cell(g))
        .collect())
}

/// A port adapter only: memory, stepping and energy stay in fs-vfit.
#[derive(Clone, Copy, Debug)]
pub(super) struct WindowRadiation {
    pub model: WindowRadiationLoad,
    state: ParallelCellLoad,
}

impl WindowRadiation {
    fn new(model: WindowRadiationLoad) -> Result<Self, FlueError> {
        let state = ParallelCellLoad::new(model.load, model.time_step_s).map_err(physics)?;
        Ok(Self { model, state })
    }
    /// Step-mean pressure = direct * CURRENT mouth flow + retained history.
    pub fn affine_pressure(&self) -> (f64, f64) {
        self.state.affine_pressure()
    }
    pub fn candidate(&self, flow: f64) -> Result<(Self, ParallelCellFrame), FlueError> {
        let mut next = *self;
        let frame = next.state.step(flow).map_err(physics)?;
        Ok((next, frame))
    }
}

impl FlueVoice {
    /// Replace the authored window loss with the window's radiation load
    /// before any accepted sample. It needs the computed window mass
    /// ([`FlueVoice::with_window_mass`]) first: the series inertance becomes
    /// that mass less the cells' static mass, so the static mass from the
    /// junction to the far field stays the computed one.
    ///
    /// # Errors
    /// Late or repeated attachment, no window mass, a compact or joint mouth
    /// reaction, a nonzero `mouth_loss_rel`, a different clock or exterior
    /// pressure, or a radiation mass not below the window's.
    pub fn with_window_radiation(mut self, model: WindowRadiationLoad) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || !self.window_mass || self.window_radiation.is_some()
            || self.mouth_radiation.is_some()
            || self.dg().is_some_and(|duct| duct.load().radiation_record()
                .is_some_and(|record| record.driven_mouth().is_some()))
        {
            return Err(invalid(
                "window radiation requires one cold attachment after the window mass and no other mouth reaction",
            ));
        }
        if self.geometry.mouth_loss_rel != 0.0 {
            return Err(invalid(
                "the window's radiation replaces the authored loss; declare mouth_loss_rel zero",
            ));
        }
        if model.time_step_s.to_bits() != self.dt.to_bits()
            || model.exterior.pressure.to_bits() != self.gas.pressure.to_bits()
        {
            return Err(invalid(
                "window radiation must share the voice's clock and static pressure",
            ));
        }
        let inertance = self.inertance - model.load.static_inertance_pa_s2_m3();
        if !(inertance.is_finite() && inertance > 0.0) {
            return Err(invalid(
                "the window's radiation mass must lie below its computed mass",
            ));
        }
        self.window_radiation = Some(WindowRadiation::new(model)?);
        self.inertance = inertance;
        Ok(self)
    }

    /// The window's radiation load, when attached.
    #[must_use]
    pub fn window_radiation_load(&self) -> Option<&WindowRadiationLoad> {
        self.window_radiation.as_ref().map(|load| &load.model)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::GasSpec;

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
    }

    fn piston(arc: f64, length: f64) -> [OuterFaceCell; 1] {
        [OuterFaceCell { theta: 0.0, arc, axial: 0.0, length, flux: 1.0 }]
    }

    #[test]
    fn a_piston_matches_the_independent_series_and_the_monopole() {
        // cylfast.py (scipy jv/yv, numpy Gauss-Legendre): A 13.5 mm, piston
        // 12 x 4.6 mm, rho 1.2, c 343.
        let (a, w, l) = (0.0135, 0.012, 0.0046);
        let opening = Opening::new(a, &piston(w / a, l)).expect("opening");
        let mut quadrature = Quadrature::default();
        for (f, expected) in [
            (200.0, 438.302_181_355_394_5),
            (2000.0, 48_993.220_151_327_19),
            (8000.0, 1_022_477.099_945_981_4),
            (24000.0, 5_317_664.952_413_325),
            (60000.0, 7_591_379.205_054_756),
        ] {
            let r = opening.resistance(TAU * f, 1.2, 343.0, &mut quadrature).expect("series");
            assert!((r / expected - 1.0).abs() < 1e-6, "{f} Hz: {r} vs {expected}");
        }
        let low = opening.resistance(TAU * 2.0, 1.2, 343.0, &mut quadrature).expect("series");
        let monopole = 1.2 * (TAU * 2.0) * (TAU * 2.0) / (4.0 * PI * 343.0);
        assert!((low / monopole - 1.0).abs() < 1e-4, "{low} vs monopole {monopole}");
    }

    #[test]
    fn splitting_a_piston_into_cells_leaves_its_resistance() {
        let (a, w, l) = (0.0085, 0.008, 0.0053);
        let mut cells = Vec::new();
        for i in 0..3 {
            for j in 0..2 {
                cells.push(OuterFaceCell {
                    theta: (i as f64 - 1.0) * w / (3.0 * a),
                    arc: w / (3.0 * a),
                    axial: (j as f64 - 0.5) * 0.5 * l,
                    length: 0.5 * l,
                    flux: 1.0 / 6.0,
                });
            }
        }
        let whole = Opening::new(a, &piston(w / a, l)).expect("whole");
        let split = Opening::new(a, &cells).expect("split");
        let mut quadrature = Quadrature::default();
        for f in [300.0, 5000.0, 30000.0] {
            let (p, q) = (
                whole.resistance(TAU * f, 1.2, 343.0, &mut quadrature).expect("whole"),
                split.resistance(TAU * f, 1.2, 343.0, &mut quadrature).expect("split"),
            );
            assert!((p / q - 1.0).abs() < 1e-11, "{f} Hz: {p} vs {q}");
        }
        assert!((split.concentration * w * l - 1.0).abs() < 1e-12);
    }

    #[test]
    fn the_realization_holds_the_band_and_the_static_mass() {
        // cylfast.static_mass for the same piston at rho 1.2: 62.104999595
        // kg/m^4. M0 is rho times an integral over k, so it scales with the
        // density alone.
        let (a, w, l) = (0.0135, 0.012, 0.0046);
        let gas = air();
        let load = WindowRadiationLoad::new(a, &piston(w / a, l), WindowFlux::CutUpPiston, gas,
            1.0 / 288_000.0, &CancelGate::new_clock_free())
            .expect("fit");
        let m0 = load.static_mass_pa_s2_m3();
        let expected = 62.104_999_595_372_97 * gas.density / 1.2;
        assert!((m0 / expected - 1.0).abs() < 1e-5, "M0 {m0} vs {expected}");
        let realized = load.load().static_inertance_pa_s2_m3();
        assert!((realized / m0 - 1.0).abs() < 2e-3, "realized mass {realized} vs M0 {m0}");
        assert!(load.load().cells().len() <= MAX_PARALLEL_CELLS);
        assert!(load.band_error() <= 0.05, "band error {}", load.band_error());
        eprintln!(
            "cells {}, band {:.4}, discrete {:.4}, tail {:.4} to {:.0} Hz, tail share {:.4}, kA {:?}",
            load.load().cells().len(), load.band_error(), load.discrete_band_error(),
            load.tail_error(), load.tail_top_hz(), load.static_mass_tail_share(), load.ka_band()
        );
    }
}
