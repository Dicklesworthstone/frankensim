//! Finite-interval Joule heat from one original three-dimensional AC response.
//!
//! For the native peak convention, E(t) = Re(E exp(-i(theta + omega*t)))
//! and likewise for J. The interval starts at the authored phase theta;
//! no absolute timestamp is divided by frequency or rounded into a midpoint.
//! The shared analytic temporal Gram kernel retains the twice-frequency term
//! and its positive small-window limit. The consistent tetrahedral P1 mass
//! form integrates the original conductor-frame E dot J in space.
//!
//! This is an observation of a settled, frozen electromagnetic response. It
//! supplies resolved conductor heat, not source work, magnetic storage change,
//! or mechanical work. It does not evolve temperature, include an operating DC
//! field, add stranded-winding resistance, or resolve switch-on transients.
//! A thermal receiver must retain the original material correspondence, admit
//! its frozen electrical-law approximation and deposit this heat exactly once.

use super::{SpatialHarmonicInductionCell, SpatialHarmonicInductionResponse};
use crate::spatial::{SpatialMagneticError, finite, invalid, poll};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_math::{c64::C64, dd::Dd, det};

type Result<T> = core::result::Result<T, SpatialMagneticError>;

/// Additional observation bounds; the original field is borrowed, not copied.
#[derive(Clone, Copy, Debug)]
pub struct SpatialHarmonicHeatLimits {
    pub maximum_cells: usize,
    /// Logical result/header and cell-record bytes, excluding allocator metadata
    /// and the original response. Temporary work is constant-size stack data.
    pub maximum_output_bytes: usize,
}

/// Heat in one original conducting tetrahedron, in original response order.
#[derive(Clone, Copy, Debug)]
pub struct SpatialHarmonicHeatCell<'r, C = ResolvedConductorState> {
    pub source: &'r C,
    pub cell: usize,
    pub component: usize,
    pub volume_m3: f64,
    /// Physical interval mean, which need not equal the cycle mean [W].
    pub mean_power_w: f64,
    pub heat_j: f64,
}

/// Unweighted physical heat over one finite interval of the original AC field.
///
/// The immutable response borrow and per-cell source references preserve the
/// actual conducting regions and constitutive state used by the field solve.
/// The result is a Joule observation, not a complete electromagnetic budget.
pub struct SpatialHarmonicHeatExposure<
    'r,
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    response: &'r SpatialHarmonicInductionResponse<'s, 'm, C>,
    duration_s: f64,
    phase_at_start_rad: f64,
    cells: Vec<SpatialHarmonicHeatCell<'r, C>>,
    mean_power_w: f64,
    heat_j: f64,
}

impl<'r, 's, 'm, C: SpatialConductorMaterial> SpatialHarmonicHeatExposure<'r, 's, 'm, C> {
    pub fn response(&self) -> &'r SpatialHarmonicInductionResponse<'s, 'm, C> {
        self.response
    }
    pub fn cells(&self) -> &[SpatialHarmonicHeatCell<'r, C>] {
        &self.cells
    }
    pub fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Physical field phase advances as phase_at_start_rad + omega*t.
    pub fn phase_at_start_rad(&self) -> f64 {
        self.phase_at_start_rad
    }
    pub fn mean_power_w(&self) -> f64 {
        self.mean_power_w
    }
    pub fn heat_j(&self) -> f64 {
        self.heat_j
    }
}

impl<'s, 'm, C: SpatialConductorMaterial> SpatialHarmonicInductionResponse<'s, 'm, C> {
    /// Integrate original conductor-frame E(t) dot J(t) on [0, duration_s].
    ///
    /// Phase is authored at the interval start in radians. The physical field
    /// is Re(phasor * exp(-i*(phase_at_start_rad + omega*t))); positive phase
    /// therefore multiplies the phasor's imaginary component by positive sine.
    /// One shared temporal kernel and O(conducting cells) P1 contractions serve
    /// the complete observation without field solves or temporal sampling.
    /// Cancellation/refusal publishes no partial exposure or field mutation.
    pub fn observe_joule_heat<'r>(
        &'r self,
        duration_s: f64,
        phase_at_start_rad: f64,
        limits: SpatialHarmonicHeatLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicHeatExposure<'r, 's, 'm, C>> {
        poll(cancelled)?;
        if !phase_at_start_rad.is_finite() {
            return Err(invalid("harmonic heat needs a finite phase at interval start"));
        }
        let window = HarmonicWindow::new(0.0, duration_s)
            .map_err(|_| invalid("harmonic heat needs a finite positive duration"))?;
        let bytes = self
            .cells
            .len()
            .checked_mul(core::mem::size_of::<SpatialHarmonicHeatCell<'r, C>>())
            .and_then(|n| {
                n.checked_add(core::mem::size_of::<SpatialHarmonicHeatExposure<'r, 's, 'm, C>>())
            })
            .ok_or_else(|| invalid("harmonic heat output size overflow"))?;
        if self.cells.len() > limits.maximum_cells || bytes > limits.maximum_output_bytes {
            return Err(invalid("harmonic heat observation exceeds its output bounds"));
        }
        let kernel = window
            .product_kernel(self.omega_rad_s, self.omega_rad_s)
            .map_err(|_| invalid("harmonic heat temporal moment is not representable"))?;
        let gram = kernel.centered_gram_dd();
        // Rotate by the start phase separately from the half-width. DD retains
        // both cancellations and a tiny interval at a large phase. The same
        // normalized temporal owner supplies the positive centered weights.
        let (sine, cosine) = det::sin_cos_dd(Dd::from_f64(phase_at_start_rad));
        let duration = Dd::from_f64(duration_s);
        let omega = Dd::from_f64(self.omega_rad_s);
        let half_width = if self.omega_rad_s >= duration_s {
            checked((omega * Dd::from_f64(0.5)) * duration)?
        } else {
            checked(omega * Dd::from_f64(0.5 * duration_s))?
        };
        let (half_sine, half_cosine) = det::sin_cos_dd(half_width);
        let rotation = [
            checked(cosine * half_cosine - sine * half_sine)?,
            checked(-(cosine * half_sine + sine * half_cosine))?,
        ];
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(self.cells.len())
            .map_err(|_| invalid("harmonic heat output allocation failed"))?;
        let mut total = Dd::ZERO;
        for cell in &self.cells {
            poll(cancelled)?;
            let source = self
                .conductors()
                .get(cell.cell)
                .and_then(Option::as_ref)
                .ok_or_else(|| invalid("harmonic heat lost its original conducting source"))?;
            let power = cell_power(cell, source.conductivity().matrix_s_m(), gram, rotation)?;
            let mean_power_w = nonnegative(power)?;
            let heat_j = nonnegative(checked(power * duration)?)?;
            total = checked(total + power)?;
            cells.push(SpatialHarmonicHeatCell {
                source,
                cell: cell.cell,
                component: cell.component,
                volume_m3: cell.volume_m3,
                mean_power_w,
                heat_j,
            });
        }
        let mean_power_w = nonnegative(total)?;
        let heat_j = nonnegative(checked(total * duration)?)?;
        poll(cancelled)?;
        Ok(SpatialHarmonicHeatExposure {
            response: self,
            duration_s,
            phase_at_start_rad,
            cells,
            mean_power_w,
            heat_j,
        })
    }
}

fn checked(value: Dd) -> Result<Dd> {
    finite(value.hi)?;
    finite(value.lo)?;
    Ok(value)
}

fn nonnegative(value: Dd) -> Result<f64> {
    if value.lt(Dd::ZERO) {
        return Err(invalid("negative finite-interval induction Joule heat"));
    }
    finite(value.to_f64())
}

fn rotate(value: C64, [cosine, minus_sine]: [Dd; 2]) -> Result<[Dd; 2]> {
    let real = Dd::from_f64(value.re);
    let imaginary = Dd::from_f64(value.im);
    Ok([
        checked(real * cosine - imaginary * minus_sine)?,
        checked(real * minus_sine + imaginary * cosine)?,
    ])
}

fn mean_quadratic(e: [[Dd; 2]; 3], sigma: [[f64; 3]; 3], gram: [Dd; 2]) -> Result<Dd> {
    let mut power = Dd::ZERO;
    for part in 0..2 {
        let mut quadratic = Dd::ZERO;
        for i in 0..3 {
            let mut current = Dd::ZERO;
            for j in 0..3 {
                current = checked(current + Dd::from_f64(sigma[i][j]) * e[j][part])?;
            }
            quadratic = checked(quadratic + e[i][part] * current)?;
        }
        power = checked(power + gram[part] * quadratic)?;
    }
    Ok(power)
}

fn cell_power(
    cell: &SpatialHarmonicInductionCell,
    sigma: [[f64; 3]; 3],
    gram: [Dd; 2],
    rotation: [Dd; 2],
) -> Result<Dd> {
    let mut sum_e = [[Dd::ZERO; 2]; 3];
    let mut power = Dd::ZERO;
    // Integral N_i N_j = V*(1+delta_ij)/20. The factored form uses
    // sum_i E_i.J_i + (sum_i E_i).(sum_i J_i), only five vector products
    // instead of sixteen, and keeps the positive spatial mass form intact.
    // Apply the ORIGINAL sigma to the rotated E in DD. Independently rounded
    // E/J phasors can have different rounded zero crossings and yield negative
    // E.J in tiny windows at a current zero; evaluating their common source
    // law after rotation preserves the physical quadratic without clipping.
    for node in 0..4 {
        let mut electric = [[Dd::ZERO; 2]; 3];
        for axis in 0..3 {
            let e = rotate(cell.conductor_electric_field_v_m[node][axis], rotation)?;
            electric[axis] = e;
            for part in 0..2 {
                sum_e[axis][part] = checked(sum_e[axis][part] + e[part])?;
            }
        }
        power = checked(power + mean_quadratic(electric, sigma, gram)?)?;
    }
    power = checked(power + mean_quadratic(sum_e, sigma, gram)?)?;
    checked(power * Dd::from_f64(cell.volume_m3 / 20.0))
}
