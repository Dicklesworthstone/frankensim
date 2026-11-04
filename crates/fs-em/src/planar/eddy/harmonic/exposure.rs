//! Finite-time winding heat and magnetic work from the original AC fields.
//!
//! For real perturbations `a(t)=Re(sum a_j exp(-i omega_j t))`, the SAME
//! differential field operator gives `U= a^T K_t a/2 + sum L_end i^2/2`.
//! Its instantaneous power identity is `sum v i = integral sigma E^2 + dU/dt`.
//! Each frequency retains its original skin/proximity current distribution;
//! adding independent cycle-mean losses would discard both beat and carrier
//! terms. Exact P1 spatial products and the shared analytic temporal Gram
//! kernel integrate these terms without a time-stepping or fitted-loss model.
//!
//! This observes the existing fixed-bias, linearized AC system. It does not
//! add the operating DC heat, DC/AC cross work, switch-on transients, finite
//! saturation excursions, moving geometry or hysteresis. A thermal consumer
//! must deposit the returned spatial heat once, not add the cycle-mean heat
//! from the original frequency responses as a second source.

use super::*;
use fs_fft::harmonic_window::{HarmonicProduct, HarmonicWindow, HarmonicWindowProfile};
use fs_math::dd::Dd;

/// Bounded source grouping; pair and workspace limits also apply before allocation.
pub const MAX_EXPOSURE_LINES: usize = 256;

/// One authored peak-current tone on an original frequency-specific field owner.
#[derive(Clone, Copy)]
pub struct HarmonicWindingTone<'f, 'p, 'b, 'm> {
    pub frequency: &'f HarmonicWindingFrequency<'p, 'b, 'm>,
    /// One peak current per original winding, with phase at source time zero.
    pub currents_a: &'f [C64],
}

/// Additional spectral work/storage bounds and the physical power allowance.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicWindingExposureLimits {
    pub maximum_lines: usize,
    /// Conservative scalar temporal contractions per observation. Preparation
    /// uses the authored line count even when equal frequencies later combine.
    pub maximum_pair_products: usize,
    /// Additional bank, grouping, temporal-kernel and output payload. Borrowed
    /// source/frequency owners and each original response's separately bounded
    /// temporary payload are excluded; no native field solve is added here.
    pub maximum_workspace_bytes: usize,
    /// Absolute defect between terminal power, spatial heat and storage rate [W].
    pub power_tolerance_w: f64,
}

/// The source's electrical role; driven solid current includes its local
/// uniform-drive/circulating-current cross term before spatial integration.
#[derive(Clone, Copy, Debug)]
pub enum HarmonicWindingHeatSource<'s> {
    Floating(&'s EddyConductor),
    Solid {
        source: &'s EddyConductor,
        binding: SolidWindingRegion,
    },
    Stranded {
        source: &'s WindingConductor,
        winding: usize,
    },
}

/// Original spatial material source and its physical AC Joule deposit.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicWindingCellExposure<'s> {
    pub source: HarmonicWindingHeatSource<'s>,
    pub volume_m3: f64,
    pub mean_power_w: f64,
    pub heat_j: f64,
}

/// An unweighted physical interval of the settled, coherent AC perturbation.
#[derive(Clone, Debug)]
pub struct HarmonicWindingExposureResult<'s> {
    pub window: HarmonicWindow,
    /// Floating/solid cells in original conductor order, then stranded cells.
    pub cells: Vec<HarmonicWindingCellExposure<'s>>,
    pub mean_magnetic_energy_j: f64,
    pub mean_magnetic_storage_rate_w: f64,
    /// Integral of dU/dt, evaluated without subtracting rounded endpoint energies.
    pub magnetic_energy_change_j: f64,
    pub mean_dissipated_power_w: f64,
    pub dissipated_heat_j: f64,
    /// Real instantaneous terminal v(t)i(t) averaged over this same window;
    /// it can be negative while stored magnetic energy returns to the sources.
    pub mean_supplied_power_w: f64,
    pub supplied_work_j: f64,
    /// Heat plus storage rate minus terminal power.
    pub power_residual_w: f64,
    pub energy_residual_j: f64,
}

struct SpectralLine {
    omega_rad_s: f64,
    currents_a: Vec<C64>,
    voltages_v: Vec<C64>,
    potential_wb_m: Vec<C64>,
    magnetic_rhs_a_m: Vec<C64>,
    /// Total constitutive E, in original projected-conductor order.
    electric_v_m: Vec<[C64; 3]>,
}

/// Compact original field responses reused across arbitrarily many windows.
/// Queries require O((F^2+1)*(vertices+eddy_cells+windings)+stranded_cells) work
/// for F distinct frequencies, without field solves or per-cell trigonometry.
pub struct HarmonicWindingExposure<'p, 'b, 'm> {
    model: &'p PlanarHarmonicWindings<'b, 'm>,
    lines: Vec<SpectralLine>,
    limits: HarmonicWindingExposureLimits,
}

fn bad(message: &'static str) -> HarmonicWindingError {
    HarmonicWindingError::Invalid(message)
}
fn moment(result: core::result::Result<Dd, fs_fft::harmonic_window::WindowError>) -> Result<Dd> {
    result.map_err(|_| bad("winding temporal moment is not representable"))
}
fn real(value: Dd) -> Result<f64> {
    finite(value.hi)?;
    finite(value.lo)?;
    finite(value.to_f64())
}
fn nonnegative(value: Dd) -> Result<f64> {
    let value = real(value)?;
    if value < 0.0 {
        Err(HarmonicWindingError::Balance(
            "negative winding exposure heat or perturbation energy",
        ))
    } else {
        Ok(value)
    }
}
fn derivative(value: C64, omega: f64) -> Result<C64> {
    phasor(C64::new(omega * value.im, -omega * value.re))
}

impl<'b, 'm> PlanarHarmonicWindings<'b, 'm> {
    /// Bind original frequency responses and authored current phases to one
    /// unchanged material/geometry/bias owner. Equal-frequency currents combine
    /// with compensated arithmetic BEFORE the native field response is formed.
    /// An empty list observes zero AC perturbation, not the DC operating field.
    pub fn prepare_exposure<'p>(
        &'p self,
        tones: &[HarmonicWindingTone<'_, '_, 'b, 'm>],
        limits: HarmonicWindingExposureLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicWindingExposure<'p, 'b, 'm>> {
        poll(cancelled)?;
        let l = tones.len();
        let n = self.tangent.nrows();
        let m = self.basis.winding_count();
        let e = self.eddy_conductors.len();
        let c = self.winding_conductors.len();
        if l > MAX_EXPOSURE_LINES || l > limits.maximum_lines || !positive(limits.power_tolerance_w)
        {
            return Err(bad("invalid winding exposure limits or too many lines"));
        }
        // Four P1 square contractions per conductor (three vertices and their
        // sum), plus source power and both magnetic storage contractions.
        let products = l.checked_mul(l).and_then(|pairs| {
            n.checked_mul(2)
                .and_then(|v| m.checked_mul(3).and_then(|w| v.checked_add(w)))
                .and_then(|v| e.checked_mul(4).and_then(|w| v.checked_add(w)))
                .and_then(|v| pairs.checked_mul(v))
        });
        let bytes = n
            .checked_mul(4)
            .and_then(|v| m.checked_mul(4).and_then(|w| v.checked_add(w)))
            .and_then(|v| e.checked_mul(6).and_then(|w| v.checked_add(w)))
            .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
            .and_then(|v| v.checked_add(core::mem::size_of::<SpectralLine>()))
            .and_then(|v| v.checked_mul(l))
            .and_then(|v| {
                l.checked_mul(l)
                    .and_then(|w| w.checked_mul(core::mem::size_of::<HarmonicProduct>()))
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| {
                e.checked_add(c)
                    .and_then(|w| {
                        w.checked_mul(core::mem::size_of::<HarmonicWindingCellExposure<'_>>())
                    })
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| {
                l.checked_mul(core::mem::size_of::<usize>())
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| {
                m.checked_mul(core::mem::size_of::<[Dd; 3]>() + core::mem::size_of::<C64>())
                    .and_then(|w| v.checked_add(w))
            });
        if products.is_none_or(|v| v > limits.maximum_pair_products)
            || bytes.is_none_or(|v| v > limits.maximum_workspace_bytes)
        {
            return Err(bad("winding exposure exceeds its pair or workspace budget"));
        }
        for tone in tones {
            poll(cancelled)?;
            if !core::ptr::eq(tone.frequency.prepared, self)
                || tone.currents_a.len() != m
                || tone
                    .currents_a
                    .iter()
                    .any(|i| !i.re.is_finite() || !i.im.is_finite())
            {
                return Err(bad(
                    "exposure needs original frequency owners and one finite current per winding",
                ));
            }
        }
        let mut order: Vec<_> = (0..l).collect();
        order.sort_unstable_by(|&a, &b| {
            tones[a]
                .frequency
                .omega_rad_s
                .total_cmp(&tones[b].frequency.omega_rad_s)
                .then_with(|| a.cmp(&b))
        });
        let mut lines = Vec::with_capacity(l);
        let mut offset = 0;
        while offset < l {
            poll(cancelled)?;
            let first = tones[order[offset]];
            let omega_rad_s = first.frequency.omega_rad_s;
            let mut currents = vec![[Dd::ZERO; 2]; m];
            while offset < l && tones[order[offset]].frequency.omega_rad_s == omega_rad_s {
                for (sum, &current) in currents.iter_mut().zip(tones[order[offset]].currents_a) {
                    sum[0] = sum[0] + Dd::from_f64(current.re);
                    sum[1] = sum[1] + Dd::from_f64(current.im);
                }
                offset += 1;
            }
            let currents: Vec<_> = currents
                .into_iter()
                .map(|[re, im]| Ok(C64::new(real(re)?, real(im)?)))
                .collect::<Result<_>>()?;
            let response = first.frequency.response(&currents, cancelled)?;
            let mut magnetic_rhs_a_m = Vec::with_capacity(n);
            for row in 0..n {
                poll(cancelled)?;
                let (columns, values) = self.tangent.row(row);
                let mut sum = [Dd::ZERO; 2];
                for (&column, &value) in columns.iter().zip(values) {
                    let a = response.potential_wb_m[column];
                    sum[0] = sum[0] + Dd::from_f64(value) * Dd::from_f64(a.re);
                    sum[1] = sum[1] + Dd::from_f64(value) * Dd::from_f64(a.im);
                }
                magnetic_rhs_a_m.push(C64::new(real(sum[0])?, real(sum[1])?));
            }
            let mut electric_v_m = vec![[C64::ZERO; 3]; e];
            for (source, field) in response
                .eddy_cells
                .iter()
                .map(|cell| (&cell.source, cell.electric_field_v_m))
                .chain(
                    response
                        .solid_cells
                        .iter()
                        .map(|cell| (&cell.source, cell.electric_field_v_m)),
                )
            {
                poll(cancelled)?;
                let index = self
                    .eddy_conductors
                    .binary_search_by_key(&source.triangle, |s| s.triangle)
                    .map_err(|_| bad("original exposure conductor was not recovered"))?;
                electric_v_m[index] = field;
            }
            lines.push(SpectralLine {
                omega_rad_s,
                currents_a: response.currents_a,
                voltages_v: response.terminal_voltages_v,
                potential_wb_m: response.potential_wb_m,
                magnetic_rhs_a_m,
                electric_v_m,
            });
        }
        poll(cancelled)?;
        Ok(HarmonicWindingExposure {
            model: self,
            lines,
            limits,
        })
    }
}

impl<'p, 'b, 'm> HarmonicWindingExposure<'p, 'b, 'm> {
    /// Original material, geometry and DC-bias owner used by every spectral line.
    pub fn model(&self) -> &'p PlanarHarmonicWindings<'b, 'm> {
        self.model
    }

    pub fn distinct_frequencies(&self) -> usize {
        self.lines.len()
    }

    /// Analytic real-field observation on the physical source clock. Uniform
    /// windows are required so duration times storage rate is the actual energy
    /// change. Pair kernels are cached once; all cell sums retain DD arithmetic.
    pub fn observe(
        &self,
        window: HarmonicWindow,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicWindingExposureResult<'_>> {
        poll(cancelled)?;
        if window.profile() != HarmonicWindowProfile::Uniform {
            return Err(bad(
                "winding work and heat require an unweighted physical interval",
            ));
        }
        let f = self.lines.len();
        let mut kernels = Vec::with_capacity(f * f);
        for a in &self.lines {
            for b in &self.lines {
                poll(cancelled)?;
                kernels.push(
                    window
                        .product_kernel(a.omega_rad_s, b.omega_rad_s)
                        .map_err(|_| bad("winding temporal kernel is not representable"))?,
                );
            }
        }
        let mut energy = Dd::ZERO;
        let mut rate = Dd::ZERO;
        let mut supplied = Dd::ZERO;
        let mut current_squares = vec![Dd::ZERO; self.model.basis.winding_count()];
        for (i, a) in self.lines.iter().enumerate() {
            for (j, b) in self.lines.iter().enumerate() {
                let kernel = kernels[i * f + j];
                for (k, &end) in self.model.basis.end_inductance_h().iter().enumerate() {
                    poll(cancelled)?;
                    supplied = supplied
                        + moment(kernel.mean_real_product_dd(a.voltages_v[k], b.currents_a[k]))?;
                    let square =
                        moment(kernel.mean_real_product_dd(a.currents_a[k], b.currents_a[k]))?;
                    current_squares[k] = current_squares[k] + square;
                    energy = energy + Dd::from_f64(0.5 * end) * square;
                    rate = rate
                        + Dd::from_f64(end)
                            * moment(kernel.mean_real_product_dd(
                                a.currents_a[k],
                                derivative(b.currents_a[k], b.omega_rad_s)?,
                            ))?;
                }
                for (&ka, &b_potential) in a.magnetic_rhs_a_m.iter().zip(&b.potential_wb_m) {
                    poll(cancelled)?;
                    energy = energy
                        + Dd::from_f64(0.5) * moment(kernel.mean_real_product_dd(ka, b_potential))?;
                    rate = rate
                        + moment(
                            kernel
                                .mean_real_product_dd(ka, derivative(b_potential, b.omega_rad_s)?),
                        )?;
                }
            }
        }
        let mut cells = Vec::with_capacity(
            self.model.eddy_conductors.len() + self.model.winding_conductors.len(),
        );
        let mut heat = Dd::ZERO;
        for (cell_index, (source, cell)) in self
            .model
            .eddy_conductors
            .iter()
            .zip(&self.model.conduction.cells)
            .enumerate()
        {
            poll(cancelled)?;
            let mut integral = Dd::ZERO;
            for (i, a) in self.lines.iter().enumerate() {
                for (j, b) in self.lines.iter().enumerate() {
                    poll(cancelled)?;
                    let ea = a.electric_v_m[cell_index];
                    let eb = b.electric_v_m[cell_index];
                    let kernel = kernels[i * f + j];
                    // Exact P1 mass: volume/12 * (sum(Ea)*sum(Eb) + Ea.Eb).
                    let sum = |e: [C64; 3]| {
                        e.into_iter().fold([Dd::ZERO; 2], |s, v| {
                            [s[0] + Dd::from_f64(v.re), s[1] + Dd::from_f64(v.im)]
                        })
                    };
                    integral =
                        integral + moment(kernel.mean_real_product_parts_dd(sum(ea), sum(eb)))?;
                    for k in 0..3 {
                        integral = integral + moment(kernel.mean_real_product_dd(ea[k], eb[k]))?;
                    }
                }
            }
            let watts = integral * Dd::from_f64(cell.sigma_volume / 12.0);
            let power = nonnegative(watts)?;
            heat = heat + watts;
            let role = self
                .model
                .solid
                .binding_for_group(cell.group)
                .map_or(HarmonicWindingHeatSource::Floating(source), |binding| {
                    HarmonicWindingHeatSource::Solid { source, binding }
                });
            cells.push(HarmonicWindingCellExposure {
                source: role,
                volume_m3: cell.volume_m3,
                mean_power_w: power,
                heat_j: finite(power * window.duration_s())?,
            });
        }
        let mesh = self.model.basis.operating_field().problem().mesh();
        for (source, &winding) in self
            .model
            .winding_conductors
            .iter()
            .zip(&self.model.conductor_windings)
        {
            poll(cancelled)?;
            let profile = self.model.basis.profiles_per_m2()[winding][source.triangle];
            let volume = mesh
                .face_measure(source.triangle)
                .ok_or_else(|| bad("original stranded volume missing"))?;
            let watts = current_squares[winding]
                * Dd::from_f64(profile)
                * Dd::from_f64(profile)
                * Dd::from_f64(volume)
                * Dd::from_f64(source.conductor.resistivity_ohm_m());
            let power = nonnegative(watts)?;
            heat = heat + watts;
            cells.push(HarmonicWindingCellExposure {
                source: HarmonicWindingHeatSource::Stranded { source, winding },
                volume_m3: volume,
                mean_power_w: power,
                heat_j: finite(power * window.duration_s())?,
            });
        }
        let mean_magnetic_energy_j = nonnegative(energy)?;
        let mean_magnetic_storage_rate_w = real(rate)?;
        let mean_dissipated_power_w = nonnegative(heat)?;
        let mean_supplied_power_w = real(supplied)?;
        let power_residual_w = real(heat + rate - supplied)?;
        if power_residual_w.abs() > self.limits.power_tolerance_w {
            return Err(HarmonicWindingError::Balance(
                "finite-window winding power, heat and magnetic storage disagree",
            ));
        }
        poll(cancelled)?;
        Ok(HarmonicWindingExposureResult {
            window,
            cells,
            mean_magnetic_energy_j,
            mean_magnetic_storage_rate_w,
            magnetic_energy_change_j: finite(mean_magnetic_storage_rate_w * window.duration_s())?,
            mean_dissipated_power_w,
            dissipated_heat_j: finite(mean_dissipated_power_w * window.duration_s())?,
            mean_supplied_power_w,
            supplied_work_j: finite(mean_supplied_power_w * window.duration_s())?,
            power_residual_w,
            energy_residual_j: finite(power_residual_w * window.duration_s())?,
        })
    }
}
