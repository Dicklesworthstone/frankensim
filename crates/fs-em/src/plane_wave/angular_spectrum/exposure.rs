//! Finite-time interference of original frequency-resolved Maxwell fields.
//!
//! For authored spectral lines, the physical fields are
//! E(t)=Re(sum_j E_j exp(+i omega_j t)), and likewise H. The shared harmonic
//! window integrates their REAL products, retaining both sum and difference
//! frequencies. An ultrashort shutter therefore does not acquire an implicit
//! optical-cycle average. Source amplitudes carry their phase at t=0; they are
//! peak line amplitudes, not a power spectral density or frequency-bin weight.
//!
//! The rectangular time-average definition and real/complex field distinction
//! follow TU Delft, Interference and Coherence, equations 295 and 303–306:
//! <https://interactivetextbooks.tudelft.nl/tn2421/content/Chap5_Interference/InterferenceCoherence_2022_01Clean.html>.
//! No stationary-source, random-phase or slowly varying envelope assumption is
//! inferred. Each frequency uses its original material/film/diffraction solve.

use super::*;
use fs_fft::harmonic_window::{HarmonicProduct, HarmonicWindow, ProperHarmonicProduct};
use fs_math::dd::Dd;

mod sources;
use sources::{ExposureLine, ExposureLines, ExposureSide};

mod stress;
pub use stress::{
    AngularSpectrumExposureStressField, AngularSpectrumExposureStressLimits,
    AngularSpectrumExposureStressPoint,
};

mod momentum;
pub use momentum::{
    AngularSpectrumExposureMomentumField, AngularSpectrumExposureMomentumLimits,
    AngularSpectrumExposureMomentumPoint,
};

mod stress_momentum;
pub use stress_momentum::{
    AngularSpectrumExposureStressMomentumField, AngularSpectrumExposureStressMomentumLimits,
};

mod ensemble;
pub use ensemble::{
    AngularSpectrumExposureEnsemble, AngularSpectrumExposureEnsembleLimits,
    AngularSpectrumFactorStatistics,
};

mod layer;
pub use layer::{
    AngularSpectrumExposureHeatChannel, AngularSpectrumExposureHeatSpectrumLimits,
    AngularSpectrumExposureLayerDissipation, AngularSpectrumExposureLayerField,
    AngularSpectrumExposureLayerHeatLoadSpectrum, AngularSpectrumExposureLayerHeatSpectrum,
    AngularSpectrumExposureLayerLimits, AngularSpectrumHeatAreaQuadrature,
    AngularSpectrumHeatDepthSample, AngularSpectrumHeatVolumeQuadrature,
    AngularSpectrumIntegratedHeatLimits, PreparedAngularSpectrumExposureLayer,
};

/// Whole-line-bank bounds for a query. Pair work counts N*F² for N spatial
/// points and F distinct frequencies. Logical numerical storage includes one
/// original beam's complete workspace, all retained frequency fields, temporal
/// kernels, a compensated accumulation plane and the output. Other borrowed
/// beams, earlier outputs and allocator overhead are outside this bound.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureLimits {
    pub maximum_lines: usize,
    pub maximum_point_line_pairs: usize,
    pub maximum_point_frequency_pairs: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

/// Rectangular-exposure averages of the real physical fields. Electric and
/// magnetic squared components are mean squares in the named SI units.
/// Poynting flux is <E(t) cross H(t)> in W/m² and retains its physical sign;
/// no monochromatic complex/reactive flux is assigned to a finite exposure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposurePoint {
    pub electric_mean_squared_v2_m2: [f64; 3],
    pub magnetic_mean_squared_a2_m2: [f64; 3],
    pub poynting_w_m2: [f64; 3],
}

/// Exposure-averaged observables on the original periodic transverse cell.
#[derive(Debug)]
pub struct AngularSpectrumExposureField {
    grid: AngularSpectrumGrid,
    points: Vec<AngularSpectrumExposurePoint>,
    normal_power_w: f64,
}
impl AngularSpectrumExposureField {
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.grid
    }
    pub fn points(&self) -> &[AngularSpectrumExposurePoint] {
        &self.points
    }
    /// Mean signed normal power integrated over the continuous periodic cell,
    /// in watts. Modal orthogonality integrates the represented band-limited
    /// fields exactly in space, without aliasing their quadratic products on
    /// the original grid. This can differ from summing the point samples.
    pub fn normal_power_w(&self) -> f64 {
        self.normal_power_w
    }
}

/// Phase-related spectral lines through one physical stack. Lines may differ
/// in frequency, polarization and spatial profile, but must retain the same
/// original incident/exit/film material owners, layer order and thicknesses.
/// Equal-frequency lines are summed as fields before forming any products.
pub struct AngularSpectrumExposure<'b, 's, 'm> {
    anchor: &'b PlanarAngularSpectrum<'s, 'm>,
    lines: ExposureLines<'b, 's, 'm>,
    frequencies: Vec<f64>,
    source_indices: Vec<usize>,
    group_offsets: Vec<usize>,
    workspace_bytes: usize,
    field_fft_work: u64,
    point_frequency_pairs: usize,
}

impl<'b, 's, 'm> AngularSpectrumExposure<'b, 's, 'm> {
    /// The anchor fixes the physical stack/grid and adds no illumination unless
    /// it is listed as a line. An empty bank is darkness. Each supplied line
    /// includes its physical amplitude; no normalization or spectrum is inferred.
    /// O(L log L) preparation groups source indices by increasing frequency,
    /// preserving the original line order within every frequency.
    pub fn try_new(
        anchor: &'b PlanarAngularSpectrum<'s, 'm>,
        lines: &'b [PlanarAngularSpectrum<'s, 'm>],
        limits: AngularSpectrumExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        Self::new(
            anchor,
            ExposureLines::OneSided(lines),
            limits,
            &mut cancelled,
        )
    }

    /// Observe jointly phased inputs at both original physical faces.
    ///
    /// Each prepared illumination retains its left AND right peak amplitudes
    /// on one source clock. A proper covariance factor may therefore correlate
    /// faces, polarizations, spatial modes and frequencies together. The native
    /// two-port caches are borrowed and never rebuilt by this constructor.
    pub fn try_new_illumination(
        anchor: &'b PlanarAngularSpectrum<'s, 'm>,
        lines: &'b [PreparedAngularSpectrumIllumination<'b, 's, 'm>],
        limits: AngularSpectrumExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        Self::new(
            anchor,
            ExposureLines::TwoSided(lines),
            limits,
            &mut cancelled,
        )
    }

    fn new(
        anchor: &'b PlanarAngularSpectrum<'s, 'm>,
        lines: ExposureLines<'b, 's, 'm>,
        limits: AngularSpectrumExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if lines.len() > 4096 || lines.len() > limits.maximum_lines {
            return Err(bad("exposure line count exceeds its budget"));
        }
        // Reserve complete grouping capacity before any allocation. Offsets
        // retain L+1 slots even when several lines share a frequency.
        let retained_group_bytes = lines
            .len()
            .checked_mul(core::mem::size_of::<f64>() + 2 * core::mem::size_of::<usize>())
            .and_then(|v| v.checked_add(core::mem::size_of::<usize>()))
            .ok_or_else(|| bad("exposure workspace overflow"))?;
        let minimum_workspace = retained_group_bytes
            .checked_add(anchor.workspace_bytes())
            .ok_or_else(|| bad("exposure workspace overflow"))?;
        if minimum_workspace > limits.maximum_workspace_bytes {
            return Err(bad("exposure exceeds its workspace budget"));
        }
        let mut frequencies = Vec::new();
        frequencies
            .try_reserve_exact(lines.len())
            .map_err(|_| bad("exposure frequency allocation failed"))?;
        let mut source_indices = Vec::new();
        source_indices
            .try_reserve_exact(lines.len())
            .map_err(|_| bad("exposure source index allocation failed"))?;
        let mut group_offsets = Vec::new();
        group_offsets
            .try_reserve_exact(lines.len() + 1)
            .map_err(|_| bad("exposure group offset allocation failed"))?;
        let mut original_workspace = anchor.workspace_bytes();
        for index in 0..lines.len() {
            poll(&mut cancelled)?;
            let original = lines.get(index);
            let line = original.beam();
            if line.grid() != anchor.grid()
                || !line.source().same_sources_and_geometry(anchor.source())
            {
                return Err(bad(
                    "exposure lines require the same original materials, layers and grid",
                ));
            }
            source_indices.push(index);
            original_workspace = original_workspace.max(original.workspace_bytes());
        }
        // In-place sorting needs no extra index allocation. The secondary
        // index gives stable within-frequency accumulation without relying on
        // sorting stability. Admission caps this bounded sort at 4096 lines.
        source_indices.sort_unstable_by(|left, right| {
            lines
                .get(*left)
                .beam()
                .source()
                .angular_frequency_rad_s()
                .total_cmp(&lines.get(*right).beam().source().angular_frequency_rad_s())
                .then_with(|| left.cmp(right))
        });
        poll(&mut cancelled)?;
        for (offset, &index) in source_indices.iter().enumerate() {
            if offset % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let omega = lines.get(index).beam().source().angular_frequency_rad_s();
            if frequencies.last().copied() != Some(omega) {
                frequencies.push(omega);
                group_offsets.push(offset);
            }
        }
        group_offsets.push(lines.len());
        let n = anchor.modes.len();
        let f = frequencies.len();
        let line_pairs = n
            .checked_mul(lines.len())
            .ok_or_else(|| bad("exposure line work overflow"))?;
        let field_fft_work = anchor
            .field_fft_work()
            .checked_mul(f as u64)
            .ok_or_else(|| bad("exposure FFT work overflow"))?;
        let pairs = f
            .checked_mul(f)
            .ok_or_else(|| bad("exposure pair overflow"))?;
        let point_frequency_pairs = n
            .checked_mul(pairs)
            .ok_or_else(|| bad("exposure work overflow"))?;
        let workspace_bytes = n
            .checked_mul(core::mem::size_of::<[C64; 6]>())
            .and_then(|v| v.checked_mul(f))
            .and_then(|v| v.checked_add(f.checked_mul(core::mem::size_of::<Vec<[C64; 6]>>())?))
            .and_then(|v| {
                v.checked_add(pairs.checked_mul(core::mem::size_of::<HarmonicProduct>())?)
            })
            .and_then(|v| {
                v.checked_add(n.checked_mul(
                    core::mem::size_of::<[Dd; 12]>()
                        + core::mem::size_of::<AngularSpectrumExposurePoint>(),
                )?)
            })
            .and_then(|v| v.checked_add(retained_group_bytes))
            .and_then(|v| v.checked_add(original_workspace))
            .ok_or_else(|| bad("exposure workspace overflow"))?;
        if line_pairs > limits.maximum_point_line_pairs
            || point_frequency_pairs > limits.maximum_point_frequency_pairs
            || field_fft_work > limits.maximum_fft_work
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad("exposure exceeds its work or workspace budget"));
        }
        Ok(Self {
            anchor,
            lines,
            frequencies,
            source_indices,
            group_offsets,
            workspace_bytes,
            field_fft_work,
            point_frequency_pairs,
        })
    }

    pub fn line_count(&self) -> usize {
        self.lines.len()
    }
    pub fn frequency_count(&self) -> usize {
        self.frequencies.len()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn field_fft_work(&self) -> u64 {
        self.field_fft_work
    }
    pub fn point_frequency_pairs(&self) -> usize {
        self.point_frequency_pairs
    }

    /// Reconstruct each original line once, cache O(F²) temporal kernels, and
    /// integrate all vector cross terms. O(L*N + F*N*log N + F²*N) work; trigonometry
    /// for the temporal window is O(F²), independent of the grid size. There is
    /// no time stepping, time aliasing or repeated constitutive solve.
    pub fn field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureField, AngularSpectrumError> {
        self.field_on(
            ExposureSide::OneSided(side),
            distance_from_face_m,
            window,
            &mut cancelled,
        )
    }

    /// Finite-window fields of incoming, outgoing or total illumination at
    /// either original face. The right-face normal flux keeps its Cartesian
    /// sign; incoming right waves carry negative z flux. For a one-sided bank,
    /// incoming-right is darkness and total-right is the transmitted field.
    pub fn illumination_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureField, AngularSpectrumError> {
        self.field_on(
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            window,
            &mut cancelled,
        )
    }

    fn field_on(
        &self,
        side: ExposureSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !distance_from_face_m.is_finite() || distance_from_face_m < 0.0 {
            return Err(bad("field distance must be finite and nonnegative"));
        }
        let f = self.frequencies.len();
        let mut kernels = Vec::new();
        kernels
            .try_reserve_exact(f * f)
            .map_err(|_| bad("exposure kernel allocation failed"))?;
        for &a in &self.frequencies {
            for &b in &self.frequencies {
                poll(&mut cancelled)?;
                kernels.push(
                    window
                        .product_kernel(a, b)
                        .map_err(|_| bad("exposure temporal kernel is unrepresentable"))?,
                );
            }
        }
        let mut fields = Vec::new();
        fields
            .try_reserve_exact(f)
            .map_err(|_| bad("exposure field allocation failed"))?;
        for group in 0..f {
            fields.push(self.frequency_modal_field(
                group,
                side,
                distance_from_face_m,
                &mut cancelled,
            )?);
        }
        let normal_power_w = integrate_cell_power(
            self.anchor.grid(),
            &fields,
            &kernels,
            HarmonicProduct::mean_real_product_parts_dd,
            &mut cancelled,
        )?
        .to_f64();
        for (group, field) in fields.iter_mut().enumerate() {
            self.synthesize_frequency_field(group, field, &mut cancelled)?;
        }
        let n = self.anchor.modes.len();
        let mut points = Vec::new();
        points
            .try_reserve_exact(n)
            .map_err(|_| bad("exposure result allocation failed"))?;
        for point in 0..n {
            poll(&mut cancelled)?;
            let values = integrate_point(f, &kernels, |i| fields[i][point], &mut cancelled)
                .map_err(|error| match error {
                    MomentError::Cancelled => AngularSpectrumError::Cancelled,
                    MomentError::Product(_) => bad("exposure field product is unrepresentable"),
                    MomentError::Unrepresentable => {
                        bad("exposure second moment is unrepresentable")
                    }
                })?;
            points.push(AngularSpectrumExposurePoint {
                electric_mean_squared_v2_m2: [values[0], values[1], values[2]],
                magnetic_mean_squared_a2_m2: [values[3], values[4], values[5]],
                poynting_w_m2: [values[6], values[7], values[8]],
            });
        }
        let grid = self.anchor.grid();
        if !normal_power_w.is_finite() {
            return Err(bad("integrated exposure power is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(AngularSpectrumExposureField {
            grid,
            points,
            normal_power_w,
        })
    }

    // Sum equal-frequency sources BEFORE products, preserving complete nulls
    // and small residual fields instead of subtracting large source powers.
    fn frequency_field(
        &self,
        group: usize,
        side: ExposureSide,
        distance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[C64; 6]>, AngularSpectrumError> {
        self.frequency_field_with(group, cancelled, |line, cancelled| {
            line.modal_field(side, distance, cancelled)
        })
    }

    fn frequency_modal_field(
        &self,
        group: usize,
        side: ExposureSide,
        distance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[C64; 6]>, AngularSpectrumError> {
        self.frequency_modal_field_with(group, cancelled, |line, cancelled| {
            line.modal_field(side, distance, cancelled)
        })
    }

    fn frequency_field_with<F: FnMut() -> bool>(
        &self,
        group: usize,
        cancelled: &mut F,
        modal_field: impl FnMut(
            ExposureLine<'b, 's, 'm>,
            &mut F,
        ) -> Result<(Vec<[C64; 3]>, Vec<[C64; 3]>), AngularSpectrumError>,
    ) -> Result<Vec<[C64; 6]>, AngularSpectrumError> {
        let mut field = self.frequency_modal_field_with(group, cancelled, modal_field)?;
        self.synthesize_frequency_field(group, &mut field, cancelled)?;
        Ok(field)
    }

    fn frequency_modal_field_with<F: FnMut() -> bool>(
        &self,
        group: usize,
        cancelled: &mut F,
        mut modal_field: impl FnMut(
            ExposureLine<'b, 's, 'm>,
            &mut F,
        )
            -> Result<(Vec<[C64; 3]>, Vec<[C64; 3]>), AngularSpectrumError>,
    ) -> Result<Vec<[C64; 6]>, AngularSpectrumError> {
        let n = self.anchor.modes.len();
        let mut sums = Vec::new();
        sums.try_reserve_exact(n)
            .map_err(|_| bad("exposure sum allocation failed"))?;
        for i in 0..n {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            sums.push([Dd::ZERO; 12]);
        }
        let indices =
            &self.source_indices[self.group_offsets[group]..self.group_offsets[group + 1]];
        for &index in indices {
            poll(cancelled)?;
            let line = self.lines.get(index);
            let (electric, magnetic) = modal_field(line, &mut *cancelled)?;
            for (i, ((sum, e), h)) in sums.iter_mut().zip(&electric).zip(&magnetic).enumerate() {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                for (c, z) in e.iter().chain(h).enumerate() {
                    add(&mut sum[2 * c], z.re);
                    add(&mut sum[2 * c + 1], z.im);
                }
            }
        }
        let mut field = Vec::new();
        field
            .try_reserve_exact(n)
            .map_err(|_| bad("exposure group allocation failed"))?;
        for (i, sum) in sums.into_iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let value: [C64; 6] =
                std::array::from_fn(|c| C64::new(sum[2 * c].to_f64(), sum[2 * c + 1].to_f64()));
            if value.iter().any(|&z| !finite(z)) {
                return Err(bad("coherent exposure field is unrepresentable"));
            }
            field.push(value);
        }
        Ok(field)
    }

    // Replace modal vectors with their original samples in the same retained
    // frequency plane. The exact cell contraction runs first, so it needs no
    // second FFT or second bank of frequency fields.
    fn synthesize_frequency_field(
        &self,
        group: usize,
        field: &mut [[C64; 6]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        let mut electric = Vec::new();
        let mut magnetic = Vec::new();
        electric
            .try_reserve_exact(field.len())
            .map_err(|_| bad("exposure group allocation failed"))?;
        magnetic
            .try_reserve_exact(field.len())
            .map_err(|_| bad("exposure group allocation failed"))?;
        for (i, value) in field.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            electric.push([value[0], value[1], value[2]]);
            magnetic.push([value[3], value[4], value[5]]);
        }
        let source = self
            .lines
            .get(self.source_indices[self.group_offsets[group]])
            .beam();
        let synthesized = source.synthesize(electric, magnetic, &mut *cancelled)?;
        for (i, ((value, e), h)) in field
            .iter_mut()
            .zip(synthesized.electric_v_per_m)
            .zip(synthesized.magnetic_a_per_m)
            .enumerate()
        {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            *value = [e[0], e[1], e[2], h[0], h[1], h[2]];
        }
        Ok(())
    }
}

// One cached temporal kernel and one retained original six-component field
// per frequency. Stress, momentum and their joint observer share this owner.
// Modal contractions must precede synthesize(): the same storage then holds
// original grid samples, with exactly one FFT synthesis per frequency.
struct ExposureObservationBank<P> {
    kernels: Vec<P>,
    fields: Vec<Vec<[C64; 6]>>,
}

impl<P> ExposureObservationBank<P> {
    fn prepare(
        bank: &AngularSpectrumExposure<'_, '_, '_>,
        side: ExposureSide,
        distance: f64,
        kernel: impl Fn(f64, f64) -> Result<P, fs_fft::harmonic_window::WindowError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        let f = bank.frequency_count();
        let mut kernels = Vec::new();
        kernels
            .try_reserve_exact(f * f)
            .map_err(|_| bad("exposure observation kernel allocation failed"))?;
        for &a in &bank.frequencies {
            for &b in &bank.frequencies {
                poll(cancelled)?;
                kernels.push(
                    kernel(a, b).map_err(|_| {
                        bad("exposure observation temporal kernel is unrepresentable")
                    })?,
                );
            }
        }
        let mut fields = Vec::new();
        fields
            .try_reserve_exact(f)
            .map_err(|_| bad("exposure observation field allocation failed"))?;
        for group in 0..f {
            fields.push(bank.frequency_modal_field(group, side, distance, cancelled)?);
        }
        Ok(Self { kernels, fields })
    }

    fn synthesize(
        &mut self,
        bank: &AngularSpectrumExposure<'_, '_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        for (group, field) in self.fields.iter_mut().enumerate() {
            bank.synthesize_frequency_field(group, field, cancelled)?;
        }
        Ok(())
    }
}

// Integrate over the continuous periodic cell BEFORE sampling the quadratic
// product. The represented Fourier bins are signed integers, not integers
// modulo the FFT length: a negative Nyquist bin has no positive partner.
// For every nonzero pair k,-k, write the real field as C(t) cos(k.r) +
// S(t) sin(k.r), with complex C=A_k+A_-k and S=-i(A_k-A_-k). Orthogonality
// gives (mean(C_E*C_H)+mean(S_E*S_H))/2. This uses the original stable
// temporal Gram kernel, retaining subcycle nulls without subtracting two
// separately rounded sum/difference-frequency powers. Proper factors use
// the same spatial algebra with their original proper temporal kernel.
fn integrate_cell_power<P: Copy>(
    grid: AngularSpectrumGrid,
    fields: &[Vec<[C64; 6]>],
    kernels: &[P],
    product: impl Fn(P, [Dd; 2], [Dd; 2]) -> Result<Dd, fs_fft::harmonic_window::WindowError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, AngularSpectrumError> {
    let f = fields.len();
    let mut sum = Dd::ZERO;
    for mode in 0..grid.width * grid.height {
        if mode % 32 == 0 {
            poll(cancelled)?;
        }
        let opposite = opposite_mode(grid, mode);
        if opposite.is_some_and(|other| other < mode) {
            continue;
        }
        let constant = opposite == Some(mode);
        let partner = opposite.filter(|_| !constant);
        for i in 0..f {
            let a = spatial_quadratures(fields[i][mode], partner.map(|other| fields[i][other]));
            for j in 0..f {
                if j % 64 == 0 {
                    poll(cancelled)?;
                }
                let b = spatial_quadratures(fields[j][mode], partner.map(|other| fields[j][other]));
                let kernel = kernels[i * f + j];
                let mut flux = Dd::ZERO;
                for quadrature in 0..if constant { 1 } else { 2 } {
                    flux = flux
                        + product(kernel, a[quadrature][0], b[quadrature][4])
                            .map_err(|_| bad("modal cell power product is unrepresentable"))?
                        - product(kernel, a[quadrature][1], b[quadrature][3])
                            .map_err(|_| bad("modal cell power product is unrepresentable"))?;
                }
                sum = sum + flux * Dd::from_f64(if constant { 1.0 } else { 0.5 });
            }
        }
    }
    let area = Dd::from_f64(grid.width as f64 * grid.spacing_m[0])
        * Dd::from_f64(grid.height as f64 * grid.spacing_m[1]);
    let power = sum * area;
    if !power.hi.is_finite() || !power.lo.is_finite() {
        return Err(bad("integrated modal cell power is unrepresentable"));
    }
    poll(cancelled)?;
    Ok(power)
}

fn opposite_mode(grid: AngularSpectrumGrid, index: usize) -> Option<usize> {
    let opposite = |bin: usize, count: usize| {
        if bin == 0 {
            Some(0)
        } else if bin == count / 2 {
            None
        } else {
            Some(count - bin)
        }
    };
    Some(
        opposite(index / grid.width, grid.height)? * grid.width
            + opposite(index % grid.width, grid.width)?,
    )
}

fn spatial_quadratures(a: [C64; 6], b: Option<[C64; 6]>) -> [[[Dd; 2]; 6]; 2] {
    let b = b.unwrap_or([C64::ZERO; 6]);
    // Conjugate both amplitudes for the shared negative-time convention.
    let cosine = core::array::from_fn(|c| {
        [
            Dd::from_f64(a[c].re) + Dd::from_f64(b[c].re),
            -Dd::from_f64(a[c].im) - Dd::from_f64(b[c].im),
        ]
    });
    let sine = core::array::from_fn(|c| {
        [
            Dd::from_f64(a[c].im) - Dd::from_f64(b[c].im),
            Dd::from_f64(a[c].re) - Dd::from_f64(b[c].re),
        ]
    });
    [cosine, sine]
}

// Pure real-vector contraction, shared by the original-source exposure path
// and the analytic Maxwell-wave regressions below.
#[derive(Debug)]
enum MomentError {
    Cancelled,
    Product(fs_fft::harmonic_window::WindowError),
    Unrepresentable,
}

fn integrate_point(
    frequencies: usize,
    kernels: &[HarmonicProduct],
    field_at: impl Fn(usize) -> [C64; 6],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 9], MomentError> {
    integrate_point_dd(frequencies, kernels, field_at, product, cancelled)
        .map(|sum| sum.map(Dd::to_f64))
}

// Both deterministic means and proper covariance factors use exactly the
// same vector contraction. Keep the pair products and their sum in extended
// precision through the complete exposure/ensemble accumulation.
fn integrate_point_dd<P: Copy>(
    frequencies: usize,
    kernels: &[P],
    field_at: impl Fn(usize) -> [C64; 6],
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 9], MomentError> {
    let mut sum = [Dd::ZERO; 9];
    for i in 0..frequencies {
        let a = field_at(i);
        for j in 0..frequencies {
            if j % 64 == 0 && cancelled() {
                return Err(MomentError::Cancelled);
            }
            let b = field_at(j);
            let kernel = kernels[i * frequencies + j];
            for c in 0..3 {
                sum[c] = sum[c] + product(kernel, a[c], b[c])?;
                sum[c + 3] = sum[c + 3] + product(kernel, a[c + 3], b[c + 3])?;
                let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                sum[c + 6] = sum[c + 6] + product(kernel, a[u], b[v + 3])?;
                sum[c + 6] = sum[c + 6] - product(kernel, a[v], b[u + 3])?;
            }
        }
    }
    let values = sum.map(Dd::to_f64);
    if values.iter().any(|x| !x.is_finite()) || values[..6].iter().any(|x| *x < 0.0) {
        return Err(MomentError::Unrepresentable);
    }
    Ok(sum)
}

fn add(sum: &mut Dd, value: f64) {
    *sum = *sum + Dd::from_f64(value);
}
fn product(kernel: HarmonicProduct, a: C64, b: C64) -> Result<Dd, MomentError> {
    // The shared window uses exp(-i omega t); EM uses exp(+i omega t).
    kernel
        .mean_real_product_dd(a.conj(), b.conj())
        .map_err(MomentError::Product)
}

fn proper_product(kernel: ProperHarmonicProduct, a: C64, b: C64) -> Result<Dd, MomentError> {
    // Proper covariance removes sum-frequency terms by its authored zero
    // pseudocovariance, not by an optical-cycle or long-window approximation.
    kernel
        .mean_real_product_dd(a.conj(), b.conj())
        .map_err(MomentError::Product)
}

#[cfg(test)]
mod tests {
    use super::super::modes::fields;
    use super::*;

    fn wave(omega: f64, backwards: bool, transverse: [C64; 2]) -> [C64; 6] {
        let (e, h) = fields(
            0.0,
            0.0,
            C64::from_re(if backwards { -omega } else { omega }),
            omega,
            transverse,
        )
        .unwrap();
        [e[0], e[1], e[2], h[0], h[1], h[2]]
    }
    fn evaluate(waves: &[[C64; 6]], frequencies: &[f64], start: f64, duration: f64) -> [f64; 9] {
        let window = HarmonicWindow::new(start, duration).unwrap();
        let kernels = frequencies
            .iter()
            .flat_map(|&a| {
                frequencies
                    .iter()
                    .map(move |&b| window.product_kernel(a, b).unwrap())
            })
            .collect::<Vec<_>>();
        integrate_point(frequencies.len(), &kernels, |i| waves[i], &mut || false).unwrap()
    }

    fn evaluate_proper(
        waves: &[[C64; 6]],
        frequencies: &[f64],
        start: f64,
        duration: f64,
    ) -> [Dd; 9] {
        let window = HarmonicWindow::new(start, duration).unwrap();
        let kernels = frequencies
            .iter()
            .flat_map(|&a| {
                frequencies
                    .iter()
                    .map(move |&b| window.proper_product_kernel(a, b).unwrap())
            })
            .collect::<Vec<_>>();
        integrate_point_dd(
            frequencies.len(),
            &kernels,
            |i| waves[i],
            proper_product,
            &mut || false,
        )
        .unwrap()
    }

    #[test]
    fn g1_proper_vector_covariance_matches_authored_phase_ensemble_quadrature() {
        let frequencies = [4.0, 5.3];
        let waves = [
            wave(4.0, false, [C64::new(0.7, -0.3), C64::new(0.2, 0.5)]),
            wave(5.3, true, [C64::new(-0.4, 0.2), C64::new(0.1, -0.6)]),
        ];
        let (start, duration) = (0.17, 0.83);
        let actual = evaluate_proper(&waves, &frequencies, start, duration).map(Dd::to_f64);
        let mut expected = [0.0; 9];
        let n = 4096;
        // This explicitly authored four-state distribution has E[z]=E[z²]=0
        // and E[|z|²]=1. It is an independent second-moment oracle, not sampled
        // noise or an assumption about the production source distribution.
        for z in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
            for k in 0..=n {
                let t = start + duration * k as f64 / n as f64;
                let mut real = [0.0; 6];
                for (wave, omega) in waves.iter().zip(frequencies) {
                    let (s, c) = (omega * t).sin_cos();
                    for j in 0..6 {
                        let a = wave[j] * z;
                        real[j] += a.re * c - a.im * s;
                    }
                }
                let weight = if k == 0 || k == n {
                    1.0
                } else if k % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                for j in 0..3 {
                    expected[j] += weight * real[j] * real[j];
                    expected[j + 3] += weight * real[j + 3] * real[j + 3];
                    let (u, v) = ((j + 1) % 3, (j + 2) % 3);
                    expected[j + 6] += weight * (real[u] * real[v + 3] - real[v] * real[u + 3]);
                }
            }
        }
        for (actual, expected) in actual
            .into_iter()
            .zip(expected.map(|v| v / (12.0 * n as f64)))
        {
            assert!((actual - expected).abs() < 4e-13, "{actual} != {expected}");
        }
    }

    #[test]
    fn g1_proper_near_null_covariance_keeps_tiny_distinct_frequency_beats() {
        let frequencies = [10.0, 10.000000001];
        let waves = [
            wave(frequencies[0], false, [C64::from_re(1e6), C64::ZERO]),
            wave(frequencies[1], false, [C64::from_re(-1e6), C64::ZERO]),
        ];
        let actual = evaluate_proper(&waves, &frequencies, 0.3, 0.7)[0].to_f64();
        let mut expected = 0.0;
        let n = 1024;
        for k in 0..=n {
            let t = 0.3 + 0.7 * k as f64 / n as f64;
            // 1/2 |A exp(i wa t)-A exp(i wb t)|², evaluated without
            // subtracting near-identical powers or cosines.
            let sine = (0.5 * (frequencies[1] - frequencies[0]) * t).sin();
            let weight = if k == 0 || k == n {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            expected += weight * 2e12 * sine * sine;
        }
        expected /= 3.0 * n as f64;
        assert!(actual > 0.0);
        assert!(
            (actual / expected - 1.0).abs() < 2e-12,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g1_proper_covariance_is_invariant_to_unitary_factor_mixing() {
        let frequencies = [3.0, 5.0];
        let a = [
            wave(3.0, false, [C64::new(0.8, 0.3), C64::new(-0.2, 0.5)]),
            wave(5.0, true, [C64::new(-0.3, 0.1), C64::new(0.4, -0.7)]),
        ];
        let b = [
            wave(3.0, true, [C64::new(0.2, -0.1), C64::new(0.5, 0.3)]),
            wave(5.0, false, [C64::new(0.7, 0.6), C64::new(-0.1, 0.2)]),
        ];
        let i = C64::new(0.0, 1.0);
        let mixed_a: [[C64; 6]; 2] = std::array::from_fn(|r| {
            std::array::from_fn(|c| (a[r][c] + i * b[r][c]).scale(core::f64::consts::FRAC_1_SQRT_2))
        });
        let mixed_b: [[C64; 6]; 2] = std::array::from_fn(|r| {
            std::array::from_fn(|c| (i * a[r][c] + b[r][c]).scale(core::f64::consts::FRAC_1_SQRT_2))
        });
        let original_a = evaluate_proper(&a, &frequencies, 0.13, 0.71);
        let original_b = evaluate_proper(&b, &frequencies, 0.13, 0.71);
        let mixed_a = evaluate_proper(&mixed_a, &frequencies, 0.13, 0.71);
        let mixed_b = evaluate_proper(&mixed_b, &frequencies, 0.13, 0.71);
        for c in 0..9 {
            let original = (original_a[c] + original_b[c]).to_f64();
            let mixed = (mixed_a[c] + mixed_b[c]).to_f64();
            assert!((original - mixed).abs() < 5e-15, "{original} != {mixed}");
        }
        // A proper monochromatic factor has no carrier-phase dependence,
        // even when the shutter duration is far below one cycle.
        let sine = wave(10.0, true, [C64::new(0.0, 1.0), C64::ZERO]);
        let short = evaluate_proper(&[sine], &[10.0], 0.0, 1e-10).map(Dd::to_f64);
        assert_eq!(short[0], 0.5);
        assert_eq!(short[8], -0.5);
    }

    #[test]
    fn g1_real_maxwell_exposure_matches_direct_simpson_integration() {
        let frequencies = [4.0, 5.3];
        let waves = [
            wave(4.0, false, [C64::new(0.7, -0.3), C64::new(0.2, 0.5)]),
            wave(5.3, true, [C64::new(-0.4, 0.2), C64::new(0.1, -0.6)]),
        ];
        let (start, duration) = (0.17, 0.83);
        let actual = evaluate(&waves, &frequencies, start, duration);
        let mut expected = [0.0; 9];
        let n = 4096;
        for k in 0..=n {
            let t = start + duration * k as f64 / n as f64;
            let mut real = [0.0; 6];
            for (wave, omega) in waves.iter().zip(frequencies) {
                let (s, c) = (omega * t).sin_cos();
                for j in 0..6 {
                    real[j] += wave[j].re * c - wave[j].im * s;
                }
            }
            let weight = if k == 0 || k == n {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for j in 0..3 {
                expected[j] += weight * real[j] * real[j];
                expected[j + 3] += weight * real[j + 3] * real[j + 3];
                let (u, v) = ((j + 1) % 3, (j + 2) % 3);
                expected[j + 6] += weight * (real[u] * real[v + 3] - real[v] * real[u + 3]);
            }
        }
        for (actual, expected) in actual
            .into_iter()
            .zip(expected.map(|v| v / (3.0 * n as f64)))
        {
            assert!((actual - expected).abs() < 3e-13, "{actual} != {expected}");
        }
    }

    #[test]
    fn g1_distinct_frequency_cancellation_retains_weak_beating_power() {
        let frequencies = [10.0, 10.0001];
        let waves = [
            wave(frequencies[0], false, [C64::from_re(1e3), C64::ZERO]),
            wave(frequencies[1], false, [C64::from_re(-1e3), C64::ZERO]),
        ];
        let actual = evaluate(&waves, &frequencies, 0.3, 0.7);
        let mut expected = 0.0;
        let n = 8192;
        for k in 0..=n {
            let t = 0.3 + 0.7 * k as f64 / n as f64;
            // Independent stable identity for cos(wa*t)-cos(wb*t).
            let e = 2e3
                * (((frequencies[0] + frequencies[1]) * 0.5) * t).sin()
                * (((frequencies[1] - frequencies[0]) * 0.5) * t).sin();
            let weight = if k == 0 || k == n {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            expected += weight * e * e;
        }
        expected /= 3.0 * n as f64;
        assert!(expected > 0.0);
        assert!(
            (actual[0] / expected - 1.0).abs() < 2e-9,
            "{} != {expected}",
            actual[0]
        );
        // The original Maxwell reconstruction obtains H via k cross E / omega;
        // its rounded multiplication/division need not equal E bit for bit.
        for component in [4, 8] {
            assert!((actual[component] / actual[0] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn g1_short_exposures_keep_phase_polarization_and_reflected_flux_sign() {
        let sine = wave(10.0, false, [C64::new(0.0, 1.0), C64::ZERO]);
        let tiny = evaluate(&[sine], &[10.0], 0.0, 1e-10);
        assert!((tiny[0] / (1e-18 / 3.0) - 1.0).abs() < 1e-12);
        let reflected = [
            wave(10.0, true, [C64::ONE, C64::ZERO]),
            wave(11.0, true, [C64::ZERO, C64::ONE]),
        ];
        let long = evaluate(&reflected, &[10.0, 11.0], 0.0, core::f64::consts::TAU);
        assert!((long[8] + 1.0).abs() < 1e-14);
        let mut cancelled = || true;
        let kernel = HarmonicWindow::new(0.0, 1.0)
            .unwrap()
            .product_kernel(10.0, 10.0)
            .unwrap();
        assert!(matches!(
            integrate_point(1, &[kernel], |_| sine, &mut cancelled),
            Err(MomentError::Cancelled)
        ));
    }
}
