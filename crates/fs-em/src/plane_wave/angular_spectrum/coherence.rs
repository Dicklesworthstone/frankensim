//! Spatial partial coherence through the original vector diffraction operator.
//!
//! Each borrowed beam is one mutually uncorrelated peak-field factor B_m of
//! W(r,r') = sum_m B_m(r) B_m(r')*. Its amplitude includes its square-root
//! statistical weight. Within a factor, Fourier components interfere before
//! any quadratic observable is formed. Factors are then added as second
//! moments, never as fields. This is one monochromatic ensemble, not a PSD,
//! an inferred source bandwidth, or a random-phase approximation.

use super::*;
use crate::plane_wave::layers::{FIELD_PREPARATION_WORKSPACE_BYTES, PreparedLayerField};
use fs_math::dd::Dd;

/// Whole-factor-bank limits for one query. The conservative workspace includes
/// one original beam's complete workspace, accumulation/output planes and,
/// for an internal query, ONE shared layer cache. Other borrowed beam storage,
/// retained previous results and allocator overhead are outside this bound.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumCoherenceLimits {
    pub maximum_modes: usize,
    pub maximum_point_mode_pairs: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

/// Local time/ensemble second moments. The peak-phasor factor 1/2 is included
/// exactly once. Electric coherency is Hermitian: entry (i,j) is <E_i E_j*>;
/// its diagonal is electric RMS squared, not irradiance. Physical energy flux
/// is the real part of the complete complex Poynting vector <E cross H*>.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumPointStatistics {
    pub electric_coherency_v2_m2: [[C64; 3]; 3],
    pub magnetic_rms_squared_a2_m2: [f64; 3],
    pub complex_poynting_w_m2: [C64; 3],
}

/// Ensemble observables on the original periodic cell; no mean complex field
/// is invented for a partially coherent source.
#[derive(Debug)]
pub struct AngularSpectrumCoherenceField {
    grid: AngularSpectrumGrid,
    points: Vec<AngularSpectrumPointStatistics>,
    normal_power_w: f64,
}
impl AngularSpectrumCoherenceField {
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.grid
    }
    pub fn points(&self) -> &[AngularSpectrumPointStatistics] {
        &self.points
    }
    pub fn normal_power_w(&self) -> f64 {
        self.normal_power_w
    }
}

/// The same ensemble inside one authored film. Local positive heat is W/m³;
/// its transverse integral is W/m of depth. Neither is an extra thermal debit
/// in addition to the prepared layer's integrated absorbed watts.
#[derive(Debug)]
pub struct AngularSpectrumCoherenceLayerField {
    field: AngularSpectrumCoherenceField,
    absorbed_power_density_w_m3: Vec<f64>,
    absorbed_power_per_depth_w_m: f64,
}
impl AngularSpectrumCoherenceLayerField {
    pub fn field(&self) -> &AngularSpectrumCoherenceField {
        &self.field
    }
    pub fn absorbed_power_density_w_m3(&self) -> &[f64] {
        &self.absorbed_power_density_w_m3
    }
    pub fn absorbed_power_per_depth_w_m(&self) -> f64 {
        self.absorbed_power_per_depth_w_m
    }
}

/// Positive spatial-coherence factors sharing the exact original source and
/// grid. Preparation borrows existing beams without copying their mode arrays.
/// A field query streams one reconstructed factor at a time, O(N) additional
/// storage and O(M*N*log N) FFT work for M factors on N points.
pub struct AngularSpectrumCoherence<'b, 's, 'm> {
    anchor: &'b PlanarAngularSpectrum<'s, 'm>,
    factors: &'b [PlanarAngularSpectrum<'s, 'm>],
    power: AngularSpectrumPower,
    workspace_bytes: usize,
    field_fft_work: u64,
}

impl<'b, 's, 'm> AngularSpectrumCoherence<'b, 's, 'm> {
    /// Factors need not be orthogonal. Their covariance, not an eigenbasis,
    /// defines the source. The anchor fixes the original source/grid but adds
    /// no power unless listed among the factors. Empty factors mean darkness.
    pub fn try_new(
        anchor: &'b PlanarAngularSpectrum<'s, 'm>,
        factors: &'b [PlanarAngularSpectrum<'s, 'm>],
        limits: AngularSpectrumCoherenceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(&mut cancelled)?;
        let first = anchor;
        if factors.len() > 4096 || factors.len() > limits.maximum_modes {
            return Err(bad("coherence factor count exceeds its budget"));
        }
        let n = first.modes.len();
        let pairs = n
            .checked_mul(factors.len())
            .ok_or_else(|| bad("coherence work overflow"))?;
        let field_fft_work = first
            .field_fft_work
            .checked_mul(factors.len() as u64)
            .ok_or_else(|| bad("coherence FFT work overflow"))?;
        let per_point = core::mem::size_of::<PointSum>()
            + core::mem::size_of::<AngularSpectrumPointStatistics>()
            + 2 * core::mem::size_of::<f64>();
        let workspace_bytes = n
            .checked_mul(per_point)
            .and_then(|v| v.checked_add(first.workspace_bytes))
            .ok_or_else(|| bad("coherence workspace overflow"))?;
        if pairs > limits.maximum_point_mode_pairs
            || field_fft_work > limits.maximum_fft_work
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad("coherence exceeds its work or workspace budget"));
        }
        let mut sums = [Dd::ZERO; 4];
        for factor in factors {
            poll(&mut cancelled)?;
            if !core::ptr::eq(factor.source(), first.source()) || factor.grid() != first.grid() {
                return Err(bad(
                    "coherence factors require the same original source and grid",
                ));
            }
            let p = factor.power();
            for (sum, value) in sums.iter_mut().zip([
                p.incident_w,
                p.reflected_w,
                p.transmitted_w,
                p.film_absorbed_w,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let [incident_w, reflected_w, transmitted_w, film_absorbed_w] = sums.map(Dd::to_f64);
        if [incident_w, reflected_w, transmitted_w, film_absorbed_w]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(bad("coherence power is unrepresentable"));
        }
        Ok(Self {
            anchor,
            factors,
            power: AngularSpectrumPower {
                incident_w,
                reflected_w,
                transmitted_w,
                film_absorbed_w,
            },
            workspace_bytes,
            field_fft_work,
        })
    }

    pub fn source(&self) -> &'s PlanarDielectricLayers<'m> {
        self.anchor.source()
    }
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.anchor.grid()
    }
    pub fn factor_count(&self) -> usize {
        self.factors.len()
    }
    pub fn power(&self) -> AngularSpectrumPower {
        self.power
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn field_fft_work(&self) -> u64 {
        self.field_fft_work
    }

    /// Net ensemble Maxwell force [N] in the original beam coordinates.
    /// Original coherent factors interfere internally, then their forces add
    /// as second moments. The common source/grid and preparation work limits
    /// are unchanged. O(factors*points), allocation-free, without any FFT.
    /// See `PlanarAngularSpectrum::radiation_force` for the exterior stress
    /// admission, periodic-cell geometry and mechanical no-claim boundaries.
    pub fn radiation_force(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3], AngularSpectrumError> {
        poll(&mut cancelled)?;
        self.source().admit_refracting_momentum_exterior()?;
        let mut sums = [Dd::ZERO; 3];
        for factor in self.factors {
            let force = factor.radiation_force(&mut cancelled)?;
            for (sum, value) in sums.iter_mut().zip(force) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let force = sums.map(Dd::to_f64);
        if force.iter().any(|f| !f.is_finite()) {
            return Err(bad("coherence radiation force is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(force)
    }

    /// Complex fields interfere within each source factor through the existing
    /// six-component Maxwell reconstruction. Their second moments then add.
    pub fn field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumCoherenceField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !distance_from_face_m.is_finite() || distance_from_face_m < 0.0 {
            return Err(bad("field distance must be finite and nonnegative"));
        }
        let mut sums = point_sums(self.anchor.modes.len(), &mut cancelled)?;
        for factor in self.factors {
            let field = factor.field(side, distance_from_face_m, &mut cancelled)?;
            accumulate(&mut sums, &field, None, &mut cancelled)?;
        }
        Ok(finish(self.grid(), sums, false, &mut cancelled)?.field)
    }

    /// Prepare identical TE/TM layer responses ONCE for the entire bank. The
    /// layer/work limits cover that shared preparation and the whole ensemble
    /// field query; no per-factor constitutive solve or cache is added.
    pub fn prepare_layer(
        &self,
        layer_index: usize,
        limits: AngularSpectrumLayerLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumCoherenceLayer<'_, 'b, 's, 'm>, AngularSpectrumError> {
        poll(&mut cancelled)?;
        let n = self.anchor.modes.len();
        let workspace_bytes = n
            .checked_mul(core::mem::size_of::<[PreparedLayerField<'_, '_>; 2]>())
            .and_then(|v| v.checked_add(self.workspace_bytes))
            .and_then(|v| v.checked_add(FIELD_PREPARATION_WORKSPACE_BYTES))
            .ok_or_else(|| bad("coherence layer workspace overflow"))?;
        if workspace_bytes > limits.maximum_workspace_bytes
            || self.field_fft_work > limits.maximum_fft_work
        {
            return Err(bad("coherence layer exceeds its work or workspace budget"));
        }
        let beam: &'b PlanarAngularSpectrum<'s, 'm> = self.anchor;
        let layer = beam.prepare_layer(layer_index, limits, &mut cancelled)?;
        let mut absorbed = Dd::ZERO;
        for factor in self.factors {
            poll(&mut cancelled)?;
            absorbed =
                absorbed + Dd::from_f64(layer.absorbed_power_for_beam(factor, &mut cancelled)?);
        }
        let absorbed_power_w = absorbed.to_f64();
        if !absorbed_power_w.is_finite() || absorbed_power_w < 0.0 {
            return Err(bad("coherence film heat is unrepresentable"));
        }
        Ok(PreparedAngularSpectrumCoherenceLayer {
            bank: self,
            layer,
            absorbed_power_w,
            workspace_bytes,
        })
    }
}

pub struct PreparedAngularSpectrumCoherenceLayer<'c, 'b, 's, 'm> {
    bank: &'c AngularSpectrumCoherence<'b, 's, 'm>,
    layer: PreparedAngularSpectrumLayer<'b, 's, 'm>,
    absorbed_power_w: f64,
    workspace_bytes: usize,
}
impl PreparedAngularSpectrumCoherenceLayer<'_, '_, '_, '_> {
    pub fn layer_index(&self) -> usize {
        self.layer.layer_index()
    }
    pub fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn preparation_layer_mode_pairs(&self) -> usize {
        self.layer.preparation_layer_mode_pairs()
    }
    pub fn field(
        &self,
        depth_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumCoherenceLayerField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        let thickness = self
            .bank
            .source()
            .layer(self.layer_index())
            .unwrap()
            .thickness_m;
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > thickness {
            return Err(bad("field depth is outside the selected layer"));
        }
        let mut sums = point_sums(self.bank.anchor.modes.len(), &mut cancelled)?;
        for factor in self.bank.factors {
            let field = self.layer.field_for_beam(factor, depth_m, &mut cancelled)?;
            accumulate(
                &mut sums,
                field.field(),
                Some(field.absorbed_power_density_w_m3()),
                &mut cancelled,
            )?;
        }
        finish(self.bank.grid(), sums, true, &mut cancelled)
    }
}

#[derive(Clone, Copy)]
struct PointSum {
    electric_diagonal: [Dd; 3],
    electric_cross: [[Dd; 2]; 3],
    magnetic_diagonal: [Dd; 3],
    poynting: [[Dd; 2]; 3],
    heat: Dd,
}
impl PointSum {
    const ZERO: Self = Self {
        electric_diagonal: [Dd::ZERO; 3],
        electric_cross: [[Dd::ZERO; 2]; 3],
        magnetic_diagonal: [Dd::ZERO; 3],
        poynting: [[Dd::ZERO; 2]; 3],
        heat: Dd::ZERO,
    };
}
const PAIRS: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];

fn point_sums(
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<PointSum>, AngularSpectrumError> {
    let mut sums = Vec::new();
    sums.try_reserve_exact(n)
        .map_err(|_| bad("coherence accumulation allocation failed"))?;
    for i in 0..n {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        sums.push(PointSum::ZERO);
    }
    Ok(sums)
}
fn add(sum: &mut Dd, value: f64) -> Result<(), AngularSpectrumError> {
    if !value.is_finite() {
        return Err(bad("coherence second moment is unrepresentable"));
    }
    *sum = *sum + Dd::from_f64(value);
    Ok(())
}
fn add_complex(sum: &mut [Dd; 2], value: C64) -> Result<(), AngularSpectrumError> {
    add(&mut sum[0], value.re)?;
    add(&mut sum[1], value.im)
}
fn accumulate(
    sums: &mut [PointSum],
    field: &AngularSpectrumField,
    heat: Option<&[f64]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    for (i, ((sum, &e), &h)) in sums
        .iter_mut()
        .zip(field.electric_v_per_m())
        .zip(field.magnetic_a_per_m())
        .enumerate()
    {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        for c in 0..3 {
            add(&mut sum.electric_diagonal[c], 0.5 * e[c].norm_sq())?;
            add(&mut sum.magnetic_diagonal[c], 0.5 * h[c].norm_sq())?;
            let a = (c + 1) % 3;
            let b = (c + 2) % 3;
            add_complex(
                &mut sum.poynting[c],
                e[a].scale(0.5) * h[b].conj() - e[b].scale(0.5) * h[a].conj(),
            )?;
        }
        for (slot, (a, b)) in PAIRS.into_iter().enumerate() {
            add_complex(&mut sum.electric_cross[slot], e[a].scale(0.5) * e[b].conj())?;
        }
        if let Some(heat) = heat {
            add(&mut sum.heat, heat[i])?;
        }
    }
    Ok(())
}
fn finish(
    grid: AngularSpectrumGrid,
    sums: Vec<PointSum>,
    with_heat: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumCoherenceLayerField, AngularSpectrumError> {
    let mut points = Vec::new();
    points
        .try_reserve_exact(sums.len())
        .map_err(|_| bad("coherence result allocation failed"))?;
    let mut heat = Vec::new();
    if with_heat {
        heat.try_reserve_exact(sums.len())
            .map_err(|_| bad("coherence heat allocation failed"))?;
    }
    let mut power = Dd::ZERO;
    let mut total_heat = Dd::ZERO;
    for (i, sum) in sums.into_iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        let mut electric = [[C64::ZERO; 3]; 3];
        for c in 0..3 {
            electric[c][c] = C64::from_re(sum.electric_diagonal[c].to_f64());
        }
        for (slot, (a, b)) in PAIRS.into_iter().enumerate() {
            electric[a][b] = C64::new(
                sum.electric_cross[slot][0].to_f64(),
                sum.electric_cross[slot][1].to_f64(),
            );
            electric[b][a] = electric[a][b].conj();
        }
        let magnetic = sum.magnetic_diagonal.map(Dd::to_f64);
        let poynting = sum.poynting.map(|z| C64::new(z[0].to_f64(), z[1].to_f64()));
        if electric
            .iter()
            .flatten()
            .chain(&poynting)
            .any(|&z| !finite(z))
            || magnetic.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(bad("coherence result is unrepresentable"));
        }
        add(&mut power, poynting[2].re)?;
        if with_heat {
            let q = sum.heat.to_f64();
            if q < 0.0 {
                return Err(bad("coherence heat must remain nonnegative"));
            }
            add(&mut total_heat, q)?;
            heat.push(q);
        }
        points.push(AngularSpectrumPointStatistics {
            electric_coherency_v2_m2: electric,
            magnetic_rms_squared_a2_m2: magnetic,
            complex_poynting_w_m2: poynting,
        });
    }
    let pixel_area = grid.spacing_m[0] * grid.spacing_m[1];
    let normal_power_w = power.to_f64() * pixel_area;
    let absorbed_power_per_depth_w_m = total_heat.to_f64() * pixel_area;
    if !normal_power_w.is_finite() || !absorbed_power_per_depth_w_m.is_finite() {
        return Err(bad("integrated coherence result is unrepresentable"));
    }
    poll(cancelled)?;
    Ok(AngularSpectrumCoherenceLayerField {
        field: AngularSpectrumCoherenceField {
            grid,
            points,
            normal_power_w,
        },
        absorbed_power_density_w_m3: heat,
        absorbed_power_per_depth_w_m,
    })
}
