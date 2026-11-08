//! Direct volume integral of the four coherent crystal eigenwaves.
//!
//! Every pair uses a decaying exponential integral anchored at the endpoint
//! where its magnitude is largest. Forward/backward interference therefore
//! needs neither a growing transfer factor nor a subtraction of port powers.
//! Short intervals sum the field series before squaring, retaining heat near
//! standing-wave nodes. The original positive dielectric-loss channels are
//! scaled independently so a strong lossless polarization cannot erase a weak
//! dissipative component.

use super::super::kernel::{
    absorption::{
        attenuated, cross_integral, decay_integral, dot, multiply, physical_power, product,
        short_mean,
    },
    slab_volume_m3,
};
use super::{
    AbsorbingSlabResponse, AbsorbingSlabWave, OrientedUniaxialSlab, PreparedAbsorbingSlabField,
    ReflectionError,
};
use fs_math::{c64::C64, dd::Dd, decay, det};

type Complex = [Dd; 2];
const ZERO: Complex = [Dd::ZERO; 2];
type Result<T> = core::result::Result<T, &'static str>;

/// Dielectric heat deposited in a finite prism inside the original film.
/// The area is measured parallel to the two faces. Both incident ports share
/// one transverse wavevector, so cycle-mean heat is uniform in that plane;
/// all four internal waves still interfere through depth. Material and illumination
/// are frozen at the prepared source snapshot, not evolved thermally here.
pub struct AbsorbingSlabVolumeAbsorption<'w, 'f, 's, 'm> {
    prepared_field: &'w PreparedAbsorbingSlabField<'f, 's, 'm>,
    area_m2: f64,
    start_depth_m: f64,
    end_depth_m: f64,
    volume_m3: f64,
    absorbed_power_w: f64,
}

impl<'w, 'f, 's, 'm> AbsorbingSlabVolumeAbsorption<'w, 'f, 's, 'm> {
    pub const fn prepared_field(&self) -> &'w PreparedAbsorbingSlabField<'f, 's, 'm> {
        self.prepared_field
    }
    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
        self.prepared_field.source()
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn start_depth_m(&self) -> f64 {
        self.start_depth_m
    }
    pub const fn end_depth_m(&self) -> f64 {
        self.end_depth_m
    }
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Integral of omega/2 E* [-Im(epsilon)] E over the declared prism [W],
    /// including every ordinary/extraordinary and forward/backward cross term.
    pub const fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
}

impl<'f, 's, 'm> PreparedAbsorbingSlabField<'f, 's, 'm> {
    /// Integrate the original tensor dielectric work over positive face area
    /// and `0 <= start_depth_m < end_depth_m <= thickness_m`. Fixed-size
    /// analytic exponential integrals retain weak loss and coherent nodes
    /// without depth quadrature, another boundary solve or material queries.
    pub fn slab_absorption(
        &self,
        area_m2: f64,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> core::result::Result<AbsorbingSlabVolumeAbsorption<'_, 'f, 's, 'm>, ReflectionError> {
        let volume_m3 = admit_volume(&self.response, area_m2, start_depth_m, end_depth_m)
            .map_err(ReflectionError::Invalid)?;
        let absorbed_power_w = integrated_heat(&self.response, area_m2, start_depth_m, end_depth_m)
            .map_err(ReflectionError::Invalid)?;
        Ok(AbsorbingSlabVolumeAbsorption {
            prepared_field: self,
            area_m2,
            start_depth_m,
            end_depth_m,
            volume_m3,
            absorbed_power_w,
        })
    }
}

fn admit_volume(response: &AbsorbingSlabResponse, area: f64, start: f64, end: f64) -> Result<f64> {
    if end > response.thickness_m {
        return Err("absorbing film control volume extends beyond the original crystal");
    }
    slab_volume_m3(area, start, end)
}

#[derive(Clone, Copy)]
struct Geometry {
    normal: C64,
    anchor: Dd,
    orientation: f64,
    alpha: f64,
}

/// Sibling consumers must retain the same prepared source and admit volume.
pub(super) fn integrated_heat(
    response: &AbsorbingSlabResponse,
    area: f64,
    start: f64,
    end: f64,
) -> Result<f64> {
    admit_volume(response, area, start, end)?;
    integrated_layer_heat(
        LayerWaves {
            thickness_m: response.thickness_m,
            omega: response.omega,
            axis: response.axis,
            loss_f_m: response.loss_f_m,
            forward: &response.forward,
            backward: &response.backward,
        },
        area,
        start,
        end,
    )
}

/// Shared physical layer view: all four waves retain their two face anchors.
/// A stack layer and a finite film use exactly the same tensor-work integral.
pub(in crate::plane_wave::uniaxial::oriented::absorbing) struct LayerWaves<'a> {
    pub thickness_m: f64,
    pub omega: f64,
    pub axis: [f64; 3],
    pub loss_f_m: [f64; 2],
    pub forward: &'a [AbsorbingSlabWave; 2],
    pub backward: &'a [AbsorbingSlabWave; 2],
}

pub(in crate::plane_wave::uniaxial::oriented::absorbing) fn integrated_layer_heat(
    layer: LayerWaves<'_>,
    area: f64,
    start: f64,
    end: f64,
) -> Result<f64> {
    slab_volume_m3(area, start, end)?;
    if end > layer.thickness_m {
        return Err("absorbing control volume extends beyond the original layer");
    }
    if layer.loss_f_m == [0.0; 2] {
        return Ok(0.0);
    }
    let start = Dd::from_f64(start);
    let end = Dd::from_f64(end);
    let waves = [
        layer.forward[0],
        layer.forward[1],
        layer.backward[0],
        layer.backward[1],
    ];
    let geometry: [Geometry; 4] = core::array::from_fn(|m| {
        let orientation = if m < 2 { 1.0 } else { -1.0 };
        Geometry {
            normal: waves[m].normal_wavenumber_rad_m,
            anchor: if m < 2 {
                Dd::ZERO
            } else {
                Dd::from_f64(layer.thickness_m)
            },
            orientation,
            alpha: -orientation * waves[m].normal_wavenumber_rad_m.im,
        }
    });
    if geometry
        .iter()
        .any(|g| !g.alpha.is_finite() || g.alpha < 0.0 || !g.normal.re.is_finite())
    {
        return Err("absorbing film volume requires resolved passive eigenwaves");
    }
    let mut electric = core::array::from_fn::<_, 4, _>(|m| {
        let g = geometry[m];
        let near = if g.orientation > 0.0 { start } else { end };
        attenuated(
            waves[m].electric_v_m,
            g.alpha,
            (near - g.anchor) * Dd::from_f64(g.orientation),
        )
    });
    // Identical anchored modes share their complete exponential. Sum them
    // before projection/integration instead of cancelling separate powers.
    for m in 0..4 {
        for j in 0..m {
            if geometry[m].normal == geometry[j].normal && geometry[m].anchor == geometry[j].anchor
            {
                for k in 0..3 {
                    for part in 0..2 {
                        electric[j][k][part] = electric[j][k][part] + electric[m][k][part];
                    }
                }
                electric[m] = [ZERO; 3];
                break;
            }
        }
    }
    let axis = layer.axis.map(Dd::from_f64);
    let norm = axis.iter().fold(Dd::ZERO, |s, &v| s + v * v);
    let projected = electric.map(|e| {
        let parallel = [0, 1].map(|part| (0..3).fold(Dd::ZERO, |s, j| s + e[j][part] * axis[j]));
        core::array::from_fn::<_, 4, _>(|j| {
            if j == 3 {
                parallel
            } else {
                [0, 1].map(|part| e[j][part] - parallel[part] * axis[j] / norm)
            }
        })
    });
    let mut power = Dd::ZERO;
    for channel in 0..4 {
        let loss = if channel == 3 {
            layer.loss_f_m[1] / norm.to_f64()
        } else {
            layer.loss_f_m[0]
        };
        if loss == 0.0 {
            continue;
        }
        let mut modes = projected.map(|e| e[channel]);
        let field_scale = modes
            .iter()
            .flatten()
            .fold(0.0_f64, |s, v| s.max(v.to_f64().abs()));
        if !field_scale.is_finite()
            || modes
                .iter()
                .flatten()
                .any(|v| !v.hi.is_finite() || !v.lo.is_finite())
        {
            return Err("absorbing film loss-channel field is unrepresentable");
        }
        if field_scale == 0.0 {
            continue;
        }
        let scale = Dd::from_f64(field_scale);
        for mode in &mut modes {
            *mode = mode.map(|v| v / scale);
        }
        let (integral, length) = channel_integral(modes, geometry, start, end)?;
        let integral = integral.to_f64();
        if !integral.is_finite() || integral < 0.0 {
            return Err("coherent film absorption integral is unresolved");
        }
        power = power
            + Dd::from_f64(physical_power(
                integral,
                length,
                loss,
                layer.omega,
                area,
                field_scale,
            ));
    }
    let power = power.to_f64();
    if !power.is_finite() || power < 0.0 {
        return Err("absorbing film volume power is unrepresentable");
    }
    Ok(power)
}

fn channel_integral(
    mut modes: [Complex; 4],
    geometry: [Geometry; 4],
    start: Dd,
    end: Dd,
) -> Result<(Dd, f64)> {
    let width = end - start;
    let active = modes.map(|v| v != ZERO);
    let mut reference = 0;
    for m in 1..4 {
        if (dot(modes[m], modes[m]) - dot(modes[reference], modes[reference])).to_f64() > 0.0 {
            reference = m;
        }
    }
    if active.iter().filter(|&&v| v).count() == 1 {
        return Ok((
            dot(modes[reference], modes[reference])
                * decay_integral(geometry[reference].alpha, width),
            1.0,
        ));
    }
    let exposures = geometry.map(|g| product(g.alpha, width));
    let beats =
        geometry.map(|g| Dd::from_f64(geometry[reference].normal.re) - Dd::from_f64(g.normal.re));
    let short = (0..4).all(|m| {
        !active[m] || {
            let phase = beats[m] * width;
            exposures[m] < 0.125
                && phase.hi.is_finite()
                && phase.lo.is_finite()
                && phase.to_f64().abs() < 0.125
        }
    });
    if short {
        let mut rates = [ZERO; 4];
        for m in 0..4 {
            if !active[m] {
                continue;
            }
            // A backward wave grows only over this small local interval;
            // start it at its correctly attenuated entrance-side value.
            let remaining = if geometry[m].orientation < 0.0 {
                Dd::ONE + Dd::from_f64(det::expm1(-exposures[m]))
            } else {
                Dd::ONE
            };
            let phase = relative_phase(geometry[reference], geometry[m], start)?;
            modes[m] = multiply(modes[m], phase).map(|v| v * remaining);
            rates[m] = [
                -Dd::from_f64(geometry[m].orientation * geometry[m].alpha) * width,
                beats[m] * width,
            ];
        }
        return Ok((short_mean(modes, rates), width.to_f64()));
    }
    let mut integral = Dd::ZERO;
    for m in 0..4 {
        if !active[m] {
            continue;
        }
        integral = integral + dot(modes[m], modes[m]) * decay_integral(geometry[m].alpha, width);
        for j in 0..m {
            if !active[j] {
                continue;
            }
            let a = geometry[j];
            let b = geometry[m];
            let (orientation, rates, exposure) = if a.orientation == b.orientation {
                (a.orientation, [a.alpha, b.alpha], 0.0)
            } else {
                let orientation = if a.alpha >= b.alpha {
                    a.orientation
                } else {
                    b.orientation
                };
                (
                    orientation,
                    [(a.alpha - b.alpha).abs(), 0.0],
                    product(a.alpha.min(b.alpha), width),
                )
            };
            let remaining = decay::remaining_fraction(exposure);
            if remaining == 0.0 {
                continue;
            }
            let factor = if exposure < 0.5 {
                Dd::ONE + Dd::from_f64(det::expm1(-exposure))
            } else {
                Dd::from_f64(remaining)
            };
            let endpoint = if orientation > 0.0 { start } else { end };
            let phase = relative_phase(a, b, endpoint)?;
            let overlap = [
                dot(modes[j], modes[m]),
                modes[j][0] * modes[m][1] - modes[j][1] * modes[m][0],
            ];
            let beat =
                (Dd::from_f64(a.normal.re) - Dd::from_f64(b.normal.re)) * Dd::from_f64(orientation);
            let weight = cross_integral(rates, beat, width)?;
            let cross = multiply(multiply(overlap, phase), weight)[0];
            integral = integral + cross * factor * Dd::from_f64(2.0);
        }
    }
    Ok((integral, 1.0))
}

fn relative_phase(a: Geometry, b: Geometry, position: Dd) -> Result<Complex> {
    let beat = Dd::from_f64(a.normal.re) - Dd::from_f64(b.normal.re);
    let da = position - a.anchor;
    let db = position - b.anchor;
    let phase = if a.anchor == b.anchor {
        beat * da
    } else if da.to_f64().abs() <= db.to_f64().abs() {
        beat * da + Dd::from_f64(b.normal.re) * (b.anchor - a.anchor)
    } else {
        beat * db + Dd::from_f64(a.normal.re) * (b.anchor - a.anchor)
    };
    if !phase.hi.is_finite() || !phase.lo.is_finite() {
        return Err("absorbing film relative anchored phase is unrepresentable");
    }
    let (sine, cosine) = det::sin_cos_dd(phase);
    Ok([cosine, sine])
}

#[cfg(test)]
mod tests;
