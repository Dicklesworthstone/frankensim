//! Positive-factor electromagnetic statistics and original dielectric work.
use super::{AbsorbingCrystalResponse, FieldSum, OrientedCrystalFieldSide};
use fs_math::{c64::C64, dd::Dd};

#[derive(Debug)]
pub(in super::super) enum IlluminationError {
    Cancelled,
    Electromagnetic(&'static str),
}
impl From<&'static str> for IlluminationError {
    fn from(message: &'static str) -> Self {
        Self::Electromagnetic(message)
    }
}
type Result<T> = core::result::Result<T, IlluminationError>;

/// Absolute normal powers, never normalized fractions of an empty source.
/// Reflected power is the nonnegative outward magnitude. Transmitted power
/// enters the crystal at z=0, before the subsequent volumetric absorption.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingCrystalIlluminationPower {
    pub incident_w_m2: f64,
    pub reflected_w_m2: f64,
    pub transmitted_w_m2: f64,
}

/// Cartesian cycle/ensemble second moments. Electric entry (i,j) is
/// sum_m E_m,i conj(E_m,j)/2: its diagonal is RMS electric field squared,
/// not irradiance. The peak-phasor factor 1/2 is included exactly once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingCrystalStatistics {
    pub electric_coherency_v2_m2: [[C64; 3]; 3],
    pub magnetic_rms_squared_a2_m2: [f64; 3],
    pub poynting_w_m2: [f64; 3],
    pub absorbed_power_density_w_m3: f64,
}

pub(in super::super) fn poll_illumination(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(IlluminationError::Cancelled)
    } else {
        Ok(())
    }
}

pub(in super::super) fn illumination_power<'a>(
    responses: impl IntoIterator<Item = &'a AbsorbingCrystalResponse>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AbsorbingCrystalIlluminationPower> {
    poll_illumination(&mut cancelled)?;
    let mut sums = [Dd::ZERO; 3];
    for response in responses {
        poll_illumination(&mut cancelled)?;
        let incoming = response.incident_normal_power_w_m2;
        for (sum, value) in sums.iter_mut().zip([
            incoming,
            incoming * response.reflected_power_fraction,
            incoming * response.transmitted_power_fraction,
        ]) {
            *sum = *sum + Dd::from_f64(value);
        }
    }
    let [incident_w_m2, reflected_w_m2, transmitted_w_m2] = sums.map(Dd::to_f64);
    if [incident_w_m2, reflected_w_m2, transmitted_w_m2]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err("crystal illumination boundary power is unrepresentable".into());
    }
    poll_illumination(&mut cancelled)?;
    Ok(AbsorbingCrystalIlluminationPower {
        incident_w_m2,
        reflected_w_m2,
        transmitted_w_m2,
    })
}

pub(in super::super) fn illumination_field_statistics<'a>(
    responses: impl IntoIterator<Item = &'a AbsorbingCrystalResponse>,
    side: OrientedCrystalFieldSide,
    point: [f64; 3],
    mut cancelled: impl FnMut() -> bool,
) -> Result<AbsorbingCrystalStatistics> {
    poll_illumination(&mut cancelled)?;
    FieldSum::new(side, point)?;
    let mut electric = [[[Dd::ZERO; 2]; 3]; 3];
    let mut magnetic = [Dd::ZERO; 3];
    let mut poynting = [Dd::ZERO; 3];
    let mut heat = Dd::ZERO;
    for response in responses {
        poll_illumination(&mut cancelled)?;
        let value = response.field(side, point)?;
        for i in 0..3 {
            for j in i..3 {
                let product =
                    half_product(value.field.electric_v_m[i], value.field.electric_v_m[j]);
                for part in 0..2 {
                    electric[i][j][part] = electric[i][j][part] + product[part];
                }
            }
            let h = value.field.magnetic_a_m[i];
            magnetic[i] = magnetic[i] + half_product(h, h)[0];
            poynting[i] = poynting[i] + Dd::from_f64(value.field.poynting_w_m2[i]);
        }
        heat = heat + Dd::from_f64(value.absorbed_power_density_w_m3);
    }
    let mut electric_coherency_v2_m2 = [[C64::ZERO; 3]; 3];
    for i in 0..3 {
        for j in i..3 {
            let [re, im] = electric[i][j].map(Dd::to_f64);
            let value = C64::new(re, if i == j { 0.0 } else { im });
            if !re.is_finite() || !im.is_finite() || i == j && re < 0.0 {
                return Err("crystal electric coherency is unrepresentable".into());
            }
            electric_coherency_v2_m2[i][j] = value;
            if i != j {
                electric_coherency_v2_m2[j][i] = value.conj();
            }
        }
    }
    let magnetic_rms_squared_a2_m2 = magnetic.map(Dd::to_f64);
    let poynting_w_m2 = poynting.map(Dd::to_f64);
    let absorbed_power_density_w_m3 = heat.to_f64();
    if magnetic_rms_squared_a2_m2
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        || poynting_w_m2.iter().any(|v| !v.is_finite())
        || !absorbed_power_density_w_m3.is_finite()
        || absorbed_power_density_w_m3 < 0.0
    {
        return Err("crystal field statistics are unrepresentable".into());
    }
    poll_illumination(&mut cancelled)?;
    Ok(AbsorbingCrystalStatistics {
        electric_coherency_v2_m2,
        magnetic_rms_squared_a2_m2,
        poynting_w_m2,
        absorbed_power_density_w_m3,
    })
}

pub(in super::super) fn illumination_poynting<'a>(
    responses: impl IntoIterator<Item = &'a AbsorbingCrystalResponse>,
    side: OrientedCrystalFieldSide,
    point: [f64; 3],
    mut cancelled: impl FnMut() -> bool,
) -> Result<[f64; 3]> {
    poll_illumination(&mut cancelled)?;
    FieldSum::new(side, point)?;
    let mut sum = [Dd::ZERO; 3];
    for response in responses {
        poll_illumination(&mut cancelled)?;
        let value = response.electromagnetic_field(side, point)?.poynting_w_m2;
        for (sum, value) in sum.iter_mut().zip(value) {
            *sum = *sum + Dd::from_f64(value);
        }
    }
    let sum = sum.map(Dd::to_f64);
    if sum.iter().any(|v| !v.is_finite()) {
        return Err("crystal illumination Poynting vector is unrepresentable".into());
    }
    poll_illumination(&mut cancelled)?;
    Ok(sum)
}

pub(in super::super) fn illumination_absorption<'a>(
    responses: impl IntoIterator<Item = &'a AbsorbingCrystalResponse>,
    area: f64,
    start: f64,
    end: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(f64, f64)> {
    poll_illumination(&mut cancelled)?;
    let volume = super::slab_volume_m3(area, start, end)?;
    let mut heat = Dd::ZERO;
    for response in responses {
        poll_illumination(&mut cancelled)?;
        heat = heat + Dd::from_f64(response.slab_absorbed_power_w(area, start, end)?);
    }
    let heat = heat.to_f64();
    if !heat.is_finite() || heat < 0.0 {
        return Err("crystal illumination absorbed power is unrepresentable".into());
    }
    poll_illumination(&mut cancelled)?;
    Ok((volume, heat))
}

fn half_product(a: C64, b: C64) -> [Dd; 2] {
    let [mut ar, mut ai, mut br, mut bi] = [a.re, a.im, b.re, b.im].map(Dd::from_f64);
    let half = Dd::from_f64(0.5);
    // Half before the quadratic product avoids overflowing a representable
    // RMS moment. Halve the larger field so a tiny partner is not erased.
    if a.re.abs().max(a.im.abs()) >= b.re.abs().max(b.im.abs()) {
        ar = ar * half;
        ai = ai * half;
    } else {
        br = br * half;
        bi = bi * half;
    }
    [ar * br + ai * bi, ai * br - ar * bi]
}
