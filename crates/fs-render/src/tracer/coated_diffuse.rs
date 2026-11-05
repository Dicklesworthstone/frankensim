//! Mixed-measure transport for the original scalar tracer's smooth coated base.
//! The shared coating owns physical reflection. This adapter only samples its
//! mirror mass and continuous lobe, with matched forward/reverse MIS densities.

use super::{
    DielectricGlass, TracerError, Vec3, cosine_sample, fresnel_dielectric,
    medium_ior, smooth_conductor_direction,
};
use crate::coated_diffuse::SmoothDiffuseCoating;

/// Both branches reflect into the same physical exterior. A constant coat
/// cannot silently replace a dispersive surrounding phase index.
pub(super) fn admit_medium(
    coating: SmoothDiffuseCoating,
    incident: Option<DielectricGlass>,
    wavelength_nm: f64,
) -> Result<(), TracerError> {
    let actual = medium_ior(incident, wavelength_nm)?;
    let declared = coating.exterior_index();
    if (actual - declared).abs() > 8.0 * f64::EPSILON * actual.max(declared) {
        return Err(TracerError::InvalidInput);
    }
    Ok(())
}

/// Unconditional solid-angle density, including the probability of reaching
/// the diffuse branch. Its specular delta contributes no continuous density.
pub(super) fn pdf(
    coating: SmoothDiffuseCoating,
    normal: Vec3,
    outgoing: Vec3,
    incoming: Vec3,
) -> f64 {
    let cos_o = normal.dot(outgoing);
    let cos_i = normal.dot(incoming);
    if cos_o <= 0.0 || cos_i <= 0.0 {
        return 0.0;
    }
    coating.branch_probabilities(cos_o.clamp(0.0, 1.0))
        .map_or(0.0, |mass| mass[1] * cos_i / core::f64::consts::PI)
}

pub(super) struct Sample {
    pub direction: Vec3,
    pub delta: bool,
    pub pdf_fwd: f64,
    pub pdf_rev: f64,
    // Mirror R / proposal mass, or entry T / diffuse mass times exit T.
    coefficient: f64,
}

impl Sample {
    /// The substrate spectrum is evaluated at the actual path wavelength by
    /// the caller. No RGB projection or wavelength-independent albedo enters.
    pub fn weight(
        &self,
        coating: SmoothDiffuseCoating,
        substrate_reflectance: f64,
    ) -> Result<f64, TracerError> {
        if self.delta {
            Ok(self.coefficient)
        } else {
            Ok(self.coefficient * coating.diffuse_weight(substrate_reflectance)?)
        }
    }
}

/// Reuse the existing two BSDF uniforms: the first selects a branch and is
/// remapped within the diffuse interval. The second supplies its azimuth.
/// No new random dimension or change to legacy material paths is required.
pub(super) fn sample(
    coating: SmoothDiffuseCoating,
    normal: Vec3,
    outgoing: Vec3,
    u1: f64,
    u2: f64,
) -> Result<Option<Sample>, TracerError> {
    if !(0.0..1.0).contains(&u1) || !(0.0..1.0).contains(&u2) {
        return Err(TracerError::InvalidInput);
    }
    let cosine = normal.dot(outgoing);
    if cosine <= 0.0 {
        return Ok(None);
    }
    let interface = fresnel_dielectric(
        cosine.clamp(0.0, 1.0), coating.exterior_index(), coating.coating_index(),
    )?;
    let [mirror_mass, diffuse_mass] = interface.event_probabilities();
    // Place the smaller branch beside zero. Subtracting a mirror mass rounded
    // to one would erase still-representable grazing transmission.
    let mirror = if mirror_mass <= diffuse_mass {
        u1 < mirror_mass
    } else {
        u1 >= diffuse_mass
    };
    if mirror || diffuse_mass == 0.0 {
        return Ok(Some(Sample {
            direction: smooth_conductor_direction(normal, outgoing)?,
            delta: true,
            pdf_fwd: 0.0,
            pdf_rev: 0.0,
            coefficient: interface.reflectance / mirror_mass,
        }));
    }
    // Only the remapping's last-ULP roundoff is bounded; physical parameters
    // and energy are never clipped. The authored uniforms remain half-open.
    let conditional = if mirror_mass <= diffuse_mass {
        (u1 - mirror_mass) / diffuse_mass
    } else {
        u1 / diffuse_mass
    };
    let remapped = conditional.clamp(0.0, 1.0_f64.next_down());
    let (direction, conditional_pdf) = cosine_sample(normal, remapped, u2);
    if conditional_pdf <= 0.0 {
        return Ok(None);
    }
    let exit = fresnel_dielectric(
        normal.dot(direction).clamp(0.0, 1.0),
        coating.exterior_index(), coating.coating_index(),
    )?;
    Ok(Some(Sample {
        direction,
        delta: false,
        pdf_fwd: diffuse_mass * normal.dot(direction) / core::f64::consts::PI,
        pdf_rev: exit.event_probabilities()[1] * cosine / core::f64::consts::PI,
        coefficient: (interface.transmittance / diffuse_mass) * exit.transmittance,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectral::SpectralReflectance;
    use crate::tracer::{Material, opaque_bsdf_eval, bsdf_pdf};

    fn material(rho: f64, coating: SmoothDiffuseCoating) -> Material {
        Material::CoatedDiffuse {
            coating,
            reflectance: SpectralReflectance::constant(rho).unwrap(),
        }
    }

    fn mirror_draw(mass: [f64; 2]) -> f64 {
        if mass[0] <= mass[1] { 0.5 * mass[0] } else { mass[1] + 0.5 * mass[0] }
    }

    fn diffuse_draw(mass: [f64; 2], unit: f64) -> f64 {
        if mass[0] <= mass[1] { mass[0] + mass[1] * unit } else { mass[1] * unit }
    }

    #[test]
    fn weak_grazing_diffuse_branch_survives_unit_rounded_mirror_mass() {
        let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let outgoing = Vec3::new(1.0, 0.0, 1e-20);
        let mass = coating.branch_probabilities(outgoing.z).unwrap();
        assert_eq!(mass[0], 1.0);
        assert!(mass[1] > 0.0);
        let event = sample(coating, normal, outgoing, 0.5 * mass[1], 0.3)
            .unwrap().unwrap();
        assert!(!event.delta);
        assert!(event.pdf_fwd > 0.0);
        assert!(event.weight(coating, 0.8).unwrap() > 0.0);
    }

    #[test]
    fn scalar_mixed_sample_and_connection_pdf_share_the_same_physical_brdf() {
        let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        for cosine in [0.03, 0.4, 1.0] {
            let outgoing = Vec3::new((1.0_f64 - cosine * cosine).sqrt(), 0.0, cosine);
            let mass = coating.branch_probabilities(cosine).unwrap();
            let mirror = sample(coating, normal, outgoing, mirror_draw(mass), 0.3)
                .unwrap().unwrap();
            assert!(mirror.delta);
            assert_eq!((mirror.pdf_fwd, mirror.pdf_rev), (0.0, 0.0));
            assert!((mirror.direction.z - cosine).abs() < 1e-14);
            for rho in [0.0, 0.2, 0.8, 1.0] {
                let mat = material(rho, coating);
                for unit in [0.0, 0.2, 0.8, 0.99] {
                    let draw = diffuse_draw(mass, unit);
                    let event = sample(coating, normal, outgoing, draw, 0.31)
                        .unwrap().unwrap();
                    assert!(!event.delta);
                    let value = opaque_bsdf_eval(
                        &mat, normal, outgoing, event.direction, 550.0, None, None,
                    ).unwrap();
                    let reverse = opaque_bsdf_eval(
                        &mat, normal, event.direction, outgoing, 550.0, None, None,
                    ).unwrap();
                    let expected = value * normal.dot(event.direction) / event.pdf_fwd;
                    assert!((event.weight(coating, rho).unwrap() - expected).abs() < 2e-13);
                    assert!((value - reverse).abs() < 2e-13);
                    assert_eq!(event.pdf_fwd, bsdf_pdf(&mat, normal, outgoing, event.direction, None));
                    assert_eq!(event.pdf_rev, bsdf_pdf(&mat, normal, event.direction, outgoing, None));
                }
            }
        }
    }

    #[test]
    fn scalar_white_furnace_and_black_base_retain_both_measures() {
        const N: usize = 32768;
        let normal = Vec3::new(0.0, 0.0, 1.0);
        for cosine in [0.1, 0.5, 1.0] {
            let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
            let outgoing = Vec3::new((1.0_f64 - cosine * cosine).sqrt(), 0.0, cosine);
            let mass = coating.branch_probabilities(cosine).unwrap();
            let mirror = sample(coating, normal, outgoing, mirror_draw(mass), 0.3)
                .unwrap().unwrap();
            for rho in [0.0, 0.3, 1.0] {
                let mut diffuse = 0.0;
                for index in 0..N {
                    let conditional = (index as f64 + 0.5) / N as f64;
                    let event = sample(
                        coating, normal, outgoing,
                        diffuse_draw(mass, conditional), 0.37,
                    ).unwrap().unwrap();
                    diffuse += event.weight(coating, rho).unwrap();
                }
                let actual = mass[0] * mirror.weight(coating, rho).unwrap()
                    + mass[1] * diffuse / N as f64;
                let expected = coating.directional_reflectance(cosine, rho).unwrap();
                assert!((actual - expected).abs() < 4e-7, "{actual} != {expected}");
            }
        }
    }

    #[test]
    fn scalar_coating_binds_medium_wavelength_and_material_identity() {
        let matched = SmoothDiffuseCoating::try_new(1.0, 1.0).unwrap();
        let coated = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
        let immersed = SmoothDiffuseCoating::try_new(1.33, 1.5).unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        assert!(admit_medium(immersed, None, 550.0).is_err());
        for draw in [0.0, 0.3, 0.999] {
            let event = sample(matched, normal, normal, draw, 0.2).unwrap().unwrap();
            assert!(!event.delta);
            assert_eq!(event.weight(matched, 0.6).unwrap(), 0.6);
        }
        assert_ne!(material(0.6, matched).content_identity(), material(0.6, coated).content_identity());
        assert_ne!(material(0.6, coated).content_identity(), material(0.7, coated).content_identity());
        for draw in [-0.1, 1.0, f64::NAN] {
            assert!(sample(coated, normal, normal, draw, 0.2).is_err());
        }
    }
}
