//! Wood moisture: the equilibrium moisture content a declared environment
//! gives, and the clear-wood stiffness and density at that moisture content,
//! in the printed forms of the FPL Wood Handbook (FPL-GTR-282). Every
//! coefficient is a claim of a material card; this module holds only the forms.
//!
//! - Equation (4-5), the Hailwood-Horrobin EMC fitted by Simpson (1973), with
//!   W, K, K1 and K2 quadratic in the Celsius temperature
//!   (pack `wood-emc-hailwood-horrobin-fpl-gtr282`, valid over Table 4-2).
//! - Equation (5-3), P = P12 (P12/Pg)^((12 - M)/(Mp - 12)) below Mp, applied
//!   with the bending modulus's 12 % and green values.
//! - Equation (4-9), shrinkage S_x = S_0 (1 - x/MCfs) below the fiber
//!   saturation point, which with the (1 + M) moisture mass gives density.
//!
//! Moisture contents are fractions (0.12 is 12 %), as the pack compiler
//! stores a `%` claim. The elastic RATIOS of a reference-ratio card are not
//! moisture laws in the source: a caller that scales every modulus by
//! [`WoodMoistureLaw::modulus_factor`] holds them at their 12 % values, an
//! assumption it must state.
use crate::MaterialError;
use fs_matdb::{MaterialCard, QueryPoint, SelectionPolicy};
use fs_math::det;
use fs_qty::Dims;

/// The FPL tables' reference moisture content, 12 %.
pub const REFERENCE_MOISTURE_CONTENT: f64 = 0.12;
const PER_KELVIN: Dims = Dims([0, 0, 0, -1, 0, 0]);
const PER_KELVIN_SQUARED: Dims = Dims([0, 0, 0, -2, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn claim(
    card: &MaterialCard,
    name: &str,
    point: &QueryPoint,
    dims: Dims,
) -> Result<f64, MaterialError> {
    let answer = card
        .claims()
        .query(name, point, SelectionPolicy::SingleClaimOnly)
        .map_err(|e| MaterialError::Parameters {
            what: format!("{name}: {e:?}"),
        })?;
    let sample = &answer.evidence.value;
    if sample.dims != dims || !sample.value.is_finite() {
        return Err(MaterialError::Parameters {
            what: format!("{name}: dims {:?}, value {}", sample.dims.0, sample.value),
        });
    }
    Ok(sample.value)
}

/// Equation (4-5): the equilibrium moisture content (fraction) of wood at
/// `temperature_k` and `relative_humidity` (0..1), from a card carrying the
/// Hailwood-Horrobin coefficients. Outside the card's validity it refuses.
///
/// # Errors
/// [`MaterialError::Parameters`] for a missing or mis-dimensioned
/// coefficient, a state outside the card's validity, or a nonphysical result.
pub fn equilibrium_moisture_content(
    card: &MaterialCard,
    temperature_k: f64,
    relative_humidity: f64,
) -> Result<f64, MaterialError> {
    let point = QueryPoint::new()
        .with("temperature", temperature_k)
        .and_then(|p| p.with("relative_humidity", relative_humidity))
        .map_err(|e| MaterialError::Parameters {
            what: format!("EMC state: {e:?}"),
        })?;
    let quadratic = |name: &str| -> Result<[f64; 3], MaterialError> {
        let stem = format!("wood_emc_hailwood_horrobin_{name}");
        Ok([
            claim(card, &format!("{stem}_constant"), &point, Dims::NONE)?,
            claim(card, &format!("{stem}_per_celsius"), &point, PER_KELVIN)?,
            claim(
                card,
                &format!("{stem}_per_celsius_squared"),
                &point,
                PER_KELVIN_SQUARED,
            )?,
        ])
    };
    let numerator = claim(
        card,
        "wood_emc_hailwood_horrobin_numerator",
        &point,
        Dims::NONE,
    )?;
    // The printed polynomials take the Celsius temperature.
    let t = temperature_k - 273.15;
    let at = |c: [f64; 3]| c[0] + c[1] * t + c[2] * t * t;
    let (w, k, k1, k2) = (
        at(quadratic("w")?),
        at(quadratic("k")?),
        at(quadratic("k1")?),
        at(quadratic("k2")?),
    );
    let kh = k * relative_humidity;
    let percent = numerator / w
        * (kh / (1.0 - kh)
            + (k1 * kh + 2.0 * k1 * k2 * kh * kh) / (1.0 + k1 * kh + k1 * k2 * kh * kh));
    if !(percent.is_finite() && percent >= 0.0 && kh < 1.0) {
        return Err(MaterialError::Parameters {
            what: format!(
                "EMC at {temperature_k} K, RH {relative_humidity} is not physical ({percent} %)"
            ),
        });
    }
    Ok(percent / 100.0)
}

/// One species' moisture laws, from its card's claims at a moisture content.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WoodMoistureLaw {
    /// Bending modulus at 12 % [Pa], the P12 of Equation (5-3).
    pub bending_modulus_12_pa: f64,
    /// Green bending modulus [Pa], the Pg of Equation (5-3).
    pub bending_modulus_green_pa: f64,
    /// Mp, where the green and dry branches of Equation (5-3) meet (fraction).
    pub intersection_moisture_content: f64,
    /// S_0, volumetric shrinkage green to ovendry (fraction of green volume).
    pub volumetric_shrinkage: f64,
    /// MCfs, the fiber saturation point of Equation (4-9) (fraction).
    pub fiber_saturation: f64,
    /// The moisture content the laws were admitted at (fraction).
    pub moisture_content: f64,
}

impl WoodMoistureLaw {
    /// The laws of `card` at `moisture_content`. The 12 % bending modulus is
    /// read at 12 %; the law parameters at `moisture_content`, so a moisture
    /// content outside their stated domain (below Mp, below MCfs) refuses.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for a missing claim, a moisture content
    /// outside a claim's validity, or inconsistent parameters.
    pub fn from_card(card: &MaterialCard, moisture_content: f64) -> Result<Self, MaterialError> {
        let at = |m: f64| {
            QueryPoint::new()
                .with("moisture_content", m)
                .map_err(|e| MaterialError::Parameters {
                    what: format!("moisture state: {e:?}"),
                })
        };
        let (reference, here) = (at(REFERENCE_MOISTURE_CONTENT)?, at(moisture_content)?);
        let law = Self {
            bending_modulus_12_pa: claim(
                card,
                "modulus_of_elasticity_bending",
                &reference,
                PRESSURE,
            )?,
            bending_modulus_green_pa: claim(
                card,
                "moisture_law_green_bending_modulus",
                &here,
                PRESSURE,
            )?,
            intersection_moisture_content: claim(
                card,
                "moisture_law_intersection_moisture_content",
                &here,
                Dims::NONE,
            )?,
            volumetric_shrinkage: claim(
                card,
                "volumetric_shrinkage_green_to_ovendry",
                &here,
                Dims::NONE,
            )?,
            fiber_saturation: claim(card, "fiber_saturation_moisture_content", &here, Dims::NONE)?,
            moisture_content,
        };
        let consistent = law.bending_modulus_green_pa > 0.0
            && law.bending_modulus_12_pa > 0.0
            && law.intersection_moisture_content > REFERENCE_MOISTURE_CONTENT
            && (0.0..1.0).contains(&law.volumetric_shrinkage)
            && law.fiber_saturation > REFERENCE_MOISTURE_CONTENT
            && (0.0..law.intersection_moisture_content.min(law.fiber_saturation))
                .contains(&moisture_content);
        if !consistent {
            return Err(MaterialError::Parameters {
                what: format!("inconsistent wood moisture laws {law:?}"),
            });
        }
        Ok(law)
    }

    /// Equation (5-3) as a factor on a 12 % property: (P12/Pg)^((12 - M)/(Mp - 12)).
    #[must_use]
    pub fn modulus_factor(&self) -> f64 {
        let exponent = (REFERENCE_MOISTURE_CONTENT - self.moisture_content)
            / (self.intersection_moisture_content - REFERENCE_MOISTURE_CONTENT);
        det::pow(
            self.bending_modulus_12_pa / self.bending_modulus_green_pa,
            exponent,
        )
    }

    /// Density at M over density at 12 %: moisture mass (1 + M) over the
    /// Equation (4-9) volume, V(M) = V_green (1 - S_0 (1 - M/MCfs)).
    #[must_use]
    pub fn density_factor(&self) -> f64 {
        let volume = |m: f64| 1.0 - self.volumetric_shrinkage * (1.0 - m / self.fiber_saturation);
        (1.0 + self.moisture_content) / (1.0 + REFERENCE_MOISTURE_CONTENT)
            * volume(REFERENCE_MOISTURE_CONTENT)
            / volume(self.moisture_content)
    }
}

/// Equation (5-3) on explicit values (moisture contents as fractions), for a
/// property with no card: P12 (P12/Pg)^((0.12 - M)/(Mp - 0.12)).
#[must_use]
pub fn moisture_adjusted_property(
    p12: f64,
    green: f64,
    intersection: f64,
    moisture_content: f64,
) -> f64 {
    p12 * det::pow(
        p12 / green,
        (REFERENCE_MOISTURE_CONTENT - moisture_content)
            / (intersection - REFERENCE_MOISTURE_CONTENT),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// FPL's own worked example under Equation (5-3): white ash modulus of
    /// rupture at 8 %, P12 = 103,000 kPa, Pg = 66,000 kPa, Mp = 24 %:
    /// "P8 = 103,000 (103,000/66,000)^(4/12) = 119,500 kPa".
    #[test]
    fn equation_5_3_reproduces_the_handbook_example() {
        let p8 = moisture_adjusted_property(103_000.0, 66_000.0, 0.24, 0.08);
        assert!((p8 - 119_500.0).abs() < 50.0, "{p8}");
        assert_eq!(
            moisture_adjusted_property(103_000.0, 66_000.0, 0.24, 0.12),
            103_000.0
        );
    }
}
