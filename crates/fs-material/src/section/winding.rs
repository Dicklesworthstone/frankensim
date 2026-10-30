//! A close-coiled round-wire winding is not a bonded cylindrical sleeve.
//!
//! The frozen helix has centre-line radius `R`, wire radius `a`, and axial
//! pitch `p`. Its length per member length is `hypot(2*pi*R, p)/p`, so its
//! mass is the wire's actual density times area times that length. The
//! close-coiled spring law gives `EA_wrap = G*(2a)^4*p/(8*(2R)^3)`.
//! See J. A. Kemp (2017), PLOS ONE 12:e0184803, equations 8--19,
//! <https://doi.org/10.1371/journal.pone.0184803>.
//!
//! This is an explicit leading-order sliding-winding reduction: the core
//! supplies flexure and its selected constitutive loss; the winding contributes mass and
//! elastic axial compliance. Core and winding share end displacements while
//! remaining free to slip along their interface; a free-ended winding requires
//! a different axial boundary condition. It has no winding/core friction, winding loss,
//! torsion, wrapping prestress or resolved three-dimensional contact. It is
//! not an acoustically qualified replacement for a measured wound specimen.
//! Both materials and geometry survive compilation. Small wire/coil ratio,
//! small helix angle and wavelengths long compared with the pitch are required.

use super::{ConcentricElasticSection, ConcentricKelvinVoigtSection, invalid, positive};
use super::relaxation::ConcentricPronySection;
use crate::state_point::{DENSITY_PROPERTY, MaterialStatePointError, ResolvedMaterialStatePoint};
use fs_qty::{Density, Pressure, QuantitySpec};

/// Single, non-interpenetrating round-wire helix on a circular core.
/// Dimensions describe the already loaded, frozen specimen in metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RoundWireWinding {
    /// Radius of the individual winding wire, not the finished member.
    pub wire_radius_m: f64,
    /// Axial advance of one complete turn.
    pub pitch_m: f64,
}

/// The original core's selected bending law. The winding does not inherit it.
/// A Prony spectrum retains each layer's independent material history; it is
/// never converted to a frequency-sampled Kelvin-Voigt viscosity.
#[derive(Clone, Debug, PartialEq)]
pub enum WoundCoreSection {
    KelvinVoigt(ConcentricKelvinVoigtSection),
    Prony(ConcentricPronySection),
}
impl WoundCoreSection {
    /// Common source geometry and equilibrium section integrals.
    #[must_use]
    pub const fn elastic(&self) -> &ConcentricElasticSection {
        match self { Self::KelvinVoigt(s) => s.elastic(), Self::Prony(s) => s.elastic() }
    }
    /// Complete support of the selected core law, before wrapping coefficients.
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        match self { Self::KelvinVoigt(s) => s.omega_band_rad_s(), Self::Prony(s) => s.omega_band_rad_s() }
    }
}
impl From<ConcentricKelvinVoigtSection> for WoundCoreSection {
    fn from(section: ConcentricKelvinVoigtSection) -> Self { Self::KelvinVoigt(section) }
}
impl From<ConcentricPronySection> for WoundCoreSection {
    fn from(section: ConcentricPronySection) -> Self { Self::Prony(section) }
}

/// Source-retaining close-coiled spring winding and core-dominated flexure.
#[derive(Clone, Debug, PartialEq)]
pub struct CloseWoundSection {
    core: WoundCoreSection,
    winding: RoundWireWinding,
    material: ResolvedMaterialStatePoint,
    shear_property: String,
    helix_radius_m: f64,
    wire_length_per_member_length: f64,
    winding_line_density_kg_m: f64,
    winding_axial_stiffness_n: f64,
    line_density_kg_m: f64,
    axial_stiffness_n: f64,
    omega_band_rad_s: (f64, f64),
}

impl CloseWoundSection {
    /// Bind actual winding density and shear modulus to the original core.
    /// The wire/coil radius ratio must be at most 1/4, and pitch/circumference
    /// at most 0.15. These are conservative model-admission assumptions, not
    /// experimental error bounds. Both materials must share the same physical
    /// query coordinates. The winding remains lossless in this reduction;
    /// its missing friction is never replaced with the core's viscosity.
    pub fn try_new(
        core: impl Into<WoundCoreSection>,
        winding: RoundWireWinding,
        material: ResolvedMaterialStatePoint,
        shear_property: &str,
    ) -> Result<Self, MaterialStatePointError> {
        let core = core.into();
        let a = positive(winding.wire_radius_m, "positive winding wire radius")?;
        let pitch = positive(winding.pitch_m, "positive winding pitch")?;
        let radius = positive(core.elastic().outer_radius_m() + a, "winding centre-line radius")?;
        let circumference = positive(std::f64::consts::TAU * radius, "winding circumference")?;
        if a / radius > 0.25 || pitch / circumference > 0.15 {
            return Err(invalid("close-coiled reduction requires wire/coil radius <= 1/4 and pitch/circumference <= 0.15"));
        }
        // A pitch of exactly one wire diameter slightly overlaps circular
        // tubes: closest points occur just before a complete turn. Solve that
        // geometric minimum, rather than checking only equal azimuths.
        let advance = pitch / std::f64::consts::TAU;
        let mut lo = 1.5 * std::f64::consts::PI;
        let mut hi = std::f64::consts::TAU;
        for _ in 0..64 {
            let angle = 0.5 * (lo + hi);
            if angle.sin() + (advance / radius).powi(2) * angle < 0.0 {
                lo = angle;
            } else {
                hi = angle;
            }
        }
        let turn = 0.5 * (lo + hi);
        let separation = (2.0 * radius * (0.5 * turn).sin()).hypot(advance * turn);
        if !separation.is_finite() || separation < 2.0 * a {
            return Err(invalid("round-wire winding turns interpenetrate at the declared pitch"));
        }
        let mut band = core.omega_band_rad_s();
        for layer in core.elastic().layers() {
            if layer.material.query_point() != material.query_point()
                || layer.material.axis_quantities() != material.axis_quantities() {
                return Err(invalid("core and winding must resolve at identical physical coordinates"));
            }
        }
        let mut coefficients = [0.0; 2];
        for (i, (name, dims)) in [(DENSITY_PROPERTY, Density::DIMS), (shear_property, Pressure::DIMS)]
            .into_iter().enumerate()
        {
            let property = material.property(name)
                .ok_or_else(|| invalid("winding needs unambiguous density and shear modulus"))?;
            if property.requirement().quantity() != QuantitySpec::dimensional(dims)
                || !property.is_constant_at_fixed_environment() {
                return Err(invalid("winding coefficients require dimensioned laws at a held typed environment"));
            }
            coefficients[i] = positive(property.value_si(), "positive winding coefficient")?;
            if let Some((lower, upper)) = property.angular_frequency_band_rad_s()? {
                band = (band.0.max(lower), band.1.min(upper));
            }
        }
        if band.0 > band.1 {
            return Err(invalid("core and winding frequency domains do not intersect"));
        }
        let wire_length_per_member_length = positive(circumference.hypot(pitch) / pitch, "helix length ratio")?;
        let area = positive(std::f64::consts::PI * a * a, "winding wire area")?;
        let winding_line_density_kg_m = positive(coefficients[0] * area * wire_length_per_member_length,
            "winding mass per member length")?;
        // Equivalent to G d^4 p/(8 D^3), evaluated without tiny fourth powers.
        let winding_axial_stiffness_n = positive(coefficients[1] * a * (a / radius).powi(3) * pitch / 4.0,
            "close-coiled winding axial stiffness")?;
        let line_density_kg_m = positive(core.elastic().line_density_kg_m() + winding_line_density_kg_m,
            "core plus winding mass per length")?;
        let axial_stiffness_n = positive(core.elastic().axial_stiffness_n() + winding_axial_stiffness_n,
            "core plus winding axial stiffness")?;
        Ok(Self { core, winding, material, shear_property: shear_property.into(), helix_radius_m: radius,
            wire_length_per_member_length, winding_line_density_kg_m, winding_axial_stiffness_n,
            line_density_kg_m, axial_stiffness_n, omega_band_rad_s: band })
    }

    /// Admit a spatial basis only when `k_max * pitch <= 0.2`. This is an
    /// explicit scale-separation assumption, not a numerical convergence test.
    pub fn admit_wavenumber(&self, maximum_wavenumber_rad_m: f64) -> Result<(), MaterialStatePointError> {
        positive(maximum_wavenumber_rad_m, "positive wound-member wavenumber")?;
        if maximum_wavenumber_rad_m * self.winding.pitch_m > 0.2 {
            return Err(invalid("wound-member wavelength does not separate from the winding pitch"));
        }
        Ok(())
    }

    /// Original core, including every material and selected loss-law receipt.
    #[must_use] pub const fn core(&self) -> &WoundCoreSection { &self.core }
    /// Actual winding-wire material, also the exposed contact material.
    #[must_use] pub const fn winding_material(&self) -> &ResolvedMaterialStatePoint { &self.material }
    /// Wire and turn geometry of the frozen specimen.
    #[must_use] pub const fn winding(&self) -> RoundWireWinding { self.winding }
    /// Name of the selected shear-modulus property.
    #[must_use] pub fn shear_property(&self) -> &str { &self.shear_property }
    /// Centre-line radius of the winding helix [m].
    #[must_use] pub const fn helix_radius_m(&self) -> f64 { self.helix_radius_m }
    /// Physical finished radius, including both sides of the winding wire [m].
    #[must_use] pub fn outer_radius_m(&self) -> f64 { self.helix_radius_m + self.winding.wire_radius_m }
    /// Actual wire length per member length, a dimensionless geometric ratio.
    #[must_use] pub const fn wire_length_per_member_length(&self) -> f64 { self.wire_length_per_member_length }
    /// Winding-only line density [kg/m].
    #[must_use] pub const fn winding_line_density_kg_m(&self) -> f64 { self.winding_line_density_kg_m }
    /// Winding spring's axial rigidity [N].
    #[must_use] pub const fn winding_axial_stiffness_n(&self) -> f64 { self.winding_axial_stiffness_n }
    /// Total line density [kg/m], without an invented equivalent bulk material.
    #[must_use] pub const fn line_density_kg_m(&self) -> f64 { self.line_density_kg_m }
    /// Core plus winding axial rigidity [N].
    #[must_use] pub const fn axial_stiffness_n(&self) -> f64 { self.axial_stiffness_n }
    /// Intersection of all participating coefficients' support [rad/s].
    #[must_use] pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::section::{ConcentricElasticSection, ConcentricMaterialLayer};
    use crate::state_point::{MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point};
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_qty::DynViscosity;

    fn material(rho: f64, shear: f64, temperature: f64, band: f64) -> ResolvedMaterialStatePoint {
        let mut claims = ClaimSet::new();
        let mut required = Vec::new();
        for (name, value, dims) in [("density", rho, Density::DIMS), ("young_modulus", 200e9, Pressure::DIMS),
            ("shear_modulus", shear, Pressure::DIMS), ("eta", 10.0, DynViscosity::DIMS)] {
            claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, band),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                provenance: Provenance { source: "synthetic winding mechanics".into(), license: "test-only".into(), artifact: None },
            }).unwrap();
            required.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
            phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
        resolve_material_state_point(&card, &QueryPoint::new().with("T", temperature).unwrap()
            .with("omega", 1000.0).unwrap(), &required, MaterialPropertySelection::SingleClaimOnly).unwrap()
    }
    fn core() -> ConcentricKelvinVoigtSection {
        ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
            outer_radius_m: 0.0002, material: material(7800.0, 80e9, 293.15, 200_000.0),
        }]).unwrap().with_kelvin_voigt_bending("eta").unwrap()
    }
    fn section(rho: f64, shear: f64) -> CloseWoundSection {
        CloseWoundSection::try_new(core(), RoundWireWinding { wire_radius_m: 0.00005, pitch_m: 0.000101 },
            material(rho, shear, 293.15, 100_000.0), "shear_modulus").unwrap()
    }
    #[test]
    fn g1_helix_mass_matches_resolved_wire_polyline_and_is_not_a_solid_sleeve() {
        let winding = section(8900.0, 76e9);
        let r = winding.helix_radius_m();
        let p = winding.winding().pitch_m;
        let steps = 16_384;
        let mut length = 0.0;
        let mut previous = [r, 0.0, 0.0];
        for i in 1..=steps {
            let t = i as f64 / steps as f64;
            let angle = std::f64::consts::TAU * t;
            let point = [r * angle.cos(), r * angle.sin(), p * t];
            length += ((point[0] - previous[0]).powi(2) + (point[1] - previous[1]).powi(2)
                + (point[2] - previous[2]).powi(2)).sqrt();
            previous = point;
        }
        let resolved_mass = 8900.0 * std::f64::consts::PI * winding.winding().wire_radius_m.powi(2) * length / p;
        assert!((resolved_mass / winding.winding_line_density_kg_m() - 1.0).abs() < 1e-8);
        let solid_sleeve_mass = 8900.0 * std::f64::consts::PI * (0.0003_f64.powi(2) - 0.0002_f64.powi(2));
        assert!(winding.winding_line_density_kg_m() < solid_sleeve_mass);
        assert_eq!(winding.omega_band_rad_s(), (0.0, 100_000.0));
    }
    #[test]
    fn g1_spring_energy_and_material_swap_keep_density_and_modulus_independent() {
        let a = section(8900.0, 76e9);
        let b = section(17800.0, 38e9);
        assert_eq!(b.winding_line_density_kg_m(), 2.0 * a.winding_line_density_kg_m());
        assert_eq!(b.winding_axial_stiffness_n(), 0.5 * a.winding_axial_stiffness_n());
        assert_eq!(a.core(), b.core());
        assert_ne!(a.winding_material().identity(), b.winding_material().identity());
        // Integrate torsional strain energy of N close-coiled turns under F:
        // U = (F R)^2 (2 pi R N)/(2 G J), so F^2 L/(2 U) = EA.
        let force = 0.2_f64;
        let turns = 40.0;
        let length = turns * a.winding().pitch_m;
        let polar_area = std::f64::consts::PI * a.winding().wire_radius_m.powi(4) / 2.0;
        let energy = (force * a.helix_radius_m()).powi(2)
            * (std::f64::consts::TAU * a.helix_radius_m() * turns) / (2.0 * 76e9 * polar_area);
        let energy_rigidity = force * force * length / (2.0 * energy);
        assert!((energy_rigidity / a.winding_axial_stiffness_n() - 1.0).abs() < 1e-14);
        assert!(a.winding_axial_stiffness_n() < 0.01 * a.core().elastic().axial_stiffness_n());
    }
    #[test]
    fn g0_overlapping_turns_state_mismatch_and_short_waves_refuse() {
        let source = material(8900.0, 76e9, 293.15, 100_000.0);
        for (radius, pitch) in [(0.00005, 0.0001), (0.0001, 0.00021), (0.00005, 0.001)] {
            assert!(CloseWoundSection::try_new(core(), RoundWireWinding { wire_radius_m: radius, pitch_m: pitch },
                source.clone(), "shear_modulus").is_err());
        }
        assert!(CloseWoundSection::try_new(core(), RoundWireWinding { wire_radius_m: 0.00005, pitch_m: 0.000101 },
            material(8900.0, 76e9, 300.0, 100_000.0), "shear_modulus").is_err());
        let winding = section(8900.0, 76e9);
        assert!(winding.admit_wavenumber(10.0).is_ok());
        assert!(winding.admit_wavenumber(2000.0).is_err());
        assert!(winding.admit_wavenumber(f64::NAN).is_err());
    }
}
