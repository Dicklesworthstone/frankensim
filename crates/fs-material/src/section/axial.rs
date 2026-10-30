//! Source-bound linear axial Kelvin--Voigt loss for bonded circular members.
//!
//! Common axial strain gives N = EA*epsilon + C_a*epsilon_dot, with
//! C_a = sum(eta_i*A_i) [N s]. The elastic section already owns A_i, E_i,
//! density and their source identities. No bending/shear viscosity, fitted Q,
//! material name or fluid property is substituted for an explicit axial claim.
//!
//! This is frozen-state SMALL-STRAIN loss. Applying it to longitudinal modes
//! dissipates integral(C_a*u_xt^2 dx); it does not model the additional
//! finite-slope rate y_x*y_xt, creep, heating, slip or evolving material memory.

use super::{ConcentricElasticSection, MaterialStatePointError, invalid};
use crate::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY};
use fs_qty::{DynViscosity, QuantitySpec};

/// Uniaxial solid viscosity, also used by the isotropic membrane-wall binding.
/// It is distinct from the bending and torsional-shear viscosity properties.
pub const YOUNG_VISCOSITY_PROPERTY: &str = "kelvin_voigt_young_viscosity";

/// Linear axial loss with the same ordered solid states and section geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricAxialViscosity {
    elastic: ConcentricElasticSection,
    viscosity_property: String,
    layer_viscous_rigidity_n_s: Vec<f64>,
    viscous_rigidity_n_s: f64,
    omega_band_rad_s: (f64, f64),
}

impl ConcentricElasticSection {
    /// Bind explicit uniaxial viscosity [Pa s] from every original layer.
    ///
    /// Density, E and viscosity must be constant at the held typed environment.
    /// Every viscosity supplies a frequency band; all coefficient bands are
    /// intersected. Explicit zero is lossless, but absent loss data is an error.
    /// The section is perfectly bonded: no averaging of material identities.
    ///
    /// # Errors
    /// Missing/wrong-quantity/negative data, frequency-dependent coefficients,
    /// disjoint source support or unrepresentable integrals refuse atomically.
    pub fn with_axial_viscosity(
        self,
        property_name: &str,
    ) -> Result<ConcentricAxialViscosity, MaterialStatePointError> {
        let mut band = (0.0_f64, f64::MAX);
        let mut total = 0.0;
        let mut layers = Vec::with_capacity(self.layers().len());
        for (layer, geometry) in self.layers().iter().zip(self.layer_integrals()) {
            let eta = layer.material.property(property_name)
                .ok_or_else(|| invalid("every axial layer needs an explicit uniaxial viscosity"))?;
            if eta.requirement().quantity() != QuantitySpec::dimensional(DynViscosity::DIMS)
                || !eta.value_si().is_finite() || eta.value_si() < 0.0 {
                return Err(invalid("axial viscosity must be nonnegative dimension-only Pa s"));
            }
            let interval = eta.angular_frequency_band_rad_s()?
                .ok_or_else(|| invalid("axial viscosity needs an explicit frequency validity band"))?;
            band = (band.0.max(interval.0), band.1.min(interval.1));
            for key in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, property_name] {
                let coefficient = layer.material.property(key)
                    .ok_or_else(|| invalid("axial loss lost an original section coefficient"))?;
                if !coefficient.is_constant_at_fixed_environment() {
                    return Err(invalid("axial loss requires fixed-state constant density, modulus and viscosity"));
                }
                if let Some((lo, hi)) = coefficient.angular_frequency_band_rad_s()? {
                    band = (band.0.max(lo), band.1.min(hi));
                }
            }
            let value = eta.value_si() * geometry.area_m2();
            if !value.is_finite() || (eta.value_si() > 0.0 && value == 0.0) {
                return Err(invalid("layer axial viscous rigidity is unrepresentable"));
            }
            total += value;
            if !total.is_finite() { return Err(invalid("axial viscous rigidity sum overflow")); }
            layers.push(value);
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid("axial coefficient frequency bands do not intersect"));
        }
        Ok(ConcentricAxialViscosity {
            elastic: self, viscosity_property: property_name.into(),
            layer_viscous_rigidity_n_s: layers, viscous_rigidity_n_s: total,
            omega_band_rad_s: band,
        })
    }
}

impl ConcentricAxialViscosity {
    /// Original elasticity, mass, geometry and every material receipt.
    #[must_use]
    pub const fn elastic(&self) -> &ConcentricElasticSection { &self.elastic }
    /// Explicit axial selector, not a bending-loss alias.
    #[must_use]
    pub fn viscosity_property(&self) -> &str { &self.viscosity_property }
    /// Each eta_i*A_i [N s], in the original spatial layer order.
    #[must_use]
    pub fn layer_viscous_rigidity_n_s(&self) -> &[f64] { &self.layer_viscous_rigidity_n_s }
    /// C_a [N s] multiplying axial strain rate.
    #[must_use]
    pub const fn viscous_rigidity_n_s(&self) -> f64 { self.viscous_rigidity_n_s }
    /// Intersected fixed-state source applicability [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }

    /// Damping of a fixed-end longitudinal eigenmode of this SAME section.
    /// For k=n*pi/L, omega^2=(EA/mu)*k^2 and 2*zeta*omega=(C_a/mu)*k^2.
    /// The consumer supplies the undamped elastic omega, not a transverse or
    /// gas-loaded frequency. No solver clock or nonlinear-harmonic admission is
    /// implied; those remain the receiving dynamics owner's responsibility.
    pub fn longitudinal_damping_ratio(&self, omega: f64) -> Result<f64, MaterialStatePointError> {
        let (lo, hi) = self.omega_band_rad_s;
        if !omega.is_finite() || omega <= 0.0 || omega < lo || omega > hi {
            return Err(invalid("longitudinal mode exceeds axial material-loss frequency support"));
        }
        let zeta = 0.5 * (self.viscous_rigidity_n_s / self.elastic.axial_stiffness_n()) * omega;
        if !zeta.is_finite() || (self.viscous_rigidity_n_s > 0.0 && zeta == 0.0) {
            return Err(invalid("longitudinal material damping is unrepresentable"));
        }
        Ok(zeta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::section::ConcentricMaterialLayer;
    use crate::state_point::{ResolvedMaterialStatePoint, MaterialPropertySelection,
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_qty::{Density, Dims, Pressure};

    fn material(eta: f64, band: Option<(f64, f64)>, sampled: bool) -> ResolvedMaterialStatePoint {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, value, dims) in [
            (DENSITY_PROPERTY, 7800.0, Density::DIMS),
            (YOUNG_MODULUS_PROPERTY, 200e9, Pressure::DIMS),
            (YOUNG_VISCOSITY_PROPERTY, eta, DynViscosity::DIMS),
        ] {
            let mut validity = ValidityDomain::unconstrained();
            if let Some((lo, hi)) = band { validity = validity.with("omega", lo, hi); }
            let curve = sampled && name == YOUNG_VISCOSITY_PROPERTY;
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if curve { PropertyValue::Curve {
                    abscissa: "omega".into(), abscissa_dims: Dims([0,0,-1,0,0,0]),
                    knots: vec![(1.0,eta),(20000.0,2.0*eta)], dims,
                } } else { PropertyValue::Scalar { value, dims } },
                validity, uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve { InterpolationPolicy::LinearInside }
                    else { InterpolationPolicy::ConstantWithinValidity },
                observations: vec![], provenance: Provenance {
                    source: "synthetic axial viscosity".into(), license: "test-only".into(), artifact: None,
                },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: "synthetic".into(), phase: "solid".into(), process: "test".into(), revision: 0,
        }, claims, vec![]).unwrap();
        let query = QueryPoint::new().with("omega", band.map_or(100.0, |(lo, hi)| 0.5*(lo+hi))).unwrap();
        resolve_material_state_point(&card, &query, &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
    }
    fn section(a: ResolvedMaterialStatePoint, b: ResolvedMaterialStatePoint) -> ConcentricElasticSection {
        ConcentricElasticSection::try_new(0.0, vec![
            ConcentricMaterialLayer { outer_radius_m: 0.001, material: a },
            ConcentricMaterialLayer { outer_radius_m: 0.002, material: b },
        ]).unwrap()
    }
    fn close(a: f64, b: f64) { assert!((a-b).abs() < 3e-13*a.abs().max(b.abs()).max(1e-30), "{a} != {b}"); }

    #[test]
    fn g1_layer_integrals_and_modal_decay_use_actual_axial_sources() {
        let a = material(1e7, Some((1.0,20000.0)), false);
        let b = material(4e7, Some((10.0,10000.0)), false);
        let elastic = section(a,b);
        let law = elastic.clone().with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY).unwrap();
        assert_eq!(law.elastic(), &elastic);
        assert_eq!(law.omega_band_rad_s(), (10.0,10000.0));
        let expected = [1e7*elastic.layer_integrals()[0].area_m2(), 4e7*elastic.layer_integrals()[1].area_m2()];
        assert_eq!(law.layer_viscous_rigidity_n_s(), &expected);
        close(law.viscous_rigidity_n_s(), expected.iter().sum());
        for omega in [10.0,1000.0,10000.0] {
            let k2 = omega*omega*elastic.line_density_kg_m()/elastic.axial_stiffness_n();
            close(2.0*law.longitudinal_damping_ratio(omega).unwrap()*omega,
                law.viscous_rigidity_n_s()/elastic.line_density_kg_m()*k2);
        }
        for omega in [0.0,9.0,10001.0,f64::NAN] {
            assert!(law.longitudinal_damping_ratio(omega).is_err());
        }
    }

    #[test]
    fn g3_homogeneous_subdivision_does_not_change_axial_decay() {
        let m = material(2e6, Some((1.0,20000.0)), false);
        let whole = ConcentricElasticSection::try_new(0.0, vec![
            ConcentricMaterialLayer { outer_radius_m: 0.002, material: m.clone() },
        ]).unwrap().with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY).unwrap();
        let split = section(m.clone(),m).with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY).unwrap();
        close(whole.viscous_rigidity_n_s(), split.viscous_rigidity_n_s());
        close(whole.longitudinal_damping_ratio(1000.0).unwrap(), split.longitudinal_damping_ratio(1000.0).unwrap());
    }

    #[test]
    fn g0_missing_negative_sampled_or_disjoint_loss_never_becomes_a_preset() {
        let m = material(0.0, Some((1.0,20000.0)), false);
        let elastic = section(m.clone(),m);
        assert_eq!(elastic.clone().with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY).unwrap()
            .longitudinal_damping_ratio(1000.0).unwrap(), 0.0);
        for key in ["missing", DENSITY_PROPERTY] {
            assert!(elastic.clone().with_axial_viscosity(key).is_err());
        }
        for m in [material(-1.0,Some((1.0,20000.0)),false),
            material(1e7,None,false), material(1e7,Some((1.0,20000.0)),true)] {
            assert!(section(m.clone(),m).with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY).is_err());
        }
        assert!(section(material(1e7,Some((1.0,100.0)),false),
            material(1e7,Some((101.0,200.0)),false))
            .with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY).is_err());
    }
}
