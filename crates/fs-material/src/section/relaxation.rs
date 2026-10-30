//! Material-bound relaxation of perfectly bonded concentric sections.
//!
//! Layers share curvature or twist gradient, so their moments add. A bending
//! arm contributes `I_i delta_E_ij`; a circular torsional arm contributes
//! `J_i delta_G_ij`, with `J_i = 2 I_i`. Every source relaxation time survives.
//! Averaging time constants would change both storage and dissipation.
//! This is the parallel-arm model used by `visco::GeneralizedMaxwell`
//! (Roylance, MIT 3.11, "Linear Viscoelasticity", Wiechert model).
//!
//! The section is frozen at the resolved material states. Axial creep, thermal
//! history, interlayer slip and helical winding are not modeled by this binding.

use super::{ConcentricElasticSection, MaterialStatePointError, invalid};
use crate::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY};
use crate::visco::GeneralizedMaxwell;
use fs_qty::{Density, Pressure, QuantitySpec, Time};

/// Cold-path bound on the total declared arms across all layers.
pub const MAX_SECTION_RELAXATION_BRANCHES: usize = 128;

/// Explicit source roles for one layer's selected strain channel.
#[derive(Clone, Debug, PartialEq)]
pub struct LayerRelaxationProperties {
    /// Long-time modulus [Pa], matching the channel's elastic Young/shear modulus.
    pub equilibrium_modulus: String,
    /// Ordered (modulus increment [Pa], relaxation time [s]) property names.
    /// An empty list declares an elastic channel; the equilibrium claim must
    /// supply explicit frequency support. Missing source data is never zero loss.
    pub branches: Vec<(String, String)>,
}

/// Bending selections use the same source-role data as other section channels.
pub use self::LayerRelaxationProperties as LayerBendingRelaxationProperties;

/// One original layer arm integrated over the actual bending/polar area moment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerRelaxationBranch {
    /// Index into the retained section's spatially ordered layers.
    pub layer_index: usize,
    /// Nonnegative `delta_E * I` or `delta_G * J` [N m2].
    pub relaxing_stiffness_n_m2: f64,
    /// Original source time [s], not averaged with another layer or arm.
    pub relaxation_time_s: f64,
}

/// The bending-channel spelling of the shared integrated-arm representation.
pub use self::LayerRelaxationBranch as LayerBendingRelaxationBranch;

/// Distinct strain channels select their actual section integral.
#[derive(Clone, Copy, Debug)]
pub(super) enum SectionRelaxationAxis { Bending, Torsion }

/// Shared coefficient admission; physical section wrappers keep channel identity.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct SectionPronySpectrum {
    pub(super) properties: Vec<LayerRelaxationProperties>,
    pub(super) branches: Vec<LayerRelaxationBranch>,
    pub(super) spectra: Vec<GeneralizedMaxwell>,
    pub(super) omega_band_rad_s: (f64, f64),
}

/// A passive bending spectrum on the same source-bound elastic section.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricPronySection {
    elastic: ConcentricElasticSection,
    spectrum: SectionPronySpectrum,
}

impl ConcentricElasticSection {
    /// Resolve every layer's independent generalized-Maxwell bending spectrum.
    ///
    /// Density, Young's modulus, equilibrium modulus and all branch coefficients
    /// must be constants at the held, explicitly typed environment. Each
    /// relaxation time declares a frequency band, intersected with every other
    /// coefficient's band. For an empty layer spectrum its equilibrium modulus
    /// declares the band instead. All original receipts and selectors survive.
    ///
    /// # Errors
    /// Refuses missing/ambiguous/wrong-quantity properties, a non-equilibrium
    /// elastic modulus, repeated source pairs, unsupported source dependence,
    /// nonpositive times, negative increments, disjoint support and arithmetic
    /// overflow/underflow. Every layer needs a binding and at most 128 total arms.
    pub fn with_prony_bending(
        self,
        properties: &[LayerRelaxationProperties],
    ) -> Result<ConcentricPronySection, MaterialStatePointError> {
        let spectrum = SectionPronySpectrum::resolve(&self, YOUNG_MODULUS_PROPERTY,
            self.bending_stiffness_n_m2(), SectionRelaxationAxis::Bending, properties)?;
        Ok(ConcentricPronySection { elastic: self, spectrum })
    }
}

impl SectionPronySpectrum {
    pub(super) fn resolve(
        elastic: &ConcentricElasticSection,
        elastic_modulus_property: &str,
        equilibrium_stiffness_n_m2: f64,
        axis: SectionRelaxationAxis,
        properties: &[LayerRelaxationProperties],
    ) -> Result<Self, MaterialStatePointError> {
        if properties.len() != elastic.layers().len() {
            return Err(invalid("section relaxation needs one source binding per layer"));
        }
        let mut count = 0_usize;
        for binding in properties {
            if binding.branches.len() > MAX_SECTION_RELAXATION_BRANCHES - count {
                return Err(invalid("section relaxation exceeds 128 total arms"));
            }
            count += binding.branches.len();
        }
        let mut branches = Vec::with_capacity(count);
        let mut spectra = Vec::with_capacity(properties.len());
        let mut band = (0.0_f64, f64::MAX);
        let mut instantaneous = equilibrium_stiffness_n_m2;
        for (layer_index, ((layer, geometry), binding)) in elastic.layers().iter()
            .zip(elastic.layer_integrals()).zip(properties).enumerate()
        {
            let mut coefficient = |key: &str, quantity: QuantitySpec, needs_band: bool| {
                let property = layer.material.property(key)
                    .ok_or_else(|| invalid("section relaxation is missing an unambiguous source coefficient"))?;
                if property.requirement().quantity() != quantity
                    || !property.is_constant_at_fixed_environment()
                {
                    return Err(invalid("section relaxation needs held-environment coefficients with the declared SI quantity"));
                }
                match property.angular_frequency_band_rad_s()? {
                    Some((lo, hi)) => band = (band.0.max(lo), band.1.min(hi)),
                    None if needs_band => return Err(invalid("each layer relaxation law needs explicit frequency support")),
                    None => (),
                }
                Ok(property.value_si())
            };
            coefficient(DENSITY_PROPERTY, QuantitySpec::dimensional(Density::DIMS), false)?;
            let elastic_modulus = coefficient(elastic_modulus_property,
                QuantitySpec::dimensional(Pressure::DIMS), false)?;
            let equilibrium = coefficient(&binding.equilibrium_modulus,
                QuantitySpec::dimensional(Pressure::DIMS), binding.branches.is_empty())?;
            if !equilibrium.is_finite() || equilibrium <= 0.0
                || equilibrium.to_bits() != elastic_modulus.to_bits()
            {
                return Err(invalid("each elastic layer must use its explicit equilibrium modulus"));
            }
            let moment = match axis {
                SectionRelaxationAxis::Bending => geometry.second_moment_m4(),
                SectionRelaxationAxis::Torsion => 2.0 * geometry.second_moment_m4(),
            };
            let mut terms = Vec::with_capacity(binding.branches.len());
            for (index, pair) in binding.branches.iter().enumerate() {
                if binding.branches[..index].contains(pair) {
                    return Err(invalid("layer relaxation repeats a source-property pair"));
                }
                let delta = coefficient(&pair.0, QuantitySpec::dimensional(Pressure::DIMS), false)?;
                let tau = coefficient(&pair.1, QuantitySpec::dimensional(Time::DIMS), true)?;
                if !delta.is_finite() || delta < 0.0 || !tau.is_finite() || tau <= 0.0 {
                    return Err(invalid("layer relaxation needs nonnegative modulus increments and positive times"));
                }
                let stiffness = delta * moment;
                instantaneous += stiffness;
                if !stiffness.is_finite() || (delta > 0.0 && stiffness == 0.0)
                    || !instantaneous.is_finite()
                {
                    return Err(invalid("section relaxation stiffness is unrepresentable"));
                }
                branches.push(LayerRelaxationBranch {
                    layer_index, relaxing_stiffness_n_m2: stiffness, relaxation_time_s: tau,
                });
                terms.push((delta, tau));
            }
            spectra.push(GeneralizedMaxwell::new(equilibrium, terms)
                .map_err(|_| invalid("layer relaxation spectrum is outside its constitutive domain"))?);
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid("layer relaxation frequency bands do not intersect"));
        }
        Ok(Self { properties: properties.to_vec(), branches, spectra, omega_band_rad_s: band })
    }
}

impl ConcentricPronySection {
    /// Original elastic geometry, mass, equilibrium stiffness and material receipts.
    #[must_use]
    pub const fn elastic(&self) -> &ConcentricElasticSection { &self.elastic }
    /// Source selectors in the original spatial layer order.
    #[must_use]
    pub fn properties(&self) -> &[LayerRelaxationProperties] { &self.spectrum.properties }
    /// Every source arm, ordered by layer and then declaration within that layer.
    #[must_use]
    pub fn branches(&self) -> &[LayerRelaxationBranch] { &self.spectrum.branches }
    /// Original constitutive spectra [Pa, s], one per spatial layer.
    #[must_use]
    pub fn spectra(&self) -> &[GeneralizedMaxwell] { &self.spectrum.spectra }
    /// Intersection of all contributing source frequency supports [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.spectrum.omega_band_rad_s }
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

    fn material(eq: f64, arms: &[(&str, &str, f64, f64)], band: Option<(f64, f64)>)
        -> ResolvedMaterialStatePoint {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        let mut values = vec![(DENSITY_PROPERTY, Density::DIMS, 7800.0),
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 2e9), ("equilibrium", Pressure::DIMS, eq)];
        for &(modulus, time, delta, tau) in arms {
            values.extend([(modulus, Pressure::DIMS, delta), (time, Time::DIMS, tau)]);
        }
        for (key, dims, value) in values {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
                validity: band.map_or_else(ValidityDomain::unconstrained,
                    |(lo, hi)| ValidityDomain::unconstrained().with("omega", lo, hi)),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                provenance: Provenance { source: "synthetic section spectrum".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(key, dims,
                ScalarAdmissibility::Finite).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
            phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
        let query = QueryPoint::new().with("omega", band.map_or(100.0, |(lo, hi)| 0.5*(lo+hi))).unwrap();
        resolve_material_state_point(&card, &query, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap()
    }
    fn selection(pairs: &[(&str, &str)]) -> LayerBendingRelaxationProperties {
        LayerBendingRelaxationProperties { equilibrium_modulus: "equilibrium".into(),
            branches: pairs.iter().map(|(a,b)| ((*a).into(),(*b).into())).collect() }
    }
    fn section(a: ResolvedMaterialStatePoint, b: ResolvedMaterialStatePoint) -> ConcentricElasticSection {
        ConcentricElasticSection::try_new(0.0, vec![
            ConcentricMaterialLayer { outer_radius_m: 0.001, material: a },
            ConcentricMaterialLayer { outer_radius_m: 0.002, material: b },
        ]).unwrap()
    }
    fn response(law: &ConcentricPronySection, omega: f64) -> (f64, f64) {
        GeneralizedMaxwell::new(law.elastic().bending_stiffness_n_m2(), law.branches().iter()
            .map(|arm| (arm.relaxing_stiffness_n_m2, arm.relaxation_time_s)).collect())
            .unwrap().modulus(omega)
    }
    fn near(a: f64, b: f64) {
        assert!((a-b).abs() < 2e-13*a.abs().max(b.abs()).max(1e-30), "{a} != {b}");
    }

    #[test]
    fn g1_distinct_layer_spectra_preserve_complex_storage_and_loss() {
        let a = material(2e9, &[("core_delta", "core_time", 4e9, 0.001)], Some((0.0, 20000.0)));
        let b = material(2e9, &[("fast", "fast_time", 0.5e9, 0.003),
            ("slow", "slow_time", 3e9, 0.02)], Some((0.0, 10000.0)));
        let elastic = section(a,b);
        let law = elastic.clone().with_prony_bending(&[
            selection(&[("core_delta", "core_time")]),
            selection(&[("fast", "fast_time"), ("slow", "slow_time")]),
        ]).unwrap();
        assert_eq!(law.elastic(), &elastic);
        assert_eq!(law.omega_band_rad_s(), (0.0,10000.0));
        assert_eq!(law.branches().iter().map(|a|a.layer_index).collect::<Vec<_>>(), vec![0,1,1]);
        assert_eq!(law.branches().iter().map(|a|a.relaxation_time_s).collect::<Vec<_>>(), vec![0.001,0.003,0.02]);
        let i0 = core::f64::consts::PI*0.001_f64.powi(4)/4.0;
        let i1 = core::f64::consts::PI*(0.002_f64.powi(4)-0.001_f64.powi(4))/4.0;
        for omega in [0.0, 5.0, 50.0, 1000.0, 10000.0] {
            let mut storage = 2e9*(i0+i1); let mut loss = 0.0;
            for (stiffness,tau) in [(i0*4e9,0.001),(i1*0.5e9,0.003),(i1*3e9,0.02)] {
                let x = omega*tau;
                storage += stiffness*x*x/(1.0+x*x);
                loss += stiffness*x/(1.0+x*x);
            }
            let actual = response(&law,omega);
            near(actual.0,storage); near(actual.1,loss);
        }
    }

    #[test]
    fn g3_homogeneous_subdivision_preserves_the_entire_response() {
        let m = material(2e9, &[("delta", "time", 4e9, 0.001)], Some((0.0,20000.0)));
        let bind = selection(&[("delta", "time")]);
        let single = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
            outer_radius_m: 0.002, material: m.clone(),
        }]).unwrap().with_prony_bending(&[bind.clone()]).unwrap();
        let split = section(m.clone(),m).with_prony_bending(&[bind.clone(),bind]).unwrap();
        for omega in [0.0,10.0,1000.0,20000.0] {
            let a = response(&single,omega); let b = response(&split,omega);
            near(a.0,b.0); near(a.1,b.1);
        }
    }

    #[test]
    fn g0_empty_is_explicit_elastic_but_absent_or_mismatched_sources_refuse() {
        let elastic = material(2e9,&[],Some((0.0,20000.0)));
        let lossy = material(2e9,&[("delta","time",4e9,0.001)],Some((0.0,10000.0)));
        let plain = selection(&[]); let loss = selection(&[("delta","time")]);
        let base = section(elastic.clone(),lossy.clone());
        let law = base.clone().with_prony_bending(&[plain.clone(),loss.clone()]).unwrap();
        assert_eq!(law.branches().len(),1);
        assert_eq!(law.branches()[0].layer_index,1);
        assert!(base.clone().with_prony_bending(&[loss.clone(),loss.clone()]).is_err());
        assert!(base.clone().with_prony_bending(&[plain.clone()]).is_err());
        assert!(base.with_prony_bending(&[plain.clone(),selection(&[("delta","time"),("delta","time")])]).is_err());
        assert!(section(material(3e9,&[],Some((0.0,20000.0))),lossy.clone())
            .with_prony_bending(&[plain.clone(),loss.clone()]).is_err());
        assert!(section(material(2e9,&[],None),lossy)
            .with_prony_bending(&[plain.clone(),loss]).is_err());
        assert!(section(elastic.clone(),material(2e9,&[("delta","time",1.0,0.001)],Some((0.0,20000.0))))
            .with_prony_bending(&[plain.clone(),selection(&[("time","delta")])]).is_err());
        for (delta,tau) in [(-1.0,0.001),(1.0,0.0),(1.0,-1.0)] {
            assert!(section(elastic.clone(), material(2e9,&[("delta","time",delta,tau)],Some((0.0,20000.0))))
                .with_prony_bending(&[plain.clone(),selection(&[("delta","time")])]).is_err());
        }
        let no_overlap = section(material(2e9,&[],Some((0.0,10.0))),
            material(2e9,&[],Some((20.0,100.0))));
        assert!(no_overlap.with_prony_bending(&[plain.clone(),plain]).is_err());
    }

    #[test]
    fn g0_section_arm_budget_is_checked_before_publication() {
        let source = material(2e9,&[("delta","time",1e9,0.001)],Some((0.0,20000.0)));
        let section = ConcentricElasticSection::try_new(0.0,(0..MAX_SECTION_RELAXATION_BRANCHES)
            .map(|index|ConcentricMaterialLayer {
                outer_radius_m: 0.001*(index+1) as f64, material: source.clone(),
            }).collect()).unwrap();
        let mut bindings = vec![selection(&[("delta","time")]);MAX_SECTION_RELAXATION_BRANCHES];
        assert_eq!(section.clone().with_prony_bending(&bindings).unwrap().branches().len(),
            MAX_SECTION_RELAXATION_BRANCHES);
        bindings[0].branches.push(("other_delta".into(),"other_time".into()));
        assert!(matches!(section.with_prony_bending(&bindings),
            Err(MaterialStatePointError::InvalidDerived { quantity: "section relaxation exceeds 128 total arms" })));
    }
}
