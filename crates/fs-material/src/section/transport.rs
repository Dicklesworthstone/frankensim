//! Axial transport on the existing source-bound concentric section.
//!
//! Mechanics and transport use the same layer ordering, geometry, density and
//! resolved material bundles. The local parallel-path reduction assumes a
//! common axial temperature gradient. It does not infer radial contact
//! conductance, a thermal expansion law, or a homogeneous replacement material.
//!
//! Coefficients describe the retained material states only. Re-resolve the
//! original layer bundles and rebuild the section when those states change;
//! these local rates are not a temperature integrator or a melting model.

use fs_qty::QuantitySpec;

use crate::state_point::{
    MaterialStatePointError, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};

use super::{ConcentricElasticSection, invalid, positive};

/// Thermal coefficients integrated over one existing mechanical layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalLayerIntegrals {
    /// `rho A cp` [J/(m K)]; uses the original mechanical line mass.
    pub heat_capacity_per_length_j_m_k: f64,
    /// `k A` [W m/K]; multiply by minus the axial temperature gradient for W.
    pub conductivity_area_w_m_k: f64,
}

/// Frozen-state thermal reduction retaining the exact mechanical section.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricThermalTransport<'a> {
    section: &'a ConcentricElasticSection,
    layers: Vec<ThermalLayerIntegrals>,
    heat_capacity_per_length_j_m_k: f64,
    conductivity_area_w_m_k: f64,
}

impl ConcentricElasticSection {
    /// Bind heat capacity and axial conduction from EVERY layer's existing
    /// resolved bundle. No second density query or averaged material is used.
    ///
    /// Explicit zero conductivity is insulating, not missing data. Each heat
    /// capacity must be positive. Layers may be at different retained states;
    /// their summed heat capacity describes a common incremental temperature
    /// change, not an assertion that the layers are in thermal equilibrium.
    ///
    /// # Errors
    /// Missing/ambiguous properties, incompatible quantities, negative
    /// conductivity, nonpositive heat capacity and unrepresentable integrals
    /// refuse without changing the mechanical section.
    pub fn thermal_transport(
        &self,
    ) -> Result<ConcentricThermalTransport<'_>, MaterialStatePointError> {
        let mut layers = Vec::with_capacity(self.layers().len());
        let mut capacity = 0.0;
        let mut conductivity_area = 0.0;
        for (layer, geometry) in self.layers().iter().zip(self.layer_integrals()) {
            let cp = layer_scalar(
                &layer.material,
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            )?;
            positive(cp, "section heat capacity must be positive")?;
            let k = layer_scalar(
                &layer.material,
                THERMAL_CONDUCTIVITY_PROPERTY,
                QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS),
            )?;
            let c = product(geometry.line_density_kg_m(), cp)?;
            let ka = product(geometry.area_m2(), k)?;
            capacity = positive(capacity + c, "section line heat capacity")?;
            conductivity_area = nonnegative(conductivity_area + ka)?;
            layers.push(ThermalLayerIntegrals {
                heat_capacity_per_length_j_m_k: c,
                conductivity_area_w_m_k: ka,
            });
        }
        Ok(ConcentricThermalTransport {
            section: self,
            layers,
            heat_capacity_per_length_j_m_k: capacity,
            conductivity_area_w_m_k: conductivity_area,
        })
    }
}

impl<'a> ConcentricThermalTransport<'a> {
    /// Original geometry, mechanics and all material-property receipts.
    #[must_use]
    pub const fn section(&self) -> &'a ConcentricElasticSection {
        self.section
    }

    /// Thermal integrals in the exact order of the mechanical layers.
    #[must_use]
    pub fn layers(&self) -> &[ThermalLayerIntegrals] {
        &self.layers
    }

    /// Sum of `rho_i A_i cp_i` [J/(m K)].
    #[must_use]
    pub const fn heat_capacity_per_length_j_m_k(&self) -> f64 {
        self.heat_capacity_per_length_j_m_k
    }

    /// Sum of `k_i A_i` [W m/K], not an interface conductance.
    #[must_use]
    pub const fn conductivity_area_w_m_k(&self) -> f64 {
        self.conductivity_area_w_m_k
    }

    /// Fourier heat flow [W], positive along the section's positive axis.
    /// This evaluates one common LOCAL axial temperature gradient [K/m].
    ///
    /// # Errors
    /// Nonfinite gradients and unrepresentable heat flows refuse.
    pub fn axial_heat_flow_w(
        &self,
        temperature_gradient_k_m: f64,
    ) -> Result<f64, MaterialStatePointError> {
        product(-self.conductivity_area_w_m_k, temperature_gradient_k_m)
    }
}

/// Frozen-state thermal and parallel DC paths on one mechanical section.
/// All layers conduct axially between common ideal equipotential terminals.
/// Contact resistance, radial currents, skin effect and thermoelectric effects
/// are not inferred. This is a constitutive reduction, not a circuit solver.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricElectrothermalTransport<'a> {
    thermal: ConcentricThermalTransport<'a>,
    conductors: Vec<crate::conductor::ResolvedConductor>,
    conductance_s: f64,
    length_m: f64,
}

/// Electrical dissipation and caloric coefficient of one original layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrothermalLayerResponse {
    /// Signed axial current [A].
    pub current_a: f64,
    /// Nonnegative heat source [W], ready for the caller's thermal balance.
    pub joule_power_w: f64,
    /// `rho A L cp` [J/K], at this layer's retained material state.
    pub heat_capacity_j_k: f64,
}

/// Signed terminal response and heat sources without advancing physical state.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectrothermalSectionResponse {
    /// Potential drop between the common terminals [V].
    pub voltage_drop_v: f64,
    /// Sum of the layer currents [A].
    pub current_a: f64,
    /// Sum of the layer heat sources [W].
    pub joule_power_w: f64,
    /// Axial Fourier heat flow [W], not a net heating power for the member.
    /// A thermal balance needs the DIFFERENCE of the end heat flows.
    pub axial_heat_flow_w: f64,
    /// Responses in the original section's spatial layer order.
    pub layers: Vec<ElectrothermalLayerResponse>,
}

impl<'a> ConcentricThermalTransport<'a> {
    /// Bind finite, positive DC resistivity from the SAME complete layer bundles.
    /// Reuses the conductor owner's admission, geometry and Joule-power laws.
    /// Every layer must be an admitted Ohmic conductor; missing resistivity
    /// does not silently become an insulating layer or a standard metal.
    ///
    /// # Errors
    /// Invalid length, incomplete/incompatible resistivity data and
    /// unrepresentable resistance/conductance refuse atomically.
    pub fn with_ohmic_transport(
        self,
        length_m: f64,
    ) -> Result<ConcentricElectrothermalTransport<'a>, crate::conductor::ConductorError> {
        let mut conductors = Vec::with_capacity(self.layers.len());
        let mut conductance_s = 0.0;
        for (layer, geometry) in self.section.layers().iter().zip(self.section.layer_integrals()) {
            let conductor = crate::conductor::ResolvedConductorState::try_from_material(
                layer.material.clone(),
            )?.into_uniform_conductor(length_m, geometry.area_m2())?;
            let conductance = positive(1.0 / conductor.resistance_ohm(), "layer conductance")?;
            conductance_s = positive(conductance_s + conductance, "section conductance")?;
            conductors.push(conductor);
        }
        Ok(ConcentricElectrothermalTransport { thermal: self, conductors, conductance_s, length_m })
    }
}

impl<'a> ConcentricElectrothermalTransport<'a> {
    /// Thermal response and original source-bound mechanical section.
    #[must_use]
    pub const fn thermal(&self) -> &ConcentricThermalTransport<'a> {
        &self.thermal
    }

    /// Equivalent parallel conductance [S] for the declared length.
    #[must_use]
    pub const fn conductance_s(&self) -> f64 {
        self.conductance_s
    }

    /// Evaluate current, layer-resolved Joule heat and local axial heat flow.
    /// The current drive is TOTAL member current, not current in each layer.
    /// No heat is spread by area or mass: each layer receives its own `I_i² R_i`.
    /// Re-resolve and rebind after temperature, geometry or material changes.
    ///
    /// # Errors
    /// Invalid drives/gradients or unrepresentable derived values refuse without
    /// mutating the input section or publishing partial heat sources.
    pub fn evaluate(
        &self,
        drive: crate::conductor::OhmicDrive,
        temperature_gradient_k_m: f64,
    ) -> Result<ElectrothermalSectionResponse, crate::conductor::ConductorError> {
        use crate::conductor::{ConductorError, OhmicDrive};
        let voltage = match drive {
            OhmicDrive::Voltage(value) => value,
            OhmicDrive::Current(value) => {
                let voltage = value / self.conductance_s;
                if value != 0.0 && voltage == 0.0 {
                    return Err(ConductorError::InvalidDerived { quantity: "section voltage drop" });
                }
                voltage
            }
        };
        let axial_heat_flow_w = self.thermal.axial_heat_flow_w(temperature_gradient_k_m)?;
        let mut layers = Vec::with_capacity(self.conductors.len());
        let mut current_a = 0.0;
        let mut joule_power_w = 0.0;
        for (conductor, thermal) in self.conductors.iter().zip(&self.thermal.layers) {
            let current = conductor.current_a_for_drive(OhmicDrive::Voltage(voltage))?;
            let power = conductor.joule_power_w(current)?;
            current_a += current;
            if !current_a.is_finite() {
                return Err(ConductorError::InvalidDerived { quantity: "total section current" });
            }
            joule_power_w = nonnegative(joule_power_w + power)?;
            layers.push(ElectrothermalLayerResponse {
                current_a: current,
                joule_power_w: power,
                heat_capacity_j_k: product(thermal.heat_capacity_per_length_j_m_k, self.length_m)?,
            });
        }
        Ok(ElectrothermalSectionResponse {
            voltage_drop_v: voltage, current_a, joule_power_w, axial_heat_flow_w, layers,
        })
    }
}

fn layer_scalar(
    material: &crate::state_point::ResolvedMaterialStatePoint,
    name: &str,
    quantity: QuantitySpec,
) -> Result<f64, MaterialStatePointError> {
    let property = material
        .property(name)
        .ok_or_else(|| invalid("every section layer needs the requested transport property"))?;
    if property.requirement().quantity() != quantity {
        return Err(invalid("section transport property has incompatible quantity semantics"));
    }
    nonnegative(property.value_si())
}

fn nonnegative(value: f64) -> Result<f64, MaterialStatePointError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(invalid("section transport coefficient must be finite and nonnegative"))
    }
}

fn product(a: f64, b: f64) -> Result<f64, MaterialStatePointError> {
    let result = a * b;
    if !a.is_finite() || !b.is_finite() || !result.is_finite()
        || (a != 0.0 && b != 0.0 && result == 0.0)
    {
        Err(invalid("section transport product is unrepresentable"))
    } else {
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::section::ConcentricMaterialLayer;
    use crate::state_point::{
        DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint,
        ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
        resolve_material_state_point,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_qty::{Density, Dims, Pressure};

    // Synthetic analytical inputs, not claims about any real metal or alloy.
    fn source(rho: f64, cp: f64, k: f64, extra: &[(&str, f64, Dims)]) -> ResolvedMaterialStatePoint {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        let mut entries = vec![
            (DENSITY_PROPERTY, rho, Density::DIMS),
            (YOUNG_MODULUS_PROPERTY, 2e11, Pressure::DIMS),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, cp, SPECIFIC_HEAT_CAPACITY_DIMS),
            (THERMAL_CONDUCTIVITY_PROPERTY, k, THERMAL_CONDUCTIVITY_DIMS),
        ];
        entries.extend_from_slice(extra);
        for (name, value, dims) in entries {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic section transport test".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(
                name, dims, ScalarAdmissibility::Finite,
            ).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: "synthetic".into(), phase: "solid".into(),
            process: "test".into(), revision: 0,
        }, claims, vec![]).unwrap();
        resolve_material_state_point(&card, &QueryPoint::new(), &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap()
    }

    fn section(a: ResolvedMaterialStatePoint, b: ResolvedMaterialStatePoint) -> ConcentricElasticSection {
        ConcentricElasticSection::try_new(0.0, vec![
            ConcentricMaterialLayer { outer_radius_m: 0.001, material: a },
            ConcentricMaterialLayer { outer_radius_m: 0.002, material: b },
        ]).unwrap()
    }

    fn near(a: f64, b: f64) {
        assert!((a - b).abs() <= 2e-13 * a.abs().max(b.abs()).max(f64::MIN_POSITIVE),
            "{a} != {b}");
    }

    #[test]
    fn same_layer_geometry_and_density_drive_mechanics_and_heat() {
        let elastic = section(source(8000.0, 400.0, 20.0, &[]),
            source(2700.0, 900.0, 200.0, &[]));
        let thermal = elastic.thermal_transport().unwrap();
        let layers = elastic.layer_integrals();
        assert!(core::ptr::eq(thermal.section(), &elastic));
        near(thermal.heat_capacity_per_length_j_m_k(),
            400.0 * layers[0].line_density_kg_m() + 900.0 * layers[1].line_density_kg_m());
        let ka = 20.0 * layers[0].area_m2() + 200.0 * layers[1].area_m2();
        near(thermal.axial_heat_flow_w(12.0).unwrap(), -12.0 * ka);
        near(thermal.axial_heat_flow_w(-12.0).unwrap(), 12.0 * ka);
    }

    #[test]
    fn homogeneous_subdivision_does_not_change_thermal_response() {
        let material = source(8000.0, 400.0, 20.0, &[]);
        for bore in [0.0, 0.0005] {
            let make = |radii: &[f64]| ConcentricElasticSection::try_new(bore,
                radii.iter().map(|&outer_radius_m| ConcentricMaterialLayer {
                    outer_radius_m, material: material.clone(),
                }).collect()).unwrap();
            let whole = make(&[0.002]);
            let split = make(&[0.001, 0.002]);
            let a = whole.thermal_transport().unwrap();
            let b = split.thermal_transport().unwrap();
            near(a.heat_capacity_per_length_j_m_k(), b.heat_capacity_per_length_j_m_k());
            near(a.conductivity_area_w_m_k(), b.conductivity_area_w_m_k());
        }
    }

    #[test]
    fn insulator_and_invalid_thermal_data_do_not_create_defaults() {
        let make = |cp, k| section(source(8000.0, cp, k, &[]), source(8000.0, cp, k, &[]));
        let insulator = make(400.0, 0.0);
        assert_eq!(insulator.thermal_transport().unwrap().axial_heat_flow_w(10.0).unwrap(), 0.0);
        for (cp, k) in [(0.0, 1.0), (-1.0, 1.0), (400.0, -1.0)] {
            assert!(make(cp, k).thermal_transport().is_err());
        }
        assert!(insulator.thermal_transport().unwrap().axial_heat_flow_w(f64::NAN).is_err());
    }

    fn electrical_source(rho: f64, cp: f64, k: f64, resistivity: f64) -> ResolvedMaterialStatePoint {
        source(rho, cp, k, &[(crate::conductor::ELECTRICAL_RESISTIVITY_PROPERTY,
            resistivity, crate::conductor::ELECTRICAL_RESISTIVITY_DIMS)])
    }

    #[test]
    fn shared_material_bundle_survives_the_conductor_bridge() {
        let material = electrical_source(8000.0, 400.0, 20.0, 1e-7);
        let conductor = crate::conductor::ResolvedConductorState::try_from_material(material.clone()).unwrap();
        assert_eq!(conductor.material(), &material);
        assert!(conductor.material().property(SPECIFIC_HEAT_CAPACITY_PROPERTY).is_some());
        assert!(conductor.material().property(YOUNG_MODULUS_PROPERTY).is_some());
        near(conductor.resistivity_ohm_m(), 1e-7);
    }

    #[test]
    fn parallel_current_and_joule_heat_follow_each_layers_resistivity() {
        use crate::conductor::OhmicDrive;
        let elastic = section(electrical_source(8000.0, 400.0, 20.0, 1e-7),
            electrical_source(2700.0, 900.0, 200.0, 4e-7));
        let coupled = elastic.thermal_transport().unwrap().with_ohmic_transport(2.0).unwrap();
        let g = elastic.layer_integrals();
        let expected_g = g[0].area_m2() / (2.0 * 1e-7) + g[1].area_m2() / (2.0 * 4e-7);
        near(coupled.conductance_s(), expected_g);
        let positive = coupled.evaluate(OhmicDrive::Voltage(3.0), 5.0).unwrap();
        let negative = coupled.evaluate(OhmicDrive::Voltage(-3.0), 5.0).unwrap();
        near(positive.current_a, 3.0 * expected_g);
        near(positive.joule_power_w, 3.0 * positive.current_a);
        near(positive.layers[0].joule_power_w / positive.layers[1].joule_power_w, 4.0 / 3.0);
        near(positive.joule_power_w, negative.joule_power_w);
        near(positive.current_a, -negative.current_a);
        near(positive.layers[0].heat_capacity_j_k, 2.0 * 400.0 * g[0].line_density_kg_m());
        let prescribed_current = coupled.evaluate(OhmicDrive::Current(positive.current_a), 5.0).unwrap();
        near(prescribed_current.voltage_drop_v, 3.0);
        near(prescribed_current.joule_power_w, positive.joule_power_w);
        assert!(core::ptr::eq(coupled.thermal().section(), &elastic));
    }

    #[test]
    fn data_only_substitution_changes_mechanics_heat_and_current_together() {
        use crate::conductor::OhmicDrive;
        let first = section(electrical_source(8000.0, 400.0, 20.0, 1e-7),
            electrical_source(2700.0, 900.0, 200.0, 4e-7));
        let swapped = section(electrical_source(16000.0, 700.0, 80.0, 2e-7),
            electrical_source(2700.0, 900.0, 200.0, 4e-7));
        let a = first.thermal_transport().unwrap().with_ohmic_transport(1.0).unwrap();
        let b = swapped.thermal_transport().unwrap().with_ohmic_transport(1.0).unwrap();
        assert_eq!(first.area_m2(), swapped.area_m2());
        assert_ne!(first.line_density_kg_m(), swapped.line_density_kg_m());
        assert_ne!(a.thermal().heat_capacity_per_length_j_m_k(), b.thermal().heat_capacity_per_length_j_m_k());
        let ra = a.evaluate(OhmicDrive::Voltage(1.0), 10.0).unwrap();
        let rb = b.evaluate(OhmicDrive::Voltage(1.0), 10.0).unwrap();
        assert_ne!(ra.axial_heat_flow_w, rb.axial_heat_flow_w);
        assert_ne!(ra.joule_power_w, rb.joule_power_w);
        assert_eq!(ra.layers[1], rb.layers[1]);
    }

    #[test]
    fn missing_wrong_quantity_and_invalid_resistivity_refuse() {
        use crate::conductor::{ELECTRICAL_RESISTIVITY_PROPERTY, OhmicDrive};
        let missing = source(8000.0, 400.0, 20.0, &[]);
        assert!(crate::conductor::ResolvedConductorState::try_from_material(missing).is_err());
        let wrong = source(8000.0, 400.0, 20.0,
            &[(ELECTRICAL_RESISTIVITY_PROPERTY, 1e-7, Pressure::DIMS)]);
        assert!(crate::conductor::ResolvedConductorState::try_from_material(wrong).is_err());
        for resistivity in [0.0, -1.0] {
            assert!(crate::conductor::ResolvedConductorState::try_from_material(
                electrical_source(8000.0, 400.0, 20.0, resistivity)).is_err());
        }
        let material = electrical_source(8000.0, 400.0, 20.0, 1e-7);
        let elastic = section(material.clone(), material);
        for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(elastic.thermal_transport().unwrap().with_ohmic_transport(length).is_err());
        }
        let coupled = elastic.thermal_transport().unwrap().with_ohmic_transport(1.0).unwrap();
        for drive in [OhmicDrive::Voltage(f64::NAN), OhmicDrive::Current(f64::INFINITY)] {
            assert!(coupled.evaluate(drive, 0.0).is_err());
        }
    }
}
