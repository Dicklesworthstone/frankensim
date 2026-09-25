//! Stationary piecewise-uniform gas on the existing acoustic graph (MR33).
//!
//! Each section derives Z=rho*c/A and transit=L/c from its own supplied state.
//! The existing parallel adaptor preserves pressure and volume-flow continuity;
//! it is not a second source or a pitch correction. See J. O. Smith, Physical
//! Audio Signal Processing, "Lossless Scattering" (parallel acoustic tubes).
//!
//! This is linear acoustics about a prescribed, common-pressure background.
//! Temperature/composition gradients are frozen. No mixing, heat evolution,
//! mean flow, density-advection or moving thermodynamic interface is implied.
//! GasState's phase/transport validity limitations still apply. A caller-built
//! state is not promoted to measured or material-card authority.
use super::*;

impl TubeNetworkSpec {
    /// The unique physical section attached to the degree-one inlet.
    ///
    /// # Errors
    /// Missing/repeated inlet or degree other than one.
    pub fn inlet_section_index(&self) -> Result<usize, AcousticRealizeError> {
        let mut inlets = self.nodes.iter().enumerate()
            .filter(|(_, node)| matches!(node, NetworkNode::Inlet));
        let inlet = inlets.next().ok_or_else(|| invalid("tube network requires one inlet"))?.0;
        if inlets.next().is_some() { return Err(invalid("tube network requires one inlet")); }
        let mut found = None;
        for (i, section) in self.sections.iter().enumerate() {
            for &node in &section.nodes {
                if node == inlet && found.replace(i).is_some() {
                    return Err(invalid("tube network inlet must meet exactly one section"));
                }
            }
        }
        found.ok_or_else(|| invalid("tube network inlet is disconnected"))
    }

    /// Check a complete ordered mapping before allocating or advancing physics.
    /// Equal static pressure excludes an unmodeled steady pressure-driven flow.
    /// The legacy speed field must name the actual inlet, never an ignored value.
    /// This validates finite inputs, not their empirical accuracy or gas phase.
    ///
    /// # Errors
    /// Incomplete/invalid states, unequal static pressures or inlet-speed mismatch.
    pub fn validate_section_gases(&self, gases: &[GasState]) -> Result<(), AcousticRealizeError> {
        if gases.len() != self.sections.len() || gases.is_empty() {
            return Err(invalid("regional duct requires exactly one gas state per physical section"));
        }
        for gas in gases {
            if ![gas.temperature, gas.pressure, gas.density, gas.sound_speed,
                gas.dynamic_viscosity, gas.thermal_conductivity, gas.specific_gas_constant,
                gas.specific_heat_cp, gas.prandtl, gas.characteristic_impedance]
                .iter().all(|v| v.is_finite() && *v > 0.0)
                || !gas.gamma.is_finite() || gas.gamma <= 1.0
                || !gas.water_mole_fraction.is_finite() || !(0.0..1.0).contains(&gas.water_mole_fraction)
                || gas.pressure.to_bits() != gases[0].pressure.to_bits()
            { return Err(invalid("regional gas states must be finite and share one static pressure")); }
        }
        if self.sound_speed_m_s.to_bits() != gases[self.inlet_section_index()?].sound_speed.to_bits() {
            return Err(invalid("regional network speed declaration must equal its inlet gas speed"));
        }
        Ok(())
    }

    /// Actual inlet impedance for constructing the valve's characteristic port.
    /// The Bernoulli density must also be this inlet region's density.
    ///
    /// # Errors
    /// Same regional admission or invalid inlet geometry.
    pub fn inlet_impedance_with_gases(&self, gases: &[GasState]) -> Result<f64, AcousticRealizeError> {
        self.validate_section_gases(gases)?;
        let i = self.inlet_section_index()?;
        self.sections[i].uniform(gases[i].sound_speed).characteristic_impedance(gases[i].density)
    }
}

impl ApertureNetwork {
    /// Bind a complete frozen gas field to the existing reciprocal runtime.
    /// All scattering, storage, contact and step publication stay in their owners.
    /// Stored gas payload is included in the existing network memory allowance.
    /// No running valve can acquire a replacement gas field or empty waves.
    ///
    /// # Errors
    /// Invalid mapping, inlet density/impedance, geometry, topology or budget.
    pub fn with_section_gases(aperture: DynamicAperture, spec: TubeNetworkSpec,
        gases: Vec<GasState>) -> Result<Self, AcousticRealizeError>
    {
        spec.validate_section_gases(&gases)?;
        Self::build(aperture, spec, Some(gases))
    }

    /// Original complete regional states, absent on the legacy uniform path.
    /// They cannot mutate behind the propagation, receiver or material histories.
    #[must_use]
    pub fn section_gases(&self) -> Option<&[GasState]> { self.section_gases.as_deref() }

    /// Actual (density [kg/m3], sound speed [m/s]) at a section on either path.
    /// Observers must not substitute the inlet medium for a different outlet.
    #[must_use]
    pub fn section_medium(&self, section: usize) -> Option<(f64, f64)> {
        self.spec.sections.get(section)?;
        Some(self.section_gases.as_ref().map_or(
            (self.aperture.spec().density_kg_m3, self.spec.sound_speed_m_s),
            |gases| (gases[section].density, gases[section].sound_speed)))
    }
}
