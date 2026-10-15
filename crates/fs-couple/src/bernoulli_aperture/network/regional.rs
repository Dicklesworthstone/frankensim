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
    /// This validates coherent inputs, not their empirical accuracy or gas phase.
    ///
    /// # Errors
    /// Incomplete/invalid states, unequal static pressures or inlet-speed mismatch.
    pub fn validate_section_gases(&self, gases: &[GasState]) -> Result<(), AcousticRealizeError> {
        if gases.len() != self.sections.len() || gases.is_empty() {
            return Err(invalid("regional duct requires exactly one gas state per physical section"));
        }
        for gas in gases {
            gas.validate_coherence()
                .map_err(|_| invalid("regional duct needs coherent material gas states"))?;
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

    /// Realize this physical graph without requiring a reed/valve source.
    /// Jets, prescribed sources and valves use the same section geometry and
    /// fs-vfit scattering/load owner. Every section retains its own gas-derived
    /// represented length and impedance. Initial waves and load states are zero.
    ///
    /// The payload budget covers the network, returned geometry and temporary
    /// segment records, not caller-owned graph/gas copies or process RSS.
    /// Cancellation is checked around the bounded network constructor and once
    /// per section; no partially constructed network is returned.
    ///
    /// # Errors
    /// Gas, topology, geometry, clock, memory, allocation or cancellation refusal.
    pub fn realize_waveguide(
        &self, gases: &[GasState], dt: f64, gate: &CancelGate,
    ) -> Result<(WaveguideNetwork, Vec<SectionRealization>), AcousticRealizeError> {
        let poll = || if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) };
        poll()?;
        self.validate_section_gases(gases)?;
        let bytes = self.sections.len().checked_mul(
            core::mem::size_of::<NetworkSegment>() + core::mem::size_of::<SectionRealization>())
            .ok_or_else(|| invalid("network geometry allocation size overflow"))?;
        let budget = self.max_wave_memory_bytes.checked_sub(bytes)
            .ok_or_else(|| invalid("network geometry payload exceeds memory allowance"))?;
        let mut represented = Vec::new();
        let mut segments = Vec::new();
        represented.try_reserve_exact(self.sections.len())
            .map_err(|_| invalid("network geometry allocation failed"))?;
        segments.try_reserve_exact(self.sections.len())
            .map_err(|_| invalid("network geometry allocation failed"))?;
        for (section, gas) in self.sections.iter().zip(gases) {
            poll()?;
            let actual = section.realize(gas.sound_speed, gas.density, dt)?;
            segments.push(NetworkSegment { nodes: section.nodes,
                one_way_samples: actual.one_way_samples,
                impedance_pa_s_m3: actual.impedance_pa_s_m3 });
            represented.push(actual);
        }
        let network = WaveguideNetwork::new(&self.nodes, &segments, dt, budget)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        poll()?;
        Ok((network, represented))
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

#[cfg(test)]
mod source_independent_tests {
    use super::*;
    use fs_material::gas::GasSpec;

    fn gas(t: f64, p: f64) -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), t, p).unwrap()
    }
    fn graph(g: &GasState) -> TubeNetworkSpec {
        TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: 0.0 }],
            sections: vec![TubeSection { nodes: [0,1], length_m: 0.3,
                radius_m: 0.006, max_length_error_m: 0.005 }],
            sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 100_000,
        }
    }
    #[test]
    fn same_geometry_owner_sets_impedance_and_transit_without_a_valve() {
        let g = gas(293.15, 101325.0);
        let spec = graph(&g);
        let gate = CancelGate::new_clock_free();
        let dt = 1.0/48000.0;
        let (mut net, sections) = spec.realize_waveguide(&[g], dt, &gate).unwrap();
        let z = spec.inlet_impedance(g.density).unwrap();
        assert_eq!(net.inlet_impedance_pa_s_m3(), z);
        assert_eq!(sections[0].one_way_samples, (0.3/(g.sound_speed*dt)).round() as usize);
        assert_eq!(sections[0].represented_length_m, sections[0].one_way_samples as f64*g.sound_speed*dt);
        for i in 0..=sections[0].one_way_samples {
            let f = net.step(if i == 0 {1.0} else {0.0}).unwrap();
            assert!(f.balance_residual_j().abs() < 1.0e-15);
        }
        assert!((net.node_frame(1).unwrap().net_flow_into_node_m3_s*z-1.0).abs() < 1.0e-12);
    }
    #[test]
    fn gas_changes_real_load_and_delay_instead_of_only_an_output_gain() {
        let gate = CancelGate::new_clock_free();
        let a = gas(280.0, 101325.0);
        let b = gas(320.0, 80000.0);
        let (_, ra) = graph(&a).realize_waveguide(&[a], 1.0/48000.0, &gate).unwrap();
        let (_, rb) = graph(&b).realize_waveguide(&[b], 1.0/48000.0, &gate).unwrap();
        assert_ne!(ra[0].impedance_pa_s_m3, rb[0].impedance_pa_s_m3);
        assert_ne!(ra[0].one_way_samples, rb[0].one_way_samples);
    }
    #[test]
    fn stale_gas_geometry_and_memory_refuse() {
        let gate = CancelGate::new_clock_free();
        let g = gas(293.15, 101325.0);
        let mut stale = g;
        stale.density *= 2.0;
        assert!(graph(&g).realize_waveguide(&[stale], 1.0/48000.0, &gate).is_err());
        let mut spec = graph(&g);
        spec.max_wave_memory_bytes = 0;
        assert!(spec.realize_waveguide(&[g], 1.0/48000.0, &gate).is_err());
        spec = graph(&g);
        spec.sections[0].max_length_error_m = 0.0;
        assert!(spec.realize_waveguide(&[g], 1.0/48000.0, &gate).is_err());
    }
}
