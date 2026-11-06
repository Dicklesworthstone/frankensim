//! Molecular liquid stress power uses the same strain estimator as fixed-domain
//! thermal flow. It is an explicit constitutive source, not a correction from
//! the difference of two numerical kinetic energies.

use super::*;

/// One accepted initial-state molecular power sample, before heat or advection.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceViscousCell {
    pub coordinates: [usize; 3],
    /// Original source tuple at the accepted caloric temperature.
    pub material: ResolvedLiquidState,
    /// Force-corrected sqrt(2 dev(D):dev(D)) [1/s], using this source's tau.
    pub shear_rate_per_s: f64,
    /// Absolute divergence [1/s], reconstructed with the actual bulk relaxation.
    pub divergence_magnitude_per_s: f64,
    /// Applied deviatoric molecular stress power [W/m3].
    pub shear_power_w_m3: f64,
    /// Applied volumetric molecular stress power [W/m3].
    pub bulk_power_w_m3: f64,
    /// Actual bulk coefficient [Pa s] from this same material response, or the
    /// selected uniform coefficient / default `2mu/3` tie.
    pub bulk_viscosity_pa_s: f64,
    /// Molecular trace relaxation used to reconstruct this sample.
    pub bulk_collision_tau: f64,
    /// Located carrier mass divided by retained reference density [m3].
    /// Zero-mass interface ghosts and numerical carry get no cell heat.
    pub liquid_volume_m3: f64,
    /// Total `2 mu dev(D):dev(D) + zeta (div u)^2` [W/m3], held over one interval.
    pub power_w_m3: f64,
    /// The same power times liquid volume and physical timestep [J].
    pub heat_j: f64,
}

/// Molecular source samples and the single accepted flow/enthalpy step.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceViscousReport {
    /// Authored and molecular heat remain distinct in the existing balance.
    pub flow: FreeSurfaceEnthalpyReport,
    pub applied: Vec<FreeSurfaceViscousCell>,
}

pub(super) struct PreparedViscousHeat {
    pub heated: BTreeMap<MassLocation3, f64>,
    pub external: f64,
    pub absolute_heat: f64,
    pub viscous: f64,
    pub applied: Vec<FreeSurfaceViscousCell>,
}

/// Conductivity, molecular heating and momentum on one common interval.
#[cfg(feature = "enthalpy-coupling")]
#[derive(Clone, Debug)]
pub struct FreeSurfaceViscousConductionReport {
    /// The thermal audit sees both inputs; its enclosing flow audit separates
    /// external and molecular energy.
    pub conduction: conduction::FreeSurfaceConductionReport,
    pub applied: Vec<FreeSurfaceViscousCell>,
}

// Both thermal integration paths use the same accepted source coefficients.
// `dilation` is sqrt(2/3)|div u|; preserve tied-model arithmetic exactly.
pub(super) fn molecular_power(
    shear: f64,
    dilation: f64,
    mu: f64,
    independent_bulk: Option<f64>,
) -> [f64; 2] {
    let bulk = match independent_bulk {
        None => (mu * dilation) * dilation,
        Some(zeta) => ((zeta * dilation) * dilation) * 1.5,
    };
    [(mu * shear) * shear, bulk]
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    pub(super) fn viscous_sources(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<FreeSurfaceViscousCell>, Error> {
        poll(gate)?;
        let masses = self.flow.lattice.mass_snapshot();
        admit_masses(&masses, self.controls.max_cells)?;
        let resolved = self.resolve_states(&self.flow.lattice, &self.h, &masses, gate)?;
        let mass_unit = self.flow.scaling.mass_unit_kg();
        let rho0 = self.flow.scaling.material().density_kg_m3();
        let dt = self.flow.scaling.timestep_s();
        let mut samples = Vec::new();
        for (&location, &native_mass) in &masses {
            poll(gate)?;
            let MassLocation3::Cell(key, lane) = location else {
                continue;
            };
            if native_mass == 0.0 {
                continue;
            }
            let material = resolved.states[&location].clone();
            let tau = resolved.tau[&(key, lane)];
            let independent_bulk = self.flow.scaling.bulk_viscosity_at(material.material())?;
            let bulk_rate = self.flow.scaling.bulk_relaxation_rate_at(material.material())?;
            let numerators = self.flow.lattice.transport_stress_numerators(location)?;
            let lattice_rate = numerators[0] / tau;
            let lattice_dilation = match bulk_rate {
                Some(rate) => numerators[1] * rate,
                None => numerators[1] / tau,
            };
            let rate = lattice_rate / dt;
            let dilation = lattice_dilation / dt;
            let mass = native_mass * mass_unit;
            let volume = mass / rho0;
            let mu = material.dynamic_viscosity_pa_s();
            let [shear_power, bulk_power] = molecular_power(
                rate, dilation, mu, independent_bulk,
            );
            let divergence = dilation * 1.5_f64.sqrt();
            let power = shear_power + bulk_power;
            let heat = (power * volume) * dt;
            if !rate.is_finite()
                || rate < 0.0
                || (lattice_rate != 0.0 && rate == 0.0)
                || !dilation.is_finite()
                || !divergence.is_finite()
                || (lattice_dilation != 0.0 && dilation == 0.0)
                || !mass.is_finite()
                || mass <= 0.0
                || !volume.is_finite()
                || volume <= 0.0
                || !power.is_finite()
                || power < 0.0
                || (rate != 0.0 && shear_power == 0.0)
                || (dilation != 0.0 && bulk_power == 0.0)
                || !heat.is_finite()
                || heat < 0.0
                || (power != 0.0 && heat == 0.0)
            {
                return Err(Error::Invalid(
                    "liquid molecular stress power is unrepresentable",
                ));
            }
            samples.push(FreeSurfaceViscousCell {
                coordinates: FreeSurface3::transport_coordinates(location).expect("cell"),
                material,
                shear_rate_per_s: rate,
                divergence_magnitude_per_s: divergence,
                shear_power_w_m3: shear_power,
                bulk_power_w_m3: bulk_power,
                bulk_viscosity_pa_s: independent_bulk.unwrap_or((2.0 / 3.0) * mu),
                bulk_collision_tau: bulk_rate.map_or(tau, |rate| 1.0 / rate),
                liquid_volume_m3: volume,
                power_w_m3: power,
                heat_j: heat,
            });
        }
        Ok(samples)
    }

    /// Convert molecular stress work into the actual moving liquid's enthalpy.
    ///
    /// Native force-corrected non-equilibrium stress supplies the accepted
    /// initial deviatoric and volumetric strain rates. Its initial source shear
    /// viscosity and selected physical bulk viscosity determine stress power.
    /// This source is held over one interval, applied with authored heat, then
    /// the original owner refreshes viscosity/buoyancy and advances momentum
    /// and VOF enthalpy transport once. Heating can change collision viscosity;
    /// these are different epochs of a first-order source split. Strong
    /// thermal feedback requires refinement.
    ///
    /// Authored heat excludes molecular work owned here. The default collision's
    /// `zeta=2mu/3` retains bulk heat; an explicitly bound independent coefficient
    /// uses its own collision rate and `zeta (div u)^2` exactly once.
    /// Pressure work, surface energy and numerical damping are not thermalized.
    /// This closes constitutive heat and extensive enthalpy, without asserting
    /// exact discrete kinetic-plus-thermal balance. Non-Newtonian memory and
    /// phase/volume change remain unsupported. Any refusal publishes nothing.
    pub fn step_with_viscous_heating(
        &mut self,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceViscousReport, Error> {
        let prepared = self.prepare_viscous_heat(heat, gate)?;
        let flow = self.advance_prepared_enthalpy(
            prepared.heated,
            prepared.external,
            prepared.viscous,
            0.0,
            prepared.absolute_heat,
            None,
            gate,
        )?;
        Ok(FreeSurfaceViscousReport {
            flow,
            applied: prepared.applied,
        })
    }

    // The same explicit molecular source preparation feeds stationary or
    // moving-wall transport. Neither path thermalizes boundary work directly.
    pub(super) fn prepare_viscous_heat(
        &self,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<PreparedViscousHeat, Error> {
        let applied = self.viscous_sources(gate)?;
        let (mut heated, external, absolute_heat) = self.prepare_heat(heat, gate)?;
        let masses = self.flow.lattice.mass_snapshot();
        let mass_unit = self.flow.scaling.mass_unit_kg();
        let mut viscous = 0.0;
        for source in &applied {
            poll(gate)?;
            if source.heat_j != 0.0 {
                let location = self.location_at(source.coordinates).expect("sampled cell");
                *heated.get_mut(&location).expect("sampled enthalpy") +=
                    source.heat_j / (masses[&location] * mass_unit);
            }
            viscous += source.heat_j;
        }
        if !viscous.is_finite() {
            return Err(Error::Invalid("liquid molecular heat sum overflows"));
        }
        Ok(PreparedViscousHeat {
            heated,
            external,
            viscous,
            absolute_heat: absolute_heat + viscous,
            applied,
        })
    }

    /// Conductivity and authored heat share thermal substeps with the held
    /// initial molecular stress power. Source k(T), collision mu(T) and optional
    /// buoyancy refresh through their original owners; flow advances once.
    /// The same source-split and physical limits as viscous heating apply.
    /// Use [`Self::step_with_conduction_and_viscous_feedback`] when molecular
    /// power must also refresh at each thermal temperature/rheology query.
    #[cfg(feature = "enthalpy-coupling")]
    pub fn step_with_conduction_and_viscous_heating(
        &mut self,
        faces: &[conduction::FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: conduction::FreeSurfaceConductionControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceViscousConductionReport, Error> {
        let applied = self.viscous_sources(gate)?;
        let conduction = self.advance_with_conduction(faces, heat, controls, &applied, gate)?;
        Ok(FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }
}
