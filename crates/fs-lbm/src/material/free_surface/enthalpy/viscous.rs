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
    /// Force-corrected sqrt(2 D:D) [1/s], using this source's collision tau.
    pub shear_rate_per_s: f64,
    /// Located carrier mass divided by retained reference density [m3].
    /// Zero-mass interface ghosts and numerical carry get no cell heat.
    pub liquid_volume_m3: f64,
    /// mu * shear_rate^2 [W/m3], held over one physical flow interval.
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

/// Conductivity, molecular heating and momentum on one common interval.
#[cfg(feature = "enthalpy-coupling")]
#[derive(Clone, Debug)]
pub struct FreeSurfaceViscousConductionReport {
    /// The thermal audit sees both inputs; its enclosing flow audit separates
    /// external and molecular energy.
    pub conduction: conduction::FreeSurfaceConductionReport,
    pub applied: Vec<FreeSurfaceViscousCell>,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    pub(super) fn viscous_sources(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<FreeSurfaceViscousCell>, Error> {
        poll(gate)?;
        let masses = self.flow.lattice.mass_snapshot();
        admit_masses(&masses, self.controls.max_cells)?;
        let (states, tau, _) = self.resolve_states(&self.h, &masses, gate)?;
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
            let material = states[&location].clone();
            let lattice_rate = self
                .flow
                .lattice
                .transport_shear_rate(location, tau[&(key, lane)])?;
            let rate = lattice_rate / dt;
            let mass = native_mass * mass_unit;
            let volume = mass / rho0;
            let power = (material.dynamic_viscosity_pa_s() * rate) * rate;
            let heat = (power * volume) * dt;
            if !rate.is_finite()
                || rate < 0.0
                || (lattice_rate != 0.0 && rate == 0.0)
                || !mass.is_finite()
                || mass <= 0.0
                || !volume.is_finite()
                || volume <= 0.0
                || !power.is_finite()
                || power < 0.0
                || (rate != 0.0 && power == 0.0)
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
    /// initial strain rate. Its initial source viscosity gives 2 mu D:D.
    /// This source is held over one interval, applied with authored heat, then
    /// the original owner refreshes viscosity/buoyancy and advances momentum
    /// and VOF enthalpy transport once. Heating can change collision viscosity;
    /// these are different epochs of a first-order source split. Strong
    /// thermal feedback requires refinement.
    ///
    /// Authored heat must exclude molecular work owned here. No pressure work,
    /// bulk-viscosity loss, surface energy or numerical damping is thermalized.
    /// This closes constitutive heat and extensive enthalpy, without asserting
    /// exact discrete kinetic-plus-thermal balance. Non-Newtonian memory and
    /// phase/volume change remain unsupported. Any refusal publishes nothing.
    pub fn step_with_viscous_heating(
        &mut self,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceViscousReport, Error> {
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
        let flow = self.advance_prepared_enthalpy(
            heated,
            external,
            viscous,
            0.0,
            absolute_heat + viscous,
            gate,
        )?;
        Ok(FreeSurfaceViscousReport { flow, applied })
    }

    /// Conductivity and authored heat share thermal substeps with the held
    /// initial molecular stress power. Source k(T), collision mu(T) and optional
    /// buoyancy refresh through their original owners; flow advances once.
    /// The same source-split and physical limits as viscous heating apply.
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
