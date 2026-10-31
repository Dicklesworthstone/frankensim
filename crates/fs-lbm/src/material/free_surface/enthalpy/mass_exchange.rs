//! Conservative material transfer on the original liquid inventory.
//!
//! A phase-change owner supplies signed mass, its liquid enthalpy in the SAME
//! caloric datum, and separately partitioned latent/sensible heat. This module
//! changes actual VOF fill and extensive enthalpy together; it never represents
//! condensation as heat with no added liquid, or silently treats a source
//! formation enthalpy as the native liquid's caloric reference.

use super::*;

/// One integrated interfacial transfer, positive into the liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceMassExchange {
    /// An already occupied native interface; no droplet is nucleated in gas.
    pub coordinates: [usize; 3],
    /// Signed requested liquid mass [kg]. Evaporation is negative.
    pub mass_kg: f64,
    /// Carried material enthalpy [J/kg] in the original liquid caloric datum.
    /// A condensation owner must explicitly align any gas formation reference.
    pub specific_enthalpy_j_kg: f64,
    /// Separately supplied heat [J], including the declared liquid share of
    /// latent heat. Carried enthalpy is not counted again here.
    pub heat_j: f64,
}

/// Absolute representation/balance bounds, not phase-transfer model accuracy.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceMassExchangeControls {
    /// Sum of absolute per-interface requested/represented mass defects [kg].
    pub mass_tolerance_kg: f64,
    /// Sum of absolute per-interface enthalpy residuals [J].
    pub energy_tolerance_j: f64,
    /// Sum of absolute component momentum defects [kg m/s].
    pub momentum_tolerance_kg_m_s: f64,
}

/// Actual original native-cell inventory and source changes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceMassExchangeSample {
    /// Original native Cartesian coordinate.
    pub coordinates: [usize; 3],
    /// Signed requested transfer [kg], retained separately from representation.
    pub requested_mass_kg: f64,
    /// Actual native mass change; the external owner must debit this amount.
    pub mass_change_kg: f64,
    /// Actual represented transfer times the supplied liquid enthalpy [J].
    pub carried_enthalpy_j: f64,
    /// Separately supplied heat, counted once [J].
    pub heat_j: f64,
    /// Accepted [initial, final] liquid specific enthalpy [J/kg].
    pub specific_enthalpy_j_kg: [f64; 2],
}

/// Closed original-liquid mass, enthalpy and carried-momentum accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceMassExchangeReport {
    /// The accepted clock is unchanged; the outer coupling owner advances the
    /// physical interval once through the normal native thermal/flow step.
    pub physical_time_s: f64,
    /// Actual represented liquid mass change, including every interface [kg].
    pub mass_change_kg: f64,
    /// Sum of absolute requested-versus-represented mass discrepancies [kg].
    pub absolute_mass_defect_kg: f64,
    /// Actual represented transfers times their incoming/outgoing enthalpies [J].
    pub carried_enthalpy_j: f64,
    /// Sum of separately supplied heat [J].
    pub heat_j: f64,
    /// Independent endpoint change, evaluated incrementally to retain small
    /// transfers on a large existing inventory or arbitrary caloric datum.
    pub enthalpy_change_j: f64,
    /// Sum of absolute local caloric balance defects, with no datum in the scale [J].
    pub absolute_energy_defect_j: f64,
    /// Momentum carried at the accepted local liquid velocity, not a force.
    pub carried_momentum_kg_m_s: [f64; 3],
    /// Actual endpoint change after refreshing the original force frame [kg m/s].
    pub momentum_change_kg_m_s: [f64; 3],
    /// Sum of absolute local component discrepancies from carried momentum [kg m/s].
    pub absolute_momentum_defect_kg_m_s: f64,
    /// Stable native-coordinate order, independently of caller enumeration.
    pub samples: Vec<FreeSurfaceMassExchangeSample>,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Apply one integrated phase transfer to real native interface inventory.
    ///
    /// This transaction does not advance flow or time. An outer coupling owner
    /// computes the physical flux/partition, evaluates this on a cloned liquid,
    /// debits the actual represented gas inventory and accepts all participants
    /// together. It then advances its thermal/flow interval exactly once.
    /// Neither a transfer rate, surface-temperature solve, gas energy nor a
    /// saturation law is inferred here. Only the liquid caloric reference is
    /// admitted: callers explicitly transform source formation enthalpies.
    ///
    /// Material crosses at the accepted local liquid velocity. Native density,
    /// pressure and nonequilibrium populations are retained; the reported dm*u
    /// is material momentum, not a recoil force. Relative gas/liquid velocity,
    /// nucleation, disappearance and transfer beyond current cell capacity need
    /// their own interface/momentum/remap closure. Such geometry refuses here.
    /// Thermocapillary geometry, buoyancy, rheology and original material sources
    /// are refreshed on the full candidate. A preexisting solved-surface stress
    /// or dynamic contact history cannot be bypassed by an instantaneous update.
    pub fn exchange_interface_mass(
        &mut self,
        exchanges: &[FreeSurfaceMassExchange],
        controls: FreeSurfaceMassExchangeControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceMassExchangeReport, Error> {
        poll(gate)?;
        if exchanges.len() > self.controls.max_cells
            || [
                controls.mass_tolerance_kg,
                controls.energy_tolerance_j,
                controls.momentum_tolerance_kg_m_s,
            ]
            .iter()
            .any(|t| !t.is_finite() || *t < 0.0)
        {
            return Err(Error::Invalid("invalid interface-transfer budget"));
        }
        if self.surface_temperature_forces || self.wetting.is_some() {
            return Err(Error::Invalid(
                "phase transfer needs the coupled surface-temperature/contact-history owner",
            ));
        }
        let mut ordered = BTreeMap::new();
        for transfer in exchanges {
            poll(gate)?;
            let location = self
                .location_at(transfer.coordinates)
                .ok_or(Error::Invalid(
                    "phase transfer needs an existing liquid interface",
                ))?;
            if !transfer.mass_kg.is_finite()
                || !transfer.specific_enthalpy_j_kg.is_finite()
                || !transfer.heat_j.is_finite()
                || ordered.insert(location, *transfer).is_some()
            {
                return Err(Error::Invalid(
                    "interface transfers must be finite and unique",
                ));
            }
        }
        let masses = self.flow.lattice.mass_snapshot();
        let mass_unit = self.flow.scaling.mass_unit_kg();
        let mut lattice = self.flow.lattice.clone();
        let mut h = self.h.clone();
        let mut report = FreeSurfaceMassExchangeReport {
            physical_time_s: self.physical_time_s(),
            mass_change_kg: 0.0,
            absolute_mass_defect_kg: 0.0,
            carried_enthalpy_j: 0.0,
            heat_j: 0.0,
            enthalpy_change_j: 0.0,
            absolute_energy_defect_j: 0.0,
            carried_momentum_kg_m_s: [0.0; 3],
            momentum_change_kg_m_s: [0.0; 3],
            absolute_momentum_defect_kg_m_s: 0.0,
            samples: Vec::with_capacity(ordered.len()),
        };
        for (&location, transfer) in &ordered {
            poll(gate)?;
            let old_mass = masses[&location];
            let new_mass = old_mass + transfer.mass_kg / mass_unit;
            let dm = (new_mass - old_mass) * mass_unit;
            if !dm.is_finite() || (transfer.mass_kg != 0.0 && dm == 0.0) {
                return Err(Error::Invalid(
                    "interface mass increment is unrepresentable",
                ));
            }
            lattice.transport_set_interface_mass(location, new_mass, gate)?;
            let new_mass_kg = new_mass * mass_unit;
            let old_h = self.h[&location];
            // The actual native dm owns the transfer. Center the update on the
            // old datum so a common caloric offset cannot manufacture heat.
            let relative_supply = dm * (transfer.specific_enthalpy_j_kg - old_h) + transfer.heat_j;
            let new_h = old_h + relative_supply / new_mass_kg;
            if !new_mass_kg.is_finite()
                || new_mass_kg <= 0.0
                || !new_h.is_finite()
                || !relative_supply.is_finite()
                || (relative_supply != 0.0 && new_h == old_h)
            {
                return Err(Error::Invalid(
                    "interface enthalpy increment is unrepresentable",
                ));
            }
            h.insert(location, new_h);
            report.mass_change_kg += dm;
            report.absolute_mass_defect_kg += (dm - transfer.mass_kg).abs();
            report.carried_enthalpy_j += dm * transfer.specific_enthalpy_j_kg;
            report.heat_j += transfer.heat_j;
            report.samples.push(FreeSurfaceMassExchangeSample {
                coordinates: transfer.coordinates,
                requested_mass_kg: transfer.mass_kg,
                mass_change_kg: dm,
                carried_enthalpy_j: dm * transfer.specific_enthalpy_j_kg,
                heat_j: transfer.heat_j,
                specific_enthalpy_j_kg: [old_h, new_h],
            });
        }
        let endpoint_masses = lattice.mass_snapshot();
        let resolved = self.resolve_states(&lattice, &h, &endpoint_masses, gate)?;
        if let Some(tensions) = resolved.surface_tensions {
            lattice.rebind_surface_tensions(tensions, gate)?;
        }
        if let Some(forces) = resolved.forces {
            lattice.rebind_force_field(forces, gate)?;
        }
        #[cfg(feature = "enthalpy-coupling")]
        if let Some(core) = &self.electrical {
            core.admit(&lattice, &resolved.states, gate)?;
        }
        // Read back every actual endpoint, including cells whose force frame
        // changed after a material refresh. No unchanged massive inventory is
        // subtracted from another large total to hide a lost small transfer.
        let endpoint_masses = lattice.mass_snapshot();
        for (&location, &old_mass) in &masses {
            poll(gate)?;
            let Some([x, y, z]) = FreeSurface3::transport_coordinates(location) else {
                continue;
            };
            let new_mass = endpoint_masses[&location] * mass_unit;
            let dm = (endpoint_masses[&location] - old_mass) * mass_unit;
            let dh = h[&location] - self.h[&location];
            let transfer = ordered.get(&location);
            let supplied_relative = transfer.map_or(0.0, |t| {
                dm * (t.specific_enthalpy_j_kg - self.h[&location]) + t.heat_j
            });
            report.enthalpy_change_j += dm * self.h[&location] + new_mass * dh;
            report.absolute_energy_defect_j += (new_mass * dh - supplied_relative).abs();
            let before = self
                .flow
                .lattice
                .sample_at(x, y, z)
                .expect("accepted liquid");
            let after = lattice.sample_at(x, y, z).expect("retained interface");
            for axis in 0..3 {
                let u = self.flow.scaling.velocity_si(before.velocity[axis])?;
                let du = self
                    .flow
                    .scaling
                    .velocity_si(after.velocity[axis] - before.velocity[axis])?;
                let carried = dm * u;
                let change = carried + new_mass * du;
                report.carried_momentum_kg_m_s[axis] += carried;
                report.momentum_change_kg_m_s[axis] += change;
                report.absolute_momentum_defect_kg_m_s += (change - carried).abs();
            }
        }
        if [
            report.mass_change_kg,
            report.absolute_mass_defect_kg,
            report.carried_enthalpy_j,
            report.heat_j,
            report.enthalpy_change_j,
            report.absolute_energy_defect_j,
            report.absolute_momentum_defect_kg_m_s,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || report
                .carried_momentum_kg_m_s
                .iter()
                .chain(&report.momentum_change_kg_m_s)
                .any(|v| !v.is_finite())
            || report.absolute_mass_defect_kg > controls.mass_tolerance_kg
            || report.absolute_energy_defect_j > controls.energy_tolerance_j
            || report.absolute_momentum_defect_kg_m_s > controls.momentum_tolerance_kg_m_s
        {
            return Err(Error::Invalid(
                "interface mass, enthalpy or momentum balance exceeded budget",
            ));
        }
        poll(gate)?;
        self.flow.lattice = lattice;
        self.h = h;
        self.states = resolved.states;
        self.rheology_samples = resolved.rheology;
        Ok(report)
    }
}
