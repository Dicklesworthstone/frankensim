//! Property-temperature queries retaining the original real-fluid constitutive law.
use super::*;

impl ThermochemicalPengRobinsonModel {
    /// Resolve a property state at the bulk pressure and declared homogeneous
    /// branch, retaining the bulk's active composition and exact source laws.
    ///
    /// This bank may contain additional species or different initial fractions.
    /// Both axes of its interaction matrix are restricted to the actual active
    /// species before comparing with the law retained by the original snapshot.
    /// A numerical match at one pressure/temperature is not source identity.
    /// Transport claims and fixed material coordinates use the existing source
    /// pinning owner. The result has no independent gas mass or stored energy.
    ///
    /// `phase` must explicitly choose Vapor or Liquid and reproduce the accepted
    /// bulk density. Stable-root selection is refused because it can change
    /// branches as property temperature changes. A unique bulk root admits
    /// either declared branch; the caller supplies its physical interpretation.
    /// The original mechanically/calorically stable homogeneous model still
    /// admits metastability and does not establish phase equilibrium or latent
    /// heat across the bulk-to-property temperature interval.
    pub fn evaluate_temperature_from_snapshot(
        &self,
        bulk: &ResolvedThermochemicalRealGasState,
        temperature: Temperature,
        phase: PengRobinsonPhaseV1,
    ) -> Result<ResolvedThermochemicalRealGasState, MaterialError> {
        if !temperature.value().is_finite() || temperature.value() <= 0.0 {
            return Err(invalid(
                "real-gas property temperature must be finite and positive",
            ));
        }
        if phase == PengRobinsonPhaseV1::Stable {
            return Err(invalid(
                "real-gas property queries require an explicit vapor or liquid branch",
            ));
        }
        let source = self.source.pinned_to_snapshot(bulk.source())?;
        // Ordinary coupled trials use the original bank and share this Arc.
        // Source/claim admission above still runs, but no PR vectors or pair
        // matrices are rebuilt for that common path.
        if !Arc::ptr_eq(&self.mechanics, &bulk.law) {
            let bank = self.source.thermodynamics().components();
            let active = source
                .thermodynamics()
                .components()
                .iter()
                .map(|component| {
                    bank.iter()
                        .position(|original| original.species() == component.species())
                        .ok_or_else(|| invalid("real-gas property bank lacks an active species"))
                })
                .collect::<Result<Vec<_>, _>>()?;
            let retained = bulk.constitutive_law();
            if active.len() != retained.components().len()
                || source.thermodynamics().mole_composition().fractions()
                    != retained.mole_fractions()
                || active
                    .iter()
                    .zip(retained.components())
                    .any(|(&i, original)| &self.mechanics.components()[i] != original)
            {
                return Err(invalid(
                    "real-gas property bank changes the original active constitutive parameters",
                ));
            }
            for (i, &row) in active.iter().enumerate() {
                for (j, &column) in active.iter().enumerate() {
                    if self.mechanics.binary_interactions()[row * bank.len() + column]
                        != retained.binary_interactions()[i * active.len() + j]
                    {
                        return Err(invalid(
                            "real-gas property bank changes the original pair interactions",
                        ));
                    }
                }
            }
        }
        let selected = bulk.mechanics().select(phase);
        if (selected.mass_density_kg_per_m3() / bulk.density_kg_m3() - 1.0).abs()
            > 16_384.0 * f64::EPSILON
        {
            return Err(invalid(
                "declared real-gas property branch differs from the accepted bulk root",
            ));
        }
        let point = bulk.source().components()[0]
            .query()
            .with("temperature", temperature.value())
            .map_err(|error| failure("real-gas property temperature coordinate", error))?;
        let source = source.evaluate(&point)?;
        let mechanics = bulk
            .law
            .evaluate_pt(temperature, Pressure::new(bulk.pressure_pa()))
            .map_err(|error| failure("real-gas property mechanics", error))?;
        let root = mechanics.select(phase).clone();
        // The admission above proves this is the original immutable law. Share
        // it with the new snapshot instead of allocating another retained copy.
        Self::resolve_real_gas(source, Arc::clone(&bulk.law), mechanics, root)
    }
}

#[cfg(test)]
mod tests;
