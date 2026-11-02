//! Pressure-flow dissipation located in a retained through-gap viscosity profile.
//!
//! The original film pressure solve supplies each channel's total work. The
//! exact profile supplies the positive fractions integral((z-z0)^2/mu dz)/D.
//! No wall heat split, layer mass, temperature law or mixing is inferred here.
//! Each complete layer explicitly names a finite caloric destination, which
//! may be shared with another layer or with the existing thermal network.
//!
//! Profiles are prescribed and frozen. Their hydraulic equivalent must match
//! the actual channel coefficient at construction and before every deposit.
//! A common native Barus factor preserves the layer work fractions. Uniform
//! or channel-only liquid feedback cannot replace this retained layer model.

use fs_tribo::resistive_film::ResistiveFilm;
pub use fs_tribo::resistive_film::stratified::{GapViscosityProfile, ViscosityLayer};

use super::{CaloricReservoir, ImpactError, invalid};

/// Complete destinations for one original passage's retained viscosity layers.
///
/// The mechanical film must already carry this profile's hydraulic-equivalent
/// viscosity, for example through `StratifiedFilm::film()`. Destination order
/// is exactly the original lower-wall to upper-wall layer order. The same
/// store may receive several layers; every layer's work is deposited once.
#[derive(Clone, Debug, PartialEq)]
pub struct SqueezeFilmLayerHeatRoute {
    /// Existing pressure-island insertion index.
    pub film: usize,
    /// Existing passage insertion index within that island.
    pub channel: usize,
    /// Original Newtonian layer interfaces and zero-gauge viscosities.
    pub profile: GapViscosityProfile,
    /// One existing finite reservoir key per original profile layer.
    pub layer_destinations: Vec<String>,
}

#[derive(Clone)]
pub(super) struct LayerHeatRoute {
    pub(super) profile: GapViscosityProfile,
    pub(super) destinations: Vec<usize>,
    largest: usize,
}

impl LayerHeatRoute {
    pub(super) fn compile(
        film: &ResistiveFilm,
        route: SqueezeFilmLayerHeatRoute,
        stores: &[CaloricReservoir<'_>],
    ) -> Result<Self, ImpactError> {
        if route.layer_destinations.len() != route.profile.layers().len() {
            return Err(invalid("layer heat needs exactly one destination per original layer"));
        }
        let mut destinations = Vec::with_capacity(route.layer_destinations.len());
        for key in &route.layer_destinations {
            destinations.push(stores.iter().position(|store|store.key()==key.as_str())
                .ok_or_else(||invalid("layer heat names an absent finite reservoir"))?);
        }
        let fractions=route.profile.layer_pressure_work_fractions();
        let mut largest=0;
        for i in 1..fractions.len() {
            if fractions[i]>fractions[largest] {largest=i;}
        }
        let result=Self {profile:route.profile,destinations,largest};
        result.admit_operator(film,route.channel)?;
        Ok(result)
    }

    fn admit_operator(&self, film:&ResistiveFilm, channel:usize) -> Result<(),ImpactError> {
        if film.channel_viscosities_pa_s().get(channel).copied()
            != Some(self.profile.hydraulic_equivalent_viscosity_pa_s()) {
            return Err(invalid("layer heat profile does not match the mechanical channel viscosity"));
        }
        Ok(())
    }

    /// Stage only into the caller's scratch heat vector. The enclosing native
    /// accepted-step callback owns all publication and cancellation decisions.
    pub(super) fn deposit(&self, film:&ResistiveFilm, channel:usize,
        heat_j:f64, stores_j:&mut[f64]) -> Result<(),ImpactError> {
        self.admit_operator(film,channel)?;
        let fractions=self.profile.layer_pressure_work_fractions();
        let mut assigned=0.0;
        let mut correction=0.0;
        for (layer, (&fraction,&to)) in fractions.iter().zip(&self.destinations).enumerate() {
            if layer==self.largest {continue;}
            let amount=heat_j*fraction;
            if !amount.is_finite() || amount<0.0 || (heat_j>0.0 && amount==0.0) {
                return Err(invalid("positive layer heat is unrepresentable"));
            }
            stores_j[to]+=amount;
            let increment=amount-correction;
            let next=assigned+increment;
            correction=(next-assigned)-increment;
            assigned=next;
        }
        // Put only the floating-point partition remainder in the largest
        // physical contributor. This preserves tiny positive layer deposits
        // and avoids erasing a thin final layer by subtracting nearly equal
        // whole-channel energies. The retained physical fractions stay exact.
        let remainder=heat_j-assigned;
        if !remainder.is_finite() || remainder<0.0 || (heat_j>0.0 && remainder==0.0) {
            return Err(invalid("conservative layer heat partition is unrepresentable"));
        }
        stores_j[self.destinations[self.largest]]+=remainder;
        Ok(())
    }
}
