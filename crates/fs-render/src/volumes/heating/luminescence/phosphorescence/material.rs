//! Original material rates refreshed at the owning thermal cell's temperature.
//! All other physical coordinates and exact selected claims remain fixed.

use super::{
    HC_J_NM, LuminescenceSource, OpticalHeatingError, RateLaw, SingletTripletDecay,
    SingletTripletLimits, SingletTripletSource, invalid,
};
use crate::volumes::heating::fluorescence::FluorescenceSource;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::MaterialCard;
use fs_material::phosphorescence::{
    PHOSPHORESCENCE_CROSSING_RATE_PROPERTIES, PHOSPHORESCENCE_ENERGY_PROPERTIES,
    PHOSPHORESCENCE_NONRADIATIVE_RATE_PROPERTIES, PHOSPHORESCENCE_RADIATIVE_RATE_PROPERTIES,
    PhosphorescenceStatePoint,
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::sync::Arc;

#[derive(Debug)]
pub(super) struct MaterialTemperaturePhosphor {
    card: Arc<MaterialCard>,
    pub(super) original: PhosphorescenceStatePoint,
    temperature_axis: String,
    initial_temperature_k: f64,
    pub(super) identity: ContentHash,
}

impl<const M: usize> SingletTripletSource<M> {
    /// Bind original native energies and six competing rates directly, without
    /// converting them through an intrinsic lifetime or a scalar quantum yield.
    /// Exactly zero rates retain truly stable excited states.
    ///
    /// The two fluorescence inputs supply only source identities, validity
    /// bounds and conditional photon spectra; their scalar yields are ignored.
    /// Both use a common wavelength grid, with independent photon weights.
    /// The returned singlet descriptor must be used unchanged by the optical
    /// heating owner before binding the returned coupling.
    ///
    /// The named query axis must declare static absolute temperature in kelvin.
    /// Each thermal interval refreshes the original material claims at the
    /// actual cell temperature, retaining all other coordinates (including any
    /// quencher concentration). No pressure/concentration conversion or extra
    /// quenching channel is inferred. Both state energies must be independent
    /// of the owned temperature: this owner has no parameter-work channel.
    pub fn try_material_temperature(
        card: Arc<MaterialCard>,
        material: PhosphorescenceStatePoint,
        temperature_axis: &str,
        singlet_spectrum: FluorescenceSource<M>,
        triplet_spectrum: FluorescenceSource<M>,
        limits: SingletTripletLimits,
    ) -> Result<(LuminescenceSource<M>, Self), OpticalHeatingError> {
        let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if material.resolved().axis_quantities().get(temperature_axis)
            != Some(&absolute_temperature)
        {
            return Err(invalid(
                "material phosphorescence requires an explicitly typed absolute-temperature axis",
            ));
        }
        let initial_temperature_k = material
            .resolved()
            .state_coordinate(temperature_axis)
            .filter(|temperature| temperature.is_finite() && *temperature > 0.0)
            .ok_or_else(|| {
                invalid("material phosphorescence requires positive absolute temperature")
            })?;
        for name in [
            PHOSPHORESCENCE_ENERGY_PROPERTIES,
            PHOSPHORESCENCE_RADIATIVE_RATE_PROPERTIES,
            PHOSPHORESCENCE_NONRADIATIVE_RATE_PROPERTIES,
            PHOSPHORESCENCE_CROSSING_RATE_PROPERTIES,
        ]
        .into_iter()
        .flatten()
        {
            let property = material
                .resolved()
                .property(name)
                .ok_or_else(|| invalid("missing original phosphorescence rate or energy"))?;
            if !property.has_only_fixed_external_conditions(Some(temperature_axis), false) {
                return Err(invalid(
                    "native phosphorescence requires one thermal temperature owner and fixed external conditions",
                ));
            }
        }
        for name in PHOSPHORESCENCE_ENERGY_PROPERTIES {
            let property = material
                .resolved()
                .property(name)
                .ok_or_else(|| invalid("missing original phosphorescence state energy"))?;
            if property
                .source_axes()
                .iter()
                .any(|(axis, _)| axis == temperature_axis)
            {
                return Err(invalid(
                    "temperature-dependent state energies require an explicit parameter-work owner",
                ));
            }
        }
        // Refreshing at the original point verifies the immutable card and
        // exact claim pins before any initial-value fast path may bypass queries.
        material
            .at_query_point(&card, &material.resolved().query())
            .map_err(|_| {
                invalid("material phosphorescence requires its original source card and claims")
            })?;
        material
            .decay()
            .prepare(0.0, limits, || false)
            .map_err(|_| invalid("invalid singlet/triplet evolution limits"))?;
        let energies = material.decay().excitation_energies_j();
        let singlet =
            LuminescenceSource::try_native_spectrum(singlet_spectrum, HC_J_NM / energies[0])?;
        let triplet =
            LuminescenceSource::try_native_spectrum(triplet_spectrum, HC_J_NM / energies[1])?;
        if singlet.fluorescence().wavelengths_nm() != triplet.fluorescence().wavelengths_nm() {
            return Err(invalid(
                "material singlet/triplet emission requires a common wavelength grid",
            ));
        }
        singlet.check_temperature(initial_temperature_k)?;
        triplet.check_temperature(initial_temperature_k)?;
        let mut hash =
            DomainHasher::new("org.frankensim.render.material-temperature-phosphorescence.v1");
        hash.update(material.identity().as_bytes());
        hash.update(&(temperature_axis.len() as u64).to_le_bytes());
        hash.update(temperature_axis.as_bytes());
        let source = MaterialTemperaturePhosphor {
            card,
            original: material,
            temperature_axis: temperature_axis.to_owned(),
            initial_temperature_k,
            identity: hash.finalize(),
        };
        Ok((
            singlet,
            Self {
                triplet,
                rates: RateLaw::Material(Arc::new(source)),
                expected_singlet: Some(singlet),
                limits,
            },
        ))
    }

    /// Complete original snapshot, including fixed external coordinates and
    /// selected source claims. The thermal owner refreshes it without mutation.
    pub fn material_source(&self) -> Option<&PhosphorescenceStatePoint> {
        match &self.rates {
            RateLaw::Material(source) => Some(&source.original),
            RateLaw::Intrinsic { .. } => None,
        }
    }
}

impl MaterialTemperaturePhosphor {
    pub(super) fn model_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<SingletTripletDecay, OpticalHeatingError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid(
                "material phosphorescence requires positive absolute temperature",
            ));
        }
        if temperature_k.to_bits() == self.initial_temperature_k.to_bits() {
            return Ok(*self.original.decay());
        }
        let point = self
            .original
            .resolved()
            .query()
            .with(self.temperature_axis.as_str(), temperature_k)
            .map_err(|_| invalid("invalid material phosphorescence temperature refresh"))?;
        let refreshed = self
            .original
            .at_query_point(&self.card, &point)
            .map_err(|_| {
                invalid("original phosphorescence claims do not support the thermal state")
            })?;
        let energies = self.original.decay().excitation_energies_j();
        let refreshed_energies = refreshed.decay().excitation_energies_j();
        if energies
            .iter()
            .zip(refreshed_energies)
            .any(|(old, new)| old.to_bits() != new.to_bits())
        {
            return Err(invalid(
                "changing phosphorescence state energy requires an explicit parameter-work owner",
            ));
        }
        Ok(*refreshed.decay())
    }
}
