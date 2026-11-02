//! Coherent film and substrate heat on the existing thermal-cell timeline.
//!
//! Preparation resolves each spectral lane once and retains the original
//! optical owner's positive loss integrals. Each authored film has one explicit
//! receiver. Opaque stacks additionally deposit substrate work; transmitting
//! stacks retain outgoing radiation without an implicit substrate heat cell.
//! Receiver cells may coincide; their powers add without duplicating incident,
//! reflected or transmitted energy.
//!
//! Each compiled response freezes its coefficients and geometry. Fixed sources
//! retain them over their declared temperature domains; stateful coating sources
//! replace the packet from the original native state before each needed interval.
//! One cell represents each whole film's absorbed power, not its depth profile.
//! Callers justify that thermal lumping and the optical/material association;
//! an optical source identity alone does not certify that association. Cooling,
//! phase change and conduction remain with the original thermal owners.
//! Absorbed power is assumed to thermalize within the resolved thermal interval;
//! fluorescence, carrier extraction and delayed nonthermal energy stores are
//! outside this heating model.

use super::response::{
    RESPONSE_ROUNDOFF, hash_geometry, hash_values, mix, polarized_power, response_hash,
    validate_packet,
};
use super::{
    CheckedSum, DepositionError, DepositionSimulation, DomainHasher, MaterialId,
    PreparedSurfaceSource, Receiver, SpectralPowerSchedule, SurfaceHeatingBudget,
    SurfaceHeatingError, SurfaceHeatingTolerance, SurfaceOpticalHeating, SurfacePowerResponse,
    invalid,
};
use crate::conductor::{CoherentCoatingStack, ConductorOptics};
use crate::dielectric::DielectricCoatingStack;
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::CellId;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Frozen, polarization-weighted spectral heat in authored film order.
#[derive(Clone, Debug, PartialEq)]
pub struct CoatedSurfacePowerResponse<const N: usize> {
    total: SurfacePowerResponse<N>,
    films: Box<[[f64; N]]>,
    substrate: [f64; N],
    identity: ContentHash,
}

impl<const N: usize> CoatedSurfacePowerResponse<N> {
    /// Compile the original conductor stack at one resolved incidence geometry.
    /// As with the native optical owner, this explicit stack replaces any
    /// coating attached to `substrate`. Fractions multiply intercepted watts;
    /// no area, cosine projection, radiance factor or camera exposure is added.
    pub fn from_coated_conductor(
        wavelengths_nm: [f64; N],
        substrate: &ConductorOptics,
        stack: CoherentCoatingStack<'_>,
        incident_eta: [f64; N],
        incidence_cosine: f64,
        s_power_fraction: [f64; N],
    ) -> Result<Self, SurfaceHeatingError> {
        validate_packet(&wavelengths_nm, &s_power_fraction)?;
        let mut hash = response_hash(b"resolved-coating-heat", substrate.content_identity())?;
        hash.update(stack.content_identity().as_bytes());
        hash_geometry(
            &mut hash,
            incidence_cosine,
            &incident_eta,
            &s_power_fraction,
        );
        let mut films = vec![[0.0; N]; stack.layers().len()].into_boxed_slice();
        let mut substrate_heat = [0.0; N];
        let mut reflected = [0.0; N];
        let mut absorbed = [0.0; N];
        for lane in 0..N {
            let power = stack.layer_power_balance(
                substrate,
                wavelengths_nm[lane],
                incident_eta[lane],
                incidence_cosine,
            )?;
            reflected[lane] = mix(power.total.reflected, s_power_fraction[lane]);
            absorbed[lane] = mix(power.total.absorbed, s_power_fraction[lane]);
            substrate_heat[lane] = mix(power.substrate_absorbed, s_power_fraction[lane]);
            let mut sum = CheckedSum::default();
            sum.add(substrate_heat[lane])?;
            for (film, fraction) in films.iter_mut().zip(power.layer_absorbed()) {
                film[lane] = mix(*fraction, s_power_fraction[lane]);
                sum.add(film[lane])?;
            }
            if (sum.value()? - absorbed[lane]).abs() > RESPONSE_ROUNDOFF {
                return Err(invalid(
                    "film and substrate heat must close the original optical balance",
                ));
            }
        }
        let total =
            SurfacePowerResponse::finish(hash, wavelengths_nm, reflected, [0.0; N], absorbed)?;
        let mut hash = DomainHasher::new("org.frankensim.render.coated-surface-heat-response.v1");
        hash.update(total.identity().as_bytes());
        hash.update(&(films.len() as u64).to_le_bytes());
        for fraction in films.iter().chain(core::iter::once(&substrate_heat)) {
            for value in fraction {
                if !value.is_finite() || !(0.0..=1.0).contains(value) {
                    return Err(invalid(
                        "finite nonnegative coating absorption fractions required",
                    ));
                }
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        Ok(Self {
            total,
            films,
            substrate: substrate_heat,
            identity: hash.finalize(),
        })
    }

    /// Original total reflection and absorption, counted once per source.
    pub fn total(&self) -> &SurfacePowerResponse<N> {
        &self.total
    }
    /// Fractions of incident power in incident-to-substrate film order.
    pub fn film_absorbed_fractions(&self) -> &[[f64; N]] {
        &self.films
    }
    /// Fraction of incident power absorbed in the opaque half space.
    pub fn substrate_absorbed_fraction(&self) -> &[f64; N] {
        &self.substrate
    }
    /// Original optical sources, incidence and complete power partition.
    pub fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Reflected/transmitted radiation and positive heat in finite dielectric films.
/// The lossless surrounding media receive no thermal absorption assignment.
#[derive(Clone, Debug, PartialEq)]
pub struct TransmittingSurfacePowerResponse<const N: usize> {
    total: SurfacePowerResponse<N>,
    films: Box<[[f64; N]]>,
    identity: ContentHash,
}

impl<const N: usize> TransmittingSurfacePowerResponse<N> {
    /// Resolve original coherent film work once per wavelength. Film order is
    /// the authored outside-to-inside order, including for reverse incidence.
    /// The caller reverses incident/exit indices and supplies the corresponding
    /// incidence angle; no camera radiance eta-squared factor enters these watts.
    ///
    /// Exactly lossless, zero-thickness and empty stacks deposit no heat.
    /// Weak absorption remains independent of the rounded R+T residual. The
    /// explicit source identity binds the surrounding media/incidence declaration.
    #[allow(clippy::too_many_arguments)]
    pub fn from_dielectric_stack(
        wavelengths_nm: [f64; N],
        layers: &[crate::conductor::CoatingLayer<'_>],
        incident_eta: [f64; N],
        transmitted_eta: [f64; N],
        incidence_cosine: f64,
        reverse: bool,
        s_power_fraction: [f64; N],
        source_identity: ContentHash,
    ) -> Result<Self, SurfaceHeatingError> {
        validate_packet(&wavelengths_nm, &s_power_fraction)?;
        // The supplied exterior indices describe these exact packet nodes;
        // they do not extend any original film's wavelength support.
        let stack = DielectricCoatingStack::try_mixed_in_band(
            layers,
            (wavelengths_nm[0], wavelengths_nm[N - 1]),
        )?;
        let mut hash = response_hash(b"transmitting-coating-stack", source_identity)?;
        hash.update(
            CoherentCoatingStack::try_new(layers)?
                .content_identity()
                .as_bytes(),
        );
        hash.update(&[u8::from(reverse)]);
        hash_geometry(
            &mut hash,
            incidence_cosine,
            &incident_eta,
            &s_power_fraction,
        );
        hash_values(&mut hash, transmitted_eta);
        let mut films = vec![[0.0; N]; layers.len()].into_boxed_slice();
        let mut reflected = [0.0; N];
        let mut transmitted = [0.0; N];
        let mut absorbed = [0.0; N];
        for lane in 0..N {
            let work = stack.layer_power_balance(
                wavelengths_nm[lane],
                reverse,
                incidence_cosine,
                incident_eta[lane],
                transmitted_eta[lane],
            )?;
            reflected[lane] = polarized_power(work.fresnel.reflection, s_power_fraction[lane])?;
            transmitted[lane] = polarized_power(work.fresnel.transmission, s_power_fraction[lane])?;
            absorbed[lane] = mix(work.absorbed, s_power_fraction[lane]);
            let mut sum = CheckedSum::default();
            for (film, fraction) in films.iter_mut().zip(work.layer_absorbed()) {
                film[lane] = mix(*fraction, s_power_fraction[lane]);
                sum.add(film[lane])?;
            }
            if (sum.value()? - absorbed[lane]).abs() > RESPONSE_ROUNDOFF {
                return Err(invalid(
                    "finite-film heat must close the original dielectric absorption",
                ));
            }
        }
        let total =
            SurfacePowerResponse::finish(hash, wavelengths_nm, reflected, transmitted, absorbed)?;
        let mut hash =
            DomainHasher::new("org.frankensim.render.transmitting-surface-heat-response.v1");
        hash.update(total.identity().as_bytes());
        hash.update(&(films.len() as u64).to_le_bytes());
        for fraction in &films {
            for value in fraction {
                if !value.is_finite() || !(0.0..=1.0).contains(value) {
                    return Err(invalid(
                        "finite nonnegative coating absorption fractions required",
                    ));
                }
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        Ok(Self {
            total,
            films,
            identity: hash.finalize(),
        })
    }

    /// Original reflection, transmission and finite-film absorption.
    pub fn total(&self) -> &SurfacePowerResponse<N> {
        &self.total
    }

    /// Fractions in authored physical film order, never incidence traversal order.
    pub fn film_absorbed_fractions(&self) -> &[[f64; N]] {
        &self.films
    }

    pub fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Explicit native cell and applicability of its frozen optical coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceHeatingReceiver {
    /// Original finite thermal cell; no optical cell is created.
    pub thermal_cell: CellId,
    /// Material identity already bound to that native cell.
    pub material_id: MaterialId,
    /// Supported positive Kelvin interval for the frozen optical response.
    pub temperature_bounds_k: (f64, f64),
}

/// One physical illumination with separate film and substrate heat destinations.
#[derive(Clone, Debug)]
pub struct CoatedSurfaceIlluminationSource<const N: usize> {
    /// Stable identity; sources reduce in ascending ID order.
    pub source_id: u64,
    /// Exactly one receiver per authored film, including lossless/zero films.
    pub film_receivers: Vec<SurfaceHeatingReceiver>,
    /// Native cell receiving absorption in the opaque substrate.
    pub substrate_receiver: SurfaceHeatingReceiver,
    /// Original optical sources and compiled coherent power fractions.
    pub response: CoatedSurfacePowerResponse<N>,
    /// Spectral watts already intercepted by the surface patch.
    pub power: SpectralPowerSchedule<N>,
}

/// One transmitting coating illumination with a receiver for each finite film.
/// Outgoing transmitted radiation remains with its original transport owner;
/// there is no absorbing substrate or implicit thermal receiver.
#[derive(Clone, Debug)]
pub struct TransmittingSurfaceIlluminationSource<const N: usize> {
    pub source_id: u64,
    /// Authored outside-to-inside order, including for reverse illumination.
    pub film_receivers: Vec<SurfaceHeatingReceiver>,
    pub response: TransmittingSurfacePowerResponse<N>,
    pub power: SpectralPowerSchedule<N>,
}

pub(super) struct PreparedCoatingSource<const N: usize> {
    pub(super) source_id: u64,
    pub(super) receiver: Option<SurfaceHeatingReceiver>,
    pub(super) bindings: Vec<(SurfaceHeatingReceiver, [f64; N])>,
    pub(super) response: SurfacePowerResponse<N>,
    pub(super) power: SpectralPowerSchedule<N>,
}

impl<const N: usize> PreparedCoatingSource<N> {
    pub(super) fn opaque(
        source: &CoatedSurfaceIlluminationSource<N>,
    ) -> Result<Self, SurfaceHeatingError> {
        if source.film_receivers.len() != source.response.films.len() {
            return Err(invalid(
                "exactly one thermal receiver per authored coating film required",
            ));
        }
        Ok(Self {
            source_id: source.source_id,
            receiver: Some(source.substrate_receiver),
            bindings: source
                .film_receivers
                .iter()
                .copied()
                .chain(core::iter::once(source.substrate_receiver))
                .zip(
                    source
                        .response
                        .films
                        .iter()
                        .copied()
                        .chain(core::iter::once(source.response.substrate)),
                )
                .collect(),
            response: source.response.total,
            power: source.power.clone(),
        })
    }

    pub(super) fn transmitting(
        source: &TransmittingSurfaceIlluminationSource<N>,
    ) -> Result<Self, SurfaceHeatingError> {
        if source.film_receivers.len() != source.response.films.len() {
            return Err(invalid(
                "exactly one thermal receiver per authored coating film required",
            ));
        }
        Ok(Self {
            source_id: source.source_id,
            receiver: None,
            bindings: source
                .film_receivers
                .iter()
                .copied()
                .zip(source.response.films.iter().copied())
                .collect(),
            response: source.response.total,
            power: source.power.clone(),
        })
    }
}

#[derive(Clone, Debug)]
pub(super) struct AbsorptionReceiver<const N: usize> {
    pub(super) cell: CellId,
    pub(super) absorbed: [f64; N],
}

impl<const N: usize> SurfaceOpticalHeating<N> {
    /// Reuse the native surface observer with explicit film/substrate receivers.
    /// Spectral preparation is outside the time loop. Native scheduling, total
    /// optical accounting, thermal sources, cancellation and rollback are shared
    /// with `try_new`; a layer partition never creates a second light source.
    pub fn try_new_coated(
        sources: &[CoatedSurfaceIlluminationSource<N>],
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        if sources.len() > budget.max_sources {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(sources.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for source in sources {
            checkpoint()?;
            retained.push(PreparedCoatingSource::opaque(source)?);
        }
        Self::try_new_partitioned(retained, thermal, budget, tolerance, checkpoint)
    }

    /// Deposit only finite-film absorption from transmitting interfaces. Each
    /// film keeps its authored receiver even under reverse incidence. A bare
    /// interface has no receivers, so reflected/transmitted light can be
    /// accounted without creating a fictitious thermal cell.
    pub fn try_new_transmitting(
        sources: &[TransmittingSurfaceIlluminationSource<N>],
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        if sources.len() > budget.max_sources {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(sources.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for source in sources {
            checkpoint()?;
            retained.push(PreparedCoatingSource::transmitting(source)?);
        }
        Self::try_new_partitioned(retained, thermal, budget, tolerance, checkpoint)
    }

    pub(super) fn try_new_partitioned(
        sources: Vec<PreparedCoatingSource<N>>,
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        if sources.len() > budget.max_sources {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(sources.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        let mut partitions = BTreeMap::new();
        let mut receivers = BTreeMap::<CellId, Receiver>::new();
        for source in sources {
            checkpoint()?;
            let mut partition = Vec::new();
            partition
                .try_reserve_exact(source.bindings.len())
                .map_err(|_| SurfaceHeatingError::Budget)?;
            for (receiver, fraction) in &source.bindings {
                checkpoint()?;
                let (lo, hi) = receiver.temperature_bounds_k;
                if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi < lo {
                    return Err(invalid(
                        "positive finite coating receiver temperature support required",
                    ));
                }
                let cell = thermal.cells().get(&receiver.thermal_cell).ok_or_else(|| {
                    invalid("every coating receiver must be an active native cell")
                })?;
                if cell.material_id != receiver.material_id {
                    return Err(invalid(
                        "coating receiver must retain its declared native material",
                    ));
                }
                let admitted = receivers.entry(receiver.thermal_cell).or_insert(Receiver {
                    material: receiver.material_id,
                    minimum_k: lo,
                    maximum_k: hi,
                });
                admitted.minimum_k = admitted.minimum_k.max(lo);
                admitted.maximum_k = admitted.maximum_k.min(hi);
                partition.push(AbsorptionReceiver {
                    cell: receiver.thermal_cell,
                    absorbed: *fraction,
                });
            }
            if partitions
                .insert(source.source_id, partition.into_boxed_slice())
                .is_some()
            {
                return Err(invalid("surface illumination source IDs must be unique"));
            }
            prepared.push(PreparedSurfaceSource {
                source_id: source.source_id,
                receiver: source.receiver,
                response: source.response,
                power: source.power,
            });
        }
        let mut result =
            Self::try_new_prepared(prepared, thermal, budget, tolerance, &mut checkpoint)?;
        let model = Arc::get_mut(&mut result.model).expect("new surface model is uniquely owned");
        let mut hash = DomainHasher::new("org.frankensim.render.coated-surface-optical-heating.v1");
        hash.update(model.identity.as_bytes());
        for (source_id, partition) in &partitions {
            checkpoint()?;
            hash.update(&source_id.to_le_bytes());
            hash.update(&(partition.len() as u64).to_le_bytes());
            for receiver in partition {
                hash.update(&receiver.cell.0.to_le_bytes());
                for value in receiver.absorbed {
                    hash.update(&value.to_bits().to_le_bytes());
                }
            }
        }
        for (id, receiver) in &receivers {
            hash.update(&id.0.to_le_bytes());
            hash.update(&receiver.material.0.to_le_bytes());
            hash.update(&receiver.minimum_k.to_bits().to_le_bytes());
            hash.update(&receiver.maximum_k.to_bits().to_le_bytes());
        }
        model.identity = hash.finalize();
        model.receivers = receivers;
        model.coating_receivers = partitions;
        result.validate_receivers(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }
}
