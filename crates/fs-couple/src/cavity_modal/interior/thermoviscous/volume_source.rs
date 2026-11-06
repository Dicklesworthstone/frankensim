//! Original material-owned force and heat mechanisms in resolved acoustics.
//!
//! The same immutable fluid, molecular capacities and finite-solid sources
//! supply both frequency and time responses. Spatial source profiles remain
//! attached to their original FE owner; only their authored amplitudes vary.

use super::*;
use fs_feec::acoustics::thermoviscous::body_force::{
    MAX_BODY_FORCE_SOURCES, ThermoviscousBodyForceSource,
};
use fs_feec::acoustics::thermoviscous::time::PreparedThermoviscousVolumeTime;

impl<'c, S> ThermoviscousInterior<'c, S> {
    /// Integrate an original vector force profile [m^-3] once. A drive [N]
    /// multiplies this signed, unnormalized profile. Curved physical geometry,
    /// P2 velocity shapes and the original moving-wall lifting come from the FE.
    pub fn volume_force_source<'a>(
        &'a self,
        quadrature_order: usize,
        density_per_m3: &mut impl FnMut(usize, [f64; 3]) -> Result<[f64; 3], ThermoviscousError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousBodyForceSource<'a, 'c>, ThermoviscousSourceError> {
        Ok(self
            .fe
            .volume_force_source(quadrature_order, density_per_m3, cancelled)?)
    }

    /// Simultaneous physical forces, heat release and boundary drives, all in
    /// exp(+i omega t). The original source evaluates molecular Cp once before
    /// the same mixed solve. Signed force and heat work are counted once.
    #[allow(clippy::too_many_arguments)]
    pub fn response_with_volume_sources(
        &self,
        frequency_hz: f64,
        traction_pa: &[C64],
        wall_forces_n_per_sqrt_kg: &[C64],
        force_sources: &[(&ThermoviscousBodyForceSource<'_, '_>, C64)],
        heat_sources: &[(&ThermoviscousHeatSource<'_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousField, ThermoviscousSourceError> {
        Ok(self
            .prepare_response(frequency_hz, cancelled)?
            .response_with_volume_sources(
                traction_pa,
                wall_forces_n_per_sqrt_kg,
                force_sources,
                heat_sources,
                cancelled,
            )?)
    }

    /// Prepare held physical forces and heating on the complete original
    /// trajectory, retaining the source's molecular reservoirs and band.
    /// Finite conducting walls require explicit depth elements through
    /// `prepare_time_with_wall_diffusion_and_volume_sources`.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_time_with_volume_sources<'a>(
        &'a self,
        duration_s: f64,
        force_sources: &[&ThermoviscousBodyForceSource<'_, '_>],
        heat_sources: &[&ThermoviscousHeatSource<'_, '_>],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousVolumeTime<'a, 'c>, ThermoviscousSourceError> {
        self.prepare_volume_time_fields(
            duration_s,
            None,
            force_sources,
            heat_sources,
            limits,
            cancelled,
        )
    }

    /// Force/heat excitation with every original fluid, molecular, structural
    /// and finite-solid field. Depth counts follow `wall_sources()`; attached
    /// in-plane material conductivity is retained with its original identity.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_time_with_wall_diffusion_and_volume_sources<'a>(
        &'a self,
        duration_s: f64,
        normal_elements: &[usize],
        force_sources: &[&ThermoviscousBodyForceSource<'_, '_>],
        heat_sources: &[&ThermoviscousHeatSource<'_, '_>],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousVolumeTime<'a, 'c>, ThermoviscousSourceError> {
        self.prepare_volume_time_fields(
            duration_s,
            Some(normal_elements),
            force_sources,
            heat_sources,
            limits,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_volume_time_fields<'a>(
        &'a self,
        duration_s: f64,
        normal_elements: Option<&[usize]>,
        force_sources: &[&ThermoviscousBodyForceSource<'_, '_>],
        heat_sources: &[&ThermoviscousHeatSource<'_, '_>],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousVolumeTime<'a, 'c>, ThermoviscousSourceError> {
        let (frozen_cp, reservoirs) = self.time_calorics(duration_s, cancelled)?;
        let walls = normal_elements.map(|normal_elements| ThermoviscousWallDiscretization {
            normal_elements,
            maximum_frequency_hz: self.band_hz[1],
        });
        let transport = (!self.wall_tangential_sources.is_empty()).then(|| {
            self.wall_tangential_sources
                .iter()
                .map(|source| TangentialWallConductivity {
                    isotropic_w_m_k: source.conductivity_w_m_k(),
                })
                .collect::<Vec<_>>()
        });
        Ok(self.fe.prepare_time_with_volume_sources(
            duration_s,
            frozen_cp,
            &reservoirs,
            walls,
            transport.as_deref(),
            force_sources,
            heat_sources,
            limits,
            cancelled,
        )?)
    }

    /// Affine power scattering `b = S a + b_source` in the DUCT convention,
    /// exp(-i omega t). Returns (row-major S, outgoing source waves), retaining
    /// every original port and the caller's positive real reference resistance.
    ///
    /// Wall-force [N/sqrt(kg)], volume-force [N] and heat [W] phasors supplied
    /// here all use that same duct convention. Source amplitudes, material Cp,
    /// scattering and outgoing waves receive their explicit conjugations at
    /// the native FE boundary. The geometry and real spatial profiles do not.
    /// One loaded factor supplies S and the internal-source contribution.
    ///
    /// This is frequency-local linear composition; it does not infer a causal
    /// source feedback mechanism or a realization outside the source band.
    #[allow(clippy::too_many_arguments)]
    pub fn duct_scattering_with_volume_sources(
        &self,
        frequency_hz: f64,
        reference_resistance: &[f64],
        wall_forces_n_per_sqrt_kg: &[C64],
        force_sources: &[(&ThermoviscousBodyForceSource<'_, '_>, C64)],
        heat_sources: &[(&ThermoviscousHeatSource<'_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<C64>, Vec<C64>), ThermoviscousSourceError> {
        if cancelled() {
            return Err(ThermoviscousError::Cancelled.into());
        }
        if self.fe.port_mode_count() == 0 {
            return Err(invalid(
                "a closed cavity has no duct scattering coordinates",
            ));
        }
        if force_sources.len() > MAX_BODY_FORCE_SOURCES {
            return Err(ThermoviscousError::Budget {
                what: "simultaneous volume force sources",
                requested: force_sources.len(),
                limit: MAX_BODY_FORCE_SOURCES,
            }
            .into());
        }
        let wall_forces: Vec<_> = wall_forces_n_per_sqrt_kg
            .iter()
            .map(|force| force.conj())
            .collect();
        let mut forces = Vec::with_capacity(force_sources.len());
        for &(source, drive) in force_sources {
            if cancelled() {
                return Err(ThermoviscousError::Cancelled.into());
            }
            forces.push((source, drive.conj()));
        }
        let mut heat = Vec::with_capacity(heat_sources.len());
        for &(source, drive) in heat_sources {
            if cancelled() {
                return Err(ThermoviscousError::Cancelled.into());
            }
            heat.push((source, drive.conj()));
        }
        let prepared = self.prepare_scattering(frequency_hz, reference_resistance, cancelled)?;
        let mut outgoing = prepared.source_outgoing_with_volume_sources(
            &wall_forces,
            &forces,
            &heat,
            cancelled,
        )?;
        let mut scattering = prepared.scattering(cancelled)?;
        for value in scattering.iter_mut().chain(&mut outgoing) {
            *value = value.conj();
        }
        if cancelled() {
            return Err(ThermoviscousError::Cancelled.into());
        }
        Ok((scattering, outgoing))
    }
}
