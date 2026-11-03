//! Molecular loads from the original accepted phase and reflection streams.
//!
//! Surfaces and the local Maxwellian gas have zero bulk velocity in one caller-
//! declared frame. Incident, emitted and reflected molecules all contribute to
//! normal traction, even at zero NET phase transfer. A fixed support absorbs
//! the solid impulse; this module does not add a velocity to `ClosedGasParcel`.
//! Its reciprocal gas impulse is a boundary transfer, not a claim that the
//! original internal-energy parcel evolved bulk momentum or kinetic energy.
//!
//! These are COMPLETE molecular loads, including carrier pressure. Do not add
//! another ambient-pressure load for the same molecules. Stationary supports
//! do zero work, and molecular kinetic energy is already in the original
//! h(T)-R*T/2 energy streams. Moving surfaces, Knudsen-layer drift, continuum
//! stress and non-Maxwellian scattering need their own momentum/energy owners.

/// Error-controlled complete-window impulses on the original adaptive chamber.
pub mod history;

use super::accommodation::{ReflectedHeatStep, ThermalAccommodation};
use super::chamber::{
    GasChamberCandidate, GasChamberSolve, GasChamberSurface, advance_chamber_observed,
    mixed::{CondensedPhaseStep, SharpMeltingCaloric, phase_reports},
};
use super::{
    ClosedGasParcel, GasSurfaceError, SpeciesId, SurfaceSnapshot, ThermochemicalGasModel, invalid,
    nonnegative, poll,
};
use fs_thermal::sublimation::molecular::momentum::{
    flux_weighted_normal_speed, maxwell_reflected_normal_pressure,
};
use std::collections::BTreeMap;

/// Explicit Maxwell scattering: each fraction is diffusely reemitted at the
/// surface temperature; the remainder is specular at the incident temperature.
/// This declaration also gives the existing energy-accommodation coefficient.
/// An arbitrary measured thermal coefficient alone does NOT imply this model.
#[derive(Clone, Debug, PartialEq)]
pub struct MaxwellReflection {
    thermal: ThermalAccommodation,
}
impl MaxwellReflection {
    pub fn from_diffuse_fractions(
        catalog: &ThermochemicalGasModel,
        fractions: &[(SpeciesId, f64)],
        temperature_range_k: [f64; 2],
    ) -> Result<Self, GasSurfaceError> {
        Ok(Self {
            thermal: ThermalAccommodation::new(catalog, fractions, temperature_range_k)?,
        })
    }

    /// The SAME explicitly declared Maxwell fractions, for the original heat
    /// owner. No independent thermal/scattering coefficient can be substituted.
    pub fn thermal_accommodation(&self) -> &ThermalAccommodation {
        &self.thermal
    }
}

/// Original exposed area comes from the corresponding `GasChamberSurface`.
/// This describes one planar patch, or an area-equivalent patch with the same
/// uniform traction resultant and moment. Coordinates share one fixed SI frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularSurfaceFrame {
    pub surface_id: u64,
    /// Unit normal directed FROM condensed material INTO the gas.
    pub outward_normal: [f64; 3],
    /// Centroid of this patch's uniformly exposed area [m].
    pub centroid_m: [f64; 3],
}

/// Integrated loads from exactly one accepted chamber interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularSurfaceImpulse {
    pub outgoing_normal_impulse_n_s: f64,
    pub sticking_normal_impulse_n_s: f64,
    /// Includes incident AND returning nonsticking molecules of every species.
    pub reflected_normal_impulse_n_s: f64,
    pub total_normal_impulse_n_s: f64,
    pub impulse_on_surface_n_s: [f64; 3],
    /// Moment about the candidate's explicit common origin [N m s].
    pub angular_impulse_on_surface_n_m_s: [f64; 3],
}
impl MolecularSurfaceImpulse {
    /// Reciprocal molecular boundary transfer; no bulk gas state is updated.
    pub fn impulse_into_gas_n_s(&self) -> [f64; 3] {
        self.impulse_on_surface_n_s.map(|x| -x)
    }

    /// Reaction required to keep this surface stationary, equal to the
    /// reciprocal molecular impulse and carrying exactly zero support work.
    pub fn fixed_support_impulse_n_s(&self) -> [f64; 3] {
        self.impulse_into_gas_n_s()
    }

    pub fn angular_impulse_into_gas_n_m_s(&self) -> [f64; 3] {
        self.angular_impulse_on_surface_n_m_s.map(|x| -x)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularChamberImpulse {
    pub duration_s: f64,
    pub origin_m: [f64; 3],
    pub surfaces: BTreeMap<u64, MolecularSurfaceImpulse>,
    pub impulse_on_surfaces_n_s: [f64; 3],
    pub angular_impulse_on_surfaces_n_m_s: [f64; 3],
}

/// Mass, thermal storage, optional fusion, reflected heat and mechanical loads
/// are published together after the ORIGINAL chamber's endpoint gates pass.
#[derive(Clone, Debug, PartialEq)]
pub struct MomentumChamberCandidate<'catalog> {
    exchange: GasChamberCandidate<'catalog>,
    reflection: BTreeMap<u64, ReflectedHeatStep>,
    phases: BTreeMap<u64, CondensedPhaseStep>,
    momentum: MolecularChamberImpulse,
}
impl<'catalog> MomentumChamberCandidate<'catalog> {
    pub fn exchange(&self) -> &GasChamberCandidate<'catalog> {
        &self.exchange
    }
    pub fn reflection(&self) -> &BTreeMap<u64, ReflectedHeatStep> {
        &self.reflection
    }
    pub fn phases(&self) -> &BTreeMap<u64, CondensedPhaseStep> {
        &self.phases
    }
    pub fn momentum(&self) -> &MolecularChamberImpulse {
        &self.momentum
    }
    #[allow(clippy::type_complexity)]
    pub fn into_parts(
        self,
    ) -> (
        GasChamberCandidate<'catalog>,
        BTreeMap<u64, ReflectedHeatStep>,
        BTreeMap<u64, CondensedPhaseStep>,
        MolecularChamberImpulse,
    ) {
        (self.exchange, self.reflection, self.phases, self.momentum)
    }
}

/// Joint phase exchange and stationary-surface force/torque impulse.
/// `None` retains the original specular, no-reflected-heat model. An explicit
/// Maxwell map must cover every selected surface and all retained gas species.
/// The original timestep, inventory, caloric and cancellation limits apply.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber_with_momentum<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    frames: &[MolecularSurfaceFrame],
    origin_m: [f64; 3],
    scattering: Option<&BTreeMap<u64, MaxwellReflection>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<MomentumChamberCandidate<'catalog>, GasSurfaceError> {
    advance(
        surfaces, frames, origin_m, scattering, None, gas, catalog, duration_s, controls, cancelled,
    )
}

/// The same original molecular loads with explicitly bound sharp fusion charts.
/// Area still means admitted exposed SOLID area. No load is inferred by
/// multiplying total surface area by a bulk solid fraction.
#[allow(clippy::too_many_arguments)]
pub fn advance_partitioned_chamber_with_momentum<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    calorics: &BTreeMap<u64, SharpMeltingCaloric<'_>>,
    frames: &[MolecularSurfaceFrame],
    origin_m: [f64; 3],
    scattering: Option<&BTreeMap<u64, MaxwellReflection>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<MomentumChamberCandidate<'catalog>, GasSurfaceError> {
    advance(
        surfaces,
        frames,
        origin_m,
        scattering,
        Some(calorics),
        gas,
        catalog,
        duration_s,
        controls,
        cancelled,
    )
}

// Compensated canonical sums retain opposing surface wrenches and weak carrier
// contributions. This does not renew any mass/energy tolerance of the chamber.
#[derive(Default)]
struct Sum(f64, f64);
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), GasSurfaceError> {
        let (sum, correction) = fs_math::eft::two_sum(self.0, term);
        self.0 = finite(sum)?;
        self.1 = finite(self.1 + correction)?;
        Ok(())
    }
    fn value(&self) -> Result<f64, GasSurfaceError> {
        finite(self.0 + self.1)
    }
}

fn finite(value: f64) -> Result<f64, GasSurfaceError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("molecular surface impulse is unrepresentable"))
    }
}

fn product(a: f64, b: f64) -> Result<f64, GasSurfaceError> {
    let value = finite(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        return Err(invalid("nonzero molecular impulse underflowed"));
    }
    Ok(value)
}

fn surface_impulse(
    surface: &GasChamberSurface<'_>,
    snapshot: &SurfaceSnapshot,
    frame: MolecularSurfaceFrame,
    origin: [f64; 3],
    incident: &BTreeMap<SpeciesId, f64>,
    scattering: Option<&MaxwellReflection>,
    duration: f64,
) -> Result<MolecularSurfaceImpulse, GasSurfaceError> {
    let emitted = product(
        snapshot.outgoing_mass_kg,
        flux_weighted_normal_speed(snapshot.gas.specific_gas_constant, snapshot.temperature_k)?,
    )?;
    let sticking = product(
        snapshot.incoming_mass_kg,
        flux_weighted_normal_speed(
            snapshot.gas.specific_gas_constant,
            snapshot.gas.temperature_k,
        )?,
    )?;
    let exposure = product(surface.solid.area_m2, duration)?;
    let mut reflected = Sum::default();
    for (species, pressure) in incident {
        let beta = if species == &surface.binding.species {
            surface.binding.accommodation_coefficient
        } else {
            0.0
        };
        let diffuse = match scattering {
            None => 0.0,
            Some(law) => *law
                .thermal
                .coefficients()
                .get(species)
                .ok_or(invalid("Maxwell scattering lacks an active gas species"))?,
        };
        let pressure = maxwell_reflected_normal_pressure(
            *pressure,
            snapshot.temperature_k,
            snapshot.gas.temperature_k,
            beta,
            diffuse,
        )?;
        reflected.add(product(pressure, exposure)?)?;
    }
    let reflected = reflected.value()?;
    let mut total = Sum::default();
    for term in [emitted, sticking, reflected] {
        total.add(term)?;
    }
    let total = total.value()?;
    let mut impulse = [0.0; 3];
    let mut arm = [0.0; 3];
    for axis in 0..3 {
        impulse[axis] = product(-frame.outward_normal[axis], total)?;
        arm[axis] = finite(frame.centroid_m[axis] - origin[axis])?;
    }
    let mut angular = [0.0; 3];
    for (axis, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
        angular[axis] = finite(arm[a].mul_add(impulse[b], -product(arm[b], impulse[a])?))?;
    }
    Ok(MolecularSurfaceImpulse {
        outgoing_normal_impulse_n_s: emitted,
        sticking_normal_impulse_n_s: sticking,
        reflected_normal_impulse_n_s: reflected,
        total_normal_impulse_n_s: total,
        impulse_on_surface_n_s: impulse,
        angular_impulse_on_surface_n_m_s: angular,
    })
}

#[allow(clippy::too_many_arguments)]
fn advance<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    frames: &[MolecularSurfaceFrame],
    origin: [f64; 3],
    scattering: Option<&BTreeMap<u64, MaxwellReflection>>,
    calorics: Option<&BTreeMap<u64, SharpMeltingCaloric<'_>>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration: f64,
    controls: GasChamberSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MomentumChamberCandidate<'catalog>, GasSurfaceError> {
    poll(&mut cancelled)?;
    if surfaces.len() > controls.maximum_surfaces
        || frames.len() != surfaces.len()
        || origin.iter().any(|x| !x.is_finite())
    {
        return Err(invalid(
            "bounded complete molecular surface frames and finite origin required",
        ));
    }
    let mut geometry = BTreeMap::new();
    for frame in frames {
        poll(&mut cancelled)?;
        let n = frame.outward_normal;
        if frame.centroid_m.iter().chain(&n).any(|x| !x.is_finite())
            || n.iter().any(|x| x.abs() > 1.0 + 128.0 * f64::EPSILON)
            || ((n[0] * n[0] + n[1] * n[1]) + n[2] * n[2] - 1.0).abs() > 128.0 * f64::EPSILON
            || geometry.insert(frame.surface_id, *frame).is_some()
        {
            return Err(invalid(
                "unique molecular surface identities and unit outward normals required",
            ));
        }
    }
    if surfaces
        .iter()
        .any(|surface| !geometry.contains_key(&surface.id))
    {
        return Err(invalid(
            "molecular geometry roster differs from selected surfaces",
        ));
    }
    let mut accommodation = BTreeMap::new();
    if let Some(scattering) = scattering {
        if scattering.len() != surfaces.len() {
            return Err(invalid(
                "Maxwell scattering must cover exactly the selected surfaces",
            ));
        }
        for (id, law) in scattering {
            poll(&mut cancelled)?;
            accommodation.insert(*id, law.thermal.clone());
        }
    }
    // Original active-gas composition is sampled ONCE, without any caloric
    // evaluation or inferred carrier defaults. Depleted species exert no load.
    let mut incident = BTreeMap::new();
    let model = gas.model().thermodynamics();
    for (component, fraction) in model
        .components()
        .iter()
        .zip(model.mole_composition().fractions())
    {
        poll(&mut cancelled)?;
        let pressure = *fraction * gas.state().state().pressure;
        if !nonnegative(pressure) {
            return Err(invalid(
                "original molecular partial pressure is unrepresentable",
            ));
        }
        incident.insert(component.species().clone(), pressure);
    }
    let mut loads = BTreeMap::new();
    let (exchange, reflection) = advance_chamber_observed(
        surfaces,
        scattering.map(|_| &accommodation),
        calorics,
        gas,
        catalog,
        duration,
        controls,
        |surface, snapshot| {
            let scattering = scattering.map(|laws| &laws[&surface.id]);
            let load = surface_impulse(
                surface,
                snapshot,
                geometry[&surface.id],
                origin,
                &incident,
                scattering,
                duration,
            )?;
            loads.insert(surface.id, load);
            Ok(())
        },
        &mut cancelled,
    )?;
    let phases = match calorics {
        Some(laws) => phase_reports(
            surfaces,
            laws,
            &exchange,
            controls.surface.gas.mass_tolerance_kg,
            &mut cancelled,
        )?,
        None => BTreeMap::new(),
    };
    let mut impulse: [Sum; 3] = core::array::from_fn(|_| Sum::default());
    let mut angular: [Sum; 3] = core::array::from_fn(|_| Sum::default());
    for load in loads.values() {
        poll(&mut cancelled)?;
        for axis in 0..3 {
            impulse[axis].add(load.impulse_on_surface_n_s[axis])?;
            angular[axis].add(load.angular_impulse_on_surface_n_m_s[axis])?;
        }
    }
    let mut impulse_value = [0.0; 3];
    let mut angular_value = [0.0; 3];
    for axis in 0..3 {
        impulse_value[axis] = impulse[axis].value()?;
        angular_value[axis] = angular[axis].value()?;
    }
    poll(&mut cancelled)?;
    Ok(MomentumChamberCandidate {
        exchange,
        reflection,
        phases,
        momentum: MolecularChamberImpulse {
            duration_s: duration,
            origin_m: origin,
            surfaces: loads,
            impulse_on_surfaces_n_s: impulse_value,
            angular_impulse_on_surfaces_n_m_s: angular_value,
        },
    })
}

#[cfg(test)]
mod tests;
