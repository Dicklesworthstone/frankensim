//! The existing lumped Helmholtz model in the reciprocal modal network.
//!
//! This is a coordinate/port lowering, not another cavity constitutive law or
//! time integrator. fs-phs owns compliance, neck end correction and the optional
//! resonance-pinned radiation resistance. The modal network owns all reactions,
//! stepping, dissipation and state publication. A panel, vented enclosure and
//! string-instrument body use the same element.
//!
//! Positive plate displacement and neck flow are OUTWARD from the cavity.
//! With modal outward areas A and neck area S, volume change is A.q + S.x.
//! Compression energy is (A.q + S.x)^2/(2 C). Its gradient gives BOTH panel
//! pressure forces and neck acceleration. It must not be added as a one-way
//! pressure drive or as independent springs on individual plate modes.
//!
//! Scope: one uniform small-signal pressure, one circular lumped neck, frozen
//! coherent gas and fixed geometry. No internal standing waves, mean flow,
//! thermal history or f-hole shape correction. `shared::SharedCavity` can attach
//! a supplied passive exterior multiport to the same coupled time solve.

use crate::acoustic_realize::AcousticRealizeError;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingError};
use fs_duct::LossModel;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_scenario::{HelmholtzCavity, NeckEnd};

/// Explicit loss choice, separate from the cavity's reversible storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CavityNeckLoss {
    /// Constant acoustic resistance [Pa s/m^3]. Zero is the lossless reference.
    /// Caller-authored data, admitted only for named test fixtures (BR-04);
    /// product paths use [`CavityNeckLoss::RadiationAtResonance`].
    Resistance(f64),
    /// The declared outer end's fs-phs compact radiation resistance plus the
    /// neck's exact fs-duct Bessel viscothermal resistance, both frozen at the
    /// isolated Helmholtz frequency. It is NOT broadband: M2.4-W3b retires
    /// the frozen value from acceptance. End corrections stay the Helmholtz
    /// owner's, counted once.
    RadiationAtResonance,
    /// The outer end opens onto a computed exterior load on the actual mouth
    /// (a port on the body's exterior skin), which owns that end's mass and
    /// radiation resistance. The neck keeps its length, its inner end
    /// correction and the same Bessel neck loss as
    /// [`CavityNeckLoss::RadiationAtResonance`], frozen at the same isolated
    /// frequency of the declared ends. The declared outer end enters only
    /// that frequency. Without the exterior load the neck has no outer mass.
    ExteriorMouth,
}

/// fs-phs terminations of a declared cavity neck (DEF-09).
#[must_use]
pub fn neck_ends(geometry: &HelmholtzCavity) -> fs_phs::NeckEnds {
    let flange = |end| match end {
        NeckEnd::Flanged => fs_phs::MouthFlange::Flanged,
        NeckEnd::Unflanged => fs_phs::MouthFlange::Unflanged,
    };
    fs_phs::NeckEnds { inner: flange(geometry.inner_end), outer: flange(geometry.outer_end) }
}

/// Series viscothermal resistance [Pa s/m^3] of a lumped neck at `omega`:
/// `Re(-i k Zc) L` of fs-duct's exact Bessel wave (`e^{-iωt}`, mass-like
/// `Im < 0`). It tends to Poiseuille `8 μ L / (π a⁴)` at small shear number.
///
/// # Errors
/// The duct owner's refusal, or a negative/non-finite result.
pub fn neck_viscothermal_resistance(gas: &GasState, radius_m: f64, length_m: f64, omega: f64)
    -> Result<f64, AcousticRealizeError>
{
    let wave = fs_duct::segment_wave(gas, radius_m, omega, LossModel::Bessel)
        .map_err(AcousticRealizeError::Duct)?;
    let resistance = (C64::new(0.0, -1.0) * wave.wavenumber * wave.characteristic_impedance).re
        * length_m;
    if !(resistance >= 0.0 && resistance.is_finite()) {
        return Err(invalid("neck viscothermal resistance must be finite and nonnegative"));
    }
    Ok(resistance)
}

/// Total neck damper of a declared cavity at its isolated resonance: the
/// outer end's compact radiation resistance plus the Bessel neck loss. A
/// mouth beyond the compact fit's `ka = 1` refuses; nothing is swallowed.
///
/// # Errors
/// Helmholtz admission, the radiation fit's refusal or the duct owner's.
pub fn resonance_neck_resistance(geometry: &HelmholtzCavity, gas: &GasState)
    -> Result<f64, AcousticRealizeError>
{
    let ends = neck_ends(geometry);
    let phs = |e: fs_phs::PhsError| AcousticRealizeError::Nonlinear(e.to_string());
    let omega = fs_phs::helmholtz_natural_frequency(geometry.volume_m3, geometry.neck_radius_m,
        geometry.neck_length_m, gas.density, gas.sound_speed, ends).map_err(phs)?;
    let (radiation, _) = fs_phs::compact_radiation_impedance(gas.density, gas.sound_speed,
        geometry.neck_radius_m, omega, ends.outer).map_err(phs)?;
    Ok(radiation + neck_viscothermal_resistance(gas, geometry.neck_radius_m,
        geometry.neck_length_m, omega)?)
}

/// Acoustic resistance [Pa s/m^3] of a sharp-edged cylindrical orifice in a rigid
/// plate at `omega`, with thin Stokes layers: Brandao and Schnitzer, J. Fluid Mech.
/// 2020 (arXiv:2001.04904), eq. 1.8, `Re Z = (rho omega / a) eps Theta / sqrt 2`,
/// `eps = sqrt(mu / (rho omega)) / a`, aspect ratio `h = L / (2 a)`. `Theta` is the
/// larger of the paper's two asymptotes, eq. 1.9 `(4 / pi) (h + 0.91)` and eq. 1.10
/// `(ln(pi / h) / pi + 1) / 2`, which cross near `h = 0.037`. It counts the hole's
/// own walls as well as the plate around the edge, so it stands in place of the
/// neck's Bessel wall resistance, not beside it.
///
/// # Errors
/// Non-positive inputs, or a hole outside the paper's limit (`eps < 1`, `h > eps`).
pub fn orifice_viscous_resistance(gas: &GasState, radius_m: f64, length_m: f64, omega: f64)
    -> Result<f64, AcousticRealizeError>
{
    if [radius_m, length_m, omega].iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err(invalid("orifice resistance needs finite positive radius, length and frequency"));
    }
    let eps = det::sqrt(gas.dynamic_viscosity / (gas.density * omega)) / radius_m;
    let h = length_m / (2.0 * radius_m);
    if !(eps.is_finite() && eps < 1.0 && h > eps) {
        return Err(invalid("orifice resistance is outside the thin Stokes layer limit"));
    }
    let pi = core::f64::consts::PI;
    let long = 4.0 / pi * (h + 0.91);
    let short = 0.5 * (det::ln(pi / h) / pi + 1.0);
    let resistance = gas.density * omega / radius_m * eps * long.max(short)
        / core::f64::consts::SQRT_2;
    if !(resistance.is_finite() && resistance > 0.0) {
        return Err(invalid("orifice resistance must be finite and positive"));
    }
    Ok(resistance)
}

/// Original physical inputs and their energy-equivalent modal coefficients.
#[derive(Clone, Debug)]
pub struct HelmholtzModalCavity {
    geometry: HelmholtzCavity,
    gas: GasState,
    loss: CavityNeckLoss,
    area: f64,
    mass: f64,
    stiffness: f64,
    resistance: f64,
}
fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn network_invalid(what: &'static str) -> ModalCouplingError {
    ModalCouplingError::Invalid { what }
}

/// Admit distinct frozen interior and exterior media around one unloaded
/// acoustic equilibrium. Composition and temperature may differ; reference
/// pressure must agree to roundoff. A finite pressure difference needs its own
/// static load/prestress solve, which this small-signal cavity does not supply.
/// This does not model species diffusion or mean heat/mass exchange at a mouth.
pub fn validate_cavity_media(interior: &GasState, exterior: &GasState)
    -> Result<(), AcousticRealizeError>
{
    for gas in [interior, exterior] {
        gas.validate_coherence().map_err(|e| AcousticRealizeError::Ambient(e.to_string()))?;
    }
    let pressure_scale = interior.pressure.max(exterior.pressure);
    if (interior.pressure-exterior.pressure).abs() > 64.0*f64::EPSILON*pressure_scale {
        return Err(invalid("cavity and exterior reference pressures must agree without an explicit static pressure-load solve"));
    }
    Ok(())
}

impl HelmholtzModalCavity {
    /// Lower the original linear pHS, preserving its exact selected parameters.
    /// The two unit-basis effort evaluations recover its quadratic coefficients
    /// algebraically; they are neither finite differences nor a fitted response.
    /// Gas source receipts remain with its producer; GasState alone adds none.
    pub fn new(geometry: HelmholtzCavity, gas: GasState, loss: CavityNeckLoss)
        -> Result<Self, AcousticRealizeError>
    {
        gas.validate_coherence().map_err(|e| AcousticRealizeError::Ambient(e.to_string()))?;
        if [geometry.volume_m3, geometry.neck_radius_m, geometry.neck_length_m]
            .iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(invalid("modal cavity requires finite positive volume and neck geometry"));
        }
        let ends = neck_ends(&geometry);
        let source = match loss {
            CavityNeckLoss::Resistance(r) => {
                if !r.is_finite() || r < 0.0 {
                    return Err(invalid("cavity acoustic resistance must be finite and nonnegative"));
                }
                fs_phs::helmholtz_resonator(geometry.volume_m3, geometry.neck_radius_m,
                    geometry.neck_length_m, gas.density, gas.sound_speed, r, ends)
            }
            CavityNeckLoss::RadiationAtResonance => {
                let omega = fs_phs::helmholtz_natural_frequency(geometry.volume_m3,
                    geometry.neck_radius_m, geometry.neck_length_m, gas.density,
                    gas.sound_speed, ends).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
                let neck = neck_viscothermal_resistance(&gas, geometry.neck_radius_m,
                    geometry.neck_length_m, omega)?;
                fs_phs::helmholtz_resonator_radiating(geometry.volume_m3, geometry.neck_radius_m,
                    geometry.neck_length_m, gas.density, gas.sound_speed, ends, neck)
            }
            CavityNeckLoss::ExteriorMouth => {
                let omega = fs_phs::helmholtz_natural_frequency(geometry.volume_m3,
                    geometry.neck_radius_m, geometry.neck_length_m, gas.density,
                    gas.sound_speed, ends).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
                let neck = neck_viscothermal_resistance(&gas, geometry.neck_radius_m,
                    geometry.neck_length_m, omega)?;
                fs_phs::helmholtz_resonator_exterior_mouth(geometry.volume_m3, geometry.neck_radius_m,
                    geometry.neck_length_m, gas.density, gas.sound_speed, ends.inner, neck)
            }
        }.map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let q = source.effort(&[1.0, 0.0]);
        let p = source.effort(&[0.0, 1.0]);
        let (j, r, g) = source.structure();
        if source.state_dim() != 2 || q[1] != 0.0 || p[0] != 0.0
            || j != [0.0, 1.0, -1.0, 0.0] || g != [0.0, 1.0]
            || r[..3] != [0.0, 0.0, 0.0] {
            return Err(invalid("Helmholtz source no longer has its documented volume/momentum basis"));
        }
        let area = core::f64::consts::PI * geometry.neck_radius_m * geometry.neck_radius_m;
        // Physical neck translation x has volume coordinate S*x and momentum
        // P/S. Thus M = M_ac*S^2, K = S^2/C, and D = R_ac*S^2.
        let area_squared = area * area;
        let mass = area_squared / p[1];
        let stiffness = area_squared * q[0];
        let resistance = area_squared * r[3];
        if [area, mass, stiffness, mass.recip(), det::sqrt(mass).recip()]
            .iter().any(|x| !x.is_finite() || *x <= 0.0)
            || !resistance.is_finite() || resistance < 0.0
            || (r[3] > 0.0 && resistance == 0.0) {
            return Err(invalid("modal cavity normalization is not representable"));
        }
        Ok(Self { geometry, gas, loss, area, mass, stiffness, resistance })
    }

    /// Replace an [`CavityNeckLoss::ExteriorMouth`] neck's Bessel wall resistance by
    /// [`orifice_viscous_resistance`] frozen at `freeze_hz`. The neck's mass and
    /// compressibility stay as they are: eq. 1.8's equal reactance is left out, about
    /// 0.6 % of a guitar neck's mass, because that mass enters a computed exterior
    /// load's aperture shape.
    ///
    /// # Errors
    /// Any other loss selection, or the resistance's own refusals.
    pub fn with_orifice_resistance(mut self, freeze_hz: f64) -> Result<Self, AcousticRealizeError> {
        if self.loss != CavityNeckLoss::ExteriorMouth {
            return Err(invalid("the orifice resistance belongs to an exterior-mouth neck"));
        }
        let r = orifice_viscous_resistance(&self.gas, self.geometry.neck_radius_m,
            self.geometry.neck_length_m, 2.0 * core::f64::consts::PI * freeze_hz)?;
        let resistance = self.area * self.area * r;
        if !(resistance.is_finite() && resistance > 0.0) {
            return Err(invalid("orifice resistance is not representable on the neck"));
        }
        self.resistance = resistance;
        Ok(self)
    }

    /// [`Self::with_orifice_resistance`] for a mouth of several measured holes in one
    /// plate, each `(area_m2, bore_m)`, such as a violin's two f-holes. Each hole is the
    /// circle of its own area through its own bore, with [`orifice_viscous_resistance`]
    /// at `freeze_hz`. The neck's volume flow divides among the holes as their compact
    /// conductances, `S / (L + 2 delta)` with [`fs_phs::neck_end_correction`] flanged,
    /// and each hole dissipates its share squared times its own resistance. The neck
    /// keeps its declared mass and compressibility: the holes set only the loss. The
    /// paper's evaluated law is for a circle, so a hole of another shape enters as the
    /// circle of its area; the example README compares a slender slot from the same
    /// quadrature, eq. 5.15.
    ///
    /// # Errors
    /// Any other loss selection, no holes, a non-positive or non-finite area or bore,
    /// or a hole's own refusals.
    pub fn with_orifice_holes(mut self, freeze_hz: f64, holes: &[(f64, f64)])
        -> Result<Self, AcousticRealizeError>
    {
        if self.loss != CavityNeckLoss::ExteriorMouth {
            return Err(invalid("the orifice resistance belongs to an exterior-mouth neck"));
        }
        if holes.is_empty() {
            return Err(invalid("orifice holes need at least one measured hole"));
        }
        let omega = 2.0 * core::f64::consts::PI * freeze_hz;
        let mut parts = Vec::with_capacity(holes.len());
        for &(area_m2, bore_m) in holes {
            if !(area_m2.is_finite() && area_m2 > 0.0 && bore_m.is_finite() && bore_m > 0.0) {
                return Err(invalid("each orifice hole needs a finite positive area and bore"));
            }
            let radius = det::sqrt(area_m2 / core::f64::consts::PI);
            let end = fs_phs::neck_end_correction(radius, fs_phs::MouthFlange::Flanged)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            let conductance = area_m2 / (bore_m + 2.0 * end);
            parts.push((conductance, orifice_viscous_resistance(&self.gas, radius, bore_m, omega)?));
        }
        let total: f64 = parts.iter().map(|part| part.0).sum();
        let r: f64 = parts.iter().map(|&(k, r)| (k / total) * (k / total) * r).sum();
        let resistance = self.area * self.area * r;
        if !(resistance.is_finite() && resistance > 0.0) {
            return Err(invalid("orifice resistance is not representable on the neck"));
        }
        self.resistance = resistance;
        Ok(self)
    }

    /// Original fixed cavity geometry.
    #[must_use] pub const fn geometry(&self) -> &HelmholtzCavity { &self.geometry }
    /// The same gas used for neck inertia, cavity compliance and selected loss.
    #[must_use] pub const fn gas(&self) -> &GasState { &self.gas }
    /// Check the medium on the other side of the actual mouth. Legacy compact
    /// end corrections include exterior inertia and loss evaluated in `gas`,
    /// so they cannot be rebound to a different atmosphere. An exterior-owned
    /// mouth retains only internal inertia/loss and hands the complete outer
    /// load to its explicitly bound exterior owner.
    ///
    /// The existing frozen Bessel resistance keeps its declared interior-gas
    /// reference frequency; this is not a broadband neck-loss model.
    pub fn validate_exterior_gas(&self, exterior: &GasState)
        -> Result<(), AcousticRealizeError>
    {
        validate_cavity_media(&self.gas, exterior)?;
        if &self.gas != exterior && self.loss != CavityNeckLoss::ExteriorMouth {
            return Err(invalid("distinct cavity media require an exterior-owned mouth; the compact neck end law owns one gas"));
        }
        Ok(())
    }
    /// Original loss-model selection; not an experimental validation claim.
    #[must_use] pub const fn loss(&self) -> CavityNeckLoss { self.loss }
    /// Physical neck mass [kg], including the original owner's end correction.
    #[must_use] pub const fn neck_mass_kg(&self) -> f64 { self.mass }
    /// Physical neck area [m^2].
    #[must_use] pub const fn neck_area_m2(&self) -> f64 { self.area }
    /// Stiffness on volume-change/neck-area [N/m].
    #[must_use] pub const fn compression_stiffness_n_m(&self) -> f64 { self.stiffness }
    /// Viscous resistance on physical neck translation [N s/m].
    #[must_use] pub const fn neck_resistance_n_s_m(&self) -> f64 { self.resistance }
    /// Isolated Helmholtz angular frequency [rad/s], not the coupled body modes.
    #[must_use] pub fn isolated_omega_rad_s(&self) -> f64 { det::sqrt(self.stiffness / self.mass) }
    /// Outward volume-displacement/velocity projection of the neck mode.
    #[must_use] pub fn neck_volume_shape_m2_per_sqrt_kg(&self) -> f64 { self.area / det::sqrt(self.mass) }

    /// Inertia ONLY: the compression spring and loss belong to the connections.
    /// A zero natural frequency is the shared owner's real translational mode,
    /// not an artificial low-frequency oscillator. Attach the returned links
    /// before stepping to obtain the complete cavity.
    pub fn neck_component(&self, sample_rate_hz: u32, budget: ModalAcousticTimeBudget)
        -> Result<ModalAcousticTimeModel, ModalCouplingError>
    {
        ModalAcousticTimeModel::try_new_with_free_coordinates(sample_rate_hz, vec![ModalAcousticMode {
            angular_frequency_rad_s: 0.0, damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        }], budget).map_err(|source| ModalCouplingError::Component { component: 0, source })
    }

    /// One collective compression spring and, when selected, one neck dashpot.
    /// Component indices address the caller's original modal bank. `areas` are
    /// signed OUTWARD volume projections [m^2/sqrt(kg)] of the body modes, not
    /// point-displacement shapes. Opposite modal signs retain their cancellation.
    /// The body's own supports/stiffness must already be in its component.
    /// A zero-volume mode is allowed: it does not excite this uniform pressure.
    pub fn connections(&self, body_component: usize, neck_component: usize,
        areas: &[f64]) -> Result<Vec<ModalConnection>, ModalCouplingError>
    {
        if body_component == neck_component || areas.is_empty() || areas.len() > 4096
            || areas.iter().any(|x| !x.is_finite()) {
            return Err(network_invalid("cavity needs distinct body/neck components and finite signed modal areas"));
        }
        let shapes: Vec<f64> = areas.iter().map(|a| a / self.area).collect();
        if shapes.iter().any(|x| !x.is_finite()) {
            return Err(network_invalid("cavity body volume projection overflows"));
        }
        let inv_root = det::sqrt(self.mass).recip();
        let mut links = vec![ModalConnection {
            left: ModalAttachment { component: body_component, shapes },
            // ModalConnection is left MINUS right. The negative neck shape
            // makes both outward displacements increase the SAME volume.
            right: ModalAttachment { component: neck_component, shapes: vec![-inv_root] },
            stiffness_n_m: self.stiffness, damping_n_s_m: 0.0, rest_extension_m: 0.0,
        }];
        if self.resistance != 0.0 {
            links.push(ModalConnection {
                left: ModalAttachment { component: neck_component, shapes: vec![inv_root] },
                right: ModalAttachment { component: neck_component, shapes: vec![0.0] },
                stiffness_n_m: 0.0, damping_n_s_m: self.resistance, rest_extension_m: 0.0,
            });
        }
        Ok(links)
    }
}

#[cfg(test)]
mod tests;

/// One shared gas volume for multiple independently retained moving surfaces.
pub mod shared;
/// Source-gas binding for the original spatial Helmholtz pressure solver.
pub mod interior;

