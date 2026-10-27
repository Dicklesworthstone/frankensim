//! The existing generalized-Maxwell material law in pressure/volume-flow units.
//! With f=h/(r*r) (including axial restraint) and wetted area A,
//! Z(s)=rho*h*s/A + E_inf*f/(A*s) + sum E_j*f/[A*(s+1/tau_j)].
//! Each last term is the shared parallel-RC impedance, NOT gas RL diffusion.
use super::{BoundMaterialWall, WallMaterialLoss, WallPatch, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use fs_blake3::ContentHash;
use fs_material::visco::GeneralizedMaxwell;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_matdb::ClaimId;
use fs_qty::{Dims, Pressure};
use fs_vfit::relaxation::{ComplianceRelaxationTerm, RelaxationImpedanceSpec, MAX_RELAXATION_BRANCHES};
use fs_vfit::waveguide::network::{NetworkNode, WaveguideNetwork};

mod normal;

/// Immutable source spectrum at the original wall's fixed material state.
/// This retains a supplied material law, not a fit from the desired sound.
#[derive(Clone, Debug)]
pub struct WallRelaxation {
    material: GeneralizedMaxwell,
    provenance: String,
    source_identity: Option<ContentHash>,
    source_arms: Option<Vec<SourceArm>>,
    normal_history: Option<normal::NormalPlaneHistory>,
    band_hz: [f64; 2],
    max_dt_over_tau: f64,
    max_angular_step: f64,
}

/// The identity of an internal variable comes from its original property pair
/// and selected claims, never from the current coefficient values or ordering.
#[derive(Clone, Debug, PartialEq, Eq)]
struct SourceArm {
    modulus_property: String,
    time_property: String,
    modulus_claim: ClaimId,
    time_claim: ClaimId,
}
impl PartialEq for WallRelaxation {
    fn eq(&self, other: &Self) -> bool {
        self.material.e_inf == other.material.e_inf && self.material.terms == other.material.terms
            && self.provenance == other.provenance && self.source_identity == other.source_identity
            && self.source_arms == other.source_arms
            && self.normal_history == other.normal_history
            && self.band_hz == other.band_hz
            && self.max_dt_over_tau == other.max_dt_over_tau && self.max_angular_step == other.max_angular_step
    }
}
impl WallRelaxation {
    /// Resolve a spectrum from the SAME material snapshot as the wall mass and
    /// equilibrium stiffness. Each pair names a modulus [Pa] and time [s].
    /// `young_modulus` is the equilibrium modulus, never the instantaneous one.
    /// Every selected coefficient, including a zero-modulus arm's time, must
    /// be constant at the declared temperature and cover the complete use band.
    /// The state identity retains the original claims, query and provenance.
    ///
    /// # Errors
    /// Missing, incorrectly dimensioned, nonconstant or out-of-band properties;
    /// duplicate branches, more than eight arms, or the original law/clock limits.
    pub fn from_state_point(material: &ResolvedMaterialStatePoint,
        properties: &[(String, String)], band_hz: [f64; 2],
        max_dt_over_tau: f64, max_angular_step: f64) -> Result<Self, AcousticRealizeError>
    {
        if properties.is_empty() || properties.len() > MAX_RELAXATION_BRANCHES {
            return Err(invalid("material wall needs 1..=8 sourced relaxation branches"));
        }
        let equilibrium = super::property(material, "young_modulus", Pressure::DIMS, band_hz)?;
        let mut terms = Vec::with_capacity(properties.len());
        let mut source_arms = Vec::with_capacity(properties.len());
        for (index, (modulus, time)) in properties.iter().enumerate() {
            if properties[..index].iter().any(|(e, tau)| e == modulus && tau == time) {
                return Err(invalid("material wall must not duplicate a sourced relaxation branch"));
            }
            terms.push((super::property(material, modulus, Pressure::DIMS, band_hz)?,
                super::property(material, time, Dims([0, 0, 1, 0, 0, 0]), band_hz)?));
            let selected = |name: &str| material.property(name)
                .map(|property| property.answer().receipt.selected)
                .ok_or_else(|| invalid("wall relaxation requires its original selected source claims"));
            source_arms.push(SourceArm {
                modulus_property: modulus.clone(), time_property: time.clone(),
                modulus_claim: selected(modulus)?, time_claim: selected(time)?,
            });
        }
        let law = GeneralizedMaxwell::new(equilibrium, terms)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let mut spectrum = Self::new(law,
            format!("material-state:{}; branches:{properties:?}", material.identity()),
            band_hz, max_dt_over_tau, max_angular_step)?;
        spectrum.source_identity = Some(material.identity());
        spectrum.source_arms = Some(source_arms);
        Ok(spectrum)
    }

    /// Admit an existing `GeneralizedMaxwell` law and retain its source/conditions.
    /// The wall card's Young modulus must be its EQUILIBRIUM modulus, with one
    /// constant Poisson ratio shared by all branches. No instantaneous modulus
    /// or separately authored modal damping is substituted or counted twice.
    /// `band_hz` is material applicability, not a claim that transient forcing
    /// has no out-of-band content. These source records do not update temperature.
    ///
    /// `max_dt_over_tau` is in (0,2], and `max_angular_step` in (0,1]. They bound
    /// midpoint memory resolution and the local instantaneous membrane frequency.
    /// Applications bind their actual clock through `validate_time_step`; the
    /// flute's material-network constructor does this before propagation setup.
    ///
    /// # Errors
    /// Invalid original model, more than eight supplied arms (including zero
    /// arms), missing source, invalid band or temporal allowances.
    pub fn new(material: GeneralizedMaxwell, provenance: String, band_hz: [f64; 2],
        max_dt_over_tau: f64, max_angular_step: f64) -> Result<Self, AcousticRealizeError>
    {
        if material.terms.len() > MAX_RELAXATION_BRANCHES || provenance.trim().is_empty()
            || !band_hz.iter().all(|x| x.is_finite() && *x >= 0.0) || band_hz[0] >= band_hz[1]
            || !(core::f64::consts::TAU * band_hz[1]).is_finite()
            || !max_dt_over_tau.is_finite() || !(0.0 < max_dt_over_tau && max_dt_over_tau <= 2.0)
            || !max_angular_step.is_finite() || !(0.0 < max_angular_step && max_angular_step <= 1.0) {
            return Err(invalid("wall relaxation needs a bounded original spectrum, source, band and explicit temporal limits"));
        }
        let material = GeneralizedMaxwell::new(material.e_inf, material.terms)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        Ok(Self { material, provenance, source_identity: None, source_arms: None,
            normal_history: None, band_hz, max_dt_over_tau, max_angular_step })
    }
    /// Exact original material law; zero arms remain in the source record.
    #[must_use]
    pub const fn material(&self) -> &GeneralizedMaxwell { &self.material }
    /// Source and conditions, or the explicit label of an authored model.
    #[must_use]
    pub fn provenance(&self) -> &str { &self.provenance }
    /// Exact source snapshot for a resolved spectrum; None for an authored law.
    #[must_use]
    pub const fn source_identity(&self) -> Option<ContentHash> { self.source_identity }
    /// Original material applicability [Hz].
    #[must_use]
    pub const fn band_hz(&self) -> [f64; 2] { self.band_hz }

    /// Exact source lineage on the same material card is checked by the wall
    /// owner. Here every declared arm must survive, even a zero arm; only the
    /// positive histories produce native storage-transfer entries.
    pub(super) fn history_correspondence(&self, replacement: &Self)
        -> Result<Vec<(usize, usize)>, AcousticRealizeError>
    {
        let (Some(old), Some(new)) = (&self.source_arms, &replacement.source_arms) else {
            return Err(invalid("live Maxwell history requires source-resolved arm identities"));
        };
        if old.len() != new.len() {
            return Err(invalid("live Maxwell update cannot create or discard source arms"));
        }
        let mut result = Vec::with_capacity(old.len());
        for (previous, identity) in old.iter().enumerate() {
            let next = new.iter().position(|candidate| candidate == identity)
                .ok_or_else(|| invalid("live Maxwell history must retain each original property pair and selected claims"))?;
            let positive = self.material.terms[previous].0 > 0.0;
            if positive != (replacement.material.terms[next].0 > 0.0) {
                return Err(invalid("live Maxwell update cannot activate or extinguish an arm"));
            }
            if positive { result.push((previous, next)); }
        }
        Ok(result)
    }

    pub(super) fn validate_source(&self, material: &ResolvedMaterialStatePoint)
        -> Result<(), AcousticRealizeError> {
        if self.source_identity.is_some_and(|identity| identity != material.identity()) {
            return Err(invalid("wall relaxation and inertia must use the same original material snapshot"));
        }
        Ok(())
    }

    pub(super) fn validate_material(&self, young: f64, band: [f64; 2])
        -> Result<(), AcousticRealizeError> {
        if self.material.e_inf.to_bits() != young.to_bits()
            || band[0] < self.band_hz[0] || band[1] > self.band_hz[1] {
            return Err(invalid("wall spectrum must match the original equilibrium modulus and cover the entire requested band"));
        }
        Ok(())
    }

    pub(super) fn load(&self, patch: WallPatch, young: f64, band: [f64; 2])
        -> Result<RelaxationImpedanceSpec, AcousticRealizeError>
    {
        self.validate_material(young, band)?;
        let mut terms = [ComplianceRelaxationTerm::default(); MAX_RELAXATION_BRANCHES];
        let mut count = 0;
        for &(modulus, tau) in &self.material.terms {
            if modulus == 0.0 { continue; } // Exact zero stores/dissipates nothing; positive arms never vanish.
            let stiffness = (patch.wall.stiffness_per_area / young) * modulus;
            let resistance = (stiffness / patch.area_m2) * tau;
            let rate = 1.0 / tau;
            if ![stiffness, resistance, rate].iter().all(|x| x.is_finite() && *x > 0.0) {
                return Err(invalid("wall relaxation geometry map lost a positive physical coefficient"));
            }
            terms[count] = ComplianceRelaxationTerm { resistance_pa_s_m3: resistance, rate_per_s: rate };
            count += 1;
        }
        RelaxationImpedanceSpec::new(patch.impedance()?, &[])
            .and_then(|load| load.with_compliance_terms(&terms[..count]))
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
    }
}

impl BoundMaterialWall {
    /// Complete physical shunt, including every positive hereditary branch.
    /// `patch()` alone describes ONLY the base inertial/equilibrium/viscous law.
    #[must_use]
    pub const fn shunt(&self) -> NetworkNode { NetworkNode::Shunt { load: self.load } }
    /// Full original-coordinate pressure/volume-flow impedance declaration.
    #[must_use]
    pub const fn load(&self) -> RelaxationImpedanceSpec { self.load }

    /// Validate the actual propagation clock against the source's time limits.
    /// This is a declared local accuracy admission, not a global error bound.
    ///
    /// # Errors
    /// Invalid dt, unresolved relaxation time or local instantaneous resonance.
    pub fn validate_time_step(&self, dt: f64) -> Result<(), AcousticRealizeError> {
        if !dt.is_finite() || dt <= 0.0 { return Err(invalid("material wall clock must be positive and finite")); }
        if let WallMaterialLoss::Maxwell(source) = &self.source.loss {
            let mut instantaneous = self.patch.wall.stiffness_per_area;
            let factor = instantaneous / source.material.e_inf;
            for &(modulus, tau) in &source.material.terms {
                if modulus == 0.0 { continue; }
                if dt / tau > source.max_dt_over_tau {
                    return Err(invalid("material wall clock does not resolve a retained relaxation time"));
                }
                instantaneous += factor * modulus;
            }
            let angular_step = dt * fs_math::det::sqrt(instantaneous / self.patch.wall.surface_density);
            if !angular_step.is_finite() || angular_step > source.max_angular_step {
                return Err(invalid("material wall clock does not resolve its instantaneous membrane frequency"));
            }
        }
        Ok(())
    }

    fn check_network(&self, network: &WaveguideNetwork) -> Result<(), AcousticRealizeError> {
        if network.node(self.source.node) != Some(self.shunt()) {
            return Err(invalid("material wall observation requires its complete original hereditary load"));
        }
        Ok(())
    }
    /// Actual accepted motion and BASE mechanical storage, excluding Maxwell arms.
    /// Total wall storage adds `relaxation_stored_energy_j`; both are already
    /// counted once by the common acoustic network's energy balance.
    ///
    /// # Errors
    /// Mismatched load, unavailable state or nonfinite motion.
    pub fn observe_waveguide(&self, network: &WaveguideNetwork) -> Result<super::super::WallState, AcousticRealizeError> {
        self.check_network(network)?;
        let state = network.terminal_state(self.source.node)
            .ok_or_else(|| invalid("material wall state is unavailable"))?;
        self.patch.motion(state)
    }
    /// Actual accepted energy in every stress-relaxation spring [J].
    /// This observes the shared pressure history; it does not step another law.
    ///
    /// # Errors
    /// Mismatched load, missing pressure states or nonfinite energy.
    pub fn relaxation_stored_energy_j(&self, network: &WaveguideNetwork) -> Result<f64, AcousticRealizeError> {
        self.check_network(network)?;
        let pressures = network.load_compliance_pressures(self.source.node)
            .ok_or_else(|| invalid("material wall compliance memory is unavailable"))?;
        let mut energy = 0.0;
        for (&p, term) in pressures.iter().zip(self.load.compliance_terms()) {
            let compliance = (1.0 / term.resistance_pa_s_m3) / term.rate_per_s;
            let scaled = fs_math::det::sqrt(compliance) * p;
            energy += 0.5 * scaled * scaled;
        }
        if !energy.is_finite() { return Err(invalid("material wall memory energy left the finite set")); }
        Ok(energy)
    }
}
