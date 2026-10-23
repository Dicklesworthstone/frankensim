//! Per-key piano hammers from a felt pack, not from code.
//!
//! The source of record is Stulov's measured set of 88 new, unvoiced Abel grand
//! hammers (Stulov 2005, Acta Acust. united Ac. 91). Each hammer was fitted with
//! the four-parameter hereditary law
//!   F = F0 [u^p - (gamma/tau0) int_0^t u^p(xi) exp((xi - t)/tau0) dxi],
//! and the per-key regressions of that set are printed in Kartofelev & Stulov,
//! Acta Mech. 225 (2014) eqs (4)-(7), with the masses in Kartofelev's thesis
//! eq (1.9). The pack `packs/hammer-felt-abel-stulov2005.fsmatpk` transcribes
//! their coefficients; this module evaluates the printed forms and nothing
//! else. The printed fits that contradict themselves are named in the pack's
//! notes and are not used.
//!
//! The contact solvers play the slow-loading form K u^p + R d(u^p)/dt, with
//! K = F0 (1 - gamma) (Stulov & Kartofelev 2009 eq (4)) and R = F0 gamma tau0:
//! the hereditary law's first two terms in tau0. Against the exact law, for a
//! half-sine compression it is off by 0.15 % of the peak force for hammer 1 at
//! 4 ms, 0.48 % for hammer 40 at 2 ms and 1.9 % for hammer 85 at 0.8 ms
//! (`hammer_felt_tests`). The exact one-memory contact is not built yet.
//!
//! The set was measured as whole hammers against a rigid surface, so a unison
//! shares one law: each string takes K/unison and R/unison, and the forces sum
//! to the measured hammer when the strings compress alike.
use super::geometry::Course;
use fs_matdb::{MaterialCard, MaterialStateId, NormalizedPack, QueryPoint, SelectionPolicy};
use fs_material::{WoolFelt, visco::GeneralizedMaxwell};
use fs_math::det;
use fs_qty::Dims;

const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);
const MASS: Dims = Dims([0, 1, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

/// Reference strain of the per-course felt card. The stress law is the
/// hammer's force law spread over its felt patch, so the choice cancels.
const REFERENCE_STRAIN: f64 = 0.2;
/// Densification bound of the felt card, as for the other cards.
const DENSIFICATION_STRAIN: f64 = 0.8;

const STULOV_ABEL: &[u8] = include_bytes!("packs/hammer-felt-abel-stulov2005.fsmatpk");

/// One hammer of a set, in SI.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hammer {
    /// Hammer number n, 1 (A0) to 88 (C8).
    pub number: u8,
    pub mass_kg: f64,
    /// Compliance nonlinearity exponent p.
    pub exponent: f64,
    /// Instantaneous stiffness F0 [N/m^p].
    pub instantaneous_n_m_p: f64,
    /// Hereditary amplitude gamma, 0 < gamma < 1.
    pub hereditary_amplitude: f64,
    /// Relaxation time tau0 [s].
    pub relaxation_time_s: f64,
}

/// A felt pack's per-key hammer regressions, read through a material card.
pub struct FeltCard {
    card: MaterialCard,
    pack_id: String,
    content_hash: String,
}

/// Admit compiled pack bytes; truncated or corrupt bytes refuse.
pub fn load_pack(bytes: &[u8]) -> Result<NormalizedPack, String> {
    NormalizedPack::from_bytes(bytes).map_err(|e| format!("felt pack: {e}"))
}

/// The source of record: Stulov's 88 measured Abel hammers.
pub fn stulov_abel() -> Result<FeltCard, String> {
    FeltCard::from_pack(&load_pack(STULOV_ABEL)?)
}

impl FeltCard {
    pub fn from_pack(pack: &NormalizedPack) -> Result<Self, String> {
        let id = pack.pack_id();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("source-pack:{id}"),
                phase: "piano hammer felt".into(),
                process: "frozen source-pack realization".into(),
                revision: 0,
            },
            pack.claims().clone(),
            vec![],
        )
        .map_err(|e| format!("felt pack {id}: {e:?}"))?;
        Ok(Self {
            card,
            pack_id: id.to_owned(),
            content_hash: pack.content_hash().to_hex(),
        })
    }

    pub fn pack_id(&self) -> &str {
        &self.pack_id
    }

    pub fn content_hash(&self) -> &str {
        &self.content_hash
    }

    /// The hammer of key `midi`. Each coefficient is queried at its hammer
    /// number, so a key outside the pack's validity refuses, as does a missing
    /// coefficient or an evaluated law outside the model's admissible region.
    pub fn hammer(&self, midi: u8) -> Result<Hammer, String> {
        let id = &self.pack_id;
        if !(21..=108).contains(&midi) {
            return Err(format!("felt pack {id}: key {midi} outside A0..C8"));
        }
        let number = midi - 20;
        let n = f64::from(number);
        let point = QueryPoint::new()
            .with("hammer_number", n)
            .map_err(|e| format!("felt pack {id}: {e:?}"))?;
        let value = |name: &str, dims: Dims| -> Result<f64, String> {
            let answer = self
                .card
                .claims()
                .query(name, &point, SelectionPolicy::SingleClaimOnly)
                .map_err(|e| format!("felt pack {id}: {name} at hammer {number}: {e:?}"))?;
            let sample = &answer.evidence.value;
            if sample.dims != dims || !sample.value.is_finite() {
                return Err(format!(
                    "felt pack {id}: {name} has dims {:?}, value {}",
                    sample.dims.0, sample.value
                ));
            }
            Ok(sample.value)
        };
        let one = Dims::NONE;
        // Acta Mech. 225 eqs (4)-(7) and thesis eq (1.9), as printed.
        let force_at_1mm = value(
            "hammer_instantaneous_force_at_1mm_regression_prefactor",
            FORCE,
        )? * det::exp(
            value(
                "hammer_instantaneous_force_regression_exponential_rate",
                one,
            )? * n,
        );
        let exponent = value("hammer_compression_exponent_regression_constant", one)?
            + value(
                "hammer_compression_exponent_regression_linear_coefficient",
                one,
            )? * n;
        let hereditary_amplitude = value("hammer_hereditary_amplitude_regression_constant", one)?
            + value(
                "hammer_hereditary_amplitude_regression_linear_coefficient",
                one,
            )? * n;
        let relaxation_time_s = value("hammer_relaxation_time_regression_constant", TIME)?
            - value(
                "hammer_relaxation_time_regression_linear_coefficient_magnitude_negative",
                TIME,
            )? * n
            + value(
                "hammer_relaxation_time_regression_quadratic_coefficient",
                TIME,
            )? * n
                * n;
        let mass_kg = value("hammer_mass_regression_constant", MASS)?
            - value(
                "hammer_mass_regression_linear_coefficient_magnitude_negative",
                MASS,
            )? * n
            + value("hammer_mass_regression_quadratic_coefficient", MASS)? * n * n;
        // F = F0 u^p with u in mm, so F0 newtons is the force at 1 mm.
        let instantaneous_n_m_p = force_at_1mm * det::pow(1.0e3, exponent);
        let admissible = exponent > 1.0
            && hereditary_amplitude > 0.0
            && hereditary_amplitude < 1.0
            && relaxation_time_s > 0.0
            && mass_kg > 0.0
            && instantaneous_n_m_p > 0.0
            && instantaneous_n_m_p.is_finite();
        if !admissible {
            return Err(format!(
                "felt pack {id}: hammer {number} evaluates outside the hereditary model \
                 (p {exponent}, gamma {hereditary_amplitude}, tau0 {relaxation_time_s} s, \
                 mass {mass_kg} kg, F0 {instantaneous_n_m_p} N/m^p)"
            ));
        }
        Ok(Hammer {
            number,
            mass_kg,
            exponent,
            instantaneous_n_m_p,
            hereditary_amplitude,
            relaxation_time_s,
        })
    }
}

impl Hammer {
    /// Relaxed stiffness K = F0 (1 - gamma) [N/m^p].
    pub fn static_n_m_p(&self) -> f64 {
        self.instantaneous_n_m_p * (1.0 - self.hereditary_amplitude)
    }

    /// Rate coefficient R = F0 gamma tau0 [N s/m^p] of the whole hammer.
    pub fn rate_n_s_m_p(&self) -> f64 {
        self.instantaneous_n_m_p * self.hereditary_amplitude * self.relaxation_time_s
    }

    /// The course's conservative felt card: K u^p of the whole hammer spread
    /// over the course's felt patch, with no crush and no Prony terms. Its loss
    /// is `rate_per_string`, applied by the contact solver.
    pub fn material(&self, c: &Course) -> Result<(WoolFelt, GeneralizedMaxwell), String> {
        c.validate()?;
        if c.midi != self.number + 20 {
            return Err(format!(
                "hammer {} does not belong to key {}",
                self.number, c.midi
            ));
        }
        let p = self.exponent;
        let stress = self.static_n_m_p() * det::pow(REFERENCE_STRAIN * c.felt_thickness_m, p)
            / c.felt_area_m2;
        let law = WoolFelt::new(stress, REFERENCE_STRAIN, p, p, 0.0, DENSIFICATION_STRAIN)
            .map_err(|e| e.to_string())?;
        let elastic = GeneralizedMaxwell::new(p * stress / REFERENCE_STRAIN, vec![])
            .map_err(|e| e.to_string())?;
        Ok((law, elastic))
    }

    /// R on each string of the course [N s/m^p].
    pub fn rate_per_string(&self, c: &Course) -> f64 {
        self.rate_n_s_m_p() / c.unison as f64
    }
}

#[cfg(test)]
#[path = "hammer_felt_tests.rs"]
mod tests;
