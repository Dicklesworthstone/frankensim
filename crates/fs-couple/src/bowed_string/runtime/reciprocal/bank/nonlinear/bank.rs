//! Material-derived stretching in the original bowed network and scheduler.
//!
//! The original exact linear owner keeps EVERY string, torsion, body, bridge,
//! cavity and cylinder-memory coordinate. This adapter adds only the existing
//! fs-nlmodal geometric stress channel, in that SAME normalized transverse
//! basis. No body or fluid component is rebuilt, omitted or stepped separately.
use super::super::{MaterialStringBankSource, MaterialStringSource, ReciprocalStringBank};
use super::super::multiple::MultiBowStringBank;
use super::{BodyPotential, ImpactError, StringStretching, WireSpan, invalid};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use crate::render::plate::impact::string::{StringObservation, StringPotential};
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::friction::nonlinear::{
    ModalPotential, ModalPotentialBinding, NonlinearModalConfig,
};
use fs_exec::CancelGate;
use fs_phs::{ConstitutiveDiscreteGradient, Storage};

/// Source-defined axial and winding-shear memory on the same mean strain.
pub mod relaxation;

/// Frozen source-derived geometric laws in original string declaration order.
/// Full linear-plus-quartic potentials are retained for diagnostics; only their
/// quartic remainders are added to the actual bank.
#[derive(Clone, Debug)]
pub struct StringBankStretching {
    potentials: Vec<StringPotential>,
    config: NonlinearModalConfig,
    material_relaxation: Option<relaxation::StringBankStretchingRelaxation>,
}
impl StringBankStretching {
    /// Original material stress channels, retaining the loaded mass basis.
    pub fn potentials(&self) -> &[StringPotential] {
        &self.potentials
    }
    /// Explicit bounded constitutive iteration selected before playback.
    pub const fn config(&self) -> NonlinearModalConfig {
        self.config
    }
    pub const fn material_relaxation(&self) -> Option<&relaxation::StringBankStretchingRelaxation> {
        self.material_relaxation.as_ref()
    }
}
impl MaterialStringBankSource {
    /// Source-owned stretching retained by this complete playable bank.
    pub const fn geometric_stretching(&self) -> Option<&StringBankStretching> {
        self.geometric_stretching.as_ref()
    }
}
struct GeometricStorage(StringPotential);
// The owner's own discrete gradient, never fs-phs's endpoint-subtraction
// fallback: from a static equilibrium a held step moves q by far less than q,
// and H(b) - H(a) then carries eps*|H| of rounding divided by |b - a|. The
// exact rule runs when the storage admits it; the bound remainder (zeroed
// frequencies) does not, so Gonzalez uses the owner's analytic remainder.
impl Storage for GeometricStorage {
    fn hamiltonian(&self, q: &[f64]) -> f64 {
        self.0.potential(q)
    }
    fn gradient(&self, q: &[f64], out: &mut [f64]) {
        self.0.gradient(q, out);
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        self.0.has_discrete_gradient_rule().then_some(self as &dyn ConstitutiveDiscreteGradient)
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], midpoint: &[f64]) -> Option<f64> {
        self.0.midpoint_energy_remainder(a, b, midpoint)
    }
}
impl ConstitutiveDiscreteGradient for GeometricStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], midpoint: &[f64], out: &mut [f64]) {
        self.0.discrete_gradient(a, b, midpoint, out);
    }
}
impl ModalPotential for GeometricStorage {
    fn coordinate_count(&self) -> usize {
        self.0.mode_count()
    }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        self.0.observe(q).map(|_| ()).map_err(|_| ModalCouplingError::Invalid {
            what: "source-derived string motion exceeds its finite continuous-span slope envelope",
        })
    }
}
impl ReciprocalStringBank {
    /// Add one shared transverse stretching energy per original material string.
    ///
    /// The material/section supplies EA. The already admitted card supplies the
    /// ACTUAL mass normalization, including instantaneous gas added inertia if
    /// selected. Both transverse planes contribute to one strain; torsion keeps
    /// its original solid polar inertia and force/torque coupling. Every memory
    /// branch, original body state, complete tensor loss, connection, pressure
    /// observer and material interface stays in the SAME FrictionModalSystem.
    ///
    /// Select gas loading and all geometric connections first, then bind this
    /// law once before playback. Existing initial coordinates are preserved;
    /// adding their geometric energy is an explicit initial condition. The
    /// original scheduled renderer executes the combined nonlinear bow equation.
    /// Additional bows and compliant stops may then be attached. Prescribed or
    /// contact-derived bow loads and normal support translation at fixed string
    /// stations share this same constitutive iteration and complete network.
    /// With released bows, static force preload solves this same strain potential
    /// and any compliant normal contacts against the complete connected network.
    ///
    /// Moderate-slope fixed-end Kirchhoff--Carrier mechanics; no axial wave
    /// inertia, dynamic temperature, helical winding slip or full rod. The
    /// bounded DG/linear iteration refuses nonconvergence; refine clock/modes.
    pub fn with_geometric_stretching(
        mut self,
        maximum_slope: f64,
        config: NonlinearModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if self.system.has_nonlinear_potentials() || self.system.samples_rendered() != 0 {
            return Err(invalid(
                "geometric stretching requires one source-bound attachment before playback",
            ));
        }
        let mut potentials = Vec::with_capacity(self.cards.len());
        let mut bindings = Vec::with_capacity(self.cards.len());
        for (index, (member, card)) in self.source.strings.iter().zip(&self.cards).enumerate() {
            if gate.is_requested() {
                return Err(ImpactError::Cancelled);
            }
            let planes = if matches!(&member.source, MaterialStringSource::WithPolarizations(_)) {
                2
            } else {
                1
            };
            let transverse = card
                .mode_count
                .checked_mul(planes)
                .ok_or_else(|| invalid("nonlinear material transverse prefix overflow"))?;
            let component = &self.system.components()[index];
            if component.states().len() != member.source.modes()
                || transverse > component.states().len()
            {
                return Err(invalid(
                    "nonlinear material string lost its original complete component layout",
                ));
            }
            let wire = WireSpan {
                endpoints_m: [[0.0, 0.0], [card.length_m, 0.0]],
                linear_density_kg_m: card.linear_density_kg_m,
                tension_n: card.tension_n,
                bending_stiffness_n_m2: card.bending_stiffness_n_m2,
                damping_per_s: (0..card.mode_count)
                    .map(|k| 2.0 * card.zetas[k] * card.mode_omega_rad_s(k + 1))
                    .collect(),
            };
            let stretching = StringStretching {
                axial_rigidity_n: member.source.axial_rigidity_n()?,
                maximum_slope,
            };
            let initial = component.states()[..transverse].to_vec();
            let body = if planes == 2 {
                wire.stretching_body_two_planes(initial, stretching)?
            } else {
                wire.stretching_body(initial, stretching)?
            };
            let BodyPotential::String(potential) = body.potential else {
                return Err(invalid(
                    "source-derived string changed its nonlinear storage owner",
                ));
            };
            if potential
                .omegas()
                .iter()
                .zip(component.modes())
                .any(|(a, b)| {
                    (a - b.angular_frequency_rad_s).abs()
                        > 32.0 * f64::EPSILON * a.abs().max(b.angular_frequency_rad_s.abs())
                })
            {
                return Err(invalid(
                    "geometric storage and retained loaded string use different linear bases",
                ));
            }
            bindings.push(ModalPotentialBinding {
                component: index,
                coordinate_start: 0,
                potential: Box::new(GeometricStorage(potential.geometric_remainder()?)),
            });
            potentials.push(potential);
        }
        self.system = self
            .system
            .with_nonlinear_potentials(bindings, config, gate)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        self.source.geometric_stretching = Some(StringBankStretching { potentials, config, material_relaxation: None });
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        Ok(self)
    }
    /// Original stress channel at accepted string coordinates; no state advance.
    pub fn stretching_observation(&self, string: usize) -> Result<StringObservation, ImpactError> {
        observe(&self.source, self.system.components(), &self.cards,
            &self.system.nonlinear_material_memory(), string)
    }
}
impl MultiBowStringBank {
    /// The same source-derived strain, tension and energy at the jointly
    /// accepted coordinates of this string, including every prescribed bow.
    pub fn stretching_observation(&self, string: usize) -> Result<StringObservation, ImpactError> {
        observe(self.source().bank(), self.system().components(), self.cards(),
            &self.system().nonlinear_material_memory(), string)
    }
}
fn observe(source: &MaterialStringBankSource, components: &[ModalAcousticTimeModel],
    cards: &[crate::bowed_string::BowedStringCard], memory: &[&[f64]], string: usize)
    -> Result<StringObservation, ImpactError>
{
    let potential = source.geometric_stretching.as_ref()
        .and_then(|binding| binding.potentials.get(string))
        .ok_or_else(|| invalid("string has no bound material geometric stretching law"))?;
    let q: Vec<_> = components[string].states()[..potential.mode_count()]
        .iter().map(|state| state.displacement_m_sqrt_kg).collect();
    let mut observation = potential.observe(&q)?;
    if let Some(memory_source) = source.geometric_stretching.as_ref()
        .and_then(|binding| binding.material_relaxation()) {
        memory_source.add_observation(string, cards[string].length_m, memory, &mut observation)?;
    }
    Ok(observation)
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use crate::modal_acoustic_time::ModalAcousticState;

    /// A fixed splitmix64 stream in [-1, 1), so every run checks the same pairs.
    pub(in super::super) fn signed(state: &mut u64) -> f64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        2.0 * ((z ^ (z >> 31)) >> 11) as f64 / (1_u64 << 53) as f64 - 1.0
    }

    /// Worst ratios, each admitted at <= 1, of the discrete gradient that the
    /// nonlinear owner evaluates (`fs_phs::discrete_gradient`) at `trials`
    /// points: (energy identity, continuity). For d = b - a, n coordinates and
    /// r = max|d| / max|a|, from 1 ulp to 1e-8 of max|a| (log-spaced), plus
    /// three large steps for the identity alone:
    ///   identity   |H(b) - H(a) - g.d| <= 8 (n+4) eps (|H(a)| + |H(b)| + sum|g_i d_i|);
    ///   continuity |g - gradH(m)| <= 64 (n+4) eps max|gradH(m)| + 16 r^2 max|gradH(m)|.
    pub(in super::super) fn discrete_gradient_ratios(storage: &dyn Storage, trials: usize,
        point: &mut dyn FnMut(&mut u64) -> Vec<f64>) -> (f64, f64)
    {
        let mut state = 0x5EED_0F_D15C_0E7E_u64;
        let (mut identity, mut continuity) = (0.0_f64, 0.0_f64);
        for _ in 0..trials {
            let a = point(&mut state);
            let n = a.len();
            let direction: Vec<f64> = (0..n).map(|_| signed(&mut state)).collect();
            let scale = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            let small = 17;
            let steps: Vec<f64> = (0..small).map(|k| 1e-8 * (f64::EPSILON / 1e-8).powf(k as f64 / (small - 1) as f64))
                .chain([1e-6, 1e-4, 1e-2]).collect();
            for (k, &step) in steps.iter().enumerate() {
                let b: Vec<f64> = a.iter().zip(&direction).map(|(a, d)| a + step * scale * d).collect();
                let d: Vec<f64> = a.iter().zip(&b).map(|(a, b)| b - a).collect();
                let r = d.iter().fold(0.0_f64, |m, v| m.max(v.abs())) / scale;
                if r == 0.0 { continue; }
                let g = fs_phs::discrete_gradient(storage, &a, &b);
                let (ha, hb) = (storage.hamiltonian(&a), storage.hamiltonian(&b));
                let work: f64 = g.iter().zip(&d).map(|(g, d)| g * d).sum();
                let size: f64 = g.iter().zip(&d).map(|(g, d)| (g * d).abs()).sum();
                let bound = 8.0 * (n + 4) as f64 * f64::EPSILON * (ha.abs() + hb.abs() + size);
                let ratio = (hb - ha - work).abs() / bound;
                assert!(ratio.is_finite(), "identity ratio is not finite at r {r:e}");
                identity = identity.max(ratio);
                if k < small {
                    let m: Vec<f64> = a.iter().zip(&b).map(|(&a, &b)| f64::midpoint(a, b)).collect();
                    let mut gm = vec![0.0; n];
                    storage.gradient(&m, &mut gm);
                    let top = gm.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
                    let gap = g.iter().zip(&gm).fold(0.0_f64, |s, (g, m)| s.max((g - m).abs()));
                    let bound = 64.0 * (n + 4) as f64 * f64::EPSILON * top + 16.0 * r * r * top;
                    assert!(top > 0.0, "the test point has no force");
                    continuity = continuity.max(gap / bound);
                }
            }
        }
        (identity, continuity)
    }

    const MODES: usize = 16;
    // A violin-like steel span (16 sine modes) and its full linear-plus-quartic potential.
    fn violin_span() -> StringPotential {
        let wire = WireSpan {
            endpoints_m: [[0.0, 0.0], [0.33, 0.0]],
            linear_density_kg_m: 6.0e-4,
            tension_n: 70.0,
            bending_stiffness_n_m2: 2.0e-5,
            damping_per_s: vec![0.0; MODES],
        };
        let body = wire.stretching_body(vec![ModalAcousticState::default(); MODES],
            StringStretching { axial_rigidity_n: 1.0e4, maximum_slope: 0.2 }).unwrap();
        let BodyPotential::String(potential) = body.potential else { panic!("the span lost its string storage") };
        potential
    }
    // About 0.3 mm of the fundamental, falling as 1/k over the modes.
    fn span_ratios(storage: &GeometricStorage) -> (f64, f64) {
        let mass = 0.5 * 6.0e-4 * 0.33_f64;
        discrete_gradient_ratios(storage, 24, &mut |state| {
            (0..MODES).map(|k| 3.0e-4 * mass.sqrt() * signed(state) / (k + 1) as f64).collect()
        })
    }

    // The bank binds the quartic remainder, with the linear frequencies left to the
    // exact linear owner. Its zeroed frequencies leave the owner's exact rule out, so
    // this is Gonzalez's form on the owner's analytic remainder.
    #[test]
    fn the_bound_geometric_remainder_has_an_exact_and_continuous_discrete_gradient() {
        let storage = GeometricStorage(violin_span().geometric_remainder().unwrap());
        assert!(storage.discrete_gradient_rule().is_none(), "the zero-frequency remainder took the exact rule");
        let (identity, continuity) = span_ratios(&storage);
        assert!(identity <= 1.0, "energy identity ratio {identity:e}");
        assert!(continuity <= 1.0, "continuity ratio {continuity:e}");
    }

    // The full potential keeps its positive frequencies, so the owner's exact rule
    // runs, through the adapter's [q, p] layout.
    #[test]
    fn the_full_string_potential_takes_the_owners_exact_discrete_gradient() {
        let storage = GeometricStorage(violin_span());
        assert!(storage.discrete_gradient_rule().is_some(), "the owner's exact rule was not forwarded");
        let (identity, continuity) = span_ratios(&storage);
        assert!(identity <= 1.0, "energy identity ratio {identity:e}");
        assert!(continuity <= 1.0, "continuity ratio {continuity:e}");
    }
}
