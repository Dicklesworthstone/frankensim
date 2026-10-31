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
use fs_phs::Storage;

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
impl Storage for GeometricStorage {
    fn hamiltonian(&self, q: &[f64]) -> f64 {
        self.0.potential(q)
    }
    fn gradient(&self, q: &[f64], out: &mut [f64]) {
        self.0.gradient(q, out);
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
                potential: Box::new(GeometricStorage(potential.geometric_remainder())),
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
