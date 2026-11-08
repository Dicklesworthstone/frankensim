//! Prescribed reservoir inflow into an existing rigid phase-equilibrium tank.
//!
//! Tank storage is original U; the supplied stream carries original H and its
//! species/entropy once. The final composition is rebound through the existing
//! source catalogs, then the native U,V owner resolves the complete endpoint.
//! This supplies no mass-flow law, empty-tank model or resolved filling path.

use super::super::{
    ThermochemicalFlashError, ThermochemicalPengRobinsonModel,
    equilibrium::{
        ResolvedThermochemicalPhaseState,
        mixing::{PhaseMixingFeed, ideal_composition_entropy},
        volume::PhaseEnergyVolumeControls,
    },
    failure,
};
use super::{IsochoricMixtureEquilibriumParcel, admit_inventory, source_specific_volume};
use crate::thermochemical_gas::energy::{
    GasEnergyError,
    parcel::transfer::{inventory, mixed_model_from_portions},
};
use crate::thermochemical_gas::nonideal::energy::parcel::transfer::mixing::model_for_source;
use fs_math::eft::{two_prod, two_sum};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// A prescribed finite incoming portion from one held equilibrium source state.
#[derive(Clone, Copy, Debug)]
pub struct MixtureEquilibriumInflow<'m> {
    /// Original species calorics, overall composition and PR interactions.
    pub model: &'m ThermochemicalPengRobinsonModel,
    /// Positive source temperature, held over the incoming portion.
    pub temperature: Temperature,
    /// Positive source pressure, held over the incoming portion.
    pub pressure: Pressure,
    /// Prescribed nonnegative mass [kg]. Zero is an exact no-op.
    pub mass_kg: f64,
}

/// Original U,V policy and independent physical mass/entropy allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureEquilibriumInflowControl {
    /// Original U,V domains and residual tolerances. Its PT limit covers the
    /// incoming selection, mixed reference and every nested U,V trial together.
    pub equilibrium: PhaseEnergyVolumeControls,
    /// Optional common source-supported caloric reference; None uses the
    /// accepted tank temperature. The reference pressure is the tank pressure.
    /// This changes numerical coordinates, not the physical target.
    pub reference_temperature: Option<Temperature>,
    /// Absolute total and individual-species mass discrepancy allowance [kg].
    pub absolute_mass_tolerance_kg: f64,
    /// Relative mass allowance, scaled separately for each species and total.
    pub relative_mass_tolerance: f64,
    /// Numerical allowance for negative adiabatic entropy production [J/K].
    pub entropy_tolerance_j_k: f64,
}
impl MixtureEquilibriumInflowControl {
    /// Retain the native U,V controls and use the current tank as the reference.
    #[must_use]
    pub fn new(equilibrium: PhaseEnergyVolumeControls) -> Self {
        Self {
            equilibrium,
            reference_temperature: None,
            absolute_mass_tolerance_kg: 1e-12,
            relative_mass_tolerance: 1e-12,
            entropy_tolerance_j_k: 1e-7,
        }
    }
    /// Admit controls without changing the inventory.
    pub fn validate(self) -> Result<(), MixtureEquilibriumInflowError> {
        self.equilibrium
            .validate()
            .map_err(|e| phase_error(2, 0, e))?;
        if self
            .reference_temperature
            .is_some_and(|t| !positive(t.value()))
            || !positive(self.absolute_mass_tolerance_kg)
            || !self.relative_mass_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_mass_tolerance)
            || !self.entropy_tolerance_j_k.is_finite()
            || self.entropy_tolerance_j_k < 0.0
        {
            return Err(MixtureEquilibriumInflowError::Invalid(
                "equilibrium inflow controls",
            ));
        }
        Ok(())
    }
    fn mass_tolerance(self, scale: f64) -> f64 {
        self.absolute_mass_tolerance_kg + self.relative_mass_tolerance * scale
    }
}

/// No partial composition, mass, energy or phase state is published on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureEquilibriumInflowError {
    /// Invalid input or unrepresentable positive inventory.
    Invalid(&'static str),
    /// Original source/card/species/PR catalog or inventory admission refused.
    Catalog(GasEnergyError),
    /// Native equilibrium refusal. Stages0/1/2 are the inlet, mixed reference
    /// and U,V solve;3/4 are the old/incoming ideal calorics. Completed PT
    /// work precedes the failing nested query.
    Equilibrium {
        stage: usize,
        completed_pt_evaluations: usize,
        source: ThermochemicalFlashError,
    },
    /// Shared PT work exhausted before another source selection.
    WorkLimit { pt_evaluations: usize },
    /// Source pressure is below an initial or accepted endpoint tank pressure.
    SupplyPressure {
        supply_pressure_pa: f64,
        tank_pressure_pa: f64,
    },
    /// Actual total/species inventory failed its separate physical allowance.
    MassBalance {
        species: Option<SpeciesId>,
        residual_kg: f64,
        tolerance_kg: f64,
    },
    /// Original U balance or retained physical vessel geometry failed closure.
    Balance {
        energy_residual_j: f64,
        volume_residual_m3: f64,
    },
    /// Original tank entropy minus transported inlet entropy is too negative.
    Entropy {
        production_j_k: f64,
        tolerance_j_k: f64,
    },
    /// Cancellation before atomic endpoint publication.
    Cancelled,
}
impl core::fmt::Display for MixtureEquilibriumInflowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid equilibrium tank inflow: {self:?}")
    }
}
impl std::error::Error for MixtureEquilibriumInflowError {}
type Error = MixtureEquilibriumInflowError;

/// Original species inventory, in canonical source order.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureEquilibriumSpeciesInflow {
    /// Conserved original species.
    pub species: SpeciesId,
    /// Initial tank species mass [kg].
    pub initial_mass_kg: f64,
    /// Prescribed supplied species mass [kg].
    pub incoming_mass_kg: f64,
    /// Actual accepted tank species mass [kg].
    pub final_mass_kg: f64,
    /// Accepted minus initial minus incoming mass [kg].
    pub residual_kg: f64,
}

/// Accepted original-source intake and its reference-safe physical accounts.
///
/// Absolute transported/stored energy and entropy reports retain original
/// species datums and can round away a small difference. Use the independent
/// residual/production fields below, not subtraction of the absolute reports.
/// Species-reference cancellation follows the intended conserved species
/// ledger; actual floating-point species amounts have their separate audit.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureEquilibriumInflowStep {
    /// Prescribed incoming mass [kg].
    pub incoming_mass_kg: f64,
    /// Actual final-minus-initial total mass minus incoming mass [kg].
    pub mass_residual_kg: f64,
    /// Canonical original-species audit; empty for the exact zero-portion path.
    pub species: Vec<MixtureEquilibriumSpeciesInflow>,
    /// Retained original inlet selection; absent for zero incoming mass.
    pub inlet_state: Option<ResolvedThermochemicalPhaseState>,
    /// Transported original flow enthalpy, dm*h_in [J].
    pub incoming_enthalpy_j: f64,
    /// Included inlet pressure-flow work dm*p_in*v_in [J]. This is already
    /// inside incoming_enthalpy_j and is never added a second time.
    pub incoming_pressure_flow_work_j: f64,
    /// Original species-reference and source-integrated tank U increment [J].
    pub stored_internal_energy_change_j: f64,
    /// Datum-free stored U increment minus incoming H [J].
    pub energy_residual_j: f64,
    /// Source-predicted final volume minus the unchanged vessel volume [m3].
    pub volume_residual_m3: f64,
    /// Transported original entropy, dm*s_in [J/K].
    pub incoming_entropy_j_k: f64,
    /// Original species-reference and source-integrated tank entropy change [J/K].
    pub stored_entropy_change_j_k: f64,
    /// Original tank entropy change minus transported inlet entropy [J/K].
    pub entropy_production_j_k: f64,
    /// Ideal overall-composition mixing contribution, included once above [J/K].
    pub ideal_mixing_entropy_j_k: f64,
    /// Actual native pressure/enthalpy inversions in the U,V closure.
    pub pressure_evaluations: usize,
    /// Whole-operation PT selections, including inlet and mixed reference.
    pub pt_evaluations: usize,
    /// Native EOS work ceilings, not a measured EOS call count.
    pub eos_evaluation_upper_bound: usize,
}

impl IsochoricMixtureEquilibriumParcel {
    /// Receive a prescribed equilibrium stream portion into this fixed volume.
    ///
    /// The supplied source holds its declared T,p,overall composition and h over
    /// the supplied portion; it is an immutable reservoir boundary, not another
    /// finite vessel whose mass or energy is silently changed. The existing tank
    /// remains nonempty. Its exact physical volume is retained, while source
    /// species are joined and the native U,V owner solves
    /// `U_final = U_initial + dm*h_in`, `v_final = V/(m_initial+dm)`.
    ///
    /// No heat, shaft work, reaction, moving wall, selective phase entrainment,
    /// resolved jet kinetic energy, time or flow rate is inferred. Complete bulk
    /// source enthalpy includes pressure-flow work once. Original source/phase
    /// domains and the existing accepted phase policy remain authoritative.
    ///
    /// Supply pressure must be at least the initial AND accepted final tank
    /// pressure. These are necessary endpoint checks, not proof that every
    /// intermediate state permits passive filling. A flow/path owner must admit
    /// that process separately. The caller may subdivide the supplied portion.
    ///
    /// All candidates are private. Source, budget, pressure, mass, energy,
    /// entropy or final-cancellation refusal leaves the complete tank unchanged.
    /// Zero mass preserves every stored bit and does no PT/source selection.
    pub fn receive_stream(
        &mut self,
        incoming: MixtureEquilibriumInflow<'_>,
        interaction_catalog: Option<&ThermochemicalPengRobinsonModel>,
        controls: MixtureEquilibriumInflowControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MixtureEquilibriumInflowStep, Error> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !incoming.mass_kg.is_finite()
            || incoming.mass_kg < 0.0
            || !positive(incoming.temperature.value())
            || !positive(incoming.pressure.value())
        {
            return Err(Error::Invalid(
                "finite nonnegative incoming mass and positive source T,p required",
            ));
        }
        let initial_volume_residual =
            self.mass_kg * (source_specific_volume(&self.state) - self.specific_volume_m3_kg);
        let initial_volume_tolerance =
            self.mass_kg * controls.equilibrium.specific_volume_tolerance_m3_kg;
        if !positive(initial_volume_tolerance)
            || !initial_volume_residual.is_finite()
            || initial_volume_residual.abs() > initial_volume_tolerance
        {
            return Err(Error::Balance {
                energy_residual_j: 0.0,
                volume_residual_m3: initial_volume_residual,
            });
        }
        if incoming.mass_kg == 0.0 {
            poll(&mut cancelled)?;
            return Ok(MixtureEquilibriumInflowStep {
                incoming_mass_kg: 0.0,
                mass_residual_kg: 0.0,
                species: Vec::new(),
                inlet_state: None,
                incoming_enthalpy_j: 0.0,
                incoming_pressure_flow_work_j: 0.0,
                stored_internal_energy_change_j: 0.0,
                energy_residual_j: 0.0,
                volume_residual_m3: initial_volume_residual,
                incoming_entropy_j_k: 0.0,
                stored_entropy_change_j_k: 0.0,
                entropy_production_j_k: 0.0,
                ideal_mixing_entropy_j_k: 0.0,
                pressure_evaluations: 0,
                pt_evaluations: 0,
                eos_evaluation_upper_bound: 0,
            });
        }
        if controls.equilibrium.phase.phase_options() != self.state.phase().options() {
            return Err(Error::Invalid(
                "continuing inflow requires the accepted original phase policy",
            ));
        }
        if incoming.pressure.value() < self.state.pressure_pa() {
            return Err(Error::SupplyPressure {
                supply_pressure_pa: incoming.pressure.value(),
                tank_pressure_pa: self.state.pressure_pa(),
            });
        }
        let final_mass = self.mass_kg + incoming.mass_kg;
        if !positive(final_mass) || final_mass <= self.mass_kg {
            return Err(Error::Invalid(
                "positive incoming mass is lost or overflows the tank inventory",
            ));
        }
        let mass_residual = sum([final_mass, -self.mass_kg, -incoming.mass_kg]);
        admit_mass(None, mass_residual, controls.mass_tolerance(final_mass))?;
        let final_specific_volume =
            admit_inventory(final_mass, self.volume_m3).map_err(|e| phase_error(2, 0, e))?;
        let energy_tolerance = final_mass * controls.equilibrium.internal_energy_tolerance_j_kg;
        let volume_tolerance = final_mass * controls.equilibrium.specific_volume_tolerance_m3_kg;
        if !positive(energy_tolerance) || !positive(volume_tolerance) {
            return Err(Error::Invalid("unrepresentable extensive U,V allowances"));
        }
        let joined = if self.model.source() == incoming.model.source() {
            self.model.source().clone()
        } else {
            mixed_model_from_portions(
                &[
                    (self.model.source(), self.mass_kg),
                    (incoming.model.source(), incoming.mass_kg),
                ],
                final_mass,
            )
            .map_err(Error::Catalog)?
        };
        let mut catalogs = vec![&self.model, incoming.model];
        if let Some(catalog) = interaction_catalog {
            catalogs.push(catalog);
        }
        let model = model_for_source(&catalogs, joined).map_err(Error::Catalog)?;
        let mut initial_species =
            inventory(self.model.source(), self.mass_kg).map_err(Error::Catalog)?;
        let mut incoming_species =
            inventory(incoming.model.source(), incoming.mass_kg).map_err(Error::Catalog)?;
        let final_species = inventory(model.source(), final_mass).map_err(Error::Catalog)?;
        let mut species = Vec::with_capacity(final_species.len());
        for (id, after) in final_species {
            let before = initial_species.remove(&id).unwrap_or(0.0);
            let received = incoming_species.remove(&id).unwrap_or(0.0);
            let residual = sum([after, -before, -received]);
            admit_mass(
                Some(id.clone()),
                residual,
                controls.mass_tolerance(after.max(before + received)),
            )?;
            species.push(MixtureEquilibriumSpeciesInflow {
                species: id,
                initial_mass_kg: before,
                incoming_mass_kg: received,
                final_mass_kg: after,
                residual_kg: residual,
            });
        }
        if !initial_species.is_empty() || !incoming_species.is_empty() {
            return Err(Error::Invalid("joined source lost a conserved species"));
        }
        poll(&mut cancelled)?;
        let inlet_state = incoming
            .model
            .evaluate_phase_equilibrium(
                incoming.temperature,
                incoming.pressure,
                controls.equilibrium.phase.phase_options(),
                &mut cancelled,
            )
            .map_err(|e| phase_error(0, 0, e))?;
        let mut pt = 1;
        let mut eos = inlet_state.caloric_sample().eos_evaluation_upper_bound;
        let remaining = controls.equilibrium.maximum_pt_evaluations - pt;
        if remaining == 0 {
            return Err(Error::WorkLimit { pt_evaluations: pt });
        }
        let reference_temperature = controls
            .reference_temperature
            .unwrap_or(Temperature::new(self.state.temperature_k()));
        let reference_pressure = Pressure::new(self.state.pressure_pa());
        let mut native_controls = controls.equilibrium;
        native_controls.maximum_pt_evaluations =
            native_controls.maximum_pt_evaluations.min(remaining);
        let path = model
            .prepare_phase_energy_volume_path(
                reference_temperature,
                reference_pressure,
                native_controls,
                &mut cancelled,
            )
            .map_err(|e| phase_error(1, pt, e))?;
        pt += path.preparation_pt_evaluations();
        eos = eos
            .checked_add(path.preparation_eos_evaluation_upper_bound())
            .ok_or(Error::Invalid("inflow source-work counter overflow"))?;
        let anchor = path.reference_sample();
        let anchor_mass = model.mechanics().molar_mass().value();
        let anchor_u = anchor.residual_u / anchor_mass;
        let anchor_s = anchor.entropy_correction / anchor_mass;
        let old = self.state.caloric_sample();
        let old_mass = self.model.mechanics().molar_mass().value();
        let old_u = old.residual_u / old_mass;
        let old_s = old.entropy_correction / old_mass;
        let feed = inlet_state.caloric_sample();
        let feed_mass = incoming.model.mechanics().molar_mass().value();
        let feed_h = feed.residual_h / feed_mass;
        let feed_s = feed.entropy_correction / feed_mass;
        let sources = [
            self.model.source().thermodynamics(),
            incoming.model.source().thermodynamics(),
        ];
        let temperatures = [
            Temperature::new(self.state.temperature_k()),
            incoming.temperature,
        ];
        let pressures = [reference_pressure, incoming.pressure];
        let mut ideal_u = [0.0; 2];
        let mut ideal_s = [0.0; 2];
        for i in 0..2 {
            poll(&mut cancelled)?;
            ideal_u[i] = sources[i]
                .specific_internal_energy_change(reference_temperature, temperatures[i])
                .map_err(|e| {
                    phase_error(
                        3 + i,
                        pt,
                        failure("inflow original ideal U increment", e).into(),
                    )
                })?
                .value();
            ideal_s[i] = sources[i]
                .specific_entropy_change(
                    reference_temperature,
                    reference_pressure,
                    temperatures[i],
                    pressures[i],
                )
                .map_err(|e| {
                    phase_error(
                        3 + i,
                        pt,
                        failure("inflow original ideal entropy increment", e).into(),
                    )
                })?
                .value();
        }
        let old_target = sum([ideal_u[0], old_u, -anchor_u]);
        // h_ideal(Tin) = u_ideal(Tref) + Delta u_ideal + R*Tin/M.
        // The ORIGINAL H departure completes the incoming flow enthalpy.
        let feed_target = sum([
            ideal_u[1],
            R * incoming.temperature.value() / feed_mass,
            feed_h,
            -anchor_u,
        ]);
        let target_energy = products([(self.mass_kg, old_target), (incoming.mass_kg, feed_target)]);
        let target_increment = target_energy / final_mass;
        if !target_energy.is_finite()
            || !target_increment.is_finite()
            || target_energy != 0.0 && target_increment == 0.0
        {
            return Err(Error::Invalid(
                "unrepresentable reference-free inflow energy target",
            ));
        }
        let weights = [self.mass_kg / final_mass, incoming.mass_kg / final_mass];
        if weights.iter().any(|w| !positive(*w)) {
            return Err(Error::Invalid(
                "positive species mixture weight is unrepresentable",
            ));
        }
        let mixing = ideal_composition_entropy(
            &[
                PhaseMixingFeed {
                    model: &self.model,
                    temperature: temperatures[0],
                    weight: self.mass_kg,
                },
                PhaseMixingFeed {
                    model: incoming.model,
                    temperature: incoming.temperature,
                    weight: incoming.mass_kg,
                },
            ],
            weights,
            &model,
        );
        let solved = path
            .solve_increment(
                final_specific_volume,
                target_increment,
                controls.equilibrium.maximum_pressure_evaluations,
                controls.equilibrium.maximum_pt_evaluations - pt,
                &mut cancelled,
            )
            .map_err(|e| phase_error(2, pt, e))?;
        pt += solved.pt_evaluations();
        eos = eos
            .checked_add(solved.eos_evaluation_upper_bound())
            .ok_or(Error::Invalid("inflow source-work counter overflow"))?;
        let pressure_evaluations = solved.pressure_evaluations();
        let increment = solved
            .accepted_increment_j_kg()
            .ok_or(Error::Invalid("missing original U,V energy increment"))?;
        let entropy_increment = solved
            .entropy_change_j_kg_k()
            .ok_or(Error::Invalid("missing original U,V entropy increment"))?;
        let energy_residual = products([
            (final_mass, increment),
            (-self.mass_kg, old_target),
            (-incoming.mass_kg, feed_target),
        ]);
        let volume_residual =
            final_mass * (source_specific_volume(solved.state()) - final_specific_volume);
        if !energy_residual.is_finite()
            || energy_residual.abs() > energy_tolerance
            || !volume_residual.is_finite()
            || volume_residual.abs() > volume_tolerance
        {
            return Err(Error::Balance {
                energy_residual_j: energy_residual,
                volume_residual_m3: volume_residual,
            });
        }
        if incoming.pressure.value() < solved.state().pressure_pa() {
            return Err(Error::SupplyPressure {
                supply_pressure_pa: incoming.pressure.value(),
                tank_pressure_pa: solved.state().pressure_pa(),
            });
        }
        let entropy = products([
            (final_mass, entropy_increment),
            (final_mass, mixing),
            (self.mass_kg, sum([anchor_s, -old_s, -ideal_s[0]])),
            (incoming.mass_kg, sum([anchor_s, -feed_s, -ideal_s[1]])),
        ]);
        if !entropy.is_finite() || entropy < -controls.entropy_tolerance_j_k {
            return Err(Error::Entropy {
                production_j_k: entropy,
                tolerance_j_k: controls.entropy_tolerance_j_k,
            });
        }
        // Absolute reports keep the original incoming species reference. The
        // independent balances above were evaluated BEFORE adding those datums.
        let reference = incoming
            .model
            .source()
            .thermodynamics()
            .evaluate_properties(reference_temperature, reference_pressure)
            .map_err(|e| {
                phase_error(4, pt, failure("inflow source reference receipt", e).into())
            })?;
        let incoming_enthalpy = incoming.mass_kg * inlet_state.specific_enthalpy_j_kg();
        let incoming_entropy = incoming.mass_kg * inlet_state.specific_entropy_j_kg_k();
        let incoming_flow_work =
            incoming.mass_kg * incoming.pressure.value() * source_specific_volume(&inlet_state);
        let stored_energy = products([
            (incoming.mass_kg, reference.u().value() / feed_mass),
            (final_mass, anchor_u),
            (-self.mass_kg, old_u),
            (final_mass, increment),
            (-self.mass_kg, ideal_u[0]),
        ]);
        let stored_entropy = products([
            (incoming.mass_kg, reference.s().value() / feed_mass),
            (final_mass, anchor_s),
            (-self.mass_kg, old_s),
            (-self.mass_kg, ideal_s[0]),
            (final_mass, entropy_increment),
            (final_mass, mixing),
        ]);
        let mixing_entropy = final_mass * mixing;
        if [
            incoming_enthalpy,
            incoming_entropy,
            incoming_flow_work,
            stored_energy,
            stored_entropy,
            mixing_entropy,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(Error::Invalid(
                "unrepresentable original extensive inflow account",
            ));
        }
        let candidate = Self::checked(model, solved.into_state(), final_mass, self.volume_m3)
            .map_err(|e| phase_error(2, pt, e))?;
        let report = MixtureEquilibriumInflowStep {
            incoming_mass_kg: incoming.mass_kg,
            mass_residual_kg: mass_residual,
            species,
            inlet_state: Some(inlet_state),
            incoming_enthalpy_j: incoming_enthalpy,
            incoming_pressure_flow_work_j: incoming_flow_work,
            stored_internal_energy_change_j: stored_energy,
            energy_residual_j: energy_residual,
            volume_residual_m3: volume_residual,
            incoming_entropy_j_k: incoming_entropy,
            stored_entropy_change_j_k: stored_entropy,
            entropy_production_j_k: entropy,
            ideal_mixing_entropy_j_k: mixing_entropy,
            pressure_evaluations,
            pt_evaluations: pt,
            eos_evaluation_upper_bound: eos,
        };
        poll(&mut cancelled)?;
        *self = candidate;
        Ok(report)
    }
}

fn admit_mass(species: Option<SpeciesId>, residual: f64, tolerance: f64) -> Result<(), Error> {
    if !positive(tolerance) || !residual.is_finite() || residual.abs() > tolerance {
        Err(Error::MassBalance {
            species,
            residual_kg: residual,
            tolerance_kg: tolerance,
        })
    } else {
        Ok(())
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut total = 0.0;
    let mut correction = 0.0;
    for value in values {
        let (next, error) = two_sum(total, value);
        total = next;
        correction += error;
    }
    total + correction
}
fn products(values: impl IntoIterator<Item = (f64, f64)>) -> f64 {
    sum(values.into_iter().flat_map(|(a, b)| {
        let (product, error) = two_prod(a, b);
        [product, error]
    }))
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn phase_error(
    stage: usize,
    completed_pt_evaluations: usize,
    source: ThermochemicalFlashError,
) -> Error {
    if source == ThermochemicalFlashError::Cancelled {
        Error::Cancelled
    } else {
        Error::Equilibrium {
            stage,
            completed_pt_evaluations,
            source,
        }
    }
}
