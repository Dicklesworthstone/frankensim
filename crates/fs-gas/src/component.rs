//! A closed, spatially uniform ideal-gas cavity on the composition contract.
//!
//! The original material parcel owns calorics, pressure, source admission and
//! the discrete Gibbs heat/work solve. This adapter exposes its work pressure
//! and entropy-conjugate temperature as ports; it supplies no second equation
//! of state, constant-capacity approximation or independent pressure quadrature.
//!
//! Pressure port `p` receives a compression flow [m3/s], positive when the
//! cavity shrinks, and returns absolute pressure [Pa]. Heat port `heat`
//! receives heat rate [W] into the gas and returns its thermodynamic mean
//! temperature [K]. Both are `FlowIn`. A piston-area transformer can join `p`
//! to an inertia's effort-input port without adding piston state to the gas.
//!
//! [`GasCavity::insulated`] is the same cavity without its heat port. Its
//! one relation, work pressure against compression flow, is scalar and
//! rising, and it says so with its class.
//!
//! [`GasCavity::streaming`] adds a stream port: gas of the reference's
//! species may enter and leave. Through it the cavity receives the mass
//! rate, the enthalpy rate in the species' own datum, and the mass rate of
//! each species, and returns its pressure, specific enthalpy and species
//! mass fractions as they were when the step began. The step is first what
//! crossed, at the old volume, by the parcel's boundary exchange, and then
//! the heat and wall work on what the cavity then holds. The stream side is
//! therefore explicit and first order in the step.
//!
//! [`GasCavity::publishing`] adds scalar signals for the pressure, temperature,
//! mass and volume at the accepted start of a step. They carry no energy and
//! introduce no new gas state; feedback through them is held over the step.
//!
//! This component has no reactions, nonideal equation of state or internal
//! wave field, and without `streaming` no open streams. It does not
//! know a connected piston's mass or area and therefore cannot enforce the
//! native piston owner's coupled acoustic-frequency timestep bound. Callers
//! must resolve their assembled dynamics; a closed ledger is not an accuracy
//! certificate or a claim to replace the complete native piston solver.

use core::num::NonZeroUsize;
use std::fmt::Write;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_material::thermochemical_gas::energy::parcel::discrete_gradient::GasDiscreteGradientCandidate;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::boundary::{
    GasSpeciesBoundarySolve, GasSpeciesBoundaryStep, SpeciesBoundaryTransfer,
};
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, STREAM_CARRIED,
    STREAM_ENTHALPY_RATE, STREAM_PRESSURE, STREAM_SPECIFIC_ENTHALPY, StreamBasis, StreamQuantity,
};
use fs_ports::{PortKind, StableId};
use fs_qty::MolarMass;
use fs_thermochem::SpeciesId;

pub use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
pub use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};

fn invalid(message: &'static str) -> GasEnergyError {
    GasEnergyError::Invalid(message)
}

fn poll(cx: &Cx<'_>) -> Result<(), GasEnergyError> {
    cx.checkpoint().map_err(|_| GasEnergyError::Cancelled)
}

fn pair(values: &[f64]) -> Result<[f64; 2], GasEnergyError> {
    match values {
        [a, b] if a.is_finite() && b.is_finite() => Ok([*a, *b]),
        _ => Err(invalid("gas cavity requires two finite scalar port values")),
    }
}

/// Source-backed closed gas with conservative pressure-work and heat ports.
///
/// Its [`Component::State`] is the original [`ClosedGasParcel`]. The immutable
/// reference fixes only an energy datum and the admitted inventory/source
/// identity. Measured energy is the original caloric increment from that
/// reference, so large formation-energy constants cannot erase small work.
#[derive(Clone, Debug)]
pub struct GasCavity {
    reference: ClosedGasParcel,
    controls: GasEnergySolve,
    ports: Vec<PortSpec>,
    /// Set by [`Self::streaming`].
    stream: Option<StreamPlan>,
    /// Set by [`Self::insulated`]: there is no heat port, and no heat.
    insulated: bool,
    /// Set by [`Self::publishing`]: accepted-state observations follow the physical ports.
    published: bool,
}

/// What a streaming cavity needs besides its parcel: the species of its
/// reference, in the order its stream port carries them, and the largest
/// mass error its boundary exchange may leave.
#[derive(Clone, Debug)]
struct StreamPlan {
    species: Vec<SpeciesId>,
    mass_tolerance_kg: f64,
}

/// What a streaming cavity answers a tangent request with.
const NO_TANGENT: &str =
    "a streaming gas cavity offers no port tangent; add it to a system without one";

/// Internal energy [J] of `parcel` in its species' own datum.
fn absolute_energy_j(parcel: &ClosedGasParcel) -> f64 {
    parcel.mass_kg() * parcel.state().specific_internal_energy_j_kg()
}

/// Hex keeps arbitrary species/convention text unambiguous inside a StableId.
fn identity_bytes(identity: &mut String, bytes: &[u8]) {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    for &byte in bytes {
        identity.push(char::from(HEX[usize::from(byte >> 4)]));
        identity.push(char::from(HEX[usize::from(byte & 15)]));
    }
}

/// Bind the complete native source catalog, never live composition or state.
fn stream_basis(parcel: &ClosedGasParcel) -> StreamBasis {
    let mut catalog = Vec::new();
    parcel.model().catalog_identity_words(&mut catalog);
    let mut medium = String::from("native-gas-catalog-v1/");
    for word in catalog {
        identity_bytes(&mut medium, &word.to_le_bytes());
    }
    let mut datum = String::from("nasa9-caloric-reference-v1");
    let mut quantities = Vec::new();
    for model in parcel.model().thermodynamics().components() {
        quantities.push(StreamQuantity::SpeciesMass(model.species().clone()));
        datum.push('/');
        identity_bytes(&mut datum, model.species().as_str().as_bytes());
        if let Some(reference) = model.caloric_reference() {
            write!(
                &mut datum,
                "/continuous/{}/{:016x}/{}/{:016x}/{:016x}",
                reference.evaluator_version(),
                reference.anchor_temperature().value().to_bits(),
                reference.anchor_region(),
                reference.maximum_enthalpy_correction().value().to_bits(),
                reference.maximum_entropy_correction().value().to_bits(),
            )
            .expect("writing to a String");
            for region in reference.regions() {
                datum.push_str("/region/");
                identity_bytes(&mut datum, region.card_identity());
                write!(
                    &mut datum,
                    "/{:016x}/{:016x}",
                    region.enthalpy_offset().value().to_bits(),
                    region.entropy_offset().value().to_bits(),
                )
                .expect("writing to a String");
            }
        } else {
            datum.push_str("/raw");
        }
    }
    StreamBasis::new(
        StableId::new(medium).expect("hex catalog identity is a stable id"),
        StableId::new(datum).expect("hex caloric identity is a stable id"),
        quantities,
    )
    .expect("the admitted mixture has distinct species")
}

impl StreamPlan {
    /// The mass fraction of each species of the plan in `parcel`; zero for
    /// one it holds none of.
    fn fractions(&self, parcel: &ClosedGasParcel) -> Result<Vec<f64>, GasEnergyError> {
        let mixture = parcel.model().thermodynamics();
        let molar_masses: Vec<MolarMass> = mixture
            .components()
            .iter()
            .map(fs_thermochem::Nasa9StandardStateModelV1::molar_mass)
            .collect();
        let by_mass = mixture
            .mole_composition()
            .to_mass_fractions(&molar_masses)
            .map_err(|_| invalid("gas cavity composition has no mass fractions"))?;
        Ok(self
            .species
            .iter()
            .map(|wanted| {
                mixture
                    .components()
                    .iter()
                    .position(|held| held.species() == wanted)
                    .map_or(0.0, |index| by_mass.fractions()[index])
            })
            .collect())
    }
}

impl GasCavity {
    /// Bind the existing parcel and its original caloric solve controls.
    /// Pass a clone of `reference` as the initial state when adding the component.
    ///
    /// # Errors
    /// Original source, temperature-branch, derivative or control refusal at
    /// the unchanged state. No material parameters or source ranges are inferred.
    pub fn new(reference: ClosedGasParcel, controls: GasEnergySolve) -> Result<Self, GasEnergyError> {
        reference.discrete_gradient_heat_volume(reference.volume_m3(), 0.0, controls, || false)?;
        Ok(Self {
            reference,
            controls,
            ports: vec![
                PortSpec::power(
                    StableId::new("p").expect("static port name"),
                    PortKind::FluidPressureFlux,
                    PortRole::FlowIn,
                    NonZeroUsize::MIN,
                )
                .expect("pressure-volume power port"),
                PortSpec::heat(
                    StableId::new("heat").expect("static port name"),
                    PortRole::FlowIn,
                )
                .expect("heat input port"),
            ],
            stream: None,
            insulated: false,
            published: false,
        })
    }

    /// Publish accepted pressure [Pa], temperature [K], mass [kg] and volume
    /// [m3] as scalar signals `pressure`, `temperature`, `mass` and `volume`.
    ///
    /// These observations are held at the start of each proposed step, not
    /// the endpoint or the work-conjugate mean returned by a physical port.
    /// They carry no energy. Their input derivatives are zero at fixed accepted
    /// state; control accuracy and stability remain properties of the assembled
    /// system and its step size. Publishing is idempotent and may precede or
    /// follow [`Self::insulated`] or [`Self::streaming`]. Physical ports always
    /// precede these signals, and publishing changes no native state codec.
    #[must_use]
    pub fn publishing(mut self) -> Self {
        if !self.published {
            for name in ["pressure", "temperature", "mass", "volume"] {
                self.ports.push(PortSpec::signal_out(
                    StableId::new(name).expect("static gas observation name"),
                    NonZeroUsize::MIN,
                ));
            }
            self.published = true;
        }
        self
    }

    /// The same cavity behind an insulating wall: its one physical port is `p`, and
    /// a step is the two-port cavity's with a heat rate of exactly zero.
    ///
    /// The work pressure it returns rises with the compression flow it
    /// receives, so without published signals its class is [`Class::Monotone1d`].
    /// With a heat port or published observations the map is vector-valued
    /// and stays [`Class::Smooth`].
    ///
    /// Two things follow from the class. A kernel that differences only
    /// smooth components (fs-system today) needs this cavity's tangent: in
    /// a loop, add it with the tangent it offers, or the step is refused as
    /// unsupported. And "rises" holds down to the rounding of the owner's
    /// energy solve, about 1e-15 of the pressure: between two flows closer
    /// than that resolves, the returned pressure can fall by as much, so a
    /// solve that brackets it must stop above that.
    ///
    /// # Errors
    /// A streaming cavity, which keeps its heat port.
    pub fn insulated(mut self) -> Result<Self, GasEnergyError> {
        if self.stream.is_some() {
            return Err(invalid("a streaming gas cavity keeps its heat port"));
        }
        if !self.insulated {
            self.ports.remove(1);
            self.insulated = true;
        }
        Ok(self)
    }

    /// The same cavity with a stream port `stream` after `p` and `heat`:
    /// gas of the reference's species may enter and leave.
    ///
    /// The port carries one quantity per species of the reference, in its
    /// order. It receives mass rate [kg/s], enthalpy rate [W] in the
    /// species' own datum and each species' mass rate [kg/s], and returns
    /// the pressure [Pa], specific enthalpy [J/kg] and species mass
    /// fractions the cavity had when the step began. The species rates are
    /// what changes the inventory. The component exposes the native total
    /// and species masses for the kernel's independent endpoint audits, so
    /// inconsistent total/species transport cannot publish an accepted step.
    /// Proposals themselves allow independent port probes for the interface
    /// Jacobian. A species increment too small to change what is held is
    /// dropped within the shared audit's rounding allowance. A step that
    /// would draw more of a species than is held, or all of the gas, is
    /// shortened.
    ///
    /// The stream declares its ordered species, source cards, molar masses
    /// and caloric conventions. Connections refuse a different catalog or
    /// energy reference, even at equal width. Composition and live state do
    /// not enter the identity, so compatible mixtures may exchange gas.
    /// A transport element must retain this basis on both ports.
    ///
    /// The stored energy of a streaming cavity is its internal energy in
    /// that datum, not the increment from the reference, so a large
    /// formation energy costs it rounding the closed cavity does not pay.
    /// It offers no port tangent: add it to a system without one.
    ///
    /// # Errors
    /// A mass tolerance that is not positive and finite, or an insulated
    /// cavity.
    pub fn streaming(mut self, mass_tolerance_kg: f64) -> Result<Self, GasEnergyError> {
        if !(mass_tolerance_kg.is_finite() && mass_tolerance_kg > 0.0) {
            return Err(invalid("a streaming gas cavity requires a positive mass tolerance"));
        }
        if self.insulated {
            return Err(invalid("an insulated gas cavity takes no stream port"));
        }
        let species: Vec<SpeciesId> = self
            .reference
            .model()
            .thermodynamics()
            .components()
            .iter()
            .map(|species| species.species().clone())
            .collect();
        if self.stream.is_none() {
            self.ports.insert(
                2,
                PortSpec::identified_stream(
                    StableId::new("stream").expect("static port name"),
                    PortRole::FlowIn,
                    stream_basis(&self.reference),
                )
                .expect("stream storage port"),
            );
        }
        self.stream = Some(StreamPlan {
            species,
            mass_tolerance_kg,
        });
        Ok(self)
    }

    /// Original parcel defining the source/inventory identity and energy datum.
    #[must_use]
    pub const fn reference(&self) -> &ClosedGasParcel {
        &self.reference
    }

    /// Original source-temperature bracket, absolute energy budget and work cap.
    #[must_use]
    pub const fn controls(&self) -> GasEnergySolve {
        self.controls
    }

    fn output_width(&self) -> usize {
        let physical = self.stream.as_ref().map_or_else(
            || if self.insulated { 1 } else { 2 },
            |stream| 2 + STREAM_CARRIED + stream.species.len(),
        );
        physical + if self.published { 4 } else { 0 }
    }

    /// Called only after admission and the final cancellation poll. Every
    /// observation is a native accepted-state value, independent of this trial.
    fn publish(&self, old: &ClosedGasParcel, outputs: &mut [f64]) {
        if self.published {
            let start = outputs.len() - 4;
            outputs[start..].copy_from_slice(&[
                old.state().state().pressure,
                old.state().state().temperature,
                old.mass_kg(),
                old.volume_m3(),
            ]);
        }
    }

    /// Depletion may remove species; it may not replace a retained source.
    fn admit_stream_state(&self, state: &ClosedGasParcel) -> Result<(), GasEnergyError> {
        self.reference.admit_inventory_state(state)
    }

    /// Compression flow and heat rate of a step: as received, or for an
    /// insulated cavity its one input and no heat.
    fn held(&self, inputs: &[f64]) -> Result<[f64; 2], GasEnergyError> {
        if !self.insulated {
            return pair(inputs);
        }
        match inputs {
            [flow] if flow.is_finite() => Ok([*flow, 0.0]),
            _ => Err(invalid("an insulated gas cavity requires one finite port value")),
        }
    }

    fn trial(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(GasDiscreteGradientCandidate, f64), GasEnergyError> {
        poll(cx)?;
        let [compression, heat_rate] = self.held(inputs)?;
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(invalid("gas cavity requires a positive finite timestep"));
        }
        // This also checks unchanged mass, composition, complete source model
        // and auxiliary query coordinates before the native proposal.
        let stored_old = self.reference.stored_energy_change_to(old, || cx.is_cancel_requested())?;
        let displaced = dt_s * compression;
        let heat = dt_s * heat_rate;
        let volume = old.volume_m3() - displaced;
        if !(volume.is_finite() && volume > 0.0)
            || !heat.is_finite()
            || (compression != 0.0 && displaced == 0.0)
            || (heat_rate != 0.0 && heat == 0.0)
        {
            return Err(invalid("gas cavity volume or heat increment is unrepresentable"));
        }
        let trial = old.discrete_gradient_heat_volume(volume, heat, self.controls, || {
            cx.is_cancel_requested()
        })?;
        Ok((trial, stored_old))
    }

    /// What `old` holds once `entered_j` of enthalpy and the species masses
    /// `dt_s * rates` have crossed at its volume, and the owner's account
    /// of that; `old` itself when nothing crossed. `Ok(None)` when the step
    /// would draw more of a species than the cavity holds, or all of its
    /// gas: a shorter step is wanted.
    ///
    /// `fractions` are `old`'s species mass fractions, formed as the owner
    /// forms its inventory, so the two tests below are the owner's own.
    #[allow(clippy::too_many_arguments)]
    fn crossed(
        &self,
        stream: &StreamPlan,
        old: &ClosedGasParcel,
        dt_s: f64,
        entered_j: f64,
        rates: &[f64],
        fractions: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Option<(ClosedGasParcel, Option<GasSpeciesBoundaryStep>)>, GasEnergyError> {
        let mut rows = Vec::with_capacity(rates.len());
        let mut entering_kg = 0.0;
        for (index, ((species, rate), fraction)) in
            stream.species.iter().zip(rates).zip(fractions).enumerate()
        {
            let held_kg = old.mass_kg() * fraction;
            let mut mass_kg = dt_s * rate;
            // An increment too small to change the inventory is none. The
            // owner refuses one, and a flow that dies away would otherwise
            // end the run when it reached the rounding of what is held.
            #[allow(clippy::float_cmp)]
            if held_kg + mass_kg == held_kg {
                mass_kg = 0.0;
            }
            if held_kg + mass_kg < 0.0 {
                return Ok(None);
            }
            entering_kg += mass_kg;
            // The owner sums the rows' energies, so the whole of the
            // stream's enthalpy goes on the first one.
            rows.push(SpeciesBoundaryTransfer {
                species: species.clone(),
                mass_into_gas_kg: mass_kg,
                energy_into_gas_j: if index == 0 { entered_j } else { 0.0 },
            });
        }
        if entered_j == 0.0 && rows.iter().all(|row| row.mass_into_gas_kg == 0.0) {
            return Ok(Some((old.clone(), None)));
        }
        if !(old.mass_kg() + entering_kg > 0.0) {
            return Ok(None);
        }
        // The owner's inversion is gated per kilogram of what the cavity
        // holds afterwards.
        let solve = GasSpeciesBoundarySolve {
            gas: self.controls,
            mass_tolerance_kg: stream.mass_tolerance_kg,
            energy_tolerance_j: old.mass_kg().max(old.mass_kg() + entering_kg)
                * self.controls.energy_tolerance_j_kg,
        };
        let candidate = old.exchange_species_boundary(self.reference.model(), &rows, solve, || {
            cx.is_cancel_requested()
        })?;
        let (gas, _, step) = candidate.into_parts();
        Ok(Some((gas, Some(step))))
    }

    /// One step of a streaming cavity: what crossed, at the old volume,
    /// then heat and wall work on what the cavity then holds.
    #[allow(clippy::too_many_arguments)]
    fn propose_streaming(
        &self,
        stream: &StreamPlan,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        new: &mut ClosedGasParcel,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, GasEnergyError> {
        poll(cx)?;
        let width = 2 + STREAM_CARRIED + stream.species.len();
        if inputs.len() != width
            || outputs.len() != self.output_width()
            || inputs.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "a streaming gas cavity requires finite values on its three ports",
            ));
        }
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(invalid("gas cavity requires a positive finite timestep"));
        }
        let (compression, heat_rate) = (inputs[0], inputs[1]);
        self.admit_stream_state(old)?;
        let crossing = &inputs[2..];
        let rates = &crossing[STREAM_CARRIED..];
        // What the cavity holds is like when the step begins.
        let pressure_old = old.state().state().pressure;
        let enthalpy_old = old.state().specific_enthalpy_j_kg();
        let fractions = stream.fractions(old)?;
        let stored_old = absolute_energy_j(old);

        // What crossed.
        let entered_j = dt_s * crossing[STREAM_ENTHALPY_RATE];
        let Some((held, exchange)) =
            self.crossed(stream, old, dt_s, entered_j, rates, &fractions, cx)?
        else {
            return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
        };
        let budget_j = old.mass_kg().max(held.mass_kg()) * self.controls.energy_tolerance_j_kg;

        // Heat and wall work on what it then holds.
        let displaced = dt_s * compression;
        let heat = dt_s * heat_rate;
        let volume = held.volume_m3() - displaced;
        if !(volume.is_finite() && volume > 0.0 && heat.is_finite())
            || (compression != 0.0 && displaced == 0.0)
            || (heat_rate != 0.0 && heat == 0.0)
        {
            return Err(invalid("gas cavity volume or heat increment is unrepresentable"));
        }
        let trial = held.discrete_gradient_heat_volume(volume, heat, self.controls, || {
            cx.is_cancel_requested()
        })?;
        let step = trial.step();
        let pressure = step.work_pressure_pa;
        let theta = trial.thermodynamic_mean().state().temperature;
        let stored_new = absolute_energy_j(trial.endpoint());
        let (crossed_change, crossed_residual) = exchange.as_ref().map_or((0.0, 0.0), |crossed| {
            (crossed.stored_energy_change_j, crossed.energy_residual_j)
        });
        let rounding = 64.0 * f64::EPSILON;
        let allowance = budget_j
            + held.mass_kg() * self.controls.energy_tolerance_j_kg
            + rounding
                * (stored_old.abs()
                    + stored_new.abs()
                    + crossed_change.abs()
                    + entered_j.abs()
                    + step.stored_energy_change_j.abs()
                    + step.supplied_heat_j.abs()
                    + step.work_by_gas_j.abs()
                    + (pressure * held.volume_m3()).abs()
                    + (pressure * step.final_volume_m3).abs());
        let residual = step
            .energy_residual_j
            .abs()
            .max(step.gibbs_residual_j.abs())
            .max(theta * step.entropy_residual_j_k.abs())
            .max(crossed_residual.abs());
        if !(allowance.is_finite() && residual.is_finite() && stored_new.is_finite()) {
            return Err(invalid("gas cavity energy account is unrepresentable"));
        }
        let receipt = Receipt {
            stored_change_j: crossed_change + step.stored_energy_change_j,
            dissipated_j: 0.0,
            numerical_loss_j: -(crossed_residual + step.energy_residual_j),
            // Entropy carried in by the stream and made by mixing it is not
            // accounted (the stream port has no entropy on the wire).
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: allowance,
            residual,
            tolerance: (budget_j + held.mass_kg() * self.controls.energy_tolerance_j_kg)
                .next_up()
                .min(f64::MAX),
        };
        poll(cx)?;
        let (endpoint, _) = trial.into_parts();
        *new = endpoint;
        outputs[0] = pressure;
        outputs[1] = theta;
        outputs[2 + STREAM_PRESSURE] = pressure_old;
        outputs[2 + STREAM_SPECIFIC_ENTHALPY] = enthalpy_old;
        outputs[2 + STREAM_CARRIED..width].copy_from_slice(&fractions);
        self.publish(old, outputs);
        Ok(Proposal::Step(receipt))
    }

    fn jacobian(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<[[f64; 2]; 2], GasEnergyError> {
        let (trial, _) = self.trial(old, dt_s, inputs, cx)?;
        let step = trial.step();
        let theta = trial.thermodynamic_mean().state().temperature;
        let pressure = step.work_pressure_pa;
        let pressure_q = step.pressure_heat_derivative_pa_j;
        // The native ideal-gas residual is e(T1)+R ln(V1/V0) theta(T1)-Q/m.
        // With p=m R a theta, a=ln(V1/V0)/(V1-V0), implicit differentiation
        // gives theta_Q=(theta/p)p_Q and theta_V=-(m R theta/V1)theta_Q.
        // The owner's p_Q already includes its original caloric mean derivative,
        // including the coincident-state limit and source-region quadrature.
        // These identities are specific to ClosedGasParcel's ideal-gas law.
        let theta_q = (theta / pressure) * pressure_q;
        let theta_v = -(old.mass_kg() * old.state().state().specific_gas_constant * theta
            / step.final_volume_m3)
            * theta_q;
        let matrix = [
            [-dt_s * step.pressure_final_volume_derivative_pa_m3, dt_s * pressure_q],
            [-dt_s * theta_v, dt_s * theta_q],
        ];
        if matrix.iter().flatten().any(|entry| !entry.is_finite()) {
            return Err(invalid("gas cavity port tangent is unrepresentable"));
        }
        poll(cx)?;
        Ok(matrix)
    }

    /// Tangent of an insulated cavity, in either direction: the slope of
    /// its work pressure in its compression flow, the two-port cavity's
    /// entry for those two at no heat, times `seed`, into `slot`.
    fn insulated_slope(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        seed: &[f64],
        slot: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), GasEnergyError> {
        let ([seed], [slot]) = (seed, slot) else {
            return Err(invalid("an insulated gas cavity has one tangent value a side"));
        };
        let value = self.jacobian(old, dt_s, inputs, cx)?[0][0] * seed;
        if !value.is_finite() {
            return Err(invalid("gas cavity port tangent is unrepresentable"));
        }
        poll(cx)?;
        *slot = value;
        Ok(())
    }
}

impl Component for GasCavity {
    type State = ClosedGasParcel;
    type Error = GasEnergyError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        if self.insulated && !self.published {
            // One input, and the work pressure rises with the compression.
            return Class::Monotone1d;
        }
        // A heat port or published observations make the map vector-valued;
        // the scalar monotone pressure solve does not classify that map.
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &ClosedGasParcel) -> Result<f64, GasEnergyError> {
        if self.stream.is_some() {
            self.admit_stream_state(state)?;
            // The inventory changes, so there is no one reference to
            // measure an increment from.
            let energy = absolute_energy_j(state);
            return if energy.is_finite() {
                Ok(energy)
            } else {
                Err(invalid("gas cavity internal energy is unrepresentable"))
            };
        }
        self.reference.stored_energy_change_to(state, || false)
    }

    fn stored_mass_kg(&self, state: &ClosedGasParcel) -> Result<Option<f64>, GasEnergyError> {
        if self.stream.is_some() {
            self.admit_stream_state(state)?;
        } else {
            // The native endpoint check retains the closed parcel's mass,
            // complete source model and auxiliary material coordinates.
            self.reference.stored_energy_change_to(state, || false)?;
        }
        Ok(Some(state.mass_kg()))
    }

    fn stored_carried_quantity(
        &self,
        state: &ClosedGasParcel,
        carried: usize,
    ) -> Result<Option<f64>, GasEnergyError> {
        let Some(species) = self
            .stream
            .as_ref()
            .and_then(|stream| stream.species.get(carried))
        else {
            return Ok(None);
        };
        self.admit_stream_state(state)?;
        // The stream keeps its reference catalog order even after depletion
        // removes a species from the active mixture. The native inventory
        // owner measures that species as zero until it enters again.
        state.species_mass_kg(species).map(Some)
    }

    fn propose(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        new: &mut ClosedGasParcel,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, GasEnergyError> {
        if let Some(stream) = &self.stream {
            return self.propose_streaming(stream, old, dt_s, inputs, new, outputs, cx);
        }
        if outputs.len() != self.output_width() {
            return Err(invalid("gas cavity output extent must match its physical and signal ports"));
        }
        let (trial, stored_old) = self.trial(old, dt_s, inputs, cx)?;
        let step = trial.step();
        let pressure = step.work_pressure_pa;
        let theta = trial.thermodynamic_mean().state().temperature;
        let stored_new = self.reference.stored_energy_change_to(trial.endpoint(), || {
            cx.is_cancel_requested()
        })?;
        let budget = old.mass_kg() * self.controls.energy_tolerance_j_kg;
        // Include p*V: forming V1=V0-dt*f rounds the represented displacement,
        // while the kernel independently books the prescribed dt*f on the wire.
        // That representation defect stays visible in the ledger, not in heat.
        let rounding = 64.0 * f64::EPSILON;
        let allowance = budget
            + rounding * stored_old.abs()
            + rounding * stored_new.abs()
            + rounding * step.stored_energy_change_j.abs()
            + rounding * step.supplied_heat_j.abs()
            + rounding * step.work_by_gas_j.abs()
            + rounding * (pressure * old.volume_m3()).abs()
            + rounding * (pressure * step.final_volume_m3).abs();
        let residual = step.energy_residual_j.abs()
            .max(step.gibbs_residual_j.abs())
            .max(theta * step.entropy_residual_j_k.abs());
        if !(allowance.is_finite() && residual.is_finite()) {
            return Err(invalid("gas cavity energy account is unrepresentable"));
        }
        let receipt = Receipt {
            stored_change_j: step.stored_energy_change_j,
            dissipated_j: 0.0,
            numerical_loss_j: -step.energy_residual_j,
            // The gas-side heat entropy is Q/theta. The native Gibbs owner
            // separately bounds Delta S-Q/theta; a numerical entropy residual
            // is not physical entropy generation.
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: allowance,
            residual,
            // Converting the owner's entropy bound budget/theta back to
            // joules can round upward once, although that native gate passed.
            tolerance: budget.next_up().min(f64::MAX),
        };
        poll(cx)?;
        let (endpoint, _) = trial.into_parts();
        *new = endpoint;
        outputs[0] = pressure;
        // An insulated cavity has no heat port to return its mean
        // temperature on.
        if !self.insulated {
            outputs[1] = theta;
        }
        self.publish(old, outputs);
        Ok(Proposal::Step(receipt))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-gas::component::GasCavity",
        state_family: "gas-cavity-parcel-state",
        schema_grammar: "closed cavity: parcel volume, temperature and pressure on the cavity's reference; streaming cavity: version tag, mass, volume, temperature, pressure, fraction basis and N reference species fractions",
        codec_grammar: "closed: three u64 words (f64 bits of volume, temperature, pressure); streaming: 6+N u64 words (version tag, f64 bits of mass, volume, temperature, pressure, basis code, and N f64 fraction bits)",
        codec_version: 1,
    });

    /// A closed cavity retains its three-word codec. An open cavity delegates
    /// accepted mass and composition to the native inventory codec. Both
    /// refuse foreign source models or auxiliary axes without writing anything.
    fn encode_state(&self, state: &ClosedGasParcel, words: &mut Vec<u64>) -> bool {
        if self.stream.is_some() {
            return self.reference.encode_inventory_state(state, words).is_ok();
        }
        let on_reference = self.reference.stored_energy_change_to(state, || false);
        if on_reference.is_err() {
            return false;
        }
        state.state_words(words);
        true
    }

    /// The parcel the three words stand for: the reference itself if they
    /// are its own (a step with no heat and no motion keeps its parcel to
    /// the bit), and otherwise the reference at that volume and temperature
    /// ([`ClosedGasParcel::at_volume_temperature`]), which computes the
    /// pressure as a step does. Either is returned only if writing it again
    /// gives all three words back, so a pressure that is not the one the
    /// owner computes for the other two is refused, and only if the cavity
    /// can measure its energy from the reference.
    fn decode_state(&self, words: &[u64]) -> Option<ClosedGasParcel> {
        if self.stream.is_some() {
            return self.reference.decode_inventory_state(words).ok();
        }
        let &[volume, temperature, _] = words else {
            return None;
        };
        let mut rebuilt = Vec::with_capacity(words.len());
        self.reference.state_words(&mut rebuilt);
        let state = if rebuilt == words {
            self.reference.clone()
        } else {
            let [volume_m3, temperature_k] = [volume, temperature].map(f64::from_bits);
            self.reference
                .at_volume_temperature(volume_m3, temperature_k)
                .ok()?
        };
        rebuilt.clear();
        state.state_words(&mut rebuilt);
        let held = self.reference.stored_energy_change_to(&state, || false);
        (rebuilt == words && held.is_ok()).then_some(state)
    }

    /// The reference, solve controls and port mode bind durable state.
    /// Unpublished identities remain unchanged; published ones add a version tag.
    fn identity_words(&self, words: &mut Vec<u64>) {
        let Self {
            reference,
            controls,
            // Derived from `insulated`, `stream` and `published`.
            ports: _,
            stream,
            insulated,
            published,
        } = self;
        let GasEnergySolve {
            minimum_temperature_k,
            maximum_temperature_k,
            energy_tolerance_j_kg,
            maximum_iterations,
        } = controls;
        words.push(u64::from_le_bytes(if stream.is_some() {
            *b"gasopen1"
        } else {
            *b"gascavty"
        }));
        reference.identity_words(words);
        words.extend([minimum_temperature_k, maximum_temperature_k].map(|x| x.to_bits()));
        words.push(energy_tolerance_j_kg.to_bits());
        words.extend([*maximum_iterations as u64, u64::from(*insulated)]);
        if let Some(stream) = stream {
            words.push(stream.mass_tolerance_kg.to_bits());
        }
        if *published {
            words.push(u64::from_le_bytes(*b"gaspub01"));
        }
    }
}

impl PortTangent for GasCavity {
    fn jvp(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), GasEnergyError> {
        if self.stream.is_some() {
            return Err(invalid(NO_TANGENT));
        }
        if d_outputs.len() != self.output_width() {
            return Err(invalid("gas cavity tangent output extent must match its ports"));
        }
        if self.insulated {
            self.insulated_slope(old, dt_s, inputs, d_inputs, &mut d_outputs[..1], cx)?;
            d_outputs[1..].fill(0.0);
            return Ok(());
        }
        let direction = pair(d_inputs)?;
        let j = self.jacobian(old, dt_s, inputs, cx)?;
        let result = [
            j[0][0] * direction[0] + j[0][1] * direction[1],
            j[1][0] * direction[0] + j[1][1] * direction[1],
        ];
        pair(&result)?;
        poll(cx)?;
        d_outputs[..2].copy_from_slice(&result);
        d_outputs[2..].fill(0.0);
        Ok(())
    }

    fn vjp(
        &self,
        old: &ClosedGasParcel,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), GasEnergyError> {
        if self.stream.is_some() {
            return Err(invalid(NO_TANGENT));
        }
        if bar_outputs.len() != self.output_width()
            || bar_outputs.iter().any(|value| !value.is_finite())
        {
            return Err(invalid("gas cavity requires one finite tangent weight per output"));
        }
        if self.insulated {
            // One by one: the transpose is the slope itself.
            return self.insulated_slope(old, dt_s, inputs, &bar_outputs[..1], bar_inputs, cx);
        }
        let weights = pair(&bar_outputs[..2])?;
        if bar_inputs.len() != 2 {
            return Err(invalid("gas cavity requires two scalar tangent inputs"));
        }
        let j = self.jacobian(old, dt_s, inputs, cx)?;
        let result = [
            j[0][0] * weights[0] + j[1][0] * weights[1],
            j[0][1] * weights[0] + j[1][1] * weights[1],
        ];
        pair(&result)?;
        poll(cx)?;
        bar_inputs.copy_from_slice(&result);
        Ok(())
    }
}
