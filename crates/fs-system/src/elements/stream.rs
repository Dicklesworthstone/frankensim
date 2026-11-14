//! Streams: a volume of perfect gas and a linear restriction between volumes.
//!
//! [`GasVolume`] is a storage. Through its stream port it receives what
//! crosses into it (mass, enthalpy and the carried quantities) and returns
//! what it holds is like (pressure, specific enthalpy, each carried quantity
//! per unit mass). It also has a heat port and, if asked, a moving wall.
//! [`Restriction`] is a transport element between two storages: the mass
//! rate follows the pressure difference, and the enthalpy and carried
//! quantities that go with it are those of the side the mass leaves. Built
//! for a gas ([`Restriction::for_gas`]), it reports the entropy its pressure
//! drop makes, from that gas's own state function.
//!
//! The gas has constant heat capacities. A gas whose properties come from
//! sources belongs to its physics owner (`fs_gas::component::GasCavity`).

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::component::{
    Carrier, Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, STREAM_CARRIED,
    STREAM_ENTHALPY_RATE, STREAM_MASS_RATE, STREAM_PRESSURE, STREAM_SPECIFIC_ENTHALPY, StreamBasis,
    StreamQuantity,
};
use fs_ports::{PortKind, StableId};

use super::{ElementError, finite, heat, inputs_finite, port, positive, roundoff, temperature};

fn stream(name: &str, role: PortRole, carried: usize) -> PortSpec {
    let id = StableId::new(name).expect("element port names are valid stable ids");
    PortSpec::stream(id, role, carried).expect("element stream ports have a power role")
}

/// What a [`GasVolume`] holds.
#[derive(Debug, Clone, PartialEq)]
pub struct GasContents {
    /// Mass [kg].
    pub mass_kg: f64,
    /// Internal energy [J], zero at absolute zero.
    pub energy_j: f64,
    /// Volume [m^3].
    pub volume_m3: f64,
    /// Mass [kg] of each carried quantity, a part of `mass_kg`.
    pub carried_kg: Vec<f64>,
}

/// A well-mixed volume of perfect gas: `p V = m R T` and `U = m c_v T` with
/// `c_v = R / (ratio - 1)`.
///
/// Ports, in this order:
/// - `stream` (a storage's stream port, with `carried` quantities): receives
///   mass rate [kg/s], enthalpy rate [W] and the rate of each carried
///   quantity [kg/s]; returns pressure [Pa], specific enthalpy `c_p T`
///   [J/kg] and the mass fraction of each carried quantity.
/// - `heat`: receives the heat rate into the gas [W], returns its
///   temperature [K].
/// - `p`, only [`Self::with_wall`]: receives the rate at which a wall takes
///   volume from the gas [m^3/s] and returns the pressure on it [Pa].
///
/// One step is the trapezoidal rule. Mass and carried masses move by what
/// crossed. The energy moves by the enthalpy and heat that entered and the
/// work of the mean of the old and new pressure, which is the pressure the
/// ports return, so the first law holds for any step. The temperature and
/// the fractions returned are the means of their old and new values.
/// Entropy made by mixing a stream into the volume is not reported, so the
/// entropy a system of volumes shows generated is a lower bound of what it
/// made.
///
/// A step that would empty the volume, close it, take more energy than it
/// holds or draw a carried quantity below nothing is shortened. So a small
/// volume between restrictions sets the step: more than about twice its
/// mass cannot pass through it in one. A prescribed draw of a quantity the
/// volume holds none of is never satisfied and ends in step underflow.
///
/// [`Self::identified`] binds gas constants, the energy datum and named
/// carried masses to the port, so incompatible streams fail admission.
/// [`Self::new`] leaves these identities unspecified: its caller must ensure
/// every joined volume holds the same gas and agrees on coordinate meanings.
#[derive(Debug, Clone)]
pub struct GasVolume {
    gas_constant_j_kg_k: f64,
    ratio: f64,
    carried: usize,
    ports: Vec<PortSpec>,
}

/// One trapezoidal step of a [`GasVolume`], before it is published.
struct Trial {
    mass_kg: f64,
    energy_j: f64,
    volume_m3: f64,
    /// `1 - (ratio - 1) dV / (2 V1)`, the pivot of the energy equation.
    pivot: f64,
    pressure_pa: f64,
    temperature_k: f64,
    change_j: f64,
    scale_j: f64,
}

impl GasVolume {
    /// The same perfect-gas law with identified carried masses. Gas constants
    /// and the zero-at-absolute-zero energy datum are bound automatically;
    /// the supplied quantities retain their order. Composition and live
    /// temperature do not change the basis. Named quantities must be masses.
    ///
    /// To connect it through a restriction, construct that restriction with
    /// [`Restriction::for_stream`] from this volume's `stream` port. A volume
    /// with different gas constants or quantity meanings cannot share it.
    ///
    /// # Errors
    /// The refusals of [`Self::new`], duplicate quantity identities, or a
    /// carried inventory whose dimensions are not mass.
    pub fn identified(
        gas_constant_j_kg_k: f64,
        heat_capacity_ratio: f64,
        quantities: Vec<StreamQuantity>,
    ) -> Result<Self, ElementError> {
        let mut volume = Self::new(gas_constant_j_kg_k, heat_capacity_ratio, quantities.len())?;
        if let Some(index) = quantities
            .iter()
            .position(|q| q.dimensions().0 != [0, 1, 0, 0, 0, 0])
        {
            return Err(ElementError::Parameter {
                name: "stream quantity must be a mass [kg]",
                value: index as f64,
            });
        }
        let medium = StableId::new(format!(
            "perfect-gas-v1/{:016x}/{:016x}",
            gas_constant_j_kg_k.to_bits(),
            heat_capacity_ratio.to_bits()
        ))
        .expect("finite gas constants encoded as hexadecimal");
        let basis = StreamBasis::new(
            medium,
            StableId::new("ideal-gas-u-zero-at-absolute-zero-v1").expect("static reference"),
            quantities,
        )
        .map_err(ElementError::StreamBasis)?;
        volume.ports[0] = PortSpec::identified_stream(
            StableId::new("stream").expect("static port name"),
            PortRole::FlowIn,
            basis,
        )
        .map_err(ElementError::StreamBasis)?;
        Ok(volume)
    }

    /// A volume of gas with specific gas constant `gas_constant_j_kg_k`
    /// [J/(kg K)] and heat capacity ratio `heat_capacity_ratio`, whose stream
    /// port carries `carried` quantities besides mass and enthalpy.
    ///
    /// # Errors
    /// A gas constant that is not positive and finite, or a ratio that is not
    /// finite and above one.
    pub fn new(
        gas_constant_j_kg_k: f64,
        heat_capacity_ratio: f64,
        carried: usize,
    ) -> Result<Self, ElementError> {
        let ratio = finite("heat_capacity_ratio", heat_capacity_ratio)?;
        if ratio <= 1.0 {
            return Err(ElementError::Parameter {
                name: "heat_capacity_ratio",
                value: ratio,
            });
        }
        Ok(Self {
            gas_constant_j_kg_k: positive("gas_constant_j_kg_k", gas_constant_j_kg_k)?,
            ratio,
            carried,
            ports: vec![
                stream("stream", PortRole::FlowIn, carried),
                heat("heat", PortRole::FlowIn),
            ],
        })
    }

    /// The same gas behind a moving wall: adds port `p`.
    #[must_use]
    pub fn with_wall(mut self) -> Self {
        if self.ports.len() == 2 {
            self.ports.push(port(
                "p",
                Some(PortKind::FluidPressureFlux),
                PortRole::FlowIn,
            ));
        }
        self
    }

    /// The contents at `pressure_pa`, `temperature_k` and `volume_m3`, with
    /// the given mass fraction of each carried quantity.
    ///
    /// # Errors
    /// A pressure, temperature or volume that is not positive and finite; a
    /// number of fractions other than the carried quantities; a fraction
    /// outside `[0, 1]`.
    pub fn contents(
        &self,
        pressure_pa: f64,
        temperature_k: f64,
        volume_m3: f64,
        carried_fractions: &[f64],
    ) -> Result<GasContents, ElementError> {
        let pressure = positive("pressure_pa", pressure_pa)?;
        let kelvin = temperature(temperature_k)?;
        let volume = positive("volume_m3", volume_m3)?;
        if carried_fractions.len() != self.carried {
            return Err(ElementError::Parameter {
                name: "carried fractions (one per carried quantity)",
                value: carried_fractions.len() as f64,
            });
        }
        if let Some(&value) = carried_fractions
            .iter()
            .find(|fraction| !(0.0..=1.0).contains(*fraction))
        {
            return Err(ElementError::Parameter {
                name: "carried fraction",
                value,
            });
        }
        let mass_kg = pressure * volume / (self.gas_constant_j_kg_k * kelvin);
        Ok(GasContents {
            mass_kg: positive("mass_kg", mass_kg)?,
            energy_j: positive("energy_j", pressure * volume / (self.ratio - 1.0))?,
            volume_m3: volume,
            carried_kg: carried_fractions
                .iter()
                .map(|fraction| fraction * mass_kg)
                .collect(),
        })
    }

    /// Pressure [Pa] of `contents`.
    #[must_use]
    pub fn pressure_pa(&self, contents: &GasContents) -> f64 {
        (self.ratio - 1.0) * contents.energy_j / contents.volume_m3
    }

    /// Temperature [K] of `contents`.
    #[must_use]
    pub fn temperature_k(&self, contents: &GasContents) -> f64 {
        contents.energy_j / (contents.mass_kg * self.heat_capacity_j_kg_k())
    }

    /// Heat capacity at constant volume [J/(kg K)].
    #[must_use]
    pub fn heat_capacity_j_kg_k(&self) -> f64 {
        self.gas_constant_j_kg_k / (self.ratio - 1.0)
    }

    /// Specific entropy [J/(kg K)] of this gas at specific enthalpy `to_j_kg`
    /// [J/kg] and pressure `to_pa` [Pa], less that at `from_j_kg` and
    /// `from_pa`: `c_p ln(h_to / h_from) - R ln(p_to / p_from)`, since
    /// `h = c_p T`. Each logarithm is taken as `ln_1p` of the relative
    /// change, so equal values give exactly zero and each term has the sign
    /// its change gives it, however small the change.
    #[must_use]
    pub fn specific_entropy_change_j_kg_k(
        &self,
        from_j_kg: f64,
        from_pa: f64,
        to_j_kg: f64,
        to_pa: f64,
    ) -> f64 {
        let capacity = self.ratio * self.heat_capacity_j_kg_k();
        capacity * ((to_j_kg - from_j_kg) / from_j_kg).ln_1p()
            - self.gas_constant_j_kg_k * ((to_pa - from_pa) / from_pa).ln_1p()
    }

    fn has_wall(&self) -> bool {
        self.ports.len() == 3
    }

    fn admit(&self, contents: &GasContents) -> Result<(), ElementError> {
        positive("mass_kg", contents.mass_kg)?;
        positive("energy_j", contents.energy_j)?;
        positive("volume_m3", contents.volume_m3)?;
        if contents.carried_kg.len() != self.carried {
            return Err(ElementError::Parameter {
                name: "carried masses (one per carried quantity)",
                value: contents.carried_kg.len() as f64,
            });
        }
        match contents
            .carried_kg
            .iter()
            .find(|mass| !(mass.is_finite() && **mass >= 0.0))
        {
            Some(&value) => Err(ElementError::Parameter {
                name: "carried mass [kg]",
                value,
            }),
            None => Ok(()),
        }
    }

    /// Volume the wall takes from the gas per second, or zero without a wall.
    fn compression(&self, inputs: &[f64]) -> f64 {
        if self.has_wall() {
            inputs[STREAM_CARRIED + self.carried + 1]
        } else {
            0.0
        }
    }

    /// The step, or `None` when it is too long for what it carries: it would
    /// empty the volume, close it, or take more work than it holds.
    fn trial(&self, old: &GasContents, dt_s: f64, inputs: &[f64]) -> Option<Trial> {
        let heat_w = inputs[STREAM_CARRIED + self.carried];
        let mass_kg = old.mass_kg + dt_s * inputs[STREAM_MASS_RATE];
        let displaced = dt_s * self.compression(inputs);
        let volume_m3 = old.volume_m3 - displaced;
        let half = 0.5 * (self.ratio - 1.0) * displaced;
        let pivot = 1.0 - half / volume_m3;
        if !(mass_kg > 0.0 && volume_m3 > 0.0 && pivot > 0.0) {
            return None;
        }
        // U1 - U0 = E + (p0 + p1) / 2 * dV with p = (ratio - 1) U / V is
        // linear in U1.
        let enthalpy_w = inputs[STREAM_ENTHALPY_RATE];
        let supplied = dt_s * (enthalpy_w + heat_w);
        let solved = (old.energy_j * (1.0 + half / old.volume_m3) + supplied) / pivot;
        let pressure_pa =
            0.5 * (self.ratio - 1.0) * (old.energy_j / old.volume_m3 + solved / volume_m3);
        // The change is formed from the terms the ports book, so the first
        // law is not left to the rounding of the division above.
        let work_j = pressure_pa * displaced;
        let change_j = supplied + work_j;
        let energy_j = old.energy_j + change_j;
        if !(energy_j > 0.0 && energy_j.is_finite() && mass_kg.is_finite()) {
            return None;
        }
        let capacity = self.heat_capacity_j_kg_k();
        Some(Trial {
            mass_kg,
            energy_j,
            volume_m3,
            pivot,
            pressure_pa,
            temperature_k: 0.5
                * (old.energy_j / (old.mass_kg * capacity) + energy_j / (mass_kg * capacity)),
            change_j,
            // The ports book the enthalpy and the heat separately, so the
            // allowance is sized by each and not by their sum, which can
            // cancel.
            scale_j: old.energy_j
                + energy_j
                + dt_s * (enthalpy_w.abs() + heat_w.abs())
                + work_j.abs(),
        })
    }

    /// The derivative of every output with respect to every input, row-major.
    fn jacobian(
        &self,
        old: &GasContents,
        dt_s: f64,
        inputs: &[f64],
    ) -> Result<Vec<f64>, ElementError> {
        inputs_finite(inputs)?;
        self.admit(old)?;
        let Some(trial) = self.trial(old, dt_s, inputs) else {
            return Err(ElementError::Parameter {
                name: "step_s (too long for what the gas volume receives)",
                value: dt_s,
            });
        };
        let k = self.carried;
        let n = self.ports.iter().map(PortSpec::width).sum::<usize>();
        let (heat_at, wall_at) = (STREAM_CARRIED + k, STREAM_CARRIED + k + 1);
        let capacity = self.heat_capacity_j_kg_k();
        let gamma_1 = self.ratio - 1.0;
        // Derivative of the new energy with respect to each input.
        let mut energy = vec![0.0; n];
        energy[STREAM_ENTHALPY_RATE] = dt_s / trial.pivot;
        energy[heat_at] = dt_s / trial.pivot;
        // And of the new pressure.
        let mut pressure = vec![0.0; n];
        if self.has_wall() {
            let half = 0.5 * gamma_1 * dt_s * self.compression(inputs);
            let d_half = 0.5 * gamma_1 * dt_s;
            let d_pivot =
                -(d_half / trial.volume_m3 + half * dt_s / (trial.volume_m3 * trial.volume_m3));
            energy[wall_at] =
                (old.energy_j * d_half / old.volume_m3 - trial.energy_j * d_pivot) / trial.pivot;
            pressure[wall_at] =
                gamma_1 * trial.energy_j * dt_s / (trial.volume_m3 * trial.volume_m3);
        }
        for (slope, energy) in pressure.iter_mut().zip(&energy) {
            *slope += gamma_1 * energy / trial.volume_m3;
        }
        let mut matrix = vec![0.0; n * n];
        for column in 0..n {
            let mut kelvin = 0.5 * energy[column] / (trial.mass_kg * capacity);
            if column == STREAM_MASS_RATE {
                kelvin -= 0.5 * trial.energy_j * dt_s / (trial.mass_kg * trial.mass_kg * capacity);
            }
            matrix[STREAM_PRESSURE * n + column] = 0.5 * pressure[column];
            matrix[STREAM_SPECIFIC_ENTHALPY * n + column] = self.ratio * capacity * kelvin;
            matrix[heat_at * n + column] = kelvin;
            if self.has_wall() {
                matrix[wall_at * n + column] = 0.5 * pressure[column];
            }
        }
        for (index, held) in old.carried_kg.iter().enumerate() {
            let row = (STREAM_CARRIED + index) * n;
            let new = held + dt_s * inputs[STREAM_CARRIED + index];
            matrix[row + STREAM_CARRIED + index] = 0.5 * dt_s / trial.mass_kg;
            matrix[row + STREAM_MASS_RATE] = -0.5 * new * dt_s / (trial.mass_kg * trial.mass_kg);
        }
        Ok(matrix)
    }
}

impl Component for GasVolume {
    type State = GasContents;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, contents: &GasContents) -> Result<f64, ElementError> {
        Ok(contents.energy_j)
    }

    fn stored_mass_kg(&self, contents: &GasContents) -> Result<Option<f64>, ElementError> {
        self.admit(contents)?;
        Ok(Some(contents.mass_kg))
    }

    fn stored_carried_quantity(
        &self,
        contents: &GasContents,
        carried: usize,
    ) -> Result<Option<f64>, ElementError> {
        Ok(contents.carried_kg.get(carried).copied())
    }

    fn propose(
        &self,
        old: &GasContents,
        dt_s: f64,
        inputs: &[f64],
        new: &mut GasContents,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        self.admit(old)?;
        let shorter = Proposal::Shorten { dt_s: 0.5 * dt_s };
        let Some(trial) = self.trial(old, dt_s, inputs) else {
            return Ok(shorter);
        };
        // A carried quantity cannot be drawn below nothing. Its rounding
        // scale is its own inventory and transfer, never the bulk gas mass:
        // otherwise an entire trace constituent can disappear as roundoff.
        let rates = &inputs[STREAM_CARRIED..STREAM_CARRIED + self.carried];
        if old
            .carried_kg
            .iter()
            .zip(rates)
            .any(|(held, rate)| {
                let transfer = dt_s * rate;
                let remaining = held + transfer;
                let rounding = 64.0 * f64::EPSILON * held.abs().max(transfer.abs());
                !(transfer.is_finite() && remaining.is_finite() && remaining >= -rounding)
            })
        {
            return Ok(shorter);
        }
        new.mass_kg = trial.mass_kg;
        new.energy_j = trial.energy_j;
        new.volume_m3 = trial.volume_m3;
        new.carried_kg.clear();
        for (index, (held, rate)) in old.carried_kg.iter().zip(rates).enumerate() {
            let carried = (held + dt_s * rate).max(0.0);
            new.carried_kg.push(carried);
            outputs[STREAM_CARRIED + index] = 0.5 * (held / old.mass_kg + carried / trial.mass_kg);
        }
        outputs[STREAM_PRESSURE] = trial.pressure_pa;
        outputs[STREAM_SPECIFIC_ENTHALPY] =
            self.ratio * self.heat_capacity_j_kg_k() * trial.temperature_k;
        outputs[STREAM_CARRIED + self.carried] = trial.temperature_k;
        if self.has_wall() {
            outputs[STREAM_CARRIED + self.carried + 1] = trial.pressure_pa;
        }
        Ok(Proposal::Step(Receipt::lossless(
            trial.change_j,
            roundoff(trial.scale_j),
        )))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::GasVolume",
        state_family: "gas-contents",
        schema_grammar: "mass_kg f64; energy_j f64; volume_m3 f64; carried_kg: one f64 per \
            carried quantity of the volume",
        codec_grammar: "u64 words, the bits of mass_kg, energy_j, volume_m3 and each carried \
            mass in that order: 3 + carried words",
        codec_version: 1,
    });

    fn encode_state(&self, contents: &GasContents, words: &mut Vec<u64>) -> bool {
        let held = [contents.mass_kg, contents.energy_j, contents.volume_m3];
        words.extend(
            held.iter()
                .chain(&contents.carried_kg)
                .map(|value| value.to_bits()),
        );
        true
    }

    /// Refuses another count than `3 + carried`, and contents no volume
    /// holds (the checks of every step's old state).
    fn decode_state(&self, words: &[u64]) -> Option<GasContents> {
        if words.len() != 3 + self.carried {
            return None;
        }
        let value = |index: usize| f64::from_bits(words[index]);
        let contents = GasContents {
            mass_kg: value(0),
            energy_j: value(1),
            volume_m3: value(2),
            carried_kg: words[3..]
                .iter()
                .map(|&word| f64::from_bits(word))
                .collect(),
        };
        self.admit(&contents).ok().map(|()| contents)
    }

    /// The gas constant and heat capacity ratio, as bits, and the number of
    /// carried quantities.
    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            self.gas_constant_j_kg_k.to_bits(),
            self.ratio.to_bits(),
            self.carried as u64,
        ]);
    }
}

impl PortTangent for GasVolume {
    fn jvp(
        &self,
        old: &GasContents,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let matrix = self.jacobian(old, dt_s, inputs)?;
        let n = d_inputs.len();
        for (row, slot) in d_outputs.iter_mut().enumerate() {
            *slot = (0..n)
                .map(|column| matrix[row * n + column] * d_inputs[column])
                .sum();
        }
        Ok(())
    }

    fn vjp(
        &self,
        old: &GasContents,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let matrix = self.jacobian(old, dt_s, inputs)?;
        let n = bar_outputs.len();
        for (column, slot) in bar_inputs.iter_mut().enumerate() {
            *slot = (0..n)
                .map(|row| matrix[row * n + column] * bar_outputs[row])
                .sum();
        }
        Ok(())
    }
}

/// A linear restriction between two storages, with stream ports `a` and `b`.
///
/// The mass rate from `a` to `b` is `conductance (p_a - p_b)`. The enthalpy
/// and the carried quantities that cross are that mass rate times the values
/// of the side the mass leaves, so what one storage loses the other gains,
/// and a flow that turns round passes through zero without a jump. It stores
/// nothing. The port map has a kink at zero flow, where the side that
/// supplies the carried values changes.
///
/// Built by [`Self::for_gas`], it reports the entropy the pressure drop
/// makes. Nothing stores or exchanges heat on the way, so what crosses
/// leaves the upstream side at that side's specific enthalpy and pressure
/// and arrives at the downstream pressure with the same specific enthalpy.
/// The entropy of a step is `dt |m_dot|` times the gas's
/// [`GasVolume::specific_entropy_change_j_kg_k`] between those two states.
/// It is exactly zero when nothing crosses and the same whichever way the
/// mass goes. Built otherwise, it reports none, and that entropy appears
/// only in the receiving storage's state. The entropy made as the stream
/// mixes into the receiving [`GasVolume`] is not reported by either, so
/// what the ledger shows generated is a lower bound.
#[derive(Debug, Clone)]
pub struct Restriction {
    conductance_kg_s_pa: f64,
    carried: usize,
    ports: [PortSpec; 2],
    /// The gas whose state function gives the entropy of the pressure drop.
    gas: Option<GasVolume>,
}

impl Restriction {
    /// A restriction that preserves the exact stream basis and width of
    /// `port` at both ends. This works for identified source-backed streams
    /// and for unidentified streams; it never erases an existing identity.
    ///
    /// # Errors
    /// A non-stream port, or the conductance refusal of [`Self::new`].
    pub fn for_stream(conductance_kg_s_pa: f64, port: &PortSpec) -> Result<Self, ElementError> {
        if port.carrier() != Carrier::Stream {
            return Err(ElementError::Parameter {
                name: "restriction requires a stream port",
                value: port.width() as f64,
            });
        }
        let mut restriction = Self::new(conductance_kg_s_pa, port.carried_quantities())?;
        if let Some(basis) = port.stream_basis() {
            for (spec, name) in restriction.ports.iter_mut().zip(["a", "b"]) {
                *spec = PortSpec::identified_stream(
                    StableId::new(name).expect("static port name"),
                    PortRole::EffortIn,
                    basis.clone(),
                )
                .map_err(ElementError::StreamBasis)?;
            }
        }
        Ok(restriction)
    }

    /// A restriction between volumes of `gas` that reports the entropy its
    /// pressure drop makes. Its ports are built from `gas`'s stream port as
    /// by [`Self::for_stream`].
    ///
    /// # Errors
    /// The refusals of [`Self::for_stream`].
    pub fn for_gas(conductance_kg_s_pa: f64, gas: &GasVolume) -> Result<Self, ElementError> {
        let mut restriction = Self::for_stream(conductance_kg_s_pa, &gas.ports[0])?;
        restriction.gas = Some(gas.clone());
        Ok(restriction)
    }

    /// A restriction that passes `conductance_kg_s_pa` [kg/(s Pa)] between
    /// stream ports with `carried` quantities.
    ///
    /// # Errors
    /// A conductance that is not positive and finite.
    pub fn new(conductance_kg_s_pa: f64, carried: usize) -> Result<Self, ElementError> {
        Ok(Self {
            conductance_kg_s_pa: positive("conductance_kg_s_pa", conductance_kg_s_pa)?,
            carried,
            ports: [
                stream("a", PortRole::EffortIn, carried),
                stream("b", PortRole::EffortIn, carried),
            ],
            gas: None,
        })
    }

    fn width(&self) -> usize {
        STREAM_CARRIED + self.carried
    }

    /// Mass rate from `a` to `b` [kg/s] and where the side it leaves starts
    /// in the input vector.
    fn flow(&self, inputs: &[f64]) -> (f64, usize) {
        let n = self.width();
        let rate =
            self.conductance_kg_s_pa * (inputs[STREAM_PRESSURE] - inputs[n + STREAM_PRESSURE]);
        (rate, if rate >= 0.0 { 0 } else { n })
    }

    /// Entropy [J/K] the pressure drop makes in a step of `dt_s` at mass
    /// rate `rate` leaving the side that starts at `upstream`.
    ///
    /// # Errors
    /// A pressure or specific enthalpy that is not positive and finite on a
    /// side the state function reads, while mass crosses.
    fn throttled(
        gas: &GasVolume,
        dt_s: f64,
        rate: f64,
        upstream: usize,
        downstream: usize,
        inputs: &[f64],
    ) -> Result<f64, ElementError> {
        if rate == 0.0 {
            return Ok(0.0);
        }
        let enthalpy = positive(
            "specific_enthalpy_j_kg (upstream)",
            inputs[upstream + STREAM_SPECIFIC_ENTHALPY],
        )?;
        let from_pa = positive("pressure_pa (upstream)", inputs[upstream + STREAM_PRESSURE])?;
        let to_pa = positive(
            "pressure_pa (downstream)",
            inputs[downstream + STREAM_PRESSURE],
        )?;
        let gained = gas.specific_entropy_change_j_kg_k(enthalpy, from_pa, enthalpy, to_pa);
        Ok(dt_s * rate.abs() * gained)
    }
}

impl Component for Restriction {
    type State = ();
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn stored_mass_kg(&self, (): &()) -> Result<Option<f64>, ElementError> {
        Ok(Some(0.0))
    }

    fn stored_carried_quantity(
        &self,
        (): &(),
        carried: usize,
    ) -> Result<Option<f64>, ElementError> {
        Ok((carried < self.ports[0].carried_quantities()).then_some(0.0))
    }

    fn propose(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        (): &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        let n = self.width();
        let (rate, upstream) = self.flow(inputs);
        outputs[STREAM_MASS_RATE] = rate;
        outputs[n + STREAM_MASS_RATE] = -rate;
        // Position 1 is the specific enthalpy and the enthalpy rate; the
        // carried values and their rates follow in the same places.
        for index in STREAM_SPECIFIC_ENTHALPY..n {
            let crossing = rate * inputs[upstream + index];
            outputs[index] = crossing;
            outputs[n + index] = -crossing;
        }
        let scale = dt_s * outputs[STREAM_ENTHALPY_RATE].abs();
        let mut receipt = Receipt::lossless(0.0, roundoff(scale));
        if let Some(gas) = &self.gas {
            receipt.entropy_generated_j_per_k =
                Self::throttled(gas, dt_s, rate, upstream, n - upstream, inputs)?;
        }
        Ok(Proposal::Step(receipt))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::Restriction",
        state_family: "no-state",
        schema_grammar: "(): a restriction holds nothing",
        codec_grammar: "no words",
        codec_version: 1,
    });

    fn encode_state(&self, (): &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    /// The conductance, as bits, and the number of carried quantities;
    /// then 0, or 1 and the words of the gas whose entropy it reports.
    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([self.conductance_kg_s_pa.to_bits(), self.carried as u64]);
        match &self.gas {
            None => words.push(0),
            Some(gas) => {
                words.push(1);
                gas.identity_words(words);
            }
        }
    }
}

impl PortTangent for Restriction {
    fn jvp(
        &self,
        (): &(),
        _dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let n = self.width();
        let (rate, upstream) = self.flow(inputs);
        let d_rate =
            self.conductance_kg_s_pa * (d_inputs[STREAM_PRESSURE] - d_inputs[n + STREAM_PRESSURE]);
        d_outputs[STREAM_MASS_RATE] = d_rate;
        d_outputs[n + STREAM_MASS_RATE] = -d_rate;
        for index in STREAM_SPECIFIC_ENTHALPY..n {
            let crossing = inputs[upstream + index] * d_rate + rate * d_inputs[upstream + index];
            d_outputs[index] = crossing;
            d_outputs[n + index] = -crossing;
        }
        Ok(())
    }

    fn vjp(
        &self,
        (): &(),
        _dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let n = self.width();
        let (rate, upstream) = self.flow(inputs);
        bar_inputs.fill(0.0);
        // Each crossing rate enters `a` and leaves `b`.
        let mut bar_rate = bar_outputs[STREAM_MASS_RATE] - bar_outputs[n + STREAM_MASS_RATE];
        for index in STREAM_SPECIFIC_ENTHALPY..n {
            let bar_crossing = bar_outputs[index] - bar_outputs[n + index];
            bar_rate += inputs[upstream + index] * bar_crossing;
            bar_inputs[upstream + index] = rate * bar_crossing;
        }
        bar_inputs[STREAM_PRESSURE] += self.conductance_kg_s_pa * bar_rate;
        bar_inputs[n + STREAM_PRESSURE] -= self.conductance_kg_s_pa * bar_rate;
        Ok(())
    }
}
