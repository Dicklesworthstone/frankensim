//! Streams: a volume of perfect gas and a linear restriction between volumes.
//!
//! [`GasVolume`] is a storage. Through its stream port it receives what
//! crosses into it (mass, enthalpy and the carried quantities) and returns
//! what it holds is like (pressure, specific enthalpy, each carried quantity
//! per unit mass). It also has a heat port and, if asked, a moving wall.
//! [`Restriction`] is a transport element between two storages: the mass
//! rate follows the pressure difference, and the enthalpy and carried
//! quantities that go with it are those of the side the mass leaves.
//!
//! The gas has constant heat capacities. A gas whose properties come from
//! sources belongs to its physics owner (`fs_gas::component::GasCavity`).

use fs_exec::Cx;
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, STREAM_CARRIED,
    STREAM_ENTHALPY_RATE, STREAM_MASS_RATE, STREAM_PRESSURE, STREAM_SPECIFIC_ENTHALPY,
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
/// Entropy made by mixing a stream into the volume is not reported.
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
        let supplied = dt_s * (inputs[STREAM_ENTHALPY_RATE] + heat_w);
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
            scale_j: old.energy_j + energy_j + supplied.abs() + work_j.abs(),
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
                kelvin -=
                    0.5 * trial.energy_j * dt_s / (trial.mass_kg * trial.mass_kg * capacity);
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
            matrix[row + STREAM_MASS_RATE] =
                -0.5 * new * dt_s / (trial.mass_kg * trial.mass_kg);
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
        // A carried quantity cannot be drawn below nothing.
        let rates = &inputs[STREAM_CARRIED..STREAM_CARRIED + self.carried];
        if old
            .carried_kg
            .iter()
            .zip(rates)
            .any(|(held, rate)| !(held + dt_s * rate >= 0.0))
        {
            return Ok(shorter);
        }
        new.mass_kg = trial.mass_kg;
        new.energy_j = trial.energy_j;
        new.volume_m3 = trial.volume_m3;
        new.carried_kg.clear();
        for (index, (held, rate)) in old.carried_kg.iter().zip(rates).enumerate() {
            let carried = held + dt_s * rate;
            new.carried_kg.push(carried);
            outputs[STREAM_CARRIED + index] =
                0.5 * (held / old.mass_kg + carried / trial.mass_kg);
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
            *slot = (0..n).map(|column| matrix[row * n + column] * d_inputs[column]).sum();
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
            *slot = (0..n).map(|row| matrix[row * n + column] * bar_outputs[row]).sum();
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
/// supplies the carried values changes. Entropy made by the pressure drop
/// appears in the receiving storage's state and is not reported here.
#[derive(Debug, Clone)]
pub struct Restriction {
    conductance_kg_s_pa: f64,
    carried: usize,
    ports: [PortSpec; 2],
}

impl Restriction {
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
        Ok(Proposal::Step(Receipt::lossless(0.0, roundoff(scale))))
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
