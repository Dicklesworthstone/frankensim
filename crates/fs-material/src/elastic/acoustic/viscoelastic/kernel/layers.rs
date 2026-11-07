//! Scalar scattering composition with positive work in each original slab.
//!
//! A real reference impedance at a join is only a choice of wave coordinates:
//! p = a+b, v = (a-b)/Z_ref. Connecting the two wave ports therefore enforces
//! traction and velocity continuity at the SAME plane, with zero propagation.
//! No fluid gap, loss law or extra material is inserted. The suffix scattering
//! recursion eliminates reflected waves without multiplying growing transfer
//! matrices; see Tan (2003), Ultrasonics 41, 229–236,
//! <https://doi.org/10.1016/S0041-624X(02)00447-X>.

use super::{C64, Dd, Drive, Field, Slab, finite};

#[derive(Clone, Copy, Debug)]
pub(crate) struct Layer {
    pub(crate) wave: Slab,
    pub(crate) maximum_strain: f64,
}

#[derive(Clone, Copy, Debug)]
struct PreparedLayer {
    layer: Layer,
    // The layer's two incident reference-port pressures from global inputs.
    inputs: [[C64; 2]; 2],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Error {
    Invalid(&'static str),
    OutsideSourceSupport,
    Cancelled,
}

type Result<T> = core::result::Result<T, Error>;

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

fn round_trip(a: C64, b: C64) -> C64 {
    // Retain the residual of 1-a*b before rounding. Passivity is not imposed
    // by clipping this denominator or by dropping weak reflected channels.
    let (ar, ai, br, bi) = (
        Dd::from_f64(a.re),
        Dd::from_f64(a.im),
        Dd::from_f64(b.re),
        Dd::from_f64(b.im),
    );
    C64::new(
        (Dd::ONE - ar * br + ai * bi).to_f64(),
        -(ar * bi + ai * br).to_f64(),
    )
}

fn linear(row: [C64; 2], input: [C64; 2]) -> C64 {
    let mut real = Dd::ZERO;
    let mut imaginary = Dd::ZERO;
    for (a, b) in row.into_iter().zip(input) {
        let (ar, ai, br, bi) = (
            Dd::from_f64(a.re),
            Dd::from_f64(a.im),
            Dd::from_f64(b.re),
            Dd::from_f64(b.im),
        );
        real = real + ar * br - ai * bi;
        imaginary = imaginary + ar * bi + ai * br;
    }
    C64::new(real.to_f64(), imaginary.to_f64())
}

fn checked(values: impl IntoIterator<Item = C64>) -> Result<()> {
    if values.into_iter().all(|v| finite(v) && v.abs().is_finite()) {
        Ok(())
    } else {
        Err(Error::Invalid(
            "layered acoustic coefficients are unrepresentable",
        ))
    }
}

pub(crate) struct Stack {
    layers: Vec<PreparedLayer>,
    pub(crate) scattering: [C64; 4],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Power {
    pub(crate) incident: [f64; 2],
    pub(crate) outgoing: [f64; 2],
    pub(crate) outgoing_pressure: [C64; 2],
    pub(crate) heat: [f64; 2],
    pub(crate) residual: f64,
}

impl Stack {
    pub(crate) fn new(layers: Vec<Layer>, cancelled: &mut impl FnMut() -> bool) -> Result<Self> {
        poll(cancelled)?;
        if layers.is_empty() {
            return Err(Error::Invalid("a solid stack needs at least one layer"));
        }
        for (i, layer) in layers.iter().enumerate() {
            poll(cancelled)?;
            if !layer.maximum_strain.is_finite()
                || layer.maximum_strain <= 0.0
                || (i > 0 && layers[i - 1].wave.exterior[1] != layer.wave.exterior[0])
                || layer.wave.omega != layers[0].wave.omega
            {
                return Err(Error::Invalid(
                    "layer supports, frequencies and joining coordinates must agree",
                ));
            }
        }
        let mut layers: Vec<_> = layers
            .into_iter()
            .map(|layer| PreparedLayer {
                layer,
                inputs: [[C64::ZERO; 2]; 2],
            })
            .collect();
        let mut suffix = [C64::ZERO, C64::ONE, C64::ONE, C64::ZERO];
        for layer in layers.iter_mut().rev() {
            poll(cancelled)?;
            let s = layer.layer.wave.scattering;
            let denominator = round_trip(s[3], suffix[0]);
            if !finite(denominator) || denominator == C64::ZERO {
                return Err(Error::Invalid("layered acoustic round trip is unresolved"));
            }
            let right = [suffix[0] * s[2] / denominator, suffix[1] / denominator];
            let composed = [
                linear([s[0], s[1]], [C64::ONE, right[0]]),
                s[1] * right[1],
                suffix[2] * (s[2] / denominator),
                linear([suffix[3], suffix[2] * s[3]], [C64::ONE, right[1]]),
            ];
            checked(right.into_iter().chain(composed))?;
            layer.inputs[1] = right;
            suffix = composed;
        }
        // Once prepared, every layer's original owner can be driven directly
        // from the same two exterior inputs, without a per-drive allocation.
        let mut left = [C64::ONE, C64::ZERO];
        for layer in &mut layers {
            poll(cancelled)?;
            let r = layer.inputs[1];
            let right = [r[0] * left[0], linear(r, [left[1], C64::ONE])];
            layer.inputs = [left, right];
            let s = layer.layer.wave.scattering;
            left = [
                linear([s[2], s[3]], [left[0], right[0]]),
                linear([s[2], s[3]], [left[1], right[1]]),
            ];
            checked(layer.inputs.into_iter().flatten().chain(left))?;
        }
        poll(cancelled)?;
        Ok(Self {
            layers,
            scattering: suffix,
        })
    }

    pub(crate) fn layer(&self, index: usize) -> Option<&Slab> {
        self.layers.get(index).map(|p| &p.layer.wave)
    }

    pub(crate) fn layer_drive(&self, index: usize, input: [C64; 2], area: f64) -> Result<Drive> {
        let layer = self
            .layers
            .get(index)
            .ok_or(Error::Invalid("solid layer index is outside the stack"))?;
        let incoming = layer.inputs.map(|row| linear(row, input));
        let drive = layer
            .layer
            .wave
            .drive(incoming, area)
            .map_err(Error::Invalid)?;
        if drive.strain_bound > layer.layer.maximum_strain {
            return Err(Error::OutsideSourceSupport);
        }
        Ok(drive)
    }

    pub(crate) fn field(&self, index: usize, input: [C64; 2], depth: f64) -> Result<Field> {
        let layer = self
            .layers
            .get(index)
            .ok_or(Error::Invalid("solid layer index is outside the stack"))?;
        let incoming = layer.inputs.map(|row| linear(row, input));
        let inward = layer
            .layer
            .wave
            .inward
            .map(|row| row[0] * incoming[0] + row[1] * incoming[1]);
        let parity = layer.layer.wave.parity(incoming);
        layer
            .layer
            .wave
            .field(inward, parity, depth)
            .map_err(Error::Invalid)
    }

    pub(crate) fn drive(
        &self,
        input: [C64; 2],
        area: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Power> {
        poll(cancelled)?;
        let mut incident = [0.0; 2];
        let mut outgoing = [0.0; 2];
        let mut outgoing_pressure = [C64::ZERO; 2];
        let mut heat = [Dd::ZERO; 2];
        for i in 0..self.layers.len() {
            poll(cancelled)?;
            let local = self.layer_drive(i, input, area)?;
            if i == 0 {
                incident[0] = local.incident_power[0];
                outgoing[0] = local.outgoing_power[0];
                outgoing_pressure[0] = local.outgoing[0];
            }
            if i + 1 == self.layers.len() {
                incident[1] = local.incident_power[1];
                outgoing[1] = local.outgoing_power[1];
                outgoing_pressure[1] = local.outgoing[1];
            }
            for (sum, value) in heat.iter_mut().zip(local.heat) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let heat = heat.map(Dd::to_f64);
        let sum = |values: &[f64]| values.iter().fold(Dd::ZERO, |s, x| s + Dd::from_f64(*x));
        let supplied = sum(&incident);
        let spent = sum(&outgoing) + sum(&heat);
        let residual = (supplied - spent).to_f64();
        let scale = (supplied + spent).to_f64();
        if !scale.is_finite()
            || !residual.is_finite()
            || residual.abs()
                > 4096.0 * self.layers.len() as f64 * f64::EPSILON * scale.max(f64::MIN_POSITIVE)
        {
            return Err(Error::Invalid(
                "solid stack boundary powers and positive local heat disagree",
            ));
        }
        poll(cancelled)?;
        Ok(Power {
            incident,
            outgoing,
            outgoing_pressure,
            heat,
            residual,
        })
    }
}

#[cfg(test)]
#[path = "layers/tests.rs"]
mod tests;
