//! Coherent coatings joined by explicitly incoherent, lossless propagation.
//!
//! The caller declares complete decorrelation between distinct reflected
//! traversals of every connection. Thin-film interference remains inside each
//! native block. No coherence length, separation thickness, bandwidth or
//! random phase is inferred from missing data. For a single connection this
//! is also the uniform round-trip-phase average of the coherent solution.
//!
//! Intensity multiple-reflection sums: S. J. Byrnes, "Multilayer optical
//! calculations", section 6, <https://arxiv.org/abs/1603.02720>. This owner
//! restricts connecting ports to the original positive lossless media so
//! incident/reflected powers are additive. Absorption in coherent blocks uses
//! their native positive loss modes, including longitudinal TM fields.

use super::{LayerIlluminationPower, PreparedLayerIllumination};
use super::super::{MAX_LAYERS, bad, same_response_source};
use crate::plane_wave::{Polarization, ReflectionError};
use fs_math::{c64::C64, dd::Dd, det};

const MAX_BLOCKS: usize = MAX_LAYERS + 1;

/// Prepared mixed coherent/incoherent planar power transport.
///
/// Blocks run left to right. Every adjacent pair must retain the same
/// original connecting medium, frequency, polarization and tangential wave
/// number. Both exterior modes of each native block must propagate without
/// loss, as already required by `prepare_illumination`. Evanescent/critical
/// connecting modes are not incoherent power channels; keep their tunnelling
/// and interference inside one coherent block instead.
///
/// Coherent films may absorb. The connections do not absorb, scatter or add
/// interfaces. This model describes an authored fully incoherent propagation
/// limit, not arbitrary partial temporal coherence or ensemble averaging of
/// several static correlated film thicknesses. Independent exterior sources
/// are supplied as normal power densities, never as invented complex fields.
///
/// Preparation is O(blocks + films), bounded to 65 blocks and 64 total films.
/// Evaluation is O(films), allocation-free, and reuses all original Maxwell
/// solves. Borrowed blocks preserve every original source and layer geometry.
pub struct IncoherentLayerCascade<'p, 's, 'm> {
    blocks: Vec<&'p PreparedLayerIllumination<'s, 'm>>,
    scattering: [[f64; 2]; 2],
    // Native block order, then native layer order; columns are exterior inputs.
    absorption: Vec<[f64; 2]>,
}

#[derive(Clone, Copy, Default)]
struct PowerBlock {
    reflected: [f64; 2],
    transmitted: [f64; 2],
    absorbed: [f64; 2],
}

#[derive(Clone, Copy, Default)]
struct Connection {
    // F_right = forward[0] * F_left + forward[1] * exterior_right.
    forward: [f64; 2],
    // B_right = backward[0] * F_left + backward[1] * exterior_right.
    backward: [f64; 2],
}

fn sum(values: &[f64]) -> f64 {
    values.iter().fold(Dd::ZERO, |s, &x| s + Dd::from_f64(x)).to_f64()
}

fn weighted_sum(coefficients: [f64; 2], input: [f64; 2]) -> f64 {
    (Dd::from_f64(coefficients[0]) * Dd::from_f64(input[0])
        + Dd::from_f64(coefficients[1]) * Dd::from_f64(input[1]))
    .to_f64()
}

fn layer_absorption(block: &PreparedLayerIllumination<'_, '_>, modes: [[C64; 2]; 2]) -> [f64; 2] {
    std::array::from_fn(|source| {
        let normalization = core::f64::consts::SQRT_2 / det::sqrt(block.admittance_s[source]);
        sum(&[
            modes[0][source].scale(normalization).norm_sq(),
            modes[1][source].scale(normalization).norm_sq(),
        ])
    })
}

fn native_power(block: &PreparedLayerIllumination<'_, '_>) -> Result<PowerBlock, ReflectionError> {
    let root_y = block.admittance_s.map(det::sqrt);
    let transmitted = std::array::from_fn(|source| {
        let field = block.scattering[1 - source][source].scale(root_y[1 - source]);
        C64::new(field.re / root_y[source], field.im / root_y[source]).norm_sq()
    });
    let mut absorbed = [Dd::ZERO; 2];
    for &modes in &block.loss_modes {
        for (total, a) in absorbed.iter_mut().zip(layer_absorption(block, modes)) {
            *total = *total + Dd::from_f64(a);
        }
    }
    let result = PowerBlock {
        reflected: [block.scattering[0][0].norm_sq(), block.scattering[1][1].norm_sq()],
        transmitted,
        absorbed: absorbed.map(Dd::to_f64),
    };
    for side in 0..2 {
        check_column(
            [result.reflected[side], result.transmitted[side]],
            result.absorbed[side],
            block.loss_modes.len() + 1,
        )?;
    }
    Ok(result)
}

fn check_column(output: [f64; 2], absorbed: f64, count: usize) -> Result<(), ReflectionError> {
    if output.into_iter().chain([absorbed]).any(|x| !x.is_finite() || x < 0.0)
        || (sum(&[output[0], output[1], absorbed]) - 1.0).abs()
            > 2048.0 * f64::EPSILON * count as f64
    {
        return Err(bad("incoherent cascade power balance failed"));
    }
    Ok(())
}

fn compose(left: PowerBlock, right: PowerBlock) -> Result<(PowerBlock, Connection), ReflectionError> {
    // 1 - R_left,right * R_right,left, expressed entirely through positive
    // escape channels. Subtracting rounded unit reflectances would erase
    // weak transmission/absorption or create a false singular cavity.
    let denominator = (Dd::from_f64(left.transmitted[1]) + Dd::from_f64(left.absorbed[1])
        + Dd::from_f64(left.reflected[1])
            * (Dd::from_f64(right.transmitted[0]) + Dd::from_f64(right.absorbed[0])))
    .to_f64();
    if !denominator.is_finite() || denominator < 0.0 {
        return Err(bad("incoherent round-trip escape is unrepresentable"));
    }
    let connection = if denominator == 0.0 {
        // An exactly closed, undriven cavity has no incoming steady power.
        // Native opaque limits on both sides remain separate passive blocks;
        // no arbitrary trapped initial field is introduced.
        if left.transmitted[0] != 0.0 || right.transmitted[1] != 0.0 {
            return Err(bad("incoherent cavity has no representable escape"));
        }
        Connection::default()
    } else {
        let forward = left.transmitted[0] / denominator;
        let backward = right.transmitted[1] / denominator;
        Connection {
            forward: [forward, left.reflected[1] * backward],
            backward: [right.reflected[0] * forward, backward],
        }
    };
    if connection.forward.into_iter().chain(connection.backward).any(|x| !x.is_finite() || x < 0.0) {
        return Err(bad("incoherent internal power is unrepresentable"));
    }
    let f = connection.forward;
    let b = connection.backward;
    Ok((PowerBlock {
        reflected: [
            sum(&[left.reflected[0], left.transmitted[1] * b[0]]),
            sum(&[right.reflected[1], right.transmitted[0] * f[1]]),
        ],
        transmitted: [right.transmitted[0] * f[0], left.transmitted[1] * b[1]],
        absorbed: [
            sum(&[left.absorbed[0], left.absorbed[1] * b[0], right.absorbed[0] * f[0]]),
            sum(&[right.absorbed[1], left.absorbed[1] * b[1], right.absorbed[0] * f[1]]),
        ],
    }, connection))
}

impl<'p, 's, 'm> IncoherentLayerCascade<'p, 's, 'm> {
    /// Declare complete inter-block incoherence while preserving each native
    /// block's coherent response. A single block is the existing independent
    /// two-face illumination operator, expressed in physical input powers.
    pub fn try_new(blocks: &[&'p PreparedLayerIllumination<'s, 'm>]) -> Result<Self, ReflectionError> {
        if blocks.is_empty() || blocks.len() > MAX_BLOCKS {
            return Err(bad("incoherent cascade requires 1..=65 native coherent blocks"));
        }
        let layers = blocks.iter().map(|b| b.loss_modes.len()).sum::<usize>();
        if layers > MAX_LAYERS {
            return Err(bad("incoherent cascade exceeds 64 original films"));
        }
        for pair in blocks.windows(2) {
            let (left, right) = (pair[0], pair[1]);
            if left.stack.omega != right.stack.omega
                || left.polarization != right.polarization
                || !same_response_source(left.stack.exit, right.stack.incident_response)
            {
                return Err(bad("incoherent connections require identical original medium, frequency and polarization"));
            }
            let exit_normal = left.stack.exit_wave.normal(left.incidence_cosine)?;
            let exit_cosine = exit_normal.re / det::sqrt(left.stack.exit_wave.index_squared_ratio.re);
            if (exit_cosine - right.incidence_cosine).abs()
                > 128.0 * f64::EPSILON * exit_cosine.max(right.incidence_cosine)
            {
                return Err(bad("incoherent connections must conserve the tangential wavevector"));
            }
        }

        let n = blocks.len();
        let mut connections = [Connection::default(); MAX_BLOCKS];
        connections[n - 1].backward = [0.0, 1.0];
        let mut total = native_power(blocks[n - 1])?;
        for i in (0..n - 1).rev() {
            (total, connections[i]) = compose(native_power(blocks[i])?, total)?;
        }

        let mut absorption = Vec::with_capacity(layers);
        let mut forward = [1.0, 0.0];
        for (&block, connection) in blocks.iter().zip(connections) {
            let backward = [
                connection.backward[0] * forward[0],
                sum(&[connection.backward[0] * forward[1], connection.backward[1]]),
            ];
            for &modes in &block.loss_modes {
                let loss = layer_absorption(block, modes);
                absorption.push(std::array::from_fn(|s| {
                    weighted_sum(loss, [forward[s], backward[s]])
                }));
            }
            forward = [
                connection.forward[0] * forward[0],
                sum(&[connection.forward[0] * forward[1], connection.forward[1]]),
            ];
        }
        let scattering = [
            [total.reflected[0], total.transmitted[1]],
            [total.transmitted[0], total.reflected[1]],
        ];
        for side in 0..2 {
            let absorbed = absorption.iter().fold(Dd::ZERO, |s, a| s + Dd::from_f64(a[side])).to_f64();
            check_column([scattering[0][side], scattering[1][side]], absorbed, layers + n)?;
        }
        Ok(Self { blocks: blocks.to_vec(), scattering, absorption })
    }

    /// Original coherent blocks, in physical left-to-right order.
    pub fn blocks(&self) -> &[&'p PreparedLayerIllumination<'s, 'm>] {
        &self.blocks
    }

    pub fn polarization(&self) -> Polarization {
        self.blocks[0].polarization
    }

    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.blocks[0].stack.omega
    }

    /// Outgoing-face major, incoming-face minor, in normal power units.
    /// Diagonal entries reflect; off-diagonal entries transmit.
    pub const fn power_scattering(&self) -> [[f64; 2]; 2] {
        self.scattering
    }

    /// Absorbed fractions per original film, each with left/right input
    /// columns. Ordering is block-major, then native layer order. Lossless
    /// connecting propagation has no film heat and adds no synthetic entry.
    pub fn layer_absorption_fractions(&self) -> &[[f64; 2]] {
        &self.absorption
    }

    /// Native film indices in the flattened absorption/output array.
    pub fn layer_range(&self, block: usize) -> Option<core::ops::Range<usize>> {
        let selected = self.blocks.get(block)?;
        let start = self.blocks[..block].iter().map(|b| b.loss_modes.len()).sum();
        Some(start..start + selected.loss_modes.len())
    }

    /// Retained payload, excluding borrowed native fields and allocator overhead.
    pub fn workspace_bytes(&self) -> usize {
        core::mem::size_of::<Self>()
            + self.blocks.capacity() * core::mem::size_of::<&PreparedLayerIllumination<'s, 'm>>()
            + self.absorption.capacity() * core::mem::size_of::<[f64; 2]>()
    }

    /// Apply independent physical normal input powers [W/m²], left then right.
    /// Layer heat is deposited only by the downstream thermal owner. No fields,
    /// optical solve, material lookup, sampling or trigonometry is performed.
    /// Output remains unchanged if admission or power accounting fails.
    pub fn evaluate(
        &self,
        incoming_w_m2: [f64; 2],
        layer_heat_w_m2: &mut [f64],
    ) -> Result<LayerIlluminationPower, ReflectionError> {
        if layer_heat_w_m2.len() != self.absorption.len()
            || incoming_w_m2.iter().any(|p| !p.is_finite() || *p < 0.0)
        {
            return Err(bad("incoherent illumination requires finite nonnegative powers and one output per film"));
        }
        let outgoing = self.scattering.map(|row| weighted_sum(row, incoming_w_m2));
        let mut heat = [0.0; MAX_LAYERS];
        let mut absorbed = Dd::ZERO;
        for (out, coefficients) in heat.iter_mut().zip(&self.absorption) {
            *out = weighted_sum(*coefficients, incoming_w_m2);
            absorbed = absorbed + Dd::from_f64(*out);
        }
        let absorbed = absorbed.to_f64();
        if outgoing.into_iter().chain([absorbed]).any(|x| !x.is_finite() || x < 0.0) {
            return Err(bad("incoherent output power is unrepresentable"));
        }
        let scale = incoming_w_m2.into_iter().chain(outgoing).chain([absorbed]).fold(0.0_f64, f64::max);
        if scale > 0.0
            && sum(&[
                incoming_w_m2[0] / scale,
                incoming_w_m2[1] / scale,
                -outgoing[0] / scale,
                -outgoing[1] / scale,
                -absorbed / scale,
            ]).abs() > 2048.0 * f64::EPSILON * (self.blocks.len() + self.absorption.len()) as f64
        {
            return Err(bad("incoherent illuminated power balance failed"));
        }
        layer_heat_w_m2.copy_from_slice(&heat[..self.absorption.len()]);
        Ok(LayerIlluminationPower { incoming_w_m2, outgoing_w_m2: outgoing, absorbed_w_m2: absorbed })
    }
}
