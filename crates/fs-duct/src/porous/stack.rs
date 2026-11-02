//! Coherent fields and heat in heterogeneous finite porous stacks.
//!
//! Every physical layer retains its original material, saturating fluid and
//! JCA(L) response. Native real-reference scattering composes pressure and
//! superficial volume flow at open-pore interfaces. No extra interface loss,
//! effective-medium average, pore-area correction or rigid backing is inferred.
//! Layer order is physical, from the left incident face to the right face.

use super::{
    DuctError, PorousLayer, PorousLayerDissipation, PorousLayerExcitation, PorousLayerField,
    PorousLayerPower, bad, finite, power_for_port,
};
use crate::compose::{ComposeOptions, TwoPort, feedback};
use fs_math::{c64::C64, dd::Dd};

/// Original ordered layers and one prepared external-to-local incident map.
/// Preparation is O(layers), bounded by the native composition work/conditioning
/// policy. Only native decaying scattering blocks are composed: reconstruction
/// never divides by a transmission coefficient or propagates a growing wave.
#[derive(Debug)]
pub struct PorousLayerStack<'l, 'm> {
    layers: &'l [PorousLayer<'m>],
    depths_m: Vec<f64>,
    incident_maps: Vec<[C64; 4]>,
    two_port: TwoPort,
    condition_estimate: f64,
}

/// A single coherent two-sided drive on the entire physical stack.
/// Each local field uses the existing `PorousLayerExcitation` owner, prepared
/// once with its solved incoming waves. Repeated local queries allocate nothing.
#[derive(Debug, Clone)]
pub struct PorousStackExcitation<'s, 'l, 'm> {
    stack: &'s PorousLayerStack<'l, 'm>,
    layers: Vec<PorousLayerExcitation<'l, 'm>>,
    boundary: PorousLayerPower,
    dissipation: PorousLayerDissipation,
    power_residual_w: f64,
}

/// A global depth query with its original local material association.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PorousStackField {
    /// Original zero-based physical layer index.
    pub layer_index: usize,
    /// Depth from that layer's left face [m].
    pub local_depth_m: f64,
    /// Original pressure, superficial flow and local mechanism losses.
    pub field: PorousLayerField,
}

fn poll(cancelled: &mut impl FnMut() -> bool, completed_steps: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps })
    } else {
        Ok(())
    }
}

impl<'l, 'm> PorousLayerStack<'l, 'm> {
    /// Bind an admitted nonempty physical stack with one frequency, full face
    /// area and positive real wave reference. Every original layer is retained.
    /// Native `maximum_blocks` and `maximum_condition` bound work/feedback;
    /// reconstruction uses the native Redheffer star in every composition mode.
    /// Cancellation or refusal publishes no partial stack.
    pub fn try_new(
        layers: &'l [PorousLayer<'m>],
        options: ComposeOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(&mut cancelled, 0)?;
        options.validate()?;
        let first = layers.first().ok_or(DuctError::EmptyDuct)?;
        if layers.len() > options.maximum_blocks {
            return Err(bad("porous stack composition block work cap exceeded"));
        }
        let n = layers.len();
        let reference = first.two_port().reference_impedance;
        let area = first.area_m2();
        let frequency = first.material_response().angular_frequency_rad_s();
        let mut depths_m = Vec::with_capacity(n + 1);
        depths_m.push(0.0);
        for (i, layer) in layers.iter().enumerate() {
            poll(&mut cancelled, i)?;
            if layer.area_m2() != area
                || layer.two_port().reference_impedance != reference
                || layer.material_response().angular_frequency_rad_s() != frequency
            {
                return Err(bad(
                    "porous stack needs one full area, frequency and real wave reference",
                ));
            }
            let depth = depths_m[i] + layer.thickness_m();
            if !depth.is_finite() || depth <= depths_m[i] {
                return Err(bad("porous stack physical depth is unrepresentable"));
            }
            depths_m.push(depth);
        }
        // At a cut, a=P21*A+P22*b and b=Q11*a+Q12*B. P and Q
        // are the native prefix/suffix blocks, A/B the external incoming waves.
        // The original feedback owner admits 1-P22*Q11. No boundary impedance,
        // inlet-flow normalization or inverse attenuated transmission is used.
        let through = TwoPort::through(reference)?;
        let mut suffix = vec![through; n + 1];
        let mut condition = 1.0_f64;
        for i in (0..n).rev() {
            poll(&mut cancelled, n + n - i - 1)?;
            let (composed, local_condition) = layers[i]
                .two_port()
                .star(suffix[i + 1], options.maximum_condition)?;
            suffix[i] = composed;
            condition = condition.max(local_condition);
        }
        let mut incident_maps = vec![[C64::ZERO; 4]; n];
        let mut prefix = through;
        for i in 0..=n {
            poll(&mut cancelled, 2 * n + i)?;
            let p = prefix.scattering;
            let q = suffix[i].scattering;
            let (denominator, local_condition) = feedback(p[3], q[0]);
            condition = condition.max(local_condition);
            if !condition.is_finite() || condition > options.maximum_condition {
                return Err(DuctError::Singular);
            }
            let from_left = p[2] / denominator;
            let from_right = q[1] / denominator;
            let forward = [from_left, p[3] * from_right];
            let backward = [q[0] * from_left, from_right];
            if !forward.into_iter().chain(backward).all(finite) {
                return Err(DuctError::Singular);
            }
            if i > 0 {
                incident_maps[i - 1][2..].copy_from_slice(&backward);
            }
            if i < n {
                incident_maps[i][..2].copy_from_slice(&forward);
                let (composed, local_condition) =
                    prefix.star(layers[i].two_port(), options.maximum_condition)?;
                prefix = composed;
                condition = condition.max(local_condition);
            }
        }
        poll(&mut cancelled, 3 * n + 1)?;
        Ok(Self {
            layers,
            depths_m,
            incident_maps,
            two_port: prefix,
            condition_estimate: condition,
        })
    }

    /// Original physical materials/layers, left incident face first.
    pub const fn layers(&self) -> &'l [PorousLayer<'m>] {
        self.layers
    }
    /// Global face depths [m], including zero and the complete thickness.
    pub fn layer_depths_m(&self) -> &[f64] {
        &self.depths_m
    }
    /// Complete normal thickness [m].
    pub fn thickness_m(&self) -> f64 {
        self.depths_m[self.layers.len()]
    }
    /// Common full physical face area [m²], never the open-pore area.
    pub fn area_m2(&self) -> f64 {
        self.layers[0].area_m2()
    }
    /// Common angular frequency of the coherent drive [rad/s].
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.layers[0].material_response().angular_frequency_rad_s()
    }
    /// Native composed exterior pressure/flow scattering block.
    pub const fn two_port(&self) -> TwoPort {
        self.two_port
    }
    /// Largest admitted native composition/internal-feedback condition estimate.
    pub const fn condition_estimate(&self) -> f64 {
        self.condition_estimate
    }

    /// Prepare arbitrary coherent peak pressure waves at BOTH exterior faces.
    /// Local fields combine both source amplitudes before any power is squared.
    /// Every original per-layer constitutive/field refusal remains authoritative.
    pub fn prepare_excitation(
        &self,
        incident_pressure_pa: [C64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PorousStackExcitation<'_, 'l, 'm>, DuctError> {
        poll(&mut cancelled, 0)?;
        let boundary = power_for_port(self.two_port, incident_pressure_pa)?;
        let mut layers = Vec::with_capacity(self.layers.len());
        let mut sums = [Dd::ZERO; 3];
        for (i, (layer, map)) in self.layers.iter().zip(&self.incident_maps).enumerate() {
            poll(&mut cancelled, i)?;
            let local_incident = [
                map[0] * incident_pressure_pa[0] + map[1] * incident_pressure_pa[1],
                map[2] * incident_pressure_pa[0] + map[3] * incident_pressure_pa[1],
            ];
            let excitation = layer.prepare_excitation(local_incident)?;
            let work = excitation.dissipation_between(0.0, layer.thickness_m())?;
            for (sum, value) in sums.iter_mut().zip([
                work.viscous_power_w,
                work.thermal_power_w,
                work.absorbed_power_w,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
            layers.push(excitation);
        }
        let [viscous_power_w, thermal_power_w, absorbed_power_w] = sums.map(Dd::to_f64);
        let power_residual_w = boundary.absorbed_power_w - absorbed_power_w;
        let scale = boundary
            .incident_power_w
            .max(boundary.outgoing_power_w)
            .max(absorbed_power_w);
        if ![
            viscous_power_w,
            thermal_power_w,
            absorbed_power_w,
            power_residual_w,
        ]
        .into_iter()
        .all(f64::is_finite)
            || power_residual_w.abs()
                > 512.0 * f64::EPSILON * (self.layers.len() + 1) as f64 * scale
        {
            return Err(bad(
                "porous stack physical cell heat does not close exterior power",
            ));
        }
        poll(&mut cancelled, self.layers.len())?;
        Ok(PorousStackExcitation {
            stack: self,
            layers,
            boundary,
            dissipation: PorousLayerDissipation {
                viscous_power_w,
                thermal_power_w,
                absorbed_power_w,
            },
            power_residual_w,
        })
    }
}

impl<'s, 'l, 'm> PorousStackExcitation<'s, 'l, 'm> {
    /// Original source-owned physical stack.
    pub const fn stack(&self) -> &'s PorousLayerStack<'l, 'm> {
        self.stack
    }
    /// Existing constant-work local field/cell APIs, in original physical order.
    pub fn layer_excitations(&self) -> &[PorousLayerExcitation<'l, 'm>] {
        &self.layers
    }
    /// External incoming/outgoing acoustic pressure, flow and power only.
    /// Internal interface powers are never counted as additional source power.
    pub const fn boundary_power(&self) -> &PorousLayerPower {
        &self.boundary
    }
    /// Sum of original analytic viscous/thermal work across every layer [W].
    pub const fn dissipation(&self) -> PorousLayerDissipation {
        self.dissipation
    }
    /// External incoming minus outgoing minus physical layer heat [W].
    pub const fn power_residual_w(&self) -> f64 {
        self.power_residual_w
    }

    fn admit_depth(&self, depth_m: f64) -> Result<(), DuctError> {
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.stack.thickness_m() {
            Err(bad("porous stack depth must lie inside the physical stack"))
        } else {
            Ok(())
        }
    }
    fn layer_index(&self, depth_m: f64) -> usize {
        self.stack
            .depths_m
            .partition_point(|&face| face <= depth_m)
            .saturating_sub(1)
            .min(self.layers.len() - 1)
    }
    fn local_depth(&self, index: usize, depth_m: f64) -> f64 {
        if depth_m == self.stack.depths_m[index + 1] {
            self.stack.layers[index].thickness_m()
        } else {
            // Cumulative global face coordinates can differ from a local
            // thickness by an addition-rounding ulp. Never extrapolate a layer.
            (depth_m - self.stack.depths_m[index]).min(self.stack.layers[index].thickness_m())
        }
    }
    /// Global depth query, O(log layers), with no allocation. At an internal
    /// interface choose the right-hand material's local loss density; pressure
    /// and superficial flow are continuous. The final face uses the last layer.
    pub fn at_depth(&self, depth_m: f64) -> Result<PorousStackField, DuctError> {
        self.admit_depth(depth_m)?;
        let layer_index = self.layer_index(depth_m);
        let local_depth_m = self.local_depth(layer_index, depth_m);
        Ok(PorousStackField {
            layer_index,
            local_depth_m,
            field: self.layers[layer_index].at_depth(local_depth_m)?,
        })
    }
    /// Exact physical heat in a global interval, allocating nothing and taking
    /// O(intersected layers) work. Every interface splits the interval onto its
    /// own original analytic loss owner; no heat is assigned to an interface.
    pub fn dissipation_between(
        &self,
        start_m: f64,
        end_m: f64,
    ) -> Result<PorousLayerDissipation, DuctError> {
        self.admit_depth(start_m)?;
        self.admit_depth(end_m)?;
        if end_m < start_m {
            return Err(bad("porous stack heat interval must be ordered"));
        }
        let mut sums = [Dd::ZERO; 3];
        if start_m < end_m {
            let mut index = self.layer_index(start_m);
            let mut start = start_m;
            while start < end_m {
                let end = end_m.min(self.stack.depths_m[index + 1]);
                let work = self.layers[index].dissipation_between(
                    self.local_depth(index, start),
                    self.local_depth(index, end),
                )?;
                for (sum, value) in sums.iter_mut().zip([
                    work.viscous_power_w,
                    work.thermal_power_w,
                    work.absorbed_power_w,
                ]) {
                    *sum = *sum + Dd::from_f64(value);
                }
                start = end;
                index += 1;
            }
        }
        let [viscous_power_w, thermal_power_w, absorbed_power_w] = sums.map(Dd::to_f64);
        if ![viscous_power_w, thermal_power_w, absorbed_power_w]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(bad("porous stack heat integral is unrepresentable"));
        }
        Ok(PorousLayerDissipation {
            viscous_power_w,
            thermal_power_w,
            absorbed_power_w,
        })
    }
}

#[cfg(test)]
mod tests;
