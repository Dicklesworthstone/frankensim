//! Simultaneous coherent illumination of the original finite-film operator.
//!
//! Incoming sources share frequency, transverse wavevector, polarization and
//! the exp(+i omega t) clock. Their peak tangential E amplitudes are referenced
//! to the two stack faces, on the same fixed electric-field axis. Both exterior
//! modes must propagate in lossless media; a lossy/evanescent exterior has no
//! two-channel additive incoming-power decomposition in this model.
//!
//! Two bounded sweeps of the SAME slab operators prepare the complete linear
//! field response. Slab absorption is contracted through the existing positive
//! even/odd dissipative modes, retaining interference before taking powers.
//! No total or layer heat is recovered by subtracting nearly equal fluxes.

use super::*;
use fs_math::dd::Dd;

/// A partial-coherence request is admitted before any factor propagation.
#[derive(Clone, Copy, Debug)]
pub struct IlluminationBudget {
    pub maximum_modes: usize,
    pub maximum_layer_mode_pairs: usize,
}

/// Time/ensemble-averaged normal Poynting power at the two physical faces.
/// All values are W/m². Incoming and outgoing port powers are nonnegative;
/// their signed difference is the net flux through that face.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LayerIlluminationPower {
    pub incoming_w_m2: [f64; 2],
    pub outgoing_w_m2: [f64; 2],
    /// Sum of physical dielectric loss in finite films only.
    pub absorbed_w_m2: f64,
}

/// Peak fields at both faces, including simultaneous-source interference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentLayerIllumination {
    /// Outgoing tangential E [V/m], left then right, on the fixed source axis.
    pub outgoing_tangential_e_v_m: [C64; 2],
    /// Total tangential E [V/m] at the two faces.
    pub boundary_tangential_e_v_m: [C64; 2],
    /// The power-conjugate tangential H [A/m], signed for positive z flux:
    /// -H_x for TE (E_y), H_y for TM (E_x). S_z=Re(E conj(H))/2.
    pub boundary_conjugate_h_a_m: [C64; 2],
    pub power: LayerIlluminationPower,
}

/// Immutable scattering and positive loss transfers, tied to original sources.
/// Preparation costs O(layers); evaluation costs O(layers) and allocates
/// nothing. Per-layer output follows the original authored layer order.
pub struct PreparedLayerIllumination<'s, 'm> {
    stack: &'s PlanarDielectricLayers<'m>,
    polarization: Polarization,
    admittance_s: [f64; 2],
    // Outgoing-face major, incoming-face minor, in physical tangential E units.
    scattering: [[C64; 2]; 2],
    // layer -> even/odd loss mode -> left/right input; each contracted mode
    // is sqrt(W/m²), including the peak-phasor time-average factor 1/2.
    loss_modes: Vec<[[C64; 2]; 2]>,
}

fn reflection(state: [C64; 2]) -> C64 {
    (state[0] - state[1]) / (state[0] + state[1])
}

fn dot(a: [C64; 2], b: [C64; 2]) -> C64 {
    let products = [a[0] * b[0], a[1] * b[1]];
    C64::new(
        (Dd::from_f64(products[0].re) + Dd::from_f64(products[1].re)).to_f64(),
        (Dd::from_f64(products[0].im) + Dd::from_f64(products[1].im)).to_f64(),
    )
}

impl<'m> PlanarDielectricLayers<'m> {
    /// Prepare physical sources at both faces, sharing the declared transverse
    /// wavevector. The original one-sided API still admits lossy/TIR exits;
    /// this two-source API requires two lossless propagating exterior ports.
    pub fn prepare_illumination(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
    ) -> Result<PreparedLayerIllumination<'_, 'm>, ReflectionError> {
        let c = positive(incidence_cosine)?;
        if c > 1.0 || self.exit_wave.epsilon_f_m.im != 0.0 {
            return Err(bad(
                "two-sided illumination requires lossless exterior media and cosine in (0,1]",
            ));
        }
        let q = self.exit_wave.normal(c)?;
        if q.im != 0.0 || q.re <= 0.0 {
            return Err(bad("two-sided illumination needs a propagating exit mode"));
        }
        let exit_cosine = positive(q.re / det::sqrt(self.exit_wave.index_squared_ratio.re))?;
        let admittance = |epsilon: f64, mu: f64, cosine: f64| {
            let normal = positive(det::sqrt(epsilon) / det::sqrt(mu))?;
            positive(match polarization {
                Polarization::Te => normal * cosine,
                Polarization::Tm => normal / cosine,
            })
        };
        let admittance_s = [
            admittance(
                linear_electric_permittivity(self.incident)?,
                linear_magnetic_permeability(self.incident)?,
                c,
            )?,
            admittance(
                self.exit_wave.epsilon_f_m.re,
                linear_magnetic_permeability(self.exit.material())?,
                exit_cosine,
            )?,
        ];
        let n = self.layers.len();
        let mut slabs = [Slab::IDENTITY; MAX_LAYERS];
        for (i, layer) in self.layers.iter().enumerate() {
            slabs[i] = Slab::new(
                layer.wave,
                self.wave_number,
                c,
                polarization,
                layer.source.thickness_m,
            )?;
        }
        let mut right = [[C64::ZERO; 2]; MAX_LAYERS + 1];
        let mut left = [[C64::ZERO; 2]; MAX_LAYERS + 1];
        let mut forward = [C64::ONE; MAX_LAYERS];
        let mut reverse = [C64::ONE; MAX_LAYERS];
        let exit = admittance_state(self.exit_wave, c, polarization)?;
        right[n] = exit;
        left[0] = [C64::ONE; 2];
        for i in (0..n).rev() {
            (right[i], forward[i]) = if self.layers[i].source.thickness_m == 0.0 {
                (right[i + 1], C64::ONE)
            } else {
                propagate_state(slabs[i], right[i + 1])?
            };
        }
        for i in 0..n {
            (left[i + 1], reverse[i]) = if self.layers[i].source.thickness_m == 0.0 {
                (left[i], C64::ONE)
            } else {
                propagate_state(slabs[i], left[i])?
            };
        }
        // exit=[a,b], b/a=Y_right/Y_left. Keep the bounded pair rather than
        // forming that possibly extreme admittance ratio. Match a physical
        // unit E incident from the right to the common left-reference ports.
        let denominator = exit[1] * left[n][0] + exit[0] * left[n][1];
        if denominator == C64::ZERO || !finite(denominator) {
            return Err(bad("right incident field is unrepresentable"));
        }
        let right_input = exit[1] * (left[n][0] + left[n][1]) / denominator;
        let right_reflection = (exit[1] * left[n][0] - exit[0] * left[n][1]) / denominator;
        // Each local slab has incident reference waves a_left and b_right.
        // Columns retain both independently applied physical face sources.
        let mut inputs = vec![[[C64::ZERO; 2]; 2]; n];
        let mut incoming = C64::ONE;
        for i in 0..n {
            let outgoing = incoming * forward[i];
            inputs[i][0][0] = incoming;
            inputs[i][1][0] = reflection(right[i + 1]) * outgoing;
            incoming = outgoing;
        }
        let transmitted_right = (exit[0].scale(2.0) / (exit[0] + exit[1])) * incoming;
        incoming = right_input;
        for i in (0..n).rev() {
            let outgoing = incoming * reverse[i];
            inputs[i][0][1] = reflection(left[i]) * outgoing;
            inputs[i][1][1] = incoming;
            incoming = outgoing;
        }
        let scattering = [
            [reflection(right[0]), incoming],
            [transmitted_right, right_reflection],
        ];
        let root_reference = det::sqrt(admittance_s[0]);
        let mut loss_modes = Vec::with_capacity(n);
        for i in 0..n {
            let even = 0.5 * root_reference * det::sqrt(slabs[i].even_loss);
            let odd = 0.5 * root_reference * det::sqrt(slabs[i].odd_loss);
            loss_modes.push([
                std::array::from_fn(|s| (inputs[i][0][s] + inputs[i][1][s]).scale(even)),
                std::array::from_fn(|s| (inputs[i][0][s] - inputs[i][1][s]).scale(odd)),
            ]);
        }
        if scattering.iter().flatten().any(|&z| !finite(z))
            || loss_modes.iter().flatten().flatten().any(|&z| !finite(z))
        {
            return Err(bad("two-sided field or heat transfer is unrepresentable"));
        }
        let prepared = PreparedLayerIllumination {
            stack: self,
            polarization,
            admittance_s,
            scattering,
            loss_modes,
        };
        // Test the full two-source power form, including both real and imaginary
        // interference. Unit incident powers prevent a large port admittance
        // from hiding a defective weak-channel coefficient.
        let unit = admittance_s.map(|y| 1.0 / det::sqrt(y));
        let mut heat = [0.0; MAX_LAYERS];
        for source in [
            [C64::from_re(unit[0]), C64::ZERO],
            [C64::ZERO, C64::from_re(unit[1])],
            [C64::from_re(unit[0]), C64::from_re(unit[1])],
            [C64::from_re(unit[0]), C64::new(0.0, unit[1])],
        ] {
            prepared.evaluate(source, &mut heat[..n])?;
        }
        Ok(prepared)
    }
}

impl<'s, 'm> PreparedLayerIllumination<'s, 'm> {
    pub fn stack(&self) -> &'s PlanarDielectricLayers<'m> {
        self.stack
    }
    pub fn polarization(&self) -> Polarization {
        self.polarization
    }
    /// The two real physical modal admittances [A/V], including TE/TM angle.
    pub fn port_admittance_s(&self) -> [f64; 2] {
        self.admittance_s
    }
    /// Outgoing-face major, incoming-face minor. Physical tangential E units;
    /// power reciprocity includes the two admittances, not ordinary symmetry.
    pub fn tangential_e_scattering(&self) -> [[C64; 2]; 2] {
        self.scattering
    }

    fn admit_output(&self, n: usize) -> Result<(), ReflectionError> {
        if n != self.loss_modes.len() {
            return Err(bad("one heat output is required per authored film"));
        }
        Ok(())
    }
    fn power(
        &self,
        input: [C64; 2],
        output: [C64; 2],
        heat: &[f64],
    ) -> Result<LayerIlluminationPower, ReflectionError> {
        let flux = |field: [C64; 2]| {
            std::array::from_fn(|i| {
                field[i]
                    .scale(det::sqrt(self.admittance_s[i]) * core::f64::consts::FRAC_1_SQRT_2)
                    .norm_sq()
            })
        };
        let result = LayerIlluminationPower {
            incoming_w_m2: flux(input),
            outgoing_w_m2: flux(output),
            absorbed_w_m2: heat
                .iter()
                .fold(Dd::ZERO, |s, &x| s + Dd::from_f64(x))
                .to_f64(),
        };
        self.check_balance(result)?;
        Ok(result)
    }
    fn check_balance(&self, power: LayerIlluminationPower) -> Result<(), ReflectionError> {
        let incoming = power.incoming_w_m2[0] + power.incoming_w_m2[1];
        let outgoing = power.outgoing_w_m2[0] + power.outgoing_w_m2[1];
        let scale = incoming.max(outgoing).max(power.absorbed_w_m2);
        if [incoming, outgoing, power.absorbed_w_m2]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || (scale > 0.0
                && (incoming / scale - outgoing / scale - power.absorbed_w_m2 / scale).abs()
                    > 2048.0 * f64::EPSILON * (self.loss_modes.len() + 1) as f64)
        {
            return Err(bad("coherent finite-film Poynting balance failed"));
        }
        Ok(())
    }

    /// Apply simultaneous peak tangential E inputs [V/m], left then right.
    /// Heat outputs [W/m²] retain cross-source interference within each film.
    /// This does not debit a thermal state; consumers own that transfer once.
    /// Allocation-free. Discard the heat buffer if evaluation refuses.
    pub fn evaluate(
        &self,
        incoming_tangential_e_v_m: [C64; 2],
        layer_heat_w_m2: &mut [f64],
    ) -> Result<CoherentLayerIllumination, ReflectionError> {
        self.admit_output(layer_heat_w_m2.len())?;
        let input = incoming_tangential_e_v_m;
        if input.iter().any(|&z| !finite(z)) {
            return Err(bad("nonfinite illumination source"));
        }
        let outgoing = self.scattering.map(|row| dot(row, input));
        for (heat, modes) in layer_heat_w_m2.iter_mut().zip(&self.loss_modes) {
            *heat = dot(modes[0], input).norm_sq() + dot(modes[1], input).norm_sq();
            if !heat.is_finite() {
                return Err(bad("film heat overflow"));
            }
        }
        let boundary_e = [input[0] + outgoing[0], outgoing[1] + input[1]];
        let boundary_h = [
            (input[0] - outgoing[0]).scale(self.admittance_s[0]),
            (outgoing[1] - input[1]).scale(self.admittance_s[1]),
        ];
        if outgoing
            .iter()
            .chain(&boundary_e)
            .chain(&boundary_h)
            .any(|&z| !finite(z))
        {
            return Err(bad("illuminated boundary field overflow"));
        }
        Ok(CoherentLayerIllumination {
            outgoing_tangential_e_v_m: outgoing,
            boundary_tangential_e_v_m: boundary_e,
            boundary_conjugate_h_a_m: boundary_h,
            power: self.power(input, outgoing, layer_heat_w_m2)?,
        })
    }

    /// Explicit peak-E second-moment factors C=B B*: each entry is one coherent
    /// two-source factor in V/m. Factors are mutually uncorrelated. Returned
    /// quantities are spectral-line averages, not PSD/Hz. Different frequencies
    /// or transverse wavevectors require separate prepared operators.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[[C64; 2]],
        budget: IlluminationBudget,
        layer_heat_w_m2: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayerIlluminationPower, ReflectionError> {
        self.admit_output(layer_heat_w_m2.len())?;
        if factors.len() > 4096
            || factors.len() > budget.maximum_modes
            || factors
                .len()
                .checked_mul(self.loss_modes.len())
                .is_none_or(|n| n > budget.maximum_layer_mode_pairs)
        {
            return Err(bad("illumination coherence exceeds work budget"));
        }
        if factors.iter().flatten().any(|&z| !finite(z)) {
            return Err(bad("nonfinite illumination factor"));
        }
        if cancelled() {
            return Err(bad("illumination cancelled"));
        }
        let mut sums = [Dd::ZERO; 5];
        let mut heat_sums = [Dd::ZERO; MAX_LAYERS];
        let mut heat = [0.0; MAX_LAYERS];
        for &factor in factors {
            if cancelled() {
                return Err(bad("illumination cancelled"));
            }
            let power = self
                .evaluate(factor, &mut heat[..self.loss_modes.len()])?
                .power;
            for (s, x) in sums.iter_mut().zip([
                power.incoming_w_m2[0],
                power.incoming_w_m2[1],
                power.outgoing_w_m2[0],
                power.outgoing_w_m2[1],
                power.absorbed_w_m2,
            ]) {
                *s = *s + Dd::from_f64(x);
            }
            for (sum, &x) in heat_sums.iter_mut().zip(&heat[..self.loss_modes.len()]) {
                *sum = *sum + Dd::from_f64(x);
            }
        }
        let total = sums.map(Dd::to_f64);
        let power = LayerIlluminationPower {
            incoming_w_m2: [total[0], total[1]],
            outgoing_w_m2: [total[2], total[3]],
            absorbed_w_m2: total[4],
        };
        self.check_balance(power)?;
        for (out, sum) in layer_heat_w_m2.iter_mut().zip(heat_sums) {
            *out = sum.to_f64();
        }
        Ok(power)
    }

    /// Independent face illumination, without inferring that assumption from
    /// unknown phase. Peak input amplitudes still use V/m, not input powers.
    pub fn evaluate_independent_sources(
        &self,
        peak_e_v_m: [C64; 2],
        layer_heat_w_m2: &mut [f64],
    ) -> Result<LayerIlluminationPower, ReflectionError> {
        self.evaluate_coherence_modes(
            &[[peak_e_v_m[0], C64::ZERO], [C64::ZERO, peak_e_v_m[1]]],
            IlluminationBudget {
                maximum_modes: 2,
                maximum_layer_mode_pairs: 2 * MAX_LAYERS,
            },
            layer_heat_w_m2,
            || false,
        )
    }
}
