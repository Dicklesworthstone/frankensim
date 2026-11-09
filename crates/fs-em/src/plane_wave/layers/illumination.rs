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
//! Prepared depth intervals restrict this same positive operator to thermal
//! cells, with exact homogeneous-film integration instead of depth quadrature.

use super::*;
use fs_math::dd::Dd;

mod interior;
pub use interior::{LayerPointField, LayerPointStatistics, PreparedLayerPoint};

mod incoherent;
pub use incoherent::IncoherentLayerCascade;

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
    incidence_cosine: f64,
    admittance_s: [f64; 2],
    // Original exterior momentum coefficients, prepared once per angle.
    momentum_normal_over_omega: [f64; 2],
    momentum_normal_difference: f64,
    momentum_parallel_over_omega: f64,
    // Outgoing-face major, incoming-face minor, in physical tangential E units.
    scattering: [[C64; 2]; 2],
    // layer -> even/odd loss mode -> left/right input; each contracted mode
    // is sqrt(W/m²), including the peak-phasor time-average factor 1/2.
    loss_modes: Vec<[[C64; 2]; 2]>,
    // Local left/right incoming reference waves, each with two physical source
    // columns. Retained from the original sweeps for O(1) interior queries.
    interior_inputs: Vec<[[C64; 2]; 2]>,
}

/// Positive dissipative source transfers for one finite depth interval in an
/// original film. Preparation and coherent evaluation are allocation-free and
/// O(1); repeated source/coherence evaluations reuse the same two loss modes.
/// These cycle-mean spectral-line powers do not include a transient reservoir
/// or transfer energy into a thermal state by themselves.
pub struct PreparedLayerDepthInterval<'p, 's, 'm> {
    illumination: &'p PreparedLayerIllumination<'s, 'm>,
    layer_index: usize,
    start_depth_m: f64,
    end_depth_m: f64,
    loss_modes: [[C64; 2]; 2],
}

fn reflection(state: [C64; 2]) -> C64 {
    (state[0] - state[1]) / (state[0] + state[1])
}

pub(super) fn dot(a: [C64; 2], b: [C64; 2]) -> C64 {
    // A dark port cancels products, not merely their rounded f64 values.
    // Preserve each multiplication residual through the coherent sum, before
    // any field, positive loss mode or outgoing power is rounded/squared.
    let product = |a: C64, b: C64| {
        let (ar, ai) = (Dd::from_f64(a.re), Dd::from_f64(a.im));
        let (br, bi) = (Dd::from_f64(b.re), Dd::from_f64(b.im));
        [ar * br - ai * bi, ar * bi + ai * br]
    };
    let products = [product(a[0], b[0]), product(a[1], b[1])];
    C64::new(
        (products[0][0] + products[1][0]).to_f64(),
        (products[0][1] + products[1][1]).to_f64(),
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
            admittance(self.incident_epsilon_f_m, self.incident_mu_h_m, c)?,
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
        let inverse_speed = det::sqrt(self.incident_epsilon_f_m) * det::sqrt(self.incident_mu_h_m);
        let momentum_normal_over_omega = [inverse_speed * c, inverse_speed * q.re];
        let momentum_normal_difference = super::super::momentum::normal_difference(
            inverse_speed,
            c,
            q,
            self.exit_wave.epsilon_ratio.re,
            self.exit_wave.mu_ratio,
        );
        let momentum_parallel_over_omega = inverse_speed * det::sqrt((1.0 - c) * (1.0 + c));
        let prepared = PreparedLayerIllumination {
            stack: self,
            polarization,
            incidence_cosine: c,
            admittance_s,
            momentum_normal_over_omega,
            momentum_normal_difference,
            momentum_parallel_over_omega,
            scattering,
            loss_modes,
            interior_inputs: inputs,
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

    /// Prepare exact dielectric heat in a depth interval of one original film.
    /// Depths are in metres from its left face, with
    /// 0 <= start <= end <= thickness. Empty intervals (including a zero-
    /// thickness identity) have exactly zero heat. The whole-film interval
    /// reuses the original loss transfer without recomputing it.
    ///
    /// Boundary traces from the existing interior solution feed the same
    /// positive even/odd slab dissipation modes over the selected width. This
    /// retains standing-wave interference, longitudinal TM electric loss and
    /// weak/opaque limits without subtracting boundary powers or sampling a
    /// depth mesh. No source law is reevaluated and the stack is not solved
    /// again. A thermal consumer applies area/time and transfers the heat once.
    pub fn prepare_depth_interval(
        &self,
        layer_index: usize,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<PreparedLayerDepthInterval<'_, 's, 'm>, ReflectionError> {
        let layer = self
            .stack
            .layers
            .get(layer_index)
            .ok_or_else(|| bad("depth-interval film index is out of range"))?;
        let width =
            fields::admit_depth_interval(layer.source.thickness_m, start_depth_m, end_depth_m)?;
        let loss_modes = if width == 0.0 || layer.wave.epsilon_f_m.im == 0.0 {
            [[C64::ZERO; 2]; 2]
        } else if start_depth_m == 0.0 && end_depth_m == layer.source.thickness_m {
            self.loss_modes[layer_index]
        } else {
            let slab = Slab::new_with_phase_thickness(
                layer.wave,
                self.stack.wave_number,
                self.incidence_cosine,
                self.polarization,
                width,
                Dd::from_f64(end_depth_m) - Dd::from_f64(start_depth_m),
            )?;
            let left = self
                .prepare_point(layer_index, start_depth_m)?
                .tangential_transfer();
            let right = self
                .prepare_point(layer_index, end_depth_m)?
                .tangential_transfer();
            let mut modes = [[C64::ZERO; 2]; 2];
            for source in 0..2 {
                let transfer = fields::interval_loss_transfer(
                    slab,
                    self.admittance_s[0],
                    [left[0][source], left[1][source]],
                    [right[0][source], right[1][source]],
                )?;
                modes[0][source] = transfer[0];
                modes[1][source] = transfer[1];
            }
            modes
        };
        Ok(PreparedLayerDepthInterval {
            illumination: self,
            layer_index,
            start_depth_m,
            end_depth_m,
            loss_modes,
        })
    }

    /// Cycle-mean stack/interface Maxwell stress resultant [Pa] under simultaneous
    /// coherent illumination of both faces. The original x/z incidence plane
    /// has k_x >= 0; both sources share this tangential wavevector while their
    /// normal wavevectors have opposite signs. Source phase can change the
    /// force, including its normal sign, even at fixed incoming powers.
    ///
    /// Each face uses its own original lossless exterior Maxwell stress.
    /// Refraction contributes even for a transparent, impedance-matched
    /// interface. The load excludes material striction/pressure, transients,
    /// finite edges and mechanical work feedback. It is not another heat input.
    /// Allocation-free, O(layers), with the original positive layer-loss modes.
    pub fn radiation_pressure(
        &self,
        incoming_tangential_e_v_m: [C64; 2],
    ) -> Result<[f64; 3], ReflectionError> {
        self.stack.admit_refracting_momentum_exterior()?;
        let mut heat = [0.0; MAX_LAYERS];
        let response = self.evaluate(
            incoming_tangential_e_v_m,
            &mut heat[..self.loss_modes.len()],
        )?;
        let mut diagonal = [Dd::ZERO; 2];
        for modes in &self.loss_modes {
            for (s, sum) in diagonal.iter_mut().enumerate() {
                *sum = *sum
                    + Dd::from_f64(modes[0][s].norm_sq())
                    + Dd::from_f64(modes[1][s].norm_sq());
            }
        }
        let normal_pressure = super::super::momentum::two_sided_normal(
            incoming_tangential_e_v_m,
            self.scattering,
            self.admittance_s,
            self.momentum_normal_over_omega,
            self.momentum_normal_difference,
            diagonal.map(Dd::to_f64),
        );
        let pressure = [
            self.momentum_parallel_over_omega * response.power.absorbed_w_m2,
            0.0,
            normal_pressure,
        ];
        if pressure.iter().any(|f| !f.is_finite()) {
            return Err(bad("coherent radiation pressure is unrepresentable"));
        }
        Ok(pressure)
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

impl<'p, 's, 'm> PreparedLayerDepthInterval<'p, 's, 'm> {
    pub fn illumination(&self) -> &'p PreparedLayerIllumination<'s, 'm> {
        self.illumination
    }
    pub fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub fn start_depth_m(&self) -> f64 {
        self.start_depth_m
    }
    pub fn end_depth_m(&self) -> f64 {
        self.end_depth_m
    }

    /// Cycle-mean absorbed power per transverse area [W/m²] for simultaneous
    /// coherent peak tangential E sources [V/m], left then right. Interference
    /// is retained before squaring. The result is nonnegative and contains
    /// only dissipation within this interval of the original film.
    pub fn evaluate(&self, input: [C64; 2]) -> Result<f64, ReflectionError> {
        if input.iter().any(|&z| !finite(z)) {
            return Err(bad("nonfinite depth-interval illumination source"));
        }
        let heat =
            dot(self.loss_modes[0], input).norm_sq() + dot(self.loss_modes[1], input).norm_sq();
        if !heat.is_finite() {
            return Err(bad("depth-interval absorbed power overflow"));
        }
        Ok(heat)
    }

    /// Explicit peak-E second-moment factors C=B B*, each in V/m. Independent
    /// factors add heat, never phasors. This interval counts as one layer in
    /// maximum_layer_mode_pairs. Allocation-free and transactional; empty
    /// factors mean exact silence. These powers are spectral-line averages,
    /// not PSD/Hz or finite-window beating between different frequencies.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[[C64; 2]],
        budget: IlluminationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, ReflectionError> {
        if factors.len() > 4096
            || factors.len() > budget.maximum_modes
            || factors.len() > budget.maximum_layer_mode_pairs
        {
            return Err(bad("depth-interval coherence exceeds work budget"));
        }
        if factors.iter().flatten().any(|&z| !finite(z)) {
            return Err(bad("nonfinite depth-interval illumination factor"));
        }
        if cancelled() {
            return Err(bad("illumination cancelled"));
        }
        let mut heat = Dd::ZERO;
        for &factor in factors {
            if cancelled() {
                return Err(bad("illumination cancelled"));
            }
            heat = heat + Dd::from_f64(self.evaluate(factor)?);
        }
        let heat = heat.to_f64();
        if !heat.is_finite() {
            return Err(bad("depth-interval statistical power overflow"));
        }
        Ok(heat)
    }

    /// Independent face illumination is an explicit diagonal source covariance.
    pub fn evaluate_independent_sources(&self, input: [C64; 2]) -> Result<f64, ReflectionError> {
        self.evaluate_coherence_modes(
            &[[input[0], C64::ZERO], [C64::ZERO, input[1]]],
            IlluminationBudget {
                maximum_modes: 2,
                maximum_layer_mode_pairs: 2,
            },
            || false,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// G0/G3: exact dyadic cancellation, including both complex quadratures.
    #[test]
    fn coherent_products_preserve_a_residual_below_one_product_ulp() {
        let epsilon = 2.0_f64.powi(-27);
        for exponent in [-450, -200, 0, 200, 450] {
            let scale = 2.0_f64.powi(exponent);
            let a = [
                C64::new(1.0 + epsilon, 1.0 + epsilon).scale(scale),
                C64::new(1.0, 1.0).scale(scale),
            ];
            let b = [C64::from_re((1.0 - epsilon) * scale), C64::from_re(-scale)];
            let residual = -2.0_f64.powi(2 * exponent - 54);
            assert_eq!(dot(a, b), C64::new(residual, residual));
            assert_eq!(dot([a[1], a[0]], [b[1], b[0]]), dot(a, b));
            assert_eq!(dot(a, [a[1], -a[0]]), C64::ZERO);
        }
    }
}
