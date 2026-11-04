//! Coherent P/SV/SH fields in a finite, perfectly bonded isotropic solid slab.
//!
//! The first face is z=0 and the second z=h in the parent's common frame.
//! Incident amplitudes are referenced to their own face, with the SAME peak
//! phasor convention `exp(i omega (p x + q z - t))`. All six velocity/traction
//! continuities are imposed at each face using the original elastic tensors.
//! P/SV share one 8x8 factorization; SH uses one 4x4 factorization. Repeated
//! coherent drives and interior queries allocate nothing.
//!
//! Ordinary interior waves are anchored at opposite faces: all evanescent
//! exponentials decay into the slab. Near a critical internal mode, a regular
//! even/odd cosine/sinc basis retains the second solution at q=0. This avoids
//! both growing transfer matrices and a singular pair of coincident waves.
//! The basis change is algebraic, not a physical cutoff or thin-plate model.
//! Interior power is reconstructed from the TOTAL velocity and traction;
//! interfering evanescent waves can transmit power without any absorption.
//!
//! This extends the parent's bonded solid boundary equations. The distinction
//! between evanescence and finite-plate transmission is also discussed for
//! fluid-loaded plates by Almeida et al., <https://arxiv.org/abs/2302.08826>.
//! That fluid/solid boundary law remains owned by `porous::layers::elastic`.
//! No damping, prestress, anisotropy, compliant contact or finite lateral edge
//! is inferred here. Exact undriven guided-mode singularities may refuse.

use super::*;

/// Three original material states and a nonnegative slab thickness [m].
#[derive(Debug)]
pub struct IsotropicElasticLayer<'a> {
    exterior: IsotropicElasticInterface<'a>,
    layer: &'a ElasticAcousticState,
    speeds: [f64; 2],
    impedances: [f64; 2],
    thickness_m: f64,
}

/// One frequency and conserved tangential slowness, with all propagating
/// external incident channels prepared. Exterior incoming evanescence is
/// refused just as for the original halfspace interface.
#[derive(Debug)]
pub struct PreparedElasticLayer<'a> {
    exterior: PreparedElasticInterface<'a>,
    layer: &'a ElasticAcousticState,
    thickness_m: f64,
    frequency_hz: f64,
    basis: Option<[ModeBasis; 3]>,
    // Two interior coefficients per mode, each linear in six incident waves.
    interior: [[C64; 6]; 6],
}

/// Actual coherent field at one depth in the layer, with common +z traction.
#[derive(Clone, Copy, Debug)]
pub struct ElasticLayerField {
    pub velocity_m_s: [C64; 3],
    pub traction_pa: [C64; 3],
    pub strain_mandel: [C64; 6],
    /// `-Re(v* . traction)/2` [W/m^2], including all interference terms.
    pub normal_flux_w_m2: f64,
}

/// A finite admitted drive, borrowing the exact prepared layer for interior
/// queries. Neither source histories nor material states are reconstructed.
#[derive(Debug)]
pub struct ElasticLayerResponse<'p, 'a> {
    prepared: &'p PreparedElasticLayer<'a>,
    coefficients: [C64; 6],
    /// Exterior fields are evaluated at z=0 and z=h respectively. The same
    /// original halfspace power and strain checks own these quantities.
    pub exterior: ElasticInterfaceResponse,
    /// Interior fields at the two faces. An exactly absent layer has no field.
    pub layer_faces: Option<[ElasticLayerField; 2]>,
    /// Conservative Mandel norm bound over the ENTIRE slab, not only its faces.
    pub layer_strain_norm_bound: f64,
}

#[derive(Clone, Copy, Debug)]
struct FieldTrace([C64; 12]);

impl FieldTrace {
    fn wave(wave: ElasticWaveTrace) -> Self {
        Self(core::array::from_fn(|i| match i {
            0..3 => wave.velocity[i],
            3..6 => wave.traction_pa_s_m[i - 3],
            _ => wave.strain_mandel_s_m[i - 6],
        }))
    }

    fn combine(first: Self, a: C64, second: Self, b: C64) -> Self {
        Self(core::array::from_fn(|i| first.0[i] * a + second.0[i] * b))
    }

    fn strain_norm(self) -> f64 {
        self.0[6..]
            .iter()
            .fold(0.0, |norm, x| det::hypot(norm, x.abs()))
    }
}

#[derive(Clone, Copy, Debug)]
struct ModeBasis {
    waves: [ElasticWaveTrace; 2],
    cosine: C64,
    k0: f64,
    regularized: bool,
    // R contains the even-q polynomial; L is the opposite parity divided by
    // the normal cosine. For SV the parity roles are interchanged. L is
    // obtained from the original tensor at unit normal cosine, so it remains
    // finite at EXACT q=0 without differencing nearby incidence angles.
    regular: [FieldTrace; 2],
}

impl<'a> IsotropicElasticLayer<'a> {
    pub fn new(
        first: &'a ElasticAcousticState,
        layer: &'a ElasticAcousticState,
        second: &'a ElasticAcousticState,
        thickness_m: f64,
    ) -> Result<Self, MaterialError> {
        if !layer.is_explicitly_isotropic() {
            return Err(invalid(
                "elastic layer requires an explicitly isotropic state",
            ));
        }
        if !thickness_m.is_finite() || thickness_m < 0.0 {
            return Err(invalid(
                "elastic layer thickness must be finite and nonnegative",
            ));
        }
        let exterior = IsotropicElasticInterface::new(first, second)?;
        let modes = layer.directional_waves([0.0, 0.0, 1.0])?.modes;
        Ok(Self {
            exterior,
            layer,
            speeds: [modes[2].phase_speed_m_s, modes[0].phase_speed_m_s],
            impedances: [modes[2].impedance_pa_s_m, modes[0].impedance_pa_s_m],
            thickness_m,
        })
    }

    pub const fn layer_state(&self) -> &'a ElasticAcousticState {
        self.layer
    }

    pub fn exterior_state(&self, side: ElasticSide) -> &'a ElasticAcousticState {
        self.exterior.state(side)
    }

    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    pub fn prepare_incidence(
        &self,
        side: ElasticSide,
        mode: ElasticMode,
        cosine: f64,
        frequency_hz: f64,
    ) -> Result<PreparedElasticLayer<'a>, MaterialError> {
        self.prepare_incidence_cancellable(side, mode, cosine, frequency_hz, &mut || false)
    }

    pub fn prepare_incidence_cancellable(
        &self,
        side: ElasticSide,
        mode: ElasticMode,
        cosine: f64,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedElasticLayer<'a>, MaterialError> {
        poll(cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 || !omega.is_finite() {
            return Err(invalid(
                "elastic layer frequency must be finite and positive",
            ));
        }
        if self.thickness_m == 0.0 {
            return Ok(PreparedElasticLayer {
                exterior: self
                    .exterior
                    .prepare_incidence_cancellable(side, mode, cosine, cancelled)?,
                layer: self.layer,
                thickness_m: 0.0,
                frequency_hz,
                basis: None,
                interior: [[C64::ZERO; 6]; 6],
            });
        }
        let exterior = self
            .exterior
            .prepare_wave_traces(side, mode, cosine, cancelled)?;
        let incident_speed = self.exterior.speeds[side.index()][mode.index().min(1)];
        let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine));
        let p = exterior.tangential_slowness_s_m;
        let mut modes = [None; 3];
        for kind in 0..3 {
            poll(cancelled)?;
            let speed = self.speeds[kind.min(1)];
            let normal = normal_ratio(incident_speed, speed, cosine)?;
            let wave_cosine = normal.scale(speed / incident_speed);
            let wave_sine = sine * (speed / incident_speed);
            let normal_slowness = divided(normal, incident_speed);
            let k0 = omega / speed;
            let full_phase = wave_cosine.scale(k0 * self.thickness_m);
            if !k0.is_finite()
                || k0 == 0.0
                || !(k0 * self.thickness_m).is_finite()
                || !finite(wave_cosine)
                || !wave_sine.is_finite()
                || !finite(full_phase)
                || (normal != C64::ZERO && normal_slowness == C64::ZERO)
            {
                return Err(invalid("elastic layer normal phase is unrepresentable"));
            }
            let mut waves = [empty_trace(); 2];
            let mut unit = [empty_trace(); 2];
            for direction in 0..2 {
                let sign = if direction == 0 { 1.0 } else { -1.0 };
                waves[direction] = trace(
                    self.layer,
                    self.impedances,
                    kind,
                    sign,
                    wave_sine,
                    wave_cosine,
                    p,
                    normal_slowness,
                )?;
                // Only the q-linear parity of this algebraic trace is used.
                // No unit-cosine state is substituted for the actual wave.
                unit[direction] = trace(
                    self.layer,
                    self.impedances,
                    kind,
                    sign,
                    wave_sine,
                    C64::ONE,
                    p,
                    C64::from_re(1.0 / speed),
                )?;
            }
            let parity = if kind == 1 { -1.0 } else { 1.0 };
            let regular = [
                FieldTrace::combine(
                    FieldTrace::wave(waves[0]),
                    C64::from_re(0.5),
                    FieldTrace::wave(waves[1]),
                    C64::from_re(0.5 * parity),
                ),
                FieldTrace::combine(
                    FieldTrace::wave(unit[0]),
                    C64::from_re(0.5),
                    FieldTrace::wave(unit[1]),
                    C64::from_re(-0.5 * parity),
                ),
            ];
            modes[kind] = Some(ModeBasis {
                waves,
                cosine: wave_cosine,
                k0,
                regular,
                regularized: full_phase.abs() <= 0.5,
            });
        }
        let mut prepared = PreparedElasticLayer {
            exterior,
            layer: self.layer,
            thickness_m: self.thickness_m,
            frequency_hz,
            basis: Some(modes.map(Option::unwrap)),
            interior: [[C64::ZERO; 6]; 6],
        };
        if self.exterior.states.iter().all(|state| {
            state.density_kg_m3() == self.layer.density_kg_m3()
                && state.stiffness_mandel_pa() == self.layer.stiffness_mandel_pa()
        }) {
            prepared.prepare_homogeneous()?;
        } else {
            prepared.prepare_block(true, cancelled)?;
            prepared.prepare_block(false, cancelled)?;
        }
        poll(cancelled)?;
        Ok(prepared)
    }
}

impl<'a> PreparedElasticLayer<'a> {
    pub fn material_state_identities(&self) -> [ContentHash; 3] {
        [
            self.exterior.states[0].material_state_identity(),
            self.layer.material_state_identity(),
            self.exterior.states[1].material_state_identity(),
        ]
    }

    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn tangential_slowness_s_m(&self) -> f64 {
        self.exterior.tangential_slowness_s_m
    }

    pub fn layer_normal_slownesses_s_m(&self) -> Option<[C64; 3]> {
        self.basis
            .map(|b| b.map(|mode| mode.waves[0].normal_slowness_s_m))
    }

    pub fn incoming_mode(&self, side: ElasticSide, mode: ElasticMode) -> Option<ElasticWaveTrace> {
        self.exterior.incoming_mode(side, mode)
    }

    pub fn outgoing_mode(&self, side: ElasticSide, mode: ElasticMode) -> ElasticWaveTrace {
        self.exterior.outgoing_mode(side, mode)
    }

    /// Infinitesimal outgoing modal velocity coefficients, referenced to the
    /// respective faces. A finite drive additionally admits all three strains.
    pub fn velocity_coefficients(
        &self,
        side: ElasticSide,
        mode: ElasticMode,
    ) -> Option<[[C64; 3]; 2]> {
        self.exterior.velocity_coefficients(side, mode)
    }

    pub fn drive(
        &self,
        incoming: [[C64; 3]; 2],
    ) -> Result<ElasticLayerResponse<'_, 'a>, MaterialError> {
        let exterior = self.exterior.drive(incoming)?;
        let amplitudes: [C64; 6] = core::array::from_fn(|i| incoming[i / 3][i % 3]);
        let coefficients: [C64; 6] = core::array::from_fn(|i| {
            complex_sum((0..6).map(|j| self.interior[i][j] * amplitudes[j]))
        });
        if coefficients.iter().any(|&c| !finite(c)) {
            return Err(invalid(
                "elastic layer interior amplitudes are unrepresentable",
            ));
        }
        let mut bound = 0.0;
        let mut faces = None;
        if let Some(basis) = self.basis {
            let bounds = basis.map(|mode| mode.strain_bounds(self.thickness_m));
            bound = real_sum((0..6).map(|i| bounds[i / 2][i % 2] * coefficients[i].abs()));
            if !bound.is_finite() || bound > self.layer.maximum_absolute_mandel_strain() {
                return Err(invalid(
                    "elastic layer drive exceeds original interior strain domain",
                ));
            }
            let fields = [
                self.field_at(coefficients, 0.0)?,
                self.field_at(coefficients, self.thickness_m)?,
            ];
            for side in 0..2 {
                if !balanced(
                    fields[side].normal_flux_w_m2,
                    exterior.boundary_normal_flux_w_m2[side],
                    exterior.total_incoming_power_w_m2,
                ) {
                    return Err(invalid("elastic layer interior power does not close"));
                }
            }
            faces = Some(fields);
        }
        Ok(ElasticLayerResponse {
            prepared: self,
            coefficients,
            exterior,
            layer_faces: faces,
            layer_strain_norm_bound: bound,
        })
    }

    fn field_at(
        &self,
        coefficients: [C64; 6],
        depth_m: f64,
    ) -> Result<ElasticLayerField, MaterialError> {
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.thickness_m {
            return Err(invalid(
                "elastic field depth must lie inside the original slab",
            ));
        }
        let basis = self
            .basis
            .ok_or_else(|| invalid("an absent elastic layer has no interior field"))?;
        let traces = basis.map(|mode| mode.at(depth_m, self.thickness_m));
        let values: [C64; 12] = core::array::from_fn(|i| {
            complex_sum((0..6).map(|j| traces[j / 2][j % 2].0[i] * coefficients[j]))
        });
        let velocity_m_s = core::array::from_fn(|i| values[i]);
        let traction_pa = core::array::from_fn(|i| values[3 + i]);
        let normal_flux_w_m2 =
            -0.5 * real_sum((0..3).map(|i| (velocity_m_s[i].conj() * traction_pa[i]).re));
        if values.iter().any(|&v| !finite(v)) || !normal_flux_w_m2.is_finite() {
            return Err(invalid("elastic layer field is unrepresentable"));
        }
        Ok(ElasticLayerField {
            velocity_m_s,
            traction_pa,
            strain_mandel: core::array::from_fn(|i| values[6 + i]),
            normal_flux_w_m2,
        })
    }

    fn prepare_homogeneous(&mut self) -> Result<(), MaterialError> {
        let basis = self.basis.expect("a positive layer has a basis");
        for side in 0..2 {
            for mode in 0..3 {
                let channel = 3 * side + mode;
                if self.exterior.incoming[channel].normal_power_factor_pa_s_m == 0.0 {
                    continue;
                }
                let b = basis[mode];
                self.exterior.scattering[3 * (1 - side) + mode][channel] =
                    b.phase(self.thickness_m);
                if b.regularized {
                    let phase = b.phase(0.5 * self.thickness_m);
                    self.interior[2 * mode][channel] =
                        phase.scale(if side == 1 && mode == 1 { -1.0 } else { 1.0 });
                    self.interior[2 * mode + 1][channel] =
                        (phase * b.cosine).scale(if side == 1 && mode != 1 { -1.0 } else { 1.0 });
                } else {
                    self.interior[2 * mode + side][channel] = C64::ONE;
                }
            }
        }
        Ok(())
    }

    fn prepare_block(
        &mut self,
        psv: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MaterialError> {
        let mode_indices: &[usize] = if psv { &[0, 1] } else { &[2] };
        let components: &[usize] = if psv { &[0, 2, 3, 5] } else { &[1, 4] };
        let m = mode_indices.len();
        let n = 4 * m;
        if (0..2).all(|s| {
            mode_indices
                .iter()
                .all(|&mode| self.exterior.incoming[3 * s + mode].normal_power_factor_pa_s_m == 0.0)
        }) {
            return Ok(());
        }
        let basis = self.basis.expect("a positive layer has a basis");
        let faces = [
            basis.map(|b| b.at(0.0, self.thickness_m)),
            basis.map(|b| b.at(self.thickness_m, self.thickness_m)),
        ];
        let mut matrix = [C64::ZERO; 64];
        let mut scales = [0.0_f64; 8];
        for face in 0..2 {
            for (r, &component) in components.iter().enumerate() {
                poll(cancelled)?;
                let row = 2 * m * face + r;
                for (j, &mode) in mode_indices.iter().enumerate() {
                    matrix[n * row + 3 * m * face + j] =
                        FieldTrace::wave(self.exterior.outgoing[3 * face + mode]).0[component]
                            .scale(if face == 0 { 1.0 } else { -1.0 });
                    for parity in 0..2 {
                        matrix[n * row + m + 2 * j + parity] = faces[face][mode][parity].0
                            [component]
                            .scale(if face == 0 { -1.0 } else { 1.0 });
                    }
                }
                scales[row] = matrix[n * row..n * (row + 1)]
                    .iter()
                    .fold(0.0, |scale, x| scale.max(x.abs()));
                if !scales[row].is_finite() || scales[row] == 0.0 {
                    return Err(invalid("elastic layer boundary row is unrepresentable"));
                }
                for col in 0..n {
                    matrix[n * row + col] = divided(matrix[n * row + col], scales[row]);
                }
            }
        }
        let factor =
            lu_complex_with_cancel(&matrix[..n * n], n, &mut *cancelled).map_err(solve_error)?;
        for side in 0..2 {
            for &mode in mode_indices {
                let channel = 3 * side + mode;
                if self.exterior.incoming[channel].normal_power_factor_pa_s_m == 0.0 {
                    continue;
                }
                let incoming = FieldTrace::wave(self.exterior.incoming[channel]);
                let mut rhs = [C64::ZERO; 8];
                for (r, &component) in components.iter().enumerate() {
                    let row = 2 * m * side + r;
                    rhs[row] = divided(incoming.0[component], scales[row]).scale(if side == 0 {
                        -1.0
                    } else {
                        1.0
                    });
                }
                let mut solution = rhs;
                factor
                    .solve_with_cancel(&mut solution[..n], &mut *cancelled)
                    .map_err(solve_error)?;
                if solution[..n].iter().any(|&x| !finite(x)) {
                    return Err(invalid(
                        "elastic layer boundary solution is unrepresentable",
                    ));
                }
                // Small critical tractions can lose relative accuracy in an
                // otherwise backward-stable elimination. Reuse this SAME LU
                // for at most two residual corrections; do not relax that
                // physical row's tolerance or snap its critical mode away.
                for refinement in 0..=2 {
                    let mut correction = [C64::ZERO; 8];
                    let mut converged = true;
                    for row in 0..n {
                        correction[row] = rhs[row]
                            - complex_sum((0..n).map(|j| matrix[n * row + j] * solution[j]));
                        let residual = correction[row].abs();
                        let scale =
                            real_sum((0..n).map(|j| (matrix[n * row + j] * solution[j]).abs()))
                                + rhs[row].abs();
                        if !residual.is_finite() || (scale > 0.0 && residual / scale > 1e-11) {
                            converged = false;
                        }
                    }
                    if converged {
                        break;
                    }
                    if refinement == 2 {
                        return Err(invalid("elastic layer continuity solve did not converge"));
                    }
                    factor
                        .solve_with_cancel(&mut correction[..n], &mut *cancelled)
                        .map_err(solve_error)?;
                    for i in 0..n {
                        solution[i] = complex_sum([solution[i], correction[i]]);
                    }
                }
                for (j, &out_mode) in mode_indices.iter().enumerate() {
                    self.exterior.scattering[out_mode][channel] = solution[j];
                    self.exterior.scattering[3 + out_mode][channel] = solution[3 * m + j];
                    for parity in 0..2 {
                        self.interior[2 * out_mode + parity][channel] =
                            solution[m + 2 * j + parity];
                    }
                }
            }
        }
        Ok(())
    }
}

impl ElasticLayerResponse<'_, '_> {
    /// Field at the physical depth, on the same source phase and material state.
    pub fn field_at(&self, depth_m: f64) -> Result<ElasticLayerField, MaterialError> {
        self.prepared.field_at(self.coefficients, depth_m)
    }
}

impl ModeBasis {
    fn phase(self, distance: f64) -> C64 {
        let phase = self.cosine.scale(self.k0 * distance);
        if phase.im > 0.0 {
            C64::from_re(det::exp(-phase.im))
        } else {
            let (sine, cosine) = det::sin_cos(phase.re);
            C64::new(cosine, sine)
        }
    }

    fn at(self, z: f64, h: f64) -> [FieldTrace; 2] {
        if !self.regularized {
            return [
                FieldTrace::combine(
                    FieldTrace::wave(self.waves[0]),
                    self.phase(z),
                    self.regular[0],
                    C64::ZERO,
                ),
                FieldTrace::combine(
                    FieldTrace::wave(self.waves[1]),
                    self.phase(h - z),
                    self.regular[0],
                    C64::ZERO,
                ),
            ];
        }
        let theta = self.k0 * (z - 0.5 * h);
        let phase = self.cosine.scale(theta);
        let squared = phase.re * phase.re - phase.im * phase.im;
        // The regular basis is selected only for |k_z h| <= 1/2, so these
        // even series have |argument| <= 1/4. The omitted terms are <3e-27.
        let mut cosine = 1.0;
        let mut sinc = 1.0;
        let mut c_term = 1.0;
        let mut s_term = 1.0;
        for j in 1..=8 {
            c_term *= -squared / ((2 * j - 1) * (2 * j)) as f64;
            s_term *= -squared / ((2 * j) * (2 * j + 1)) as f64;
            cosine += c_term;
            sinc += s_term;
        }
        let i = C64::new(0.0, 1.0);
        let sine = phase.scale(sinc);
        [
            FieldTrace::combine(
                self.regular[0],
                C64::from_re(cosine),
                self.regular[1],
                i * self.cosine * sine,
            ),
            FieldTrace::combine(
                self.regular[1],
                C64::from_re(cosine),
                self.regular[0],
                i.scale(theta * sinc),
            ),
        ]
    }

    fn strain_bounds(self, h: f64) -> [f64; 2] {
        if !self.regularized {
            return self.waves.map(|wave| FieldTrace::wave(wave).strain_norm());
        }
        let theta = 0.5 * (self.k0 * h);
        let phase = self.cosine.abs() * theta;
        // cosh(|phase|) and sinh(|phase|)/|phase| bound both trigonometric and
        // hyperbolic cases over the whole symmetric interval. Here phase<=1/4.
        let upper_cosine = (det::exp(phase) + det::exp(-phase)) * 0.5;
        let upper_sinc = if phase == 0.0 {
            1.0
        } else {
            (det::expm1(phase) - det::expm1(-phase)) / (2.0 * phase)
        };
        let normal = self.cosine.abs();
        let r = self.regular[0].strain_norm();
        let l = self.regular[1].strain_norm();
        [
            r * upper_cosine + l * normal * phase * upper_sinc,
            l * upper_cosine + r * theta * upper_sinc,
        ]
    }
}

#[cfg(test)]
mod tests;
