//! Resolved finite elastic coatings backed by an original local support law.
//!
//! The incident solid occupies z<0, the coating 0<z<h, and the support acts at
//! z=h. The same original tensors, peak `exp(-i omega t)` phase and +z traction
//! convention own every field. The coating is lossless: its interference and
//! evanescence transport work, while only the authored backing can absorb it.
//!
//! One row-scaled 12x12 native solve retains three reflected amplitudes, six
//! original coating-basis coefficients and three actual rear velocities.
//! This preserves weak motion at a stiff support without subtracting nearly
//! opposite waves. The parent slab's decaying/critical basis is reused exactly;
//! no exterior evanescent channel is treated as an incoming power source.

use super::super::surface::{
    ElasticSurfaceCondition, ElasticSurfaceImpedanceTensor, ElasticSurfaceResponse,
    IsotropicElasticSurface, PreparedElasticSurface,
};
use super::*;

/// Original incident/coating materials, physical thickness and rear support.
#[derive(Debug)]
pub struct IsotropicElasticCoating<'a> {
    interface: IsotropicElasticInterface<'a>,
    thickness_m: f64,
    backing: ElasticSurfaceCondition,
}

/// One common frequency and tangential slowness, with fixed-size drive maps.
#[derive(Debug)]
pub struct PreparedElasticCoating<'a> {
    incident: &'a ElasticAcousticState,
    coating: &'a ElasticAcousticState,
    thickness_m: f64,
    frequency_hz: f64,
    backing: ElasticSurfaceCondition,
    tangential_slowness_s_m: f64,
    direct: Option<PreparedElasticSurface<'a>>,
    incoming: [ElasticWaveTrace; 3],
    reflected: [ElasticWaveTrace; 3],
    scattering: [[C64; 3]; 3],
    basis: Option<[ModeBasis; 3]>,
    interior: [[C64; 3]; 6],
    rear_velocity: [[C64; 3]; 3],
    impedance: Option<ElasticSurfaceImpedanceTensor>,
}

/// Coherent reflection, resolved coating fields and actual backing motion.
#[derive(Debug)]
pub struct ElasticCoatingResponse<'p, 'a> {
    prepared: &'p PreparedElasticCoating<'a>,
    coefficients: [C64; 6],
    /// Incident-solid waves and traces at z=0. Absorption is the sole positive
    /// work at the rear support; the lossless coating contributes no heat.
    pub reflection: ElasticSurfaceResponse,
    /// Coating traces at z=0 and z=h. An exactly absent coating has no field.
    /// These retain the original wave reconstruction; near a stiff backing,
    /// use the explicit support fields for weak motion and mechanical work.
    pub coating_faces: Option<[ElasticLayerField; 2]>,
    /// Conservative Mandel-strain bound throughout the original coating.
    pub coating_strain_norm_bound: f64,
    /// Actual velocity at z=h, retained by the coupled solve [m/s].
    pub support_velocity_m_s: [C64; 3],
    /// Traction sigma.e_z at z=h, in the same +z frame [Pa].
    pub support_traction_pa: [C64; 3],
}

impl<'a> IsotropicElasticCoating<'a> {
    /// Bind two original isotropic materials and a nonnegative physical thickness.
    /// Zero thickness recovers the original supported half-space exactly.
    pub fn new(
        solid: &'a ElasticAcousticState,
        coating: &'a ElasticAcousticState,
        thickness_m: f64,
        backing: ElasticSurfaceCondition,
    ) -> Result<Self, MaterialError> {
        if !thickness_m.is_finite() || thickness_m < 0.0 {
            return Err(invalid(
                "elastic coating thickness must be finite and nonnegative",
            ));
        }
        Ok(Self {
            interface: IsotropicElasticInterface::new(solid, coating)?,
            thickness_m,
            backing,
        })
    }

    pub fn prepare_incidence(
        &self,
        mode: ElasticMode,
        cosine: f64,
        frequency_hz: f64,
    ) -> Result<PreparedElasticCoating<'a>, MaterialError> {
        self.prepare_incidence_cancellable(mode, cosine, frequency_hz, &mut || false)
    }

    /// Prepare all propagating incident channels; cancellation publishes nothing.
    /// Exact singular guided responses and unrepresentable systems refuse.
    pub fn prepare_incidence_cancellable(
        &self,
        mode: ElasticMode,
        cosine: f64,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedElasticCoating<'a>, MaterialError> {
        poll(cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 || !omega.is_finite() {
            return Err(invalid(
                "elastic coating frequency must be finite and positive",
            ));
        }
        let mut prepared = PreparedElasticCoating {
            incident: self.interface.states[0],
            coating: self.interface.states[1],
            thickness_m: self.thickness_m,
            frequency_hz,
            backing: self.backing,
            tangential_slowness_s_m: 0.0,
            direct: None,
            incoming: [empty_trace(); 3],
            reflected: [empty_trace(); 3],
            scattering: [[C64::ZERO; 3]; 3],
            basis: None,
            interior: [[C64::ZERO; 3]; 6],
            rear_velocity: [[C64::ZERO; 3]; 3],
            impedance: None,
        };
        if self.thickness_m == 0.0 {
            let direct = IsotropicElasticSurface::new(prepared.incident, self.backing)?
                .prepare_incidence_cancellable(mode, cosine, cancelled)?;
            prepared.tangential_slowness_s_m = direct.tangential_slowness_s_m();
            prepared.direct = Some(direct);
            return Ok(prepared);
        }
        // This is the actual front interface's native wave preparation only;
        // no fictitious second half-space boundary problem is solved.
        let traces =
            self.interface
                .prepare_wave_traces(ElasticSide::First, mode, cosine, cancelled)?;
        prepared.tangential_slowness_s_m = traces.tangential_slowness_s_m;
        prepared.incoming = core::array::from_fn(|i| traces.incoming[i]);
        prepared.reflected = core::array::from_fn(|i| traces.outgoing[i]);
        prepared.basis = Some(prepare_basis(
            prepared.coating,
            self.interface.speeds[1],
            self.interface.impedances[1],
            self.interface.speeds[0][mode.index().min(1)],
            cosine,
            prepared.tangential_slowness_s_m,
            omega,
            self.thickness_m,
            cancelled,
        )?);
        prepared.impedance = match self.backing {
            ElasticSurfaceCondition::Impedance(z) => {
                let values = [
                    z.tangential_pa_s_m(),
                    z.tangential_pa_s_m(),
                    z.normal_pa_s_m(),
                ];
                Some(ElasticSurfaceImpedanceTensor::new(core::array::from_fn(
                    |i| core::array::from_fn(|j| if i == j { values[i] } else { C64::ZERO }),
                ))?)
            }
            ElasticSurfaceCondition::TensorImpedance(z) => Some(z),
            _ => None,
        };
        prepared.prepare_boundary(cancelled)?;
        poll(cancelled)?;
        Ok(prepared)
    }
}

impl<'a> PreparedElasticCoating<'a> {
    pub const fn incident_state(&self) -> &'a ElasticAcousticState {
        self.incident
    }
    pub const fn coating_state(&self) -> &'a ElasticAcousticState {
        self.coating
    }
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn backing(&self) -> ElasticSurfaceCondition {
        self.backing
    }
    pub const fn tangential_slowness_s_m(&self) -> f64 {
        self.tangential_slowness_s_m
    }

    pub fn material_state_identities(&self) -> [ContentHash; 2] {
        [
            self.incident.material_state_identity(),
            self.coating.material_state_identity(),
        ]
    }

    pub fn incoming_mode(&self, mode: ElasticMode) -> Option<ElasticWaveTrace> {
        if let Some(direct) = &self.direct {
            return direct.incoming_mode(mode);
        }
        let wave = self.incoming[mode.index()];
        (wave.normal_power_factor_pa_s_m > 0.0).then_some(wave)
    }

    pub fn reflected_mode(&self, mode: ElasticMode) -> ElasticWaveTrace {
        if let Some(direct) = &self.direct {
            return direct.reflected_mode(mode);
        }
        self.reflected[mode.index()]
    }

    pub fn velocity_coefficients(&self, mode: ElasticMode) -> Option<[C64; 3]> {
        if let Some(direct) = &self.direct {
            return direct.velocity_coefficients(mode);
        }
        self.incoming_mode(mode)?;
        Some(core::array::from_fn(|i| self.scattering[i][mode.index()]))
    }

    /// Apply coherent incident modal velocities without allocation or a new solve.
    /// Both original material strain domains and physical power balances apply.
    pub fn drive(
        &self,
        incoming: [C64; 3],
    ) -> Result<ElasticCoatingResponse<'_, 'a>, MaterialError> {
        if let Some(direct) = &self.direct {
            let reflection = direct.drive(incoming)?;
            return Ok(ElasticCoatingResponse {
                prepared: self,
                coefficients: [C64::ZERO; 6],
                support_velocity_m_s: reflection.boundary_velocity_m_s,
                support_traction_pa: reflection.boundary_traction_pa,
                reflection,
                coating_faces: None,
                coating_strain_norm_bound: 0.0,
            });
        }
        for (i, &value) in incoming.iter().enumerate() {
            if !finite(value)
                || (value != C64::ZERO && self.incoming[i].normal_power_factor_pa_s_m == 0.0)
            {
                return Err(invalid(
                    "elastic coating drive needs finite propagating amplitudes",
                ));
            }
        }
        let reflected: [C64; 3] = core::array::from_fn(|i| {
            complex_sum((0..3).map(|j| self.scattering[i][j] * incoming[j]))
        });
        let coefficients: [C64; 6] = core::array::from_fn(|i| {
            complex_sum((0..3).map(|j| self.interior[i][j] * incoming[j]))
        });
        if reflected.iter().chain(&coefficients).any(|&v| !finite(v)) {
            return Err(invalid("elastic coating amplitudes are unrepresentable"));
        }
        let mut input_power = [0.0; 3];
        let mut reflected_power = [0.0; 3];
        let mut strain = Dd::ZERO;
        for mode in 0..3 {
            input_power[mode] = power(self.incoming[mode], incoming[mode])?;
            reflected_power[mode] = power(self.reflected[mode], reflected[mode])?;
            for (wave, amplitude) in [
                (self.incoming[mode], incoming[mode]),
                (self.reflected[mode], reflected[mode]),
            ] {
                strain =
                    strain + Dd::from_f64(FieldTrace::wave(wave).strain_norm() * amplitude.abs());
            }
        }
        let strain_norm_bound = strain.to_f64();
        let basis = self
            .basis
            .expect("a positive coating has its original basis");
        let bounds = basis.map(|b| b.strain_bounds(self.thickness_m));
        let coating_strain_norm_bound =
            real_sum((0..6).map(|i| bounds[i / 2][i % 2] * coefficients[i].abs()));
        if !strain_norm_bound.is_finite()
            || strain_norm_bound > self.incident.maximum_absolute_mandel_strain()
            || !coating_strain_norm_bound.is_finite()
            || coating_strain_norm_bound > self.coating.maximum_absolute_mandel_strain()
        {
            return Err(invalid(
                "elastic coating drive exceeds an original strain domain",
            ));
        }
        let front: [C64; 6] = core::array::from_fn(|i| {
            complex_sum((0..3).flat_map(|m| {
                [
                    FieldTrace::wave(self.incoming[m]).0[i] * incoming[m],
                    FieldTrace::wave(self.reflected[m]).0[i] * reflected[m],
                ]
            }))
        });
        let boundary_velocity_m_s: [C64; 3] = core::array::from_fn(|i| front[i]);
        let boundary_traction_pa: [C64; 3] = core::array::from_fn(|i| front[3 + i]);
        let boundary_normal_flux_w_m2 = -0.5
            * real_sum(
                (0..3).map(|i| (boundary_velocity_m_s[i].conj() * boundary_traction_pa[i]).re),
            );
        let faces = [
            layer_field(basis, coefficients, 0.0, self.thickness_m)?,
            layer_field(basis, coefficients, self.thickness_m, self.thickness_m)?,
        ];
        let mut support_velocity_m_s = core::array::from_fn(|i| {
            complex_sum((0..3).map(|j| self.rear_velocity[i][j] * incoming[j]))
        });
        let (support_traction_pa, absorbed) = if let Some(z) = self.impedance {
            (
                z.traction(support_velocity_m_s),
                z.absorbed_power(support_velocity_m_s)?,
            )
        } else if self.backing == ElasticSurfaceCondition::Clamped {
            support_velocity_m_s = [C64::ZERO; 3];
            (faces[1].traction_pa, 0.0)
        } else {
            ([C64::ZERO; 3], 0.0)
        };
        let total_incoming = real_sum(input_power);
        let total_reflected = real_sum(reflected_power);
        if front
            .iter()
            .chain(&support_velocity_m_s)
            .chain(&support_traction_pa)
            .any(|&v| !finite(v))
            || ![
                total_incoming,
                total_reflected,
                absorbed,
                boundary_normal_flux_w_m2,
            ]
            .into_iter()
            .all(f64::is_finite)
            || !balanced(total_incoming, total_reflected + absorbed, total_incoming)
            || ![
                boundary_normal_flux_w_m2,
                faces[0].normal_flux_w_m2,
                faces[1].normal_flux_w_m2,
            ]
            .into_iter()
            .all(|flux| balanced(flux, absorbed, total_incoming))
        {
            return Err(invalid(
                "elastic coating power is unrepresentable or does not close",
            ));
        }
        Ok(ElasticCoatingResponse {
            prepared: self,
            coefficients,
            reflection: ElasticSurfaceResponse {
                material_state_identity: self.incident.material_state_identity(),
                incoming_velocity_m_s: incoming,
                reflected_velocity_m_s: reflected,
                incoming_normal_power_w_m2: input_power,
                reflected_normal_power_w_m2: reflected_power,
                total_incoming_power_w_m2: total_incoming,
                total_reflected_power_w_m2: total_reflected,
                absorbed_normal_power_w_m2: absorbed,
                boundary_velocity_m_s,
                boundary_traction_pa,
                boundary_normal_flux_w_m2,
                strain_norm_bound,
            },
            coating_faces: Some(faces),
            coating_strain_norm_bound,
            support_velocity_m_s,
            support_traction_pa,
        })
    }

    fn prepare_boundary(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MaterialError> {
        let basis = self
            .basis
            .expect("a positive coating has its original basis");
        let faces = [
            basis.map(|b| b.at(0.0, self.thickness_m)),
            basis.map(|b| b.at(self.thickness_m, self.thickness_m)),
        ];
        let mut matrix = [C64::ZERO; 144];
        for row in 0..6 {
            poll(cancelled)?;
            for mode in 0..3 {
                matrix[12 * row + mode] = FieldTrace::wave(self.reflected[mode]).0[row];
            }
            for j in 0..6 {
                matrix[12 * row + 3 + j] = -faces[0][j / 2][j % 2].0[row];
            }
        }
        for axis in 0..3 {
            poll(cancelled)?;
            for j in 0..6 {
                matrix[12 * (6 + axis) + 3 + j] = faces[1][j / 2][j % 2].0[axis];
            }
            matrix[12 * (6 + axis) + 9 + axis] = -C64::ONE;
            if self.backing == ElasticSurfaceCondition::Clamped {
                matrix[12 * (9 + axis) + 9 + axis] = C64::ONE;
            } else {
                for j in 0..6 {
                    matrix[12 * (9 + axis) + 3 + j] = faces[1][j / 2][j % 2].0[3 + axis];
                }
                if let Some(z) = self.impedance {
                    for j in 0..3 {
                        matrix[12 * (9 + axis) + 9 + j] = z.matrix_pa_s_m()[axis][j];
                    }
                }
            }
        }
        let mut scales = [0.0_f64; 12];
        for row in 0..12 {
            let entries = &matrix[12 * row..12 * (row + 1)];
            scales[row] = entries
                .iter()
                .fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
            if scales[row] == 0.0 || entries.iter().any(|&z| !finite(z)) {
                return Err(invalid("elastic coating boundary row is unrepresentable"));
            }
            for j in 0..12 {
                matrix[12 * row + j] = divided(matrix[12 * row + j], scales[row]);
            }
        }
        let factor = lu_complex_with_cancel(&matrix, 12, &mut *cancelled).map_err(solve_error)?;
        for channel in 0..3 {
            if self.incoming[channel].normal_power_factor_pa_s_m == 0.0 {
                continue;
            }
            let trace = FieldTrace::wave(self.incoming[channel]);
            let rhs: [C64; 12] = core::array::from_fn(|row| {
                if row < 6 {
                    -divided(trace.0[row], scales[row])
                } else {
                    C64::ZERO
                }
            });
            let mut solution = rhs;
            factor
                .solve_with_cancel(&mut solution, &mut *cancelled)
                .map_err(solve_error)?;
            // Reuse the same factor for at most two corrections, as the original
            // slab does for weak critical traction rows; never relax that row.
            for refinement in 0..=2 {
                let mut correction = [C64::ZERO; 12];
                let mut converged = solution.iter().all(|&v| finite(v));
                for row in 0..12 {
                    let terms: [C64; 12] =
                        core::array::from_fn(|j| matrix[12 * row + j] * solution[j]);
                    correction[row] = rhs[row] - complex_sum(terms);
                    let residual = correction[row].abs();
                    let scale = real_sum(terms.map(C64::abs)) + rhs[row].abs();
                    if !residual.is_finite() || (scale > 0.0 && residual / scale > 1e-11) {
                        converged = false;
                    }
                }
                if converged {
                    break;
                }
                if refinement == 2 {
                    return Err(invalid("elastic coating boundary solve did not converge"));
                }
                factor
                    .solve_with_cancel(&mut correction, &mut *cancelled)
                    .map_err(solve_error)?;
                for i in 0..12 {
                    solution[i] = complex_sum([solution[i], correction[i]]);
                }
            }
            for i in 0..3 {
                self.scattering[i][channel] = solution[i];
                self.rear_velocity[i][channel] = solution[9 + i];
            }
            for i in 0..6 {
                self.interior[i][channel] = solution[3 + i];
            }
        }
        Ok(())
    }
}

impl ElasticCoatingResponse<'_, '_> {
    /// Sample the actual coherent coating field at physical depth z in [0,h].
    /// The unmodified wave basis supplies strain and traction everywhere. At
    /// a stiff endpoint, the explicit support velocity retains better relative
    /// accuracy than a sum of nearly cancelling interior waves.
    pub fn field_at(&self, depth_m: f64) -> Result<ElasticLayerField, MaterialError> {
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.prepared.thickness_m {
            return Err(invalid(
                "elastic coating field depth must lie inside the original coating",
            ));
        }
        let basis = self
            .prepared
            .basis
            .ok_or_else(|| invalid("an absent elastic coating has no interior field"))?;
        layer_field(basis, self.coefficients, depth_m, self.prepared.thickness_m)
    }
}

#[cfg(test)]
mod tests;
