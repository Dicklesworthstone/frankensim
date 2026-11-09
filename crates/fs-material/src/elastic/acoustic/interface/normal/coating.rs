//! Resolved normal-incidence anisotropic coatings with a passive backing.
//!
//! The incident solid occupies z<0 and the coating 0<z<h. Each original
//! characteristic maps the SAME physical surface axes into its own material
//! axes. Positive normal is from the incident solid toward the support.
//! Phase waves are normal to the plane; their anisotropic energy velocity may
//! walk off laterally. Infinite homogeneous planes are assumed.
//!
//! The original three Christoffel polarizations and their opposite-direction
//! traces propagate through the lossless coating. One native 12x12 solve
//! retains reflection, six interior amplitudes and the actual rear velocity.
//! No isotropic replacement, new eigensolver or transfer-matrix inversion is
//! used. Only the original passive backing absorbs work.
//!
//! Source strain bounds apply over the entire incident half-space and entire
//! coating, not just their sampled faces. This owner adds no oblique branches,
//! bulk attenuation, prestress, finite lateral geometry or measured calibration.

use super::*;
use crate::elastic::acoustic::interface::surface::{
    ElasticSurfaceCondition, ElasticSurfaceImpedanceTensor,
};

/// Two original normal-wave characteristics and a finite physical coating.
/// Their surface axes describe one common physical frame, although each
/// surface-to-material matrix can differ.
#[derive(Debug)]
pub struct NormalElasticCoating<'a> {
    incident: NormalElasticCharacteristic<'a>,
    coating: NormalElasticCharacteristic<'a>,
    thickness_m: f64,
    backing: ElasticSurfaceCondition,
}

/// Fixed-frequency coating response maps. Preparation uses bounded native LU;
/// coherent drives and interior field queries allocate nothing.
#[derive(Debug)]
pub struct PreparedNormalElasticCoating<'a> {
    incident: NormalElasticCharacteristic<'a>,
    coating: NormalElasticCharacteristic<'a>,
    thickness_m: f64,
    frequency_hz: f64,
    backing: ElasticSurfaceCondition,
    impedance: Option<ElasticSurfaceImpedanceTensor>,
    direct: Option<PreparedNormalElasticSurface<'a>>,
    wavenumbers_rad_m: [Dd; 3],
    phase_through: [C64; 3],
    scattering: [[C64; 3]; 3],
    // For each mode: +normal amplitude at z=0, -normal amplitude at z=h.
    interior: [[C64; 3]; 6],
    rear_velocity: [[C64; 3]; 3],
}

/// Coherent physical field at one depth. Cartesian vectors use surface axes;
/// Mandel strain remains in the ORIGINAL COATING MATERIAL axes.
#[derive(Clone, Copy, Debug)]
pub struct NormalElasticCoatingField {
    pub velocity_m_s: [C64; 3],
    /// Traction sigma.n on the common positive normal [Pa].
    pub traction_pa: [C64; 3],
    pub strain_mandel: [C64; 6],
    /// Independently reconstructed -Re(v*.traction)/2 [W/m²]. Near a stiff
    /// backing, cancellation can limit relative accuracy of this small flux.
    /// The response's positive support work remains its authoritative loss.
    pub normal_flux_w_m2: f64,
}

/// Admitted coherent reflection and resolved interior/support fields.
#[derive(Debug)]
pub struct NormalElasticCoatingResponse<'p, 'a> {
    prepared: &'p PreparedNormalElasticCoating<'a>,
    coefficients: [C64; 6],
    /// Native incident quasi-mode amplitudes and original front fields.
    /// Absorption is the sole positive work at the rear support.
    pub reflection: NormalElasticSurfaceResponse,
    /// Original-wave fields at z=0 and z=h; an exactly absent coating has none.
    /// Use the explicit support fields for weak rear motion beside stiff loads.
    pub coating_faces: Option<[NormalElasticCoatingField; 2]>,
    /// Conservative Mandel norm bound over the COMPLETE original coating.
    pub coating_strain_norm_bound: f64,
    /// Actual rear velocity retained as an independent solve unknown [m/s].
    pub support_velocity_m_s: [C64; 3],
    /// Rear traction on the common +normal [Pa].
    pub support_traction_pa: [C64; 3],
}

impl<'a> NormalElasticCoating<'a> {
    /// Bind original incident/coating characteristics, nonnegative thickness
    /// and one passive support in their common physical surface frame.
    /// Zero thickness delegates exactly to the incident supported half-space.
    pub fn new(
        incident: &NormalElasticCharacteristic<'a>,
        coating: &NormalElasticCharacteristic<'a>,
        thickness_m: f64,
        backing: ElasticSurfaceCondition,
    ) -> Result<Self, MaterialError> {
        if !thickness_m.is_finite() || thickness_m < 0.0 {
            return Err(invalid("normal elastic coating thickness must be finite and nonnegative"));
        }
        Ok(Self {
            incident: incident.clone(),
            coating: coating.clone(),
            thickness_m,
            backing,
        })
    }

    pub fn prepare(
        &self,
        frequency_hz: f64,
    ) -> Result<PreparedNormalElasticCoating<'a>, MaterialError> {
        self.prepare_cancellable(frequency_hz, &mut || false)
    }

    /// All three original incident quasi-modes share one frequency and phase
    /// reference. Cancellation or a singular/unresolved solve publishes nothing.
    pub fn prepare_cancellable(
        &self,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedNormalElasticCoating<'a>, MaterialError> {
        poll(cancelled)?;
        let omega = Dd::from_f64(core::f64::consts::TAU) * Dd::from_f64(frequency_hz);
        if !frequency_hz.is_finite() || frequency_hz <= 0.0
            || !omega.hi.is_finite() || !omega.lo.is_finite() || omega.to_f64() <= 0.0
        {
            return Err(invalid("normal elastic coating frequency must be finite and positive"));
        }
        let mut result = PreparedNormalElasticCoating {
            incident: self.incident.clone(),
            coating: self.coating.clone(),
            thickness_m: self.thickness_m,
            frequency_hz,
            backing: self.backing,
            impedance: None,
            direct: None,
            wavenumbers_rad_m: [Dd::ZERO; 3],
            phase_through: [C64::ONE; 3],
            scattering: [[C64::ZERO; 3]; 3],
            interior: [[C64::ZERO; 3]; 6],
            rear_velocity: [[C64::ZERO; 3]; 3],
        };
        if self.thickness_m == 0.0 {
            result.direct = Some(self.incident.prepare_surface_cancellable(self.backing, cancelled)?);
            poll(cancelled)?;
            return Ok(result);
        }
        for m in 0..3 {
            poll(cancelled)?;
            let k = omega / Dd::from_f64(self.coating.waves.modes[m].phase_speed_m_s);
            if !k.hi.is_finite() || !k.lo.is_finite() || k.to_f64() <= 0.0 {
                return Err(invalid("normal elastic coating wavenumber is unrepresentable"));
            }
            result.wavenumbers_rad_m[m] = k;
            result.phase_through[m] = propagation_phase(k, Dd::from_f64(self.thickness_m))?;
        }
        result.impedance = match self.backing {
            ElasticSurfaceCondition::TensorImpedance(z) => Some(z),
            ElasticSurfaceCondition::Impedance(z) => {
                let mut matrix = [[C64::ZERO; 3]; 3];
                matrix[0][0] = z.tangential_pa_s_m();
                matrix[1][1] = z.tangential_pa_s_m();
                matrix[2][2] = z.normal_pa_s_m();
                Some(ElasticSurfaceImpedanceTensor::new(matrix)?)
            }
            ElasticSurfaceCondition::TractionFree | ElasticSurfaceCondition::Clamped => None,
        };
        result.prepare_boundary(cancelled)?;
        poll(cancelled)?;
        Ok(result)
    }
}

impl<'a> PreparedNormalElasticCoating<'a> {
    pub const fn incident_characteristic(&self) -> &NormalElasticCharacteristic<'a> {
        &self.incident
    }

    pub const fn coating_characteristic(&self) -> &NormalElasticCharacteristic<'a> {
        &self.coating
    }

    pub const fn thickness_m(&self) -> f64 { self.thickness_m }
    pub const fn frequency_hz(&self) -> f64 { self.frequency_hz }
    pub const fn backing(&self) -> ElasticSurfaceCondition { self.backing }

    pub fn material_state_identities(&self) -> [ContentHash; 2] {
        [self.incident.state.material_state_identity(), self.coating.state.material_state_identity()]
    }

    /// Reflected native quasi-mode velocities per unit incident native mode.
    pub fn velocity_coefficients(&self, index: usize) -> Option<[C64; 3]> {
        if let Some(direct) = &self.direct { return direct.velocity_coefficients(index); }
        (index < 3).then(|| core::array::from_fn(|i| self.scattering[i][index]))
    }

    /// Incoming physical velocity in the common surface Cartesian axes.
    pub fn drive_cartesian(
        &self,
        incoming: [C64; 3],
    ) -> Result<NormalElasticCoatingResponse<'_, 'a>, MaterialError> {
        self.drive(self.incident.surface_velocity_to_modal(incoming))
    }

    /// Coherent peak velocities in the original incident ascending-speed basis.
    /// No allocation, constitutive source query, eigensolve or boundary solve.
    pub fn drive(
        &self,
        incoming: [C64; 3],
    ) -> Result<NormalElasticCoatingResponse<'_, 'a>, MaterialError> {
        if let Some(direct) = &self.direct {
            let reflection = direct.drive(incoming)?;
            return Ok(NormalElasticCoatingResponse {
                prepared: self,
                coefficients: [C64::ZERO; 6],
                support_velocity_m_s: reflection.boundary_velocity_m_s,
                support_traction_pa: reflection.boundary_traction_pa,
                reflection,
                coating_faces: None,
                coating_strain_norm_bound: 0.0,
            });
        }
        let input_power = self.incident.powers(incoming)?;
        let reflected = core::array::from_fn(|i| {
            complex_dot((0..3).map(|j| (self.scattering[i][j], incoming[j])))
        });
        let coefficients: [C64; 6] = core::array::from_fn(|i| {
            complex_dot((0..3).map(|j| (self.interior[i][j], incoming[j])))
        });
        if coefficients.iter().any(|&v| !finite(v)) {
            return Err(invalid("normal elastic coating amplitudes are unrepresentable"));
        }
        let output_power = self.incident.powers(reflected)?;
        let front = self.incident.field(incoming, reflected)?;
        // All three interior phase waves propagate without loss. Their
        // amplitudes have constant magnitude across the whole layer, giving
        // this single conservative bound even when both faces are strain nodes.
        let coating_strain_norm_bound = real_sum((0..3).map(|m| {
            let norm = self.coating.incoming[m].strain_mandel_s_m.iter()
                .fold(0.0, |n, e| det::hypot(n, e.abs()));
            norm * (coefficients[2 * m].abs() + coefficients[2 * m + 1].abs())
        }));
        if !coating_strain_norm_bound.is_finite()
            || coating_strain_norm_bound > self.coating.state.maximum_absolute_mandel_strain()
        {
            return Err(invalid("normal elastic coating exceeds its original whole-layer strain domain"));
        }
        let faces = [
            self.field_at(coefficients, 0.0)?,
            self.field_at(coefficients, self.thickness_m)?,
        ];
        let mut velocity = core::array::from_fn(|i| {
            complex_dot((0..3).map(|j| (self.rear_velocity[i][j], incoming[j])))
        });
        let (traction, absorbed) = if let Some(z) = self.impedance {
            (z.traction(velocity), z.absorbed_power(velocity)?)
        } else if self.backing == ElasticSurfaceCondition::Clamped {
            velocity = [C64::ZERO; 3];
            (faces[1].traction_pa, 0.0)
        } else {
            ([C64::ZERO; 3], 0.0)
        };
        let total_input = real_sum(input_power);
        let total_output = real_sum(output_power);
        if velocity.iter().chain(&traction).any(|&v| !finite(v))
            || !total_input.is_finite() || !total_output.is_finite()
            || !balanced(total_input, total_output + absorbed, total_input)
            || ![front.flux, faces[0].normal_flux_w_m2, faces[1].normal_flux_w_m2]
                .into_iter().all(|flux| balanced(flux, absorbed, total_input))
        {
            return Err(invalid("normal elastic coating power is unrepresentable or does not close"));
        }
        let velocity_scale = real_sum(incoming.into_iter().chain(reflected).chain(coefficients).map(C64::abs));
        let traction_scale = real_sum((0..3).map(|m| {
            self.incident.impedances[m] * (incoming[m].abs() + reflected[m].abs())
                + self.coating.impedances[m] * (coefficients[2 * m].abs() + coefficients[2 * m + 1].abs())
        }));
        check_continuity(front.velocity, faces[0].velocity_m_s, velocity_scale)?;
        check_continuity(front.traction, faces[0].traction_pa, traction_scale)?;
        check_continuity(faces[1].velocity_m_s, velocity, velocity_scale)?;
        check_continuity(faces[1].traction_pa, traction, traction_scale)?;
        Ok(NormalElasticCoatingResponse {
            prepared: self,
            coefficients,
            reflection: NormalElasticSurfaceResponse {
                material_state_identity: self.incident.state.material_state_identity(),
                incoming_velocity_m_s: incoming,
                reflected_velocity_m_s: reflected,
                incoming_normal_power_w_m2: input_power,
                reflected_normal_power_w_m2: output_power,
                total_incoming_power_w_m2: total_input,
                total_reflected_power_w_m2: total_output,
                absorbed_normal_power_w_m2: absorbed,
                boundary_velocity_m_s: front.velocity,
                boundary_traction_pa: front.traction,
                boundary_normal_flux_w_m2: front.flux,
                boundary_strain_mandel: front.strain,
                strain_norm_bound: front.strain_bound,
            },
            coating_faces: Some(faces),
            coating_strain_norm_bound,
            support_velocity_m_s: velocity,
            support_traction_pa: traction,
        })
    }

    fn prepare_boundary(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialError> {
        let front: [[C64; 12]; 6] = core::array::from_fn(|j| {
            self.wave_trace(j, if j % 2 == 0 { C64::ONE } else { self.phase_through[j / 2] })
        });
        let rear: [[C64; 12]; 6] = core::array::from_fn(|j| {
            self.wave_trace(j, if j % 2 == 0 { self.phase_through[j / 2] } else { C64::ONE })
        });
        let mut matrix = [C64::ZERO; 144];
        for row in 0..6 {
            poll(cancelled)?;
            for m in 0..3 {
                let wave = self.incident.reflected_mode(m).expect("native mode index");
                matrix[12 * row + m] = trace_entry(wave, row);
            }
            for j in 0..6 { matrix[12 * row + 3 + j] = -front[j][row]; }
        }
        for axis in 0..3 {
            for j in 0..6 { matrix[12 * (6 + axis) + 3 + j] = rear[j][axis]; }
            matrix[12 * (6 + axis) + 9 + axis] = -C64::ONE;
            if self.backing == ElasticSurfaceCondition::Clamped {
                matrix[12 * (9 + axis) + 9 + axis] = C64::ONE;
            } else {
                for j in 0..6 { matrix[12 * (9 + axis) + 3 + j] = rear[j][3 + axis]; }
                if let Some(z) = self.impedance {
                    for j in 0..3 { matrix[12 * (9 + axis) + 9 + j] = z.matrix_pa_s_m()[axis][j]; }
                }
            }
        }
        let mut scales = [0.0_f64; 12];
        for row in 0..12 {
            poll(cancelled)?;
            let entries = &matrix[12 * row..12 * (row + 1)];
            scales[row] = entries.iter().fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
            if scales[row] == 0.0 || !scales[row].is_finite() || entries.iter().any(|&v| !finite(v)) {
                return Err(invalid("normal elastic coating boundary row is unrepresentable"));
            }
            for j in 0..12 { matrix[12 * row + j] = divided(matrix[12 * row + j], scales[row]); }
        }
        let factor = lu_complex_with_cancel(&matrix, 12, &mut *cancelled).map_err(solve_error)?;
        for incident in 0..3 {
            let rhs: [C64; 12] = core::array::from_fn(|row| {
                if row < 6 { -divided(trace_entry(self.incident.incoming[incident], row), scales[row]) }
                else { C64::ZERO }
            });
            let mut solution = rhs;
            factor.solve_with_cancel(&mut solution, &mut *cancelled).map_err(solve_error)?;
            // Bounded residual correction reuses the original factor. Actual
            // rear motion remains an unknown throughout, including stiff loads.
            for _ in 0..2 {
                if check_solution(&matrix, &solution, &rhs).is_ok() { break; }
                let mut correction: [C64; 12] = core::array::from_fn(|i| {
                    rhs[i] - complex_dot((0..12).map(|j| (matrix[12 * i + j], solution[j])))
                });
                factor.solve_with_cancel(&mut correction, &mut *cancelled).map_err(solve_error)?;
                for i in 0..12 { solution[i] = complex_sum([solution[i], correction[i]]); }
            }
            check_solution(&matrix, &solution, &rhs)?;
            for i in 0..3 {
                self.scattering[i][incident] = solution[i];
                self.rear_velocity[i][incident] = solution[9 + i];
            }
            for i in 0..6 { self.interior[i][incident] = solution[3 + i]; }
        }
        Ok(())
    }

    fn wave_trace(&self, index: usize, phase: C64) -> [C64; 12] {
        let wave = if index % 2 == 0 { self.coating.incoming[index / 2] }
            else { self.coating.reflected_mode(index / 2).expect("native mode index") };
        core::array::from_fn(|row| trace_entry(wave, row) * phase)
    }

    fn field_at(
        &self,
        coefficients: [C64; 6],
        depth_m: f64,
    ) -> Result<NormalElasticCoatingField, MaterialError> {
        if self.thickness_m == 0.0 || !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.thickness_m {
            return Err(invalid("normal elastic coating field needs a depth in a positive-thickness layer"));
        }
        let mut positive = [C64::ZERO; 3];
        let mut negative = [C64::ZERO; 3];
        for m in 0..3 {
            let forward = if depth_m == 0.0 { C64::ONE }
                else if depth_m == self.thickness_m { self.phase_through[m] }
                else { propagation_phase(self.wavenumbers_rad_m[m], Dd::from_f64(depth_m))? };
            let backward = if depth_m == self.thickness_m { C64::ONE }
                else if depth_m == 0.0 { self.phase_through[m] }
                else { propagation_phase(self.wavenumbers_rad_m[m], Dd::from_f64(self.thickness_m) - Dd::from_f64(depth_m))? };
            positive[m] = coefficients[2 * m] * forward;
            negative[m] = coefficients[2 * m + 1] * backward;
        }
        let field = self.coating.field(positive, negative)?;
        Ok(NormalElasticCoatingField {
            velocity_m_s: field.velocity,
            traction_pa: field.traction,
            strain_mandel: field.strain,
            normal_flux_w_m2: field.flux,
        })
    }
}

impl<'p, 'a> NormalElasticCoatingResponse<'p, 'a> {
    /// Query the admitted coating at physical depth 0<=z<=h. The wave
    /// amplitudes and original material/frame remain those of this response.
    /// An exactly absent layer has no interior field, including at z=0.
    pub fn field_at(&self, depth_m: f64) -> Result<NormalElasticCoatingField, MaterialError> {
        self.prepared.field_at(self.coefficients, depth_m)
    }
}

fn trace_entry(wave: ElasticWaveTrace, row: usize) -> C64 {
    match row {
        0..3 => wave.velocity[row],
        3..6 => wave.traction_pa_s_m[row - 3],
        _ => wave.strain_mandel_s_m[row - 6],
    }
}

fn propagation_phase(wavenumber: Dd, distance: Dd) -> Result<C64, MaterialError> {
    let angle = wavenumber * distance;
    if !angle.hi.is_finite() || !angle.lo.is_finite()
        || (distance != Dd::ZERO && angle.hi == 0.0 && angle.lo == 0.0)
    {
        return Err(invalid("normal elastic coating propagation phase is unrepresentable"));
    }
    // The existing DD range-reduction owner retains phase-product and rear
    // distance-subtraction low bits. No custom trigonometric implementation.
    let (sine, cosine) = det::sin_cos_dd(angle);
    let result = C64::new(cosine.to_f64(), sine.to_f64());
    if !finite(result) { return Err(invalid("normal elastic coating phase is unresolved")); }
    Ok(result)
}

fn check_continuity(first: [C64; 3], second: [C64; 3], scale: f64) -> Result<(), MaterialError> {
    let defect = (0..3).fold(0.0, |n, i| det::hypot(n, (first[i] - second[i]).abs()));
    if !scale.is_finite() || !defect.is_finite() || !balanced(defect, 0.0, scale) {
        return Err(invalid("normal elastic coating physical fields are discontinuous"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
