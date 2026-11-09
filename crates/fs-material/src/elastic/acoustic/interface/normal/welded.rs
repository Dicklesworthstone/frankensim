//! Six-channel welded contact of two original anisotropic normal-wave owners.

use super::*;

/// Prepared welded contact. Both characteristics' surface axes describe the
/// SAME physical frame; +normal points from the first solid into the second.
/// Their frame matrices can differ because the original material axes differ.
/// A native 6x6 LU enforces all three velocity and traction continuities.
#[derive(Debug)]
pub struct PreparedNormalElasticInterface<'a> {
    characteristics: [NormalElasticCharacteristic<'a>; 2],
    scattering: [[C64; 6]; 6],
}

/// Coherent two-sided normal reflection/transmission. Modal triples retain
/// each original material's ascending-speed basis and same-polarization sign.
#[derive(Clone, Debug)]
pub struct NormalElasticInterfaceResponse {
    pub material_state_identities: [ContentHash; 2],
    pub incoming_velocity_m_s: [[C64; 3]; 2],
    pub outgoing_velocity_m_s: [[C64; 3]; 2],
    pub incoming_normal_power_w_m2: [[f64; 3]; 2],
    pub outgoing_normal_power_w_m2: [[f64; 3]; 2],
    pub total_incoming_power_w_m2: f64,
    pub total_outgoing_power_w_m2: f64,
    /// Both physical velocities use the common surface Cartesian frame.
    pub boundary_velocity_m_s: [[C64; 3]; 2],
    /// Both tractions act on the SAME +normal, not opposite outward normals.
    pub boundary_traction_pa: [[C64; 3]; 2],
    /// Original coherent -Re(v*.traction)/2, on the same positive normal.
    pub boundary_normal_flux_w_m2: [f64; 2],
    /// Each strain is in that side's ORIGINAL MATERIAL frame.
    pub boundary_strain_mandel: [[C64; 6]; 2],
    /// Conservative original strain-domain bound over each whole half-space.
    pub strain_norm_bound: [f64; 2],
}

impl<'a> PreparedNormalElasticInterface<'a> {
    pub fn new(
        first: &NormalElasticCharacteristic<'a>,
        second: &NormalElasticCharacteristic<'a>,
    ) -> Result<Self, MaterialError> {
        Self::new_cancellable(first, second, &mut || false)
    }

    /// Cache all six propagating incident channels with bounded cancellation.
    pub fn new_cancellable(
        first: &NormalElasticCharacteristic<'a>,
        second: &NormalElasticCharacteristic<'a>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MaterialError> {
        poll(cancelled)?;
        let mut result = Self {
            characteristics: [first.clone(), second.clone()],
            scattering: [[C64::ZERO; 6]; 6],
        };
        if first.polarizations == second.polarizations && first.impedances == second.impedances {
            for j in 0..6 {
                result.scattering[(j + 3) % 6][j] = C64::ONE;
            }
            poll(cancelled)?;
            return Ok(result);
        }
        let outgoing: [ElasticWaveTrace; 6] = core::array::from_fn(|j| {
            if j < 3 {
                first.reflected_mode(j).expect("native mode index")
            } else {
                second.incoming[j - 3]
            }
        });
        let incoming: [ElasticWaveTrace; 6] = core::array::from_fn(|j| {
            if j < 3 {
                first.incoming[j]
            } else {
                second.reflected_mode(j - 3).expect("native mode index")
            }
        });
        let mut matrix = [C64::ZERO; 36];
        let mut scales = [0.0_f64; 6];
        for row in 0..6 {
            poll(cancelled)?;
            for j in 0..6 {
                matrix[6 * row + j] = entry(outgoing[j], row)
                    .scale(if j < 3 { 1.0 } else { -1.0 });
                scales[row] = scales[row].max(matrix[6 * row + j].re.abs());
            }
            if scales[row] <= 0.0 || !scales[row].is_finite() {
                return Err(invalid("normal welded boundary row is unrepresentable"));
            }
            for j in 0..6 {
                matrix[6 * row + j] = divided(matrix[6 * row + j], scales[row]);
            }
        }
        let factor = lu_complex_with_cancel(&matrix, 6, &mut *cancelled).map_err(solve_error)?;
        for j in 0..6 {
            let rhs: [C64; 6] = core::array::from_fn(|row| {
                divided(entry(incoming[j], row), scales[row])
                    .scale(if j < 3 { -1.0 } else { 1.0 })
            });
            let mut solution = rhs;
            factor.solve_with_cancel(&mut solution, &mut *cancelled).map_err(solve_error)?;
            check_solution(&matrix, &solution, &rhs)?;
            for i in 0..6 {
                result.scattering[i][j] = solution[i];
            }
        }
        poll(cancelled)?;
        Ok(result)
    }

    /// One original characteristic and its original material/frame.
    pub const fn characteristic(&self, side: ElasticSide) -> &NormalElasticCharacteristic<'a> {
        &self.characteristics[side.index()]
    }

    /// Outgoing triples on both sides per incident original mode on side.
    pub fn velocity_coefficients(
        &self,
        side: ElasticSide,
        index: usize,
    ) -> Option<[[C64; 3]; 2]> {
        (index < 3).then(|| core::array::from_fn(|s| {
            core::array::from_fn(|i| self.scattering[3 * s + i][3 * side.index() + index])
        }))
    }

    /// Incoming Cartesian velocities on both sides in their COMMON surface
    /// axes. The second side travels toward -normal with the same polarization.
    pub fn drive_cartesian(
        &self,
        incoming: [[C64; 3]; 2],
    ) -> Result<NormalElasticInterfaceResponse, MaterialError> {
        self.drive(core::array::from_fn(|s| {
            self.characteristics[s].surface_velocity_to_modal(incoming[s])
        }))
    }

    /// Coherent native quasi-mode velocities, with no allocation or new solve.
    pub fn drive(
        &self,
        incoming: [[C64; 3]; 2],
    ) -> Result<NormalElasticInterfaceResponse, MaterialError> {
        let input_power = [
            self.characteristics[0].powers(incoming[0])?,
            self.characteristics[1].powers(incoming[1])?,
        ];
        let outgoing: [[C64; 3]; 2] = core::array::from_fn(|s| core::array::from_fn(|i| {
            complex_dot((0..6).map(|j| {
                (self.scattering[3 * s + i][j], incoming[j / 3][j % 3])
            }))
        }));
        let output_power = [
            self.characteristics[0].powers(outgoing[0])?,
            self.characteristics[1].powers(outgoing[1])?,
        ];
        let fields = [
            self.characteristics[0].field(incoming[0], outgoing[0])?,
            self.characteristics[1].field(outgoing[1], incoming[1])?,
        ];
        let total_input = real_sum(input_power.into_iter().flatten());
        let total_output = real_sum(output_power.into_iter().flatten());
        let first_flux = real_sum(input_power[0]) - real_sum(output_power[0]);
        let second_flux = real_sum(output_power[1]) - real_sum(input_power[1]);
        if !total_input.is_finite() || !total_output.is_finite()
            || !balanced(total_input, total_output, total_input)
            || !balanced(fields[0].flux, first_flux, total_input)
            || !balanced(fields[1].flux, second_flux, total_input)
            || !balanced(fields[0].flux, fields[1].flux, total_input)
        {
            return Err(invalid("normal welded elastic power is unrepresentable or does not close"));
        }
        // Check the physical traces independently of the cached solve. Scale
        // each field by incident/outgoing wave magnitudes, allowing coherent
        // cancellation without inventing a boundary-field absolute floor.
        for traction in [false, true] {
            let mut scale = 0.0_f64;
            for s in 0..2 {
                for m in 0..3 {
                    let wave = self.characteristics[s].incoming[m];
                    let components = if traction { wave.traction_pa_s_m } else { wave.velocity };
                    let norm = components.iter().fold(0.0, |n, v| det::hypot(n, v.abs()));
                    scale += norm * (incoming[s][m].abs() + outgoing[s][m].abs());
                }
            }
            let a = if traction { fields[0].traction } else { fields[0].velocity };
            let b = if traction { fields[1].traction } else { fields[1].velocity };
            let defect = (0..3).fold(0.0, |n, i| det::hypot(n, (a[i] - b[i]).abs()));
            if !scale.is_finite() || !balanced(defect, 0.0, scale) {
                return Err(invalid("normal welded elastic fields are discontinuous"));
            }
        }
        Ok(NormalElasticInterfaceResponse {
            material_state_identities: self.characteristics.each_ref()
                .map(|c| c.state.material_state_identity()),
            incoming_velocity_m_s: incoming,
            outgoing_velocity_m_s: outgoing,
            incoming_normal_power_w_m2: input_power,
            outgoing_normal_power_w_m2: output_power,
            total_incoming_power_w_m2: total_input,
            total_outgoing_power_w_m2: total_output,
            boundary_velocity_m_s: fields.map(|f| f.velocity),
            boundary_traction_pa: fields.map(|f| f.traction),
            boundary_normal_flux_w_m2: fields.map(|f| f.flux),
            boundary_strain_mandel: fields.map(|f| f.strain),
            strain_norm_bound: fields.map(|f| f.strain_bound),
        })
    }
}

fn entry(wave: ElasticWaveTrace, row: usize) -> C64 {
    if row < 3 { wave.velocity[row] } else { wave.traction_pa_s_m[row - 3] }
}
