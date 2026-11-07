//! Joint mass/entropy ports on the original two-fluid wave operator.

use super::boundary::endpoint_inverse;
use super::*;

/// Exact reciprocal endpoint response when BOTH component currents may pass.
///
/// At each end the outward flows are `[A*v, A*(v_n-v_s)]` in m³/s, where
/// `v=f_n*v_n+f_s*v_s`. Their conjugate efforts are `[p', rho_s*s*T']` in Pa.
/// The dot product is outward wave-availability power in W. In particular,
/// the second flow is neither entropy current nor physical heat by itself.
///
/// Flattened matrix order is `[left bulk, left relative, right bulk, right
/// relative]`. The exact affine relation is `e = e_reference - R*(q-q_reference)`.
/// Cross-component entries retain thermomechanical coupling; they may have
/// either sign. The whole matrix is symmetric and positive semidefinite.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermomechanicalEndsResponse {
    /// Midpoint effort pairs about the selected equilibrium, in Pa.
    pub reference_effort_pa: [[f64; 2]; 2],
    /// Outward volume-current pairs at the algebraic reference, in m³/s.
    pub reference_flow_outward_m3_s: [[f64; 2]; 2],
    /// Reciprocal Schur response, in Pa s/m³.
    pub impedance_pa_s_m3: [[f64; 4]; 4],
}

impl TwoFluidChannel {
    /// Convert outward barycentric/relative volume currents to physical ports.
    ///
    /// No assumption of zero mass or zero entropy flow is made. Original
    /// amplitude and source-domain limits are enforced by the eventual step.
    pub fn thermomechanical_boundary_velocities(
        &self,
        flow_outward_m3_s: [[f64; 2]; 2],
    ) -> Result<BoundaryVelocities, ChannelError> {
        let mut ends = [EndVelocities::default(); 2];
        for (end, [bulk, relative]) in flow_outward_m3_s.into_iter().enumerate() {
            let sign = if end == 0 { -1.0 } else { 1.0 };
            let v = sign * bulk / self.area_m2;
            let w = sign * relative / self.area_m2;
            let normal = v + self.material.superfluid_fraction() * w;
            let superfluid = v - self.material.normal_fraction() * w;
            if ![bulk, relative, v, w, normal, superfluid]
                .into_iter()
                .all(f64::is_finite)
            {
                return Err(ChannelError::Domain {
                    field: "thermomechanical boundary flow",
                    index: end,
                });
            }
            ends[end] = EndVelocities {
                normal_m_s: normal,
                superfluid_m_s: superfluid,
            };
        }
        Ok(BoundaryVelocities {
            left: ends[0],
            right: ends[1],
        })
    }

    /// Query coupled pressure/temperature efforts without publishing any state.
    ///
    /// Uses one original proposal and its cached endpoint inverse columns, with
    /// O(cells) work and no allocation. Both sound modes, thermal expansion,
    /// configured normal-wall memory and the selected friction law remain on
    /// their original paths. The reference may exceed amplitude limits, but it
    /// must pass finite, residual and storage checks. The final physical step
    /// must independently pass all original limits.
    pub fn thermomechanical_ends_response(
        &mut self,
        dt_s: f64,
        reference_flow_outward_m3_s: [[f64; 2]; 2],
        force: BodyForce<'_>,
        friction: MutualFriction,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermomechanicalEndsResponse, ChannelError> {
        let boundary = self.thermomechanical_boundary_velocities(reference_flow_outward_m3_s)?;
        self.stage_step_with_limits(dt_s, boundary, friction, force, false, &mut cancelled)?;
        let reference_effort_pa = self.thermomechanical_midpoint_efforts()?;
        let rotation = self.wave_cell_rotation(dt_s, friction)?;
        let bulk = self.modes.face_to_modes(1.0, 0.0);
        let relative = self.modes.face_to_modes(0.0, 1.0);
        let columns = [bulk, relative].map(|column| {
            let original = [
                self.modes.speed[0] * column[0],
                self.modes.speed[1] * column[1],
            ];
            rotation.map(|row| row[0] * original[0] + row[1] * original[1])
        });
        let factor = 0.5 * dt_s / self.cell_volume;
        let mut impedance = [[0.0; 4]; 4];
        let last = self.cells() - 1;
        for mode in 0..2 {
            for (end, cell) in [0, last].into_iter().enumerate() {
                let inverse =
                    endpoint_inverse(&mut self.workspace[mode], end, cell, &mut cancelled)?;
                for other in 0..2 {
                    for row in 0..2 {
                        for column in 0..2 {
                            impedance[2 * other + row][2 * end + column] += (factor
                                * columns[row][mode])
                                * columns[column][mode]
                                * inverse[other];
                        }
                    }
                }
            }
        }
        for row in 0..4 {
            if !positive(impedance[row][row]) {
                return Err(ChannelError::Invalid(
                    "thermomechanical response diagonal is unresolved",
                ));
            }
            for column in 0..row {
                let a = impedance[row][column];
                let b = impedance[column][row];
                let scale = impedance[row][row].sqrt() * impedance[column][column].sqrt();
                if ![a, b, scale].into_iter().all(f64::is_finite)
                    || (a - b).abs() > 8192.0 * f64::EPSILON * scale
                {
                    return Err(ChannelError::Invalid(
                        "thermomechanical endpoint reciprocity is unresolved",
                    ));
                }
                let mean = 0.5 * a + 0.5 * b;
                impedance[row][column] = mean;
                impedance[column][row] = mean;
            }
        }
        poll(&mut cancelled)?;
        Ok(ThermomechanicalEndsResponse {
            reference_effort_pa,
            reference_flow_outward_m3_s,
            impedance_pa_s_m3: impedance,
        })
    }

    pub(super) fn thermomechanical_midpoint_efforts(&self) -> Result<[[f64; 2]; 2], ChannelError> {
        let entropy_coefficient =
            self.material.superfluid_density_kg_m3() * self.material.entropy_j_kg_k();
        let mut efforts = [[0.0; 2]; 2];
        for (end, cell) in [0, self.cells() - 1].into_iter().enumerate() {
            let old = self.modes.cell(
                &self.material,
                [self.accepted[0].cells[cell], self.accepted[1].cells[cell]],
            );
            let new = self.modes.cell(
                &self.material,
                [self.pending[0].cells[cell], self.pending[1].cells[cell]],
            );
            efforts[end] = [
                0.5 * old.pressure_perturbation_pa + 0.5 * new.pressure_perturbation_pa,
                entropy_coefficient
                    * (0.5 * old.temperature_perturbation_k + 0.5 * new.temperature_perturbation_k),
            ];
        }
        if !positive(entropy_coefficient) || !efforts.into_iter().flatten().all(f64::is_finite) {
            return Err(ChannelError::Invalid(
                "thermomechanical endpoint efforts are not representable",
            ));
        }
        Ok(efforts)
    }
}
