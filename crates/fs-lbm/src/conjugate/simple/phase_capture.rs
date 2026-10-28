//! Explicit inelastic capture by a stationary support, without pressure projection.
use super::*;

/// Native extensive transfers, before any thermal or support publication.
#[derive(Debug, Clone)]
pub(in crate::conjugate) struct NativePhaseCapture {
    pub initial_kinetic_energy_j: f64,
    pub final_kinetic_energy_j: f64,
    pub heat_j: f64,
    pub cell_heat_j: Vec<f64>,
    pub cell_support_impulse_n_s: Vec<[f64; 3]>,
    pub captured_faces: usize,
}

impl NativePhaseCapture {
    pub fn stationary(cells: usize) -> Self {
        Self {
            initial_kinetic_energy_j: 0.0,
            final_kinetic_energy_j: 0.0,
            heat_j: 0.0,
            cell_heat_j: vec![0.0; cells],
            cell_support_impulse_n_s: vec![[0.0; 3]; cells],
            captured_faces: 0,
        }
    }
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "supported_phase_capture",
        reason: reason.into(),
    }
}

impl Solver<'_> {
    /// Each removed DOF loses p=M*u and K=M*u²/2 to the declared stationary
    /// support. Its two original half-cell volumes partition impulse and heat.
    /// A live face joining surviving liquid cannot be zeroed here: that needs
    /// a pressure/impulse redistribution solve, whose loss is not capture heat.
    pub(in crate::conjugate) fn supported_phase_capture(
        &self,
        target: &VoxelDomain,
        supported: &[bool],
        gate: &CancelGate,
    ) -> Result<(MomentumState, NativePhaseCapture), ChtError> {
        poll(gate)?;
        if target.dims() != self.n || supported.len() != self.domain.cell_count() {
            return Err(invalid(
                "capture must retain the original staggered grid and support map",
            ));
        }
        for axis in 0..3 {
            for face in 0..=self.n[axis] {
                if face % 256 == 0 {
                    poll(gate)?;
                }
                if target.face_coord(axis, face).to_bits()
                    != self.domain.face_coord(axis, face).to_bits()
                {
                    return Err(invalid("capture cannot change native reference geometry"));
                }
            }
        }
        let mut next = self.momentum_state();
        let mut report = NativePhaseCapture::stationary(self.domain.cell_count());
        for axis in 0..3 {
            for &index in &self.comps[axis].unknowns {
                poll(gate)?;
                let plus = self.comps[axis].coords(index);
                if plus[axis] == 0 || plus[axis] == self.n[axis] {
                    return Err(invalid(
                        "supported capture requires closed stationary-normal boundaries",
                    ));
                }
                let mut minus = plus;
                minus[axis] -= 1;
                let cells = [cell_of(self.n, minus), cell_of(self.n, plus)];
                let speed = self.vel[axis][index];
                let mass = self.interior_face_mass_kg(axis, index);
                finite_positive("capture.face_mass", mass)?;
                let energy = 0.5 * mass * speed * speed;
                finite("capture.face_energy", energy)?;
                report.initial_kinetic_energy_j += energy;
                let mobile = cells.map(|cell| target.is_fluid(cell));
                if mobile == [true, true] {
                    report.final_kinetic_energy_j += energy;
                    continue;
                }
                if speed == 0.0 {
                    continue;
                }
                if mobile != [false, false] {
                    return Err(invalid(
                        "a moving liquid/solid interface needs pressure-impulse redistribution; refine the phase event",
                    ));
                }
                if cells.iter().any(|&cell| !supported[cell]) {
                    return Err(invalid(
                        "moving frozen material lacks its explicitly stationary support",
                    ));
                }
                let width_minus = self.domain.width(axis, minus[axis]);
                let fraction = width_minus / (width_minus + self.domain.width(axis, plus[axis]));
                let heat_minus = fraction * energy;
                let impulse = mass * speed;
                let impulse_minus = fraction * impulse;
                for (cell, heat, impulse) in [
                    (cells[0], heat_minus, impulse_minus),
                    (cells[1], energy - heat_minus, impulse - impulse_minus),
                ] {
                    report.cell_heat_j[cell] += heat;
                    report.cell_support_impulse_n_s[cell][axis] += impulse;
                    finite("capture.cell_heat", report.cell_heat_j[cell])?;
                    finite(
                        "capture.cell_impulse",
                        report.cell_support_impulse_n_s[cell][axis],
                    )?;
                }
                report.heat_j += energy;
                report.captured_faces += 1;
                next.velocity[axis][index] = 0.0;
            }
        }
        for value in [
            report.initial_kinetic_energy_j,
            report.final_kinetic_energy_j,
            report.heat_j,
        ] {
            finite("capture.energy", value)?;
        }
        next.admit_phase_domain(target, gate)?;
        next.clear_solid_pressure(target);
        poll(gate)?;
        Ok((next, report))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjugate::Voxel;

    #[test]
    fn graded_half_cell_capture_preserves_a_separate_moving_component() {
        let mut source =
            VoxelDomain::graded([vec![1.0, 2.0, 1.0, 2.0, 1.0], vec![1.5, 2.5], vec![2.0]])
                .unwrap();
        for y in 0..2 {
            source.set(2, y, 0, Voxel::Solid(0));
        }
        let mut target = source.clone();
        let mut supported = vec![false; source.cell_count()];
        for y in 0..2 {
            for x in 0..2 {
                target.set(x, y, 0, Voxel::Solid(0));
                supported[cell_of(source.dims(), [x, y, 0])] = true;
            }
        }
        let fluid = FluidProperties {
            density_kg_m3: 2.5,
            specific_heat_j_kg_k: 2.0,
            conductivity_w_m_k: 0.1,
            kinematic_viscosity_m2_s: 0.1,
        };
        let config = SimpleConfig::new([FvBoundary::wall(); 6]);
        let mut original = Solver::new(&source, &fluid, &config);
        for (x, q) in [(0, 1.2), (3, 0.7)] {
            for (axis, at, flux) in [
                (0, [x + 1, 0, 0], -q),
                (0, [x + 1, 1, 0], q),
                (1, [x, 1, 0], q),
                (1, [x + 1, 1, 0], -q),
            ] {
                let i = original.comps[axis].index(at);
                original.vel[axis][i] = flux / original.comp_face_area(axis, i);
            }
        }
        let before = original.momentum_state();
        let (after, report) = original
            .supported_phase_capture(&target, &supported, &CancelGate::new())
            .unwrap();
        // Independent face masses: 11.25,18.75,10,20 kg; velocities -.4,.24,.6,-.3.
        let expected = 0.5
            * (11.25 * 0.4_f64.powi(2)
                + 18.75 * 0.24_f64.powi(2)
                + 10.0 * 0.6_f64.powi(2)
                + 20.0 * 0.3_f64.powi(2));
        assert!((report.heat_j - expected).abs() < 1e-13);
        assert!(
            (report.initial_kinetic_energy_j - report.final_kinetic_energy_j - expected).abs()
                < 1e-13
        );
        assert_eq!(report.captured_faces, 4);
        for (cell, heat) in [(0, 0.975), (1, 0.9375), (5, 1.305), (6, 0.9225)] {
            assert!((report.cell_heat_j[cell] - heat).abs() < 1e-13);
        }
        for axis in 0..3 {
            for (i, &speed) in before.velocity[axis].iter().enumerate() {
                let at = original.comps[axis].coords(i);
                if at[0] >= 3 {
                    assert_eq!(after.velocity[axis][i].to_bits(), speed.to_bits());
                }
            }
        }
        let mut remaining = Solver::new(&target, &fluid, &config);
        remaining.restore_momentum(&after).unwrap();
        let flow = remaining.finish_projected(&fluid, 0, (0.0, 0.0));
        assert!(flow.field.max_divergence(&target) < 1e-14);
        assert!(report.final_kinetic_energy_j > 0.0);
        let mut impulse = [0.0; 3];
        let mut angular_z = 0.0;
        for (cell, p) in report.cell_support_impulse_n_s.iter().enumerate() {
            let [x, y, z] = source.coords(cell);
            let r = source.center(x, y, z);
            for axis in 0..3 {
                impulse[axis] += p[axis];
            }
            angular_z += r[0] * p[1] - r[1] * p[0];
        }
        assert!(impulse.iter().all(|v| v.abs() < 1e-13));
        assert!((angular_z + 18.0).abs() < 1e-13);
    }
}
