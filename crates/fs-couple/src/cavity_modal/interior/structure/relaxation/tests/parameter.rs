use super::*;
use crate::cavity_modal::interior::structure::resonance::{
    StructuralInteriorParameterDirection, StructuralInteriorResponseField,
    StructuralInteriorSourceParameterJet, StructuralInteriorSources,
};
use fs_modal::harmonic::resonance::{
    HarmonicPowerParameterJet, HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet,
};
use fs_phs::mechanical::MechanicalRelaxationDirection;

// Independent five-order scalar calculus of the original physical two-equation
// uniform-pressure law. No native matrix derivative or work helper is reused.
#[derive(Clone, Copy)]
struct Jet([C64; 5]);
impl Jet {
    fn law(value: f64, parameter: f64) -> Self {
        Self([
            C64::from_re(value),
            C64::ZERO,
            C64::ZERO,
            C64::from_re(parameter),
            C64::ZERO,
        ])
    }
    fn scale(self, c: f64) -> Self {
        Self(self.0.map(|z| z.scale(c)))
    }
    fn conj(self) -> Self {
        Self(self.0.map(C64::conj))
    }
    fn inverse(self) -> Self {
        let a = C64::ONE / self.0[0];
        Self([
            a,
            -(self.0[1] * a * a),
            (self.0[1] * self.0[1] * a * a * a).scale(2.0) - self.0[2] * a * a,
            -(self.0[3] * a * a),
            (self.0[1] * self.0[3] * a * a * a).scale(2.0) - self.0[4] * a * a,
        ])
    }
}
impl core::ops::Add for Jet {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self(core::array::from_fn(|i| self.0[i] + b.0[i]))
    }
}
impl core::ops::Sub for Jet {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + b.scale(-1.0)
    }
}
impl core::ops::Mul for Jet {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self([
            self.0[0] * b.0[0],
            self.0[1] * b.0[0] + self.0[0] * b.0[1],
            self.0[2] * b.0[0] + (self.0[1] * b.0[1]).scale(2.0) + self.0[0] * b.0[2],
            self.0[3] * b.0[0] + self.0[0] * b.0[3],
            self.0[4] * b.0[0] + self.0[1] * b.0[3] + self.0[3] * b.0[1] + self.0[0] * b.0[4],
        ])
    }
}
impl core::ops::Div for Jet {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        self * b.inverse()
    }
}
fn work(a: Jet, b: Jet) -> Jet {
    (a.conj() * b).scale(0.5)
}
fn check_power(actual: HarmonicPowerParameterJet, expected: Jet) {
    for (a, b) in [
        actual.sample.average_power_w,
        actual.sample.angular_derivative,
        actual.angular_second_derivative,
        actual.parameter_derivative,
        actual.angular_parameter_derivative,
    ]
    .into_iter()
    .zip(expected.0)
    {
        assert!(
            (a - b.re).abs() < 8e-9 * b.re.abs().max(1.0),
            "{a:e} != {:e}",
            b.re
        );
    }
}
fn check_field(actual: HarmonicTransferFieldParameterJet<'_>, index: usize, expected: Jet) {
    for (a, b) in [
        actual.values[index],
        actual.angular_derivative[index],
        actual.angular_second_derivative[index],
        actual.parameter_derivative[index],
        actual.angular_parameter_derivative[index],
    ]
    .into_iter()
    .zip(expected.0)
    {
        close(a, b, 8e-9);
    }
}
fn check_scalar(actual: HarmonicTransferParameterJet, expected: Jet) {
    for (a, b) in [
        actual.sample.amplitude,
        actual.sample.angular_derivative,
        actual.angular_second_derivative,
        actual.parameter_derivative,
        actual.angular_parameter_derivative,
    ]
    .into_iter()
    .zip(expected.0)
    {
        close(a, b, 8e-9);
    }
}
fn balances(r: &StructuralInteriorRelaxationFrequencyParameterResponse) {
    for p in [
        r.power.structural_balance(),
        r.power.acoustic_balance(),
        r.power.relaxation_balance(),
        r.power.total_balance(),
    ] {
        check_power(p.unwrap(), Jet::law(0.0, 0.0));
    }
    assert!(r.relative_parameter_derivative_residual < 1e-9);
    assert!(r.relative_mixed_derivative_residual < 1e-9);
}
fn solve(
    p: &Prepared<'_>,
    d: &PreparedDirection<'_, '_>,
    omega: f64,
    f: Jet,
    q: Jet,
) -> StructuralInteriorRelaxationFrequencyParameterResponse {
    let forces = f.0.map(|z| [z]);
    let flows = q.0.map(|z| [z.scale(0.25); 4]);
    let source = |i: usize| StructuralInteriorSources {
        structural_force: &forces[i],
        inward_flow_load: &flows[i],
    };
    p.factor(omega, &mut || false)
        .unwrap()
        .solve_with_frequency_parameter_derivative(
            StructuralInteriorSourceParameterJet {
                values: source(0),
                angular_derivative: Some(source(1)),
                angular_second_derivative: Some(source(2)),
                parameter_derivative: Some(source(3)),
                angular_parameter_derivative: Some(source(4)),
            },
            d,
            &mut || false,
        )
        .unwrap()
}

#[test]
fn mixed_material_cavity_tuning_matches_original_law_with_coherent_source_jets() {
    let mut cavity = Cavity::new();
    cavity.damping[0] = 0.13;
    for i in 1..5 {
        cavity.damping[i * 5 + i] = 0.02;
    }
    let branches = [MaxwellDampingBranch {
        strain: vec![1.0, 0.05, 0.05, 0.05, 0.05],
        stiffness: 1.25,
        relaxation_time_s: 0.8,
    }];
    let prepared =
        Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let dm: Vec<_> = (0..4)
        .flat_map(|i| (0..4).map(move |j| (i, j)))
        .map(|(i, j)| cavity.mass[(i + 1) * 5 + j + 1] * (-0.07 / (5.0 / 6.0)))
        .collect();
    let mut dw = [0.0; 16];
    for i in 0..4 {
        dw[i * 4 + i] = 0.015 / 4.0;
    }
    let db = [0.11, -0.03 / 4.0, -0.03 / 4.0, -0.03 / 4.0, -0.03 / 4.0];
    let branch_d = [MechanicalRelaxationDirection {
        projection: Some(&db),
        stiffness: 0.31,
        relaxation_time_s: 0.16,
    }];
    let direction = prepared
        .prepare_parameter_direction(
            StructuralInteriorRelaxationParameterDirection {
                spatial: StructuralInteriorParameterDirection {
                    structural_mass: &[0.21],
                    structural_damping: &[0.03],
                    structural_stiffness: &[-0.12],
                    acoustic_mass: &dm,
                    acoustic_stiffness: &[0.0; 16],
                    wall_damping: &dw,
                    coupling: &[0.01; 4],
                },
                branches: Some(&branch_d),
            },
            &mut || false,
        )
        .unwrap();
    let f = Jet([
        C64::new(0.7, -0.2),
        C64::new(0.03, 0.08),
        C64::new(-0.05, 0.02),
        C64::new(0.1, -0.04),
        C64::new(-0.02, 0.06),
    ]);
    let q = Jet([
        C64::new(-0.02, 0.05),
        C64::new(0.01, -0.03),
        C64::new(0.02, 0.01),
        C64::new(-0.01, 0.03),
        C64::new(0.04, -0.02),
    ]);
    let (m, c, k, ma, wall, b, bs, ba, arm, tau) = (
        Jet::law(1.0, 0.21),
        Jet::law(0.13, 0.03),
        Jet::law(0.7, -0.12),
        Jet::law(5.0 / 6.0, -0.07),
        Jet::law(0.08, 0.015),
        Jet::law(0.5, 0.04),
        Jet::law(1.0, 0.11),
        Jet::law(0.2, -0.03),
        Jet::law(1.25, 0.31),
        Jet::law(0.8, 0.16),
    );
    for omega in [0.2, 0.7, 1.1, 1.4, 2.0, 7.0] {
        let s = Jet([
            C64::new(0.0, -omega),
            C64::new(0.0, -1.0),
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ]);
        let arm_transfer = arm / (s + tau.inverse());
        let aa = m * s * s + c * s + k + arm_transfer * bs * bs * s;
        let ab = b.scale(-1.0) + arm_transfer * bs * ba;
        let ba_matrix = (b + arm_transfer * ba * bs) * s;
        let bb = ma * s + wall + arm_transfer * ba * ba;
        let det = aa * bb - ab * ba_matrix;
        let u = (f * bb - ab * q) / det;
        let pressure = (aa * q - ba_matrix * f) / det;
        let velocity = s * u;
        let stress = arm_transfer * (bs * velocity + ba * pressure);
        let result = solve(&prepared, &direction, omega, f, q);
        check_field(
            result.field_parameter_jet(StructuralInteriorResponseField::StructuralDisplacement),
            0,
            u,
        );
        check_field(
            result.field_parameter_jet(StructuralInteriorResponseField::StructuralVelocity),
            0,
            velocity,
        );
        check_field(
            result.field_parameter_jet(StructuralInteriorResponseField::FluidForce),
            0,
            b * pressure,
        );
        for i in 0..4 {
            check_field(
                result.field_parameter_jet(StructuralInteriorResponseField::Pressure),
                i,
                pressure,
            );
            check_field(
                result.field_parameter_jet(StructuralInteriorResponseField::PressurePotential),
                i,
                pressure / s,
            );
            check_field(
                result.field_parameter_jet(StructuralInteriorResponseField::DisplacedFlow),
                i,
                (b * velocity).scale(0.25),
            );
        }
        check_scalar(result.branch_stress_parameter_jet(0).unwrap(), stress);
        check_power(result.power.structural_input, work(velocity, f));
        check_power(result.power.acoustic_input, work(pressure, q));
        check_power(
            result.power.fluid_to_structure,
            work(velocity, b * pressure),
        );
        check_power(
            result.power.structural_dissipation,
            work(velocity, velocity) * c,
        );
        check_power(
            result.power.wall_dissipation,
            work(pressure, pressure) * wall,
        );
        let heat = work(stress, stress) / (arm * tau);
        let stored = work(stress, stress) / arm.scale(2.0);
        check_power(result.memory[0].dissipation, heat);
        check_power(result.power.relaxation_dissipation, heat);
        check_power(
            result.memory[0].structural_work,
            work(bs * velocity, stress),
        );
        check_power(result.memory[0].acoustic_work, work(ba * pressure, stress));
        near(
            result.memory[0].mean_stored_energy_parameter_derivative,
            stored.0[3].re,
        );
        near(
            result.memory[0].mean_stored_energy_mixed_derivative,
            stored.0[4].re,
        );
        let energy =
            (work(velocity, velocity) * m + work(u, u) * k + work(pressure, pressure) * ma)
                .scale(0.5)
                + stored;
        near(
            result.mean_stored_energy_parameter_derivative,
            energy.0[3].re,
        );
        near(result.mean_stored_energy_mixed_derivative, energy.0[4].re);
        balances(&result);

        // A changing common propagation delay must not manufacture a heat
        // sensitivity or erase its material/source-amplitude derivative.
        let (delay, delay_theta) = (1.7, 0.23);
        let phase = C64::new((omega * delay).cos(), (omega * delay).sin());
        let phase = Jet([
            phase,
            C64::new(0.0, delay) * phase,
            phase.scale(-delay * delay),
            C64::new(0.0, omega * delay_theta) * phase,
            C64::new(-omega * delay * delay_theta, delay_theta) * phase,
        ]);
        let delayed = solve(&prepared, &direction, omega, f * phase, q * phase);
        check_field(
            delayed.field_parameter_jet(StructuralInteriorResponseField::StructuralDisplacement),
            0,
            u * phase,
        );
        check_scalar(
            delayed.branch_stress_parameter_jet(0).unwrap(),
            stress * phase,
        );
        check_power(delayed.power.relaxation_dissipation, heat);
        check_power(delayed.power.structural_input, work(velocity, f));
        check_power(delayed.power.acoustic_input, work(pressure, q));
        balances(&delayed);
    }
}

fn elastic_direction() -> StructuralInteriorRelaxationParameterDirection<'static> {
    StructuralInteriorRelaxationParameterDirection {
        spatial: StructuralInteriorParameterDirection {
            structural_mass: &[0.0],
            structural_damping: &[0.0],
            structural_stiffness: &[0.3],
            acoustic_mass: &[0.0; 16],
            acoustic_stiffness: &[0.0; 16],
            wall_damping: &[0.0; 16],
            coupling: &[0.0; 4],
        },
        branches: None,
    }
}

#[test]
fn zero_memory_parameter_limit_and_refusals_preserve_reusable_original_factor() {
    let cavity = Cavity::new();
    let prepared = Prepared::new(cavity.source(), &[], Limits::default(), &mut || false).unwrap();
    let direction = prepared
        .prepare_parameter_direction(elastic_direction(), &mut || false)
        .unwrap();
    let omega = 0.6;
    let factor = prepared.factor(omega, &mut || false).unwrap();
    let source = StructuralInteriorSourceParameterJet::fixed(StructuralInteriorSources {
        structural_force: &[C64::ONE],
        inward_flow_load: &[C64::ZERO; 4],
    });
    let result = factor
        .solve_with_frequency_parameter_derivative(source, &direction, &mut || false)
        .unwrap();
    let s = Jet([
        C64::new(0.0, -omega),
        C64::new(0.0, -1.0),
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
    ]);
    let u = Jet::law(1.0, 0.0) / (s * s + Jet::law(1.0, 0.3));
    check_field(
        result.field_parameter_jet(StructuralInteriorResponseField::StructuralDisplacement),
        0,
        u,
    );
    check_field(
        result.field_parameter_jet(StructuralInteriorResponseField::Pressure),
        0,
        u.scale(-0.6),
    );
    assert!(result.memory.is_empty());
    assert!(result.branch_stress_parameter_jet(0).is_err());
    check_power(result.power.relaxation_dissipation, Jet::law(0.0, 0.0));
    balances(&result);
    assert!(matches!(
        factor.solve_with_frequency_parameter_derivative(source, &direction, &mut || true),
        Err(Error::Cancelled)
    ));
    let bad = StructuralInteriorSourceParameterJet {
        angular_parameter_derivative: Some(StructuralInteriorSources {
            structural_force: &[],
            inward_flow_load: &[C64::ZERO; 4],
        }),
        ..source
    };
    assert!(
        factor
            .solve_with_frequency_parameter_derivative(bad, &direction, &mut || false)
            .is_err()
    );
    let other = Prepared::new(cavity.source(), &[], Limits::default(), &mut || false).unwrap();
    assert!(
        other
            .factor(omega, &mut || false)
            .unwrap()
            .solve_with_frequency_parameter_derivative(source, &direction, &mut || false)
            .is_err()
    );
    let replay = factor
        .solve_with_frequency_parameter_derivative(source, &direction, &mut || false)
        .unwrap();
    assert_eq!(
        result.parameter_derivative.structural_displacement,
        replay.parameter_derivative.structural_displacement
    );
    let bounded = Prepared::new(
        cavity.source(),
        &[],
        Limits {
            maximum_adapter_bytes: 40_000,
            ..Limits::default()
        },
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        bounded.prepare_parameter_direction(elastic_direction(), &mut || false),
        Err(Error::Limit { .. })
    ));
    assert!(matches!(
        prepared.prepare_parameter_direction(elastic_direction(), &mut || true),
        Err(Error::Cancelled)
    ));
}
