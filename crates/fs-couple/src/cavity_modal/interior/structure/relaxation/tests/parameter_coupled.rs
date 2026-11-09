//! Manufactured original fields, differentiated independently of the native solve.
use super::super::parameter::StructuralInteriorRelaxationParameterDirection;
use super::*;
use crate::cavity_modal::interior::structure::resonance::{
    StructuralInteriorParameterDirection, StructuralInteriorResponseField,
    StructuralInteriorSourceParameterJet, StructuralInteriorSources,
};
use fs_modal::harmonic::resonance::HarmonicPowerParameterJet;
use fs_phs::mechanical::MechanicalRelaxationDirection;

// Test-only scalar product rules, in value/omega/omega-omega/theta/omega-theta
// order. The source below is manufactured from physical q and stress laws;
// no production J/R/Q/G derivative or solver participates in the oracle.
#[derive(Clone, Copy, Default)]
struct Jet([C64; 5]);
impl Jet {
    fn real(value: f64, derivative: f64) -> Self {
        Self([
            C64::from_re(value),
            C64::ZERO,
            C64::ZERO,
            C64::from_re(derivative),
            C64::ZERO,
        ])
    }
    fn conj(self) -> Self {
        Self(self.0.map(C64::conj))
    }
    fn inverse(self) -> Self {
        let [a, w, ww, t, wt] = self.0;
        let v = C64::ONE / a;
        Self([
            v,
            -w * v * v,
            (w * w * v * v * v).scale(2.0) - ww * v * v,
            -t * v * v,
            (w * t * v * v * v).scale(2.0) - wt * v * v,
        ])
    }
}
impl core::ops::Add for Jet {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self(core::array::from_fn(|i| self.0[i] + b.0[i]))
    }
}
impl core::ops::Mul for Jet {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let [a, w, ww, t, wt] = self.0;
        let [b, bw, bww, bt, bwt] = b.0;
        Self([
            a * b,
            w * b + a * bw,
            ww * b + (w * bw).scale(2.0) + a * bww,
            t * b + a * bt,
            wt * b + w * bt + t * bw + a * bwt,
        ])
    }
}
impl core::ops::Div for Jet {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        self * b.inverse()
    }
}
fn sum(values: impl Iterator<Item = Jet>) -> Jet {
    values.fold(Jet::default(), |a, b| a + b)
}
fn power(actual: HarmonicPowerParameterJet, expected: Jet) {
    for (a, e) in [
        actual.sample.average_power_w,
        actual.sample.angular_derivative,
        actual.angular_second_derivative,
        actual.parameter_derivative,
        actual.angular_parameter_derivative,
    ]
    .into_iter()
    .zip(expected.0)
    {
        near(a, e.re);
    }
}

#[test]
fn manufactured_coupled_material_and_interface_jets_keep_every_original_field_and_heat_channel() {
    let ms = [1.2, 0.1, 0.1, 0.9];
    let ma = [0.8, 0.15, 0.15, 0.65];
    let cs = [0.4, 0.7, -0.5, 0.3];
    let ca = [0.08, -0.01, -0.01, 0.06];
    let ks = [5.0, 0.4, 0.4, 3.0];
    let ka = [2.0, -2.0, -2.0, 2.0];
    let b = [0.18, 0.04, -0.06, 0.12];
    let dms = [0.13, -0.02, -0.02, -0.08];
    let dma = [0.04, 0.01, 0.01, -0.03];
    let dcs = [-0.02, 0.05, -0.07, 0.03];
    let dca = [0.01, -0.004, -0.004, -0.02];
    let dks = [0.2, -0.08, -0.08, 0.1];
    let dka = [0.07, -0.07, -0.07, 0.07];
    let db = [0.03, -0.02, 0.04, 0.01];
    let mut matrices = [[Jet::default(); 16]; 3];
    for (matrix, offset, base, derivative) in [
        (0, 0, ms, dms),
        (0, 2, ma, dma),
        (1, 0, cs, dcs),
        (1, 2, ca, dca),
        (2, 0, ks, dks),
        (2, 2, ka, dka),
    ] {
        for i in 0..2 {
            for j in 0..2 {
                matrices[matrix][(offset + i) * 4 + offset + j] =
                    Jet::real(base[i * 2 + j], derivative[i * 2 + j]);
            }
        }
    }
    let interface: [Jet; 4] = core::array::from_fn(|i| Jet::real(b[i], db[i]));
    for a in 0..2 {
        for i in 0..2 {
            matrices[1][(a + 2) * 4 + i] = interface[a * 2 + i];
            matrices[1][i * 4 + a + 2] = Jet::real(-1.0, 0.0) * interface[a * 2 + i];
        }
    }
    let [mass, damping, stiffness] = matrices.map(|m| m.map(|v| v.0[0].re));
    let branches = [
        MaxwellDampingBranch {
            strain: vec![0.8, -0.3, 0.2, -0.1],
            stiffness: 1.7,
            relaxation_time_s: 0.6,
        },
        MaxwellDampingBranch {
            strain: vec![-0.2, 0.7, -0.15, 0.35],
            stiffness: 0.8,
            relaxation_time_s: 0.6,
        },
        MaxwellDampingBranch {
            strain: vec![0.0; 4],
            stiffness: 2.5,
            relaxation_time_s: 0.15,
        },
    ];
    let branch_directions = [
        MechanicalRelaxationDirection {
            projection: Some(&[-0.05, 0.07, 0.02, -0.03]),
            stiffness: 0.21,
            relaxation_time_s: -0.04,
        },
        MechanicalRelaxationDirection {
            projection: Some(&[0.03, -0.04, 0.08, 0.01]),
            stiffness: -0.13,
            relaxation_time_s: 0.06,
        },
        MechanicalRelaxationDirection {
            projection: Some(&[0.2, -0.1, 0.15, 0.05]),
            stiffness: 0.4,
            relaxation_time_s: 0.02,
        },
    ];
    let prepared = Prepared::new(
        Source {
            structural: 2,
            coordinates: 4,
            mass: &mass,
            damping: &damping,
            stiffness: &stiffness,
            coupling: &b,
        },
        &branches,
        Limits::default(),
        &mut || false,
    )
    .unwrap();
    let direction = prepared
        .prepare_parameter_direction(
            StructuralInteriorRelaxationParameterDirection {
                spatial: StructuralInteriorParameterDirection {
                    structural_mass: &dms,
                    structural_damping: &dcs,
                    structural_stiffness: &dks,
                    acoustic_mass: &dma,
                    acoustic_stiffness: &dka,
                    wall_damping: &dca,
                    coupling: &db,
                },
                branches: Some(&branch_directions),
            },
            &mut || false,
        )
        .unwrap();
    for omega in [0.4, 1.7, 3.2] {
        let s = Jet([
            C64::new(0.0, -omega),
            C64::new(0.0, -1.0),
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ]);
        let q: [Jet; 4] = core::array::from_fn(|i| {
            let i = i as f64;
            Jet([
                C64::new(
                    0.2 - 0.04 * i + 0.01 * omega * omega,
                    -0.08 + 0.03 * i + 0.02 * omega,
                ),
                C64::new(0.02 * omega, 0.02),
                C64::from_re(0.02),
                C64::new(
                    0.04 - 0.01 * i + 0.03 * omega,
                    -0.02 + 0.01 * i - 0.015 * omega,
                ),
                C64::new(0.03, -0.015),
            ])
        });
        let velocity = q.map(|q| s * q);
        // Constant acoustic potential is retained, not projected away.
        assert!((q[2].0[0] + q[3].0[0]).abs() > 0.1);
        let mut forcing: [Jet; 4] = core::array::from_fn(|i| {
            sum((0..4).map(|j| {
                (matrices[2][i * 4 + j]
                    + s * matrices[1][i * 4 + j]
                    + s * s * matrices[0][i * 4 + j])
                    * q[j]
            }))
        });
        let mut stresses = Vec::new();
        let mut memory = Vec::new();
        for (branch, d) in branches.iter().zip(branch_directions) {
            let k = Jet::real(branch.stiffness, d.stiffness);
            let tau = Jet::real(branch.relaxation_time_s, d.relaxation_time_s);
            let projection: [Jet; 4] =
                core::array::from_fn(|i| Jet::real(branch.strain[i], d.projection.unwrap()[i]));
            let solid = sum((0..2).map(|i| projection[i] * velocity[i]));
            let acoustic = sum((2..4).map(|i| projection[i] * velocity[i]));
            let stress = k * (solid + acoustic) / (s + tau.inverse());
            for i in 0..4 {
                forcing[i] = forcing[i] + projection[i] * stress;
            }
            let stored = Jet::real(0.25, 0.0) * stress.conj() * stress / k;
            let heat = Jet::real(0.5, 0.0) * stress.conj() * stress / (k * tau);
            memory.push((
                stored,
                heat,
                Jet::real(0.5, 0.0) * stress.conj() * solid,
                Jet::real(0.5, 0.0) * stress.conj() * acoustic,
            ));
            stresses.push(stress);
        }
        let loads: [Vec<C64>; 5] =
            core::array::from_fn(|order| forcing.iter().map(|j| j.0[order]).collect());
        let load = |order: usize| StructuralInteriorSources {
            structural_force: &loads[order][..2],
            inward_flow_load: &loads[order][2..],
        };
        let r = prepared
            .factor(omega, &mut || false)
            .unwrap()
            .solve_with_frequency_parameter_derivative(
                StructuralInteriorSourceParameterJet {
                    values: load(0),
                    angular_derivative: Some(load(1)),
                    angular_second_derivative: Some(load(2)),
                    parameter_derivative: Some(load(3)),
                    angular_parameter_derivative: Some(load(4)),
                },
                &direction,
                &mut || false,
            )
            .unwrap();
        let force: [Jet; 2] =
            core::array::from_fn(|i| sum((0..2).map(|a| interface[a * 2 + i] * velocity[a + 2])));
        let flow: [Jet; 2] =
            core::array::from_fn(|a| sum((0..2).map(|i| interface[a * 2 + i] * velocity[i])));
        for (field, expected) in [
            (
                StructuralInteriorResponseField::StructuralDisplacement,
                &q[..2],
            ),
            (
                StructuralInteriorResponseField::StructuralVelocity,
                &velocity[..2],
            ),
            (StructuralInteriorResponseField::PressurePotential, &q[2..]),
            (StructuralInteriorResponseField::Pressure, &velocity[2..]),
            (
                StructuralInteriorResponseField::FluidForce,
                force.as_slice(),
            ),
            (
                StructuralInteriorResponseField::DisplacedFlow,
                flow.as_slice(),
            ),
        ] {
            let f = r.field_parameter_jet(field);
            for (order, values) in [
                f.values,
                f.angular_derivative,
                f.angular_second_derivative,
                f.parameter_derivative,
                f.angular_parameter_derivative,
            ]
            .into_iter()
            .enumerate()
            {
                for (actual, expected) in values.iter().zip(expected) {
                    close(*actual, expected.0[order], 3e-10);
                }
            }
        }
        for (j, &expected) in stresses.iter().enumerate() {
            let actual = r.branch_stress_parameter_jet(j).unwrap();
            for (a, e) in [
                actual.sample.amplitude,
                actual.sample.angular_derivative,
                actual.angular_second_derivative,
                actual.parameter_derivative,
                actual.angular_parameter_derivative,
            ]
            .into_iter()
            .zip(expected.0)
            {
                close(a, e, 3e-10);
            }
            let (stored, heat, solid, acoustic) = memory[j];
            power(r.memory[j].dissipation, heat);
            power(r.memory[j].structural_work, solid);
            power(r.memory[j].acoustic_work, acoustic);
            let first = &r.frequency_response.memory[j];
            for (a, e) in [
                first.mean_stored_energy_j,
                first.mean_stored_energy_derivative,
                first.mean_stored_energy_second_derivative,
                r.memory[j].mean_stored_energy_parameter_derivative,
                r.memory[j].mean_stored_energy_mixed_derivative,
            ]
            .into_iter()
            .zip(stored.0)
            {
                near(a, e.re);
            }
        }
        let half = Jet::real(0.5, 0.0);
        let heat = |offset: usize| {
            half * sum((offset..offset + 2).flat_map(|i| {
                (offset..offset + 2)
                    .map(move |j| velocity[i].conj() * matrices[1][i * 4 + j] * velocity[j])
            }))
        };
        for (actual, expected) in [
            (
                r.power.structural_input,
                half * sum((0..2).map(|i| velocity[i].conj() * forcing[i])),
            ),
            (
                r.power.acoustic_input,
                half * sum((2..4).map(|i| velocity[i].conj() * forcing[i])),
            ),
            (
                r.power.fluid_to_structure,
                half * sum((0..2).map(|i| velocity[i].conj() * force[i])),
            ),
            (r.power.structural_dissipation, heat(0)),
            (r.power.wall_dissipation, heat(2)),
            (
                r.power.relaxation_dissipation,
                sum(memory.iter().map(|m| m.1)),
            ),
            (
                r.power.structural_to_relaxation,
                sum(memory.iter().map(|m| m.2)),
            ),
            (
                r.power.acoustic_to_relaxation,
                sum(memory.iter().map(|m| m.3)),
            ),
        ] {
            power(actual, expected);
        }
        let energy = Jet::real(0.25, 0.0)
            * sum((0..4).flat_map(|i| {
                (0..4).map(move |j| {
                    q[i].conj() * matrices[2][i * 4 + j] * q[j]
                        + velocity[i].conj() * matrices[0][i * 4 + j] * velocity[j]
                })
            }))
            + sum(memory.iter().map(|m| m.0));
        let c = &r.frequency_response;
        for (a, e) in [
            c.frequency_response.response.mean_stored_energy_j(),
            c.mean_stored_energy_derivative,
            c.mean_stored_energy_second_derivative,
            r.mean_stored_energy_parameter_derivative,
            r.mean_stored_energy_mixed_derivative,
        ]
        .into_iter()
        .zip(energy.0)
        {
            near(a, e.re);
        }
        for balance in [
            r.power.structural_balance(),
            r.power.acoustic_balance(),
            r.power.relaxation_balance(),
            r.power.total_balance(),
        ] {
            power(balance.unwrap(), Jet::default());
        }
        let dark = r.branch_stress_parameter_jet(2).unwrap();
        assert_eq!(dark.sample.amplitude, C64::ZERO);
        assert!(dark.parameter_derivative.abs() > 1e-5);
        assert_eq!(r.memory[2].dissipation.parameter_derivative, 0.0);
        assert_eq!(r.memory[2].dissipation.angular_parameter_derivative, 0.0);
        assert!(r.original_relative_parameter_derivative_residual < 1e-10);
        assert!(r.original_relative_mixed_derivative_residual < 1e-10);
    }
}
