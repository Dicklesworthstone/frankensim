//! Independent endpoint equations for synthetic reversible isomers. These
//! check implementation and coupling, not empirical combustion kinetics.

use super::*;
use crate::cavity_phs::piston::PistonSupport;
use crate::cavity_phs::piston::reactive::coupled::{
    CoupledReactivePistonError, CoupledReactivePistonSolve,
};
use crate::cavity_phs::piston::reactive::implicit::coupled::ImplicitMechanismSolve;
use crate::cavity_phs::piston::reactive::material_rate::{
    ACTIVITY_RATE_SCALE_DIMS, ACTIVITY_RATE_SCALE_PROPERTY, MaterialRateScale,
};

pub(super) fn temperature_reaction(release: f64, scale: f64) -> GasReaction {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, scale / 3.0), (1800.0, 3.0 * scale)],
                dims: ACTIVITY_RATE_SCALE_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                200.0,
                1800.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-endpoint-isomer-kinetics".into(),
            phase: "gas".into(),
            process: "linear-temperature-scale".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let law = reaction(release, 3.5, scale)
        .with_material_rate(
            MaterialRateScale::new(
                material,
                MaterialPropertySelection::SingleClaimOnly,
                Pressure::new(100000.0),
            )
            .unwrap(),
        )
        .unwrap();
    let policy = ImplicitMechanismSolve {
        extent_brackets_mol: law
            .channels()
            .map(|(kinetics, _)| {
                (
                    kinetics
                        .reverse_rate_closure()
                        .equilibrium()
                        .reaction()
                        .clone(),
                    [-1.0, 1.0],
                )
            })
            .collect(),
        extent_tolerance_mol: 1e-12,
        interior_fraction: 0.05,
        jacobian_relative_step: 1e-5,
        maximum_iterations: 32,
        maximum_evaluations: 1024,
        newton: fs_solver::NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            max_linear_cycles: 8,
            forcing_maximum: 1e-6,
            ..Default::default()
        },
    };
    law.with_implicit_mechanism(policy).unwrap()
}

pub(super) fn coupled_controls() -> CoupledReactivePistonSolve {
    CoupledReactivePistonSolve {
        piston: mechanical_controls(),
        pressure_interval_pa: [50000.0, 300000.0],
        pressure_tolerance_pa: 1e-6,
    }
}

fn chemistry() -> GasReactionSolve {
    GasReactionSolve {
        maximum_consumed_fraction: 0.8,
        ..chemical_controls()
    }
}

fn capacity(amounts: &[f64]) -> f64 {
    R * (2.5 * amounts[0] + 2.5 * amounts[1] + 1.5 * amounts[2])
}

#[test]
fn fixed_wall_heat_changes_the_same_endpoint_kinetics_and_support_reaction() {
    let scale = 100000.0;
    let law = temperature_reaction(0.0, scale);
    let piston = GasPiston::new(gas(&law, [0.75, 0.0, 0.25]), 0.1, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let original = piston.clone();
    let input = PistonInput {
        heat_j: 30.0,
        ..drive()
    };
    let controls = CoupledReactivePistonSolve {
        piston: PistonSolve {
            volume_bracket_m3: [piston.gas().volume_m3(); 2],
            ..mechanical_controls()
        },
        ..coupled_controls()
    };
    let result = piston
        .coupled_reactive_candidate(&law, input, chemistry(), controls, || false)
        .unwrap();
    let before = law.inventory(piston.gas()).unwrap();
    let after = law.inventory(result.piston().gas()).unwrap();
    let initial_t = piston.gas().state().state().temperature;
    let expected_t = initial_t + input.heat_j / capacity(&before);
    let z = input.duration_s * scale * (expected_t / 600.0) * R * expected_t / 100000.0;
    let extent = z * (before[0] - before[1]) / (1.0 + 2.0 * z);
    let expected_p = before.iter().sum::<f64>() * R * expected_t / piston.gas().volume_m3();
    near(
        result.piston().gas().state().state().temperature,
        expected_t,
        1e-7,
    );
    near(
        result.piston().gas().state().state().pressure,
        expected_p,
        1e-6,
    );
    near(after[0], before[0] - extent, 2e-12);
    near(after[1], before[1] + extent, 2e-12);
    near(after[2], before[2], 1e-14);
    let report = result.step();
    near(report.reaction.reaction.extent_mol, extent, 2e-12);
    let source = report.reaction.reaction.rate_material.as_ref().unwrap();
    near(
        source
            .property(ACTIVITY_RATE_SCALE_PROPERTY)
            .unwrap()
            .value_si(),
        scale * expected_t / 600.0,
        1e-5,
    );
    near(
        *source.query().axes().get("temperature").unwrap(),
        expected_t,
        1e-7,
    );
    near(
        report.support_impulse_n_s,
        -input.duration_s * piston.area_m2() * (expected_p - input.exterior_pressure_pa),
        1e-10,
    );
    assert_eq!(result.piston().support(), PistonSupport::Fixed);
    assert_eq!(result.piston().gas().volume_m3(), piston.gas().volume_m3());
    assert_eq!(result.piston().gas().mass_kg(), piston.gas().mass_kg());
    assert_eq!(report.displacement_m, 0.0);
    assert_eq!(report.support_work_j, 0.0);
    assert_eq!(report.kinetic_energy_change_j, 0.0);
    assert!(report.total_energy_residual_j.abs() <= controls.piston.energy_tolerance_j);
    assert!(report.pressure_residual_pa.abs() <= controls.pressure_tolerance_pa);
    let staged = piston
        .reactive_candidate(&law, input, chemistry(), controls.piston, || false)
        .unwrap();
    assert!(report.reaction.reaction.extent_mol > 1.05 * staged.step().reaction.extent_mol);
    assert_eq!(piston, original);

    // An equilibrium chemical residual is zero even while the thermal state
    // changes. The endpoint owner must still apply the heat and resolve k(T).
    let equilibrium = GasPiston::new(gas(&law, [0.5, 0.5, 0.0]), 0.1, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let equilibrium_amounts = law.inventory(equilibrium.gas()).unwrap();
    let heated = equilibrium
        .coupled_reactive_candidate(&law, input, chemistry(), controls, || false)
        .unwrap();
    near(heated.step().reaction.reaction.extent_mol, 0.0, 1e-14);
    for (actual, expected) in law
        .inventory(heated.piston().gas())
        .unwrap()
        .iter()
        .zip(&equilibrium_amounts)
    {
        near(*actual, *expected, 1e-14);
    }
    near(
        heated.piston().gas().state().state().temperature,
        600.0 + input.heat_j / capacity(&equilibrium_amounts),
        1e-7,
    );
}

#[test]
fn exothermic_compression_continues_with_one_pressure_energy_and_endpoint_rate() {
    let release = 4000.0;
    let scale = 100000.0;
    let law = temperature_reaction(release, scale);
    let mut piston = GasPiston::new(gas(&law, [0.75, 0.0, 0.25]), 0.1, 0.01, -10.0).unwrap();
    let input = PistonInput {
        heat_j: 10.0,
        external_force_n: 5.0,
        ..drive()
    };
    let controls = coupled_controls();
    for _ in 0..2 {
        let original = piston.clone();
        let before = law.inventory(piston.gas()).unwrap();
        let candidate = piston
            .coupled_reactive_candidate(&law, input, chemistry(), controls, || false)
            .unwrap();
        let next = candidate.piston();
        let after = law.inventory(next.gas()).unwrap();
        let report = candidate.step();
        let temperature = next.gas().state().state().temperature;
        let pressure = next.gas().state().state().pressure;
        let extent = after[1] - before[1];
        let dv = next.gas().volume_m3() - piston.gas().volume_m3();
        let dx = dv / piston.area_m2();
        // NASA-9 fixtures have equal cp and entropy offsets. Their independent
        // detailed-balance factor is K(T)=exp(release/(R*T)).
        let endpoint_rate_extent =
            input.duration_s * scale * (temperature / 600.0) * R * temperature / 100000.0
                * (after[0] - after[1] * (-release / (R * temperature)).exp());
        near(extent, endpoint_rate_extent, 2e-12);
        near(after[0], before[0] - extent, 1e-14);
        near(after[2], before[2], 1e-14);
        near(
            pressure,
            after.iter().sum::<f64>() * R * temperature / next.gas().volume_m3(),
            1e-6,
        );
        let u0 = capacity(&before) * piston.gas().state().state().temperature - release * before[1];
        let u1 = capacity(&after) * temperature - release * after[1];
        near(u1 - u0 + pressure * dv, input.heat_j, 2e-7);
        near(
            piston.mass_kg() * (next.velocity_m_s() - piston.velocity_m_s()),
            input.duration_s
                * (piston.area_m2() * (pressure - input.exterior_pressure_pa)
                    + input.external_force_n),
            1e-10,
        );
        near(
            dx,
            0.5 * input.duration_s * (piston.velocity_m_s() + next.velocity_m_s()),
            1e-14,
        );
        let dk =
            0.5 * piston.mass_kg() * (next.velocity_m_s().powi(2) - piston.velocity_m_s().powi(2));
        near(
            u1 - u0 + dk + input.exterior_pressure_pa * dv - input.external_force_n * dx,
            input.heat_j,
            2e-7,
        );
        near(report.reaction.boundary.work_by_gas_j, pressure * dv, 1e-10);
        near(report.kinetic_energy_change_j, dk, 1e-12);
        near(report.reaction.reaction.extent_mol, extent, 1e-14);
        assert!(report.reaction.boundary.work_by_gas_j < 0.0);
        assert_eq!(report.support_impulse_n_s, 0.0);
        assert_eq!(next.gas().mass_kg(), piston.gas().mass_kg());
        assert!(report.total_energy_residual_j.abs() <= controls.piston.energy_tolerance_j);
        assert_eq!(piston, original);
        piston = candidate.into_parts().0;
    }
}

#[test]
fn pressure_domain_and_late_cancellation_leave_every_accepted_owner_unchanged() {
    let law = temperature_reaction(4000.0, 100000.0);
    let piston = GasPiston::new(gas(&law, [0.75, 0.0, 0.25]), 0.1, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let original = piston.clone();
    let input = PistonInput {
        heat_j: 20.0,
        ..drive()
    };
    let controls = coupled_controls();
    let mut calls = 0;
    let accepted = piston
        .coupled_reactive_candidate(&law, input, chemistry(), controls, || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(accepted.step().reaction.reaction.extent_mol > 0.0);
    assert!(calls > 1);
    for cancel_at in [1, calls] {
        let mut observed = 0;
        assert!(matches!(
            piston.coupled_reactive_candidate(&law, input, chemistry(), controls, || {
                observed += 1;
                observed == cancel_at
            },),
            Err(CoupledReactivePistonError::Cancelled)
        ));
        assert_eq!(piston, original);
        assert_eq!(piston.ports.next_tick(), original.ports.next_tick());
    }
    let narrow = CoupledReactivePistonSolve {
        pressure_interval_pa: [99000.0, 101000.0],
        ..controls
    };
    assert!(matches!(
        piston.coupled_reactive_candidate(&law, input, chemistry(), narrow, || false,),
        Err(CoupledReactivePistonError::OutsidePressureBracket { .. })
    ));
    let restricted = GasReactionSolve {
        gas: GasEnergySolve {
            maximum_temperature_k: 600.01,
            ..chemistry().gas
        },
        ..chemistry()
    };
    let restricted_controls = CoupledReactivePistonSolve {
        piston: PistonSolve {
            gas: restricted.gas,
            ..controls.piston
        },
        ..controls
    };
    assert!(matches!(
        piston.coupled_reactive_candidate(&law, input, restricted, restricted_controls, || false,),
        Err(CoupledReactivePistonError::Reaction(GasReactionError::Gas(
            _
        )))
    ));
    assert_eq!(piston, original);
    assert_eq!(piston.ports.next_tick(), original.ports.next_tick());
}
