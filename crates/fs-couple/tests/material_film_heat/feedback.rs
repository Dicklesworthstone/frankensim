use super::*;
use fs_couple::render::plate::impact::material::heating::liquid_film::{
    LiquidFilmBinding, LiquidFilmImpactSystem,
};
use fs_couple::render::plate::impact::material::heating::transport::CoupledThermalStep;
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::UniformEnthalpyStep;
use fs_qty::Temperature;

fn card(cp: f64, hot_mu: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims, quantity) in [
        (
            "density",
            1000.0,
            Dims([-3, 1, 0, 0, 0, 0]),
            QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
        ),
        (
            "specific_heat_capacity",
            cp,
            Dims([2, 0, -2, -1, 0, 0]),
            QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0])),
        ),
        (
            "dynamic_viscosity",
            0.01,
            Dims([-1, 1, -1, 0, 0, 0]),
            QuantitySpec::dimensional(Dims([-1, 1, -1, 0, 0, 0])),
        ),
        (
            "thermal_conductivity",
            0.6,
            Dims([1, 1, -3, -1, 0, 0]),
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
        ),
    ] {
        let (value, interpolation) = if name == "dynamic_viscosity" {
            (
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, value), (300.001, hot_mu)],
                    dims,
                },
                InterpolationPolicy::LinearInside,
            )
        } else {
            (
                PropertyValue::Scalar { value, dims },
                InterpolationPolicy::ConstantWithinValidity,
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value,
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 300.001)
                    .with("p", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic temperature-sensitive film; not a measured lubricant".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("p", 100000.0)
            .unwrap(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn coupled<'m>(
    card: &'m MaterialCard,
    curve: &'m EquilibriumEnthalpyPhaseCurve,
) -> LiquidFilmImpactSystem<'m> {
    let law = law(card);
    LiquidFilmImpactSystem::new(
        heated(law.reference_state(), curve, 0.0),
        vec![LiquidFilmBinding {
            film: 0,
            reservoir: "fluid".into(),
            law,
            maximum_relative_viscosity_change: 0.1,
        }],
    )
    .unwrap()
}

#[test]
fn dissipative_work_heats_the_liquid_and_changes_the_next_actual_joint_force() {
    let card = card(0.2, 0.001);
    let law = law(&card);
    let c = curve(card.content_hash(), 0.2, 1.0);
    for analytic in [false, true] {
        let mut live = coupled(&card, &c);
        live.set_analytic_newton(analytic);
        let mut frozen = heated(law.reference_state(), &c, 0.0);
        frozen.set_analytic_newton(analytic);
        let initial = live.system().mechanics().stored_energy_j();
        for _ in 0..128 {
            let old = live.system().mechanics().state().to_vec();
            let mu = live.sources()[0].dynamic_viscosity_pa_s();
            let frame = live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            frozen
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            let next = live.system().mechanics().state();
            let closure =
                0.5 * ((old[0] + next[0]) / 0.2_f64.sqrt() - (old[2] + next[2]) / 0.4_f64.sqrt());
            let speed =
                0.5 * ((old[1] + next[1]) / 0.2_f64.sqrt() - (old[3] + next[3]) / 0.4_f64.sqrt());
            let conductance = WIDTH * (GAP - closure).powi(3) / (12.0 * mu * LENGTH);
            // The completed step used the OLD accepted coefficient. Repricing
            // its heat with the hotter endpoint viscosity would violate this.
            near(
                frame.thermal.impact.deposited_heat_j,
                DT * (AREA * speed).powi(2) / conductance,
                1e-20,
            );
            let temperature = live.system().reservoirs()[0].state().temperature_k();
            let expected = law.at_temperature(Temperature::new(temperature)).unwrap();
            assert_eq!(live.sources()[0], expected);
            assert_eq!(
                live.system()
                    .mechanics()
                    .squeeze_film(0)
                    .unwrap()
                    .viscosity_pa_s(),
                expected.dynamic_viscosity_pa_s()
            );
        }
        assert!(live.sources()[0].dynamic_viscosity_pa_s() < 0.009);
        assert!(live.system().mechanics().state()[3] < frozen.mechanics().state()[3]);
        let motion = live.system().mechanics();
        near(
            0.2_f64.sqrt() * motion.state()[1] + 0.4_f64.sqrt() * motion.state()[3],
            0.0002,
            1e-12,
        );
        near(
            motion.stored_energy_j()
                + live
                    .system()
                    .reservoirs()
                    .iter()
                    .map(CaloricReservoir::storage_change_j)
                    .sum::<f64>(),
            initial,
            1e-13,
        );
    }
}

#[test]
fn heat_capacity_is_data_and_feeds_back_into_motion_through_viscosity() {
    let a = card(0.2, 0.001);
    let b = card(0.4, 0.001);
    let ca = curve(a.content_hash(), 0.2, 1.0);
    let cb = curve(b.content_hash(), 0.4, 1.0);
    let mut sa = coupled(&a, &ca);
    let mut sb = coupled(&b, &cb);
    for _ in 0..128 {
        sa.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        sb.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
    }
    assert!(
        sa.system().reservoirs()[0].state().temperature_k()
            > sb.system().reservoirs()[0].state().temperature_k()
    );
    assert!(sa.sources()[0].dynamic_viscosity_pa_s() < sb.sources()[0].dynamic_viscosity_pa_s());
    assert!(sa.system().mechanics().state()[3] < sb.system().mechanics().state()[3]);
}

#[test]
fn live_material_feedback_composes_with_the_existing_spatial_heat_network() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let mut live = coupled(&card, &c);
    let mut thermal = network(&c);
    let mut observer = ();
    let bindings = [
        ReservoirCell {
            reservoir: "fluid".into(),
            cell: CellId(1),
        },
        ReservoirCell {
            reservoir: "sink".into(),
            cell: CellId(2),
        },
    ];
    let initial = live.system().mechanics().stored_energy_j();
    for _ in 0..128 {
        let frame = live
            .step_with_thermal_network(
                &mut thermal,
                &mut observer,
                &bindings,
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-13,
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None),
            )
            .unwrap();
        assert!(frame.network.substeps >= 4);
        near(frame.external_heat_j, 0.0, 1e-19);
        for (i, b) in bindings.iter().enumerate() {
            assert_eq!(
                live.system().reservoirs()[i]
                    .state()
                    .specific_enthalpy_j_kg(),
                thermal.cells()[&b.cell].specific_h_j_kg
            );
        }
    }
    assert!(live.system().reservoirs()[1].state().temperature_k() > 300.0);
    assert_eq!(live.system().reservoirs()[1].received_heat_j(), 0.0);
    assert!(live.sources()[0].dynamic_viscosity_pa_s() < 0.01);
    near(
        live.system().mechanics().stored_energy_j()
            + live
                .system()
                .reservoirs()
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>(),
        initial,
        1e-13,
    );
}

#[test]
fn endpoint_source_refusal_and_late_cancellation_preserve_every_accepted_owner() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    for cancel in [false, true] {
        let mut live = coupled(&card, &c);
        let motion = live.system().mechanics().state().to_vec();
        let stores = live.system().reservoirs().to_vec();
        let source = live.sources().to_vec();
        let mu = live
            .system()
            .mechanics()
            .squeeze_film(0)
            .unwrap()
            .viscosity_pa_s();
        let gate = CancelGate::new_clock_free();
        let result = live.step_with_coupled_transport(
            &[0.0; 2],
            &gate,
            |_, _| Ok(None),
            |trial| {
                let transfers = trial
                    .reservoirs()
                    .iter()
                    .map(|r| {
                        let i = r.thermal;
                        let heat = if !cancel && r.key == "fluid" {
                            i.mass_kg * 0.2 * 0.002
                        } else {
                            0.0
                        };
                        Ok(UniformEnthalpyStep {
                            state: i
                                .curve
                                .advance_specific_energy(
                                    i.initial,
                                    (i.internal_heat_j + heat) / i.mass_kg,
                                )
                                .unwrap(),
                            external_heat_j: heat,
                            energy_residual_tolerance_j: 1e-13,
                            report: (),
                        })
                    })
                    .collect::<Result<Vec<_>, fs_couple::render::plate::impact::ImpactError>>()?;
                let boundary_heat_j = transfers.iter().map(|t| t.external_heat_j).sum();
                if cancel {
                    gate.request();
                }
                Ok(CoupledThermalStep {
                    transfers,
                    boundary_heat_j,
                })
            },
        );
        assert!(result.is_err());
        assert_eq!(live.system().mechanics().state(), motion);
        assert_eq!(live.system().reservoirs(), stores);
        assert_eq!(live.sources(), source);
        assert_eq!(live.system().mechanics().samples(), 0);
        assert_eq!(
            live.system()
                .mechanics()
                .squeeze_film(0)
                .unwrap()
                .viscosity_pa_s(),
            mu
        );
    }
}

#[test]
fn mechanical_heat_and_conservative_advection_share_the_original_flux_and_cfl() {
    use fs_thermal::enthalpy::advection::{Advected, MassBoundary, MassFlow};
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let mut live = coupled(&card, &c);
    let mut thermal = network(&c);
    // Fast prescribed circulation forces the ADVECTIVE timestep below the
    // authored mechanical period. Inflow and outflow balance each fixed mass;
    // exported enthalpy belongs to the outside reservoir, never a heat repair.
    let rate = MASS / (DT / 8.0);
    let mut observer = Advected::new(
        (),
        vec![MassFlow {
            id: 10,
            a: CellId(1),
            b: CellId(2),
            kg_s: rate,
        }],
        3,
    )
    .unwrap()
    .with_boundaries(
        vec![
            MassBoundary::Inlet {
                id: 11,
                cell: CellId(1),
                kg_s: rate,
                material_id: MaterialId(1),
                specific_h_j_kg: 0.0,
            },
            MassBoundary::Outlet {
                id: 12,
                cell: CellId(2),
                kg_s: rate,
            },
        ],
        3,
    )
    .unwrap();
    let bindings = [
        ReservoirCell {
            reservoir: "fluid".into(),
            cell: CellId(1),
        },
        ReservoirCell {
            reservoir: "sink".into(),
            cell: CellId(2),
        },
    ];
    let initial = live.system().mechanics().stored_energy_j();
    let mut external = 0.0;
    for _ in 0..32 {
        let frame = live
            .step_with_thermal_network(
                &mut thermal,
                &mut observer,
                &bindings,
                StepConfig {
                    max_dt_s: DT,
                    energy_tolerance_j: 1e-13,
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None),
            )
            .unwrap();
        let transport = observer.audit();
        assert_eq!(transport.substeps, frame.network.substeps);
        assert!(transport.substeps >= 10);
        assert!(transport.maximum_courant <= 0.8);
        assert!(transport.absolute_enthalpy_transfer_j > 0.0);
        assert!(transport.external_enthalpy_j < 0.0);
        near(transport.inlet_mass_kg, rate * DT, 1e-16);
        near(transport.outlet_mass_kg, transport.inlet_mass_kg, 1e-16);
        near(frame.external_heat_j, transport.external_enthalpy_j, 1e-18);
        external += frame.external_heat_j;
        for (i, binding) in bindings.iter().enumerate() {
            assert_eq!(
                live.system().reservoirs()[i]
                    .state()
                    .specific_enthalpy_j_kg(),
                thermal.cells()[&binding.cell].specific_h_j_kg
            );
        }
    }
    near(
        live.system().mechanics().stored_energy_j()
            + live
                .system()
                .reservoirs()
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>()
            - external,
        initial,
        1e-13,
    );
}
