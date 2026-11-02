use super::*;
use fs_couple::render::plate::impact::material::heating::SqueezeFilmChannelHeatRoute;
use fs_couple::render::plate::impact::material::heating::liquid_film::{
    LiquidFilmBinding, LiquidFilmChannelBinding, LiquidFilmImpactSystem,
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

fn channel_mechanics(source: &ResolvedLiquidState) -> ImpactSystem {
    let (a, wa) = ImpactBody::free_mass(0.2, 0.0, 0.001).unwrap();
    let (b, wb) = ImpactBody::free_mass(0.4, 0.0, 0.0).unwrap();
    let gap = GapPort {
        reference_m: GAP,
        closure: vec![wa, -wb],
    };
    let sources = [source, source];
    let film = MaterialChannelSqueezeFilm::new(
        &sources,
        vec![FilmCell {
            area_m2: AREA,
            gap: gap.clone(),
        }],
        [LENGTH, 2.0 * LENGTH]
            .into_iter()
            .map(|length_m| FilmChannel {
                from: 0,
                to: None,
                width_m: WIDTH,
                length_m,
                gap: gap.clone(),
            })
            .collect(),
        2,
        FilmLimits {
            minimum_cell_gap_m: 1e-5,
            maximum_gap_m: 5e-4,
            maximum_pressure_pa: 10000.0,
        },
    )
    .unwrap();
    let mechanics = ImpactSystem::new(
        vec![a, b],
        vec![],
        vec![],
        vec![],
        ImpactConfig {
            dt_s: DT,
            max_steps: 512,
            maximum_energy_j: 1.0,
            energy_absolute_tolerance_j: 1e-14,
            energy_relative_tolerance: 1e-7,
            maximum_generalized_force: 1.0,
        },
    )
    .unwrap();
    let (mechanics, _) = film.attach(mechanics).unwrap();
    mechanics
}

fn channel_heat_routes() -> Vec<SqueezeFilmChannelHeatRoute> {
    ["fluid", "sink"]
        .into_iter()
        .enumerate()
        .map(|(channel, key)| SqueezeFilmChannelHeatRoute {
            film: 0,
            channel,
            to_a: key.into(),
            to_b: key.into(),
            fraction_to_a: 1.0,
        })
        .collect()
}

fn channel_heated<'m>(
    source: &ResolvedLiquidState,
    curve: &'m EquilibriumEnthalpyPhaseCurve,
) -> HeatedImpactSystem<'m> {
    HeatedImpactSystem::new_with_channel_losses(
        channel_mechanics(source).prepare_analytic().unwrap(),
        stores(curve),
        vec![],
        vec![],
        channel_heat_routes(),
        vec![],
        1e-13,
    )
    .unwrap()
}

fn channel_bindings(card: &MaterialCard) -> Vec<LiquidFilmChannelBinding<'_>> {
    // Reversed authoring order must not alter native channel identities.
    [(1, "sink"), (0, "fluid")]
        .into_iter()
        .map(|(channel, reservoir)| LiquidFilmChannelBinding {
            film: 0,
            channel,
            reservoir: reservoir.into(),
            law: law(card),
            maximum_relative_viscosity_change: 0.1,
        })
        .collect()
}

fn channel_coupled<'m>(
    card: &'m MaterialCard,
    curve: &'m EquilibriumEnthalpyPhaseCurve,
) -> LiquidFilmImpactSystem<'m> {
    LiquidFilmImpactSystem::new_with_channels(
        channel_heated(law(card).reference_state(), curve),
        vec![],
        channel_bindings(card),
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
                Some(expected.dynamic_viscosity_pa_s())
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

#[test]
fn passage_work_heats_its_local_liquid_and_changes_the_joint_pressure_force() {
    // G2: independent parallel Poiseuille flow/work; G3: unequal local heating.
    let card = card(0.2, 0.001);
    let law = law(&card);
    let c = curve(card.content_hash(), 0.2, 1.0);
    for analytic in [false, true] {
        let mut live = channel_coupled(&card, &c);
        let mut frozen = channel_heated(law.reference_state(), &c);
        live.set_analytic_newton(analytic);
        frozen.set_analytic_newton(analytic);
        assert!(live.bindings().is_empty());
        assert!(live.sources().is_empty());
        assert_eq!(live.channel_bindings()[0].channel, 0);
        assert_eq!(live.channel_bindings()[1].channel, 1);
        let initial = live.system().mechanics().stored_energy_j();
        for _ in 0..128 {
            let old = live.system().mechanics().state().to_vec();
            let mu: Vec<_> = live
                .channel_sources()
                .iter()
                .map(ResolvedLiquidState::dynamic_viscosity_pa_s)
                .collect();
            let received: Vec<_> = live
                .system()
                .reservoirs()
                .iter()
                .map(CaloricReservoir::received_heat_j)
                .collect();
            let frame = live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            frozen
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            let next = live.system().mechanics().state();
            let closure =
                0.5 * ((old[0] + next[0]) / 0.2_f64.sqrt() - (old[2] + next[2]) / 0.4_f64.sqrt());
            let speed =
                0.5 * ((old[1] + next[1]) / 0.2_f64.sqrt() - (old[3] + next[3]) / 0.4_f64.sqrt());
            let conductance: Vec<_> = [LENGTH, 2.0 * LENGTH]
                .iter()
                .zip(&mu)
                .map(|(length, mu)| WIDTH * (GAP - closure).powi(3) / (12.0 * mu * length))
                .collect();
            let pressure = AREA * speed / conductance.iter().sum::<f64>();
            let mut expected_heat = 0.0;
            for i in 0..2 {
                // Each passage receives its actual Q_i * delta_p_i, using the
                // OLD accepted source. Endpoint viscosity never reprices heat.
                let heat = DT * conductance[i] * pressure.powi(2);
                expected_heat += heat;
                let store = &live.system().reservoirs()[i];
                near(store.received_heat_j() - received[i], heat, 2e-20);
                let expected = law
                    .at_temperature(Temperature::new(store.state().temperature_k()))
                    .unwrap();
                assert_eq!(live.channel_sources()[i], expected);
                assert_eq!(
                    live.system()
                        .mechanics()
                        .squeeze_film(0)
                        .unwrap()
                        .channel_viscosities_pa_s()[i],
                    expected.dynamic_viscosity_pa_s()
                );
            }
            near(frame.thermal.impact.deposited_heat_j, expected_heat, 2e-20);
            near(frame.thermal.impact.remaining_mechanical_loss_j, 0.0, 1e-20);
        }
        let local = live.channel_sources();
        assert!(local[0].dynamic_viscosity_pa_s() < local[1].dynamic_viscosity_pa_s());
        assert!(local[1].dynamic_viscosity_pa_s() < 0.01);
        assert!(
            live.system().reservoirs()[0].state().temperature_k()
                > live.system().reservoirs()[1].state().temperature_k()
        );
        assert_eq!(
            live.system()
                .mechanics()
                .squeeze_film(0)
                .unwrap()
                .viscosity_pa_s(),
            None
        );
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
fn equal_channel_temperatures_recover_uniform_feedback_and_rebinding_keeps_time() {
    // G3: splitting one isothermal source into explicit channel rows changes
    // neither the hydraulic solve nor its accepted work, motion or heat history.
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let law = law(&card);
    let mut uniform = LiquidFilmImpactSystem::new(
        channel_heated(law.reference_state(), &c),
        vec![LiquidFilmBinding {
            film: 0,
            reservoir: "fluid".into(),
            law: law.clone(),
            maximum_relative_viscosity_change: 0.1,
        }],
    )
    .unwrap();
    let mut bindings = channel_bindings(&card);
    for binding in &mut bindings {
        binding.reservoir = "fluid".into();
    }
    let mut local = LiquidFilmImpactSystem::new_with_channels(
        channel_heated(law.reference_state(), &c),
        vec![],
        bindings,
    )
    .unwrap();
    for _ in 0..32 {
        uniform
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        local
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(
            uniform.system().mechanics().state(),
            local.system().mechanics().state()
        );
        assert_eq!(uniform.system().reservoirs(), local.system().reservoirs());
        for source in local.channel_sources() {
            assert_eq!(source, &uniform.sources()[0]);
        }
    }
    let before = local.channel_sources().to_vec();
    let (inner, uniforms, channels) = local.into_all_parts();
    let rebound = LiquidFilmImpactSystem::new_with_channels(inner, uniforms, channels).unwrap();
    assert_eq!(rebound.system().mechanics().samples(), 32);
    assert_eq!(
        rebound.system().mechanics().state(),
        uniform.system().mechanics().state()
    );
    assert_eq!(rebound.system().reservoirs(), uniform.system().reservoirs());
    assert_eq!(rebound.channel_sources(), before);
}

fn local_network_bindings() -> [ReservoirCell; 2] {
    [
        ReservoirCell {
            reservoir: "fluid".into(),
            cell: CellId(1),
        },
        ReservoirCell {
            reservoir: "sink".into(),
            cell: CellId(2),
        },
    ]
}

#[test]
fn each_passage_resolves_at_its_actual_conducting_thermal_node_on_the_shared_clock() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let law = law(&card);
    let mut live = channel_coupled(&card, &c);
    let mut thermal = network(&c);
    let mut observer = ();
    let bindings = local_network_bindings();
    let initial = live.system().mechanics().stored_energy_j();
    for _ in 0..64 {
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
        assert_eq!(
            thermal.time_s(),
            live.system().mechanics().samples() as f64 * DT
        );
        for (i, binding) in bindings.iter().enumerate() {
            let h = thermal.cells()[&binding.cell].specific_h_j_kg;
            let state = c.state_at_specific_enthalpy(h).unwrap();
            assert_eq!(live.system().reservoirs()[i].state(), state);
            assert_eq!(
                live.channel_sources()[i],
                law.at_temperature(Temperature::new(state.temperature_k()))
                    .unwrap()
            );
        }
    }
    assert!(
        live.channel_sources()[0].dynamic_viscosity_pa_s()
            < live.channel_sources()[1].dynamic_viscosity_pa_s()
    );
    let stores = live.system().reservoirs();
    // The hotter, shorter passage also conducts heat into its neighbor. Local
    // viscosity follows the full network state, not just its own deposited work.
    assert!(stores[0].storage_change_j() < stores[0].received_heat_j());
    assert!(stores[1].storage_change_j() > stores[1].received_heat_j());
    near(
        live.system().mechanics().stored_energy_j()
            + stores
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>(),
        initial,
        1e-13,
    );
}

#[derive(Clone)]
struct LocalHeater<'g> {
    power_w: f64,
    steps: usize,
    cancel: Option<&'g CancelGate>,
}

impl fs_thermal::enthalpy::DepositionObserver for LocalHeater<'_> {
    type Error = fs_thermal::enthalpy::DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), Self::Error>,
    ) -> Result<std::collections::BTreeMap<CellId, f64>, Self::Error> {
        Ok([(CellId(2), self.power_w)].into_iter().collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &std::collections::BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), Self::Error>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if let Some(gate) = self.cancel {
            gate.request();
        }
        Ok(())
    }
}

#[test]
fn late_local_source_refusal_or_cancellation_rolls_back_every_network_and_mechanical_owner() {
    // G4: the first passage is admissible, but the SECOND fails source support
    // or temporal refinement. No staged source/heat/observer/clock may publish.
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    for (rise_k, cancel) in [(0.002, false), (0.0002, false), (0.0, true)] {
        let mut live = channel_coupled(&card, &c);
        let mut thermal = network(&c);
        let gate = CancelGate::new_clock_free();
        let mut observer = LocalHeater {
            power_w: MASS * 0.2 * rise_k / DT,
            steps: 0,
            cancel: cancel.then_some(&gate),
        };
        let motion = live.system().mechanics().state().to_vec();
        let stores = live.system().reservoirs().to_vec();
        let sources = live.channel_sources().to_vec();
        let cells = thermal.cells().clone();
        let result = live.step_with_thermal_network(
            &mut thermal,
            &mut observer,
            &local_network_bindings(),
            StepConfig {
                max_dt_s: DT / 4.0,
                energy_tolerance_j: 1e-13,
            },
            &[0.0; 2],
            &gate,
            |_, _| Ok(None),
        );
        let error = result.unwrap_err();
        use fs_couple::render::plate::impact::ImpactError;
        if cancel {
            assert!(matches!(error, ImpactError::Cancelled));
        } else if rise_k > 0.001 {
            assert!(
                matches!(error, ImpactError::Owner(ref source) if source.starts_with("Query {")),
                "unexpected refusal: {error}"
            );
        } else {
            assert!(
                matches!(error, ImpactError::Invalid(reason)
                if reason.contains("viscosity changed too far")),
                "unexpected refusal: {error}"
            );
        }
        assert_eq!(live.system().mechanics().state(), motion);
        assert_eq!(live.system().mechanics().samples(), 0);
        assert_eq!(live.system().reservoirs(), stores);
        assert_eq!(live.channel_sources(), sources);
        assert_eq!(
            live.system()
                .mechanics()
                .squeeze_film(0)
                .unwrap()
                .channel_viscosities_pa_s(),
            &[0.01, 0.01]
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(observer.steps, 0);
        // A rejected trial leaves no hidden state: retry agrees with a fresh
        // owner when the authored heater/cancellation is removed.
        let mut fresh = channel_coupled(&card, &c);
        let mut fresh_thermal = network(&c);
        observer.power_w = 0.0;
        observer.cancel = None;
        live.step_with_thermal_network(
            &mut thermal,
            &mut observer,
            &local_network_bindings(),
            StepConfig {
                max_dt_s: DT / 4.0,
                energy_tolerance_j: 1e-13,
            },
            &[0.0; 2],
            &CancelGate::new_clock_free(),
            |_, _| Ok(None),
        )
        .unwrap();
        fresh
            .step_with_thermal_network(
                &mut fresh_thermal,
                &mut (),
                &local_network_bindings(),
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-13,
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None),
            )
            .unwrap();
        assert_eq!(
            live.system().mechanics().state(),
            fresh.system().mechanics().state()
        );
        assert_eq!(live.system().reservoirs(), fresh.system().reservoirs());
        assert_eq!(live.channel_sources(), fresh.channel_sources());
        assert_eq!(thermal.time_s(), fresh_thermal.time_s());
        assert_eq!(thermal.cells(), fresh_thermal.cells());
        assert!(observer.steps >= 4);
    }
}

#[test]
fn local_film_sources_require_complete_unique_original_channel_coverage() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let original = channel_bindings(&card);
    let mut outside = original.clone();
    outside[0].channel = 2;
    for channels in [
        vec![],
        vec![original[0].clone()],
        vec![original[0].clone(), original[0].clone()],
        outside,
    ] {
        assert!(
            LiquidFilmImpactSystem::new_with_channels(
                channel_heated(law(&card).reference_state(), &c),
                vec![],
                channels,
            )
            .is_err()
        );
    }
    assert!(
        LiquidFilmImpactSystem::new_with_channels(
            channel_heated(law(&card).reference_state(), &c),
            vec![LiquidFilmBinding {
                film: 0,
                reservoir: "fluid".into(),
                law: law(&card),
                maximum_relative_viscosity_change: 0.1
            }],
            original,
        )
        .is_err()
    );
}

#[test]
fn channel_heat_routes_refuse_double_charging_and_invalid_destinations() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let law = law(&card);
    let selected = channel_heat_routes()[0].clone();
    let refuses = |films, channels| {
        HeatedImpactSystem::new_with_channel_losses(
            channel_mechanics(law.reference_state())
                .prepare_analytic()
                .unwrap(),
            stores(&c),
            vec![],
            films,
            channels,
            vec![],
            1e-13,
        )
        .is_err()
    };
    assert!(refuses(vec![], vec![selected.clone(), selected.clone()]));
    assert!(refuses(vec![route()], vec![selected.clone()]));
    for invalid in [
        SqueezeFilmChannelHeatRoute {
            channel: 2,
            ..selected.clone()
        },
        SqueezeFilmChannelHeatRoute {
            film: 1,
            ..selected.clone()
        },
        SqueezeFilmChannelHeatRoute {
            to_a: "absent".into(),
            ..selected.clone()
        },
        SqueezeFilmChannelHeatRoute {
            to_b: "absent".into(),
            ..selected.clone()
        },
        SqueezeFilmChannelHeatRoute {
            fraction_to_a: -0.1,
            ..selected.clone()
        },
        SqueezeFilmChannelHeatRoute {
            fraction_to_a: 1.1,
            ..selected.clone()
        },
        SqueezeFilmChannelHeatRoute {
            fraction_to_a: f64::NAN,
            ..selected
        },
    ] {
        assert!(refuses(vec![], vec![invalid]));
    }
}

#[test]
fn selective_channel_heat_retains_unrouted_work_and_identical_mechanics() {
    // G2/G3: thermal routing changes the energy boundary, never the hydraulic
    // pressure solve or mechanical trajectory. The omitted passage still resists.
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let law = law(&card);
    let source = law.reference_state();
    let mut all = channel_heated(source, &c);
    let mut selective = HeatedImpactSystem::new_with_channel_losses(
        channel_mechanics(source).prepare_analytic().unwrap(),
        stores(&c),
        vec![],
        vec![],
        vec![channel_heat_routes()[0].clone()],
        vec![],
        1e-13,
    )
    .unwrap();
    let mut plain = channel_mechanics(source).prepare_analytic().unwrap();
    let initial = plain.stored_energy_j();
    let mut retained = 0.0;
    let mut outside = 0.0;
    for _ in 0..64 {
        let old = plain.state().to_vec();
        let partial = selective
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        let complete = all.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        let mechanical = plain
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        let next = plain.state();
        assert_eq!(selective.mechanics().state(), next);
        assert_eq!(all.mechanics().state(), next);
        let closure =
            0.5 * ((old[0] + next[0]) / 0.2_f64.sqrt() - (old[2] + next[2]) / 0.4_f64.sqrt());
        let speed =
            0.5 * ((old[1] + next[1]) / 0.2_f64.sqrt() - (old[3] + next[3]) / 0.4_f64.sqrt());
        let g0 =
            WIDTH * (GAP - closure).powi(3) / (12.0 * source.dynamic_viscosity_pa_s() * LENGTH);
        let g1 = 0.5 * g0;
        let pressure = AREA * speed / (g0 + g1);
        let local_heat = DT * g0 * pressure.powi(2);
        let unrouted_heat = DT * g1 * pressure.powi(2);
        near(partial.deposited_heat_j, local_heat, 2e-20);
        near(partial.remaining_mechanical_loss_j, unrouted_heat, 2e-20);
        near(complete.deposited_heat_j, local_heat + unrouted_heat, 2e-20);
        near(complete.remaining_mechanical_loss_j, 0.0, 2e-20);
        near(
            mechanical.dissipated_energy_j,
            local_heat + unrouted_heat,
            2e-20,
        );
        retained += local_heat;
        outside += partial.remaining_mechanical_loss_j;
    }
    assert!(outside > 0.0);
    near(selective.reservoirs()[0].received_heat_j(), retained, 2e-20);
    assert_eq!(selective.reservoirs()[1].received_heat_j(), 0.0);
    assert_eq!(selective.reservoirs()[1].storage_change_j(), 0.0);
    near(
        selective.mechanics().stored_energy_j()
            + selective
                .reservoirs()
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>()
            + outside,
        initial,
        1e-13,
    );
    near(
        all.mechanics().stored_energy_j()
            + all
                .reservoirs()
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>(),
        initial,
        1e-13,
    );
}
