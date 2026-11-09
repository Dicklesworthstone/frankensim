use super::*;
use crate::lattice3::open_boundaries::sgs::{OpenThermalSgsConfig3, OpenThermalSgsError3};
use crate::sgs3::{SgsConfig, SgsLattice3, SgsModel, SgsUse};

fn closure() -> OpenThermalSgsConfig3 {
    OpenThermalSgsConfig3 {
        closure: SgsConfig {
            model: SgsModel::Smagorinsky { cs: 0.24 },
            filter_width_lu: 1.3,
            maximum_eddy_to_molecular: 20.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.73,
        maximum_diffusivity_to_molecular: 20.0,
    }
}

fn source_and_native_sgs_oracle<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
) {
    let mut actual = make::<L, Q>(
        collision,
        [6, 4, 3],
        native_channel(),
        channel_faces(),
        0.5,
        BETA,
        0.02,
    );
    let mut explicit = make::<L, Q>(
        collision,
        [6, 4, 3],
        native_channel(),
        channel_faces(),
        0.5,
        BETA,
        0.02,
    );
    let cells = actual.thermal().populations().len();
    let heat: Vec<_> = (0..cells)
        .map(|cell| 0.0004 * (1 + cell % 4) as f64)
        .collect();
    for owner in [&mut actual, &mut explicit] {
        advance(owner, &heat);
    }
    let reference = material_at(310.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.0).unwrap();
    let mut material_workspace = actual
        .material_transport_workspace(&scaling, material_config())
        .unwrap();
    let mut sgs_workspace = actual.thermal_sgs_workspace(1024 * 1024).unwrap();
    let mut native = explicit.flow().sgs_workspace(1024 * 1024).unwrap();
    let config = closure();
    let reads = Arc::new(AtomicUsize::new(0));
    let mut previous_molecular = Vec::<f64>::new();
    let mut maximum_molecular_change = 0.0_f64;
    for epoch in 0..3 {
        let opening_k: Vec<_> = (0..cells)
            .map(|cell| 300.0 + 20.0 * explicit.thermal().temperature(cell).unwrap())
            .collect();
        let mut collisions = Vec::new();
        let mut molecular_tau = Vec::new();
        for &temperature in &opening_k {
            let source = material_at(temperature);
            let units = (DT / DX) / DX;
            // Independent source-to-native conversion uses fixed storage, not
            // the source's admitted density/Cp variations or last SGS output.
            let tau = 0.5 + (source.dynamic_viscosity_pa_s() / PHYSICAL_RHO) * units / (1.0 / 3.0);
            collisions.push(L::with_shear_relaxation(collision, tau).unwrap());
            molecular_tau.push(
                0.5 + ((source.thermal_conductivity_w_m_k() / PHYSICAL_RHO) / PHYSICAL_CP) * units
                    / 0.25,
            );
        }
        let fit = explicit
            .flow
            .inner
            .prepare_open_sgs(
                config.closure,
                Some(&collisions),
                &mut native,
                Some(&explicit.flow.forces),
                || Ok(()),
            )
            .unwrap();
        // Only the original momentum closure is used as oracle. Thermal
        // augmentation is independently formed from its returned nu_t.
        let effective_tau: Vec<_> = molecular_tau
            .iter()
            .zip(native.eddy_viscosities())
            .map(|(&tau, &eddy)| {
                if eddy == 0.0 {
                    tau
                } else {
                    tau + 4.0 * (eddy / config.turbulent_prandtl)
                }
            })
            .collect();
        let mut background: Vec<_> = (0..cells)
            .map(|cell| {
                [
                    BACKGROUND[0] + (epoch + 1) as f64 * 1e-5 * (1 + cell % 3) as f64,
                    BACKGROUND[1] - 2e-5 * (1 + cell % 4) as f64,
                    BACKGROUND[2] + 1e-5 * (1 + cell % 5) as f64,
                ]
            })
            .collect();
        let mut reference_background = background.clone();
        let expected = explicit
            .advance_with_diffusion_poll(
                Some(&mut reference_background),
                Some(native.collisions()),
                &heat,
                Some(&effective_tau),
                None,
                [0.0; 3],
                &mut || Ok(()),
            )
            .unwrap();
        let mut ordinal = 0;
        reads.store(0, Ordering::SeqCst);
        let report = actual
            .advance_material_sgs_with_poll(
                Some(&mut background),
                config,
                None,
                &heat,
                &mut material_workspace,
                &mut sgs_workspace,
                [0.0; 3],
                &mut |cell, temperature_k, retained| {
                    assert_eq!(*retained, reference);
                    let source = TracedMaterial {
                        cell,
                        temperature_k,
                        ordinal,
                        state: material_at(temperature_k),
                        reads: Arc::clone(&reads),
                    };
                    ordinal += 1;
                    Ok(source)
                },
                &mut || Ok(()),
            )
            .unwrap();
        assert_eq!(ordinal, 2 * cells);
        assert_eq!(reads.load(Ordering::SeqCst), 2 * cells);
        assert_eq!(report.material.transport, expected);
        assert_eq!(report.sgs, fit);
        assert_eq!(accepted(&actual), accepted(&explicit));
        assert_eq!(background, reference_background);
        assert_eq!(material_workspace.thermal_relaxation(), molecular_tau);
        assert_eq!(sgs_workspace.thermal_relaxation(), effective_tau);
        assert!(report.sgs.maximum_eddy_viscosity_lu > 0.0);
        assert!(effective_tau.iter().zip(&molecular_tau).any(|(a, b)| a > b));
        near(
            report.maximum_eddy_diffusivity_lu,
            fit.maximum_eddy_viscosity_lu / config.turbulent_prandtl,
            2e-16,
        );
        for (query_epoch, samples) in [&report.material.applied, &report.material.endpoint]
            .into_iter()
            .enumerate()
        {
            assert_eq!(samples.len(), cells);
            for sample in samples {
                let source = sample.material.material();
                assert_eq!(source.ordinal, query_epoch * cells + sample.cell);
                assert_eq!(source.cell, sample.cell);
                let expected_k = if query_epoch == 0 {
                    opening_k[sample.cell]
                } else {
                    300.0 + 20.0 * actual.thermal().temperature(sample.cell).unwrap()
                };
                assert_eq!(sample.temperature_k.to_bits(), expected_k.to_bits());
                assert_eq!(source.temperature_k.to_bits(), expected_k.to_bits());
                if query_epoch == 0 {
                    assert_eq!(
                        sample.thermal_relaxation.to_bits(),
                        molecular_tau[sample.cell].to_bits()
                    );
                }
            }
        }
        for (&before, &after) in previous_molecular.iter().zip(&molecular_tau) {
            maximum_molecular_change = maximum_molecular_change.max((after - before).abs());
        }
        previous_molecular = molecular_tau;
    }
    assert!(
        maximum_molecular_change > 1e-8,
        "{maximum_molecular_change:e}"
    );
}

#[test]
fn refreshed_material_and_sgs_match_separate_native_rates_across_repeated_steps() {
    source_and_native_sgs_oracle::<D3q19, 19>(CollisionModel3::Trt {
        tau: 0.9,
        magic: 0.19,
    });
    source_and_native_sgs_oracle::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Trt {
            tau: 0.9,
            magic: 0.17,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn material_sgs_late_refusal_and_cancel_preserve_wall_energy_force_epochs_and_exact_retry() {
    let mut actual = channel(BETA);
    let mut reference = channel(BETA);
    let heat = vec![0.0004; 72];
    for owner in [&mut actual, &mut reference] {
        advance(owner, &heat);
        owner.set_finite_wall_heat(2, &vec![0.0003; 18]).unwrap();
    }
    let physical_reference = material_at(310.0);
    let scaling = plan_scaling_from_material(&physical_reference, DX, DT, 0.0).unwrap();
    let mut material = actual
        .material_transport_workspace(&scaling, material_config())
        .unwrap();
    let mut expected_material = reference
        .material_transport_workspace(&scaling, material_config())
        .unwrap();
    let mut sgs = actual.thermal_sgs_workspace(1024 * 1024).unwrap();
    let mut expected_sgs = reference.thermal_sgs_workspace(1024 * 1024).unwrap();
    let before = accepted(&actual);
    let energy_before = total_energy(&actual);
    let opening_k: Vec<_> = (0..72)
        .map(|cell| 300.0 + 20.0 * actual.thermal().temperature(cell).unwrap())
        .collect();
    let mut background: Vec<_> = (0..72)
        .map(|cell| {
            [
                BACKGROUND[0] + 1e-5 * (1 + cell % 3) as f64,
                BACKGROUND[1] - 2e-5,
                BACKGROUND[2],
            ]
        })
        .collect();
    let saved_background = background.clone();
    let mut expected_background = background.clone();
    let mut polls = 0;
    let expected = reference
        .advance_material_sgs_with_poll(
            Some(&mut expected_background),
            closure(),
            None,
            &heat,
            &mut expected_material,
            &mut expected_sgs,
            [0.0; 3],
            &mut |_, t, _| Ok(material_at(t)),
            &mut || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let changed = expected
        .material
        .endpoint
        .iter()
        .rev()
        .find(|sample| (sample.temperature_k - opening_k[sample.cell]).abs() > 1e-6)
        .expect("conducting wall and heating change actual endpoint temperature");
    let refusal_cell = changed.cell;
    let allowed_change = 0.5 * (changed.temperature_k - opening_k[refusal_cell]).abs();
    let refused = actual.advance_material_sgs_with_poll(
        Some(&mut background),
        closure(),
        None,
        &heat,
        &mut material,
        &mut sgs,
        [0.0; 3],
        &mut |cell, t, _| {
            if cell == refusal_cell && (t - opening_k[cell]).abs() > allowed_change {
                Err(MaterialError::State {
                    what: "outside source temperature support".into(),
                })
            } else {
                Ok(material_at(t))
            }
        },
        &mut || Ok(()),
    );
    assert!(
        matches!(refused, Err(BoussinesqOpenMaterialSgsError3::Material(
        BoussinesqOpenMaterialError3::Material(crate::material::transport::MaterialTransportError::State { cell, .. })
    )) if cell == refusal_cell)
    );
    assert_eq!(accepted(&actual), before);
    assert_eq!(background, saved_background);
    let mut limited = closure();
    limited.maximum_diffusivity_to_molecular = 0.0;
    let refused = actual.advance_material_sgs_with_poll(
        Some(&mut background),
        limited,
        None,
        &heat,
        &mut material,
        &mut sgs,
        [0.0; 3],
        &mut |_, t, _| Ok(material_at(t)),
        &mut || Ok(()),
    );
    assert!(matches!(
        refused,
        Err(BoussinesqOpenMaterialSgsError3::Sgs(
            OpenThermalSgsError3::DiffusivityRatio { .. }
        ))
    ));
    assert_eq!(accepted(&actual), before);
    assert_eq!(background, saved_background);
    let mut calls = 0;
    let cancelled = actual.advance_material_sgs_with_poll(
        Some(&mut background),
        closure(),
        None,
        &heat,
        &mut material,
        &mut sgs,
        [0.0; 3],
        &mut |_, t, _| Ok(material_at(t)),
        &mut || {
            calls += 1;
            if calls == polls {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(BoussinesqOpenMaterialSgsError3::Material(
            BoussinesqOpenMaterialError3::Transport(BoussinesqOpenError3::Flow(
                GridError3::Cancelled
            ))
        ))
    ));
    assert_eq!(calls, polls);
    assert_eq!(accepted(&actual), before);
    assert_eq!(background, saved_background);
    let retry = actual
        .advance_material_sgs_with_poll(
            Some(&mut background),
            closure(),
            None,
            &heat,
            &mut material,
            &mut sgs,
            [0.0; 3],
            &mut |_, t, _| Ok(material_at(t)),
            &mut || Ok(()),
        )
        .unwrap();
    assert_eq!(retry.material.transport, expected.material.transport);
    assert_eq!(retry.material.applied, expected.material.applied);
    assert_eq!(retry.material.endpoint, expected.material.endpoint);
    assert_eq!(retry.sgs, expected.sgs);
    assert_eq!(
        retry.maximum_eddy_diffusivity_lu,
        expected.maximum_eddy_diffusivity_lu
    );
    assert_eq!(
        retry.maximum_diffusivity_to_molecular,
        expected.maximum_diffusivity_to_molecular
    );
    assert_eq!(accepted(&actual), accepted(&reference));
    assert_eq!(background, expected_background);
    let receipt = retry.material.transport.heat;
    let wall_source: f64 = receipt
        .walls
        .into_iter()
        .flatten()
        .filter_map(|wall| wall.solid)
        .map(|solid| solid.source)
        .sum();
    near(
        total_energy(&actual) - energy_before,
        heat.iter().sum::<f64>()
            + wall_source
            + receipt
                .reservoirs
                .iter()
                .map(|reservoir| reservoir.heat_in)
                .sum::<f64>(),
        2e-11,
    );
}
