//! Original material density, molecular transport and native source admission.
use super::*;
use crate::lattice3::open_boundaries::sgs::OpenThermalSgsConfig3;
use crate::sgs3::{SgsConfig, SgsModel, SgsUse};
use crate::thermal::open::buoyancy::BoussinesqDensityLawError3;

const PHYSICAL_GRAVITY: [f64; 3] = [0.0, -0.02, 0.0];

fn nonlinear_material(temperature_k: f64) -> FluidTransportState {
    let original = material_at(temperature_k);
    let difference = temperature_k - 310.0;
    FluidTransportState::try_new(
        PHYSICAL_RHO * (1.0 - 0.00002 * difference * difference),
        original.specific_heat_j_kg_k(),
        original.thermal_conductivity_w_m_k(),
        original.dynamic_viscosity_pa_s(),
    )
    .unwrap()
}

fn resolve_material(
    temperature_k: f64,
    _: &FluidTransportState,
) -> Result<FluidTransportState, MaterialError> {
    Ok(nonlinear_material(temperature_k))
}

#[test]
fn source_density_uses_original_reference_kelvin_and_si_force_units() {
    let reference = material_at(310.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.0).unwrap();
    for config in [
        material_config(),
        BoussinesqOpenMaterialConfig3 {
            temperature_at_zero_k: 320.0,
            temperature_scale_k: -20.0,
            ..material_config()
        },
    ] {
        let law =
            material_density_law3(&scaling, config, PHYSICAL_GRAVITY, resolve_material).unwrap();
        assert!(core::ptr::eq(law.source().reference_source(), &reference));
        assert_eq!(law.source().reference_transport(), &reference);
        assert_eq!(law.source().config(), config);
        assert_eq!(law.reference_density(), 1.0);
        for theta in [0.0, 0.5, 1.0] {
            let kelvin = config.temperature_at_zero_k + config.temperature_scale_k * theta;
            let difference = kelvin - 310.0;
            let departure = -0.00002 * difference * difference;
            let actual = law.sample(theta, BACKGROUND).unwrap();
            near(actual.relative_density_departure, departure, 2e-16);
            for axis in 0..3 {
                let expected = departure * PHYSICAL_GRAVITY[axis] * DT * DT / DX;
                near(actual.buoyancy[axis], expected, 2e-17);
                near(actual.total[axis], BACKGROUND[axis] + expected, 2e-17);
            }
            if theta == 0.5 {
                assert_eq!(actual.total.map(f64::to_bits), BACKGROUND.map(f64::to_bits));
            } else {
                assert!(actual.buoyancy[1] > 0.0);
            }
        }
    }
    let law = material_density_law3(
        &scaling,
        material_config(),
        PHYSICAL_GRAVITY,
        |kelvin, _| {
            if kelvin > 315.0 {
                return Err(MaterialError::State {
                    what: "original density support".into(),
                });
            }
            resolve_material(kelvin, &reference)
        },
    )
    .unwrap();
    assert!(matches!(
        law.sample(1.0, BACKGROUND),
        Err(BoussinesqDensityLawError3::Source(MaterialError::State { what }))
            if what == "original density support"
    ));
    assert!(matches!(
        law.sample(-20.0, BACKGROUND),
        Err(BoussinesqDensityLawError3::Source(_))
    ));
    let bounded_cp =
        material_density_law3(&scaling, material_config(), PHYSICAL_GRAVITY, |_, _| {
            FluidTransportState::try_new(PHYSICAL_RHO, 2.0 * PHYSICAL_CP, 2.0, 0.002)
        })
        .unwrap();
    assert!(matches!(
        bounded_cp.sample(0.5, BACKGROUND),
        Err(BoussinesqDensityLawError3::Source(_))
    ));
}

#[test]
fn nonlinear_density_composes_with_original_molecular_and_sgs_transport() {
    let reference = material_at(310.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.0).unwrap();
    let config = material_config();
    let collision = CollisionModel3::Bgk { tau: 0.9 };
    let make = || {
        make_with_policy::<D3q19, 19, _>(
            collision,
            [6, 4, 3],
            native_channel(),
            channel_faces(),
            0.5,
            material_density_law3(&scaling, config, PHYSICAL_GRAVITY, resolve_material).unwrap(),
        )
    };
    let mut actual = make();
    let mut explicit = make();
    let cells = actual.thermal().populations().len();
    let heat: Vec<_> = (0..cells)
        .map(|cell| 0.001 * (1 + cell % 4) as f64)
        .collect();
    advance(&mut actual, &heat);
    advance(&mut explicit, &heat);
    let mut material_workspace = actual
        .material_transport_workspace(&scaling, config)
        .unwrap();
    let mut sgs_workspace = actual.thermal_sgs_workspace(1024 * 1024).unwrap();
    let mut native_sgs = explicit.flow().sgs_workspace(1024 * 1024).unwrap();
    let closure = OpenThermalSgsConfig3 {
        closure: SgsConfig {
            model: SgsModel::Smagorinsky { cs: 0.24 },
            filter_width_lu: 1.3,
            maximum_eddy_to_molecular: 20.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.73,
        maximum_diffusivity_to_molecular: 20.0,
    };
    let reads = Arc::new(AtomicUsize::new(0));
    let mut observed_buoyancy = 0.0_f64;
    let mut observed_eddy = 0.0_f64;
    for epoch in 0..3 {
        let mut collisions = Vec::new();
        let mut molecular_tau = Vec::new();
        for cell in 0..cells {
            let kelvin = 300.0 + 20.0 * explicit.thermal().temperature(cell).unwrap();
            let source = nonlinear_material(kelvin);
            let units = (DT / DX) / DX;
            // Independent SI conversion preserves original reference rho/Cp.
            let tau = 0.5 + source.dynamic_viscosity_pa_s() / PHYSICAL_RHO * units / (1.0 / 3.0);
            collisions.push(D3q19::with_shear_relaxation(collision, tau).unwrap());
            molecular_tau.push(
                0.5 + source.thermal_conductivity_w_m_k() / PHYSICAL_RHO / PHYSICAL_CP * units
                    / 0.25,
            );
        }
        let fit = explicit
            .flow
            .inner
            .prepare_open_sgs(
                closure.closure,
                Some(&collisions),
                &mut native_sgs,
                Some(&explicit.flow.forces),
                || Ok(()),
            )
            .unwrap();
        let effective_tau: Vec<_> = molecular_tau
            .iter()
            .zip(native_sgs.eddy_viscosities())
            .map(|(&tau, &eddy)| {
                if eddy == 0.0 {
                    tau
                } else {
                    tau + 4.0 * (eddy / closure.turbulent_prandtl)
                }
            })
            .collect();
        let mut background: Vec<_> = (0..cells)
            .map(|cell| {
                [
                    BACKGROUND[0] + (epoch + 1) as f64 * 1e-5 * (1 + cell % 3) as f64,
                    BACKGROUND[1],
                    BACKGROUND[2],
                ]
            })
            .collect();
        let desired_background = background.clone();
        let mut expected_background = background.clone();
        let expected = explicit
            .advance_with_diffusion_poll(
                Some(&mut expected_background),
                Some(native_sgs.collisions()),
                &heat,
                Some(&effective_tau),
                None,
                [0.0; 3],
                &mut || Ok(()),
            )
            .unwrap();
        reads.store(0, Ordering::SeqCst);
        let energy_before = total_energy(&actual);
        let mut ordinal = 0;
        let report = actual
            .advance_material_sgs_with_poll(
                Some(&mut background),
                closure,
                None,
                &heat,
                &mut material_workspace,
                &mut sgs_workspace,
                [0.0; 3],
                &mut |cell, temperature_k, retained| {
                    assert_eq!(*retained, reference);
                    let original = TracedMaterial {
                        cell,
                        temperature_k,
                        ordinal,
                        state: nonlinear_material(temperature_k),
                        reads: Arc::clone(&reads),
                    };
                    ordinal += 1;
                    Ok(original)
                },
                &mut || Ok(()),
            )
            .unwrap();
        assert_eq!(report.material.transport, expected);
        assert_eq!(report.sgs, fit);
        assert_eq!(accepted(&actual), accepted(&explicit));
        assert_eq!(background, expected_background);
        assert_eq!(material_workspace.thermal_relaxation(), molecular_tau);
        assert_eq!(sgs_workspace.thermal_relaxation(), effective_tau);
        assert_eq!(reads.load(Ordering::SeqCst), 2 * cells);
        for samples in [&report.material.applied, &report.material.endpoint] {
            for sample in samples {
                assert_eq!(sample.material.material().cell, sample.cell);
                assert_eq!(
                    sample.material.material().temperature_k,
                    sample.temperature_k
                );
                assert_eq!(sample.transport.density_kg_m3(), PHYSICAL_RHO);
                assert_eq!(sample.transport.specific_heat_j_kg_k(), PHYSICAL_CP);
            }
        }
        let thermal = report.material.transport.heat;
        near(
            total_energy(&actual) - energy_before,
            heat.iter().sum::<f64>() + thermal.reservoirs.iter().map(|r| r.heat_in).sum::<f64>(),
            2e-10,
        );
        near(thermal.combined_defect, 0.0, 2e-10);
        observed_eddy = observed_eddy.max(report.sgs.maximum_eddy_viscosity_lu);
        for cell in 0..cells {
            let kelvin = 300.0 + 20.0 * actual.thermal().temperature(cell).unwrap();
            let rho = PHYSICAL_RHO * (1.0 - 0.00002 * (kelvin - 310.0).powi(2));
            for axis in 0..3 {
                let buoyancy =
                    (rho - PHYSICAL_RHO) / PHYSICAL_RHO * PHYSICAL_GRAVITY[axis] * DT * DT / DX;
                near(
                    actual.flow().forces()[cell][axis],
                    desired_background[cell][axis] + buoyancy,
                    2e-12,
                );
                observed_buoyancy = observed_buoyancy.max(buoyancy.abs());
            }
        }
    }
    assert!(observed_buoyancy > 1e-7);
    assert!(observed_eddy > 0.0);
}

#[test]
fn source_density_refusal_and_final_cancellation_preserve_joint_history_then_retry() {
    let reference = material_at(310.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.0).unwrap();
    let resolve = |kelvin: f64, _: &FluidTransportState| {
        if kelvin > 305.3 {
            Err(MaterialError::State {
                what: "finite original density domain".into(),
            })
        } else {
            Ok(nonlinear_material(kelvin))
        }
    };
    let make = || {
        make_with_policy::<D3q19, 19, _>(
            CollisionModel3::Bgk { tau: 0.9 },
            [4, 3, 2],
            BoundarySpec3::periodic(),
            [ThermalOpenFace3::Periodic; 6],
            0.25,
            material_density_law3(&scaling, material_config(), PHYSICAL_GRAVITY, resolve).unwrap(),
        )
    };
    let mut actual = make();
    let mut clean = make();
    let cells = actual.thermal().populations().len();
    let mut workspace = actual
        .material_transport_workspace(&scaling, material_config())
        .unwrap();
    let mut clean_workspace = clean
        .material_transport_workspace(&scaling, material_config())
        .unwrap();
    let saved = accepted(&actual);
    let mut background = vec![[4e-6, 1e-6, -2e-6]; cells];
    let proposed = background.clone();
    let source = |_: usize, kelvin: f64, _: &FluidTransportState| Ok(nonlinear_material(kelvin));
    let result = actual.advance_material_with_poll(
        Some(&mut background),
        None,
        &vec![0.02; cells],
        &mut workspace,
        [0.0; 3],
        &mut source.clone(),
        &mut || Ok(()),
    );
    assert!(matches!(
        result,
        Err(BoussinesqOpenMaterialError3::Transport(BoussinesqOpenError3::Law {
            source: BoussinesqDensityLawError3::Source(MaterialError::State { what }), ..
        })) if what == "finite original density domain"
    ));
    assert_eq!(accepted(&actual), saved);
    assert_eq!(background, proposed);
    let heat = vec![0.002; cells];
    let mut control_background = proposed.clone();
    let expected = clean
        .advance_material_with_poll(
            Some(&mut control_background),
            None,
            &heat,
            &mut clean_workspace,
            [0.0; 3],
            &mut source.clone(),
            &mut || Ok(()),
        )
        .unwrap();
    let retried = actual
        .advance_material_with_poll(
            Some(&mut background),
            None,
            &heat,
            &mut workspace,
            [0.0; 3],
            &mut source.clone(),
            &mut || Ok(()),
        )
        .unwrap();
    assert_eq!(retried.transport, expected.transport);
    assert_eq!(accepted(&actual), accepted(&clean));
    let saved = accepted(&actual);
    let mut proposed = vec![[5e-6, 2e-6, -3e-6]; cells];
    let unchanged = proposed.clone();
    let mut control_proposed = proposed.clone();
    let mut polls = 0;
    let expected = clean
        .advance_material_with_poll(
            Some(&mut control_proposed),
            None,
            &heat,
            &mut clean_workspace,
            [0.0; 3],
            &mut source.clone(),
            &mut || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut seen = 0;
    let cancelled = actual.advance_material_with_poll(
        Some(&mut proposed),
        None,
        &heat,
        &mut workspace,
        [0.0; 3],
        &mut source.clone(),
        &mut || {
            seen += 1;
            if seen == polls {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(BoussinesqOpenMaterialError3::Transport(
            BoussinesqOpenError3::Flow(GridError3::Cancelled)
        ))
    ));
    assert_eq!(accepted(&actual), saved);
    assert_eq!(proposed, unchanged);
    let retried = actual
        .advance_material_with_poll(
            Some(&mut proposed),
            None,
            &heat,
            &mut workspace,
            [0.0; 3],
            &mut source.clone(),
            &mut || Ok(()),
        )
        .unwrap();
    assert_eq!(retried.transport, expected.transport);
    assert_eq!(accepted(&actual), accepted(&clean));
    assert_eq!(proposed, control_proposed);
    near(retried.transport.heat.combined_defect, 0.0, 1e-10);
    for cell in 0..cells {
        near(actual.flow().macroscopic(cell).unwrap().0, RHO, 1e-11);
    }
}
