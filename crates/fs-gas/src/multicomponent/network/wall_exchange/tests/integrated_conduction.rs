//! Manufactured source-integrated material profiles in the real gas/solid owner.

use super::conduction::source;
use super::*;
use crate::multicomponent::network::wall_exchange::{
    MaterialWallConduction, SharedWall, SharedWallControl, SharedWallGas, SharedWallLink,
    SharedWallPort, SharedWallStorage, WallConductionEvaluation,
};
use crate::multicomponent::wall_exchange::MaterialWallCaloric;
use fs_material::thermal::conduction::ReservoirLayerGeometry;
use fs_thermal::enthalpy::DepositionError;

fn network_at(g: &MixtureGas, temperature: f64) -> MixtureNetwork1d {
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        0.4,
        vec![g.conserved(&[1.0], 0.0, temperature).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let plenum = Plenum::new(0.1, g.conserved(&[2.0], 0.0, 2.0 * temperature).unwrap()).unwrap();
    MixtureNetwork1d::new(vec![branch], vec![plenum], vec![], || false).unwrap()
}

fn wall(
    n: &MixtureNetwork1d,
    caloric: MaterialWallCaloric<'static>,
    temperature: f64,
    film: f64,
    layer: MaterialWallConduction,
) -> Result<NetworkWallOperator<'static>, WallExchangeError> {
    uniform_operator(
        n,
        caloric.contact(temperature, film, WallDrag::None).unwrap(),
    )
    .with_material_heat_storage(
        vec![vec![Some(caloric.clone()); 2]],
        vec![Some(caloric)],
        || false,
    )?
    .with_material_conduction(
        vec![vec![Some(layer.clone()); 2]],
        vec![Some(layer)],
        || false,
    )
}

#[test]
fn integrated_profile_matches_nonlinear_storage_and_material_face_calorimetry() {
    let g = gas(1, false, 0.0);
    let geometries = [
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
        ReservoirLayerGeometry::cylindrical_face(0.1, 0.2, 0.1).unwrap(),
    ];
    for geometry in geometries {
        for (tw, tg, heat, span) in [
            (5.0_f64, 3.0_f64, 0.5_f64, [3.5, 10.0]),
            (3.0, 5.0, -0.5, [1.0, 10.0]),
        ] {
            // Original k(T)=T and Cwall(T)=1+0.2T, one explicit kg. With an
            // independently manufactured material face Ts=4, Fourier gives
            // P=C*(Tw²-Ts²)/2 and the film gives P=Gfilm*(Ts-Tg).
            let face = 4.0_f64;
            let power = geometry.conductance_multiplier_m() * 0.5 * (tw * tw - face * face);
            let film = power / (face - tg);
            let dt = heat / power;
            let initial_energy = tw + 0.1 * tw * tw + heat;
            let tw0 = 2.0 * initial_energy / (1.0 + (1.0 + 0.4 * initial_energy).sqrt());
            let tg0 = tg - heat / 0.5;
            let (transport, caloric) = source([0.0, 1.0], [1.0, 0.2], span);
            let material = MaterialWallConduction::from_geometry(transport.clone(), geometry)
                .unwrap()
                .with_integrated_profile()
                .unwrap();
            assert!(material.has_integrated_profile());
            assert_eq!(
                material.evaluation(),
                WallConductionEvaluation::ImplicitEndpoint
            );
            let mut n = network_at(&g, tg0);
            let mut operator = wall(&n, caloric, tw0, film, material).unwrap();
            let original_gas = n.totals().absolute_energy(&g).unwrap();
            let report = operator
                .advance_for(&mut n, Time::new(dt), controls(), || false)
                .unwrap();
            for contact in operator
                .branch_contacts(0)
                .unwrap()
                .iter()
                .chain(operator.plenum_contacts())
            {
                close(wall_temperature(*contact), tw, 3e-10);
                let WallThermal::Lumped {
                    conductance_w_per_k,
                    ..
                } = contact.thermal
                else {
                    panic!("original material temperature storage must remain in the contact");
                };
                close(conductance_w_per_k, power / (tw - tg), 3e-10);
            }
            close(
                g.primitive(n.branches()[0].cell(0).unwrap(), || false)
                    .unwrap()
                    .temperature,
                tg,
                3e-10,
            );
            close(
                g.primitive(n.plenums()[0].state(), || false)
                    .unwrap()
                    .temperature,
                tg,
                3e-10,
            );
            close(report.heat_to_gas_j, 3.0 * heat, 3e-10);
            close(report.wall_energy_change_j, -3.0 * heat, 3e-10);
            close(
                n.totals().absolute_energy(&g).unwrap() - original_gas,
                3.0 * heat,
                3e-10,
            );
            close(report.energy_residual_j, 0.0, 3e-10);
            close(report.maximum_thermal_residual_j, 0.0, 3e-10);
            assert!(
                report.maximum_iterations <= 12,
                "the exact material power tangent must reach the scalar Newton solve"
            );
            assert_eq!(
                operator.branch_conduction(0, 0).unwrap().source(),
                Some(transport.as_ref())
            );
            if heat > 0.0 {
                // The gas is outside the solid's chart; its ACTUAL face is
                // supported. Extending material k all the way to Tg is wrong.
                assert!(transport.at_temperature(tg).is_err());
                assert!(transport.at_temperature(face).is_ok());
            }
            // A wall-endpoint sample predicts a different physical power.
            let pointwise = 1.0 / (1.0 / film + 1.0 / (geometry.conductance_multiplier_m() * tw));
            assert!((dt * pointwise * (tw - tg) - heat).abs() > 1e-3);
        }
    }
}

#[test]
fn integrated_profile_refuses_an_unsupported_material_face_without_state_publication() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let (transport, caloric) = source([0.0, 1.0], [1.0, 0.0], [4.1, 10.0]);
    let layer = MaterialWallConduction::from_geometry(
        transport,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
    )
    .unwrap()
    .with_integrated_profile()
    .unwrap();
    let mut operator = wall(&n, caloric, 5.5, 0.45, layer).unwrap();
    let before = fields(&n);
    let contacts = operator.branch_contacts(0).unwrap().to_vec();
    let plenums = operator.plenum_contacts().to_vec();
    assert!(matches!(
        operator.advance_for(&mut n, Time::new(10.0 / 9.0), controls(), || false),
        Err(WallExchangeError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(operator.branch_contacts(0).unwrap(), contacts);
    assert_eq!(operator.plenum_contacts(), plenums);
    assert!(operator.last_report().is_none());
}

#[test]
fn integrated_profile_cannot_silently_become_lagged_or_infer_a_spatial_phase_state() {
    let g = gas(1, false, 0.0);
    let n = network(&g, &[0.0], true);
    let (transport, caloric) = source([0.0, 1.0], [1.0, 0.0], [1.0, 10.0]);
    let layer = MaterialWallConduction::from_geometry(
        transport,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
    )
    .unwrap()
    .with_integrated_profile()
    .unwrap()
    .with_evaluation(WallConductionEvaluation::InitialState);
    assert!(matches!(
        wall(&n, caloric, 5.5, 0.45, layer),
        Err(WallExchangeError::Invalid(
            "integrated wall profile requires implicit endpoint evaluation"
        ))
    ));
    let phase = super::phase_conduction::source(1.0, 1.0, 5.0);
    let layer = MaterialWallConduction::from_phase_source_with_geometry(
        phase,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        layer.with_integrated_profile(),
        Err(WallExchangeError::Invalid(
            "integrated conductivity requires a temperature-only source; spatial phase state is absent"
        ))
    ));
}

#[test]
fn integrated_profile_refuses_wall_links_without_a_two_layer_endpoint_law() {
    let g = gas(1, false, 0.0);
    let n = network(&g, &[0.0], true);
    let (transport, caloric) = source([0.0, 1.0], [1.0, 0.0], [1.0, 10.0]);
    let layer = MaterialWallConduction::from_geometry(
        transport,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
    )
    .unwrap()
    .with_integrated_profile()
    .unwrap();
    let shared = |cell| SharedWall {
        storage: SharedWallStorage {
            thermal: caloric.contact(5.0, 0.0, WallDrag::None).unwrap().thermal,
            heat_storage: Some(caloric.clone()),
        },
        ports: vec![SharedWallPort {
            gas: SharedWallGas::Branch { branch: 0, cell },
            conductance_w_per_k: 0.45,
        }],
    };
    let build = || {
        uniform_operator(
            &n,
            WallContact {
                drag: WallDrag::None,
                thermal: WallThermal::Adiabatic,
            },
        )
        .with_shared_walls(
            &n,
            vec![shared(0), shared(1)],
            SharedWallControl::default(),
            || false,
        )
        .unwrap()
    };
    build()
        .with_shared_wall_conduction(vec![vec![Some(layer.clone())], vec![None]], || false)
        .unwrap();
    assert!(matches!(
        build()
            .with_shared_wall_links(
                vec![SharedWallLink {
                    walls: [0, 1],
                    conductance_w_per_k: 0.7
                }],
                || false,
            )
            .unwrap()
            .with_shared_wall_link_conduction(vec![[Some(layer), None]], || false,),
        Err(WallExchangeError::Invalid(
            "integrated wall link profiles require a two-layer endpoint law"
        ))
    ));
}

#[test]
fn integrated_profile_disabled_film_remains_insulating_outside_the_material_chart() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let (transport, caloric) = source([0.0, 1.0], [1.0, 0.0], [4.1, 10.0]);
    let layer = MaterialWallConduction::from_geometry(
        transport,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
    )
    .unwrap()
    .with_integrated_profile()
    .unwrap();
    let mut operator = wall(&n, caloric, 5.5, 0.0, layer).unwrap();
    let before = fields(&n);
    let report = operator
        .advance_for(&mut n, Time::new(1.0), controls(), || false)
        .unwrap();
    assert_eq!(fields(&n), before);
    assert_eq!(report.heat_to_gas_j, 0.0);
    assert_eq!(report.wall_energy_change_j, 0.0);
    close(
        wall_temperature(operator.branch_contacts(0).unwrap()[0]),
        5.5,
        0.0,
    );
}

#[test]
fn integrated_profile_brackets_supported_small_and_large_steps_without_projecting_a_face() {
    let g = gas(1, false, 0.0);
    for dt in [0.1_f64, 100.0, 200.0] {
        let mut n = network(&g, &[0.0], true);
        let (transport, caloric) = source([0.0, 1.0], [1.0, 0.0], [4.5, 10.0]);
        let layer = MaterialWallConduction::from_geometry(
            transport,
            ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
        )
        .unwrap()
        .with_integrated_profile()
        .unwrap();
        let mut operator = wall(&n, caloric, 5.5, 0.01, layer).unwrap();
        let original = fields(&n);
        let original_contacts = operator.branch_contacts(0).unwrap().to_vec();
        let original_plenums = operator.plenum_contacts().to_vec();
        // Exact original Fourier, film and storage equations eliminate Ts:
        // Ts=2+(2+100/dt)*Q, Tw=5.5-Q, Tg=2+2Q, and
        // ((2+100/dt)^2-1)*Q^2 + (19+420/dt)*Q - 26.25 = 0.
        let a = (2.0 + 100.0 / dt).powi(2) - 1.0;
        let b = 19.0 + 420.0 / dt;
        let heat = 52.5 / (b + (b * b + 105.0 * a).sqrt());
        let face = 2.0 + (2.0 + 100.0 / dt) * heat;
        let result = operator.advance_for(&mut n, Time::new(dt), controls(), || false);
        if dt == 200.0 {
            // The real root has a supported wall but an unsupported face. A
            // source-aware bracket must refuse it, never clip Ts to 4.5 K.
            assert!(5.5 - heat > 4.5 && face < 4.5);
            assert!(result.is_err());
            assert_eq!(fields(&n), original);
            assert_eq!(operator.branch_contacts(0).unwrap(), original_contacts);
            assert_eq!(operator.plenum_contacts(), original_plenums);
            assert!(operator.last_report().is_none());
        } else {
            assert!(face > 4.5);
            let report = result.unwrap();
            for &contact in operator
                .branch_contacts(0)
                .unwrap()
                .iter()
                .chain(operator.plenum_contacts())
            {
                close(wall_temperature(contact), 5.5 - heat, 3e-10);
            }
            close(
                g.primitive(n.branches()[0].cell(0).unwrap(), || false)
                    .unwrap()
                    .temperature,
                2.0 + 2.0 * heat,
                3e-10,
            );
            close(report.heat_to_gas_j, 3.0 * heat, 3e-10);
            close(report.wall_energy_change_j, -3.0 * heat, 3e-10);
            if dt == 100.0 {
                // The storage bound Q=1 is unsupported; its first supported
                // half Q=.5 is still BELOW the root and cannot close the bracket.
                assert!(heat > 0.5);
            }
        }
    }
}

#[test]
fn integrated_profile_gas_inverse_separates_physical_heat_from_a_global_slope_bound() {
    use crate::multicomponent::wall_exchange::{
        WallGasCaloric, WallGasState, refine_caloric_temperature_with_allowance,
    };
    use fs_material::thermal::conduction::{ConductivityCurve, LayerConductionConfig};
    // An explicit unit-capacity caloric owner isolates the inverse-error budget;
    // the real original material kernel below evaluates the sharply varying k.
    struct UnitGas(crate::dynamics::CaloricDomain);
    impl WallGasCaloric for UnitGas {
        fn domain(&self) -> crate::dynamics::CaloricDomain {
            self.0
        }
        fn primitive(
            &self,
            momentum: f64,
            energy: f64,
            _: &mut dyn FnMut() -> bool,
        ) -> Result<WallGasState, WallExchangeError> {
            Ok(WallGasState {
                density: 1.0,
                velocity: momentum,
                temperature: energy - 0.5 * momentum * momentum,
            })
        }
        fn caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError> {
            Ok((temperature, 1.0))
        }
    }
    let original_gas = gas(1, false, 0.0);
    let caloric = UnitGas(original_gas.domain());
    let curve = ConductivityCurve::new(vec![
        [1.0, 1.0],
        [1.99, 1.0],
        [2.0, 50.0],
        [2.01, 1.0],
        [100.0, 1.0],
    ])
    .unwrap();
    let geometry = ReservoirLayerGeometry::planar(1.0, 1.0).unwrap();
    let film = 98.245_f64;
    let bound = 1.0 / (1.0 / film + 1.0 / 50.0);
    let control = WallControl {
        absolute_energy_tolerance_j: 1e-14,
        relative_energy_tolerance: 1e-6,
        maximum_iterations: 80,
    };
    let physical_budget = control.tolerance(98.245);
    let coarse = 1.0 + 1e-5;
    let refined = refine_caloric_temperature_with_allowance(
        &caloric,
        coarse,
        1.0,
        1.0,
        bound,
        physical_budget,
        control,
        &mut || false,
    )
    .unwrap();
    let power = |temperature| {
        geometry
            .evaluate_material_to_bath(
                &curve,
                [100.0, temperature],
                1.0 / film,
                LayerConductionConfig {
                    max_iterations: 96,
                    absolute_heat_tolerance_w: 1e-12,
                    relative_heat_tolerance: 1e-12,
                },
                || Ok(()),
            )
            .unwrap()
            .power_to_reservoir_w
    };
    let exact = power(1.0);
    close(exact, 98.245, 1e-9);
    assert!((power(coarse) - exact).abs() > physical_budget);
    assert!((power(refined) - exact).abs() <= 0.125 * physical_budget);
}
