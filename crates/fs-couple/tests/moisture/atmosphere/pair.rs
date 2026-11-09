//! G1/G3/G4: two finite chambers, original diffusion, and atomic refusal.
//! All material and gas data are authored verification fixtures.
#[path = "pair/adaptive.rs"]
mod adaptive;

use super::*;
use fs_blake3::ContentHash;
use fs_couple::moisture::{
    atmosphere::pair::*,
    cartesian::{MoistureCartesian, MoistureCartesianGeometry},
};

fn pair_controls() -> PairedAtmosphereMoistureSolve {
    PairedAtmosphereMoistureSolve {
        atmospheres: [controls(); 2],
        maximum_evaluations: 10_000,
        mass_tolerance_kg: 5e-12,
    }
}
fn port(atmosphere: MoistureAtmospherePort, slope: f64) -> Option<AtmosphereMoistureFace> {
    Some(AtmosphereMoistureFace {
        atmosphere,
        sorption: sorption(slope),
    })
}
fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

/// Independently eliminate the two vapor inventories from the one-cell
/// backward-Euler equations. alpha[i] sums 2*dt*D/length^2 over that port's faces.
/// q_i = dry*alpha_i*(s_i*r_i - M), C_i*(r_i-r_i0) = -q_i.
fn linear_endpoint(
    dry: f64,
    old: f64,
    rh: [f64; 2],
    volume: [f64; 2],
    slope: [f64; 2],
    alpha: [f64; 2],
) -> (f64, [f64; 2], [f64; 2]) {
    let capacity = volume.map(|v| saturation_pressure() * v * WATER_MOLAR_MASS / (R * TEMPERATURE));
    let weight: [f64; 2] = std::array::from_fn(|i| {
        dry * alpha[i] / (1.0 + dry * alpha[i] * slope[i] / capacity[i])
    });
    let moisture = (dry * old + weight[0] * slope[0] * rh[0] + weight[1] * slope[1] * rh[1])
        / (dry + weight[0] + weight[1]);
    let transfer: [f64; 2] =
        std::array::from_fn(|i| weight[i] * (slope[i] * rh[i] - moisture));
    let endpoint = std::array::from_fn(|i| rh[i] - transfer[i] / capacity[i]);
    (moisture, endpoint, transfer)
}

#[test]
fn g1_two_independent_inventories_match_the_one_cell_balance_and_continue() {
    with_cx(|cx| {
        let diffusion = 0.01;
        let dt = 2.0;
        let volume = [0.7, 1.3];
        let slope = [0.22, 0.28];
        let faces = [
            port(MoistureAtmospherePort::First, slope[0]),
            port(MoistureAtmospherePort::Second, slope[1]),
        ];
        let original = slab(card([diffusion; 2], 0.3), 1.0, vec![0.1], cx);
        for mut moisture in [
            original.clone(),
            original.clone().with_implicit_diffusivity(nonlinear()).unwrap(),
            original.clone().with_kirchhoff_diffusivity(nonlinear()).unwrap(),
        ] {
            let mut gas = [atmosphere(0.85, volume[0]), atmosphere(0.15, volume[1])];
            let carrier = gas.each_ref().map(|g| carrier_mass(g));
            let identity = moisture.definition_identity();
            for interval in 0..2 {
                let rh = gas.each_ref().map(|g| g.relative_humidity().unwrap());
                let expected = linear_endpoint(
                    moisture.dry_mass_kg(), moisture.moisture()[0],
                    rh, volume, slope, [2.0 * dt * diffusion; 2],
                );
                let before_mass = gas[0].gas().mass_kg() + gas[1].gas().mass_kg()
                    + moisture.wet_mass_kg();
                let candidate = gas[0].advance_paired(
                    &gas[1], &moisture, dt, &faces, config(), pair_controls(), cx,
                ).unwrap();
                near(candidate.moisture().moisture()[0], expected.0, 3e-12);
                let step = candidate.step();
                assert!(step.field_evaluations > 1);
                assert!(step.field_evaluations <= pair_controls().maximum_evaluations);
                near(step.water_balance_residual_kg, 0.0, pair_controls().mass_tolerance_kg);
                near(step.total_mass_residual_kg, 0.0, pair_controls().mass_tolerance_kg);
                near(step.water_change_kg,
                    moisture.dry_mass_kg() * (expected.0 - moisture.moisture()[0]), 2e-12);
                for i in 0..2 {
                    let next = &candidate.atmospheres()[i];
                    let account = &step.atmospheres[i];
                    near(next.relative_humidity().unwrap(), expected.1[i], 3e-11);
                    near(account.water_to_solids_kg, expected.2[i], 3e-13);
                    near(step.transport.boundary_water_kg[i], account.water_to_solids_kg, 0.0);
                    near(carrier_mass(next), carrier[i], 3e-14);
                    near(account.water_balance_residual_kg, 0.0, controls().mass_tolerance_kg);
                    near(next.gas().state().state().pressure,
                        gas[i].gas().state().state().pressure
                            + saturation_pressure() * (expected.1[i] - rh[i]), 2e-7);
                    assert_eq!(next.source_catalog(), gas[i].source_catalog());
                    assert_eq!(account.interval_s, [interval as f64 * dt, (interval + 1) as f64 * dt]);
                    let source = step.transport.surface_sorption[i].as_ref().unwrap();
                    near(source.property(SORPTION_PROPERTY).unwrap().value_si(),
                        slope[i] * expected.1[i], 1e-11);
                }
                near(candidate.atmospheres()[0].gas().mass_kg()
                    + candidate.atmospheres()[1].gas().mass_kg()
                    + candidate.moisture().wet_mass_kg(), before_mass, 5e-12);
                assert_eq!(candidate.moisture().definition_identity(), identity);
                let (next_gas, next_field, _) = candidate.into_parts();
                gas = next_gas;
                moisture = next_field;
            }
        }
    });
}

#[test]
fn g1_chamber_to_chamber_transfer_survives_zero_net_solid_uptake_and_port_reversal() {
    with_cx(|cx| {
        let gas = [atmosphere(0.8, 1.0), atmosphere(0.2, 1.0)];
        let moisture = slab(card([0.01; 2], 0.3), 1.0, vec![0.125], cx);
        let expected = linear_endpoint(
            moisture.dry_mass_kg(), 0.125, [0.8, 0.2], [1.0; 2], [0.25; 2], [0.04; 2],
        );
        for reversed in [false, true] {
            let assignments = if reversed {
                [MoistureAtmospherePort::Second, MoistureAtmospherePort::First]
            } else {
                [MoistureAtmospherePort::First, MoistureAtmospherePort::Second]
            };
            let faces = assignments.map(|p| port(p, 0.25));
            let candidate = gas[0].advance_paired(
                &gas[1], &moisture, 2.0, &faces, config(), pair_controls(), cx,
            ).unwrap();
            near(candidate.moisture().moisture()[0], 0.125, 2e-12);
            near(candidate.step().water_change_kg, 0.0, 1e-12);
            assert!(candidate.step().atmospheres[0].water_to_solids_kg > 1e-4);
            assert!(candidate.step().atmospheres[1].water_to_solids_kg < -1e-4);
            for i in 0..2 {
                near(candidate.atmospheres()[i].relative_humidity().unwrap(), expected.1[i], 3e-11);
                near(candidate.step().atmospheres[i].water_to_solids_kg, expected.2[i], 3e-13);
                let face = if reversed { 1 - i } else { i };
                near(candidate.step().transport.boundary_water_kg[face], expected.2[i], 3e-13);
            }
            assert_eq!(candidate.step().face_atmospheres, assignments.map(Some));
        }
        // Swap the ordered gas owners as well as the face labels: the physical
        // result must be independent of which chamber is the nested inner solve.
        let reversed = gas[1].advance_paired(
            &gas[0], &moisture, 2.0,
            &[port(MoistureAtmospherePort::Second, 0.25),
              port(MoistureAtmospherePort::First, 0.25)],
            config(), pair_controls(), cx,
        ).unwrap();
        near(reversed.atmospheres()[1].relative_humidity().unwrap(), expected.1[0], 3e-11);
        near(reversed.atmospheres()[0].relative_humidity().unwrap(), expected.1[1], 3e-11);
    });
}

#[test]
fn g1_cartesian_six_faces_debit_only_their_assigned_chamber() {
    with_cx(|cx| {
        let length = 0.1_f64;
        let dt = 2.0;
        let diffusion = 0.01;
        let field = MoistureCartesian::new(
            card([diffusion; 2], 0.3), point(0.1), ["moisture_diffusivity"; 3],
            MaterialPropertySelection::SingleClaimOnly,
            MoistureCartesianGeometry {
                lengths_m: [length; 3], cells: [1; 3], frame: ContentHash([7; 32]),
                maximum_moisture: 0.3, maximum_cells: 1,
            },
            vec![0.1], cx,
        ).unwrap();
        let gas = [atmosphere(0.8, 0.8), atmosphere(0.2, 1.2)];
        let assignments = [
            MoistureAtmospherePort::Second, MoistureAtmospherePort::First,
            MoistureAtmospherePort::First, MoistureAtmospherePort::Second,
            MoistureAtmospherePort::Second, MoistureAtmospherePort::First,
        ];
        let faces = assignments.map(|p| port(p, 0.25));
        let expected = linear_endpoint(
            field.dry_mass_kg(), 0.1, [0.8, 0.2], [0.8, 1.2], [0.25; 2],
            [3.0 * 2.0 * dt * diffusion / length.powi(2); 2],
        );
        let candidate = gas[0].advance_paired_cartesian(
            &gas[1], &field, dt, &faces, config(), nonlinear(), pair_controls(), cx,
        ).unwrap();
        near(candidate.moisture().moisture()[0], expected.0, 3e-12);
        for i in 0..2 {
            near(candidate.atmospheres()[i].relative_humidity().unwrap(), expected.1[i], 3e-11);
            near(candidate.step().atmospheres[i].water_to_solids_kg, expected.2[i], 4e-12);
        }
        for (side, assignment) in assignments.into_iter().enumerate() {
            let i = usize::from(assignment == MoistureAtmospherePort::Second);
            near(candidate.step().transport.boundary_water_kg[side], expected.2[i] / 3.0, 2e-12);
        }
        near(candidate.step().water_balance_residual_kg, 0.0, pair_controls().mass_tolerance_kg);
        near(candidate.step().total_mass_residual_kg, 0.0, pair_controls().mass_tolerance_kg);
        assert_eq!(candidate.moisture().geometry().cells, [1; 3]);
    });
}

#[test]
fn g3_single_chamber_and_sealed_limits_reuse_the_original_transport() {
    with_cx(|cx| {
        let first = atmosphere(0.8, 1.0);
        let second = atmosphere(0.2, 0.7);
        let field = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 4], cx);
        let single = first.advance(
            &field, 2.0, &[Some(sorption(0.25)), None], config(), controls(), cx,
        ).unwrap();
        let paired = first.advance_paired(
            &second, &field, 2.0, &[port(MoistureAtmospherePort::First, 0.25), None],
            config(), pair_controls(), cx,
        ).unwrap();
        for (&a, &b) in paired.moisture().moisture().iter().zip(single.moisture().moisture()) {
            near(a, b, 3e-12);
        }
        near(paired.atmospheres()[0].relative_humidity().unwrap(),
            single.atmosphere().relative_humidity().unwrap(), 3e-11);
        assert_eq!(&paired.atmospheres()[1], &second);
        assert_eq!(paired.step().atmospheres[1].water_to_solids_kg, 0.0);
        let nonuniform = slab(card([0.01; 2], 0.3), 1.0, vec![0.05, 0.15, 0.2, 0.1], cx);
        let native = nonuniform.advance(
            2.0, &[MoistureFace::Sealed, MoistureFace::Sealed], config(), cx,
        ).unwrap();
        let sealed = first.advance_paired(
            &second, &nonuniform, 2.0, &[None, None], config(), pair_controls(), cx,
        ).unwrap();
        assert_eq!(sealed.moisture().moisture(), native.0.moisture());
        assert_eq!(sealed.atmospheres(), &[first, second]);
        assert_eq!(sealed.step().field_evaluations, 1);
        near(sealed.step().water_change_kg, 0.0, 2e-14);
    });
}

#[test]
fn g4_work_source_and_global_balance_refusals_leave_all_three_inputs_unchanged() {
    with_cx(|cx| {
        let gas = [atmosphere(0.8, 1.0), atmosphere(0.2, 1.0)];
        let saved = gas.clone();
        let field = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let faces = [
            port(MoistureAtmospherePort::First, 0.25),
            port(MoistureAtmospherePort::Second, 0.25),
        ];
        let mut bounded = pair_controls();
        bounded.maximum_evaluations = 1;
        assert!(matches!(
            gas[0].advance_paired(&gas[1], &field, 2.0, &faces, config(), bounded, cx),
            Err(AtmosphereMoistureError::EvaluationBudget { evaluations: 1, maximum: 1 })
        ));
        // Individually tolerated reservoirs cannot spend the whole-system
        // budget independently. Initial RH is deliberately accepted by each
        // loose gas gate despite a nonzero combined solid uptake.
        let mut loose = pair_controls();
        loose.atmospheres[0].mass_tolerance_kg = 0.1;
        loose.atmospheres[1].mass_tolerance_kg = 0.1;
        assert!(matches!(
            gas[0].advance_paired(&gas[1], &field, 2.0, &faces, config(), loose, cx),
            Err(AtmosphereMoistureError::MassBalance { .. })
        ));
        let unsupported = DryBasisSorption::new(
            sorption_card(0.25, MoistureBasis::DryMass, 0.6), sorption_point(0.5),
            SORPTION_PROPERTY, MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        let bad_faces = [
            Some(AtmosphereMoistureFace {
                atmosphere: MoistureAtmospherePort::First, sorption: unsupported,
            }),
            faces[1].clone(),
        ];
        assert!(gas[0].advance_paired(
            &gas[1], &field, 2.0, &bad_faces, config(), pair_controls(), cx,
        ).is_err());
        assert_eq!(gas, saved);
        assert_eq!(field.moisture(), &[0.1]);
        assert_eq!(field.time_s(), 0.0);
        // An exhausted nested budget leaves the accepted input usable.
        assert!(gas[0].advance_paired(
            &gas[1], &field, 2.0, &faces, config(), pair_controls(), cx,
        ).is_ok());
    });
}

#[test]
fn g4_cancelled_pair_publishes_no_gas_or_field_endpoint() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 3, kernel_id: 17, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        );
        let gas = [atmosphere(0.8, 1.0), atmosphere(0.2, 1.0)];
        let saved = gas.clone();
        let field = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], &cx);
        gate.request();
        assert!(matches!(
            gas[0].advance_paired(
                &gas[1], &field, 2.0,
                &[port(MoistureAtmospherePort::First, 0.25),
                  port(MoistureAtmospherePort::Second, 0.25)],
                config(), pair_controls(), &cx,
            ),
            Err(AtmosphereMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(gas, saved);
        assert_eq!(field.moisture(), &[0.1]);
        assert_eq!(field.time_s(), 0.0);
    });
}
