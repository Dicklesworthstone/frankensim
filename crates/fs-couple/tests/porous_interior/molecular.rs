//! G1/G3: original molecular reservoirs throughout complete pressure volumes.
//! Authored reservoir data exercise equations, not empirical calibration.
use super::*;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::MechanicalBoundary;
use fs_material::porous::{
    PorousFluidResponse,
    layers::{Layer, RigidBackedLayers},
};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_phs::{WallPin, acoustic_outlet::AcousticWallRelaxation};
use fs_thermochem::SpeciesId;

fn molecular(count: usize) -> RelaxationState {
    let gas = gas(293.15);
    let species = SpeciesId::new("synthetic").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: (0..count)
            .map(|i| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * i as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * [0.0002, 0.003][i],
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap();
    law.evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
        .unwrap()
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?} at {tolerance}"
    );
}

#[test]
fn molecular_coefficients_and_volume_heat_follow_the_original_reservoir_equation() {
    let state = molecular(2);
    let unused = molecular(0);
    let (complex, mut positions) = single_tet();
    for p in &mut positions {
        *p = p.map(|x| x * 0.03);
    }
    let geometry = HelmholtzFe::new(&complex, &positions, 1, medium(state.gas())).unwrap();
    for frequency in [20.0, 200.0, 2000.0] {
        let omega = core::f64::consts::TAU * frequency;
        let model = PorousRegionsInterior::new(
            &geometry,
            frequency,
            &[
                AcousticVolumeRegion::MolecularGas(&state),
                AcousticVolumeRegion::MolecularGas(&unused),
            ],
            &[0],
            &[],
            &[],
            EquivalentFluidControls::default(),
            &mut || false,
        )
        .unwrap();
        assert!(
            matches!(model.sources()[0], AcousticVolumeRegion::MolecularGas(original)
            if core::ptr::eq(original, &state))
        );
        assert!(model.region_state(1).is_none());
        let snapshot = match model.region_state(0).unwrap() {
            AcousticVolumeRegionState::MolecularGas(response) => response,
            _ => panic!("free molecular gas was replaced"),
        };
        assert_eq!(
            snapshot,
            PorousFluidResponse::relaxing_fluid(&state, omega).unwrap()
        );
        assert_eq!(snapshot.angular_frequency_rad_s(), omega);
        assert_eq!(snapshot.flow_resistivity_pa_s_m2(), 0.0);

        // Independent +i omega Debye sum; no attenuation coefficient, fitted
        // modulus or production Debye evaluator enters this constitutive oracle.
        let cv = state.processes().iter().fold(
            C64::from_re(state.frozen_cv_j_kg_k()),
            |sum, process| {
                sum + C64::from_re(process.heat_capacity_j_kg_k)
                    / C64::new(1.0, omega * process.relaxation_time_s)
            },
        );
        let cp = cv + C64::from_re(state.gas().specific_gas_constant);
        let bulk = (cp / cv).scale(state.gas().pressure);
        let coefficients = model.region_state(0).unwrap().coefficients();
        assert_eq!(
            coefficients.dynamic_density_kg_m3,
            C64::from_re(state.gas().density)
        );
        near(coefficients.dynamic_bulk_modulus_pa, bulk, 2e-14);
        assert!(coefficients.dynamic_bulk_modulus_pa.im > 0.0);
        near(
            snapshot.characteristic().1,
            state.relaxing_wavenumber(omega).unwrap().conj(),
            2e-14,
        );

        // A spatially uniform compression still dissipates through the
        // molecular reservoirs. Local heat is not inferred from propagation
        // attenuation or a boundary-power difference.
        let pressure = C64::new(30.0, -7.0);
        let powers = model
            .operator()
            .element_dissipation(&vec![pressure; geometry.ndof()], &mut || false)
            .unwrap();
        assert_eq!(powers.len(), 1);
        assert_eq!(powers[0].density_dissipation_w, 0.0);
        let expected = -0.5 * omega * bulk.recip().im * pressure.norm_sq() * 0.03_f64.powi(3) / 6.0;
        assert!(expected > 0.0);
        assert!((powers[0].compressional_dissipation_w / expected - 1.0).abs() < 3e-13);
    }

    let lossless = PorousRegionsInterior::new(
        &geometry,
        800.0,
        &[AcousticVolumeRegion::MolecularGas(&unused)],
        &[0],
        &[],
        &[],
        EquivalentFluidControls::default(),
        &mut || false,
    )
    .unwrap();
    let heat = lossless
        .operator()
        .element_dissipation(&vec![C64::ONE; geometry.ndof()], &mut || false)
        .unwrap();
    assert_eq!(heat[0].total_dissipation_w(), 0.0);
    near(
        lossless
            .region_state(0)
            .unwrap()
            .coefficients()
            .dynamic_bulk_modulus_pa,
        C64::from_re(
            unused
                .gas()
                .acoustic_state()
                .unwrap()
                .adiabatic_bulk_modulus_pa(),
        ),
        2e-14,
    );
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(PorousFluidResponse::relaxing_fluid(&state, omega).is_err());
    }
}

#[test]
fn molecular_gas_and_relaxing_porous_volumes_converge_to_original_layer_solution() {
    let state = molecular(2);
    let porous = law();
    let frequency = 800.0;
    let omega = core::f64::consts::TAU * frequency;
    let layers = RigidBackedLayers::new(
        state.gas(),
        omega,
        &[
            Layer::RelaxingFluid {
                relaxation: &state,
                thickness_m: 0.05,
            },
            Layer::RelaxingPorous {
                law: &porous,
                relaxation: &state,
                thickness_m: 0.05,
            },
        ],
    )
    .unwrap();
    let [p, v] = layers
        .reflection(1.0)
        .unwrap()
        .normalized_surface_coordinates();
    let expected = (p / v).scale(state.gas().characteristic_impedance);
    let sources = [
        AcousticVolumeRegion::MolecularGas(&state),
        AcousticVolumeRegion::RelaxingPorous {
            law: &porous,
            relaxation: &state,
        },
    ];
    let mut errors = Vec::new();
    for nx in [4, 8, 16] {
        let (complex, positions) = brick(nx, 0.1);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, medium(state.gas())).unwrap();
        let regions: Vec<_> = complex
            .tets
            .iter()
            .map(|tet| {
                let x = tet.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 4.0;
                usize::from(x > 0.05)
            })
            .collect();
        let model = PorousRegionsInterior::new(
            &geometry,
            frequency,
            &sources,
            &regions,
            &[],
            &[],
            EquivalentFluidControls::default(),
            &mut || false,
        )
        .unwrap();
        let faces = inlet(&geometry);
        let result = model
            .operator()
            .solve(&drive(&faces), &mut || false)
            .unwrap();
        let actual = mean_pressure(&geometry, &faces, &result.pressure_pa);
        errors.push((actual - expected).abs() / expected.abs());
        for (&region, power) in regions.iter().zip(&result.element_power) {
            assert!(power.compressional_dissipation_w > 0.0);
            if region == 0 {
                assert_eq!(power.density_dissipation_w, 0.0);
            } else {
                assert!(power.density_dissipation_w > 0.0);
            }
        }
        assert_eq!(result.boundary_dissipation_w, 0.0);
        assert!(result.power_balance_residual_w().abs() < 3e-10 * result.input_power_w);
        assert_eq!(
            model
                .region_state(0)
                .unwrap()
                .coefficients()
                .dynamic_bulk_modulus_pa,
            layers
                .material_response(0)
                .unwrap()
                .dynamic_bulk_modulus_pa()
        );
        assert_eq!(
            model
                .region_state(1)
                .unwrap()
                .coefficients()
                .dynamic_bulk_modulus_pa,
            layers
                .material_response(1)
                .unwrap()
                .dynamic_bulk_modulus_pa()
        );
    }
    assert!(
        errors[1] < errors[0] / 4.0 && errors[2] < errors[1] / 4.0,
        "{errors:?}"
    );
    assert!(errors[2] < 5e-4, "{errors:?}");
}

#[test]
fn molecular_compression_and_mechanical_wall_heat_share_one_complete_pressure_solve() {
    let state = molecular(2);
    let (complex, mut positions) = single_tet();
    for p in &mut positions {
        *p = p.map(|x| x * 0.03);
    }
    let geometry = HelmholtzFe::new(&complex, &positions, 2, medium(state.gas())).unwrap();
    let faces = geometry.boundary_faces();
    let wall = MechanicalBoundary {
        faces: vec![faces[0]],
        wall: WallPin {
            surface_density: 0.35,
            stiffness_per_area: 2e8,
            resistance: 500.0,
        },
        relaxation: vec![AcousticWallRelaxation {
            stiffness_per_area: 5e7,
            relaxation_time_s: 0.002,
        }],
    };
    let model = PorousRegionsInterior::new_with_coupled_mechanical_walls(
        &geometry,
        800.0,
        &[AcousticVolumeRegion::MolecularGas(&state)],
        &[0],
        &[],
        &[wall],
        &[],
        EquivalentFluidControls::default(),
        &mut || false,
    )
    .unwrap();
    let forcing = [(faces[1], [C64::from_re(0.001); 3])];
    let result = model.operator().solve(&forcing, &mut || false).unwrap();
    assert_eq!(result.walls.len(), 1);
    assert_eq!(result.density_dissipation_w, 0.0);
    assert!(result.compressional_dissipation_w > 0.0);
    assert!(result.wall_accounting.resistance_dissipation_w > 0.0);
    assert!(result.wall_accounting.relaxation_dissipation_w > 0.0);
    assert_eq!(result.boundary_dissipation_w, 0.0);
    assert!(result.power_balance_residual_w().abs() < 1e-8 * result.input_power_w);
    assert!(result.fluid_power_balance_residual_w().abs() < 1e-8 * result.input_power_w);
    let volume = model
        .operator()
        .element_dissipation(&result.pressure_pa, &mut || false)
        .unwrap();
    assert_eq!(volume, result.element_power);
    assert_eq!(
        model
            .operator()
            .solve(&forcing, &mut || false)
            .unwrap()
            .pressure_pa,
        result.pressure_pa
    );
    assert!(
        matches!(model.sources()[0], AcousticVolumeRegion::MolecularGas(original)
        if core::ptr::eq(original, &state))
    );
    assert!(matches!(
        model.operator().solve(&forcing, &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
}
