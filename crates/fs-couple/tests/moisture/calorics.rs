//! Native finite-vapor moisture and independent caloric reference balances.
//! All constitutive values below are authored verification data.
use super::atmosphere as air;
use super::*;
use fs_blake3::ContentHash;
use fs_couple::moisture::calorics::*;
use fs_couple::moisture::cartesian::{MoistureCartesian, MoistureCartesianGeometry};
use fs_material::state_point::{
    MOISTURE_ISOSTERIC_HEAT_DIMS as ENERGY, MOISTURE_ISOSTERIC_HEAT_PROPERTY as NAME,
    MoistureSorptionEnergy, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn source(q: [f64; 2]) -> MaterialCard {
    let base = card([0.01; 2], 0.3);
    let mut claims = base.claims().clone();
    let mut heat = base.claims_for("moisture_diffusivity")[0].1.clone();
    heat.key = PropertyKey::new(NAME, ENERGY);
    heat.value = PropertyValue::Curve {
        abscissa: "moisture_content".into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(0.0, q[0]), (0.3, q[1])],
        dims: ENERGY,
    };
    heat.provenance.source = "authored total differential isosteric heat; J/kg water".into();
    claims.insert_claim(heat).unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}
fn law(c: &MaterialCard, bounds: [f64; 2]) -> MoistureSorptionEnergy {
    let bound = resolve_material_state_point(
        c,
        &point(0.1),
        &[ScalarPropertyRequirement::try_new(NAME, ENERGY, ScalarAdmissibility::Finite).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MoistureSorptionEnergy::from_binding(c, &bound, NAME, bounds).unwrap()
}
fn energy_controls() -> SorptionEnergyControls {
    SorptionEnergyControls {
        maximum_source_segments: 256,
        energy_balance_tolerance_j: 1e-6,
    }
}
fn rt(atmosphere: &fs_couple::moisture::atmosphere::FiniteMoistureAtmosphere) -> f64 {
    let catalog = atmosphere.source_catalog();
    let id = catalog.water_species().unwrap();
    let mass = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|s| s.species() == id)
        .unwrap()
        .molar_mass()
        .value();
    R * atmosphere.gas().state().state().temperature / mass
}

#[test]
fn native_uptake_and_desorption_have_correct_rigid_system_thermostat_duty() {
    with_cx(|cx| {
        let q = 2.5e6;
        let material = source([q; 2]);
        let energy = law(&material, [0.0, 0.3]);
        let field = slab(material, 1.0, vec![0.1], cx);
        let boundaries = [Some(air::sorption(0.25)), Some(air::sorption(0.25))];
        for humidity in [0.8, 0.2] {
            let atmosphere = air::atmosphere(humidity, 1.0);
            let native = atmosphere
                .advance(&field, 2.0, &boundaries, config(), air::controls(), cx)
                .unwrap();
            let result = atmosphere
                .advance_with_sorption_energy(
                    &field,
                    2.0,
                    &boundaries,
                    config(),
                    air::controls(),
                    &energy,
                    energy_controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(
                result.native().moisture().moisture(),
                native.moisture().moisture()
            );
            assert_eq!(result.native().atmosphere(), native.atmosphere());
            let dm = field.dry_mass_kg() * (native.moisture().moisture()[0] - 0.1);
            let expected = (rt(&atmosphere) - q) * dm;
            let calorics = result.energy();
            near(calorics.thermostat_heat_j, expected, 2e-9);
            near(calorics.mean_thermostat_power_w, expected / 2.0, 1e-9);
            assert_eq!(calorics.thermostat_heat_j.signum(), -dm.signum());
            // Native gas energy is retained rather than replaced by a presumed
            // exact water cancellation. It closes against bound storage.
            assert_eq!(
                calorics.gas_internal_energy_change_j,
                native.step().gas_internal_energy_change_j
            );
            near(
                calorics.gas_internal_energy_change_j + calorics.bound_internal_energy_change_j
                    - calorics.thermostat_heat_j,
                0.0,
                1e-6,
            );
            assert!(calorics.energy_balance_residual_j.abs() < 1e-6);
            assert_eq!(
                calorics.vapor_calorics.species(),
                atmosphere.source_catalog().water_species().unwrap()
            );
            assert_eq!(calorics.sorption_energy_source, energy.identity());
            assert_eq!(calorics.source_segments, 1);
            assert_eq!(field.time_s(), 0.0);
            assert_eq!(field.moisture(), &[0.1]);
        }
    });
}

#[test]
fn sealed_internal_redistribution_changes_energy_at_zero_net_vapor_uptake() {
    with_cx(|cx| {
        // q(M)=3e6-2e6*M. The bound potential has the nonlinear term +1e6*M².
        let material = source([3e6, 2.4e6]);
        let energy = law(&material, [0.0, 0.3]);
        let field = slab(material, 1.0, vec![0.05, 0.25], cx);
        let atmosphere = air::atmosphere(0.5, 1.0);
        let result = atmosphere
            .advance_with_sorption_energy(
                &field,
                1.0,
                &[None, None],
                config(),
                air::controls(),
                &energy,
                energy_controls(),
                cx,
            )
            .unwrap();
        let decay = 1.0 / 1.08; // Independent two-cell Neumann backward Euler.
        let expected = [0.15 - 0.1 * decay, 0.15 + 0.1 * decay];
        for (&actual, target) in result.native().moisture().moisture().iter().zip(expected) {
            near(actual, target, 1e-14);
        }
        let mass = field.dry_mass_kg() / 2.0;
        let expected_heat = mass
            * 1e6
            * (expected[0] * expected[0] + expected[1] * expected[1] - 0.05 * 0.05 - 0.25 * 0.25);
        near(result.energy().thermostat_heat_j, expected_heat, 2e-9);
        near(result.native().step().water_to_solid_kg, 0.0, 1e-14);
        near(
            result.native().step().gas_internal_energy_change_j,
            0.0,
            1e-7,
        );
        assert!(result.energy().thermostat_heat_j < -100.0);
        assert_eq!(result.energy().cell_binding_energy_change_j.len(), 2);
        // A single mean-moisture model would miss this entire heat duty.
        near(
            result.native().moisture().water_mass_kg(),
            field.water_mass_kg(),
            1e-14,
        );
    });
}

#[test]
fn six_face_native_diffusion_integrates_each_real_cartesian_cell() {
    with_cx(|cx| {
        let material = source([3e6, 2.4e6]);
        let energy = law(&material, [0.0, 0.3]);
        let field = MoistureCartesian::new(
            material,
            point(0.1),
            ["moisture_diffusivity"; 3],
            MaterialPropertySelection::SingleClaimOnly,
            MoistureCartesianGeometry {
                lengths_m: [0.02, 0.03, 0.04],
                cells: [2; 3],
                frame: ContentHash([4; 32]),
                maximum_moisture: 0.3,
                maximum_cells: 8,
            },
            vec![0.1; 8],
            cx,
        )
        .unwrap();
        let atmosphere = air::atmosphere(0.8, 0.01);
        let boundaries = [
            Some(air::sorption(0.25)),
            None,
            None,
            None,
            None,
            Some(air::sorption(0.25)),
        ];
        let nonlinear = MoistureNonlinearConfig {
            relative_tolerance: 1e-12,
            maximum_iterations: 128,
            relaxation: 1.0,
        };
        let native = atmosphere
            .advance_cartesian(
                &field,
                0.002,
                &boundaries,
                config(),
                nonlinear,
                air::controls(),
                cx,
            )
            .unwrap();
        let result = atmosphere
            .advance_cartesian_with_sorption_energy(
                &field,
                0.002,
                &boundaries,
                config(),
                nonlinear,
                air::controls(),
                &energy,
                energy_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            result.native().moisture().moisture(),
            native.moisture().moisture()
        );
        assert_eq!(result.native().atmosphere(), native.atmosphere());
        let mass = field.dry_mass_kg() / 8.0;
        let mut expected = 0.0;
        for (i, (&a, &b)) in field
            .moisture()
            .iter()
            .zip(native.moisture().moisture())
            .enumerate()
        {
            // Independent primitive of R*T-q(M), with real per-cell loading.
            let increment = mass * ((rt(&atmosphere) - 3e6) * (b - a) + 1e6 * (b * b - a * a));
            near(
                result.energy().cell_binding_energy_change_j[i],
                increment,
                1e-10,
            );
            expected += increment;
        }
        near(result.energy().thermostat_heat_j, expected, 1e-9);
        assert_eq!(result.energy().source_segments, 8);
        assert!(
            native
                .moisture()
                .moisture()
                .windows(2)
                .any(|pair| pair[0] != pair[1])
        );
        assert_eq!(field.time_s(), 0.0);
    });
}

#[test]
fn late_source_and_work_refusals_and_cancellation_publish_no_transport_endpoint() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let material = source([2.5e6; 2]);
        let energy = law(&material, [0.0, 0.3]);
        let narrow = law(&material, [0.0, 0.10000001]);
        let field = slab(material, 1.0, vec![0.1], &cx);
        let atmosphere = air::atmosphere(0.8, 1.0);
        let saved_gas = atmosphere.clone();
        let boundaries = [Some(air::sorption(0.25)), Some(air::sorption(0.25))];
        let advance = |law, controls| {
            atmosphere.advance_with_sorption_energy(
                &field,
                2.0,
                &boundaries,
                config(),
                air::controls(),
                law,
                controls,
                &cx,
            )
        };
        assert!(matches!(
            advance(&narrow, energy_controls()),
            Err(SorptionEnergyError::Moisture(MoistureError::Material(_)))
        ));
        assert!(matches!(
            advance(
                &energy,
                SorptionEnergyControls {
                    maximum_source_segments: 0,
                    ..energy_controls()
                }
            ),
            Err(SorptionEnergyError::SourceBudget {
                required: 1,
                maximum: 0
            })
        ));
        let foreign = law(&source([2.6e6; 2]), [0.0, 0.3]);
        assert!(matches!(
            advance(&foreign, energy_controls()),
            Err(SorptionEnergyError::Invalid(_))
        ));
        gate.request();
        assert!(matches!(
            advance(&energy, energy_controls()),
            Err(SorptionEnergyError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(atmosphere, saved_gas);
        assert_eq!(field.time_s(), 0.0);
        assert_eq!(field.moisture(), &[0.1]);
    });
}
