//! Source/solver integration, not experimental qualification of water flow.

use super::common_material_acquisition::{compile, point};
use fs_lbm::material::{MaterialBuoyancyBoundary, plan_scaling_from_material};
use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
use fs_matdb_store::{CatalogPack, MaterialStore};
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::liquid::{VOLUMETRIC_EXPANSION_PROPERTY, resolve_liquid_state_with_expansion};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::state_point::MaterialPropertySelection::SingleClaimOnly;
use fs_qty::Temperature;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0e-20),
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={relative}"
    );
}

#[test]
fn g1_g3_iapws_liquid_expansion_reaches_native_thermal_flow() {
    let (pack, _) = compile("water-liquid-iapws-sr6-08");
    let original = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "IAPWS ordinary water".into(),
            phase: "liquid".into(),
            process: "SR6-08(2011), frozen single-phase properties at 0.1 MPa".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let store = MaterialStore::open(":memory:").unwrap();
    store
        .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
        .unwrap();
    store.seal_corpus().unwrap();
    let CatalogPack::MaterialCard(loaded) = store.load_catalog_pack(original.pack_id()).unwrap()
    else {
        panic!("water pack changed family during storage")
    };
    assert_eq!(loaded, original);
    let card = loaded.card();
    let density = card.claims().claims_for("density")[0].1;
    let reference = point(density, &[("temperature", 293.15)]);
    let law = LiquidTemperatureLaw::try_new_with_expansion(
        card,
        &reference,
        "temperature",
        LiquidViscosityBasis::Dynamic,
        SingleClaimOnly,
    )
    .unwrap();

    // Independent printed IAPWS Table 8 check, with its limited decimal precision.
    // v_T/v = rho*v_T; these are the official 298.15 K verification ordinates.
    let verification = law.at_temperature(Temperature::new(298.15)).unwrap();
    close(verification.density_kg_m3(), 997.047_013, 1.0e-9);
    close(
        verification.volumetric_expansion_per_k().unwrap(),
        997.047_013 * 2.580_541_78e-7,
        1.0e-8,
    );

    let (dx, dt, gravity) = (0.001, 0.01, 9.81);
    let mut outcomes = Vec::new();
    for temperature in [293.15, 333.15] {
        let liquid = law.at_temperature(Temperature::new(temperature)).unwrap();
        let thermal = liquid.thermal_transport_state().unwrap();
        let transport = liquid.transport_state().unwrap();
        let beta = thermal.volumetric_expansion_per_k().unwrap();
        assert_eq!(thermal.transport(), transport);
        assert_eq!(thermal.temperature_k(), temperature);
        assert_eq!(liquid.material().card_identity(), card.content_hash());
        assert_eq!(liquid.material().properties().len(), 5);
        for property in liquid.material().properties() {
            let receipt = &property.answer().receipt;
            card.claims().verify_receipt(receipt).unwrap();
            assert_eq!(
                receipt.selected,
                law.reference_state()
                    .material()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert_eq!(receipt.query_point, liquid.material().query_point());
        }
        assert_eq!(
            liquid
                .material()
                .property(VOLUMETRIC_EXPANSION_PROPERTY)
                .unwrap()
                .value_si(),
            beta
        );

        let scaling = plan_scaling_from_material(&liquid, dx, dt, 0.001).unwrap();
        let walls = MaterialBuoyancyBoundary {
            bottom_temperature_k: temperature + 1.0,
            top_temperature_k: temperature - 1.0,
            gravity_m_s2: gravity,
            max_relative_density_change: 0.01,
        };
        let mut two = scaling
            .thermal_slab_from_material([4, 8], walls, 40)
            .unwrap();
        let mut three = scaling
            .thermal_slab3_from_material([4, 2, 8], walls, 80)
            .unwrap();
        let expected_ra = gravity * beta * 2.0 * (8.0_f64 * dx).powi(3)
            / (transport.kinematic_viscosity_m2_s() * transport.thermal_diffusivity_m2_s());
        close(two.lattice.rayleigh(), expected_ra, 1.0e-10);
        let nu_lattice = (scaling.flow().tau - 0.5) / 3.0;
        let alpha3_lattice = (three.lattice.tau_g - 0.5) / 4.0;
        close(
            three.lattice.gbeta * 8.0_f64.powi(3) / (nu_lattice * alpha3_lattice),
            expected_ra,
            1.0e-10,
        );

        // Identical dimensionless perturbation; the sourced beta reaches actual
        // momentum forcing, while the sourced transport controls heat evolution.
        let i2 = two.lattice.grid.idx(0, 4);
        let i3 = three.lattice.grid.idx(0, 0, 4);
        two.lattice.gpop[i2][0] += 0.05;
        three.lattice.gpop[i3][0] += 0.05;
        let delta2 = two.temperature_si(two.lattice.temperature(0, 4)).unwrap() - temperature;
        let delta3 = three
            .temperature_si(three.lattice.temperature(0, 0, 4))
            .unwrap()
            - temperature;
        two.lattice.step();
        three.lattice.step();
        close(
            two.lattice.grid.fext[i2][1],
            gravity * beta * delta2 * dt * dt / dx,
            1.0e-10,
        );
        close(
            three.lattice.grid.fext[i3][2],
            gravity * beta * delta3 * dt * dt / dx,
            1.0e-10,
        );
        let first_force = two.lattice.grid.fext[i2][1];
        for _ in 1..12 {
            two.lattice.step();
            three.lattice.step();
        }
        let response = (
            two.lattice.temperature(0, 4),
            three.lattice.temperature(0, 0, 4),
            first_force,
        );
        assert!(response.0.is_finite() && response.1.is_finite());
        outcomes.push(response);
    }
    assert!((outcomes[1].0 - outcomes[0].0).abs() > 1.0e-8);
    assert!((outcomes[1].1 - outcomes[0].1).abs() > 1.0e-8);
    assert!(outcomes[1].2 > outcomes[0].2 * 2.0);

    for (axis, unsupported) in [
        ("pressure", 100_001.0),
        ("phase-liquid", 0.0),
        ("temperature", 353.16),
    ] {
        let outside = reference.clone().with(axis, unsupported).unwrap();
        assert!(
            resolve_liquid_state_with_expansion(
                card,
                &outside,
                LiquidViscosityBasis::Dynamic,
                SingleClaimOnly,
            )
            .is_err()
        );
        assert!(
            law.reference_state()
                .at_query_point(card, &outside)
                .is_err()
        );
    }
    assert_eq!(
        &law.at_temperature(Temperature::new(293.15)).unwrap(),
        law.reference_state()
    );
}
