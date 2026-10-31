use super::*;

#[test]
fn original_gas_source_reaches_causal_field_work_without_caloric_shortcuts() {
    let (mesh, points0) = kuhn_cube(1);
    let xyz = points(&points0);
    let gas = gas();
    let original =
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, BAND).unwrap();
    let prepared = original
        .prepare_time(
            1.0 / 48_000.0,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        prepared.inner().medium(),
        source_medium(original.gas(), 0.0).unwrap()
    );
    let step = prepared
        .step(&prepared.zero_state(), &[1.0, 0.0], &[], &mut || false)
        .unwrap();
    assert!(step.work.pressure_work_j > 0.0 && step.work.final_energy_j > 0.0);
    assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.pressure_work_j);
    assert!(
        original
            .prepare_time(
                1.0 / BAND[1],
                ThermoviscousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
    let molecular = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, true),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let prepared = molecular
        .prepare_time(
            1.0 / 48_000.0,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let source = molecular.relaxation().unwrap();
    assert_eq!(
        prepared.frozen_cp_j_kg_k(),
        source.frozen_cv_j_kg_k() + gas.specific_gas_constant
    );
    assert_eq!(
        prepared.caloric_reservoirs().len(),
        source.processes().len()
    );
    for (state, original) in prepared.caloric_reservoirs().iter().zip(source.processes()) {
        assert_eq!(
            state.heat_capacity_j_kg_k.to_bits(),
            original.heat_capacity_j_kg_k.to_bits()
        );
        assert_eq!(
            state.relaxation_time_s.to_bits(),
            original.relaxation_time_s.to_bits()
        );
    }
    // The same material law, including the explicit phasor convention change.
    let frequency = 400.0;
    let mut cp = C64::from_re(prepared.frozen_cp_j_kg_k());
    for state in prepared.caloric_reservoirs() {
        cp = cp
            + C64::new(
                1.0,
                std::f64::consts::TAU * frequency * state.relaxation_time_s,
            )
            .recip()
            .scale(state.heat_capacity_j_kg_k);
    }
    close(
        cp,
        molecular.heat_capacity(frequency, &mut || false).unwrap(),
        cp.abs(),
        4e-15,
    );
    let mut initial = prepared.zero_state();
    // Mesh vertex zero is an isothermal wall coefficient. Its molecular state
    // must be admitted independently of the constrained translational field.
    assert_eq!(initial.temperature_k[0], 0.0);
    initial.molecular_temperature_k[0][0] = 1.0;
    let step = prepared
        .step(&initial, &[0.0, 0.0], &[], &mut || false)
        .unwrap();
    assert!(step.work.molecular_relaxation_heat_j > 0.0);
    assert!(step.work.final_energy_j < step.work.initial_energy_j);
    assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.initial_energy_j);
}

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSourceQuery};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn wall_transport(rho: f64, cp: f64, k: f64) -> PreparedThermalTransport {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let keys = [
        PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
        PropertyKey::new("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0])),
    ];
    let mut claims = ClaimSet::new();
    for (key, value) in keys.iter().zip([rho, cp, k]) {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Scalar {
                    value,
                    dims: key.quantity().dims(),
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic fixed-geometry radial thermal source".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic wall".into(),
            phase: "solid".into(),
            process: "authored numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let lower = QueryPoint::new()
        .with_quantity("T", temperature, 280.0)
        .unwrap();
    let upper = QueryPoint::new()
        .with_quantity("T", temperature, 360.0)
        .unwrap();
    PreparedThermalTransport::try_new(
        &card,
        ThermalTransportSourceQuery {
            density: &keys[0],
            specific_heat: &keys[1],
            conductivity: &keys[2],
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
        },
    )
    .unwrap()
}

#[test]
fn source_wall_material_changes_actual_shared_time_fields_and_physical_heat() {
    use super::super::wall::MaterialWallRegion;
    use fs_feec::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
    use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
        ThermalBacking, ThermalWallDepth,
    };
    let (mesh, points0) = kuhn_cube(1);
    let xyz = points(&points0);
    let gas = gas();
    let geometry = HelmholtzFe::new(
        &mesh,
        &xyz,
        2,
        Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    let face = geometry
        .boundary_faces()
        .into_iter()
        .find(|f| f.iter().all(|&v| xyz[v as usize][1] == 0.0))
        .unwrap();
    let mut heats = Vec::new();
    let mut identities = Vec::new();
    for k in [0.003, 0.012] {
        let wall = wall_transport(1000.0, 2.0, k);
        let sample = wall.at_temperature(gas.temperature).unwrap();
        let source = ThermoviscousInterior::from_relaxation(
            build(&mesh, &xyz, &gas, false),
            source(gas, true, 1.0),
            BAND,
        )
        .unwrap()
        .with_conductive_walls(
            &[MaterialWallRegion {
                faces: &[face],
                sample: &sample,
                depth: ThermalWallDepth::Slab {
                    thickness_m: 0.00002,
                    backing: ThermalBacking::Adiabatic,
                },
                minimum_tangential_scale_m: 0.004,
                geometry: BoundaryLayerGeometry {
                    minimum_curvature_radius_m: f64::INFINITY,
                    minimum_half_gap_m: 0.0001,
                    maximum_layer_ratio: 0.05,
                },
            }],
            &mut || false,
        )
        .unwrap();
        assert_eq!(source.wall_sources(), &[sample]);
        assert!(
            source
                .prepare_time(1e-7, ThermoviscousTimeLimits::default(), &mut || false)
                .is_err()
        );
        let prepared = source
            .prepare_time_with_wall_diffusion(
                1e-7,
                &[2],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(prepared.maximum_wall_frequency_hz(), Some(BAND[1]));
        assert_eq!(
            prepared.inner().conductive_walls()[0]
                .wall
                .transport
                .conductivity_w_m_k,
            k
        );
        assert_eq!(
            prepared.inner().conductive_walls()[0]
                .wall
                .transport
                .volumetric_heat_capacity_j_m3_k,
            sample.volumetric_heat_capacity_j_m3_k()
        );
        assert_eq!(
            prepared.caloric_reservoirs().len(),
            source.relaxation().unwrap().processes().len()
        );
        let mut initial = prepared.zero_state();
        initial.conductive_wall_temperature_k[0][0] = 1.0;
        let step = prepared
            .step(&initial, &[0.0, 0.0], &[], &mut || false)
            .unwrap();
        assert!(step.work.solid_conductive_heat_j > 0.0);
        assert!(step.state.temperature_k.iter().any(|t| *t != 0.0));
        assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.initial_energy_j);
        heats.push(step.work.solid_conductive_heat_j);
        identities.push(sample.path_identity());
    }
    assert_ne!(identities[0], identities[1]);
    assert!(heats[1] > 3.8 * heats[0]);
}
