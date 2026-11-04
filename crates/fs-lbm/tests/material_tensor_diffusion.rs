//! Synthetic, explicitly oriented material claims on the actual thermal kernels.
//! These establish numerical/source coupling, not measured material validation.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::material::transport::tensor::{MaterialThermalTensorError, ThermalMaterialField};
use fs_lbm::material::{
    MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3,
    plan_scaling_from_material,
};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::thermal::transport::{
    PreparedThermalTransport, ThermalTransportSourceQuery,
    orthotropic::PreparedOrthotropicThermalTransport,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
}
fn paths<const D: usize>(
    conductivity: [[f64; 2]; D],
    span: [f64; 2],
    capacity_excursion: bool,
) -> [PreparedThermalTransport; D] {
    let density = PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0]));
    let cp = PropertyKey::new("specific_heat", Dims([2, 0, -2, -1, 0, 0]));
    let k: [_; D] = core::array::from_fn(|axis| {
        PropertyKey::new(
            format!("principal_conductivity_{axis}"),
            Dims([1, 1, -3, -1, 0, 0]),
        )
    });
    let mut claims = ClaimSet::new();
    let mut add = |key: &PropertyKey, knots: Vec<(f64, f64)>| {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims: key.quantity().dims(),
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity(),
                    span[0],
                    span[1],
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored orthotropic transport numerical fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(&density, vec![(span[0], 1000.0), (span[1], 1000.0)]);
    let capacities = if capacity_excursion {
        vec![
            (span[0], 2000.0),
            (315.0, 2000.0),
            (320.0, 2400.0),
            (325.0, 2000.0),
            (span[1], 2000.0),
        ]
    } else {
        vec![(span[0], 2000.0), (span[1], 2000.0)]
    };
    add(&cp, capacities);
    for axis in 0..D {
        add(
            &k[axis],
            vec![
                (span[0], conductivity[axis][0]),
                (span[1], conductivity[axis][1]),
            ],
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic orthotropic medium".into(),
            phase: "declared medium".into(),
            process: "numerical source fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    core::array::from_fn(|axis| {
        PreparedThermalTransport::try_new(
            &card,
            ThermalTransportSourceQuery {
                density: &density,
                specific_heat: &cp,
                conductivity: &k[axis],
                temperature_axis: "T",
                lower: &point(span[0]),
                upper: &point(span[1]),
                selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
            },
        )
        .unwrap()
    })
}
fn identity<const D: usize>() -> [[f64; D]; D] {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }))
}
fn rotated2() -> [[f64; 2]; 2] {
    [[0.6, -0.8], [0.8, 0.6]]
}
fn limits() -> HeatSourceLimits {
    HeatSourceLimits {
        max_cells: 64,
        temperature_range_k: [280.0, 340.0],
        balance_tolerance_j: 1e-8,
    }
}
fn slab2() -> MaterialSlab {
    let material = FluidTransportState::try_new(1000.0, 2000.0, 1000.0, 1.0).unwrap();
    plan_scaling_from_material(&material, 0.01, 0.01, 0.1)
        .unwrap()
        .thermal_slab(MaterialThermalSlab {
            nx: 4,
            ny: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: 32,
        })
        .unwrap()
}
fn slab3() -> MaterialSlab3 {
    let material = FluidTransportState::try_new(1000.0, 2000.0, 1000.0, 1.0).unwrap();
    plan_scaling_from_material(&material, 0.01, 0.01, 0.1)
        .unwrap()
        .thermal_slab3(MaterialThermalSlab3 {
            nx: 2,
            ny: 2,
            nz: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: 32,
        })
        .unwrap()
}
fn power(value: f64) -> Vec<f64> {
    let mut values = vec![0.0; 32];
    values[4..28].fill(value);
    values
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn rotated_source_tensor_preserves_directional_values_receipts_and_support() {
    let original = paths([[1000.0, 1600.0], [3000.0, 1800.0]], [280.0, 340.0], false);
    let law = PreparedOrthotropicThermalTransport::try_new(original.clone(), rotated2()).unwrap();
    assert_eq!(law.principal_paths(), &original);
    assert_eq!(law.constant_storage(), Some([1000.0, 2000.0]));
    let sample = law.at_temperature(310.0).unwrap();
    assert_eq!(sample.path_identity(), law.identity());
    let tensor = sample.conductivity_w_m_k();
    // Independently rotated principal values k1=1300, k2=2400 W/(m K).
    close(tensor[0][0], 0.36 * 1300.0 + 0.64 * 2400.0, 1e-12);
    close(tensor[1][1], 0.64 * 1300.0 + 0.36 * 2400.0, 1e-12);
    close(tensor[0][1], 0.48 * (1300.0 - 2400.0), 1e-12);
    assert_eq!(tensor[1][0], tensor[0][1]);
    let gradient = [2.0, -3.0];
    let heat_flux = tensor.map(|row| -(row[0] * gradient[0] + row[1] * gradient[1]));
    assert!(heat_flux[0] * gradient[0] + heat_flux[1] * gradient[1] < 0.0);
    for i in 0..2 {
        for j in 0..2 {
            close(sample.diffusivity_m2_s()[i][j], tensor[i][j] / 2e6, 1e-18);
        }
    }
    assert!(law.at_temperature(279.0).is_err());
    assert!(law.at_temperature(341.0).is_err());
    assert!(law.at_temperature(f64::NAN).is_err());
    let unrotated =
        PreparedOrthotropicThermalTransport::try_new(original.clone(), identity()).unwrap();
    assert_ne!(unrotated.identity(), law.identity());
    assert!(
        PreparedOrthotropicThermalTransport::try_new(original.clone(), [[1.0, 0.1], [0.0, 1.0]])
            .is_err()
    );
    assert!(
        PreparedOrthotropicThermalTransport::try_new(original.clone(), [[1.0, 0.0], [0.0, -1.0]])
            .is_err()
    );
    let other = paths([[1000.0, 1600.0], [3100.0, 1800.0]], [280.0, 340.0], false);
    assert!(
        PreparedOrthotropicThermalTransport::try_new(
            [original[0].clone(), other[1].clone()],
            identity(),
        )
        .is_err()
    );
}

#[test]
fn physical_isotropic_limit_matches_existing_heat_steps_in_both_dimensions() {
    let gate = CancelGate::new();
    let p = power(1e6);
    let source2 = PreparedOrthotropicThermalTransport::try_new(
        paths([[1000.0; 2]; 2], [280.0, 340.0], false),
        identity(),
    )
    .unwrap();
    let source3 = PreparedOrthotropicThermalTransport::try_new(
        paths([[1000.0; 2]; 3], [280.0, 340.0], false),
        identity(),
    )
    .unwrap();
    let (mut tensor2, mut scalar2) = (slab2(), slab2());
    let (mut tensor3, mut scalar3) = (slab3(), slab3());
    for _ in 0..5 {
        let scalar = scalar2
            .step_with_power_density(&p, &[[1.0, 0.0]; 32], 0.01, limits())
            .unwrap();
        let actual = tensor2
            .step_with_material_tensor_diffusion(
                &p,
                &[[1.0, 0.0]; 32],
                0.01,
                ThermalMaterialField::Uniform(&source2),
                limits(),
                &gate,
            )
            .unwrap();
        close(actual.heat.source_j, 0.24, 1e-15);
        close(actual.heat.stored_change_j, scalar.stored_change_j, 1e-11);
        assert!(actual.heat.defect_j.abs() < 1e-10);
        for (a, b) in tensor2
            .lattice
            .gpop
            .iter()
            .flatten()
            .zip(scalar2.lattice.gpop.iter().flatten())
        {
            close(*a, *b, 2e-14);
        }
        let scalar = scalar3
            .step_with_power_density(&p, &[[1.0, 0.0, 0.0]; 32], limits())
            .unwrap();
        let actual = tensor3
            .step_with_material_tensor_diffusion(
                &p,
                &[[1.0, 0.0, 0.0]; 32],
                ThermalMaterialField::Uniform(&source3),
                limits(),
                &gate,
            )
            .unwrap();
        close(actual.heat.stored_change_j, scalar.stored_change_j, 1e-11);
        assert!(actual.heat.defect_j.abs() < 1e-10);
        for (a, b) in tensor3
            .lattice
            .gpop
            .iter()
            .flatten()
            .zip(scalar3.lattice.gpop.iter().flatten())
        {
            close(*a, *b, 2e-14);
        }
    }
}

#[test]
fn temperature_and_cell_assignments_reach_the_real_tensor_collision_and_si_balance() {
    let original = paths([[1000.0, 2200.0], [3000.0, 1800.0]], [280.0, 340.0], false);
    let a = PreparedOrthotropicThermalTransport::try_new(original.clone(), rotated2()).unwrap();
    let b = PreparedOrthotropicThermalTransport::try_new(original, identity()).unwrap();
    let assigned: Vec<_> = (0..32)
        .map(|cell| if cell % 4 < 2 { &a } else { &b })
        .collect();
    let (mut slab, mut isotropic) = (slab2(), slab2());
    slab.lattice.gpop[13][0] += 0.1;
    isotropic.lattice.gpop[13][0] += 0.1;
    for _ in 0..4 {
        let temperatures: Vec<_> = slab
            .lattice
            .gpop
            .iter()
            .map(|p| slab.temperature_si(p.iter().sum()).unwrap())
            .collect();
        let result = slab
            .step_with_material_tensor_diffusion(
                &power(1e6),
                &[[0.0; 2]; 32],
                0.02,
                ThermalMaterialField::PerCell(&assigned),
                limits(),
                &CancelGate::new(),
            )
            .unwrap();
        close(result.heat.source_j, 0.48, 2e-15);
        assert!(result.heat.defect_j.abs() < 1e-10);
        for sample in &result.applied {
            let t = temperatures[sample.cell];
            close(sample.material.temperature_k(), t, 0.0);
            assert_eq!(
                sample.material.path_identity(),
                assigned[sample.cell].identity()
            );
            let k1 = 1000.0 + 20.0 * (t - 280.0);
            let k2 = 3000.0 - 20.0 * (t - 280.0);
            let expected_xy = if sample.cell % 4 < 2 {
                0.48 * (k1 - k2)
            } else {
                0.0
            };
            close(
                sample.material.conductivity_w_m_k()[0][1],
                expected_xy,
                1e-11,
            );
            close(
                sample.diffusion.lattice_diffusivity()[0][1],
                expected_xy / 2e6 * 100.0,
                1e-15,
            );
        }
        for sample in &result.endpoint {
            close(
                sample.material.temperature_k(),
                slab.temperature_si(slab.lattice.gpop[sample.cell].iter().sum())
                    .unwrap(),
                0.0,
            );
        }
        isotropic
            .step_with_power_density(&power(1e6), &[[0.0; 2]; 32], 0.02, limits())
            .unwrap();
    }
    assert!(
        slab.lattice
            .gpop
            .iter()
            .flatten()
            .zip(isotropic.lattice.gpop.iter().flatten())
            .any(|(a, b)| (a - b).abs() > 1e-6)
    );
}

#[test]
fn interior_storage_excursions_endpoint_source_refusal_and_cancellation_are_atomic() {
    let constant = PreparedOrthotropicThermalTransport::try_new(
        paths([[1000.0; 2]; 2], [300.0, 310.0], false),
        identity(),
    )
    .unwrap();
    let variable = PreparedOrthotropicThermalTransport::try_new(
        paths([[1000.0; 2]; 2], [280.0, 340.0], true),
        identity(),
    )
    .unwrap();
    assert_eq!(variable.constant_storage(), None);
    // All initially sampled temperatures have Cp=2000; an unsupported future
    // interior excursion must still reject the frozen-storage model now.
    for t in [300.0, 305.0, 310.0, 340.0] {
        assert_eq!(
            variable.at_temperature(t).unwrap().specific_heat_j_kg_k(),
            2000.0
        );
    }
    let (mut slab, mut replay) = (slab2(), slab2());
    let before = (
        slab.lattice.grid.f.clone(),
        slab.lattice.gpop.clone(),
        slab.lattice.grid.fext.clone(),
        slab.lattice.grid.tau.clone(),
    );
    assert!(matches!(
        slab.step_with_material_tensor_diffusion(
            &power(0.0),
            &[[1.0, 0.0]; 32],
            0.01,
            ThermalMaterialField::Uniform(&variable),
            limits(),
            &CancelGate::new(),
        ),
        Err(MaterialThermalTensorError::Invalid(_))
    ));
    assert!(matches!(
        slab.step_with_material_tensor_diffusion(
            &power(2e9),
            &[[1.0, 0.0]; 32],
            0.01,
            ThermalMaterialField::Uniform(&constant),
            limits(),
            &CancelGate::new(),
        ),
        Err(MaterialThermalTensorError::Source { .. })
    ));
    let cancel = CancelGate::new();
    cancel.request();
    assert!(matches!(
        slab.step_with_material_tensor_diffusion(
            &power(0.0),
            &[[1.0, 0.0]; 32],
            0.01,
            ThermalMaterialField::Uniform(&constant),
            limits(),
            &cancel,
        ),
        Err(MaterialThermalTensorError::Cancelled)
    ));
    assert_eq!(
        (
            &slab.lattice.grid.f,
            &slab.lattice.gpop,
            &slab.lattice.grid.fext,
            &slab.lattice.grid.tau
        ),
        (&before.0, &before.1, &before.2, &before.3)
    );
    for target in [&mut slab, &mut replay] {
        target
            .step_with_material_tensor_diffusion(
                &power(0.0),
                &[[1.0, 0.0]; 32],
                0.01,
                ThermalMaterialField::Uniform(&constant),
                limits(),
                &CancelGate::new(),
            )
            .unwrap();
    }
    assert_eq!(slab.lattice.gpop, replay.lattice.gpop);
    assert_eq!(slab.lattice.grid.f, replay.lattice.grid.f);
}

#[test]
fn rotated_three_dimensional_conductivity_retains_cross_plane_heat_transport() {
    let rotation = [[0.6, 0.0, -0.8], [0.0, 1.0, 0.0], [0.8, 0.0, 0.6]];
    let law = PreparedOrthotropicThermalTransport::try_new(
        paths(
            [[1000.0; 2], [2000.0; 2], [3000.0; 2]],
            [280.0, 340.0],
            false,
        ),
        rotation,
    )
    .unwrap();
    let k = law.at_temperature(305.0).unwrap().conductivity_w_m_k();
    close(k[0][2], -960.0, 1e-12);
    assert_eq!(k[2][0], k[0][2]);
    assert_eq!(k[1][0], 0.0);
    let (mut tensor, mut isotropic) = (slab3(), slab3());
    tensor.lattice.gpop[12][0] += 0.1;
    isotropic.lattice.gpop[12][0] += 0.1;
    for _ in 0..4 {
        let result = tensor
            .step_with_material_tensor_diffusion(
                &power(1e6),
                &[[0.0; 3]; 32],
                ThermalMaterialField::Uniform(&law),
                limits(),
                &CancelGate::new(),
            )
            .unwrap();
        close(result.heat.source_j, 0.24, 1e-15);
        assert!(result.heat.defect_j.abs() < 1e-10);
        close(
            result.applied[0].diffusion.lattice_diffusivity()[0][2],
            -0.048,
            1e-16,
        );
        isotropic
            .step_with_power_density(&power(1e6), &[[0.0; 3]; 32], limits())
            .unwrap();
    }
    assert!(
        tensor
            .lattice
            .gpop
            .iter()
            .flatten()
            .zip(isotropic.lattice.gpop.iter().flatten())
            .any(|(a, b)| (a - b).abs() > 1e-6)
    );
}
