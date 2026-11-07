use fs_blake3::ContentHash;
use fs_material::elastic::acoustic::viscoelastic::{
    ViscoelasticSlab, ViscoelasticSlabError, ViscoelasticSlabSource,
    layers::{ViscoelasticLayer, ViscoelasticLayerBudget, ViscoelasticLayers},
};
use fs_material::visco::GeneralizedMaxwell;
use fs_material::{
    MaterialError,
    fluid_acoustics::{FluidAcousticState, FluidAcoustics},
};
use fs_math::c64::C64;
use std::cell::Cell;

struct CountedFluid(Cell<usize>);
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.0.set(self.0.get() + 1);
        FluidAcousticState::try_new(2.0, 10.0)
    }
}
fn source<'a>(
    bulk: &'a GeneralizedMaxwell,
    shear: &'a GeneralizedMaxwell,
    id: u8,
    density: f64,
) -> ViscoelasticSlabSource<'a> {
    ViscoelasticSlabSource {
        bulk,
        shear,
        density_kg_m3: density,
        frequency_band_rad_s: (10.0, 100.0),
        maximum_absolute_strain: 0.1,
        source_identity: ContentHash([id; 32]),
    }
}
fn budget(n: usize) -> ViscoelasticLayerBudget {
    ViscoelasticLayerBudget {
        maximum_layers: n,
        maximum_total_branches: 2 * n,
    }
}
fn near(a: C64, b: C64) {
    assert!((a - b).abs() < 3e-13 * b.abs().max(1e-10), "{a:?} != {b:?}");
}

#[test]
fn g1_original_heterogeneous_laws_own_interface_fields_and_local_heat() {
    let bulk = [
        GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap(),
        GeneralizedMaxwell::new(130.0, vec![(15.0, 0.02)]).unwrap(),
    ];
    let shear = [
        GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap(),
        GeneralizedMaxwell::new(45.0, vec![(27.0, 0.03)]).unwrap(),
    ];
    let layers = [
        ViscoelasticLayer {
            source: source(&bulk[0], &shear[0], 3, 4.0),
            thickness_m: 0.2,
        },
        ViscoelasticLayer {
            source: source(&bulk[1], &shear[1], 8, 7.0),
            thickness_m: 0.3,
        },
    ];
    let fluid = CountedFluid(Cell::new(0));
    let stack =
        ViscoelasticLayers::try_new(&layers, [&fluid; 2], 20.0, budget(2), &mut || false).unwrap();
    assert_eq!(fluid.0.get(), 1);
    assert_eq!(stack.layer_count(), 2);
    assert_eq!(stack.angular_frequency_rad_s(), 20.0);
    assert!(core::ptr::eq(
        stack.exterior_sources()[0],
        &fluid as &dyn FluidAcoustics
    ));
    assert_eq!(
        stack.exterior_states()[0].characteristic_impedance_pa_s_m(),
        20.0
    );
    for i in 0..2 {
        let retained = stack.layer(i).unwrap();
        assert!(core::ptr::eq(retained.source.bulk, &bulk[i]));
        assert!(core::ptr::eq(retained.source.shear, &shear[i]));
        for (actual, law) in stack
            .layer_moduli_pa(i)
            .unwrap()
            .into_iter()
            .zip([&bulk[i], &shear[i]])
        {
            let (real, imaginary) = law.modulus(20.0);
            assert_eq!(actual, C64::new(real, imaginary));
        }
    }
    let input = [C64::new(1.0, 0.2), C64::new(-0.1, 0.3)];
    let drive = stack.drive(input, 0.7, &mut || false).unwrap();
    let left = drive.field(0, 0.2, &mut || false).unwrap();
    let right = drive.field(1, 0.0, &mut || false).unwrap();
    near(left.compressive_traction_pa, right.compressive_traction_pa);
    near(left.normal_velocity_m_s, right.normal_velocity_m_s);
    assert!((left.longitudinal_strain - right.longitudinal_strain).abs() > 1e-4);
    let mut heat = [0.0; 2];
    for (i, id) in [3_u8, 8].into_iter().enumerate() {
        let absorption = drive.layer_absorption(i, &mut || false).unwrap();
        assert_eq!(absorption.source_identity, ContentHash([id; 32]));
        assert!(absorption.maximum_absolute_strain_bound <= 0.1);
        assert!(absorption.bulk_absorbed_power_w > 0.0 && absorption.shear_absorbed_power_w > 0.0);
        heat[0] += absorption.bulk_absorbed_power_w;
        heat[1] += absorption.shear_absorbed_power_w;
    }
    let p = drive.power();
    assert_eq!(heat, [p.bulk_absorbed_power_w, p.shear_absorbed_power_w]);
    assert!(p.balance_residual_w.abs() < 1e-13 * p.incident_power_w.iter().sum::<f64>());
    let s = stack.pressure_scattering();
    let outgoing = drive.outgoing_pressure_pa();
    near(outgoing[0], s[0] * input[0] + s[1] * input[1]);
    near(outgoing[1], s[2] * input[0] + s[3] * input[1]);
    assert_eq!(fluid.0.get(), 1);
    assert!(stack.layer(2).is_none() && stack.layer_moduli_pa(2).is_none());
    assert!(drive.layer_absorption(2, &mut || false).is_err());
    assert!(drive.field(0, 0.21, &mut || false).is_err());
}

#[test]
fn g0_g3_g4_single_layer_identity_and_each_original_domain_remain_binding() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let original = source(&bulk, &shear, 4, 4.0);
    let layer = ViscoelasticLayer {
        source: original,
        thickness_m: 0.3,
    };
    let slab =
        ViscoelasticSlab::try_new(original, [&fluid; 2], 20.0, 0.3, 2, &mut || false).unwrap();
    let stack =
        ViscoelasticLayers::try_new(&[layer], [&fluid; 2], 20.0, budget(1), &mut || false).unwrap();
    let input = [C64::new(0.4, 0.1), C64::new(-0.1, 0.2)];
    let single = slab.drive(input, 0.7, &mut || false).unwrap();
    let full = stack.drive(input, 0.7, &mut || false).unwrap();
    assert_eq!(single.outgoing_pressure_pa(), full.outgoing_pressure_pa());
    assert_eq!(
        single.power().bulk_absorbed_power_w,
        full.power().bulk_absorbed_power_w
    );
    assert_eq!(
        single.field(0.1, &mut || false).unwrap(),
        full.field(0, 0.1, &mut || false).unwrap()
    );
    assert!(ViscoelasticLayers::try_new(&[], [&fluid; 2], 20.0, budget(1), &mut || false).is_err());
    assert!(matches!(
        ViscoelasticLayers::try_new(&[layer; 2], [&fluid; 2], 20.0, budget(1), &mut || false),
        Err(ViscoelasticSlabError::Budget)
    ));
    assert!(matches!(
        ViscoelasticLayers::try_new(
            &[layer; 2],
            [&fluid; 2],
            20.0,
            ViscoelasticLayerBudget {
                maximum_layers: 2,
                maximum_total_branches: 3
            },
            &mut || false
        ),
        Err(ViscoelasticSlabError::Budget)
    ));
    let mut second = layer;
    second.source.frequency_band_rad_s = (25.0, 100.0);
    assert!(matches!(
        ViscoelasticLayers::try_new(&[layer, second], [&fluid; 2], 20.0, budget(2), &mut || {
            false
        }),
        Err(ViscoelasticSlabError::OutsideSourceSupport)
    ));
    second = layer;
    second.source.source_identity = ContentHash([0; 32]);
    assert!(
        ViscoelasticLayers::try_new(&[layer, second], [&fluid; 2], 20.0, budget(2), &mut || {
            false
        })
        .is_err()
    );
    second = layer;
    second.source.maximum_absolute_strain = 1e-12;
    let limited =
        ViscoelasticLayers::try_new(&[layer, second], [&fluid; 2], 20.0, budget(2), &mut || {
            false
        })
        .unwrap();
    assert!(matches!(
        limited.drive(input, 0.7, &mut || false),
        Err(ViscoelasticSlabError::OutsideSourceSupport)
    ));
    let stack =
        ViscoelasticLayers::try_new(&[layer; 2], [&fluid; 2], 20.0, budget(2), &mut || false)
            .unwrap();
    let mut calls = 0;
    assert!(matches!(
        stack.drive(input, 0.7, &mut || {
            calls += 1;
            calls == 4
        }),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    let good = stack.drive(input, 0.7, &mut || false).unwrap();
    assert!(matches!(
        good.layer_absorption(1, &mut || true),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    assert!(matches!(
        good.field(1, 0.1, &mut || true),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    assert!(good.field(1, 0.1, &mut || false).is_ok());
}
