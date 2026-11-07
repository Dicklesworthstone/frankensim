use fs_blake3::ContentHash;
use fs_material::MaterialError;
use fs_material::elastic::acoustic::viscoelastic::{
    ViscoelasticSlab, ViscoelasticSlabError, ViscoelasticSlabSource,
};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use std::cell::Cell;

struct CountedFluid {
    calls: Cell<usize>,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.calls.set(self.calls.get() + 1);
        FluidAcousticState::try_new(2.0, 10.0)
    }
}
fn source<'a>(
    bulk: &'a GeneralizedMaxwell,
    shear: &'a GeneralizedMaxwell,
) -> ViscoelasticSlabSource<'a> {
    ViscoelasticSlabSource {
        bulk,
        shear,
        density_kg_m3: 4.0,
        frequency_band_rad_s: (10.0, 100.0),
        maximum_absolute_strain: 0.1,
        source_identity: ContentHash([9; 32]),
    }
}

#[test]
fn g1_original_maxwell_sources_drive_coherent_solid_heat_without_requery() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = CountedFluid {
        calls: Cell::new(0),
    };
    let slab = ViscoelasticSlab::try_new(
        source(&bulk, &shear),
        [&fluid; 2],
        20.0,
        0.3,
        2,
        &mut || false,
    )
    .unwrap();
    assert_eq!(fluid.calls.get(), 1);
    assert!(core::ptr::eq(slab.source().bulk, &bulk));
    assert!(core::ptr::eq(slab.source().shear, &shear));
    assert_eq!(slab.source().source_identity, ContentHash([9; 32]));
    for (actual, original) in slab.moduli_pa().into_iter().zip([&bulk, &shear]) {
        let (storage, loss) = original.modulus(20.0);
        assert_eq!(actual, C64::new(storage, loss));
    }
    let wave = slab
        .drive([C64::new(1.0, 0.2), C64::new(-0.1, 0.3)], 0.7, &mut || {
            false
        })
        .unwrap();
    let power = wave.power();
    assert!(power.bulk_absorbed_power_w > 0.0 && power.shear_absorbed_power_w > 0.0);
    assert!(power.maximum_absolute_strain_bound <= slab.source().maximum_absolute_strain);
    let field = wave.field(0.1, &mut || false).unwrap();
    let modulus = slab.moduli_pa()[0] + slab.moduli_pa()[1].scale(4.0 / 3.0);
    assert!((field.longitudinal_strain + field.compressive_traction_pa / modulus).abs() < 1e-16);
    assert!(field.bulk_heat_density_w_m3 > 0.0 && field.shear_heat_density_w_m3 > 0.0);
    assert_eq!(fluid.calls.get(), 1);
}

#[test]
fn g0_g4_source_domain_strain_budget_and_cancellation_are_atomic() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let original = source(&bulk, &shear);
    for frequency in [0.0, 9.9, 100.1, f64::NAN] {
        assert!(matches!(
            ViscoelasticSlab::try_new(original, [&fluid; 2], frequency, 0.3, 2, &mut || false),
            Err(ViscoelasticSlabError::OutsideSourceSupport)
        ));
    }
    assert!(matches!(
        ViscoelasticSlab::try_new(original, [&fluid; 2], 20.0, 0.3, 1, &mut || false),
        Err(ViscoelasticSlabError::Budget)
    ));
    let mut missing_identity = original;
    missing_identity.source_identity = ContentHash([0; 32]);
    assert!(
        ViscoelasticSlab::try_new(missing_identity, [&fluid; 2], 20.0, 0.3, 2, &mut || false)
            .is_err()
    );
    for frequency in [10.0, 100.0] {
        assert!(
            ViscoelasticSlab::try_new(original, [&fluid; 2], frequency, 0.3, 2, &mut || false)
                .is_ok()
        );
    }
    let slab =
        ViscoelasticSlab::try_new(original, [&fluid; 2], 20.0, 0.3, 2, &mut || false).unwrap();
    assert!(matches!(
        slab.drive([C64::from_re(1000.0), C64::ZERO], 1.0, &mut || false),
        Err(ViscoelasticSlabError::OutsideSourceSupport)
    ));
    let mut checkpoints = 0;
    assert!(matches!(
        slab.drive([C64::ONE, C64::ZERO], 1.0, &mut || {
            checkpoints += 1;
            checkpoints == 2
        }),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    let good = slab
        .drive([C64::ONE, C64::ZERO], 1.0, &mut || false)
        .unwrap();
    assert!(matches!(
        good.field(0.1, &mut || true),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    assert!(good.field(0.1, &mut || false).is_ok());
    let mut invalid_bulk = bulk.clone();
    invalid_bulk.terms[0].1 = -1.0;
    assert!(
        ViscoelasticSlab::try_new(
            source(&invalid_bulk, &shear),
            [&fluid; 2],
            20.0,
            0.3,
            2,
            &mut || false
        )
        .is_err()
    );
}
