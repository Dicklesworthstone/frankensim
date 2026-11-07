//! Discovered planes keep exact thermal sources; geometry is not history transfer.
use super::{curve_card, near, query, resolve};
use fs_material::SmallStrainLaw;
use fs_material::brittle::initiation::TensileInitiationSettings;
use fs_material::brittle::source::ResolvedBrittleCrackBand;
use fs_material::brittle::source::expansion::dynamics::resolve_thermal_brittle_dynamic_material;
use fs_material::state_point::{DENSITY_PROPERTY, MaterialPropertySelection};
use fs_qty::Density;

fn uniaxial_y(stress: f64) -> [f64; 6] {
    [
        -0.2 * stress / 1e9,
        stress / 1e9,
        -0.2 * stress / 1e9,
        0.0,
        0.0,
        0.0,
    ]
}

#[test]
fn discovered_thermal_plane_retains_reference_receipts_and_native_mechanical_history() {
    let source = curve_card();
    let warm = resolve(&source, 300.0, 350.0);
    let start = uniaxial_y(0.2e6);
    let end = uniaxial_y(2e6);
    let event = warm
        .band()
        .tensile_initiator(TensileInitiationSettings {
            maximum_strain: 0.05,
            relative_stress_tolerance: 1e-10,
            fraction_tolerance: 1e-12,
            maximum_bisections: 64,
        })
        .unwrap()
        .between(&start, &end, None)
        .unwrap()
        .unwrap();
    let receiving = warm
        .with_resolved_band(event.resolve_band(0.008).unwrap())
        .unwrap();
    assert_eq!(receiving.band().material(), warm.band().material());
    assert_eq!(receiving.expansion(), warm.expansion());
    assert_eq!(
        receiving.band().law().parameters(),
        warm.band().law().parameters()
    );
    assert!(receiving.band().law().normal()[1].abs() > 1.0 - 1e-14);
    near(receiving.band().law().band_width_m(), 0.008);
    near(receiving.free_linear_strain(), 0.001);

    // Only the original tensile owner may seed the OLD mechanical endpoint.
    // Thermal preparation alone charges no load, eigenstrain work or fracture.
    let seed = event.seed(receiving.band()).unwrap();
    assert_eq!(seed.state.strain(), &start);
    assert_eq!(seed.state.damage(), 0.0);
    assert_eq!(seed.state.dissipated_j_m3(), 0.0);
    let total_start = receiving.total_strain(&seed.state).unwrap();
    for i in 0..6 {
        near(total_start[i], start[i] + if i < 3 { 0.001 } else { 0.0 });
    }
    let total_end = core::array::from_fn(|i| end[i] + if i < 3 { 0.001 } else { 0.0 });
    let fracture = receiving.trial_total(&total_end, &seed.state).unwrap();
    assert!(fracture.state.damage() > 0.0);
    assert!(fracture.state.dissipated_j_m3() > 0.0);
    let native = receiving.band().law().trial(&end, &seed.state).unwrap();
    for i in 0..6 {
        near(fracture.stress_pa[i], native.stress_pa[i]);
    }
    near(
        fracture.state.stored_energy_j_m3(),
        native.state.stored_energy_j_m3(),
    );

    let hot = receiving.at_temperature(&source, 400.0).unwrap();
    let original_hot = warm.at_temperature(&source, 400.0).unwrap();
    assert_eq!(hot.expansion(), original_hot.expansion());
    assert_eq!(hot.band().law().normal(), receiving.band().law().normal());
    assert_eq!(
        hot.band().law().band_width_m(),
        receiving.band().law().band_width_m()
    );
    receiving.require_temperature_successor(&hot).unwrap();
    assert_eq!(warm.band().law().normal(), [1.0, 0.0, 0.0]);
    assert_eq!(warm.band().law().band_width_m(), 0.01);

    // Equal sources still have distinct native owners; preparation is not a
    // route around mechanical-history ownership or irreversible state transfer.
    let foreign = warm.band().law().initial_state();
    assert!(receiving.trial_total(&total_start, &foreign).is_err());
}

#[test]
fn thermal_plane_preparation_refuses_other_sources_contexts_and_controls() {
    let source = curve_card();
    let warm = resolve(&source, 300.0, 350.0);
    let hot = warm.at_temperature(&source, 400.0).unwrap();
    assert_eq!(
        hot.band().law().parameters(),
        warm.band().law().parameters()
    );
    assert_eq!(hot.band().law().elastic(), warm.band().law().elastic());
    assert!(warm.with_resolved_band(hot.band().clone()).is_err());

    let mut claims = super::super::properties();
    claims[0].provenance.source = "independent elastic measurement with identical values".into();
    claims.push(super::super::scalar(
        fs_material::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        fs_material::state_point::INVERSE_TEMPERATURE_DIMS,
        2e-5,
    ));
    let foreign = resolve(&super::super::card(claims), 300.0, 350.0);
    assert_eq!(foreign.band().law().elastic(), warm.band().law().elastic());
    near(foreign.free_linear_strain(), warm.free_linear_strain());
    assert!(warm.with_resolved_band(foreign.band().clone()).is_err());

    for change in [0, 1] {
        let mut controls = warm.band().law().settings();
        if change == 0 {
            controls.maximum_strain = 0.04;
        } else {
            controls.energy_tolerance_j_m3 = 1e-5;
        }
        let other =
            ResolvedBrittleCrackBand::try_from_material(warm.band().material().clone(), controls)
                .unwrap();
        assert!(warm.with_resolved_band(other).is_err());
    }
}

#[test]
fn thermal_plane_preparation_preserves_extra_reference_density_receipts() {
    let mut claims = super::super::properties();
    claims.push(super::super::scalar(
        DENSITY_PROPERTY,
        Density::DIMS,
        2300.0,
    ));
    claims.push(super::super::scalar(
        fs_material::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        fs_material::state_point::INVERSE_TEMPERATURE_DIMS,
        1e-5,
    ));
    let source = super::super::card(claims);
    let warm = resolve_thermal_brittle_dynamic_material(
        &source,
        &query(300.0),
        &query(350.0),
        MaterialPropertySelection::SingleClaimOnly,
        super::super::settings(0.01),
    )
    .unwrap();
    let mut geometry = warm.band().law().settings();
    geometry.normal = [0.0, 0.0, 1.0];
    geometry.band_width_m = 0.008;
    let band =
        ResolvedBrittleCrackBand::try_from_material(warm.band().material().clone(), geometry)
            .unwrap();
    let receiving = warm.with_resolved_band(band).unwrap();
    assert_eq!(receiving.band().material(), warm.band().material());
    assert_eq!(receiving.band().material().properties().len(), 9);
    assert_eq!(
        receiving.band().material().property(DENSITY_PROPERTY),
        warm.band().material().property(DENSITY_PROPERTY),
    );
    let virgin = receiving.band().law().initial_state();
    let free = receiving.free_linear_strain();
    let trial = receiving
        .trial_total(&[free, free, free, 0.0, 0.0, 0.0], &virgin)
        .unwrap();
    assert_eq!(trial.stress_pa, [0.0; 6]);
    assert_eq!(trial.state.damage(), 0.0);
    assert_eq!(
        receiving.band().law().elastic().stress(&[0.0; 6], &()),
        [0.0; 6]
    );
}
