//! Original liquid identity and native adaptive schedule parity.
use super::*;
use fs_lbm::film::evolution::adaptive::{
    ThinFilmAdaptiveControls, ThinFilmDrive, ThinFilmDriveQuery,
};
use fs_lbm::film::evolution::{
    ThinFilmBoundary, ThinFilmControls, ThinFilmGeometry, ThinFilmLimits, ThinFilmLine,
    ThinFilmOpenEnd, ThinFilmOpenEnds,
};

#[test]
fn adaptive_original_liquid_preserves_source_receipts_and_native_signed_transport() {
    let geometry = ThinFilmGeometry {
        spacing_m: 0.01,
        normal_acceleration_m_s2: 1.0,
        boundary: ThinFilmBoundary::Open,
    };
    let limits = ThinFilmLimits {
        minimum_thickness_m: 5e-5,
        maximum_thickness_m: 2e-4,
        characteristic_length_m: 0.04,
        maximum_height_to_length: 0.01,
        maximum_slope: 0.01,
        maximum_reduced_reynolds: 0.1,
    };
    let controls = ThinFilmControls {
        height_scale_m: 1e-4,
        ..ThinFilmControls::default()
    };
    let adaptive = ThinFilmAdaptiveControls {
        initial_dt_s: 0.01,
        minimum_dt_s: 1e-7,
        maximum_dt_s: 0.01,
        transport_absolute_tolerance_m2: 1e-12,
        energy_absolute_tolerance_j_per_m: 1e-12,
        ..ThinFilmAdaptiveControls::default()
    };
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let card = card(basis, true, false);
        let source = resolve(&card, basis, LiquidCapillarity::Sourced);
        let liquid = MaterialFilmLiquid::try_new(&source).unwrap();
        let mut bound = liquid.thin_film_line(geometry, vec![1e-4; 4], 0.0).unwrap();
        let mut native =
            ThinFilmLine::newtonian(liquid.coefficients(), geometry, vec![1e-4; 4], 0.0).unwrap();
        let resolve = |query: ThinFilmDriveQuery| {
            let end = ThinFilmOpenEnd {
                pressure_pa: 0.1,
                surface_shear_pa: 0.0,
            };
            Ok::<_, std::convert::Infallible>(ThinFilmDrive {
                acceleration_m_s2: if query.segment == 0 { 2.0 } else { -1.0 },
                surface_shear_pa: vec![0.0; 3],
                open_ends: Some(ThinFilmOpenEnds {
                    left: end,
                    right: end,
                }),
                face_temperatures_k: None,
            })
        };
        let actual = bound
            .advance_adaptive(0.02, &[0.01], limits, controls, adaptive, resolve, || false)
            .unwrap();
        let expected = native
            .advance_adaptive(0.02, &[0.01], limits, controls, adaptive, resolve, || false)
            .unwrap();
        assert_eq!(bound.flow().heights_m(), native.heights_m());
        assert_eq!(bound.flow().accepted_time_s(), 0.02);
        assert_eq!(actual.face_integrals, expected.face_integrals);
        assert_eq!(actual.viscous_heat_j_per_m, expected.viscous_heat_j_per_m);
        assert_eq!(actual.endpoint_evaluations, expected.endpoint_evaluations);
        assert_eq!(actual.events_reached, 1);
        assert!(actual.face_integrals[0].transported_volume_m2 > 0.0);
        assert!(core::ptr::eq(bound.liquid().material(), &source));
        assert_eq!(bound.liquid().material().viscosity_basis(), basis);
        assert_eq!(
            bound.liquid().material().material().card_identity(),
            card.content_hash()
        );
    }
}
