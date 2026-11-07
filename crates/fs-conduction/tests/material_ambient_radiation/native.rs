//! Original material queries and phase tangents reach native gray heat owners.

use super::*;
use fs_conduction::MaterialGraySurfaceLaw;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    GrayAmbientBoundary, GrayAmbientControls, ImplicitBoundaryGeometry, ImplicitGraySurfaceLaw,
};
use fs_thermal::enthalpy::conduction::implicit::link::{GrayRadiativeControls, GrayRadiativeLink};
use fs_thermal::enthalpy::conduction::implicit::{ImplicitConductionLaw, ImplicitConductionState};
use fs_thermal::enthalpy::conduction::{ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{CellId, ContactPatch, DepositionError, ThermalLaw, ThermalLink};

#[test]
fn source_scalars_temperature_and_phase_curves_keep_exact_calorics_and_receipts() {
    for axis in [None, Some("T"), Some("melt")] {
        for datum in [0.0, 1e6] {
            let card = card(axis);
            let curve = curve(&card, datum);
            let context = context();
            let material = binding(&card, &context);
            let bounds = [datum - 100.0, datum + LATENT + 200.0];
            let source = MaterialGraySurfaceLaw::new(
                "native",
                &curve,
                material,
                MaterialPropertySelection::SingleClaimOnly,
                bounds,
            )
            .unwrap();
            assert!(std::ptr::addr_eq(
                source.caloric_law() as *const dyn ThermalLaw,
                &curve as &dyn ThermalLaw as *const dyn ThermalLaw
            ));
            assert_eq!(source.enthalpy_bounds_j_kg(), bounds);
            for (offset, temperature, dt_dh, fraction, df_dh) in [
                (-50.0, 350.0, 1.0, 0.0, 0.0),
                (250.0, 400.0, 0.0, 0.25, 1.0 / LATENT),
                (1100.0, 450.0, 0.5, 1.0, 0.0),
            ] {
                let h = datum + offset;
                let (state, sample) = source.evaluate(h).unwrap();
                let (epsilon, slope) = match axis {
                    None => (0.5, 0.0),
                    Some("T") => (0.2 + 0.003 * (temperature - 300.0), 0.003 * dt_dh),
                    _ => (0.2 + 0.6 * fraction, 0.6 * df_dh),
                };
                close(state.emissivity, epsilon, 2e-15);
                close(state.emissivity_enthalpy_derivative_kg_j, slope, 2e-15);
                assert_eq!(state.convection_w_m2_k, 0.0);
                assert_eq!(state.convection_enthalpy_derivative_kg_w_m2_k_j, 0.0);
                assert_eq!(sample.phase.phase_curve_identity(), curve.identity());
                assert_eq!(sample.material.card_identity(), card.content_hash());
                assert_eq!(source.state_at_specific_enthalpy(h).unwrap(), state);
                let expected = SurfaceEmissivity::from_thermal_state(
                    "native",
                    material,
                    curve.state_at_specific_enthalpy(h).unwrap(),
                    MaterialPropertySelection::PinnedByProperty(vec![(
                        SURFACE_EMISSIVITY_PROPERTY.into(),
                        source.reference_emissivity().receipt().selected,
                    )]),
                )
                .unwrap();
                assert_eq!(sample.emissivity.receipt(), expected.0.receipt());
                assert_eq!(sample.material, expected.1);
            }
            assert!(source.evaluate(bounds[0] - 1.0).is_err());
            assert_eq!(
                source.state_at_specific_enthalpy(bounds[1] + 1.0),
                Err(DepositionError::MaterialDomain)
            );
        }
    }
}

#[test]
fn joint_source_surface_keeps_both_original_axis_orders_and_phase_tangent() {
    for reverse in [false, true] {
        let card = super::surface::rectangle(reverse);
        let curve = curve(&card, 0.0);
        let context = context();
        let source = MaterialGraySurfaceLaw::new(
            "joint",
            &curve,
            binding(&card, &context),
            MaterialPropertySelection::SingleClaimOnly,
            [-100.0, 1200.0],
        )
        .unwrap();
        for (h, t, dt, fraction, df) in [
            (-50.0, 350.0, 1.0, 0.0, 0.0),
            (250.0, 400.0, 0.0, 0.25, 1.0 / LATENT),
            (1100.0, 450.0, 0.5, 1.0, 0.0),
        ] {
            let (state, sample) = source.evaluate(h).unwrap();
            close(
                state.emissivity,
                0.1875 + 0.000625 * (t - 300.0) + 0.5 * fraction,
                2e-15,
            );
            close(
                state.emissivity_enthalpy_derivative_kg_j,
                0.000625 * dt + 0.5 * df,
                2e-15,
            );
            assert_eq!(sample.phase.temperature_k(), t);
            assert_eq!(sample.phase.liquid_mass_fraction(), fraction);
        }
    }
}

// Explicit constant-k test leg; caloric storage is the ORIGINAL material chart.
struct Conductivity<'a>(&'a EquilibriumEnthalpyPhaseCurve);
impl ConductivityLaw for Conductivity<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.0.temperature_k(h)?;
        Ok(2.0)
    }
}
impl ImplicitConductionLaw for Conductivity<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.0
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [-100.0, 1200.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.0.temperature_k(h)?,
            temperature_slope_kg_k_j: self
                .0
                .temperature_derivative_at_specific_enthalpy(h)
                .map_err(|_| DepositionError::MaterialDomain)?,
            conductivity_w_m_k: 2.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

#[test]
fn source_card_changes_native_boundary_and_finite_link_on_latent_plateau() {
    let card = card(Some("melt"));
    let curve = curve(&card, 0.0);
    let conductivity = Conductivity(&curve);
    let context = context();
    let source = MaterialGraySurfaceLaw::new(
        "native",
        &curve,
        binding(&card, &context),
        MaterialPropertySelection::SingleClaimOnly,
        [-100.0, 1200.0],
    )
    .unwrap();
    let boundary =
        GrayAmbientBoundary::new(&conductivity, 0.0, 0.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&source)
            .unwrap();
    let face = ImplicitBoundaryGeometry {
        area_m2: 0.7,
        half_distance_m: 0.02,
        surface_resistance_m2_k_w: 0.01,
    };
    let first = boundary
        .evaluate(250.0, 300.0, face, &mut || Ok(()))
        .unwrap();
    let second = boundary
        .evaluate(750.0, 300.0, face, &mut || Ok(()))
        .unwrap();
    assert!(second.outward_power_w > first.outward_power_w);
    assert!(first.boundary.power_enthalpy_derivative_w_kg_j > 0.0);
    let boundary_fd = (boundary
        .evaluate(250.001, 300.0, face, &mut || Ok(()))
        .unwrap()
        .outward_power_w
        - boundary
            .evaluate(249.999, 300.0, face, &mut || Ok(()))
            .unwrap()
            .outward_power_w)
        / 0.002;
    close(
        first.boundary.power_enthalpy_derivative_w_kg_j,
        boundary_fd,
        2e-8,
    );

    let link = GrayRadiativeLink::new(
        [&conductivity; 2],
        [0.0; 2],
        0.0,
        [0.01; 2],
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .with_surface_laws([&source; 2])
    .unwrap();
    let geometry = ConductiveLink {
        link_id: 1,
        half_distances_m: [0.02; 2],
        contact_resistance_m2_k_w: 0.0,
    };
    let patch = ContactPatch::symmetric(ThermalLink {
        id: 1,
        a: CellId(0),
        b: CellId(1),
        area_m2: 0.7,
        conductance_w_k: 0.0,
    });
    let first = link
        .evaluate([250.0, -50.0], geometry, patch, &mut || Ok(()))
        .unwrap();
    let second = link
        .evaluate([750.0, -50.0], geometry, patch, &mut || Ok(()))
        .unwrap();
    assert!(second.power_w > first.power_w);
    let fixed = GrayRadiativeLink::new(
        [&conductivity; 2],
        [0.35, 0.2],
        0.0,
        [0.01; 2],
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .evaluate([250.0, -50.0], geometry, patch, &mut || Ok(()))
    .unwrap();
    close(first.power_w, fixed.power_w, 2e-12);
    assert_eq!(fixed.link.power_enthalpy_derivatives_w_kg_j[0], 0.0);
    assert!(first.link.power_enthalpy_derivatives_w_kg_j[0] > 0.0);
}

#[test]
fn native_source_refuses_wrong_card_missing_conditions_and_unsupported_bounds() {
    let card = card(Some("melt"));
    let curve = curve(&card, 0.0);
    let context = context();
    for bounds in [
        [-101.0, 1200.0],
        [-100.0, 1201.0],
        [0.0, 0.0],
        [f64::NAN, 100.0],
    ] {
        assert!(
            MaterialGraySurfaceLaw::new(
                "native",
                &curve,
                binding(&card, &context),
                MaterialPropertySelection::SingleClaimOnly,
                bounds
            )
            .is_err()
        );
    }
    let foreign = super::card(Some("T"));
    assert!(
        MaterialGraySurfaceLaw::new(
            "native",
            &curve,
            binding(&foreign, &context),
            MaterialPropertySelection::SingleClaimOnly,
            [-100.0, 1200.0]
        )
        .is_err()
    );
    let missing = QueryPoint::new()
        .with("T", 400.0)
        .unwrap()
        .with("melt", 0.5)
        .unwrap();
    assert!(
        MaterialGraySurfaceLaw::new(
            "native",
            &curve,
            binding(&card, &missing),
            MaterialPropertySelection::SingleClaimOnly,
            [-100.0, 1200.0]
        )
        .is_err()
    );
    let source = MaterialGraySurfaceLaw::new(
        "native",
        &curve,
        binding(&card, &context),
        MaterialPropertySelection::SingleClaimOnly,
        [-100.0, 1200.0],
    )
    .unwrap();
    let same_content_different_object = curve.clone();
    let foreign_leg = Conductivity(&same_content_different_object);
    assert!(
        GrayAmbientBoundary::new(&foreign_leg, 0.0, 0.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&source)
            .is_err()
    );
}
