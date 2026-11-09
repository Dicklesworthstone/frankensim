//! Original liquid cards in the existing pressure/LNS modal interface.
//!
//! These regressions retain the adiabatic pressure-region approximation and
//! its explicitly authored thin-wall law. They test original-owner composition,
//! not continuum convergence or a fully resolved liquid field in that region.

use super::*;
use crate::cavity_modal::interior::thermoviscous::source_fluid_medium;
use crate::cavity_modal::interior::thermoviscous::tests::fluid::liquid;
use fs_feec::acoustics::hybrid::{HybridInterior as NativeHybrid, HybridResponse};
use fs_material::fluid_acoustics::thermoviscous::ResolvedThermoviscousLiquidState;

fn liquid_medium(source: &ResolvedThermoviscousLiquidState) -> Medium {
    let state = source.state();
    let acoustic = state.acoustic();
    Medium {
        density_kg_m3: acoustic.density_kg_m3(),
        sound_speed_m_s: acoustic.sound_speed_m_s(),
    }
}

fn resolved_liquid<'c>(
    mesh: &'c TetComplex,
    points: &[[f64; 3]],
    source: &ResolvedThermoviscousLiquidState,
) -> ThermoviscousFe<'c> {
    let (medium, thermodynamics) = source_fluid_medium(source, 101325.0, 0.002).unwrap();
    let fe = HelmholtzFe::new(mesh, points, 2, liquid_medium(source)).unwrap();
    let ports = ports(&fe, H_LENGTH, L_LENGTH);
    let walls = fe
        .boundary_faces()
        .into_iter()
        .filter(|face| !ports.iter().any(|p| p.traces.iter().any(|(f, _)| f == face)))
        .collect::<Vec<_>>();
    ThermoviscousFe::new_with_thermodynamics(
        mesh,
        points,
        medium,
        thermodynamics,
        &ThermoviscousBoundary {
            no_slip_faces: walls.clone(),
            isothermal_faces: walls,
            ports,
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap()
}

fn same_phasors(actual: &[C64], expected: &[C64], sign: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|z| z.abs()).fold(1e-30, f64::max);
    for (&actual, &expected) in actual.iter().zip(expected) {
        assert!((actual - expected.scale(sign)).abs() < 3e-9 * scale);
    }
}

fn same_fields(actual: &HybridResponse, expected: &HybridResponse, temperature_sign: f64) {
    same_phasors(&actual.outgoing, &expected.outgoing, 1.0);
    same_phasors(&actual.pressure_pa, &expected.pressure_pa, 1.0);
    same_phasors(&actual.inward_flow_m3_s, &expected.inward_flow_m3_s, 1.0);
    same_phasors(&actual.helmholtz.amplitudes, &expected.helmholtz.amplitudes, 1.0);
    same_phasors(
        &actual.thermoviscous.pressure_pa,
        &expected.thermoviscous.pressure_pa,
        1.0,
    );
    same_phasors(
        &actual.thermoviscous.velocity_m_s.iter().flatten().copied().collect::<Vec<_>>(),
        &expected.thermoviscous.velocity_m_s.iter().flatten().copied().collect::<Vec<_>>(),
        1.0,
    );
    same_phasors(
        &actual.thermoviscous.temperature_k,
        &expected.thermoviscous.temperature_k,
        temperature_sign,
    );
    let scale = expected.power.external_input_w.abs().max(expected.power.dissipated_w());
    assert!((actual.power.external_input_w - expected.power.external_input_w).abs() < 3e-9 * scale);
    assert!((actual.power.dissipated_w() - expected.power.dissipated_w()).abs() < 3e-9 * scale);
}

#[test]
fn original_liquid_eos_fields_work_and_duct_convention_match_native_hybrid() {
    let (hm, cube) = kuhn_cube(2);
    let hp = positions(&cube, false);
    let (lm, cube) = kuhn_cube(2);
    let lp = positions(&cube, true);
    let frequency = 400.0;
    let incident = [C64::new(0.002, 0.0006), C64::new(-0.0004, 0.0002)];
    let mut signed_fields = Vec::new();
    for alpha in [2e-4, -2e-4] {
        let (card, fluid) = liquid(Some(alpha), 0.004, 1.0);
        let h = pressure_model(&hm, &hp, liquid_medium(&fluid), false, None);
        let source = ThermoviscousInterior::from_material(
            resolved_liquid(&lm, &lp, &fluid),
            fluid.clone(),
            BAND,
        )
        .unwrap();
        let hybrid: HybridInterior<'_, '_, ResolvedThermoviscousLiquidState> =
            HybridInterior::new(&h, &source, &[1], &[0], &|| false).unwrap();
        assert!(core::ptr::eq(hybrid.source(), &source));
        assert_eq!(hybrid.source().material_source().identity(), fluid.identity());
        assert_eq!(
            hybrid.source().material_source().material().card_identity(),
            card.content_hash(),
        );
        let (medium, thermodynamics) = source_fluid_medium(&fluid, 101325.0, 0.002).unwrap();
        assert_eq!(hybrid.source().inner().medium(), medium);
        assert_eq!(hybrid.source().inner().thermodynamics(), thermodynamics);
        assert_eq!(thermodynamics.isobaric_expansion_per_k.to_bits(), alpha.to_bits());
        assert_eq!(hybrid.inner().external_port_count(), 2);
        assert!(hybrid.source().relaxation().is_none());

        let reference = fluid.state().acoustic().characteristic_impedance_pa_s_m()
            / (HEIGHT * WIDTH);
        let prepared = hybrid
            .prepare(frequency, &[reference; 2], &[reference; 2], &|| false)
            .unwrap();
        let actual = prepared.response(&incident, &|| false).unwrap();
        // Independently construct the original spatial owner and supply the
        // card's actual Cp. No adapter response or inferred ideal-gas Cp is used.
        let native = NativeHybrid::new(&h, source.inner(), &[1], &[0], &|| false).unwrap();
        let native_prepared = native
            .prepare(
                frequency,
                C64::from_re(fluid.state().transport().specific_heat_j_kg_k()),
                &[reference; 2],
                &[reference; 2],
                &|| false,
            )
            .unwrap();
        let expected = native_prepared.response(&incident, &|| false).unwrap();
        same_fields(&actual, &expected, 1.0);
        assert!(actual.thermoviscous.temperature_k.iter().any(|t| t.abs() > 0.0));
        assert!(actual.power.thermoviscous.viscous_w > 0.0);
        assert!(actual.power.thermoviscous.conductive_w > 0.0);
        assert_eq!(actual.power.thermoviscous.caloric_w, 0.0);
        assert_eq!(actual.power.helmholtz_bulk_w, 0.0);
        assert_eq!(actual.power.source_input_w, 0.0);
        assert!(actual.power.external_input_w > 0.0);
        assert!(
            actual.power.balance_residual_w().abs()
                < 2e-8 * actual.power.external_input_w
        );
        let interface_scale = actual.power.helmholtz_interface_w.abs()
            + actual.power.thermoviscous_interface_w.abs();
        assert!(actual.power.interface_residual_w().abs() < 2e-10 * interface_scale);

        let plus = prepared.scattering(&|| false).unwrap();
        let minus = hybrid
            .duct_scattering(frequency, &[reference; 2], &[reference; 2], &|| false)
            .unwrap();
        for (a, b) in minus.iter().zip(&plus) {
            assert_eq!(a.re.to_bits(), b.re.to_bits());
            assert_eq!(a.im.to_bits(), (-b.im).to_bits());
        }
        let block = hybrid.duct_two_port(frequency, reference, &|| false).unwrap();
        assert_eq!(block.scattering.as_slice(), minus.as_slice());
        signed_fields.push(actual);
    }
    // With no thermal source, alpha -> -alpha maps T' -> -T' in the ORIGINAL
    // coupled EOS. Pressure, velocity, port waves and positive losses are even.
    same_fields(&signed_fields[1], &signed_fields[0], -1.0);
}

#[test]
fn liquid_density_and_original_wall_transport_cannot_be_substituted() {
    let (_, fluid) = liquid(Some(-2e-4), 0.004, 1.0);
    let (hm, cube) = kuhn_cube(1);
    let hp = positions(&cube, false);
    let (lm, cube) = kuhn_cube(1);
    let lp = positions(&cube, true);
    let medium = liquid_medium(&fluid);
    let source = ThermoviscousInterior::from_material(
        resolved_liquid(&lm, &lp, &fluid),
        fluid.clone(),
        BAND,
    )
    .unwrap();
    for wrong in [
        Medium { density_kg_m3: medium.density_kg_m3 * (1.0 + 1e-12), ..medium },
        Medium { sound_speed_m_s: medium.sound_speed_m_s * (1.0 + 1e-12), ..medium },
    ] {
        let h = pressure_model(&hm, &hp, wrong, false, None);
        assert!(matches!(
            HybridInterior::new(&h, &source, &[1], &[0], &|| false),
            Err(HybridSourceError::InvalidInput(_))
        ));
    }
    let state = fluid.state();
    let transport = state.transport();
    let wall = BoundaryLayerTransport {
        dynamic_viscosity_pa_s: transport.dynamic_viscosity_pa_s(),
        thermal_conductivity_w_m_k: transport.thermal_conductivity_w_m_k(),
        cp_j_kg_k: transport.specific_heat_j_kg_k(),
        heat_capacity_ratio: state.heat_capacity_ratio(),
    };
    for wrong in [
        BoundaryLayerTransport {
            dynamic_viscosity_pa_s: wall.dynamic_viscosity_pa_s * 1.01,
            ..wall
        },
        BoundaryLayerTransport {
            thermal_conductivity_w_m_k: wall.thermal_conductivity_w_m_k * 1.01,
            ..wall
        },
        BoundaryLayerTransport {
            cp_j_kg_k: wall.cp_j_kg_k * 1.01,
            ..wall
        },
        BoundaryLayerTransport {
            heat_capacity_ratio: wall.heat_capacity_ratio * 1.01,
            ..wall
        },
    ] {
        let h = pressure_model(&hm, &hp, medium, false, Some(wrong));
        assert!(matches!(
            HybridInterior::new(&h, &source, &[1], &[0], &|| false),
            Err(HybridSourceError::InvalidInput(_))
        ));
    }
    let h = pressure_model(&hm, &hp, medium, false, Some(wall));
    let hybrid = HybridInterior::new(&h, &source, &[1], &[0], &|| false).unwrap();
    let reference = state.acoustic().characteristic_impedance_pa_s_m() / (HEIGHT * WIDTH);
    // At 1100 Hz the original declared thin-layer geometry is admitted.
    // This checks the authored approximation, not boundary-layer convergence.
    let field = hybrid
        .prepare(1100.0, &[reference; 2], &[reference; 2], &|| false)
        .unwrap()
        .response(&[C64::new(0.001, 0.0002), C64::ZERO], &|| false)
        .unwrap();
    assert!(field.power.helmholtz_wall.viscous_w > 0.0);
    assert!(field.power.balance_residual_w().abs() < 2e-8 * field.power.external_input_w);

    let (_, changed) = liquid(Some(-2e-4), 0.008, 1.0);
    let changed = ThermoviscousInterior::from_material(
        resolved_liquid(&lm, &lp, &changed),
        changed,
        BAND,
    )
    .unwrap();
    assert!(matches!(
        HybridInterior::new(&h, &changed, &[1], &[0], &|| false),
        Err(HybridSourceError::InvalidInput(_))
    ));
    let adiabatic = pressure_model(&hm, &hp, medium, false, None);
    assert!(HybridInterior::new(&adiabatic, &changed, &[1], &[0], &|| false).is_ok());
    assert!(matches!(
        hybrid.prepare(BAND[1] + 1.0, &[reference; 2], &[reference; 2], &|| false),
        Err(HybridSourceError::InvalidInput(_))
    ));
    assert!(matches!(
        hybrid.prepare(1100.0, &[reference; 2], &[reference; 2], &|| true),
        Err(HybridSourceError::Hybrid(HybridError::Cancelled))
    ));
}
