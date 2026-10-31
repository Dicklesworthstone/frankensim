//! Native heat -> original material sources -> compatible MITC deformation.
//! Authored numerical cards and independent analytical oracles, not measurements.

mod laminate_thermal_common;
use fs_couple::laminate_thermal::{
    LaminateShellError, LaminateShellLimits, LaminateShellSnapshot, LayeredThermalShell,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::thermal::transport::stack::transient::implicit::{
    LayeredThermalSolveControls, LayeredThermalWorkspace,
};
use laminate_thermal_common::*;
use std::sync::Arc;

// Independent scalar CLT: nu=0 decouples the two identical normal directions.
// Integrate each original cell's constant E and free strain analytically.
fn clt(thicknesses: &[f64], temperatures: &[f64]) -> (f64, f64, f64) {
    let (mut a, mut b, mut d, mut n, mut m) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut z = -0.5 * thicknesses.iter().sum::<f64>();
    for (&h, &t) in thicknesses.iter().zip(temperatures) {
        let top = z + h;
        let e = young(t);
        let moment = (top * top - z * z) / 2.0;
        a += e * h;
        b += e * moment;
        d += e * (top.powi(3) - z.powi(3)) / 3.0;
        n += e * free(t) * h;
        m += e * free(t) * moment;
        z = top;
    }
    let determinant = a * d - b * b;
    let eps = (d * n - b * m) / determinant;
    let kappa = (a * m - b * n) / determinant;
    let mut energy = 0.0;
    z = -0.5 * thicknesses.iter().sum::<f64>();
    for (&h, &t) in thicknesses.iter().zip(temperatures) {
        let top = z + h;
        let r = eps - free(t);
        // Two normal directions cancel the 1/2 in each quadratic energy.
        energy += young(t)
            * (r * r * h
                + r * kappa * (top * top - z * z)
                + kappa * kappa * (top.powi(3) - z.powi(3)) / 3.0);
        z = top;
    }
    (eps, kappa, AREA * energy)
}
fn assert_motion(state: &LaminateShellSnapshot, thicknesses: &[f64], temperatures: &[f64]) {
    let (eps, kappa, energy) = clt(thicknesses, temperatures);
    let mesh = discretization().mesh;
    for (node, p) in mesh.nodes.iter().enumerate() {
        let expected = [
            eps * p[0],
            eps * p[1],
            -0.5 * kappa * (p[0] * p[0] + p[1] * p[1]),
            -kappa * p[1],
            kappa * p[0],
            0.0,
        ];
        for (i, want) in expected.into_iter().enumerate() {
            close(
                state.response().displacement_6n[6 * node + i],
                want,
                2e-8 * want.abs().max(0.001),
            );
        }
    }
    close(
        state.response().strain_energy_j,
        energy,
        2e-8 * energy.abs().max(1e-8),
    );
}

#[test]
fn g1_implicit_heat_exchange_drives_original_temperature_dependent_bimetal() {
    let cards = [card("lower"), card("upper")];
    let body = prepare_body(&cards, &[1, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[300.0, 340.0]);
    let old_cells = initial.cells().clone();
    let old_links = initial.links().clone();
    let owner = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let dt = 0.5;
    let step = owner
        .advance_to(
            &initial,
            dt,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    // Two half-layer distances in series: G=A*k/(h/2), C=rho*A*(h/2)*cp.
    let g = AREA * K / (H / 2.0);
    let c = RHO * AREA * (H / 2.0) * CP;
    let difference = 40.0 / (1.0 + 2.0 * g * dt / c);
    let expected = [320.0 - difference / 2.0, 320.0 + difference / 2.0];
    for (i, &t) in expected.iter().enumerate() {
        close(step.shell.temperatures_k()[i], t, 5e-10);
        let layer = &step.shell.material().layers()[i];
        close(layer.elasticity().state_coordinate("T").unwrap(), t, 5e-10);
        close(
            layer
                .elasticity()
                .property("young_modulus")
                .unwrap()
                .value_si(),
            young(t),
            1e-5,
        );
        close(
            step.shell.material().prepared().profiles()[i].faces()[0][0],
            free(t),
            1e-14,
        );
        assert_eq!(layer.elasticity().card_identity(), cards[i].content_hash());
    }
    assert_motion(&step.shell, &[H / 2.0; 2], &expected);
    assert!(step.shell.response().strain_energy_j > 0.0);
    close(step.thermal.mass_kg(), RHO * H * AREA, 1e-12);
    let section = step.shell.material().prepared().section();
    let reference = owner.source().reference_section();
    assert_eq!(section.i0_kg_m2().to_bits(), reference.i0_kg_m2().to_bits());
    assert_eq!(section.i1_kg_m().to_bits(), reference.i1_kg_m().to_bits());
    assert_eq!(section.i2_kg().to_bits(), reference.i2_kg().to_bits());
    // Heat conservation remains the caloric owner's separate account; stored
    // mechanical energy is measured above and is not inserted into this audit.
    close(step.thermal.energy_j(), initial.energy_j(), 1e-8);
    close(step.heat_report.energy.stored_energy_change_j, 0.0, 1e-8);
    close(step.heat_report.energy.external_heat_j, 0.0, 0.0);
    close(step.heat_report.energy.residual_j, 0.0, 1e-8);
    assert_eq!(step.heat_report.energy.substeps, 1);
    assert_eq!(step.shell.time_s(), dt);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &old_cells);
    assert_eq!(initial.links(), &old_links);
}

#[test]
fn g1_physical_layers_retain_every_unequal_thermal_cell_without_averaging() {
    let cards = [card("two-cell material"), card("three-cell material")];
    let body = prepare_body(&cards, &[2, 3], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let temperatures = [300.0, 340.0, 305.0, 325.0, 350.0];
    let initial = initial(&bound, &temperatures);
    let owner = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let state = owner.resolve(&initial, &[0.0; 24], || false).unwrap();
    let thicknesses = [H / 4.0, H / 4.0, H / 6.0, H / 6.0, H / 6.0];
    assert_eq!(owner.source().layers().len(), 5);
    for (i, source) in owner.source().layers().iter().enumerate() {
        let physical = usize::from(i >= 2);
        assert!(Arc::ptr_eq(&source.input().card, &cards[physical]));
        close(source.input().thickness_m, thicknesses[i], 1e-18);
        close(state.temperatures_k()[i], temperatures[i], 5e-10);
        close(
            state.material().prepared().profiles()[i].faces()[0][0],
            free(temperatures[i]),
            1e-14,
        );
    }
    assert_motion(&state, &thicknesses, &temperatures);
    let exact = clt(&thicknesses, &temperatures);
    let averaged = clt(&[H / 2.0; 2], &[320.0, (305.0 + 325.0 + 350.0) / 3.0]);
    assert!(
        (exact.1 - averaged.1).abs() > 1e-3,
        "fixture must distinguish actual sublayers from temperature averages"
    );
    close(
        AREA * state.material().prepared().section().i0_kg_m2(),
        initial.mass_kg(),
        1e-12,
    );
}

#[test]
fn g1_uniform_free_expansion_has_near_zero_elastic_strain_bound_and_storage() {
    let original = card("uniform");
    let cards = [Arc::clone(&original), original];
    let body = prepare_body(&cards, &[2, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[330.0; 3]);
    let owner = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let state = owner.resolve(&initial, &[0.0; 24], || false).unwrap();
    assert_motion(&state, &[H / 4.0, H / 4.0, H / 2.0], &[330.0; 3]);
    assert!(state.reduced_elastic_strain_upper_bound() < 1e-10);
    assert!(state.response().strain_energy_j < 1e-17);
    close(state.response().displacement_6n[6], free(330.0), 1e-11);
    // The card's density falls with temperature, but the original inventory
    // remains exactly rho(300 K) times reference thickness and patch area.
    close(
        AREA * state.material().prepared().section().i0_kg_m2(),
        RHO * H * AREA,
        1e-12,
    );
}

#[test]
fn g0_binding_envelope_and_late_cancellation_leave_original_timeline_unchanged() {
    let cards = [card("binding lower"), card("binding upper")];
    let body = prepare_body(&cards, &[1, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[300.0, 340.0]);
    let old_cells = initial.cells().clone();
    let old_links = initial.links().clone();
    let mut wrong_card = mechanical(&cards);
    wrong_card[0].card = card("different immutable source");
    assert!(matches!(
        LayeredThermalShell::try_new(
            &bound,
            wrong_card,
            ReferenceSurface::Midsurface,
            discretization(),
            limits(),
            || false
        ),
        Err(LaminateShellError::Invalid(_))
    ));
    let mut wrong_thickness = mechanical(&cards);
    wrong_thickness[0].thickness_m *= 1.1;
    assert!(matches!(
        LayeredThermalShell::try_new(
            &bound,
            wrong_thickness,
            ReferenceSurface::Midsurface,
            discretization(),
            limits(),
            || false
        ),
        Err(LaminateShellError::Invalid(_))
    ));
    let wrong_mass_body = prepare_body(&cards, &[1, 1], 1.1);
    let wrong_mass_laws = wrong_mass_body.thermal_laws(false).unwrap();
    let wrong_mass_bound = wrong_mass_body.bind(&wrong_mass_laws, ids()).unwrap();
    assert!(matches!(
        LayeredThermalShell::try_new(
            &wrong_mass_bound,
            mechanical(&cards),
            ReferenceSurface::Midsurface,
            discretization(),
            limits(),
            || false
        ),
        Err(LaminateShellError::Invalid(
            "thermal and mechanical reference cell mass differ"
        ))
    ));
    let restrictive = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        LaminateShellLimits {
            maximum_translation_m: 1e-8,
            ..limits()
        },
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    assert!(matches!(
        restrictive.advance_to(
            &initial,
            0.5,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false
        ),
        Err(LaminateShellError::Envelope {
            quantity: "translation",
            ..
        })
    ));
    let owner = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let mut checkpoints = 0;
    let accepted = owner
        .advance_to(
            &initial,
            0.5,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(accepted.thermal.time_s(), 0.5);
    let mut retry = 0;
    // Cancel at the final checkpoint, after both physical candidates were
    // computed: neither candidate may alter the input owner's accepted state.
    assert!(matches!(
        owner.advance_to(
            &initial,
            0.5,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || {
                retry += 1;
                retry == checkpoints
            }
        ),
        Err(LaminateShellError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &old_cells);
    assert_eq!(initial.links(), &old_links);
}
