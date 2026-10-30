use super::*;
use crate::multicomponent::dissipation::{evolve, face as kernel};
use crate::multicomponent::{NormalCell, NormalPrimitive};

mod spatial;

fn prepared(formation: f64) -> PlanarTransport<ConstantTransport> {
    let g = MixtureGas::new(
        bank(vec![
            source("A", R, 3.5, 0.0, formation),
            source("B", R, 5.0, 0.0, -2.0 * formation),
        ]),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let flow = field(g, 2, Boundaries::Periodic, |_| (vec![0.5, 0.5], 0.0, 2.0));
    PlanarTransport::new(
        flow.gas().clone(),
        law(&flow, 0.07, 0.02, 0.09, 0.03),
        || false,
    )
    .unwrap()
}

fn input(densities: &[f64], velocity: [f64; 2], pressure: f64) -> NormalPrimitive<'_> {
    NormalPrimitive {
        partial_densities: densities,
        normal_velocity: velocity[0],
        transverse_velocity: velocity[1],
        pressure,
    }
}

fn cell(q: &[f64; 5]) -> NormalCell<'_> {
    NormalCell {
        partial_densities: &q[..2],
        normal_momentum: q[2],
        transverse_momentum: q[3],
        energy: q[4],
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn planar_binary_diffusion_full_stress_and_referenced_enthalpy_match_independent_flux() {
    let mut op = prepared(1e6);
    let left = input(&[0.8, 0.2], [0.3, -0.4], 1.5);
    let right = input(&[0.2, 0.8], [1.3, 1.6], 2.5);
    let mut ql = [0.0; 5];
    let mut qr = [0.0; 5];
    op.gas().normal_conserved_into(left, &mut ql).unwrap();
    op.gas().normal_conserved_into(right, &mut qr).unwrap();
    let before = (ql, qr, op.workspace_bytes());
    let mut flux = [0.0; 5];
    let report = op
        .face_from_primitives(
            left,
            right,
            0.5,
            [6.0, -3.0],
            MaxwellStefanSolveV1::default(),
            &mut flux,
            || false,
        )
        .unwrap();

    // Equal molar masses make the mass-frame binary law Fick's law and
    // eliminate barodiffusion even though this face has a pressure gradient.
    let j_a = -0.03 * (0.2 - 0.8) / 0.5;
    let tau_nn = (4.0 / 3.0 * 0.07 + 0.02) * 2.0 + (0.02 - 2.0 / 3.0 * 0.07) * -3.0;
    let tau_tn = 0.07 * (4.0 + 6.0);
    let h_a = 2.5 * (2.0 - 0.001) + 2.0;
    let h_b = 4.0 * (2.0 - 0.001) + 2.0;
    let expected = [
        j_a,
        -j_a,
        -tau_nn,
        -tau_tn,
        j_a * (h_a - h_b) - 0.09 * 2.0 - 0.8 * tau_nn - 0.6 * tau_tn,
    ];
    for (a, b) in flux.into_iter().zip(expected) {
        close(a, b);
    }
    assert!(report.maximum_gradient_residual_per_m < 1e-13);
    assert!(report.mass_frame_residual_kg_m2_s.abs() < 1e-14);
    let mut conservative = [0.0; 5];
    op.face(
        cell(&ql),
        cell(&qr),
        0.5,
        [6.0, -3.0],
        MaxwellStefanSolveV1::default(),
        &mut conservative,
        || false,
    )
    .unwrap();
    for (a, b) in conservative.into_iter().zip(flux) {
        close(a, b);
    }
    assert_eq!((ql, qr, op.workspace_bytes()), before);
}

#[test]
fn planar_flux_obeys_tangential_frame_and_formation_gauge_covariance() {
    let mut original = prepared(0.0);
    let mut shifted = prepared(1e6);
    let mut flux = [0.0; 5];
    original
        .face_from_primitives(
            input(&[0.8, 0.2], [0.3, -0.4], 1.5),
            input(&[0.2, 0.8], [1.3, 1.6], 2.5),
            0.5,
            [6.0, -3.0],
            MaxwellStefanSolveV1::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    let mut gauge = [0.0; 5];
    shifted
        .face_from_primitives(
            input(&[0.8, 0.2], [0.3, -0.4], 1.5),
            input(&[0.2, 0.8], [1.3, 1.6], 2.5),
            0.5,
            [6.0, -3.0],
            MaxwellStefanSolveV1::default(),
            &mut gauge,
            || false,
        )
        .unwrap();
    for (a, b) in gauge.into_iter().zip(flux) {
        close(a, b);
    }
    let boost = 2.0;
    let mut boosted = [0.0; 5];
    original
        .face_from_primitives(
            input(&[0.8, 0.2], [0.3, -0.4 + boost], 1.5),
            input(&[0.2, 0.8], [1.3, 1.6 + boost], 2.5),
            0.5,
            [6.0, -3.0],
            MaxwellStefanSolveV1::default(),
            &mut boosted,
            || false,
        )
        .unwrap();
    for i in 0..4 {
        close(boosted[i], flux[i]);
    }
    close(boosted[4], flux[4] + boost * flux[3]);
    let mut reversed = [0.0; 5];
    original
        .face_from_primitives(
            input(&[0.8, 0.2], [0.3, 0.4], 1.5),
            input(&[0.2, 0.8], [1.3, -1.6], 2.5),
            0.5,
            [-6.0, -3.0],
            MaxwellStefanSolveV1::default(),
            &mut reversed,
            || false,
        )
        .unwrap();
    for i in [0, 1, 2, 4] {
        close(reversed[i], flux[i]);
    }
    close(reversed[3], -flux[3]);
}

#[test]
fn planar_zero_transverse_limit_reuses_original_axial_flux_exactly() {
    let mut op = prepared(0.0);
    let g = op.gas().clone();
    let left = input(&[0.8, 0.2], [0.3, 0.0], 1.5);
    let right = input(&[0.2, 0.8], [1.3, 0.0], 2.5);
    let mut ql = [0.0; 5];
    let mut qr = [0.0; 5];
    let pl = g.normal_conserved_into(left, &mut ql).unwrap();
    let pr = g.normal_conserved_into(right, &mut qr).unwrap();
    let cvl = evolve::heat_capacity_density(&g, &ql[..2], pl).unwrap();
    let cvr = evolve::heat_capacity_density(&g, &qr[..2], pr).unwrap();
    let flow = field(g.clone(), 2, Boundaries::Periodic, |_| {
        (vec![0.5, 0.5], 0.0, 2.0)
    });
    let mut axial_law = law(&flow, 0.07, 0.02, 0.09, 0.03);
    let (diffusion, mass_ratio) = kernel::prepare_diffusion(&g).unwrap();
    let mut workspace = diffusion.workspace().unwrap();
    let mut stage = evolve::Stage {
        left: [0.0; MAX_COMPONENTS],
        right: [0.0; MAX_COMPONENTS],
        rate: 0.0,
        gradient_residual: 0.0,
        mass_frame_residual: 0.0,
    };
    let mut axial = [0.0; 4];
    evolve::interior_face(
        &g,
        &mut axial_law,
        &diffusion,
        &mut workspace,
        &ql[..2],
        &qr[..2],
        pl,
        pr,
        cvl,
        cvr,
        0.5,
        0.5,
        mass_ratio,
        &mut [0.0; 4],
        &mut axial,
        MaxwellStefanSolveV1::default(),
        &mut stage,
        &mut || false,
    )
    .unwrap();
    let mut planar = [0.0; 5];
    op.face_from_primitives(
        left,
        right,
        0.5,
        [0.0; 2],
        MaxwellStefanSolveV1::default(),
        &mut planar,
        || false,
    )
    .unwrap();
    assert_eq!([planar[0], planar[1], planar[2], planar[4]], axial);
    assert_eq!(planar[3], 0.0);
}

#[test]
fn planar_slip_wall_carries_full_normal_traction_without_mass_or_energy() {
    let mut op = prepared(0.0);
    let trace = input(&[1.0, 0.0], [0.3, 1.6], 2.0);
    let mut q = [0.0; 5];
    op.gas().normal_conserved_into(trace, &mut q).unwrap();
    let mut out = [f64::NAN; 5];
    let first = op
        .wall_from_primitive(trace, 2.0, -3.0, &mut out, || false)
        .unwrap();
    close(
        out[2],
        -((4.0 / 3.0 * 0.07 + 0.02) * 2.0 + (0.02 - 2.0 / 3.0 * 0.07) * -3.0),
    );
    for i in [0, 1, 3, 4] {
        assert_eq!(out[i], 0.0);
    }
    let mut conservative = [0.0; 5];
    let second = op
        .wall(cell(&q), 2.0, -3.0, &mut conservative, || false)
        .unwrap();
    for (a, b) in conservative.into_iter().zip(out) {
        close(a, b);
    }
    close(first.diffusivity_m2_s, second.diffusivity_m2_s);
    op.face_from_primitives(
        trace,
        trace,
        0.5,
        [0.0; 2],
        MaxwellStefanSolveV1::default(),
        &mut out,
        || false,
    )
    .unwrap();
    assert!(out.iter().all(|v| *v == 0.0));
}

#[test]
fn planar_faces_refuse_foreign_sources_invalid_axes_and_final_cancellation_then_replay() {
    let mut op = prepared(0.0);
    let other = prepared(1e6);
    let flow = field(other.gas().clone(), 2, Boundaries::Periodic, |_| {
        (vec![0.5, 0.5], 0.0, 2.0)
    });
    assert!(
        PlanarTransport::new(op.gas().clone(), law(&flow, 0.07, 0.02, 0.09, 0.03), || {
            false
        })
        .is_err()
    );
    let left = input(&[0.8, 0.2], [0.3, -0.4], 1.5);
    let right = input(&[0.2, 0.8], [1.3, 1.6], 2.5);
    let mut out = [0.0; 5];
    for (distance, gradient) in [(0.0, [0.0; 2]), (0.5, [f64::NAN, 0.0])] {
        assert!(
            op.face_from_primitives(
                left,
                right,
                distance,
                gradient,
                MaxwellStefanSolveV1::default(),
                &mut out,
                || false
            )
            .is_err()
        );
    }
    assert!(
        op.face_from_primitives(
            left,
            right,
            0.5,
            [0.0; 2],
            MaxwellStefanSolveV1::default(),
            &mut out[..4],
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    op.face_from_primitives(
        left,
        right,
        0.5,
        [6.0, -3.0],
        MaxwellStefanSolveV1::default(),
        &mut out,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let expected = out;
    let mut count = 0;
    assert_eq!(
        op.face_from_primitives(
            left,
            right,
            0.5,
            [6.0, -3.0],
            MaxwellStefanSolveV1::default(),
            &mut out,
            || {
                count += 1;
                count == polls
            }
        ),
        Err(DissipationError::Cancelled)
    );
    op.face_from_primitives(
        left,
        right,
        0.5,
        [6.0, -3.0],
        MaxwellStefanSolveV1::default(),
        &mut out,
        || false,
    )
    .unwrap();
    assert_eq!(out, expected);
}


#[test]
fn planar_no_slip_wall_carries_shear_heat_and_prescribed_wall_work() {
    let mut op = prepared(0.0);
    let trace = input(&[1.0, 0.0], [0.3, 1.6], 2.0);
    let mut q = [0.0; 5];
    let primitive = op.gas().normal_conserved_into(trace, &mut q).unwrap();
    let mut out = [f64::NAN; 5];
    op.no_slip_wall_from_primitive(
        trace,
        0.25,
        0.4,
        WallThermalCondition::Isothermal(Temperature::new(primitive.temperature - 0.5)),
        2.0,
        -3.0,
        &mut out,
        || false,
    )
    .unwrap();
    let cross = 0.02 - 2.0 / 3.0 * 0.07;
    let normal = (4.0 / 3.0 * 0.07 + 0.02) * 2.0 + cross * -3.0;
    let shear = 0.07 * (0.4 - 1.6) / 0.25;
    close(out[2], -normal);
    close(out[3], -shear);
    close(out[4], 0.4 * (-shear) - 0.09 * (-0.5 / 0.25));
    assert_eq!(out[0], 0.0);
    assert_eq!(out[1], 0.0);

    let mut stationary = [0.0; 5];
    op.no_slip_wall_from_primitive(
        trace,
        0.25,
        0.0,
        WallThermalCondition::Adiabatic,
        2.0,
        -3.0,
        &mut stationary,
        || false,
    )
    .unwrap();
    assert_ne!(stationary[3], 0.0);
    assert_eq!(stationary[4], 0.0);
}
