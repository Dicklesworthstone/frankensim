use super::*;
use crate::impulse::{ImpulseReport, WrenchImpulse2, WrenchImpulseReport2};

fn controls() -> MaterialFreeSurfaceWrenchImpulseControls2 {
    MaterialFreeSurfaceWrenchImpulseControls2 {
        impulse: MaterialFreeSurfaceImpulseControls2 {
            max_cells: 64,
            maximum_speed_m_s: 1.0,
            maximum_momentum_defect_kg_s: 1e-12,
            maximum_work_defect_j_m: 1e-12,
        },
        maximum_angular_momentum_defect_kg_m_s: 1e-12,
    }
}

fn native() -> WrenchImpulseReport2 {
    WrenchImpulseReport2 {
        requested: WrenchImpulse2 {
            origin: [2.0, -1.0],
            linear: [2.0, -3.0],
            angular: 4.0,
        },
        impulse: ImpulseReport {
            applied_impulse: [2.0, -3.0],
            momentum_change: [2.0, -3.0],
            mass_change: 0.0,
            work: -0.125,
            kinetic_energy_change: -0.125,
            absolute_momentum_defect: 0.0,
            absolute_work_defect: 0.0,
            affected_cells: 3,
        },
        stencil_mass: 2.5,
        stencil_center_of_mass: [1.5, 2.5],
        applied_angular_impulse: 4.0,
        angular_momentum_change: 4.0,
        absolute_linear_momentum_defect: 0.0,
        absolute_angular_momentum_defect: 0.0,
    }
}

fn requested() -> MaterialFreeSurfaceWrenchImpulse2 {
    MaterialFreeSurfaceWrenchImpulse2 {
        origin_m: [0.5, -0.25],
        linear_impulse_kg_s: [12.0, -18.0],
        angular_impulse_kg_m_s: 6.0,
    }
}

#[test]
fn si_receipt_uses_distinct_mass_linear_angular_and_work_units() {
    let r = receipt(
        native(),
        requested(),
        controls(),
        [2.0, 6.0, 18.0],
        0.25,
        1.5,
        [0.0; 2],
    )
    .unwrap();
    assert_eq!(r.stencil_mass_kg_m, 5.0);
    assert_eq!(r.stencil_center_of_mass_m, [0.375, 0.625]);
    assert_eq!(
        r.impulse.applied_impulse_kg_s,
        requested().linear_impulse_kg_s
    );
    assert_eq!(
        r.applied_angular_impulse_kg_m_s,
        requested().angular_impulse_kg_m_s
    );
    assert_eq!(r.impulse.work_j_m, -2.25);
    assert_eq!(r.impulse.kinetic_energy_change_j_m, -2.25);
}

#[test]
fn authored_si_origin_and_request_are_admitted_after_conversion() {
    let shift = [0.125, -0.25];
    let torque = cross(shift, requested().linear_impulse_kg_s);
    let mut shifted = requested();
    for a in 0..2 {
        shifted.origin_m[a] -= shift[a];
    }
    shifted.angular_impulse_kg_m_s += torque;
    let r = receipt(
        native(),
        shifted,
        controls(),
        [2.0, 6.0, 18.0],
        0.25,
        1.5,
        shift,
    )
    .unwrap();
    assert_eq!(
        r.angular_momentum_change_kg_m_s,
        shifted.angular_impulse_kg_m_s
    );
    assert_eq!(r.absolute_angular_momentum_defect_kg_m_s, 0.0);
    // Reporting the representable lattice-origin torque as though it were
    // about the user's original origin would conceal this finite discrepancy.
    assert!(
        receipt(
            native(),
            requested(),
            controls(),
            [2.0, 6.0, 18.0],
            0.25,
            1.5,
            shift
        )
        .is_err()
    );
    let mut wrong_linear = requested();
    wrong_linear.linear_impulse_kg_s[0] += 0.01;
    assert!(
        receipt(
            native(),
            wrong_linear,
            controls(),
            [2.0, 6.0, 18.0],
            0.25,
            1.5,
            [0.0; 2]
        )
        .is_err()
    );
}

#[test]
fn unrepresentable_si_mass_or_angular_momentum_refuses_the_dependent_receipt() {
    assert!(
        receipt(
            native(),
            requested(),
            controls(),
            [f64::MAX, 6.0, 18.0],
            0.25,
            1.5,
            [0.0; 2]
        )
        .is_err()
    );
    assert!(
        receipt(
            native(),
            requested(),
            controls(),
            [2.0, 6.0, 18.0],
            0.25,
            f64::MAX,
            [0.0; 2]
        )
        .is_err()
    );
}
