//! An impulse changes original motion while retaining the stress source owner.

use super::*;
use fs_lbm::material::free_surface::MaterialFreeSurfaceImpulseControls3;

#[test]
fn source_rheology_accepts_si_impulse_without_a_hidden_collision() {
    let gate = CancelGate::new_clock_free();
    let source = card(0.003, 1000.0, 0.001);
    let binding = law(&source);
    let cells = initial(true);
    let gravity = [0.03, -0.02, 0.01];
    with_cx(&gate, |cx| {
        let mut flow = MaterialFreeSurfaceRheology3::from_cells(
            &binding, DX, DT, 0.1, domain(gravity), &cells, solve(), cx, &gate,
        ).unwrap();
        flow.step_with_gate(&gate).unwrap();
        let before = flow.clone();
        let [x, y, z] = CELL;
        let old = before.sample_at_si(x, y, z).unwrap().unwrap();
        let index = (z * 8 + y) * 8 + x;
        let mut field = vec![[0.0; 3]; 512];
        field[index] = [2.0, 0.0, 0.0];
        let controls = MaterialFreeSurfaceImpulseControls3 {
            max_cells: 512,
            maximum_speed_m_s: 0.1,
            maximum_momentum_defect_kg_m_s: 1e-16,
            maximum_work_defect_j: 1e-17,
        };
        let report = flow.apply_impulse_density(&field, controls, &gate).unwrap();
        let new = flow.sample_at_si(x, y, z).unwrap().unwrap();
        assert!((new.velocity_m_s[0] - old.velocity_m_s[0] - 2.0 / old.density_kg_m3).abs() < 2e-13);
        assert_eq!(report.affected_cells, 1);
        assert_eq!(report.mass_change_kg, 0.0);
        assert_eq!(flow.lattice().grid().steps(), before.lattice().grid().steps());
        assert_eq!(flow.ledger_mass_kg().unwrap(), before.ledger_mass_kg().unwrap());
        assert_eq!(new.density_kg_m3.to_bits(), old.density_kg_m3.to_bits());
        assert_eq!(new.fill_fraction, old.fill_fraction);
        assert!((numerator(&flow, gravity, CELL) - numerator(&before, gravity, CELL)).abs() < 1e-12);
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..8 {
                    let p = [x, y, z];
                    assert_eq!(flow.material_at(p), before.material_at(p));
                    assert_eq!(flow.rheology_at(p), before.rheology_at(p));
                }
            }
        }
        let accepted = flow.lattice().grid().state_bits();
        field[0] = [0.5, 0.0, 0.0];
        field[index] = [1e9, 0.0, 0.0];
        assert!(flow.apply_impulse_density(&field, controls, &gate).is_err());
        assert_eq!(flow.lattice().grid().state_bits(), accepted);
        assert_eq!(flow.rheology_at(CELL), before.rheology_at(CELL));
    });
}
