//! G0/G1/G3 numerical fixtures; synthetic coefficients, not measured materials.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ElectricMotionInput, ElectricSetup, ElectricalContact,
    Electrode, ElectrodeDrive, MovingConductorSolution, Voxel, VoxelDomain,
};

fn rig() -> (VoxelDomain, ElectricSetup) {
    let domain = VoxelDomain::graded_from_fn(
        [vec![1.0, 0.5, 2.0, 1.5, 1.0], vec![2.0], vec![3.0]],
        |p| if p[0] < 1.5 { Voxel::Solid(7) } else { Voxel::Fluid },
    ).unwrap();
    let setup = ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo: [0,0,0], hi: [1,1,1], drive: ElectrodeDrive::Voltage(0.0) },
        Electrode { lo: [4,0,0], hi: [5,1,1], drive: ElectrodeDrive::Voltage(9.0) },
    ] };
    (domain, setup)
}
fn contact() -> ElectricalContact {
    ElectricalContact { cells: [1,2], resistance_ohm_m2: 3.0, heat_fraction_to_first: 0.25 }
}
fn solve(domain: &VoxelDomain, setup: &ElectricSetup, u: &[[f64;3]], b: &[[f64;3]], contacts: &[ElectricalContact], scale: f64) -> MovingConductorSolution {
    setup.solve_with_cell_motion(domain, &[300.0;5], contacts,
        ElectricMotionInput { velocity_m_s: u, magnetic_flux_density_t: b, max_cells: 5, power_tolerance_w: 1e-8 },
        |c, _, _| Ok(Some(scale * [2.0,4.0,8.0][c-1])), &CancelGate::new()).unwrap()
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-8, "{a} versus {b}"); }

#[test]
fn g1_motional_voltage_and_contact_heat_match_a_graded_series_circuit() {
    let (domain, mut setup) = rig();
    let u = [[0.0,1.0,0.0];5];
    let b = [[0.0,0.0,2.0];5];
    // Bulk R = (2*.5 + 4*2 + 8*1.5)/6 = 3.5 ohm; contact = .5.
    // EMF(left -> right) = 2*(.5+2+1.5) = 8 V.
    for drive in [ElectrodeDrive::Voltage(0.0), ElectrodeDrive::Voltage(9.0),
                  ElectrodeDrive::Current(2.0), ElectrodeDrive::Current(-1.0)] {
        setup.electrodes[1].drive = drive;
        let result = solve(&domain, &setup, &u, &b, &[contact()], 1.0);
        let i = match drive { ElectrodeDrive::Voltage(v) => (8.0-v)/4.0, ElectrodeDrive::Current(i) => -i };
        near(result.electric.electrodes[0].current_a, i);
        near(result.electric.electrodes[1].current_a, -i);
        near(result.electric.electrodes[1].voltage_v, 8.0 - i*4.0);
        for (c, resistance) in [(1,1.0/6.0+0.125),(2,8.0/6.0+0.375),(3,2.0)] {
            near(result.electric.joule_w[c], i*i*resistance);
            near(result.lorentz_force_density_n_m3[c][1], -i*2.0/6.0);
        }
        near(result.electric.total_joule_w, i*i*4.0);
        near(result.mechanical_power_w, -8.0*i);
        near(result.power_defect_w, 0.0);
        near(result.transfer_defect_w, 0.0);
        assert_eq!(result.lorentz_force_density_n_m3[0], [0.0;3]);
        assert_eq!(result.lorentz_force_density_n_m3[4], [0.0;3]);
    }
    // Mechanical occupancy is never inferred from conductivity.
    assert_eq!(domain.voxel_at(1), Voxel::Solid(7));
    assert_eq!(domain.voxel_at(2), Voxel::Fluid);
}

#[test]
fn g3_short_circuit_braking_turns_mechanical_work_into_heat_and_scales_with_material() {
    let (domain, mut setup) = rig();
    setup.electrodes[1].drive = ElectrodeDrive::Voltage(0.0);
    for sign in [-1.0,1.0] {
        let u = [[0.0,sign,0.0];5];
        let b = [[0.0,0.0,2.0];5];
        let a = solve(&domain, &setup, &u, &b, &[], 1.0);
        let c = solve(&domain, &setup, &u, &b, &[], 2.0);
        assert!(a.mechanical_power_w < 0.0);
        near(a.electric.delivered_w, 0.0);
        near(a.mechanical_power_w + a.electric.total_joule_w, 0.0);
        near(c.electric.total_joule_w, 0.5*a.electric.total_joule_w);
        near(c.mechanical_power_w, 0.5*a.mechanical_power_w);
        for cell in 1..4 { near(c.lorentz_force_density_n_m3[cell][1], 0.5*a.lorentz_force_density_n_m3[cell][1]); }
    }
}

#[test]
fn g5_zero_magnetic_field_reuses_the_stationary_kernel_and_contact_partition() {
    let (domain, setup) = rig();
    let u = [[1.0,2.0,3.0];5];
    let b = [[0.0;3];5];
    let moving = solve(&domain, &setup, &u, &b, &[contact()], 1.0);
    let original = setup.solve_with_cell_contacts(&domain, &[300.0;5], &[contact()],
        |c,_,_| Ok(Some([2.0,4.0,8.0][c-1])), &CancelGate::new()).unwrap();
    assert_eq!(moving.electric, original);
    assert_eq!(moving.mechanical_power_w, 0.0);
    assert!(moving.lorentz_force_density_n_m3.iter().flatten().all(|f| *f == 0.0));
}

#[test]
fn g0_three_dimensional_work_transpose_handles_nonuniform_fields_and_axis_permutations() {
    // Cyclic permutation is a proper rotation: preserve cross-product handedness.
    for axis in 0..3 {
        let mut widths = [vec![2.0],vec![2.0],vec![2.0]];
        widths[axis] = vec![1.0,0.5,2.0,1.5,1.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mut hi0 = [1;3]; hi0[axis] = 1;
        let mut lo1 = [0;3]; lo1[axis] = 4;
        let mut hi1 = [1;3]; hi1[axis] = 5;
        let setup = ElectricSetup { conductors: vec![], electrodes: vec![
            Electrode { lo: [0;3], hi: hi0, drive: ElectrodeDrive::Voltage(0.0) },
            Electrode { lo: lo1, hi: hi1, drive: ElectrodeDrive::Voltage(3.0) },
        ] };
        let rotate = |v: [f64;3]| { let mut out=[0.0;3]; for j in 0..3 {out[(j+axis)%3]=v[j];} out };
        let u: Vec<_> = (0..5).map(|c| rotate([0.3,0.2*c as f64,-0.1])).collect();
        let b: Vec<_> = (0..5).map(|c| rotate([0.7,-0.4,1.0+c as f64])).collect();
        let result = solve(&domain, &setup, &u, &b, &[], 1.0);
        let mut work=0.0;
        for c in 0..5 { for j in 0..3 { work += result.lorentz_force_density_n_m3[c][j]*u[c][j]*domain.volume(c); } }
        near(work, result.mechanical_power_w);
        near(result.electric.delivered_w, result.electric.total_joule_w+work);
        near(result.transfer_defect_w, 0.0);
    }
}

#[test]
fn g0_motion_admission_and_missing_reference_refuse_without_hiding_currents() {
    let (domain, mut setup) = rig();
    let u = [[0.0,1.0,0.0];5]; let b = [[0.0,0.0,2.0];5];
    let good = ElectricMotionInput {velocity_m_s:&u, magnetic_flux_density_t:&b,max_cells:5,power_tolerance_w:1e-8};
    for input in [ElectricMotionInput {max_cells:4,..good}, ElectricMotionInput {velocity_m_s:&u[..4],..good},
                  ElectricMotionInput {power_tolerance_w:f64::NAN,..good}] {
        assert!(setup.solve_with_cell_motion(&domain,&[300.0;5],&[],input,
            |_,_,_| panic!("admission must precede material queries"),&CancelGate::new()).is_err());
    }
    setup.electrodes.pop(); // right island has no reference after cell 2 is insulated
    assert!(matches!(setup.solve_with_cell_motion(&domain,&[300.0;5],&[],good,
        |c,_,_| Ok(if c==2 {None} else {Some(2.0)}),&CancelGate::new()),
        Err(ChtError::InvalidInput {field:"electric.motion",..})));
    let gate=CancelGate::new();gate.request();
    assert_eq!(setup.solve_with_cell_motion(&domain,&[300.0;5],&[],good,
        |_,_,_| panic!("cancel before callback"),&gate),Err(ChtError::Cancelled));
}
