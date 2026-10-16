use fs_vfit::impedance::{ImpedanceFrame, SeriesImpedance, SeriesImpedanceSpec};
use fs_vfit::relaxation::{RelaxationImpedance, RelaxationImpedanceSpec, RelaxationTerm, MAX_RELAXATION_BRANCHES};

fn base() -> SeriesImpedanceSpec {
    SeriesImpedanceSpec { resistance_pa_s_m3: 0.4, inertance_pa_s2_m3: 0.3, compliance_m3_pa: Some(0.7) }
}
fn terms() -> [RelaxationTerm; 2] {
    [RelaxationTerm { resistance_pa_s_m3: 0.6, rate_per_s: 2.0 },
     RelaxationTerm { resistance_pa_s_m3: 0.3, rate_per_s: 40.0 }]
}
fn load() -> RelaxationImpedance {
    RelaxationImpedance::new(RelaxationImpedanceSpec::new(base(), &terms()).unwrap(), 1.0, 0.02).unwrap()
}
fn bits(f: ImpedanceFrame) -> Vec<u64> {
    [f.state.inertive_flow_m3_s, f.state.compliance_pressure_pa, f.reflected_pressure_pa,
     f.pressure_pa, f.flow_m3_s, f.stored_energy_j, f.storage_change_j,
     f.dissipated_energy_j, f.supplied_work_j].map(f64::to_bits).to_vec()
}

#[test]
fn empty_relaxation_reproduces_the_existing_series_load_exactly() {
    for spec in [base(), SeriesImpedanceSpec { resistance_pa_s_m3: 0.0,
        inertance_pa_s2_m3: 0.0, compliance_m3_pa: None }] {
        let mut a = SeriesImpedance::new(spec, 1.0, 0.02).unwrap();
        let mut b = RelaxationImpedance::new(RelaxationImpedanceSpec::new(spec, &[]).unwrap(), 1.0, 0.02).unwrap();
        for n in 0..512 {
            let x = f64::from(n).sin();
            assert_eq!(bits(a.step(x).unwrap()), bits(b.step(x).unwrap().port));
        }
    }
}

#[test]
fn branch_states_satisfy_independent_constitutive_and_energy_equations() {
    let mut model = load();
    let mut old_base = fs_vfit::impedance::ImpedanceState::default();
    let mut released = false;
    let mut omitted_storage_defect = 0.0_f64;
    for n in 0..2000 {
        let old = model.branch_flows_m3_s().to_vec();
        let before = model.stored_energy_j();
        let f = model.step(if n < 300 { f64::from(n % 13) - 6.0 } else { 0.0 }).unwrap();
        let p = f.port;
        let s = p.state;
        let mut energy = 0.5 * 0.3 * s.inertive_flow_m3_s.powi(2)
            + 0.5 * 0.7 * s.compliance_pressure_pa.powi(2);
        let mut pressure = 0.4 * p.flow_m3_s
            + 0.3 * (s.inertive_flow_m3_s - old_base.inertive_flow_m3_s) / 0.02
            + 0.5 * (s.compliance_pressure_pa + old_base.compliance_pressure_pa);
        let mut old_branch_energy = 0.0;
        for (i, term) in terms().iter().enumerate() {
            let q = f.branch_flows_m3_s()[i];
            let mid = 0.5 * (old[i] + q);
            let rhs = term.rate_per_s * (p.flow_m3_s - mid);
            assert!(((q - old[i]) / 0.02 - rhs).abs() <= 1e-11 * (1.0 + rhs.abs()));
            pressure += term.resistance_pa_s_m3 * (p.flow_m3_s - mid);
            energy += 0.5 * term.resistance_pa_s_m3 / term.rate_per_s * q * q;
            old_branch_energy += 0.5 * term.resistance_pa_s_m3 / term.rate_per_s * old[i] * old[i];
        }
        let scale = (before + energy + p.supplied_work_j.abs() + p.dissipated_energy_j).max(f64::MIN_POSITIVE);
        assert!((energy - p.stored_energy_j).abs() < 1e-12 * scale);
        assert!(p.balance_residual_j().abs() < 1e-12 * scale);
        assert!((pressure - p.pressure_pa).abs() < 1e-11 * (1.0 + pressure.abs()));
        assert!(p.dissipated_energy_j >= 0.0);
        omitted_storage_defect = omitted_storage_defect.max(
            (p.balance_residual_j() - (f.relaxation_stored_energy_j - old_branch_energy)).abs());
        if n >= 300 {
            assert!(p.stored_energy_j <= before + 1e-12 * scale);
            released |= p.supplied_work_j < 0.0;
        }
        old_base = s;
    }
    assert!(released, "internal memory must be able to return energy");
    assert!(omitted_storage_defect > 1e-4, "omitting branch energy must fail");
}

#[test]
fn reflection_matches_independent_complex_impedance_including_bilinear_warp() {
    // Independent real/imaginary Z(jw), not the runtime's elimination constants.
    for k in [3, 13, 53, 151, 401, 701] {
        let theta = 2.0 * core::f64::consts::PI * f64::from(k) / 2048.0;
        let mut model = load();
        let (mut re, mut im) = (0.0, 0.0);
        for n in 0..12288 {
            let phase = theta * f64::from(n);
            let wave = model.step(phase.cos()).unwrap().port.reflected_pressure_pa;
            if n >= 10240 { re += wave * phase.cos() / 1024.0; im -= wave * phase.sin() / 1024.0; }
        }
        let w = 100.0 * (theta * 0.5).tan();
        let (mut zr, mut zi) = (0.4, 0.3 * w - 1.0 / (0.7 * w));
        for t in terms() {
            let d = w * w + t.rate_per_s * t.rate_per_s;
            zr += t.resistance_pa_s_m3 * w * w / d;
            zi += t.resistance_pa_s_m3 * w * t.rate_per_s / d;
        }
        let den = (zr + 1.0).powi(2) + zi * zi;
        let expected_re = (zr * zr - 1.0 + zi * zi) / den;
        let expected_im = 2.0 * zi / den;
        assert!((re - expected_re).hypot(im - expected_im) < 1e-10);
    }
}

#[test]
fn preview_and_refusal_preserve_every_internal_history() {
    let mut a = load();
    a.step(20.0).unwrap();
    let mut b = a;
    let before = a.branch_flows_m3_s().to_vec();
    let energy = a.stored_energy_j();
    let expected = a.preview_step(7.0).unwrap();
    for x in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(a.step(x).is_err());
        assert_eq!(a.branch_flows_m3_s(), before);
        assert_eq!(a.stored_energy_j().to_bits(), energy.to_bits());
    }
    assert_eq!(a.step(7.0).unwrap(), expected);
    assert_eq!(expected, b.step(7.0).unwrap());
}

#[test]
fn branch_admission_is_bounded_and_never_repairs_negative_residues() {
    let t = terms()[0];
    assert!(RelaxationImpedanceSpec::new(base(), &[t; MAX_RELAXATION_BRANCHES]).is_ok());
    assert!(RelaxationImpedanceSpec::new(base(), &[t; MAX_RELAXATION_BRANCHES + 1]).is_err());
    for invalid in [-1.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(RelaxationImpedanceSpec::new(base(), &[RelaxationTerm { resistance_pa_s_m3: invalid, ..t }]).is_err());
        assert!(RelaxationImpedanceSpec::new(base(), &[RelaxationTerm { rate_per_s: invalid, ..t }]).is_err());
    }
    let huge = RelaxationImpedanceSpec::new(base(), &[RelaxationTerm {
        resistance_pa_s_m3: 1e300, rate_per_s: 1e-300 }]).unwrap();
    assert!(RelaxationImpedance::new(huge, 1.0, 0.02).is_err());
}
