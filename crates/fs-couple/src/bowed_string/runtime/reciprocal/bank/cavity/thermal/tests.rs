use super::*;
use crate::vibroacoustic::CavityModes;
use fs_feec::acoustics::helmholtz::Medium;
use fs_material::gas::GasSpec;

fn gas(pressure:f64)->GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,pressure).unwrap()
}
fn with_volume<R>(side:f64, gas:GasState, f:impl FnOnce(&HelmholtzFe<'_>)->R)->R {
    let p=vec![[0.0,0.0,0.0],[side,0.0,0.0],[side,side,0.0],[0.0,side,0.0],
        [0.0,0.0,side],[side,0.0,side],[side,side,side],[0.0,side,side],
        [side/2.0,side/2.0,side/2.0]];
    let faces=[[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
        [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]];
    let complex=fs_rep_mesh::TetComplex::from_tets(p.len(),faces.into_iter()
        .map(|[a,b,c]|[8,a,b,c]).collect());
    let fe=HelmholtzFe::new(&complex,&p,1,Medium {
        density_kg_m3:gas.density,sound_speed_m_s:gas.sound_speed,
    }).unwrap();
    f(&fe)
}
fn options()->ShellCavityThermalOptions {
    ShellCavityThermalOptions { maximum_decay_rate_per_s:10.0,
        maximum_thermal_modes:1,maximum_thermal_dofs:16,
        maximum_projection_terms:1_000_000,maximum_unresolved_fraction:0.376 }
}
fn prepare(fe:&HelmholtzFe<'_>, gas:GasState, side:f64)
    ->(CavityCoupling,ShellCavityThermalLoss)
{
    let volume=side.powi(3);
    let shape=vec![1.0/det::sqrt(volume);fe.ndof()];
    let basis=CavityModes { omegas:vec![0.0],lambdas:vec![1.0],
        interface:vec![Vec::new()],loss_factor:0.0,rho0:gas.density,c0:gas.sound_speed };
    let coupling=CavityCoupling::new(&basis,2,&[0.03,-0.02],&[0.0]).unwrap();
    let loss=ShellCavityThermalLoss::prepare(fe,&gas,&[shape],&coupling,&[],options(),
        &SliceOptions::default(),&CancelGate::new_clock_free()).unwrap();
    (coupling,loss)
}
fn close(a:f64,b:f64) {
    assert!((a-b).abs()<=1e-9*a.abs().max(b.abs()).max(1e-15),"{a} != {b}");
}

#[test]
fn exact_star_tet_heat_mode_and_source_gas_determine_the_relaxation() {
    let side=0.02;
    with_volume(side,gas(101325.0),|fe| {
        let (_,loss)=prepare(fe,gas(101325.0),side);
        // Analytic P1 star mode: int(phi_c)=V/4, int(phi_c²)=V/10,
        // int(|grad phi_c|²)=4V/L². No production quadrature in this oracle.
        close(loss.diffusion_modes().eigenvalues_per_m2[0],40.0/(side*side));
        close(loss.normalized_pressure_overlaps()[0][0].powi(2),5.0/8.0);
        close(loss.unresolved_fraction_bound(),3.0/8.0);
        let air=gas(101325.0);
        close(loss.decay_rates_per_s()[0],
            air.thermal_conductivity/(air.density*(air.specific_heat_cp/air.gamma))*40.0/(side*side));
        assert_eq!(loss.memory_states(),1);
    });
    let ordinary=with_volume(side,gas(101325.0),|fe|prepare(fe,gas(101325.0),side).1);
    let pressure=with_volume(side,gas(202650.0),|fe|prepare(fe,gas(202650.0),side).1);
    close(pressure.decay_rates_per_s()[0]/ordinary.decay_rates_per_s()[0],0.5);
    close(pressure.memory_rows[0][0]/ordinary.memory_rows[0][0],det::sqrt(2.0));
    let larger=with_volume(2.0*side,gas(101325.0),|fe|prepare(fe,gas(101325.0),2.0*side).1);
    close(larger.decay_rates_per_s()[0]/ordinary.decay_rates_per_s()[0],0.25);
}

#[test]
fn pressure_force_storage_and_thermal_loss_are_one_reciprocal_law() {
    with_volume(0.02,gas(101325.0),|fe| {
        let (mut coupling,loss)=prepare(fe,gas(101325.0),0.02);
        let q=[1.0e-5,-2.0e-5];
        let r=&loss.memory_rows[0];
        let rq=r.iter().zip(q).map(|(r,q)|r*q).sum::<f64>();
        let rate=loss.rates[0];
        let mut last_pressure=None;
        for time in [0.0,0.1,1.0,100.0] {
            let z=rq*(1.0-det::exp(-rate*time));
            let state=[q[0],0.0,q[1],0.0];
            let mut pressure=[0.0];
            coupling.adiabatic_pressures_into(&state,&mut pressure).unwrap();
            let adiabatic=pressure[0];
            loss.correct_pressure(&[z],&mut pressure).unwrap();
            let ratio=1.0-(1.0-1.0/gas(101325.0).gamma)*(5.0/8.0)
                *(1.0-det::exp(-rate*time));
            close(pressure[0],adiabatic*ratio);
            if let Some(previous)=last_pressure { assert!(pressure[0].abs()<=previous); }
            last_pressure=Some(pressure[0].abs());
            let mut force=[0.0;2];
            for row in &loss.equilibrium_rows {
                let extension=row.iter().zip(q).map(|(r,q)|r*q).sum::<f64>();
                for (force,r) in force.iter_mut().zip(row) { *force-=r*extension; }
            }
            for (force,r) in force.iter_mut().zip(r) { *force-=r*(rq-z); }
            for i in 0..2 { close(force[i],coupling.compression_storage()[0].areas[i]*pressure[0]); }
            let zdot=rate*(rq-z);
            let energy_rate=(z-rq)*zdot;
            close(energy_rate,-rate*(rq-z).powi(2));
            assert!(energy_rate<=0.0);
        }
        coupling.require_compression_memory();
        let mut pressure=[123.0];
        assert!(coupling.pressures_into(&[q[0],0.0,q[1],0.0],&mut pressure).is_err());
        assert_eq!(pressure,[123.0]);
        assert!(coupling.pressure_at(&[q[0],0.0,q[1],0.0],&[1.0]).is_err());
    });
}

#[test]
fn inadequate_heat_space_and_budgets_refuse_without_repair() {
    with_volume(0.02,gas(101325.0),|fe| {
        let (coupling,_)=prepare(fe,gas(101325.0),0.02);
        let shape=vec![1.0/det::sqrt(0.02_f64.powi(3));fe.ndof()];
        for choice in 0..4 {
            let mut o=options();
            match choice { 0=>o.maximum_unresolved_fraction=0.374,
                1=>o.maximum_decay_rate_per_s=1e-9,
                2=>o.maximum_thermal_dofs=1,_=>o.maximum_projection_terms=1 }
            assert!(ShellCavityThermalLoss::prepare(fe,&gas(101325.0),&[shape.clone()],
                &coupling,&[],o,&SliceOptions::default(),&CancelGate::new_clock_free()).is_err());
        }
        let gate=CancelGate::new_clock_free();gate.request();
        assert!(ShellCavityThermalLoss::prepare(fe,&gas(101325.0),&[shape],&coupling,&[],
            options(),&SliceOptions::default(),&gate).is_err());
    });
}
