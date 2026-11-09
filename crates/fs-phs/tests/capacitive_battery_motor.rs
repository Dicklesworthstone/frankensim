//! G1/G3/G4: finite chemical energy through four stored device capacitors to
//! two moving shafts. The independent reference integrates SI voltages,
//! currents and speeds, not the production generalized charge coordinates.

use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    CapacitiveBatteryMotorCycleStep, PreparedCapacitiveBatteryMotorCycle,
};
use fs_phs::circuit::converter::commutation::capacitance::CapacitiveCommutationParameters;
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutationInterval, CommutationParameters, DeviceGateState,
    DeviceOperatingBounds, SwitchDevice,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::switching::SwitchCycleControls;
use fs_phs::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::{PreparedStepError, SolveBudget};
use std::cell::Cell;

const CAPACITY: f64 = 80.0;
const LEAD: f64 = 0.08;
// Deliberately macroscopic synthetic capacitances keep the independent RK4
// reference resolved. The standalone analytic test separately covers nF/ns.
const COSS: [f64; 4] = [0.02, 0.03, 0.025, 0.015];
const BATTERY: BatteryParameters = BatteryParameters {
    series_resistance_ohm: 0.1,
    polarization_resistance_ohm: 0.2,
    polarization_capacitance_f: 0.5,
    bus_capacitance_f: 0.03,
};
const CONVERTER: BuckBoostParameters = BuckBoostParameters {
    inductance_h: 0.04,
    bus_capacitance_f: [0.04, 0.05],
    winding_resistance_ohm: 0.07,
    input_on_resistance_ohm: [0.031, 0.047],
    output_on_resistance_ohm: [0.053, 0.071],
};
const MOTOR: DcMotorParameters = DcMotorParameters {
    inductance_h: 0.02,
    winding_resistance_ohm: 0.2,
    rotor_inertia_kg_m2: 0.08,
    stator_inertia_kg_m2: 0.4,
    coupling_nm_per_a: 0.4,
    bearing_drag_nm_s_per_rad: 0.01,
};
const DEVICE: SwitchDevice = SwitchDevice {
    on_resistance_ohm: 0.08,
    off_resistance_ohm: 500.0,
    body_diode: BodyDiode {
        forward_voltage_v: 0.45,
        differential_resistance_ohm: 0.12,
    },
    bounds: DeviceOperatingBounds {
        maximum_reverse_voltage_v: 100.0,
        maximum_channel_current_a: 200.0,
        maximum_diode_current_a: 100.0,
    },
};
const HIGH: [bool; 4] = [true, false, true, false];

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}");
}
fn model(constant_ocv: bool) -> BatteryMotorDrive {
    let knots = if constant_ocv { [(0.0,12.0),(1.0,12.0)] } else { [(0.0,10.0),(1.0,14.0)] };
    AveragedBuckBoost::new(CONVERTER, BridgeModulation {input_high_fraction:0.2, output_high_fraction:0.8}).unwrap()
        .with_battery_source(BatterySource::new(BatteryOcvCurve::new(CAPACITY,&knots).unwrap(),BATTERY).unwrap(),LEAD).unwrap()
        .with_dc_motor(DcMotor::new(MOTOR).unwrap()).unwrap()
}
fn parameters() -> CapacitiveCommutationParameters {
    CapacitiveCommutationParameters {
        converter: CommutationParameters {
            inductance_h: CONVERTER.inductance_h,
            bus_capacitance_f: CONVERTER.bus_capacitance_f,
            winding_resistance_ohm: CONVERTER.winding_resistance_ohm,
            devices: [DEVICE;4],
        },
        output_capacitance_f: COSS,
    }
}
fn electrical() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {relative:1e-13,absolute:1e-13,max_iterations:40},
        backtracking_halvings:10,
        absolute_energy_tolerance_j:1e-10,
        relative_energy_tolerance:1e-11,
    }
}
fn controls() -> SwitchCycleControls {
    SwitchCycleControls {absolute_energy_tolerance_j:4e-9,relative_energy_tolerance:2e-11,..SwitchCycleControls::default()}
}
fn interval(on: [bool;4], duration_s:f64) -> CommutationInterval {
    CommutationInterval {gates:DeviceGateState{on},duration_s}
}
fn prepare(deck:&[CommutationInterval], subdivisions:usize) -> PreparedCapacitiveBatteryMotorCycle {
    PreparedCapacitiveBatteryMotorCycle::new_controlled(model(false),[DEVICE;4],COSS,deck,subdivisions,electrical(),controls(),||false).unwrap()
}
fn command(power:f64) -> BatteryMotorSwitchCommand {
    BatteryMotorSwitchCommand {
        shaft_torque_nm:[-0.2,0.03],
        auxiliary_load:(power>0.0).then(||RegulatedPowerLoad::try_new(power,3.0,0.5,50.0).unwrap()),
    }
}
fn initial_physical() -> [f64;11] {
    // [chemical C, Vpol, Vbattery, Iconv, Vin, Vout, Vleft, Vright, Im, wr, ws]
    [0.65*CAPACITY,0.1,12.1,0.5,11.8,6.4,10.0,4.0,0.1,4.0,-0.4]
}
fn state(native:&PreparedCapacitiveBatteryMotorCycle, y:&[f64;11]) -> [f64;11] {
    native.state(
        &[y[0],y[1]*BATTERY.polarization_capacitance_f,y[2]*BATTERY.bus_capacitance_f],
        &parameters().state(y[3],[y[4],y[5]],[y[6],y[7]]).unwrap(),
        &[y[8]*MOTOR.inductance_h,y[9]*MOTOR.rotor_inertia_kg_m2,y[10]*MOTOR.stator_inertia_kg_m2],
    ).unwrap()
}
fn physical(x:&[f64;11]) -> [f64;11] {
    let mut y=[0.0;11];
    y[0]=x[0]; y[1]=x[1]/BATTERY.polarization_capacitance_f; y[2]=x[2]/BATTERY.bus_capacitance_f;
    y[3]=x[3]/CONVERTER.inductance_h;
    for leg in 0..2 {
        let cb=CONVERTER.bus_capacitance_f[leg]; let ch=COSS[2*leg]; let cl=COSS[2*leg+1];
        let det=cb*(ch+cl)+ch*cl;
        y[4+leg]=((ch+cl)*x[4+leg]+ch*x[6+leg])/det;
        y[6+leg]=(ch*x[4+leg]+(cb+ch)*x[6+leg])/det;
    }
    y[8]=x[8]/MOTOR.inductance_h; y[9]=x[9]/MOTOR.rotor_inertia_kg_m2; y[10]=x[10]/MOTOR.stator_inertia_kg_m2;
    y
}
fn load_current(load:Option<RegulatedPowerLoad>, voltage:f64) -> f64 {
    let Some(load)=load else{return 0.0;};
    assert!(voltage.abs()<=load.maximum_voltage_v());
    let magnitude=if voltage.abs()<=load.foldback_voltage_v(){load.requested_power_w()*voltage.abs()/load.foldback_voltage_v().powi(2)}else{load.requested_power_w()/voltage.abs()};
    magnitude.min(load.maximum_current_a()).copysign(voltage)
}

// Eleven PHYSICAL coordinates, fourteen heat integrals, eight signed work
// integrals. No production device, capacitance, port, storage or solver helper
// participates in this SI reference derivative.
fn derivative(x:&[f64;33], on:[bool;4], cmd:BatteryMotorSwitchCommand) -> [f64;33] {
    let (q,pol,vb,ic,vi,vo,left,right,im,wr,ws)=(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10]);
    let cell=(10.0+4.0*q/CAPACITY-pol-vb)/BATTERY.series_resistance_ohm;
    let lead=(vb-vi)/LEAD;
    let load=load_current(cmd.auxiliary_load,vo);
    let relative=wr-ws;
    let torque=MOTOR.coupling_nm_per_a*im-MOTOR.bearing_drag_nm_s_per_rad*relative;
    let voltage=[vi-left,left,vo-right,right];
    let channel: [f64;4]=core::array::from_fn(|i|voltage[i]/if on[i]{DEVICE.on_resistance_ohm}else{DEVICE.off_resistance_ohm});
    let diode=voltage.map(|v|((-v-DEVICE.body_diode.forward_voltage_v)/DEVICE.body_diode.differential_resistance_ohm).max(0.0));
    let current: [f64;4]=core::array::from_fn(|i|channel[i]-diode[i]);
    let mut out=[0.0;33];
    out[0]=-cell;
    out[1]=(cell-pol/BATTERY.polarization_resistance_ohm)/BATTERY.polarization_capacitance_f;
    out[2]=(cell-lead)/BATTERY.bus_capacitance_f;
    out[3]=(left-right-CONVERTER.winding_resistance_ohm*ic)/CONVERTER.inductance_h;
    for leg in 0..2 {
        let rhs_bus=if leg==0{lead-current[0]}else{-im-load-current[2]};
        let rhs_node=current[2*leg]-current[2*leg+1]+if leg==0{-ic}else{ic};
        let cb=CONVERTER.bus_capacitance_f[leg]; let ch=COSS[2*leg]; let cl=COSS[2*leg+1];
        let det=cb*(ch+cl)+ch*cl;
        out[4+leg]=((ch+cl)*rhs_bus+ch*rhs_node)/det;
        out[6+leg]=(ch*rhs_bus+(cb+ch)*rhs_node)/det;
    }
    out[8]=(vo-MOTOR.winding_resistance_ohm*im-MOTOR.coupling_nm_per_a*relative)/MOTOR.inductance_h;
    out[9]=(torque+cmd.shaft_torque_nm[0])/MOTOR.rotor_inertia_kg_m2;
    out[10]=(-torque+cmd.shaft_torque_nm[1])/MOTOR.stator_inertia_kg_m2;
    out[11]=CONVERTER.winding_resistance_ohm*ic*ic;
    for i in 0..4 {out[12+i]=voltage[i]*channel[i];out[21+i]=-voltage[i]*diode[i];}
    out[16]=BATTERY.series_resistance_ohm*cell*cell;
    out[17]=pol*pol/BATTERY.polarization_resistance_ohm;
    out[18]=LEAD*lead*lead;
    out[19]=MOTOR.winding_resistance_ohm*im*im;
    out[20]=MOTOR.bearing_drag_nm_s_per_rad*relative*relative;
    out[25..].copy_from_slice(&[-vb*lead,vi*lead,-vo*(im+load),vo*im,MOTOR.coupling_nm_per_a*relative*im,cmd.shaft_torque_nm[0]*wr,cmd.shaft_torque_nm[1]*ws,vo*load]);
    out
}
fn reference(deck:&[CommutationInterval], cycles:usize, subdivisions:usize, cmd:BatteryMotorSwitchCommand) -> [f64;33] {
    let mut x=[0.0;33];x[..11].copy_from_slice(&initial_physical());
    for _ in 0..cycles {for interval in deck {
        let h=interval.duration_s/subdivisions as f64;
        for _ in 0..subdivisions {
            let a=derivative(&x,interval.gates.on,cmd);
            let b=derivative(&core::array::from_fn(|i|x[i]+0.5*h*a[i]),interval.gates.on,cmd);
            let c=derivative(&core::array::from_fn(|i|x[i]+0.5*h*b[i]),interval.gates.on,cmd);
            let d=derivative(&core::array::from_fn(|i|x[i]+h*c[i]),interval.gates.on,cmd);
            x=core::array::from_fn(|i|x[i]+h*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0);
        }
    }}
    x
}
fn converter_storage(y:&[f64;11]) -> [f64;7] {
    [0.5*CONVERTER.inductance_h*y[3].powi(2),
     0.5*CONVERTER.bus_capacitance_f[0]*y[4].powi(2),0.5*CONVERTER.bus_capacitance_f[1]*y[5].powi(2),
     0.5*COSS[0]*(y[4]-y[6]).powi(2),0.5*COSS[1]*y[6].powi(2),
     0.5*COSS[2]*(y[5]-y[7]).powi(2),0.5*COSS[3]*y[7].powi(2)]
}
fn storage(x:&[f64;11], constant_ocv:bool) -> [f64;7] {
    let y=physical(x);
    [if constant_ocv{12.0*y[0]}else{10.0*y[0]+2.0*y[0].powi(2)/CAPACITY},
     0.5*BATTERY.polarization_capacitance_f*y[1].powi(2),0.5*BATTERY.bus_capacitance_f*y[2].powi(2),
     converter_storage(&y).iter().sum(),0.5*MOTOR.inductance_h*y[8].powi(2),
     0.5*MOTOR.rotor_inertia_kg_m2*y[9].powi(2),0.5*MOTOR.stator_inertia_kg_m2*y[10].powi(2)]
}
fn works(step:&CapacitiveBatteryMotorCycleStep) -> [f64;8] {
    [step.battery_terminal_work_j,step.converter.terminal_work_j[0],step.converter.terminal_work_j[1],step.motor_electrical_work_j,step.electromagnetic_work_j,step.shaft_work_j[0],step.shaft_work_j[1],step.load_work_j]
}
fn account(before:&[f64;11], after:&[f64;11], step:&CapacitiveBatteryMotorCycleStep, constant_ocv:bool) {
    let a=storage(before,constant_ocv);let b=storage(after,constant_ocv);
    let change: [f64;7]=core::array::from_fn(|i|b[i]-a[i]);
    let ca=converter_storage(&physical(before));let cb=converter_storage(&physical(after));
    close(step.converter.inductor_energy_change_j,cb[0]-ca[0],3e-9);
    for i in 0..2 {close(step.converter.bus_capacitor_energy_change_j[i],cb[1+i]-ca[1+i],3e-9);}
    for i in 0..4 {close(step.converter.device_capacitor_energy_change_j[i],cb[3+i]-ca[3+i],3e-9);}
    let q=step.heat_joules();let w=works(step);
    for i in 0..7 {close(step.storage_change_j[i],change[i],3e-9);}
    assert!(q.iter().all(|v|v.is_finite()&&*v>=0.0));
    close(change[..3].iter().sum::<f64>()+q[5]+q[6],w[0],3e-9);
    close(w[0]+q[7]+w[1],0.0,3e-9);
    close(change[3]+q[..5].iter().sum::<f64>()+q[10..].iter().sum::<f64>(),w[1]+w[2],3e-9);
    close(w[2]+w[3]+w[7],0.0,3e-9);
    close(change[4]+q[8]+w[4],w[3],3e-9);
    close(change[5]+change[6]+q[9],w[4]+w[5]+w[6],3e-9);
    close(step.native.delta_h,change.iter().sum(),3e-9);
    close(step.native.dissipated,q.iter().sum(),3e-9);
    close(step.native.supplied,w[5]+w[6]-w[7],3e-9);
    close(change.iter().sum::<f64>()+q.iter().sum::<f64>(),w[5]+w[6]-w[7],3e-9);
    close(step.balance_residual_j(),0.0,3e-9);
    close(step.angular_momentum_residual_nms,0.0,3e-11);
}

#[test]
fn eleven_physical_states_fourteen_heats_and_eight_works_refine_to_independent_node_kcl_ode() {
    let deck=[interval(HIGH,0.0015),interval([false;4],0.0002),interval([false,true,false,true],0.0007),interval([false,true,true,false],0.0004),interval([true;4],0.0001),interval(HIGH,0.0011)];
    let cmd=command(2.0);
    let oracle=reference(&deck,6,512,cmd);let fine=reference(&deck,6,256,cmd);
    for i in 0..33 {close(oracle[i],fine[i],1e-8);}
    let mut errors=Vec::new();
    for n in [4,8,16] {
        let mut native=prepare(&deck,n);let mut x=state(&native,&initial_physical());
        let mut heat=[0.0;14];let mut work=[0.0;8];
        for _ in 0..6 {
            let mut next=[0.0;11];let step=native.step_controlled(&x,cmd,&mut next,||false).unwrap();
            account(&x,&next,&step,false);
            assert_eq!(step.integration_steps,deck.len()*n);
            for (q,value) in heat.iter_mut().zip(step.heat_joules()){*q+=value;}
            for (w,value) in work.iter_mut().zip(works(&step)){*w+=value;}
            x=next;
        }
        let actual=physical(&x);
        let state_error=actual.iter().zip(&oracle[..11]).map(|(a,b)|(a-b).abs()).fold(0.0,f64::max);
        let power_error=heat.iter().chain(&work).zip(&oracle[11..]).map(|(a,b)|(a-b).abs()).fold(0.0,f64::max);
        errors.push([state_error,power_error]);
        if n==16 {
            assert!(state_error<5e-4&&power_error<5e-4,"{errors:?}");
            assert!(x[0]<initial_physical()[0]);
            assert!(heat[..10].iter().all(|v|*v>0.0));
        }
    }
    for metric in 0..2 {for pair in errors.windows(2) {
        let ratio=pair[0][metric]/pair[1][metric];
        assert!((3.2..4.8).contains(&ratio),"second-order physical state/work refinement: {errors:?}");
    }}
}

#[test]
fn finite_node_storage_does_not_interrupt_exact_regeneration_into_chemical_charge() {
    // Actual finite OFF leakage distinguishes cell, converter and motor currents.
    // The four explicitly initialized device capacitors retain constant energy.
    let cell=-2.0;
    let pol=BATTERY.polarization_resistance_ohm*cell;
    let vb=12.0-pol-BATTERY.series_resistance_ohm*cell;
    let vi=vb-LEAD*cell;
    let left=vi-DEVICE.on_resistance_ohm*cell;
    let ic=cell-left/DEVICE.off_resistance_ohm;
    let right=left-CONVERTER.winding_resistance_ohm*ic;
    let high=-ic+right/DEVICE.off_resistance_ohm;
    let vo=right+DEVICE.on_resistance_ohm*high;
    let im=-high;
    let relative=(vo-MOTOR.winding_resistance_ohm*im)/MOTOR.coupling_nm_per_a;
    let reaction=-MOTOR.coupling_nm_per_a*im+MOTOR.bearing_drag_nm_s_per_rad*relative;
    let ws=-2.0;let wr=ws+relative;
    let y=[CAPACITY/2.0,pol,vb,ic,vi,vo,left,right,im,wr,ws];
    let deck=[interval(HIGH,1.0/256.0);3];
    let mut native=PreparedCapacitiveBatteryMotorCycle::new_controlled(model(true),[DEVICE;4],COSS,&deck,4,electrical(),controls(),||false).unwrap();
    let before=state(&native,&y);let dt=native.period_s();let mut next=[0.0;11];
    let step=native.step_controlled(&before,BatteryMotorSwitchCommand{shaft_torque_nm:[reaction,-reaction],auxiliary_load:None},&mut next,||false).unwrap();
    close(next[0],before[0]-cell*dt,2e-11);
    for i in 1..11 {close(next[i],before[i],2e-11);}
    account(&before,&next,&step,true);
    let expected=[CONVERTER.winding_resistance_ohm*ic*ic,
        DEVICE.on_resistance_ohm*cell*cell,left*left/DEVICE.off_resistance_ohm,
        DEVICE.on_resistance_ohm*high*high,right*right/DEVICE.off_resistance_ohm,
        BATTERY.series_resistance_ohm*cell*cell,pol*pol/BATTERY.polarization_resistance_ohm,
        LEAD*cell*cell,MOTOR.winding_resistance_ohm*im*im,MOTOR.bearing_drag_nm_s_per_rad*relative*relative,
        0.0,0.0,0.0,0.0];
    for (got,power) in step.heat_joules().into_iter().zip(expected){close(got,dt*power,2e-10);}
    let expected_work=[-vb*cell,vi*cell,-vo*im,vo*im,MOTOR.coupling_nm_per_a*relative*im,reaction*wr,-reaction*ws,0.0];
    for (got,power) in works(&step).into_iter().zip(expected_work){close(got,dt*power,2e-10);}
    for change in step.converter.device_capacitor_energy_change_j {close(change,0.0,2e-12);}
    close(step.storage_change_j[0],-12.0*cell*dt,2e-10);
    assert!(step.motor_electrical_work_j<0.0&&step.storage_change_j[0]>0.0);
    assert!(step.shaft_work_j.iter().all(|w|*w>0.0));
}

#[test]
fn both_dead_time_directions_keep_dynamic_node_kcl_and_each_physical_heat_at_the_same_effort() {
    let dt=1.0/65536.0;
    for sign in [1.0,-1.0] {
        let mut y=initial_physical();y[3]=2.0*sign;
        if sign>0.0 {y[6]=-0.8;y[7]=y[5]+0.8;}else{y[6]=y[4]+0.8;y[7]=-0.8;}
        let mut native=prepare(&[interval([false;4],dt)],1);
        let before=state(&native,&y);let mut next=[0.0;11];let cmd=command(0.0);
        let step=native.step_controlled(&before,cmd,&mut next,||false).unwrap();
        account(&before,&next,&step,false);
        let after=physical(&next);let mut midpoint=[0.0;33];
        for i in 0..11 {midpoint[i]=f64::midpoint(y[i],after[i]);}
        let expected=derivative(&midpoint,[false;4],cmd);
        for i in 0..11 {close(after[i]-y[i],dt*expected[i],3e-10);}
        for (got,power) in step.heat_joules().iter().zip(&expected[11..25]) {close(*got,dt*power,2e-10);}
        for (got,power) in works(&step).iter().zip(&expected[25..]) {close(*got,dt*power,2e-10);}
        let volts=[midpoint[4]-midpoint[6],midpoint[6],midpoint[5]-midpoint[7],midpoint[7]];
        for (i,v) in volts.into_iter().enumerate() {
            close(step.converter.last_devices[i].drain_source_voltage_v,v,2e-11);
            let diode=((-v-DEVICE.body_diode.forward_voltage_v)/DEVICE.body_diode.differential_resistance_ohm).max(0.0);
            close(step.converter.last_devices[i].diode_current_a,diode,2e-10);
        }
        let active=if sign>0.0{[1,2]}else{[0,3]};
        for i in 0..4 {
            if active.contains(&i){assert!(step.heat_joules()[10+i]>0.0);}else{assert_eq!(step.heat_joules()[10+i],0.0);}
        }
    }
}

#[test]
fn late_gate_support_load_soc_and_final_cancellation_preserve_all_eleven_coordinates() {
    let dt=1.0/65536.0;let deck=[interval(HIGH,dt);2];let cmd=command(0.0);
    let mut native=prepare(&deck,3);let before=state(&native,&initial_physical());
    let mut expected=[0.0;11];let polls=Cell::new(0usize);
    let good=native.step_controlled(&before,cmd,&mut expected,||{polls.set(polls.get()+1);false}).unwrap();
    let stop=polls.get()-1;polls.set(0);let sentinel=[77.0;11];let mut next=sentinel;
    assert!(matches!(native.step_controlled(&before,cmd,&mut next,||{polls.set(polls.get()+1);polls.get()>=stop}),Err(PreparedStepError::Cancelled)));
    assert_eq!(next,sentinel);
    let bad=BatteryMotorSwitchCommand{shaft_torque_nm:[0.0;2],auxiliary_load:Some(RegulatedPowerLoad::try_new(1.0,1.0,0.25,0.5).unwrap())};
    assert!(native.step_interval_commands_controlled(&before,&[cmd,bad],&mut next,||false).is_err());
    assert_eq!(next,sentinel);
    let mut depleted=[0.0;11];depleted[0]=1e-8;
    assert!(native.step_controlled(&depleted,cmd,&mut next,||false).is_err());assert_eq!(next,sentinel);
    let retry=native.step_controlled(&before,cmd,&mut next,||false).unwrap();
    assert_eq!(next.map(f64::to_bits),expected.map(f64::to_bits));assert_eq!(retry.heat_joules(),good.heat_joules());
    let mut limited=[DEVICE;4];for d in &mut limited{d.bounds.maximum_channel_current_a=100.0;}
    let mut hidden=PreparedCapacitiveBatteryMotorCycle::new_controlled(model(false),limited,COSS,&[interval(HIGH,dt),interval([true;4],dt)],1,electrical(),controls(),||false).unwrap();
    next=sentinel;
    assert!(hidden.step_controlled(&before,cmd,&mut next,||false).is_err());assert_eq!(next,sentinel);
    for n in [0,257] {
        assert!(PreparedCapacitiveBatteryMotorCycle::new_controlled(model(false),[DEVICE;4],COSS,&deck,n,electrical(),controls(),||false).is_err());
    }
    for (actual, expected) in physical(&before).into_iter().zip(initial_physical()) {
        close(actual, expected, 2e-13);
    }
}
