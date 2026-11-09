//! Original voltage/charge curves in coupled bus and switching-node storage.
//!
//! Each device retains its complete monotone piecewise-affine V(q) curve.
//! Cdiff is dq/dV, while stored energy is the original integral of V(q).
//! Neither a sampled C(V) nor half-C(V)-V-squared replaces that law.
//!
//! A grounded leg's inverse is a scalar monotone equation in its high-device
//! charge. Original knot searches locate its affine cell; no nonlinear dynamic
//! solver or curve resampling is added. The storage-owned discrete gradient is
//! the exact average-vector-field integral along the ELECTRODE charge chord.
//! Affine cells are traced with explicit work limits, and their constant
//! Hessians also give the analytic end-state Jacobian integral.
//!
//! Entire chord support is required: supported endpoints alone need not admit
//! a chord in this nonlinear electrode domain. On a chord lying on a kink,
//! directional cell selection supplies a one-sided Newton representative;
//! continuous tangents refuse genuine kinks and source support boundaries.

use super::{CommutationParameters, invalid, number, positive, quadratic, total};
use crate::circuit::capacitance::CapacitorChargeCurve;
use crate::{ConstitutiveDiscreteGradient, ContinuousStorageTangent, PhsError, Storage};
use fs_math::dd::Dd;

/// A hard arithmetic-work limit, independent of the original source knot cap.
pub const MAX_CHARGE_CELLS_PER_CHORD: usize = 65_536;

#[derive(Clone,Debug,PartialEq)]
pub struct ChargeCommutationParameters {
    pub converter:CommutationParameters,
    /// Input high, input low, output high, output low; original C/V source laws.
    pub output_charge_curves:[CapacitorChargeCurve;4],
    /// Maximum active cells per leg and per charge-chord evaluation, including
    /// the initial cell. Exhaustion refuses without changing physical state.
    pub maximum_cells_per_chord:usize,
}
impl ChargeCommutationParameters {
    pub fn validate(&self)->Result<(),PhsError>{ChargeCommutationStorage::new(self).map(|_|())}
    pub fn state(&self,current_a:f64,bus_voltage_v:[f64;2],node_voltage_v:[f64;2])->Result<[f64;5],PhsError>{
        ChargeCommutationStorage::new(self)?.state(current_a,bus_voltage_v,node_voltage_v)
    }
    pub fn stored_energy_j(&self,state:&[f64;5])->Result<f64,PhsError>{
        ChargeCommutationStorage::new(self)?.energy(state)
    }
}
fn dd(x:f64)->Dd{Dd::from_f64(x)}
fn sign(x:Dd)->i8{
    if x.hi>0.0||(x.hi==0.0&&x.lo>0.0){1}
    else if x.hi<0.0||(x.hi==0.0&&x.lo<0.0){-1}else{0}
}
fn less(a:Dd,b:Dd)->bool{sign(a-b)<0}
fn min(a:Dd,b:Dd)->Dd{if less(a,b){a}else{b}}
fn max(a:Dd,b:Dd)->Dd{if less(a,b){b}else{a}}
fn slope(curve:&CapacitorChargeCurve,i:usize)->f64{
    let k=curve.knots();(k[i+1].1-k[i].1)/(k[i+1].0-k[i].0)
}
#[derive(Clone,Copy)]
struct Evaluation { voltage:Dd,slope:f64,index:usize }
fn evaluate(curve:&CapacitorChargeCurve,q:Dd)->Result<Evaluation,PhsError>{
    number(q)?;
    let knots=curve.knots();
    if less(q,dd(knots[0].0))||less(dd(knots[knots.len()-1].0),q){
        return Err(invalid("device charge leaves the original curve support"));
    }
    let index=knots.partition_point(|&(x,_)|!less(q,dd(x)))
        .saturating_sub(1).min(knots.len()-2);
    let anchor=if sign(q)<0{index+1}else{index};
    let derivative=slope(curve,index);
    let voltage=dd(knots[anchor].1)+dd(derivative)*(q-dd(knots[anchor].0));
    number(voltage)?;
    Ok(Evaluation{voltage,slope:derivative,index})
}
fn curve_error()->PhsError{invalid("original device charge/energy curve refused")}
fn tail(curve:&CapacitorChargeCurve,q:Dd)->Result<(f64,Dd,Dd,f64),PhsError>{
    let rounded=number(q)?;
    let voltage=curve.voltage_v(rounded).map_err(|_|curve_error())?;
    let actual=evaluate(curve,q)?;
    Ok((rounded,q-dd(rounded),dd(voltage),actual.slope))
}
fn energy_curve(curve:&CapacitorChargeCurve,q:Dd)->Result<Dd,PhsError>{
    let(r,t,v,s)=tail(curve,q)?;
    let base=curve.stored_energy_j(r).map_err(|_|curve_error())?;
    Ok(dd(base)+t*(v+dd(0.5*s)*t))
}
fn change_curve(curve:&CapacitorChargeCurve,a:Dd,b:Dd)->Result<Dd,PhsError>{
    let(ar,at,av,aslope)=tail(curve,a)?;
    let(br,bt,bv,bslope)=tail(curve,b)?;
    let base=curve.energy_change_j(ar,br).map_err(|_|curve_error())?;
    Ok(dd(base)+bt*(bv+dd(0.5*bslope)*bt)-at*(av+dd(0.5*aslope)*at))
}
fn remainder_curve(curve:&CapacitorChargeCurve,a:Dd,b:Dd,m:Dd,reference_v:Dd)->Result<Dd,PhsError>{
    let(ar,at,av,aslope)=tail(curve,a)?;
    let(br,bt,bv,bslope)=tail(curve,b)?;
    let mr=number(m)?;
    let mv=dd(curve.voltage_v(mr).map_err(|_|curve_error())?);
    let base=curve.energy_change_remainder_j(ar,br,mr).map_err(|_|curve_error())?;
    Ok(dd(base)+(mv-reference_v)*(dd(br)-dd(ar))
        +bt*(bv-reference_v+dd(0.5*bslope)*bt)-at*(av-reference_v+dd(0.5*aslope)*at))
}

#[derive(Clone,Debug)]
struct Leg { bus_f:f64,high:CapacitorChargeCurve,low:CapacitorChargeCurve }
#[derive(Clone,Copy)]
struct Point { charge:[Dd;2], voltage:[Dd;2] }
#[derive(Clone,Copy)]
struct Cell { index:[usize;2], derivative:[Dd;2], hessian:[f64;3] }
struct Chord { mean:[Dd;2], end_hessian:[Dd;3] }
impl Leg {
    fn indices(&self,index:[usize;2],direction:[Dd;2])->Result<Cell,PhsError>{
        let beta=dd(1.0)/dd(self.bus_f);
        let ah=slope(&self.high,index[0]);let al=slope(&self.low,index[1]);
        let scale=number(beta)?.max(ah).max(al);
        let bn=beta/dd(scale);let hn=dd(ah)/dd(scale);let ln=dd(al)/dd(scale);
        let denominator=bn+hn+ln;
        let dh=(bn*direction[0]-ln*direction[1])/denominator;
        let dl=direction[1]+dh;
        let qbb=number(beta*(hn+ln)/denominator)?;
        let qbn=number(beta*ln/denominator)?;
        let qnn=number(dd(al)*(bn+hn)/denominator)?;
        if !positive(qbb)||!positive(qbn)||!positive(qnn){return Err(invalid("coupled charge Hessian is unrepresentable"));}
        number(dh)?;number(dl)?;
        Ok(Cell{index,derivative:[dh,dl],hessian:[qbb,qbn,qnn]})
    }
    fn equation(&self,qb:Dd,qn:Dd,qh:Dd)->Result<Dd,PhsError>{
        let h=evaluate(&self.high,qh)?;let l=evaluate(&self.low,qn+qh)?;
        let beta=dd(1.0)/dd(self.bus_f);let scale=number(beta)?.max(h.slope).max(l.slope);
        let value=h.voltage/dd(scale)+l.voltage/dd(scale)+(beta/dd(scale))*(qh-qb);
        number(value)?;Ok(value)
    }
    fn point(&self,q:[Dd;2])->Result<Point,PhsError>{
        number(q[0])?;number(q[1])?;
        let hk=self.high.knots();let lk=self.low.knots();
        let lo=max(dd(hk[0].0),dd(lk[0].0)-q[1]);
        let hi=min(dd(hk[hk.len()-1].0),dd(lk[lk.len()-1].0)-q[1]);
        if less(hi,lo){return Err(invalid("coupled electrode charge has no supported device-charge inverse"));}
        let flo=self.equation(q[0],q[1],lo)?;let fhi=self.equation(q[0],q[1],hi)?;
        if sign(flo)>0||sign(fhi)<0{return Err(invalid("electrode charge leaves the original coupled curve domain"));}
        let root=if sign(flo)==0{lo}else if sign(fhi)==0{hi}else{
            let locate=|curve:&CapacitorChargeCurve,shift:Dd|->Result<usize,PhsError>{
                let mut a=0usize;let mut b=curve.knots().len();
                while a<b{
                    let m=a+(b-a)/2;let x=dd(curve.knots()[m].0)-shift;
                    let nonpositive=if less(x,lo){true}else if less(hi,x){false}
                        else{sign(self.equation(q[0],q[1],x)?)<=0};
                    if nonpositive{a=m+1;}else{b=m;}
                }
                Ok(a.saturating_sub(1).min(curve.knots().len()-2))
            };
            let index=[locate(&self.high,Dd::ZERO)?,locate(&self.low,q[1])?];
            let lower=max(lo,max(dd(hk[index[0]].0),dd(lk[index[1]].0)-q[1]));
            let upper=min(hi,min(dd(hk[index[0]+1].0),dd(lk[index[1]+1].0)-q[1]));
            if less(upper,lower){return Err(invalid("original charge-knot searches disagree"));}
            let x=dd(0.5)*lower+dd(0.5)*upper;
            let h=evaluate(&self.high,x)?;let l=evaluate(&self.low,q[1]+x)?;
            let beta=dd(1.0)/dd(self.bus_f);
            let scale=number(beta)?.max(h.slope).max(l.slope);
            let denominator=beta/dd(scale)+dd(h.slope)/dd(scale)+dd(l.slope)/dd(scale);
            x-self.equation(q[0],q[1],x)?/denominator
        };
        let charge=[root,q[1]+root];
        let h=evaluate(&self.high,charge[0])?;let l=evaluate(&self.low,charge[1])?;
        let voltage=[(q[0]-root)/dd(self.bus_f),l.voltage];
        let error=number(voltage[0]-h.voltage-l.voltage)?;
        let magnitude=total([number(voltage[0])?.abs(),number(h.voltage)?.abs(),number(l.voltage)?.abs()]);
        if !magnitude.is_finite()||error.abs()>128.0*f64::EPSILON*magnitude{
            return Err(invalid("original coupled charge/voltage inverse is numerically unresolved"));
        }
        Ok(Point{charge,voltage})
    }
    fn directional_cell(&self,charge:[Dd;2],direction:[Dd;2])->Result<Cell,PhsError>{
        let curves=[&self.high,&self.low];
        let ev=[evaluate(curves[0],charge[0])?,evaluate(curves[1],charge[1])?];
        let mut options=[[None;2];2];
        for k in 0..2{
            let i=ev[k].index;options[k][0]=Some(i);
            if i>0&&sign(charge[k]-dd(curves[k].knots()[i].0))==0{options[k][1]=Some(i-1);}
        }
        for hi in options[0].into_iter().flatten(){for li in options[1].into_iter().flatten(){
            let cell=self.indices([hi,li],direction)?;
            let mut supported=true;
            for k in 0..2{
                let bounds=curves[k].knots();let i=cell.index[k];
                if (sign(charge[k]-dd(bounds[i].0))==0&&sign(cell.derivative[k])<0)
                    ||(sign(charge[k]-dd(bounds[i+1].0))==0&&sign(cell.derivative[k])>0){supported=false;}
            }
            if supported{return Ok(cell);}
        }}
        Err(invalid("charge chord leaves support at a device knot"))
    }
    fn chord(&self,a:[Dd;2],b:[Dd;2],budget:usize)->Result<Chord,PhsError>{
        let first=self.point(a)?;self.point(b)?;
        let direction=[b[0]-a[0],b[1]-a[1]];
        let mut charge=first.charge;
        let mut position=Dd::ZERO;
        let mut mean=[Dd::ZERO;2];let mut end_hessian=[Dd::ZERO;3];
        for _ in 0..budget{
            let cell=self.directional_cell(charge,direction)?;
            let remaining=dd(1.0)-position;
            let mut width=remaining;
            let mut hits=[None;2];
            for k in 0..2{
                let curve=if k==0{&self.high}else{&self.low};
                let d=cell.derivative[k];let ds=sign(d);
                if ds==0{continue;}
                let boundary=dd(curve.knots()[cell.index[k]+usize::from(ds>0)].0);
                let endpoint=charge[k]+d*remaining;
                if (ds>0&&!less(endpoint,boundary))||(ds<0&&!less(boundary,endpoint)){
                    let w=(boundary-charge[k])/d;
                    if sign(w)<=0{return Err(invalid("unrepresentable positive charge-cell span"));}
                    if less(w,width){width=w;}
                    hits[k]=Some((w,boundary));
                }
            }
            // The final span is defined by the authored endpoint, not by
            // a rounded reconstruction of position + (1 - position). That DD
            // sum can retain a nonzero low word on either side of one.
            let terminal=sign(width-remaining)==0;
            let next=if terminal{dd(1.0)}else{position+width};
            if sign(width)<=0||!less(position,next)||(!terminal&&!less(next,dd(1.0))){
                return Err(invalid("charge-cell traversal cannot advance representably"));
            }
            let qb=a[0]+direction[0]*position;
            let vbus=(qb-charge[0])/dd(self.bus_f);
            let vlo=evaluate(&self.low,charge[1])?.voltage;
            let dv=[(direction[0]-cell.derivative[0])/dd(self.bus_f),
                dd(slope(&self.low,cell.index[1]))*cell.derivative[1]];
            for k in 0..2{mean[k]=mean[k]+width*([vbus,vlo][k]+dd(0.5)*width*dv[k]);}
            let weight=width*(position+dd(0.5)*width);
            for k in 0..3{end_hessian[k]=end_hessian[k]+weight*dd(cell.hessian[k]);}
            position=next;
            if terminal{
                for v in mean.into_iter().chain(end_hessian){number(v)?;}
                return Ok(Chord{mean,end_hessian});
            }
            // Snap only an analytically encountered source-knot coordinate.
            // The other coordinate retains q_low-q_high=q_node on the chord.
            let node=a[1]+direction[1]*position;
            let mut high_hit=None;let mut low_hit=None;
            if let Some((w,boundary))=hits[0]{if sign(w-width)==0{high_hit=Some(boundary);}}
            if let Some((w,boundary))=hits[1]{if sign(w-width)==0{low_hit=Some(boundary);}}
            if let (Some(high),Some(low))=(high_hit,low_hit){
                if sign((low-high)-node)!=0{
                    return Err(invalid("simultaneous charge-knot intersection is numerically unresolved"));
                }
                charge=[high,low];
            }
            else if let Some(high)=high_hit{charge=[high,node+high];}
            else if let Some(low)=low_hit{charge=[low-node,low];}
            else{return Err(invalid("charge-cell traversal lost its original knot boundary"));}
        }
        Err(invalid("explicit charge-chord active-cell budget exhausted"))
    }
    fn smooth_hessian(&self,q:[Dd;2])->Result<[f64;3],PhsError>{
        let point=self.point(q)?;
        for (curve,charge)in[(&self.high,point.charge[0]),(&self.low,point.charge[1])]{
            let knots=curve.knots();let e=evaluate(curve,charge)?;
            if sign(charge-dd(knots[0].0))==0||sign(charge-dd(knots[knots.len()-1].0))==0
                ||(e.index>0&&sign(charge-dd(knots[e.index].0))==0&&slope(curve,e.index-1)!=e.slope){
                return Err(invalid("continuous coupled charge tangent needs smooth interior source support"));
            }
        }
        Ok(self.directional_cell(point.charge,[Dd::ZERO;2])?.hessian)
    }
}

/// Prepared immutable five-coordinate storage. Curve clones share original
/// source tables; gradient, AVF and tangent queries allocate no scratch.
#[derive(Clone,Debug)]
pub struct ChargeCommutationStorage {
    inductance_h:f64,
    legs:[Leg;2],
    maximum_cells:usize,
    electrode_scales:[f64;4],
}
impl ChargeCommutationStorage {
    pub fn new(parameters:&ChargeCommutationParameters)->Result<Self,PhsError>{
        parameters.converter.validate()?;
        if !(1..=MAX_CHARGE_CELLS_PER_CHORD).contains(&parameters.maximum_cells_per_chord){
            return Err(invalid("bounded positive charge-chord cell budget required"));
        }
        let curves=&parameters.output_charge_curves;
        let mut electrode_scales=[0.0_f64;4];
        for k in 0..2{
            let high=curves[2*k].knots();let low=curves[2*k+1].knots();
            let high_ends=[high[0],high[high.len()-1]];
            let low_ends=[low[0],low[low.len()-1]];
            for (qh,vh)in high_ends{for (ql,vl)in low_ends{
                let qb=number(dd(parameters.converter.bus_capacitance_f[k])*(dd(vh)+dd(vl))+dd(qh))?;
                let qn=number(dd(ql)-dd(qh))?;
                electrode_scales[k]=electrode_scales[k].max(qb.abs());
                electrode_scales[2+k]=electrode_scales[2+k].max(qn.abs());
            }}
        }
        if electrode_scales.iter().any(|&v|!positive(v)){
            return Err(invalid("coupled electrode charge support is unrepresentable"));
        }
        Ok(Self{inductance_h:parameters.converter.inductance_h,
            legs:core::array::from_fn(|k|Leg{bus_f:parameters.converter.bus_capacitance_f[k],
                high:curves[2*k].clone(),low:curves[2*k+1].clone()}),
            maximum_cells:parameters.maximum_cells_per_chord,electrode_scales})
    }
    pub(super) fn maximum_cells_per_chord(&self)->usize{self.maximum_cells}
    pub(super) fn state(&self,current:f64,bus:[f64;2],node:[f64;2])->Result<[f64;5],PhsError>{
        if !current.is_finite()||bus.iter().chain(&node).any(|v|!v.is_finite()){
            return Err(invalid("finite authored nonlinear capacitor initialization required"));
        }
        let mut x=[0.0;5];x[0]=number(dd(self.inductance_h)*dd(current))?;
        for k in 0..2{
            let high=self.legs[k].high.charge_c(number(dd(bus[k])-dd(node[k]))?).map_err(|_|curve_error())?;
            let low=self.legs[k].low.charge_c(node[k]).map_err(|_|curve_error())?;
            x[1+k]=number(dd(self.legs[k].bus_f)*dd(bus[k])+dd(high))?;
            x[3+k]=number(dd(low)-dd(high))?;
        }
        self.energy(&x)?;Ok(x)
    }
    pub(super) fn energy(&self,x:&[f64;5])->Result<f64,PhsError>{
        if x.iter().any(|v|!v.is_finite()){return Err(invalid("finite original nonlinear charge state required"));}
        let mut energy=dd(0.5)*dd(x[0])*dd(x[0])/dd(self.inductance_h);
        for k in 0..2{
            let p=self.legs[k].point([dd(x[1+k]),dd(x[3+k])])?;
            energy=energy+dd(0.5*self.legs[k].bus_f)*p.voltage[0]*p.voltage[0]
                +energy_curve(&self.legs[k].high,p.charge[0])?+energy_curve(&self.legs[k].low,p.charge[1])?;
        }
        let h=number(energy)?;
        if h<0.0||(h==0.0&&x.iter().any(|v|*v!=0.0)){return Err(invalid("nonlinear capacitor stored energy is unrepresentable"));}
        Ok(h)
    }
    pub(super) fn parts(&self,a:&[f64;5],b:&[f64;5])->Result<[f64;7],PhsError>{
        let mut parts=[0.0;7];parts[0]=quadratic(a[0],b[0],self.inductance_h);
        for k in 0..2{
            let leg=&self.legs[k];
            let pa=leg.point([dd(a[1+k]),dd(a[3+k])])?;
            let pb=leg.point([dd(b[1+k]),dd(b[3+k])])?;
            parts[1+k]=number(dd(leg.bus_f)*(pb.voltage[0]-pa.voltage[0])
                *(dd(0.5)*pa.voltage[0]+dd(0.5)*pb.voltage[0]))?;
            parts[3+2*k]=number(change_curve(&leg.high,pa.charge[0],pb.charge[0])?)?;
            parts[4+2*k]=number(change_curve(&leg.low,pa.charge[1],pb.charge[1])?)?;
        }
        Ok(parts)
    }
    fn efforts(&self,x:&[f64;5])->Result<[f64;5],PhsError>{
        let mut e=[0.0;5];e[0]=number(dd(x[0])/dd(self.inductance_h))?;
        for k in 0..2{
            let p=self.legs[k].point([dd(x[1+k]),dd(x[3+k])])?;
            e[1+k]=number(p.voltage[0])?;e[3+k]=number(p.voltage[1])?;
        }
        Ok(e)
    }
    fn chord(&self,a:&[f64;5],b:&[f64;5])->Result<([f64;5],[[Dd;3];2]),PhsError>{
        let mut e=[0.0;5];let mut jacobian=[[Dd::ZERO;3];2];
        e[0]=number((dd(0.5)*dd(a[0])+dd(0.5)*dd(b[0]))/dd(self.inductance_h))?;
        for k in 0..2{
            let c=self.legs[k].chord([dd(a[1+k]),dd(a[3+k])],[dd(b[1+k]),dd(b[3+k])],self.maximum_cells)?;
            e[1+k]=number(c.mean[0])?;e[3+k]=number(c.mean[1])?;jacobian[k]=c.end_hessian;
        }
        Ok((e,jacobian))
    }
    fn remainder(&self,a:&[f64;5],b:&[f64;5],m:&[f64;5])->Result<f64,PhsError>{
        let mut r=(dd(b[0])-dd(a[0]))*dd(0.5)/dd(self.inductance_h)
            *((dd(a[0])-dd(m[0]))+(dd(b[0])-dd(m[0])));
        for k in 0..2{
            let leg=&self.legs[k];
            let pa=leg.point([dd(a[1+k]),dd(a[3+k])])?;
            let pb=leg.point([dd(b[1+k]),dd(b[3+k])])?;
            let pm=leg.point([dd(m[1+k]),dd(m[3+k])])?;
            let da=pa.voltage[0]-pm.voltage[0];let db=pb.voltage[0]-pm.voltage[0];
            r=r+dd(0.5*leg.bus_f)*(db-da)*(db+da)
                +remainder_curve(&leg.high,pa.charge[0],pb.charge[0],pm.charge[0],pm.voltage[0]-pm.voltage[1])?
                +remainder_curve(&leg.low,pa.charge[1],pb.charge[1],pm.charge[1],pm.voltage[1])?;
        }
        number(r)
    }
}
impl Storage for ChargeCommutationStorage{
    fn shared_clone(&self)->Option<Box<dyn Storage>>{Some(Box::new(self.clone()))}
    fn continuous_tangent(&self)->Option<&dyn ContinuousStorageTangent>{Some(self)}
    fn discrete_gradient_rule(&self)->Option<&dyn ConstitutiveDiscreteGradient>{Some(self)}
    fn hamiltonian(&self,x:&[f64])->f64{x.try_into().ok().and_then(|x|self.energy(x).ok()).unwrap_or(f64::NAN)}
    fn gradient(&self,x:&[f64],out:&mut[f64]){
        let e=x.try_into().ok().and_then(|x|self.efforts(x).ok());
        if out.len()==5{if let Some(e)=e{out.copy_from_slice(&e);return;}}out.fill(f64::NAN);
    }
    fn midpoint_energy_remainder(&self,a:&[f64],b:&[f64],m:&[f64])->Option<f64>{
        let result=a.try_into().ok().zip(b.try_into().ok()).zip(m.try_into().ok())
            .and_then(|((a,b),m)|self.remainder(a,b,m).ok());
        Some(result.unwrap_or(f64::NAN))
    }
    fn coordinate_scale(&self,coordinate:usize)->Option<f64>{
        coordinate.checked_sub(1).and_then(|k|self.electrode_scales.get(k).copied())
    }
    fn quadratic_part(&self)->Option<crate::linear_exact::QuadraticPart>{
        let mut matrix=vec![0.0;25];matrix[0]=1.0/self.inductance_h;
        Some(crate::linear_exact::QuadraticPart{matrix,exact:false})
    }
}
impl ConstitutiveDiscreteGradient for ChargeCommutationStorage{
    fn evaluate(&self,a:&[f64],b:&[f64],_:&[f64],out:&mut[f64]){
        let mean=a.try_into().ok().zip(b.try_into().ok()).and_then(|(a,b)|self.chord(a,b).ok()).map(|c|c.0);
        if out.len()==5{if let Some(e)=mean{out.copy_from_slice(&e);return;}}out.fill(f64::NAN);
    }
    fn end_jacobian_action(&self,a:&[f64],b:&[f64],direction:&[f64],out:&mut[f64])->bool{
        if direction.len()!=5||out.len()!=5||direction.iter().any(|v|!v.is_finite()){return false;}
        let Some((_,jac))=a.try_into().ok().zip(b.try_into().ok()).and_then(|(a,b)|self.chord(a,b).ok())else{return false;};
        out.fill(0.0);out[0]=0.5*direction[0]/self.inductance_h;
        for k in 0..2{
            out[1+k]=(jac[k][0]*dd(direction[1+k])+jac[k][1]*dd(direction[3+k])).to_f64();
            out[3+k]=(jac[k][1]*dd(direction[1+k])+jac[k][2]*dd(direction[3+k])).to_f64();
        }
        out.iter().all(|v|v.is_finite())
    }
}
impl ContinuousStorageTangent for ChargeCommutationStorage{
    fn apply(&self,state:&[f64],direction:&[f64],out:&mut[f64])->bool{
        if state.len()!=5||direction.len()!=5||out.len()!=5||state.iter().chain(direction).any(|v|!v.is_finite()){return false;}
        out.fill(0.0);out[0]=direction[0]/self.inductance_h;
        for k in 0..2{
            let Ok(h)=self.legs[k].smooth_hessian([dd(state[1+k]),dd(state[3+k])])else{return false;};
            out[1+k]=(dd(h[0])*dd(direction[1+k])+dd(h[1])*dd(direction[3+k])).to_f64();
            out[3+k]=(dd(h[1])*dd(direction[1+k])+dd(h[2])*dd(direction[3+k])).to_f64();
        }
        out.iter().all(|v|v.is_finite())
    }
}
