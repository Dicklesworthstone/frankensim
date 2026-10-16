//! Axisymmetric rigid striker properties from a supplied longitudinal profile.
//! A small rotation about an explicit grip pivot maps to m_eff=I/lever^2.
//! This does not predict finger mechanics, wood flexure or contact elasticity.
use super::{ImpactError,invalid};

/// One cross-section of a real measured or explicitly estimated stick/mallet.
#[derive(Debug,Clone,Copy)]
pub struct RadiusStation {
    /// Axial position [m], strictly increasing.
    pub position_m:f64,
    /// Nonnegative radius [m], linearly interpolated to the next station.
    pub radius_m:f64,
}
/// Mass properties and contact-coordinate inertia of the declared rigid shape.
#[derive(Debug,Clone,Copy)]
pub struct StrikerProperties {
    /// Full body mass [kg].
    pub mass_kg:f64,
    /// Axial center of mass [m].
    pub center_of_mass_m:f64,
    /// Transverse moment about the declared pivot [kg m^2].
    pub pivot_inertia_kg_m2:f64,
    /// Equivalent mass at the contact point [kg], not the full stick mass.
    pub contact_effective_mass_kg:f64,
}
impl StrikerProperties {
    /// Integrate actual radius/density, using degree-five-exact Gauss quadrature
    /// for the piecewise-linear radius. Includes each disk's transverse inertia.
    /// # Errors
    /// Invalid profile, density, pivot/contact positions or finite mass.
    pub fn from_profile(stations:&[RadiusStation],density_kg_m3:f64,pivot_m:f64,contact_m:f64)
        ->Result<Self,ImpactError> {
        if stations.len()<2 || stations.len()>4096 || !density_kg_m3.is_finite() || density_kg_m3<=0.0
            || !pivot_m.is_finite() || !contact_m.is_finite() || pivot_m==contact_m
            || stations.iter().any(|s|!s.position_m.is_finite() || !s.radius_m.is_finite() || s.radius_m<0.0)
            || stations.windows(2).any(|w|w[1].position_m<=w[0].position_m)
            || pivot_m<stations[0].position_m || pivot_m>stations[stations.len()-1].position_m
            || contact_m<stations[0].position_m || contact_m>stations[stations.len()-1].position_m {
            return Err(invalid("striker needs finite ordered radius stations, density and distinct in-profile pivot/contact"));
        }
        let(mut mass,mut moment,mut inertia)=(0.0,0.0,0.0);
        let root=(3.0_f64/5.0).sqrt();
        for w in stations.windows(2) {
            let length=w[1].position_m-w[0].position_m;
            for (point,weight) in [(-root,5.0/9.0),(0.0,8.0/9.0),(root,5.0/9.0)] {
                let t=0.5*(point+1.0);
                let x=w[0].position_m+t*length;
                let r=w[0].radius_m+t*(w[1].radius_m-w[0].radius_m);
                let dm=0.5*length*weight*density_kg_m3*core::f64::consts::PI*r*r;
                mass+=dm;moment+=dm*x;inertia+=dm*((x-pivot_m).powi(2)+0.25*r*r);
            }
        }
        let result=Self{mass_kg:mass,center_of_mass_m:moment/mass,pivot_inertia_kg_m2:inertia,
            contact_effective_mass_kg:inertia/(contact_m-pivot_m).powi(2)};
        if ![mass,inertia,result.contact_effective_mass_kg].iter().all(|v|v.is_finite() && *v>0.0)
            || !result.center_of_mass_m.is_finite() {return Err(invalid("striker inertia is not representable"));}
        Ok(result)
    }
}

/// Geometry-derived flexible shaft storage and reciprocal tip/hand ports.
pub mod flexible;
