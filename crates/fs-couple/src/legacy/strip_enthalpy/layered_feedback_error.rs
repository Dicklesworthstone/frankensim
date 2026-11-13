// Typed refusals for layered thermal feedback.
use fs_material::thermovisco::source::SourceError;
use fs_material::thermovisco::thermal::FeedbackError;
use fs_solid::deposition;
use fs_solver::{NewtonError,NewtonReport,SolverRunProgress};
use fs_thermal::enthalpy::DepositionError;

#[derive(Debug)]
pub enum StripEnthalpyError<E>{
    Invalid(&'static str),
    Source(E),
    Binding(SourceError),
    CaloricAdmission(FeedbackError),
    Mechanical(deposition::Error),
    Thermal(DepositionError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    StructuralBudget,
    Cancelled,
}
impl<E> From<DepositionError> for StripEnthalpyError<E>{
    fn from(e:DepositionError)->Self{
        if e==DepositionError::Cancelled{Self::Cancelled}else{Self::Thermal(e)}
    }
}
impl<E> From<deposition::Error> for StripEnthalpyError<E>{
    fn from(e:deposition::Error)->Self{
        if e==deposition::Error::Cancelled{Self::Cancelled}else{Self::Mechanical(e)}
    }
}
impl<E> From<NewtonError> for StripEnthalpyError<E>{
    fn from(e:NewtonError)->Self{Self::Newton(e)}
}
impl<E:core::fmt::Debug> core::fmt::Display for StripEnthalpyError<E>{
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result{
        write!(f,"strip/enthalpy coupling: {self:?}")
    }
}
impl<E:std::error::Error+'static> std::error::Error for StripEnthalpyError<E>{}
