# Synthetic uniform shaft for a reproducible bending experiment, NOT a measured stick.
frankensim-flexible-stick-v1
# E_Pa, density_kg_m3, modal damping ratio. Rigid rotation is never damped.
material,12000000000,800,0.001
# Pin, contact, hand stations along the same axial coordinate [m].
support,0.1,0.39,0.16
# Subdivisions per profile segment, maximum frequency [Hz], total mode ceiling.
basis,8,3000,17
station,0,0.005
station,0.4,0.005
