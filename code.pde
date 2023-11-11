TITLE
'T6 Part 1'
SELECT
errlim=1e-4
ngrid=7
!spectral_colors
COORDINATES
cartesian3
VARIABLES
u !Displacement in x
v !Displacement in y
w !Displacement in z
DEFINITIONS
!mag = 10
mag = .3*globalmax(magnitude(x,y,z))/globalmax(magnitude(u,v,w))
Lx=2
Ly=.3
Lz=.2
nu=0.3
E=600e6
G=E/(2*(1+nu))
Az = Lx*Ly
p=3e3*x
Mz = 5e3
C11 =E*(1-nu)/(1+nu)/(1-2*nu)
C22 = C11
C33 = C11
C12 = E*nu/(1+nu)/(1-2*nu)
C13 = C12
C21 = C12
C23 = C12
C31 = C12
C32 = C12
!! Strain
!Axial Strain
ex=dx(u)
ey=dy(v)
ez=dz(w)
!Engineering Shear Strain
gxy=(dx(v)+dy(u))
gyz=(dy(w)+dz(v))
gxz=(dz(u)+dx(w))
!!Stress via Hooke's law
!Axial Stress
sx = C11*ex+C12*ey+C13*ez
sy = C21*ex+C22*ey+C23*ez
sz = C31*ex+C32*ey+C33*ez
!Shear stress
sxy=G*gxy
sxz=G*gxz
syz=G*gyz
Ax = Ly*Lz
Ay = Lx*Lz
EQUATIONS
!FNet = 0
u: dx(sx)+dy(sxy)+dz(sxz)=0
v: dx(sxy)+dy(sy)+dz(syz)=0
w: dx(sxz)+dy(syz)+dz(sz)=0
EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz
BOUNDARIES
surface 'bottom'
load(u)=0
load(v)=0
load(w)=0
surface 'top'
load(u)=0
load(v)=0
load(w)=-p/Ly
REGION 1
START(0,0) !y=0 surface:
load(u)=0
load(v)=0
load(w)=0!!-tauzy
LINE TO (Lx,0) !x=Lx surface
load(u)=0!40e3*(z/Lz-.5)
load(v)=0
load(w)=-33333.333 !tauzx
LINE TO (Lx,Ly) !y=Ly surface
load(u)=0
load(v)=0
load(w)=0!tauzy
LINE TO (0,Ly) !x=0 surface
value(u)=0
value(v)=0
value(w)=0!-tauzx
LINE TO CLOSE
PLOTS
grid(x+u*mag, y+v*mag, z+w*mag)
elevation(w) from (0,Ly/2,Lz/2) to (Lx,0,0)
summary
!report val(u,Lx,Ly/2,Lz/2)
!report val(v,Lx,Ly/2,Lz/2)
report val(w,Lx,Ly/2,Lz/2)
end