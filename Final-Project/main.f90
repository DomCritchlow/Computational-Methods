Program soccer
implicit none
double precision::u, Fu,boundx , boundy,Drag, Magnus, rho, DragC, MagnusC, m, A, vMag, wMag,wxv_x, wxv_y, wxv_z, h, g, ax,ay, az
double precision:: vMag2, wMag2,wxv_x2, wxv_y2, wxv_z2, ax2,ay2, az2,eps,tempx,tempy
integer :: i, n,b 
double precision, allocatable,dimension(:):: Vx,Vy,Vz,rx,ry,rz,wx,wy,wz
double precision, allocatable,dimension(:):: Vx2,Vy2,Vz2,rx2,ry2,rz2,wx2,wy2,wz2
n = 5000
rho=1.2041 
allocate(Vx(n), Vy(n), Vz(n), rx(n), ry(n), rz(n), wx(n),wy(n),wz(n)) 
allocate(Vx2(n), Vy2(n), Vz2(n), rx2(n), ry2(n), rz2(n), wx2(n),wy2(n),wz2(n)) 
DragC=.250d0 
MagnusC=17.0
m= .165
A= (.00615**2)*3.1415
h=.001
g=9.81
u=.0015
boundx = 1.35d0
boundy = 2.7d0

Vx(1)=1.94d0 
Vy(1)=1.94d0
Vz(1)=0   ! 7.2
rx(1)=.675
ry(1)=.635
rz(1)=0
wx(1)=0
wy(1)=0
wz(1)=0.0d0

Vx2(1)=0.597d0
Vy2(1)=.597d0
Vz2(1)=0   ! 7.2
rx2(1)=.675
ry2(1)=1.35
rz2(1)=0
wx2(1)=0
wy2(1)=0
wz2(1)=0

eps = 8.0E-3
Drag = ((rho*A)/2*m)*DragC
Magnus = ((rho*A)/(2*m))*MagnusC

b =0
do i=1,n-1
if (((abs(rx(i)/rx2(i))-abs(rx(i))).lt.eps).and.((abs(ry(i)/ry2(i))-abs(rx(i))).lt.eps).and.b ==0)then 
b=1
print *, vx(i),vy(i),vx2(i),vy2(i)
tempx = (2*m*vx2(i))/(2*m)
tempy = (2*m*vy2(i))/(2*m)

vx2(i) = (2*m*vx(i))/(2*m)
vy2(i) = (2*m*vy(i))/(2*m)
vx(i) = tempx
vy(i) = tempy

print *,"collision"
end if 

vMag = sqrt((Vx(i)**2)+(Vy(i)**2)+(Vz(i)**2))
vMag2 = sqrt((Vx2(i)**2)+(Vy2(i)**2)+(Vz2(i)**2))
!! wMag = max(1E-16, sqrt((wx(1)**2)+(wy(1)**2)+(wz(1)**2)))
wxv_x = (wy(1)*vz(i)-wz(1)*vy(i))
wxv_y = -(wx(1)*vz(i)-wz(1)*vx(i))
wxv_z = (wx(1)*vy(i)-wy(1)*vx(i))

wxv_x2 = (wy2(1)*vz2(i)-wz2(1)*vy2(i))
wxv_y2 = -(wx2(1)*vz2(i)-wz2(1)*vx2(i))
wxv_z2 = (wx2(1)*vy2(i)-wy2(1)*vx2(i))

if ((rx(i) > boundx).or.(rx(i) <0.0d0))then
        vx(i) = vx(i)*(-1)
        
end if
if ((ry(i) > boundy).or.(ry(i) <0.0d0))then
        vy(i) = vy(i)*(-1)
end if

if ((rx2(i) > boundx).or.(rx2(i) <0.0d0))then
        vx2(i) = vx2(i)*(-1)
end if
if ((ry2(i) > boundy).or.(ry2(i) <0.0d0))then
        vy2(i) = vy2(i)*(-1)
end if


ax = -(Drag*VMag*Vx(i))+ ((Magnus*vMag)*(wxv_x))-m*g*u*vx(i)/.00615
ay = -(Drag*VMag*Vy(i))+ ((Magnus*vMag)*(wxv_y))-m*g*u*vy(i)/.00615
az = -(Drag*VMag*Vz(i))- g + ((Magnus*vMag)*(wxv_z))

ax2 = -(Drag*VMag2*Vx2(i))+ ((Magnus*vMag2)*(wxv_x2))-m*g*u*vx2(i)/.00615
ay2 = -(Drag*VMag2*Vy2(i))+ ((Magnus*vMag2)*(wxv_y2))-m*g*u*vy2(i)/.00615
az2 = -(Drag*VMag2*Vz2(i))- g + ((Magnus*vMag2)*(wxv_z2))

vx(i+1) = vx(i) + ax*h
vy(i+1) = vy(i) + ay*h
vz(i+1) = vz(i) + az*h

vx2(i+1) = vx2(i) + ax2*h
vy2(i+1) = vy2(i) + ay2*h
vz2(i+1) = vz2(i) + az2*h

rx(i+1) = rx(i)+ h*(vx(i+1)+vx(i))/2.0
ry(i+1) = ry(i)+ h*(vy(i+1)+vy(i))/2.0
rz(i+1) = rz(i)+ h*(vz(i+1)+vz(i))/2.0
rz = 0 

rx2(i+1) = rx2(i)+ h*(vx2(i+1)+vx2(i))/2.0
ry2(i+1) = ry2(i)+ h*(vy2(i+1)+vy2(i))/2.0
rz2(i+1) = rz2(i)+ h*(vz2(i+1)+vz2(i))/2.0
rz2 = 0 

end do 


open(unit=12, file= 'out1.csv')
do i=1,n-1
write(12,*) rx(i)," , ",ry(i)," , ", rz(i)

end do 
close(unit =12)

open(unit=13, file= 'out2.csv')
do i=1,n-1
write(13,*) rx2(i)," , ",ry2(i)," , ", rz2(i)

end do 
close(unit =13)





open(unit=14, file= 'out1_v.csv')
do i=1,n-1
write(14,*) vx(i)," , ",vy(i)," , ", vz(i)

end do
close(unit =14)

open(unit=15, file= 'out2_v.csv')
do i=1,n-1
write(15,*) vx2(i)," , ",vy2(i)," , ", vz2(i)

end do
close(unit =15)





deallocate(Vx,Vy,Vz,rx,ry,rz,wx,wy,wz)
deallocate(Vx2,Vy2,Vz2,rx2,ry2,rz2,wx2,wy2,wz2)

end Program soccer
