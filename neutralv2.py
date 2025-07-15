import shapely
from shapely.geometry import Polygon, Point, LineString
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


def circular_arc(center, radius, start_angle, end_angle, num_points):
    """Generate points along a circular arc."""
    angles = np.linspace(start_angle, end_angle, num_points)
    return [
        (
            center[0] + radius * np.cos(angle),
            center[1] + radius * np.sin(angle)
        )
        for angle in angles
    ]

def cubiroot(first,last,offset,amp,num_points):
    points=np.linspace(first[1],last[1],num_points) 
    return [
        (
            amp*(point-offset[1])**3+offset[0],point
        )
        for point in points
    ]

inside1 = [(-12,12.1),(-13,12.1),(-14,12.1),(-15,12.1),(-16,12.1),(-17,12.1),(-18,12.1)]#First straight
outer1 = [(-12,14.9),(-13,14.9),(-14,14.9),(-15,14.9),(-16,14.9),(-17,14.9),(-18,14.9)]#First straight

inside2=circular_arc((-18,6.5), 5.6, np.radians(90), np.radians(180), 15)#First bend
outer2=circular_arc((-18,6.5), 5.6+2.8, np.radians(90), np.radians(180), 15)#First bend

inside3=[(-23.6,5.5),(-23.6,4.5),(-23.6,3.5),(-23.6,2.5),(-23.6,1.5),(-23.6,0.5),(-23.6,-0.5),(-23.6,-1.5),(-23.6,-2.5),(-23.6,-3.5)]#Second straight
outer3=[(-26.4,5.5),(-26.4,4.5),(-26.4,3.5),(-26.4,2.5),(-26.4,1.5),(-26.4,0.5),(-26.4,-0.5),(-26.4,-1.5),(-26.4,-2.5),(-26.4,-3.5)]#Second straight

inside4=circular_arc((-18,-3.5), 5.6, np.radians(180), np.radians(270), 15)#Second bend
outer4=circular_arc((-18,-3.5), 5.6+2.8, np.radians(180), np.radians(270), 15)#Second bend

inside5=[(-18,-9.1),(5,-9.1)]
outer5=[(-18,-11.9),(6,-11.9)]

inside6=cubiroot((6,-8.5),(27,5.5),(22,-0.5),0.03,20)
outer6=cubiroot((8,-10.9),(27.9,2.7),(25,-3),0.03,20)

#From top left of the tunnel around in a counter clockwise pattern
outer = outer1 + outer2 + outer3 + outer4 + outer5 +outer6
inside = inside1 + inside2 + inside3 + inside4 + inside5 +inside6
insideline = LineString(inside)
c_detector = Polygon(outer + inside[::-1])

#test = Polygon([(0,0),(30*np.cos(np.radians(-135)),30*np.sin(np.radians(-135))),(30*np.cos(np.radians(-90)),30*np.sin(np.radians(-90)))])

#c_detector_sect = shapely.intersection(c_detector,test) 


origin = Point(0,0)
#dist = shapely.distance(origin, c_detector_sect)
#print(dist)


#plotting
#fig, ax = plt.subplots(figsize=(10, 10))
#x, y = c_detector.exterior.xy
#ax.plot(*origin.xy, 'o', color="Red")
#ax.fill(x, y, facecolor='green', edgecolor='black', linewidth=2, label="Detector")
#ax.set_aspect('equal')
#ax.set_xlim(-40, 40)
#ax.set_ylim(-40, 40)
#plt.legend(loc="upper right")
#plt.title("Detector")
#plt.grid(True)
#plt.show()

#slicing the geometry 
#getting the L1 distance per angle section
phi = np.linspace(-235, 15, 100)
deltaphi = np.array([])
for j in range(len(phi)-1):
    #phi[j] = np.radians(phi[j])
    diff = phi[j] - phi[j+1]
    diff= abs(phi)
    deltaphi = np.append(deltaphi,diff)
xarr=np.array([])
yarr=np.array([])
L1=np.array([])


for i in range(len(phi)):
    r=30
    x=r*np.cos(np.radians(phi[i]))
    y=r*np.sin(np.radians(phi[i]))
   
    xarr=np.append(xarr,x)
    yarr=np.append(yarr,y)
    if i!=0 :
        temp= Polygon([(0,0),(xarr[i-1],yarr[i-1]),(xarr[i],yarr[i])])
        slice = shapely.intersection(c_detector,temp) 
        dist = shapely.distance(origin, slice)
        L1= np.append(L1,dist)
      

Estimate=np.array([])
d=15

for k in range(len(phi)-1):
    Probestimate = (1/(4*np.pi))*deltaphi[k]*np.exp(-L1[k]/d)*2.8/d
    Estimate=np.append(Estimate,Probestimate)
#print(Estimate)
plt.figure(figsize=(10, 4))
plt.plot(phi[:-1],Estimate, color='skyblue')
plt.xlabel('Phi (degrees)')
plt.ylabel('Probability')
plt.title('Estimated probability d=15')
plt.grid(True)
plt.tight_layout()
plt.show()

def solidangle(a):
    return a
def depth(L,d):
    return (1/d)*np.exp(-L/d)

Exact=np.array([])
for k in range(len(phi)-1):
    a=deltaphi[k]
    angle =quad(solidangle,phi[k],phi[k+1])
    d=15
    detector =quad(depth,L1[k],L1[k]+2.8,args=(d))
    Prob = (1/(4*np.pi))*-1*angle[0]*detector[0]
    Exact=np.append(Exact,Prob)
#print(Estimate)
plt.figure(figsize=(10, 4))
plt.plot(phi[:-1],Exact, color='Red')
plt.xlabel('Phi (degrees)')
plt.ylabel('Probability')
plt.title('Probability d=15')
plt.grid(True)
plt.tight_layout()
plt.show()