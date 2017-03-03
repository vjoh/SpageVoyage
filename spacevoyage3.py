from __future__ import division
from visual import *
from visual.graph import *
scene.y=400
scene.width=1024
scene.height=760


# CONSTANTS
G = 6.7e-11
mEarth = 6e24
mcraft = 15e3
mMoon = 0.000000000000001
deltat = 60
t = 0
pscale = 0.5
fscale = 2000
dpscale = 40
netscale = 1000


#OBJECTS AND INITIAL VALUES
Earth = sphere(pos=vector(0,0,0), radius=6.4e6, color=color.cyan)
scene.range=75*Earth.radius
# Add a radius for the spacecraft. It should be BIG, so it can be seen.
Moon = sphere(pos=vector(4e8,0,0), radius=1.75e6, color=color.white)
craft = sphere(pos=vector(-65920000,1536000,0), radius=1.5e6, color=color.yellow)
vcraft = vector(508,2847,0)
pcraft = mcraft*vcraft

trail = curve(color=craft.color)    # This creates a trail for the spacecraft
trail_Moon = curve(color=Moon.color)
scene.autoscale = 0                 # And this prevents zooming in or out

pArrow = arrow(color=color.green)
fArrow = arrow(color=color.cyan)
dpArrow = arrow(color=color.red)
Fnet_tangent_arrow = arrow(color=color.yellow)
Fnet_perp_arrow = arrow(color=color.magenta)


print("p=", pcraft)

U_graph = gcurve(color=color.blue)
K_graph = gcurve(color=color.yellow)
Energy_graph = gcurve(color=color.green)

Moon.v = sqrt((G * mEarth) / mag(Moon.pos - Earth.pos)) * vector(0,1,0)
momentum_Moon = mMoon * Moon.v

# CALCULATIONS
# old time: 194760
while t < 500000:
    rate(2000)   # This slows down the animation (runs faster with bigger number)

    # Add statements here for the iterative update of gravitational
    # force, momentum, and position.
    rE = craft.pos - Earth.pos
    rmagEarth = sqrt(rE.x**2 + rE.y**2 + rE.z**2)
    rhatEarth = rE/rmagEarth

    rM = craft.pos - Moon.pos
    rmagMoon = sqrt(rM.x**2 + rM.y**2 + rM.z**2)
    rhatMoon = rM/rmagMoon

    
    FmagEarth = (-G*mEarth*mcraft)/(rmagEarth**2)
    FnetE = FmagEarth*rhatEarth

    FmagMoon = (-G*mMoon*mcraft)/(rmagMoon**2)
    FnetM = FmagMoon*rhatMoon

    Fnet = FnetE + FnetM
    
    pcraft_i = pcraft + vector(0,0,0)

    p_init = mag(pcraft)
    pcraft = pcraft + Fnet*deltat
    p_final = mag(pcraft)
    Fnet_tangent = ((p_final - p_init)/deltat)*norm(pcraft)
    Fnet_perp = Fnet - Fnet_tangent
    
    #Moon Stuff
    r_EarthMoon = Moon.pos - Earth.pos
    r_craftMoon = Moon.pos - craft.pos
    Force_EarthMoon = ((-G * mEarth * mMoon) / (mag(r_EarthMoon) ** 2)) * norm(r_EarthMoon)
    Force_craftMoon = ((-G * mcraft *mMoon) / (mag(r_craftMoon) ** 2)) * norm(r_craftMoon)
    Fnet_Moon = Force_EarthMoon + Force_craftMoon
    momentum_Moon = momentum_Moon + Fnet_Moon * deltat
    Moon.pos = Moon.pos + (momentum_Moon/mMoon)*deltat

    craft.pos = craft.pos + (pcraft*deltat)/mcraft
    Fnet_tangent_arrow.pos = craft.pos
    Fnet_tangent_arrow.axis = Fnet_tangent*netscale
    Fnet_perp_arrow.pos = craft.pos
    Fnet_perp_arrow.axis = Fnet_perp*netscale

    pArrow.pos = craft.pos
    pArrow.axis = pcraft*pscale
    fArrow.pos = craft.pos
    fArrow.axis = Fnet*fscale
    deltap = pcraft - pcraft_i
    dpArrow.pos = craft.pos
    dpArrow.axis = deltap*dpscale

    # Uncomment these two lines to exit the loop if
    # the spacecraft crashes onto the Earth.
    if rmagEarth < Earth.radius: 
        break
    if rmagMoon < Moon.radius:
        break

    trail.append(pos=craft.pos)  
    trail_Moon.append(pos=Moon.pos)
    t = t+deltat

    #print("Mag SpaceCraft momentum: = ", mag(pcraft))
    #print("Mag SpaceCraft velocity: = ", mag(pcraft/mcraft))
    #print("Mag Perpendicular component: = ", mag(Fnet_perp))
    #print("Mag Seperation vector: = ", (mag(pcraft)*mag(pcraft/mcraft))/mag(Fnet_perp))
    K_craft = 0.5 * mcraft * ((mag(pcraft/mcraft))**2)
    U_craft_Earth = (-G * mcraft * mEarth) / rmagEarth
    U_craft_Moon = (-G * mcraft * mMoon) / rmagMoon
    E = K_craft + U_craft_Earth + U_craft_Moon
    U_graph.plot(pos=(t,U_craft_Earth+U_craft_Moon))
    K_graph.plot(pos=(t,K_craft))
    Energy_graph.plot(pos=(t,E))
    
     
print("Calculations finished after ",t, "seconds")
print("Fnet =", Fnet)
print("Craft Position = ", craft.pos)
print("Craft Velocity = ", pcraft/mcraft) 
print("parallel ", mag(Fnet_tangent))
print("perpendicular ", mag(Fnet_perp))
