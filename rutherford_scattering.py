import numpy as np
from numpy import pi
from matplotlib import pyplot as plt
import random

class Particle():

    def __init__(self, mass, location, charge, velocity):
        self.mass = mass
        self.location = location
        self.charge = charge
        self.velocity = velocity

    #calculates distance from alpha particle to target nuclei (origin)
    def calculate_distance(self):
        return np.linalg.norm(self.location)

    def coloumb_force(self, target_particle_charge, r):

        #coulombs constant
        k = 8.99e9

        #calculating force in Newtons
        force = self.charge*target_particle_charge*k/r**2
        return force
#time step in seconds
TIME_STEP = 1e-23

#b = 13e-15 # initialized later

#mass in kg
alpha_particle_mass = 6.6446573357e-27
            
#charge in coulombs
alpha_particle_charge = 2 * 1.6e-19

#no need
target_particle_mass = None
#197 * 1.6605390666e-27

target_particle_charge = 79 * 1.6e-19

#target nuclei at origin
target_location = np.array([0,0])

#it's at rest
target_velocity = np.array([0,0])

#intilializing the target
target = Particle(target_particle_mass, target_location, target_particle_charge, target_velocity)

def simulation_loop(alpha): # must pass the object now
  
  #stores total history of alpha particles movement
  history = []

  #stores closest distance reached by alpha particle to target nuclei
  closest_distance = 1

  #stores the x,y coordinates of that closest point
  closest_location = np.array([])



  while (alpha.location[0] >= -2e-14 and alpha.location[0] < 2e-14):
      
      #adds location of alpha particle to history
      history.append(alpha.location)
      
      distance = alpha.calculate_distance()
      if distance < closest_distance:
          closest_distance = distance
          closest_location = alpha.location
      


      #calculating psialpha particle - target nuclei - X axis
      #adding a '+' to x coordinate because it is negative and it must be positive to calculate angle
      angle_psi = np.arctan(alpha.location[1]/-alpha.location[0])


      #calculating force between two particle in Newtons
      force = alpha.coloumb_force(target.charge, distance)


      #calculating the velocity acting on alpha particle due to force - this is the hypoenuse
      acting_velocity = TIME_STEP * (force / alpha.mass)



      #calculating x,y components of velocity
      velocity_x = acting_velocity * np.cos(angle_psi)


      #cheking if alpha particle is behind or infront of target nuclei, because that affects if acting x veloxity is '+' or '-'
      if alpha.location[0] < 0:
          velocity_x = velocity_x * -1

      velocity_y = acting_velocity * np.sin(angle_psi)
      acting_velocity = np.array([velocity_x, velocity_y])

      #vector addition with acting velocity and initial velocity
      alpha.velocity = alpha.velocity + acting_velocity

      #calculating the movement of the alpha particle within the timestep for x,y components
      movement_x = alpha.velocity[0] * TIME_STEP
      movement_y = alpha.velocity[1] * TIME_STEP
      movement = np.array([movement_x, movement_y])

      #moving the alpha particle
      alpha.location = alpha.location + movement

      #moving the velocity the same amount so relative origin stays same for next step
      alpha.velocity = alpha.velocity + movement

  #print(closest_location, alpha.location, history)
  return history, closest_location

def calculating_scattering_angle(closest_approach,alpha, b):

  #we split the closest approach into x and y components
  closest_y = closest_approach[1]

  #if the x component is negative we make it positive
  if closest_approach[0] < 0:
    closest_x = -1 * closest_approach[0]
  else:
    closest_x = closest_approach[0]

  #we split the final location of the alpha particle into x and y components
  end_loc_y = alpha.location[1]

  #if the x component is negative we make it positive
  if alpha.location[0] < 0:
    end_loc_x = -1 * alpha.location[0]
  else:
    end_loc_x = alpha.location[0]



  #calculates angle psi of the closest approach of alpha particle - returns radians
  psi_of_closest_approach = np.arctan(closest_y/closest_x)


  #finding the point where the line between closest approach and target nuclei intercepts with b
  intersection_y = b
  #np.tan() takes in rad
  intersection_x = intersection_y/np.tan(psi_of_closest_approach)
  intersection_location = np.array([intersection_x, intersection_y])

  #since we can't divide by 0, if it bounces straight back we hardcode the scattering angle to 0
  if end_loc_y == 0:
    scattering_angle = 180
  else:
    #finds difference between intersection point and final alpha location
    difference_x = end_loc_x - intersection_location[0]
    difference_y = end_loc_y - intersection_location[1]

    #calculates scattering angle
    scattering_angle = np.arctan(difference_y/difference_x)

    #converts to degrees from radians
    scattering_angle = np.rad2deg(scattering_angle)

    #if the alpha particle bounces back it returns 180 - scattering angle 
    if alpha.location[0] < closest_approach[0]:
      scattering_angle = 180 - scattering_angle 
  
  #for case when end_loc x is negative as well as closest approach x but end_lox x is bigger 
  if scattering_angle < 0:
    scattering_angle = scattering_angle *-1

  return scattering_angle 

def plot_prepare_trajectory(history):
    x_loc = []
    y_loc = []

    for i in history:
      x_loc.append(i[0])
      y_loc.append(i[1])
    
    return x_loc, y_loc

def simulater(alpha, b): #must pass the model now
  history, closest_approach = simulation_loop(alpha)
  scattering_angle = calculating_scattering_angle(closest_approach, alpha, b)
  return scattering_angle, history, closest_approach

#converts MeV to Joules
def mev_to_vel(mev, mass):
  joule = mev * 1.6022E-13
  velocity = np.sqrt(joule/(0.5*mass))
  return velocity 

def myround(x, base=15):
    return base * round(x/base)

def main_trajectory(number_of_alpha_particles=10, b_max=1e-14, target_proton_num=79, mev=7.7):

  start_velocity = mev_to_vel(mev, alpha_particle_mass)

  target.charge = target_proton_num * 1.6e-19


  plt.figure(figsize=(10,10))
  plt.style.use('ggplot')

  #plotting target nucleus
  plt.plot(0,0, marker = 'o')

  for i in range(number_of_alpha_particles):

    #selecting y distance (impact parameter) range to randomly initialize particles from
    b = np.random.uniform(0,5)*(b_max)

    #start location in meters
    alpha_start_location = np.array([-2e-14,b])


    alpha_start_velocity = np.array([start_velocity- alpha_start_location[0], b])

    alpha = Particle(alpha_particle_mass, alpha_start_location, alpha_particle_charge, alpha_start_velocity)

    #first alpha particle
    scattering_angle, history, closest_approach = simulater(alpha, b) # passing the first particle

    x_loc, y_loc = plot_prepare_trajectory(history)


    plt.plot(x_loc, y_loc, '--',label=f'{i+1}')



  plt.legend()

  plt.title('Alpha Particle Trajectory')
  plt.ylabel('Y distance - impact parameter - (m e-14)')
  plt.xlabel('X distance (m e-14)')

  plt.tight_layout()

  ax = plt.gca()
  ax.set_ylim([-5e-14,5e-14])
  ax.set_xlim([-2.1e-14, 2e-14])
    
  plt.grid(True)
  plt.show()

def main_r_formula(number_of_alpha_particles=1000,  b_max=1e-13, target_proton_num=79, mev=7.7, label = 'Data points'):
  global alpha_particle_mass

  #calcualtes the start x velocity of alpha particle
  start_velocity = mev_to_vel(mev, alpha_particle_mass)

  #calculates the charge of the target in joules
  target.charge = target_proton_num * 1.6e-19


  plt.figure(figsize=(10,10))
  plt.style.use('ggplot')

  #counts for the scattering angle rounded to nearest 15 degrees
  scattering_angle_count = {0: 0, 15: 0, 30: 0, 45: 0, 60: 0, 75: 0, 90: 0, 105: 0, 120: 0, 135: 0, 150: 0, 165: 0, 180: 0}

  #counts for the scattering angle rounded to nearest 60 degrees - this is used later for quantitative analysis
  scattering_angle_count_grouped = {0: 0,60: 0,120: 0, 180: 0}


  for i in range(number_of_alpha_particles):
    b = np.random.uniform(0,1)*(b_max)
    #start location in meters
    alpha_start_location = np.array([-2e-14,b])

    #velocity at the start of the experiment
    alpha_start_velocity = np.array([start_velocity - alpha_start_location[0], b])

    #initializing the particle
    alpha = Particle(alpha_particle_mass, alpha_start_location, alpha_particle_charge, alpha_start_velocity)

    #first alpha particle
    scattering_angle, history, closest_approach = simulater(alpha, b) # passing the first particle

    #rounding to nearest 15 degrees
    rounded_scattering_angle = myround(scattering_angle)

    higher_rounded_scattering_angle = myround(scattering_angle, base=60)

    #adding to count
    scattering_angle_count[rounded_scattering_angle] = scattering_angle_count[rounded_scattering_angle] + 1

    #adding to count
    scattering_angle_count_grouped[higher_rounded_scattering_angle] = scattering_angle_count_grouped[higher_rounded_scattering_angle] + 1

    #shows the user how many particles have been sent
    print(f'{i+1} particles fired')

  




  #gets values from scattering angle dictionary
  angles = list(scattering_angle_count.keys())
  count = list(scattering_angle_count.values())

  #gets count from dictionary rounded to nearest 60 degrees
  rounded_count = list(scattering_angle_count_grouped.values())

  #plots number of particles against scattering angle
  plt.plot(angles, count, marker = 'o', label=label)

  plt.title("Number of alpha particles scattered against scattering angle (Rutherford's Formula)")
  plt.ylabel('Number of alpha particles')
  plt.xlabel('Scattering angles (rounded to nearest 15)')

  plt.legend()

  plt.tight_layout()
    
  plt.grid(True)
  plt.show()

  return rounded_count

#compares rutherford's formula for alpha particles with different energies
def compare(number_of_alpha_particles=1000,  b_max=1e-13, target_proton_num=79, mev=7.7, label = '1', number_of_alpha_particles1=1000,  b_max1=1e-13, target_proton_num1=79, mev1=7.7, label1 ='2', title='Comparison'):
  angles = [0,15,30,45,60,75,90,105,120,135,150,165,180]

  one = main_r_formula(number_of_alpha_particles,  b_max, target_proton_num, mev)

  two = main_r_formula(number_of_alpha_particles1,  b_max1, target_proton_num1, mev1)

  plt.plot(angles, one, marker = 'o', label = label)
  plt.plot(angles, two, marker = 'o', label = label1)
  plt.title(title)
  plt.ylabel('Number of alpha particles')
  plt.xlabel('Scattering angles (rounded to nearest 15)')

  plt.tight_layout()
    
  plt.grid(True)
  plt.legend()
  plt.show()

def test(number_of_alpha_particles=1000, target_proton_num=79, mev=7.7):
  #atoms per unit volume gold
  n = 1

  #thichness of target (we only have one nucleus) in meters 
  L = 1.2e-15 * 5.8

  #electron charge
  e = -1.60217663e-19

  #coulombs constant
  k = 8.99e9

  #atomic number of target 
  Z = target_proton_num

  #kinetic energy of particle in mev
  KE = mev

  #total fired
  Ni = number_of_alpha_particles

  #distance to detector
  r = 2e-14

  scattering_counts = []

  #goes through angles 0-180 in steps of 60
  for i in range(0,240,60):
    
    #converts the degrees into radians
    scattering_angle = np.deg2rad(i)

    #since the formula won't work with 0, and the angle never is truely 0
    if i==0:
      N = (Ni * n * L * Z**2 * k**2 * e**4) / (4 * r**2 * KE**2 * np.sin(0.5/2)**4)

    #we plug in the values here: all SI units
    else:
      N = (Ni * n * L * Z**2 * k**2 * e**4) / (4 * r**2 * KE**2 * np.sin(scattering_angle/2)**4)

    scattering_counts.append([i, N])

  total = 0

  #adds solutions to total
  for i in scattering_counts:
   total+=i[1]

  #creates a factor so that the total = Ni
  factor = total/Ni

  #uses factor and we get new results
  for i in scattering_counts:
     i[1] = i[1]/factor
    
    #rounded to nearest whole number
     i[1] = round(i[1])

  counts = []
  #we add the factored counts to the counts array
  for i in scattering_counts:
    counts.append(i[1])
  
  angles = [0,60,120,180]

  simulation_count = main_r_formula(number_of_alpha_particles=1000,  b_max=3.1e-13, target_proton_num=79, mev=7.7, label = 'Simulation results')

  #width of bar
  w = 0.4

  #this part is done so we can have two bars side by side.
  bar1 = np.arange(len(angles))
  bar2 = [i+w for i in bar1]
  plt.bar(bar1, counts, width = w, label = "Rutherford's formula",align="center")
  plt.bar(bar2, simulation_count, width = w, label = "Simulation results",align="center")
  plt.xticks(bar1 + w/2, angles)

  plt.legend()
  plt.xlabel('Angles (degrees)')
  plt.ylabel('Number of particles')

  print(f'Expected: {counts}')
  print(f'Observed: {simulation_count}')

  plt.title('Comparison between expected results and observed results')
  plt.legend()
  plt.tight_layout()
  plt.grid(True)
  plt.show()